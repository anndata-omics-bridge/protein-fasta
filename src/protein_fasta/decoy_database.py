"""Generate one search database from an existing biological inventory."""

from __future__ import annotations

import importlib.metadata
import os
from dataclasses import dataclass, replace
from pathlib import Path
from typing import cast

import polars as pl
from loguru import logger

from protein_fasta.analytics.hashing import sequence_hash
from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.database.models import (
    DecoyInventoryEntry,
    ProteinInventoryEntry,
    SearchDatabase,
)
from protein_fasta.decoy_compile import make_decoy_generation
from protein_fasta.inventory import (
    protein_inventory_entries,
    read_protein_inventory,
    search_database_frame,
)
from protein_fasta.reading.parser import FastaRecord
from protein_fasta.reading.writer import write_records
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.schema.decoy import (
    DecoyCountsDocument,
    DecoyGenerationEvidenceDocument,
    DecoyRequestDocument,
    DecoyResultDocument,
    DecoySummaryDocument,
    EffectiveDecoyRequestDocument,
)
from protein_fasta.summary import FastaSummary, summarize_sequences

_DECOY_SOURCE_KINDS = {
    EntryKind.TARGET.value,
    EntryKind.CONTAMINANT.value,
    EntryKind.ENTRAPMENT.value,
}


@dataclass(frozen=True, slots=True)
class DecoyOverrides:
    """Common explicitly supplied values taking precedence over one request."""

    output_fasta: Path | None = None
    decoy_prefix: str | None = None


DEFAULT_DECOY_OVERRIDES = DecoyOverrides()


@dataclass(frozen=True, slots=True)
class DecoyExecution:
    """Search artifacts and evidence from one completed decoy operation."""

    document: DecoyResultDocument
    database: SearchDatabase
    search_fasta_path: Path
    search_inventory_path: Path
    effective_request_path: Path
    result_path: Path


def resolve_decoy_request(
    request: DecoyRequestDocument,
    /,
    *,
    request_base: Path,
    overrides: DecoyOverrides = DEFAULT_DECOY_OVERRIDES,
) -> EffectiveDecoyRequestDocument:
    """Resolve request-relative output and explicitly supplied common overrides."""
    output = overrides.output_fasta or request.output_fasta
    if not output.is_absolute():
        output = request_base / output
    return EffectiveDecoyRequestDocument(
        output_fasta=output.resolve(),
        decoy_prefix=overrides.decoy_prefix or request.decoy_prefix,
        strategy=request.strategy,
    )


def run_decoy_generation(
    biological_inventory_path: Path,
    effective: EffectiveDecoyRequestDocument,
    /,
) -> DecoyExecution:
    """Generate and atomically publish search FASTA, inventory, and evidence."""
    biological_path = biological_inventory_path.resolve()
    biological = read_protein_inventory(biological_path)
    _reject_decoy_rows(biological, biological_path)
    biological_entries = protein_inventory_entries(biological)
    output = effective.output_fasta
    search_inventory_path = output.with_suffix(f"{output.suffix}.search-inventory.parquet")
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, search_inventory_path, result_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(effective_path, effective.model_dump(mode="json"), replace_existing=True)

    generation = make_decoy_generation(effective.strategy)
    source_entries_runtime = [
        entry for entry in biological_entries if entry.kind in _DECOY_SOURCE_KINDS
    ]
    source_entries = tuple((entry.raw_header, entry.sequence) for entry in source_entries_runtime)
    batch = generation.generate(source_entries, prefix=effective.decoy_prefix)
    source_by_header = {entry.raw_header: entry for entry in source_entries_runtime}
    search_biological_entries = list(biological_entries)
    if search_biological_entries:
        first = search_biological_entries[0]
        search_biological_entries[0] = replace(
            first,
            raw_header=_search_sentinel(
                first.raw_header,
                generation.annotation(
                    initial_collisions=batch.initial_collisions,
                    dropped_peptides=batch.dropped_peptides,
                ),
                generation.mode.value,
            ),
        )
    decoy_entries = [
        _decoy_row(
            len(search_biological_entries) + index,
            header,
            sequence,
            source_by_header,
            effective.decoy_prefix,
            generation.mode.value,
        )
        for index, (header, sequence) in enumerate(batch.entries)
    ]
    database = SearchDatabase((*search_biological_entries, *decoy_entries))
    search = search_database_frame(database)
    summary = summarize_sequences(
        entry.sequence
        for entry in database.entries
        if entry.kind in {*_DECOY_SOURCE_KINDS, EntryKind.DECOY.value}
    )

    with temporary_sibling(output) as staged_fasta:
        write_records(
            (FastaRecord(entry.raw_header, entry.sequence) for entry in database.entries),
            staged_fasta,
        )
        with temporary_sibling(search_inventory_path) as staged_inventory:
            search.write_parquet(staged_inventory)
            document = DecoyResultDocument(
                protein_fasta_version=importlib.metadata.version("protein-fasta"),
                effective_request=effective,
                biological_inventory=artifact_document(
                    biological_path,
                    recorded_path=Path(os.path.relpath(biological_path, start=output.parent)),
                    schema_name="protein-inventory",
                    schema_version="1",
                    row_count=biological.height,
                ),
                search_fasta=artifact_document(
                    staged_fasta,
                    recorded_path=Path(output.name),
                    schema_name="search-fasta",
                    schema_version="1",
                    row_count=search.height,
                ),
                search_inventory=artifact_document(
                    staged_inventory,
                    recorded_path=Path(search_inventory_path.name),
                    schema_name="search-inventory",
                    schema_version="1",
                    row_count=search.height,
                ),
                counts=DecoyCountsDocument(
                    biological=biological.height,
                    decoy=len(decoy_entries),
                    total=search.height,
                ),
                summary=_summary_document(summary),
                generation=DecoyGenerationEvidenceDocument(
                    strategy=generation.mode.value,
                    seed=generation.seed,
                    parameters=batch.parameters,
                    initial_collisions=batch.initial_collisions,
                    unresolved_collisions=batch.unresolved_collisions,
                    dropped_peptides=batch.dropped_peptides,
                    omitted_decoys=batch.omitted_decoys,
                ),
            )
            publish_exclusive(staged_fasta, output)
            try:
                publish_exclusive(staged_inventory, search_inventory_path)
                write_json_atomic(
                    result_path,
                    document.model_dump(mode="json"),
                    replace_existing=False,
                )
            except BaseException:
                output.unlink(missing_ok=True)
                search_inventory_path.unlink(missing_ok=True)
                raise
    logger.info("Generated {} decoys -> {}", len(decoy_entries), output)
    return DecoyExecution(
        document,
        database,
        output,
        search_inventory_path,
        effective_path,
        result_path,
    )


def _reject_decoy_rows(frame: pl.DataFrame, path: Path) -> None:
    kinds = cast("list[str]", frame.get_column("kind").to_list())
    if EntryKind.DECOY.value in kinds:
        raise ValueError(f"biological inventory {path} already contains decoy rows")


def _decoy_row(
    final_order: int,
    raw_header: str,
    sequence: str,
    source_by_header: dict[str, ProteinInventoryEntry],
    prefix: str,
    strategy: str,
) -> DecoyInventoryEntry:
    source_header = raw_header.removeprefix(prefix)
    source = source_by_header.get(source_header)
    if source is None:
        raise ValueError(f"generated decoy {raw_header!r} does not identify a biological source")
    source_id = source.identifier
    if source.kind == "sentinel":
        raise ValueError(f"generated decoy {raw_header!r} identifies a metadata source")
    return DecoyInventoryEntry(
        final_order=final_order,
        raw_header=raw_header,
        identifier=f"{prefix}{source_id}",
        description=source.description,
        sequence=sequence,
        contaminant_group=source.contaminant_group,
        sequence_hash=sequence_hash(sequence).hex(),
        entrapment_strategy=source.entrapment_strategy,
        decoy_strategy=strategy,
        decoy_source_order=source.final_order,
        decoy_source_id=source_id,
        decoy_source_kind=source.kind,
    )


def _search_sentinel(header: str, annotation: str, strategy: str) -> str:
    note = annotation or f"decoys {strategy}"
    return f"{header}; {note}"


def _summary_document(summary: FastaSummary) -> DecoySummaryDocument:
    frequencies = (
        {
            residue: count / summary.total_residues
            for residue, count in summary.aa_frequencies.items()
        }
        if summary.total_residues
        else {}
    )
    return DecoySummaryDocument(
        n_sequences=summary.n_sequences,
        length_min=summary.length_min,
        length_max=summary.length_max,
        length_mean=summary.length_mean,
        length_q1=summary.length_q1,
        length_median=summary.length_median,
        length_q3=summary.length_q3,
        total_residues=summary.total_residues,
        aa_counts=summary.aa_frequencies,
        aa_frequencies=frequencies,
    )


def _refuse_existing(*paths: Path) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing decoy artifacts: {names}")
