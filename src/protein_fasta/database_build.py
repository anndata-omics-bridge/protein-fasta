"""Assemble reproducible protein FASTA databases.

The assembly operation is a Python port of prozor's ``create_fgcz_fasta_db``.

Assembly order mirrors prozor:

    [db sentinel] + targets + [contaminants] + [entrapment] + decoys(all of those)

with two FGCZ-convention additions from the newer collection: an ``aa|`` database
sentinel as the first entry, and an ``aa|Cont_...`` section marker before each
contaminant set's block. Decoys are generated for targets **and** contaminants
**and** entrapment records, but NOT for the metadata entries. Entrapment records
belong to the target space, so they are generated before the decoys rather than
after: decoying only the biological entries would leave a target space the decoy
space does not cover.

Every input sequence is normalized first -- upper-cased, one trailing stop removed
-- so what is written, what is decoyed, and what deduplication compares are the same
content. Duplicate entries are then removed by id, the first occurrence winning as
in prozor's ``!duplicated(names(...))``, but only when their sequences agree: an id
carrying two different sequences is a disagreement between sources and raises
rather than silently losing one of them.
"""

from __future__ import annotations

import datetime
import importlib.metadata
import json
import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass, replace
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

from loguru import logger

from protein_fasta.analytics.hashing import FILE_CHECKSUM_VERSION, file_checksum, sequence_hash
from protein_fasta.build.generation.decoy import (
    DEFAULT_DECOY_SPEC,
    make_decoy_generation,
)
from protein_fasta.build.generation.decoy_types import DecoyBatch, DecoyGeneration
from protein_fasta.build.metadata import build_section_marker_header, build_sentinel_header
from protein_fasta.build.naming import build_dbname, build_fasta_name
from protein_fasta.diagnostics.messages import describe_illegal_residues
from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import FastaRecord
from protein_fasta.registry.classification import ContaminantBlockState, classify_record
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.registry.rules import RegistryDiagnosticRules, load_registry_diagnostics
from protein_fasta.schema.build import (
    BuildArtifactDocument,
    DatabaseBuildCountsDocument,
    DatabaseBuildDecoyEvidenceDocument,
    DatabaseBuildDocument,
    DatabaseBuildEntrapmentEvidenceDocument,
    DatabaseBuildNormalizationDocument,
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
    DatabaseBuildResultDocument,
    DatabaseBuildSummaryDocument,
    DecoyDocument,
    DecoyMode,
    EffectiveDatabaseBuildDocument,
    EntrapmentDocument,
    MetadataDocument,
    NamingDocument,
)
from protein_fasta.summary import FastaSummary, summarize_sequences
from protein_fasta.validation.sequence import normalize_sequence
from protein_fasta.writing import write_records

if TYPE_CHECKING:
    from fdr_benchmark.models import PeptidePairRecord

    from protein_fasta.build.generation.entrapment import (
        EntrapmentBatch,
        EntrapmentGeneration,
    )

Entry = tuple[str, str]

#: Conflicting ids named in the error before it summarizes the rest. Enough to act
#: on, few enough to read.
_REPORTED_CONFLICTS = 5


def _make_entrapment_generation(spec: EntrapmentDocument) -> EntrapmentGeneration:
    """Load the optional entrapment adapter only when a build selects it."""
    try:
        from protein_fasta.build.generation.entrapment import make_entrapment_generation
    except ModuleNotFoundError as error:
        missing_name = error.name or ""
        if missing_name == "fdr_benchmark" or missing_name.startswith("fdr_benchmark."):
            message = "entrapment generation requires the 'protein-fasta[generation]' extra"
            raise RuntimeError(message) from error
        raise
    return make_entrapment_generation(spec)


def _format_entrapment_pairs(pairs: tuple[PeptidePairRecord, ...]) -> str:
    """Format mappings through the package that generated them."""
    from fdr_benchmark.provenance import format_peptide_pairs

    return format_peptide_pairs(pairs)


@dataclass(frozen=True, slots=True)
class ContaminantBlock:
    """One resolved contaminant module supplied to the build pipeline."""

    name: str
    description: str
    entries: tuple[Entry, ...]


@dataclass(frozen=True)
class DecoyBuildEvidence:
    """Decoy algorithm identity, parameters, and collision outcomes."""

    mode: str
    seed: int | None
    parameters: dict[str, Any]
    initial_collisions: int
    unresolved_collisions: int
    dropped_peptides: int
    omitted_decoys: int


@dataclass(frozen=True)
class EntrapmentBuildEvidence:
    """Entrapment algorithm identity and achieved multiplicity."""

    strategy: str
    seed: int
    requested_fold: int
    achieved_fold: int
    failures: int
    proteins_affected: int
    source_proteins: int


@dataclass(frozen=True)
class PipelineResult:
    """Outcome of a database build."""

    path: Path
    dbname: str
    n_target: int
    n_contaminant: int
    n_decoy: int
    n_total: int
    contaminant_sets: list[str]
    summary: FastaSummary
    decoy: DecoyBuildEvidence | None = None
    n_entrapment: int = 0
    entrapment: EntrapmentBuildEvidence | None = None
    entrapment_pairs_path: Path | None = None
    # What normalization and deduplication changed on the way in, so a produced
    # file can account for the difference from its sources.
    upper_cased_entries: int = 0
    stop_stripped_entries: int = 0
    duplicates_dropped: int = 0


class BuildDecoyOverride(StrEnum):
    """Typed CLI/API override for the configured decoy stage."""

    NONE = "none"
    REVERSE = DecoyMode.REVERSE.value
    SHUFFLE = DecoyMode.SHUFFLE.value
    DECOYPYRAT = DecoyMode.DECOYPYRAT.value


@dataclass(frozen=True, slots=True)
class DatabaseBuildOverrides:
    """Common explicitly supplied values that take precedence over JSON."""

    date: datetime.date | None = None
    decoy: BuildDecoyOverride | None = None


DEFAULT_BUILD_OVERRIDES = DatabaseBuildOverrides()


@dataclass(frozen=True, slots=True)
class DatabaseBuildExecution:
    """Runtime result and its durable JSON evidence paths."""

    result: PipelineResult
    document: DatabaseBuildResultDocument
    effective_request_path: Path
    result_path: Path


def _decoy_build_evidence(
    generation: DecoyGeneration | None,
    batch: DecoyBatch,
) -> DecoyBuildEvidence | None:
    if generation is None:
        return None
    return DecoyBuildEvidence(
        mode=generation.mode.value,
        seed=generation.seed,
        parameters=batch.parameters,
        initial_collisions=batch.initial_collisions,
        unresolved_collisions=batch.unresolved_collisions,
        dropped_peptides=batch.dropped_peptides,
        omitted_decoys=batch.omitted_decoys,
    )


def _entrapment_build_evidence(
    generation: EntrapmentGeneration | None,
    batch: EntrapmentBatch | None,
) -> EntrapmentBuildEvidence | None:
    if generation is None or batch is None:
        return None
    return EntrapmentBuildEvidence(
        strategy=generation.strategy.value,
        seed=generation.seed,
        requested_fold=batch.requested_fold,
        achieved_fold=batch.achieved_fold,
        failures=batch.failures,
        proteins_affected=batch.proteins_affected,
        source_proteins=batch.source_proteins,
    )


def resolve_database_build(
    profile: DatabaseBuildProfileDocument,
    request: DatabaseBuildRequestDocument,
    *,
    profile_base: Path,
    request_base: Path,
    overrides: DatabaseBuildOverrides = DEFAULT_BUILD_OVERRIDES,
) -> EffectiveDatabaseBuildDocument:
    """Resolve profile, request, and explicit overrides into replayable paths."""
    naming = request.naming or profile.naming
    metadata = request.metadata or profile.metadata
    diagnostics = _resolved_profile_or_request_path(
        profile.diagnostics,
        request.diagnostics,
        request_field_is_set="diagnostics" in request.model_fields_set,
        profile_base=profile_base,
        request_base=request_base,
    )
    decoy = _resolved_decoy(profile, request, overrides.decoy)
    entrapment = (
        request.entrapment
        if "entrapment" in request.model_fields_set
        else profile.default_entrapment
    )
    return EffectiveDatabaseBuildDocument(
        targets=tuple(_resolved_path(path, request_base) for path in request.targets),
        output_dir=_resolved_path(request.output_dir, request_base),
        date=overrides.date or request.date,
        name_fields=request.name_fields,
        template=request.template or naming.default_dbname,
        naming=naming,
        metadata=metadata,
        diagnostics=diagnostics,
        contaminant_blocks=tuple(
            block.model_copy(update={"path": _resolved_path(block.path, request_base)})
            for block in request.contaminant_blocks
        ),
        decoy=decoy,
        entrapment=entrapment,
        foreign_sources=tuple(
            _resolved_path(path, request_base) for path in request.foreign_sources
        ),
        annotation=request.annotation,
        installer=request.installer,
    )


def run_database_build(effective: EffectiveDatabaseBuildDocument) -> DatabaseBuildExecution:
    """Read one effective request, build its FASTA, and persist typed evidence."""
    effective.output_dir.mkdir(parents=True, exist_ok=True)
    expected_name = build_fasta_name(
        config=effective.naming,
        template=effective.template,
        date=effective.date,
        decoy=effective.decoy is not None,
        entrapment=effective.entrapment is not None,
        **effective.name_fields,
    )
    expected_path = effective.output_dir / expected_name
    effective_path = expected_path.with_suffix(f"{expected_path.suffix}.effective.json")
    _write_json_atomic(effective_path, effective.model_dump(mode="json"))

    target_entries = _entries_from_paths(effective.targets)
    foreign_entries = _entries_from_paths(effective.foreign_sources)
    blocks = tuple(
        ContaminantBlock(
            name=block.name,
            description=block.description,
            entries=tuple(_entries_from_paths((block.path,))),
        )
        for block in effective.contaminant_blocks
    )

    diagnostics = load_registry_diagnostics(effective.diagnostics)
    result = build_database(
        targets=target_entries,
        name_fields=effective.name_fields,
        naming=effective.naming,
        metadata=effective.metadata,
        diagnostics=diagnostics,
        output_dir=effective.output_dir,
        date=effective.date,
        template=effective.template,
        contaminant_blocks=blocks,
        decoy_spec=effective.decoy,
        entrapment_spec=effective.entrapment,
        foreign_entries=foreign_entries,
        annotation=effective.annotation,
        installer=effective.installer,
    )
    input_paths = (
        *effective.targets,
        *(block.path for block in effective.contaminant_blocks),
        *effective.foreign_sources,
    )
    inventory_path = result.path.with_suffix(f"{result.path.suffix}.protein-inventory.parquet")
    inventory_rows = _write_protein_inventory(result, diagnostics, inventory_path)
    result_document = _result_document(
        effective,
        result,
        effective_path,
        input_paths,
        inventory_path,
        inventory_rows,
    )
    result_path = result.path.with_suffix(f"{result.path.suffix}.result.json")
    _write_json_atomic(result_path, result_document.model_dump(mode="json"))
    return DatabaseBuildExecution(result, result_document, effective_path, result_path)


def _write_protein_inventory(
    result: PipelineResult,
    diagnostics: RegistryDiagnosticRules,
    path: Path,
) -> int:
    """Write the final FASTA order and operational provenance as Parquet."""
    try:
        import polars as pl
    except ModuleNotFoundError as error:
        if error.name == "polars":
            message = "protein inventory output requires the 'protein-fasta[frame]' extra"
            raise RuntimeError(message) from error
        raise
    rows: list[dict[str, object]] = []
    block_state: ContaminantBlockState | None = None
    for order, record in enumerate(_read_existing_records(result.path)):
        header = parse_header(record.raw_header)
        _, classifications = diagnostics.rules.diagnose_identifier(header.id)
        kind, contaminant_group, block_state = classify_record(
            record.raw_header,
            classifications,
            block_state,
            diagnostics.decoy_prefix,
        )
        rows.append(
            {
                "final_order": order,
                "raw_header": record.raw_header,
                "id": header.id,
                "description": header.description,
                "sequence": record.sequence,
                "kind": kind.value,
                "contaminant_group": contaminant_group,
                "sequence_hash": sequence_hash(record.sequence).hex(),
                "decoy_mode": (
                    result.decoy.mode
                    if kind is EntryKind.DECOY and result.decoy is not None
                    else None
                ),
                "entrapment_strategy": (
                    result.entrapment.strategy
                    if kind is EntryKind.ENTRAPMENT and result.entrapment is not None
                    else None
                ),
            }
        )
    pl.DataFrame(rows).write_parquet(path)
    return len(rows)


def _resolved_path(path: Path, base: Path) -> Path:
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _resolved_profile_or_request_path(
    profile_path: Path | None,
    request_path: Path | None,
    *,
    request_field_is_set: bool,
    profile_base: Path,
    request_base: Path,
) -> Path | None:
    if request_field_is_set:
        return None if request_path is None else _resolved_path(request_path, request_base)
    return None if profile_path is None else _resolved_path(profile_path, profile_base)


def _resolved_decoy(
    profile: DatabaseBuildProfileDocument,
    request: DatabaseBuildRequestDocument,
    override: BuildDecoyOverride | None,
) -> DecoyDocument | None:
    configured = request.decoy if "decoy" in request.model_fields_set else profile.default_decoy
    if override is None:
        return configured
    if override is BuildDecoyOverride.NONE:
        return None
    base = configured or DecoyDocument()
    return base.model_copy(update={"mode": DecoyMode(override.value)})


def _entries_from_paths(paths: tuple[Path, ...]) -> list[Entry]:
    return [
        (record.raw_header, record.sequence)
        for path in paths
        for record in _read_existing_records(path)
    ]


def _read_existing_records(path: Path) -> Iterable[FastaRecord]:
    from protein_fasta.reading.parser import read_records

    try:
        return read_records(path)
    except OSError as error:
        raise ValueError(f"cannot read database-build source {path}: {error}") from error


def _artifact(
    path: Path,
    *,
    relative_to: Path,
    schema_name: str,
    schema_version: str,
    row_count: int | None = None,
) -> BuildArtifactDocument:
    return BuildArtifactDocument(
        schema_name=schema_name,
        schema_version=schema_version,
        path=Path(os.path.relpath(path, start=relative_to)),
        checksum_version=FILE_CHECKSUM_VERSION,
        checksum=file_checksum(path),
        byte_count=path.stat().st_size,
        row_count=row_count,
    )


def _result_document(
    effective: EffectiveDatabaseBuildDocument,
    result: PipelineResult,
    effective_path: Path,
    input_paths: tuple[Path, ...],
    inventory_path: Path,
    inventory_rows: int,
) -> DatabaseBuildResultDocument:
    summary = result.summary
    frequencies = (
        {
            residue: count / summary.total_residues
            for residue, count in summary.aa_frequencies.items()
        }
        if summary.total_residues
        else {}
    )
    artifacts = [
        _artifact(
            effective_path,
            relative_to=effective.output_dir,
            schema_name="effective-database-build",
            schema_version=effective.schema_version,
        ),
        *(
            _artifact(
                path,
                relative_to=effective.output_dir,
                schema_name="protein-fasta-input",
                schema_version="1",
            )
            for path in input_paths
        ),
        _artifact(
            result.path,
            relative_to=effective.output_dir,
            schema_name="protein-fasta",
            schema_version="1",
            row_count=result.n_total,
        ),
        _artifact(
            inventory_path,
            relative_to=effective.output_dir,
            schema_name="protein-inventory",
            schema_version="1",
            row_count=inventory_rows,
        ),
    ]
    if result.entrapment_pairs_path is not None:
        artifacts.append(
            _artifact(
                result.entrapment_pairs_path,
                relative_to=effective.output_dir,
                schema_name="entrapment-pairs-tsv",
                schema_version="legacy-1",
            )
        )
    decoy = None
    if result.decoy is not None:
        decoy = DatabaseBuildDecoyEvidenceDocument(
            mode=result.decoy.mode,
            seed=result.decoy.seed,
            initial_collisions=result.decoy.initial_collisions,
            unresolved_collisions=result.decoy.unresolved_collisions,
            dropped_peptides=result.decoy.dropped_peptides,
            omitted_decoys=result.decoy.omitted_decoys,
        )
    entrapment = None
    if result.entrapment is not None:
        entrapment = DatabaseBuildEntrapmentEvidenceDocument(
            strategy=result.entrapment.strategy,
            seed=result.entrapment.seed,
            requested_fold=result.entrapment.requested_fold,
            achieved_fold=result.entrapment.achieved_fold,
            failures=result.entrapment.failures,
            proteins_affected=result.entrapment.proteins_affected,
            source_proteins=result.entrapment.source_proteins,
        )
    return DatabaseBuildResultDocument(
        protein_fasta_version=importlib.metadata.version("protein_fasta"),
        effective_request=effective,
        artifacts=tuple(artifacts),
        counts=DatabaseBuildCountsDocument(
            target=result.n_target,
            contaminant=result.n_contaminant,
            entrapment=result.n_entrapment,
            decoy=result.n_decoy,
            total=result.n_total,
        ),
        normalization=DatabaseBuildNormalizationDocument(
            upper_cased=result.upper_cased_entries,
            terminal_stops_stripped=result.stop_stripped_entries,
            duplicates_dropped=result.duplicates_dropped,
        ),
        summary=DatabaseBuildSummaryDocument(
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
        ),
        decoy=decoy,
        entrapment=entrapment,
    )


def _write_json_atomic(path: Path, payload: object) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        temporary_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def _normalize_entries(
    entries: Iterable[Entry],
    rules: DiagnosticRules,
    *,
    source: str,
) -> tuple[list[Entry], int, int]:
    """Normalize build input, refusing what normalization cannot make valid.

    Returns the normalized entries plus how many were upper-cased and how many lost
    a trailing stop. Runs before deduplication, so two entries that differ only by
    case or a terminator are recognised as the duplicate they are rather than
    reported as a conflict.
    """
    normalized_entries: list[Entry] = []
    upper_cased = 0
    stop_stripped = 0
    for description, sequence in entries:
        normalized = normalize_sequence(sequence)
        upper_cased += normalized.upper_cased
        stop_stripped += normalized.stop_stripped
        entry_id = parse_header(description).id
        if not normalized.sequence:
            raise ValueError(f"{source}: entry {entry_id} has a header but no sequence.")
        illegal = rules.illegal_residues(normalized.sequence)
        if illegal:
            raise ValueError(
                f"{source}: entry {entry_id} has {describe_illegal_residues(illegal)}."
            )
        normalized_entries.append((description, normalized.sequence))
    return normalized_entries, upper_cased, stop_stripped


class ConflictingIdError(ValueError):
    """Two entries share an id but not a sequence, so one of them is wrong.

    prozor -- and this pipeline until now -- kept whichever came first and dropped
    the other silently. Dropping half of a genuine disagreement is a decision no
    build should make on an operator's behalf.
    """

    def __init__(self, entry_ids: list[str]) -> None:
        self.entry_ids = entry_ids
        shown = ", ".join(entry_ids[:_REPORTED_CONFLICTS])
        remainder = len(entry_ids) - _REPORTED_CONFLICTS
        if remainder > 0:
            shown += f", and {remainder:,} more"
        super().__init__(
            f"{len(entry_ids):,} identifier(s) appear more than once with different sequences: {shown}. "
            "Resolve the sources before building."
        )


def _deduplicate_by_id(
    groups: Iterable[Iterable[Entry]],
) -> tuple[list[Entry], list[int], int]:
    """Keep the first entry per id, refusing ids whose sequences disagree.

    Returns the survivors, the survivor count per input group, and how many
    identical duplicates were dropped. An identical duplicate is prozor's case and
    stays a silent-but-counted drop; a conflicting one raises.
    """
    seen: dict[str, str] = {}
    entries: list[Entry] = []
    survivor_counts: list[int] = []
    duplicates_dropped = 0
    conflicts: list[str] = []
    for group in groups:
        kept = 0
        for description, sequence in group:
            entry_id = parse_header(description).id
            if entry_id in seen:
                if seen[entry_id] != sequence:
                    conflicts.append(entry_id)
                else:
                    duplicates_dropped += 1
                continue
            seen[entry_id] = sequence
            entries.append((description, sequence))
            kept += 1
        survivor_counts.append(kept)
    if conflicts:
        raise ConflictingIdError(conflicts)
    return entries, survivor_counts, duplicates_dropped


def build_database(
    *,
    targets: Iterable[Entry],
    name_fields: dict[str, Any],
    naming: NamingDocument,
    metadata: MetadataDocument,
    diagnostics: RegistryDiagnosticRules,
    output_dir: Path,
    date: datetime.date,
    template: str | None = None,
    contaminant_blocks: Iterable[ContaminantBlock] = (),
    decoy_spec: DecoyDocument | None = DEFAULT_DECOY_SPEC,
    entrapment_spec: EntrapmentDocument | None = None,
    foreign_entries: Iterable[Entry] = (),
    annotation: str = "",
    installer: str | None = None,
) -> PipelineResult:
    """Build a FASTA database from resolved sources and compiled policy.

    ``name_fields`` are substituted into the chosen naming ``template`` (e.g.
    ``{"project": 42261, "dbn": 1, "description": "human"}``). ``targets`` is
    consumed once; pass a list if you need to reuse it.
    """
    resolved_blocks = list(contaminant_blocks)
    decoy_generation = None if decoy_spec is None else make_decoy_generation(decoy_spec)
    entrapment_generation: EntrapmentGeneration | None = (
        None if entrapment_spec is None else _make_entrapment_generation(entrapment_spec)
    )
    sentinel_cfg = metadata
    if installer is not None:
        sentinel_cfg = sentinel_cfg.model_copy(update={"installer": installer})

    filename = build_fasta_name(
        config=naming,
        template=template,
        date=date,
        decoy=decoy_generation is not None,
        entrapment=entrapment_spec is not None,
        **name_fields,
    )
    # Sentinel/database identity is the name without the date or decoy marker,
    # matching the collection (e.g. file p42261_db1_..._20260605.fasta -> aa|p42261_db1_...).
    dbname = build_dbname(config=naming, template=template, **name_fields)

    # Normalize every input before anything reads a sequence: what is written, what
    # is decoyed, and what deduplication compares are then all the same content.
    target_list, upper_cased, stop_stripped = _normalize_entries(
        targets,
        diagnostics.rules,
        source="Selected targets",
    )
    normalized_blocks: list[ContaminantBlock] = []
    for block in resolved_blocks:
        block_entries, block_upper, block_stops = _normalize_entries(
            block.entries,
            diagnostics.rules,
            source=f"Contaminant set {block.name}",
        )
        upper_cased += block_upper
        stop_stripped += block_stops
        normalized_blocks.append(replace(block, entries=tuple(block_entries)))
    resolved_blocks = normalized_blocks
    normalized_foreign, foreign_upper, foreign_stops = _normalize_entries(
        foreign_entries,
        diagnostics.rules,
        source="Foreign entrapment source",
    )
    upper_cased += foreign_upper
    stop_stripped += foreign_stops
    foreign_entries = normalized_foreign
    if entrapment_generation is not None:
        target_list = list(entrapment_generation.normalize(tuple(target_list)))
        resolved_blocks = [
            replace(block, entries=entrapment_generation.normalize(block.entries))
            for block in resolved_blocks
        ]
        foreign_entries = entrapment_generation.normalize(tuple(foreign_entries))
    contaminant_entries = [entry for block in resolved_blocks for entry in block.entries]

    biological_source, _, _ = _deduplicate_by_id([target_list, contaminant_entries])
    entrapment_batch: EntrapmentBatch | None = None
    entrapment_entries: tuple[Entry, ...] = ()
    if entrapment_generation is not None:
        entrapment_batch = entrapment_generation.generate(
            tuple(biological_source),
            foreign_entries=tuple(foreign_entries),
        )
        entrapment_entries = entrapment_batch.entries
        entrapment_note = entrapment_generation.annotation(entrapment_batch)
        annotation = f"{annotation}; {entrapment_note}" if annotation else entrapment_note

    decoy_source, _, _ = _deduplicate_by_id(
        [target_list, contaminant_entries, list(entrapment_entries)]
    )
    decoy_batch = DecoyBatch(
        (),
        {} if decoy_generation is None else decoy_generation.parameters(),
    )
    if decoy_generation is not None:
        prefix = diagnostics.decoy_prefix
        decoy_batch = decoy_generation.generate(tuple(decoy_source), prefix=prefix)
    decoys = decoy_batch.entries

    if decoy_generation is not None:
        decoy_note = decoy_generation.annotation(
            initial_collisions=decoy_batch.initial_collisions,
            dropped_peptides=decoy_batch.dropped_peptides,
        )
        if decoy_note:
            annotation = f"{annotation}; {decoy_note}" if annotation else decoy_note

    # Assemble in prozor order with marker-delimited contaminant blocks, deduplicating
    # by id across the whole file (first occurrence wins, matching prozor's
    # !duplicated). The category counts are the survivors of dedup, so they always
    # reconcile with the written file — important for merged proteomes that share
    # accessions (the multi-select feature's core case).
    sentinel_header = build_sentinel_header(dbname, annotation, date, sentinel_cfg)
    groups: list[Iterable[Entry]] = [
        [(sentinel_header, sentinel_cfg.body_sequence)],
        target_list,
    ]
    contaminant_group_indexes: list[int] = []
    for block in resolved_blocks:
        marker_header = build_section_marker_header(
            f"Cont_{block.name}",
            block.description,
            sentinel_cfg,
        )
        groups.append([(marker_header, sentinel_cfg.marker_body_sequence)])
        groups.append(block.entries)
        contaminant_group_indexes.append(len(groups) - 1)
    groups.append(entrapment_entries)
    entrapment_group_index = len(groups) - 1
    groups.append(decoys)

    final, survivor_counts, duplicates_dropped = _deduplicate_by_id(groups)
    n_target = survivor_counts[1]
    n_contaminant = sum(survivor_counts[index] for index in contaminant_group_indexes)
    n_entrapment = survivor_counts[entrapment_group_index]
    n_decoy = survivor_counts[-1]

    out_path = output_dir / filename
    write_records(
        (FastaRecord(raw_header=header, sequence=sequence) for header, sequence in final),
        out_path,
    )

    pairs_path: Path | None = None
    if entrapment_batch is not None and entrapment_batch.peptide_pairs:
        # Beside the FASTA and named for it: the later FDP evaluation has nothing
        # to join search results on without this mapping, and a sidecar travels
        # with the database the way the build manifest already does.
        #
        # Only when there are pairs. Foreign-species selection produces none,
        # because a foreign protein is not peptide-paired with the target it
        # entraps, and a header-only file would claim a pairing that does not
        # exist -- which a paired FDP estimate would then silently be run against.
        pairs_path = out_path.with_suffix(f"{out_path.suffix}.entrapment_pairs.tsv")
        pairs_path.write_text(
            _format_entrapment_pairs(entrapment_batch.peptide_pairs), encoding="utf-8"
        )

    summary = summarize_sequences(seq for _, seq in final)
    logger.info(
        "built {}: {} target + {} contaminant + {} entrapment + {} decoy = {} entries",
        out_path.name,
        n_target,
        n_contaminant,
        n_entrapment,
        n_decoy,
        len(final),
    )
    decoy_evidence = _decoy_build_evidence(decoy_generation, decoy_batch)
    entrapment_evidence = _entrapment_build_evidence(entrapment_generation, entrapment_batch)
    return PipelineResult(
        path=out_path,
        dbname=dbname,
        n_target=n_target,
        n_contaminant=n_contaminant,
        n_decoy=n_decoy,
        n_total=len(final),
        contaminant_sets=[block.name for block in resolved_blocks],
        summary=summary,
        decoy=decoy_evidence,
        n_entrapment=n_entrapment,
        entrapment=entrapment_evidence,
        entrapment_pairs_path=pairs_path,
        upper_cased_entries=upper_cased,
        stop_stripped_entries=stop_stripped,
        duplicates_dropped=duplicates_dropped,
    )


def build_manifest(
    request: DatabaseBuildDocument,
    result: PipelineResult,
    input_paths: tuple[Path, ...],
) -> dict[str, object]:
    """Return reproducibility evidence for one completed database build."""
    return {
        "schema_version": "0.1",
        "protein_fasta_version": importlib.metadata.version("protein_fasta"),
        "request": request.model_dump(mode="json"),
        "inputs": [
            {
                "path": str(path),
                "checksum_version": FILE_CHECKSUM_VERSION,
                "checksum": file_checksum(path),
            }
            for path in input_paths
        ],
        "output": {
            "path": str(result.path),
            "checksum_version": FILE_CHECKSUM_VERSION,
            "checksum": file_checksum(result.path),
        },
        "counts": {
            "target": result.n_target,
            "contaminant": result.n_contaminant,
            "entrapment": result.n_entrapment,
            "decoy": result.n_decoy,
            "total": result.n_total,
        },
        "normalization": {
            "upper_cased": result.upper_cased_entries,
            "terminal_stops_stripped": result.stop_stripped_entries,
            "duplicates_dropped": result.duplicates_dropped,
        },
        "decoy": {} if result.decoy is None else result.decoy.parameters,
        "entrapment": None
        if result.entrapment is None
        else {
            "strategy": result.entrapment.strategy,
            "seed": result.entrapment.seed,
            "requested_fold": result.entrapment.requested_fold,
            "achieved_fold": result.entrapment.achieved_fold,
            "failures": result.entrapment.failures,
        },
    }


def write_build_manifest(
    request: DatabaseBuildDocument,
    result: PipelineResult,
    input_paths: tuple[Path, ...],
) -> Path:
    """Atomically write a reproducibility manifest beside a completed FASTA."""
    manifest_path = result.path.with_suffix(f"{result.path.suffix}.manifest.json")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{manifest_path.name}.",
        suffix=".tmp",
        dir=manifest_path.parent,
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        temporary_path.write_text(
            json.dumps(
                build_manifest(request, result, input_paths),
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        os.replace(temporary_path, manifest_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return manifest_path
