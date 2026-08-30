"""Prepare ordered FASTA sources as one canonical Parquet build handoff."""

from __future__ import annotations

import importlib.metadata
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl
from loguru import logger

from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.database.models import DecoyInventoryEntry, ProteinInventoryEntry
from protein_fasta.inventory import (
    PROTEIN_INPUT_SCHEMA,
    PROTEIN_INVENTORY_SCHEMA,
    SEARCH_INVENTORY_SCHEMA,
    protein_inventory_entries,
    search_inventory_entries,
)
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import read_records
from protein_fasta.schema.protein_input import (
    ContaminantProteinSourceDocument,
    DerivedProteinInputCountsDocument,
    DerivedProteinInputRequestDocument,
    DerivedProteinInputResultDocument,
    DerivedProteinInputSourceEvidenceDocument,
    ProteinInputNormalizationDocument,
    ProteinInputRequestDocument,
    ProteinInputResultDocument,
    ProteinInputSourceEvidenceDocument,
    ProteinSourceDocument,
)
from protein_fasta.validation.sequence import normalize_sequence


@dataclass(frozen=True, slots=True)
class ProteinSourceEntry:
    """One exact source header and sequence supplied to input preparation."""

    raw_header: str
    sequence: str


@dataclass(frozen=True, slots=True)
class TargetProteinSource:
    """One ordered source of biological target proteins."""

    source_id: str
    entries: tuple[ProteinSourceEntry, ...]


@dataclass(frozen=True, slots=True)
class ContaminantProteinBlock:
    """One named contaminant block with exact source entries."""

    source_id: str
    block_name: str
    block_description: str
    entries: tuple[ProteinSourceEntry, ...]


@dataclass(frozen=True, slots=True)
class ForeignProteinBlock:
    """One foreign source available to biological entrapment generation."""

    source_id: str
    entries: tuple[ProteinSourceEntry, ...]


type ProteinInputSource = TargetProteinSource | ContaminantProteinBlock | ForeignProteinBlock


@dataclass(frozen=True, slots=True)
class PreparedProteinInput:
    """Canonical frame and normalization evidence from in-memory sources."""

    frame: pl.DataFrame
    source_row_counts: tuple[int, ...]
    upper_cased: int
    terminal_stops_stripped: int


@dataclass(frozen=True, slots=True)
class ProteinInputExecution:
    """Canonical frame and durable evidence from one source preparation."""

    frame: pl.DataFrame
    document: ProteinInputResultDocument
    protein_input_path: Path
    effective_request_path: Path
    result_path: Path


@dataclass(frozen=True, slots=True)
class DerivedProteinInputExecution:
    """Canonical source rows and evidence derived from database inventories."""

    frame: pl.DataFrame
    document: DerivedProteinInputResultDocument
    protein_input_path: Path
    effective_request_path: Path
    result_path: Path


def resolve_protein_input_request(
    request: ProteinInputRequestDocument,
    /,
    *,
    request_base: Path,
) -> ProteinInputRequestDocument:
    """Resolve every request-relative path without reading source bytes."""
    sources = tuple(
        source.model_copy(update={"path": _resolved_path(source.path, request_base)})
        for source in request.sources
    )
    return request.model_copy(
        update={
            "sources": sources,
            "output_parquet": _resolved_path(request.output_parquet, request_base),
        }
    )


def resolve_derived_protein_input_request(
    request: DerivedProteinInputRequestDocument,
    /,
    *,
    request_base: Path,
) -> DerivedProteinInputRequestDocument:
    """Resolve request-relative inventory and output paths without reading them."""
    source_inventory = _resolved_path(request.source_inventory, request_base)
    foreign_inventory = (
        None
        if request.foreign_inventory is None
        else _resolved_path(request.foreign_inventory, request_base)
    )
    if foreign_inventory == source_inventory:
        raise ValueError("the foreign inventory must differ from the source inventory")
    return request.model_copy(
        update={
            "source_inventory": source_inventory,
            "foreign_inventory": foreign_inventory,
            "output_parquet": _resolved_path(request.output_parquet, request_base),
        }
    )


def run_protein_input_preparation(
    effective: ProteinInputRequestDocument,
    /,
) -> ProteinInputExecution:
    """Normalize sources once and atomically publish their canonical Parquet."""
    output = effective.output_parquet
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, result_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(
        effective_path,
        effective.model_dump(mode="json"),
        replace_existing=True,
    )
    runtime_sources = tuple(_runtime_source(source) for source in effective.sources)
    prepared = _prepare_protein_sources(runtime_sources)
    source_evidence: list[ProteinInputSourceEvidenceDocument] = []
    for source_order, (source, row_count) in enumerate(
        zip(effective.sources, prepared.source_row_counts, strict=True)
    ):
        source_evidence.append(
            ProteinInputSourceEvidenceDocument(
                source_id=source.source_id,
                source_order=source_order,
                role=source.type,
                artifact=artifact_document(
                    source.path,
                    recorded_path=Path(os.path.relpath(source.path, start=output.parent)),
                    schema_name="protein-fasta-input",
                    schema_version="1",
                    row_count=row_count,
                ),
            )
        )
    frame = prepared.frame
    with temporary_sibling(output) as staged:
        frame.write_parquet(staged)
        input_artifact = artifact_document(
            staged,
            recorded_path=Path(output.name),
            schema_name="protein-input",
            schema_version="1",
            row_count=frame.height,
        )
        document = ProteinInputResultDocument(
            protein_fasta_version=importlib.metadata.version("protein-fasta"),
            effective_request=effective,
            protein_input=input_artifact,
            sources=tuple(source_evidence),
            normalization=ProteinInputNormalizationDocument(
                upper_cased=prepared.upper_cased,
                terminal_stops_stripped=prepared.terminal_stops_stripped,
            ),
        )
        publish_exclusive(staged, output)
        try:
            write_json_atomic(
                result_path,
                document.model_dump(mode="json"),
                replace_existing=False,
            )
        except BaseException:
            output.unlink(missing_ok=True)
            raise
    logger.info(
        "Prepared {} protein rows from {} sources -> {}",
        frame.height,
        len(effective.sources),
        output,
    )
    return ProteinInputExecution(frame, document, output, effective_path, result_path)


def prepare_protein_input_frame(
    target_sources: tuple[TargetProteinSource, ...],
    contaminant_blocks: tuple[ContaminantProteinBlock, ...] = (),
    foreign_blocks: tuple[ForeignProteinBlock, ...] = (),
    /,
) -> PreparedProteinInput:
    """Prepare canonical rows directly from ordered in-memory protein sources.

    Source groups are emitted in target, contaminant, then foreign order. Callers
    requiring arbitrary interleaving should use the request/file adapter, which
    preserves the authored source order while invoking the same computation.
    """
    if not target_sources:
        raise ValueError("protein input requires at least one target source")
    return _prepare_protein_sources((*target_sources, *contaminant_blocks, *foreign_blocks))


def run_derived_protein_input_preparation(
    effective: DerivedProteinInputRequestDocument,
    /,
) -> DerivedProteinInputExecution:
    """Select clean target, contaminant, and optional foreign rows from inventories."""
    output = effective.output_parquet
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, result_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(
        effective_path,
        effective.model_dump(mode="json"),
        replace_existing=True,
    )

    source_entries, source_schema = _inventory_entries(effective.source_inventory)
    source_rows, source_counts = _derived_rows(
        source_entries,
        source_order=0,
        source_id=effective.source_id,
        foreign=False,
    )
    rows = list(source_rows)
    sources = [
        _derived_source_evidence(
            effective.source_inventory,
            output=output,
            source_id=effective.source_id,
            source_order=0,
            purpose="biological",
            schema_name=source_schema,
            row_count=len(source_entries),
        )
    ]
    totals = source_counts

    if effective.foreign_inventory is not None and effective.foreign_source_id is not None:
        foreign_entries, foreign_schema = _inventory_entries(effective.foreign_inventory)
        foreign_rows, foreign_counts = _derived_rows(
            foreign_entries,
            source_order=1,
            source_id=effective.foreign_source_id,
            foreign=True,
        )
        rows.extend(foreign_rows)
        totals = _add_derived_counts(totals, foreign_counts)
        sources.append(
            _derived_source_evidence(
                effective.foreign_inventory,
                output=output,
                source_id=effective.foreign_source_id,
                source_order=1,
                purpose="foreign",
                schema_name=foreign_schema,
                row_count=len(foreign_entries),
            )
        )

    if totals.target == 0:
        raise ValueError(f"source inventory {effective.source_inventory} contains no target rows")
    frame = pl.DataFrame(rows, schema=PROTEIN_INPUT_SCHEMA)
    with temporary_sibling(output) as staged:
        frame.write_parquet(staged)
        input_artifact = artifact_document(
            staged,
            recorded_path=Path(output.name),
            schema_name="protein-input",
            schema_version="1",
            row_count=frame.height,
        )
        document = DerivedProteinInputResultDocument(
            protein_fasta_version=importlib.metadata.version("protein-fasta"),
            effective_request=effective,
            protein_input=input_artifact,
            sources=tuple(sources),
            counts=totals,
        )
        publish_exclusive(staged, output)
        try:
            write_json_atomic(
                result_path,
                document.model_dump(mode="json"),
                replace_existing=False,
            )
        except BaseException:
            output.unlink(missing_ok=True)
            raise
    logger.info(
        "Derived {} target, {} contaminant, and {} foreign rows -> {}",
        totals.target,
        totals.contaminant,
        totals.foreign,
        output,
    )
    return DerivedProteinInputExecution(frame, document, output, effective_path, result_path)


def _inventory_entries(
    path: Path,
) -> tuple[tuple[ProteinInventoryEntry | DecoyInventoryEntry, ...], str]:
    """Read one canonical database inventory and retain its precise schema identity."""
    frame = pl.read_parquet(path)
    if frame.schema == PROTEIN_INVENTORY_SCHEMA:
        return protein_inventory_entries(frame), "protein-inventory"
    if frame.schema == SEARCH_INVENTORY_SCHEMA:
        return search_inventory_entries(frame), "search-inventory"
    raise ValueError(f"invalid biological/search inventory schema in {path}: {frame.schema!r}")


def _derived_rows(
    entries: tuple[ProteinInventoryEntry | DecoyInventoryEntry, ...],
    *,
    source_order: int,
    source_id: str,
    foreign: bool,
) -> tuple[list[dict[str, object]], DerivedProteinInputCountsDocument]:
    """Project the same source categories retained by the former application workflow."""
    rows: list[dict[str, object]] = []
    target = 0
    contaminant = 0
    foreign_count = 0
    skipped_sentinel = 0
    skipped_entrapment = 0
    skipped_decoy = 0
    for entry in entries:
        if isinstance(entry, DecoyInventoryEntry):
            skipped_decoy += 1
            continue
        if entry.kind == "sentinel":
            skipped_sentinel += 1
            continue
        if entry.kind == "entrapment":
            skipped_entrapment += 1
            continue
        if foreign:
            role = "foreign"
            block_name = None
            block_description = None
            foreign_count += 1
        elif entry.kind == "contaminant":
            role = "contaminant"
            block_name = entry.contaminant_group or "derived-contaminants"
            block_description = f"derived from {source_id}"
            contaminant += 1
        else:
            role = "target"
            block_name = None
            block_description = None
            target += 1
        rows.append(
            {
                "source_order": source_order,
                "record_order": entry.final_order,
                "source_id": source_id,
                "role": role,
                "block_name": block_name,
                "block_description": block_description,
                "raw_header": entry.raw_header,
                "id": entry.identifier,
                "description": entry.description,
                "sequence": entry.sequence,
                "upper_cased": False,
                "terminal_stop_stripped": False,
            }
        )
    return rows, DerivedProteinInputCountsDocument(
        target=target,
        contaminant=contaminant,
        foreign=foreign_count,
        skipped_sentinel=skipped_sentinel,
        skipped_entrapment=skipped_entrapment,
        skipped_decoy=skipped_decoy,
    )


def _add_derived_counts(
    left: DerivedProteinInputCountsDocument,
    right: DerivedProteinInputCountsDocument,
) -> DerivedProteinInputCountsDocument:
    """Combine independent source-accounting evidence."""
    return DerivedProteinInputCountsDocument(
        target=left.target + right.target,
        contaminant=left.contaminant + right.contaminant,
        foreign=left.foreign + right.foreign,
        skipped_sentinel=left.skipped_sentinel + right.skipped_sentinel,
        skipped_entrapment=left.skipped_entrapment + right.skipped_entrapment,
        skipped_decoy=left.skipped_decoy + right.skipped_decoy,
    )


def _derived_source_evidence(
    path: Path,
    *,
    output: Path,
    source_id: str,
    source_order: int,
    purpose: Literal["biological", "foreign"],
    schema_name: str,
    row_count: int,
) -> DerivedProteinInputSourceEvidenceDocument:
    """Record one exact database inventory consumed by derivation."""
    return DerivedProteinInputSourceEvidenceDocument(
        source_id=source_id,
        source_order=source_order,
        purpose=purpose,
        artifact=artifact_document(
            path,
            recorded_path=Path(os.path.relpath(path, start=output.parent)),
            schema_name=schema_name,
            schema_version="1",
            row_count=row_count,
        ),
    )


def _prepare_protein_sources(
    sources: tuple[ProteinInputSource, ...],
) -> PreparedProteinInput:
    rows: list[dict[str, object]] = []
    row_counts: list[int] = []
    upper_cased = 0
    stops_stripped = 0
    for source_order, source in enumerate(sources):
        source_rows, source_upper, source_stops = _source_rows(source, source_order)
        rows.extend(source_rows)
        row_counts.append(len(source_rows))
        upper_cased += source_upper
        stops_stripped += source_stops
    return PreparedProteinInput(
        pl.DataFrame(rows, schema=PROTEIN_INPUT_SCHEMA),
        tuple(row_counts),
        upper_cased,
        stops_stripped,
    )


def _runtime_source(source: ProteinSourceDocument) -> ProteinInputSource:
    entries = tuple(
        ProteinSourceEntry(record.raw_header, record.sequence)
        for record in read_records(source.path)
    )
    if isinstance(source, ContaminantProteinSourceDocument):
        return ContaminantProteinBlock(
            source.source_id,
            source.block_name,
            source.block_description,
            entries,
        )
    if source.type == "foreign":
        return ForeignProteinBlock(source.source_id, entries)
    return TargetProteinSource(source.source_id, entries)


def _source_rows(
    source: ProteinInputSource,
    source_order: int,
) -> tuple[list[dict[str, object]], int, int]:
    rows: list[dict[str, object]] = []
    upper_cased = 0
    stops_stripped = 0
    block_name: str | None = None
    block_description: str | None = None
    role: Literal["target", "contaminant", "foreign"] = "target"
    if isinstance(source, ContaminantProteinBlock):
        role = "contaminant"
        block_name = source.block_name
        block_description = source.block_description
    elif isinstance(source, ForeignProteinBlock):
        role = "foreign"
    for record_order, entry in enumerate(source.entries):
        normalized = normalize_sequence(entry.sequence)
        if not normalized.sequence:
            raise ValueError(
                f"source {source.source_id!r}: entry {entry.raw_header!r} has no sequence"
            )
        header = parse_header(entry.raw_header)
        upper_cased += int(normalized.upper_cased)
        stops_stripped += int(normalized.stop_stripped)
        rows.append(
            {
                "source_order": source_order,
                "record_order": record_order,
                "source_id": source.source_id,
                "role": role,
                "block_name": block_name,
                "block_description": block_description,
                "raw_header": entry.raw_header,
                "id": header.id,
                "description": header.description,
                "sequence": normalized.sequence,
                "upper_cased": normalized.upper_cased,
                "terminal_stop_stripped": normalized.stop_stripped,
            }
        )
    return rows, upper_cased, stops_stripped


def _resolved_path(path: Path, base: Path) -> Path:
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _refuse_existing(output: Path, result_path: Path) -> None:
    existing = [path for path in (output, result_path) if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing protein-input artifacts: {names}")
