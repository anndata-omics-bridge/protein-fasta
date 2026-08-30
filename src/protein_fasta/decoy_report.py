"""Peptide-level comparison of decoy methods from one biological inventory."""

from __future__ import annotations

import importlib.metadata
import os
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from loguru import logger

from protein_fasta.analytics.decoy_diagnostics import (
    DecoyMethodStats,
    PeptidePopulationStats,
    compare_decoy_population,
    peptide_population,
)
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.decoy_compile import make_decoy_generation
from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.inventory import protein_inventory_entries, read_protein_inventory
from protein_fasta.record import iter_protein_diagnostics
from protein_fasta.registry.classification import ContaminantBlockState, classify_record
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.decoy import DecoyStrategyDocument
from protein_fasta.schema.decoy_report import (
    DecoyReportRequestDocument,
    DecoyReportResultDocument,
    EffectiveDecoyReportDocument,
)

DECOY_REPORT_SCHEMA = pl.Schema(
    {
        "method": pl.String,
        "proteins": pl.Int64,
        "peptides": pl.Int64,
        "unique_peptides": pl.Int64,
        "unique_ratio": pl.Float64,
        "length_p10": pl.Float64,
        "length_median": pl.Float64,
        "length_p90": pl.Float64,
        "shared_with_targets": pl.Int64,
        "repeated_peptides": pl.Int64,
        "composition_overlap": pl.Float64,
        "initial_collisions": pl.Int64,
        "unresolved_collisions": pl.Int64,
        "dropped_peptides": pl.Int64,
        "omitted_decoys": pl.Int64,
    }
)


@dataclass(frozen=True, slots=True)
class DecoyReport:
    """Fixed biological baseline and every requested decoy-method statistic."""

    target: PeptidePopulationStats
    methods: tuple[DecoyMethodStats, ...]


@dataclass(frozen=True, slots=True)
class DecoyReportExecution:
    """Runtime report plus its durable artifacts."""

    report: DecoyReport
    frame: pl.DataFrame
    document: DecoyReportResultDocument
    comparison_path: Path
    effective_request_path: Path
    result_path: Path


def resolve_decoy_report_request(
    request: DecoyReportRequestDocument,
    /,
    *,
    request_base: Path,
) -> EffectiveDecoyReportDocument:
    """Resolve one report destination beside its request."""
    output = request.output_parquet
    if not output.is_absolute():
        output = request_base / output
    return EffectiveDecoyReportDocument(
        output_parquet=output.resolve(),
        decoy_prefix=request.decoy_prefix,
        digestion=request.digestion,
        strategies=request.strategies,
    )


def run_decoy_report(
    biological_inventory_path: Path,
    effective: EffectiveDecoyReportDocument,
    /,
) -> DecoyReportExecution:
    """Generate, digest, compare, and atomically publish requested methods."""
    inventory_path = biological_inventory_path.resolve()
    output = effective.output_parquet
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, result_path)
    write_json_atomic(effective_path, effective.model_dump(mode="json"), replace_existing=True)
    entries = protein_inventory_entries(read_protein_inventory(inventory_path))
    sources = tuple(
        (entry.raw_header, entry.sequence)
        for entry in entries
        if entry.kind in {"target", "contaminant", "entrapment"}
    )
    report = compare_decoy_methods(
        sources,
        prefix=effective.decoy_prefix,
        digestion_document=effective.digestion,
        strategies=effective.strategies,
    )
    target = report.target
    frame = decoy_report_frame(report)
    with temporary_sibling(output) as staged:
        frame.write_parquet(staged)
        document = DecoyReportResultDocument(
            protein_fasta_version=importlib.metadata.version("protein-fasta"),
            effective_request=effective,
            biological_inventory=artifact_document(
                inventory_path,
                recorded_path=Path(os.path.relpath(inventory_path, start=output.parent)),
                schema_name="protein-inventory",
                schema_version="1",
                row_count=len(entries),
            ),
            comparison=artifact_document(
                staged,
                recorded_path=Path(output.name),
                schema_name="decoy-method-comparison",
                schema_version="1",
                row_count=frame.height,
            ),
            target_proteins=target.proteins,
            target_peptides=target.peptides,
            target_unique_peptides=target.unique_peptides,
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
    logger.info("Compared {} decoy methods -> {}", len(report.methods), output)
    return DecoyReportExecution(report, frame, document, output, effective_path, result_path)


def compare_decoy_methods(
    sources: tuple[tuple[str, str], ...],
    /,
    *,
    prefix: str,
    digestion_document: DigestionDocument,
    strategies: tuple[DecoyStrategyDocument, ...],
) -> DecoyReport:
    """Compare decoy strategies in memory against one fixed biological source."""
    digestion = make_digestion(digestion_document)
    target, target_peptides = peptide_population((sequence for _, sequence in sources), digestion)
    methods: list[DecoyMethodStats] = []
    for strategy in strategies:
        generation = make_decoy_generation(strategy)
        batch = generation.generate(sources, prefix=prefix)
        methods.append(
            compare_decoy_population(
                method=generation.mode.value,
                target_peptides=target_peptides,
                decoy_sequences=(sequence for _, sequence in batch.entries),
                digestion=digestion,
                initial_collisions=batch.initial_collisions,
                unresolved_collisions=batch.unresolved_collisions,
                dropped_peptides=batch.dropped_peptides,
                omitted_decoys=batch.omitted_decoys,
            )
        )
    return DecoyReport(target, tuple(methods))


def read_decoy_sources(
    path: Path,
    rules: DiagnosticRules,
    decoy_prefix: str,
    /,
) -> tuple[tuple[str, str], ...]:
    """Adapt a legacy FASTA to biological entries used by decoy comparison."""
    entries: list[tuple[str, str]] = []
    block_state: ContaminantBlockState | None = None
    for record in iter_protein_diagnostics(path, rules):
        kind, _, block_state = classify_record(
            record.raw_header,
            record.classifications,
            block_state,
            decoy_prefix,
        )
        if kind in {EntryKind.TARGET, EntryKind.CONTAMINANT, EntryKind.ENTRAPMENT}:
            entries.append((record.raw_header, record.protein.sequence))
    return tuple(entries)


def decoy_report_frame(report: DecoyReport, /) -> pl.DataFrame:
    """Project one report to a canonical baseline-plus-method table."""
    target = report.target
    rows: list[dict[str, object]] = [
        {
            "method": "biological",
            "proteins": target.proteins,
            "peptides": target.peptides,
            "unique_peptides": target.unique_peptides,
            "unique_ratio": 1.0 if target.unique_peptides else 0.0,
            "length_p10": target.length_p10,
            "length_median": target.length_median,
            "length_p90": target.length_p90,
            "shared_with_targets": target.unique_peptides,
            "repeated_peptides": target.repeated_peptides,
            "composition_overlap": 1.0 if target.unique_peptides else 0.0,
            "initial_collisions": 0,
            "unresolved_collisions": 0,
            "dropped_peptides": 0,
            "omitted_decoys": 0,
        }
    ]
    rows.extend(_method_row(method) for method in report.methods)
    return pl.DataFrame(rows, schema=DECOY_REPORT_SCHEMA)


def read_decoy_report(path: Path, /) -> pl.DataFrame:
    """Read and validate one canonical decoy-method comparison."""
    frame = pl.read_parquet(path)
    if frame.schema != DECOY_REPORT_SCHEMA:
        raise ValueError(f"invalid decoy-method comparison schema in {path}: {frame.schema!r}")
    return frame


def _method_row(method: DecoyMethodStats) -> dict[str, object]:
    population = method.population
    return {
        "method": method.method,
        "proteins": population.proteins,
        "peptides": population.peptides,
        "unique_peptides": population.unique_peptides,
        "unique_ratio": method.unique_ratio,
        "length_p10": population.length_p10,
        "length_median": population.length_median,
        "length_p90": population.length_p90,
        "shared_with_targets": method.shared_with_targets,
        "repeated_peptides": population.repeated_peptides,
        "composition_overlap": method.composition_overlap,
        "initial_collisions": method.initial_collisions,
        "unresolved_collisions": method.unresolved_collisions,
        "dropped_peptides": method.dropped_peptides,
        "omitted_decoys": method.omitted_decoys,
    }


def _refuse_existing(*paths: Path) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing decoy-report artifacts: {names}")
