"""Artifact workflows for theoretical peptide construction and comparison."""

from __future__ import annotations

import importlib.metadata
import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import polars as pl
from loguru import logger

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.inventory import read_database_inventory
from protein_fasta.peptide.computation import (
    DigestedProtein,
    PartitionDigestResult,
    peptide_database_from_partitions,
)
from protein_fasta.peptide.executors import (
    MemoryPeptideExecutor,
    WorkspacePeptideExecutor,
)
from protein_fasta.peptide.models import (
    PeptideDatabase,
    PeptideExecutor,
    PeptideProtein,
    PeptideProteinKind,
)
from protein_fasta.reading.parser import FastaRecord
from protein_fasta.reading.writer import write_records
from protein_fasta.registry.backend import factory
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.peptide import (
    DuckDBPeptideExecutionDocument,
    EffectivePeptideBuildDocument,
    EffectivePeptideComparisonDocument,
    MemoryPeptideExecutionDocument,
    PeptideBuildCountsDocument,
    PeptideBuildRequestDocument,
    PeptideBuildResultDocument,
    PeptideComparisonRequestDocument,
    PeptideComparisonResultDocument,
)

PEPTIDE_SCHEMA = pl.Schema(
    {
        "peptide_id": pl.String,
        "sequence": pl.String,
        "length": pl.Int64,
        "missed_cleavages": pl.Int64,
        "mapping_count": pl.Int64,
        "protein_count": pl.Int64,
        "target_count": pl.Int64,
        "contaminant_count": pl.Int64,
        "entrapment_count": pl.Int64,
        "decoy_count": pl.Int64,
    }
)

PROTEIN_PEPTIDE_MAP_SCHEMA = pl.Schema(
    {
        "peptide_id": pl.String,
        "peptide_sequence": pl.String,
        "protein_order": pl.Int64,
        "protein_id": pl.String,
        "protein_kind": pl.String,
        "missed_cleavages": pl.Int64,
    }
)

PEPTIDE_COMPARISON_SCHEMA = pl.Schema(
    {
        "population": pl.String,
        "distinct_peptides_a": pl.Int64,
        "distinct_peptides_b": pl.Int64,
        "shared_peptides": pl.Int64,
        "jaccard": pl.Float64,
        "coverage_a_by_b": pl.Float64,
        "coverage_b_by_a": pl.Float64,
        "containment": pl.Float64,
    }
)


@dataclass(frozen=True, slots=True)
class PeptideBuildExecution:
    """Runtime peptide database plus its durable artifacts."""

    database: PeptideDatabase
    document: PeptideBuildResultDocument
    peptides_path: Path
    mapping_path: Path
    fasta_path: Path
    effective_request_path: Path
    result_path: Path


@dataclass(frozen=True, slots=True)
class PeptideComparisonExecution:
    """Completed exact peptide-set comparison."""

    frame: pl.DataFrame
    document: PeptideComparisonResultDocument
    comparison_path: Path
    effective_request_path: Path
    result_path: Path


def resolve_peptide_build_request(
    request: PeptideBuildRequestDocument,
    /,
    *,
    request_base: Path,
) -> EffectivePeptideBuildDocument:
    """Resolve every peptide artifact path beside its request."""
    return EffectivePeptideBuildDocument(
        peptides_parquet=_resolve_path(request.peptides_parquet, request_base),
        mapping_parquet=_resolve_path(request.mapping_parquet, request_base),
        peptide_fasta=_resolve_path(request.peptide_fasta, request_base),
        digestion=request.digestion,
        execution=request.execution,
    )


def run_peptide_build(
    protein_inventory_path: Path,
    effective: EffectivePeptideBuildDocument,
    /,
) -> PeptideBuildExecution:
    """Digest one canonical protein inventory and atomically publish products."""
    inventory_path = protein_inventory_path.resolve()
    effective_path = effective.peptides_parquet.with_suffix(
        f"{effective.peptides_parquet.suffix}.effective.json"
    )
    result_path = effective.peptides_parquet.with_suffix(
        f"{effective.peptides_parquet.suffix}.result.json"
    )
    _refuse_existing(
        effective.peptides_parquet,
        effective.mapping_parquet,
        effective.peptide_fasta,
        result_path,
    )
    write_json_atomic(effective_path, effective.model_dump(mode="json"), replace_existing=True)
    proteins = _peptide_proteins(inventory_path)
    executor = _make_executor(effective)
    database = executor.execute(proteins, make_digestion(effective.digestion))
    peptides = peptide_frame(database)
    mappings = protein_peptide_map_frame(database)
    with temporary_sibling(effective.peptides_parquet) as staged_peptides:
        peptides.write_parquet(staged_peptides)
        with temporary_sibling(effective.mapping_parquet) as staged_mapping:
            mappings.write_parquet(staged_mapping)
            with temporary_sibling(effective.peptide_fasta) as staged_fasta:
                write_records(
                    (
                        FastaRecord(
                            f"{peptide.peptide_id} mappings={peptide.mapping_count}",
                            peptide.sequence,
                        )
                        for peptide in database.peptides
                    ),
                    staged_fasta,
                )
                document = _peptide_result_document(
                    database,
                    proteins,
                    effective,
                    inventory_path,
                    staged_peptides,
                    staged_mapping,
                    staged_fasta,
                )
                published: list[Path] = []
                try:
                    for staged, destination in (
                        (staged_peptides, effective.peptides_parquet),
                        (staged_mapping, effective.mapping_parquet),
                        (staged_fasta, effective.peptide_fasta),
                    ):
                        publish_exclusive(staged, destination)
                        published.append(destination)
                    write_json_atomic(
                        result_path,
                        document.model_dump(mode="json"),
                        replace_existing=False,
                    )
                except BaseException:
                    for path in published:
                        path.unlink(missing_ok=True)
                    raise
    logger.info(
        "Built {} unique peptides and {} mappings with {} -> {}",
        len(database.peptides),
        len(database.mappings),
        executor.name,
        effective.peptides_parquet,
    )
    return PeptideBuildExecution(
        database,
        document,
        effective.peptides_parquet,
        effective.mapping_parquet,
        effective.peptide_fasta,
        effective_path,
        result_path,
    )


def peptide_frame(database: PeptideDatabase, /) -> pl.DataFrame:
    """Project one peptide database to its canonical inventory frame."""
    return pl.DataFrame(
        [
            {
                "peptide_id": peptide.peptide_id,
                "sequence": peptide.sequence,
                "length": peptide.length,
                "missed_cleavages": peptide.missed_cleavages,
                "mapping_count": peptide.mapping_count,
                "protein_count": peptide.protein_count,
                "target_count": peptide.target_count,
                "contaminant_count": peptide.contaminant_count,
                "entrapment_count": peptide.entrapment_count,
                "decoy_count": peptide.decoy_count,
            }
            for peptide in database.peptides
        ],
        schema=PEPTIDE_SCHEMA,
    )


def protein_peptide_map_frame(database: PeptideDatabase, /) -> pl.DataFrame:
    """Project one peptide database to its canonical protein mapping frame."""
    return pl.DataFrame(
        [
            {
                "peptide_id": mapping.peptide_id,
                "peptide_sequence": mapping.peptide_sequence,
                "protein_order": mapping.protein_order,
                "protein_id": mapping.protein_id,
                "protein_kind": mapping.protein_kind,
                "missed_cleavages": mapping.missed_cleavages,
            }
            for mapping in database.mappings
        ],
        schema=PROTEIN_PEPTIDE_MAP_SCHEMA,
    )


def read_peptides(path: Path, /) -> pl.DataFrame:
    """Read and validate one canonical peptide inventory."""
    return _read_frame(path, PEPTIDE_SCHEMA, "peptide inventory")


def read_protein_peptide_map(path: Path, /) -> pl.DataFrame:
    """Read and validate one canonical protein-peptide mapping."""
    return _read_frame(path, PROTEIN_PEPTIDE_MAP_SCHEMA, "protein-peptide mapping")


def resolve_peptide_comparison_request(
    request: PeptideComparisonRequestDocument,
    /,
    *,
    request_base: Path,
) -> EffectivePeptideComparisonDocument:
    """Resolve one peptide-comparison destination."""
    return EffectivePeptideComparisonDocument(
        output_parquet=_resolve_path(request.output_parquet, request_base)
    )


def run_peptide_comparison(
    peptides_a_path: Path,
    peptides_b_path: Path,
    effective: EffectivePeptideComparisonDocument,
    /,
) -> PeptideComparisonExecution:
    """Compare two canonical peptide populations exactly and atomically."""
    path_a = peptides_a_path.resolve()
    path_b = peptides_b_path.resolve()
    output = effective.output_parquet
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, result_path)
    write_json_atomic(effective_path, effective.model_dump(mode="json"), replace_existing=True)
    frame_a = read_peptides(path_a)
    frame_b = read_peptides(path_b)
    comparison = peptide_comparison_frame(frame_a, frame_b)
    with temporary_sibling(output) as staged:
        comparison.write_parquet(staged)
        document = PeptideComparisonResultDocument(
            protein_fasta_version=importlib.metadata.version("protein-fasta"),
            effective_request=effective,
            peptides_a=_artifact(path_a, output.parent, "peptides", "1", frame_a.height),
            peptides_b=_artifact(path_b, output.parent, "peptides", "1", frame_b.height),
            comparison=artifact_document(
                staged,
                recorded_path=Path(output.name),
                schema_name="peptide-comparisons",
                schema_version="1",
                row_count=comparison.height,
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
    return PeptideComparisonExecution(comparison, document, output, effective_path, result_path)


def peptide_comparison_frame(a: pl.DataFrame, b: pl.DataFrame, /) -> pl.DataFrame:
    """Return exact set statistics for all and each protein-kind population."""
    rows: list[dict[str, object]] = []
    for population, column in (
        ("all", None),
        ("target", "target_count"),
        ("contaminant", "contaminant_count"),
        ("entrapment", "entrapment_count"),
        ("decoy", "decoy_count"),
    ):
        sequences_a = _sequence_set(a, column)
        sequences_b = _sequence_set(b, column)
        shared = len(sequences_a & sequences_b)
        union = len(sequences_a | sequences_b)
        smaller = min(len(sequences_a), len(sequences_b))
        rows.append(
            {
                "population": population,
                "distinct_peptides_a": len(sequences_a),
                "distinct_peptides_b": len(sequences_b),
                "shared_peptides": shared,
                "jaccard": shared / union if union else 0.0,
                "coverage_a_by_b": shared / len(sequences_a) if sequences_a else 0.0,
                "coverage_b_by_a": shared / len(sequences_b) if sequences_b else 0.0,
                "containment": shared / smaller if smaller else 0.0,
            }
        )
    return pl.DataFrame(rows, schema=PEPTIDE_COMPARISON_SCHEMA)


def read_peptide_comparisons(path: Path, /) -> pl.DataFrame:
    """Read and validate one canonical peptide comparison."""
    return _read_frame(path, PEPTIDE_COMPARISON_SCHEMA, "peptide comparison")


def _peptide_proteins(path: Path) -> tuple[PeptideProtein, ...]:
    proteins: list[PeptideProtein] = []
    for entry in read_database_inventory(path):
        if entry.kind == "sentinel":
            continue
        proteins.append(
            PeptideProtein(
                ordinal=entry.final_order,
                identifier=entry.identifier,
                kind=entry.kind,
                sequence=entry.sequence,
            )
        )
    return tuple(proteins)


def _make_executor(effective: EffectivePeptideBuildDocument) -> PeptideExecutor:
    execution = effective.execution
    if isinstance(execution, MemoryPeptideExecutionDocument):
        return MemoryPeptideExecutor(execution.workers, execution.partition_size)
    backend = "duckdb" if isinstance(execution, DuckDBPeptideExecutionDocument) else "sqlite"
    return WorkspacePeptideExecutor(
        execution.workers,
        execution.partition_size,
        backend,
        _RegistryPartitionWorkspace(backend),
    )


@dataclass(frozen=True, slots=True)
class _RegistryPartitionWorkspace:
    """Temporary SQLite or DuckDB storage adapted to the peptide capability."""

    backend: str

    def merge(
        self,
        partitions: Iterable[PartitionDigestResult],
        proteins: tuple[PeptideProtein, ...],
    ) -> PeptideDatabase:
        """Persist completed partitions, then merge them in stable index order."""
        protein_by_order = {protein.ordinal: protein for protein in proteins}
        with tempfile.NamedTemporaryFile(
            prefix="protein-fasta-peptides-",
            suffix=factory.suffix_for(self.backend),
            delete=False,
        ) as handle:
            path = Path(handle.name)
        path.unlink(missing_ok=True)
        try:
            with factory.connect(path, backend=self.backend) as connection:
                connection.execute(
                    "CREATE TABLE partition_rows ("
                    "partition_index INTEGER NOT NULL, protein_order INTEGER NOT NULL, "
                    "protein_id TEXT NOT NULL, protein_kind TEXT NOT NULL, "
                    "sequence TEXT NOT NULL, missed_cleavages INTEGER NOT NULL)"
                )
                for partition in partitions:
                    rows = [
                        (
                            partition.index,
                            digested.protein.ordinal,
                            digested.protein.identifier,
                            digested.protein.kind,
                            sequence,
                            missed_cleavages,
                        )
                        for digested in partition.proteins
                        for sequence, missed_cleavages in digested.peptides
                    ]
                    connection.executemany(
                        "INSERT INTO partition_rows VALUES (?, ?, ?, ?, ?, ?)", rows
                    )
                records = connection.execute(
                    "SELECT partition_index, protein_order, protein_id, protein_kind, "
                    "sequence, missed_cleavages FROM partition_rows "
                    "ORDER BY partition_index, protein_order, sequence"
                ).fetchall()
            grouped: dict[int, dict[tuple[int, str, str], list[tuple[str, int]]]] = {}
            for row in records:
                partition_index = int(row[0])
                protein_key = (int(row[1]), str(row[2]), str(row[3]))
                grouped.setdefault(partition_index, {}).setdefault(protein_key, []).append(
                    (str(row[4]), int(row[5]))
                )
            results = tuple(
                PartitionDigestResult(
                    index,
                    tuple(
                        DigestedProtein(
                            PeptideProtein(
                                order,
                                protein_id,
                                cast("PeptideProteinKind", protein_kind),
                                protein_by_order[order].sequence,
                            ),
                            tuple(peptides),
                        )
                        for (order, protein_id, protein_kind), peptides in sorted(group.items())
                    ),
                )
                for index, group in sorted(grouped.items())
            )
            return peptide_database_from_partitions(results)
        finally:
            path.unlink(missing_ok=True)


def _peptide_result_document(
    database: PeptideDatabase,
    proteins: tuple[PeptideProtein, ...],
    effective: EffectivePeptideBuildDocument,
    inventory_path: Path,
    staged_peptides: Path,
    staged_mapping: Path,
    staged_fasta: Path,
) -> PeptideBuildResultDocument:
    relative_to = effective.peptides_parquet.parent
    return PeptideBuildResultDocument(
        protein_fasta_version=importlib.metadata.version("protein-fasta"),
        effective_request=effective,
        protein_inventory=_artifact(
            inventory_path,
            relative_to,
            "protein-database-inventory",
            "1",
            len(proteins),
        ),
        peptides=artifact_document(
            staged_peptides,
            recorded_path=Path(os.path.relpath(effective.peptides_parquet, start=relative_to)),
            schema_name="peptides",
            schema_version="1",
            row_count=len(database.peptides),
        ),
        protein_peptide_mapping=artifact_document(
            staged_mapping,
            recorded_path=Path(os.path.relpath(effective.mapping_parquet, start=relative_to)),
            schema_name="protein-peptide-map",
            schema_version="1",
            row_count=len(database.mappings),
        ),
        peptide_fasta=artifact_document(
            staged_fasta,
            recorded_path=Path(os.path.relpath(effective.peptide_fasta, start=relative_to)),
            schema_name="unique-peptide-fasta",
            schema_version="1",
            row_count=len(database.peptides),
        ),
        counts=PeptideBuildCountsDocument(
            input_proteins=len(proteins),
            peptides=len(database.peptides),
            mappings=len(database.mappings),
            target_peptides=sum(peptide.target_count > 0 for peptide in database.peptides),
            contaminant_peptides=sum(
                peptide.contaminant_count > 0 for peptide in database.peptides
            ),
            entrapment_peptides=sum(peptide.entrapment_count > 0 for peptide in database.peptides),
            decoy_peptides=sum(peptide.decoy_count > 0 for peptide in database.peptides),
        ),
    )


def _artifact(
    path: Path,
    relative_to: Path,
    schema_name: str,
    schema_version: str,
    row_count: int,
) -> ArtifactDocument:
    return artifact_document(
        path,
        recorded_path=Path(os.path.relpath(path, start=relative_to)),
        schema_name=schema_name,
        schema_version=schema_version,
        row_count=row_count,
    )


def _sequence_set(frame: pl.DataFrame, count_column: str | None) -> set[str]:
    sequences = cast("list[str]", frame.get_column("sequence").to_list())
    if count_column is None:
        return set(sequences)
    counts = cast("list[int]", frame.get_column(count_column).to_list())
    return {sequence for sequence, count in zip(sequences, counts, strict=True) if count > 0}


def _read_frame(path: Path, schema: pl.Schema, label: str) -> pl.DataFrame:
    frame = pl.read_parquet(path)
    if frame.schema != schema:
        raise ValueError(f"invalid {label} schema in {path}: {frame.schema!r}")
    return frame


def _resolve_path(path: Path, base: Path) -> Path:
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _refuse_existing(*paths: Path) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing peptide artifacts: {names}")
