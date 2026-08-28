"""Read-only review of one canonical candidate inventory against a registry."""

from __future__ import annotations

import importlib.metadata
import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from loguru import logger

from protein_fasta.analytics.clustering import ClusteringMetric, DatabaseClustering
from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.inventory import read_database_inventory
from protein_fasta.registry.clustering import cluster_candidate_with_similar_databases
from protein_fasta.registry.comparisons import DatabaseComparison, compare_candidate
from protein_fasta.registry.indexing import (
    SCHEMA_VERSION,
    RegistryRecord,
    RegistrySchemaError,
    RegistrySettings,
    connect_registry,
    open_or_create_registry,
    populate_candidate_entries,
    populate_candidate_files,
)
from protein_fasta.registry.kinds import DetailLevel, EntryKind
from protein_fasta.registry_workflow import make_registry_settings
from protein_fasta.schema.candidate import (
    CandidateCountsDocument,
    CandidateNeighbourhoodDocument,
    CandidateRequestDocument,
    CandidateResultDocument,
    EffectiveCandidateRequestDocument,
)
from protein_fasta.schema.registry import RegistryDocument

CANDIDATE_COMPARISON_SCHEMA = pl.Schema(
    {
        "database_id": pl.Int64,
        "filename": pl.String,
        "annotation": pl.String,
        "kind": pl.String,
        "relationship": pl.String,
        "selected_ids": pl.Int64,
        "other_ids": pl.Int64,
        "shared_ids": pl.Int64,
        "selected_coverage": pl.Float64,
        "other_coverage": pl.Float64,
        "containment": pl.Float64,
        "id_jaccard": pl.Float64,
        "selected_sequences": pl.Int64,
        "other_sequences": pl.Int64,
        "shared_sequences": pl.Int64,
        "sequence_jaccard": pl.Float64,
        "shared_descriptions": pl.Int64,
        "description_jaccard": pl.Float64,
        "shared_exact_pairs": pl.Int64,
        "changed_shared_ids": pl.Int64,
        "selected_only_ids": pl.Int64,
        "other_only_ids": pl.Int64,
        "exact_id_set": pl.Boolean,
        "exact_sequence_set": pl.Boolean,
        "exact_description_set": pl.Boolean,
        "exact_content": pl.Boolean,
    }
)


@dataclass(frozen=True, slots=True)
class CandidateAnalysis:
    """Runtime candidate facts, comparisons, and bounded neighbourhood."""

    candidate: RegistryRecord
    target_comparisons: tuple[DatabaseComparison, ...]
    contaminant_comparisons: tuple[DatabaseComparison, ...]
    checked_database_count: int
    excluded_metadata_database_count: int
    neighbourhood: DatabaseClustering


@dataclass(frozen=True, slots=True)
class CandidateExecution:
    """Completed candidate analysis and its durable artifacts."""

    analysis: CandidateAnalysis
    document: CandidateResultDocument
    comparison_path: Path
    effective_request_path: Path
    result_path: Path


@dataclass(frozen=True, slots=True)
class FastaCandidateAnalysis:
    """Within- and between-database facts for transient FASTA sources."""

    per_file: tuple[RegistryRecord, ...]
    combined: RegistryRecord
    target_comparisons: tuple[DatabaseComparison, ...]
    contaminant_comparisons: tuple[DatabaseComparison, ...]
    checked_database_count: int
    excluded_metadata_database_count: int
    neighbourhood: DatabaseClustering | None = None


def inspect_fasta_candidates(
    paths: list[Path],
    settings: RegistrySettings,
    labels: list[str],
    combined_label: str,
    *,
    kind_override: EntryKind | None = None,
    contaminant_groups: list[str] | None = None,
    on_progress: Callable[[str], None] | None = None,
    strict: bool = True,
    clustering_metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
    neighbour_limit: int = 50,
) -> FastaCandidateAnalysis:
    """Adapt transient FASTA files to the shared candidate computation."""
    with open_or_create_registry(settings) as connection:
        per_file, combined = populate_candidate_files(
            connection,
            paths,
            settings,
            labels=labels,
            combined_label=combined_label,
            kind_override=kind_override,
            contaminant_groups=contaminant_groups,
            on_progress=on_progress,
            strict=strict,
        )
        if on_progress is not None:
            on_progress("Comparing targets against the registered databases ...")
        target_comparisons = tuple(
            compare_candidate(connection, settings.overlap_threshold, kind=EntryKind.TARGET)
        )
        if on_progress is not None:
            on_progress(f"Compared targets against {len(target_comparisons):,} databases")
            on_progress("Comparing contaminants against the registered databases ...")
        contaminant_comparisons = tuple(
            compare_candidate(
                connection,
                settings.overlap_threshold,
                kind=EntryKind.CONTAMINANT,
            )
        )
        if on_progress is not None:
            on_progress(f"Compared contaminants against {len(contaminant_comparisons):,} databases")
        availability = {
            DetailLevel(str(row[0])): int(row[1])
            for row in connection.execute(
                "SELECT detail_level, COUNT(*) FROM databases GROUP BY detail_level"
            ).fetchall()
        }
        neighbourhood = cluster_candidate_with_similar_databases(
            connection,
            combined_label,
            target_comparisons,
            candidate_count=combined.distinct_target_ids or 0,
            metric=clustering_metric,
            limit=neighbour_limit,
        )
    return FastaCandidateAnalysis(
        per_file=tuple(per_file),
        combined=combined,
        target_comparisons=target_comparisons,
        contaminant_comparisons=contaminant_comparisons,
        checked_database_count=availability.get(DetailLevel.FULL, 0),
        excluded_metadata_database_count=availability.get(DetailLevel.METADATA_ONLY, 0),
        neighbourhood=neighbourhood,
    )


def resolve_candidate_request(
    request: CandidateRequestDocument,
    /,
    *,
    request_base: Path,
) -> EffectiveCandidateRequestDocument:
    """Resolve the candidate comparison destination beside its request."""
    output = request.output_parquet
    if not output.is_absolute():
        output = request_base / output
    return EffectiveCandidateRequestDocument(
        output_parquet=output.resolve(),
        overlap_threshold=request.overlap_threshold,
        clustering_metric=request.clustering_metric,
        neighbour_limit=request.neighbour_limit,
    )


def run_candidate_analysis(
    inventory_path: Path,
    registry_path: Path,
    effective: EffectiveCandidateRequestDocument,
    registry_document: RegistryDocument,
    /,
) -> CandidateExecution:
    """Compare one inventory without adding it to or mutating the registry."""
    inventory_path = inventory_path.resolve()
    registry_path = registry_path.resolve()
    output = effective.output_parquet
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    _refuse_existing(output, result_path)
    write_json_atomic(effective_path, effective.model_dump(mode="json"), replace_existing=True)
    entries = read_database_inventory(inventory_path)
    settings = make_registry_settings(
        registry_document,
        fasta_root=inventory_path.parent,
        registry_path=registry_path,
    )
    metric = ClusteringMetric(effective.clustering_metric)
    label = f"{inventory_path.name} [candidate]"
    with connect_registry(registry_path, read_only=True) as connection:
        schema_version = connection.schema_version()
        if schema_version != SCHEMA_VERSION:
            raise RegistrySchemaError(schema_version, path=registry_path)
        candidate = populate_candidate_entries(connection, entries, settings, label=label)
        target_comparisons = tuple(
            compare_candidate(
                connection,
                effective.overlap_threshold,
                kind=EntryKind.TARGET,
            )
        )
        contaminant_comparisons = tuple(
            compare_candidate(
                connection,
                effective.overlap_threshold,
                kind=EntryKind.CONTAMINANT,
            )
        )
        availability = {
            DetailLevel(str(row[0])): int(row[1])
            for row in connection.execute(
                "SELECT detail_level, COUNT(*) FROM databases GROUP BY detail_level"
            ).fetchall()
        }
        neighbourhood = cluster_candidate_with_similar_databases(
            connection,
            label,
            target_comparisons,
            candidate_count=candidate.distinct_target_ids or 0,
            metric=metric,
            limit=effective.neighbour_limit,
        )
    analysis = CandidateAnalysis(
        candidate=candidate,
        target_comparisons=target_comparisons,
        contaminant_comparisons=contaminant_comparisons,
        checked_database_count=availability.get(DetailLevel.FULL, 0),
        excluded_metadata_database_count=availability.get(DetailLevel.METADATA_ONLY, 0),
        neighbourhood=neighbourhood,
    )
    frame = candidate_comparison_frame(analysis)
    with temporary_sibling(output) as staged:
        frame.write_parquet(staged)
        document = _result_document(
            analysis,
            effective,
            inventory_path,
            registry_path,
            staged,
            output,
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
        "Compared candidate with {} full-detail databases -> {}",
        analysis.checked_database_count,
        output,
    )
    return CandidateExecution(analysis, document, output, effective_path, result_path)


def candidate_comparison_frame(analysis: CandidateAnalysis, /) -> pl.DataFrame:
    """Project one runtime candidate analysis to canonical comparison rows."""
    rows = [
        _comparison_row(comparison)
        for comparison in (*analysis.target_comparisons, *analysis.contaminant_comparisons)
    ]
    return pl.DataFrame(rows, schema=CANDIDATE_COMPARISON_SCHEMA)


def read_candidate_comparisons(path: Path, /) -> pl.DataFrame:
    """Read and validate one candidate-comparison Parquet."""
    frame = pl.read_parquet(path)
    if frame.schema != CANDIDATE_COMPARISON_SCHEMA:
        raise ValueError(f"invalid candidate-comparison schema in {path}: {frame.schema!r}")
    return frame


def _comparison_row(comparison: DatabaseComparison) -> dict[str, object]:
    return {
        "database_id": comparison.database.database_id,
        "filename": comparison.database.filename,
        "annotation": comparison.database.annotation,
        "kind": comparison.kind.value,
        "relationship": comparison.relationship.value,
        "selected_ids": comparison.selected_ids,
        "other_ids": comparison.other_ids,
        "shared_ids": comparison.shared_ids,
        "selected_coverage": comparison.selected_coverage,
        "other_coverage": comparison.other_coverage,
        "containment": comparison.containment,
        "id_jaccard": comparison.id_jaccard,
        "selected_sequences": comparison.selected_sequences,
        "other_sequences": comparison.other_sequences,
        "shared_sequences": comparison.shared_sequence_checksums,
        "sequence_jaccard": comparison.sequence_jaccard,
        "shared_descriptions": comparison.shared_descriptions,
        "description_jaccard": comparison.description_jaccard,
        "shared_exact_pairs": comparison.shared_exact_pairs,
        "changed_shared_ids": comparison.changed_shared_ids,
        "selected_only_ids": comparison.selected_only_ids,
        "other_only_ids": comparison.other_only_ids,
        "exact_id_set": comparison.exact_id_set,
        "exact_sequence_set": comparison.exact_sequence_set,
        "exact_description_set": comparison.exact_description_set,
        "exact_content": comparison.exact_content,
    }


def _result_document(
    analysis: CandidateAnalysis,
    effective: EffectiveCandidateRequestDocument,
    inventory_path: Path,
    registry_path: Path,
    staged_output: Path,
    output: Path,
) -> CandidateResultDocument:
    relative_to = output.parent
    neighbourhood = analysis.neighbourhood
    return CandidateResultDocument(
        protein_fasta_version=importlib.metadata.version("protein-fasta"),
        effective_request=effective,
        candidate_inventory=artifact_document(
            inventory_path,
            recorded_path=Path(os.path.relpath(inventory_path, start=relative_to)),
            schema_name="protein-database-inventory",
            schema_version="1",
            row_count=analysis.candidate.entry_count,
        ),
        registry=artifact_document(
            registry_path,
            recorded_path=Path(os.path.relpath(registry_path, start=relative_to)),
            schema_name="protein-database-registry",
            schema_version=str(SCHEMA_VERSION),
            row_count=(analysis.checked_database_count + analysis.excluded_metadata_database_count),
        ),
        comparisons=artifact_document(
            staged_output,
            recorded_path=Path(output.name),
            schema_name="candidate-comparisons",
            schema_version="1",
            row_count=(len(analysis.target_comparisons) + len(analysis.contaminant_comparisons)),
        ),
        counts=CandidateCountsDocument(
            candidate_records=analysis.candidate.entry_count,
            candidate_targets=analysis.candidate.target_count,
            candidate_contaminants=analysis.candidate.contaminant_count,
            checked_databases=analysis.checked_database_count,
            excluded_metadata_databases=analysis.excluded_metadata_database_count,
            target_comparisons=len(analysis.target_comparisons),
            contaminant_comparisons=len(analysis.contaminant_comparisons),
        ),
        neighbourhood=CandidateNeighbourhoodDocument(
            metric=neighbourhood.metric.value,
            relative_paths=neighbourhood.relative_paths,
            excluded_empty_paths=neighbourhood.excluded_empty_paths,
            leaf_order=neighbourhood.leaf_order,
            omitted_metadata_paths=neighbourhood.omitted_metadata_paths,
            merge_count=len(neighbourhood.merges),
        ),
    )


def _refuse_existing(*paths: Path) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing candidate artifacts: {names}")
