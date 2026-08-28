"""Passive documents for read-only candidate review."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import Field

from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class CandidateRequestDocument(DocumentBase):
    """Candidate-comparison policy and explicit Parquet destination."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path
    overlap_threshold: float = Field(default=0.99, ge=0.0, le=1.0)
    clustering_metric: Literal["target_ids", "target_sequences"] = "target_ids"
    neighbour_limit: int = Field(default=50, ge=1)


class EffectiveCandidateRequestDocument(DocumentBase):
    """Fully resolved candidate review request."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path
    overlap_threshold: float = Field(ge=0.0, le=1.0)
    clustering_metric: Literal["target_ids", "target_sequences"]
    neighbour_limit: int = Field(ge=1)


class CandidateCountsDocument(DocumentBase):
    """Candidate size and registry comparison availability."""

    candidate_records: int = Field(ge=0)
    candidate_targets: int = Field(ge=0)
    candidate_contaminants: int = Field(ge=0)
    checked_databases: int = Field(ge=0)
    excluded_metadata_databases: int = Field(ge=0)
    target_comparisons: int = Field(ge=0)
    contaminant_comparisons: int = Field(ge=0)


class CandidateNeighbourhoodDocument(DocumentBase):
    """Deterministic candidate-neighbour clustering evidence."""

    metric: Literal["target_ids", "target_sequences"]
    relative_paths: tuple[str, ...]
    excluded_empty_paths: tuple[str, ...]
    leaf_order: tuple[int, ...]
    omitted_metadata_paths: tuple[str, ...]
    merge_count: int = Field(ge=0)


class CandidateResultDocument(DocumentBase):
    """Evidence for one completed, non-mutating candidate comparison."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectiveCandidateRequestDocument
    candidate_inventory: ArtifactDocument
    registry: ArtifactDocument
    comparisons: ArtifactDocument
    counts: CandidateCountsDocument
    neighbourhood: CandidateNeighbourhoodDocument
