"""Passive documents for peptide construction and comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import Field

from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class MemoryPeptideExecutionDocument(DocumentBase):
    """Request in-memory deterministic partition merging."""

    type: Literal["memory"] = "memory"
    workers: int = Field(default=1, ge=1)
    partition_size: int = Field(default=500, ge=1)


class SQLitePeptideExecutionDocument(DocumentBase):
    """Request a temporary SQLite merge workspace."""

    type: Literal["sqlite"] = "sqlite"
    workers: int = Field(default=1, ge=1)
    partition_size: int = Field(default=500, ge=1)


class DuckDBPeptideExecutionDocument(DocumentBase):
    """Request a temporary DuckDB merge workspace."""

    type: Literal["duckdb"] = "duckdb"
    workers: int = Field(default=1, ge=1)
    partition_size: int = Field(default=500, ge=1)


PeptideExecutionDocument = Annotated[
    MemoryPeptideExecutionDocument
    | SQLitePeptideExecutionDocument
    | DuckDBPeptideExecutionDocument,
    Field(discriminator="type"),
]


class PeptideBuildRequestDocument(DocumentBase):
    """Digestion policy, execution behavior, and explicit artifact paths."""

    schema_version: Literal["0.1"] = "0.1"
    peptides_parquet: Path
    mapping_parquet: Path
    peptide_fasta: Path
    digestion: DigestionDocument = Field(default_factory=DigestionDocument)
    execution: PeptideExecutionDocument = Field(default_factory=MemoryPeptideExecutionDocument)


class EffectivePeptideBuildDocument(DocumentBase):
    """Fully resolved and directly replayable peptide request."""

    schema_version: Literal["0.1"] = "0.1"
    peptides_parquet: Path
    mapping_parquet: Path
    peptide_fasta: Path
    digestion: DigestionDocument
    execution: PeptideExecutionDocument


class PeptideBuildCountsDocument(DocumentBase):
    """Unique peptide and protein-mapping counts."""

    input_proteins: int = Field(ge=0)
    peptides: int = Field(ge=0)
    mappings: int = Field(ge=0)
    target_peptides: int = Field(ge=0)
    contaminant_peptides: int = Field(ge=0)
    entrapment_peptides: int = Field(ge=0)
    decoy_peptides: int = Field(ge=0)


class PeptideBuildResultDocument(DocumentBase):
    """Checksummed evidence for one protein-to-peptide build."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectivePeptideBuildDocument
    protein_inventory: ArtifactDocument
    peptides: ArtifactDocument
    protein_peptide_mapping: ArtifactDocument
    peptide_fasta: ArtifactDocument
    counts: PeptideBuildCountsDocument


class PeptideComparisonRequestDocument(DocumentBase):
    """Explicit destination for one exact peptide-set comparison."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path


class EffectivePeptideComparisonDocument(DocumentBase):
    """Fully resolved peptide comparison destination."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path


class PeptideComparisonResultDocument(DocumentBase):
    """Checksummed evidence for one exact peptide comparison."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectivePeptideComparisonDocument
    peptides_a: ArtifactDocument
    peptides_b: ArtifactDocument
    comparison: ArtifactDocument
