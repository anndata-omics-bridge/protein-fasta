"""Passive documents for peptide-level decoy-method comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Self

from pydantic import Field, model_validator

from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase
from protein_fasta.schema.decoy import DecoyStrategyDocument


class DecoyReportRequestDocument(DocumentBase):
    """Selected methods, digestion policy, prefix, and report destination."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path
    decoy_prefix: str = Field(default="REV_", min_length=1)
    digestion: DigestionDocument = Field(default_factory=DigestionDocument)
    strategies: tuple[DecoyStrategyDocument, ...]

    @model_validator(mode="after")
    def validate_strategies(self) -> Self:
        """Require a non-empty set of independently named methods."""
        if not self.strategies:
            raise ValueError("strategies must contain at least one decoy method")
        names = [strategy.type for strategy in self.strategies]
        if len(names) != len(set(names)):
            raise ValueError("strategies must not repeat a decoy method")
        return self


class EffectiveDecoyReportDocument(DocumentBase):
    """Fully resolved and replayable decoy-method report request."""

    schema_version: Literal["0.1"] = "0.1"
    output_parquet: Path
    decoy_prefix: str = Field(min_length=1)
    digestion: DigestionDocument
    strategies: tuple[DecoyStrategyDocument, ...]


class DecoyReportResultDocument(DocumentBase):
    """Checksummed evidence for one completed decoy-method comparison."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectiveDecoyReportDocument
    biological_inventory: ArtifactDocument
    comparison: ArtifactDocument
    target_proteins: int = Field(ge=0)
    target_peptides: int = Field(ge=0)
    target_unique_peptides: int = Field(ge=0)
