"""Passive documents for inventory-to-search-database decoy generation."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import Field

from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class ReverseDecoyDocument(DocumentBase):
    """Request deterministic whole-protein reversal."""

    type: Literal["reverse"] = "reverse"


class ShuffleDecoyDocument(DocumentBase):
    """Request deterministic seeded whole-protein shuffling."""

    type: Literal["shuffle"] = "shuffle"
    seed: int


class DecoyPyratDocument(DocumentBase):
    """Request digestion-aware DecoyPYrat generation."""

    type: Literal["decoypyrat"] = "decoypyrat"
    seed: int
    digestion: DigestionDocument = Field(default_factory=DigestionDocument)


DecoyStrategyDocument = Annotated[
    ReverseDecoyDocument | ShuffleDecoyDocument | DecoyPyratDocument,
    Field(discriminator="type"),
]


class DecoyRequestDocument(DocumentBase):
    """One required decoy strategy and explicit search-database destination."""

    schema_version: Literal["0.1"] = "0.1"
    output_fasta: Path
    decoy_prefix: str = Field(default="REV_", min_length=1)
    strategy: DecoyStrategyDocument


class EffectiveDecoyRequestDocument(DocumentBase):
    """Fully resolved and directly replayable decoy request."""

    schema_version: Literal["0.1"] = "0.1"
    output_fasta: Path
    decoy_prefix: str = Field(min_length=1)
    strategy: DecoyStrategyDocument


class DecoyGenerationEvidenceDocument(DocumentBase):
    """Algorithm identity and observed collision or omission outcomes."""

    strategy: Literal["reverse", "shuffle", "decoypyrat"]
    seed: int | None = None
    parameters: dict[str, object]
    initial_collisions: int = Field(ge=0)
    unresolved_collisions: int = Field(ge=0)
    dropped_peptides: int = Field(ge=0)
    omitted_decoys: int = Field(ge=0)


class DecoyCountsDocument(DocumentBase):
    """Physical biological and generated record counts."""

    biological: int = Field(ge=0)
    decoy: int = Field(ge=0)
    total: int = Field(ge=0)


class DecoySummaryDocument(DocumentBase):
    """Scientific-sequence summary for one completed search database."""

    n_sequences: int = Field(ge=0)
    length_min: int | None = Field(default=None, ge=0)
    length_max: int | None = Field(default=None, ge=0)
    length_mean: float | None = Field(default=None, ge=0)
    length_q1: float | None = Field(default=None, ge=0)
    length_median: float | None = Field(default=None, ge=0)
    length_q3: float | None = Field(default=None, ge=0)
    total_residues: int = Field(ge=0)
    aa_counts: dict[str, int]
    aa_frequencies: dict[str, float]


class DecoyResultDocument(DocumentBase):
    """Evidence for one completed biological-inventory decoy operation."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectiveDecoyRequestDocument
    biological_inventory: ArtifactDocument
    search_fasta: ArtifactDocument
    search_inventory: ArtifactDocument
    counts: DecoyCountsDocument
    summary: DecoySummaryDocument
    generation: DecoyGenerationEvidenceDocument
    warnings: tuple[str, ...] = ()
