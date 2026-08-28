"""Passive documents for ordered FASTA source preparation."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import Field, model_validator

from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class TargetProteinSourceDocument(DocumentBase):
    """One ordered biological target FASTA source."""

    type: Literal["target"] = "target"
    source_id: str = Field(min_length=1)
    path: Path


class ContaminantProteinSourceDocument(DocumentBase):
    """One ordered named contaminant FASTA source."""

    type: Literal["contaminant"] = "contaminant"
    source_id: str = Field(min_length=1)
    path: Path
    block_name: str = Field(min_length=1)
    block_description: str = ""


class ForeignProteinSourceDocument(DocumentBase):
    """One ordered foreign source available to entrapment generation."""

    type: Literal["foreign"] = "foreign"
    source_id: str = Field(min_length=1)
    path: Path


ProteinSourceDocument = Annotated[
    TargetProteinSourceDocument | ContaminantProteinSourceDocument | ForeignProteinSourceDocument,
    Field(discriminator="type"),
]


class ProteinInputRequestDocument(DocumentBase):
    """Ordered FASTA sources and destination for one canonical protein input."""

    schema_version: Literal["0.1"] = "0.1"
    sources: tuple[ProteinSourceDocument, ...] = Field(min_length=1)
    output_parquet: Path


class ProteinInputSourceEvidenceDocument(DocumentBase):
    """Exact source bytes and rows represented in one prepared input."""

    source_id: str
    source_order: int = Field(ge=0)
    role: Literal["target", "contaminant", "foreign"]
    artifact: ArtifactDocument


class ProteinInputNormalizationDocument(DocumentBase):
    """Aggregate sequence changes made while preparing source rows."""

    upper_cased: int = Field(ge=0)
    terminal_stops_stripped: int = Field(ge=0)


class ProteinInputResultDocument(DocumentBase):
    """Evidence for one successfully prepared canonical protein input."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: ProteinInputRequestDocument
    protein_input: ArtifactDocument
    sources: tuple[ProteinInputSourceEvidenceDocument, ...]
    normalization: ProteinInputNormalizationDocument
    warnings: tuple[str, ...] = ()


class DerivedProteinInputRequestDocument(DocumentBase):
    """Select clean biological rows from existing database inventories."""

    schema_version: Literal["0.1"] = "0.1"
    source_inventory: Path
    source_id: str = Field(min_length=1)
    foreign_inventory: Path | None = None
    foreign_source_id: str | None = None
    output_parquet: Path

    @model_validator(mode="after")
    def validate_foreign_source(self) -> DerivedProteinInputRequestDocument:
        """Require the foreign source identity and artifact together."""
        if (self.foreign_inventory is None) != (self.foreign_source_id is None):
            raise ValueError(
                "foreign_inventory and foreign_source_id must either both be supplied or both be absent"
            )
        return self


class DerivedProteinInputSourceEvidenceDocument(DocumentBase):
    """One canonical inventory consumed while deriving clean source rows."""

    source_id: str
    source_order: int = Field(ge=0)
    purpose: Literal["biological", "foreign"]
    artifact: ArtifactDocument


class DerivedProteinInputCountsDocument(DocumentBase):
    """Selected and deliberately excluded rows from the source inventories."""

    target: int = Field(ge=0)
    contaminant: int = Field(ge=0)
    foreign: int = Field(ge=0)
    skipped_sentinel: int = Field(ge=0)
    skipped_entrapment: int = Field(ge=0)
    skipped_decoy: int = Field(ge=0)


class DerivedProteinInputResultDocument(DocumentBase):
    """Evidence for one inventory-to-protein-input derivation."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: DerivedProteinInputRequestDocument
    protein_input: ArtifactDocument
    sources: tuple[DerivedProteinInputSourceEvidenceDocument, ...]
    counts: DerivedProteinInputCountsDocument
    warnings: tuple[str, ...] = ()
