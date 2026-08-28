"""Pydantic documents for enzyme and digestion configuration."""

from __future__ import annotations

from typing import Literal, Self

from pydantic import Field, model_validator

from protein_fasta.schema.base import DocumentBase


class EnzymeDocument(DocumentBase):
    """One named cleavage rule authored as a regular expression."""

    schema_version: Literal["0.1"] = "0.1"
    file_version: str = Field(min_length=1)
    name: str = Field(min_length=1)
    cleavage_pattern: str = Field(min_length=1)


class DigestionDocument(DocumentBase):
    """Peptide-length and missed-cleavage limits for one enzyme."""

    enzyme: str = "trypsin"
    min_length: int = Field(default=7, ge=1)
    max_length: int = Field(default=50, ge=1)
    missed_cleavages: int = Field(default=0, ge=0, le=2)

    @model_validator(mode="after")
    def validate_length_range(self) -> Self:
        if self.min_length > self.max_length:
            raise ValueError("min_length must not exceed max_length")
        return self
