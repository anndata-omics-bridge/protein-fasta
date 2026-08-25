"""Pydantic documents for protein-sequence and identifier validation."""

from __future__ import annotations

import re

from pydantic import Field, field_validator

from protein_fasta.schema.base import PolicyDocument


class IdentifierNamespaceDocument(PolicyDocument):
    """One named protein-identifier namespace regular expression."""

    name: str
    pattern: str

    @field_validator("pattern")
    @classmethod
    def _valid_pattern(cls, pattern: str) -> str:
        try:
            re.compile(pattern)
        except re.error as error:
            raise ValueError(f"invalid identifier namespace regex {pattern!r}: {error}") from error
        return pattern


def _default_namespaces() -> tuple[IdentifierNamespaceDocument, ...]:
    return (
        IdentifierNamespaceDocument(name="fgcz_sentinel", pattern=r"^aa\|"),
        IdentifierNamespaceDocument(
            name="fgcz_contaminant",
            pattern=r"^(?:sp|tr)\|Cont_|^zh\|C[0-9]{4}_",
        ),
        IdentifierNamespaceDocument(
            name="uniprot",
            pattern=r"^(?:sp|tr)\|[A-Z0-9]{6,10}(?:-[0-9]+)?(?:\.[0-9]+)?\|\S+$",
        ),
        IdentifierNamespaceDocument(
            name="uniprot_bare",
            pattern=(
                r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
                r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$"
            ),
        ),
        IdentifierNamespaceDocument(name="pdb", pattern=r"^pdb\|[0-9A-Za-z]{4}\|[A-Za-z0-9]$"),
        IdentifierNamespaceDocument(
            name="refseq",
            pattern=r"^(?:AP|NP|XP|YP|ZP|WP)_[0-9]{6,}(?:\.[0-9]+)?$",
        ),
        IdentifierNamespaceDocument(
            name="genbank",
            pattern=(
                r"^[A-Z]{3}[0-9]{5}(?:\.[0-9]+)?$|"
                r"^[A-Z]{2}[0-9]{6}(?:\.[0-9]+)?$"
            ),
        ),
        IdentifierNamespaceDocument(
            name="ensembl",
            pattern=r"^ENS[A-Z]{0,4}[GTP][0-9]{6,}(?:\.[0-9]+)?$",
        ),
        IdentifierNamespaceDocument(name="nextprot", pattern=r"^NX_[A-Z0-9]+$"),
    )


class SequenceValidationDocument(PolicyDocument):
    """Configured protein alphabet, normalization, and namespace reporting."""

    standard_residues: str = "ACDEFGHIKLMNPQRSTVWY"
    tolerated_residues: str = "XBZUOJ"
    strip_trailing_stop: bool = True
    max_reported_id_namespaces: int = Field(default=32, gt=0)
    id_namespaces: tuple[IdentifierNamespaceDocument, ...] = Field(
        default_factory=_default_namespaces
    )
