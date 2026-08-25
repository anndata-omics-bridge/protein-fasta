"""Pydantic documents for protein-identifier classification policies."""

from __future__ import annotations

import re
from typing import Annotated, Literal

from pydantic import Field, field_validator

from protein_fasta.schema.base import PolicyDocument

DEFAULT_DECOY_CANDIDATES = (
    r"^REV_",
    r"^rev_",
    r"^DECOY_",
    r"^decoy_",
    r"^XXX_",
    r"^reverse_",
)
DEFAULT_CONTAMINANT_CANDIDATES = (
    r"^zz(?:\||_)",
    r"^CON__",
    r"^CON_",
    r"^Cont_",
    r"^contam_",
)


def _validated_patterns(patterns: tuple[str, ...]) -> tuple[str, ...]:
    for pattern in patterns:
        try:
            re.compile(pattern)
        except re.error as error:
            raise ValueError(f"invalid FASTA identifier regex {pattern!r}: {error}") from error
    return patterns


class InferredPatternDocument(PolicyDocument):
    """Infer effective patterns from those observed in one input."""

    mode: Literal["infer"] = "infer"
    candidates: tuple[str, ...]

    @field_validator("candidates")
    @classmethod
    def _valid_candidates(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        return _validated_patterns(patterns)


class ExplicitPatternDocument(PolicyDocument):
    """Use exactly the supplied patterns, including an empty set."""

    mode: Literal["explicit"] = "explicit"
    patterns: tuple[str, ...]

    @field_validator("patterns")
    @classmethod
    def _valid_patterns(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        return _validated_patterns(patterns)


type PatternDocument = Annotated[
    InferredPatternDocument | ExplicitPatternDocument,
    Field(discriminator="mode"),
]


class IdentifierPatternDocument(PolicyDocument):
    """Independent decoy and contaminant annotation-pattern policies."""

    decoy: PatternDocument = InferredPatternDocument(candidates=DEFAULT_DECOY_CANDIDATES)
    contaminant: PatternDocument = InferredPatternDocument(
        candidates=DEFAULT_CONTAMINANT_CANDIDATES
    )


class PatternMatchCountDocument(PolicyDocument):
    """Serialized match count for one tested pattern."""

    pattern: str
    count: int = Field(ge=0)


class ResolvedPatternDocument(PolicyDocument):
    """Serialized effective patterns and count provenance."""

    patterns: tuple[str, ...]
    source: Literal["explicit", "inferred", "none"]
    match_counts: tuple[PatternMatchCountDocument, ...]


class ResolvedIdentifierPatternsDocument(PolicyDocument):
    """Serialized decoy/contaminant pattern resolution for one input."""

    schema_version: Literal["0.1"] = "0.1"
    decoy: ResolvedPatternDocument
    contaminant: ResolvedPatternDocument
    n_identifiers: int = Field(ge=0)


class IdentifierClassificationDocument(PolicyDocument):
    """First-match primary entry classification used by database workflows."""

    sentinel_regexes: tuple[str, ...] = (r"^(REV_)?aa\|",)
    decoy_prefix: str = Field(default="REV_", min_length=1)
    entrapment_regexes: tuple[str, ...] = (r"_p_target$",)
    contaminant_regexes: tuple[str, ...] = (
        r"^(sp|tr)\|Cont_",
        r"^zh\|C[0-9]{4}_",
    )

    @field_validator("sentinel_regexes", "entrapment_regexes", "contaminant_regexes")
    @classmethod
    def _valid_regexes(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        return _validated_patterns(patterns)
