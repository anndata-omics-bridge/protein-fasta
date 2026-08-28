"""Pydantic documents for protein diagnostic and classification rules."""

from __future__ import annotations

import re
from typing import Literal, Self

from pydantic import Field, field_validator, model_validator

from protein_fasta.schema.base import DocumentBase


def _compile_pattern(pattern: str, *, purpose: str) -> str:
    try:
        re.compile(pattern)
    except re.error as error:
        raise ValueError(f"invalid {purpose} regex {pattern!r}: {error}") from error
    return pattern


class NamedPatternDocument(DocumentBase):
    """One named identifier namespace expression."""

    name: str = Field(min_length=1)
    pattern: str

    @field_validator("pattern")
    @classmethod
    def _valid_pattern(cls, pattern: str) -> str:
        return _compile_pattern(pattern, purpose="identifier namespace")


class DiagnosticDocument(DocumentBase):
    """Versioned alphabet and identifier-namespace rules."""

    schema_version: Literal["0.1"] = "0.1"
    file_version: str = Field(min_length=1)
    allowed_residues: str = Field(min_length=1)
    identifier_namespaces: tuple[NamedPatternDocument, ...]

    @field_validator("allowed_residues")
    @classmethod
    def _unique_upper_residues(cls, residues: str) -> str:
        if residues != residues.upper():
            raise ValueError("allowed_residues must be upper case")
        if len(set(residues)) != len(residues):
            raise ValueError("allowed_residues must not contain duplicates")
        return residues

    @model_validator(mode="after")
    def _unique_namespace_names(self) -> Self:
        names = [rule.name for rule in self.identifier_namespaces]
        if len(names) != len(set(names)):
            raise ValueError("identifier namespace names must be unique")
        return self


class EntryClassifierDocument(DocumentBase):
    """One independent label and its identifier decorations."""

    name: str = Field(min_length=1)
    output_column: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    match_patterns: tuple[str, ...] = ()
    removable_prefix_patterns: tuple[str, ...] = ()
    removable_suffix_patterns: tuple[str, ...] = ()

    @field_validator("match_patterns")
    @classmethod
    def _valid_match_patterns(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        return tuple(_compile_pattern(pattern, purpose="classifier") for pattern in patterns)

    @field_validator("removable_prefix_patterns")
    @classmethod
    def _valid_prefix_patterns(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        for pattern in patterns:
            compiled = re.compile(_compile_pattern(pattern, purpose="removable prefix"))
            if not pattern.startswith("^"):
                raise ValueError(f"removable prefix regex must start with '^': {pattern!r}")
            match = compiled.match("")
            if match is not None:
                raise ValueError(f"removable prefix regex matches empty text: {pattern!r}")
        return patterns

    @field_validator("removable_suffix_patterns")
    @classmethod
    def _valid_suffix_patterns(cls, patterns: tuple[str, ...]) -> tuple[str, ...]:
        for pattern in patterns:
            compiled = re.compile(_compile_pattern(pattern, purpose="removable suffix"))
            if not pattern.endswith("$"):
                raise ValueError(f"removable suffix regex must end with '$': {pattern!r}")
            match = compiled.search("")
            if match is not None:
                raise ValueError(f"removable suffix regex matches empty text: {pattern!r}")
        return patterns

    @model_validator(mode="after")
    def _has_pattern(self) -> Self:
        if not (
            self.match_patterns or self.removable_prefix_patterns or self.removable_suffix_patterns
        ):
            raise ValueError("an entry classifier must declare at least one pattern")
        return self


class EntryClassifierCatalogDocument(DocumentBase):
    """Versioned independent entry classifiers in authored peeling order."""

    schema_version: Literal["0.1"] = "0.1"
    file_version: str = Field(min_length=1)
    classifiers: tuple[EntryClassifierDocument, ...]

    @model_validator(mode="after")
    def _unique_classifier_names(self) -> Self:
        names = [classifier.name for classifier in self.classifiers]
        columns = [classifier.output_column for classifier in self.classifiers]
        if len(names) != len(set(names)):
            raise ValueError("entry classifier names must be unique")
        if len(columns) != len(set(columns)):
            raise ValueError("entry classifier output columns must be unique")
        return self
