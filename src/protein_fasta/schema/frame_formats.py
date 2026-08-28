"""Pydantic documents for configured FASTA frame enrichment."""

from __future__ import annotations

import re
from typing import Literal, Self

from pydantic import Field, field_validator, model_validator

from protein_fasta.schema.base import DocumentBase

type ColumnType = Literal["string", "integer", "number", "boolean"]
type ColumnValue = str | int | float | bool

_BASE_COLUMNS = frozenset({"id", "description", "sequence"})


def _compile_pattern(pattern: str, *, purpose: str) -> re.Pattern[str]:
    try:
        return re.compile(pattern)
    except re.error as error:
        raise ValueError(f"invalid {purpose} regex {pattern!r}: {error}") from error


class HeaderColumnDocument(DocumentBase):
    """One literal or regex-extracted output column."""

    name: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    type: ColumnType = "string"
    pattern: str | None = None
    value: ColumnValue | None = None
    values: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _valid_source_and_type(self) -> Self:
        has_pattern = self.pattern is not None
        has_value = self.value is not None
        if has_pattern == has_value:
            raise ValueError("a header column must declare exactly one of pattern or value")
        if self.name in _BASE_COLUMNS:
            raise ValueError(f"configured columns cannot replace base column {self.name!r}")
        if self.pattern is not None:
            compiled = _compile_pattern(self.pattern, purpose="header extraction")
            if compiled.groups != 1:
                raise ValueError("a header extraction regex must contain exactly one capture group")
        if self.values and (self.pattern is None or self.type != "string"):
            raise ValueError("values mappings require a regex-extracted string column")
        if self.value is not None and not _value_matches_type(self.value, self.type):
            raise ValueError(f"literal value does not match declared type {self.type!r}")
        return self


class HeaderColumnsDocument(DocumentBase):
    """Required and optional columns in output order."""

    required: tuple[HeaderColumnDocument, ...] = ()
    optional: tuple[HeaderColumnDocument, ...] = ()

    @model_validator(mode="after")
    def _unique_names(self) -> Self:
        names = [column.name for column in (*self.required, *self.optional)]
        if len(names) != len(set(names)):
            raise ValueError("header output column names must be unique")
        return self


class HeaderFormatDocument(DocumentBase):
    """Versioned recognition and extraction rules for one database."""

    schema_version: Literal["0.1"] = "0.1"
    file_version: str = Field(min_length=1)
    format: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    detection_pattern: str
    columns: HeaderColumnsDocument

    @field_validator("detection_pattern")
    @classmethod
    def _valid_detection_pattern(cls, pattern: str) -> str:
        _compile_pattern(pattern, purpose="header detection")
        return pattern


class HeaderFormatCatalogDocument(DocumentBase):
    """Validated in-memory collection of database format documents."""

    formats: tuple[HeaderFormatDocument, ...]

    @model_validator(mode="after")
    def _unique_formats(self) -> Self:
        names = [document.format for document in self.formats]
        if len(names) != len(set(names)):
            raise ValueError("header format names must be unique")
        return self


def _value_matches_type(value: ColumnValue, column_type: ColumnType) -> bool:
    if column_type == "string":
        return isinstance(value, str)
    if column_type == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if column_type == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, bool)
