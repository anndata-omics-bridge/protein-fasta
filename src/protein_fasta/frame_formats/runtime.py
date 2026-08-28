"""Compiled runtime values for Polars FASTA enrichment."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol, cast

import polars as pl

type RuntimeColumnType = str

_DTYPES: dict[str, type[pl.DataType]] = {
    "string": pl.String,
    "integer": pl.Int64,
    "number": pl.Float64,
    "boolean": pl.Boolean,
}


class _SelectColumns(Protocol):
    def __call__(self, columns: Sequence[str], /) -> pl.DataFrame: ...


class _WithColumns(Protocol):
    def __call__(self, expressions: Iterable[pl.Expr | pl.Series], /) -> pl.DataFrame: ...


class _FilterRows(Protocol):
    def __call__(self, predicate: pl.Series, /) -> pl.DataFrame: ...


class _SortRows(Protocol):
    def __call__(self, column: str, /) -> pl.DataFrame: ...


class _ReplaceStrict(Protocol):
    def __call__(
        self,
        old: Mapping[str, str | bool],
        /,
        *,
        default: None,
    ) -> pl.Expr: ...


class _ConcatStrings(Protocol):
    def __call__(
        self,
        expressions: Sequence[pl.Expr],
        /,
        *,
        separator: str,
        ignore_nulls: bool,
    ) -> pl.Expr: ...


class FrameColumn(Protocol):
    """Produce one configured Polars output expression."""

    @property
    def name(self) -> str:
        """Return the stable output column name."""
        ...

    @property
    def required(self) -> bool:
        """Return whether every selected row must produce a value."""
        ...

    @property
    def column_type(self) -> RuntimeColumnType:
        """Return the configured output type."""
        ...

    def expression(self, header_column: str, /) -> pl.Expr:
        """Return the complete output expression."""
        ...


@dataclass(frozen=True, slots=True)
class RegexFrameColumn:
    """Extract and cast one regex capture group."""

    name: str
    required: bool
    pattern: str
    column_type: RuntimeColumnType
    values: tuple[tuple[str, str], ...]

    def expression(self, header_column: str, /) -> pl.Expr:
        """Extract, translate, and cast this configured column."""
        expression = pl.col(header_column).str.extract(self.pattern, group_index=1)
        if self.values:
            expression = replace_strict(expression, dict(self.values))
        if self.column_type == "boolean":
            expression = replace_strict(
                expression.str.to_lowercase(),
                {"true": True, "false": False, "1": True, "0": False},
            )
        return expression.cast(_DTYPES[self.column_type], strict=True).alias(self.name)


@dataclass(frozen=True, slots=True)
class LiteralFrameColumn:
    """Append one typed format-level constant."""

    name: str
    required: bool
    value: str | int | float | bool
    column_type: RuntimeColumnType

    def expression(self, header_column: str, /) -> pl.Expr:
        """Return a typed literal output expression."""
        del header_column
        return pl.lit(self.value, dtype=_DTYPES[self.column_type]).alias(self.name)


@dataclass(frozen=True, slots=True)
class CompiledFrameParser:
    """One database detector and its configured output columns."""

    format: str
    detection_pattern: str
    columns: tuple[FrameColumn, ...]


@dataclass(frozen=True, slots=True)
class CompiledEntryClassifier:
    """One Polars label and its identifier expressions."""

    name: str
    output_column: str
    match_patterns: tuple[str, ...]
    removable_prefix_patterns: tuple[str, ...]
    removable_suffix_patterns: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CompiledFrameClassifiers:
    """Independent labels in deterministic decoration-peeling order."""

    classifiers: tuple[CompiledEntryClassifier, ...]

    @property
    def output_columns(self) -> tuple[str, ...]:
        """Return classification columns in authored order."""
        return tuple(classifier.output_column for classifier in self.classifiers)


def select_columns(frame: pl.DataFrame, columns: Sequence[str], /) -> pl.DataFrame:
    """Call Polars selection through its narrow used signature."""
    select = cast(_SelectColumns, object.__getattribute__(frame, "select"))
    return select(columns)


def with_columns(
    frame: pl.DataFrame,
    expressions: Iterable[pl.Expr | pl.Series],
    /,
) -> pl.DataFrame:
    """Call Polars column append through its narrow used signature."""
    append = cast(_WithColumns, object.__getattribute__(frame, "with_columns"))
    return append(expressions)


def filter_rows(frame: pl.DataFrame, predicate: pl.Series, /) -> pl.DataFrame:
    """Call Polars row filtering through its narrow used signature."""
    filter_frame = cast(_FilterRows, object.__getattribute__(frame, "filter"))
    return filter_frame(predicate)


def sort_rows(frame: pl.DataFrame, column: str, /) -> pl.DataFrame:
    """Call Polars sorting through its narrow used signature."""
    sort_frame = cast(_SortRows, object.__getattribute__(frame, "sort"))
    return sort_frame(column)


def replace_strict(
    expression: pl.Expr,
    values: Mapping[str, str | bool],
    /,
) -> pl.Expr:
    """Call Polars value replacement through its narrow used signature."""
    replace = cast(_ReplaceStrict, object.__getattribute__(expression, "replace_strict"))
    return replace(values, default=None)


def concat_strings(
    expressions: Sequence[pl.Expr],
    /,
    *,
    separator: str,
    ignore_nulls: bool,
) -> pl.Expr:
    """Call Polars string concatenation through its narrow used signature."""
    concatenate = cast(_ConcatStrings, object.__getattribute__(pl, "concat_str"))
    return concatenate(
        expressions,
        separator=separator,
        ignore_nulls=ignore_nulls,
    )
