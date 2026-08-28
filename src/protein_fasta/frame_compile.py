"""Compile frame documents into Polars runtime values."""

from __future__ import annotations

import polars as pl

from protein_fasta.frame_formats.runtime import (
    CompiledEntryClassifier,
    CompiledFrameClassifiers,
    CompiledFrameParser,
    FrameColumn,
    LiteralFrameColumn,
    RegexFrameColumn,
)
from protein_fasta.schema.diagnostics import EntryClassifierCatalogDocument
from protein_fasta.schema.frame_formats import (
    HeaderColumnDocument,
    HeaderFormatCatalogDocument,
)


def make_frame_parsers(
    catalog: HeaderFormatCatalogDocument,
    /,
) -> tuple[CompiledFrameParser, ...]:
    """Compile every database document into one regex-driven parser."""
    parsers: list[CompiledFrameParser] = []
    for document in catalog.formats:
        _validate_polars_pattern(document.detection_pattern, extraction=False)
        required = tuple(
            _make_column(column, required=True) for column in document.columns.required
        )
        optional = tuple(
            _make_column(column, required=False) for column in document.columns.optional
        )
        parsers.append(
            CompiledFrameParser(
                format=document.format,
                detection_pattern=document.detection_pattern,
                columns=(*required, *optional),
            )
        )
    return tuple(parsers)


def make_frame_classifiers(
    document: EntryClassifierCatalogDocument,
    /,
) -> CompiledFrameClassifiers:
    """Compile classifier documents for vectorized frame execution."""
    classifiers: list[CompiledEntryClassifier] = []
    for classifier in document.classifiers:
        patterns = (
            *classifier.match_patterns,
            *classifier.removable_prefix_patterns,
            *classifier.removable_suffix_patterns,
        )
        for pattern in patterns:
            _validate_polars_pattern(pattern, extraction=False)
        classifiers.append(
            CompiledEntryClassifier(
                name=classifier.name,
                output_column=classifier.output_column,
                match_patterns=classifier.match_patterns,
                removable_prefix_patterns=classifier.removable_prefix_patterns,
                removable_suffix_patterns=classifier.removable_suffix_patterns,
            )
        )
    return CompiledFrameClassifiers(tuple(classifiers))


def _make_column(document: HeaderColumnDocument, *, required: bool) -> FrameColumn:
    if document.pattern is not None:
        _validate_polars_pattern(document.pattern, extraction=True)
        return RegexFrameColumn(
            name=document.name,
            required=required,
            pattern=document.pattern,
            column_type=document.type,
            values=tuple(document.values.items()),
        )
    if document.value is None:
        raise ValueError(f"column {document.name!r} has no runtime source")
    return LiteralFrameColumn(
        name=document.name,
        required=required,
        value=document.value,
        column_type=document.type,
    )


def _validate_polars_pattern(pattern: str, *, extraction: bool) -> None:
    try:
        series = pl.Series([""])
        if extraction:
            series.str.extract(pattern, group_index=1)
        else:
            series.str.contains(pattern)
    except pl.exceptions.PolarsError as error:
        raise ValueError(f"regex is not supported by Polars: {pattern!r}: {error}") from error
