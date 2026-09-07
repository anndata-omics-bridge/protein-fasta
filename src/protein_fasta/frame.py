"""Polars protein-FASTA frame APIs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from protein_fasta.analytics.hashing import file_checksum
from protein_fasta.documents import (
    load_builtin_entry_classifier_document,
    load_builtin_header_format_catalog,
)
from protein_fasta.frame_compile import make_frame_classifiers, make_frame_parsers
from protein_fasta.frame_formats.classification import (
    WORKING_HEADER,
    append_classifications,
)
from protein_fasta.frame_formats.detection import format_matches, parser_for
from protein_fasta.frame_formats.extraction import append_parser_columns
from protein_fasta.frame_formats.runtime import (
    CompiledFrameClassifiers,
    CompiledFrameParser,
    filter_rows,
    select_columns,
    sort_rows,
    with_columns,
)
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import read_records
from protein_fasta.schema.diagnostics import EntryClassifierCatalogDocument
from protein_fasta.schema.frame_formats import HeaderFormatCatalogDocument
from protein_fasta.validation.sequence import normalize_sequence

_RAW_HEADER = "__raw_header"
_ROW_INDEX = "__row_index"
_BASE_COLUMNS = ("id", "description", "sequence")
_PROVENANCE_SCHEMA: dict[str, type[pl.DataType]] = {
    "fasta_source_path": pl.String,
    "fasta_source_checksum": pl.String,
    "fasta_source_ordinal": pl.UInt32,
    "fasta_record_ordinal": pl.UInt64,
}
_COLUMN_TYPES: dict[str, type[pl.DataType]] = {
    "string": pl.String,
    "integer": pl.Int64,
    "number": pl.Float64,
    "boolean": pl.Boolean,
}


@dataclass(frozen=True, slots=True)
class ProteinFormat:
    """A packaged FASTA header format selectable by callers."""

    name: str


uniprotkb = ProteinFormat("uniprotkb")
"""The packaged UniProtKB header format."""

refseq = ProteinFormat("refseq")
"""The packaged RefSeq header format."""


class ProteinDatabase:
    """Parse ordered FASTA sources into one configured Polars frame."""

    def __init__(self, *formats: ProteinFormat) -> None:
        """Configure how FASTA headers become protein-frame columns.

        Args:
            *formats: Packaged FASTA formats in recognition and output order.

        Raises:
            ValueError: No formats were supplied, a format is duplicated, or a format is unknown.
        """
        if not formats:
            raise ValueError("ProteinDatabase requires at least one protein format")
        names = tuple(format_.name for format_ in formats)
        if len(names) != len(set(names)):
            raise ValueError("ProteinDatabase formats must be unique")
        available = {
            document.format: document for document in load_builtin_header_format_catalog().formats
        }
        unknown = tuple(name for name in names if name not in available)
        if unknown:
            rendered = ", ".join(repr(name) for name in unknown)
            raise ValueError(f"unknown packaged protein format: {rendered}")
        catalog = HeaderFormatCatalogDocument(
            formats=tuple(available[name] for name in names),
        )
        self._parsers = make_frame_parsers(catalog)
        self._classifiers = make_frame_classifiers(
            load_builtin_entry_classifier_document(),
        )
        _validate_output_names(self._parsers, self._classifiers)
        self._schema = _database_schema(self._parsers, self._classifiers)

    def parse(self, paths: tuple[Path, ...], /) -> pl.DataFrame:
        """Return one source-aware frame for the supplied FASTA paths."""
        frames = tuple(
            self._parse_source(path, source_ordinal) for source_ordinal, path in enumerate(paths)
        )
        if not frames:
            return pl.DataFrame(schema=self._schema)
        return pl.concat(frames, how="vertical")

    def _parse_source(self, path: Path, source_ordinal: int) -> pl.DataFrame:
        frame = _read_with_runtime(path, self._parsers, self._classifiers)
        missing = tuple(name for name in self._schema if name not in frame.columns)
        if missing:
            frame = with_columns(
                frame,
                [pl.lit(None, dtype=self._schema[name]).alias(name) for name in missing],
            )
        frame = select_columns(frame, tuple(self._schema)[:-4])
        return with_columns(
            frame,
            [
                pl.lit(str(path), dtype=pl.String).alias("fasta_source_path"),
                pl.lit(file_checksum(path), dtype=pl.String).alias("fasta_source_checksum"),
                pl.lit(source_ordinal, dtype=pl.UInt32).alias("fasta_source_ordinal"),
                pl.int_range(0, frame.height, dtype=pl.UInt64).alias("fasta_record_ordinal"),
            ],
        )


def read_basic_protein_frame(path: Path, /) -> pl.DataFrame:
    """Return exactly the normalized base protein columns."""
    return select_columns(_read_internal_frame(path), _BASE_COLUMNS)


def read_protein_frame(path: Path, /) -> pl.DataFrame:
    """Enrich each row accepted by exactly one packaged parser."""
    return _read_with_runtime(
        path,
        make_frame_parsers(load_builtin_header_format_catalog()),
        make_frame_classifiers(load_builtin_entry_classifier_document()),
    )


def read_strict_protein_frame(path: Path, /) -> pl.DataFrame:
    """Enrich only when exactly one packaged parser accepts the complete file."""
    return _read_strict_with_runtime(
        path,
        make_frame_parsers(load_builtin_header_format_catalog()),
        make_frame_classifiers(load_builtin_entry_classifier_document()),
    )


def read_configured_protein_frame(
    path: Path,
    catalog: HeaderFormatCatalogDocument,
    classifiers: EntryClassifierCatalogDocument,
    /,
) -> pl.DataFrame:
    """Enrich each row accepted by exactly one explicit parser."""
    return _read_with_runtime(
        path,
        make_frame_parsers(catalog),
        make_frame_classifiers(classifiers),
    )


def read_strict_configured_protein_frame(
    path: Path,
    catalog: HeaderFormatCatalogDocument,
    classifiers: EntryClassifierCatalogDocument,
    /,
) -> pl.DataFrame:
    """Enrich only when one explicit parser accepts the complete file."""
    return _read_strict_with_runtime(
        path,
        make_frame_parsers(catalog),
        make_frame_classifiers(classifiers),
    )


def read_header_format_diagnostics_frame(
    path: Path,
    catalog: HeaderFormatCatalogDocument,
    classifiers: EntryClassifierCatalogDocument,
    /,
) -> pl.DataFrame:
    """Explain whole-frame format recognition without enriching proteins."""
    parsers = make_frame_parsers(catalog)
    classified = append_classifications(
        _read_internal_frame(path),
        make_frame_classifiers(classifiers),
    )
    matches = format_matches(classified[WORKING_HEADER], parsers)
    selected = parser_for(classified[WORKING_HEADER], parsers)
    rows = [
        {
            "format": match.parser.format,
            "matched_rows": match.matched_rows,
            "total_rows": match.total_rows,
            "status": _match_status(match.matched_rows, match.total_rows, match.parser, selected),
        }
        for match in matches
    ]
    return pl.DataFrame(
        rows,
        schema={
            "format": pl.String,
            "matched_rows": pl.Int64,
            "total_rows": pl.Int64,
            "status": pl.String,
        },
    )


def _read_with_runtime(
    path: Path,
    parsers: tuple[CompiledFrameParser, ...],
    classifiers: CompiledFrameClassifiers,
) -> pl.DataFrame:
    internal = _read_internal_frame(path)
    base = select_columns(internal, _BASE_COLUMNS)
    if internal.is_empty():
        return base
    _validate_output_names(parsers, classifiers)
    classified = append_classifications(internal, classifiers)
    indexed = with_columns(
        classified,
        [pl.Series(_ROW_INDEX, range(classified.height), dtype=pl.UInt64)],
    )
    headers = classified.get_column(WORKING_HEADER)
    matches = tuple(headers.str.contains(parser.detection_pattern) for parser in parsers)
    match_counts = pl.Series([0] * classified.height, dtype=pl.UInt32)
    for match in matches:
        match_counts = match_counts + match.cast(pl.UInt32)

    selected = tuple(
        (parser, (match_counts == 1) & match)
        for parser, match in zip(parsers, matches, strict=True)
    )
    used = tuple((parser, mask) for parser, mask in selected if mask.any())
    if not used:
        return base
    parts = [
        append_parser_columns(
            filter_rows(indexed, mask),
            parser,
            WORKING_HEADER,
        )
        for parser, mask in used
    ]
    parts.append(filter_rows(indexed, match_counts != 1))
    enriched = sort_rows(pl.concat(parts, how="diagonal"), _ROW_INDEX)
    output_columns = (
        *_BASE_COLUMNS,
        *classifiers.output_columns,
        *_parser_output_columns(tuple(parser for parser, _ in used)),
    )
    return select_columns(enriched, output_columns)


def _read_strict_with_runtime(
    path: Path,
    parsers: tuple[CompiledFrameParser, ...],
    classifiers: CompiledFrameClassifiers,
) -> pl.DataFrame:
    internal = _read_internal_frame(path)
    base = select_columns(internal, _BASE_COLUMNS)
    if internal.is_empty():
        return base
    _validate_output_names(parsers, classifiers)
    classified = append_classifications(internal, classifiers)
    parser = parser_for(classified[WORKING_HEADER], parsers)
    if parser is None:
        return base
    enriched = append_parser_columns(classified, parser, WORKING_HEADER)
    output_columns = (
        *_BASE_COLUMNS,
        *classifiers.output_columns,
        *(column.name for column in parser.columns),
    )
    return select_columns(enriched, output_columns)


def _parser_output_columns(
    parsers: tuple[CompiledFrameParser, ...],
) -> tuple[str, ...]:
    names: list[str] = []
    for parser in parsers:
        for column in parser.columns:
            if column.name not in names:
                names.append(column.name)
    return tuple(names)


def _database_schema(
    parsers: tuple[CompiledFrameParser, ...],
    classifiers: CompiledFrameClassifiers,
) -> dict[str, type[pl.DataType]]:
    schema: dict[str, type[pl.DataType]] = {
        "id": pl.String,
        "description": pl.String,
        "sequence": pl.String,
    }
    schema.update(dict.fromkeys(classifiers.output_columns, pl.Boolean))
    for parser in parsers:
        for column in parser.columns:
            schema[column.name] = _COLUMN_TYPES[column.column_type]
    schema.update(_PROVENANCE_SCHEMA)
    return schema


def _read_internal_frame(path: Path) -> pl.DataFrame:
    raw_headers: list[str] = []
    identifiers: list[str] = []
    descriptions: list[str | None] = []
    sequences: list[str] = []
    for lexical in read_records(path):
        parsed = parse_header(lexical.raw_header)
        normalized = normalize_sequence(lexical.sequence)
        raw_headers.append(lexical.raw_header)
        identifiers.append(parsed.id)
        descriptions.append(parsed.description)
        sequences.append(normalized.sequence)
    return pl.DataFrame(
        {
            _RAW_HEADER: raw_headers,
            "id": identifiers,
            "description": descriptions,
            "sequence": sequences,
        },
        schema={
            _RAW_HEADER: pl.String,
            "id": pl.String,
            "description": pl.String,
            "sequence": pl.String,
        },
    )


def _validate_output_names(
    parsers: tuple[CompiledFrameParser, ...],
    classifiers: CompiledFrameClassifiers,
) -> None:
    classifier_names = set(classifiers.output_columns)
    parser_types: dict[str, str] = {}
    for parser in parsers:
        overlaps = classifier_names.intersection(column.name for column in parser.columns)
        if overlaps:
            rendered = ", ".join(sorted(overlaps))
            raise ValueError(
                f"format {parser.format!r} replaces classification columns: {rendered}"
            )
        for column in parser.columns:
            previous = parser_types.setdefault(column.name, column.column_type)
            if previous != column.column_type:
                raise ValueError(
                    f"configured column {column.name!r} has incompatible types: "
                    f"{previous!r} and {column.column_type!r}"
                )


def _match_status(
    matched_rows: int,
    total_rows: int,
    parser: CompiledFrameParser,
    selected: CompiledFrameParser | None,
) -> str:
    if total_rows == 0 or matched_rows == 0:
        return "no_match"
    if selected is parser:
        return "selected"
    if matched_rows == total_rows:
        return "ambiguous"
    return "partial"
