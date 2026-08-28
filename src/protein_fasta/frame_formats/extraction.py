"""Native Polars extraction for one selected database parser."""

from __future__ import annotations

import polars as pl

from protein_fasta.frame_formats.runtime import CompiledFrameParser, with_columns


class FrameExtractionError(ValueError):
    """A selected format failed its configured required-column contract."""


def append_parser_columns(
    frame: pl.DataFrame,
    parser: CompiledFrameParser,
    header_column: str,
    /,
) -> pl.DataFrame:
    """Append every configured column and enforce required extraction."""
    enriched = with_columns(
        frame,
        (column.expression(header_column) for column in parser.columns),
    )
    missing = [
        column.name
        for column in parser.columns
        if column.required and enriched.get_column(column.name).null_count()
    ]
    if missing:
        rendered = ", ".join(missing)
        raise FrameExtractionError(
            f"format {parser.format!r} did not extract required columns: {rendered}"
        )
    return enriched
