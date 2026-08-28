"""Whole-frame database format recognition."""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from protein_fasta.frame_formats.runtime import CompiledFrameParser


@dataclass(frozen=True, slots=True)
class FormatMatch:
    """Matched row count for one configured format."""

    parser: CompiledFrameParser
    matched_rows: int
    total_rows: int


def format_matches(
    headers: pl.Series,
    parsers: tuple[CompiledFrameParser, ...],
    /,
) -> tuple[FormatMatch, ...]:
    """Count complete-header matches for every configured parser."""
    return tuple(
        FormatMatch(
            parser=parser,
            matched_rows=int(headers.str.contains(parser.detection_pattern).sum() or 0),
            total_rows=headers.len(),
        )
        for parser in parsers
    )


def parser_for(
    headers: pl.Series,
    parsers: tuple[CompiledFrameParser, ...],
    /,
) -> CompiledFrameParser | None:
    """Select the sole parser matching every nonempty frame row."""
    if headers.len() == 0:
        return None
    complete = [
        match.parser
        for match in format_matches(headers, parsers)
        if match.matched_rows == match.total_rows
    ]
    return complete[0] if len(complete) == 1 else None
