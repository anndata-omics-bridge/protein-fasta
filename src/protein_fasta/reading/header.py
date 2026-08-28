"""Generic protein-FASTA header splitting."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ParsedHeader:
    """Identifier token and normalized optional description."""

    id: str
    description: str | None


def parse_header(raw_header: str, /) -> ParsedHeader:
    """Split a raw header into its first nonempty token and description."""
    fields = raw_header.removeprefix(">").split(maxsplit=1)
    if not fields:
        return ParsedHeader(id="", description=None)
    description = " ".join(fields[1].split()) if len(fields) == 2 else None
    return ParsedHeader(id=fields[0], description=description or None)
