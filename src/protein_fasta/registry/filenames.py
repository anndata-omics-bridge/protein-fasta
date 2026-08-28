"""Interpret configured database filenames for indexing and comparison."""

from __future__ import annotations

import datetime
import re
from dataclasses import dataclass

from protein_fasta.schema.build import NamingDocument

_EXTENSIONS = (".fasta.gz", ".fasta.bz2", ".fasta", ".fas", ".fna", ".fa")
_EIGHT_DIGITS = re.compile(r"^\d{8}$")


@dataclass(frozen=True, slots=True)
class ParsedName:
    """Facts recoverable from one existing database filename."""

    dbname: str
    date: datetime.date | None
    is_decoy: bool
    reverse_prefix: bool
    is_entrapment: bool = False


def _strip_extension(filename: str, config: NamingDocument) -> str:
    configured = f".{config.extension}" if config.extension else ""
    candidates = (configured, *_EXTENSIONS) if configured else _EXTENSIONS
    for extension in candidates:
        if extension and filename.endswith(extension):
            return filename[: -len(extension)]
    return filename


def _try_date(token: str) -> datetime.date | None:
    if not _EIGHT_DIGITS.match(token):
        return None
    try:
        return datetime.datetime.strptime(token, "%Y%m%d").date()
    except ValueError:
        return None


def parse_filename(filename: str, config: NamingDocument) -> ParsedName:
    """Recover database name, build date, and decoration flags from a filename."""
    stem = _strip_extension(filename, config)
    reverse_prefix = stem.startswith("R_")
    if reverse_prefix:
        stem = stem.removeprefix("R_")

    tokens = stem.split(config.separator)
    date: datetime.date | None = None
    date_index: int | None = None
    for index, token in enumerate(tokens):
        parsed = _try_date(token)
        if parsed is not None:
            date, date_index = parsed, index

    remaining = [token for index, token in enumerate(tokens) if index != date_index]
    decoy = reverse_prefix
    entrapment = False
    dbname_tokens: list[str] = []
    for token in remaining:
        if token == config.decoy_token:
            decoy = True
        elif token == config.entrapment_token:
            entrapment = True
        else:
            dbname_tokens.append(token)

    return ParsedName(
        dbname=config.separator.join(dbname_tokens),
        date=date,
        is_decoy=decoy,
        reverse_prefix=reverse_prefix,
        is_entrapment=entrapment,
    )
