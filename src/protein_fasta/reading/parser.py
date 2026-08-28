"""Streaming protein-FASTA parsing for text, paths, and compressed paths."""

from __future__ import annotations

import bz2
import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import IO


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """One FASTA record with an uninterpreted header and lexical sequence."""

    raw_header: str
    sequence: str


class FastaReadError(ValueError):
    """A source cannot be decoded or parsed safely as FASTA."""

    def __init__(
        self,
        source_name: str,
        reason: str,
        *,
        line_number: int | None = None,
    ) -> None:
        """Record the source and precise parse failure."""
        self.source_name = source_name
        self.reason = reason
        self.line_number = line_number
        location = f" at line {line_number}" if line_number is not None else ""
        super().__init__(f"{source_name}: {reason}{location}")


def parse_records(
    lines: Iterable[str],
    *,
    source_name: str = "<stream>",
) -> Iterator[FastaRecord]:
    """Parse FASTA records from text lines without owning the input resource."""
    header: str | None = None
    sequence_parts: list[str] = []

    for line_number, raw_line in enumerate(lines, start=1):
        stripped = raw_line.strip()
        if not stripped:
            continue
        if raw_line.startswith(">"):
            if header is not None:
                yield FastaRecord(header, "".join(sequence_parts))
            header = raw_line[1:].rstrip("\r\n")
            sequence_parts = []
            continue
        if header is None:
            raise FastaReadError(
                source_name,
                "sequence content before the first FASTA header",
                line_number=line_number,
            )
        sequence_parts.append("".join(stripped.split()))

    if header is not None:
        yield FastaRecord(header, "".join(sequence_parts))


def read_records(path: Path) -> Iterator[FastaRecord]:
    """Read FASTA records from a plain, gzip, or bzip2 path."""
    try:
        with _open_text(path) as handle:
            yield from parse_records(handle, source_name=str(path))
    except FastaReadError:
        raise
    except UnicodeDecodeError as error:
        raise FastaReadError(str(path), f"file is not valid UTF-8 ({error})") from error
    except (EOFError, OSError) as error:
        raise FastaReadError(str(path), f"file cannot be read as FASTA ({error})") from error


def parse_text(text: str) -> Iterator[FastaRecord]:
    """Parse records from explicit inline FASTA text."""
    yield from parse_records(StringIO(text), source_name="<inline-fasta>")


def read_headers(path: Path) -> Iterator[str]:
    """Stream raw header text without materializing sequences."""
    try:
        with _open_text(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    yield line[1:].rstrip("\r\n")
    except UnicodeDecodeError as error:
        raise FastaReadError(str(path), f"file is not valid UTF-8 ({error})") from error
    except (EOFError, OSError) as error:
        raise FastaReadError(str(path), f"file cannot be read as FASTA ({error})") from error


def _open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    if path.name.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")
