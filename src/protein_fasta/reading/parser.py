"""Streaming protein-FASTA parsing for text, paths, and compressed paths."""

from __future__ import annotations

import bz2
import gzip
from collections.abc import Iterable, Iterator
from io import StringIO
from pathlib import Path
from typing import IO

from protein_fasta.reading.record import FastaReadError, FastaRecord


def parse_records(
    lines: Iterable[str],
    *,
    source_name: str = "<stream>",
) -> Iterator[FastaRecord]:
    """Parse FASTA records from text lines without owning the input resource."""
    header: str | None = None
    sequence_parts: list[str] = []

    try:
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
    except FastaReadError:
        raise
    except UnicodeDecodeError as error:
        raise FastaReadError(source_name, f"input is not valid UTF-8 ({error})") from error
    except (EOFError, OSError) as error:
        raise FastaReadError(source_name, f"input cannot be read as FASTA ({error})") from error

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
        raise FastaReadError(str(path), f"input is not valid UTF-8 ({error})") from error
    except (EOFError, OSError) as error:
        raise FastaReadError(str(path), f"input cannot be read as FASTA ({error})") from error


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
        raise FastaReadError(str(path), f"input is not valid UTF-8 ({error})") from error
    except (EOFError, OSError) as error:
        raise FastaReadError(str(path), f"input cannot be read as FASTA ({error})") from error


def _open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    if path.name.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")
