"""Generic FASTA serialization."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from protein_fasta.reading.parser import FastaRecord


def write_records(
    records: Iterable[FastaRecord],
    path: Path,
    /,
    *,
    line_width: int = 60,
) -> None:
    """Write lexical records, wrapping sequences at the requested width."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.raw_header}\n")
            if line_width > 0:
                for start in range(0, len(record.sequence), line_width):
                    handle.write(record.sequence[start : start + line_width])
                    handle.write("\n")
            else:
                handle.write(f"{record.sequence}\n")
