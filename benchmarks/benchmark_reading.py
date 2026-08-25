"""Simple repeatable reader benchmark against a supplied protein FASTA path."""

from __future__ import annotations

import time
from pathlib import Path

from protein_fasta.reading.parser import read_records


def benchmark(path: Path) -> tuple[int, float]:
    """Return the streamed record count and elapsed seconds."""
    started = time.perf_counter()
    count = sum(1 for _ in read_records(path))
    return count, time.perf_counter() - started
