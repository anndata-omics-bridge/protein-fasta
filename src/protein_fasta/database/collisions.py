"""Collision-safe sequence permutation for database generation.

Derived from FDR Benchmark revision ``bbf582e`` and modified for the native
``protein_fasta`` runtime model.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Set


def shuffled_candidate(
    sequence: str,
    *,
    rng: random.Random,
    forbidden: Set[str],
    reserved: Set[str] = frozenset(),
    fix_n_term: bool = False,
    fix_c_term: bool = False,
    accepts: Callable[[str], bool] | None = None,
    max_attempts: int,
) -> tuple[str | None, int]:
    """Generate one unique composition-preserving permutation.

    Fixed terminal residues never enter the shuffled core. A core containing
    only one distinct residue has no alternative arrangement and is refused
    without spending the retry budget.
    """
    start = 1 if fix_n_term else 0
    stop = len(sequence) - 1 if fix_c_term else len(sequence)
    if stop <= start:
        return None, 0
    prefix = sequence[:start]
    suffix = sequence[stop:]
    core = list(sequence[start:stop])
    if len(set(core)) == 1:
        return None, 0
    for attempt in range(1, max_attempts + 1):
        shuffled = core.copy()
        rng.shuffle(shuffled)
        candidate = prefix + "".join(shuffled) + suffix
        if candidate == sequence or candidate in forbidden or candidate in reserved:
            continue
        if accepts is None or accepts(candidate):
            return candidate, attempt
    return None, max_attempts
