"""Fixed normalization for high-level protein sequences."""

from __future__ import annotations

from dataclasses import dataclass

_STOP = "*"


@dataclass(frozen=True, slots=True)
class NormalizedSequence:
    """Normalized sequence together with an exact change audit."""

    sequence: str
    upper_cased: bool
    stop_stripped: bool


def normalize_sequence(sequence: str, /) -> NormalizedSequence:
    """Uppercase a lexical sequence and remove exactly one terminal stop."""
    upper = sequence.upper()
    upper_cased = upper != sequence
    if upper.endswith(_STOP):
        return NormalizedSequence(upper[:-1], upper_cased, stop_stripped=True)
    return NormalizedSequence(upper, upper_cased, stop_stripped=False)
