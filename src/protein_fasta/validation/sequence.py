"""Configured protein-sequence normalization and residue validation."""

from __future__ import annotations

from dataclasses import dataclass

STOP = "*"


@dataclass(frozen=True, slots=True)
class NormalizedSequence:
    """Normalized sequence together with an exact change audit."""

    sequence: str
    upper_cased: bool
    stop_stripped: bool


class ProteinSequenceValidator:
    """Normalize sequences and report residues outside the configured alphabet."""

    __slots__ = ("_delete_allowed", "_strip_trailing_stop")

    def __init__(self, *, allowed_residues: str, strip_trailing_stop: bool) -> None:
        """Compile a translation table for the configured residue alphabet."""
        self._delete_allowed = str.maketrans("", "", allowed_residues)
        self._strip_trailing_stop = strip_trailing_stop

    def normalize(self, sequence: str, /) -> NormalizedSequence:
        """Uppercase a sequence and optionally remove one trailing stop."""
        upper = sequence.upper()
        upper_cased = upper != sequence
        if self._strip_trailing_stop and upper.endswith(STOP):
            return NormalizedSequence(upper[:-1], upper_cased, stop_stripped=True)
        return NormalizedSequence(upper, upper_cased, stop_stripped=False)

    def illegal_residues(self, normalized_sequence: str, /) -> str:
        """Return every character outside the configured allowed alphabet."""
        return normalized_sequence.translate(self._delete_allowed)

    def describe_violation(self, illegal: str, /) -> str:
        """Explain an illegal-residue result."""
        distinct = sorted(set(illegal))
        if distinct == [STOP]:
            return "stop codon inside the sequence, which indicates translated nucleotide output"
        rendered = " ".join(repr(character) for character in distinct)
        return f"illegal sequence characters {rendered}"
