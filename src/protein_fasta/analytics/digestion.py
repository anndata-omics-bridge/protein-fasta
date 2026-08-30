"""Backend-free theoretical protein digestion."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from re import Pattern

_LETTER_SEGMENT = re.compile(r"[A-Z]+")


@dataclass(frozen=True, slots=True)
class Digestion:
    """Compiled enzyme behavior and peptide-selection limits."""

    enzyme: str
    cleavage: Pattern[str]
    min_length: int
    max_length: int
    missed_cleavages: int


@dataclass(frozen=True, slots=True)
class DigestedPeptide:
    """One theoretical peptide and its number of missed cleavage sites."""

    sequence: str
    missed_cleavages: int

    @property
    def length(self) -> int:
        return len(self.sequence)


def protein_residue_length(sequence: str, /) -> int:
    """Count uppercase ASCII residue segments without changing the sequence."""
    _require_normalized(sequence)
    return sum(len(match.group()) for match in _LETTER_SEGMENT.finditer(sequence))


def digest_sequence(sequence: str, digestion: Digestion, /) -> tuple[DigestedPeptide, ...]:
    """Return ordered peptide candidates from an already-normalized sequence."""
    _require_normalized(sequence)
    peptides: list[DigestedPeptide] = []
    for segment_match in _LETTER_SEGMENT.finditer(sequence):
        segment = segment_match.group()
        boundaries = [0, *(match.end() for match in digestion.cleavage.finditer(segment))]
        if boundaries[-1] != len(segment):
            boundaries.append(len(segment))
        fully_cleaved_count = len(boundaries) - 1
        for start in range(fully_cleaved_count):
            for missed_cleavages in range(digestion.missed_cleavages + 1):
                end = start + missed_cleavages + 1
                if end > fully_cleaved_count:
                    break
                peptide = segment[boundaries[start] : boundaries[end]]
                if digestion.min_length <= len(peptide) <= digestion.max_length:
                    peptides.append(DigestedPeptide(peptide, missed_cleavages))
    return tuple(peptides)


def digest_segments(sequence: str, digestion: Digestion, /) -> tuple[str, ...]:
    """Return the complete zero-missed-cleavage partition of a protein."""
    _require_normalized(sequence)
    segments: list[str] = []
    for segment_match in _LETTER_SEGMENT.finditer(sequence):
        segment = segment_match.group()
        start = 0
        for cleavage_match in digestion.cleavage.finditer(segment):
            end = cleavage_match.end()
            segments.append(segment[start:end])
            start = end
        if start < len(segment):
            segments.append(segment[start:])
    return tuple(segments)


def peptide_universe(sequences: Iterable[str], digestion: Digestion, /) -> frozenset[str]:
    """Return every selected peptide from normalized protein sequences."""
    return frozenset(
        peptide.sequence
        for sequence in sequences
        for peptide in digest_sequence(sequence, digestion)
    )


def is_selected_peptide(peptide: str, digestion: Digestion, /) -> bool:
    """Return whether a peptide participates in the configured collision universe."""
    return digestion.min_length <= len(peptide) <= digestion.max_length


def _require_normalized(sequence: str) -> None:
    if sequence != sequence.upper() or any(character.isspace() for character in sequence):
        raise ValueError("digest_sequence requires an explicitly normalized protein sequence")
