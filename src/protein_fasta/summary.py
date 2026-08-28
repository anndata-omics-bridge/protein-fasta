"""Streaming protein-sequence summary statistics."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field


def _empty_frequencies() -> dict[str, int]:
    return {}


@dataclass(frozen=True, slots=True)
class FastaSummary:
    """Summary statistics over a set of sequences."""

    n_sequences: int
    length_min: int
    length_q1: float
    length_median: float
    length_mean: float
    length_q3: float
    length_max: int
    total_residues: int
    aa_frequencies: dict[str, int] = field(default_factory=_empty_frequencies)


class SummaryAccumulator:
    """Accumulate exact length and optional residue statistics incrementally."""

    __slots__ = (
        "_aa_counts",
        "_length_counts",
        "_length_max",
        "_length_min",
        "_n_sequences",
        "_total_residues",
    )

    def __init__(self) -> None:
        """Initialize an empty accumulator."""
        self._aa_counts: Counter[str] = Counter()
        self._length_counts: Counter[int] = Counter()
        self._length_min: int | None = None
        self._length_max: int | None = None
        self._n_sequences = 0
        self._total_residues = 0

    def add(self, sequence: str, /) -> None:
        """Add one sequence and its amino-acid composition."""
        self.add_length(len(sequence))
        self._aa_counts.update(sequence)

    def add_length(self, length: int, /) -> None:
        """Add one known sequence length without amino-acid composition."""
        if length < 0:
            raise ValueError("sequence length must be nonnegative")
        self._n_sequences += 1
        self._total_residues += length
        self._length_counts[length] += 1
        if self._length_min is None or length < self._length_min:
            self._length_min = length
        if self._length_max is None or length > self._length_max:
            self._length_max = length

    def merge(self, other: SummaryAccumulator, /) -> None:
        """Merge another accumulator without revisiting its sequences."""
        if not other._n_sequences:
            return
        self._aa_counts.update(other._aa_counts)
        self._length_counts.update(other._length_counts)
        self._n_sequences += other._n_sequences
        self._total_residues += other._total_residues
        if self._length_min is None or (
            other._length_min is not None and other._length_min < self._length_min
        ):
            self._length_min = other._length_min
        if self._length_max is None or (
            other._length_max is not None and other._length_max > self._length_max
        ):
            self._length_max = other._length_max

    def summary(self) -> FastaSummary:
        """Return an immutable snapshot of the accumulated statistics."""
        if not self._n_sequences:
            return _empty_summary()
        if self._length_min is None or self._length_max is None:
            raise RuntimeError("nonempty summary has no length bounds")
        if self._n_sequences == 1:
            q1 = median = q3 = float(self._length_min)
        else:
            q1 = self._inclusive_quantile(1, 4)
            median = self._inclusive_quantile(2, 4)
            q3 = self._inclusive_quantile(3, 4)
        return FastaSummary(
            n_sequences=self._n_sequences,
            length_min=self._length_min,
            length_q1=q1,
            length_median=median,
            length_mean=self._total_residues / self._n_sequences,
            length_q3=q3,
            length_max=self._length_max,
            total_residues=self._total_residues,
            aa_frequencies=dict(sorted(self._aa_counts.items())),
        )

    def _inclusive_quantile(self, numerator: int, denominator: int) -> float:
        scaled_rank = numerator * (self._n_sequences - 1)
        lower_rank, remainder = divmod(scaled_rank, denominator)
        lower = self._length_at_rank(lower_rank)
        upper = self._length_at_rank(lower_rank + 1)
        return (lower * (denominator - remainder) + upper * remainder) / denominator

    def _length_at_rank(self, rank: int) -> int:
        entries_seen = 0
        for length, count in sorted(self._length_counts.items()):
            entries_seen += count
            if rank < entries_seen:
                return length
        raise IndexError(f"length rank {rank} is outside the accumulated values")


def summarize_sequences(sequences: Iterable[str], /) -> FastaSummary:
    """Compute sequence counts, length statistics, and residue frequencies."""
    accumulator = SummaryAccumulator()
    for sequence in sequences:
        accumulator.add(sequence)
    return accumulator.summary()


def _empty_summary() -> FastaSummary:
    return FastaSummary(
        n_sequences=0,
        length_min=0,
        length_q1=0.0,
        length_median=0.0,
        length_mean=0.0,
        length_q3=0.0,
        length_max=0,
        total_residues=0,
        aa_frequencies={},
    )
