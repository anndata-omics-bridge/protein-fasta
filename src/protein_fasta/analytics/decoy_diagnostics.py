"""Backend-free peptide-level diagnostics for decoy populations."""

from __future__ import annotations

import statistics
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass

from protein_fasta.analytics.digestion import Digestion, digest_sequence

_DECILE_COUNT = 10


@dataclass(frozen=True, slots=True)
class PeptidePopulationStats:
    """Exact summary of one digested protein population."""

    proteins: int
    peptides: int
    unique_peptides: int
    length_p10: float
    length_median: float
    length_p90: float
    repeated_peptides: int


@dataclass(frozen=True, slots=True)
class DecoyMethodStats:
    """One generated method measured against a target peptide population."""

    method: str
    population: PeptidePopulationStats
    unique_ratio: float
    shared_with_targets: int
    composition_overlap: float
    initial_collisions: int
    unresolved_collisions: int
    dropped_peptides: int
    omitted_decoys: int


def peptide_population(
    sequences: Iterable[str],
    digestion: Digestion,
) -> tuple[PeptidePopulationStats, tuple[str, ...]]:
    """Digest protein sequences and return summary plus exact peptide occurrences."""
    protein_sequences = tuple(sequences)
    peptides = tuple(
        peptide.sequence
        for sequence in protein_sequences
        for peptide in digest_sequence(sequence, digestion)
    )
    lengths = _length_deciles(peptides)
    return (
        PeptidePopulationStats(
            proteins=len(protein_sequences),
            peptides=len(peptides),
            unique_peptides=len(set(peptides)),
            length_p10=lengths[0],
            length_median=lengths[1],
            length_p90=lengths[2],
            repeated_peptides=sum(1 for count in Counter(peptides).values() if count > 1),
        ),
        peptides,
    )


def compare_decoy_population(
    *,
    method: str,
    target_peptides: tuple[str, ...],
    decoy_sequences: Iterable[str],
    digestion: Digestion,
    initial_collisions: int,
    unresolved_collisions: int,
    dropped_peptides: int,
    omitted_decoys: int,
) -> DecoyMethodStats:
    """Measure one generated decoy population against fixed target peptides."""
    population, decoy_peptides = peptide_population(decoy_sequences, digestion)
    target_unique = set(target_peptides)
    decoy_unique = set(decoy_peptides)
    target_compositions = {_composition(peptide) for peptide in target_unique}
    matched_compositions = sum(
        _composition(peptide) in target_compositions for peptide in decoy_unique
    )
    return DecoyMethodStats(
        method=method,
        population=population,
        unique_ratio=(population.unique_peptides / len(target_unique) if target_unique else 0.0),
        shared_with_targets=len(decoy_unique & target_unique),
        composition_overlap=(matched_compositions / len(decoy_unique) if decoy_unique else 0.0),
        initial_collisions=initial_collisions,
        unresolved_collisions=unresolved_collisions,
        dropped_peptides=dropped_peptides,
        omitted_decoys=omitted_decoys,
    )


def _composition(peptide: str) -> str:
    return "".join(sorted(peptide))


def _length_deciles(peptides: tuple[str, ...]) -> tuple[float, float, float]:
    if not peptides:
        return 0.0, 0.0, 0.0
    lengths = [len(peptide) for peptide in peptides]
    if len(lengths) < _DECILE_COUNT:
        return float(min(lengths)), float(statistics.median(lengths)), float(max(lengths))
    deciles = statistics.quantiles(lengths, n=_DECILE_COUNT)
    return deciles[0], deciles[4], deciles[8]
