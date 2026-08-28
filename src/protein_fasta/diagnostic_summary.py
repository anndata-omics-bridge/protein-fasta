"""Aggregate record-level protein FASTA diagnostics."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field

from protein_fasta.record import ProteinDiagnostics
from protein_fasta.summary import FastaSummary, SummaryAccumulator


def _empty_string_counts() -> dict[str, int]:
    return {}


def _empty_classification_counts() -> dict[tuple[str, ...], int]:
    return {}


@dataclass(frozen=True, slots=True)
class ProteinDiagnosticsSummary:
    """Aggregate facts from a stream of record-level diagnostics."""

    proteins: FastaSummary
    namespace_counts: dict[str, int] = field(default_factory=_empty_string_counts)
    classification_counts: dict[str, int] = field(default_factory=_empty_string_counts)
    classification_combination_counts: dict[tuple[str, ...], int] = field(
        default_factory=_empty_classification_counts
    )
    upper_cased_count: int = 0
    stop_stripped_count: int = 0
    illegal_sequence_count: int = 0
    illegal_residue_counts: dict[str, int] = field(default_factory=_empty_string_counts)


def summarize_protein_diagnostics(
    records: Iterable[ProteinDiagnostics],
    /,
) -> ProteinDiagnosticsSummary:
    """Summarize record diagnostics while preserving overlapping classifications."""
    sequences = SummaryAccumulator()
    namespaces: Counter[str] = Counter()
    classifications: Counter[str] = Counter()
    combinations: Counter[tuple[str, ...]] = Counter()
    illegal_residues: Counter[str] = Counter()
    upper_cased_count = 0
    stop_stripped_count = 0
    illegal_sequence_count = 0

    for record in records:
        sequences.add(record.protein.sequence)
        namespaces[record.identifier_namespace] += 1
        classifications.update(record.classifications)
        combinations[tuple(sorted(record.classifications))] += 1
        upper_cased_count += record.upper_cased
        stop_stripped_count += record.stop_stripped
        if record.illegal_residues:
            illegal_sequence_count += 1
            illegal_residues.update(record.illegal_residues)

    return ProteinDiagnosticsSummary(
        proteins=sequences.summary(),
        namespace_counts=dict(sorted(namespaces.items())),
        classification_counts=dict(sorted(classifications.items())),
        classification_combination_counts=dict(sorted(combinations.items())),
        upper_cased_count=upper_cased_count,
        stop_stripped_count=stop_stripped_count,
        illegal_sequence_count=illegal_sequence_count,
        illegal_residue_counts=dict(sorted(illegal_residues.items())),
    )
