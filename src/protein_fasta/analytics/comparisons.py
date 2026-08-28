"""Backend-free database comparison values and set metrics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class KindCounts:
    """Set denominators for one entry kind in one database."""

    entries: int
    distinct_ids: int
    distinct_sequences: int
    distinct_pairs: int
    distinct_descriptions: int = 0


@dataclass(frozen=True, slots=True)
class KindSimilarity:
    """One materialized kind comparison and its per-side denominators."""

    a: KindCounts
    b: KindCounts
    shared_ids: int
    shared_sequences: int
    shared_exact_pairs: int
    matching_shared_ids: int
    shared_descriptions: int = 0

    @property
    def id_jaccard(self) -> float:
        """Return Jaccard similarity between the two distinct ID sets."""
        return _jaccard(self.shared_ids, self.a.distinct_ids, self.b.distinct_ids)

    @property
    def sequence_jaccard(self) -> float:
        """Return Jaccard similarity between distinct sequence-checksum sets."""
        return _jaccard(
            self.shared_sequences,
            self.a.distinct_sequences,
            self.b.distinct_sequences,
        )

    @property
    def description_jaccard(self) -> float:
        """Return Jaccard similarity between normalized-description hash sets."""
        return _jaccard(
            self.shared_descriptions,
            self.a.distinct_descriptions,
            self.b.distinct_descriptions,
        )

    @property
    def id_coverage_a_by_b(self) -> float:
        """Return the fraction of A's distinct IDs also present in B."""
        return _coverage(self.shared_ids, self.a.distinct_ids)

    @property
    def id_coverage_b_by_a(self) -> float:
        """Return the fraction of B's distinct IDs also present in A."""
        return _coverage(self.shared_ids, self.b.distinct_ids)

    @property
    def id_containment(self) -> float:
        """Return shared IDs divided by the smaller non-empty ID set."""
        smaller = min(self.a.distinct_ids, self.b.distinct_ids)
        return self.shared_ids / smaller if smaller else 0.0

    @property
    def sequence_coverage_a_by_b(self) -> float:
        """Return the fraction of A's distinct sequences also present in B."""
        return _coverage(self.shared_sequences, self.a.distinct_sequences)

    @property
    def sequence_coverage_b_by_a(self) -> float:
        """Return the fraction of B's distinct sequences also present in A."""
        return _coverage(self.shared_sequences, self.b.distinct_sequences)

    @property
    def sequence_containment(self) -> float:
        """Return shared sequences divided by the smaller non-empty set."""
        return _containment(
            self.shared_sequences, self.a.distinct_sequences, self.b.distinct_sequences
        )

    @property
    def description_coverage_a_by_b(self) -> float:
        """Return the fraction of A's distinct descriptions also present in B."""
        return _coverage(self.shared_descriptions, self.a.distinct_descriptions)

    @property
    def description_coverage_b_by_a(self) -> float:
        """Return the fraction of B's distinct descriptions also present in A."""
        return _coverage(self.shared_descriptions, self.b.distinct_descriptions)

    @property
    def description_containment(self) -> float:
        """Return shared descriptions divided by the smaller non-empty set."""
        return _containment(
            self.shared_descriptions,
            self.a.distinct_descriptions,
            self.b.distinct_descriptions,
        )

    @property
    def a_contained_in_b(self) -> bool:
        """Return whether A is a non-empty ID subset of B."""
        return bool(self.shared_ids) and self.shared_ids == self.a.distinct_ids

    @property
    def b_contained_in_a(self) -> bool:
        """Return whether B is a non-empty ID subset of A."""
        return bool(self.shared_ids) and self.shared_ids == self.b.distinct_ids

    @property
    def sequence_a_contained_in_b(self) -> bool:
        """Return whether A is a non-empty sequence subset of B."""
        return bool(self.shared_sequences) and self.shared_sequences == self.a.distinct_sequences

    @property
    def sequence_b_contained_in_a(self) -> bool:
        """Return whether B is a non-empty sequence subset of A."""
        return bool(self.shared_sequences) and self.shared_sequences == self.b.distinct_sequences

    @property
    def description_a_contained_in_b(self) -> bool:
        """Return whether A is a non-empty description subset of B."""
        return (
            bool(self.shared_descriptions)
            and self.shared_descriptions == self.a.distinct_descriptions
        )

    @property
    def description_b_contained_in_a(self) -> bool:
        """Return whether B is a non-empty description subset of A."""
        return (
            bool(self.shared_descriptions)
            and self.shared_descriptions == self.b.distinct_descriptions
        )

    @property
    def exact_id_set(self) -> bool:
        """Return whether both databases have the same non-empty ID set."""
        return self.a_contained_in_b and self.b_contained_in_a

    @property
    def exact_sequence_set(self) -> bool:
        """Return whether both databases have the same non-empty sequence set."""
        return self.sequence_a_contained_in_b and self.sequence_b_contained_in_a

    @property
    def exact_description_set(self) -> bool:
        """Return whether both databases have the same non-empty description set."""
        return self.description_a_contained_in_b and self.description_b_contained_in_a

    @property
    def changed_shared_ids(self) -> int:
        """Return shared IDs without any exact sequence match."""
        return self.shared_ids - self.matching_shared_ids

    @property
    def exact_content(self) -> bool:
        """Return whether both databases have the same non-empty ID/sequence pairs."""
        return (
            bool(self.shared_exact_pairs)
            and self.a.distinct_pairs == self.shared_exact_pairs
            and self.b.distinct_pairs == self.shared_exact_pairs
        )


@dataclass(frozen=True, slots=True)
class SimilarityPair:
    """Long-form target and contaminant statistics for one database pair."""

    database_a_relative_path: str
    database_b_relative_path: str
    target: KindSimilarity
    contaminant: KindSimilarity


@dataclass(frozen=True, slots=True)
class SimilarityData:
    """Stable path order and every unordered pair in that order."""

    relative_paths: tuple[str, ...]
    pairs: tuple[SimilarityPair, ...]
    omitted_relative_paths: tuple[str, ...] = ()


def _coverage(shared: int, total: int) -> float:
    return shared / total if total else 0.0


def _containment(shared: int, count_a: int, count_b: int) -> float:
    smaller = min(count_a, count_b)
    return shared / smaller if smaller else 0.0


def _jaccard(shared: int, count_a: int, count_b: int) -> float:
    union = count_a + count_b - shared
    return shared / union if union else 0.0
