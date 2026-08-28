"""Backend-free deterministic clustering of database similarities."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from itertools import combinations

from protein_fasta.analytics.comparisons import SimilarityData, SimilarityPair


class ClusteringMetric(StrEnum):
    """Supported kind-scoped similarities for database clustering."""

    TARGET_IDS = "target_ids"
    TARGET_SEQUENCES = "target_sequences"


@dataclass(frozen=True, slots=True)
class ClusterMerge:
    """One average-linkage merge using SciPy-compatible cluster identifiers."""

    cluster_id: int
    left_id: int
    right_id: int
    distance: float
    leaf_count: int


@dataclass(frozen=True, slots=True)
class DatabaseClustering:
    """Average-linkage tree for one registered-database target metric."""

    metric: ClusteringMetric
    relative_paths: tuple[str, ...]
    excluded_empty_paths: tuple[str, ...]
    leaf_order: tuple[int, ...]
    merges: tuple[ClusterMerge, ...]
    omitted_metadata_paths: tuple[str, ...] = ()

    @property
    def ordered_relative_paths(self) -> tuple[str, ...]:
        """Return database paths in deterministic dendrogram leaf order."""
        return tuple(self.relative_paths[index] for index in self.leaf_order)


@dataclass(frozen=True, slots=True)
class _WorkingCluster:
    members: tuple[int, ...]
    leaf_order: tuple[int, ...]


def _pair_key(left_id: int, right_id: int) -> tuple[int, int]:
    return (left_id, right_id) if left_id < right_id else (right_id, left_id)


def _metric_value_and_counts(
    pair: SimilarityPair,
    metric: ClusteringMetric,
) -> tuple[float, int, int]:
    if metric is ClusteringMetric.TARGET_IDS:
        return (
            pair.target.id_jaccard,
            pair.target.a.distinct_ids,
            pair.target.b.distinct_ids,
        )
    return (
        pair.target.sequence_jaccard,
        pair.target.a.distinct_sequences,
        pair.target.b.distinct_sequences,
    )


def _record_count(counts: dict[str, int], path: str, count: int) -> None:
    previous = counts.setdefault(path, count)
    if previous != count:
        raise ValueError(f"Similarity data has inconsistent set sizes for {path!r}.")


def _target_distances(
    data: SimilarityData,
    metric: ClusteringMetric,
    known_counts: dict[str, int] | None = None,
) -> tuple[tuple[str, ...], tuple[str, ...], dict[tuple[int, int], float]]:
    all_paths = tuple(sorted(data.relative_paths))
    if len(all_paths) != len(set(all_paths)):
        raise ValueError("Database relative paths must be unique for clustering.")

    path_set = set(all_paths)
    similarities: dict[tuple[str, str], float] = {}
    observed_counts: dict[str, int] = {}
    for pair in data.pairs:
        path_a = pair.database_a_relative_path
        path_b = pair.database_b_relative_path
        if path_a not in path_set or path_b not in path_set:
            raise ValueError("Similarity pairs must reference registered database paths.")
        if path_a == path_b:
            raise ValueError("Similarity pairs must reference two different databases.")

        similarity, count_a, count_b = _metric_value_and_counts(pair, metric)
        _record_count(observed_counts, path_a, count_a)
        _record_count(observed_counts, path_b, count_b)
        if not math.isfinite(similarity) or not 0.0 <= similarity <= 1.0:
            raise ValueError("Target Jaccard similarities must be between zero and one.")
        key = (path_a, path_b) if path_a < path_b else (path_b, path_a)
        if key in similarities:
            raise ValueError("Each database pair must have exactly one similarity value.")
        similarities[key] = similarity

    expected_pairs = len(all_paths) * (len(all_paths) - 1) // 2
    if len(similarities) != expected_pairs:
        raise ValueError("Target similarity data is incomplete for clustering.")
    counts = known_counts if known_counts is not None else observed_counts
    if counts and set(counts) != path_set:
        raise ValueError("Target set sizes are incomplete for clustering.")

    excluded_paths = tuple(path for path in all_paths if counts.get(path, 1) == 0)
    excluded_set = set(excluded_paths)
    paths = tuple(path for path in all_paths if path not in excluded_set)
    path_indexes = {path: index for index, path in enumerate(paths)}
    distances = {
        _pair_key(path_indexes[path_a], path_indexes[path_b]): 1.0 - similarity
        for (path_a, path_b), similarity in similarities.items()
        if path_a in path_indexes and path_b in path_indexes
    }
    return paths, excluded_paths, distances


def _cluster_key(cluster: _WorkingCluster, paths: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(paths[index] for index in cluster.members)


def cluster_comparisons(
    data: SimilarityData,
    metric: ClusteringMetric,
    known_counts: dict[str, int] | None = None,
) -> DatabaseClustering:
    paths, excluded_paths, distances = _target_distances(data, metric, known_counts)
    database_count = len(paths)
    if database_count == 0:
        return DatabaseClustering(
            metric=metric,
            relative_paths=paths,
            excluded_empty_paths=excluded_paths,
            leaf_order=(),
            merges=(),
            omitted_metadata_paths=data.omitted_relative_paths,
        )
    if database_count == 1:
        return DatabaseClustering(
            metric=metric,
            relative_paths=paths,
            excluded_empty_paths=excluded_paths,
            leaf_order=(0,),
            merges=(),
            omitted_metadata_paths=data.omitted_relative_paths,
        )

    active = {
        index: _WorkingCluster(members=(index,), leaf_order=(index,))
        for index in range(database_count)
    }
    merges: list[ClusterMerge] = []

    while len(active) > 1:
        ordered_ids = sorted(active, key=lambda cluster_id: _cluster_key(active[cluster_id], paths))
        candidates = (
            (
                distances[_pair_key(left_id, right_id)],
                _cluster_key(active[left_id], paths),
                _cluster_key(active[right_id], paths),
                left_id,
                right_id,
            )
            for left_id, right_id in combinations(ordered_ids, 2)
        )
        distance, _, _, left_id, right_id = min(candidates)
        left = active.pop(left_id)
        right = active.pop(right_id)
        cluster_id = database_count + len(merges)
        merged = _WorkingCluster(
            members=tuple(sorted((*left.members, *right.members))),
            leaf_order=(*left.leaf_order, *right.leaf_order),
        )

        for other_id in active:
            left_distance = distances[_pair_key(left_id, other_id)]
            right_distance = distances[_pair_key(right_id, other_id)]
            distances[_pair_key(cluster_id, other_id)] = (
                len(left.members) * left_distance + len(right.members) * right_distance
            ) / len(merged.members)

        active[cluster_id] = merged
        merges.append(
            ClusterMerge(
                cluster_id=cluster_id,
                left_id=left_id,
                right_id=right_id,
                distance=distance,
                leaf_count=len(merged.members),
            )
        )

    root = next(iter(active.values()))
    return DatabaseClustering(
        metric=metric,
        relative_paths=paths,
        excluded_empty_paths=excluded_paths,
        leaf_order=root.leaf_order,
        merges=tuple(merges),
        omitted_metadata_paths=data.omitted_relative_paths,
    )
