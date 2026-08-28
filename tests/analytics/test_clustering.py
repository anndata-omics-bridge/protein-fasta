from __future__ import annotations

from itertools import combinations

import pytest

from protein_fasta.analytics.clustering import ClusteringMetric, ClusterMerge, cluster_comparisons
from protein_fasta.analytics.comparisons import (
    KindCounts,
    KindSimilarity,
    SimilarityData,
    SimilarityPair,
)
from protein_fasta.registry import clustering as registry_clustering
from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.clustering import cluster_registered_databases
from tests.registry_helpers import execute_script, in_memory_registry


def _kind_similarity(
    shared_ids: int,
    *,
    count_a: int = 10,
    count_b: int = 10,
    shared_sequences: int | None = None,
) -> KindSimilarity:
    counts_a = KindCounts(
        entries=count_a,
        distinct_ids=count_a,
        distinct_sequences=count_a,
        distinct_pairs=count_a,
    )
    counts_b = KindCounts(
        entries=count_b,
        distinct_ids=count_b,
        distinct_sequences=count_b,
        distinct_pairs=count_b,
    )
    sequence_overlap = shared_ids if shared_sequences is None else shared_sequences
    return KindSimilarity(
        a=counts_a,
        b=counts_b,
        shared_ids=shared_ids,
        shared_sequences=sequence_overlap,
        shared_exact_pairs=min(shared_ids, sequence_overlap),
        matching_shared_ids=min(shared_ids, sequence_overlap),
    )


def _pair(
    path_a: str,
    path_b: str,
    shared_ids: int,
    *,
    count_a: int = 10,
    count_b: int = 10,
    shared_sequences: int | None = None,
) -> SimilarityPair:
    return SimilarityPair(
        database_a_relative_path=path_a,
        database_b_relative_path=path_b,
        target=_kind_similarity(
            shared_ids,
            count_a=count_a,
            count_b=count_b,
            shared_sequences=shared_sequences,
        ),
        contaminant=_kind_similarity(0),
    )


@pytest.mark.parametrize(
    ("paths", "expected_order"),
    [
        ((), ()),
        (("only.fasta",), (0,)),
    ],
)
def testcluster_comparisons_handles_zero_and_one_database(
    paths: tuple[str, ...],
    expected_order: tuple[int, ...],
) -> None:
    result = cluster_comparisons(
        SimilarityData(relative_paths=paths, pairs=()),
        ClusteringMetric.TARGET_IDS,
        dict.fromkeys(paths, 1),
    )

    assert result.metric is ClusteringMetric.TARGET_IDS
    assert result.relative_paths == paths
    assert result.excluded_empty_paths == ()
    assert result.leaf_order == expected_order
    assert result.merges == ()


def testcluster_comparisons_handles_two_databases() -> None:
    data = SimilarityData(
        relative_paths=("b.fasta", "a.fasta"),
        pairs=(_pair("b.fasta", "a.fasta", 10),),
    )

    result = cluster_comparisons(data, ClusteringMetric.TARGET_IDS)

    assert result.relative_paths == ("a.fasta", "b.fasta")
    assert result.ordered_relative_paths == ("a.fasta", "b.fasta")
    assert result.merges == (
        ClusterMerge(
            cluster_id=2,
            left_id=0,
            right_id=1,
            distance=0.0,
            leaf_count=2,
        ),
    )


def testcluster_comparisons_uses_weighted_average_linkage() -> None:
    data = SimilarityData(
        relative_paths=("c.fasta", "a.fasta", "b.fasta"),
        pairs=(
            _pair("b.fasta", "c.fasta", 2),
            _pair("a.fasta", "c.fasta", 4),
            _pair("a.fasta", "b.fasta", 9),
        ),
    )

    result = cluster_comparisons(data, ClusteringMetric.TARGET_IDS)

    distance_ab = 1.0 - 9 / 11
    distance_ac = 1.0 - 4 / 16
    distance_bc = 1.0 - 2 / 18
    assert result.merges[0] == ClusterMerge(
        cluster_id=3,
        left_id=0,
        right_id=1,
        distance=result.merges[0].distance,
        leaf_count=2,
    )
    assert result.merges[0].distance == pytest.approx(distance_ab)
    assert result.merges[1] == ClusterMerge(
        cluster_id=4,
        left_id=3,
        right_id=2,
        distance=result.merges[1].distance,
        leaf_count=3,
    )
    assert result.merges[1].distance == pytest.approx((distance_ac + distance_bc) / 2)


def testcluster_comparisons_resolves_ties_by_relative_path() -> None:
    paths = ("d.fasta", "b.fasta", "a.fasta", "c.fasta")
    pairs = tuple(_pair(left, right, 0) for left, right in combinations(reversed(paths), 2))

    result = cluster_comparisons(
        SimilarityData(relative_paths=paths, pairs=tuple(reversed(pairs))),
        ClusteringMetric.TARGET_IDS,
    )

    assert result.relative_paths == tuple(sorted(paths))
    assert result.ordered_relative_paths == tuple(sorted(paths))
    assert [(merge.left_id, merge.right_id) for merge in result.merges] == [
        (0, 1),
        (4, 2),
        (5, 3),
    ]


def testcluster_comparisons_supports_target_sequence_jaccard() -> None:
    paths = ("a.fasta", "b.fasta", "c.fasta")
    data = SimilarityData(
        relative_paths=paths,
        pairs=(
            _pair("a.fasta", "b.fasta", 10, shared_sequences=0),
            _pair("a.fasta", "c.fasta", 0, shared_sequences=10),
            _pair("b.fasta", "c.fasta", 0, shared_sequences=0),
        ),
    )

    by_ids = cluster_comparisons(data, ClusteringMetric.TARGET_IDS)
    by_sequences = cluster_comparisons(
        data,
        ClusteringMetric.TARGET_SEQUENCES,
    )

    assert (by_ids.merges[0].left_id, by_ids.merges[0].right_id) == (0, 1)
    assert by_sequences.metric is ClusteringMetric.TARGET_SEQUENCES
    assert (by_sequences.merges[0].left_id, by_sequences.merges[0].right_id) == (0, 2)


def testcluster_comparisons_excludes_empty_metric_sets() -> None:
    data = SimilarityData(
        relative_paths=("empty.fasta", "populated.fasta"),
        pairs=(
            _pair(
                "empty.fasta",
                "populated.fasta",
                0,
                count_a=0,
                count_b=10,
                shared_sequences=0,
            ),
        ),
    )

    result = cluster_comparisons(data, ClusteringMetric.TARGET_IDS)

    assert result.relative_paths == ("populated.fasta",)
    assert result.excluded_empty_paths == ("empty.fasta",)
    assert result.leaf_order == (0,)
    assert result.merges == ()


@pytest.mark.parametrize(
    "data",
    [
        SimilarityData(relative_paths=("a", "b"), pairs=()),
        SimilarityData(
            relative_paths=("a", "b"),
            pairs=(_pair("a", "b", 1), _pair("b", "a", 1)),
        ),
        SimilarityData(
            relative_paths=("a", "b"),
            pairs=(_pair("a", "unknown", 1),),
        ),
    ],
)
def testcluster_comparisons_rejects_inconsistent_pair_data(
    data: SimilarityData,
) -> None:
    with pytest.raises(ValueError):
        cluster_comparisons(data, ClusteringMetric.TARGET_IDS)


def test_cluster_registered_databases_uses_materialized_tables_only() -> None:
    connection = in_memory_registry()
    execute_script(
        connection,
        """
        CREATE TABLE databases (
            id INTEGER PRIMARY KEY,
            relative_path TEXT NOT NULL,
            detail_level TEXT NOT NULL
        );
        CREATE TABLE database_kind_stats (
            database_id INTEGER NOT NULL,
            kind TEXT NOT NULL,
            entry_count INTEGER NOT NULL,
            distinct_ids INTEGER NOT NULL,
            distinct_sequences INTEGER NOT NULL,
            distinct_descriptions INTEGER NOT NULL,
            distinct_pairs INTEGER NOT NULL
        );
        CREATE TABLE database_pair_stats (
            database_id_low INTEGER NOT NULL,
            database_id_high INTEGER NOT NULL,
            kind TEXT NOT NULL,
            shared_ids INTEGER NOT NULL,
            shared_sequence_checksums INTEGER NOT NULL,
            shared_descriptions INTEGER NOT NULL,
            shared_exact_pairs INTEGER NOT NULL,
            matching_shared_ids INTEGER NOT NULL
        );
        INSERT INTO databases VALUES
            (1, 'a.fasta', 'full'),
            (2, 'b.fasta', 'full'),
            (3, 'large.fasta', 'metadata_only');
        INSERT INTO database_kind_stats VALUES
            (1, 'target', 2, 2, 2, 2, 2),
            (1, 'contaminant', 0, 0, 0, 0, 0),
            (2, 'target', 2, 2, 2, 2, 2),
            (2, 'contaminant', 0, 0, 0, 0, 0);
        INSERT INTO database_pair_stats VALUES
            (1, 2, 'target', 1, 1, 1, 1, 1),
            (1, 2, 'contaminant', 0, 0, 0, 0, 0);
        """,
    )

    result = cluster_registered_databases(connection)

    assert result.relative_paths == ("a.fasta", "b.fasta")
    assert result.omitted_metadata_paths == ("large.fasta",)
    assert result.merges[0].distance == pytest.approx(2 / 3)
    connection.close()


def test_cluster_registered_databases_excludes_a_single_empty_database() -> None:
    connection = in_memory_registry()
    execute_script(
        connection,
        """
        CREATE TABLE databases (
            id INTEGER PRIMARY KEY,
            relative_path TEXT NOT NULL,
            detail_level TEXT NOT NULL
        );
        CREATE TABLE database_kind_stats (
            database_id INTEGER NOT NULL,
            kind TEXT NOT NULL,
            entry_count INTEGER NOT NULL,
            distinct_ids INTEGER NOT NULL,
            distinct_sequences INTEGER NOT NULL,
            distinct_descriptions INTEGER NOT NULL,
            distinct_pairs INTEGER NOT NULL
        );
        CREATE TABLE database_pair_stats (
            database_id_low INTEGER NOT NULL,
            database_id_high INTEGER NOT NULL,
            kind TEXT NOT NULL,
            shared_ids INTEGER NOT NULL,
            shared_sequence_checksums INTEGER NOT NULL,
            shared_descriptions INTEGER NOT NULL,
            shared_exact_pairs INTEGER NOT NULL,
            matching_shared_ids INTEGER NOT NULL
        );
        INSERT INTO databases VALUES (1, 'empty.fasta', 'full');
        INSERT INTO database_kind_stats VALUES
            (1, 'target', 0, 0, 0, 0, 0),
            (1, 'contaminant', 0, 0, 0, 0, 0);
        """,
    )

    result = cluster_registered_databases(connection)

    assert result.relative_paths == ()
    assert result.excluded_empty_paths == ("empty.fasta",)
    assert result.leaf_order == ()
    assert result.merges == ()
    connection.close()


def _bounded_registry(database_count: int) -> RegistryConnection:
    """Registry with one hub plus ``database_count`` decreasingly similar databases."""
    connection = in_memory_registry()
    execute_script(
        connection,
        """
        CREATE TABLE databases (
            id INTEGER PRIMARY KEY,
            relative_path TEXT NOT NULL,
            detail_level TEXT NOT NULL
        );
        CREATE TABLE database_kind_stats (
            database_id INTEGER NOT NULL,
            kind TEXT NOT NULL,
            entry_count INTEGER NOT NULL,
            distinct_ids INTEGER NOT NULL,
            distinct_sequences INTEGER NOT NULL,
            distinct_descriptions INTEGER NOT NULL,
            distinct_pairs INTEGER NOT NULL
        );
        CREATE TABLE database_pair_stats (
            database_id_low INTEGER NOT NULL,
            database_id_high INTEGER NOT NULL,
            kind TEXT NOT NULL,
            shared_ids INTEGER NOT NULL,
            shared_sequence_checksums INTEGER NOT NULL,
            shared_descriptions INTEGER NOT NULL,
            shared_exact_pairs INTEGER NOT NULL,
            matching_shared_ids INTEGER NOT NULL
        );
        """,
    )
    for database_id in range(1, database_count + 1):
        connection.execute(
            "INSERT INTO databases VALUES (?, ?, 'full')",
            (database_id, f"db{database_id:04d}.fasta"),
        )
        connection.execute(
            "INSERT INTO database_kind_stats VALUES (?, 'target', 100, 100, 100, 100, 100)",
            (database_id,),
        )
        connection.execute(
            "INSERT INTO database_kind_stats VALUES (?, 'contaminant', 0, 0, 0, 0, 0)",
            (database_id,),
        )
    for low, high in combinations(range(1, database_count + 1), 2):
        # Database 1 is the hub: similarity to it decreases as the id grows, so
        # the expected neighbourhood is deterministic. All other pairs overlap
        # slightly, which keeps the pair cache complete.
        shared = max(0, 100 - high) if low == 1 else 1
        connection.execute(
            "INSERT INTO database_pair_stats VALUES (?, ?, 'target', ?, ?, 0, 0, 0)",
            (low, high, shared, shared),
        )
        connection.execute(
            "INSERT INTO database_pair_stats VALUES (?, ?, 'contaminant', 0, 0, 0, 0, 0)",
            (low, high),
        )
    return connection


def test_select_similar_database_ids_ranks_by_jaccard_and_bounds() -> None:
    connection = _bounded_registry(40)

    ids = registry_clustering.select_similar_database_ids(connection, "db0001.fasta", limit=5)

    # The hub itself first, then ids 2..6 — the most similar by construction.
    assert ids == (1, 2, 3, 4, 5, 6)
    connection.close()


def test_cluster_similar_databases_clusters_only_the_neighbourhood() -> None:
    connection = _bounded_registry(40)

    result = registry_clustering.cluster_similar_databases(connection, "db0001.fasta", limit=5)

    assert result.relative_paths == tuple(f"db{index:04d}.fasta" for index in range(1, 7))
    # 6 leaves always produce 5 merges, regardless of how large the registry is.
    assert len(result.merges) == 5
    connection.close()


def test_cluster_similar_databases_returns_empty_for_unknown_selection() -> None:
    connection = _bounded_registry(5)

    result = registry_clustering.cluster_similar_databases(connection, "missing.fasta")

    assert result.relative_paths == ()
    assert result.merges == ()
    connection.close()
