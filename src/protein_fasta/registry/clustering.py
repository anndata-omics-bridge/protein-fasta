"""Registry queries that feed backend-free database clustering."""

from __future__ import annotations

from collections.abc import Collection

from protein_fasta.analytics.clustering import (
    ClusteringMetric,
    DatabaseClustering,
    cluster_comparisons,
)
from protein_fasta.analytics.comparisons import (
    KindCounts,
    KindSimilarity,
    SimilarityData,
    SimilarityPair,
)
from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.comparisons import DatabaseComparison
from protein_fasta.registry.export import query_similarity_data

NEIGHBOURHOOD_LIMIT = 50
"""Neighbours shown alongside the selected database in the clustering panel."""


def _comparison_metric(
    comparison: DatabaseComparison,
    metric: ClusteringMetric,
) -> float:
    if metric is ClusteringMetric.TARGET_IDS:
        return comparison.id_jaccard
    return comparison.sequence_jaccard


def _candidate_similarity_pair(
    candidate_label: str,
    other_path: str,
    comparison: DatabaseComparison,
) -> SimilarityPair:
    candidate_counts = KindCounts(
        entries=comparison.selected_ids,
        distinct_ids=comparison.selected_ids,
        distinct_sequences=comparison.selected_sequences,
        distinct_descriptions=comparison.selected_descriptions,
        distinct_pairs=0,
    )
    other_counts = KindCounts(
        entries=comparison.other_entries,
        distinct_ids=comparison.other_ids,
        distinct_sequences=comparison.other_sequences,
        distinct_descriptions=comparison.other_descriptions,
        distinct_pairs=0,
    )
    target = KindSimilarity(
        a=candidate_counts,
        b=other_counts,
        shared_ids=comparison.shared_ids,
        shared_sequences=comparison.shared_sequence_checksums,
        shared_descriptions=comparison.shared_descriptions,
        shared_exact_pairs=comparison.shared_exact_pairs,
        matching_shared_ids=comparison.shared_ids - comparison.changed_shared_ids,
    )
    empty_counts = KindCounts(
        entries=0,
        distinct_ids=0,
        distinct_sequences=0,
        distinct_descriptions=0,
        distinct_pairs=0,
    )
    contaminant = KindSimilarity(
        a=empty_counts,
        b=empty_counts,
        shared_ids=0,
        shared_sequences=0,
        shared_descriptions=0,
        shared_exact_pairs=0,
        matching_shared_ids=0,
    )
    return SimilarityPair(
        database_a_relative_path=candidate_label,
        database_b_relative_path=other_path,
        target=target,
        contaminant=contaminant,
    )


def cluster_candidate_with_similar_databases(
    connection: RegistryConnection,
    candidate_label: str,
    target_comparisons: Collection[DatabaseComparison],
    *,
    candidate_count: int,
    metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
    limit: int = NEIGHBOURHOOD_LIMIT,
) -> DatabaseClustering:
    """Cluster a transient candidate with its nearest registered neighbours."""
    comparison_by_id = {
        comparison.database.database_id: comparison for comparison in target_comparisons
    }
    database_ids = tuple(sorted(comparison_by_id))
    path_by_id: dict[int, str] = {}
    if database_ids:
        placeholders = ", ".join("?" for _ in database_ids)
        path_by_id = {
            int(row[0]): str(row[1])
            for row in connection.execute(
                f"""
                SELECT id, relative_path
                FROM databases
                WHERE detail_level = 'full'
                  AND id IN ({placeholders})
                """,
                database_ids,
            ).fetchall()
        }
    ranked = sorted(
        (
            comparison
            for database_id, comparison in comparison_by_id.items()
            if database_id in path_by_id
        ),
        key=lambda comparison: (
            -_comparison_metric(comparison, metric),
            path_by_id[comparison.database.database_id],
            comparison.database.database_id,
        ),
    )
    selected = ranked[:limit]
    selected_ids = tuple(comparison.database.database_id for comparison in selected)
    registered_data = query_similarity_data(connection, database_ids=selected_ids)
    if candidate_label in registered_data.relative_paths:
        raise ValueError(
            f"Candidate label {candidate_label!r} conflicts with a registered database path."
        )
    candidate_pairs = tuple(
        _candidate_similarity_pair(
            candidate_label,
            path_by_id[comparison.database.database_id],
            comparison,
        )
        for comparison in selected
    )
    data = SimilarityData(
        relative_paths=(*registered_data.relative_paths, candidate_label),
        pairs=(*registered_data.pairs, *candidate_pairs),
        omitted_relative_paths=registered_data.omitted_relative_paths,
    )
    counts = target_counts(connection, metric, selected_ids)
    counts[candidate_label] = candidate_count
    return cluster_comparisons(data, metric, counts)


def _target_count_column(metric: ClusteringMetric) -> str:
    return "distinct_ids" if metric is ClusteringMetric.TARGET_IDS else "distinct_sequences"


def target_counts(
    connection: RegistryConnection,
    metric: ClusteringMetric,
    database_ids: Collection[int] | None = None,
) -> dict[str, int]:
    predicate, id_params = ("", ())
    if database_ids is not None:
        ids = tuple(sorted(database_ids))
        if not ids:
            predicate = " AND 1 = 0"
        else:
            predicate = f" AND databases.id IN ({', '.join('?' for _ in ids)})"
            id_params = ids
    return {
        str(row[0]): int(row[1])
        for row in connection.execute(
            f"""
            SELECT databases.relative_path, stats.{_target_count_column(metric)}
            FROM databases
            JOIN database_kind_stats AS stats
              ON stats.database_id = databases.id
             AND stats.kind = 'target'
            WHERE databases.detail_level = 'full'
            {predicate}
            """,
            id_params,
        ).fetchall()
    }


def select_similar_database_ids(
    connection: RegistryConnection,
    selected_relative_path: str,
    *,
    metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
    limit: int = NEIGHBOURHOOD_LIMIT,
) -> tuple[int, ...]:
    """Return the selected database plus its ``limit`` nearest neighbours.

    This reads only the selected database's own row of the materialized pair
    table, ranks it by target Jaccard in SQL, and returns identifiers. Ties are
    broken by relative path so the neighbourhood is deterministic.
    """
    shared_column = (
        "shared_ids" if metric is ClusteringMetric.TARGET_IDS else "shared_sequence_checksums"
    )
    count_column = _target_count_column(metric)
    selected_row = connection.execute(
        f"""
        SELECT databases.id, stats.{count_column}
        FROM databases
        JOIN database_kind_stats AS stats
          ON stats.database_id = databases.id AND stats.kind = 'target'
        WHERE databases.relative_path = ? AND databases.detail_level = 'full'
        """,
        (selected_relative_path,),
    ).fetchone()
    if selected_row is None:
        return ()
    selected_id, selected_count = int(selected_row[0]), int(selected_row[1])

    # Two indexed lookups, one per position of the selected id in the pair key:
    # the (low, high, kind) primary key serves the first, the
    # database_pair_stats_high index the second. Matching the id with
    # `IN (database_id_low, database_id_high)` instead would defeat both and
    # scan the entire pair table.
    rows = connection.execute(
        f"""
        WITH neighbour_pairs AS (
            SELECT database_id_high AS other_id, {shared_column} AS shared
            FROM database_pair_stats
            WHERE database_id_low = ? AND kind = 'target'
            UNION ALL
            SELECT database_id_low AS other_id, {shared_column} AS shared
            FROM database_pair_stats
            WHERE database_id_high = ? AND kind = 'target'
        )
        SELECT neighbour_pairs.other_id,
               CAST(neighbour_pairs.shared AS REAL)
                 / (? + stats.{count_column} - neighbour_pairs.shared) AS jaccard
        FROM neighbour_pairs
        JOIN databases AS neighbour
          ON neighbour.id = neighbour_pairs.other_id AND neighbour.detail_level = 'full'
        JOIN database_kind_stats AS stats
          ON stats.database_id = neighbour_pairs.other_id AND stats.kind = 'target'
        WHERE ? + stats.{count_column} - neighbour_pairs.shared > 0
        ORDER BY jaccard DESC, neighbour.relative_path{connection.binary_collation()}, neighbour.id
        LIMIT ?
        """,
        (selected_id, selected_id, selected_count, selected_count, limit),
    ).fetchall()
    return (selected_id, *(int(row[0]) for row in rows))


def cluster_similar_databases(
    connection: RegistryConnection,
    selected_relative_path: str,
    *,
    metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
    limit: int = NEIGHBOURHOOD_LIMIT,
) -> DatabaseClustering:
    """Cluster the selected database together with its nearest neighbours.

    Bounded on purpose: pair loading and the merge loop are both quadratic or
    worse in the number of databases, and a dendrogram of the whole collection
    is unreadable anyway.
    """
    database_ids = select_similar_database_ids(
        connection, selected_relative_path, metric=metric, limit=limit
    )
    if not database_ids:
        return DatabaseClustering(
            metric=metric,
            relative_paths=(),
            excluded_empty_paths=(),
            leaf_order=(),
            merges=(),
            omitted_metadata_paths=(),
        )
    data = query_similarity_data(connection, database_ids=database_ids)
    return cluster_comparisons(data, metric, target_counts(connection, metric, database_ids))


def cluster_registered_databases(
    connection: RegistryConnection,
    *,
    metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
) -> DatabaseClustering:
    """Cluster every registered database using materialized target statistics."""
    data = query_similarity_data(connection)
    return cluster_comparisons(data, metric, target_counts(connection, metric))
