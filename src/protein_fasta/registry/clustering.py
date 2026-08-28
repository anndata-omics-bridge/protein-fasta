"""Registry queries that feed backend-free database clustering."""

from __future__ import annotations

from collections.abc import Collection

from protein_fasta.analytics.clustering import (
    ClusteringMetric,
    DatabaseClustering,
    cluster_comparisons,
)
from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.export import query_similarity_data

NEIGHBOURHOOD_LIMIT = 50
"""Neighbours shown alongside the selected database in the clustering panel."""


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
