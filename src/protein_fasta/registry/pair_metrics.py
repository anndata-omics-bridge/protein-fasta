"""Query exact pair metrics from a registry connection."""

from __future__ import annotations

from dataclasses import dataclass

from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.kinds import EntryKind


@dataclass(frozen=True, slots=True)
class PairMetricCounts:
    """Raw symmetric intersections between a selection and registered databases."""

    shared_ids: dict[int, int]
    shared_sequences: dict[int, int]
    shared_descriptions: dict[int, int]
    shared_pairs: dict[int, int]
    matching_ids: dict[int, int]


@dataclass(frozen=True, slots=True)
class PairMetricSelection:
    """Trusted SQL relations describing one already-scoped entry selection."""

    ids_table: str
    sequences_table: str
    descriptions_table: str
    pairs_table: str
    where: str
    params: tuple[object, ...]
    kind_filtered: bool = False


def entry_kind_sql_literal(kind: EntryKind) -> str:
    """Return a trusted SQL literal so an engine can use kind-specific indexes."""
    return "'" + kind.value.replace("'", "''") + "'"


def pair_metric_counts(
    connection: RegistryConnection,
    *,
    selection: PairMetricSelection,
    kind: EntryKind,
    excluded_database_id: int | None,
    minimum_other_database_id: int | None = None,
    other_database_table: str | None = None,
) -> PairMetricCounts:
    """Return raw pair intersections for one trusted internal SQL selection."""
    if excluded_database_id is not None and minimum_other_database_id is not None:
        raise ValueError(
            "Pair metrics cannot exclude one database and apply a minimum database ID together."
        )
    kind_literal = entry_kind_sql_literal(kind)
    selected_kind_filter = "" if selection.kind_filtered else f"AND kind = {kind_literal}"
    query_params = selection.params
    if excluded_database_id is not None:
        other_filter = "WHERE other.database_id != ?"
        query_params = (*query_params, excluded_database_id)
    elif minimum_other_database_id is not None:
        other_filter = "WHERE other.database_id > ?"
        query_params = (*query_params, minimum_other_database_id)
    else:
        other_filter = ""
    candidate_join = (
        ""
        if other_database_table is None
        else f"JOIN {other_database_table} AS candidate ON candidate.database_id = other.database_id"
    )

    shared_ids = {
        int(row[0]): int(row[1])
        for row in connection.execute(
            f"""
            SELECT other.database_id, COUNT(DISTINCT other.sequence_id)
            FROM (
                SELECT DISTINCT sequence_id
                FROM {selection.ids_table}
                WHERE {selection.where} {selected_kind_filter}
            ) AS selected
            JOIN entries AS other{connection.index_hint("entries_kind_id_db")}
              ON other.sequence_id = selected.sequence_id
             AND other.kind = {kind_literal}
            {candidate_join}
            {other_filter}
            GROUP BY other.database_id
            """,
            query_params,
        )
    }
    shared_sequences = {
        int(row[0]): int(row[1])
        for row in connection.execute(
            f"""
            SELECT other.database_id, COUNT(DISTINCT other.sequence_hash)
            FROM (
                SELECT DISTINCT sequence_hash
                FROM {selection.sequences_table}
                WHERE {selection.where} {selected_kind_filter}
            ) AS selected
            JOIN entries AS other{connection.index_hint("entries_kind_sequence_db")}
              ON other.sequence_hash = selected.sequence_hash
             AND other.kind = {kind_literal}
            {candidate_join}
            {other_filter}
            GROUP BY other.database_id
            """,
            query_params,
        )
    }
    shared_descriptions = {
        int(row[0]): int(row[1])
        for row in connection.execute(
            f"""
            SELECT other.database_id, COUNT(DISTINCT other.description_hash)
            FROM (
                SELECT DISTINCT description_hash
                FROM {selection.descriptions_table}
                WHERE {selection.where} {selected_kind_filter}
                  AND description_hash IS NOT NULL
            ) AS selected
            JOIN entries AS other{connection.index_hint("entries_kind_description_db")}
              ON other.description_hash = selected.description_hash
             AND other.kind = {kind_literal}
            {candidate_join}
            {other_filter}
            GROUP BY other.database_id
            """,
            query_params,
        )
    }
    pair_rows = connection.execute(
        f"""
        SELECT matching.database_id,
               COUNT(*) AS shared_exact_pairs,
               COUNT(DISTINCT matching.sequence_id) AS matching_shared_ids
        FROM (
            SELECT DISTINCT other.database_id,
                            other.sequence_id,
                            other.sequence_hash
            FROM (
                SELECT DISTINCT sequence_id, sequence_hash
                FROM {selection.pairs_table}
                WHERE {selection.where} {selected_kind_filter}
            ) AS selected
            JOIN entries AS other{connection.index_hint("entries_kind_sequence_db")}
             ON other.sequence_id = selected.sequence_id
             AND other.sequence_hash = selected.sequence_hash
             AND other.kind = {kind_literal}
            {candidate_join}
            {other_filter}
        ) AS matching
        GROUP BY matching.database_id
        """,
        query_params,
    ).fetchall()
    return PairMetricCounts(
        shared_ids=shared_ids,
        shared_sequences=shared_sequences,
        shared_descriptions=shared_descriptions,
        shared_pairs={int(row[0]): int(row[1]) for row in pair_rows},
        matching_ids={int(row[0]): int(row[2]) for row in pair_rows},
    )
