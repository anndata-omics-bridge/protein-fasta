"""Indexed target and contaminant comparisons between FASTA databases.

Target and contaminant evidence stays in separate, explicitly kind-scoped
results. Decoys and sentinels are packaging/metadata rather than database
similarity inputs. Target ID-set containment remains the duplicate-warning
score, while directional coverage, ID and sequence Jaccard similarity, and
exact ID/checksum pairs make subset and changed-sequence relationships explicit.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from protein_fasta.analytics.hashing import content_fingerprint, id_set_fingerprint
from protein_fasta.registry.backend.base import (
    RegistryConnection,
    RegistryIntegrityError,
    TempTableSpec,
)
from protein_fasta.registry.indexing import (
    CANDIDATE_TABLE,
    SCHEMA_VERSION,
    RegistrySchemaError,
    RegistrySettings,
    connect_registry,
    target_id_set,
)
from protein_fasta.registry.kinds import DetailLevel, EntryKind
from protein_fasta.registry.pair_metrics import (
    PairMetricSelection,
    entry_kind_sql_literal,
    pair_metric_counts,
)

_SELECTED_IDS_TABLE = "overlap_selected_ids"
_SELECTED_DESCRIPTIONS_TABLE = "overlap_selected_descriptions"
_SELECTED_PAIRS_TABLE = "overlap_selected_pairs"
_SELECTED_SEQUENCES_TABLE = "overlap_selected_sequences"
_COMPARABLE_KINDS = frozenset({EntryKind.TARGET, EntryKind.CONTAMINANT})

_SELECTED_TABLE_SPECS: tuple[TempTableSpec, ...] = (
    TempTableSpec(
        name=_SELECTED_IDS_TABLE,
        columns=(("sequence_id", "TEXT NOT NULL"),),
        primary_key=("sequence_id",),
        deduplicating=True,
    ),
    TempTableSpec(
        name=_SELECTED_PAIRS_TABLE,
        columns=(("sequence_id", "TEXT NOT NULL"), ("sequence_hash", "BLOB NOT NULL")),
        primary_key=("sequence_id", "sequence_hash"),
        deduplicating=True,
    ),
    TempTableSpec(
        name=_SELECTED_SEQUENCES_TABLE,
        columns=(("sequence_hash", "BLOB NOT NULL"),),
        primary_key=("sequence_hash",),
        deduplicating=True,
    ),
    TempTableSpec(
        name=_SELECTED_DESCRIPTIONS_TABLE,
        columns=(("description_hash", "BLOB NOT NULL"),),
        primary_key=("description_hash",),
        deduplicating=True,
    ),
)
"""The four staging tables one comparison selects into."""


class DatabaseRelationship(StrEnum):
    """Strongest kind-scoped ID/content relationship to another database."""

    EXACT_CONTENT = "exact_content"
    EXACT_ID_SET = "exact_id_set"
    SUBSET = "subset"
    SUPERSET = "superset"
    NEAR_MATCH = "near_match"
    PARTIAL_OVERLAP = "partial_overlap"
    NO_OVERLAP = "no_overlap"


@dataclass(frozen=True, slots=True)
class DatabaseReference:
    """Small registry reference returned with overlap results."""

    database_id: int
    filename: str
    annotation: str | None


@dataclass(frozen=True, slots=True)
class DatabaseComparison:
    """Aggregate comparison for one entry kind across two databases."""

    database: DatabaseReference
    kind: EntryKind
    selected_ids: int
    other_ids: int
    other_entries: int
    shared_ids: int
    selected_coverage: float
    other_coverage: float
    containment: float
    id_jaccard: float
    selected_sequences: int
    other_sequences: int
    shared_sequence_checksums: int
    sequence_jaccard: float
    shared_exact_pairs: int
    changed_shared_ids: int
    selected_only_ids: int
    other_only_ids: int
    exact_id_set: bool
    exact_content: bool
    relationship: DatabaseRelationship
    selected_sequence_coverage: float = 0.0
    other_sequence_coverage: float = 0.0
    sequence_containment: float = 0.0
    selected_descriptions: int = 0
    other_descriptions: int = 0
    shared_descriptions: int = 0
    selected_description_coverage: float = 0.0
    other_description_coverage: float = 0.0
    description_containment: float = 0.0
    description_jaccard: float = 0.0
    exact_sequence_set: bool = False
    exact_description_set: bool = False


@dataclass(frozen=True, slots=True)
class OverlapHit:
    """The best registry database at or above the containment threshold."""

    record: DatabaseReference
    score: float
    exact: bool


@dataclass(frozen=True, slots=True)
class _DatabaseRow:
    database: DatabaseReference
    entries: int
    ids: int
    sequences: int
    descriptions: int
    content_fingerprint: str


def containment(a: set[str], b: set[str]) -> float:
    """Return ``|a intersection b| / min(|a|, |b|)``.

    Empty sets return zero rather than treating two empty databases as a useful
    duplicate match.
    """
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def _containment_from_counts(shared: int, selected: int, other: int) -> float:
    smaller = min(selected, other)
    return shared / smaller if smaller else 0.0


def _coverage(shared: int, total: int) -> float:
    return shared / total if total else 0.0


def _jaccard(shared: int, selected: int, other: int) -> float:
    union = selected + other - shared
    return shared / union if union else 0.0


def _validate_comparison_kind(kind: EntryKind) -> None:
    if kind not in _COMPARABLE_KINDS:
        raise ValueError(
            f"Only target and contaminant entries can be compared; got {kind.value!r}."
        )


@contextmanager
def _registry_connection(
    source: Path | RegistryConnection,
) -> Generator[RegistryConnection]:
    """Yield an already-open registry, or open the one at a path and check it.

    Discriminated on ``Path`` rather than on the connection's concrete class: a
    protocol has no single class to test, and the caller either handed us a
    filename or something that already behaves like a registry.
    """
    if not isinstance(source, Path):
        yield source
        return
    with connect_registry(source) as connection:
        current_version = connection.schema_version()
        if current_version != SCHEMA_VERSION:
            raise RegistrySchemaError(current_version, path=source)
        yield connection


def _prepare_selected_tables(connection: RegistryConnection) -> None:
    """Replace the four staging tables with empty ones."""
    for spec in _SELECTED_TABLE_SPECS:
        connection.create_temp_table(spec)


def _stage_from_table(
    connection: RegistryConnection,
    table: str,
    where: str,
    kind: EntryKind,
) -> None:
    _prepare_selected_tables(connection)
    kind_literal = entry_kind_sql_literal(kind)
    connection.execute(
        f"""
        INSERT OR IGNORE INTO {connection.temp(_SELECTED_IDS_TABLE)} (sequence_id)
        SELECT sequence_id FROM {table} WHERE {where} AND kind = {kind_literal}
        """
    )
    connection.execute(
        f"""
        INSERT OR IGNORE INTO {connection.temp(_SELECTED_PAIRS_TABLE)}
            (sequence_id, sequence_hash)
        SELECT sequence_id, sequence_hash
        FROM {table}
        WHERE {where} AND kind = {kind_literal}
        """
    )
    connection.execute(
        f"""
        INSERT OR IGNORE INTO {connection.temp(_SELECTED_SEQUENCES_TABLE)} (sequence_hash)
        SELECT sequence_hash FROM {table} WHERE {where} AND kind = {kind_literal}
        """
    )
    connection.execute(
        f"""
        INSERT OR IGNORE INTO {connection.temp(_SELECTED_DESCRIPTIONS_TABLE)} (description_hash)
        SELECT description_hash
        FROM {table}
        WHERE {where} AND kind = {kind_literal} AND description_hash IS NOT NULL
        """
    )
    for selected_table in (
        _SELECTED_IDS_TABLE,
        _SELECTED_PAIRS_TABLE,
        _SELECTED_SEQUENCES_TABLE,
        _SELECTED_DESCRIPTIONS_TABLE,
    ):
        connection.analyze(connection.temp(selected_table))


def _stage_target_ids(connection: RegistryConnection, ids: Iterable[str]) -> None:
    connection.create_temp_table(_SELECTED_TABLE_SPECS[0])
    connection.executemany(
        f"INSERT OR IGNORE INTO {connection.temp(_SELECTED_IDS_TABLE)} (sequence_id) VALUES (?)",
        ((entry_id,) for entry_id in ids),
    )
    connection.analyze(connection.temp(_SELECTED_IDS_TABLE))


def _selected_id_count(connection: RegistryConnection) -> int:
    row = connection.execute(
        f"SELECT COUNT(*) FROM {connection.temp(_SELECTED_IDS_TABLE)}"
    ).fetchone()
    assert row is not None
    return int(row[0])


def _selected_sequence_count(connection: RegistryConnection) -> int:
    row = connection.execute(
        f"SELECT COUNT(*) FROM {connection.temp(_SELECTED_SEQUENCES_TABLE)}"
    ).fetchone()
    assert row is not None
    return int(row[0])


def _selected_description_count(connection: RegistryConnection) -> int:
    row = connection.execute(
        f"SELECT COUNT(*) FROM {connection.temp(_SELECTED_DESCRIPTIONS_TABLE)}"
    ).fetchone()
    assert row is not None
    return int(row[0])


def _selected_content_fingerprint(connection: RegistryConnection) -> str:
    rows = connection.execute(
        f"""
        SELECT sequence_id, sequence_hash
        FROM {connection.temp(_SELECTED_PAIRS_TABLE)}
        ORDER BY sequence_id, sequence_hash
        """
    )
    return content_fingerprint((str(row[0]), bytes(row[1])) for row in rows)


def _database_rows(
    connection: RegistryConnection,
    excluded_database_id: int | None,
    kind: EntryKind,
) -> list[_DatabaseRow]:
    rows = connection.execute(
        """
        SELECT databases.id, databases.filename, databases.annotation,
               stats.entry_count, stats.distinct_ids, stats.distinct_sequences,
               stats.distinct_descriptions, stats.content_fingerprint
        FROM databases
        LEFT JOIN database_kind_stats AS stats
          ON stats.database_id = databases.id
         AND stats.kind = ?
        WHERE databases.detail_level = ?
          AND (? IS NULL OR databases.id != ?)
        ORDER BY databases.filename COLLATE NOCASE, databases.id
        """,
        (kind.value, DetailLevel.FULL.value, excluded_database_id, excluded_database_id),
    ).fetchall()
    database_rows: list[_DatabaseRow] = []
    for row in rows:
        if any(row[index] is None for index in (3, 4, 5, 6, 7)):
            raise RegistryIntegrityError(
                f"kind statistics are incomplete for database id {int(row[0])}; run a full registry reindex"
            )
        database_rows.append(
            _DatabaseRow(
                database=DatabaseReference(
                    database_id=int(row[0]),
                    filename=str(row[1]),
                    annotation=None if row[2] is None else str(row[2]),
                ),
                entries=int(row[3]),
                ids=int(row[4]),
                sequences=int(row[5]),
                descriptions=int(row[6]),
                content_fingerprint=str(row[7]),
            )
        )
    return database_rows


def _relationship(
    *,
    exact_content: bool,
    exact_id_set: bool,
    shared: int,
    selected: int,
    other: int,
    score: float,
    threshold: float,
) -> DatabaseRelationship:
    if exact_content:
        return DatabaseRelationship.EXACT_CONTENT
    if exact_id_set:
        return DatabaseRelationship.EXACT_ID_SET
    if shared and shared == selected and selected < other:
        return DatabaseRelationship.SUBSET
    if shared and shared == other and other < selected:
        return DatabaseRelationship.SUPERSET
    if shared and score >= threshold:
        return DatabaseRelationship.NEAR_MATCH
    if shared:
        return DatabaseRelationship.PARTIAL_OVERLAP
    return DatabaseRelationship.NO_OVERLAP


def _sort_key(
    comparison: DatabaseComparison,
    threshold: float,
) -> tuple[int, float, float, int, int, float, int, str, int]:
    if comparison.exact_content:
        group = 0
    elif comparison.exact_id_set:
        group = 1
    elif comparison.shared_ids and comparison.containment >= threshold:
        group = 2
    elif comparison.shared_ids:
        group = 3
    elif comparison.shared_sequence_checksums:
        group = 4
    else:
        group = 5
    return (
        group,
        -comparison.containment,
        -comparison.id_jaccard,
        -comparison.shared_ids,
        -comparison.shared_sequence_checksums,
        -comparison.description_jaccard,
        -comparison.shared_descriptions,
        comparison.database.filename.casefold(),
        comparison.database.database_id,
    )


def _compare_staged(
    connection: RegistryConnection,
    *,
    excluded_database_id: int | None,
    kind: EntryKind,
    threshold: float,
) -> list[DatabaseComparison]:
    selected_ids = _selected_id_count(connection)
    selected_sequences = _selected_sequence_count(connection)
    selected_descriptions = _selected_description_count(connection)
    selected_content_fingerprint = _selected_content_fingerprint(connection)
    database_rows = _database_rows(connection, excluded_database_id, kind)
    metrics = pair_metric_counts(
        connection,
        selection=PairMetricSelection(
            ids_table=f"{connection.temp(_SELECTED_IDS_TABLE)}",
            sequences_table=f"{connection.temp(_SELECTED_SEQUENCES_TABLE)}",
            descriptions_table=f"{connection.temp(_SELECTED_DESCRIPTIONS_TABLE)}",
            pairs_table=f"{connection.temp(_SELECTED_PAIRS_TABLE)}",
            where="1",
            params=(),
            kind_filtered=True,
        ),
        kind=kind,
        excluded_database_id=excluded_database_id,
    )

    return _comparisons_from_counts(
        kind=kind,
        database_rows=database_rows,
        selected_ids=selected_ids,
        selected_sequences=selected_sequences,
        selected_descriptions=selected_descriptions,
        selected_content_fingerprint=selected_content_fingerprint,
        shared_ids=metrics.shared_ids,
        shared_sequences=metrics.shared_sequences,
        shared_descriptions=metrics.shared_descriptions,
        shared_pairs=metrics.shared_pairs,
        matching_ids=metrics.matching_ids,
        threshold=threshold,
    )


def _comparisons_from_counts(
    *,
    kind: EntryKind,
    database_rows: list[_DatabaseRow],
    selected_ids: int,
    selected_sequences: int,
    selected_descriptions: int,
    selected_content_fingerprint: str,
    shared_ids: dict[int, int],
    shared_sequences: dict[int, int],
    shared_descriptions: dict[int, int],
    shared_pairs: dict[int, int],
    matching_ids: dict[int, int],
    threshold: float,
) -> list[DatabaseComparison]:
    """Build and rank public comparison records from aggregate pair metrics."""

    comparisons: list[DatabaseComparison] = []
    for database_row in database_rows:
        reference = database_row.database
        database_id = reference.database_id
        other_ids = database_row.ids
        other_sequences = database_row.sequences
        other_descriptions = database_row.descriptions
        shared = shared_ids.get(database_id, 0)
        shared_sequence_count = shared_sequences.get(database_id, 0)
        shared_description_count = shared_descriptions.get(database_id, 0)
        score = _containment_from_counts(shared, selected_ids, other_ids)
        exact_id_set = selected_ids == other_ids == shared and bool(shared)
        exact_sequence_set = (
            selected_sequences == other_sequences == shared_sequence_count
            and bool(shared_sequence_count)
        )
        exact_description_set = (
            selected_descriptions == other_descriptions == shared_description_count
            and bool(shared_description_count)
        )
        exact_content = (
            exact_id_set and selected_content_fingerprint == database_row.content_fingerprint
        )
        exact_pairs = shared_pairs.get(database_id, 0)
        changed_ids = shared - matching_ids.get(database_id, 0)
        relationship = _relationship(
            exact_content=exact_content,
            exact_id_set=exact_id_set,
            shared=shared,
            selected=selected_ids,
            other=other_ids,
            score=score,
            threshold=threshold,
        )
        comparisons.append(
            DatabaseComparison(
                database=reference,
                kind=kind,
                selected_ids=selected_ids,
                other_ids=other_ids,
                other_entries=database_row.entries,
                shared_ids=shared,
                selected_coverage=_coverage(shared, selected_ids),
                other_coverage=_coverage(shared, other_ids),
                containment=score,
                id_jaccard=_jaccard(shared, selected_ids, other_ids),
                selected_sequences=selected_sequences,
                other_sequences=other_sequences,
                shared_sequence_checksums=shared_sequence_count,
                selected_sequence_coverage=_coverage(shared_sequence_count, selected_sequences),
                other_sequence_coverage=_coverage(shared_sequence_count, other_sequences),
                sequence_containment=_containment_from_counts(
                    shared_sequence_count,
                    selected_sequences,
                    other_sequences,
                ),
                sequence_jaccard=_jaccard(
                    shared_sequence_count,
                    selected_sequences,
                    other_sequences,
                ),
                selected_descriptions=selected_descriptions,
                other_descriptions=other_descriptions,
                shared_descriptions=shared_description_count,
                selected_description_coverage=_coverage(
                    shared_description_count, selected_descriptions
                ),
                other_description_coverage=_coverage(shared_description_count, other_descriptions),
                description_containment=_containment_from_counts(
                    shared_description_count,
                    selected_descriptions,
                    other_descriptions,
                ),
                description_jaccard=_jaccard(
                    shared_description_count,
                    selected_descriptions,
                    other_descriptions,
                ),
                shared_exact_pairs=exact_pairs,
                changed_shared_ids=changed_ids,
                selected_only_ids=selected_ids - shared,
                other_only_ids=other_ids - shared,
                exact_id_set=exact_id_set,
                exact_sequence_set=exact_sequence_set,
                exact_description_set=exact_description_set,
                exact_content=exact_content,
                relationship=relationship,
            )
        )

    comparisons.sort(key=lambda comparison: _sort_key(comparison, threshold))
    return comparisons


def compare_database(
    source: Path | RegistryConnection,
    database_id: int,
    threshold: float = 0.99,
    *,
    kind: EntryKind = EntryKind.TARGET,
) -> list[DatabaseComparison]:
    """Compare one registered database using kind-scoped materialized statistics."""
    _validate_comparison_kind(kind)
    with _registry_connection(source) as connection:
        selected = connection.execute(
            """
            SELECT databases.detail_level, stats.distinct_ids,
                   stats.distinct_sequences, stats.distinct_descriptions,
                   stats.content_fingerprint
            FROM databases
            LEFT JOIN database_kind_stats AS stats
              ON stats.database_id = databases.id
             AND stats.kind = ?
            WHERE databases.id = ?
            """,
            (kind.value, database_id),
        ).fetchone()
        if selected is None:
            raise ValueError(f"database id {database_id} is not registered")
        if DetailLevel(str(selected[0])) is DetailLevel.METADATA_ONLY:
            raise ValueError(
                f"database id {database_id} is metadata-only because it exceeds the configured detail limit; "
                "sequence comparison is unavailable"
            )
        if any(selected[index] is None for index in (1, 2, 3, 4)):
            raise RegistryIntegrityError(
                f"kind statistics are incomplete for database id {database_id}; run a full registry reindex"
            )
        database_rows = _database_rows(connection, database_id, kind)
        pair_rows = connection.execute(
            """
            SELECT CASE
                       WHEN database_id_low = ? THEN database_id_high
                   ELSE database_id_low
                   END AS other_database_id,
                   shared_ids,
                   shared_sequence_checksums,
                   shared_descriptions,
                   shared_exact_pairs,
                   matching_shared_ids
            FROM database_pair_stats
            JOIN databases AS low ON low.id = database_pair_stats.database_id_low
            JOIN databases AS high ON high.id = database_pair_stats.database_id_high
            WHERE low.detail_level = ?
              AND high.detail_level = ?
              AND kind = ?
              AND (database_id_low = ? OR database_id_high = ?)
            """,
            (
                database_id,
                DetailLevel.FULL.value,
                DetailLevel.FULL.value,
                kind.value,
                database_id,
                database_id,
            ),
        ).fetchall()
        if len(pair_rows) != len(database_rows):
            raise RegistryIntegrityError(
                f"pair-statistics cache is incomplete for database id {database_id}; run a full registry reindex"
            )
        shared_ids = {int(row[0]): int(row[1]) for row in pair_rows}
        shared_sequences = {int(row[0]): int(row[2]) for row in pair_rows}
        shared_descriptions = {int(row[0]): int(row[3]) for row in pair_rows}
        shared_pairs = {int(row[0]): int(row[4]) for row in pair_rows}
        matching_ids = {int(row[0]): int(row[5]) for row in pair_rows}
        return _comparisons_from_counts(
            kind=kind,
            database_rows=database_rows,
            selected_ids=int(selected[1]),
            selected_sequences=int(selected[2]),
            selected_descriptions=int(selected[3]),
            selected_content_fingerprint=str(selected[4]),
            shared_ids=shared_ids,
            shared_sequences=shared_sequences,
            shared_descriptions=shared_descriptions,
            shared_pairs=shared_pairs,
            matching_ids=matching_ids,
            threshold=threshold,
        )


def compare_candidate(
    connection: RegistryConnection,
    threshold: float = 0.99,
    *,
    kind: EntryKind = EntryKind.TARGET,
) -> list[DatabaseComparison]:
    """Compare one kind in the connection-local candidate with every database.

    The caller first populates :data:`protein_fasta.registry.indexing.CANDIDATE_TABLE` using
    the registry's streaming candidate parser.  Keeping the table temporary
    avoids persisting uploads before the user chooses to build them.
    """
    _validate_comparison_kind(kind)
    _stage_from_table(connection, f"{connection.temp(CANDIDATE_TABLE)}", "database_id = 0", kind)
    return _compare_staged(
        connection,
        excluded_database_id=None,
        kind=kind,
        threshold=threshold,
    )


def find_best_overlap(
    ids: set[str],
    source: Path | RegistryConnection,
    threshold: float,
) -> OverlapHit | None:
    """Return the best SQLite registry hit at or above ``threshold``.

    This ID-only fast path preserves the pre-SQLite duplicate-check API.  Rich
    sequence comparisons are available through :func:`compare_candidate` after
    populating the temporary candidate table.
    """
    if not ids:
        return None
    with _registry_connection(source) as connection:
        selected_count = len(ids)
        fingerprint = id_set_fingerprint(sorted(ids))
        exact = connection.execute(
            """
            SELECT databases.id, databases.filename, databases.annotation
            FROM databases
            JOIN database_kind_stats AS target_stats
              ON target_stats.database_id = databases.id
             AND target_stats.kind = 'target'
            WHERE target_stats.distinct_ids = ?
              AND target_stats.id_fingerprint = ?
              AND databases.detail_level = ?
            ORDER BY databases.filename COLLATE NOCASE, databases.id
            LIMIT 1
            """,
            (selected_count, fingerprint, DetailLevel.FULL.value),
        ).fetchone()
        if exact is not None:
            return OverlapHit(
                record=DatabaseReference(
                    database_id=int(exact[0]),
                    filename=str(exact[1]),
                    annotation=None if exact[2] is None else str(exact[2]),
                ),
                score=1.0,
                exact=True,
            )

        _stage_target_ids(connection, ids)
        rows = connection.execute(
            f"""
            SELECT databases.id,
                   databases.filename,
                   databases.annotation,
                   target_stats.distinct_ids,
                   COUNT(DISTINCT entries.sequence_id) AS shared_target_ids
            FROM entries
            JOIN {connection.temp(_SELECTED_IDS_TABLE)} AS selected
              ON selected.sequence_id = entries.sequence_id
            JOIN databases ON databases.id = entries.database_id
            JOIN database_kind_stats AS target_stats
              ON target_stats.database_id = databases.id
             AND target_stats.kind = 'target'
            WHERE entries.kind = 'target'
              AND databases.detail_level = 'full'
            GROUP BY databases.id, databases.filename, databases.annotation, target_stats.distinct_ids
            """
        ).fetchall()

    candidates: list[tuple[tuple[float, float, int, str, int], OverlapHit]] = []
    for row in rows:
        database_id = int(row[0])
        filename = str(row[1])
        other_count = int(row[3])
        shared = int(row[4])
        score = _containment_from_counts(shared, selected_count, other_count)
        if score < threshold:
            continue
        hit = OverlapHit(
            record=DatabaseReference(
                database_id=database_id,
                filename=filename,
                annotation=None if row[2] is None else str(row[2]),
            ),
            score=score,
            exact=selected_count == other_count == shared,
        )
        candidates.append(
            (
                (
                    -score,
                    -_jaccard(shared, selected_count, other_count),
                    -shared,
                    filename.casefold(),
                    database_id,
                ),
                hit,
            )
        )
    return min(candidates, key=lambda candidate: candidate[0])[1] if candidates else None


def check_fasta(
    path: Path,
    source: Path | RegistryConnection,
    settings: RegistrySettings,
) -> OverlapHit | None:
    """Find the best target-ID overlap for a FASTA against the SQLite registry."""
    ids = target_id_set(path, settings)
    return find_best_overlap(ids, source, settings.overlap_threshold)
