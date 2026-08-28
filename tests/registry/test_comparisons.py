# pyright: basic, reportAttributeAccessIssue=false, reportOptionalSubscript=false
from __future__ import annotations

import hashlib
import sqlite3
from collections.abc import Iterator
from pathlib import Path

import pytest
from tests.registry_helpers import (
    execute_script,
    in_memory_registry,
    stamp_schema_version,
    stamp_sqlite_schema_version,
)

from protein_fasta.analytics.hashing import (
    content_fingerprint,
    description_set_fingerprint,
    id_set_fingerprint,
    sequence_hash,
)
from protein_fasta.registry.backend.base import RegistryConnection, RegistryIntegrityError
from protein_fasta.registry.comparisons import (
    DatabaseRelationship,
    compare_candidate,
    compare_database,
    containment,
    find_best_overlap,
)
from protein_fasta.registry.indexing import (
    CANDIDATE_TABLE,
    SCHEMA_VERSION,
    _refresh_database_pair_stats,
)
from protein_fasta.registry.kinds import EntryKind


def _checksum(sequence: str) -> bytes:
    return sequence_hash(sequence)


def _description_hash(entry_id: str) -> bytes:
    return hashlib.blake2b(f"description {entry_id}".encode(), digest_size=16).digest()


def _content_fingerprint(
    entries: list[tuple[str, bytes, str]],
    kind: EntryKind,
) -> str:
    pairs = sorted(
        {
            (entry_id, checksum)
            for entry_id, checksum, entry_kind in entries
            if entry_kind == kind.value
        }
    )
    return content_fingerprint(pairs)


@pytest.fixture
def registry_connection() -> Iterator[RegistryConnection]:
    connection = in_memory_registry()
    execute_script(
        connection,
        """
        CREATE TABLE databases (
            id INTEGER PRIMARY KEY,
            filename TEXT NOT NULL,
            annotation TEXT,
            detail_level TEXT NOT NULL DEFAULT 'full'
        );
        CREATE TABLE entries (
            database_id INTEGER NOT NULL,
            ordinal INTEGER NOT NULL,
            sequence_id TEXT NOT NULL,
            kind TEXT NOT NULL,
            sequence_length INTEGER NOT NULL,
            sequence_hash BLOB NOT NULL,
            description_hash BLOB,
            PRIMARY KEY (database_id, ordinal)
        );
        CREATE INDEX entries_kind_db_id
            ON entries(kind, database_id, sequence_id);
        CREATE INDEX entries_kind_id_db
            ON entries(kind, sequence_id, database_id);
        CREATE INDEX entries_kind_db_sequence
            ON entries(kind, database_id, sequence_hash);
        CREATE INDEX entries_kind_sequence_db
            ON entries(kind, sequence_hash, database_id);
        CREATE INDEX entries_kind_description_db
            ON entries(kind, description_hash, database_id)
            WHERE description_hash IS NOT NULL;
        CREATE TABLE database_kind_stats (
            database_id INTEGER NOT NULL,
            kind TEXT NOT NULL,
            entry_count INTEGER NOT NULL,
            distinct_ids INTEGER NOT NULL,
            distinct_sequences INTEGER NOT NULL,
            distinct_descriptions INTEGER NOT NULL,
            distinct_pairs INTEGER NOT NULL,
            id_fingerprint TEXT NOT NULL,
            description_fingerprint TEXT NOT NULL,
            content_fingerprint TEXT NOT NULL,
            PRIMARY KEY (database_id, kind)
        );
        CREATE TABLE database_pair_stats (
            database_id_low INTEGER NOT NULL,
            database_id_high INTEGER NOT NULL,
            kind TEXT NOT NULL,
            shared_ids INTEGER NOT NULL,
            shared_sequence_checksums INTEGER NOT NULL,
            shared_descriptions INTEGER NOT NULL,
            shared_exact_pairs INTEGER NOT NULL,
            matching_shared_ids INTEGER NOT NULL,
            PRIMARY KEY (database_id_low, database_id_high, kind)
        );
        """,
    )
    stamp_schema_version(connection, SCHEMA_VERSION)
    yield connection
    connection.close()


def _insert_database(
    connection: RegistryConnection,
    filename: str,
    entries: list[tuple[str, bytes, str]],
) -> int:
    database_id = int(
        connection.execute("SELECT COALESCE(MAX(id), 0) + 1 FROM databases").fetchone()[0]
    )
    connection.execute(
        "INSERT INTO databases (id, filename, annotation) VALUES (?, ?, ?)",
        (
            database_id,
            filename,
            f"annotation for {filename}",
        ),
    )
    connection.executemany(
        """
        INSERT INTO entries (
            database_id, ordinal, sequence_id, kind,
            sequence_length, sequence_hash, description_hash
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        (
            (database_id, ordinal, entry_id, kind, 10, checksum, _description_hash(entry_id))
            for ordinal, (entry_id, checksum, kind) in enumerate(entries)
        ),
    )
    for kind in EntryKind:
        kind_entries = [
            (entry_id, checksum)
            for entry_id, checksum, entry_kind in entries
            if entry_kind == kind.value
        ]
        ids = {entry_id for entry_id, _ in kind_entries}
        sequences = {checksum for _, checksum in kind_entries}
        descriptions = {_description_hash(entry_id) for entry_id, _ in kind_entries}
        pairs = set(kind_entries)
        connection.execute(
            """
            INSERT INTO database_kind_stats (
                database_id, kind, entry_count, distinct_ids,
                distinct_sequences, distinct_descriptions, distinct_pairs,
                id_fingerprint, description_fingerprint, content_fingerprint
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                database_id,
                kind.value,
                len(kind_entries),
                len(ids),
                len(sequences),
                len(descriptions),
                len(pairs),
                id_set_fingerprint(sorted(ids)),
                description_set_fingerprint(sorted(descriptions)),
                _content_fingerprint(entries, kind),
            ),
        )
    _refresh_database_pair_stats(connection, database_id)
    return database_id


def _stage_registered_database_as_candidate(
    connection: RegistryConnection,
    database_id: int,
) -> None:
    connection.execute(
        f"""
        CREATE TEMP TABLE {CANDIDATE_TABLE} (
            database_id INTEGER NOT NULL,
            ordinal INTEGER NOT NULL,
            sequence_id TEXT NOT NULL,
            kind TEXT NOT NULL,
            sequence_length INTEGER NOT NULL,
            sequence_hash BLOB NOT NULL,
            description_hash BLOB,
            PRIMARY KEY (database_id, ordinal)
        )
        """
    )
    connection.execute(
        f"""
        INSERT INTO temp.{CANDIDATE_TABLE} (
            database_id, ordinal, sequence_id, kind,
            sequence_length, sequence_hash, description_hash
        )
        SELECT 0, ordinal, sequence_id, kind,
               sequence_length, sequence_hash, description_hash
        FROM entries
        WHERE database_id = ?
        """,
        (database_id,),
    )


@pytest.fixture
def comparison_registry(
    registry_connection: RegistryConnection,
) -> tuple[RegistryConnection, int]:
    connection = registry_connection
    checksum_a = _checksum("AAAA")
    checksum_b = _checksum("BBBB")
    checksum_c = _checksum("CCCC")
    checksum_x = _checksum("XXXX")
    contaminant_a = _checksum("CONTAMINANT-A")
    contaminant_b = _checksum("CONTAMINANT-B")
    contaminant_x = _checksum("CONTAMINANT-X")

    selected_id = _insert_database(
        connection,
        "00-selected.fasta",
        [
            ("a", checksum_a, "target"),
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("cont-a", contaminant_a, "contaminant"),
            ("cont-a", contaminant_a, "contaminant"),
            ("cont-b", contaminant_b, "contaminant"),
        ],
    )
    _insert_database(
        connection,
        "10-exact-content.fasta",
        [
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("ignored-decoy", checksum_x, "decoy"),
            ("cont-a", contaminant_a, "contaminant"),
            ("cont-b", contaminant_b, "contaminant"),
            ("ignored-sentinel", checksum_x, "sentinel"),
        ],
    )
    _insert_database(
        connection,
        "20-exact-ids-changed.fasta",
        [
            ("a", checksum_x, "target"),
            ("b", checksum_b, "target"),
            ("cont-a", contaminant_x, "contaminant"),
            ("cont-b", contaminant_b, "contaminant"),
        ],
    )
    _insert_database(
        connection,
        "30-selected-is-subset.fasta",
        [
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("c", checksum_c, "target"),
        ],
    )
    _insert_database(
        connection,
        "40-selected-is-superset.fasta",
        [("a", checksum_a, "target")],
    )
    _insert_database(
        connection,
        "50-partial.fasta",
        [("b", checksum_b, "target"), ("c", checksum_c, "target")],
    )
    _insert_database(
        connection,
        "60-renamed-sequences.fasta",
        [
            ("renamed-a", checksum_a, "target"),
            ("renamed-b", checksum_b, "target"),
            ("renamed-cont-a", contaminant_a, "contaminant"),
            ("renamed-cont-b", contaminant_b, "contaminant"),
        ],
    )
    _insert_database(
        connection,
        "70-no-overlap.fasta",
        [("x", checksum_x, "target")],
    )
    _insert_database(
        connection,
        "80-near.fasta",
        [
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("x", checksum_x, "target"),
        ],
    )
    connection.commit()
    return connection, selected_id


def test_containment_edges() -> None:
    assert containment({"a", "b", "c"}, {"a", "b", "c"}) == 1.0
    assert containment({"a", "b", "c"}, {"x", "y", "z"}) == 0.0
    assert containment(set(), {"a"}) == 0.0
    assert containment({"a", "b"}, {"a", "b", "c", "d"}) == 1.0


def test_registered_comparison_reports_id_and_content_metrics(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    comparisons = compare_database(connection, selected_id)
    by_filename = {comparison.database.filename: comparison for comparison in comparisons}

    exact = by_filename["10-exact-content.fasta"]
    assert exact.kind is EntryKind.TARGET
    assert exact.selected_ids == 2
    assert exact.other_ids == 2
    assert exact.shared_ids == 2
    assert exact.selected_coverage == 1.0
    assert exact.other_coverage == 1.0
    assert exact.containment == 1.0
    assert exact.id_jaccard == 1.0
    assert exact.selected_sequences == 2
    assert exact.other_sequences == 2
    assert exact.shared_sequence_checksums == 2
    assert exact.sequence_jaccard == 1.0
    assert exact.shared_exact_pairs == 2
    assert exact.changed_shared_ids == 0
    assert exact.selected_only_ids == 0
    assert exact.other_only_ids == 0
    assert exact.exact_id_set is True
    assert exact.exact_content is True
    assert exact.relationship is DatabaseRelationship.EXACT_CONTENT

    changed = by_filename["20-exact-ids-changed.fasta"]
    assert changed.exact_id_set is True
    assert changed.exact_content is False
    assert changed.shared_sequence_checksums == 1
    assert changed.sequence_jaccard == pytest.approx(1 / 3)
    assert changed.shared_exact_pairs == 1
    assert changed.changed_shared_ids == 1
    assert changed.relationship is DatabaseRelationship.EXACT_ID_SET


@pytest.mark.parametrize("kind", [EntryKind.TARGET, EntryKind.CONTAMINANT])
def test_registered_comparison_reads_materialized_metrics_without_entry_joins(
    comparison_registry: tuple[RegistryConnection, int],
    kind: EntryKind,
) -> None:
    connection, selected_id = comparison_registry
    statements: list[str] = []
    connection.raw.set_trace_callback(statements.append)
    try:
        compare_database(connection, selected_id, kind=kind)
    finally:
        connection.raw.set_trace_callback(None)

    assert any("database_pair_stats" in statement for statement in statements)
    assert not any("JOIN entries" in statement for statement in statements)


def test_registered_comparison_reports_incomplete_pair_cache(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    connection.execute(
        "DELETE FROM database_pair_stats WHERE kind = 'target' AND database_id_low = ?",
        (selected_id,),
    )

    with pytest.raises(
        RegistryIntegrityError,
        match=r"pair-statistics cache is incomplete.*full registry reindex",
    ):
        compare_database(connection, selected_id)


@pytest.mark.parametrize("kind", [EntryKind.TARGET, EntryKind.CONTAMINANT])
def test_live_candidate_and_materialized_registered_metrics_agree(
    comparison_registry: tuple[RegistryConnection, int],
    kind: EntryKind,
) -> None:
    connection, selected_id = comparison_registry
    materialized = {
        comparison.database.database_id: comparison
        for comparison in compare_database(connection, selected_id, kind=kind)
    }
    _stage_registered_database_as_candidate(connection, selected_id)

    live = {
        comparison.database.database_id: comparison
        for comparison in compare_candidate(connection, kind=kind)
        if comparison.database.database_id != selected_id
    }

    assert live == materialized


def test_pair_aggregation_pins_the_kind_specific_join_indexes(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    statements: list[str] = []
    connection.raw.set_trace_callback(statements.append)
    try:
        _refresh_database_pair_stats(connection, selected_id)
    finally:
        connection.raw.set_trace_callback(None)

    sql = "\n".join(statements)
    assert "INDEXED BY entries_kind_id_db" in sql
    assert "INDEXED BY entries_kind_sequence_db" in sql
    assert "INDEXED BY entries_kind_description_db" in sql


def test_contaminant_comparison_keeps_id_and_sequence_evidence_separate(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    by_filename = {
        comparison.database.filename: comparison
        for comparison in compare_database(
            connection,
            selected_id,
            kind=EntryKind.CONTAMINANT,
        )
    }

    exact = by_filename["10-exact-content.fasta"]
    assert exact.kind is EntryKind.CONTAMINANT
    assert exact.selected_ids == exact.other_ids == exact.shared_ids == 2
    assert exact.selected_sequences == exact.other_sequences == 2
    assert exact.shared_sequence_checksums == 2
    assert exact.id_jaccard == exact.sequence_jaccard == 1.0
    assert exact.exact_id_set is True
    assert exact.exact_content is True
    assert exact.relationship is DatabaseRelationship.EXACT_CONTENT

    changed = by_filename["20-exact-ids-changed.fasta"]
    assert changed.shared_ids == 2
    assert changed.shared_sequence_checksums == 1
    assert changed.shared_exact_pairs == 1
    assert changed.changed_shared_ids == 1
    assert changed.id_jaccard == 1.0
    assert changed.sequence_jaccard == pytest.approx(1 / 3)
    assert changed.exact_id_set is True
    assert changed.exact_content is False

    renamed = by_filename["60-renamed-sequences.fasta"]
    assert renamed.shared_ids == 0
    assert renamed.id_jaccard == 0.0
    assert renamed.shared_sequence_checksums == 2
    assert renamed.sequence_jaccard == 1.0
    assert renamed.shared_exact_pairs == 0

    empty = by_filename["70-no-overlap.fasta"]
    assert empty.other_ids == 0
    assert empty.other_sequences == 0
    assert empty.shared_ids == 0
    assert empty.shared_sequence_checksums == 0
    assert empty.other_coverage == 0.0
    assert empty.id_jaccard == empty.sequence_jaccard == 0.0


def test_decoy_packaging_does_not_change_target_comparison(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    exact = next(
        comparison
        for comparison in compare_database(connection, selected_id)
        if comparison.database.filename == "10-exact-content.fasta"
    )

    assert exact.kind is EntryKind.TARGET
    assert exact.exact_content is True
    assert exact.selected_ids == exact.other_ids == exact.shared_ids == 2
    assert exact.selected_sequences == exact.other_sequences == 2
    assert exact.shared_sequence_checksums == 2


@pytest.mark.parametrize("kind", [EntryKind.DECOY, EntryKind.SENTINEL])
def test_generated_and_metadata_kinds_are_not_comparable(
    comparison_registry: tuple[RegistryConnection, int],
    kind: EntryKind,
) -> None:
    connection, selected_id = comparison_registry

    with pytest.raises(ValueError, match="Only target and contaminant"):
        compare_database(connection, selected_id, kind=kind)
    with pytest.raises(ValueError, match="Only target and contaminant"):
        compare_candidate(connection, kind=kind)


def test_registered_comparison_distinguishes_direction_and_renamed_sequences(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    by_filename = {
        comparison.database.filename: comparison
        for comparison in compare_database(connection, selected_id)
    }

    subset = by_filename["30-selected-is-subset.fasta"]
    assert subset.selected_coverage == 1.0
    assert subset.other_coverage == pytest.approx(2 / 3)
    assert subset.containment == 1.0
    assert subset.id_jaccard == pytest.approx(2 / 3)
    assert subset.selected_only_ids == 0
    assert subset.other_only_ids == 1
    assert subset.relationship is DatabaseRelationship.SUBSET

    superset = by_filename["40-selected-is-superset.fasta"]
    assert superset.selected_coverage == 0.5
    assert superset.other_coverage == 1.0
    assert superset.relationship is DatabaseRelationship.SUPERSET

    renamed = by_filename["60-renamed-sequences.fasta"]
    assert renamed.shared_ids == 0
    assert renamed.shared_sequence_checksums == 2
    assert renamed.sequence_jaccard == 1.0
    assert renamed.shared_exact_pairs == 0
    assert renamed.relationship is DatabaseRelationship.NO_OVERLAP


def test_comparisons_are_ranked_exact_near_nonzero_then_no_overlap(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    comparisons = compare_database(connection, selected_id, threshold=0.6)

    relationships = [comparison.relationship for comparison in comparisons]
    assert relationships[:5] == [
        DatabaseRelationship.EXACT_CONTENT,
        DatabaseRelationship.EXACT_ID_SET,
        DatabaseRelationship.SUBSET,
        DatabaseRelationship.SUBSET,
        DatabaseRelationship.SUPERSET,
    ]
    assert relationships[5] is DatabaseRelationship.PARTIAL_OVERLAP
    assert relationships[6:] == [
        DatabaseRelationship.NO_OVERLAP,
        DatabaseRelationship.NO_OVERLAP,
    ]
    assert comparisons[6].database.filename == "60-renamed-sequences.fasta"
    assert comparisons[6].shared_sequence_checksums == 2
    assert comparisons[7].database.filename == "70-no-overlap.fasta"


def test_near_match_classification_is_not_limited_to_subsets(
    registry_connection: RegistryConnection,
) -> None:
    checksum_a = _checksum("AAAA")
    checksum_b = _checksum("BBBB")
    checksum_c = _checksum("CCCC")
    checksum_x = _checksum("XXXX")
    selected_id = _insert_database(
        registry_connection,
        "selected.fasta",
        [
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("c", checksum_c, "target"),
        ],
    )
    _insert_database(
        registry_connection,
        "near.fasta",
        [
            ("a", checksum_a, "target"),
            ("b", checksum_b, "target"),
            ("x", checksum_x, "target"),
        ],
    )

    comparison = compare_database(registry_connection, selected_id, threshold=0.6)[0]
    assert comparison.containment == pytest.approx(2 / 3)
    assert comparison.relationship is DatabaseRelationship.NEAR_MATCH


def test_temporary_candidate_is_compared_without_registration(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, selected_id = comparison_registry
    _stage_registered_database_as_candidate(connection, selected_id)
    connection.execute(
        f"""
        INSERT INTO temp.{CANDIDATE_TABLE} (
            database_id, ordinal, sequence_id, kind,
            sequence_length, sequence_hash
        ) VALUES (1, 0, 'per-file-only', 'target', 4, ?)
        """,
        (_checksum("ONLY"),),
    )

    comparisons = compare_candidate(connection)
    assert all(comparison.selected_ids == 2 for comparison in comparisons)
    matches = {
        comparison.database.filename: comparison
        for comparison in comparisons
        if comparison.exact_content
    }
    assert "00-selected.fasta" in matches
    assert "10-exact-content.fasta" in matches


def test_find_best_overlap_uses_sqlite_registry(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, _ = comparison_registry
    hit = find_best_overlap({"a", "b"}, connection, threshold=0.99)

    assert hit is not None
    assert hit.exact is True
    assert hit.score == 1.0
    assert hit.record.filename == "00-selected.fasta"
    assert find_best_overlap(set(), connection, threshold=0.99) is None


def test_find_best_overlap_nonexact_path_runs_one_id_only_entry_join(
    comparison_registry: tuple[RegistryConnection, int],
) -> None:
    connection, _ = comparison_registry
    statements: list[str] = []
    connection.raw.set_trace_callback(statements.append)
    try:
        hit = find_best_overlap({"b", "novel"}, connection, threshold=0.5)
    finally:
        connection.raw.set_trace_callback(None)

    entry_selects = [
        statement
        for statement in statements
        if statement.lstrip().upper().startswith("SELECT") and "FROM entries" in statement
    ]
    assert hit is not None
    assert hit.exact is False
    assert hit.score == 0.5
    assert len(entry_selects) == 1
    assert "sequence_hash" not in entry_selects[0]


def test_find_best_overlap_accepts_registry_path(
    comparison_registry: tuple[RegistryConnection, int],
    tmp_path: Path,
) -> None:
    connection, _ = comparison_registry
    registry_path = tmp_path / "registry.sqlite3"
    target = sqlite3.connect(registry_path)
    try:
        connection.raw.backup(target)
    finally:
        target.close()

    hit = find_best_overlap({"a", "b"}, registry_path, threshold=0.99)
    assert hit is not None
    assert hit.record.filename == "00-selected.fasta"


def test_registry_path_rejects_an_outdated_schema_before_comparing(
    comparison_registry: tuple[RegistryConnection, int],
    tmp_path: Path,
) -> None:
    connection, _ = comparison_registry
    registry_path = tmp_path / "registry-v2.sqlite3"
    target = sqlite3.connect(registry_path)
    try:
        connection.raw.backup(target)
        stamp_sqlite_schema_version(target, SCHEMA_VERSION - 1)
    finally:
        target.close()

    # Still a ValueError, so callers that already handled an unreadable registry
    # keep handling this one -- but the message now names the full rebuild, which
    # is the only reindex that can migrate it.
    with pytest.raises(
        ValueError,
        match=rf"has schema version {SCHEMA_VERSION - 1}.*reindex --full",
    ):
        find_best_overlap({"a", "b"}, registry_path, threshold=0.99)


def test_unknown_registered_database_is_rejected(
    registry_connection: RegistryConnection,
) -> None:
    with pytest.raises(ValueError, match="database id 999 is not registered"):
        compare_database(registry_connection, 999)
