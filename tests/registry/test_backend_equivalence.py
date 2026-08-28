# pyright: basic
"""The two backends must produce the same registry, not merely a working one.

This is the file that makes the second implementation safe. Everything else in
the suite proves one engine behaves; these tests prove the engines agree, which
is a different claim and the one that matters: a deployment can hold registries
written by either, the GUI compares fingerprints across them, and
``registry_meta.sequence_hash_version`` promises they mean the same thing.

The fingerprints are the sharp end. ``_record_id_fingerprint`` and
``_record_content_fingerprint`` hash rows in the order the query returns them, so
the ``ORDER BY`` collation is part of the fingerprint. DuckDB refuses
``COLLATE BINARY`` outright and its default text ordering is byte-wise instead;
that is the right answer, and it is only the right answer because it is checked
here rather than assumed.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from tests.registry_support import BackendSettings, Settings

from protein_fasta.analytics.clustering import ClusteringMetric
from protein_fasta.registry.backend import duckdb as duckdb_backend
from protein_fasta.registry.backend import factory
from protein_fasta.registry.backend.schema import REGISTRY_TABLES
from protein_fasta.registry.clustering import target_counts
from protein_fasta.registry.export import query_similarity_data
from protein_fasta.registry.indexing import (
    SCHEMA_VERSION,
    RegistryRecord,
    connect_registry,
    export_registry_stats_only,
    initialize_registry,
    list_databases,
    rebuild_registry,
)
from protein_fasta.registry.kinds import EntryKind

BACKENDS = ("sqlite", "duckdb")

# Identifiers chosen to separate a byte-wise ordering from a locale-aware or
# case-folding one: capitals sort before lower case by byte, digits before both,
# and the underscore sits above the capitals.
_ORDERING_TRAPS = ("Z_last", "a_lower", "A_upper", "0_digit", "_leading", "zz", "ZZ", "Ab", "aB")


def _corpus(root: Path) -> Path:
    """Write FASTAs whose identifiers and residues exercise the awkward cases."""
    root.mkdir(parents=True, exist_ok=True)
    first = root / "p1_db1_ordering_20260101.fasta"
    first.write_text(
        "".join(
            f">sp|{name}|{name}_HUMAN description of {name}\nMKVLAA\n" for name in _ORDERING_TRAPS
        )
    )
    second = root / "p2_db1_shared_20260101.fasta"
    second.write_text(
        # Overlaps the first on some ids so pair statistics are non-trivial, and
        # carries the normalizations: lower case, a trailing stop, tolerated
        # residues, and a duplicate id whose sequence matches after stripping.
        ">sp|a_lower|a_lower_HUMAN description of a_lower\nmkvlaa*\n"
        ">sp|A_upper|A_upper_HUMAN description of A_upper\nMKVLAA\n"
        ">sp|only_here|only_here_HUMAN unique entry\nMKXUBZOJ\n"
        ">REV_sp|a_lower|a_lower_HUMAN decoy\nAALVKM\n"
        ">sp|CON_KRT1|KRT1_HUMAN Keratin\nMKVLAAG\n"
    )
    return root


def _build(
    backend: str, tmp_path: Path, corpus: Path
) -> tuple[list[RegistryRecord], list[tuple[str, ...]], int]:
    """Rebuild the corpus on one backend and return everything it determined."""
    settings = Settings(
        fasta_root=corpus,
        registry_dir=tmp_path / backend,
        registry=BackendSettings(backend=backend),
    )
    path = tmp_path / backend / f"registry{factory.suffix_for(backend)}"
    rebuild_registry(corpus, path, settings)
    with connect_registry(path) as connection:
        records = sorted(list_databases(connection), key=lambda record: record.relative_path)
        pairs = [
            tuple(str(value) for value in row)
            for row in connection.execute(
                "SELECT database_id_low, database_id_high, kind, shared_ids, shared_sequence_checksums, "
                "shared_descriptions, shared_exact_pairs, matching_shared_ids FROM database_pair_stats "
                "ORDER BY database_id_low, database_id_high, kind"
            )
        ]
        entries = connection.scalar("SELECT COUNT(*) FROM entries")
    return records, pairs, entries


@pytest.fixture(scope="module")
def built(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[str, tuple[list[RegistryRecord], list, int]]:
    """Build the same corpus once per backend, because a rebuild is not cheap."""
    root = tmp_path_factory.mktemp("equivalence")
    corpus = _corpus(root / "databases")
    return {backend: _build(backend, root / "registries", corpus) for backend in BACKENDS}


def test_every_fingerprint_is_identical_across_backends(built: dict) -> None:
    """The claim the whole design rests on.

    A difference here would mean the same FASTA has two identities depending on
    which engine indexed it, so overlap detection, duplicate warnings, and the
    "exact content" relationship would all disagree across a mixed deployment.
    """
    sqlite_records, duckdb_records = built["sqlite"][0], built["duckdb"][0]
    assert [record.relative_path for record in sqlite_records] == [
        record.relative_path for record in duckdb_records
    ]

    for expected, actual in zip(sqlite_records, duckdb_records, strict=True):
        assert actual.target_id_fingerprint == expected.target_id_fingerprint, (
            expected.relative_path
        )
        assert actual.target_content_fingerprint == expected.target_content_fingerprint, (
            expected.relative_path
        )
        assert actual.target_description_fingerprint == expected.target_description_fingerprint
        for kind in EntryKind:
            assert (
                actual.kind_stats[kind].id_fingerprint == expected.kind_stats[kind].id_fingerprint
            )
            assert (
                actual.kind_stats[kind].content_fingerprint
                == expected.kind_stats[kind].content_fingerprint
            )
            assert (
                actual.kind_stats[kind].description_fingerprint
                == expected.kind_stats[kind].description_fingerprint
            )


def test_the_whole_record_is_identical_across_backends(built: dict) -> None:
    """Counts, statistics, normalization tallies, and namespaces all agree.

    Compared field by field rather than fingerprint by fingerprint, so a column
    added later that one backend fills differently fails here without anyone
    remembering to extend this test. ``indexed_at`` is excluded because it
    records when the build ran, not what it found.
    """
    for expected, actual in zip(built["sqlite"][0], built["duckdb"][0], strict=True):
        assert replace(actual, indexed_at="") == replace(expected, indexed_at="")


def test_materialized_pair_statistics_agree_row_for_row(built: dict) -> None:
    """The comparison views read only this table, so it is what users see."""
    assert built["duckdb"][1] == built["sqlite"][1]
    assert built["duckdb"][1], "the corpus must produce at least one pair row"


def test_both_backends_store_the_same_number_of_entries(built: dict) -> None:
    assert built["duckdb"][2] == built["sqlite"][2]


@pytest.mark.parametrize("backend", BACKENDS)
def test_empty_database_selection_means_no_rows(backend: str, tmp_path: Path) -> None:
    """An empty neighbourhood is a valid set, never the invalid SQL ``IN ()``."""
    settings = Settings(
        fasta_root=tmp_path / "databases",
        registry_dir=tmp_path / backend,
        registry=BackendSettings(backend=backend),
    )
    path = tmp_path / backend / f"registry{factory.suffix_for(backend)}"
    with connect_registry(path, backend=backend) as connection:
        initialize_registry(connection, settings)

        data = query_similarity_data(connection, database_ids=())
        counts = target_counts(connection, ClusteringMetric.TARGET_IDS, ())

    assert data.relative_paths == ()
    assert data.pairs == ()
    assert counts == {}


@pytest.mark.parametrize("table", sorted(REGISTRY_TABLES))
def test_the_duckdb_rewrite_keeps_every_column_in_order(table: str) -> None:
    """Guards the cost of writing DDL per backend rather than generating it.

    ``list_databases`` does ``SELECT *`` and several queries index rows
    positionally, so a column the rewrite dropped or reordered would corrupt
    records rather than raise. The DuckDB DDL is a regex rewrite of the SQLite
    text, which stays cheap and readable precisely because this holds it to
    account.
    """
    original = REGISTRY_TABLES[table]
    rewritten = duckdb_backend._portable(original, table=table)

    assert _declared_columns(rewritten) == _declared_columns(original)


def _declared_columns(ddl: str) -> list[str]:
    """Return the column names one CREATE TABLE declares, in order.

    Splits on top-level commas rather than on lines, because a declaration can
    wrap: ``database_id_low INTEGER NOT NULL`` and its ``REFERENCES`` clause sit
    on separate lines.
    """
    body = ddl[ddl.index("(") + 1 : ddl.rindex(")")]
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for character in body:
        if character == "," and depth == 0:
            parts.append("".join(current))
            current = []
            continue
        depth += (character == "(") - (character == ")")
        current.append(character)
    parts.append("".join(current))

    names = []
    for part in parts:
        collapsed = " ".join(part.split())
        if collapsed and not collapsed.upper().startswith(
            ("PRIMARY KEY", "CHECK", "FOREIGN KEY", "UNIQUE")
        ):
            names.append(collapsed.split()[0])
    return names


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_result_survives_running_another_statement(backend: str, tmp_path: Path) -> None:
    """Iterating one query while running another must not lose rows.

    DuckDB invalidates the previous result when the next statement runs on the
    same connection, and this package does read a cursor in a loop. Whether the
    backend materializes eagerly or hands out an independent cursor is its
    choice; that this pattern works is not.
    """
    path = tmp_path / f"registry{factory.suffix_for(backend)}"
    with connect_registry(path, backend=backend) as connection:
        connection.create_tables()
        connection.execute("INSERT INTO registry_meta (key, value) VALUES ('a', '1')")
        connection.execute("INSERT INTO registry_meta (key, value) VALUES ('b', '2')")
        connection.commit()

        seen = []
        for row in connection.execute("SELECT key FROM registry_meta ORDER BY key"):
            seen.append(str(row[0]))
            connection.execute("SELECT COUNT(*) FROM registry_meta").fetchone()
        assert seen == ["a", "b"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_failed_transaction_leaves_nothing_behind(backend: str, tmp_path: Path) -> None:
    """Atomicity, asserted portably rather than through a SQLite trigger.

    The existing rollback test installs a ``BEFORE INSERT`` trigger, which DuckDB
    cannot express. A CHECK constraint both engines declare fails the same way.
    """
    path = tmp_path / f"registry{factory.suffix_for(backend)}"
    with connect_registry(path, backend=backend) as connection:
        connection.create_tables()
        connection.execute("INSERT INTO registry_meta (key, value) VALUES ('keep', 'yes')")
        connection.commit()

        with (
            pytest.raises(Exception),  # noqa: B017 - either engine's constraint error
            connection.transaction(),
        ):
            connection.execute("INSERT INTO registry_meta (key, value) VALUES ('gone', 'no')")
            # database_id_low < database_id_high is declared by both.
            connection.execute(
                "INSERT INTO database_pair_stats (database_id_low, database_id_high, kind, shared_ids, "
                "shared_sequence_checksums, shared_descriptions, shared_exact_pairs, matching_shared_ids) "
                "VALUES (9, 1, 'target', 0, 0, 0, 0, 0)"
            )

        remaining = {str(row[0]) for row in connection.execute("SELECT key FROM registry_meta")}
        assert "gone" not in remaining
        assert "keep" in remaining


@pytest.mark.parametrize("backend", BACKENDS)
def test_deleting_a_database_removes_every_row_that_belongs_to_it(
    backend: str, tmp_path: Path
) -> None:
    """One backend declares ON DELETE CASCADE and one cannot; both must clear the children."""
    corpus = _corpus(tmp_path / "databases")
    settings = Settings(
        fasta_root=corpus,
        registry_dir=tmp_path / "registry",
        registry=BackendSettings(backend=backend),
    )
    path = tmp_path / "registry" / f"registry{factory.suffix_for(backend)}"
    rebuild_registry(corpus, path, settings)

    with connect_registry(path) as connection:
        records = list_databases(connection)
        assert records
        doomed = records[0].relative_path
        with connection.transaction():
            connection.delete_database(doomed)

        assert (
            connection.scalar("SELECT COUNT(*) FROM databases WHERE relative_path = ?", (doomed,))
            == 0
        )
        for table, column in (
            ("entries", "database_id"),
            ("database_kind_stats", "database_id"),
        ):
            orphans = connection.scalar(
                f"SELECT COUNT(*) FROM {table} LEFT JOIN databases ON databases.id = {table}.{column} "
                "WHERE databases.id IS NULL"
            )
            assert orphans == 0, f"{backend}/{table}"
        connection.check_physical_integrity()


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_stats_only_copy_carries_the_read_paths_and_drops_the_entries(
    backend: str, tmp_path: Path
) -> None:
    """The copy the application serves from, on either engine.

    ``export_registry_stats_only`` is the most engine-specific function in the
    package -- it attaches the source and writes the four read-path tables -- and
    the SQLite test for it spells ``sqlite_master`` and ``PRAGMA user_version``,
    so it cannot cover the other implementation.
    """
    corpus = _corpus(tmp_path / "databases")
    settings = Settings(
        fasta_root=corpus,
        registry_dir=tmp_path / "registry",
        registry=BackendSettings(backend=backend),
    )
    suffix = factory.suffix_for(backend)
    source = tmp_path / "registry" / f"registry{suffix}"
    rebuild_registry(corpus, source, settings)
    with connect_registry(source) as connection:
        expected = [record.relative_path for record in list_databases(connection)]
        expected_pairs = connection.scalar("SELECT COUNT(*) FROM database_pair_stats")
        full_entries = connection.scalar("SELECT COUNT(*) FROM entries")
    assert full_entries, "the corpus must store entry details for the copy to be smaller"

    destination = tmp_path / "export" / f"stats{suffix}"
    size_bytes, pair_rows = export_registry_stats_only(source, destination)

    assert size_bytes == destination.stat().st_size
    assert pair_rows == expected_pairs
    with connect_registry(destination) as exported:
        assert exported.schema_version() == SCHEMA_VERSION
        assert [record.relative_path for record in list_databases(exported)] == expected
        # The point of the copy: the read paths are complete and entries is not
        # merely empty but absent, which is where the size saving comes from.
        assert not exported.has_table("entries")


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_stats_only_copy_refuses_to_change_engine(backend: str, tmp_path: Path) -> None:
    """A physical copy cannot migrate between engines, and says so."""
    other = "duckdb" if backend == "sqlite" else "sqlite"
    source = tmp_path / f"registry{factory.suffix_for(backend)}"
    with connect_registry(source, backend=backend) as connection:
        connection.create_tables()
        connection.set_schema_version(SCHEMA_VERSION)
        connection.commit()

    with pytest.raises(ValueError, match="cannot change engine"):
        export_registry_stats_only(source, tmp_path / f"stats{factory.suffix_for(other)}")


def test_a_second_writer_is_told_which_command_to_stop() -> None:
    """DuckDB allows one writer, so the refusal has to name the way out.

    Not a lock test -- acquiring one from a second process here would be a slow
    and flaky way to check a string. This asserts the message an operator sees,
    which is the part that was easy to get wrong.
    """
    message = duckdb_backend._open_failure(
        Path("/registries/fasta_registry-20260101T000000Z.duckdb"),
        False,
        duckdb_backend.duckdb.Error("Could not set lock on file: Conflicting lock is held"),
    )

    assert "another process still held it" in message
    assert "one writer per file" in message
    assert "reindex --full" in message
    # The wait is bounded and the message says how long it waited, so an operator
    # can tell a busy registry from a stuck one.
    assert "30s" in message
