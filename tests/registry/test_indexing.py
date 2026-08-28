# pyright: basic, reportArgumentType=false, reportCallIssue=false, reportOptionalMemberAccess=false, reportOptionalSubscript=false
from __future__ import annotations

import datetime
import gzip
import json
import os
import sqlite3
from collections import Counter
from collections.abc import Iterator
from io import StringIO
from pathlib import Path
from typing import Never

import pytest
from loguru import logger
from tests.registry_helpers import (
    open_test_registry,
    stamp_schema_version,
    stamp_sqlite_schema_version,
)
from tests.registry_support import Settings

import protein_fasta.registry.indexing as registry
from protein_fasta.analytics.hashing import sequence_hash
from protein_fasta.build.metadata import build_sentinel_header
from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.record import ProteinDiagnostics
from protein_fasta.registry.backend.base import (
    RegistryBackendError,
    RegistryConnection,
    RegistryIntegrityError,
)
from protein_fasta.registry.backend.schema import STATS_ONLY_TABLES
from protein_fasta.registry.comparisons import compare_candidate, compare_database
from protein_fasta.registry.indexing import (
    CANDIDATE_TABLE,
    SCHEMA_VERSION,
    RegistryRecord,
    RegistrySettings,
    connect_registry,
    get_database,
    index_fasta,
    initialize_registry,
    iter_fasta_files,
    list_databases,
    open_registry,
    populate_candidate,
    populate_candidate_files,
    read_registry,
    rebuild_registry,
    target_id_set,
    update_registry,
)
from protein_fasta.registry.kinds import DetailLevel, EntryKind
from protein_fasta.registry.rules import load_registry_diagnostic_document


def _write_database(path: Path, settings: Settings, *, changed: bool = False) -> None:
    sentinel = build_sentinel_header(
        "p1_db1_test",
        "SQLite registry test",
        datetime.date(2026, 7, 21),
        settings.sentinel,
    )
    sequence_a = "DDDD" if changed else "CCCC"
    path.write_text(
        f">{sentinel}\nCRAP\n"
        ">sp|A|A_TEST first target\nAAAA\n"
        f">sp|A|A_TEST conflicting duplicate\n{sequence_a}\n"
        ">sp|B|B_TEST repeated sequence\naaaa\n"
        ">sp|Cont_X|CONT_TEST contaminant\nAC\n"
        ">REV_sp|A|A_TEST decoy\nTTTT\n"
    )


def test_sequence_hash_uses_the_already_normalized_sequence_exactly() -> None:
    assert sequence_hash("AAAA") != sequence_hash("aaaa")
    assert sequence_hash("AAAA") != sequence_hash("AA AA")
    assert sequence_hash("AAAA") != sequence_hash("AAA-")


def test_open_registry_initializes_an_explicit_path(settings: Settings, tmp_path: Path) -> None:
    registry_path = tmp_path / "alternate" / "registry.sqlite3"

    with open_registry(settings, registry_path) as connection:
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        tables = {
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            ).fetchall()
        }

    assert version == SCHEMA_VERSION
    assert {
        "registry_meta",
        "databases",
        "database_kind_stats",
        "entries",
        "database_pair_stats",
    } <= tables


def test_index_fasta_rejects_database_column_value_contract_drift(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    monkeypatch.setattr(
        registry,
        "_DATABASE_COLUMN_NAMES",
        (*registry._DATABASE_COLUMN_NAMES, "unexpected_column"),
    )

    with (
        open_registry(settings, tmp_path / "registry.sqlite3") as connection,
        pytest.raises(RuntimeError, match="44 columns but 43 values"),
    ):
        index_fasta(connection, fasta, settings, root=tmp_path)


def test_indexing_normalizes_a_terminal_stop_and_reports_other_symbols(
    settings: Settings,
    tmp_path: Path,
) -> None:
    """A trailing stop is a translated terminator; anything else is kept and reported.

    A sweep over the installed collection never refuses a file, so the dash is
    indexed as it stands and counted as an illegal residue instead.
    """
    fasta = tmp_path / "symbols.fasta.gz"
    registry_path = tmp_path / "registry.sqlite3"
    with gzip.open(fasta, "wt") as handle:
        handle.write(
            ">sp|STAR|STAR_TEST\nAAA*\n>sp|DASH|DASH_TEST\nAAA-\n>sp|LOWER|LOWER_TEST\naaa\n"
        )

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        indexed = index_fasta(connection, fasta, settings, root=tmp_path)
        rows = connection.execute(
            "SELECT sequence_id, sequence_length, sequence_hash FROM entries ORDER BY ordinal"
        ).fetchall()

    assert indexed.entry_count == 3
    assert indexed.total_residues == 10
    assert "*" not in indexed.aa_counts
    assert indexed.aa_counts["-"] == 1
    assert indexed.aa_counts["A"] == 9
    assert indexed.stop_stripped_entries == 1
    assert indexed.upper_cased_entries == 1
    assert indexed.illegal_residue_entries == 1
    assert indexed.illegal_residues == {"-": 1}
    assert [(row["sequence_id"], row["sequence_length"]) for row in rows] == [
        ("sp|STAR|STAR_TEST", 3),
        ("sp|DASH|DASH_TEST", 4),
        ("sp|LOWER|LOWER_TEST", 3),
    ]
    assert bytes(rows[0]["sequence_hash"]) == sequence_hash("AAA")
    assert bytes(rows[1]["sequence_hash"]) == sequence_hash("AAA-")
    assert bytes(rows[2]["sequence_hash"]) == sequence_hash("AAA")


def test_index_fasta_stores_entries_and_within_database_stats(
    settings: Settings, tmp_path: Path
) -> None:
    fasta = tmp_path / "p1_db1_test_d_20260721.fasta"
    registry_path = tmp_path / "registry.sqlite3"
    _write_database(fasta, settings)

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        assert connection.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION == 11
        assert (
            connection.execute(
                "SELECT value FROM registry_meta WHERE key = 'schema_version'"
            ).fetchone()[0]
            == "11"
        )
        assert (
            connection.execute(
                "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = 'database_pair_stats'"
            ).fetchone()[0]
            == 1
        )
        indexed = index_fasta(connection, fasta, settings, root=tmp_path)
        listed = list_databases(connection)

        assert len(listed) == 1
        assert listed[0] == indexed
        assert indexed.id is not None
        assert get_database(connection, indexed.id) == indexed
        assert indexed.relative_path == fasta.name
        assert indexed.dbname == "p1_db1_test"
        assert indexed.filename_is_decoy is True
        assert indexed.is_decoy is True
        assert indexed.entry_count == 6
        assert indexed.target_count == 3
        assert indexed.decoy_count == 1
        assert indexed.contaminant_count == 1
        assert indexed.sentinel_count == 1
        assert indexed.distinct_target_ids == 2
        assert indexed.distinct_target_sequences == 2
        assert indexed.distinct_target_descriptions == 3
        assert indexed.duplicate_target_id_occurrences == 1
        assert indexed.conflicting_target_ids == 1
        assert indexed.repeated_target_sequences == 1
        assert indexed.total_residues == 22
        assert indexed.aa_sample_size == indexed.entry_count
        assert indexed.length_min == 2
        assert indexed.length_max == 4
        assert "SQLite registry test" in (indexed.annotation or "")
        assert indexed.sentinel_header and indexed.sentinel_header.startswith("aa|p1_db1_test|")
        assert indexed.contaminant_markers == []
        assert indexed.target_id_fingerprint.startswith("blake2b-128:")
        assert indexed.target_description_fingerprint.startswith("blake2b-128:")
        assert indexed.target_content_fingerprint.startswith("blake2b-128:")
        assert set(indexed.kind_stats) == set(EntryKind)
        target_stats = indexed.kind_stats[EntryKind.TARGET]
        assert target_stats.entry_count == 3
        assert target_stats.distinct_ids == 2
        assert target_stats.distinct_sequences == 2
        assert target_stats.distinct_descriptions == 3
        assert target_stats.distinct_pairs == 3
        assert target_stats.duplicate_id_occurrences == 1
        assert target_stats.conflicting_ids == 1
        assert target_stats.repeated_sequences == 1
        assert target_stats.length_min == target_stats.length_max == 4
        assert target_stats.total_residues == 12
        assert target_stats.aa_sample_size == target_stats.entry_count
        assert target_stats.id_fingerprint == indexed.target_id_fingerprint
        assert target_stats.description_fingerprint == indexed.target_description_fingerprint
        assert target_stats.content_fingerprint == indexed.target_content_fingerprint
        assert indexed.kind_stats[EntryKind.CONTAMINANT].entry_count == 1
        assert indexed.kind_stats[EntryKind.DECOY].entry_count == 1
        assert indexed.kind_stats[EntryKind.SENTINEL].entry_count == 1
        assert connection.execute(
            "SELECT COUNT(*) FROM database_kind_stats WHERE database_id = ?",
            (indexed.id,),
        ).fetchone()[0] == len(EntryKind)

        rows = connection.execute(
            "SELECT ordinal, sequence_id, kind, contaminant_group, sequence_length, "
            "length(sequence_hash) AS checksum_bytes, length(description_hash) AS description_hash_bytes "
            "FROM entries ORDER BY ordinal"
        ).fetchall()
        assert [row["ordinal"] for row in rows] == list(range(6))
        assert rows[1]["sequence_id"] == "sp|A|A_TEST"
        assert rows[1]["kind"] == "target"
        assert rows[1]["contaminant_group"] is None
        assert rows[1]["sequence_length"] == 4
        assert rows[1]["checksum_bytes"] == 16
        assert rows[1]["description_hash_bytes"] == 16
        assert rows[4]["contaminant_group"] == "unlabelled"

    assert target_id_set(fasta, settings) == {"sp|A|A_TEST", "sp|B|B_TEST"}
    assert read_registry(registry_path)[0].filename == fasta.name


def test_same_ids_and_changed_sequence_have_different_content_fingerprint(
    settings: Settings,
    tmp_path: Path,
) -> None:
    fasta = tmp_path / "test.fasta"
    registry_path = tmp_path / "registry.sqlite3"
    _write_database(fasta, settings)
    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        first = index_fasta(connection, fasta, settings, root=tmp_path)
        _write_database(fasta, settings, changed=True)
        second = index_fasta(connection, fasta, settings, root=tmp_path)

    assert first.id == second.id
    assert first.target_id_fingerprint == second.target_id_fingerprint
    assert first.target_content_fingerprint != second.target_content_fingerprint


def test_description_sets_normalize_whitespace_exclude_empty_and_detect_case_changes(
    settings: Settings,
    tmp_path: Path,
) -> None:
    first_path = tmp_path / "first.fasta"
    second_path = tmp_path / "second.fasta"
    first_path.write_text(">sp|A|A_TEST Alpha   beta\nAAAA\n>sp|B|B_TEST\nBBBB\n")
    second_path.write_text(">sp|A|A_TEST\tAlpha beta\nAAAA\n>sp|B|B_TEST   \nBBBB\n")

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        first = index_fasta(connection, first_path, settings, root=tmp_path)
        second = index_fasta(connection, second_path, settings, root=tmp_path)
        assert first.distinct_target_descriptions == 1
        assert second.distinct_target_descriptions == 1
        comparison = compare_database(connection, first.id or 0)[0]
        assert comparison.shared_descriptions == 1
        assert comparison.description_jaccard == 1.0
        assert comparison.exact_description_set is True

        second_path.write_text(">sp|A|A_TEST alpha beta\nAAAA\n>sp|B|B_TEST\nBBBB\n")
        index_fasta(connection, second_path, settings, root=tmp_path)
        changed = compare_database(connection, first.id or 0)[0]

    assert changed.shared_descriptions == 0
    assert changed.description_jaccard == 0.0
    assert changed.exact_description_set is False


def test_record_loading_rejects_incomplete_kind_stats(settings: Settings, tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        record = index_fasta(connection, fasta, settings)
        connection.execute(
            "DELETE FROM database_kind_stats WHERE database_id = ? AND kind = 'contaminant'",
            (record.id,),
        )
        with pytest.raises(
            RegistryIntegrityError,
            match=r"Kind statistics are incomplete.*full registry reindex",
        ):
            list_databases(connection)


def test_marker_blocks_store_contaminant_provenance_and_ambiguous_blocks(
    settings: Settings,
    tmp_path: Path,
) -> None:
    fasta = tmp_path / "marker-blocks.fasta"
    fasta.write_text(
        ">aa|Cont_alpha|marker\nM\n"
        ">zh|FGCZContaminants2023|target-looking contaminant\nAAAA\n"
        ">REV_sp|D|D_TEST decoy inside block\nDD\n"
        ">aa|Cont_beta|marker\nM\n"
        ">sp|B|B_TEST target-looking contaminant\nBBBB\n"
        ">aa|Cont_old_one|marker\nM\n"
        ">aa|Cont_old_two|marker\nM\n"
        ">sp|C|C_TEST ambiguous provenance\nCCCC\n"
        ">aa|main|boundary\nM\n"
        ">sp|Cont_OUT|OUT_TEST unmarked contaminant\nOO\n"
        ">sp|T|T_TEST target\nTTTT\n"
    )

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        indexed = index_fasta(connection, fasta, settings)
        rows = connection.execute(
            "SELECT sequence_id, kind, contaminant_group FROM entries ORDER BY ordinal"
        ).fetchall()

    classified = {
        str(row["sequence_id"]): (str(row["kind"]), row["contaminant_group"]) for row in rows
    }
    assert classified["zh|FGCZContaminants2023|target-looking"] == ("contaminant", "alpha")
    assert classified["sp|B|B_TEST"] == ("contaminant", "beta")
    assert classified["sp|C|C_TEST"] == ("contaminant", "unlabelled")
    assert classified["sp|Cont_OUT|OUT_TEST"] == ("contaminant", "unlabelled")
    assert classified["REV_sp|D|D_TEST"] == ("decoy", None)
    assert classified["sp|T|T_TEST"] == ("target", None)
    assert indexed.contaminant_markers == [
        "Cont_alpha",
        "Cont_beta",
        "Cont_old_one",
        "Cont_old_two",
    ]
    assert indexed.contaminant_count == 4
    assert target_id_set(fasta, settings) == {"sp|T|T_TEST"}


def test_filename_and_content_decoy_evidence_are_stored_separately(
    settings: Settings,
    tmp_path: Path,
) -> None:
    filename_decoy = tmp_path / "p1_db2_named_d_20260721.fasta"
    content_decoy = tmp_path / "plain.fasta"
    filename_decoy.write_text(">sp|A|A_TEST\nAAAA\n")
    content_decoy.write_text(">REV_sp|B|B_TEST\nBBBB\n")

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        named = index_fasta(connection, filename_decoy, settings)
        content = index_fasta(connection, content_decoy, settings)

    assert named.filename_is_decoy is True
    assert named.decoy_count == 0
    assert named.is_decoy is True
    assert content.filename_is_decoy is False
    assert content.decoy_count == 1
    assert content.is_decoy is True


def test_pair_stats_are_symmetric_and_refreshed_with_changed_database(
    settings: Settings,
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nBBBB\n>sp|Cont_C|C_TEST\nCCCC\n")
    second.write_text(
        ">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nXXXX\n>sp|C|C_TEST\nBBBB\n>sp|Cont_C|C_TEST\nCCCC\n"
    )
    registry_path = tmp_path / "registry.sqlite3"

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        first_record = index_fasta(connection, first, settings, root=tmp_path)
        second_record = index_fasta(connection, second, settings, root=tmp_path)
        assert first_record.id is not None
        assert second_record.id is not None

        pairs = connection.execute(
            """
            SELECT database_id_low, database_id_high, kind, shared_ids,
                   shared_sequence_checksums, shared_exact_pairs,
                   matching_shared_ids
            FROM database_pair_stats
            ORDER BY kind
            """
        ).fetchall()
        assert len(pairs) == 2
        assert tuple(pairs[1]) == (
            min(first_record.id, second_record.id),
            max(first_record.id, second_record.id),
            "target",
            2,
            2,
            1,
            1,
        )
        assert tuple(pairs[0]) == (
            min(first_record.id, second_record.id),
            max(first_record.id, second_record.id),
            "contaminant",
            1,
            1,
            1,
            1,
        )

        first.write_text(">sp|A|A_TEST\nZZZZ\n>sp|B|B_TEST\nBBBB\n>sp|Cont_C|C_TEST\nCCCC\n")
        replaced = index_fasta(connection, first, settings, root=tmp_path)
        refreshed = connection.execute(
            """
            SELECT shared_ids, shared_sequence_checksums,
                   shared_exact_pairs, matching_shared_ids
            FROM database_pair_stats
            WHERE kind = 'target'
            """
        ).fetchone()

    assert replaced.id == first_record.id
    assert tuple(refreshed) == (2, 1, 0, 0)


def test_pair_stats_refresh_rolls_back_with_database_replacement(
    settings: Settings,
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">sp|A|A_TEST\nAAAA\n")
    second.write_text(">sp|A|A_TEST\nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        original = index_fasta(connection, first, settings, root=tmp_path)
        index_fasta(connection, second, settings, root=tmp_path)
        original_pair = tuple(connection.execute("SELECT * FROM database_pair_stats").fetchone())
        connection.execute(
            """
            CREATE TRIGGER reject_pair_refresh
            BEFORE INSERT ON database_pair_stats
            BEGIN
                SELECT RAISE(ABORT, 'pair refresh rejected');
            END
            """
        )
        connection.commit()

        first.write_text(">sp|A|A_TEST\nXXXX\n")
        # The driver error arrives as the backend-neutral one, which is what every
        # caller outside this package now catches.
        with pytest.raises(RegistryBackendError, match="pair refresh rejected"):
            index_fasta(connection, first, settings, root=tmp_path)

        assert original.id is not None
        persisted = get_database(connection, original.id)
        assert persisted is not None
        assert persisted.target_content_fingerprint == original.target_content_fingerprint
        assert (
            tuple(connection.execute("SELECT * FROM database_pair_stats").fetchone())
            == original_pair
        )


def test_candidate_table_is_temporary_and_supports_combined_fastas(
    settings: Settings, tmp_path: Path
) -> None:
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">sp|A|A_TEST\nAAAA\n")
    second.write_text(">sp|B|B_TEST\nBBBB\n")
    registry_path = tmp_path / "registry.sqlite3"

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        candidate = populate_candidate(connection, [first, second], settings, label="combined")
        assert candidate.filename == "combined"
        assert candidate.entry_count == 2
        assert candidate.distinct_target_ids == 2
        assert list_databases(connection) == []
        assert connection.execute(f"SELECT COUNT(*) FROM temp.{CANDIDATE_TABLE}").fetchone()[0] == 2

    reopened = sqlite3.connect(registry_path)
    try:
        assert (
            reopened.execute(
                "SELECT COUNT(*) FROM sqlite_temp_master WHERE type = 'table' AND name = ?",
                (CANDIDATE_TABLE,),
            ).fetchone()[0]
            == 0
        )
    finally:
        reopened.close()


def test_candidate_files_are_read_once_and_keep_file_and_combined_stats(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    registered = tmp_path / "registered.fasta"
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    registered.write_text(">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nCC\n")
    first.write_text(">sp|A|A_TEST\nAAAA\n")
    second.write_text(">sp|B|B_TEST\nCC\n")
    registry_path = tmp_path / "registry.sqlite3"

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        index_fasta(connection, registered, settings, root=tmp_path)

        read_counts: Counter[Path] = Counter()
        original_iter_diagnostics = registry.iter_protein_diagnostics

        def counted_iter_diagnostics(
            path: Path,
            rules: DiagnosticRules,
        ) -> Iterator[ProteinDiagnostics]:
            read_counts[path] += 1
            yield from original_iter_diagnostics(path, rules)

        monkeypatch.setattr(registry, "iter_protein_diagnostics", counted_iter_diagnostics)
        per_file, combined = populate_candidate_files(
            connection,
            [first, second],
            settings,
            labels=["alpha.fasta", "beta.fasta"],
            combined_label="combined",
        )
        comparisons = compare_candidate(connection)
        row_counts = dict(
            connection.execute(
                f"SELECT database_id, COUNT(*) FROM temp.{CANDIDATE_TABLE} GROUP BY database_id"
            ).fetchall()
        )

    assert read_counts == Counter({first: 1, second: 1})
    assert [record.filename for record in per_file] == ["alpha.fasta", "beta.fasta"]
    assert [record.entry_count for record in per_file] == [1, 1]
    assert [record.total_residues for record in per_file] == [4, 2]
    assert combined.filename == "combined"
    assert combined.entry_count == 2
    assert combined.total_residues == 6
    assert combined.length_q1 == 2.5
    assert combined.length_median == 3.0
    assert combined.length_q3 == 3.5
    assert row_counts == {0: 2, 1: 1, 2: 1}
    assert comparisons[0].database.filename == registered.name
    assert comparisons[0].exact_content is True


def test_a_candidate_on_its_way_in_is_refused_but_an_installed_file_is_reported(
    settings: Settings,
    tmp_path: Path,
) -> None:
    """Ingest refuses; a sweep records the same finding and keeps going.

    An installed database has to stay indexable whatever it contains, or one bad
    file in the collection would make the registry unbuildable.
    """
    candidate = tmp_path / "orf.fasta"
    candidate.write_text(">sp|A|A_TEST orf\nMPEP*TIDEK\n")

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        with pytest.raises(registry.FastaValidationError, match="stop codon inside the sequence"):
            populate_candidate_files(connection, candidate, settings, labels=["orf.fasta"])
        per_file, _ = populate_candidate_files(
            connection,
            candidate,
            settings,
            labels=["orf.fasta"],
            strict=False,
        )

    assert per_file[0].illegal_residue_entries == 1
    assert per_file[0].illegal_residues == {"*": 1}


def test_indexing_reports_identifier_namespaces_and_bare_identifiers(
    settings: Settings,
    tmp_path: Path,
) -> None:
    fasta = tmp_path / "mixed_20260101.fasta"
    fasta.write_text(
        ">sp|P02769|ALBU_BOVIN serum albumin\nMPEPTIDEK\n"
        ">ENST00000028008.9\nMAAAAAAAK\n"
        ">pf|Pf2004_000005200|Pf2004_000005200 hypothetical\nMCCCCCCCK\n"
        ">Cluster-838.148120;orf1 assembly output\nMDDDDDDDK\n"
    )

    with connect_registry(tmp_path / "registry.sqlite3") as connection:
        initialize_registry(connection, settings)
        indexed = index_fasta(connection, fasta, settings, root=tmp_path)

    assert indexed.id_namespaces == {"ensembl": 1, "pf|": 1, "unmatched": 1, "uniprot": 1}
    assert indexed.dominant_id_namespace == "ensembl"
    assert indexed.unmatched_id_entries == 1
    # Only the bare Ensembl transcript id carries no description at all.
    assert indexed.bare_identifier_entries == 1


def test_candidate_progress_is_reported_against_labels_not_cache_paths(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A downloaded candidate's path is a cache token, so no line may carry it."""
    candidate = tmp_path / "0f9c4b2e1d.fasta"
    candidate.write_text(">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nCC\n>sp|C|C_TEST\nDD\n")
    second = tmp_path / "7b1af03c55.fasta"
    second.write_text(">sp|D|D_TEST\nEE\n")
    # Tick on every entry: the real cadence needs both a thousand entries and
    # seconds of work to produce an in-flight line.
    monkeypatch.setattr(registry, "_ENTRY_PROGRESS_CHECK_INTERVAL", 1)
    monkeypatch.setattr(registry, "_ENTRY_PROGRESS_REPORT_SECONDS", 0.0)
    reported: list[str] = []

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        populate_candidate_files(
            connection,
            [candidate, second],
            settings,
            labels=[
                "Escherichia coli · UP000000625 · Swiss-Prot (reviewed)",
                "Homo sapiens · UP000005640 · Swiss-Prot (reviewed)",
            ],
            on_progress=reported.append,
        )

    assert reported[0] == "Indexing Escherichia coli · UP000000625 · Swiss-Prot (reviewed) ..."
    assert (
        "Indexing Escherichia coli · UP000000625 · Swiss-Prot (reviewed): 1 entries read"
        in reported
    )
    assert "Indexed Escherichia coli · UP000000625 · Swiss-Prot (reviewed): 3 entries" in reported
    assert "Indexed Homo sapiens · UP000005640 · Swiss-Prot (reviewed): 1 entries" in reported
    assert reported[-1] == "Merging 2 sources into the combined selection ..."
    assert not any(candidate.name in line or second.name in line for line in reported)


def test_single_candidate_preserves_filename_decoy_evidence(
    settings: Settings,
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "upload_0_R_fgcz_9606_20260721.fasta"
    candidate.write_text(">sp|A|A_TEST\nAAAA\n>REV_sp|A|A_TEST\nAAAA\n")
    original_name = "R_fgcz_9606_20260721.fasta"

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        per_file, combined = populate_candidate_files(
            connection,
            candidate,
            settings,
            labels=[original_name],
        )

    assert per_file[0].filename_is_decoy is True
    assert per_file[0].decoy_count == 1
    assert combined.filename_is_decoy is True
    assert combined.decoy_count == 1


def test_candidate_contaminant_override_preserves_sentinels_and_decoys(
    settings: Settings,
    tmp_path: Path,
) -> None:
    first = tmp_path / "fgcz.fasta"
    second = tmp_path / "universal.fasta"
    first.write_text(
        ">aa|metadata|raw set\nM\n>zh|FGCZContaminants2023|ordinary entry\nAAAA\n>REV_sp|D|D_TEST decoy\nDDDD\n"
    )
    second.write_text(">sp|U|U_TEST ordinary entry\nUU\n")

    with open_registry(settings, tmp_path / "registry.sqlite3") as connection:
        per_file, combined = populate_candidate_files(
            connection,
            [first, second],
            settings,
            kind_override=EntryKind.CONTAMINANT,
            contaminant_groups=["fgcz2023", "universal"],
        )
        rows = connection.execute(
            f"SELECT database_id, sequence_id, kind, contaminant_group "
            f"FROM temp.{CANDIDATE_TABLE} WHERE database_id > 0 ORDER BY database_id, ordinal"
        ).fetchall()

        with pytest.raises(ValueError, match="groups must match"):
            populate_candidate_files(
                connection,
                [first, second],
                settings,
                kind_override=EntryKind.CONTAMINANT,
                contaminant_groups=["only-one"],
            )
        with pytest.raises(ValueError, match="require a contaminant kind override"):
            populate_candidate_files(
                connection,
                first,
                settings,
                contaminant_groups=["fgcz2023"],
            )

    assert [(row["kind"], row["contaminant_group"]) for row in rows] == [
        ("sentinel", None),
        ("contaminant", "fgcz2023"),
        ("decoy", None),
        ("contaminant", "universal"),
    ]
    assert set(per_file[0].kind_stats) == set(EntryKind)
    assert per_file[0].sentinel_count == 1
    assert per_file[0].decoy_count == 1
    assert per_file[0].contaminant_count == 1
    assert per_file[1].contaminant_count == 1
    assert combined.target_count == 0
    assert combined.contaminant_count == 2
    assert combined.sentinel_count == 1
    assert combined.decoy_count == 1


def test_candidate_rejects_empty_sequence_id(settings: Settings, tmp_path: Path) -> None:
    malformed = tmp_path / "malformed.fasta"
    malformed.write_text("> \nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"

    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        with pytest.raises(
            ValueError,
            match=r"malformed\.fasta: FASTA entry ordinal 0 has an empty sequence ID",
        ):
            populate_candidate(connection, malformed, settings)
        with pytest.raises(
            ValueError,
            match=r"malformed\.fasta: FASTA entry ordinal 0 has an empty sequence ID",
        ):
            target_id_set(malformed, settings)


def test_rebuild_and_incremental_update_with_explicit_prune(
    settings: Settings, tmp_path: Path
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    first = fasta_root / "first.fasta"
    second = fasta_root / "second.fasta"
    first.write_text(">sp|A|A_TEST\nAAAA\n")
    second.write_text(">sp|B|B_TEST\nBBBB\n")
    registry_path = tmp_path / "registry.sqlite3"

    rebuilt = rebuild_registry(fasta_root, registry_path, settings)
    assert {record.filename for record in rebuilt} == {"first.fasta", "second.fasta"}
    with connect_registry(registry_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM database_pair_stats").fetchone()[0] == 2

    second.unlink()
    without_prune = update_registry(fasta_root, registry_path, settings)
    assert {record.filename for record in without_prune} == {"first.fasta", "second.fasta"}

    pruned = update_registry(fasta_root, registry_path, settings, prune=True)
    assert [record.filename for record in pruned] == ["first.fasta"]
    with connect_registry(registry_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM database_pair_stats").fetchone()[0] == 0


def test_rebuild_and_update_skip_fastas_at_the_file_size_limit(
    settings: Settings,
    tmp_path: Path,
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    small = fasta_root / "small.fasta"
    growing = fasta_root / "growing.fasta"
    oversized = fasta_root / "oversized.fasta"
    small.write_text(">sp|A|A_TEST\nAAAA\n")
    growing.write_text(">sp|B|B_TEST\nBBBB\n")
    limit_bytes = int(settings.max_fasta_file_size_gib * 1024**3)
    with oversized.open("wb") as handle:
        handle.truncate(limit_bytes)
    registry_path = tmp_path / "registry.sqlite3"
    captured = StringIO()
    sink_id = logger.add(captured, level="INFO", format="{message}")

    try:
        rebuilt = rebuild_registry(fasta_root, registry_path, settings)
        with growing.open("r+b") as handle:
            handle.truncate(limit_bytes)
        updated = update_registry(fasta_root, registry_path, settings)
    finally:
        logger.remove(sink_id)

    assert {record.filename for record in rebuilt} == {"growing.fasta", "small.fasta"}
    assert [record.filename for record in updated] == ["small.fasta"]
    log = captured.getvalue()
    assert "3 FASTA files; 2 eligible below 5 GiB; 1 skipped" in log
    assert f"skipped FASTA {oversized}: 5.00 GiB" in log
    assert f"skipped FASTA {growing}: 5.00 GiB" in log


def test_full_rebuild_defers_secondary_indexes_and_pair_refresh(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    (fasta_root / "first.fasta").write_text(">sp|A|A_TEST\nAAAA\n")
    (fasta_root / "second.fasta").write_text(">sp|B|B_TEST\nBBBB\n")
    registry_path = tmp_path / "registry.sqlite3"
    original_index = registry._index_fasta
    observed_indexes: list[set[str]] = []
    observed_refresh_modes: list[bool] = []

    def _tracked_index(
        connection: RegistryConnection,
        path: Path,
        configured: RegistrySettings,
        *,
        root: Path | None = None,
        refresh_pair_stats: bool,
        log_record: bool,
        progress_label: str | None = None,
    ) -> RegistryRecord:
        observed_indexes.append(
            {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'index' AND name LIKE 'entries_%'"
                )
            }
        )
        observed_refresh_modes.append(refresh_pair_stats)
        return original_index(
            connection,
            path,
            configured,
            root=root,
            refresh_pair_stats=refresh_pair_stats,
            log_record=log_record,
            progress_label=progress_label,
        )

    monkeypatch.setattr(registry, "_index_fasta", _tracked_index)

    rebuild_registry(fasta_root, registry_path, settings)

    assert observed_indexes == [set(), set()]
    assert observed_refresh_modes == [False, False]
    with connect_registry(registry_path) as connection:
        final_indexes = {
            str(row[0])
            for row in connection.execute("SELECT name FROM sqlite_master WHERE type = 'index'")
        }
    assert {
        "entries_kind_db_id",
        "entries_kind_id_db",
        "entries_kind_db_sequence",
        "entries_kind_sequence_db",
        "database_pair_stats_high",
    } <= final_indexes


def test_full_rebuild_logs_each_database_and_long_scan_heartbeats(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    (fasta_root / "first.fasta").write_text(">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nBBBB\n")
    (fasta_root / "second.fasta").write_text(">sp|C|C_TEST\nCCCC\n")
    monkeypatch.setattr(registry, "_INSERT_BATCH_SIZE", 1)
    monkeypatch.setattr(registry, "_ENTRY_PROGRESS_LOG_SECONDS", 0.0)
    monkeypatch.setattr(registry, "_ENTRY_PROGRESS_CHECK_INTERVAL", 1)
    captured = StringIO()
    sink_id = logger.add(captured, level="INFO", format="{message}")

    try:
        rebuild_registry(fasta_root, tmp_path / "registry.sqlite3", settings)
    finally:
        logger.remove(sink_id)

    log = captured.getvalue()
    assert "file 1/2 started: first.fasta" in log
    assert "file 1/2 reading first.fasta: 1 entries" in log
    assert "file 1/2 scan complete: 2 entries" in log
    assert "file 1/2 complete: first.fasta (2 entries" in log
    assert "file 2/2 started: second.fasta" in log
    assert "file 2/2 complete: second.fasta (1 entries" in log


def test_bulk_and_incremental_pair_statistics_are_identical(
    settings: Settings, tmp_path: Path
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    duplicated_contents = (
        ">sp|A|A_TEST\nAAAA\n>sp|A|A_CONFLICT\nAAAT\n>sp|B|B_TEST\nAAAA\n>sp|Cont_X|CONT\nXX\n"
    )
    files = [
        ("alpha.fasta", duplicated_contents),
        ("alpha-copy.fasta", duplicated_contents),
        ("beta.fasta", ">sp|A|A_TEST\nAAAA\n>sp|C|C_TEST\nCCCC\n>sp|Cont_X|CONT\nXY\n"),
        ("gamma.fasta", ">sp|D|D_TEST\nBBBB\n>sp|Cont_Y|CONT\nXX\n"),
    ]
    for filename, contents in files:
        (fasta_root / filename).write_text(contents)
    bulk_path = tmp_path / "bulk.sqlite3"
    incremental_path = tmp_path / "incremental.sqlite3"

    rebuild_registry(fasta_root, bulk_path, settings)
    with open_registry(settings, incremental_path) as connection:
        for fasta_path in iter_fasta_files(fasta_root):
            index_fasta(connection, fasta_path, settings, root=fasta_root)

    query = (
        "SELECT database_id_low, database_id_high, kind, shared_ids, "
        "shared_sequence_checksums, shared_exact_pairs, matching_shared_ids "
        "FROM database_pair_stats ORDER BY database_id_low, database_id_high, kind"
    )
    with connect_registry(bulk_path) as bulk, connect_registry(incremental_path) as incremental:
        bulk_rows = [tuple(row) for row in bulk.execute(query)]
        incremental_rows = [tuple(row) for row in incremental.execute(query)]

    assert bulk_rows == incremental_rows
    assert len(bulk_rows) == 12


def test_interrupted_bulk_rebuild_preserves_registry_and_removes_temporary_file(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    (fasta_root / "database.fasta").write_text(">sp|A|A_TEST\nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"
    rebuild_registry(fasta_root, registry_path, settings)
    original_bytes = registry_path.read_bytes()

    def _interrupt(
        connection: RegistryConnection,
        path: Path,
        configured: RegistrySettings,
        *,
        root: Path | None = None,
        refresh_pair_stats: bool,
        log_record: bool,
        progress_label: str | None = None,
    ) -> Never:
        del connection, path, configured, root, refresh_pair_stats, log_record, progress_label
        raise KeyboardInterrupt

    monkeypatch.setattr(registry, "_index_fasta", _interrupt)

    with pytest.raises(KeyboardInterrupt):
        rebuild_registry(fasta_root, registry_path, settings)

    assert registry_path.read_bytes() == original_bytes
    assert list(tmp_path.glob(".registry.sqlite3.*.tmp")) == []


def test_recursive_discovery_and_prune_are_explicit(settings: Settings, tmp_path: Path) -> None:
    fasta_root = tmp_path / "fastas"
    nested_root = fasta_root / "nested"
    nested_root.mkdir(parents=True)
    top = fasta_root / "top.fasta"
    nested = nested_root / "nested.fasta"
    top.write_text(">sp|T|T_TEST\nTTTT\n")
    nested.write_text(">sp|N|N_TEST\nNNNN\n")
    registry_path = tmp_path / "registry.sqlite3"

    assert list(iter_fasta_files(fasta_root)) == [top]
    assert set(iter_fasta_files(fasta_root, recursive=True)) == {top, nested}
    rebuilt = rebuild_registry(fasta_root, registry_path, settings)
    assert [record.relative_path for record in rebuilt] == ["top.fasta"]

    recursive = update_registry(fasta_root, registry_path, settings, recursive=True)
    assert {record.relative_path for record in recursive} == {"top.fasta", "nested/nested.fasta"}
    with connect_registry(registry_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM database_pair_stats").fetchone()[0] == 2

    top_level_prune = update_registry(fasta_root, registry_path, settings, prune=True)
    assert {record.relative_path for record in top_level_prune} == {
        "top.fasta",
        "nested/nested.fasta",
    }

    nested.unlink()
    pruned = update_registry(fasta_root, registry_path, settings, recursive=True, prune=True)
    assert [record.relative_path for record in pruned] == ["top.fasta"]


def test_an_older_schema_is_refused_with_the_command_that_fixes_it(
    settings: Settings,
    tmp_path: Path,
) -> None:
    """ "Reindex required" misleads someone who is already running reindex.

    The incremental sweep opens the existing registry and hits this same error, so
    the message has to name the full rebuild specifically.
    """
    registry_path = tmp_path / "registry.sqlite3"
    with sqlite3.connect(registry_path) as connection:
        stamp_sqlite_schema_version(connection, 2)

    with (
        connect_registry(registry_path) as connection,
        pytest.raises(registry.RegistrySchemaError) as failure,
    ):
        initialize_registry(connection, settings)

    assert failure.value.found == 2
    assert failure.value.expected == 11
    assert "fasta-gen reindex --full" in str(failure.value)
    assert "make reindex-full" in str(failure.value)
    assert "incremental reindex cannot migrate it" in str(failure.value)
    # Still a ValueError, so every caller that already handled an unreadable
    # registry keeps handling this one.
    assert isinstance(failure.value, ValueError)

    with pytest.raises(registry.RegistrySchemaError, match=registry_path.name):
        read_registry(registry_path)


def test_rebuild_rejects_malformed_fasta_and_keeps_valid_files(
    settings: Settings, tmp_path: Path
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    good = fasta_root / "good.fasta"
    good.write_text(">sp|A|A_TEST\nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"
    malformed = fasta_root / "malformed.fasta"
    malformed.write_text("> \nBBBB\n")
    legacy = fasta_root / "legacy.fasta"
    legacy.write_bytes(b">sp|B|B_TEST 20\xb0 C\nBBBB\n")
    rejections = []

    rebuilt = rebuild_registry(
        fasta_root,
        registry_path,
        settings,
        rejections=rejections,
    )

    assert [record.filename for record in rebuilt] == ["good.fasta"]
    assert [record.filename for record in read_registry(registry_path)] == ["good.fasta"]
    assert [rejection.path for rejection in rejections] == [legacy, malformed]
    assert "not valid UTF-8" in rejections[0].reason
    assert "empty sequence ID" in rejections[1].reason


def test_incremental_update_removes_a_database_that_becomes_invalid(
    settings: Settings, tmp_path: Path
) -> None:
    fasta_root = tmp_path / "fastas"
    fasta_root.mkdir()
    fasta = fasta_root / "database.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"
    rebuild_registry(fasta_root, registry_path, settings)
    fasta.write_bytes(b">sp|A|A_TEST 20\xb0 C\nAAAA\n")
    rejections = []

    records = update_registry(
        fasta_root,
        registry_path,
        settings,
        rejections=rejections,
    )

    assert records == []
    assert [rejection.path for rejection in rejections] == [fasta]


def test_registry_rejects_changed_classifier_without_reindex(
    settings: Settings, tmp_path: Path
) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    registry_path = tmp_path / "registry.sqlite3"
    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        index_fasta(connection, fasta, settings)

    document = load_registry_diagnostic_document().model_dump(mode="json")
    document["decoy_prefix"] = "DECOY_"
    document["classifiers"][0]["removable_prefix_patterns"] = ["^DECOY_"]
    changed_document = tmp_path / "changed-diagnostics.json"
    changed_document.write_text(json.dumps(document), encoding="utf-8")
    changed = settings.model_copy(update={"registry_diagnostics_path": changed_document})
    with connect_registry(registry_path) as connection:
        try:
            initialize_registry(connection, changed)
        except ValueError as error:
            assert "fasta_diagnostics_fingerprint changed" in str(error)
        else:
            raise AssertionError("Changed classification rules should require reindexing")


def test_registered_detail_limit_keeps_exact_metadata_without_large_entry_rows(
    settings: Settings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    limited = settings.model_copy(update={"max_detailed_entries": 2})
    full = tmp_path / "full.fasta"
    full.write_text(">sp|A|A_TEST\nAA\n>sp|B|B_TEST\nCCCC\n")
    large = tmp_path / "large.fasta"
    large.write_text(">sp|C|C_TEST\nA\n>sp|D|D_TEST\nCCC\n>sp|E|E_TEST\nDD\n")
    checksum_calls = 0
    original_hash = registry.sequence_hash

    def counted_hash(sequence: str) -> bytes:
        nonlocal checksum_calls
        checksum_calls += 1
        return original_hash(sequence)

    monkeypatch.setattr(registry, "sequence_hash", counted_hash)
    with open_registry(limited, tmp_path / "registry.sqlite3") as connection:
        full_record = index_fasta(connection, full, limited, root=tmp_path)
        calls_after_full = checksum_calls
        large_record = index_fasta(connection, large, limited, root=tmp_path)
        stored_counts = {
            int(row[0]): int(row[1])
            for row in connection.execute(
                "SELECT database_id, COUNT(*) FROM entries GROUP BY database_id"
            ).fetchall()
        }
        pair_count = int(
            connection.execute("SELECT COUNT(*) FROM database_pair_stats").fetchone()[0]
        )

    assert full_record.detail_level is DetailLevel.FULL
    assert stored_counts[full_record.id] == 2
    assert large_record.detail_level is DetailLevel.METADATA_ONLY
    assert large_record.id not in stored_counts
    assert checksum_calls - calls_after_full == limited.max_detailed_entries
    assert large_record.entry_count == large_record.target_count == 3
    assert large_record.total_residues == 6
    assert large_record.length_min == 1
    assert large_record.length_median == 2
    assert large_record.length_max == 3
    assert large_record.aa_counts == {"A": 1, "C": 3, "D": 2}
    assert large_record.distinct_target_ids is None
    assert large_record.target_id_fingerprint is None
    assert all(stats.distinct_ids is None for stats in large_record.kind_stats.values())
    assert all(stats.id_fingerprint is None for stats in large_record.kind_stats.values())
    assert pair_count == 0


def test_metadata_only_amino_acids_use_a_uniform_bounded_sample(
    settings: Settings,
    tmp_path: Path,
) -> None:
    limited = settings.model_copy(update={"max_detailed_entries": 2, "metadata_aa_sample_size": 2})
    fasta = tmp_path / "large.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n>sp|B|B_TEST\nCC\n>sp|C|C_TEST\nDDD\n>sp|D|D_TEST\nE\n")

    with open_registry(limited, tmp_path / "registry.sqlite3") as connection:
        record = index_fasta(connection, fasta, limited, root=tmp_path)
        assert record.id is not None
        loaded = get_database(connection, record.id)

    assert record.detail_level is DetailLevel.METADATA_ONLY
    assert loaded == record
    assert record.entry_count == 4
    assert record.total_residues == 10
    assert record.length_mean == 2.5
    assert record.aa_sample_size == 2
    assert record.aa_counts == {"A": 4, "D": 3}
    assert record.kind_stats[EntryKind.TARGET].aa_sample_size == 2
    assert record.kind_stats[EntryKind.TARGET].aa_counts == record.aa_counts


def test_reindexing_can_cross_the_detail_limit_in_both_directions(
    settings: Settings,
    tmp_path: Path,
) -> None:
    limited = settings.model_copy(update={"max_detailed_entries": 2})
    fasta = tmp_path / "changing.fasta"
    fasta.write_text(">sp|A|A_TEST\nAA\n>sp|B|B_TEST\nBB\n")

    with open_registry(limited, tmp_path / "registry.sqlite3") as connection:
        first = index_fasta(connection, fasta, limited, root=tmp_path)
        fasta.write_text(">sp|A|A_TEST\nAA\n>sp|B|B_TEST\nBB\n>sp|C|C_TEST\nCC\n")
        large = index_fasta(connection, fasta, limited, root=tmp_path)
        large_entry_count = int(
            connection.execute(
                "SELECT COUNT(*) FROM entries WHERE database_id = ?", (large.id,)
            ).fetchone()[0]
        )
        fasta.write_text(">sp|A|A_TEST\nAA\n>sp|B|B_TEST\nBB\n")
        restored = index_fasta(connection, fasta, limited, root=tmp_path)
        restored_entry_count = int(
            connection.execute(
                "SELECT COUNT(*) FROM entries WHERE database_id = ?", (restored.id,)
            ).fetchone()[0]
        )

    assert first.detail_level is DetailLevel.FULL
    assert large.detail_level is DetailLevel.METADATA_ONLY
    assert large_entry_count == 0
    assert restored.detail_level is DetailLevel.FULL
    assert restored_entry_count == 2


def test_metadata_only_database_is_excluded_from_sequence_comparisons(
    settings: Settings,
    tmp_path: Path,
) -> None:
    limited = settings.model_copy(update={"max_detailed_entries": 1})
    full = tmp_path / "full.fasta"
    full.write_text(">sp|A|A_TEST\nAA\n")
    large = tmp_path / "large.fasta"
    large.write_text(">sp|A|A_TEST\nAA\n>sp|B|B_TEST\nBB\n")

    with open_registry(limited, tmp_path / "registry.sqlite3") as connection:
        full_record = index_fasta(connection, full, limited, root=tmp_path)
        large_record = index_fasta(connection, large, limited, root=tmp_path)
        populate_candidate(connection, full, limited)
        candidate_comparisons = compare_candidate(connection)
        assert large_record.id is not None
        with pytest.raises(ValueError, match=r"metadata-only.*comparison is unavailable"):
            compare_database(connection, large_record.id)

    assert [comparison.database.database_id for comparison in candidate_comparisons] == [
        full_record.id
    ]


def test_registry_rejects_changed_detail_limit_without_full_reindex(
    settings: Settings,
    tmp_path: Path,
) -> None:
    registry_path = tmp_path / "registry.sqlite3"
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    limited = settings.model_copy(update={"max_detailed_entries": 2})
    with open_registry(limited, registry_path) as connection:
        index_fasta(connection, fasta, limited)

    changed = settings.model_copy(update={"max_detailed_entries": 3})
    with (
        connect_registry(registry_path) as connection,
        pytest.raises(
            ValueError,
            match=r"max_detailed_entries changed from '2' to '3'.*fasta-gen reindex --full",
        ),
    ):
        initialize_registry(connection, changed)


def test_registry_rejects_changed_file_size_limit_without_full_reindex(
    settings: Settings,
    tmp_path: Path,
) -> None:
    registry_path = tmp_path / "registry.sqlite3"
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    with open_registry(settings, registry_path) as connection:
        index_fasta(connection, fasta, settings)

    changed = settings.model_copy(update={"max_fasta_file_size_gib": 6.0})
    with (
        connect_registry(registry_path) as connection,
        pytest.raises(
            ValueError,
            match=r"max_fasta_file_size_gib changed from '5.0' to '6.0'.*fasta-gen reindex --full",
        ),
    ):
        initialize_registry(connection, changed)


def test_registry_rejects_changed_metadata_aa_sample_size_without_full_reindex(
    settings: Settings,
    tmp_path: Path,
) -> None:
    registry_path = tmp_path / "registry.sqlite3"
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|A|A_TEST\nAAAA\n")
    sampled = settings.model_copy(update={"metadata_aa_sample_size": 2})
    with open_registry(sampled, registry_path) as connection:
        index_fasta(connection, fasta, sampled)

    changed = settings.model_copy(update={"metadata_aa_sample_size": 3})
    with (
        connect_registry(registry_path) as connection,
        pytest.raises(
            ValueError,
            match=r"metadata_aa_sample_size changed from '2' to '3'.*fasta-gen reindex --full",
        ),
    ):
        initialize_registry(connection, changed)


def test_export_registry_stats_only_drops_entries_and_keeps_read_paths(
    settings: Settings,
    tmp_path: Path,
) -> None:
    source = tmp_path / "registry.sqlite3"
    fasta_path = tmp_path / "p1_db1_test_20260721.fasta"
    _write_database(fasta_path, settings)
    with open_registry(settings, source) as connection:
        index_fasta(connection, fasta_path, settings, root=tmp_path)

    destination = tmp_path / "export" / "stats.sqlite3"
    size_bytes, pair_rows = registry.export_registry_stats_only(source, destination)

    assert size_bytes == destination.stat().st_size
    with connect_registry(destination) as exported:
        tables = {
            row[0]
            for row in exported.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        }
        assert "entries" not in tables
        assert STATS_ONLY_TABLES[0] in tables
        assert exported.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION
        assert [record.relative_path for record in list_databases(exported)] == [fasta_path.name]
        assert (
            pair_rows == exported.execute("SELECT count(*) FROM database_pair_stats").fetchone()[0]
        )


def test_export_registry_stats_only_refuses_to_overwrite(
    settings: Settings, tmp_path: Path
) -> None:
    source = tmp_path / "registry.sqlite3"
    with open_registry(settings, source):
        pass
    destination = tmp_path / "stats.sqlite3"
    destination.write_bytes(b"existing")

    with pytest.raises(FileExistsError):
        registry.export_registry_stats_only(source, destination)

    assert destination.read_bytes() == b"existing"


def test_export_registry_stats_only_rejects_a_missing_source(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        registry.export_registry_stats_only(tmp_path / "absent.sqlite3", tmp_path / "stats.sqlite3")


def test_export_registry_stats_only_leaves_no_partial_file(
    settings: Settings, tmp_path: Path
) -> None:
    source = tmp_path / "registry.sqlite3"
    with open_registry(settings, source):
        pass
    # A schema version the exporter must refuse, so the failure happens after
    # the destination path has been accepted.
    with connect_registry(source) as connection:
        stamp_schema_version(connection, 3)
    destination = tmp_path / "stats.sqlite3"

    with pytest.raises(ValueError, match="schema version 3"):
        registry.export_registry_stats_only(source, destination)

    assert not destination.exists()


def _dated_settings(settings: Settings, root: Path, cutoff: datetime.date | None) -> Settings:
    return settings.model_copy(update={"fasta_root": root, "min_fasta_date": cutoff})


def test_update_registry_skips_fastas_built_before_min_fasta_date(
    settings: Settings,
    tmp_path: Path,
) -> None:
    old_path = tmp_path / "p1_db1_test_20090520.fasta"
    new_path = tmp_path / "p1_db1_test_20240115.fasta"
    undated_recent = tmp_path / "recent_no_date.fasta"
    undated_old = tmp_path / "old_no_date.fasta"
    for path in (old_path, new_path, undated_recent, undated_old):
        _write_database(path, settings)
    # No date in the name, so the filter has to fall back to the filesystem.
    old_stamp = datetime.datetime(2009, 5, 20, tzinfo=datetime.UTC).timestamp()
    os.utime(undated_old, (old_stamp, old_stamp))
    dated = _dated_settings(settings, tmp_path, datetime.date(2023, 1, 1))

    records = update_registry(tmp_path, tmp_path / "registry.sqlite3", dated)

    assert sorted(record.relative_path for record in records) == sorted(
        [new_path.name, undated_recent.name]
    )


def test_update_registry_prunes_records_that_the_date_filter_now_excludes(
    settings: Settings,
    tmp_path: Path,
) -> None:
    old_path = tmp_path / "p1_db1_test_20090520.fasta"
    _write_database(old_path, settings)
    registry_path = tmp_path / "registry.sqlite3"
    unfiltered = _dated_settings(settings, tmp_path, None)

    assert [
        record.relative_path for record in update_registry(tmp_path, registry_path, unfiltered)
    ] == [old_path.name]

    filtered = _dated_settings(settings, tmp_path, datetime.date(2023, 1, 1))
    assert update_registry(tmp_path, registry_path, filtered) == []


def test_rebuild_registry_applies_min_fasta_date(settings: Settings, tmp_path: Path) -> None:
    old_path = tmp_path / "p1_db1_test_20090520.fasta"
    new_path = tmp_path / "p1_db1_test_20240115.fasta"
    for path in (old_path, new_path):
        _write_database(path, settings)
    dated = _dated_settings(settings, tmp_path, datetime.date(2023, 1, 1))

    records = rebuild_registry(tmp_path, tmp_path / "registry.sqlite3", dated)

    assert [record.relative_path for record in records] == [new_path.name]


def test_min_fasta_date_none_indexes_every_date(settings: Settings, tmp_path: Path) -> None:
    old_path = tmp_path / "p1_db1_test_20090520.fasta"
    _write_database(old_path, settings)

    records = update_registry(
        tmp_path, tmp_path / "registry.sqlite3", _dated_settings(settings, tmp_path, None)
    )

    assert [record.relative_path for record in records] == [old_path.name]


def test_registry_counts_entrapment_records_separately_from_targets(tmp_path: Path) -> None:
    """They were classified as targets, so entrapment_count was always zero."""
    settings = Settings(fasta_root=tmp_path / "databases", registry_dir=tmp_path / "registry")
    settings.fasta_root.mkdir(parents=True)
    path = settings.fasta_root / "bench_e_20260813.fasta"
    path.write_text(
        ">sp|P1|ONE first\nMPEPTIDEKTESTSEQUENCEK\n"
        ">sp|P1|ONE_p_target entrapment of sp|P1|ONE\nMEDIPEPTKTEEEESQNSUTCK\n"
        ">REV_sp|P1|ONE first\nKECNEUQESTSETKEDITPEPM\n"
    )
    with open_test_registry(settings) as connection:
        index_fasta(connection, path, settings, root=settings.fasta_root)
        record = list_databases(connection)[0]

    assert record.target_count == 1
    assert record.entrapment_count == 1
    assert record.decoy_count == 1
