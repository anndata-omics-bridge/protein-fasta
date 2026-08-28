"""Dated registries of two engines sharing one directory.

The suffix is what keeps them apart. These tests exist because getting it wrong
is not a cosmetic bug: an unscoped retention count would count another engine's
live registries as superseded copies of this one's and delete them.
"""

# pyright: basic
from __future__ import annotations

import datetime
from pathlib import Path

import pytest
from tests.registry_support import BackendSettings, Settings

from protein_fasta.registry.backend import factory
from protein_fasta.registry.snapshots import (
    latest_snapshot,
    list_snapshots,
    new_snapshot_path,
    require_latest_snapshot,
    snapshot_name,
    superseded_snapshots,
)

SQLITE = factory.suffix_for("sqlite")
DUCKDB = factory.suffix_for("duckdb")


def _publish(directory: Path, *, suffix: str, count: int) -> list[Path]:
    """Write ``count`` dated registries a second apart, oldest first."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for offset in range(count):
        moment = datetime.datetime(2026, 1, 1, tzinfo=datetime.UTC) + datetime.timedelta(
            seconds=offset
        )
        path = directory / snapshot_name(moment, suffix=suffix)
        path.write_bytes(b"")
        paths.append(path)
    return paths


def test_each_engine_sees_only_its_own_registries(tmp_path: Path) -> None:
    """Two backends publish into one directory without seeing each other's files."""
    sqlite_files = _publish(tmp_path, suffix=SQLITE, count=2)
    duckdb_files = _publish(tmp_path, suffix=DUCKDB, count=3)

    assert list_snapshots(tmp_path, suffix=SQLITE) == sqlite_files
    assert list_snapshots(tmp_path, suffix=DUCKDB) == duckdb_files
    assert latest_snapshot(tmp_path, suffix=SQLITE) == sqlite_files[-1]
    assert latest_snapshot(tmp_path, suffix=DUCKDB) == duckdb_files[-1]


def test_retention_never_counts_the_other_engines_registries(tmp_path: Path) -> None:
    """The bug this guards: pruning must not delete a live registry of the other engine.

    Three of each, keeping two. Unscoped, the six would be ranked together and
    four "superseded" files named -- including DuckDB registries that are the
    newest of their own kind.
    """
    sqlite_files = _publish(tmp_path, suffix=SQLITE, count=3)
    duckdb_files = _publish(tmp_path, suffix=DUCKDB, count=3)

    assert superseded_snapshots(tmp_path, keep=2, suffix=SQLITE) == sqlite_files[:1]
    assert superseded_snapshots(tmp_path, keep=2, suffix=DUCKDB) == duckdb_files[:1]


def test_sidecars_are_never_mistaken_for_registries(tmp_path: Path) -> None:
    """Logs, journals, and write-ahead files sit beside a registry, not among them."""
    published = _publish(tmp_path, suffix=SQLITE, count=1)[0]
    for sidecar in (
        f"{published.name}.log",
        f"{published.name}-journal",
        f".{published.name}.tmp",
    ):
        (tmp_path / sidecar).write_bytes(b"")
    duckdb = _publish(tmp_path, suffix=DUCKDB, count=1)[0]
    (tmp_path / f"{duckdb.name}.wal").write_bytes(b"")

    assert list_snapshots(tmp_path, suffix=SQLITE) == [published]
    assert list_snapshots(tmp_path, suffix=DUCKDB) == [duckdb]


def test_a_new_path_never_lands_on_an_existing_registry(tmp_path: Path) -> None:
    """Two rebuilds in the same second get distinct names, per engine."""
    moment = datetime.datetime(2026, 1, 1, tzinfo=datetime.UTC)
    first = new_snapshot_path(tmp_path, suffix=SQLITE, moment=moment)
    first.parent.mkdir(parents=True, exist_ok=True)
    first.write_bytes(b"")
    second = new_snapshot_path(tmp_path, suffix=SQLITE, moment=moment)

    assert second != first
    # A DuckDB rebuild in the same second is not blocked by the SQLite file.
    assert new_snapshot_path(tmp_path, suffix=DUCKDB, moment=moment).name.endswith(DUCKDB)


def test_switching_backend_before_a_rebuild_reports_no_registry(tmp_path: Path) -> None:
    """Configuring an engine that has published nothing yet must say so.

    Reading the other engine's file instead would be worse: it would either fail
    deep inside a driver or, if it happened to parse, report the wrong registry.
    """
    _publish(tmp_path, suffix=SQLITE, count=1)
    settings = Settings(
        fasta_root=tmp_path / "databases",
        registry_dir=tmp_path,
        registry=BackendSettings(backend="duckdb"),
    )

    with pytest.raises(FileNotFoundError, match="reindex --full"):
        require_latest_snapshot(settings.registry_dir, backend=settings.registry.backend)


def test_the_engine_that_wrote_a_file_is_read_off_its_name(tmp_path: Path) -> None:
    """A path names its engine, so a registry stays readable when the default changes."""
    assert (
        factory.backend_for_path(tmp_path / f"fasta_registry-20260101T000000Z{SQLITE}") == "sqlite"
    )
    assert (
        factory.backend_for_path(tmp_path / f"fasta_registry-20260101T000000Z{DUCKDB}") == "duckdb"
    )

    with pytest.raises(ValueError, match="Cannot tell which engine"):
        factory.backend_for_path(tmp_path / "registry.db")
