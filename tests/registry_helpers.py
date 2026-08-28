"""Helpers for tests that work with one registry inside a registry directory."""

from __future__ import annotations

import datetime
import sqlite3
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

from protein_fasta.registry.backend import factory
from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.backend.sqlite import SqliteRegistryConnection
from protein_fasta.registry.indexing import connect_registry, initialize_registry
from protein_fasta.registry.snapshots import latest_snapshot, snapshot_name
from tests.registry_support import Settings

# One fixed instant, so a test's registry has a stable name it can reopen. Production
# stamps each rebuild with the real time; a test only needs the newest to be findable.
_FIXED_MOMENT = datetime.datetime(2026, 1, 1, tzinfo=datetime.UTC)


def registry_file(settings: Settings) -> Path:
    """Return the dated registry a test works with, creating its directory.

    The newest published one when the code under test already wrote a registry --
    that carries the real build instant, not this fixed one -- and otherwise a name
    to create. Either way it is a name ``open_registry`` discovers, so tests keep
    exercising discovery rather than bypassing it.
    """
    settings.registry_dir.mkdir(parents=True, exist_ok=True)
    suffix = factory.suffix_for(settings.registry.backend)
    return latest_snapshot(
        settings.registry_dir, suffix=suffix
    ) or settings.registry_dir / snapshot_name(_FIXED_MOMENT, suffix=suffix)


def in_memory_registry() -> RegistryConnection:
    """Return an empty in-memory registry for tests that build rows by hand.

    These tests predate the backend seam and used a bare driver connection. They
    need a registry connection now, but not a file: the point of them is the SQL,
    and an in-memory database keeps the suite fast.
    """
    return SqliteRegistryConnection(Path(":memory:"))


def execute_script(connection: RegistryConnection, script: str) -> None:
    """Run a semicolon-separated DDL script one statement at a time.

    The registry connection deliberately has no ``executescript``: nothing in the
    package needed multi-statement execution, and losing it is what let the schema
    be inspected and reused a table at a time. Hand-built test fixtures still find
    a script the clearest way to spell a throwaway schema.
    """
    for statement in script.split(";"):
        if statement.strip():
            connection.execute(statement)


def stamp_schema_version(connection: RegistryConnection, version: int) -> None:
    """Record a schema version the way the registry itself records it.

    Tests that simulate a registry from another release used to set
    ``PRAGMA user_version`` alone. The registry reads the portable
    ``registry_meta`` row now, so a fixture writing only the pragma stamps
    nothing the code under test will look at -- and the pragma does not exist on
    every engine. Going through the connection keeps the fixture faithful to a
    real file on either.
    """
    connection.create_tables(("registry_meta",))
    connection.set_schema_version(version)
    connection.commit()


def stamp_sqlite_schema_version(connection: sqlite3.Connection, version: int) -> None:
    """Stamp a version into a hand-built SQLite file, driver connection and all.

    For the tests that simulate a registry left behind by another release by
    writing the file themselves. Those are SQLite-specific by construction -- a
    ``.sqlite3`` on disk is what they are asserting about -- so they stay on the
    driver rather than pretending to be backend-neutral.
    """
    connection.execute(
        "CREATE TABLE IF NOT EXISTS registry_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID"
    )
    connection.execute(
        "INSERT INTO registry_meta (key, value) VALUES ('schema_version', ?) "
        "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
        (str(version),),
    )
    connection.execute(f"PRAGMA user_version = {version}")
    connection.commit()


@contextmanager
def open_test_registry(settings: Settings) -> Generator[RegistryConnection]:
    """Open a test's registry, creating it when it does not exist yet.

    ``open_registry`` deliberately refuses to create one, so that a deployment with
    no registry reports the fact instead of serving an empty database list. A test
    that is seeding the fixture needs the opposite, and gets it here.
    """
    with connect_registry(registry_file(settings)) as connection:
        initialize_registry(connection, settings)
        yield connection
