"""The SQLite registry, which is what every published registry is today.

Everything here is the behaviour the registry already had; the value of naming it
as one implementation is that the next one has something to be identical to.
"""

from __future__ import annotations

import contextlib
import sqlite3
from collections.abc import Generator, Iterable, Sequence
from contextlib import contextmanager
from pathlib import Path

from protein_fasta.registry.backend import schema
from protein_fasta.registry.backend.base import (
    Cursor,
    RegistryBackendError,
    RegistryIntegrityError,
    Row,
    TempTableSpec,
)

NAME = "sqlite"
SUFFIX = ".sqlite3"


class SqliteRegistryConnection:
    """One open SQLite registry.

    Rows come back as :class:`sqlite3.Row`, which already answers to both a
    column position and a column name, so the record builders need no adapter.
    """

    entry_batch_size = 2_000
    """Rows per executemany. Small on purpose: this engine binds row by row."""

    __slots__ = ("_connection", "path")

    def __init__(self, path: Path, *, read_only: bool = False, uri: bool = False) -> None:
        self.path = path
        try:
            if read_only:
                connection = sqlite3.connect(
                    f"file:{path.resolve()}?mode=ro", uri=True, timeout=30.0
                )
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                connection = sqlite3.connect(path, timeout=30.0, uri=uri)
            connection.row_factory = sqlite3.Row
            connection.execute("PRAGMA foreign_keys = ON")
            connection.execute("PRAGMA busy_timeout = 30000")
        except sqlite3.Error as error:
            raise RegistryBackendError(f"Cannot open registry {path}: {error}") from error
        self._connection = connection

    @property
    def raw(self) -> sqlite3.Connection:
        """The driver connection, for the few places that still need the real thing."""
        return self._connection

    # --- statements ---------------------------------------------------------
    def execute(self, sql: str, parameters: Sequence[object] = ()) -> Cursor:
        try:
            return self._connection.execute(sql, parameters)
        except sqlite3.Error as error:
            raise RegistryBackendError(str(error)) from error

    def executemany(self, sql: str, parameters: Iterable[Sequence[object]]) -> None:
        try:
            self._connection.executemany(sql, parameters)
        except sqlite3.Error as error:
            raise RegistryBackendError(str(error)) from error

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> int:
        row = self.execute(sql, parameters).fetchone()
        assert row is not None
        return int(row[0])

    @contextmanager
    def transaction(self) -> Generator[None]:
        try:
            yield
        except BaseException:
            self.rollback()
            raise
        self.commit()

    def commit(self) -> None:
        try:
            self._connection.commit()
        except sqlite3.Error as error:
            raise RegistryBackendError(str(error)) from error

    def rollback(self) -> None:
        self._connection.rollback()

    def close(self) -> None:
        self._connection.close()

    # --- schema -------------------------------------------------------------
    def create_tables(self, tables: Iterable[str] | None = None) -> None:
        for name in tables if tables is not None else schema.REGISTRY_TABLES:
            self.execute(self.portable_ddl(schema.REGISTRY_TABLES[name], table=name))

    def portable_ddl(self, statement: str, *, table: str = "") -> str:
        """Return the statement unchanged: it is already written in this dialect."""
        return statement

    def create_entry_indexes(self) -> None:
        for statement in schema.ENTRY_LOOKUP_INDEXES:
            self.execute(statement)

    def create_pair_indexes(self) -> None:
        for statement in schema.PAIR_LOOKUP_INDEXES:
            self.execute(statement)

    def create_stats_only_indexes(self) -> None:
        for statement in schema.STATS_ONLY_INDEXES:
            self.execute(statement)

    def has_table(self, name: str) -> bool:
        return (
            self.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?", (name,)
            ).fetchone()
            is not None
        )

    def schema_version(self) -> int:
        if not self.has_table("registry_meta"):
            return 0
        row = self.execute(
            "SELECT value FROM registry_meta WHERE key = 'schema_version'"
        ).fetchone()
        return int(row[0]) if row is not None else 0

    def set_schema_version(self, version: int) -> None:
        self.upsert_meta("schema_version", str(version))
        # Written as well as the row so a binary from an earlier release, which
        # reads only the pragma, still reports a clean schema error.
        self.execute(f"PRAGMA user_version = {version}")

    def upsert_meta(self, key: str, value: str) -> None:
        self.execute(
            "INSERT INTO registry_meta (key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
            (key, value),
        )

    # --- tuning -------------------------------------------------------------
    def configure_bulk(self) -> None:
        self.execute("PRAGMA journal_mode = MEMORY")
        self.execute("PRAGMA synchronous = OFF")
        self.execute("PRAGMA locking_mode = EXCLUSIVE")

    def configure_low_memory(self) -> None:
        self.execute("PRAGMA journal_mode = OFF")
        self.execute("PRAGMA synchronous = OFF")
        self.execute("PRAGMA temp_store = FILE")

    def analyze(self, table: str | None = None) -> None:
        self.execute("ANALYZE" if table is None else f"ANALYZE {table}")

    def index_hint(self, index: str) -> str:
        return f" INDEXED BY {index}"

    def binary_collation(self) -> str:
        return " COLLATE BINARY"

    # --- staging ------------------------------------------------------------
    def temp(self, table: str) -> str:
        return f"temp.{table}"

    def create_temp_table(self, spec: TempTableSpec) -> None:
        self.drop_temp_table(spec.name)
        columns = ",\n            ".join(
            f"{name} {declaration}" for name, declaration in spec.columns
        )
        key = (
            f",\n            PRIMARY KEY ({', '.join(spec.primary_key)})"
            if spec.primary_key
            else ""
        )
        without_rowid = " WITHOUT ROWID" if spec.primary_key else ""
        self.execute(
            f"CREATE TEMP TABLE {spec.name} (\n            {columns}{key}\n        ){without_rowid}"
        )
        for index in spec.indexes:
            name = f"{spec.name}_{'_'.join(index)}"
            partial = (
                f" WHERE {spec.partial_index_column} IS NOT NULL"
                if spec.partial_index_column is not None and spec.partial_index_column in index
                else ""
            )
            self.execute(f"CREATE INDEX {name} ON {spec.name}({', '.join(index)}){partial}")

    def drop_temp_table(self, table: str) -> None:
        self.execute(f"DROP TABLE IF EXISTS {self.temp(table)}")

    # --- writes -------------------------------------------------------------
    def insert_database(self, columns: str, values: Sequence[object]) -> int:
        placeholders = ",".join("?" for _ in values)
        row = self.execute(
            f"INSERT INTO databases ({columns}) VALUES ({placeholders}) RETURNING id",
            values,
        ).fetchone()
        if row is None:
            raise RegistryBackendError("The registry did not return an ID for a new database row.")
        return int(row[0])

    def insert_entries(self, table: str, rows: Sequence[tuple[object, ...]]) -> None:
        self.executemany(
            f"INSERT INTO {table} "
            "(database_id, ordinal, sequence_id, kind, contaminant_group, sequence_length, "
            "sequence_hash, description_hash) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )

    def delete_database(self, relative_path: str) -> None:
        row = self.execute(
            "SELECT id FROM databases WHERE relative_path = ?", (relative_path,)
        ).fetchone()
        if row is None:
            return
        database_id = int(row[0])
        self.execute(
            "DELETE FROM database_pair_stats WHERE database_id_low = ? OR database_id_high = ?",
            (database_id, database_id),
        )
        self.execute("DELETE FROM database_kind_stats WHERE database_id = ?", (database_id,))
        self.execute("DELETE FROM entries WHERE database_id = ?", (database_id,))
        self.execute("DELETE FROM databases WHERE id = ?", (database_id,))

    # --- integrity ----------------------------------------------------------
    def check_physical_integrity(self) -> None:
        integrity = [str(row[0]) for row in self.execute("PRAGMA integrity_check")]
        if integrity != ["ok"]:
            raise RegistryIntegrityError(f"Registry integrity check failed: {'; '.join(integrity)}")
        foreign_key_errors = self.execute("PRAGMA foreign_key_check").fetchall()
        if foreign_key_errors:
            raise RegistryIntegrityError(
                f"Registry foreign-key check failed with {len(foreign_key_errors)} error(s)."
            )
        _check_relational_integrity(self, RegistryIntegrityError)


def _check_relational_integrity(
    connection: SqliteRegistryConnection, error_type: type[Exception]
) -> None:
    """Assert as queries the guarantees a backend may decline to declare."""
    for table, column in schema.CHILD_TABLES:
        orphans = connection.scalar(
            f"SELECT COUNT(*) FROM {table} LEFT JOIN databases ON databases.id = {table}.{column} "
            "WHERE databases.id IS NULL"
        )
        if orphans:
            raise error_type(f"{table}.{column} has {orphans} row(s) with no database.")
    for table, key in schema.UNIQUE_KEYS:
        duplicates = connection.scalar(
            f"SELECT COUNT(*) FROM (SELECT {key} FROM {table} GROUP BY {key} HAVING COUNT(*) > 1)"
        )
        if duplicates:
            raise error_type(f"{table} has {duplicates} duplicated key(s) on ({key}).")


__all__ = ["NAME", "SUFFIX", "Cursor", "Row", "SqliteRegistryConnection"]


def copy_stats_only(source_path: Path, destination: Path, schema_version: int) -> int:
    """Copy the read-path tables into a new registry, returning the pair-row count.

    The destination's tables come from this module's own DDL rather than from DDL
    harvested out of the source file, so the copy states the schema it is writing
    instead of trusting whatever the source happens to say.
    """
    source_uri = f"file:{source_path.resolve()}?mode=ro"
    # ATTACH resolves a URI only on a connection opened for URIs.
    target = SqliteRegistryConnection(destination, uri=True)
    with contextlib.closing(target.raw):
        target.execute("ATTACH DATABASE ? AS source", (source_uri,))
        target.create_tables(schema.STATS_ONLY_TABLES)
        target.create_stats_only_indexes()
        for table_name in schema.STATS_ONLY_TABLES:
            target.execute(f"INSERT INTO main.{table_name} SELECT * FROM source.{table_name}")
        target.commit()
        target.execute("DETACH DATABASE source")
        target.set_schema_version(schema_version)
        target.analyze()
        target.commit()
        return target.scalar("SELECT count(*) FROM database_pair_stats")
