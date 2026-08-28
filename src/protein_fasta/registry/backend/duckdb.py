"""The DuckDB registry: same logical schema, different physics.

Three measurements shaped everything here, taken on the 456,826-entry dev
registry (see ``TODO/PLAN_registry_backend_protocol.md`` §2):

* Row-wise ``executemany`` runs at 1.8k rows/s against SQLite's 586k. It is
  never used for ``entries``; a columnar batch is, at 0.9-1.9M rows/s.
* A primary key on ``entries`` costs 18x on ingest, because it builds an index
  over every row. There is none, and the guarantee moves to the integrity check.
* The four pair-metric self-joins run 23.5x faster with no secondary indexes at
  all, so this backend creates none. That also deletes the entry-index phase of a
  full rebuild, which was 272.7s on the production collection.

What is *not* different: the tables, their columns, the CHECK constraints that
carry meaning, and every fingerprint computed from them. Cross-backend
fingerprint equality is asserted by the test suite, because the fingerprints hash
rows in query order and a different collation would silently change them.
"""

from __future__ import annotations

import re
import time
from collections.abc import Generator, Iterable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

import duckdb
import polars as pl

from protein_fasta.registry.backend import schema
from protein_fasta.registry.backend.base import (
    Cursor,
    RegistryBackendError,
    RegistryIntegrityError,
    Row,
    TempTableSpec,
)

NAME = "duckdb"
SUFFIX = ".duckdb"

_ENTRY_COLUMNS = (
    "database_id",
    "ordinal",
    "sequence_id",
    "kind",
    "contaminant_group",
    "sequence_length",
    "sequence_hash",
    "description_hash",
)

# Declarations SQLite understands and DuckDB does not, or does not need.
_WITHOUT_ROWID = re.compile(r"\s*\)\s*WITHOUT ROWID\s*$", re.IGNORECASE)
_COLLATE_BINARY = re.compile(r"\s+COLLATE BINARY\b", re.IGNORECASE)
_CASCADE = re.compile(r"\s+REFERENCES\s+databases\(id\)\s+ON DELETE CASCADE", re.IGNORECASE)
_ENTRIES_PRIMARY_KEY = re.compile(r",\s*\n\s*PRIMARY KEY \(database_id, ordinal\)", re.IGNORECASE)
# Both numeric types are narrower here than in SQLite; widening keeps the columns
# interchangeable and the statistics comparable.
_INTEGER = re.compile(r"\bINTEGER\b")
# REAL is single precision here and double in SQLite, which showed up as a mean
# length of 6.599999904632568 where the other engine stored 6.6.
_REAL = re.compile(r"\bREAL\b")
# length() is text-only here; the blob checks mean byte count either way.
_BLOB_LENGTH = re.compile(r"\blength\((sequence_hash|description_hash)\)", re.IGNORECASE)


def _portable(statement: str, *, table: str = "") -> str:
    """Rewrite one table's DDL into what DuckDB accepts, and only that.

    Every removal here is a physical decision, never a semantic one:

    ``WITHOUT ROWID`` has no equivalent. The numeric types are both narrower than
    SQLite's: ``INTEGER`` is 32 bits against 64, which ``mtime_ns`` overflows, and
    ``REAL`` is single precision against double, which cost a length mean its last
    digits. Both widen, so the columns stay interchangeable.
    ``length()`` is text-only, so the two
    checksum checks say ``octet_length`` and mean exactly the same byte count.
    ``COLLATE BINARY`` is refused outright
    by the parser, and DuckDB's default text ordering is already byte-wise, which
    is what the clause asks for. ``ON DELETE CASCADE`` is not implemented, and
    declaring the foreign key without it would put an index on the referencing
    column of the largest table -- so the cascade is done explicitly in
    ``delete_database`` and the reference itself is checked as a query. The
    ``entries`` primary key goes for the same reason, its uniqueness holding by
    construction because ``ordinal`` comes from a counter.
    """
    portable = _WITHOUT_ROWID.sub(")", statement.strip())
    portable = _COLLATE_BINARY.sub("", portable)
    portable = _CASCADE.sub("", portable)
    portable = _BLOB_LENGTH.sub(r"octet_length(\1)", portable)
    portable = _INTEGER.sub("BIGINT", portable)
    portable = _REAL.sub("DOUBLE", portable)
    if table == "entries":
        portable = _ENTRIES_PRIMARY_KEY.sub("", portable)
    return portable


class DuckdbRegistryConnection:
    """One open DuckDB registry.

    Not thread-safe, like every registry connection: the package opens one per
    operation and closes it, which is also what DuckDB's process-wide file lock
    requires.
    """

    entry_batch_size = 1_000_000
    """Effectively one batch per database, because building the column vectors is
    the per-batch cost. At SQLite's 2,000 the same load measured 8x slower than
    SQLite; in one batch it is faster."""

    __slots__ = ("_connection", "path")

    def __init__(self, path: Path, *, read_only: bool = False) -> None:
        self.path = path
        if not read_only:
            path.parent.mkdir(parents=True, exist_ok=True)
        self._connection = _connect_waiting_for_the_lock(path, read_only=read_only)

    @property
    def raw(self) -> duckdb.DuckDBPyConnection:
        """The driver connection, for the few places that need the real thing."""
        return self._connection

    # --- statements ---------------------------------------------------------
    def execute(self, sql: str, parameters: Sequence[object] = ()) -> Cursor:
        """Run one statement on this connection and materialize its rows.

        Both halves matter. It must be *this* connection, because a DuckDB cursor
        is a separate connection with its own temporary schema and its own
        transaction, so a staging table or a ``BEGIN`` issued on one would vanish
        with it. And the rows must be read out now, because the next statement on
        a connection invalidates the previous result -- which is exactly what
        ``for row in connection.execute(...)`` around another query would hit.
        """
        try:
            self._connection.execute(sql, list(parameters))
            description = self._connection.description
            if not description:
                return _DuckdbCursor([], {})
            columns = {column[0]: index for index, column in enumerate(description)}
            return _DuckdbCursor(self._connection.fetchall(), columns)
        except duckdb.Error as error:
            raise RegistryBackendError(str(error)) from error

    def executemany(self, sql: str, parameters: Iterable[Sequence[object]]) -> None:
        rows = [list(row) for row in parameters]
        if not rows:
            return
        try:
            self._connection.executemany(sql, rows)
        except duckdb.Error as error:
            raise RegistryBackendError(str(error)) from error

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> int:
        row = self.execute(sql, parameters).fetchone()
        assert row is not None
        return int(row[0])

    @contextmanager
    def transaction(self) -> Generator[None]:
        """Open one explicitly, because this engine otherwise autocommits.

        Without the ``BEGIN`` each statement would already be durable, so a
        failure part-way through would leave exactly the half-written registry
        the transaction exists to prevent.
        """
        self.execute("BEGIN TRANSACTION")
        try:
            yield
        except BaseException:
            self.rollback()
            raise
        self.commit()

    def commit(self) -> None:
        try:
            self._connection.commit()
        except duckdb.Error as error:
            # Autocommit means callers reach here with nothing open, having
            # already made their work durable. That is a no-op, not a failure.
            if "no transaction is active" in str(error):
                return
            raise RegistryBackendError(str(error)) from error

    def rollback(self) -> None:
        try:
            self._connection.rollback()
        except duckdb.Error as error:
            if "no transaction is active" not in str(error):
                raise RegistryBackendError(str(error)) from error

    def close(self) -> None:
        self._connection.close()

    # --- schema -------------------------------------------------------------
    def portable_ddl(self, statement: str, *, table: str = "") -> str:
        """Return one SQLite-dialect CREATE statement as this engine accepts it."""
        return _portable(statement, table=table)

    def create_tables(self, tables: Iterable[str] | None = None) -> None:
        for name in tables if tables is not None else schema.REGISTRY_TABLES:
            self.execute(self.portable_ddl(schema.REGISTRY_TABLES[name], table=name))
        if tables is None or "databases" in tables:
            # databases.id is a rowid alias on SQLite. Here it needs a sequence,
            # and insert_database reads the value back with RETURNING.
            self.execute("CREATE SEQUENCE IF NOT EXISTS registry_databases_id START 1")

    def create_entry_indexes(self) -> None:
        """Create nothing, deliberately.

        The queries these indexes served on SQLite -- the pair-metric self-joins
        and the per-kind COUNT(DISTINCT) battery -- measured 23.5x faster here
        without them, because they become hash joins over compressed columns
        rather than B-tree probes. Every index would also slow the load that
        precedes this call and inflate the file it writes.
        """

    def create_pair_indexes(self) -> None:
        """Create nothing: the pair table is small and read by its primary key."""

    def create_stats_only_indexes(self) -> None:
        """Create nothing, for the same reason as the other two."""

    def has_table(self, name: str) -> bool:
        return (
            self.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_name = ? AND table_schema = 'main'",
                (name,),
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

    def upsert_meta(self, key: str, value: str) -> None:
        self.execute(
            "INSERT INTO registry_meta (key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
            (key, value),
        )

    # --- tuning -------------------------------------------------------------
    def configure_bulk(self) -> None:
        # Insertion order is not meaningful: every query orders explicitly, and
        # the fingerprints order by value. Spill beside the registry rather than
        # in the system temp directory, where tens of gigabytes would not fit.
        self.execute("SET preserve_insertion_order = false")
        self.execute(f"SET temp_directory = '{self.path.parent}'")

    def configure_low_memory(self) -> None:
        # Spill beside the workspace file rather than into the system temp
        # directory, which is where a digest of a large proteome would not fit.
        self.execute("SET preserve_insertion_order = false")
        self.execute(f"SET temp_directory = '{self.path.parent}'")

    def analyze(self, table: str | None = None) -> None:
        """Do nothing: DuckDB keeps its own statistics as data is written."""

    def index_hint(self, index: str) -> str:
        """Return nothing: there are no indexes to pin and no syntax to pin them."""
        return ""

    def binary_collation(self) -> str:
        """Return nothing: the parser refuses COLLATE BINARY and the default is byte-wise.

        Load-bearing, and verified rather than assumed: the identifier and content
        fingerprints hash rows in query order, and the suite asserts they match
        SQLite's for every example database.
        """
        return ""

    # --- staging ------------------------------------------------------------
    def temp(self, table: str) -> str:
        return f"temp.main.{table}"

    def create_temp_table(self, spec: TempTableSpec) -> None:
        self.drop_temp_table(spec.name)
        columns = ", ".join(f"{name} {declaration}" for name, declaration in spec.columns)
        # The key is declared only where it means something: INSERT OR IGNORE
        # needs one to ignore against. Where rows are unique by construction it
        # would cost more than the single scan it could accelerate, and the
        # indexes never earn their keep on a table written and read once.
        key = (
            f", PRIMARY KEY ({', '.join(spec.primary_key)})"
            if spec.deduplicating and spec.primary_key
            else ""
        )
        self.execute(f"CREATE TEMP TABLE {spec.name} ({columns}{key})")

    def drop_temp_table(self, table: str) -> None:
        self.execute(f"DROP TABLE IF EXISTS {self.temp(table)}")

    # --- writes -------------------------------------------------------------
    def insert_database(self, columns: str, values: Sequence[object]) -> int:
        placeholders = ",".join("?" for _ in values)
        row = self.execute(
            f"INSERT INTO databases (id, {columns}) "
            f"VALUES (nextval('registry_databases_id'), {placeholders}) RETURNING id",
            values,
        ).fetchone()
        if row is None:
            raise RegistryBackendError("The registry did not return an ID for a new database row.")
        return int(row[0])

    def insert_entries(self, table: str, rows: Sequence[tuple[object, ...]]) -> None:
        """Append one batch of entry rows as a column batch, never row by row.

        The transpose is the point: DuckDB binds Python rows one value at a time,
        which measured 330x slower than SQLite, while a columnar batch named in
        the statement is 3x faster than SQLite. The frame is referenced by name
        from the SQL below -- DuckDB resolves it from the enclosing scope.
        """
        if not rows:
            return
        entry_batch = pl.DataFrame(
            {name: [row[index] for row in rows] for index, name in enumerate(_ENTRY_COLUMNS)},
            strict=False,
        )
        try:
            self._connection.execute(
                f"INSERT INTO {table} ({', '.join(_ENTRY_COLUMNS)}) SELECT * FROM entry_batch"
            )
        except duckdb.Error as error:
            raise RegistryBackendError(str(error)) from error
        del entry_batch

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
        """Check what this engine can be asked, then check the rest as queries.

        There is no page-level equivalent of ``PRAGMA integrity_check``: DuckDB
        validates blocks as it reads them, so a corrupt file fails the scan below
        rather than a dedicated check. Weaker than SQLite's, and said plainly
        rather than implied. The relational guarantees are stronger here than a
        foreign-key pragma, because this backend declares neither the references
        nor the ``entries`` key and therefore checks the data itself.
        """
        self.scalar("SELECT COUNT(*) FROM entries")
        for table, column in schema.CHILD_TABLES:
            orphans = self.scalar(
                f"SELECT COUNT(*) FROM {table} LEFT JOIN databases ON databases.id = {table}.{column} "
                "WHERE databases.id IS NULL"
            )
            if orphans:
                raise RegistryIntegrityError(
                    f"{table}.{column} has {orphans} row(s) with no database."
                )
        for table, key in schema.UNIQUE_KEYS:
            duplicates = self.scalar(
                f"SELECT COUNT(*) FROM (SELECT {key} FROM {table} GROUP BY {key} HAVING COUNT(*) > 1)"
            )
            if duplicates:
                raise RegistryIntegrityError(
                    f"{table} has {duplicates} duplicated key(s) on ({key})."
                )


class _DuckdbCursor:
    """One statement's rows, already read out and addressable by position or name."""

    __slots__ = ("_columns", "_rows")

    def __init__(self, rows: list[tuple[object, ...]], columns: dict[str, int]) -> None:
        self._rows = rows
        self._columns = columns

    def __iter__(self) -> Iterator[Row]:
        return (_DuckdbRow(row, self._columns) for row in self._rows)

    def fetchone(self) -> Row | None:
        return _DuckdbRow(self._rows[0], self._columns) if self._rows else None

    def fetchall(self) -> list[Row]:
        return [_DuckdbRow(row, self._columns) for row in self._rows]


class _DuckdbRow:
    """A row that answers to a column position and to a column name.

    Both are used in this package, sometimes on the same query: the aggregate
    queries index positionally while the record builders read by name.
    """

    __slots__ = ("_columns", "_values")

    def __init__(self, values: tuple[object, ...], columns: dict[str, int]) -> None:
        self._values = values
        self._columns = columns

    def __getitem__(self, key: int | str) -> object:
        return self._values[key] if isinstance(key, int) else self._values[self._columns[key]]

    def __len__(self) -> int:
        return len(self._values)

    def __iter__(self) -> Iterator[object]:
        return iter(self._values)

    def keys(self) -> list[str]:
        return list(self._columns)


_LOCK_TIMEOUT_SECONDS = 30.0
"""How long to wait for another process to release the file.

Stands in for SQLite's ``PRAGMA busy_timeout = 30000``, which is what this engine
has no equivalent of. The operations that hold the file are short -- a page render
reads a few megabytes -- so a sweep that starts mid-render should wait rather than
fail, while a sweep started against a registry someone is rebuilding should give up
and say so.
"""
_LOCK_RETRY_SECONDS = 0.25


def _is_lock_conflict(error: BaseException) -> bool:
    """Report whether an open failed because another process holds the file."""
    return "lock" in str(error).lower()


def _connect_waiting_for_the_lock(path: Path, *, read_only: bool) -> duckdb.DuckDBPyConnection:
    """Open the file, waiting a bounded time for a conflicting lock to clear."""
    deadline = time.monotonic() + _LOCK_TIMEOUT_SECONDS
    while True:
        try:
            return duckdb.connect(str(path), read_only=read_only)
        except duckdb.Error as error:
            if not _is_lock_conflict(error) or time.monotonic() >= deadline:
                raise RegistryBackendError(_open_failure(path, read_only, error)) from error
            time.sleep(_LOCK_RETRY_SECONDS)


def _open_failure(path: Path, read_only: bool, error: duckdb.Error) -> str:
    """Explain an open failure, naming the lock when that is what it is."""
    if _is_lock_conflict(error):
        return (
            f"Cannot open registry {path.name}: another process still held it after "
            f"{_LOCK_TIMEOUT_SECONDS:.0f}s. This engine allows one writer per file, so an "
            "incremental 'fasta-gen reindex' cannot run while something else is writing; "
            "'reindex --full' publishes a new file and is unaffected."
        )
    return f"Cannot open registry {path}{' read-only' if read_only else ''}: {error}"


def copy_stats_only(source_path: Path, destination: Path, schema_version: int) -> int:
    """Copy the read-path tables into a new registry, returning the pair-row count."""
    target = DuckdbRegistryConnection(destination)
    try:
        target.execute(f"ATTACH '{source_path.resolve()}' AS source (READ_ONLY)")
        target.create_tables(schema.STATS_ONLY_TABLES)
        for table_name in schema.STATS_ONLY_TABLES:
            target.execute(f"INSERT INTO main.{table_name} SELECT * FROM source.main.{table_name}")
        target.commit()
        target.execute("DETACH source")
        target.set_schema_version(schema_version)
        target.commit()
        return target.scalar("SELECT count(*) FROM database_pair_stats")
    finally:
        target.close()


__all__ = ["NAME", "SUFFIX", "DuckdbRegistryConnection", "copy_stats_only"]
