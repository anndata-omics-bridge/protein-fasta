"""Backend-neutral access to one registry file.

The registry's queries are portable; its schema and session handling are not.
This protocol is drawn on exactly that line: every statement both engines accept
verbatim stays as shared SQL in the calling module, and only what cannot be
written once becomes a method here -- DDL, session tuning, the schema version,
identity, transaction boundaries, cascade deletes, staging tables, bulk insert,
and physical integrity.

The alternative was a wide repository of domain operations, one implementation
per engine. That would put the four pair-metric self-joins in two places, and
those counts are what every comparison view reports; two copies would eventually
disagree about a ``COUNT(DISTINCT description_hash)`` edge and no test would say
which was right.

The operations are methods on the connection rather than a separate dialect
object because every function in the package already takes ``connection`` first.
A second object would have to be threaded through some sixty signatures and
every GUI call site to reach the same places this reaches for free.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from contextlib import AbstractContextManager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol


class RegistryBackendError(RuntimeError):
    """The storage engine refused a statement, or could not open the file.

    Exists so nothing outside this package has to name a driver to catch a broken
    registry. Both implementations translate their driver's error hierarchy into
    this one at ``execute``, which is the only place a driver's exceptions can
    surface.
    """


class RegistryIntegrityError(RuntimeError):
    """The registry is physically damaged or has incomplete materialized data."""


class Row(Protocol):
    """One result row addressed either by column position or by column name.

    Both are load-bearing in existing code, sometimes on the same query: the
    aggregate queries read ``row[0]`` while the record builders read
    ``row["relative_path"]``, and one comparison query indexes positionally into
    a named-column result. Supporting both is what lets the row-to-model
    functions stay as they are.
    """

    def __getitem__(self, key: int | str) -> Any: ...  # noqa: ANN401

    def __len__(self) -> int: ...


class Cursor(Protocol):
    """One executed statement's results."""

    def __iter__(self) -> Iterator[Row]: ...

    def fetchone(self) -> Row | None: ...

    def fetchall(self) -> list[Row]: ...


@dataclass(frozen=True, slots=True)
class TempTableSpec:
    """One connection-local staging table, described rather than spelled in SQL.

    The columns are named by the shared code that stages into them; how they are
    declared is the backend's decision. That matters because a primary key or an
    index on a staging table is a win on a B-tree engine and a cost on one that
    hash-joins the scan it would have accelerated.
    """

    name: str
    columns: tuple[tuple[str, str], ...]
    primary_key: tuple[str, ...] = ()
    indexes: tuple[tuple[str, ...], ...] = ()
    partial_index_column: str | None = None
    deduplicating: bool = False
    """Whether the key carries meaning rather than speed.

    ``INSERT OR IGNORE`` is how the staging tables drop repeats, and it needs the
    key to exist -- so on those the key is semantics and every backend declares
    it. Where rows are unique by construction the key is an index like any other,
    and a backend that gains nothing from one may skip it.
    """


class RegistryConnection(Protocol):
    """One open registry, and the engine-specific operations it needs.

    Obtained and closed by :func:`protein_fasta.registry.backend.factory.connect`.
    Never cached: the file lock is process-wide on some engines, and the Dash
    background callbacks fork, so a connection that outlives one operation is a
    crash waiting for the right timing.
    """

    path: Path

    entry_batch_size: int
    """How many entry rows to accumulate before writing them.

    A row-binding engine wants small batches so nothing large is held twice; a
    columnar one wants the opposite, because the per-batch cost is building the
    column vectors and a small batch pays it over and over. Measured, the wrong
    choice here is an order of magnitude, not a percentage.
    """

    # --- statements ---------------------------------------------------------
    def execute(self, sql: str, parameters: Sequence[object] = ()) -> Cursor:
        """Run one portable statement, translating driver errors."""
        ...

    def executemany(self, sql: str, parameters: Iterable[Sequence[object]]) -> None:
        """Run one portable statement over many parameter rows."""
        ...

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> int:
        """Return the single integer a counting query produced."""
        ...

    def transaction(self) -> AbstractContextManager[None]:
        """Commit the enclosed work, or roll all of it back."""
        ...

    def commit(self) -> None:
        """Commit the open transaction."""
        ...

    def rollback(self) -> None:
        """Discard the open transaction."""
        ...

    def close(self) -> None:
        """Release the file, and with it any exclusive lock on it."""
        ...

    # --- schema -------------------------------------------------------------
    def create_tables(self, tables: Iterable[str] | None = None) -> None:
        """Create the registry tables, or the named subset of them."""
        ...

    def portable_ddl(self, statement: str, *, table: str = "") -> str:
        """Return one SQLite-dialect CREATE statement as this engine accepts it.

        DDL is written once, in the dialect every published registry was created
        with, and each backend says what it has to change. The alternative is a
        second copy per engine, which drifts: the column order matters, because
        ``SELECT *`` feeds record builders that index positionally.
        """
        ...

    def create_entry_indexes(self) -> None:
        """Create whatever the engine needs to look up and join entries.

        A bulk rebuild calls this after loading, so it is also where a backend
        decides to create nothing: an engine that hash-joins these queries pays
        for every index twice over, once building it and once in the load.
        """
        ...

    def create_pair_indexes(self) -> None:
        """Create whatever the engine needs to read materialized pair rows."""
        ...

    def create_stats_only_indexes(self) -> None:
        """Create the indexes belonging to the tables a stats-only copy carries."""
        ...

    def has_table(self, name: str) -> bool:
        """Report whether one registry table exists yet."""
        ...

    def schema_version(self) -> int:
        """Return the recorded schema version, or 0 when there is no schema."""
        ...

    def set_schema_version(self, version: int) -> None:
        """Record the schema version this release wrote."""
        ...

    def upsert_meta(self, key: str, value: str) -> None:
        """Insert or replace one ``registry_meta`` entry."""
        ...

    # --- tuning -------------------------------------------------------------
    def configure_bulk(self) -> None:
        """Trade durability for throughput on an unpublished temporary file."""
        ...

    def configure_low_memory(self) -> None:
        """Keep a throwaway workspace's working set on disk rather than in memory.

        The opposite trade from :meth:`configure_bulk`, and the reason the
        database-backed peptide workspaces exist at all: they are the fallback
        chosen when a digest will not fit in memory, so spilling is the point
        rather than a cost.
        """
        ...

    def analyze(self, table: str | None = None) -> None:
        """Refresh planner statistics, for engines that need to be told."""
        ...

    def index_hint(self, index: str) -> str:
        """Return the clause pinning one join to an index, or an empty string."""
        ...

    def binary_collation(self) -> str:
        """Return the clause that orders text by bytes, or an empty string."""
        ...

    # --- staging ------------------------------------------------------------
    def temp(self, table: str) -> str:
        """Return the qualified name of one connection-local staging table."""
        ...

    def create_temp_table(self, spec: TempTableSpec) -> None:
        """Replace one connection-local staging table with an empty one."""
        ...

    def drop_temp_table(self, table: str) -> None:
        """Drop one connection-local staging table if it exists."""
        ...

    # --- writes -------------------------------------------------------------
    def insert_database(self, columns: str, values: Sequence[object]) -> int:
        """Insert one database row and return the identifier it was given."""
        ...

    def insert_entries(self, table: str, rows: Sequence[tuple[object, ...]]) -> None:
        """Append one batch of entry detail rows.

        Separate from ``executemany`` because this is the only write whose
        throughput decides how long a full rebuild takes, and the fastest way to
        do it differs by orders of magnitude between engines.
        """
        ...

    def delete_database(self, relative_path: str) -> None:
        """Remove one database and every row that belongs to it."""
        ...

    # --- integrity ----------------------------------------------------------
    def check_physical_integrity(self) -> None:
        """Raise when the file is damaged or its rows contradict each other."""
        ...
