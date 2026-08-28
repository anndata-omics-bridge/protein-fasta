"""The registry of storage engines, and which one owns a given file.

A registry's engine is a property of the file, not of the invocation. A
``--registry-path`` naming a ``.sqlite3`` has to be read by SQLite even when the
configured default is something else, or a perfectly good registry gets reported
as corrupt. So reading dispatches on the suffix and only creation consults
configuration.

Engines are looked up in a table rather than selected by a chain of comparisons,
so adding one is registering it and nothing else. Everything that needs to know
what exists -- the config's accepted values, the snapshot suffixes, the error
messages -- reads that table.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from protein_fasta.registry.backend import sqlite as sqlite_backend
from protein_fasta.registry.backend.base import RegistryBackendError, RegistryConnection


class _ConnectionFactory(Protocol):
    """Open one registry file, read-write unless asked otherwise."""

    def __call__(self, path: Path, *, read_only: bool = False) -> RegistryConnection: ...


class _StatsCopier(Protocol):
    """Copy one registry's read-path tables into a new file, returning pair rows."""

    def __call__(self, source_path: Path, destination: Path, schema_version: int) -> int: ...


@dataclass(frozen=True, slots=True)
class Backend:
    """One registered storage engine.

    ``suffix`` is what makes engines coexist in a directory: it names the files a
    backend owns, so discovery, retention, and dispatch all agree without anyone
    consulting configuration.
    """

    name: str
    suffix: str
    connect: _ConnectionFactory
    copy_stats_only: _StatsCopier


def _sqlite() -> Backend:
    return Backend(
        name=sqlite_backend.NAME,
        suffix=sqlite_backend.SUFFIX,
        connect=sqlite_backend.SqliteRegistryConnection,
        copy_stats_only=sqlite_backend.copy_stats_only,
    )


def _duckdb() -> Backend:
    from protein_fasta.registry.backend import duckdb as duckdb_backend

    return Backend(
        name=duckdb_backend.NAME,
        suffix=duckdb_backend.SUFFIX,
        connect=duckdb_backend.DuckdbRegistryConnection,
        copy_stats_only=duckdb_backend.copy_stats_only,
    )


_BACKEND_LOADERS: dict[str, Callable[[], Backend]] = {
    "sqlite": _sqlite,
    "duckdb": _duckdb,
}
"""Lazy storage-engine constructors, keyed by configured name."""

SUFFIXES: dict[str, str] = {"sqlite": ".sqlite3", "duckdb": ".duckdb"}
"""Registry filename suffix per backend, so several can share one registry_dir."""


def backend_named(name: str) -> Backend:
    """Return one registered engine, or say which names exist."""
    try:
        loader = _BACKEND_LOADERS[name]
    except KeyError:
        expected = ", ".join(sorted(_BACKEND_LOADERS))
        raise ValueError(
            f"Unknown registry backend {name!r}; expected one of {expected}."
        ) from None
    return loader()


def suffix_for(backend: str) -> str:
    """Return the filename suffix a backend publishes registries under."""
    try:
        return SUFFIXES[backend]
    except KeyError:
        expected = ", ".join(sorted(SUFFIXES))
        raise ValueError(
            f"Unknown registry backend {backend!r}; expected one of {expected}."
        ) from None


def backend_for_path(path: Path) -> str:
    """Return the backend that owns one registry file, from its suffix.

    An unrecognised suffix is refused rather than guessed, because guessing wrong
    reports a readable registry as unreadable.
    """
    for backend, suffix in SUFFIXES.items():
        if path.name.endswith(suffix):
            return backend
    raise ValueError(
        f"Cannot tell which engine wrote {path.name}; a registry filename ends in "
        f"{' or '.join(sorted(SUFFIXES.values()))}."
    )


@contextmanager
def connect(
    path: Path, *, backend: str | None = None, read_only: bool = False
) -> Generator[RegistryConnection]:
    """Open one registry for an operation or callback, and close it after.

    Without ``backend`` the engine is taken from the filename, which is what
    every read path wants. Creation passes it explicitly.
    """
    resolved = backend if backend is not None else backend_for_path(path)
    connection = backend_named(resolved).connect(path, read_only=read_only)
    try:
        yield connection
    finally:
        connection.close()


def copy_stats_only(source_path: Path, destination: Path, schema_version: int) -> int:
    """Copy one registry's read-path tables into a new file of the same engine.

    Dispatched on the source, and the destination has to match it: this is a
    physical copy, not a migration between engines.
    """
    backend = backend_for_path(source_path)
    if backend != backend_for_path(destination):
        raise ValueError(
            f"A stats-only copy cannot change engine: {source_path.name} and {destination.name} must share a suffix."
        )
    return backend_named(backend).copy_stats_only(source_path, destination, schema_version)


__all__ = [
    "SUFFIXES",
    "Backend",
    "RegistryBackendError",
    "RegistryConnection",
    "backend_for_path",
    "backend_named",
    "connect",
    "copy_stats_only",
    "suffix_for",
]
