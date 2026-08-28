"""Dated registry files in one directory, with the newest one loaded.

The same shape as the UniProt proteome cache in :mod:`a dated external-data catalogue`:
``reindex --full`` publishes ``fasta_registry-<UTC stamp>.sqlite3`` and every reader
opens the newest. Overwriting one fixed path instead leaves a running application
serving the replaced, already-unlinked inode until it restarts, and destroys the
previous registry the instant the rename lands.

Unlike a proteome snapshot, a registry can be tens of gigabytes, so nothing here
deletes anything unless a retention count is configured explicitly.
"""

from __future__ import annotations

import datetime
from pathlib import Path

from protein_fasta.registry.backend import factory

PREFIX = "fasta_registry-"
_STAMP = "%Y%m%dT%H%M%SZ"


def snapshot_name(moment: datetime.datetime, *, suffix: str) -> str:
    """Return the registry filename for one UTC instant."""
    return f"{PREFIX}{moment.strftime(_STAMP)}{suffix}"


def list_snapshots(directory: Path, *, suffix: str) -> list[Path]:
    """Return dated registries oldest first (empty when none exist).

    Ordering is lexicographic on the zero-padded UTC stamp, which sorts
    chronologically as text. Matching both ends excludes everything that must not be
    mistaken for a published registry: the in-progress ``.fasta_registry-….tmp``
    dotfile, the ``….sqlite3.log`` beside each one, and any ``….sqlite3-journal``
    or ``….duckdb.wal``.

    The suffix also scopes the listing to one engine, so two backends can publish
    into one directory without either seeing the other's files as its own.
    """
    if not directory.is_dir():
        return []
    return sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.name.startswith(PREFIX) and path.name.endswith(suffix)
    )


def latest_snapshot(directory: Path, *, suffix: str) -> Path | None:
    """Return the newest dated registry, or ``None`` when the directory has none."""
    snapshots = list_snapshots(directory, suffix=suffix)
    return snapshots[-1] if snapshots else None


def new_snapshot_path(
    directory: Path, *, suffix: str, moment: datetime.datetime | None = None
) -> Path:
    """Return an unused path for a fresh full rebuild to publish to.

    The stamp resolves to the second, so two rebuilds starting within the same
    second would otherwise land on one name and the second would overwrite a file
    readers are already following -- the very thing dated registries prevent.
    Advancing a second at a time keeps the names sortable and chronological.
    """
    stamp = moment or datetime.datetime.now(datetime.UTC)
    path = directory / snapshot_name(stamp, suffix=suffix)
    while path.exists():
        stamp += datetime.timedelta(seconds=1)
        path = directory / snapshot_name(stamp, suffix=suffix)
    return path


def require_latest_snapshot(directory: Path, *, backend: str) -> Path:
    """Return the newest registry, or explain how to build one.

    Discovery legitimately finds nothing on a new deployment, so the message names
    the directory that was searched and the command that fills it rather than
    letting an empty registry be created and reported as "no databases".
    """
    snapshot = latest_snapshot(directory, suffix=factory.suffix_for(backend))
    if snapshot is None:
        raise FileNotFoundError(
            f"No registry in {directory}. Build one with 'fasta-gen reindex --full'."
        )
    return snapshot


def superseded_snapshots(directory: Path, *, keep: int, suffix: str) -> list[Path]:
    """Return the snapshots beyond the newest ``keep``, oldest first.

    Deleting one is safe while an application holds it open: on POSIX the inode
    survives until that process closes it, so a running app keeps serving the
    snapshot it started with and the space is reclaimed when it restarts.

    ``suffix`` is not optional and not cosmetic. Counting across engines in a
    directory holding snapshots of both would treat the other engine's registries
    as superseded copies of this one's and delete live files.
    """
    if keep < 1:
        raise ValueError("keep must be at least 1")
    snapshots = list_snapshots(directory, suffix=suffix)
    return snapshots[:-keep] if len(snapshots) > keep else []
