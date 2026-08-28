"""Atomic persistence and evidence for workflow artifacts."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

from protein_fasta.analytics.hashing import FILE_CHECKSUM_VERSION, file_checksum
from protein_fasta.schema.artifacts import ArtifactDocument


@contextmanager
def temporary_sibling(destination: Path, /) -> Generator[Path]:
    """Yield a same-directory temporary path and remove it on exit."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    os.close(descriptor)
    path = Path(name)
    try:
        yield path
    finally:
        path.unlink(missing_ok=True)


def publish_exclusive(staged: Path, destination: Path, /) -> None:
    """Publish a staged file atomically without replacing an existing artifact."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.link(staged, destination)
    staged.unlink()


def write_json_atomic(
    path: Path,
    payload: object,
    /,
    *,
    replace_existing: bool,
) -> None:
    """Write deterministic JSON through a same-directory atomic publication."""
    with temporary_sibling(path) as staged:
        staged.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if replace_existing:
            os.replace(staged, path)
        else:
            publish_exclusive(staged, path)


def artifact_document(
    path: Path,
    /,
    *,
    recorded_path: Path,
    schema_name: str,
    schema_version: str,
    row_count: int | None = None,
) -> ArtifactDocument:
    """Describe exact bytes at one path using a portable recorded path."""
    return ArtifactDocument(
        schema_name=schema_name,
        schema_version=schema_version,
        path=recorded_path,
        checksum_version=FILE_CHECKSUM_VERSION,
        checksum=file_checksum(path),
        byte_count=path.stat().st_size,
        row_count=row_count,
    )
