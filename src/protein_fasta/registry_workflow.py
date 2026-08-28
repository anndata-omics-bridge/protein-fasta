"""Compile portable registry policy for concrete registry operations."""

from __future__ import annotations

import datetime
from dataclasses import dataclass
from pathlib import Path

from protein_fasta.inventory import read_database_inventory
from protein_fasta.registry.backend import factory as registry_backends
from protein_fasta.registry.indexing import (
    RegistryRecord,
    connect_registry,
    index_inventory_entries,
    initialize_registry,
)
from protein_fasta.schema.build import MetadataDocument, NamingDocument
from protein_fasta.schema.registry import RegistryBackendDocument, RegistryDocument


@dataclass(frozen=True, slots=True)
class RegistryOperationSettings:
    """Exact runtime settings exercised by registry indexing and review."""

    fasta_root: Path
    registry_dir: Path
    registry: RegistryBackendDocument
    max_fasta_file_size_gib: float
    max_detailed_entries: int
    metadata_aa_sample_size: int
    min_fasta_date: datetime.date | None
    overlap_threshold: float
    naming: NamingDocument
    sentinel: MetadataDocument
    registry_diagnostics_path: Path | None = None


def make_registry_settings(
    document: RegistryDocument,
    /,
    *,
    fasta_root: Path,
    registry_path: Path,
) -> RegistryOperationSettings:
    """Compile one passive registry document for a concrete path pair."""
    return RegistryOperationSettings(
        fasta_root=fasta_root,
        registry_dir=registry_path.parent,
        registry=document.registry,
        max_fasta_file_size_gib=document.max_fasta_file_size_gib,
        max_detailed_entries=document.max_detailed_entries,
        metadata_aa_sample_size=document.metadata_aa_sample_size,
        min_fasta_date=document.min_fasta_date,
        overlap_threshold=document.overlap_threshold,
        naming=document.naming,
        sentinel=document.metadata,
    )


def index_database_inventory(
    inventory_path: Path,
    registry_path: Path,
    document: RegistryDocument,
    /,
    *,
    label: str | None = None,
) -> RegistryRecord:
    """Index a canonical biological or search inventory through either backend."""
    path_backend = registry_backends.backend_for_path(registry_path)
    if path_backend != document.registry.backend:
        raise ValueError(
            f"registry path selects {path_backend!r}, but configuration selects "
            f"{document.registry.backend!r}"
        )
    resolved_inventory = inventory_path.resolve()
    resolved_registry = registry_path.resolve()
    settings = make_registry_settings(
        document,
        fasta_root=resolved_inventory.parent,
        registry_path=resolved_registry,
    )
    entries = read_database_inventory(resolved_inventory)
    stat = resolved_inventory.stat()
    database_label = label or _database_label(resolved_inventory)
    resolved_registry.parent.mkdir(parents=True, exist_ok=True)
    with connect_registry(resolved_registry) as connection:
        initialize_registry(connection, settings)
        return index_inventory_entries(
            connection,
            entries,
            settings,
            relative_path=database_label,
            filename=Path(database_label).name,
            artifact_size_bytes=stat.st_size,
            artifact_mtime_ns=stat.st_mtime_ns,
        )


def _database_label(inventory_path: Path) -> str:
    """Recover the database filename carried by a canonical sidecar name."""
    name = inventory_path.name
    for suffix in (".protein-inventory.parquet", ".search-inventory.parquet"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name
