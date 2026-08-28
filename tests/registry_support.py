"""Shared registry configuration and fixture helpers."""

from __future__ import annotations

import datetime
from pathlib import Path

from pydantic import BaseModel, Field

from protein_fasta.schema.build import MetadataDocument, NamingDocument


class BackendSettings(BaseModel):
    """Backend selection used by registry tests."""

    backend: str = "sqlite"


class Settings(BaseModel):
    """Exact configuration capability consumed by registry operations."""

    fasta_root: Path
    registry_dir: Path
    registry: BackendSettings = Field(default_factory=BackendSettings)
    max_fasta_file_size_gib: float = 5.0
    max_detailed_entries: int = 200_000
    min_fasta_date: datetime.date | None = None
    metadata_aa_sample_size: int = 20_000
    overlap_threshold: float = 0.99
    registry_diagnostics_path: Path | None = None
    naming: NamingDocument = Field(default_factory=NamingDocument)
    sentinel: MetadataDocument = Field(default_factory=MetadataDocument)
