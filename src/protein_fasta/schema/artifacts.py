"""Passive file-artifact evidence shared by workflow result documents."""

from pathlib import Path

from pydantic import Field

from protein_fasta.schema.base import DocumentBase


class ArtifactDocument(DocumentBase):
    """One versioned, checksummed workflow artifact."""

    schema_name: str = Field(min_length=1)
    schema_version: str = Field(min_length=1)
    path: Path
    checksum_version: str = Field(min_length=1)
    checksum: str = Field(min_length=1)
    byte_count: int = Field(ge=0)
    row_count: int | None = Field(default=None, ge=0)
