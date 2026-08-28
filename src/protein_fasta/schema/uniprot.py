"""Passive JSON documents for UniProt catalog and download workflows."""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Annotated, Literal

from pydantic import Field

from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class ReferenceProteomesDocument(DocumentBase):
    """Select the UniProt reference-proteome catalog."""

    type: Literal["reference"] = "reference"


class AllProteomesDocument(DocumentBase):
    """Select every UniProt proteome explicitly."""

    type: Literal["all"] = "all"


class ProteomeQueryDocument(DocumentBase):
    """Select proteomes using one explicit UniProt provider query."""

    type: Literal["query"] = "query"
    query: str = Field(min_length=1)


CatalogSelectionDocument = Annotated[
    ReferenceProteomesDocument | AllProteomesDocument | ProteomeQueryDocument,
    Field(discriminator="type"),
]


class UniProtCatalogRequestDocument(DocumentBase):
    """Request one immutable local snapshot of the UniProt proteome catalog."""

    schema_version: Literal["0.1"] = "0.1"
    output_dir: Path
    selection: CatalogSelectionDocument = Field(default_factory=ReferenceProteomesDocument)
    timeout_seconds: float = Field(default=120.0, gt=0.0)


class TaxonomySelectionDocument(DocumentBase):
    """Resolve a UniProt proteome from one positive NCBI taxonomy identifier."""

    type: Literal["taxid"] = "taxid"
    taxid: int = Field(gt=0)


class ProteomeIdSelectionDocument(DocumentBase):
    """Resolve one explicit UniProt proteome identifier."""

    type: Literal["proteome_id"] = "proteome_id"
    proteome_id: str = Field(min_length=1, pattern=r"^UP[0-9A-Z]+$")


ProteomeSelectionDocument = Annotated[
    TaxonomySelectionDocument | ProteomeIdSelectionDocument,
    Field(discriminator="type"),
]


class ReviewedDownloadDocument(DocumentBase):
    """Acquire reviewed UniProtKB records for one proteome."""

    type: Literal["swissprot"] = "swissprot"


class CompleteDownloadDocument(DocumentBase):
    """Acquire reviewed and unreviewed UniProtKB records for one proteome."""

    type: Literal["swissprot_trembl"] = "swissprot_trembl"


class CanonicalGeneDownloadDocument(DocumentBase):
    """Acquire one canonical UniProt protein sequence per gene."""

    type: Literal["one_seq_per_gene"] = "one_seq_per_gene"


UniProtAcquisitionDocument = Annotated[
    ReviewedDownloadDocument | CompleteDownloadDocument | CanonicalGeneDownloadDocument,
    Field(discriminator="type"),
]


class UniProtDownloadRequestDocument(DocumentBase):
    """Authored request for one UniProt proteome FASTA."""

    schema_version: Literal["0.1"] = "0.1"
    selection: ProteomeSelectionDocument
    acquisition: UniProtAcquisitionDocument = Field(default_factory=ReviewedDownloadDocument)
    output_fasta: Path
    timeout_seconds: float = Field(default=120.0, gt=0.0)


class EffectiveUniProtDownloadDocument(DocumentBase):
    """Resolved and directly replayable request for one UniProt FASTA."""

    schema_version: Literal["0.1"] = "0.1"
    selection: ProteomeSelectionDocument
    acquisition: UniProtAcquisitionDocument
    output_fasta: Path
    timeout_seconds: float = Field(gt=0.0)


class ResolvedProteomeDocument(DocumentBase):
    """Resolved provider identity and the evidence used to select it."""

    proteome_id: str = Field(min_length=1)
    taxid: int | None = Field(default=None, gt=0)
    protein_count: int | None = Field(default=None, ge=0)
    gene_count: int | None = Field(default=None, ge=0)
    organism: str | None = None
    resolution_method: Literal["explicit_id", "reference_taxid", "taxonomy_fallback"]
    resolution_query: str


class UniProtDownloadResultDocument(DocumentBase):
    """Machine-readable evidence for one completed UniProt FASTA acquisition."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    effective_request: EffectiveUniProtDownloadDocument
    resolved_proteome: ResolvedProteomeDocument
    provider_query: str
    observed_releases: tuple[str, ...]
    actual_entry_count: int = Field(gt=0)
    provider_reported_counts: tuple[int, ...] = ()
    artifact: ArtifactDocument
    warnings: tuple[str, ...] = ()


class UniProtCatalogResultDocument(DocumentBase):
    """Machine-readable evidence for one completed proteome-catalog snapshot."""

    schema_version: Literal["0.1"] = "0.1"
    protein_fasta_version: str
    request: UniProtCatalogRequestDocument
    provider_query: str
    retrieved_at: datetime.datetime
    observed_releases: tuple[str, ...]
    provider_reported_counts: tuple[int, ...] = ()
    artifact: ArtifactDocument
    warnings: tuple[str, ...] = ()
