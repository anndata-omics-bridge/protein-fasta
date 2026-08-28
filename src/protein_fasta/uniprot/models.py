"""Runtime values for UniProt source acquisition."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class UniProtDownloadMode(StrEnum):
    """Stable serialized identities for the three supported acquisitions."""

    SWISSPROT = "swissprot"
    SWISSPROT_TREMBL = "swissprot_trembl"
    ONE_SEQ_PER_GENE = "one_seq_per_gene"


class ResolutionMethod(StrEnum):
    """Evidence describing how one proteome was selected."""

    EXPLICIT_ID = "explicit_id"
    REFERENCE_TAXID = "reference_taxid"
    TAXONOMY_FALLBACK = "taxonomy_fallback"


@dataclass(frozen=True, slots=True)
class ResolvedProteome:
    """One resolved UniProt proteome and its selection evidence."""

    proteome_id: str
    taxid: int | None
    protein_count: int | None
    gene_count: int | None
    organism: str | None
    resolution_method: ResolutionMethod
    resolution_query: str


@dataclass(frozen=True, slots=True)
class ProviderPage:
    """One page of provider objects and response-level evidence."""

    records: tuple[dict[str, object], ...]
    release: str | None
    reported_count: int | None


@dataclass(frozen=True, slots=True)
class ProviderTransferEvidence:
    """Counts and headers observed while acquiring provider records."""

    actual_entry_count: int
    observed_releases: tuple[str, ...]
    provider_reported_counts: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class UniProtDownloadOutcome:
    """Runtime outcome before file-artifact persistence."""

    mode: UniProtDownloadMode
    resolved_proteome: ResolvedProteome
    provider_query: str
    transfer: ProviderTransferEvidence


@dataclass(frozen=True, slots=True)
class UniProtCatalogRow:
    """One canonical row in a local UniProt proteome catalog."""

    proteome_id: str
    taxid: int | None
    organism: str | None
    proteome_type: str | None
    swissprot: int | None
    swissprot_trembl: int | None
    one_seq_per_gene: int | None
