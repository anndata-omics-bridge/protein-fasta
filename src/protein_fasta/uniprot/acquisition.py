"""Behavior-owning UniProt protein acquisition variants."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from protein_fasta.reading.parser import FastaRecord
from protein_fasta.reading.writer import write_records
from protein_fasta.uniprot.models import (
    ProviderTransferEvidence,
    ResolvedProteome,
    UniProtDownloadMode,
    UniProtDownloadOutcome,
)
from protein_fasta.uniprot.provider_rows import canonical_fasta_record
from protein_fasta.uniprot.queries import (
    canonical_gene_query,
    complete_proteome_query,
    reviewed_proteome_query,
)
from protein_fasta.uniprot.transport import UniProtTransport


class ProteinAcquisition(Protocol):
    """Acquire one selected UniProt protein record set."""

    @property
    def mode(self) -> UniProtDownloadMode:
        """Return the stable serialized acquisition identity."""
        ...

    def acquire(
        self,
        transport: UniProtTransport,
        proteome: ResolvedProteome,
        destination: Path,
        /,
    ) -> UniProtDownloadOutcome:
        """Write unpublished records and return provider evidence."""
        ...


@dataclass(frozen=True, slots=True)
class ReviewedProteins:
    """Acquire reviewed UniProtKB proteins."""

    @property
    def mode(self) -> UniProtDownloadMode:
        return UniProtDownloadMode.SWISSPROT

    def acquire(
        self,
        transport: UniProtTransport,
        proteome: ResolvedProteome,
        destination: Path,
        /,
    ) -> UniProtDownloadOutcome:
        query = reviewed_proteome_query(proteome.proteome_id)
        transfer = transport.stream_uniprotkb_fasta(query, destination)
        return UniProtDownloadOutcome(self.mode, proteome, query, transfer)


@dataclass(frozen=True, slots=True)
class CompleteProteins:
    """Acquire reviewed and unreviewed UniProtKB proteins."""

    @property
    def mode(self) -> UniProtDownloadMode:
        return UniProtDownloadMode.SWISSPROT_TREMBL

    def acquire(
        self,
        transport: UniProtTransport,
        proteome: ResolvedProteome,
        destination: Path,
        /,
    ) -> UniProtDownloadOutcome:
        query = complete_proteome_query(proteome.proteome_id)
        transfer = transport.stream_uniprotkb_fasta(query, destination)
        return UniProtDownloadOutcome(self.mode, proteome, query, transfer)


@dataclass(frozen=True, slots=True)
class CanonicalGeneProteins:
    """Acquire one canonical protein per GeneCentric record."""

    @property
    def mode(self) -> UniProtDownloadMode:
        return UniProtDownloadMode.ONE_SEQ_PER_GENE

    def acquire(
        self,
        transport: UniProtTransport,
        proteome: ResolvedProteome,
        destination: Path,
        /,
    ) -> UniProtDownloadOutcome:
        query = canonical_gene_query(proteome.proteome_id)
        pages = transport.iter_genecentric_pages(query)
        release_values: list[str] = []
        reported_counts: list[int] = []
        actual_count = 0

        def records() -> Iterator[FastaRecord]:
            nonlocal actual_count
            for page in pages:
                if page.release is not None:
                    release_values.append(page.release)
                if page.reported_count is not None:
                    reported_counts.append(page.reported_count)
                for raw in page.records:
                    record = canonical_fasta_record(raw)
                    if record is not None:
                        actual_count += 1
                        yield record

        write_records(records(), destination)
        transfer = ProviderTransferEvidence(
            actual_entry_count=actual_count,
            observed_releases=tuple(dict.fromkeys(release_values)),
            provider_reported_counts=tuple(dict.fromkeys(reported_counts)),
        )
        return UniProtDownloadOutcome(self.mode, proteome, query, transfer)
