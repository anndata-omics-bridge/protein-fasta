"""Behavior-owning UniProt proteome resolution variants."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from protein_fasta.uniprot.models import ResolutionMethod, ResolvedProteome
from protein_fasta.uniprot.provider_rows import resolved_proteome
from protein_fasta.uniprot.queries import reference_taxonomy_query, taxonomy_query
from protein_fasta.uniprot.transport import UniProtTransport


class ProteomeResolution(Protocol):
    """Resolve exactly one UniProt proteome through a provider transport."""

    def resolve(self, transport: UniProtTransport, /) -> ResolvedProteome:
        """Resolve and return one proteome with selection evidence."""
        ...


@dataclass(frozen=True, slots=True)
class ResolveProteomeId:
    """Resolve one explicit UniProt proteome identifier."""

    proteome_id: str

    def resolve(self, transport: UniProtTransport, /) -> ResolvedProteome:
        query = f"proteome_id:{self.proteome_id}"
        return resolved_proteome(
            transport.proteome(self.proteome_id),
            method=ResolutionMethod.EXPLICIT_ID,
            query=query,
        )


@dataclass(frozen=True, slots=True)
class ResolveTaxonomy:
    """Resolve one taxonomy, preferring its reference proteome."""

    taxid: int

    def resolve(self, transport: UniProtTransport, /) -> ResolvedProteome:
        reference_query = reference_taxonomy_query(self.taxid)
        raw = transport.first_proteome(reference_query)
        if raw is not None:
            return resolved_proteome(
                raw,
                method=ResolutionMethod.REFERENCE_TAXID,
                query=reference_query,
            )
        fallback_query = taxonomy_query(self.taxid)
        raw = transport.first_proteome(fallback_query)
        if raw is None:
            raise LookupError(f"no UniProt proteome found for taxid {self.taxid}")
        return resolved_proteome(
            raw,
            method=ResolutionMethod.TAXONOMY_FALLBACK,
            query=fallback_query,
        )
