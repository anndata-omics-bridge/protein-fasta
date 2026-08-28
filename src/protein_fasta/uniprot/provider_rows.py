"""Project UniProt provider objects into exact runtime records."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from protein_fasta.reading.parser import FastaRecord
from protein_fasta.uniprot.models import ResolutionMethod, ResolvedProteome, UniProtCatalogRow


def resolved_proteome(
    raw: Mapping[str, object],
    /,
    *,
    method: ResolutionMethod,
    query: str,
) -> ResolvedProteome:
    """Project one provider proteome object into a resolved runtime value."""
    proteome_id = _required_string(raw.get("id"), "proteome id")
    taxonomy = _mapping(raw.get("taxonomy"))
    return ResolvedProteome(
        proteome_id=proteome_id,
        taxid=_integer(taxonomy.get("taxonId")),
        protein_count=_integer(raw.get("proteinCount")),
        gene_count=_integer(raw.get("geneCount")),
        organism=_string(taxonomy.get("scientificName")),
        resolution_method=method,
        resolution_query=query,
    )


def catalog_row(raw: Mapping[str, object], /) -> UniProtCatalogRow:
    """Project one provider proteome object into the canonical catalog row."""
    stats = _mapping(raw.get("proteomeStatistics"))
    reviewed = _integer(stats.get("reviewedProteinCount"))
    unreviewed = _integer(stats.get("unreviewedProteinCount"))
    total = _integer(raw.get("proteinCount"))
    if total is None and (reviewed is not None or unreviewed is not None):
        total = (reviewed or 0) + (unreviewed or 0)
    taxonomy = _mapping(raw.get("taxonomy"))
    return UniProtCatalogRow(
        proteome_id=_required_string(raw.get("id"), "proteome id"),
        taxid=_integer(taxonomy.get("taxonId")),
        organism=_string(taxonomy.get("scientificName")),
        proteome_type=_string(raw.get("proteomeType")),
        swissprot=reviewed,
        swissprot_trembl=total,
        one_seq_per_gene=_integer(raw.get("geneCount")),
    )


def canonical_fasta_record(raw: Mapping[str, object], /) -> FastaRecord | None:
    """Project one GeneCentric canonical protein into a UniProt FASTA record."""
    canonical = _mapping(raw.get("canonicalProtein"))
    if not canonical:
        return None
    sequence = _string(_mapping(canonical.get("sequence")).get("value"))
    if not sequence:
        return None
    return FastaRecord(raw_header=canonical_fasta_header(canonical), sequence=sequence)


def canonical_fasta_header(protein: Mapping[str, object], /) -> str:
    """Reconstruct a standard UniProt FASTA header for one canonical protein."""
    entry_type = (_string(protein.get("entryType")) or "").lower()
    namespace = "sp" if "reviewed" in entry_type and "unreviewed" not in entry_type else "tr"
    accession = _string(protein.get("id")) or ""
    entry_name = _string(protein.get("uniProtkbId")) or accession
    protein_name = _string(protein.get("proteinName")) or ""
    parts = [f"{namespace}|{accession}|{entry_name} {protein_name}".rstrip()]

    organism = _mapping(protein.get("organism"))
    scientific_name = _string(organism.get("scientificName"))
    taxid = _integer(organism.get("taxonId"))
    if scientific_name:
        parts.append(f"OS={scientific_name}")
    if taxid is not None:
        parts.append(f"OX={taxid}")
    gene_name = _string(protein.get("geneName"))
    if gene_name:
        parts.append(f"GN={gene_name}")
    existence = _string(protein.get("proteinExistence")) or ""
    if existence[:1].isdigit():
        parts.append(f"PE={existence[0]}")
    sequence_version = _integer(protein.get("sequenceVersion"))
    if sequence_version is not None:
        parts.append(f"SV={sequence_version}")
    return " ".join(parts)


def _mapping(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        return {}
    mapping = cast("Mapping[object, object]", value)
    return {str(key): item for key, item in mapping.items()}


def _string(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _required_string(value: object, field: str) -> str:
    result = _string(value)
    if not result:
        raise ValueError(f"UniProt response is missing {field}")
    return result


def _integer(value: object) -> int | None:
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value
