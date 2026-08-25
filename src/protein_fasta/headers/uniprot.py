"""UniProt protein-FASTA header interpretation."""

from __future__ import annotations

import re

from protein_fasta.headers.generic import ProteinHeader, parse_header

_GENE_NAME = re.compile(r"(?:^|\s)GN=(?P<gene>\S+)")
_PREFIXED_ACCESSION = re.compile(r"^(?P<prefix>.*?)(?:sp|tr)\|(?P<accession>[^|]+)\|")
_MIDDLE_FIELD = re.compile(r".+\|(?P<accession>[^|]+)\|.*")


def uniprot_accession(identifier: str) -> str:
    """Extract a UniProt accession while retaining a target/decoy prefix."""
    match = _PREFIXED_ACCESSION.match(identifier)
    if match is not None:
        return f"{match.group('prefix')}{match.group('accession')}"
    match = _MIDDLE_FIELD.match(identifier)
    return match.group("accession") if match is not None else identifier


def uniprot_gene_name(raw_header: str) -> str | None:
    """Return a UniProt `GN=` value when the header declares one."""
    match = _GENE_NAME.search(raw_header)
    return match.group("gene") if match is not None else None


class UniProtHeaderInterpreter:
    """Interpret UniProt identifiers, accessions, and optional gene names."""

    def interpret(self, raw_header: str, /) -> ProteinHeader:
        """Return UniProt-aware values while retaining the original identifier."""
        parsed = parse_header(raw_header)
        return ProteinHeader(
            identifier=parsed.identifier,
            description=parsed.description,
            protein_name=uniprot_accession(parsed.identifier),
            gene_name=uniprot_gene_name(raw_header),
        )
