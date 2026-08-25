"""Generic protein-FASTA header interpretation."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True, slots=True)
class ParsedHeader:
    """Identifier token and normalized optional description."""

    identifier: str
    description: str | None


@dataclass(frozen=True, slots=True)
class ProteinHeader:
    """Application-neutral protein identity derived from one raw header."""

    identifier: str
    description: str | None
    protein_name: str
    gene_name: str | None


class HeaderInterpreter(Protocol):
    """Interpret one raw FASTA header according to a configured convention."""

    def interpret(self, raw_header: str, /) -> ProteinHeader:
        """Return semantic protein header values."""
        ...


def parse_header(raw_header: str) -> ParsedHeader:
    """Split a raw header into its first token and normalized description."""
    text = raw_header.removeprefix(">")
    if not text or text[0].isspace():
        return ParsedHeader(identifier="", description=None)
    fields = text.split(maxsplit=1)
    description = " ".join(fields[1].split()) if len(fields) == 2 else None
    return ParsedHeader(identifier=fields[0], description=description or None)


def normalized_description_hash(raw_header: str) -> bytes | None:
    """Return a stable 128-bit hash of the semantic header description."""
    description = parse_header(raw_header).description
    if description is None:
        return None
    return hashlib.blake2b(description.encode(), digest_size=16).digest()


class GenericHeaderInterpreter:
    """Use the complete identifier token as the protein name."""

    def interpret(self, raw_header: str, /) -> ProteinHeader:
        """Interpret a header without vendor-specific accession rules."""
        parsed = parse_header(raw_header)
        return ProteinHeader(
            identifier=parsed.identifier,
            description=parsed.description,
            protein_name=parsed.identifier,
            gene_name=None,
        )
