"""High-level normalized protein records and diagnostic iteration."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import read_records
from protein_fasta.validation.sequence import normalize_sequence


@dataclass(frozen=True, slots=True)
class ProteinRecord:
    """One generically parsed and normalized protein FASTA entry."""

    id: str
    description: str | None
    sequence: str


@dataclass(frozen=True, slots=True)
class ProteinDiagnostics:
    """A normalized protein composed with configured audit facts."""

    protein: ProteinRecord
    raw_header: str
    identifier_namespace: str
    classifications: frozenset[str]
    upper_cased: bool
    stop_stripped: bool
    illegal_residues: str


def iter_proteins(path: Path, /) -> Iterator[ProteinRecord]:
    """Stream normalized protein records from a FASTA path."""
    for lexical in read_records(path):
        parsed = parse_header(lexical.raw_header)
        normalized = normalize_sequence(lexical.sequence)
        yield ProteinRecord(parsed.id, parsed.description, normalized.sequence)


def iter_protein_diagnostics(
    path: Path,
    rules: DiagnosticRules,
    /,
) -> Iterator[ProteinDiagnostics]:
    """Stream normalized records with configured diagnostic facts."""
    for lexical in read_records(path):
        parsed = parse_header(lexical.raw_header)
        normalized = normalize_sequence(lexical.sequence)
        namespace, classifications = rules.diagnose_identifier(parsed.id)
        protein = ProteinRecord(parsed.id, parsed.description, normalized.sequence)
        yield ProteinDiagnostics(
            protein=protein,
            raw_header=lexical.raw_header,
            identifier_namespace=namespace,
            classifications=classifications,
            upper_cased=normalized.upper_cased,
            stop_stripped=normalized.stop_stripped,
            illegal_residues=rules.illegal_residues(normalized.sequence),
        )
