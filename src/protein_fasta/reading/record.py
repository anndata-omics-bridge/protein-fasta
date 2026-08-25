"""Plain values emitted by protein-FASTA readers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """One FASTA record with an uninterpreted header and lexical sequence."""

    raw_header: str
    sequence: str


class FastaReadError(ValueError):
    """A source cannot be decoded or parsed safely as FASTA."""

    def __init__(
        self,
        source_name: str,
        reason: str,
        *,
        line_number: int | None = None,
    ) -> None:
        """Record the source and precise parse failure."""
        self.source_name = source_name
        self.reason = reason
        self.line_number = line_number
        location = f" at line {line_number}" if line_number is not None else ""
        super().__init__(f"{source_name}: {reason}{location}")
