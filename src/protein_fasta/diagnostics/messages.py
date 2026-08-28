"""Human-readable messages for scalar FASTA diagnostic findings."""

from __future__ import annotations


def describe_illegal_residues(illegal: str, /) -> str:
    """Explain illegal residues without changing or re-diagnosing a sequence."""
    distinct = sorted(set(illegal))
    if distinct == ["*"]:
        return "stop codon inside the sequence, which means this is translated nucleotide output"
    rendered = " ".join(repr(character) for character in distinct)
    return f"illegal sequence characters {rendered}"
