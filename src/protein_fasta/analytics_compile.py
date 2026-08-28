"""Compile analytical JSON documents into schema-free runtime values."""

from __future__ import annotations

import re

from protein_fasta.analytics.digestion import Digestion
from protein_fasta.documents import load_builtin_enzyme_document
from protein_fasta.schema.analytics import DigestionDocument, EnzymeDocument


def compile_digestion(
    document: DigestionDocument,
    enzyme: EnzymeDocument,
    /,
) -> Digestion:
    """Compile one matching digestion/enzyme document pair."""
    if document.enzyme != enzyme.name:
        raise ValueError(
            f"digestion requests enzyme {document.enzyme!r}, but received {enzyme.name!r}"
        )
    try:
        cleavage = re.compile(enzyme.cleavage_pattern)
    except re.error as error:
        raise ValueError(f"invalid cleavage pattern for enzyme {enzyme.name!r}: {error}") from error
    return Digestion(
        enzyme=enzyme.name,
        cleavage=cleavage,
        min_length=document.min_length,
        max_length=document.max_length,
        missed_cleavages=document.missed_cleavages,
    )


def make_digestion(document: DigestionDocument, /) -> Digestion:
    """Resolve and compile one digestion against its packaged enzyme rule."""
    return compile_digestion(document, load_builtin_enzyme_document(document.enzyme))
