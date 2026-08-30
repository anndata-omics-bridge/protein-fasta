"""Compile passive decoy documents into runtime generation behavior."""

from __future__ import annotations

import re

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.decoy import DecoyGeneration, DecoyMode, ReverseDecoyGeneration
from protein_fasta.database.decoy_generation import (
    make_decoypyrat_generation,
    make_shuffle_decoy_generation,
)
from protein_fasta.schema.decoy import (
    DecoyPyratDocument,
    DecoyStrategyDocument,
    ReverseDecoyDocument,
    ShuffleDecoyDocument,
)

_DECOY_ANNOTATION = re.compile(r"decoys (?P<mode>[a-z]+) seed (?P<seed>\d+)")


def make_decoy_generation(document: DecoyStrategyDocument, /) -> DecoyGeneration:
    """Compile one decoy strategy document at the root composition boundary."""
    if isinstance(document, ReverseDecoyDocument):
        return ReverseDecoyGeneration()
    if isinstance(document, ShuffleDecoyDocument):
        return make_shuffle_decoy_generation(document.seed)
    digestion = make_digestion(document.digestion)
    return make_decoypyrat_generation(
        seed=document.seed,
        digestion=digestion,
    )


def parse_decoy_annotation(annotation: str | None, /) -> DecoyStrategyDocument | None:
    """Recover explicit strategy evidence from a historical sentinel annotation."""
    match = _DECOY_ANNOTATION.search(annotation or "")
    if match is None:
        return None
    try:
        mode = DecoyMode(match.group("mode"))
    except ValueError:
        return None
    seed = int(match.group("seed"))
    if mode is DecoyMode.SHUFFLE:
        return ShuffleDecoyDocument(seed=seed)
    if mode is DecoyMode.DECOYPYRAT:
        return DecoyPyratDocument(seed=seed)
    return ReverseDecoyDocument()
