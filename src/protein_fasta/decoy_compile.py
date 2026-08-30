"""Compile passive decoy documents into runtime generation behavior."""

from __future__ import annotations

import re

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.decoy import DecoyGeneration, DecoyMode, ReverseDecoyGeneration
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
    try:
        from protein_fasta.database.decoy_advanced import (
            make_decoypyrat_generation,
            make_shuffle_decoy_generation,
        )
    except ModuleNotFoundError as error:
        missing_name = error.name or ""
        if missing_name == "fdr_benchmark" or missing_name.startswith("fdr_benchmark."):
            raise RuntimeError(
                f"{document.type} decoy generation requires the 'protein-fasta[generation]' extra"
            ) from error
        raise
    if isinstance(document, ShuffleDecoyDocument):
        return make_shuffle_decoy_generation(document.seed)
    digestion = make_digestion(document.digestion)
    return make_decoypyrat_generation(
        seed=document.seed,
        enzyme=digestion.cleavage.pattern,
        minimum_length=document.digestion.min_length,
        maximum_length=document.digestion.max_length,
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
