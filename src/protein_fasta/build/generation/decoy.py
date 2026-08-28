"""Generate decoy FASTA entries.

Default is plain full-sequence reversal (matching prozor and the existing FGCZ
collection), with the decoy prefix prepended to the whole header — so a decoy of
``sp|P123|X ...`` becomes ``REV_sp|P123|X ...`` with the reversed sequence.
Selectable shuffle and DecoyPYrat modes delegate to the standalone
``fdr_benchmark`` algorithm package.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from protein_fasta.build.generation.decoy_types import (
    DecoyBatch,
    DecoyGeneration,
    Entry,
)
from protein_fasta.schema.build import DecoyDocument, DecoyMode

DEFAULT_DECOY_SPEC = DecoyDocument()

# The shape decoy_annotation writes, read back to inherit a source's settings.
_DECOY_ANNOTATION = re.compile(r"decoys (?P<mode>[a-z]+) seed (?P<seed>\d+)")


@dataclass(frozen=True, slots=True)
class ReverseDecoyGeneration:
    """Generate deterministic whole-sequence reversals."""

    @property
    def mode(self) -> DecoyMode:
        """Return the reverse algorithm label."""
        return DecoyMode.REVERSE

    @property
    def seed(self) -> None:
        """Return no seed because reversal is deterministic."""
        return None

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Reverse each source sequence and prefix its complete header."""
        generated = tuple(
            make_decoy(description, sequence, prefix=prefix) for description, sequence in entries
        )
        return DecoyBatch(generated, self.parameters())

    def parameters(self) -> dict[str, Any]:
        """Return reverse-generation provenance."""
        return {"mode": self.mode.value}

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Return no note because reversal is the documented default."""
        del initial_collisions, dropped_peptides
        return ""


def make_decoy_generation(spec: DecoyDocument, /) -> DecoyGeneration:
    """Compile one storage document into algorithm-owning runtime behavior."""
    if spec.mode == DecoyMode.REVERSE:
        return ReverseDecoyGeneration()
    return _make_advanced_decoy_generation(spec)


def _make_advanced_decoy_generation(spec: DecoyDocument) -> DecoyGeneration:
    """Load the optional advanced-generation adapter only when selected."""
    try:
        from protein_fasta.build.generation.fdr_decoy import make_fdr_decoy_generation
    except ModuleNotFoundError as error:
        missing_name = error.name or ""
        if missing_name == "fdr_benchmark" or missing_name.startswith("fdr_benchmark."):
            message = (
                f"{spec.mode.value} decoy generation requires the 'protein-fasta[generation]' extra"
            )
            raise RuntimeError(message) from error
        raise
    return make_fdr_decoy_generation(spec)


def decoy_sequence(sequence: str, mode: str = "reverse") -> str:
    """Return the compatibility whole-sequence reverse transform."""
    if mode != DecoyMode.REVERSE:
        message = "Use generate_decoys() for seeded shuffle or DecoyPYrat generation."
        raise ValueError(message)
    return sequence[::-1]


def make_decoy(
    description: str, sequence: str, *, prefix: str, mode: str = "reverse"
) -> tuple[str, str]:
    """Build a decoy ``(description, sequence)`` from a target entry."""
    return prefix + description, decoy_sequence(sequence, mode)


def generate_decoys(
    entries: tuple[Entry, ...],
    *,
    prefix: str,
    spec: DecoyDocument = DEFAULT_DECOY_SPEC,
) -> DecoyBatch:
    """Generate one decoy per unique source entry under the selected algorithm."""
    return make_decoy_generation(spec).generate(entries, prefix=prefix)


def decoy_parameters(spec: DecoyDocument) -> dict[str, Any]:
    """Return machine-readable parameters for sentinel-sidecar provenance."""
    return make_decoy_generation(spec).parameters()


def parse_decoy_annotation(annotation: str | None) -> DecoyDocument | None:
    """Recover the decoy settings a build recorded in its sentinel annotation.

    The pipeline writes this note only for a non-default mode, so a database with
    decoys and no note was built with the default reverse. ``None`` means the
    annotation says nothing about decoys at all, which a caller must not read as
    "reverse": an unannotated file may predate the convention.
    """
    match = _DECOY_ANNOTATION.search(annotation or "")
    if match is None:
        return None
    try:
        mode = DecoyMode(match.group("mode"))
    except ValueError:
        # A mode this build does not know is not a reason to guess a different one.
        return None
    return DecoyDocument(mode=mode, seed=int(match.group("seed")))


def decoy_annotation(
    spec: DecoyDocument,
    *,
    initial_collisions: int | None = None,
    dropped_peptides: int | None = None,
) -> str:
    """Render compact non-default decoy provenance for the FASTA sentinel."""
    return make_decoy_generation(spec).annotation(
        initial_collisions=initial_collisions,
        dropped_peptides=dropped_peptides,
    )
