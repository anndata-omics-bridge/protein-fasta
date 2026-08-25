"""Protein-FASTA entry classes."""

from enum import StrEnum


class EntryKind(StrEnum):
    """Mutually exclusive primary kind of one protein-FASTA entry."""

    TARGET = "target"
    DECOY = "decoy"
    CONTAMINANT = "contaminant"
    ENTRAPMENT = "entrapment"
    SENTINEL = "sentinel"
