"""Operational record kinds stored by a protein-FASTA registry."""

from enum import StrEnum


class EntryKind(StrEnum):
    """Mutually exclusive record kinds required by registry operations."""

    TARGET = "target"
    DECOY = "decoy"
    CONTAMINANT = "contaminant"
    ENTRAPMENT = "entrapment"
    SENTINEL = "sentinel"


class DetailLevel(StrEnum):
    """Availability of record-level details for one registered FASTA."""

    FULL = "full"
    METADATA_ONLY = "metadata_only"
