"""The one public protein_fasta module for other anndata_bridge packages."""

from __future__ import annotations

from protein_fasta.analytics.digestion import digest_sequence
from protein_fasta.analytics.hashing import file_checksum, peptide_hash, sequence_hash
from protein_fasta.analytics.peptide_properties import PeptideProperties
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.frame import (
    ProteinDatabase,
    ProteinFormat,
    read_basic_protein_frame,
    read_configured_protein_frame,
    read_protein_frame,
    refseq,
    uniprotkb,
)
from protein_fasta.peptide_frame import peptide_property_frame
from protein_fasta.record import iter_protein_diagnostics, iter_proteins

__all__ = [
    "PeptideProperties",
    "ProteinDatabase",
    "ProteinFormat",
    "digest_sequence",
    "file_checksum",
    "iter_protein_diagnostics",
    "iter_proteins",
    "make_digestion",
    "peptide_hash",
    "peptide_property_frame",
    "read_basic_protein_frame",
    "read_configured_protein_frame",
    "read_protein_frame",
    "refseq",
    "sequence_hash",
    "uniprotkb",
]
