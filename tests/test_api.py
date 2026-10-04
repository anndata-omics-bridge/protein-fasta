"""The public module other anndata_bridge packages import."""

from __future__ import annotations

from protein_fasta import api


def test_api_exports_exactly_the_approved_names() -> None:
    assert sorted(api.__all__) == [
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
