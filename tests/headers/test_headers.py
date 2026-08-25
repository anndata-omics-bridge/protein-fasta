"""Tests for generic and UniProt header interpretation."""

from protein_fasta.compile import make_header_interpreter
from protein_fasta.headers.generic import (
    ParsedHeader,
    ProteinHeader,
    normalized_description_hash,
    parse_header,
)
from protein_fasta.headers.uniprot import uniprot_accession, uniprot_gene_name
from protein_fasta.schema.headers import GenericHeaderDocument, UniProtHeaderDocument


def test_parse_header_normalizes_description_without_losing_identifier() -> None:
    assert parse_header(">sp|P1|ONE   Protein kinase,   alpha ") == ParsedHeader(
        identifier="sp|P1|ONE",
        description="Protein kinase, alpha",
    )
    assert parse_header("P1") == ParsedHeader("P1", None)
    assert parse_header(" description only") == ParsedHeader("", None)


def test_generic_interpreter_uses_identifier_as_protein_name() -> None:
    interpreter = make_header_interpreter(GenericHeaderDocument())
    assert interpreter.interpret("custom|P1 description") == ProteinHeader(
        identifier="custom|P1",
        description="description",
        protein_name="custom|P1",
        gene_name=None,
    )


def test_uniprot_interpreter_preserves_decoy_prefix_and_extracts_gene() -> None:
    interpreter = make_header_interpreter(UniProtHeaderDocument())
    assert interpreter.interpret(
        "REV_sp|P12345|KINASE Protein kinase OS=Human GN=MAPK1 PE=1"
    ) == ProteinHeader(
        identifier="REV_sp|P12345|KINASE",
        description="Protein kinase OS=Human GN=MAPK1 PE=1",
        protein_name="REV_P12345",
        gene_name="MAPK1",
    )


def test_uniprot_helpers_have_generic_fallbacks() -> None:
    assert uniprot_accession("db|ACC|description") == "ACC"
    assert uniprot_accession("plain") == "plain"
    assert uniprot_gene_name("sp|P1|ONE no gene") is None


def test_description_hash_uses_normalized_semantic_text() -> None:
    assert normalized_description_hash(">P1  Protein   kinase ") == normalized_description_hash(
        "P1 Protein kinase"
    )
    assert normalized_description_hash("P1") is None
