"""Tests for the two high-level Python record products."""

from pathlib import Path

from protein_fasta.compile import make_diagnostic_rules
from protein_fasta.documents import (
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
)
from protein_fasta.record import (
    ProteinDiagnostics,
    ProteinRecord,
    iter_protein_diagnostics,
    iter_proteins,
)


def test_iter_proteins_returns_only_base_normalized_meanings(tmp_path: Path) -> None:
    path = tmp_path / "proteins.fasta"
    path.write_text(">REV_sp|P12345|RL40_YEAST   Protein   name OS=Yeast\nac d*\n")

    assert list(iter_proteins(path)) == [
        ProteinRecord(
            id="REV_sp|P12345|RL40_YEAST",
            description="Protein name OS=Yeast",
            sequence="ACD",
        )
    ]


def test_iter_diagnostics_composes_record_and_audit_once(tmp_path: Path) -> None:
    path = tmp_path / "proteins.fasta"
    path.write_text(">REV_CON__sp|P1|ENTRY_ORG description\nacd-*\n")
    rules = make_diagnostic_rules(
        load_builtin_diagnostic_document(),
        load_builtin_entry_classifier_document(),
    )

    assert list(iter_protein_diagnostics(path, rules)) == [
        ProteinDiagnostics(
            protein=ProteinRecord(
                id="REV_CON__sp|P1|ENTRY_ORG",
                description="description",
                sequence="ACD-",
            ),
            raw_header="REV_CON__sp|P1|ENTRY_ORG description",
            identifier_namespace="uniprot",
            classifications=frozenset({"decoy", "contaminant"}),
            upper_cased=True,
            stop_stripped=True,
            illegal_residues="-",
        )
    ]
