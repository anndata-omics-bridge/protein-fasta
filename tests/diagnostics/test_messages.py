"""Agreement tests for shared normalization and configured FGCZ diagnostics."""

from __future__ import annotations

from pathlib import Path

import pytest

from protein_fasta.diagnostics.messages import describe_illegal_residues
from protein_fasta.diagnostics.runtime import UNMATCHED_NAMESPACE
from protein_fasta.registry.rules import RegistryDiagnosticRules, load_registry_diagnostics
from protein_fasta.validation.sequence import normalize_sequence


def _diagnostics() -> RegistryDiagnosticRules:
    return load_registry_diagnostics()


def test_normalization_upper_cases_and_removes_one_trailing_stop() -> None:
    plain = normalize_sequence("MPEPTIDEK")
    lowered = normalize_sequence("mpeptidek")
    terminated = normalize_sequence("MPEPTIDEK*")
    doubly_terminated = normalize_sequence("MPEPTIDEK**")

    assert (plain.sequence, plain.upper_cased, plain.stop_stripped) == ("MPEPTIDEK", False, False)
    assert (lowered.sequence, lowered.upper_cased, lowered.stop_stripped) == (
        "MPEPTIDEK",
        True,
        False,
    )
    assert (terminated.sequence, terminated.upper_cased, terminated.stop_stripped) == (
        "MPEPTIDEK",
        False,
        True,
    )
    assert doubly_terminated.sequence == "MPEPTIDEK*"


def test_normalization_is_fixed_instead_of_switchable() -> None:
    assert normalize_sequence("MPEPTIDEK*").sequence == "MPEPTIDEK"


def test_tolerated_residues_pass_and_everything_else_is_reported() -> None:
    rules = _diagnostics().rules

    assert rules.illegal_residues("MXBZUOJPEPTIDEK") == ""
    assert rules.illegal_residues("MPEP*TIDEK") == "*"
    assert rules.illegal_residues("MPEP-TIDE.K1") == "-.1"
    assert rules.illegal_residues("mpeptidek") == "mpeptidek"


def test_an_internal_stop_is_explained_as_what_it_means() -> None:
    assert "translated nucleotide output" in describe_illegal_residues("*")
    assert describe_illegal_residues("-.") == "illegal sequence characters '-' '.'"


@pytest.mark.parametrize(
    ("entry_id", "expected"),
    [
        ("aa|p42261_db1_human|2026-06-05", "fgcz_sentinel"),
        ("sp|Cont_ALBU|ALBU_BOVIN", "fgcz_contaminant"),
        ("zh|C0012_KERATIN", "fgcz_contaminant"),
        ("sp|P02769|ALBU_BOVIN", "uniprot"),
        ("tr|A0A0B4J2D5|A0A0B4J2D5_HUMAN", "uniprot"),
        ("sp|J7HBH4.1|GLYC_GOGV", "uniprot"),
        ("sp|P02769-2|ALBU_BOVIN", "uniprot"),
        ("P02769", "uniprot_bare"),
        ("NP_000001.1", "refseq"),
        ("AKH49314.1", "genbank"),
        ("ENST00000028008.9", "ensembl"),
        ("NX_P02769", "nextprot"),
        ("pf|Pf2004_000005200|Pf2004_000005200", "pf|"),
        ("poi|pep_0001|synthetic", "poi|"),
        ("Cluster-838.148120;orf1", UNMATCHED_NAMESPACE),
    ],
)
def test_identifier_namespaces_compose_shared_and_fgcz_rules(entry_id: str, expected: str) -> None:
    namespace, _ = _diagnostics().rules.diagnose_identifier(entry_id)

    assert namespace == expected


def test_a_decoy_is_diagnosed_as_what_it_was_made_from() -> None:
    namespace, classifications = _diagnostics().rules.diagnose_identifier(
        "REV_sp|P02769|ALBU_BOVIN"
    )

    assert namespace == "uniprot"
    assert classifications == frozenset({"decoy"})


def test_invalid_document_reports_its_path(tmp_path: Path) -> None:
    path = tmp_path / "diagnostics.json"
    path.write_text("{}", encoding="utf-8")

    with pytest.raises(ValueError, match=str(path)):
        load_registry_diagnostics(path)
