from __future__ import annotations

import pytest

from protein_fasta.reading.header import parse_header
from protein_fasta.registry.classification import (
    ContaminantBlockState,
    classify_identifier,
    classify_record,
)
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.registry.rules import RegistryDiagnosticRules, load_registry_diagnostics

CASES = [
    ("sp|P12345|FOO_HUMAN", EntryKind.TARGET),
    ("tr|B2C1U8|B2C1U8_THAPS", EntryKind.TARGET),
    ("ENST00000028008.9", EntryKind.TARGET),
    ("poi|pep_0001", EntryKind.TARGET),
    ("YP_009724390.1", EntryKind.TARGET),
    ("REV_sp|P12345|FOO_HUMAN", EntryKind.DECOY),
    ("REV_tr|B2C1U8|B2C1U8_THAPS", EntryKind.DECOY),
    ("REV_sp|Cont_CC261|IRTPROT_FGCZ", EntryKind.DECOY),
    ("sp|Cont_P00761|TRYP_PIG", EntryKind.CONTAMINANT),
    ("tr|Cont_Q29RJ0|Q29RJ0_BOVIN", EntryKind.CONTAMINANT),
    ("zh|C0001_P61626|LYSC_HUMAN", EntryKind.CONTAMINANT),
    ("aa|p42261_db1_9606wIsoforms_plus_custom", EntryKind.SENTINEL),
    ("aa|Cont_UniversalContaminants", EntryKind.SENTINEL),
    ("REV_aa|Cont_specialContaminants", EntryKind.SENTINEL),
]


def _diagnostics() -> RegistryDiagnosticRules:
    return load_registry_diagnostics()


@pytest.mark.parametrize(("entry_id", "expected"), CASES)
def test_configured_labels_resolve_to_operational_kinds(
    entry_id: str,
    expected: EntryKind,
) -> None:
    assert classify_identifier(entry_id, _diagnostics().rules) is expected


def test_stacked_decoy_contaminant_labels_are_both_preserved() -> None:
    _, classifications = _diagnostics().rules.diagnose_identifier("REV_sp|Cont_CC261|IRTPROT_FGCZ")

    assert classifications == frozenset({"contaminant", "decoy"})
    assert (
        classify_identifier("REV_sp|Cont_CC261|IRTPROT_FGCZ", _diagnostics().rules)
        is EntryKind.DECOY
    )


def test_parse_header_extracts_the_identifier_token() -> None:
    assert parse_header("sp|P12345|FOO OS=Homo sapiens").id == "sp|P12345|FOO"
    assert parse_header(">sp|P12345|FOO desc").id == "sp|P12345|FOO"
    assert parse_header("").id == ""


def test_entrapment_records_are_resolved_as_entrapment() -> None:
    diagnostics = _diagnostics()
    header = "sp|P1|ONE_p_target entrapment of sp|P1|ONE first protein"
    _, classifications = diagnostics.rules.diagnose_identifier(parse_header(header).id)

    kind, group, _ = classify_record(
        header,
        classifications,
        None,
        diagnostics.decoy_prefix,
    )

    assert kind is EntryKind.ENTRAPMENT
    assert group is None


def test_an_entrapment_after_a_contaminant_block_is_not_a_contaminant() -> None:
    diagnostics = _diagnostics()
    open_block = ContaminantBlockState("core_bottom_up", has_entries=True)
    header = "sp|P1|ONE_p_target entrapment of sp|P1|ONE"
    _, classifications = diagnostics.rules.diagnose_identifier(parse_header(header).id)

    kind, _, block_state = classify_record(
        header,
        classifications,
        open_block,
        diagnostics.decoy_prefix,
    )

    assert kind is EntryKind.ENTRAPMENT
    assert block_state is open_block


def test_a_decoy_of_an_entrapment_record_stays_a_decoy() -> None:
    diagnostics = _diagnostics()
    header = "REV_sp|P1|ONE_p_target entrapment of sp|P1|ONE"
    _, classifications = diagnostics.rules.diagnose_identifier(parse_header(header).id)

    kind, _, _ = classify_record(
        header,
        classifications,
        None,
        diagnostics.decoy_prefix,
    )

    assert classifications == frozenset({"decoy", "entrapment"})
    assert kind is EntryKind.DECOY
