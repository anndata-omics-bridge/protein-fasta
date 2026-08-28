from protein_fasta.build.generation.entrapment_types import (
    EntrapmentBatch,
    EntrapmentPeptidePair,
    format_entrapment_peptide_pairs,
)


def test_entrapment_evidence_is_dependency_free_and_serializes_stably() -> None:
    pair = EntrapmentPeptidePair(
        source_id="sp|P1|ONE",
        target_peptide="PEPTIDEK",
        generated_peptide="EDITPEPK",
        fold_index=0,
    )
    batch = EntrapmentBatch(
        entries=(("sp|P1|ONE_p_target", "EDITPEPK"),),
        peptide_pairs=(pair,),
        parameters={"strategy": "shuffled"},
        requested_fold=2,
        achieved_fold=1,
        failures=1,
        proteins_affected=1,
        source_proteins=2,
    )

    assert batch.is_complete is False
    assert batch.complete_proteins == 1
    assert format_entrapment_peptide_pairs(batch.peptide_pairs) == (
        "source_id\tfold_index\ttarget_peptide\tgenerated_peptide\n"
        "sp|P1|ONE\t0\tPEPTIDEK\tEDITPEPK\n"
    )


def test_complete_entrapment_batch_reports_complete() -> None:
    batch = EntrapmentBatch(
        entries=(),
        peptide_pairs=(),
        parameters={},
        requested_fold=1,
        achieved_fold=1,
        failures=0,
        proteins_affected=0,
        source_proteins=0,
    )

    assert batch.is_complete is True
