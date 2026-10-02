from __future__ import annotations

import pytest

from protein_fasta.analytics.peptide_properties import PeptideProperties, peptide_properties

# Reference values computed with the R package Peptides 2.4.6: mw(), pI(), hydrophobicity(scale =
# "KyteDoolittle"), instaIndex(), boman(), and charge(pH = 7, pKscale = "Sillero").
_R_PEPTIDES: dict[str, tuple[float, float, float, float, float, float]] = {
    "KQPWWR": (900.05064, 11.6517686762183, -2.55, 75.1, 3.55833333333333, 1.94039957245729),
    "AADAEAEVASLNRR": (
        1472.57784,
        4.42636699059199,
        -0.485714285714286,
        44.1357142857143,
        3.15785714285714,
        -1.05190885933289,
    ),
    "AAAFYVR": (
        796.92434,
        9.34881785998867,
        0.942857142857143,
        55.1428571428571,
        0.372857142857143,
        0.939798520202593,
    ),
    "CAAALMASR": (
        893.08694,
        8.54992124392649,
        1.12222222222222,
        44.0,
        0.281111111111111,
        0.930896531102584,
    ),
    "AADDPENGER": (1073.04094, 3.73790191433735, -2.04, -4.68, 4.806, -3.05089985843388),
    "EAAAMGPTK": (
        875.00754,
        6.40878378972325,
        -0.311111111111111,
        8.88888888888889,
        0.69,
        -0.05643811845945,
    ),
}


# Reference values from Pyteomics 4.7.5: achrom.calculate_RT(sequence, achrom.RCs_zubarev).
_PYTEOMICS_RETENTION = {
    "KQPWWR": 40.517363089835776,
    "AADAEAEVASLNRR": 39.800342365180306,
    "AAAFYVR": 42.71754169294273,
    "CAAALMASR": 39.291806905513575,
    "AADDPENGER": 23.300594882477434,
    "EAAAMGPTK": 28.35318945031029,
}


def _properties(sequence: str) -> PeptideProperties:
    properties = peptide_properties(sequence)
    assert properties is not None
    return properties


@pytest.mark.parametrize("sequence", sorted(_R_PEPTIDES))
def test_physicochemical_properties_match_r_peptides(sequence: str) -> None:
    weight, isoelectric, hydrophobicity, instability, boman, charge = _R_PEPTIDES[sequence]
    properties = _properties(sequence)

    assert properties.length == len(sequence)
    assert properties.molecular_weight == pytest.approx(weight, abs=1e-6)
    # R's root finder stops at about 1e-4; the bisection here runs to 1e-9.
    assert properties.isoelectric_point == pytest.approx(isoelectric, abs=1e-4)
    assert properties.hydrophobicity == pytest.approx(hydrophobicity, abs=1e-9)
    assert properties.instability_index == pytest.approx(instability, abs=1e-9)
    assert properties.boman_index == pytest.approx(boman, abs=1e-9)
    assert properties.charge == pytest.approx(charge, abs=1e-9)


@pytest.mark.parametrize("sequence", sorted(_PYTEOMICS_RETENTION))
def test_predicted_retention_time_matches_pyteomics(sequence: str) -> None:
    assert _properties(sequence).predicted_retention_time == pytest.approx(
        _PYTEOMICS_RETENTION[sequence], abs=1e-9
    )


def test_sequence_flags_and_counts() -> None:
    kqpwwr = _properties("KQPWWR")
    assert (kqpwwr.missed_cleavages, kqpwwr.proline_count, kqpwwr.c_terminal_residue) == (1, 1, "R")
    assert kqpwwr.contains_tryptophan
    assert not kqpwwr.contains_cysteine

    cysteine_first = _properties("CAAALMASR")
    assert cysteine_first.n_terminal_cysteine
    assert cysteine_first.contains_cysteine
    assert cysteine_first.contains_methionine

    motifs = _properties("AADDPENGER")
    assert motifs.contains_ng_motif
    assert motifs.contains_dp_motif

    glutamate_first = _properties("EAAAMGPTK")
    assert glutamate_first.n_terminal_glutamine_or_glutamate
    assert glutamate_first.c_terminal_residue == "K"


def test_missed_cleavages_ignore_the_c_terminus_and_proline() -> None:
    assert _properties("AKPAKAR").missed_cleavages == 1


def test_non_standard_residue_has_no_properties() -> None:
    assert peptide_properties("RSGASILQAGCUG") is None


@pytest.mark.parametrize("sequence", ["", "PEPtIDE", "PEP TIDE", "PEPT1DE"])
def test_malformed_sequence_is_rejected(sequence: str) -> None:
    with pytest.raises(ValueError, match="upper-case letters"):
        peptide_properties(sequence)
