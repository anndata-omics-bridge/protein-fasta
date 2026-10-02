"""Backend-free physicochemical properties of peptide sequences.

The definitions reproduce the defaults of the R package Peptides: average molecular weight, EMBOSS
isoelectric point, Kyte-Doolittle hydrophobicity, Guruprasad instability index, Boman index, and
net charge at pH 7 with the Sillero pK scale. The predicted retention time is the length-corrected
additive reversed-phase model of Goloborodko et al. (2010), as in Pyteomics ``achrom``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

_STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")

_AVERAGE_RESIDUE_MASS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
_WATER_AVERAGE_MASS = 18.01524

_KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Boman (2003) solubility scale; the index is the negated mean over the sequence.
_BOMAN = {
    "L": 4.92,
    "I": 4.92,
    "V": 4.04,
    "F": 2.98,
    "M": 2.35,
    "W": 2.33,
    "A": 1.81,
    "C": 1.28,
    "G": 0.94,
    "P": 0.0,
    "Y": -0.14,
    "T": -2.57,
    "S": -3.40,
    "H": -4.66,
    "Q": -5.54,
    "K": -5.55,
    "N": -6.64,
    "E": -6.81,
    "D": -8.72,
    "R": -14.92,
}

# Guruprasad (1990) dipeptide instability weight values (DIWV) by value; every dipeptide not
# listed weighs 1.0. Taken from Biopython's ProtParamData (BSD 3-Clause), with KQ (24.64) and
# FY (33.601) set to the 24.68 and 33.6 of the R package Peptides.
_INSTABILITY_WEIGHTS: dict[float, str] = {
    58.28: "MH RR RW",
    44.94: "AC EC HI HY IE MP MS NI QS RS SP YM",
    33.6: "CH CM CT EE FY KM KR LQ SC",
    24.68: "CW HK HN KQ LW MY NK WH WM YA YD",
    20.26: (
        "AP CD CL CP DS ED EI EP EQ ES FP IL LP LR PA PF PP PQ PS PV QD QE QP QQ RH RP RQ "
        "SE SQ SR SS TE VP"
    ),
    18.38: "PE",
    13.34: "FD GG GW IH MA RN TF WL WN YH YP YY",
    -1.88: "HP HW IP MM MT NC NP PW VK",
    -6.54: "CQ CV DF DR EH GE HT KP MQ MR NQ PC PD PM PR QC QF QV QY RY TQ VY YE",
    -7.49: "AD AH DK GA GI GK GN GT GY IK IV KG KI KL KV LK NT RG TG VG VT WV YG YT",
    -9.37: "HF HG NW WG YW",
    -14.03: "DT EW FK NF NG TN TW VD WA WT",
    -15.91: "YR",
}
_DIPEPTIDE_INSTABILITY = {
    dipeptide: weight
    for weight, dipeptides in _INSTABILITY_WEIGHTS.items()
    for dipeptide in dipeptides.split()
}

# Reversed-phase retention coefficients (minutes) of Goloborodko et al. (2010), Rapid Commun Mass
# Spectrom 24(4):454-62: Reprosil-Pur C18-AQ, 0.5% acetic acid, 0.5% water/min, carbamidomethyl
# Cys. Taken from Pyteomics 4.7.5 ``achrom.RCs_zubarev`` (Apache License 2.0).
_RETENTION_COEFFICIENT = {
    "A": 6.73,
    "C": 3.25,
    "D": 5.64,
    "E": 5.66,
    "F": 27.43,
    "G": 2.35,
    "H": -0.66,
    "I": 20.50,
    "K": -4.47,
    "L": 23.38,
    "M": 17.39,
    "N": 2.57,
    "P": 5.66,
    "Q": 2.93,
    "R": -2.55,
    "S": 3.58,
    "T": 4.88,
    "V": 13.05,
    "W": 31.27,
    "Y": 13.22,
}
_RETENTION_LENGTH_CORRECTION = -0.21
_RETENTION_SHIFT = 0.53

_MISSED_CLEAVAGE = re.compile(r"[KR](?!P)")


@dataclass(frozen=True, slots=True)
class _PkScale:
    n_terminus: float
    c_terminus: float
    lysine: float
    arginine: float
    histidine: float
    aspartate: float
    glutamate: float
    cysteine: float
    tyrosine: float


_EMBOSS = _PkScale(8.6, 3.6, 10.8, 12.5, 6.5, 3.9, 4.1, 8.5, 10.1)
_SILLERO = _PkScale(8.2, 3.2, 10.4, 12.0, 6.4, 4.0, 4.5, 9.0, 10.0)


@dataclass(frozen=True, slots=True)
class PeptideProperties:
    """Sequence-derived properties of one peptide.

    Attributes:
        length: Number of residues.
        molecular_weight: Average molecular weight in Da.
        isoelectric_point: pH of zero net charge, EMBOSS pK scale.
        hydrophobicity: Mean Kyte-Doolittle hydropathy (GRAVY).
        instability_index: Guruprasad instability index.
        boman_index: Boman protein-binding potential index.
        charge: Net charge at pH 7, Sillero pK scale.
        predicted_retention_time: Reversed-phase retention in minutes on the Goloborodko et al.
            (2010) reference gradient; comparable between peptides, not across LC setups.
        missed_cleavages: Internal K or R not followed by P.
        proline_count: Number of prolines.
        c_terminal_residue: Last residue.
        contains_cysteine: The sequence contains C.
        contains_methionine: The sequence contains M.
        contains_tryptophan: The sequence contains W.
        n_terminal_glutamine_or_glutamate: The first residue is Q or E (pyro-Glu formation).
        n_terminal_cysteine: The first residue is C.
        contains_ng_motif: The sequence contains NG (deamidation-prone).
        contains_dp_motif: The sequence contains DP (acid-labile).
    """

    length: int
    molecular_weight: float
    isoelectric_point: float
    hydrophobicity: float
    instability_index: float
    boman_index: float
    charge: float
    predicted_retention_time: float
    missed_cleavages: int
    proline_count: int
    c_terminal_residue: str
    contains_cysteine: bool
    contains_methionine: bool
    contains_tryptophan: bool
    n_terminal_glutamine_or_glutamate: bool
    n_terminal_cysteine: bool
    contains_ng_motif: bool
    contains_dp_motif: bool


def peptide_properties(sequence: str, /) -> PeptideProperties | None:
    """Compute the properties of one stripped, upper-case peptide sequence.

    Args:
        sequence: Non-empty peptide of upper-case ASCII residue letters.

    Returns:
        The properties, or ``None`` when the sequence contains a residue other than the 20
        standard amino acids, such as selenocysteine ``U``, for which they are undefined.

    Raises:
        ValueError: If the sequence is empty or contains anything but upper-case ASCII letters.
    """
    if not (sequence.isascii() and sequence.isalpha() and sequence.isupper()):
        raise ValueError(f"peptide sequence must be non-empty upper-case letters, got {sequence!r}")
    if not _STANDARD_RESIDUES.issuperset(sequence):
        return None
    length = len(sequence)
    return PeptideProperties(
        length=length,
        molecular_weight=sum(_AVERAGE_RESIDUE_MASS[residue] for residue in sequence)
        + _WATER_AVERAGE_MASS,
        isoelectric_point=_isoelectric_point(sequence, _EMBOSS),
        hydrophobicity=sum(_KYTE_DOOLITTLE[residue] for residue in sequence) / length,
        instability_index=10
        / length
        * sum(
            _DIPEPTIDE_INSTABILITY.get(sequence[index : index + 2], 1.0)
            for index in range(length - 1)
        ),
        boman_index=-sum(_BOMAN[residue] for residue in sequence) / length,
        charge=_net_charge(sequence, 7.0, _SILLERO),
        predicted_retention_time=sum(_RETENTION_COEFFICIENT[residue] for residue in sequence)
        * (1 + _RETENTION_LENGTH_CORRECTION * math.log(length))
        + _RETENTION_SHIFT,
        missed_cleavages=len(_MISSED_CLEAVAGE.findall(sequence[:-1])),
        proline_count=sequence.count("P"),
        c_terminal_residue=sequence[-1],
        contains_cysteine="C" in sequence,
        contains_methionine="M" in sequence,
        contains_tryptophan="W" in sequence,
        n_terminal_glutamine_or_glutamate=sequence[0] in "QE",
        n_terminal_cysteine=sequence[0] == "C",
        contains_ng_motif="NG" in sequence,
        contains_dp_motif="DP" in sequence,
    )


def _net_charge(sequence: str, ph: float, scale: _PkScale) -> float:
    def positive(count: int, pk: float) -> float:
        return count / (1 + 10 ** (ph - pk))

    def negative(count: int, pk: float) -> float:
        return count / (1 + 10 ** (pk - ph))

    return (
        positive(1, scale.n_terminus)
        + positive(sequence.count("K"), scale.lysine)
        + positive(sequence.count("R"), scale.arginine)
        + positive(sequence.count("H"), scale.histidine)
        - negative(1, scale.c_terminus)
        - negative(sequence.count("D"), scale.aspartate)
        - negative(sequence.count("E"), scale.glutamate)
        - negative(sequence.count("C"), scale.cysteine)
        - negative(sequence.count("Y"), scale.tyrosine)
    )


def _isoelectric_point(sequence: str, scale: _PkScale) -> float:
    # The net charge falls strictly with pH, so bisection over 0-14 finds its single root.
    low, high = 0.0, 14.0
    while high - low > 1e-9:
        middle = (low + high) / 2
        if _net_charge(sequence, middle, scale) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2
