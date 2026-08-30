from __future__ import annotations

from typing import cast

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.decoy import DecoyMode, make_decoy, reverse_sequence
from protein_fasta.decoy_compile import make_decoy_generation
from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.decoy import DecoyPyratDocument, ShuffleDecoyDocument
from protein_fasta.summary import SummaryAccumulator, summarize_sequences


def test_reverse_is_plain_reversal() -> None:
    assert reverse_sequence("PEPTIDE") == "EDITPEP"


def test_make_decoy_prefixes_description_and_reverses_sequence() -> None:
    desc, seq = make_decoy("sp|P12345|FOO_HUMAN Some protein", "PEPTIDE", prefix="REV_")
    assert desc == "REV_sp|P12345|FOO_HUMAN Some protein"
    assert seq == "EDITPEP"


def test_shuffle_mode_preserves_composition() -> None:
    original = "PEPTIDEKR"
    batch = make_decoy_generation(ShuffleDecoyDocument(seed=7)).generate(
        (("sp|P1|ONE protein", original),), prefix="REV_"
    )
    description, shuffled = batch.entries[0]
    assert description == "REV_sp|P1|ONE protein"
    assert sorted(shuffled) == sorted(original)
    assert batch.parameters["seed"] == 7


def test_compiled_shuffle_owns_generation_provenance_and_annotation() -> None:
    generation = make_decoy_generation(ShuffleDecoyDocument(seed=7))

    batch = generation.generate((("sp|P1|ONE protein", "PEPTIDEKR"),), prefix="REV_")

    assert generation.mode is DecoyMode.SHUFFLE
    assert generation.seed == 7
    assert batch.parameters == generation.parameters()
    assert "decoys shuffle seed 7 with protein_fasta" in generation.annotation()


def test_decoypyrat_preserves_headers_and_reports_collisions() -> None:
    batch = make_decoy_generation(DecoyPyratDocument(seed=3)).generate(
        (("sp|P1|ONE protein", "MIPK"),), prefix="REV_"
    )
    assert batch.entries[0][0] == "REV_sp|P1|ONE protein"
    digestion = batch.parameters["digestion"]
    assert isinstance(digestion, dict)
    assert digestion["enzyme"] == make_digestion(DigestionDocument()).cleavage.pattern
    assert batch.unresolved_collisions == 0


def test_decoypyrat_collision_universe_matches_the_application_digester() -> None:
    """The collision check must see the peptides a search of this build would.

    Pyteomics' named ``trypsin`` also cleaves W-K-P and M-R-P, which the search
    engines these databases are built for do not, so the rule is passed
    explicitly rather than by name.
    """
    config = DigestionDocument(min_length=7, max_length=50, missed_cleavages=0)
    generation = make_decoy_generation(DecoyPyratDocument(seed=2000, digestion=config))
    parameters = generation.parameters()
    digestion = parameters["digestion"]
    assert isinstance(digestion, dict)
    values = cast(dict[str, object], digestion)
    missed_cleavages = values["missed_cleavages"]
    min_length = values["min_length"]
    max_length = values["max_length"]
    assert isinstance(missed_cleavages, int)
    assert isinstance(min_length, int)
    assert isinstance(max_length, int)
    compiled = make_digestion(config)
    assert values["enzyme"] == compiled.cleavage.pattern
    assert (min_length, max_length, missed_cleavages) == (7, 50, 0)


def test_decoypyrat_digestion_follows_the_configured_length_window() -> None:
    spec = DecoyPyratDocument(
        seed=2000,
        digestion=DigestionDocument(min_length=9, max_length=40, missed_cleavages=2),
    )
    generation = make_decoy_generation(spec)
    parameters = generation.parameters()
    digestion = parameters["digestion"]
    assert isinstance(digestion, dict)

    assert digestion["min_length"] == 9
    assert digestion["max_length"] == 40
    # DecoyPYrat replaces whole fully cleaved segments, so the collision
    # universe stays fully cleaved even when peptide analysis allows misses.
    assert digestion["missed_cleavages"] == 0
    assert "length 9-40" in generation.annotation()


def test_decoypyrat_drops_a_peptide_that_admits_no_replacement() -> None:
    batch = make_decoy_generation(DecoyPyratDocument(seed=2000)).generate(
        (("p", "PEPTIDEKAAAAAAA"), ("t", "AAAAAAK")), prefix="REV_"
    )

    # The decoy of p digests into AAAAAAK and AEDITPEP. AAAAAAK is a target and
    # its residues are one repeated letter, so no rearrangement differs from
    # itself; it is removed instead of shipped as a decoy that is a target,
    # and the neighbouring peptide is untouched.
    assert batch.entries[0] == ("REV_p", "AEDITPEP")
    assert batch.dropped_peptides == 1
    assert batch.unresolved_collisions == 0
    assert batch.omitted_decoys == 0


def test_decoypyrat_omits_a_decoy_left_with_no_sequence() -> None:
    batch = make_decoy_generation(DecoyPyratDocument(seed=2000)).generate(
        (("p", "AAAAAAA"),), prefix="REV_"
    )

    # The protein is a single unreplaceable peptide, so there is no decoy to
    # emit for it and the entry count no longer matches the input.
    assert batch.entries == ()
    assert batch.omitted_decoys == 1
    assert batch.dropped_peptides == 1


def test_summarize_basic() -> None:
    summary = summarize_sequences(["AAA", "CCCC", "DD"])
    assert summary.n_sequences == 3
    assert summary.length_min == 2
    assert summary.length_max == 4
    assert summary.length_mean == 3.0
    assert summary.total_residues == 9
    assert summary.aa_frequencies == {"A": 3, "C": 4, "D": 2}


def test_summarize_empty() -> None:
    summary = summarize_sequences([])
    assert summary.n_sequences == 0
    assert summary.total_residues == 0
    assert summary.aa_frequencies == {}


def test_summarize_uses_inclusive_quartiles() -> None:
    summary = summarize_sequences("A" * length for length in range(1, 9))

    assert summary.length_q1 == 2.75
    assert summary.length_median == 4.5
    assert summary.length_q3 == 6.25
    assert summary.total_residues == 36


def test_accumulator_returns_independent_incremental_snapshots() -> None:
    accumulator = SummaryAccumulator()
    accumulator.add("AA")
    first = accumulator.summary()

    accumulator.add("CCCCCCCCCC")
    second = accumulator.summary()

    assert first.n_sequences == 1
    assert first.length_q1 == 2.0
    assert first.total_residues == 2
    assert first.aa_frequencies == {"A": 2}
    assert second.n_sequences == 2
    assert second.length_q1 == 4.0
    assert second.length_median == 6.0
    assert second.length_q3 == 8.0
    assert second.total_residues == 12
    assert second.aa_frequencies == {"A": 2, "C": 10}


def test_accumulator_merges_without_losing_distribution_information() -> None:
    first = SummaryAccumulator()
    first.add("A")
    first.add("CCC")
    second = SummaryAccumulator()
    second.add("GG")
    second.add("TTTT")

    first.merge(second)
    summary = first.summary()

    assert summary.n_sequences == 4
    assert summary.length_min == 1
    assert summary.length_q1 == 1.75
    assert summary.length_median == 2.5
    assert summary.length_q3 == 3.25
    assert summary.length_max == 4
    assert summary.total_residues == 10
    assert summary.aa_frequencies == {"A": 1, "C": 3, "G": 2, "T": 4}


def test_accumulator_can_keep_exact_lengths_without_counting_amino_acids() -> None:
    accumulator = SummaryAccumulator()
    accumulator.add_length(3)

    summary = accumulator.summary()

    assert summary.n_sequences == 1
    assert summary.length_mean == 3
    assert summary.total_residues == 3
    assert summary.aa_frequencies == {}
