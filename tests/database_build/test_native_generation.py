"""Scientific characterization of native decoy and entrapment generation."""

from __future__ import annotations

import random
from collections import Counter
from pathlib import Path

import pytest

from protein_fasta.analytics.digestion import (
    Digestion,
    digest_segments,
    digest_sequence,
    peptide_universe,
)
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.collisions import shuffled_candidate
from protein_fasta.database.decoy_generation import (
    make_decoypyrat_generation,
    make_shuffle_decoy_generation,
    reverse_and_switch,
)
from protein_fasta.database.entrapment_generation import (
    make_foreign_species_entrapment_generation,
    make_shuffled_entrapment_generation,
)
from protein_fasta.reading.parser import read_records
from protein_fasta.schema.analytics import DigestionDocument

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "fdrbench_v1_1_1"


def _digestion(
    *,
    minimum_length: int = 2,
    maximum_length: int = 50,
    missed_cleavages: int = 0,
) -> Digestion:
    return make_digestion(
        DigestionDocument(
            min_length=minimum_length,
            max_length=maximum_length,
            missed_cleavages=missed_cleavages,
        )
    )


def test_digestion_segments_and_collision_candidate_preserve_boundaries() -> None:
    digestion = _digestion()
    assert digest_segments("MPEPTIDEKTESTAR", digestion) == ("MPEPTIDEK", "TESTAR")
    assert peptide_universe(iter(("MPEPTIDEKTESTAR", "PEPTIDEK")), digestion) == frozenset(
        ("MPEPTIDEK", "TESTAR", "PEPTIDEK")
    )

    candidate, attempts = shuffled_candidate(
        "LKPPR",
        rng=random.Random(1),
        forbidden=frozenset(("LKPPR", "AAAK", "RSPKPL")),
        fix_c_term=True,
        accepts=lambda value: len(digest_segments(value, digestion)) == 1,
        max_attempts=10,
    )
    assert (candidate, attempts) == ("LPKPR", 3)


def test_collision_candidate_reports_impossible_and_exhausted_searches() -> None:
    impossible = shuffled_candidate(
        "AAAA",
        rng=random.Random(1),
        forbidden=frozenset(("AAAA",)),
        max_attempts=3,
    )
    exhausted = shuffled_candidate(
        "ABA",
        rng=random.Random(1),
        forbidden=frozenset(("AAB", "BAA")),
        max_attempts=5,
    )
    assert impossible == (None, 0)
    assert exhausted == (None, 5)


def test_whole_protein_shuffle_is_deterministic_and_isolates_random_state() -> None:
    entries = (("p1 first", "MPEPTIDEK"), ("p2", "TESTSEQUENCE"))
    generation = make_shuffle_decoy_generation(17)

    random.seed(71)
    expected_global = random.random()
    random.seed(71)
    first = generation.generate(entries, prefix="REV_")
    second = generation.generate(entries, prefix="REV_")

    assert first == second
    assert [description for description, _ in first.entries] == ["REV_p1 first", "REV_p2"]
    assert random.random() == expected_global
    for (_, source), (_, generated) in zip(entries, first.entries, strict=True):
        assert Counter(source) == Counter(generated)


def test_whole_protein_shuffle_seed_and_unique_ids_are_explicit() -> None:
    entries = (("p", "ABCDEFGHIJKLMNOP"),)
    first = make_shuffle_decoy_generation(1).generate(entries, prefix="REV_")
    second = make_shuffle_decoy_generation(2).generate(entries, prefix="REV_")
    assert first.entries != second.entries
    with pytest.raises(ValueError, match="unique"):
        make_shuffle_decoy_generation(1).generate(
            (("p first", "AAAA"), ("p second", "BBBB")),
            prefix="REV_",
        )


def test_reverse_and_switch_matches_reviewed_pgatk_behavior() -> None:
    assert reverse_and_switch("PEPTIDEK") == "KEDITPEP"
    assert reverse_and_switch("PEPTIDEK", switch=False) == "KEDITPEP"
    assert reverse_and_switch("AKR") == "KRA"
    assert reverse_and_switch("AKR", switch=False) == "RKA"


def test_decoypyrat_resolves_collisions_deterministically() -> None:
    generation = make_decoypyrat_generation(seed=1, digestion=_digestion())
    first = generation.generate((("p", "AKPA"),), prefix="DECOY_")
    second = generation.generate((("p", "AKPA"),), prefix="DECOY_")

    assert first == second
    assert first.entries == (("DECOY_p", "KPAA"),)
    assert first.initial_collisions == 1
    assert first.unresolved_collisions == 0


def test_decoypyrat_drops_unreplaceable_segments_and_empty_decoys() -> None:
    generation = make_decoypyrat_generation(seed=1, digestion=_digestion())
    partial = generation.generate(
        (("p", "PEPTIDEKAAAA"), ("t", "AAAK")),
        prefix="DECOY_",
    )
    empty = generation.generate((("p", "AAAA"),), prefix="DECOY_")

    assert partial.entries[0] == ("DECOY_p", "AEDITPEP")
    assert partial.dropped_peptides == 1
    assert partial.unresolved_collisions == 0
    assert empty.entries == ()
    assert empty.dropped_peptides == 1
    assert empty.omitted_decoys == 1


def test_decoypyrat_replacement_preserves_neighbouring_cleavage() -> None:
    generation = make_decoypyrat_generation(seed=1, digestion=_digestion())
    result = generation.generate(
        (("p", "AAAKRSPKPL"), ("t", "LKPPR")),
        prefix="DECOY_",
    )

    assert result.entries[0][1] == "LPKPRKSAAA"
    assert digest_segments("LPKPR", _digestion()) == ("LPKPR",)
    assert result.initial_collisions == 1
    assert result.unresolved_collisions == 0


def test_shuffled_entrapment_is_deterministic_paired_and_composition_preserving() -> None:
    generation = make_shuffled_entrapment_generation(
        fold=2,
        seed=1,
        digestion=_digestion(),
        normalize_i_to_l=False,
        fix_peptide_n_term=True,
        fix_peptide_c_term=True,
    )
    entries = (("p protein", "MPEPTIDEKTESTAR"),)
    first = generation.generate(entries)
    second = generation.generate(entries)

    assert first == second
    assert [header.partition(" ")[0] for header, _ in first.entries] == [
        "p_0_p_target",
        "p_1_p_target",
    ]
    assert first.achieved_fold == 2
    assert len(first.peptide_pairs) == 4
    for _, sequence in first.entries:
        assert Counter(sequence) == Counter(entries[0][1])


def test_shuffled_entrapment_caches_shared_peptides_and_reserves_candidates() -> None:
    generation = make_shuffled_entrapment_generation(
        fold=3,
        seed=1,
        digestion=_digestion(),
        normalize_i_to_l=False,
        fix_peptide_n_term=True,
        fix_peptide_c_term=True,
    )
    result = generation.generate(
        (
            ("p1", "MPEPTIDEKTESTARSAMPLERQVNTLPWR"),
            ("p2", "MSAMPLERGGGYTFDKQVNTLPWRTESTAR"),
        )
    )

    mapping = {
        (pair.target_peptide, pair.fold_index): pair.generated_peptide
        for pair in result.peptide_pairs
    }
    assert len(set(mapping.values())) == len(mapping)
    assert not set(mapping.values()) & {
        "MPEPTIDEK",
        "TESTAR",
        "SAMPLER",
        "QVNTLPWR",
        "GGGYTFDK",
    }
    assert result.achieved_fold == 3


def test_shuffled_entrapment_partial_generation_and_normalization_are_explicit() -> None:
    generation = make_shuffled_entrapment_generation(
        fold=1,
        seed=1,
        digestion=_digestion(),
        normalize_i_to_l=True,
        fix_peptide_n_term=True,
        fix_peptide_c_term=True,
    )
    normalized = generation.normalize((("p", "MILEK"),))
    exhausted = generation.generate((("p", "ABCDEFGKQRSTUVWR"),))

    assert normalized == (("p", "MLLEK"),)
    assert exhausted.entries == ()
    assert exhausted.achieved_fold == 0
    assert exhausted.failures == 1


def test_foreign_species_filters_shared_proteins_and_samples_deterministically() -> None:
    generation = make_foreign_species_entrapment_generation(
        fold=1,
        seed=3,
        digestion=_digestion(),
        normalize_i_to_l=False,
        reject_shared_foreign=True,
    )
    targets = (("target", "MPEPTIDEK"),)
    foreign = (
        ("shared", "MPEPTIDEK"),
        ("f1", "ABCDEFGK"),
        ("f2", "QRSTUVWR"),
    )
    first = generation.generate(targets, foreign_entries=foreign)
    second = generation.generate(targets, foreign_entries=foreign)

    assert first == second
    assert len(first.entries) == 1
    assert first.entries[0][0].partition(" ")[0].endswith("_p_target")
    assert first.failures == 1
    assert first.achieved_fold == 1


def test_foreign_species_can_allow_overlap_and_emit_a_partial_fold() -> None:
    allow_overlap = make_foreign_species_entrapment_generation(
        fold=1,
        seed=1,
        digestion=_digestion(),
        normalize_i_to_l=True,
        reject_shared_foreign=False,
    )
    targets = allow_overlap.normalize((("target", "MILEK"),))
    foreign = allow_overlap.normalize((("foreign protein", "MILEK"),))
    overlap = allow_overlap.generate(targets, foreign_entries=foreign)

    partial = make_foreign_species_entrapment_generation(
        fold=1,
        seed=1,
        digestion=_digestion(),
        normalize_i_to_l=False,
        reject_shared_foreign=True,
    ).generate(
        (("t1", "AAAA"), ("t2", "BBBB")),
        foreign_entries=(("f1", "CCCC"),),
    )

    assert overlap.entries[0][1] == "MLLEK"
    assert overlap.failures == 0
    assert len(partial.entries) == 1
    assert partial.achieved_fold == 0


def test_shuffled_entrapment_matches_fdrbench_semantics() -> None:
    source = tuple(
        (record.raw_header, record.sequence) for record in read_records(FIXTURE_DIR / "input.fasta")
    )
    reference = tuple(read_records(FIXTURE_DIR / "expected.fasta"))
    reference_entrapments = reference[1::2]
    digestion = _digestion(minimum_length=7, maximum_length=35)
    generation = make_shuffled_entrapment_generation(
        fold=1,
        seed=2000,
        digestion=digestion,
        normalize_i_to_l=False,
        fix_peptide_n_term=True,
        fix_peptide_c_term=True,
    )
    result = generation.generate(source)

    assert [header.partition(" ")[0] for header, _ in result.entries] == [
        record.raw_header.partition(" ")[0] for record in reference_entrapments
    ]
    target_peptides = peptide_universe((sequence for _, sequence in source), digestion)
    for (_, target), (_, generated), reference_record in zip(
        source,
        result.entries,
        reference_entrapments,
        strict=True,
    ):
        assert Counter(generated) == Counter(target)
        assert Counter(reference_record.sequence) == Counter(target)
        target_segments = digest_segments(target, digestion)
        generated_segments = digest_segments(generated, digestion)
        reference_segments = digest_segments(reference_record.sequence, digestion)
        for source_segment, generated_segment, reference_segment in zip(
            target_segments,
            generated_segments,
            reference_segments,
            strict=True,
        ):
            assert generated_segment[0] == source_segment[0]
            assert generated_segment[-1] == source_segment[-1]
            assert reference_segment[0] == source_segment[0]
            assert reference_segment[-1] == source_segment[-1]
        assert target_peptides.isdisjoint(
            peptide.sequence for peptide in digest_sequence(generated, digestion)
        )
