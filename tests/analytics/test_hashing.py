"""Tests for versioned analytical hashes."""

import hashlib
from pathlib import Path

from protein_fasta.analytics.hashing import (
    CONTENT_FINGERPRINT_VERSION,
    DESCRIPTION_SET_FINGERPRINT_VERSION,
    FILE_CHECKSUM_VERSION,
    ID_SET_FINGERPRINT_VERSION,
    PEPTIDE_HASH_VERSION,
    SEQUENCE_HASH_VERSION,
    content_fingerprint,
    description_set_fingerprint,
    file_checksum,
    id_set_fingerprint,
    peptide_hash,
    sequence_hash,
)


def test_sequence_hash_is_blake2b_128_without_hidden_normalization() -> None:
    expected = hashlib.blake2b(b"PEPTIDE*", digest_size=16).digest()

    assert sequence_hash("PEPTIDE*") == expected
    assert len(sequence_hash("PEPTIDE*")) == 16
    assert sequence_hash("PEPTIDE") != sequence_hash("peptide")
    assert SEQUENCE_HASH_VERSION == "blake2b-128-sequence-v1"


def test_peptide_hash_has_its_own_semantic_version() -> None:
    assert peptide_hash("PEPTIDE") == hashlib.blake2b(b"PEPTIDE", digest_size=16).digest()
    assert PEPTIDE_HASH_VERSION == "blake2b-128-peptide-v1"


def test_file_checksum_is_md5_over_exact_bytes(tmp_path: Path) -> None:
    path = tmp_path / "database.fasta"
    content = b">P1 mixed case\nAa*\n"
    path.write_bytes(content)

    assert file_checksum(path) == hashlib.md5(content, usedforsecurity=False).hexdigest()
    assert FILE_CHECKSUM_VERSION == "md5-file-v1"


def test_registry_fingerprints_use_blake2b_128_with_distinct_versions() -> None:
    sequence = sequence_hash("PEPTIDE")

    assert id_set_fingerprint(["P1", "P2"]).startswith("blake2b-128:")
    assert content_fingerprint([("P1", sequence)]).startswith("blake2b-128:")
    assert description_set_fingerprint([peptide_hash("description")]).startswith("blake2b-128:")
    assert ID_SET_FINGERPRINT_VERSION == "blake2b-128-id-set-v1"
    assert CONTENT_FINGERPRINT_VERSION == "blake2b-128-id-sequence-pairs-v1"
    assert DESCRIPTION_SET_FINGERPRINT_VERSION == "blake2b-128-description-set-v1"
