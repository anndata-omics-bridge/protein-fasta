"""Versioned fast hashes for protein-FASTA analytical values."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path

FILE_CHECKSUM_VERSION = "md5-file-v1"
"""Version label for exact-byte file checksums."""

SEQUENCE_HASH_VERSION = "blake2b-128-sequence-v1"
"""Version label for supplied protein-sequence hashes."""

PEPTIDE_HASH_VERSION = "blake2b-128-peptide-v1"
"""Version label for supplied peptide-sequence hashes."""

_DIGEST_SIZE = 16
_READ_SIZE = 1024 * 1024

ID_SET_FINGERPRINT_VERSION = "blake2b-128-id-set-v1"
CONTENT_FINGERPRINT_VERSION = "blake2b-128-id-sequence-pairs-v1"
DESCRIPTION_SET_FINGERPRINT_VERSION = "blake2b-128-description-set-v1"


def sequence_hash(sequence: str, /) -> bytes:
    """Hash exactly one supplied protein sequence with BLAKE2b-128.

    The caller owns normalization. This function deliberately does not change case,
    remove whitespace, or strip a terminal stop.
    """
    return hashlib.blake2b(sequence.encode("ascii"), digest_size=_DIGEST_SIZE).digest()


def peptide_hash(sequence: str, /) -> bytes:
    """Hash exactly one supplied peptide sequence with BLAKE2b-128."""
    return hashlib.blake2b(sequence.encode("ascii"), digest_size=_DIGEST_SIZE).digest()


def file_checksum(path: Path, /) -> str:
    """Return an MD5 checksum of exact file bytes for non-security provenance."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as source:
        for block in iter(lambda: source.read(_READ_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def id_set_fingerprint(ids: Iterable[str], /) -> str:
    """Fingerprint one ordered stream of distinct protein identifiers."""
    digest = hashlib.blake2b(digest_size=_DIGEST_SIZE)
    for identifier in ids:
        digest.update(identifier.encode())
        digest.update(b"\n")
    return "blake2b-128:" + digest.hexdigest()


def content_fingerprint(pairs: Iterable[tuple[str, bytes]], /) -> str:
    """Fingerprint ordered distinct identifier and sequence-hash pairs."""
    digest = hashlib.blake2b(digest_size=_DIGEST_SIZE)
    for identifier, hashed_sequence in pairs:
        digest.update(identifier.encode())
        digest.update(b"\0")
        digest.update(hashed_sequence)
        digest.update(b"\n")
    return "blake2b-128:" + digest.hexdigest()


def description_set_fingerprint(description_hashes: Iterable[bytes], /) -> str:
    """Fingerprint one ordered stream of distinct normalized-description hashes."""
    digest = hashlib.blake2b(digest_size=_DIGEST_SIZE)
    for description_hash in description_hashes:
        digest.update(description_hash)
        digest.update(b"\n")
    return "blake2b-128:" + digest.hexdigest()
