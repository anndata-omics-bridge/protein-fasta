"""SQLite index of FASTA databases and bounded sequence-level contents.

The registry stores exact counts and length metadata for every FASTA. Databases
within the configured detail limit also store each entry's ordinal, classified
kind, length, ID token, normalized sequence checksum, and exact amino-acid
composition; larger databases use bounded composition sampling and explicitly
remain metadata-only. Full ID and checksum arrays stay in SQLite and are never
serialized into the Dash UI.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import random
import tempfile
import time
from collections import Counter
from collections.abc import Callable, Generator, Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from functools import partial
from pathlib import Path
from typing import Protocol

from loguru import logger

from protein_fasta.analytics.hashing import (
    CONTENT_FINGERPRINT_VERSION,
    DESCRIPTION_SET_FINGERPRINT_VERSION,
    ID_SET_FINGERPRINT_VERSION,
    SEQUENCE_HASH_VERSION,
    content_fingerprint,
    description_set_fingerprint,
    id_set_fingerprint,
    sequence_hash,
)
from protein_fasta.diagnostics.messages import describe_illegal_residues
from protein_fasta.diagnostics.runtime import UNMATCHED_NAMESPACE
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import FastaReadError, read_headers
from protein_fasta.record import iter_protein_diagnostics
from protein_fasta.registry import snapshots
from protein_fasta.registry.backend import factory
from protein_fasta.registry.backend.base import (
    RegistryConnection,
    RegistryIntegrityError,
    Row,
    TempTableSpec,
)
from protein_fasta.registry.classification import (
    ContaminantBlockState,
    classify_record,
    normalize_contaminant_group,
    sentinel_details,
)
from protein_fasta.registry.filenames import parse_filename
from protein_fasta.registry.kinds import DetailLevel, EntryKind
from protein_fasta.registry.metadata import parse_database_metadata
from protein_fasta.registry.pair_metrics import PairMetricSelection, pair_metric_counts
from protein_fasta.registry.rules import load_registry_diagnostics
from protein_fasta.schema.build import MetadataDocument, NamingDocument
from protein_fasta.summary import SummaryAccumulator

SCHEMA_VERSION = 11
CANDIDATE_TABLE = "candidate_entries"
_DESCRIPTION_HASH_VERSION = "blake2b-128-normalized-whitespace-v1"
_CLASSIFICATION_ALGORITHM_VERSION = "diagnostic-labels-marker-block-v2"
_FASTA_GLOBS = ("*.fasta", "*.fa", "*.fas", "*.fna", "*.fasta.gz", "*.fasta.bz2")
_INSERT_BATCH_SIZE = 2_000
_ENTRY_PROGRESS_LOG_SECONDS = 60.0
# Someone watching a browser needs a line far more often than an unattended sweep
# needs one in the log.
_ENTRY_PROGRESS_REPORT_SECONDS = 2.0
_ENTRY_PROGRESS_CHECK_INTERVAL = 1_000
_GIBIBYTE = 1024**3
_PAIR_INSERT_BATCH_SIZE = 10_000
_PAIR_REPRESENTATIVES_TABLE = "bulk_pair_representatives"
_MATERIALIZED_PAIR_KINDS = (EntryKind.TARGET, EntryKind.CONTAMINANT)


class RegistryBackendSelection(Protocol):
    """Storage-engine selection required by registry snapshot operations."""

    @property
    def backend(self) -> str:
        """Configured backend name."""
        ...


class RegistrySettings(Protocol):
    """Small structural boundary for registry indexing configuration."""

    @property
    def fasta_root(self) -> Path: ...

    @property
    def registry_dir(self) -> Path: ...

    @property
    def registry(self) -> RegistryBackendSelection: ...

    @property
    def max_fasta_file_size_gib(self) -> float: ...

    @property
    def max_detailed_entries(self) -> int: ...

    @property
    def metadata_aa_sample_size(self) -> int: ...

    @property
    def min_fasta_date(self) -> datetime.date | None: ...

    @property
    def registry_diagnostics_path(self) -> Path | None: ...

    @property
    def naming(self) -> NamingDocument: ...

    @property
    def sentinel(self) -> MetadataDocument: ...

    @property
    def overlap_threshold(self) -> float: ...


class RegistrySchemaError(ValueError):
    """A registry was written by a release with a different schema.

    Carries the remedy rather than leaving it to the reader. "Reindex required" is
    actively misleading to someone who is already running ``reindex``: the
    incremental sweep opens the existing registry and so hits this same error. Only
    a full rebuild helps, and it publishes a new dated file rather than replacing
    the one that cannot be read.

    A ``ValueError`` because every caller already treats an unreadable registry as
    one; what this adds is a type the GUI can recognise as an expected, actionable
    state instead of logging a traceback for it.
    """

    def __init__(self, found: int, *, path: Path | None = None) -> None:
        self.found = found
        self.expected = SCHEMA_VERSION
        subject = f"Registry {path.name}" if path is not None else "Registry"
        super().__init__(
            f"{subject} has schema version {found}, but this release reads version {SCHEMA_VERSION}. "
            "Rebuild it with 'fasta-gen reindex --full' (Makefile: 'make reindex-full'). "
            "An incremental reindex cannot migrate it, and a full rebuild publishes a new dated "
            "registry while leaving the current one in place."
        )


class FastaValidationError(FastaReadError):
    """One entry breaks fasta_gen's strict ingest policy."""

    def __init__(self, path: Path, entry_id: str, reason: str) -> None:
        self.entry_id = entry_id
        super().__init__(str(path), f"{reason} ({entry_id})")


def _count_namespace(
    counts: Counter[str],
    namespace: str,
    max_namespaces: int,
) -> None:
    """Count one namespace while bounding one-off authority rows."""
    if namespace not in counts and len(counts) >= max_namespaces:
        counts[UNMATCHED_NAMESPACE] += 1
        return
    counts[namespace] += 1


def normalized_description_hash(raw_header: str) -> bytes | None:
    """Hash the normalized description according to the registry contract."""
    description = parse_header(raw_header).description
    if description is None:
        return None
    return hashlib.blake2b(description.encode(), digest_size=16).digest()


@dataclass(frozen=True, slots=True)
class RejectedFasta:
    """One FASTA rejected from a registry sweep without stopping the batch."""

    path: Path
    reason: str


@dataclass(frozen=True, slots=True)
class DatabaseKindStats:
    """Within-database facts for one classified entry kind."""

    kind: EntryKind
    entry_count: int = 0
    distinct_ids: int | None = 0
    distinct_sequences: int | None = 0
    distinct_descriptions: int | None = 0
    distinct_pairs: int | None = 0
    duplicate_id_occurrences: int | None = 0
    conflicting_ids: int | None = 0
    repeated_sequences: int | None = 0
    length_min: int = 0
    length_q1: float = 0.0
    length_median: float = 0.0
    length_mean: float = 0.0
    length_q3: float = 0.0
    length_max: int = 0
    total_residues: int = 0
    aa_sample_size: int = 0
    aa_counts: dict[str, int] = field(default_factory=dict[str, int])
    id_fingerprint: str | None = ""
    description_fingerprint: str | None = ""
    content_fingerprint: str | None = ""


@dataclass(frozen=True, slots=True)
class RegistryRecord:
    """Aggregate facts for one indexed or transient FASTA database."""

    id: int | None = None
    relative_path: str = ""
    filename: str = ""
    dbname: str = ""
    file_size_bytes: int = 0
    mtime_ns: int = 0
    sentinel_header: str | None = None
    annotation: str | None = None
    filename_is_decoy: bool = False
    is_decoy: bool = False
    contaminant_markers: list[str] = field(default_factory=list[str])
    indexed_at: str = ""
    detail_level: DetailLevel = DetailLevel.FULL
    entry_count: int = 0
    target_count: int = 0
    decoy_count: int = 0
    contaminant_count: int = 0
    entrapment_count: int = 0
    sentinel_count: int = 0
    distinct_target_ids: int | None = 0
    distinct_target_sequences: int | None = 0
    distinct_target_descriptions: int | None = 0
    duplicate_target_id_occurrences: int | None = 0
    conflicting_target_ids: int | None = 0
    repeated_target_sequences: int | None = 0
    length_min: int = 0
    length_q1: float = 0.0
    length_median: float = 0.0
    length_mean: float = 0.0
    length_q3: float = 0.0
    length_max: int = 0
    total_residues: int = 0
    aa_sample_size: int = 0
    aa_counts: dict[str, int] = field(default_factory=dict[str, int])
    target_id_fingerprint: str | None = ""
    target_description_fingerprint: str | None = ""
    target_content_fingerprint: str | None = ""
    # Validation evidence. The counts describe the difference between what is on
    # disk and what was indexed, so normalizing at scan time does not hide it.
    upper_cased_entries: int = 0
    stop_stripped_entries: int = 0
    illegal_residue_entries: int = 0
    illegal_residues: dict[str, int] = field(default_factory=dict[str, int])
    empty_sequence_entries: int = 0
    bare_identifier_entries: int = 0
    id_namespaces: dict[str, int] = field(default_factory=dict[str, int])
    kind_stats: dict[EntryKind, DatabaseKindStats] = field(
        default_factory=dict[EntryKind, DatabaseKindStats]
    )

    @property
    def dominant_id_namespace(self) -> str:
        """Name the namespace most of this database's identifiers belong to.

        A tie prefers a namespace that was actually recognised: "unmatched" is the
        absence of an answer, so it should never win by alphabetical accident.
        """
        if not self.id_namespaces:
            return ""
        return min(
            self.id_namespaces.items(),
            key=lambda item: (-item[1], item[0] == UNMATCHED_NAMESPACE, item[0]),
        )[0]

    @property
    def unmatched_id_entries(self) -> int:
        """Count identifiers belonging to no known namespace and no authority."""
        return self.id_namespaces.get(UNMATCHED_NAMESPACE, 0)

    @property
    def n_ids(self) -> int | None:
        """Compatibility alias for the number of distinct target IDs."""
        return self.distinct_target_ids


def _kind_summaries() -> dict[EntryKind, SummaryAccumulator]:
    return {kind: SummaryAccumulator() for kind in EntryKind}


def _kind_aa_counts() -> dict[EntryKind, Counter[str]]:
    return {kind: Counter() for kind in EntryKind}


@dataclass(slots=True)
class _AminoAcidReservoir:
    """Keep a deterministic uniform sample of sequences from one FASTA scan."""

    sample_size: int
    entries: list[tuple[EntryKind, str]] = field(default_factory=list[tuple[EntryKind, str]])
    seen_entries: int = 0
    random_source: random.Random = field(default_factory=lambda: random.Random(0))

    def add(self, kind: EntryKind, sequence: str) -> None:
        """Consider one sequence for standard reservoir sampling."""
        self.seen_entries += 1
        if len(self.entries) < self.sample_size:
            self.entries.append((kind, sequence))
            return
        replacement = self.random_source.randrange(self.seen_entries)
        if replacement < self.sample_size:
            self.entries[replacement] = (kind, sequence)

    def summarize(self) -> tuple[dict[EntryKind, Counter[str]], Counter[EntryKind]]:
        """Count sampled sequences and residues, grouped by entry kind."""
        counts = _kind_aa_counts()
        sample_sizes: Counter[EntryKind] = Counter()
        for kind, sequence in self.entries:
            sample_sizes[kind] += 1
            counts[kind].update(sequence)
        return counts, sample_sizes


@dataclass(slots=True)
class _ScanFacts:
    """Facts collected while streaming one or more FASTA files."""

    filenames: list[str] = field(default_factory=list[str])
    file_size_bytes: int = 0
    mtime_ns: int = 0
    sentinel_header: str | None = None
    annotation: str | None = None
    contaminant_markers: list[str] = field(default_factory=list[str])
    detail_level: DetailLevel = DetailLevel.FULL
    counts: Counter[EntryKind] = field(default_factory=Counter[EntryKind])
    kind_summaries: dict[EntryKind, SummaryAccumulator] = field(default_factory=_kind_summaries)
    sampled_aa_counts: dict[EntryKind, Counter[str]] | None = None
    sampled_aa_sizes: Counter[EntryKind] | None = None
    upper_cased_entries: int = 0
    stop_stripped_entries: int = 0
    illegal_residue_entries: int = 0
    illegal_residues: Counter[str] = field(default_factory=Counter[str])
    empty_sequence_entries: int = 0
    bare_identifier_entries: int = 0
    id_namespaces: Counter[str] = field(default_factory=Counter[str])

    def merge(self, other: _ScanFacts) -> None:  # noqa: C901
        """Merge facts collected from another file without reparsing it."""
        self.filenames.extend(other.filenames)
        self.upper_cased_entries += other.upper_cased_entries
        self.stop_stripped_entries += other.stop_stripped_entries
        self.illegal_residue_entries += other.illegal_residue_entries
        self.illegal_residues.update(other.illegal_residues)
        self.empty_sequence_entries += other.empty_sequence_entries
        self.bare_identifier_entries += other.bare_identifier_entries
        self.id_namespaces.update(other.id_namespaces)
        self.file_size_bytes += other.file_size_bytes
        self.mtime_ns = max(self.mtime_ns, other.mtime_ns)
        if self.sentinel_header is None:
            self.sentinel_header = other.sentinel_header
            self.annotation = other.annotation
        for marker in other.contaminant_markers:
            if marker not in self.contaminant_markers:
                self.contaminant_markers.append(marker)
        if other.detail_level is DetailLevel.METADATA_ONLY:
            self.detail_level = DetailLevel.METADATA_ONLY
        self.counts.update(other.counts)
        for kind in EntryKind:
            self.kind_summaries[kind].merge(other.kind_summaries[kind])
        if other.sampled_aa_counts is not None:
            if self.sampled_aa_counts is None:
                self.sampled_aa_counts = _kind_aa_counts()
            for kind in EntryKind:
                self.sampled_aa_counts[kind].update(other.sampled_aa_counts[kind])
        if other.sampled_aa_sizes is not None:
            if self.sampled_aa_sizes is None:
                self.sampled_aa_sizes = Counter()
            self.sampled_aa_sizes.update(other.sampled_aa_sizes)


@dataclass(frozen=True, slots=True)
class _KindContentGroup:
    """Databases with identical distinct ID-and-sequence content for one kind."""

    representative_id: int
    database_ids: tuple[int, ...]
    distinct_ids: int
    distinct_sequences: int
    distinct_descriptions: int
    distinct_pairs: int


def _classification_fingerprint(settings: RegistrySettings) -> str:
    diagnostics = load_registry_diagnostics(settings.registry_diagnostics_path)
    payload = json.dumps(
        {
            "algorithm_version": _CLASSIFICATION_ALGORITHM_VERSION,
            "diagnostics_fingerprint": diagnostics.fingerprint,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return "blake2b-128:" + hashlib.blake2b(payload.encode(), digest_size=16).hexdigest()


@contextmanager
def connect_registry(
    path: Path,
    *,
    backend: str | None = None,
    read_only: bool = False,
) -> Generator[RegistryConnection]:
    """Open one configured registry connection for an operation or callback.

    Which engine opens it comes from the filename unless ``backend`` says
    otherwise, so a path names a readable registry regardless of what the
    configured default happens to be.
    """
    with factory.connect(path, backend=backend, read_only=read_only) as connection:
        yield connection


def _create_schema(connection: RegistryConnection) -> None:
    """Create the complete schema used by ordinary incremental operations."""
    connection.create_tables()
    connection.create_entry_indexes()
    connection.create_pair_indexes()


def _expected_meta(settings: RegistrySettings) -> dict[str, str]:
    """Return the contract a registry must agree with to be readable."""
    return {
        "schema_version": str(SCHEMA_VERSION),
        "sequence_hash_version": SEQUENCE_HASH_VERSION,
        "id_set_fingerprint_version": ID_SET_FINGERPRINT_VERSION,
        "content_fingerprint_version": CONTENT_FINGERPRINT_VERSION,
        "description_set_fingerprint_version": DESCRIPTION_SET_FINGERPRINT_VERSION,
        "description_hash_version": _DESCRIPTION_HASH_VERSION,
        "fasta_diagnostics_fingerprint": _classification_fingerprint(settings),
        "max_fasta_file_size_gib": str(settings.max_fasta_file_size_gib),
        "max_detailed_entries": str(settings.max_detailed_entries),
        "metadata_aa_sample_size": str(settings.metadata_aa_sample_size),
    }


def _stored_meta(connection: RegistryConnection) -> dict[str, str]:
    """Return the contract a registry records about itself."""
    return {
        str(row["key"]): str(row["value"])
        for row in connection.execute("SELECT key, value FROM registry_meta").fetchall()
    }


def _reconcile_meta(connection: RegistryConnection, settings: RegistrySettings) -> None:
    """Validate classification metadata and record the current schema contract."""
    expected = _expected_meta(settings)
    existing = _stored_meta(connection)
    database_count = connection.scalar("SELECT COUNT(*) FROM databases")
    for key, value in expected.items():
        previous = existing.get(key)
        if previous is not None and previous != value and database_count:
            raise ValueError(
                f"Registry {key} changed from {previous!r} to {value!r}; "
                "run 'fasta-gen reindex --full' before using this database."
            )
        connection.upsert_meta(key, value)
    connection.set_schema_version(SCHEMA_VERSION)
    connection.commit()


def _verify_meta(connection: RegistryConnection, settings: RegistrySettings) -> None:
    """Check the recorded contract without writing anything.

    The read path's half of :func:`_reconcile_meta`. Same refusals, no writes, so
    a reader needs no more than read access to the file -- which on an engine
    with one writer is the difference between the application blocking every
    other command and not.
    """
    expected = _expected_meta(settings)
    existing = _stored_meta(connection)
    for key, value in expected.items():
        previous = existing.get(key)
        if previous is not None and previous != value:
            raise ValueError(
                f"Registry {key} is {previous!r} but this release expects {value!r}; "
                "run 'fasta-gen reindex --full' before using this database."
            )


def initialize_registry(connection: RegistryConnection, settings: RegistrySettings) -> None:
    """Create or validate the versioned registry schema."""
    current_version = connection.schema_version()
    if current_version not in (0, SCHEMA_VERSION):
        raise RegistrySchemaError(current_version, path=connection.path)
    _create_schema(connection)
    _reconcile_meta(connection, settings)


@contextmanager
def open_registry(
    settings: RegistrySettings,
    path: Path | None = None,
) -> Generator[RegistryConnection]:
    """Open a published registry for reading, without writing to it.

    Without an explicit ``path`` this opens the newest dated registry in
    ``settings.registry_dir`` and fails when there is none, rather than creating
    an empty one that would then report "no databases".

    Opened read-only, and it used to run the whole DDL script and commit a
    transaction on every page render. That was harmless on SQLite and disqualifying
    on an engine that allows one writer per file: a render would have locked out
    every other process, including read-only ones. A reader also has no business
    being able to write to the collection's index.
    """
    registry_path = (
        path
        if path is not None
        else snapshots.require_latest_snapshot(
            settings.registry_dir,
            backend=settings.registry.backend,
        )
    )
    if registry_path.exists():
        with connect_registry(registry_path, read_only=True) as connection:
            current_version = connection.schema_version()
            if current_version == SCHEMA_VERSION:
                _verify_meta(connection, settings)
                yield connection
                return
            if current_version != 0:
                raise RegistrySchemaError(current_version, path=registry_path)
    # Either an explicit path that does not exist yet, or a file carrying no
    # schema at all. Both are a caller asking for a registry rather than reading
    # one, so this falls back to creating it.
    with connect_registry(registry_path) as connection:
        initialize_registry(connection, settings)
        yield connection


@contextmanager
def open_or_create_registry(settings: RegistrySettings) -> Generator[RegistryConnection]:
    """Open the newest registry, creating a dated one when the directory has none.

    For write paths. Indexing a freshly built database, or inspecting an upload
    against whatever is registered, must not fail merely because no sweep has run
    yet. Read paths use :func:`open_registry`, which reports an absent registry
    rather than inventing an empty one and calling it "no databases".
    """
    suffix = factory.suffix_for(settings.registry.backend)
    path = snapshots.latest_snapshot(
        settings.registry_dir, suffix=suffix
    ) or snapshots.new_snapshot_path(settings.registry_dir, suffix=suffix)
    path.parent.mkdir(parents=True, exist_ok=True)
    with connect_registry(path) as connection:
        initialize_registry(connection, settings)
        yield connection


def _record_id_fingerprint(
    connection: RegistryConnection,
    table: str,
    database_id: int,
    kind: EntryKind,
) -> str:
    ids = (
        str(row[0])
        for row in connection.execute(
            f"SELECT DISTINCT sequence_id FROM {table} "
            f"WHERE database_id = ? AND kind = ? ORDER BY sequence_id{connection.binary_collation()}",
            (database_id, kind.value),
        )
    )
    return id_set_fingerprint(ids)


def _record_content_fingerprint(
    connection: RegistryConnection,
    table: str,
    database_id: int,
    kind: EntryKind,
) -> str:
    rows = connection.execute(
        f"SELECT DISTINCT sequence_id, sequence_hash FROM {table} "
        "WHERE database_id = ? AND kind = ? "
        f"ORDER BY sequence_id{connection.binary_collation()}, sequence_hash",
        (database_id, kind.value),
    )
    return content_fingerprint((str(row[0]), bytes(row[1])) for row in rows)


def _record_description_fingerprint(
    connection: RegistryConnection,
    table: str,
    database_id: int,
    kind: EntryKind,
) -> str:
    rows = connection.execute(
        f"SELECT DISTINCT description_hash FROM {table} "
        "WHERE database_id = ? AND kind = ? AND description_hash IS NOT NULL "
        "ORDER BY description_hash",
        (database_id, kind.value),
    )
    return description_set_fingerprint(bytes(row[0]) for row in rows)


def _relative_path(path: Path, root: Path | None) -> str:
    if root is None:
        return path.name
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.name


def _collect_sentinel(
    raw_header: str,
    facts: _ScanFacts,
    settings: RegistrySettings,
    decoy_prefix: str,
) -> None:
    core, marker_id, _ = sentinel_details(raw_header, decoy_prefix)
    if marker_id is not None:
        if marker_id not in facts.contaminant_markers:
            facts.contaminant_markers.append(marker_id)
        return
    if facts.sentinel_header is None:
        facts.sentinel_header = core
        parsed = parse_database_metadata(core, settings.sentinel)
        facts.annotation = parsed.description if parsed and parsed.description else core


def _stream_into_table(  # noqa: C901
    connection: RegistryConnection,
    table: str,
    database_id: int,
    paths: Iterable[Path],
    settings: RegistrySettings,
    *,
    kind_override: EntryKind | None = None,
    contaminant_group: str | None = None,
    progress_label: str | None = None,
    on_progress: Callable[[int], None] | None = None,
    max_detailed_entries: int | None = None,
    strict: bool = False,
) -> _ScanFacts:
    """Scan FASTAs into bounded details, exact lengths, and sampled composition.

    Every sequence is normalized first (upper-cased, one trailing stop removed), so
    what is stored and compared is comparable content rather than whatever a source
    happened to write. What changed is counted into the facts.

    ``strict`` decides what an invalid entry does. Ingest -- an upload, a download,
    a build input -- raises, so a bad file cannot become a database. A sweep over
    the installed collection records the same findings and continues, because a file
    that is already there must stay indexable.

    ``progress_label`` narrates a bulk sweep into the log; ``on_progress`` receives
    the entries read so far, for a caller narrating one scan to a browser.
    """
    diagnostics = load_registry_diagnostics(settings.registry_diagnostics_path)
    facts = _ScanFacts()
    aa_sampler = (
        _AminoAcidReservoir(settings.metadata_aa_sample_size)
        if max_detailed_entries is not None
        else None
    )
    rows: list[tuple[int, int, str, str, str | None, int, bytes, bytes | None]] = []
    ordinal = 0
    started = time.perf_counter()
    last_progress = started
    last_report = started
    batch_size = connection.entry_batch_size
    for path in paths:
        block_state: ContaminantBlockState | None = None
        stat = path.stat()
        facts.filenames.append(path.name)
        facts.file_size_bytes += stat.st_size
        facts.mtime_ns = max(facts.mtime_ns, stat.st_mtime_ns)
        for record in iter_protein_diagnostics(path, diagnostics.rules):
            sequence_id = record.protein.id
            if not sequence_id:
                raise FastaReadError(
                    str(path),
                    f"FASTA entry ordinal {ordinal} has an empty sequence ID",
                )
            sequence = record.protein.sequence
            facts.upper_cased_entries += record.upper_cased
            facts.stop_stripped_entries += record.stop_stripped
            _count_namespace(
                facts.id_namespaces,
                record.identifier_namespace,
                diagnostics.max_reported_id_namespaces,
            )
            if record.protein.description is None:
                facts.bare_identifier_entries += 1
            if not sequence:
                facts.empty_sequence_entries += 1
                if strict:
                    raise FastaValidationError(
                        path, sequence_id, "FASTA entry has a header but no sequence"
                    )
            illegal = record.illegal_residues
            if illegal:
                facts.illegal_residue_entries += 1
                facts.illegal_residues.update(illegal)
                if strict:
                    raise FastaValidationError(
                        path, sequence_id, describe_illegal_residues(illegal)
                    )
            classified_kind, entry_group, block_state = classify_record(
                record.raw_header,
                record.classifications,
                block_state,
                diagnostics.decoy_prefix,
            )
            kind = classified_kind
            if kind_override is not None and classified_kind in {
                EntryKind.TARGET,
                EntryKind.CONTAMINANT,
            }:
                kind = kind_override
                entry_group = (
                    normalize_contaminant_group(contaminant_group)
                    if kind is EntryKind.CONTAMINANT
                    else None
                )
            if (
                max_detailed_entries is not None
                and ordinal >= max_detailed_entries
                and facts.detail_level is DetailLevel.FULL
            ):
                facts.detail_level = DetailLevel.METADATA_ONLY
                rows.clear()
                if progress_label is not None:
                    logger.info(
                        "bulk load: file {} exceeded {} entries; continuing metadata-only",
                        progress_label,
                        f"{max_detailed_entries:,}",
                    )
            facts.counts[kind] += 1
            if facts.detail_level is DetailLevel.FULL:
                facts.kind_summaries[kind].add(sequence)
            else:
                facts.kind_summaries[kind].add_length(len(sequence))
            if aa_sampler is not None:
                aa_sampler.add(kind, sequence)
            if classified_kind is EntryKind.SENTINEL:
                _collect_sentinel(record.raw_header, facts, settings, diagnostics.decoy_prefix)
            if facts.detail_level is DetailLevel.FULL:
                rows.append(
                    (
                        database_id,
                        ordinal,
                        sequence_id,
                        kind.value,
                        entry_group,
                        len(sequence),
                        sequence_hash(sequence),
                        normalized_description_hash(record.raw_header),
                    )
                )
            ordinal += 1
            if max_detailed_entries is None and len(rows) >= batch_size:
                connection.insert_entries(table, rows)
                rows.clear()
            if (progress_label is not None or on_progress is not None) and (
                ordinal % _ENTRY_PROGRESS_CHECK_INTERVAL == 0
            ):
                now = time.perf_counter()
                if (
                    progress_label is not None
                    and now - last_progress >= _ENTRY_PROGRESS_LOG_SECONDS
                ):
                    logger.info(
                        "bulk load: file {} reading {}: {} entries ({:.1f}s)",
                        progress_label,
                        path.name,
                        f"{ordinal:,}",
                        now - started,
                    )
                    last_progress = now
                if on_progress is not None and now - last_report >= _ENTRY_PROGRESS_REPORT_SECONDS:
                    on_progress(ordinal)
                    last_report = now
    if facts.detail_level is DetailLevel.METADATA_ONLY:
        assert aa_sampler is not None
        facts.sampled_aa_counts, facts.sampled_aa_sizes = aa_sampler.summarize()
    if facts.detail_level is DetailLevel.FULL:
        for offset in range(0, len(rows), batch_size):
            connection.insert_entries(table, rows[offset : offset + batch_size])
    return facts


def _aggregate_kind_stats(
    connection: RegistryConnection,
    *,
    table: str,
    database_id: int,
    kind: EntryKind,
    accumulator: SummaryAccumulator,
    detail_level: DetailLevel,
    sampled_aa_counts: Counter[str] | None,
    sampled_aa_size: int | None,
) -> DatabaseKindStats:
    """Aggregate one entry kind without materializing IDs or checksums in Python."""
    summary = accumulator.summary()
    if detail_level is DetailLevel.METADATA_ONLY:
        return DatabaseKindStats(
            kind=kind,
            entry_count=summary.n_sequences,
            distinct_ids=None,
            distinct_sequences=None,
            distinct_descriptions=None,
            distinct_pairs=None,
            duplicate_id_occurrences=None,
            conflicting_ids=None,
            repeated_sequences=None,
            length_min=summary.length_min,
            length_q1=summary.length_q1,
            length_median=summary.length_median,
            length_mean=summary.length_mean,
            length_q3=summary.length_q3,
            length_max=summary.length_max,
            total_residues=summary.total_residues,
            aa_sample_size=sampled_aa_size or 0,
            aa_counts=dict(sorted((sampled_aa_counts or {}).items())),
            id_fingerprint=None,
            description_fingerprint=None,
            content_fingerprint=None,
        )
    params: tuple[object, ...] = (database_id, kind.value)
    distinct_ids = connection.scalar(
        f"SELECT COUNT(DISTINCT sequence_id) FROM {table} WHERE database_id = ? AND kind = ?",
        params,
    )
    distinct_sequences = connection.scalar(
        f"SELECT COUNT(DISTINCT sequence_hash) FROM {table} WHERE database_id = ? AND kind = ?",
        params,
    )
    distinct_descriptions = connection.scalar(
        f"SELECT COUNT(DISTINCT description_hash) FROM {table} "
        "WHERE database_id = ? AND kind = ? AND description_hash IS NOT NULL",
        params,
    )
    distinct_pairs = connection.scalar(
        "SELECT COUNT(*) FROM ("
        f"SELECT DISTINCT sequence_id, sequence_hash FROM {table} "
        "WHERE database_id = ? AND kind = ?)",
        params,
    )
    conflicting_ids = connection.scalar(
        "SELECT COUNT(*) FROM ("
        f"SELECT sequence_id FROM {table} WHERE database_id = ? AND kind = ? "
        "GROUP BY sequence_id HAVING COUNT(DISTINCT sequence_hash) > 1)",
        params,
    )
    repeated_sequences = connection.scalar(
        "SELECT COUNT(*) FROM ("
        f"SELECT sequence_hash FROM {table} WHERE database_id = ? AND kind = ? "
        "GROUP BY sequence_hash HAVING COUNT(DISTINCT sequence_id) > 1)",
        params,
    )
    return DatabaseKindStats(
        kind=kind,
        entry_count=summary.n_sequences,
        distinct_ids=distinct_ids,
        distinct_sequences=distinct_sequences,
        distinct_descriptions=distinct_descriptions,
        distinct_pairs=distinct_pairs,
        duplicate_id_occurrences=max(0, summary.n_sequences - distinct_ids),
        conflicting_ids=conflicting_ids,
        repeated_sequences=repeated_sequences,
        length_min=summary.length_min,
        length_q1=summary.length_q1,
        length_median=summary.length_median,
        length_mean=summary.length_mean,
        length_q3=summary.length_q3,
        length_max=summary.length_max,
        total_residues=summary.total_residues,
        aa_sample_size=summary.n_sequences,
        aa_counts=summary.aa_frequencies,
        id_fingerprint=_record_id_fingerprint(connection, table, database_id, kind),
        description_fingerprint=_record_description_fingerprint(
            connection, table, database_id, kind
        ),
        content_fingerprint=_record_content_fingerprint(connection, table, database_id, kind),
    )


def _aggregate_record(
    connection: RegistryConnection,
    *,
    table: str,
    database_id: int,
    relative_path: str,
    filename: str,
    facts: _ScanFacts,
    parsed_decoy: bool,
) -> RegistryRecord:
    kind_stats = {
        kind: _aggregate_kind_stats(
            connection,
            table=table,
            database_id=database_id,
            kind=kind,
            accumulator=facts.kind_summaries[kind],
            detail_level=facts.detail_level,
            sampled_aa_counts=(
                facts.sampled_aa_counts[kind] if facts.sampled_aa_counts is not None else None
            ),
            sampled_aa_size=(
                facts.sampled_aa_sizes[kind] if facts.sampled_aa_sizes is not None else None
            ),
        )
        for kind in EntryKind
    }
    summary_accumulator = SummaryAccumulator()
    for accumulator in facts.kind_summaries.values():
        summary_accumulator.merge(accumulator)
    summary = summary_accumulator.summary()
    aa_counts = summary.aa_frequencies
    if facts.sampled_aa_counts is not None:
        combined_sample = Counter[str]()
        for counts in facts.sampled_aa_counts.values():
            combined_sample.update(counts)
        aa_counts = dict(sorted(combined_sample.items()))
    aa_sample_size = (
        sum(facts.sampled_aa_sizes.values())
        if facts.sampled_aa_sizes is not None
        else summary.n_sequences
    )
    target_stats = kind_stats[EntryKind.TARGET]
    return RegistryRecord(
        id=database_id,
        relative_path=relative_path,
        filename=filename,
        file_size_bytes=facts.file_size_bytes,
        mtime_ns=facts.mtime_ns,
        sentinel_header=facts.sentinel_header,
        annotation=facts.annotation,
        filename_is_decoy=parsed_decoy,
        is_decoy=bool(facts.counts[EntryKind.DECOY]) or parsed_decoy,
        contaminant_markers=facts.contaminant_markers,
        indexed_at=datetime.datetime.now(datetime.UTC).isoformat(),
        detail_level=facts.detail_level,
        entry_count=summary.n_sequences,
        target_count=target_stats.entry_count,
        decoy_count=facts.counts[EntryKind.DECOY],
        contaminant_count=facts.counts[EntryKind.CONTAMINANT],
        entrapment_count=facts.counts[EntryKind.ENTRAPMENT],
        sentinel_count=facts.counts[EntryKind.SENTINEL],
        distinct_target_ids=target_stats.distinct_ids,
        distinct_target_sequences=target_stats.distinct_sequences,
        distinct_target_descriptions=target_stats.distinct_descriptions,
        duplicate_target_id_occurrences=target_stats.duplicate_id_occurrences,
        conflicting_target_ids=target_stats.conflicting_ids,
        repeated_target_sequences=target_stats.repeated_sequences,
        length_min=summary.length_min,
        length_q1=summary.length_q1,
        length_median=summary.length_median,
        length_mean=summary.length_mean,
        length_q3=summary.length_q3,
        length_max=summary.length_max,
        total_residues=summary.total_residues,
        aa_sample_size=aa_sample_size,
        aa_counts=aa_counts,
        target_id_fingerprint=target_stats.id_fingerprint,
        target_description_fingerprint=target_stats.description_fingerprint,
        target_content_fingerprint=target_stats.content_fingerprint,
        upper_cased_entries=facts.upper_cased_entries,
        stop_stripped_entries=facts.stop_stripped_entries,
        illegal_residue_entries=facts.illegal_residue_entries,
        illegal_residues=dict(sorted(facts.illegal_residues.items())),
        empty_sequence_entries=facts.empty_sequence_entries,
        bare_identifier_entries=facts.bare_identifier_entries,
        id_namespaces=dict(sorted(facts.id_namespaces.items())),
        kind_stats=kind_stats,
    )


_DATABASE_COLUMNS = (
    "relative_path, filename, dbname, file_size_bytes, mtime_ns, sentinel_header, annotation, "
    "filename_is_decoy, is_decoy, contaminant_markers_json, indexed_at, detail_level, entry_count, target_count, "
    "decoy_count, contaminant_count, entrapment_count, sentinel_count, distinct_target_ids, distinct_target_sequences, "
    "distinct_target_descriptions, "
    "duplicate_target_id_occurrences, conflicting_target_ids, repeated_target_sequences, length_min, length_q1, "
    "length_median, length_mean, length_q3, length_max, total_residues, aa_sample_size, aa_counts_json, "
    "target_id_fingerprint, target_description_fingerprint, target_content_fingerprint, "
    "upper_cased_entries, stop_stripped_entries, illegal_residue_entries, illegal_residues_json, "
    "empty_sequence_entries, bare_identifier_entries, id_namespaces_json"
)
_DATABASE_COLUMN_NAMES = tuple(column.strip() for column in _DATABASE_COLUMNS.split(","))
_KIND_STATS_COLUMNS = (
    "database_id, kind, entry_count, distinct_ids, distinct_sequences, distinct_descriptions, distinct_pairs, "
    "duplicate_id_occurrences, "
    "conflicting_ids, repeated_sequences, length_min, length_q1, length_median, length_mean, length_q3, "
    "length_max, total_residues, aa_sample_size, aa_counts_json, id_fingerprint, description_fingerprint, "
    "content_fingerprint"
)


def _record_values(record: RegistryRecord) -> tuple[object, ...]:
    return (
        record.relative_path,
        record.filename,
        record.dbname,
        record.file_size_bytes,
        record.mtime_ns,
        record.sentinel_header,
        record.annotation,
        int(record.filename_is_decoy),
        int(record.is_decoy),
        json.dumps(record.contaminant_markers, separators=(",", ":")),
        record.indexed_at,
        record.detail_level.value,
        record.entry_count,
        record.target_count,
        record.decoy_count,
        record.contaminant_count,
        record.entrapment_count,
        record.sentinel_count,
        record.distinct_target_ids,
        record.distinct_target_sequences,
        record.distinct_target_descriptions,
        record.duplicate_target_id_occurrences,
        record.conflicting_target_ids,
        record.repeated_target_sequences,
        record.length_min,
        record.length_q1,
        record.length_median,
        record.length_mean,
        record.length_q3,
        record.length_max,
        record.total_residues,
        record.aa_sample_size,
        json.dumps(record.aa_counts, sort_keys=True, separators=(",", ":")),
        record.target_id_fingerprint,
        record.target_description_fingerprint,
        record.target_content_fingerprint,
        record.upper_cased_entries,
        record.stop_stripped_entries,
        record.illegal_residue_entries,
        json.dumps(record.illegal_residues, sort_keys=True, separators=(",", ":")),
        record.empty_sequence_entries,
        record.bare_identifier_entries,
        json.dumps(record.id_namespaces, sort_keys=True, separators=(",", ":")),
    )


def _validated_record_values(record: RegistryRecord) -> tuple[object, ...]:
    """Return persisted values after checking the column/value contract."""
    values = _record_values(record)
    if len(values) != len(_DATABASE_COLUMN_NAMES):
        raise RuntimeError(
            f"Registry database column/value mismatch: {len(_DATABASE_COLUMN_NAMES)} columns but {len(values)} values."
        )
    return values


def _kind_stats_values(database_id: int, stats: DatabaseKindStats) -> tuple[object, ...]:
    return (
        database_id,
        stats.kind.value,
        stats.entry_count,
        stats.distinct_ids,
        stats.distinct_sequences,
        stats.distinct_descriptions,
        stats.distinct_pairs,
        stats.duplicate_id_occurrences,
        stats.conflicting_ids,
        stats.repeated_sequences,
        stats.length_min,
        stats.length_q1,
        stats.length_median,
        stats.length_mean,
        stats.length_q3,
        stats.length_max,
        stats.total_residues,
        stats.aa_sample_size,
        json.dumps(stats.aa_counts, sort_keys=True, separators=(",", ":")),
        stats.id_fingerprint,
        stats.description_fingerprint,
        stats.content_fingerprint,
    )


def _replace_database_kind_stats(
    connection: RegistryConnection,
    database_id: int,
    kind_stats: dict[EntryKind, DatabaseKindStats],
) -> None:
    """Persist exactly one normalized summary row for every entry kind."""
    missing = set(EntryKind) - set(kind_stats)
    if missing:
        raise RuntimeError(
            f"Missing database kind statistics for: {', '.join(sorted(kind.value for kind in missing))}"
        )
    connection.execute("DELETE FROM database_kind_stats WHERE database_id = ?", (database_id,))
    connection.executemany(
        f"INSERT INTO database_kind_stats ({_KIND_STATS_COLUMNS}) "
        f"VALUES ({','.join('?' for _ in _KIND_STATS_COLUMNS.split(','))})",
        (_kind_stats_values(database_id, kind_stats[kind]) for kind in EntryKind),
    )


def _refresh_database_pair_stats(
    connection: RegistryConnection,
    database_id: int,
) -> None:
    """Replace materialized comparisons involving one registered database."""
    connection.execute(
        "DELETE FROM database_pair_stats WHERE database_id_low = ? OR database_id_high = ?",
        (database_id, database_id),
    )
    detail_row = connection.execute(
        "SELECT detail_level FROM databases WHERE id = ?",
        (database_id,),
    ).fetchone()
    if detail_row is None:
        raise ValueError(f"database id {database_id} is not registered")
    if DetailLevel(str(detail_row[0])) is DetailLevel.METADATA_ONLY:
        return
    other_database_ids = [
        int(row[0])
        for row in connection.execute(
            "SELECT id FROM databases WHERE id != ? AND detail_level = ? ORDER BY id",
            (database_id, DetailLevel.FULL.value),
        )
    ]
    rows: list[tuple[int, int, str, int, int, int, int, int]] = []
    selection = PairMetricSelection(
        ids_table="entries",
        sequences_table="entries",
        descriptions_table="entries",
        pairs_table="entries",
        where="database_id = ?",
        params=(database_id,),
    )
    for kind in _MATERIALIZED_PAIR_KINDS:
        metrics = pair_metric_counts(
            connection,
            selection=selection,
            kind=kind,
            excluded_database_id=database_id,
        )
        rows.extend(
            (
                min(database_id, other_database_id),
                max(database_id, other_database_id),
                kind.value,
                metrics.shared_ids.get(other_database_id, 0),
                metrics.shared_sequences.get(other_database_id, 0),
                metrics.shared_descriptions.get(other_database_id, 0),
                metrics.shared_pairs.get(other_database_id, 0),
                metrics.matching_ids.get(other_database_id, 0),
            )
            for other_database_id in other_database_ids
        )
    connection.executemany(
        """
        INSERT INTO database_pair_stats (
            database_id_low, database_id_high, kind, shared_ids,
            shared_sequence_checksums, shared_descriptions, shared_exact_pairs,
            matching_shared_ids
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def _kind_content_groups(
    connection: RegistryConnection, kind: EntryKind
) -> list[_KindContentGroup]:
    """Group databases whose complete distinct content is interchangeable."""
    grouped: dict[tuple[str, str, int, int, int, int], list[int]] = {}
    for row in connection.execute(
        """
        SELECT database_id, content_fingerprint, description_fingerprint,
               distinct_ids, distinct_sequences, distinct_descriptions,
               distinct_pairs
        FROM database_kind_stats AS stats
        JOIN databases ON databases.id = stats.database_id
        WHERE stats.kind = ? AND databases.detail_level = ?
        ORDER BY stats.database_id
        """,
        (kind.value, DetailLevel.FULL.value),
    ):
        key = (
            str(row[1]),
            str(row[2]),
            int(row[3]),
            int(row[4]),
            int(row[5]),
            int(row[6]),
        )
        grouped.setdefault(key, []).append(int(row[0]))
    return [
        _KindContentGroup(
            representative_id=database_ids[0],
            database_ids=tuple(database_ids),
            distinct_ids=key[2],
            distinct_sequences=key[3],
            distinct_descriptions=key[4],
            distinct_pairs=key[5],
        )
        for key, database_ids in grouped.items()
    ]


def _representative_pair_metrics(
    connection: RegistryConnection,
    kind: EntryKind,
    groups: list[_KindContentGroup],
) -> dict[tuple[int, int], tuple[int, int, int, int, int]]:
    """Compute one overlap for each pair of distinct content groups."""
    representatives = [group.representative_id for group in groups]
    connection.create_temp_table(
        TempTableSpec(
            name=_PAIR_REPRESENTATIVES_TABLE,
            columns=(("database_id", "INTEGER NOT NULL"),),
            primary_key=("database_id",),
        )
    )
    connection.executemany(
        f"INSERT INTO {connection.temp(_PAIR_REPRESENTATIVES_TABLE)} (database_id) VALUES (?)",
        ((database_id,) for database_id in representatives),
    )

    total_group_pairs = len(groups) * (len(groups) - 1) // 2
    completed_group_pairs = 0
    progress_interval = max(1, len(groups) // 20)
    started = time.perf_counter()
    result: dict[tuple[int, int], tuple[int, int, int, int, int]] = {}
    for index, group in enumerate(groups[:-1]):
        representative_id = group.representative_id
        selection = PairMetricSelection(
            ids_table="entries",
            sequences_table="entries",
            descriptions_table="entries",
            pairs_table="entries",
            where="database_id = ?",
            params=(representative_id,),
        )
        metrics = pair_metric_counts(
            connection,
            selection=selection,
            kind=kind,
            excluded_database_id=None,
            minimum_other_database_id=representative_id,
            other_database_table=connection.temp(_PAIR_REPRESENTATIVES_TABLE),
        )
        for other_group in groups[index + 1 :]:
            other_id = other_group.representative_id
            result[(representative_id, other_id)] = (
                metrics.shared_ids.get(other_id, 0),
                metrics.shared_sequences.get(other_id, 0),
                metrics.shared_descriptions.get(other_id, 0),
                metrics.shared_pairs.get(other_id, 0),
                metrics.matching_ids.get(other_id, 0),
            )
        completed_group_pairs += len(groups) - index - 1
        if index % progress_interval == 0 or completed_group_pairs == total_group_pairs:
            logger.info(
                "pair statistics ({}): {}/{} unique-content comparisons ({:.1f}%, {:.1f}s)",
                kind.value,
                f"{completed_group_pairs:,}",
                f"{total_group_pairs:,}",
                completed_group_pairs / total_group_pairs * 100 if total_group_pairs else 100.0,
                time.perf_counter() - started,
            )
    connection.execute(f"DROP TABLE {connection.temp(_PAIR_REPRESENTATIVES_TABLE)}")
    return result


def _materialize_kind_pair_stats(
    connection: RegistryConnection,
    kind: EntryKind,
    database_ids: list[int],
) -> int:
    """Write all database pairs for one kind, reusing identical-content results."""
    groups = _kind_content_groups(connection, kind)
    group_by_database = {
        database_id: group for group in groups for database_id in group.database_ids
    }
    representative_metrics = _representative_pair_metrics(connection, kind, groups)
    logger.info(
        "pair statistics ({}): {} databases collapse to {} unique content set(s)",
        kind.value,
        f"{len(database_ids):,}",
        f"{len(groups):,}",
    )

    insert_sql = """
        INSERT INTO database_pair_stats (
            database_id_low, database_id_high, kind, shared_ids,
            shared_sequence_checksums, shared_descriptions, shared_exact_pairs,
            matching_shared_ids
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """
    rows: list[tuple[int, int, str, int, int, int, int, int]] = []
    written = 0
    for index, database_id_low in enumerate(database_ids[:-1]):
        low_group = group_by_database[database_id_low]
        for database_id_high in database_ids[index + 1 :]:
            high_group = group_by_database[database_id_high]
            if low_group is high_group:
                values = (
                    low_group.distinct_ids,
                    low_group.distinct_sequences,
                    low_group.distinct_descriptions,
                    low_group.distinct_pairs,
                    low_group.distinct_ids,
                )
            else:
                representative_pair = (
                    min(low_group.representative_id, high_group.representative_id),
                    max(low_group.representative_id, high_group.representative_id),
                )
                values = representative_metrics[representative_pair]
            rows.append((database_id_low, database_id_high, kind.value, *values))
            if len(rows) >= _PAIR_INSERT_BATCH_SIZE:
                connection.executemany(insert_sql, rows)
                written += len(rows)
                rows.clear()
    if rows:
        connection.executemany(insert_sql, rows)
        written += len(rows)
    return written


def _materialize_all_pair_stats(connection: RegistryConnection) -> int:
    """Materialize every unordered database pair exactly once per comparable kind."""
    database_ids = [
        int(row[0])
        for row in connection.execute(
            "SELECT id FROM databases WHERE detail_level = ? ORDER BY id",
            (DetailLevel.FULL.value,),
        )
    ]
    expected_rows = len(database_ids) * (len(database_ids) - 1) // 2 * len(_MATERIALIZED_PAIR_KINDS)
    started = time.perf_counter()
    with connection.transaction():
        connection.execute("DELETE FROM database_pair_stats")
        written = sum(
            _materialize_kind_pair_stats(connection, kind, database_ids)
            for kind in _MATERIALIZED_PAIR_KINDS
        )
    logger.info(
        "pair statistics: wrote {}/{} materialized rows ({:.1f}s)",
        f"{written:,}",
        f"{expected_rows:,}",
        time.perf_counter() - started,
    )
    return written


def _index_fasta(
    connection: RegistryConnection,
    path: Path,
    settings: RegistrySettings,
    *,
    root: Path | None = None,
    refresh_pair_stats: bool,
    log_record: bool,
    progress_label: str | None = None,
) -> RegistryRecord:
    """Index one FASTA, optionally refreshing its materialized pair rows."""
    started = time.perf_counter()
    relative_path = _relative_path(path, root)
    parsed = parse_filename(path.name, settings.naming)
    if progress_label is not None:
        logger.info(
            "bulk load: file {} started: {}",
            progress_label,
            path.name,
        )
    with connection.transaction():
        existing = connection.execute(
            "SELECT id FROM databases WHERE relative_path = ?",
            (relative_path,),
        ).fetchone()
        if existing is None:
            initial_values = _validated_record_values(
                RegistryRecord(filename=path.name, relative_path=relative_path)
            )
            # How an identifier is assigned is the backend's business: a rowid
            # alias on one engine, an explicit sequence on another. Either way it
            # comes back from the statement that assigned it.
            database_id = connection.insert_database(_DATABASE_COLUMNS, initial_values)
        else:
            database_id = int(existing[0])
            connection.execute("DELETE FROM entries WHERE database_id = ?", (database_id,))
        facts = _stream_into_table(
            connection,
            "entries",
            database_id,
            [path],
            settings,
            progress_label=progress_label,
            max_detailed_entries=settings.max_detailed_entries,
        )
        if progress_label is not None:
            logger.info(
                "bulk load: file {} scan complete: {} entries, {:.2f} MB; aggregating statistics ({:.1f}s)",
                progress_label,
                f"{sum(facts.counts.values()):,}",
                facts.file_size_bytes / 1_000_000,
                time.perf_counter() - started,
            )
        record = replace(
            _aggregate_record(
                connection,
                table="entries",
                database_id=database_id,
                relative_path=relative_path,
                filename=path.name,
                facts=facts,
                parsed_decoy=parsed.is_decoy,
            ),
            dbname=parsed.dbname,
        )
        record_values = _validated_record_values(record)
        assignments = ", ".join(f"{column} = ?" for column in _DATABASE_COLUMN_NAMES)
        connection.execute(
            f"UPDATE databases SET {assignments} WHERE id = ?",
            (*record_values, database_id),
        )
        _replace_database_kind_stats(connection, database_id, record.kind_stats)
        if refresh_pair_stats:
            _refresh_database_pair_stats(connection, database_id)
    if log_record:
        logger.info(
            "indexed {} ({} entries, {})",
            path.name,
            record.entry_count,
            record.detail_level.value,
        )
    if progress_label is not None:
        logger.info(
            "bulk load: file {} complete: {} ({} entries, {}, {:.1f}s)",
            progress_label,
            path.name,
            f"{record.entry_count:,}",
            record.detail_level.value,
            time.perf_counter() - started,
        )
    return record


def index_fasta(
    connection: RegistryConnection,
    path: Path,
    settings: RegistrySettings,
    *,
    root: Path | None = None,
) -> RegistryRecord:
    """Atomically index one FASTA and refresh all pair rows involving it."""
    return _index_fasta(
        connection,
        path,
        settings,
        root=root,
        refresh_pair_stats=True,
        log_record=True,
    )


def _candidate_paths(paths: Path | Iterable[Path]) -> list[Path]:
    return [paths] if isinstance(paths, Path) else list(paths)


def _candidate_filename_is_decoy(path: Path, label: str, settings: RegistrySettings) -> bool:
    """Preserve naming evidence when an upload was staged under a safe prefix."""
    return any(
        parse_filename(filename, settings.naming).is_decoy
        for filename in (path.name, Path(label).name)
    )


def _create_candidate_table(connection: RegistryConnection) -> None:
    """Replace the connection-local candidate table with an empty table."""
    connection.create_temp_table(
        TempTableSpec(
            name=CANDIDATE_TABLE,
            columns=(
                ("database_id", "INTEGER NOT NULL"),
                ("ordinal", "INTEGER NOT NULL"),
                ("sequence_id", "TEXT NOT NULL"),
                ("kind", "TEXT NOT NULL"),
                ("contaminant_group", "TEXT"),
                ("sequence_length", "INTEGER NOT NULL"),
                ("sequence_hash", "BLOB NOT NULL"),
                ("description_hash", "BLOB"),
            ),
            primary_key=("database_id", "ordinal"),
            indexes=(
                ("kind", "sequence_id"),
                ("kind", "sequence_hash"),
                ("kind", "description_hash"),
            ),
            partial_index_column="description_hash",
        )
    )


def _report_entries_read(on_progress: Callable[[str], None], label: str, entries: int) -> None:
    """Narrate an in-flight scan of one labelled candidate file."""
    on_progress(f"Indexing {label}: {entries:,} entries read")


def populate_candidate_files(
    connection: RegistryConnection,
    paths: Path | Iterable[Path],
    settings: RegistrySettings,
    *,
    labels: Iterable[str] | None = None,
    combined_label: str | None = None,
    kind_override: EntryKind | None = None,
    contaminant_groups: Iterable[str] | None = None,
    on_progress: Callable[[str], None] | None = None,
    strict: bool = True,
) -> tuple[list[RegistryRecord], RegistryRecord]:
    """Populate per-file and combined candidate aggregates in one FASTA pass.

    Rows for each source file use temporary database IDs starting at one. The
    combined selection is copied inside SQLite to database ID zero, which is
    the scope consumed by candidate overlap comparisons.

    ``on_progress`` receives one line per indexing phase, reported against the
    caller's labels rather than the paths, which for a download are cache tokens.

    ``strict`` defaults to refusing an invalid entry, because a candidate is
    something on its way in. Pass ``strict=False`` to inspect a file that is
    already installed, which must stay inspectable whatever it contains.
    """
    selected_paths = _candidate_paths(paths)
    if not selected_paths:
        raise ValueError("Select at least one FASTA file to inspect.")
    selected_labels = list(labels) if labels is not None else [path.name for path in selected_paths]
    if len(selected_labels) != len(selected_paths):
        raise ValueError("Candidate labels must match the number of FASTA files.")
    selected_groups: list[str | None]
    if contaminant_groups is None:
        selected_groups = [None] * len(selected_paths)
    else:
        if kind_override is not EntryKind.CONTAMINANT:
            raise ValueError("Candidate contaminant groups require a contaminant kind override.")
        selected_groups = list(contaminant_groups)
        if len(selected_groups) != len(selected_paths):
            raise ValueError("Candidate contaminant groups must match the number of FASTA files.")

    _create_candidate_table(connection)
    combined_facts = _ScanFacts()
    per_file: list[RegistryRecord] = []
    for database_id, (path, filename, contaminant_group) in enumerate(
        zip(selected_paths, selected_labels, selected_groups, strict=True),
        start=1,
    ):
        if on_progress is not None:
            on_progress(f"Indexing {filename} ...")
        facts = _stream_into_table(
            connection,
            f"{connection.temp(CANDIDATE_TABLE)}",
            database_id,
            [path],
            settings,
            kind_override=kind_override,
            contaminant_group=contaminant_group,
            on_progress=None
            if on_progress is None
            else partial(_report_entries_read, on_progress, filename),
            strict=strict,
        )
        combined_facts.merge(facts)
        if on_progress is not None:
            on_progress(f"Indexed {filename}: {sum(facts.counts.values()):,} entries")
        record = _aggregate_record(
            connection,
            table=f"{connection.temp(CANDIDATE_TABLE)}",
            database_id=database_id,
            relative_path=filename,
            filename=filename,
            facts=facts,
            parsed_decoy=_candidate_filename_is_decoy(path, filename, settings),
        )
        per_file.append(replace(record, dbname=Path(filename).stem))

    if on_progress is not None and len(selected_paths) > 1:
        on_progress(f"Merging {len(selected_paths)} sources into the combined selection ...")
    connection.execute(
        f"""
        INSERT INTO {connection.temp(CANDIDATE_TABLE)} (
            database_id, ordinal, sequence_id, kind, contaminant_group,
            sequence_length, sequence_hash, description_hash
        )
        SELECT 0,
               ROW_NUMBER() OVER (ORDER BY database_id, ordinal) - 1,
               sequence_id,
               kind,
               contaminant_group,
               sequence_length,
               sequence_hash,
               description_hash
        FROM {connection.temp(CANDIDATE_TABLE)}
        WHERE database_id > 0
        ORDER BY database_id, ordinal
        """
    )
    filename = combined_label or (
        selected_labels[0]
        if len(selected_labels) == 1
        else f"{len(selected_labels)} selected FASTAs"
    )
    combined = _aggregate_record(
        connection,
        table=f"{connection.temp(CANDIDATE_TABLE)}",
        database_id=0,
        relative_path=filename,
        filename=filename,
        facts=combined_facts,
        parsed_decoy=per_file[0].filename_is_decoy if len(per_file) == 1 else False,
    )
    return per_file, replace(combined, dbname=Path(filename).stem)


def populate_candidate(
    connection: RegistryConnection,
    paths: Path | Iterable[Path],
    settings: RegistrySettings,
    *,
    label: str | None = None,
) -> RegistryRecord:
    """Stream transient FASTA input into a connection-local comparison table."""
    selected_paths = _candidate_paths(paths)
    if not selected_paths:
        raise ValueError("Select at least one FASTA file to inspect.")
    _create_candidate_table(connection)
    facts = _stream_into_table(
        connection, f"{connection.temp(CANDIDATE_TABLE)}", 0, selected_paths, settings
    )
    filename = label or (
        selected_paths[0].name
        if len(selected_paths) == 1
        else f"{len(selected_paths)} selected FASTAs"
    )
    record = _aggregate_record(
        connection,
        table=f"{connection.temp(CANDIDATE_TABLE)}",
        database_id=0,
        relative_path=filename,
        filename=filename,
        facts=facts,
        parsed_decoy=False,
    )
    return replace(record, dbname=Path(filename).stem)


def _optional_int(value: int | str | bytes | None) -> int | None:
    """Convert one nullable SQLite integer without conflating NULL and zero."""
    return None if value is None else int(value)


def _optional_str(value: object) -> str | None:
    """Convert one nullable SQLite string without rendering NULL as text."""
    return None if value is None else str(value)


def _kind_stats_from_row(row: Row) -> DatabaseKindStats:
    return DatabaseKindStats(
        kind=EntryKind(str(row["kind"])),
        entry_count=int(row["entry_count"]),
        distinct_ids=_optional_int(row["distinct_ids"]),
        distinct_sequences=_optional_int(row["distinct_sequences"]),
        distinct_descriptions=_optional_int(row["distinct_descriptions"]),
        distinct_pairs=_optional_int(row["distinct_pairs"]),
        duplicate_id_occurrences=_optional_int(row["duplicate_id_occurrences"]),
        conflicting_ids=_optional_int(row["conflicting_ids"]),
        repeated_sequences=_optional_int(row["repeated_sequences"]),
        length_min=int(row["length_min"]),
        length_q1=float(row["length_q1"]),
        length_median=float(row["length_median"]),
        length_mean=float(row["length_mean"]),
        length_q3=float(row["length_q3"]),
        length_max=int(row["length_max"]),
        total_residues=int(row["total_residues"]),
        aa_sample_size=int(row["aa_sample_size"]),
        aa_counts=json.loads(str(row["aa_counts_json"])),
        id_fingerprint=_optional_str(row["id_fingerprint"]),
        description_fingerprint=_optional_str(row["description_fingerprint"]),
        content_fingerprint=_optional_str(row["content_fingerprint"]),
    )


def _load_kind_stats(
    connection: RegistryConnection,
    database_ids: Iterable[int],
) -> dict[int, dict[EntryKind, DatabaseKindStats]]:
    ids = list(database_ids)
    result: dict[int, dict[EntryKind, DatabaseKindStats]] = {database_id: {} for database_id in ids}
    if not ids:
        return result
    placeholders = ",".join("?" for _ in ids)
    rows = connection.execute(
        f"SELECT * FROM database_kind_stats WHERE database_id IN ({placeholders}) ORDER BY database_id, kind",
        ids,
    ).fetchall()
    for row in rows:
        database_id = int(row["database_id"])
        stats = _kind_stats_from_row(row)
        result[database_id][stats.kind] = stats
    expected_kinds = set(EntryKind)
    for database_id, database_stats in result.items():
        if set(database_stats) != expected_kinds:
            missing = ", ".join(sorted(kind.value for kind in expected_kinds - set(database_stats)))
            raise RegistryIntegrityError(
                f"Kind statistics are incomplete for database id {database_id}"
                f" (missing: {missing or 'unknown'}); run a full registry reindex."
            )
    return result


def _row_to_record(
    row: Row,
    kind_stats: dict[EntryKind, DatabaseKindStats],
) -> RegistryRecord:
    return RegistryRecord(
        id=int(row["id"]),
        relative_path=str(row["relative_path"]),
        filename=str(row["filename"]),
        dbname=str(row["dbname"]),
        file_size_bytes=int(row["file_size_bytes"]),
        mtime_ns=int(row["mtime_ns"]),
        sentinel_header=row["sentinel_header"],
        annotation=row["annotation"],
        filename_is_decoy=bool(row["filename_is_decoy"]),
        is_decoy=bool(row["is_decoy"]),
        contaminant_markers=json.loads(str(row["contaminant_markers_json"])),
        indexed_at=str(row["indexed_at"]),
        detail_level=DetailLevel(str(row["detail_level"])),
        entry_count=int(row["entry_count"]),
        target_count=int(row["target_count"]),
        decoy_count=int(row["decoy_count"]),
        contaminant_count=int(row["contaminant_count"]),
        entrapment_count=int(row["entrapment_count"]),
        sentinel_count=int(row["sentinel_count"]),
        distinct_target_ids=_optional_int(row["distinct_target_ids"]),
        distinct_target_sequences=_optional_int(row["distinct_target_sequences"]),
        distinct_target_descriptions=_optional_int(row["distinct_target_descriptions"]),
        duplicate_target_id_occurrences=_optional_int(row["duplicate_target_id_occurrences"]),
        conflicting_target_ids=_optional_int(row["conflicting_target_ids"]),
        repeated_target_sequences=_optional_int(row["repeated_target_sequences"]),
        length_min=int(row["length_min"]),
        length_q1=float(row["length_q1"]),
        length_median=float(row["length_median"]),
        length_mean=float(row["length_mean"]),
        length_q3=float(row["length_q3"]),
        length_max=int(row["length_max"]),
        total_residues=int(row["total_residues"]),
        aa_sample_size=int(row["aa_sample_size"]),
        aa_counts=json.loads(str(row["aa_counts_json"])),
        target_id_fingerprint=_optional_str(row["target_id_fingerprint"]),
        target_description_fingerprint=_optional_str(row["target_description_fingerprint"]),
        target_content_fingerprint=_optional_str(row["target_content_fingerprint"]),
        upper_cased_entries=int(row["upper_cased_entries"]),
        stop_stripped_entries=int(row["stop_stripped_entries"]),
        illegal_residue_entries=int(row["illegal_residue_entries"]),
        illegal_residues=json.loads(str(row["illegal_residues_json"])),
        empty_sequence_entries=int(row["empty_sequence_entries"]),
        bare_identifier_entries=int(row["bare_identifier_entries"]),
        id_namespaces=json.loads(str(row["id_namespaces_json"])),
        kind_stats=kind_stats,
    )


def export_registry_stats_only(source_path: Path, destination: Path) -> tuple[int, int]:
    """Copy a registry without ``entries``, returning (bytes written, pair rows).

    Everything the read paths use is materialized, so this copy serves the GUI,
    clustering, and similarity export while being a small fraction of the full
    registry. Reindexing, registry validation, and uploaded-FASTA overlap read
    ``entries`` and still need the complete file.
    """
    if not source_path.is_file():
        raise FileNotFoundError(
            f"Registry does not exist: {source_path}. Run 'fasta-gen reindex' first."
        )
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite an existing file: {destination}")

    with connect_registry(source_path, read_only=True) as source:
        version = source.schema_version()
        if version != SCHEMA_VERSION:
            raise RegistrySchemaError(version, path=source_path)

    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        pair_rows = factory.copy_stats_only(source_path, destination, SCHEMA_VERSION)
    except BaseException:
        # A partial export is worse than none: it would be indistinguishable
        # from a complete one and would block the retry on the exists check.
        destination.unlink(missing_ok=True)
        raise
    return destination.stat().st_size, pair_rows


def list_databases(connection: RegistryConnection) -> list[RegistryRecord]:
    """Return compact aggregate records for all indexed databases."""
    rows = connection.execute(
        "SELECT * FROM databases ORDER BY filename COLLATE NOCASE, id"
    ).fetchall()
    stats = _load_kind_stats(connection, (int(row["id"]) for row in rows))
    return [_row_to_record(row, stats[int(row["id"])]) for row in rows]


def get_database(connection: RegistryConnection, database_id: int) -> RegistryRecord | None:
    """Return one aggregate registry record by integer ID."""
    row = connection.execute("SELECT * FROM databases WHERE id = ?", (database_id,)).fetchone()
    if row is None:
        return None
    stats = _load_kind_stats(connection, [database_id])
    return _row_to_record(row, stats[database_id])


def iter_target_ids(path: Path, settings: RegistrySettings) -> Iterator[str]:
    """Yield target ID tokens without loading sequences."""
    diagnostics = load_registry_diagnostics(settings.registry_diagnostics_path)
    block_state: ContaminantBlockState | None = None
    for ordinal, header in enumerate(read_headers(path)):
        entry_id = parse_header(header).id
        if not entry_id:
            raise ValueError(f"{path}: FASTA entry ordinal {ordinal} has an empty sequence ID")
        _, classifications = diagnostics.rules.diagnose_identifier(entry_id)
        kind, _, block_state = classify_record(
            header,
            classifications,
            block_state,
            diagnostics.decoy_prefix,
        )
        if kind is EntryKind.TARGET:
            yield entry_id


def target_id_set(path: Path, settings: RegistrySettings) -> set[str]:
    """Return the unique target ID set for one FASTA."""
    return set(iter_target_ids(path, settings))


def iter_fasta_files(directory: Path, *, recursive: bool = False) -> Iterator[Path]:
    """Yield sorted FASTA files once across supported filename patterns."""
    matches: set[Path] = set()
    for pattern in _FASTA_GLOBS:
        discovered = directory.rglob(pattern) if recursive else directory.glob(pattern)
        matches.update(path for path in discovered if path.is_file())
    yield from sorted(matches)


def _is_oversized_fasta(path: Path, settings: RegistrySettings) -> bool:
    """Return whether a FASTA meets or exceeds the configured size limit."""
    return path.stat().st_size >= settings.max_fasta_file_size_gib * _GIBIBYTE


def _log_oversized_fasta(path: Path, settings: RegistrySettings) -> None:
    """Log one intentional pre-scan file-size exclusion."""
    logger.info(
        "skipped FASTA {}: {:.2f} GiB is at or above max_fasta_file_size_gib={:g}",
        path,
        path.stat().st_size / _GIBIBYTE,
        settings.max_fasta_file_size_gib,
    )


def _fasta_build_date(path: Path, settings: RegistrySettings) -> tuple[datetime.date, str]:
    """Return a FASTA's build date and where it came from.

    The filename's ``YYYYMMDD`` token wins because it travels with the file,
    whereas the filesystem timestamp is rewritten by copies. When the name
    carries no date, the modification time is the only evidence available; it is
    also the field the sweep already trusts for change detection.
    """
    parsed = parse_filename(path.name, settings.naming)
    if parsed.date is not None:
        return parsed.date, "filename"
    return datetime.date.fromtimestamp(path.stat().st_mtime), "mtime"


def _is_outdated_fasta(path: Path, settings: RegistrySettings) -> bool:
    """Return whether a FASTA predates the configured minimum build date."""
    if settings.min_fasta_date is None:
        return False
    build_date, _ = _fasta_build_date(path, settings)
    return build_date < settings.min_fasta_date


def _log_outdated_fasta(path: Path, settings: RegistrySettings) -> None:
    """Log one intentional pre-scan build-date exclusion."""
    build_date, source = _fasta_build_date(path, settings)
    logger.info(
        "skipped FASTA {}: build date {} (from {}) is before min_fasta_date={}",
        path,
        build_date,
        source,
        settings.min_fasta_date,
    )


def _is_excluded_fasta(path: Path, settings: RegistrySettings) -> bool:
    """Return whether a FASTA is excluded before any sequence is read."""
    return _is_oversized_fasta(path, settings) or _is_outdated_fasta(path, settings)


def _log_excluded_fasta(path: Path, settings: RegistrySettings) -> None:
    """Log the reason one FASTA was excluded before scanning."""
    if _is_oversized_fasta(path, settings):
        _log_oversized_fasta(path, settings)
    else:
        _log_outdated_fasta(path, settings)


def _record_rejection(
    error: FastaReadError,
    rejections: list[RejectedFasta] | None,
) -> None:
    rejection = RejectedFasta(path=Path(error.source_name), reason=error.reason)
    if rejections is not None:
        rejections.append(rejection)
    logger.error("rejected FASTA {}: {}", rejection.path, rejection.reason)


def _validate_complete_registry(connection: RegistryConnection) -> tuple[int, int, int]:
    """Validate integrity and return database, entry, and pair-row counts."""
    connection.check_physical_integrity()

    database_count = connection.scalar("SELECT COUNT(*) FROM databases", ())
    entry_count = connection.scalar("SELECT COUNT(*) FROM entries", ())
    full_database_count = connection.scalar(
        "SELECT COUNT(*) FROM databases WHERE detail_level = ?",
        (DetailLevel.FULL.value,),
    )
    expected_entry_count = connection.scalar(
        "SELECT COALESCE(SUM(entry_count), 0) FROM databases WHERE detail_level = ?",
        (DetailLevel.FULL.value,),
    )
    if entry_count != expected_entry_count:
        raise RegistryIntegrityError(
            f"Entry details are incomplete: expected {expected_entry_count}, found {entry_count}."
        )
    database_entry_mismatches = connection.scalar(
        """
        SELECT COUNT(*)
        FROM databases
        LEFT JOIN (
            SELECT database_id, COUNT(*) AS stored_entries
            FROM entries
            GROUP BY database_id
        ) AS stored ON stored.database_id = databases.id
        WHERE (
            databases.detail_level = 'full'
            AND COALESCE(stored.stored_entries, 0) != databases.entry_count
        ) OR (
            databases.detail_level = 'metadata_only'
            AND COALESCE(stored.stored_entries, 0) != 0
        )
        """,
        (),
    )
    if database_entry_mismatches:
        raise RegistryIntegrityError(
            f"Entry detail counts disagree for {database_entry_mismatches} databases."
        )
    metadata_entry_count = connection.scalar(
        """
        SELECT COUNT(*)
        FROM entries
        JOIN databases ON databases.id = entries.database_id
        WHERE databases.detail_level = ?
        """,
        (DetailLevel.METADATA_ONLY.value,),
    )
    if metadata_entry_count:
        raise RegistryIntegrityError(
            f"Metadata-only databases retain {metadata_entry_count} unexpected entry detail rows."
        )
    kind_stat_count = connection.scalar("SELECT COUNT(*) FROM database_kind_stats", ())
    expected_kind_stats = database_count * len(EntryKind)
    if kind_stat_count != expected_kind_stats:
        raise RegistryIntegrityError(
            f"Kind-statistics cache is incomplete: expected {expected_kind_stats}, found {kind_stat_count}."
        )
    invalid_kind_stats = connection.scalar(
        """
        SELECT COUNT(*)
        FROM database_kind_stats AS stats
        JOIN databases ON databases.id = stats.database_id
        WHERE (
            databases.detail_level = 'full'
            AND (
                stats.distinct_ids IS NULL
                OR stats.distinct_sequences IS NULL
                OR stats.distinct_descriptions IS NULL
                OR stats.distinct_pairs IS NULL
                OR stats.duplicate_id_occurrences IS NULL
                OR stats.conflicting_ids IS NULL
                OR stats.repeated_sequences IS NULL
                OR stats.id_fingerprint IS NULL
                OR stats.description_fingerprint IS NULL
                OR stats.content_fingerprint IS NULL
            )
        ) OR (
            databases.detail_level = 'metadata_only'
            AND (
                stats.distinct_ids IS NOT NULL
                OR stats.distinct_sequences IS NOT NULL
                OR stats.distinct_descriptions IS NOT NULL
                OR stats.distinct_pairs IS NOT NULL
                OR stats.duplicate_id_occurrences IS NOT NULL
                OR stats.conflicting_ids IS NOT NULL
                OR stats.repeated_sequences IS NOT NULL
                OR stats.id_fingerprint IS NOT NULL
                OR stats.description_fingerprint IS NOT NULL
                OR stats.content_fingerprint IS NOT NULL
            )
        )
        """,
        (),
    )
    if invalid_kind_stats:
        raise RegistryIntegrityError(
            f"Detail availability is inconsistent in {invalid_kind_stats} database kind-statistics rows."
        )

    pair_row_count = connection.scalar("SELECT COUNT(*) FROM database_pair_stats", ())
    expected_pair_rows = (
        full_database_count * (full_database_count - 1) // 2 * len(_MATERIALIZED_PAIR_KINDS)
    )
    if pair_row_count != expected_pair_rows:
        raise RegistryIntegrityError(
            f"Pair-statistics cache is incomplete: expected {expected_pair_rows}, found {pair_row_count}."
        )
    metadata_pair_count = connection.scalar(
        """
        SELECT COUNT(*)
        FROM database_pair_stats AS pairs
        JOIN databases AS low ON low.id = pairs.database_id_low
        JOIN databases AS high ON high.id = pairs.database_id_high
        WHERE low.detail_level = 'metadata_only' OR high.detail_level = 'metadata_only'
        """,
        (),
    )
    if metadata_pair_count:
        raise RegistryIntegrityError(
            f"Metadata-only databases participate in {metadata_pair_count} unexpected pair-statistics rows."
        )
    return database_count, entry_count, pair_row_count


def rebuild_registry(
    directory: Path,
    path: Path,
    settings: RegistrySettings,
    *,
    recursive: bool = False,
    rejections: list[RejectedFasta] | None = None,
) -> list[RegistryRecord]:
    """Bulk-build a registry, then create indexes and pair statistics once."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False
    ) as handle:
        temporary_path = Path(handle.name)
    temporary_path.unlink()
    started = time.perf_counter()
    try:
        # The staging file is named for atomic replacement, not for its engine, so
        # the backend has to be named rather than read off the suffix.
        with connect_registry(temporary_path, backend=factory.backend_for_path(path)) as connection:
            connection.configure_bulk()
            connection.create_tables()
            _reconcile_meta(connection, settings)
            discovered_paths = list(iter_fasta_files(directory, recursive=recursive))
            fasta_paths: list[Path] = []
            excluded_paths: list[Path] = []
            for discovered_path in discovered_paths:
                destination = (
                    excluded_paths if _is_excluded_fasta(discovered_path, settings) else fasta_paths
                )
                destination.append(discovered_path)
            total_files = len(fasta_paths)
            progress_interval = max(1, total_files // 20)
            accepted_files = 0
            accepted_entries = 0
            logger.info(
                "bulk load: discovered {} FASTA files; {} eligible below {:g} GiB; {} skipped",
                f"{len(discovered_paths):,}",
                f"{total_files:,}",
                settings.max_fasta_file_size_gib,
                f"{len(excluded_paths):,}",
            )
            if settings.min_fasta_date is not None:
                logger.info(
                    "bulk load: excluding FASTAs built before min_fasta_date={}",
                    settings.min_fasta_date,
                )
            for excluded_path in excluded_paths:
                _log_excluded_fasta(excluded_path, settings)
            for processed_files, fasta_path in enumerate(fasta_paths, start=1):
                try:
                    record = _index_fasta(
                        connection,
                        fasta_path,
                        settings,
                        root=directory,
                        refresh_pair_stats=False,
                        log_record=False,
                        progress_label=f"{processed_files:,}/{total_files:,}",
                    )
                    accepted_files += 1
                    accepted_entries += record.entry_count
                except FastaReadError as error:
                    _record_rejection(error, rejections)
                if processed_files % progress_interval == 0 or processed_files == total_files:
                    logger.info(
                        "bulk load: {}/{} files, {} accepted, {} rejected, {} entries ({:.1f}s)",
                        f"{processed_files:,}",
                        f"{total_files:,}",
                        f"{accepted_files:,}",
                        f"{processed_files - accepted_files:,}",
                        f"{accepted_entries:,}",
                        time.perf_counter() - started,
                    )

            index_started = time.perf_counter()
            logger.info("bulk indexes: creating entry lookup indexes")
            connection.create_entry_indexes()
            connection.execute("ANALYZE")
            logger.info(
                "bulk indexes: entry indexes ready ({:.1f}s)", time.perf_counter() - index_started
            )

            pair_started = time.perf_counter()
            logger.info("pair statistics: materializing all database pairs")
            _materialize_all_pair_stats(connection)
            connection.create_pair_indexes()
            connection.execute("ANALYZE")
            logger.info("pair statistics: complete ({:.1f}s)", time.perf_counter() - pair_started)

            validation_started = time.perf_counter()
            database_count, entry_count, pair_row_count = _validate_complete_registry(connection)
            full_database_count = connection.scalar(
                "SELECT COUNT(*) FROM databases WHERE detail_level = ?",
                (DetailLevel.FULL.value,),
            )
            scanned_entry_count = connection.scalar(
                "SELECT COALESCE(SUM(entry_count), 0) FROM databases",
                (),
            )
            logger.info(
                "registry validation: integrity ok; {} databases ({} full-detail, {} metadata-only), "
                "{} entries scanned, {} details stored, {} pair rows ({:.1f}s)",
                f"{database_count:,}",
                f"{full_database_count:,}",
                f"{database_count - full_database_count:,}",
                f"{scanned_entry_count:,}",
                f"{entry_count:,}",
                f"{pair_row_count:,}",
                time.perf_counter() - validation_started,
            )
            records = list_databases(connection)
        os.replace(temporary_path, path)
        logger.info(
            "published registry: {} ({:.2f} GB, {:.1f}s total)",
            path,
            path.stat().st_size / 1_000_000_000,
            time.perf_counter() - started,
        )
        return records
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def update_registry(
    directory: Path,
    path: Path,
    settings: RegistrySettings,
    *,
    force: bool = False,
    prune: bool = False,
    recursive: bool = False,
    rejections: list[RejectedFasta] | None = None,
) -> list[RegistryRecord]:
    """Index changed files while removing and reporting invalid FASTAs."""
    with connect_registry(path) as connection:
        initialize_registry(connection, settings)
        existing = {record.relative_path: record for record in list_databases(connection)}
        seen: set[str] = set()
        for fasta_path in iter_fasta_files(directory, recursive=recursive):
            relative_path = _relative_path(fasta_path, directory)
            seen.add(relative_path)
            if _is_excluded_fasta(fasta_path, settings):
                _log_excluded_fasta(fasta_path, settings)
                # Drop any record indexed before the exclusion applied, so the
                # registry always matches the current filters.
                with connection.transaction():
                    connection.delete_database(relative_path)
                continue
            stat = fasta_path.stat()
            previous = existing.get(relative_path)
            unchanged = (
                previous
                and previous.file_size_bytes == stat.st_size
                and previous.mtime_ns == stat.st_mtime_ns
            )
            if force or not unchanged:
                try:
                    index_fasta(connection, fasta_path, settings, root=directory)
                except FastaReadError as error:
                    _record_rejection(error, rejections)
                    with connection.transaction():
                        connection.delete_database(relative_path)
        if prune:
            # ``recursive`` controls discovery, not the meaning of missing. A
            # top-level maintenance run must not delete a still-existing nested
            # record that was indexed by an earlier recursive sweep.
            missing = {
                relative_path
                for relative_path in set(existing) - seen
                if not (directory / relative_path).is_file()
            }
            with connection.transaction():
                for relative_path in missing:
                    connection.delete_database(relative_path)
        return list_databases(connection)


def read_registry(path: Path) -> list[RegistryRecord]:
    """Compatibility reader for an initialized SQLite registry."""
    if not path.exists():
        return []
    with connect_registry(path) as connection:
        current_version = connection.schema_version()
        has_databases = connection.has_table("databases")
        if current_version == 0 and not has_databases:
            return []
        if current_version != SCHEMA_VERSION:
            raise RegistrySchemaError(current_version, path=path)
        return list_databases(connection) if has_databases else []
