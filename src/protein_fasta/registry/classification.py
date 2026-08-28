"""Resolve independent diagnostic labels into configured registry entry kinds."""

from __future__ import annotations

from collections.abc import Container
from dataclasses import dataclass

from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.registry.kinds import EntryKind

UNLABELLED_CONTAMINANT_GROUP = "unlabelled"


@dataclass(frozen=True, slots=True)
class ContaminantBlockState:
    """Marker provenance waiting for or already applied to ordinary entries."""

    group: str
    has_entries: bool = False


def normalize_contaminant_group(group: str | None) -> str:
    """Return a non-empty normalized contaminant group."""
    normalized = group.strip() if group is not None else ""
    return normalized or UNLABELLED_CONTAMINANT_GROUP


def sentinel_details(
    raw_header: str,
    decoy_prefix: str,
    /,
) -> tuple[str, str | None, str | None]:
    """Return normalized metadata text and optional contaminant marker details."""
    core = raw_header.removeprefix(decoy_prefix)
    fields = core.split("|", 2)
    record_id = fields[1] if len(fields) > 1 else ""
    if record_id.startswith("Cont_"):
        return core, record_id, normalize_contaminant_group(record_id.removeprefix("Cont_"))
    return core, None, None


def operational_kind(classifications: Container[str], /) -> EntryKind:
    """Resolve overlapping shared labels by the one FGCZ precedence rule."""
    if "sentinel" in classifications:
        return EntryKind.SENTINEL
    if "decoy" in classifications:
        return EntryKind.DECOY
    if "entrapment" in classifications:
        return EntryKind.ENTRAPMENT
    if "contaminant" in classifications:
        return EntryKind.CONTAMINANT
    return EntryKind.TARGET


def classify_identifier(identifier: str, rules: DiagnosticRules, /) -> EntryKind:
    """Diagnose one identifier and resolve its FGCZ operational kind."""
    _, classifications = rules.diagnose_identifier(identifier)
    return operational_kind(classifications)


def classify_record(
    raw_header: str,
    classifications: Container[str],
    block_state: ContaminantBlockState | None,
    decoy_prefix: str,
    /,
) -> tuple[EntryKind, str | None, ContaminantBlockState | None]:
    """Resolve one record while honoring marker-delimited contaminant blocks."""
    configured_kind = operational_kind(classifications)
    if configured_kind is EntryKind.SENTINEL:
        _, marker_id, marker_group = sentinel_details(raw_header, decoy_prefix)
        if marker_id is None:
            return EntryKind.SENTINEL, None, None
        if marker_group is None:
            raise RuntimeError("contaminant marker has no normalized group")
        if block_state is not None and not block_state.has_entries:
            marker_group = UNLABELLED_CONTAMINANT_GROUP
        return EntryKind.SENTINEL, None, ContaminantBlockState(marker_group)
    if configured_kind is EntryKind.DECOY:
        return EntryKind.DECOY, None, block_state
    if configured_kind is EntryKind.ENTRAPMENT:
        return EntryKind.ENTRAPMENT, None, block_state
    if block_state is not None:
        next_state = ContaminantBlockState(block_state.group, has_entries=True)
        return EntryKind.CONTAMINANT, block_state.group, next_state
    if configured_kind is EntryKind.CONTAMINANT:
        return EntryKind.CONTAMINANT, UNLABELLED_CONTAMINANT_GROUP, None
    return EntryKind.TARGET, None, None
