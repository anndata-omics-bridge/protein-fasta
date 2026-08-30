"""Build configured database metadata and contaminant-marker records.

The sentinel is a self-documenting pseudo-FASTA record written as the first entry
of every produced database::

    >aa|<dbname>|<date> <description>, generated w <tool> and installed by <person>, fgcz
    CRAPCRAPCRAP

Contaminant blocks are delimited by fixed section-marker records
(``aa|Cont_specialContaminants`` / ``aa|Cont_UniversalContaminants``). This module
handles construction only; registry parsing is a separate capability.
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class DatabaseMetadata:
    """Resolved runtime policy for sentinel and section-marker records."""

    id_namespace: str
    body_sequence: str
    marker_body_sequence: str
    date_format: str
    provenance_template: str
    tool: str
    installer: str
    org: str


def build_sentinel_header(
    dbname: str,
    description: str,
    date: datetime.date | None,
    config: DatabaseMetadata,
) -> str:
    """Build a database sentinel header line (without the leading ``>``)."""
    provenance = config.provenance_template.format(
        tool=config.tool, installer=config.installer, org=config.org
    )
    if date is not None:
        date_string = date.strftime(config.date_format)
        body = f"{date_string} {description}".strip()
    else:
        body = description
    return f"{config.id_namespace}|{dbname}|{body}, {provenance}"


def build_section_marker_header(
    record_id: str,
    description: str,
    config: DatabaseMetadata,
) -> str:
    """Build a contaminant section-marker header line (without the leading ``>``)."""
    return f"{config.id_namespace}|{record_id}|{description}"
