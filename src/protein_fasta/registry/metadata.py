"""Parse configured database metadata records during registry indexing."""

from __future__ import annotations

import re
from dataclasses import dataclass

from protein_fasta.schema.build import MetadataDocument

_LEADING_DATE = re.compile(r"^(?P<date>\d{4}-\d{2}-\d{2}|\d{8})\s*,?\s*")


@dataclass(frozen=True, slots=True)
class DatabaseMetadata:
    """Parsed database metadata record."""

    dbname: str
    date: str | None
    description: str
    raw_header: str


def parse_database_metadata(
    header: str,
    document: MetadataDocument,
    /,
) -> DatabaseMetadata | None:
    """Parse a configured database metadata header, excluding section markers."""
    text = header.lstrip(">")
    if not text.startswith(f"{document.id_namespace}|"):
        return None

    fields = text.split("|", 2)
    if len(fields) < 3:
        return None
    _, record_id, description = fields
    if record_id.startswith("Cont_"):
        return None

    description = description.lstrip()
    date: str | None = None
    match = _LEADING_DATE.match(description)
    if match:
        date = match.group("date")
        description = description[match.end() :]

    return DatabaseMetadata(
        dbname=record_id,
        date=date,
        description=description.strip(),
        raw_header=text,
    )
