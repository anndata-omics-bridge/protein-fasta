"""Build configured database names and FASTA filenames.

``build_fasta_name`` generates the canonical modern form
``p<project>_db<dbn>_<description>_<YYYYMMDD>.fasta`` (with a ``_d`` token before
the date for decoy databases).
"""

from __future__ import annotations

import datetime
import re
from collections.abc import Iterable
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class DescriptionRules:
    """Runtime flags used to derive compact database descriptions."""

    mode_flags: dict[str, str]


@dataclass(frozen=True, slots=True)
class DatabaseNaming:
    """Resolved runtime grammar for database and FASTA names."""

    default_dbname: str
    dbname: dict[str, str]
    filename: dict[str, str]
    description: DescriptionRules
    decoy_token: str
    entrapment_token: str
    separator: str
    date_format: str
    extension: str


class _BlankDefault(dict[str, object]):
    """A format mapping that renders missing/None fields as an empty string."""

    def __missing__(self, key: str) -> str:
        return ""


def build_dbname(*, config: DatabaseNaming, template: str | None = None, **fields: object) -> str:
    """Render the ``{dbname}`` (no date/decoy) from a named ``naming.dbname`` pattern.

    ``fields`` are substituted into the pattern (e.g. ``project``, ``dbn``,
    ``description``, ``taxid``); missing or empty fields collapse away so a pattern
    like ``p{project}_db{dbn}_{description}`` with no description yields ``p999_db1``.
    """
    name = template or config.default_dbname
    pattern = config.dbname.get(name)
    if pattern is None:
        raise KeyError(
            f"unknown naming.dbname pattern {name!r}; configured: {sorted(config.dbname)}"
        )
    cleaned = {key: ("" if value is None else str(value)) for key, value in fields.items()}
    rendered = pattern.format_map(_BlankDefault(cleaned))
    separator = re.escape(config.separator)
    rendered = re.sub(f"(?:{separator}){{2,}}", config.separator, rendered)
    return rendered.strip(config.separator)


def build_fasta_name(
    *,
    config: DatabaseNaming,
    template: str | None = None,
    date: datetime.date,
    decoy: bool,
    entrapment: bool = False,
    **fields: object,
) -> str:
    """Render a full FGCZ database filename from the configured ``naming.filename`` patterns.

    The dbname is rendered from the chosen ``template``, then substituted into the
    ``decoy`` or ``nondecoy`` filename pattern (over ``{dbname}``, ``{decoy_token}``,
    ``{date}``, ``{extension}``) — so both naming levels live in the config.
    """
    dbname = build_dbname(config=config, template=template, **fields)
    key = "decoy" if decoy else "nondecoy"
    if entrapment:
        key = f"entrapment_{key}"
    pattern = config.filename.get(key)
    if pattern is None:
        raise KeyError(f"unknown filename pattern {key!r}; configured: {sorted(config.filename)}")
    return pattern.format_map(
        _BlankDefault(
            {
                "dbname": dbname,
                "decoy_token": config.decoy_token,
                "entrapment_token": config.entrapment_token,
                "date": date.strftime(config.date_format),
                "extension": config.extension,
            }
        )
    )


def build_description(*, config: DatabaseNaming, taxids: Iterable[object], mode: str) -> str:
    """Prefill the terse ``{description}`` from a proteome selection.

    Joins the selected ``taxids`` with the naming ``separator`` and appends the
    per-mode flag from ``naming.description.mode_flags`` (empty taxids or an empty
    flag collapse away). E.g. taxids ``[9606, 4932]`` in ``swissprot`` mode ->
    ``9606_4932_sp``. The full detail goes into the sentinel annotation, not here.
    """
    parts = [str(taxid) for taxid in taxids if taxid not in (None, "")]
    flag = config.description.mode_flags.get(mode, "")
    if flag:
        parts.append(flag)
    return config.separator.join(parts)
