"""Exercise the installed package surface and bundled documents."""

from __future__ import annotations

import importlib
import json
import pkgutil
from importlib import resources
from importlib.resources.abc import Traversable
from typing import Final

import protein_fasta
from protein_fasta.documents import (
    load_builtin_database_build_profile,
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
    load_builtin_enzyme_document,
    load_builtin_header_format_catalog,
)

_OPTIONAL_GENERATION_MODULES: Final = {
    "protein_fasta.database.decoy_advanced",
    "protein_fasta.database.entrapment_advanced",
}


def main() -> None:
    """Import installed modules and load every packaged rule/schema resource."""
    discovered = {
        module.name
        for module in pkgutil.walk_packages(
            protein_fasta.__path__,
            prefix="protein_fasta.",
        )
    }
    for module_name in sorted(discovered - _OPTIONAL_GENERATION_MODULES):
        importlib.import_module(module_name)

    load_builtin_database_build_profile()
    load_builtin_diagnostic_document()
    load_builtin_entry_classifier_document()
    load_builtin_enzyme_document()
    load_builtin_header_format_catalog()

    document_root = resources.files("protein_fasta").joinpath("documents")
    json_resources = _json_resources(document_root)
    if not json_resources:
        raise AssertionError("installed wheel contains no packaged JSON documents")
    for resource in json_resources:
        json.loads(resource.read_text(encoding="utf-8"))


def _json_resources(root: Traversable) -> tuple[Traversable, ...]:
    """Return all JSON descendants of one installed resource directory."""
    found: list[Traversable] = []
    pending = [root]
    while pending:
        current = pending.pop()
        for child in current.iterdir():
            if child.is_dir():
                pending.append(child)
            elif child.name.endswith(".json"):
                found.append(child)
    return tuple(found)


if __name__ == "__main__":
    main()
