"""Architecture facts not represented directly by Import Linter."""

import ast
from pathlib import Path

_CHILDREN = {
    "analytics",
    "build",
    "diagnostics",
    "frame_formats",
    "peptide",
    "reading",
    "registry",
    "schema",
    "validation",
}
_MULTI_CHILD_COMPOSERS = {
    "analytics_compile",
    "artifact_io",
    "candidate_analysis",
    "cli",
    "compile",
    "database_build",
    "decoy_database",
    "decoy_report",
    "frame",
    "frame_compile",
    "peptide_workflow",
    "protein_input",
    "record",
    "registry_workflow",
    "uniprot_catalog",
}


def test_package_initializers_are_empty() -> None:
    root = Path(__file__).parents[1] / "src" / "protein_fasta"
    initializers = sorted(root.rglob("__init__.py"))
    assert initializers
    assert all(not path.read_text().strip() for path in initializers)


def test_only_named_root_modules_compose_multiple_children() -> None:
    root = Path(__file__).parents[1] / "src" / "protein_fasta"
    observed: set[str] = set()
    for path in root.glob("*.py"):
        children = {
            child
            for node in ast.walk(ast.parse(path.read_text()))
            if isinstance(node, ast.ImportFrom)
            and node.module is not None
            and node.module.startswith("protein_fasta.")
            for child in [node.module.split(".")[1]]
            if child in _CHILDREN
        }
        if len(children) > 1:
            observed.add(path.stem)
    assert observed == _MULTI_CHILD_COMPOSERS


def test_cli_composes_only_authorized_package_components() -> None:
    path = Path(__file__).parents[1] / "src" / "protein_fasta" / "cli.py"
    imported_modules = {
        node.module
        for node in ast.walk(ast.parse(path.read_text()))
        if isinstance(node, ast.ImportFrom)
        and node.module is not None
        and node.module.startswith("protein_fasta.")
    }
    assert {module.split(".")[1] for module in imported_modules} == {
        "analytics",
        "analytics_compile",
        "compile",
        "candidate_analysis",
        "database_build",
        "decoy_database",
        "decoy_report",
        "diagnostic_summary",
        "documents",
        "frame",
        "peptide_workflow",
        "protein_input",
        "record",
        "registry",
        "registry_workflow",
        "schema",
        "uniprot_catalog",
        "uniprot_download",
    }
