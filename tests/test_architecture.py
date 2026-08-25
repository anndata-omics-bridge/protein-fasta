"""Architecture facts not represented directly by Import Linter."""

from pathlib import Path


def test_package_initializers_are_empty() -> None:
    root = Path(__file__).parents[1] / "src" / "protein_fasta"
    initializers = sorted(root.rglob("__init__.py"))
    assert initializers
    assert all(not path.read_text().strip() for path in initializers)
