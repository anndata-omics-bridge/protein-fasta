"""Documentation structure tests for contextual JSON examples."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import cast

_ROOT = Path(__file__).parents[1]


def _json_leaf_paths(value: object, prefix: str = "") -> set[str]:
    """Collect dotted paths to every displayed scalar field."""
    if isinstance(value, dict):
        paths: set[str] = set()
        mapping = cast(dict[object, object], value)
        for raw_key, raw_child in mapping.items():
            if not isinstance(raw_key, str):
                continue
            key: str = raw_key
            child: object = raw_child
            path = f"{prefix}.{key}" if prefix else key
            child_paths = _json_leaf_paths(child, path)
            paths.update(child_paths or {path})
        return paths
    if isinstance(value, list):
        paths: set[str] = set()
        for raw_child in cast(list[object], value):
            child: object = raw_child
            paths.update(_json_leaf_paths(child, prefix))
        return paths
    return set()


def test_static_json_examples_are_contextual_and_explained() -> None:
    """Every literal JSON example has a consumer and names every displayed field."""
    for relative_path in ("docs/cli_walkthrough.md", "docs/configuration.md"):
        text = (_ROOT / relative_path).read_text(encoding="utf-8")
        for match in re.finditer(r"```json\n(?P<payload>.*?)\n```", text, re.DOTALL):
            section_start = text.rfind("\n## ", 0, match.start())
            section_end = text.find("\n## ", match.end())
            section = text[section_start : None if section_end < 0 else section_end]
            before_json = section[: match.start() - section_start]
            assert "protein-fasta " in before_json
            payload = json.loads(match.group("payload"))
            for path in _json_leaf_paths(payload):
                leaf = path.rsplit(".", maxsplit=1)[-1]
                assert f"`{path}`" in section or f"`{leaf}`" in section


def test_walkthrough_has_no_detached_request_gallery() -> None:
    """Generated requests stay inside the workflow section that consumes them."""
    text = (_ROOT / "docs/cli_walkthrough.md").read_text(encoding="utf-8")
    assert "Request documents used by the workflow" not in text
    for request_name in (
        "protein-input.parquet.request.json",
        "build.request.json",
        "p1_db1_walkthrough_d_20260827.fasta.request.json",
        "candidate-comparisons.parquet.request.json",
        "derived-protein-input.parquet.request.json",
        "peptides.request.json",
        "peptide-comparison.parquet.request.json",
        "decoy-methods.parquet.request.json",
    ):
        assert request_name in text
