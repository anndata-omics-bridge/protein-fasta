"""Fixtures for shared protein-FASTA tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.registry_support import BackendSettings, Settings


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--registry-backend",
        default="sqlite",
        choices=("sqlite", "duckdb"),
        help="Storage engine for registry tests (default: sqlite).",
    )


@pytest.fixture(scope="session")
def registry_backend(request: pytest.FixtureRequest) -> str:
    return str(request.config.getoption("--registry-backend"))


@pytest.fixture
def settings(
    registry_backend: str,
    tmp_path: Path,
) -> Settings:
    return Settings(
        fasta_root=tmp_path / "fastas",
        registry_dir=tmp_path / "registry",
        registry=BackendSettings(backend=registry_backend),
    )
