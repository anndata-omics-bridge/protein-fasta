"""CLI parity tests for UniProt workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import httpx
import polars as pl
import pytest

import protein_fasta.uniprot_catalog as catalog_module
import protein_fasta.uniprot_download as download_module
from protein_fasta.cli import app, uniprot_catalog, uniprot_download, uniprot_proteomes
from protein_fasta.uniprot.transport import UniProtTransport

_BASE_URL = "https://unit.test"


def test_uniprot_catalog_direct_authors_request_and_replays(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    proteome = {
        "id": "UP000005640",
        "proteomeType": "Reference proteome",
        "proteinCount": 1,
        "geneCount": 1,
        "proteomeStatistics": {"reviewedProteinCount": 1, "unreviewedProteinCount": 0},
        "taxonomy": {"taxonId": 9606, "scientificName": "Homo sapiens"},
    }

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            json={"results": [proteome]},
            headers={"x-uniprot-release": "2026_03", "x-total-results": "1"},
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    transport = UniProtTransport(base_url=_BASE_URL, client=client)

    def transport_for_cli(*, timeout_seconds: float) -> UniProtTransport:
        assert timeout_seconds == 120.0
        return transport

    monkeypatch.setattr(catalog_module, "UniProtTransport", transport_for_cli)
    direct = tmp_path / "direct"
    replay = tmp_path / "replay"
    saved = tmp_path / "saved"
    saved_request = tmp_path / "requests" / "catalog.json"
    try:
        uniprot_catalog(output=direct)
        request_path = direct / "uniprot-catalog.request.json"
        assert json.loads(request_path.read_text(encoding="utf-8")) == {
            "output_dir": ".",
            "schema_version": "0.1",
            "selection": {"type": "reference"},
            "timeout_seconds": 120.0,
        }
        uniprot_catalog(request=request_path, output=replay)
        uniprot_catalog(output=saved, save=saved_request)
    finally:
        client.close()

    assert len(list(direct.glob("*.parquet"))) == 1
    assert len(list(replay.glob("*.parquet"))) == 1
    assert len(list(saved.glob("*.parquet"))) == 1
    assert json.loads(saved_request.read_text(encoding="utf-8"))["output_dir"] == "../saved"


def test_uniprot_download_command_uses_the_shared_artifact_workflow(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request_path = tmp_path / "download.json"
    output = tmp_path / "human.fasta"
    request_path.write_text(
        json.dumps(
            {
                "selection": {"type": "proteome_id", "proteome_id": "UP000005640"},
                "acquisition": {"type": "swissprot"},
                "output_fasta": output.name,
            }
        ),
        encoding="utf-8",
    )
    fasta = b">sp|P1|ONE One\nMPEPTIDE\n"

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/UP000005640":
            return httpx.Response(
                200,
                json={
                    "id": "UP000005640",
                    "taxonomy": {"taxonId": 9606, "scientificName": "Homo sapiens"},
                },
            )
        return httpx.Response(
            200,
            content=fasta,
            headers={"x-uniprot-release": "2026_03", "x-total-results": "1"},
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    transport = UniProtTransport(base_url=_BASE_URL, client=client)

    def transport_for_cli(*, timeout_seconds: float) -> UniProtTransport:
        assert timeout_seconds == 120.0
        return transport

    monkeypatch.setattr(download_module, "UniProtTransport", transport_for_cli)
    try:
        uniprot_download(request=request_path)
    finally:
        client.close()

    assert output.read_bytes() == fasta
    result = json.loads(output.with_suffix(".fasta.result.json").read_text(encoding="utf-8"))
    assert result["actual_entry_count"] == 1
    assert result["provider_query"] == "(proteome:UP000005640) AND (reviewed:true)"


@pytest.mark.parametrize(
    ("mode", "expected_query", "expected_acquisition"),
    [
        ("reviewed", "(proteome:UP000005640) AND (reviewed:true)", "swissprot"),
        ("canonical", "(proteome:UP000005640)", "swissprot_trembl"),
    ],
)
def test_uniprot_download_command_accepts_direct_proteome_arguments(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: Literal["reviewed", "canonical"],
    expected_query: str,
    expected_acquisition: str,
) -> None:
    fasta = b">sp|P1|ONE One\nMPEPTIDE\n"

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/UP000005640":
            return httpx.Response(
                200,
                json={
                    "id": "UP000005640",
                    "taxonomy": {"taxonId": 9606, "scientificName": "Homo sapiens"},
                },
            )
        return httpx.Response(
            200,
            content=fasta,
            headers={"x-uniprot-release": "2026_03", "x-total-results": "1"},
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    transport = UniProtTransport(base_url=_BASE_URL, client=client)

    def transport_for_cli(*, timeout_seconds: float) -> UniProtTransport:
        assert timeout_seconds == 45.0
        return transport

    monkeypatch.setattr(download_module, "UniProtTransport", transport_for_cli)
    monkeypatch.chdir(tmp_path)
    try:
        uniprot_download("UP000005640", mode, timeout=45.0)
    finally:
        client.close()

    output = tmp_path / f"UP000005640_{mode}.fasta"
    assert output.read_bytes() == fasta
    request = json.loads(output.with_suffix(".fasta.request.json").read_text(encoding="utf-8"))
    assert request == {
        "acquisition": {"type": expected_acquisition},
        "output_fasta": output.name,
        "schema_version": "0.1",
        "selection": {"proteome_id": "UP000005640", "type": "proteome_id"},
        "timeout_seconds": 45.0,
    }
    result = json.loads(output.with_suffix(".fasta.result.json").read_text(encoding="utf-8"))
    assert result["provider_query"] == expected_query
    assert result["effective_request"]["acquisition"]["type"] == expected_acquisition


def test_uniprot_download_command_accepts_direct_one_protein_per_gene(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    canonical_record = {
        "canonicalProtein": {
            "id": "P1",
            "uniProtkbId": "ONE_HUMAN",
            "proteinName": "One",
            "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
            "geneName": "ONE",
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "sequence": {"value": "MPEPTIDE"},
        }
    }

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/UP000005640":
            return httpx.Response(
                200,
                json={
                    "id": "UP000005640",
                    "taxonomy": {"taxonId": 9606, "scientificName": "Homo sapiens"},
                },
            )
        return httpx.Response(
            200,
            json={"results": [canonical_record]},
            headers={"x-uniprot-release": "2026_03", "x-total-results": "1"},
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    transport = UniProtTransport(base_url=_BASE_URL, client=client)

    def transport_for_cli(*, timeout_seconds: float) -> UniProtTransport:
        assert timeout_seconds == 120.0
        return transport

    monkeypatch.setattr(download_module, "UniProtTransport", transport_for_cli)
    monkeypatch.chdir(tmp_path)
    try:
        with pytest.raises(SystemExit) as exit_info:
            app(["uniprot-download", "UP000005640", "opg"])
    finally:
        client.close()

    assert exit_info.value.code == 0
    output = tmp_path / "UP000005640_opg.fasta"
    assert output.read_text(encoding="utf-8").startswith(">sp|P1|ONE_HUMAN")
    request = json.loads(output.with_suffix(".fasta.request.json").read_text(encoding="utf-8"))
    assert request["selection"] == {
        "proteome_id": "UP000005640",
        "type": "proteome_id",
    }
    assert request["acquisition"] == {"type": "one_seq_per_gene"}
    result = json.loads(output.with_suffix(".fasta.result.json").read_text(encoding="utf-8"))
    assert result["provider_query"] == "(upid:UP000005640)"
    assert result["effective_request"]["acquisition"]["type"] == "one_seq_per_gene"


def test_uniprot_proteomes_command_filters_a_local_catalog(tmp_path: Path) -> None:
    catalog_path = tmp_path / "catalog.parquet"
    table_path = tmp_path / "human.csv"
    pl.DataFrame(
        {
            "proteome_id": ["UP1", "UP2"],
            "taxid": [9606, 10090],
            "organism": ["Homo sapiens", "Mus musculus"],
            "proteome_type": ["Reference proteome", "Reference proteome"],
            "swissprot": [1, 2],
            "swissprot_trembl": [3, 4],
            "one_seq_per_gene": [1, 2],
        },
        schema={
            "proteome_id": pl.String,
            "taxid": pl.Int64,
            "organism": pl.String,
            "proteome_type": pl.String,
            "swissprot": pl.Int64,
            "swissprot_trembl": pl.Int64,
            "one_seq_per_gene": pl.Int64,
        },
    ).write_parquet(catalog_path)

    uniprot_proteomes(catalog_path, table_path, query="9606")

    assert pl.read_csv(table_path).get_column("proteome_id").to_list() == ["UP1"]
