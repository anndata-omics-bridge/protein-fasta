"""CLI parity tests for UniProt workflows."""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import polars as pl
import pytest

import protein_fasta.uniprot_download as download_module
from protein_fasta.cli import uniprot_download, uniprot_proteomes
from protein_fasta.uniprot.transport import UniProtTransport

_BASE_URL = "https://unit.test"


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
        uniprot_download(request_path)
    finally:
        client.close()

    assert output.read_bytes() == fasta
    result = json.loads(output.with_suffix(".fasta.result.json").read_text(encoding="utf-8"))
    assert result["actual_entry_count"] == 1
    assert result["provider_query"] == "(proteome:UP000005640) AND (reviewed:true)"


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
