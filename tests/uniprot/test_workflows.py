"""UniProt provider, catalog, and download workflow tests."""

from __future__ import annotations

import datetime
from pathlib import Path

import httpx
import pytest

from protein_fasta.analytics.hashing import file_checksum
from protein_fasta.schema.uniprot import (
    AllProteomesDocument,
    CanonicalGeneDownloadDocument,
    CompleteDownloadDocument,
    ProteomeIdSelectionDocument,
    ProteomeQueryDocument,
    TaxonomySelectionDocument,
    UniProtCatalogRequestDocument,
    UniProtDownloadRequestDocument,
)
from protein_fasta.uniprot.provider_rows import canonical_fasta_header
from protein_fasta.uniprot.transport import UniProtTransport
from protein_fasta.uniprot_catalog import (
    catalog_query,
    filter_uniprot_catalog,
    latest_uniprot_catalog,
    read_uniprot_catalog,
    sync_uniprot_catalog,
)
from protein_fasta.uniprot_download import resolve_uniprot_download, run_uniprot_download

_BASE_URL = "https://unit.test"

_HUMAN_PROTEOME: dict[str, object] = {
    "id": "UP000005640",
    "proteomeType": "Reference proteome",
    "proteinCount": 2,
    "geneCount": 2,
    "proteomeStatistics": {
        "reviewedProteinCount": 1,
        "unreviewedProteinCount": 1,
    },
    "taxonomy": {"taxonId": 9606, "scientificName": "Homo sapiens"},
}

_MOUSE_PROTEOME: dict[str, object] = {
    "id": "UP000000589",
    "proteomeType": "Reference proteome",
    "proteinCount": 1,
    "geneCount": 1,
    "proteomeStatistics": {
        "reviewedProteinCount": 1,
        "unreviewedProteinCount": 0,
    },
    "taxonomy": {"taxonId": 10090, "scientificName": "Mus musculus"},
}


def _transport(handler: httpx.MockTransport) -> tuple[httpx.Client, UniProtTransport]:
    client = httpx.Client(transport=handler)
    return client, UniProtTransport(base_url=_BASE_URL, client=client)


def _explicit_request(tmp_path: Path) -> UniProtDownloadRequestDocument:
    return UniProtDownloadRequestDocument(
        selection=ProteomeIdSelectionDocument(proteome_id="UP000005640"),
        output_fasta=Path("human.fasta"),
    )


def test_catalog_selection_variants_compile_to_one_query(tmp_path: Path) -> None:
    assert catalog_query(UniProtCatalogRequestDocument(output_dir=tmp_path)) == "reference:true"
    assert (
        catalog_query(
            UniProtCatalogRequestDocument(
                output_dir=tmp_path,
                selection=AllProteomesDocument(),
            )
        )
        == "*"
    )
    assert (
        catalog_query(
            UniProtCatalogRequestDocument(
                output_dir=tmp_path,
                selection=ProteomeQueryDocument(query="taxonomy_id:9606"),
            )
        )
        == "taxonomy_id:9606"
    )


def test_canonical_header_reconstructs_standard_uniprot_fields() -> None:
    canonical: dict[str, object] = {
        "id": "A0A024RBG1",
        "uniProtkbId": "NUD4B_HUMAN",
        "proteinName": "Diphosphoinositol polyphosphate phosphohydrolase NUDT4B",
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "geneName": "NUDT4B",
        "proteinExistence": "1: Evidence at protein level",
        "sequenceVersion": 1,
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
    }

    assert canonical_fasta_header(canonical) == (
        "sp|A0A024RBG1|NUD4B_HUMAN "
        "Diphosphoinositol polyphosphate phosphohydrolase NUDT4B "
        "OS=Homo sapiens OX=9606 GN=NUDT4B PE=1 SV=1"
    )


def test_reviewed_download_publishes_fasta_and_exact_evidence(tmp_path: Path) -> None:
    fasta = b">sp|P1|ONE One OS=Homo sapiens OX=9606\nMPEPTIDE\n"

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/UP000005640":
            return httpx.Response(200, json=_HUMAN_PROTEOME)
        assert request.url.path == "/uniprotkb/stream"
        assert request.url.params["query"] == "(proteome:UP000005640) AND (reviewed:true)"
        return httpx.Response(
            200,
            content=fasta,
            headers={"x-uniprot-release": "2026_03", "x-total-results": "1"},
        )

    client, transport = _transport(httpx.MockTransport(handler))
    try:
        effective = resolve_uniprot_download(
            _explicit_request(tmp_path),
            request_base=tmp_path,
        )
        execution = run_uniprot_download(effective, transport=transport)
    finally:
        client.close()

    assert execution.fasta_path.read_bytes() == fasta
    assert execution.document.actual_entry_count == 1
    assert execution.document.observed_releases == ("2026_03",)
    assert execution.document.provider_reported_counts == (1,)
    assert execution.document.warnings == ()
    assert execution.document.artifact.checksum == file_checksum(execution.fasta_path)
    assert execution.effective_request_path.is_file()
    assert execution.result_path.is_file()


def test_taxonomy_resolution_prefers_reference_then_records_fallback(tmp_path: Path) -> None:
    queries: list[str] = []
    fasta = b">tr|P2|TWO Two\nPEPTIDE\n"

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/search":
            query = request.url.params["query"]
            queries.append(query)
            if "proteome_type:1" in query:
                return httpx.Response(200, json={"results": []})
            return httpx.Response(200, json={"results": [_HUMAN_PROTEOME]})
        return httpx.Response(200, content=fasta, headers={"x-total-results": "1"})

    request = UniProtDownloadRequestDocument(
        selection=TaxonomySelectionDocument(taxid=9606),
        acquisition=CompleteDownloadDocument(),
        output_fasta=Path("complete.fasta"),
    )
    client, transport = _transport(httpx.MockTransport(handler))
    try:
        execution = run_uniprot_download(
            resolve_uniprot_download(request, request_base=tmp_path),
            transport=transport,
        )
    finally:
        client.close()

    assert queries == [
        "(taxonomy_id:9606) AND (proteome_type:1)",
        "(taxonomy_id:9606)",
    ]
    assert execution.document.resolved_proteome.resolution_method == "taxonomy_fallback"
    assert execution.document.resolved_proteome.resolution_query == "(taxonomy_id:9606)"


def test_canonical_gene_download_pages_and_retains_disagreement_evidence(tmp_path: Path) -> None:
    first = {
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
    second = {
        "canonicalProtein": {
            "id": "P2",
            "uniProtkbId": "TWO_HUMAN",
            "proteinName": "Two",
            "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
            "entryType": "UniProtKB unreviewed (TrEMBL)",
            "sequence": {"value": "PEPTIDER"},
        }
    }

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/proteomes/UP000005640":
            return httpx.Response(200, json=_HUMAN_PROTEOME)
        if request.url.params.get("cursor") == "next":
            return httpx.Response(
                200,
                json={"results": [second]},
                headers={"x-uniprot-release": "2026_04", "x-total-results": "3"},
            )
        return httpx.Response(
            200,
            json={"results": [first, {"canonicalProtein": {}}]},
            headers={
                "x-uniprot-release": "2026_03",
                "x-total-results": "3",
                "link": f'<{_BASE_URL}/genecentric/search?cursor=next>; rel="next"',
            },
        )

    request = _explicit_request(tmp_path).model_copy(
        update={"acquisition": CanonicalGeneDownloadDocument()}
    )
    client, transport = _transport(httpx.MockTransport(handler))
    try:
        execution = run_uniprot_download(
            resolve_uniprot_download(request, request_base=tmp_path),
            transport=transport,
        )
    finally:
        client.close()

    content = execution.fasta_path.read_text(encoding="utf-8")
    assert content.count(">") == 2
    assert ">sp|P1|ONE_HUMAN" in content
    assert ">tr|P2|TWO_HUMAN" in content
    assert execution.document.observed_releases == ("2026_03", "2026_04")
    assert execution.document.provider_reported_counts == (3,)
    assert len(execution.document.warnings) == 2


def test_empty_download_leaves_effective_request_but_no_success_artifacts(tmp_path: Path) -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path.startswith("/proteomes/"):
            return httpx.Response(200, json=_HUMAN_PROTEOME)
        return httpx.Response(200, content=b"", headers={"x-total-results": "0"})

    effective = resolve_uniprot_download(_explicit_request(tmp_path), request_base=tmp_path)
    client, transport = _transport(httpx.MockTransport(handler))
    try:
        with pytest.raises(ValueError, match="no FASTA entries"):
            run_uniprot_download(effective, transport=transport)
    finally:
        client.close()

    assert not effective.output_fasta.exists()
    assert not effective.output_fasta.with_suffix(".fasta.result.json").exists()
    assert effective.output_fasta.with_suffix(".fasta.effective.json").is_file()
    assert not list(tmp_path.glob(".*.tmp"))


def test_catalog_sync_is_versioned_filterable_and_manifest_committed(tmp_path: Path) -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.params.get("cursor") == "next":
            return httpx.Response(
                200,
                json={"results": [_MOUSE_PROTEOME]},
                headers={"x-uniprot-release": "2026_03", "x-total-results": "2"},
            )
        return httpx.Response(
            200,
            json={"results": [_HUMAN_PROTEOME]},
            headers={
                "x-uniprot-release": "2026_03",
                "x-total-results": "2",
                "link": f'<{_BASE_URL}/proteomes/search?cursor=next>; rel="next"',
            },
        )

    request = UniProtCatalogRequestDocument(output_dir=Path("catalog"))
    moment = datetime.datetime(2026, 8, 28, 12, 0, tzinfo=datetime.UTC)
    client, transport = _transport(httpx.MockTransport(handler))
    try:
        execution = sync_uniprot_catalog(
            request,
            request_base=tmp_path,
            transport=transport,
            retrieved_at=moment,
        )
    finally:
        client.close()

    assert execution.catalog_path.name == "uniprot-proteomes-20260828T120000000000Z.parquet"
    assert execution.result_path.name == "uniprot-catalog-20260828T120000000000Z.result.json"
    assert execution.document.artifact.checksum == file_checksum(execution.catalog_path)
    assert execution.document.artifact.row_count == 2
    assert execution.document.warnings == ()
    assert latest_uniprot_catalog(execution.catalog_path.parent) == execution.catalog_path
    frame = read_uniprot_catalog(execution.catalog_path)
    assert filter_uniprot_catalog(frame, "homo").get_column("proteome_id").to_list() == [
        "UP000005640"
    ]
    assert filter_uniprot_catalog(frame, "10090").get_column("proteome_id").to_list() == [
        "UP000000589"
    ]


def test_failed_catalog_refresh_preserves_previous_readable_snapshot(tmp_path: Path) -> None:
    def success(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={"results": [_HUMAN_PROTEOME]})

    request = UniProtCatalogRequestDocument(output_dir=Path("catalog"))
    first_moment = datetime.datetime(2026, 8, 28, 12, 0, tzinfo=datetime.UTC)
    client, transport = _transport(httpx.MockTransport(success))
    try:
        first = sync_uniprot_catalog(
            request,
            request_base=tmp_path,
            transport=transport,
            retrieved_at=first_moment,
        )
    finally:
        client.close()

    def failure(request: httpx.Request) -> httpx.Response:
        raise httpx.ConnectError("connection reset", request=request)

    client, transport = _transport(httpx.MockTransport(failure))
    try:
        with pytest.raises(httpx.ConnectError, match="connection reset"):
            sync_uniprot_catalog(
                request,
                request_base=tmp_path,
                transport=transport,
                retrieved_at=first_moment + datetime.timedelta(seconds=1),
            )
    finally:
        client.close()

    assert latest_uniprot_catalog(first.catalog_path.parent) == first.catalog_path
    assert len(list(first.catalog_path.parent.glob("*.parquet"))) == 1
