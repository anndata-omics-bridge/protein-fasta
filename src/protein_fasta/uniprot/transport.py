"""Thin HTTP transport for the UniProt REST provider."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from types import TracebackType
from typing import Self, cast

import httpx

from protein_fasta.uniprot.models import ProviderPage, ProviderTransferEvidence

BASE_URL = "https://rest.uniprot.org"


class UniProtTransport:
    """Own or borrow one synchronous HTTP client for UniProt REST calls."""

    def __init__(
        self,
        *,
        base_url: str = BASE_URL,
        timeout_seconds: float = 120.0,
        client: httpx.Client | None = None,
    ) -> None:
        self._base_url = base_url.rstrip("/")
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=timeout_seconds,
            follow_redirects=True,
        )

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        """Close the HTTP client only when this transport created it."""
        if self._owns_client:
            self._client.close()

    def proteome(self, proteome_id: str, /) -> dict[str, object]:
        """Return one explicit proteome provider object."""
        response = self._client.get(f"{self._base_url}/proteomes/{proteome_id}")
        response.raise_for_status()
        return _json_object(response)

    def first_proteome(self, query: str, /) -> dict[str, object] | None:
        """Return the first proteome matching a provider query."""
        response = self._client.get(
            f"{self._base_url}/proteomes/search",
            params={"query": query, "format": "json", "size": "1"},
        )
        response.raise_for_status()
        records = _result_objects(response)
        return records[0] if records else None

    def iter_proteome_pages(self, query: str, /) -> Iterator[ProviderPage]:
        """Yield all pages from the UniProt proteomes endpoint."""
        yield from self._iter_pages(
            f"{self._base_url}/proteomes/search",
            {"query": query, "format": "json", "size": "500"},
        )

    def iter_genecentric_pages(self, query: str, /) -> Iterator[ProviderPage]:
        """Yield all pages from the UniProt GeneCentric endpoint."""
        yield from self._iter_pages(
            f"{self._base_url}/genecentric/search",
            {"query": query, "format": "json", "size": "500"},
        )

    def stream_uniprotkb_fasta(
        self,
        query: str,
        destination: Path,
        /,
    ) -> ProviderTransferEvidence:
        """Stream one UniProtKB FASTA response into an unpublished file."""
        actual_count = 0
        releases: list[str] = []
        reported_counts: list[int] = []
        with self._client.stream(
            "GET",
            f"{self._base_url}/uniprotkb/stream",
            params={"query": query, "format": "fasta"},
        ) as response:
            response.raise_for_status()
            _append_header_evidence(response, releases, reported_counts)
            with destination.open("w", encoding="utf-8", newline="\n") as handle:
                for line in response.iter_lines():
                    handle.write(line)
                    handle.write("\n")
                    if line.startswith(">"):
                        actual_count += 1
        return ProviderTransferEvidence(
            actual_entry_count=actual_count,
            observed_releases=_distinct(releases),
            provider_reported_counts=_distinct_int(reported_counts),
        )

    def _iter_pages(
        self,
        initial_url: str,
        initial_params: dict[str, str],
    ) -> Iterator[ProviderPage]:
        url: str | None = initial_url
        params: dict[str, str] | None = initial_params
        while url is not None:
            response = self._client.get(url, params=params)
            response.raise_for_status()
            yield ProviderPage(
                records=tuple(_result_objects(response)),
                release=response.headers.get("x-uniprot-release"),
                reported_count=_reported_count(response),
            )
            url = response.links.get("next", {}).get("url")
            params = None


def _json_object(response: httpx.Response) -> dict[str, object]:
    value = cast(object, response.json())
    if not isinstance(value, dict):
        raise ValueError(f"UniProt returned a non-object JSON response from {response.url}")
    mapping = cast("dict[object, object]", value)
    return {str(key): item for key, item in mapping.items()}


def _result_objects(response: httpx.Response) -> list[dict[str, object]]:
    payload = _json_object(response)
    value = payload.get("results", [])
    if not isinstance(value, list):
        raise ValueError(f"UniProt returned a non-list results field from {response.url}")
    records: list[dict[str, object]] = []
    for item in cast("list[object]", value):
        if not isinstance(item, dict):
            raise ValueError(f"UniProt returned a non-object result from {response.url}")
        mapping = cast("dict[object, object]", item)
        records.append({str(key): field for key, field in mapping.items()})
    return records


def _reported_count(response: httpx.Response) -> int | None:
    value = response.headers.get("x-total-results")
    if value is None:
        return None
    try:
        count = int(value)
    except ValueError as error:
        raise ValueError(f"UniProt returned invalid x-total-results {value!r}") from error
    if count < 0:
        raise ValueError(f"UniProt returned negative x-total-results {count}")
    return count


def _append_header_evidence(
    response: httpx.Response,
    releases: list[str],
    reported_counts: list[int],
) -> None:
    release = response.headers.get("x-uniprot-release")
    if release is not None:
        releases.append(release)
    reported_count = _reported_count(response)
    if reported_count is not None:
        reported_counts.append(reported_count)


def _distinct(values: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(values))


def _distinct_int(values: list[int]) -> tuple[int, ...]:
    return tuple(dict.fromkeys(values))
