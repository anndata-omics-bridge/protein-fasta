"""Stable tabular exports of materialized registry similarity statistics."""

from __future__ import annotations

import os
import tempfile
import unicodedata
from collections.abc import Collection
from itertools import combinations
from pathlib import Path
from typing import Protocol, cast

import polars as pl

from protein_fasta.analytics.comparisons import (
    KindCounts,
    KindSimilarity,
    SimilarityData,
    SimilarityPair,
)
from protein_fasta.registry.backend.base import RegistryConnection, RegistryIntegrityError
from protein_fasta.registry.kinds import DetailLevel, EntryKind

_EXPORTED_KINDS = (EntryKind.TARGET, EntryKind.CONTAMINANT)

LONG_FIELDNAMES: tuple[str, ...] = (
    "database_a_relative_path",
    "database_b_relative_path",
    "target_entries_a",
    "target_entries_b",
    "target_distinct_ids_a",
    "target_distinct_ids_b",
    "target_shared_ids",
    "target_id_jaccard",
    "target_id_coverage_a_by_b",
    "target_id_coverage_b_by_a",
    "target_id_containment",
    "target_a_contained_in_b",
    "target_b_contained_in_a",
    "target_exact_id_set",
    "target_distinct_sequences_a",
    "target_distinct_sequences_b",
    "target_shared_sequences",
    "target_sequence_jaccard",
    "target_sequence_coverage_a_by_b",
    "target_sequence_coverage_b_by_a",
    "target_sequence_containment",
    "target_sequence_a_contained_in_b",
    "target_sequence_b_contained_in_a",
    "target_exact_sequence_set",
    "target_distinct_descriptions_a",
    "target_distinct_descriptions_b",
    "target_shared_descriptions",
    "target_description_jaccard",
    "target_description_coverage_a_by_b",
    "target_description_coverage_b_by_a",
    "target_description_containment",
    "target_description_a_contained_in_b",
    "target_description_b_contained_in_a",
    "target_exact_description_set",
    "target_distinct_pairs_a",
    "target_distinct_pairs_b",
    "target_shared_exact_pairs",
    "target_matching_shared_ids",
    "target_changed_shared_ids",
    "target_exact_content",
    "contaminant_entries_a",
    "contaminant_entries_b",
    "contaminant_distinct_ids_a",
    "contaminant_distinct_ids_b",
    "contaminant_shared_ids",
    "contaminant_id_jaccard",
    "contaminant_id_coverage_a_by_b",
    "contaminant_id_coverage_b_by_a",
    "contaminant_id_containment",
    "contaminant_a_contained_in_b",
    "contaminant_b_contained_in_a",
    "contaminant_distinct_sequences_a",
    "contaminant_distinct_sequences_b",
    "contaminant_shared_sequences",
    "contaminant_sequence_jaccard",
    "contaminant_sequence_coverage_a_by_b",
    "contaminant_sequence_coverage_b_by_a",
    "contaminant_sequence_containment",
    "contaminant_sequence_a_contained_in_b",
    "contaminant_sequence_b_contained_in_a",
    "contaminant_exact_sequence_set",
    "contaminant_distinct_descriptions_a",
    "contaminant_distinct_descriptions_b",
    "contaminant_shared_descriptions",
    "contaminant_description_jaccard",
    "contaminant_description_coverage_a_by_b",
    "contaminant_description_coverage_b_by_a",
    "contaminant_description_containment",
    "contaminant_description_a_contained_in_b",
    "contaminant_description_b_contained_in_a",
    "contaminant_exact_description_set",
    "contaminant_distinct_pairs_a",
    "contaminant_distinct_pairs_b",
    "contaminant_shared_exact_pairs",
    "contaminant_matching_shared_ids",
    "contaminant_changed_shared_ids",
    "contaminant_exact_id_set",
    "contaminant_exact_content",
)


class _SortableFrame(Protocol):
    def sort(self, *, by: pl.Expr) -> pl.DataFrame:
        """Sort by one fully specified Polars expression."""
        ...


class _SelectableFrame(Protocol):
    def select(self, *expressions: pl.Expr) -> pl.DataFrame:
        """Select fully specified Polars expressions."""
        ...


def _id_filter(database_ids: Collection[int] | None, column: str) -> tuple[str, tuple[int, ...]]:
    """Return an extra SQL predicate and parameters restricting ``column``."""
    if database_ids is None:
        return "", ()
    ids = tuple(sorted(database_ids))
    placeholders = ", ".join("?" for _ in ids)
    return f" AND {column} IN ({placeholders})", ids


def _kind_counts(
    connection: RegistryConnection,
    database_ids: Collection[int] | None = None,
) -> dict[tuple[int, EntryKind], KindCounts]:
    predicate, id_params = _id_filter(database_ids, "stats.database_id")
    rows = connection.execute(
        f"""
        SELECT stats.database_id, stats.kind, stats.entry_count, stats.distinct_ids,
               stats.distinct_sequences, stats.distinct_descriptions,
               stats.distinct_pairs
        FROM database_kind_stats AS stats
        JOIN databases ON databases.id = stats.database_id
        WHERE databases.detail_level = ?
          AND stats.kind IN ('target', 'contaminant')
        {predicate}
        """,
        (DetailLevel.FULL.value, *id_params),
    ).fetchall()
    return {
        (int(row[0]), EntryKind(str(row[1]))): KindCounts(
            entries=int(row[2]),
            distinct_ids=int(row[3]),
            distinct_sequences=int(row[4]),
            distinct_descriptions=int(row[5]),
            distinct_pairs=int(row[6]),
        )
        for row in rows
    }


def _pair_counts(
    connection: RegistryConnection,
    database_ids: Collection[int] | None = None,
) -> dict[tuple[int, int, EntryKind], tuple[int, int, int, int, int]]:
    low_predicate, id_params = _id_filter(database_ids, "database_id_low")
    high_predicate, _ = _id_filter(database_ids, "database_id_high")
    rows = connection.execute(
        f"""
        SELECT database_id_low, database_id_high, kind, shared_ids,
               shared_sequence_checksums, shared_descriptions,
               shared_exact_pairs, matching_shared_ids
        FROM database_pair_stats
        WHERE kind IN ('target', 'contaminant')
        {low_predicate}
        {high_predicate}
        """,
        (*id_params, *id_params),
    ).fetchall()
    return {
        (int(row[0]), int(row[1]), EntryKind(str(row[2]))): (
            int(row[3]),
            int(row[4]),
            int(row[5]),
            int(row[6]),
            int(row[7]),
        )
        for row in rows
    }


def query_similarity_data(
    connection: RegistryConnection,
    *,
    database_ids: Collection[int] | None = None,
) -> SimilarityData:
    """Read materialized target/contaminant pair statistics from SQLite.

    ``database_ids`` restricts the result to those full-detail databases, which
    bounds both the pair count and the objects built here to the subset. Pair
    work is quadratic in the number of databases, so callers that only need a
    neighbourhood must pass it rather than filter afterwards.
    """
    database_predicate, database_params = _id_filter(database_ids, "id")
    databases = [
        (int(row[0]), str(row[1]))
        for row in connection.execute(
            f"""
            SELECT id, relative_path FROM databases
            WHERE detail_level = ?
            {database_predicate}
            ORDER BY relative_path{connection.binary_collation()}, id
            """,
            (DetailLevel.FULL.value, *database_params),
        ).fetchall()
    ]
    omitted_relative_paths = tuple(
        str(row[0])
        for row in connection.execute(
            f"SELECT relative_path FROM databases WHERE detail_level = ? "
            f"ORDER BY relative_path{connection.binary_collation()}, id",
            (DetailLevel.METADATA_ONLY.value,),
        ).fetchall()
    )
    kind_counts = _kind_counts(connection, database_ids)
    pair_counts = _pair_counts(connection, database_ids)

    missing_kind_stats = [
        (relative_path, kind.value)
        for database_id, relative_path in databases
        for kind in _EXPORTED_KINDS
        if (database_id, kind) not in kind_counts
    ]
    if missing_kind_stats:
        path, kind = missing_kind_stats[0]
        raise RegistryIntegrityError(
            f"kind-statistics cache is incomplete for {path!r} ({kind}); run a full registry reindex"
        )

    pairs: list[SimilarityPair] = []
    for database_a, database_b in combinations(databases, 2):
        id_a, path_a = database_a
        id_b, path_b = database_b
        low, high = sorted((id_a, id_b))
        similarities: dict[EntryKind, KindSimilarity] = {}
        for kind in _EXPORTED_KINDS:
            shared = pair_counts.get((low, high, kind))
            if shared is None:
                raise RegistryIntegrityError(
                    f"pair-statistics cache is incomplete for {path_a!r} and "
                    f"{path_b!r} ({kind.value}); run a full registry reindex"
                )
            similarities[kind] = KindSimilarity(
                a=kind_counts[(id_a, kind)],
                b=kind_counts[(id_b, kind)],
                shared_ids=shared[0],
                shared_sequences=shared[1],
                shared_descriptions=shared[2],
                shared_exact_pairs=shared[3],
                matching_shared_ids=shared[4],
            )
        pairs.append(
            SimilarityPair(
                database_a_relative_path=path_a,
                database_b_relative_path=path_b,
                target=similarities[EntryKind.TARGET],
                contaminant=similarities[EntryKind.CONTAMINANT],
            )
        )
    return SimilarityData(
        relative_paths=tuple(relative_path for _, relative_path in databases),
        pairs=tuple(pairs),
        omitted_relative_paths=omitted_relative_paths,
    )


def _format_float(value: float) -> str:
    if value in {0.0, 1.0}:
        return f"{value:.1f}"
    return format(value, ".12g")


def _format_bool(value: bool) -> str:
    return str(value).lower()


def _paths_alias(left: Path, right: Path) -> bool:
    """Conservatively identify existing and prospective output aliases."""
    left_resolved = left.resolve()
    right_resolved = right.resolve()
    if left_resolved == right_resolved:
        return True
    left_key = unicodedata.normalize("NFC", str(left_resolved)).casefold()
    right_key = unicodedata.normalize("NFC", str(right_resolved)).casefold()
    if left_key == right_key:
        return True
    return left.exists() and right.exists() and left.samefile(right)


def _kind_row(prefix: str, similarity: KindSimilarity) -> dict[str, str | int]:
    row: dict[str, str | int] = {
        f"{prefix}_entries_a": similarity.a.entries,
        f"{prefix}_entries_b": similarity.b.entries,
        f"{prefix}_distinct_ids_a": similarity.a.distinct_ids,
        f"{prefix}_distinct_ids_b": similarity.b.distinct_ids,
        f"{prefix}_shared_ids": similarity.shared_ids,
        f"{prefix}_id_jaccard": _format_float(similarity.id_jaccard),
        f"{prefix}_id_coverage_a_by_b": _format_float(similarity.id_coverage_a_by_b),
        f"{prefix}_id_coverage_b_by_a": _format_float(similarity.id_coverage_b_by_a),
        f"{prefix}_id_containment": _format_float(similarity.id_containment),
        f"{prefix}_a_contained_in_b": _format_bool(similarity.a_contained_in_b),
        f"{prefix}_b_contained_in_a": _format_bool(similarity.b_contained_in_a),
        f"{prefix}_distinct_sequences_a": similarity.a.distinct_sequences,
        f"{prefix}_distinct_sequences_b": similarity.b.distinct_sequences,
        f"{prefix}_shared_sequences": similarity.shared_sequences,
        f"{prefix}_sequence_jaccard": _format_float(similarity.sequence_jaccard),
        f"{prefix}_sequence_coverage_a_by_b": _format_float(similarity.sequence_coverage_a_by_b),
        f"{prefix}_sequence_coverage_b_by_a": _format_float(similarity.sequence_coverage_b_by_a),
        f"{prefix}_sequence_containment": _format_float(similarity.sequence_containment),
        f"{prefix}_sequence_a_contained_in_b": _format_bool(similarity.sequence_a_contained_in_b),
        f"{prefix}_sequence_b_contained_in_a": _format_bool(similarity.sequence_b_contained_in_a),
        f"{prefix}_distinct_descriptions_a": similarity.a.distinct_descriptions,
        f"{prefix}_distinct_descriptions_b": similarity.b.distinct_descriptions,
        f"{prefix}_shared_descriptions": similarity.shared_descriptions,
        f"{prefix}_description_jaccard": _format_float(similarity.description_jaccard),
        f"{prefix}_description_coverage_a_by_b": _format_float(
            similarity.description_coverage_a_by_b
        ),
        f"{prefix}_description_coverage_b_by_a": _format_float(
            similarity.description_coverage_b_by_a
        ),
        f"{prefix}_description_containment": _format_float(similarity.description_containment),
        f"{prefix}_description_a_contained_in_b": _format_bool(
            similarity.description_a_contained_in_b
        ),
        f"{prefix}_description_b_contained_in_a": _format_bool(
            similarity.description_b_contained_in_a
        ),
        f"{prefix}_distinct_pairs_a": similarity.a.distinct_pairs,
        f"{prefix}_distinct_pairs_b": similarity.b.distinct_pairs,
        f"{prefix}_shared_exact_pairs": similarity.shared_exact_pairs,
        f"{prefix}_matching_shared_ids": similarity.matching_shared_ids,
        f"{prefix}_changed_shared_ids": similarity.changed_shared_ids,
        f"{prefix}_exact_id_set": _format_bool(similarity.exact_id_set),
        f"{prefix}_exact_sequence_set": _format_bool(similarity.exact_sequence_set),
        f"{prefix}_exact_description_set": _format_bool(similarity.exact_description_set),
        f"{prefix}_exact_content": _format_bool(similarity.exact_content),
    }
    return row


def _long_row(pair: SimilarityPair) -> dict[str, str | int]:
    return {
        "database_a_relative_path": pair.database_a_relative_path,
        "database_b_relative_path": pair.database_b_relative_path,
        **_kind_row("target", pair.target),
        **_kind_row("contaminant", pair.contaminant),
    }


def _write_long(data: SimilarityData, path: Path) -> None:
    """Write one row per unordered pair, in the order the query returned them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [_long_row(pair) for pair in data.pairs]
    frame = pl.DataFrame(
        rows,
        # Named even when empty, so a registry with fewer than two databases still
        # writes a header the reader can bind to rather than an empty file.
        schema=dict.fromkeys(LONG_FIELDNAMES, pl.String),
        orient="row" if rows else None,
    )
    _write_tsv(frame, path)


def _write_matrix(
    data: SimilarityData,
    path: Path,
    *,
    sequence: bool,
) -> None:
    """Write a symmetric Jaccard matrix in the stable path order.

    The pairs arrive as one row per unordered pair, so each is mirrored and the
    diagonal added before pivoting: a matrix is what this exports, and stating it
    as a pivot beats filling cells by hand from a lookup keyed on sorted names.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    paths = list(data.relative_paths)
    if not paths:
        _write_tsv(pl.DataFrame(schema={"relative_path": pl.String}), path)
        return

    cells = [
        (a, b, pair.target.sequence_jaccard if sequence else pair.target.id_jaccard)
        for pair in data.pairs
        for a, b in (
            (pair.database_a_relative_path, pair.database_b_relative_path),
            (pair.database_b_relative_path, pair.database_a_relative_path),
        )
    ]
    cells.extend((one, one, 1.0) for one in paths)

    pivoted = pl.DataFrame(
        cells, schema={"row": pl.String, "column": pl.String, "value": pl.Float64}, orient="row"
    ).pivot(on="column", index="row", values="value")
    ordered = cast(_SortableFrame, pivoted).sort(
        by=pl.col("row").map_elements(paths.index, return_dtype=pl.Int64)
    )
    matrix = cast(_SelectableFrame, ordered).select(
        pl.col("row").alias("relative_path"),
        # Formatted here rather than by the writer: the exports have always
        # spelled an exact 0 or 1 with one decimal and everything else to
        # twelve significant digits, and a float formatter would not.
        *(pl.col(name).map_elements(_format_float, return_dtype=pl.String) for name in paths),
    )
    _write_tsv(matrix, path)


def _write_tsv(frame: pl.DataFrame, path: Path) -> None:
    """Write one frame as the tab-separated form these exports have always used."""
    frame.write_csv(path, separator="\t", line_terminator="\n", quote_style="necessary")


def _prepare_export_destinations(exports: list[tuple[Path, bool | None]]) -> None:
    """Validate distinct destinations and create their parent directories."""
    for index, (destination, _) in enumerate(exports):
        for previous, _ in exports[:index]:
            if _paths_alias(destination, previous):
                raise ValueError("Similarity export destinations must be distinct.")
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not destination.is_file():
            raise IsADirectoryError(f"Export destination is not a regular file: {destination}")


def write_similarity_exports(
    data: SimilarityData,
    output: Path,
    *,
    id_matrix: Path | None = None,
    sequence_matrix: Path | None = None,
) -> None:
    """Write each requested TSV through a same-directory atomic replacement."""
    exports: list[tuple[Path, bool | None]] = [(output, None)]
    if id_matrix is not None:
        exports.append((id_matrix, False))
    if sequence_matrix is not None:
        exports.append((sequence_matrix, True))

    _prepare_export_destinations(exports)

    staged: list[tuple[Path, Path]] = []
    try:
        for destination, sequence in exports:
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{destination.name}.",
                suffix=".tmp",
                dir=destination.parent,
            )
            os.close(descriptor)
            temporary_path = Path(temporary_name)
            staged.append((destination, temporary_path))
            if sequence is None:
                _write_long(data, temporary_path)
            else:
                _write_matrix(data, temporary_path, sequence=sequence)
            with temporary_path.open("rb") as handle:
                os.fsync(handle.fileno())

        for destination, temporary_path in staged:
            os.replace(temporary_path, destination)
    finally:
        for _, temporary_path in staged:
            temporary_path.unlink(missing_ok=True)
