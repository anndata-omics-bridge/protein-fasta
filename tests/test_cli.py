"""Tests for the installed FASTA table-export command."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile

import polars as pl
import pytest

from protein_fasta.analytics.hashing import content_fingerprint, sequence_hash
from protein_fasta.cli import (
    basic,
    build,
    cluster,
    compare,
    digest,
    formats,
    index,
    pairs,
    registry,
    table,
)

_CLI = str(Path(sys.executable).parent / "protein-fasta")


def _write_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "proteins.fasta"
    path.write_text(">P1 Protein one\nac d*\n>P2\nEFG\n")
    return path


def _write_mixed_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "mixed.fasta"
    path.write_text(
        ">sp|P12345|RL40_YEAST Protein OS=Yeast OX=1 PE=1 SV=1\nAA\n"
        ">ref|NP_123456.1| ATP synthase subunit [Homo sapiens]\nBB\n"
    )
    return path


def test_export_table_writes_csv_with_sequence_by_default(tmp_path: Path) -> None:
    output = tmp_path / "proteins.csv"

    table(_write_fasta(tmp_path), output)

    assert output.read_text() == ("id,description,sequence\nP1,Protein one,ACD\nP2,,EFG\n")


def test_export_basic_table_writes_only_base_columns(tmp_path: Path) -> None:
    output = tmp_path / "proteins.csv"

    basic(_write_fasta(tmp_path), output, sequence=False)

    assert output.read_text() == "id,description\nP1,Protein one\nP2,\n"


def test_installed_cli_can_exclude_sequence(tmp_path: Path) -> None:
    output = tmp_path / "proteins.csv"

    result = subprocess.run(
        [_CLI, "table", str(_write_fasta(tmp_path)), str(output), "--no-sequence"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output.read_text() == "id,description\nP1,Protein one\nP2,\n"
    assert "Wrote 2 FASTA records and 2 columns" in result.stderr


def test_export_table_strict_mode_returns_base_schema_for_mixed_input(
    tmp_path: Path,
) -> None:
    output = tmp_path / "proteins.csv"

    table(_write_mixed_fasta(tmp_path), output, strict=True)

    assert pl.read_csv(output).columns == ["id", "description", "sequence"]


def test_export_table_adds_per_row_checksums_when_requested(tmp_path: Path) -> None:
    output = tmp_path / "proteins.csv"

    table(_write_fasta(tmp_path), output, checksums=True)

    frame = pl.read_csv(output)
    first_hash = sequence_hash("ACD")
    second_hash = sequence_hash("EFG")
    assert frame["id"].to_list() == ["P1", "P2"]
    assert frame["sequence"].to_list() == ["ACD", "EFG"]
    assert frame["sequence_hash"].to_list() == [first_hash.hex(), second_hash.hex()]
    assert frame["id_sequence_fingerprint"].to_list() == [
        content_fingerprint((("P1", first_hash),)),
        content_fingerprint((("P2", second_hash),)),
    ]


def test_export_table_writes_tsv(tmp_path: Path) -> None:
    output = tmp_path / "proteins.tsv"

    table(_write_fasta(tmp_path), output)

    assert output.read_text().splitlines()[0] == "id\tdescription\tsequence"


def test_export_table_writes_parquet(tmp_path: Path) -> None:
    output = tmp_path / "proteins.parquet"

    table(_write_fasta(tmp_path), output)

    assert pl.read_parquet(output).to_dicts() == [
        {"id": "P1", "description": "Protein one", "sequence": "ACD"},
        {"id": "P2", "description": None, "sequence": "EFG"},
    ]


def test_export_table_writes_xlsx(tmp_path: Path) -> None:
    output = tmp_path / "proteins.xlsx"

    table(_write_fasta(tmp_path), output)

    with ZipFile(output) as workbook:
        assert "xl/worksheets/sheet1.xml" in workbook.namelist()


def test_export_table_rejects_unknown_suffix(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"unsupported table suffix '.txt'"):
        table(_write_fasta(tmp_path), tmp_path / "proteins.txt")


def test_export_format_diagnostics_writes_one_row_per_builtin_format(tmp_path: Path) -> None:
    output = tmp_path / "formats.csv"

    formats(_write_fasta(tmp_path), output)

    assert pl.read_csv(output).to_dicts() == [
        {"format": "refseq", "matched_rows": 0, "total_rows": 2, "status": "no_match"},
        {"format": "uniprotkb", "matched_rows": 0, "total_rows": 2, "status": "no_match"},
    ]


def test_installed_cli_summarizes_record_diagnostics(tmp_path: Path) -> None:
    result = subprocess.run(
        [_CLI, "diagnostics", str(_write_fasta(tmp_path))],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Records: 2; total residues: 6" in result.stderr
    assert (
        "Amino-acid frequencies: A=1 (16.67%), C=1 (16.67%), D=1 (16.67%), "
        "E=1 (16.67%), F=1 (16.67%), G=1 (16.67%)"
    ) in result.stderr
    assert "Identifier namespaces: unmatched=2" in result.stderr
    assert "Normalization changes: upper-cased=1; terminal-stop-stripped=1" in result.stderr


def test_installed_cli_accepts_explicit_format_documents(tmp_path: Path) -> None:
    output = tmp_path / "configured.csv"
    documents = Path(__file__).parents[1] / "src" / "protein_fasta" / "documents"

    result = subprocess.run(
        [
            _CLI,
            "configured",
            str(_write_fasta(tmp_path)),
            str(output),
            "--rules",
            str(documents / "frame_formats" / "uniprotkb" / "rules.json"),
            "--rules",
            str(documents / "frame_formats" / "refseq" / "rules.json"),
            "--classifiers",
            str(documents / "entry_classifiers" / "rules.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output.read_text() == ("id,description,sequence\nP1,Protein one,ACD\nP2,,EFG\n")


def test_digest_writes_peptides_with_missed_cleavage_evidence(tmp_path: Path) -> None:
    fasta = tmp_path / "digest.fasta"
    fasta.write_text(">P1\nAKBK\n")
    config = tmp_path / "digestion.json"
    config.write_text('{"enzyme":"trypsin","min_length":1,"max_length":50,"missed_cleavages":1}')
    output = tmp_path / "peptides.csv"

    digest(fasta, output, config=config)

    assert pl.read_csv(output).to_dicts() == [
        {"protein_id": "P1", "peptide": "AK", "missed_cleavages": 0},
        {"protein_id": "P1", "peptide": "AKBK", "missed_cleavages": 1},
        {"protein_id": "P1", "peptide": "BK", "missed_cleavages": 0},
    ]


def test_installed_cli_reports_versioned_file_and_content_checksums(tmp_path: Path) -> None:
    result = subprocess.run(
        [_CLI, "checksum", str(_write_fasta(tmp_path))],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "md5-file-v1" in result.stderr
    assert "blake2b-128-id-sequence-pairs-v1" in result.stderr


def _write_registry_fastas(tmp_path: Path) -> Path:
    root = tmp_path / "fastas"
    root.mkdir()
    (root / "alpha_20260101.fasta").write_text(">P1 one\nAAAA\n>P2 two\nBBBB\n")
    (root / "beta_20260101.fasta").write_text(">P1 one\nAAAA\n>P3 three\nCCCC\n")
    return root


def test_registry_cli_functions_cover_index_compare_pairs_and_cluster(tmp_path: Path) -> None:
    registry_path = tmp_path / "registry.sqlite3"
    index(_write_registry_fastas(tmp_path), registry_path)

    registry_table = tmp_path / "registry.csv"
    comparison_table = tmp_path / "comparison.csv"
    pair_table = tmp_path / "pairs.tsv"
    id_matrix = tmp_path / "ids.tsv"
    cluster_table = tmp_path / "cluster.csv"
    registry(registry_path, registry_table)
    compare(registry_path, 1, comparison_table)
    pairs(registry_path, pair_table, ids=id_matrix)
    cluster(registry_path, cluster_table)

    registry_rows = pl.read_csv(registry_table).to_dicts()
    assert [row["filename"] for row in registry_rows] == [
        "alpha_20260101.fasta",
        "beta_20260101.fasta",
    ]
    assert all(
        str(row["target_content_fingerprint"]).startswith("blake2b-128:") for row in registry_rows
    )
    assert pl.read_csv(comparison_table).to_dicts()[0]["shared_ids"] == 1
    assert len(pl.read_csv(pair_table, separator="\t")) == 1
    assert pl.read_csv(id_matrix, separator="\t").shape == (2, 3)
    assert pl.read_csv(cluster_table).to_dicts()[0]["leaf_count"] == 2


def test_index_accepts_json_registry_policy(tmp_path: Path) -> None:
    registry_path = tmp_path / "registry.duckdb"
    config = tmp_path / "registry.json"
    config.write_text(
        '{"schema_version":"0.1","registry":{"backend":"duckdb"},"max_detailed_entries":100}'
    )

    index(_write_registry_fastas(tmp_path), registry_path, config=config)

    assert registry_path.is_file()


def test_build_resolves_profile_request_and_writes_typed_result(tmp_path: Path) -> None:
    target = tmp_path / "target.fasta"
    target.write_text(">P1 one\nac d*\n")
    config = tmp_path / "build.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": "0.2",
                "targets": ["target.fasta"],
                "output_dir": "out",
                "date": "2026-08-27",
                "name_fields": {"project": 1, "dbn": 2, "description": "demo"},
                "decoy": None,
            }
        )
    )

    build(config)

    fasta = tmp_path / "out" / "p1_db2_demo_20260827.fasta"
    effective = fasta.with_suffix(".fasta.effective.json")
    result = fasta.with_suffix(".fasta.result.json")
    assert fasta.read_text().endswith(">P1 one\nACD\n")
    effective_payload = json.loads(effective.read_text())
    payload = json.loads(result.read_text())
    assert effective_payload["decoy"] is None
    assert effective_payload["targets"] == [str(target)]
    assert payload["counts"] == {
        "contaminant": 0,
        "decoy": 0,
        "entrapment": 0,
        "target": 1,
        "total": 2,
    }
    fasta_artifact = next(
        artifact for artifact in payload["artifacts"] if artifact["schema_name"] == "protein-fasta"
    )
    assert fasta_artifact["checksum_version"] == "md5-file-v1"
    assert payload["summary"]["aa_counts"] == {"A": 4, "C": 4, "D": 1, "P": 3, "R": 3}
