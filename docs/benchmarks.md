# Benchmarks

Measurements are descriptive, not pass/fail thresholds.

## Current record/frame agreement

Measured 2026-08-27 on Apple Silicon with Python 3.13.9. The input
`p39466_db2_uniprotkb_UP000000532_20250820.fasta` is 1,216,714 bytes and contains 2,626 records.
`iter_proteins()` and `read_basic_protein_frame()` agreed row-for-row on every `(id, description,
sequence)` value.

The stored file contains three FGCZ `aa|` sentinel records. RefSeq matched 0/2,626 rows and
UniProtKB matched 2,623/2,626, so automatic whole-file detection correctly returned the exact base
frame. A temporary view excluding only those three sentinels selected UniProtKB for all 2,623
remaining rows.

Each result below is the median of seven fresh measured processes after one discarded warm-up.
Temporary filtering and runtime preparation were outside their respective timed regions.

| Operation | Median |
| --- | ---: |
| `iter_proteins` | 6.502 ms |
| `iter_protein_diagnostics` with precompiled rules | 11.461 ms |
| `read_basic_protein_frame` | 6.481 ms |
| Format selection on a prepared complete-file frame | 0.532 ms |
| Homogeneous UniProtKB `read_protein_frame` | 12.907 ms |

## Initial lexical-reader comparison

Initial standalone sanity check on 2026-08-25, Apple M4 Pro, Python 3.13.9. The input was APB's
15,215,445-byte `ProteoBenchFASTA_DDAQuantification_noecoli.fasta`, containing 27,488 records. Each
reader streamed and counted the same records five times in a fresh Python process.

| Reader | Median | Five runs (seconds) |
| --- | ---: | --- |
| `protein_fasta.reading.parser.read_records` | 0.0217 s | 0.0345, 0.0211, 0.0227, 0.0217, 0.0216 |
| Current APB reader | 0.0168 s | 0.0292, 0.0167, 0.0169, 0.0168, 0.0155 |
| Current `fasta_gen` reader | 0.0169 s | 0.0169, 0.0171, 0.0169, 0.0162, 0.0165 |

The shared reader is about 5 ms slower over this 15 MB file. That is accepted for the reviewed API:
it returns named immutable `FastaRecord` values and applies the strict shared whitespace/error
contract. This is a measurement, not a speed gate. Consumer integration should re-run it on large
plain and compressed inputs and profile before changing the record contract.

The reproducible counting function is
[`benchmarks/benchmark_reading.py`](https://github.com/anndata-omics-bridge/protein-fasta/blob/main/benchmarks/benchmark_reading.py).
