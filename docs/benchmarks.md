# Reader benchmark

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
