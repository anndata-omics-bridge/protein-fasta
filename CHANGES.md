# Changes

- 2026-08-28: Added `table --checksums` with per-row normalized `sequence_hash` and versioned
  `id_sequence_fingerprint` columns, and demonstrated the concise export in the executable
  walkthrough.

- 2026-08-28: Exposed the amino-acid counts and percentages already computed by the streaming FASTA
  summary in the `diagnostics` CLI report and executable walkthrough.

- 2026-08-27: Changed automatic and configured frame parsing to enrich each row accepted by exactly
  one format while preserving unmatched rows with null parser fields. Added separately named strict
  Python readers and `--strict` CLI options for the whole-file-or-base-schema contract. Expanded the
  executable walkthrough with mixed UniProtKB/RefSeq and stacked `REV_`/`CON__` examples.

- 2026-08-27: Moved the documentation build to Zensical and added an executable shell walkthrough
  that runs every Cyclopts command against deterministic FASTA fixtures, captures its output, and
  fails the docs build on CLI drift without requiring Quarto or Jupyter.

- 2026-08-27: Moved reusable FASTA database construction and registry ownership from `fasta_gen`:
  configured tryptic digestion, BLAKE2b-128 protein/peptide hashes, non-security MD5 file
  checksums, naming, `aa|` metadata and contaminant markers, decoy/entrapment generation,
  deterministic assembly, registry schema 11, SQLite/DuckDB backends, indexing, materialized pair
  metrics, comparison fingerprints, exports, and clustering. Added Pydantic build/registry JSON
  documents and the single-word Cyclopts commands `digest`, `checksum`, `build`, `index`,
  `registry`, `compare`, `pairs`, and `cluster`. Build commands emit a reproducibility manifest.

- 2026-08-27: Added the single-word Cyclopts commands `basic`, `table`, `configured`, `formats`, and
  `diagnostics` for table export and reporting. CSV, TSV, XLSX, and Parquet protein tables support
  `--no-sequence`. Added
  `ProteinDiagnosticsSummary`, its streaming aggregation function, and an aggregate diagnostic CLI
  report. Expanded the maintained API, configuration, CLI, architecture, and benchmark
  documentation to cover the complete 0.2 design.

- 2026-08-27: Limited the packaged entry-classifier rules to broadly used decoy and contaminant
  decorations. FGCZ-specific `aa|`, `sp|Cont_`, `zh|C...`, and `_p_target` conventions remain
  solely in `fasta_gen`'s application JSON. Added direct UniProtKB/RefSeq non-overlap and mixed-file
  fallback coverage, and made wheel smoke tests bypass stale isolated-environment caches so they
  always exercise the artifact just built.

- 2026-08-26: Replaced the prototype header/entry-kind API with minimal `ProteinRecord` and
  compositional `ProteinDiagnostics` products; made normalization fixed; moved namespace,
  classification, and extraction policy into strict versioned JSON documents; added exact base and
  homogeneous UniProtKB/RefSeq Polars frames, generic writing, streaming summaries, packaged JSON
  Schemas, and directed import contracts.
- 2026-08-25: Created the shared protein-FASTA package with streaming compressed input, header
  interpretation, identifier classification, pattern resolution, and sequence validation.
