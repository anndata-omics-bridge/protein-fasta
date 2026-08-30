# Changes

- 2026-08-30: Kept the credential-free public type-check aligned with the 0.3 runtime-module
  locations, so public CI excludes only the two private advanced-generation adapters.

- 2026-08-30: Bumped the package to 0.3.0 for the breaking removal of the legacy combined-build
  types and `build/` import paths.

- 2026-08-30: Removed the peptide runtime's dependency on registry storage. Spill executors now
  exercise a peptide-owned partition-workspace capability; the root peptide workflow adapts the
  existing SQLite or DuckDB registry connection to that capability.

- 2026-08-30: Began the final database-runtime cutover. Decoy and entrapment strategy documents now
  compile once at root boundaries into schema-independent runtime behavior; database naming and
  metadata have explicit runtime values; biological build results no longer contain decoy fields;
  and in-memory target, contaminant, and foreign sources use the same canonical protein-input
  preparation as file requests. Removed the superseded `build/` implementation tree.

- 2026-08-30: Made request-driven Cyclopts workflows direct-first. Common CLI arguments now author
  a validated non-overwriting `*.request.json` before computation; `--request` explicitly replays
  it and `--output` redirects replay products. Added explicit diagnostic and enzyme-rule CLI paths,
  request/replay coverage across all workflows, and contextual executable documentation that
  introduces and explains each JSON beside its consuming command.

- 2026-08-28: Added a credential-free public Quality gate that checks every public capability while
  leaving only the private `fdr_benchmark` algorithm adapters to the full generation gate. Moved
  entrapment peptide-pair evidence and stable TSV serialization behind a dependency-free
  `protein_fasta` runtime type so biological builds and downstream consumers no longer leak the
  optional package's data model.

- 2026-08-28: Completed the staged protein and peptide workflow API. Added immutable UniProt
  catalog snapshots and all three acquisition modes; canonical protein-input, biological,
  search, peptide, mapping, comparison, entrapment-pair, candidate, and decoy-report artifacts;
  SQLite/DuckDB inventory indexing; and matching Cyclopts commands. Biological `build` and
  subsequent `decoy` are independent replayable operations. Added a validated in-memory Polars
  build adapter that persists its replay input and delegates to the artifact build path. Moved
  the reusable UniProt, candidate, peptide, decoy-diagnostic, and derived-entrapment workflows
  out of `fasta_gen`. Reworked the documentation around portrait use-case workflows and an
  executable 50-entry UniProtKB/RefSeq/mixed/contaminant walkthrough.

- 2026-08-28: Made empty registry neighbourhoods portable: an empty database-id selection now
  compiles to a false predicate instead of `IN ()`, which SQLite accepted but DuckDB rejected.
  Added cross-backend coverage for the empty selection.

- 2026-08-28: Declared `license-files`, without which uv_build shipped wheels carrying no licence
  text even though the repository had a `LICENSE` file.

- 2026-08-28: Added the build workflow API: packaged or explicit profile JSON resolves with a
  per-run request and typed CLI overrides into a replayable effective request, then
  `run_database_build()` writes the effective JSON before source reading and a typed result JSON
  after completion. The result records relative checksummed artifacts, complete length/amino-acid
  summaries, reconciled counts, normalization, and biological-generation evidence.
  `protein-fasta build` produces only the biological database; `protein-fasta decoy` consumes
  that database's inventory afterward. Every workflow build also writes a final-order
  protein inventory Parquet with raw headers, normalized fields, operational kinds, contaminant
  blocks, sequence hashes, and generation labels. Added portrait workflow and class diagrams
  before the implementation documentation.

- 2026-08-28: Corrected the recursive `build/` ignore rule that excluded the database-build
  implementation and its tests from Git. Moved the tests to the discoverable `tests/database_build/`
  capability folder, added direct shuffled/foreign entrapment assembly coverage, and closed the raw
  SQLite connection used by the old-schema test. Kept ordinary reverse-decoy builds in the base/CLI
  installation while isolating shuffle, DecoyPYrat, and entrapment adapters behind the explicit
  `generation` extra with actionable missing-extra errors.

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
