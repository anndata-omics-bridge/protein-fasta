# Workflow artifacts

Workflow stages exchange Parquet tables and FASTA files. Parameter and result documents are JSON
serializations of the same typed Pydantic values accepted by the Python API. CSV, TSV, and XLSX
are inspection or compatibility exports, not canonical handoffs.

## Common artifact evidence

Every result document describes each owned file with the same fields:

| Field | Meaning |
| --- | --- |
| `schema_name`, `schema_version` | Semantic artifact contract, independent of its filename |
| `path` | Path recorded relative to the result location where possible |
| `checksum_version`, `checksum` | Exact-file checksum algorithm and digest |
| `byte_count` | Exact published byte count |
| `row_count` | Table or FASTA record count when meaningful |

The effective parameter document is written before computation. Data files are staged and
published atomically; for a multi-file operation the result document is published last and acts
as the success marker. A failed operation does not write success evidence.

## Protein-input Parquet

`prepare` and `derive-input` produce the only canonical file input to biological construction.

| Column | Type | Meaning |
| --- | --- | --- |
| `source_order`, `record_order` | Int64 | Stable ordering coordinates |
| `source_id` | String | Caller-authored source identity |
| `role` | String | `target`, `contaminant`, or `foreign` |
| `block_name`, `block_description` | nullable String | Contaminant-block identity and label |
| `raw_header`, `id`, `description` | String / nullable String | Exact header and parsed base fields |
| `sequence` | String | Normalized sequence used by later computation |
| `upper_cased`, `terminal_stop_stripped` | Boolean | Per-record normalization audit |

## Biological and search inventories

Biological build projects one immutable `BiologicalDatabase` to both FASTA and
`protein-inventory.parquet`. Its inventory columns are:

| Column | Type | Meaning |
| --- | --- | --- |
| `final_order` | Int64 | Exact FASTA/inventory order |
| `source_order`, `record_order` | nullable Int64 | Prepared-source coordinates |
| `source_id`, `source_role` | nullable String | Prepared-source provenance |
| `raw_header`, `id`, `description`, `sequence` | String / nullable String | Written protein record |
| `kind` | String | `sentinel`, `target`, `contaminant`, or `entrapment` |
| `contaminant_group` | nullable String | Marker-derived contaminant block |
| `sequence_hash` | String | BLAKE2b-128 sequence digest in hexadecimal |
| `entrapment_strategy` | nullable String | Generation identity for entrapment rows |

`search-inventory.parquet` contains those columns unchanged and adds nullable
`decoy_strategy`, `decoy_source_order`, `decoy_source_id`, and `decoy_source_kind`. Those four
fields are required for a decoy runtime row and absent for biological rows. Decoy generation
rejects an input inventory that already contains decoys.

## Entrapment pairs

Shuffled entrapment can additionally produce `entrapment-pairs.parquet` with `source_id`,
`target_peptide`, `generated_peptide`, and `fold_index`. Foreign-species entrapment has no
target-to-generated pairing and therefore does not publish an empty or misleading pair artifact.

## Peptide artifacts

`peptides.parquet` contains one row per distinct peptide:

| Column | Type | Meaning |
| --- | --- | --- |
| `peptide_id`, `sequence` | String | Stable peptide hash identity and sequence |
| `length`, `missed_cleavages` | Int64 | Digestion evidence |
| `mapping_count`, `protein_count` | Int64 | Mapping and distinct-protein multiplicity |
| `target_count`, `contaminant_count`, `entrapment_count`, `decoy_count` | Int64 | Kind-specific protein multiplicity |

`protein-peptide-map.parquet` contains `peptide_id`, `peptide_sequence`, `protein_order`,
`protein_id`, `protein_kind`, and `missed_cleavages`. Repeated identical mappings are deduplicated;
conflicting facts for the same mapping are rejected.

`peptide-comparisons.parquet` contains one row for each of `all`, `target`, `contaminant`,
`entrapment`, and `decoy`, with distinct counts, shared count, Jaccard, both directional coverages,
and containment.

## Candidate and decoy diagnostics

Candidate comparison Parquet contains one row per registry database and scientific kind. It
records selected/other/shared ID and sequence counts, directional coverage, containment, Jaccard,
description and exact-pair overlap, changed shared IDs, one-sided IDs, and exact-set/content flags.
The result document separately records candidate/registry checksums, availability counts, and the
bounded-neighbour clustering.

Decoy-method comparison Parquet records each requested method's protein and peptide counts,
unique ratio, peptide-length quantiles, target sharing, repeated peptides, composition overlap,
initial and unresolved collisions, dropped peptides, and omitted decoys.

## UniProt artifacts

Catalog synchronization writes immutable timestamped Parquet with `proteome_id`, `taxid`,
`organism`, `proteome_type`, `swissprot`, `swissprot_trembl`, and `one_seq_per_gene`. The paired
result records the provider query, retrieval time, releases, reported counts, actual row count,
warnings, and exact artifact checksum. Only result-committed checksum-valid snapshots are returned
by catalog discovery.

UniProt download writes a provider or reconstructed FASTA plus an effective request and result.
The result records the resolved proteome, resolution method and query, download query, every
observed release, actual and provider-reported entry counts, checksum, size, and warnings.

## Result ownership

| Workflow | Result owns |
| --- | --- |
| source preparation | source FASTA evidence and protein-input Parquet |
| biological build | input Parquet, biological FASTA/inventory, effective parameters, optional entrapment pairs, normalization and complete AA/length summary |
| decoy generation | biological inventory input, search FASTA/inventory, counts, complete AA/length summary, and algorithm evidence |
| peptide build | protein/search inventory input, peptide and mapping Parquet, unique-peptide FASTA, counts, and digestion/executor parameters |
| peptide comparison | both peptide inputs and comparison Parquet |
| candidate review | candidate inventory, registry snapshot, comparison Parquet, counts, and neighbourhood |
| UniProt catalog/download | provider evidence and the published Parquet or FASTA |

Committed JSON Schemas under `protein_fasta/documents/_schema/` are the machine-readable
contracts for every parameter and result document. Canonical Parquet readers compare the complete
column order and data types and reject schema drift before constructing runtime values.
