# Protein and peptide database workflows

The workflows are the design. Configuration documents, Python classes, CLI commands, and file
formats exist to make these use cases reproducible. In every workflow figure, rounded blue nodes
are computations, green rectangles are data, cylinders are durable registries, and amber diamonds
are user or application decisions. Every figure flows from top to bottom.

## Use cases

### Acquire or select protein sources

```mermaid
flowchart TB
    A["source request parameters<br/>taxon or proteome · reviewed policy<br/>catalog filter · destination"] --> B(["resolve source selection"])
    C["local UniProt proteome catalog<br/>Parquet snapshot"] --> B
    B --> D{"download or use<br/>local source?"}
    D -->|download| E(["stream one UniProt source"])
    D -->|local| F(["validate selected FASTA"])
    E --> G["source FASTA<br/>release · query · checksum<br/>entry-count evidence"]
    F --> H["validated upload or<br/>curated contaminant FASTA"]
    G --> I(["prepare ordered source rows"])
    H --> I
    I --> J["protein-input Parquet<br/>role · block · source order<br/>normalized sequence · checksum"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    classDef decision fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#4e2600
    class B,E,F,I compute
    class A,C,G,H,J data
    class D decision
```

*Figure 1. Acquire or select inputs and prepare one canonical protein source inventory.*

The curator selects a UniProt proteome, an uploaded FASTA, or application-curated sources. A
UniProt download is one acquisition operation: it does not build, add decoys, install, or index.
`prepare` then normalizes the selected target, contaminant, and optional foreign rows once and
records their order and origin in `protein-input.parquet`.

### Build a routine search database

```mermaid
flowchart TB
    A["protein-input Parquet<br/>targets · contaminant blocks<br/>source provenance"] --> B(["assemble biological database"])
    P["biological build parameters<br/>identity · date · naming<br/>metadata · optional entrapment"] --> B
    B --> C["biological FASTA<br/>protein inventory Parquet<br/>counts · summaries · checksums"]
    C --> D{"generate a<br/>search database?"}
    Q["decoy request parameters<br/>strategy · prefix · seed<br/>collision digestion · destination"] --> E(["generate source-linked decoys"])
    D -->|yes| E
    D -->|no| F(["review biological candidate"])
    E --> G["search FASTA<br/>search inventory Parquet<br/>decoy-generation evidence"]
    G --> H(["review search candidate"])
    F --> I{"accept?"}
    H --> I
    R[("existing SQLite or<br/>DuckDB registry")] --> F
    R --> H
    I -->|yes| J(["install in site collection"])
    I -->|no| K(["discard staged artifacts"])
    J --> L(["index accepted inventory"])
    L --> R

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    classDef decision fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#4e2600
    class B,E,F,H,J,K,L compute
    class A,C,P,Q,G data
    class D,I decision
```

*Figure 2. Build biological entries once, then optionally generate and review a search database.*

`protein-fasta build` ends with a decoy-free biological database. `protein-fasta decoy` is the
subsequent operation and can be rerun with reverse, shuffle, or DecoyPYrat parameters without
repeating source preparation, contaminants, or entrapment. Candidate review is read-only.
Installation remains a `fasta_gen` site mutation, followed by an explicit inventory-indexing step.

### Derive an entrapment database

```mermaid
flowchart TB
    A["registered source database<br/>protein or search inventory Parquet"] --> B(["retain original target and<br/>contaminant proteins"])
    B --> C["prior sentinel, section marker,<br/>entrapment and decoy rows excluded"]
    F["optional foreign database<br/>protein or search inventory Parquet"] --> G(["retain foreign biological proteins"])
    C --> H(["prepare derived source rows"])
    G --> H
    H --> I["derived protein-input Parquet<br/>selected and skipped counts<br/>source checksums"]
    P["entrapment request parameters<br/>strategy · fold · seed<br/>digestion · termini policy"] --> J(["generate entrapment entries"])
    I --> J
    J --> K(["assemble biological database"])
    K --> L["biological FASTA<br/>protein inventory Parquet<br/>entrapment evidence"]
    J --> M["entrapment-pairs Parquet<br/>target peptide · generated peptide<br/>source protein · fold"]
    L --> N(["optional decoy generation"])
    D["decoy request parameters<br/>strategy · collision policy<br/>destination"] --> N
    N --> O["search FASTA<br/>search inventory Parquet<br/>decoy evidence"]
    O --> Q(["optional peptide build"])
    M --> Q
    Q --> R["peptide Parquet and FASTA<br/>protein-peptide mapping Parquet<br/>digestion evidence"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    class B,G,H,J,K,N,Q compute
    class A,C,D,F,I,L,M,O,P,R data
```

*Figure 3. Derive clean biological sources, add entrapment, and only then generate decoys.*

This is the former `fasta_gen` derived-entrapment use case at its real boundary. The source filter
keeps only existing target and contaminant proteins, exactly matching the former application
behavior. A foreign source is filtered the same way and relabeled as foreign input. New
entrapment entries and pair evidence are generated during the biological build; decoys remain a
separate later operation over the completed biological inventory.

### Build and compare peptide databases

```mermaid
flowchart TB
    A["biological or search<br/>protein inventory Parquet"] --> B(["select scientific protein rows"])
    P["digestion parameters<br/>enzyme · length bounds<br/>missed cleavages"] --> C(["compile digestion once"])
    E["execution parameters<br/>memory, SQLite, or DuckDB<br/>workers · partition size"] --> D(["select execution behavior"])
    B --> F(["partition and digest proteins"])
    C --> F
    D --> F
    F --> G(["merge unique peptides<br/>and protein mappings"])
    G --> H["peptides Parquet<br/>protein-peptide mapping Parquet<br/>unique peptide FASTA"]
    H --> I{"compare with a<br/>second peptide inventory?"}
    J["second peptides Parquet"] --> K(["compute exact population overlap"])
    I -->|yes| K
    K --> L["peptide comparison Parquet<br/>target · contaminant · entrapment<br/>decoy · all populations"]
    I -->|no| M["completed peptide products"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    classDef decision fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#4e2600
    class B,C,D,F,G,K compute
    class A,E,H,J,L,M,P data
    class I decision
```

*Figure 4. Build canonical peptide artifacts and optionally compare exact peptide populations.*

The peptide workflow preserves protein order, identity, and scientific kind. Memory, SQLite, and
DuckDB are interchangeable execution behaviors behind one result contract; they do not change
the peptide or mapping bytes. The canonical handoffs are Parquet. The FASTA is a unique-peptide
search artifact, while TSV is reserved for explicit compatibility export at an application edge.

### Compare decoy methods

```mermaid
flowchart TB
    A["biological protein<br/>inventory Parquet"] --> B(["select target, contaminant<br/>and entrapment proteins"])
    P["comparison parameters<br/>methods · prefix · seeds<br/>digestion policy"] --> C(["generate each requested<br/>decoy method"])
    B --> C
    C --> D(["digest targets and decoys"])
    D --> E(["compute method diagnostics"])
    E --> F["decoy-method comparison Parquet<br/>protein and peptide counts<br/>lengths · sharing · repetition<br/>mass overlap · omissions"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    class B,C,D,E compute
    class A,F,P data
```

*Figure 5. Compare selectable decoy methods from the same biological inventory.*

`decoy-report` reuses one biological input and reports the scientific consequences of each
requested method. It does not publish a search database, mutate a registry, or combine
entrapment-overlap analysis with decoy comparison.

## Requirements derived from the use cases

1. Source acquisition, source preparation, biological assembly, decoy generation, candidate
   review, installation, indexing, peptide construction, and comparison are separate operations.
2. Biological build never generates decoys. Decoy requests always name exactly one strategy and
   consume a completed biological inventory.
3. A canonical protein-input Parquet is the only file-based input to biological build. It carries
   source order, role, contaminant block, raw header, normalized sequence, and normalization
   evidence.
4. Biological and search databases each have one canonical frozen entry tuple. FASTA and Parquet
   are projections of that tuple rather than separately accumulated products.
5. Result manifests checksum every artifact they own. Multi-file products are staged, published
   together, and committed by publishing the result manifest last.
6. Every workflow is callable through Python and the CLI. The CLI loads serialized request
   documents and invokes the same resolver and execution function as a Python caller.
7. Parquet is the tabular workflow exchange format. JSON serializes parameters and result
   evidence. FASTA is retained where a search engine or protein-source boundary requires it.
8. Storage documents are passive Pydantic models. Frozen runtime types own computation; Polars
   frames and filesystem paths remain at adapters and workflow boundaries.

## Command boundaries

| Command | One responsibility | Primary products |
| --- | --- | --- |
| `uniprot-catalog` | Synchronize a checksummed local proteome catalog | catalog Parquet + sync evidence |
| `uniprot-proteomes` | Filter a local catalog without network access | inspection table |
| `uniprot-download` | Acquire one selected UniProt source | FASTA + acquisition evidence |
| `prepare` | Normalize ordered FASTA sources and roles | protein-input Parquet |
| `derive-input` | Filter registered biological/search inventories for a new entrapment build | protein-input Parquet + skipped counts |
| `build` | Assemble targets, contaminants, and optional entrapment | biological FASTA + protein inventory Parquet |
| `decoy` | Generate one search database from a biological inventory | search FASTA + search inventory Parquet |
| `candidate` | Compare an unregistered inventory with an existing registry | candidate comparison Parquet |
| `index-inventory` | Index one accepted canonical inventory directly | SQLite or DuckDB registry update |
| `index` | Adapt a legacy FASTA directory into the registry | SQLite or DuckDB registry update |
| `peptides` | Digest one biological or search inventory | peptides and mapping Parquet + peptide FASTA |
| `pepcompare` | Compare two canonical peptide inventories | comparison Parquet |
| `decoy-report` | Compare scientific diagnostics for requested decoy methods | decoy comparison Parquet |

## Parameter precedence

```mermaid
flowchart TB
    A["biological parameter layers<br/>packaged build policy<br/>explicit project policy<br/>run identity · date · output · entrapment<br/>explicit CLI overrides"] --> D(["resolve biological request<br/>in precedence order"])
    D --> F["effective biological request<br/>resolved paths and values"]
    F --> G(["build biological database"])
    G --> H["biological result evidence"]
    H ~~~ I["independent decoy parameter layers<br/>one required strategy<br/>explicit output · prefix overrides"]
    I --> J(["resolve decoy request"])
    J --> K["effective decoy request<br/>resolved strategy and destination"]
    K --> L(["generate search database"])
    L --> M["decoy result evidence"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    class D,G,J,L compute
    class A,F,H,I,K,M data
```

*Figure 6. Resolve biological-build and decoy parameters independently.*

Biological precedence is `packaged policy < explicit policy < run parameters < explicitly
supplied CLI values`. The effective request is written before sequence work begins. Decoy
resolution is independent: neither the build profile nor the biological request contains a decoy
default.

## Programmatic composition

The serialized documents are optional at a Python boundary. `fasta_gen` or another caller can
construct the same Pydantic request objects directly, resolve them once, and pass the returned
execution artifact to the next API:

```python
from pathlib import Path

from protein_fasta.database_build import resolve_database_build, run_database_build
from protein_fasta.decoy_database import resolve_decoy_request, run_decoy_generation
from protein_fasta.documents import (
    load_builtin_database_build_profile,
    load_database_build_request,
    load_decoy_request,
)
from protein_fasta.peptide_workflow import (
    resolve_peptide_build_request,
    run_peptide_build,
)
from protein_fasta.documents import load_peptide_build_request

protein_input = Path("protein-input.parquet")
build_path = Path("biological-build.json")
build = resolve_database_build(
    load_builtin_database_build_profile(),
    load_database_build_request(build_path),
    profile_base=build_path.parent,
    request_base=build_path.parent,
)
biological = run_database_build(protein_input, build)

decoy_path = Path("reverse-decoys.json")
decoy = resolve_decoy_request(
    load_decoy_request(decoy_path),
    request_base=decoy_path.parent,
)
search = run_decoy_generation(biological.inventory_path, decoy)

peptide_path = Path("peptides.json")
peptide = resolve_peptide_build_request(
    load_peptide_build_request(peptide_path),
    request_base=peptide_path.parent,
)
peptides = run_peptide_build(search.search_inventory_path, peptide)
```

The chain passes paths when durable replay is desired. At the computation boundary, the
equivalent values are frozen `BiologicalDatabase`, `SearchDatabase`, and `PeptideDatabase`
instances containing tuples of typed entries. Polars frames are returned by artifact-reading and
report APIs; they are not passed into backend-free sequence computation.

## API and storage types

```mermaid
classDiagram
    direction TB
    class ProteinInputRequestDocument {
      ordered FASTA sources
      source roles and blocks
      output destination
    }
    class ProteinInputExecution {
      DataFrame frame
      result evidence
      protein input path
    }
    class EffectiveDatabaseBuildDocument {
      resolved biological policy
      identity and destination
      optional entrapment request
    }
    class BiologicalDatabase {
      tuple~ProteinInventoryEntry~ entries
    }
    class DatabaseBuildExecution {
      BiologicalDatabase database
      result evidence
      inventory path
    }
    ProteinInputRequestDocument --> ProteinInputExecution : prepare
    ProteinInputExecution --> EffectiveDatabaseBuildDocument : input to
    EffectiveDatabaseBuildDocument --> DatabaseBuildExecution : execute
    DatabaseBuildExecution *-- BiologicalDatabase
```

*Figure 7. Source preparation and biological build types.*

The source document represents authored parameters. Preparation returns a checksummed Parquet
handoff. The effective build document is replayable, while `BiologicalDatabase` is the immutable
runtime product from which both FASTA and inventory are projected.

```mermaid
classDiagram
    direction TB
    class DecoyRequestDocument {
      required strategy variant
      prefix and destination
    }
    class SearchDatabase {
      tuple~ProteinInventoryEntry or DecoyInventoryEntry~ entries
    }
    class DecoyExecution {
      SearchDatabase database
      result evidence
      inventory path
    }
    class PeptideBuildRequestDocument {
      digestion parameters
      execution behavior
      artifact destinations
    }
    class PeptideDatabase {
      tuple~PeptideInventoryEntry~ peptides
      tuple~ProteinPeptideMapping~ mappings
    }
    class PeptideBuildExecution {
      PeptideDatabase database
      result evidence
      artifact paths
    }
    DecoyRequestDocument --> DecoyExecution : execute against biological inventory
    DecoyExecution *-- SearchDatabase
    DecoyExecution --> PeptideBuildRequestDocument : inventory input
    PeptideBuildRequestDocument --> PeptideBuildExecution : execute
    PeptideBuildExecution *-- PeptideDatabase
```

*Figure 8. Search-database and peptide workflow types.*

The decoy request is a discriminated storage union for reverse, shuffle, or DecoyPYrat. Resolution
constructs one behavior-owning runtime generation. Peptide execution similarly selects memory,
SQLite, or DuckDB once; every behavior returns the same canonical `PeptideDatabase`.
