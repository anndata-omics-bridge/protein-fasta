# Build workflows

The use cases come first. Configuration documents, runtime classes, CLI commands, and storage
formats exist to make these workflows reproducible from either Python or the shell.

Rounded blue nodes are computations. Green rectangles are exchanged data. Cylinders are durable
stores, and amber diamonds are decisions. All diagrams flow top to bottom for portrait reading.

## Routine protein database build

```mermaid
flowchart TB
    A(["fasta_gen source adapters<br/>UniProt · upload<br/>curated catalog"]) --> B["selected sources<br/>JSON + protein-input.parquet"]
    P["build profile JSON"] --> C(["resolve effective request"])
    R["run request JSON<br/>plus typed CLI overrides"] --> C
    subgraph BUILD["protein-fasta build"]
        direction TB
        D(["assemble biological entries"])
        N{"decoy configured?"}
        O(["generate decoys"])
        W(["render and persist database"])
        D --> N
        N -->|yes| O
        O --> W
        N -->|no| W
    end
    B --> D
    C --> D
    W --> F["protein FASTA<br/>effective JSON<br/>result JSON"]
    F --> H(["candidate comparison"])
    I[("existing SQLite/<br/>DuckDB registry")] --> H
    H --> J{"user/application accepts build?"}
    J -->|yes| K(["fasta_gen atomic install"])
    K --> L(["protein-fasta index"])
    J -->|no| M(["discard staged artifacts"])

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    classDef decision fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#4e2600
    class A,C,D,O,W,H,K,L,M compute
    class B,P,R,F,I data
    class N,J decision
```

`protein-fasta build` owns the complete database build, including configured decoy generation.
Biological assembly and decoy generation remain separate internal computations with separate
evidence and tests. Installation is intentionally outside the build because it mutates a
site-owned collection.

## Entrapment database build

```mermaid
flowchart TB
    A["registered source selection"] --> B(["remove prior sentinel,<br/>markers and decoys"])
    F["optional foreign<br/>registered source"] --> C(["compile strategy-specific runtime"])
    E["build request JSON<br/>entrapment + decoy"] --> C
    B --> C
    subgraph BUILD["protein-fasta build"]
        direction TB
        D(["assemble target + contaminant<br/>+ entrapment"])
        N{"decoy configured?"}
        G(["generate decoys over all<br/>biological entries"])
        W(["render and persist database"])
        D --> N
        N -->|yes| G
        G --> W
        N -->|no| W
    end
    C --> D
    W --> H["protein FASTA<br/>effective JSON<br/>result JSON"]
    C --> I["entrapment pairs"]
    H --> J(["optional peptide build"])
    I --> J
    J --> K["peptide FASTA<br/>mapping Parquet<br/>result JSON"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    classDef decision fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#4e2600
    class B,C,D,G,W,J compute
    class A,E,F,H,I,K data
    class N decision
```

Entrapment is a biological-source stage, not a decoy strategy. Entrapment entries join the target
space first; when decoys are configured, `build` generates counterparts for targets,
contaminants, and entrapments together.

## Configuration precedence

```mermaid
flowchart TB
    A["packaged FGCZ<br/>profile JSON"] --> D(["resolve once"])
    B["explicit profile JSON"] --> D
    C["request JSON"] --> D
    E["explicit typed<br/>CLI overrides"] --> D
    D --> F["effective request JSON"]
    F --> G(["protein-fasta build"])
    G --> H["build result JSON"]

    classDef compute fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:#0d2742
    classDef data fill:#e8f5e9,stroke:#2e7d32,stroke-width:1.5px,color:#143d1c
    class A,B,C,E,F,H data
    class D,G compute
```

Precedence is `packaged profile < explicit profile < request < explicitly supplied CLI option`.
The effective document contains resolved paths and is written before FASTA processing begins.
An explicit request value of `"decoy": null`, or CLI `--decoy none`, disables profile decoys.

## CLI workflow

The request owns per-run facts. The profile owns reusable defaults.

```bash
protein-fasta build "$WORK/request.json" \
  --profile "$WORK/fgcz.json" \
  --date 2026-08-28 \
  --decoy reverse
head -n 3 "$WORK/out/p42261_db1_human_d_20260828.fasta"
```

The command writes three primary files beside one another:

- the protein FASTA;
- `*.fasta.effective.json`, written before sequence computation; and
- `*.fasta.result.json`, containing the effective request, checksummed artifacts, counts,
  normalization, sequence and amino-acid summaries, and decoy/entrapment evidence.

## Programmatic workflow

The Python API composes the same boundaries as the CLI:

```python
from pathlib import Path

from protein_fasta.database_build import resolve_database_build, run_database_build
from protein_fasta.documents import (
    load_database_build_profile,
    load_database_build_request,
)

profile_path = Path("fgcz.json")
request_path = Path("request.json")
profile = load_database_build_profile(profile_path)
request = load_database_build_request(request_path)
effective = resolve_database_build(
    profile,
    request,
    profile_base=profile_path.parent,
    request_base=request_path.parent,
)
execution = run_database_build(effective)

print(execution.result.path)
print(execution.result_path)
```

`fasta_gen` can construct the two Pydantic documents directly instead of loading files, call the
same resolver, and pass the effective document to `run_database_build()`. It receives a frozen
`DatabaseBuildExecution`; its `document` is the same typed result serialized by the CLI.

## Build classes and documents

```mermaid
classDiagram
    direction TB
    class DatabaseBuildProfileDocument {
      <<Pydantic JSON>>
      naming
      metadata
      diagnostics
      default_decoy
      default_entrapment
    }
    class DatabaseBuildRequestDocument {
      <<Pydantic JSON>>
      targets
      output_dir
      date
      name_fields
      explicit overrides
    }
    class EffectiveDatabaseBuildDocument {
      <<Pydantic JSON>>
      complete replayable request
    }
    class DatabaseBuildExecution {
      <<frozen runtime>>
      result
      document
      effective_request_path
      result_path
    }
    class DatabaseBuildResultDocument {
      <<Pydantic JSON>>
      artifacts
      counts
      normalization
      summary
      generation evidence
    }

    DatabaseBuildProfileDocument --> EffectiveDatabaseBuildDocument : resolve_database_build
    DatabaseBuildRequestDocument --> EffectiveDatabaseBuildDocument : resolve_database_build
    EffectiveDatabaseBuildDocument --> DatabaseBuildExecution : run_database_build
    DatabaseBuildExecution o-- DatabaseBuildResultDocument
```

Pydantic classes are storage and exchange documents. They validate and serialize; sequence work
remains in the database-build computation. The low-level build receives either one decoy
specification or `None`; it has no separate boolean that can disagree with the configured stage.
