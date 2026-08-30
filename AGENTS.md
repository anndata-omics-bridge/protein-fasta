# Protein FASTA — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --frozen --group dev --extra cli --extra frame --extra duckdb --extra generation` |
| Format | `.venv/bin/ruff format src tests benchmarks scripts docs_macros.py && .venv/bin/ruff check --fix src tests benchmarks scripts docs_macros.py` |
| Lint | `.venv/bin/ruff check src tests benchmarks scripts docs_macros.py` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Architecture | `.venv/bin/lint-imports` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Documentation | `make docs` |
| Full gate | `make check` |

## Code conventions

- Fully annotate every function and method in `src/` and `tests/`, including
  private functions, callbacks, generators, fixtures, and special methods.
- Standard Pyright strict and Ruff are mandatory. Do not create baselines,
  blanket exclusions, file-wide ignores, or unqualified `# type: ignore`.
- Ruff is the sole formatter and linter. Do not add Black, isort, Flake8, mypy,
  or another overlapping formatter/type checker.
- Keep `__init__.py` empty and import public objects from their defining modules.
- Use Google-style docstrings for public APIs and the configured 100-character
  line length.

## Dependency rules

`protein_fasta` has an explicit inward dependency structure. The outer `cli.py` adapter may compose
the root products. `protein_input.py` and `database_build.py` own source preparation and biological
assembly; `decoy_database.py` owns the subsequent search-database operation. `database/` owns
schema-independent naming, metadata, and sequence-generation behavior; root
`database_compile.py` and `decoy_compile.py` translate passive documents into those runtime
values. `analytics/` owns backend-free hashing, digestion, comparison, and clustering and never
imports database construction or persistence.
`registry/` owns indexing, comparison queries, schema versioning, snapshots, and concrete
SQLite/optional DuckDB adapters. `uniprot/` and `peptide/` own their focused runtime capabilities;
root workflow modules compose them with storage documents and artifacts. `record.py`, `frame.py`,
`reading/`, and `summary.py` retain their focused products. The exhaustive contracts in
`.importlinter` are the executable definition.

Pydantic documents are passive and live in `schema/`; loading and compilation occur at root
boundaries. Polars is restricted to the optional CLI/frame and registry-export boundary. Scalar
reading, normalization, diagnostics, and analytics import neither framework nor a consuming
application. Sequence and peptide hashes are versioned BLAKE2b-128 over already-normalized
supplied values; exact-file provenance is versioned non-security MD5. Build naming and `aa|`
metadata construction do not own analytical hashes or registry persistence.

### MUST

- Declare every imported runtime dependency directly in `[project.dependencies]`.
- Put tests, linting, typing, building, and documentation tools in dependency
  groups; optional user-facing capabilities belong in extras.
- Update `pyproject.toml` and `uv.lock` together and run `make check`.

### SHOULD

- Prefer the standard library, then an existing direct dependency, then a small,
  maintained, typed dependency.
- Keep source independent of test, build, documentation, and CLI-only packages.
- Keep the executable CLI walkthrough shell-first: authored commands live in its Markdown page,
  Zensical captures them through `docs_macros.py`, and `make docs` must fail on command drift.

### MUST NOT

- Depend on unpinned branches or undeclared transitive dependencies.
- Add parallel manifests, lockfiles, formatters, type checkers, or test runners.
- Silence a dependency or typing defect instead of fixing its source.

## Workflow

1. Preserve unrelated worktree changes.
2. Add or update focused tests with each behavioral change.
3. Run the smallest relevant check while iterating.
4. Run `make check` before handoff and report its actual result.
