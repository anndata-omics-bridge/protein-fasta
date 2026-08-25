# Protein FASTA — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --frozen --group dev` |
| Format | `.venv/bin/ruff format src tests benchmarks && .venv/bin/ruff check --fix src tests benchmarks` |
| Lint | `.venv/bin/ruff check src tests benchmarks` |
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

`protein_fasta` is a strict component tree. Root module `compile.py` may compose the independent
`schema`, `headers`, `classification`, and `validation` children. Child packages never import the
root or one another. `reading` is independent of every other child. The exhaustive contracts in
`.importlinter` are the executable definition.

Pydantic is restricted to `schema/` and the root composition module. Runtime reading, header,
classification, and validation modules do not import Pydantic or either consuming application.

### MUST

- Declare every imported runtime dependency directly in `[project.dependencies]`.
- Put tests, linting, typing, building, and documentation tools in dependency
  groups; optional user-facing capabilities belong in extras.
- Update `pyproject.toml` and `uv.lock` together and run `make check`.

### SHOULD

- Prefer the standard library, then an existing direct dependency, then a small,
  maintained, typed dependency.
- Keep source independent of test, build, documentation, and CLI-only packages.

### MUST NOT

- Depend on unpinned branches or undeclared transitive dependencies.
- Add parallel manifests, lockfiles, formatters, type checkers, or test runners.
- Silence a dependency or typing defect instead of fixing its source.

## Workflow

1. Preserve unrelated worktree changes.
2. Add or update focused tests with each behavioral change.
3. Run the smallest relevant check while iterating.
4. Run `make check` before handoff and report its actual result.
