VENV_BIN := .venv/bin

.DEFAULT_GOAL := help
.PHONY: help sync schemas format format-check lint imports typecheck deps test build docs check clean

help:  ## Show developer commands
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

sync:  ## Synchronize the locked development environment
	uv sync --frozen --group dev --extra cli --extra frame --extra duckdb --extra generation

schemas:  ## Regenerate committed Pydantic JSON Schemas
	$(VENV_BIN)/python scripts/generate_json_schemas.py

format:  ## Format and autofix source and tests
	$(VENV_BIN)/ruff format src tests benchmarks scripts docs_macros.py
	$(VENV_BIN)/ruff check --fix src tests benchmarks scripts docs_macros.py

format-check:  ## Check formatting without changing files
	$(VENV_BIN)/ruff format --check src tests benchmarks scripts docs_macros.py

lint:  ## Run Ruff lint checks
	$(VENV_BIN)/ruff check src tests benchmarks scripts docs_macros.py

imports:  ## Enforce directed package dependencies
	$(VENV_BIN)/lint-imports

typecheck:  ## Run standard Pyright in strict mode
	$(VENV_BIN)/pyright

deps:  ## Validate dependency declarations
	$(VENV_BIN)/deptry .

test:  ## Run tests with branch coverage
	$(VENV_BIN)/pytest --cov --cov-branch

build:  ## Build, validate, and smoke-test source and wheel distributions
	uv build --clear
	$(VENV_BIN)/twine check dist/*
	uv run --isolated --no-project --no-cache \
		--with "$$(printf '%s\n' dist/*.whl)" \
		python -c 'import importlib.util; import protein_fasta.database_build; import protein_fasta.record; assert importlib.util.find_spec("polars") is None'
	@output="$$(uv run --isolated --no-project --no-cache \
		--with "$$(printf '%s\n' dist/*.whl)" \
		python -c 'from protein_fasta.build.generation.decoy import make_decoy_generation; from protein_fasta.schema.build import DecoyDocument, DecoyMode; make_decoy_generation(DecoyDocument(mode=DecoyMode.SHUFFLE))' 2>&1)"; \
		status=$$?; \
		test $$status -ne 0; \
		printf '%s\n' "$$output" | grep -Fq "protein-fasta[generation]"
	uv run --isolated --no-project --no-cache \
		--with 'pytest>=9,<10' \
		--with "$$(printf '%s\n' dist/*.whl)[frame]" \
		pytest -q \
		tests/frame/test_frame.py::test_uniprot_frame_peels_decorations_and_extracts_best_columns \
		tests/frame/test_frame.py::test_refseq_frame_extracts_accession_name_and_optional_organism
	uv run --isolated --no-project --no-cache \
		--with "$$(printf '%s\n' dist/*.whl)[cli]" \
		protein-fasta --help
	uv run --isolated --no-project --no-cache \
		--with 'pytest>=9,<10' \
		--with "$$(printf '%s\n' dist/*.whl)[cli,duckdb]" \
		pytest -q \
		tests/test_cli.py::test_digest_writes_peptides_with_missed_cleavage_evidence \
		tests/test_cli.py::test_registry_cli_functions_cover_index_compare_pairs_and_cluster \
		tests/test_cli.py::test_index_accepts_json_registry_policy \
		tests/test_cli.py::test_build_resolves_profile_request_and_writes_typed_result

docs:  ## Build documentation with strict warnings
	uv run --frozen --group docs zensical build --clean --strict

check:  ## Run every merge-blocking quality gate
	uv lock --check
	$(MAKE) format-check lint imports typecheck deps test build docs

clean:  ## Remove generated build and quality artifacts
	$(VENV_BIN)/python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', 'public', '.pytest_cache', '.ruff_cache')]"
