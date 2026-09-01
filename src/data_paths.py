"""Single source of truth for where the repository's files live.

Three trees, three rules (see tasks/data_restructure_plan.md, §2):

* ``data/``    -- curated inputs. Hand-edited, cited, READ-ONLY at runtime
                  (enforced in CI by scripts/ci/data_readonly_gate.py).
* ``results/`` -- generated artifacts. Regenerable, never hand-edited.
* ``docs/``    -- human-facing documents; a few are also parsed as inputs
                  (the directional claims panel and its report).

Before 2026-09-01 every module re-derived ``Path(__file__).resolve().parents[1] / "data"
/ ...`` for itself (41 ``ROOT`` definitions, ~110 literal ``"data/..."`` strings across
``src/`` and ``scripts/``), two modules resolved data files relative to the process CWD,
and a file move meant a hundred-site edit. Everything now goes through the names below.

Rules for this module:

* One constant per curated file or directory, named after what the file IS, not where
  it is. Adding a curated file means adding a constant here first.
* ``MAILLARD_DATA_ROOT`` / ``MAILLARD_RESULTS_ROOT`` relocate the two trees (tests,
  scratch checkouts, containers). Nothing else reads those variables.
* ``rel()`` is the only sanctioned way to turn a path into the repo-relative string that
  reports, payloads and tests quote. Never hand-write ``"data/lit/x.json"`` in a string.
* No I/O here. Loading, caching and validation live in ``src.data_access``.
"""
from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------- roots
REPO_ROOT: Path = Path(__file__).resolve().parents[1]
DATA_ROOT: Path = Path(os.environ.get("MAILLARD_DATA_ROOT", REPO_ROOT / "data")).resolve()
RESULTS_ROOT: Path = Path(os.environ.get("MAILLARD_RESULTS_ROOT", REPO_ROOT / "results")).resolve()
DOCS_ROOT: Path = REPO_ROOT / "docs"
TESTS_ROOT: Path = REPO_ROOT / "tests"

# --------------------------------------------------------------------------- data/ dirs
LIT_DIR: Path = DATA_ROOT / "lit"
SPECIES_DIR: Path = DATA_ROOT / "species"
BENCHMARKS_DIR: Path = DATA_ROOT / "benchmarks"
PROTOCOLS_DIR: Path = DATA_ROOT / "protocols"
TIMESERIES_DIR: Path = LIT_DIR / "timeseries"
EXTRACTION_DOSSIERS_DIR: Path = LIT_DIR / "extraction_dossiers"
# LLM research dumps kept only so citation_gate can detect laundered citations. Not
# provenance; nothing loads them as inputs.
RESEARCH_CORPUS_DIR: Path = DATA_ROOT / "Gemini_Deep_Research"
# Primary-source PDFs. Local only (~160 MB, gitignored); code cites them by name.
ARTICLES_DIR: Path = DATA_ROOT / "articles"

# --------------------------------------------------------------------------- data/lit
ARRHENIUS_PARAMS: Path = LIT_DIR / "arrhenius_params.yml"
HENRY_CONSTANTS: Path = LIT_DIR / "henry_constants.yml"
BINDING_CONSTANTS: Path = LIT_DIR / "binding_constants.yml"
LIPID_OXIDATION_CALIBRATION: Path = LIT_DIR / "lipid_oxidation_calibration.json"
COMPUTATIONAL_PRIORS: Path = LIT_DIR / "computational_priors.json"
FLAVOR_REFERENCE_PAYLOADS: Path = LIT_DIR / "flavor_reference_payloads.json"
SAFETY_REFERENCE_PAYLOADS: Path = LIT_DIR / "safety_reference_payloads.json"
RETENTION_REFERENCE_PAYLOADS: Path = LIT_DIR / "retention_reference_payloads.json"
EXTRUSION_DAMAGE_REFERENCE_PAYLOADS: Path = LIT_DIR / "extrusion_damage_reference_payloads.json"
PROCESS_STATE_CALIBRATIONS: Path = LIT_DIR / "process_state_calibrations.json"
PROTEIN_SOURCE_REGISTRY: Path = LIT_DIR / "protein_source_registry.json"
MATRIX_DECISION_PANEL: Path = LIT_DIR / "matrix_decision_panel.json"
MATRIX_FAMILY_COVERAGE_REGISTRY: Path = LIT_DIR / "matrix_family_coverage_registry.json"
CHEMISTRY_FAMILY_SCOPE_REGISTRY: Path = LIT_DIR / "chemistry_family_scope_registry.json"
FAMILY_INGESTION_PLAN: Path = LIT_DIR / "family_ingestion_plan.json"
PROCESS_GAP_REGISTRY: Path = LIT_DIR / "process_gap_registry.json"
BENCHMARK_INTAKE_REGISTRY: Path = LIT_DIR / "benchmark_intake_registry.json"
SLR_INCORPORATION_MATRIX: Path = LIT_DIR / "slr_incorporation_matrix.json"
DEEP_RESEARCH_BACKLOG: Path = LIT_DIR / "deep_research_backlog.json"
REFINEMENT_SURROGATE_PATCHES: Path = LIT_DIR / "refinement_surrogate_patches.json"

# --------------------------------------------------------------------------- data/species
DESIRABLE_TARGETS: Path = SPECIES_DIR / "desirable_targets.yml"
OFF_FLAVOUR_TARGETS: Path = SPECIES_DIR / "off_flavour_targets.yml"
TOXIC_MARKERS: Path = SPECIES_DIR / "toxic_markers.yml"
SENSORY_TAGS: Path = SPECIES_DIR / "sensory_tags.yml"
PRECURSORS: Path = SPECIES_DIR / "precursors.yml"

# --------------------------------------------------------------------------- data/ (loose)
FORMULATION_GRID: Path = DATA_ROOT / "formulation_grid.yml"
INTERVENTIONS: Path = DATA_ROOT / "interventions.yml"

# --------------------------------------------------------------------------- data/benchmarks
# Panel discovery is deliberately NON-recursive over BENCHMARKS_DIR (see
# scripts/ci/holdout_guard.py): the subdirectories below are excluded from calibration.
EXTERNAL_VALIDATION_DIR: Path = BENCHMARKS_DIR / "external_validation"
EXTERNAL_VALIDATION_INTAKE_DIR: Path = EXTERNAL_VALIDATION_DIR / "intake"
MAILLARD_PATH_HOLDOUT_DIR: Path = EXTERNAL_VALIDATION_DIR / "maillard_path"
QUARANTINED_BENCHMARKS_DIR: Path = BENCHMARKS_DIR / "quarantined"
STEP_LEVEL_UNREACHABLE_DIR: Path = BENCHMARKS_DIR / "step_level_unreachable"
# Retired evidence record, self-declared partly fabricated; kept at its historical path
# because ~35 provenance notes cite it by that path.
MAILLARD_VALIDATION_BENCHMARKS_MD: Path = BENCHMARKS_DIR / "maillard_validation_benchmarks.md"

# --------------------------------------------------------------------------- data/protocols
MATRIX_EXPERIMENT_INTAKE_SCHEMA: Path = PROTOCOLS_DIR / "matrix_experiment_intake_schema.json"
EXAMPLE_MATRIX_EXPERIMENT_INTAKE: Path = PROTOCOLS_DIR / "example_matrix_experiment_intake.yaml"
PPI_SPI_PRIMARY_BENCHMARK_CONTRACT: Path = PROTOCOLS_DIR / "ppi_spi_primary_benchmark_contract.json"
EXTRUSION_EXTERNAL_CLOSURE_CONTRACT: Path = PROTOCOLS_DIR / "extrusion_external_closure_contract.json"
DHA_LYSINOALANINE_EXTERNAL_PACKAGE_CONTRACT: Path = (
    PROTOCOLS_DIR / "dha_lysinoalanine_external_package_contract.json"
)
PEA_SOY_MIXED_EXTERNAL_PACKAGE_CONTRACT: Path = PROTOCOLS_DIR / "pea_soy_mixed_external_package_contract.json"

# --------------------------------------------------------------------------- docs/ inputs
DIRECTIONAL_CLAIMS_PANEL: Path = DOCS_ROOT / "validation" / "directional_claims_panel.yml"
DIRECTIONAL_ACCURACY_REPORT: Path = DOCS_ROOT / "validation" / "directional_accuracy_report.md"
PRIMARY_BENCHMARK_PROTOCOL_MD: Path = DOCS_ROOT / "protocols" / "PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md"
SCIENTIFIC_REFERENCE_MD: Path = DOCS_ROOT / "reference" / "SCIENTIFIC_REFERENCE.md"
EXAMPLES_DIR: Path = DOCS_ROOT / "examples"
INGEST_TEMPLATES_DIR: Path = DOCS_ROOT / "templates"

# --------------------------------------------------------------------------- results/
VALIDATION_DIR: Path = RESULTS_ROOT / "validation"
EXPERIMENT_REQUESTS_DIR: Path = VALIDATION_DIR / "experiment_requests"
CALIBRATION_HISTORY_DIR: Path = RESULTS_ROOT / "calibration_history"
# Written by src/matrix_recalibration and read back by src/matrix_calibration_registry.
# Until 2026-09-01 this feedback loop ran through data/lit/ (a generated file inside the
# curated tree); it is a result.
MATRIX_CALIBRATION_OFFSETS: Path = RESULTS_ROOT / "calibration" / "matrix_calibration_offsets.json"
SLR_FAMILY_REPORTS_DIR: Path = RESULTS_ROOT / "literature" / "slr_family_reports"
RESULTS_DB: Path = RESULTS_ROOT / "maillard_results.db"

# --------------------------------------------------------------------------- tests/
TEST_FIXTURES_DIR: Path = TESTS_ROOT / "fixtures"


def rel(path: Path | str) -> str:
    """Repo-relative POSIX string for ``path`` (``"data/lit/x.json"``).

    Reports, payloads and provenance notes quote paths in this form, and several tests
    pin the exact strings. Paths outside the repository are returned absolute.
    """
    candidate = Path(path)
    try:
        return candidate.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return candidate.as_posix()


def benchmark_path(benchmark_id: str) -> Path:
    """Path of a panel benchmark by id (``benchmark_id`` is the file stem)."""
    return BENCHMARKS_DIR / f"{benchmark_id}.json"
