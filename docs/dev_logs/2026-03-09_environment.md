# Maillard Project Lessons & Guardrails

## Test Naming Convention
- **Lesson**: Avoid generic names like `test_phase12_features.py`.
- **Guardrail**: Always use descriptive names that reflect the chemical or physical principle being tested (e.g., `test_advanced_kinetics.py`).
- **Context**: User found `test_phase12_features.py` to be a "terrible undescriptive name".

## Environment Management
- **Lesson**: Discrepancies between `conda_env` and `.venv` can lead to confusing test skips.
- **Guardrail**: Consolidate all SOTA dependencies (MACE, PySCF, Sella) into the primary `conda_env` to match `environment.yml`. Always verify with `pip list` in the target environment.
