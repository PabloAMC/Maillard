# Maillard convenience targets.
#
# Every recipe goes through ./scripts/docker_maillard.sh run "..." per the
# project rule that all execution happens inside the validated container.

DOCKER := ./scripts/docker_maillard.sh run

.PHONY: help trust-loop artifacts test test-unit test-scientific summary clean-pyc

help:
	@echo "Maillard Make targets (all run inside Docker):"
	@echo "  make artifacts        Regenerate every live artifact in dependency order (alias: trust-loop)"
	@echo "  make test-unit        Run fast unit tests"
	@echo "  make test-scientific  Run scientific regression tests"
	@echo "  make test             Run the full test suite"
	@echo "  make summary          Regenerate the benchmark summary artifact"
	@echo "  make clean-pyc        Remove __pycache__ directories"

trust-loop: artifacts

# Regenerate every live artifact in dependency order (2026-10-28: the old target called three generators
# deleted at the lane retirement; the ranking must precede the wishlist, which reads it).
artifacts:
	./scripts/docker_maillard.sh core-envelope --workers 3
	./scripts/docker_maillard.sh core-scores
	./scripts/docker_maillard.sh core-directional
	./scripts/docker_maillard.sh run "python scripts/generators/build_sulfur_flux_attribution.py"
	./scripts/docker_maillard.sh run "python scripts/generators/build_null_model_baselines.py"
	./scripts/docker_maillard.sh experiment-value-ranking
	./scripts/docker_maillard.sh wishlist
	./scripts/docker_maillard.sh run "python scripts/generators/generate_gap_heatmap.py"
	./scripts/docker_maillard.sh model-card
	./scripts/docker_maillard.sh results-readme
	./scripts/docker_maillard.sh data-readme

test-unit:
	$(DOCKER) "pytest tests/unit -q"

test-scientific:
	$(DOCKER) "pytest tests/scientific -q"

test:
	$(DOCKER) "pytest tests/ -q"

summary:
	./scripts/docker_maillard.sh summary

clean-pyc:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
