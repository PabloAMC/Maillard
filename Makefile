# Maillard convenience targets.
#
# Every recipe goes through ./scripts/docker_maillard.sh run "..." per the
# project rule that all execution happens inside the validated container.

DOCKER := ./scripts/docker_maillard.sh run

.PHONY: help trust-loop test test-unit test-scientific summary clean-pyc

help:
	@echo "Maillard Make targets (all run inside Docker):"
	@echo "  make trust-loop       Regenerate the 5 trust-dashboard artifacts"
	@echo "                        (envelope, external hold-outs, VoI ranking, LOO leverage, gap heatmap)"
	@echo "  make test-unit        Run fast unit tests"
	@echo "  make test-scientific  Run scientific regression tests"
	@echo "  make test             Run the full test suite"
	@echo "  make summary          Regenerate the benchmark summary artifact"
	@echo "  make clean-pyc        Remove __pycache__ directories"

trust-loop:
	$(DOCKER) "python scripts/generators/generate_prediction_uncertainty.py --n-samples 200 --seed 0"
	$(DOCKER) "python scripts/generators/generate_external_validation_report.py --n-samples 200 --seed 0"
	$(DOCKER) "python scripts/generators/generate_experiment_value_ranking.py"
	$(DOCKER) "python scripts/generators/generate_loo_leverage.py"
	$(DOCKER) "python scripts/generators/generate_gap_heatmap.py"
	@echo ""
	@echo "Trust-loop artifacts refreshed under results/validation/:"
	@echo "  - prediction_uncertainty.{md,json}"
	@echo "  - external_validation_report.{md,json}"
	@echo "  - experiment_value_ranking.{md,json}"
	@echo "  - loo_leverage.{md,json}"
	@echo "  - gap_heatmap.png"

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
