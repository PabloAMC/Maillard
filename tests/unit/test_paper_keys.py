"""Referential integrity of the paper key registry (data/keys/papers.yml)."""
from __future__ import annotations

import subprocess
import sys
from collections import Counter

from src import data_paths, paper_keys
from scripts.generators.build_paper_registry import gather


def test_every_doi_under_data_has_exactly_one_paper_row_and_vice_versa():
    live = {p["doi"] for p in gather()["papers"]}
    registry = {p.doi for p in paper_keys.all_papers()}
    assert live == registry, (
        f"DOIs under data/ not in the registry: {sorted(live - registry)}; "
        f"registry rows citing nothing: {sorted(registry - live)}"
    )


def test_paper_ids_are_unique_and_every_paper_is_cited_somewhere():
    ids = Counter(p.paper_id for p in paper_keys.all_papers())
    assert not [i for i, c in ids.items() if c > 1]
    assert all(p.record_ids for p in paper_keys.all_papers())


def test_doi_normalisation_handles_wiley_and_markdown():
    assert paper_keys.normalise_doi("**10.1111/joss.12567**") == "10.1111/joss.12567"
    assert paper_keys.normalise_doi("https://doi.org/10.3168/JDS.2019-17495 ✔") == "10.3168/jds.2019-17495"
    assert (
        paper_keys.normalise_doi("10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9")
        == "10.1002/1521-3803(20010601)45:3<150::aid-food150>3.0.co;2-9"
    )
    assert paper_keys.normalise_doi("(see 10.1021/jf0200826)") == "10.1021/jf0200826"
    assert paper_keys.normalise_doi("10.0000/example-pea-matrix-package") is None


def test_registry_is_current():
    completed = subprocess.run(
        [sys.executable, "scripts/generators/build_paper_registry.py", "--check"],
        cwd=data_paths.REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
