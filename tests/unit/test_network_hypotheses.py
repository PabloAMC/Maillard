"""The hypothesis layer: cited rules with controls, structures, the firewall, and a small enumeration."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from src import data_access, data_paths  # noqa: E402
from src.network_hypotheses import classify, rules as R, structures as S  # noqa: E402
from src.network_hypotheses.enumerate import enumerate_from  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def table():
    return S.load()


@pytest.fixture(scope="module")
def rules():
    return R.load()


def test_every_rule_has_a_source_dossier_on_disk(rules):
    missing = []
    for rule in rules:
        dossier = rule.source["dossier"]
        if not (data_paths.EXTRACTION_DOSSIERS_DIR / dossier).exists():
            missing.append((rule.id, dossier))
        assert rule.source.get("anchor"), f"{rule.id}: no anchor"
        assert rule.status in ("established", "net", "proposed"), rule.id
    assert not missing, missing


def test_every_rule_has_both_kinds_of_control(rules):
    for rule in rules:
        assert rule.controls.get("positive"), f"{rule.id}: no positive control"
        assert rule.controls.get("negative"), f"{rule.id}: no negative control"


def test_rule_controls_pass(rules, table):
    failures = []
    for rule in rules:
        failures += R.run_controls(rule, table)
    assert not failures, "\n".join(failures)


def test_rule_ids_are_unique(rules):
    ids = [r.id for r in rules]
    assert len(ids) == len(set(ids))


def test_the_engine_never_imports_the_hypothesis_layer():
    """The firewall: nothing under src/kinetic_core may import src.network_hypotheses."""
    offenders = []
    for path in (ROOT / "src" / "kinetic_core").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            if any("network_hypotheses" in n for n in names):
                offenders.append(str(path.relative_to(ROOT)))
    assert not offenders, offenders


def test_enumeration_from_the_reference_pot_reproduces_modelled_steps(rules, table):
    """Pentose + cysteine (the Hofmann pot): the rules must rediscover the engine's own route and
    propose nothing that breaks the carbon count."""
    known = classify.species_by_canonical(table)
    steps, pool = enumerate_from(["PENT", "Cys"], table, rules, depth=2, may_expand=lambda c: c in known)
    by_canon = classify.species_by_canonical(table)
    reactions = classify.engine_reactions()
    placements = {}
    for step in steps:
        rule = next(r for r in rules if r.id == step.rule_id)
        placements.setdefault(step.rule_id, set()).add(classify.place_step(step, rule.status, by_canon, reactions)[0])
    # cysteine -> H2S is the engine's r_cys_h2s: modelled. The thiazolidine (Zhai 2020) is charged by the
    # engine but never formed by it (only r_ttca_cys, the ring opening): the layer must report the
    # forward step as mechanism known, not modelled. The 1-deoxypentosone must be reached.
    assert placements.get("R08_cysteine_thermolysis") == {"modelled"}, placements
    assert placements.get("R15_thiazolidine") == {"mechanism_known"}, placements
    assert table["DPO"].canonical in pool, "the 1-deoxypentosone was not reached within two steps"
    # carbon is conserved by every step (each rule writes every carbon it reads)
    Chem = S._rdkit()

    def carbons(smiles):
        return sum(1 for a in Chem.MolFromSmiles(smiles).GetAtoms() if a.GetSymbol() == "C")

    for step in steps:
        rule = next(r for r in rules if r.id == step.rule_id)
        lost = sum(carbons(s) for s in step.reactants) - sum(carbons(s) for s in step.products)
        allowed = 1 if rule.klass == "strecker" or rule.id == "R08_cysteine_thermolysis" or rule.id.startswith("R17") else 0
        assert 0 <= lost <= allowed, (step.rule_id, step.reactants, step.products, lost)


def test_enumeration_is_deterministic_across_processes(rules, table):
    """The artifact is compared byte for byte by the freshness gate, so the same charge must give the
    same steps in the same order under a different hash seed."""
    import json
    import os
    import subprocess
    import sys

    code = (
        "import json; from src.network_hypotheses import rules as R, structures as S; "
        "from src.network_hypotheses.enumerate import enumerate_from; "
        "t = S.load(); steps, _ = enumerate_from(['PENT', 'Cys'], t, R.load(), depth=2); "
        "print(json.dumps([[s.rule_id, list(s.reactants), list(s.products)] for s in steps]))"
    )
    outs = []
    for seed in ("1", "2"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=str(ROOT))
        outs.append(subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True, check=True).stdout)
    assert outs[0] == outs[1]

