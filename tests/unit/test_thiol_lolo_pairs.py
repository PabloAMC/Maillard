"""H0's residual-only pair-fit runner (results/validation/kinetic_core_between_lab_tau_h0_prereg.md sec. 4.3): a pair
(L, L') trains on the data minus both laboratories; the single-laboratory case reproduces EX-056's folds; a record
written on a stub engine carries no score key, and the runner's stdout carries no residual statistic."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))

import fit_thiol_lolo as lolo  # noqa: E402
import fit_thiol_lolo_pairs as pairs  # noqa: E402

LABS = ["jiangnan_zhang_x", "nestle_blank", "dfa_hofmann_schieberle", "reading_mottram", "ochanomizu_murata",
        "dalian_polytechnic_chen_d"]


def _row(i, laboratory, ledger_id=None):
    return SimpleNamespace(ledger_id=ledger_id or f"paper{i % 4}{laboratory[:3]}_t1", laboratory=laboratory,
                           species="MFT", t_min=float(i), method="sida", table="thiol", pot=f"p{i}",
                           sigma=0.1, value=1.0, unit="ug/L", series=f"s{i % 3}", i=i)


def _rows():
    return [_row(i, LABS[i % len(LABS)]) for i in range(30)]


def test_a_pair_trains_on_the_data_minus_both_laboratories():
    rows = _rows()
    for a, b in pairs.pairs():
        train, held = pairs.split(rows, (a, b))
        assert not {lolo.lab_of(r) for r in train} & {a, b}
        assert {r.i for r in train} == {r.i for r in rows if lolo.lab_of(r) not in (a, b)}
        for lab in (a, b):
            assert [r.i for r in held[lab]] == [r.i for r in rows if lolo.lab_of(r) == lab]
        assert len(train) + len(held[a]) + len(held[b]) == len(rows)


def test_the_single_laboratory_case_reproduces_make_folds():
    rows = _rows()
    for lab, (train, held) in lolo.make_folds(rows).items():
        t, h = pairs.split(rows, (lab,))
        assert [r.i for r in t] == [r.i for r in train]
        assert [r.i for r in h[lab]] == [r.i for r in held]


def test_the_ten_pairs_are_the_scorable_laboratories_and_never_nestle():
    ps = pairs.pairs()
    assert len(ps) == 10 and len(set(ps)) == 10
    assert all(a < b for a, b in ps)
    assert not {lab for p in ps for lab in p} & set(pairs.NEVER_HELD)
    assert {lab for p in ps for lab in p} == set(pairs.SCORABLE)


def test_converged_means_status_above_zero():
    assert pairs.converged({"status": 2}) and pairs.converged({"status": 1})
    assert not pairs.converged({"status": 0}) and not pairs.converged({"status": -1})


def _walk_keys(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _walk_keys(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk_keys(v)


@pytest.fixture
def stub_engine(monkeypatch, tmp_path):
    rows = _rows()
    keys = ("k1", "k2")
    monkeypatch.setattr(pairs, "load_fit_table", lambda: (tuple(rows), {}, {}))
    monkeypatch.setattr(pairs.fj, "coordinates",
                        lambda *a, **k: (keys, ("log10k", "log10k"), np.array([-5.0, -5.0]), np.array([5.0, 5.0]),
                                         np.array([0.0, 0.0])))
    full = {"fit": {"theta": {"k1": 0.5, "k2": -0.5}, "cost": 7.5},
            "laplace": {"identified": {"k1": True, "k2": True}, "at_bound": {"k1": False, "k2": False}}}
    noise = {"series_extra_noise": {pairs.fj.series_of(r): 0.1 for r in rows}}
    monkeypatch.setattr(pairs, "_fit_record", lambda suffix: full if suffix == pairs.FULL_FIT else noise)
    raw_value = lambda theta, rs: np.array([0.01 * r.i + float(np.sum(theta)) for r in rs])  # noqa: E731
    monkeypatch.setattr(pairs.fj, "make_raw_many",
                        lambda jobs, ks: (lambda points, rs: [raw_value(np.asarray(p), rs) for p in points]))
    # the stub's whiten makes the full-data cost equal the record's, so the reproduction check passes
    model = SimpleNamespace(whiten=lambda r, rs: np.full(len(rs), np.sqrt(2 * 7.5 / len(rs))) if len(rs) == len(rows)
                            else np.asarray(r, dtype=float), nuisance_dof=lambda rs: 0)
    monkeypatch.setattr(pairs.lolo, "_model", lambda *a, **k: model)

    def fake_fit(problem, theta, budget=1200, ftol=None, **kw):
        theta = np.asarray(theta, dtype=float)
        return SimpleNamespace(theta=theta, status=2, cost=1.0, evaluations=3, seconds=0.1,
                               as_dict=lambda ks: {"theta": dict(zip(ks, map(float, theta))), "cost": 1.0,
                                                   "start_cost": 2.0, "start": dict(zip(ks, map(float, theta))),
                                                   "status": 2, "budget_exhausted": False, "evaluations": 3,
                                                   "budget": budget, "seconds": 0.1})
    monkeypatch.setattr(pairs.fitting, "fit", fake_fit)
    # k2 loses identification without the pair: the reset arm must run
    monkeypatch.setattr(pairs.fitting, "laplace",
                        lambda problem, theta, **kw: {"identified": {"k1": True, "k2": False},
                                                      "at_bound": {"k1": False, "k2": False},
                                                      "sigma": {"k1": 0.1, "k2": float("inf")}})
    monkeypatch.setattr(pairs, "out_path", lambda labs, guard=False: tmp_path / f"pair_{'__'.join(labs)}.json")
    monkeypatch.setattr(pairs, "_provenance", lambda name: {"artifact": name})
    return rows


def test_a_stub_record_carries_raw_residuals_and_no_score(stub_engine, capsys, tmp_path):
    a, b = "dfa_hofmann_schieberle", "reading_mottram"
    state = pairs.run_pair((a, b), budget=10, max_continuations=1, freeze_commit="abc")
    record = json.loads((tmp_path / f"pair_{a}__{b}.json").read_text())
    assert not set(_walk_keys(record)) & pairs.FORBIDDEN_KEYS
    assert record["complete"] and record["freeze_commit"] == "abc" and record["leak_prone"] == ["k2"]
    assert set(record["held_engine_raw"]) == {"warm", "reset"}
    for v in ("warm", "reset"):
        for lab in (a, b):
            assert len(record["held_engine_raw"][v][lab]) == len(record["row_keys"]["held"][lab])
        assert len(record["train_engine_raw"][v]) == len(record["row_keys"]["train"]) == record["training_rows"]
    assert record["warm"]["converged"] is True
    out = capsys.readouterr().out
    for word in ("tau", "omega", "coverage", "log_density", "lpd", "rms", "median", "mean", "residual"):
        assert word not in out.lower(), f"stdout mentions {word!r}: {out}"
    assert state["training_rows"] == sum(1 for r in stub_engine if lolo.lab_of(r) not in (a, b))
