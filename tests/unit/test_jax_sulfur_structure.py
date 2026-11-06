"""The JAX sulfur lane's structure against the NumPy lane (plan step A2 + H3; docs/guides/A2_MODULE_CONTRACT.md sec.
6.2.2). Core env, NO JAX: ``src/kinetic_core_jax/sulfur_structure.py`` is plain numpy, so this runs in the unit tier of
every PR. It is the half of the one-writer rule's replacement that every env can run: an edit of
``src/kinetic_core/sulfur.py`` (or of the ``parameters_sulfur`` / ``ph_state`` logic the port re-implements) that the
JAX lane does not follow fails here, with a message naming the file to update in the same PR. The numerical half
(rates, panel trajectories, random draws, the gradient) is tests/scientific/test_jax_sulfur_parity.py, in the
``maillard-jax`` env (CI job "JAX sulfur lane parity").
"""
from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

from src.kinetic_core_jax import sulfur_structure as ss

ROOT = Path(__file__).resolve().parents[2]


def test_the_numpy_sulfur_lane_has_not_moved_under_the_jax_port():
    problems = ss.structure_problems()
    assert not problems, "\n".join(problems)


def test_every_hand_ported_reaction_has_a_term_in_the_jax_lane():
    """The JAX rate function's ``derived`` table plus the held-at-zero list is exactly ``HAND_PORTED``: read from the
    module's source, so this needs no JAX."""
    tree = ast.parse((ROOT / ss.JAX_FILE).read_text(encoding="utf-8"))
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "_make_rate_constants")
    tables = [n.value for n in ast.walk(fn) if isinstance(n, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == "derived" for t in n.targets)]
    assert len(tables) == 1 and isinstance(tables[0], ast.Dict)
    coded = {k.value for k in tables[0].keys}
    assert not coded & set(ss.OFF)
    assert coded | set(ss.OFF) == set(ss.HAND_PORTED), (
        f"coded {sorted(coded)}, held at zero {sorted(ss.OFF)}, registry {sorted(ss.HAND_PORTED)} -- {ss.UPDATE}")


def test_the_registry_covers_the_engine_exactly():
    found = ss.classify()
    assert found["custom"] == set(ss.HAND_PORTED)
    assert found["zero"] == set(ss.OFF)
    for name, (_, reads) in ss.HAND_PORTED.items():
        assert tuple(sorted(reads)) == found["reads"][name], name
    derived = ss.derive()
    assert len(derived["reactions"]) == len(set(derived["reactions"]))


# --------------------------------------------------------------------------- negative controls: the check fires
def test_an_engine_term_the_port_lacks_is_reported(monkeypatch):
    registry = dict(ss.HAND_PORTED)
    registry.pop("r_amadori")
    monkeypatch.setattr(ss, "HAND_PORTED", registry)
    problems = ss.structure_problems()
    assert any("r_amadori" in p and "update src/kinetic_core_jax/sulfur_lane.py" in p for p in problems), problems


def test_a_new_custom_rate_in_the_engine_is_reported(monkeypatch):
    """A reaction whose constant ``sulfur_rate_constants_at`` overrides (here: r_schiff doubled after the generic
    loop, as a new derived term would) is found numerically, and the edit of the function trips its pin."""
    from src.kinetic_core import sulfur

    original = sulfur.sulfur_rate_constants_at

    def patched(*args, **kwargs):
        out = original(*args, **kwargs)
        out["r_schiff"] = 2.0 * out["r_schiff"]
        return out

    monkeypatch.setattr(sulfur, "sulfur_rate_constants_at", patched)
    problems = ss.structure_problems()
    assert any("r_schiff" in p and "outside the generic rule" in p for p in problems), problems


def test_an_unhandled_switch_and_a_moved_default_are_reported(monkeypatch):
    from src.kinetic_core import config

    handling = dict(ss.CONFIG_HANDLING)
    handling.pop("sulfur.ttca_formation")
    monkeypatch.setattr(ss, "CONFIG_HANDLING", handling)
    monkeypatch.setattr(ss, "HELD", tuple(k for k in ss.HELD if k != "sulfur.ttca_formation"))
    problems = ss.structure_problems(numeric=False)
    assert any("sulfur.ttca_formation" in p and "neither ports nor refuses" in p for p in problems), problems
    monkeypatch.undo()
    shipped = config.Config.shipped()
    moved = shipped.replace({"sulfur.strecker_sulfur": dataclasses.replace(shipped["sulfur.strecker_sulfur"],
                                                                           enabled=True)})
    monkeypatch.setattr(config.Config, "shipped", classmethod(lambda cls: moved))
    problems = ss.structure_problems(numeric=False)
    assert any("shipped defaults of the held keys moved" in p for p in problems), problems


def test_a_logic_edit_of_a_ported_function_is_reported(monkeypatch):
    pins = dict(ss.SOURCE_PINS)
    key = ("src/kinetic_core/sulfur.py", "sulfur_rate_constants_at")
    pins[key] = "0" * 16
    monkeypatch.setattr(ss, "SOURCE_PINS", pins)
    problems = ss.structure_problems(numeric=False)
    assert any("sulfur_rate_constants_at changed" in p and ss.UPDATE in p for p in problems), problems


def test_the_source_pin_ignores_docstrings_and_comments(tmp_path, monkeypatch):
    a = tmp_path / "a.py"
    b = tmp_path / "b.py"
    a.write_text('def f(x):\n    """doc"""\n    return x + 1\n')
    b.write_text('def f(x):\n    # a comment\n    """another doc"""\n    return x + 1\n')
    c = tmp_path / "c.py"
    c.write_text('def f(x):\n    return x + 2\n')
    monkeypatch.setattr(ss, "ROOT", tmp_path)
    assert ss.source_hash("a.py", "f") == ss.source_hash("b.py", "f")
    assert ss.source_hash("a.py", "f") != ss.source_hash("c.py", "f")


def test_the_registry_module_needs_no_jax():
    tree = ast.parse((ROOT / "src/kinetic_core_jax/sulfur_structure.py").read_text(encoding="utf-8"))
    imported = {n.names[0].name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import)} | {
        (n.module or "").split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert not imported & {"jax", "diffrax", "equinox"}


@pytest.mark.parametrize("key", sorted(ss.HELD))
def test_every_held_key_names_what_the_port_implements(key):
    how, what = ss.CONFIG_HANDLING[key]
    assert how == "held" and what
