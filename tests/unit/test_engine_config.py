"""src/kinetic_core/config.py (A3 phases 1-2, 2026-11-01): the engine's switches as one frozen Config.

Guarantees:
* COVERAGE (phase 2, replacing phase 1's pattern scan): every public UPPERCASE module-level global that a function in
  src/kinetic_core reads at call time is either a ConfigKey binding or on the reviewed allowlist of constants
  (tests/engine_constants.py); a new one, or a stale entry, fails. Every switch a script or test assigns has a key, and
  every module-level re-import of a switch is a listed binding;
* ``Config.shipped()`` is the engine as first imported, and ``from_globals()`` reads it back;
* a Config validates types, choices, bands and field bands on construction, and has a stable ``digest()``;
* build-phase keys are read-only: ``install`` / ``applied`` / ``update`` refuse a non-shipped one, naming it;
* ``applied()`` gives every binding back the very object it held, also on an exception;
* a SPAWN-started worker given a Config with every run-time switch flipped reproduces the forked worker.

Dropping ANY key fails ``test_a_config_is_frozen_hashable_picklable_and_validated`` (FLIPPED names every key).
"""
from __future__ import annotations

import ast
import dataclasses
import importlib
import multiprocessing as mp
import os
import pickle
import subprocess
import sys
import types
from pathlib import Path

import pytest

from src.kinetic_core import config, parallel
from src.kinetic_core import parameters_dicarbonyl as pd
from src.kinetic_core import parameters_lipid as pl
from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import ph_balance_rivals as pbr
from src.kinetic_core import ph_state
from tests import engine_constants, support

ROOT = Path(__file__).resolve().parents[2]
KC = ROOT / "src" / "kinetic_core"
PACKAGE = "src.kinetic_core"

#: Names a script or test assigns that are NOT engine switches, each with the reason.
NOT_SWITCHES = {
    ("parameters_furanic", "UG_PER_MMOL_EQUALS_MG_PER_MOL"): "a unit identity (ug/mmol == mg/mol), asserted, not chosen",
    ("fit_rows", "FIT_ROWS_DIR"): "a data path, redirected to tmp_path by tests",
    ("fit_rows", "LABORATORIES"): "a data path, redirected to tmp_path by tests",
}

#: A non-default value for EVERY key (the spawn test flips them all; a new key needs one here).
FLIPPED = {
    "ph.coupling": "inline",
    "ph.amine_survival": pbr.AmineSurvival(True, 0.5, 0.25),
    "ph.acid_lag": pbr.AcidLag(False, 0.3, 0.5, 90.0),   # B2 with B1 + B3 is refused: flipped in its fields
    "ph.acid_supply": pbr.AcidSupply(True, 0.2, 0.5, 0.9, (5.2, 6.2, 6.7), 0.1),
    "ph.arp_forcing": pbr.ArpForcing(True, (0.0, 20.0), (5.0, 3.2)),
    "ph.node_grid": "uniform",
    "ph.fixed_point": "fixed_sweeps",
    "ph.fixed_point.relaxation": 0.5,
    "ph.fixed_point.tolerance": 5e-4,
    "ph.fixed_point.max_sweeps": 30,
    "ph.fixed_point.flag": 2e-2,
    "ph.b51": ph_state.PhStateCorrections(False, False, False, False, False, False),
    "ph.trunk_solutes_in_balance": True,
    "ph.ionic_strength.sweeps": 25,
    "ph.ionic_strength.tol_m": 1e-10,
    "ph.trunk_terms": True,
    "ph.trunk_exponent_override": {"k_ama_odg": 0.2, "k_schiff": 0.5},
    "sulfur.flux_budget": "linear_rederived",
    "sulfur.profile_node_step_k": 1.5,
    "engine.ramp_step_k": 1.0,
    "sulfur.ttca_formation": dataclasses.replace(ps.TtcaFormation(), enabled=True),
    "sulfur.nf_hydroxide_decay": dataclasses.replace(ps.NfHydroxideDecay(), enabled=True, mgo_share=0.1),
    "sulfur.pentose_mgo_routes": ps.PentoseMgoRoutes(k_tdp_mgo_100c_ph5=1e-5, k_pent_mgo_100c_ph5=1e-6),
    "sulfur.pentose_amine_entry": ps.PentoseAmineEntry(k_entry_l_per_mmol_min_100c_ph6=1e-4,
                                                       k_dpo_aa_per_min_100c_ph6=1e-3),
    "sulfur.ddp_formation_ph": ps.DdpFormationPh(law="slope", slope=0.5),
    "sulfur.mgo_cysteine_strecker": ps.MgoCysteineStrecker(k_l_per_mmol_min_100c_ph5=1e-3),
    "sulfur.tdp_cysteine_trap": ps.TdpCysteineTrap(k_l_per_mmol_min_100c_ph5=1e-3),
    "sulfur.hexose_c5_entry": ps.HexoseC5Entry(k_per_min_100c=2e-3),
    "sulfur.strecker_sulfur": dataclasses.replace(ps.StreckerSulfur(), enabled=True, h2s_share=0.08),
    "sulfur.cys_amadori_route": dataclasses.replace(ps.CysAmadoriRoute(), enabled=True, split_law="base_saturating"),
    "sulfur.tdp_own_decay": ps.TdpOwnDecay(enabled=True),
    "sulfur.furfural_barriers": ps.FurfuralBarriers(enabled=True),
    "sulfur.single_temperature_barriers": ps.SingleTemperatureBarriers(enabled=True),
    "fit.dimer_encounter_ceiling": dataclasses.replace(ps.DimerEncounterCeiling(), enabled=True, mft_cap_mmol_l=2.0),
    "fit.mft_decay_prior": dataclasses.replace(ps.MftDecayPrior(), enabled=True, sd_dex=0.5),
    "trunk.three_deoxy_exit_ph_term": False,
    "trunk.three_deoxy_exit_ph": {"k_tdg_fa": (0.5, (0.4, 0.6), "flipped for the test"),
                                  "k_tdg_mgo": (0.63, (0.37, 0.92), "flipped for the test")},
    "dicarbonyl.aqueous_glyoxal_parameters": pd.with_aqueous_glyoxal(-2.5, -1.0),
    "lipid.q10_assumption": dataclasses.replace(pl.Q10Assumption(), default=2.0),
    "lipid.decadienal_breakage": dataclasses.replace(pl.DecadienalBreakage(), enabled=False),
    "lipid.covalent_sink": dataclasses.replace(pl.CovalentSinkCeiling(), enabled=True),
    "dicarbonyl.shipped_b39": False,
    "methionine.shipped_b22": True,
    "methionine.shipped_b22b": True,
    "proline.shipped_b24": True,
    "proline.shipped_b24b": True,
    "panel.thermal_history": "as_recorded",
    "uncertainty.b55_single_temperature_ea_prior": False,
    "uncertainty.b56_structural_discrepancy": True,
}


#: The run-time keys flipped: a Config that changes a build-phase key cannot be installed (the owner's decision 1).
RUNTIME_FLIPPED = {name: value for name, value in FLIPPED.items() if name not in config.BUILD_KEYS}


def _flipped() -> config.Config:
    return config.Config.shipped().replace(FLIPPED)


def _flipped_runtime() -> config.Config:
    return config.Config.shipped().replace(RUNTIME_FLIPPED)


def _covered():
    return {binding for key in config.KEYS for binding in key.bindings}


# --- the scans: what counts as a switch ---------------------------------------------------------------------------


def _module_names():
    return sorted(p.stem for p in KC.glob("*.py") if p.stem != "__init__")


def _engine_module(name: str):
    return importlib.import_module(f"{PACKAGE}.{name}" if name else PACKAGE)


def _is_global_name(name: str) -> bool:
    return name.isupper() and not name.startswith("_")


def _call_time_reads():
    """{(module, NAME): {modules reading it}} for every public UPPERCASE module-level global that a function (any
    def, method or lambda) in src/kinetic_core reads: by bare name in the module that holds it, or as ``alias.NAME``
    through an import of an engine module. Callables (classes, functions) and modules are not values and are left
    out; module-level code (import-time table building) is not a call-time read."""
    out = {}
    for path in sorted(KC.glob("*.py")):
        here = "" if path.stem == "__init__" else path.stem
        module = _engine_module(here)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        aliases = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and ((node.level == 1 and not node.module) or
                                                     (node.level == 0 and node.module == PACKAGE)):
                aliases.update({a.asname or a.name: a.name for a in node.names if (KC / f"{a.name}.py").exists()})
            elif isinstance(node, ast.Import):
                aliases.update({a.asname: a.name[len(PACKAGE) + 1:] for a in node.names
                                if a.name.startswith(PACKAGE + ".") and a.asname})
        for fn in ast.walk(tree):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                continue
            for n in ast.walk(fn):
                if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load) and _is_global_name(n.id):
                    owner, attr, holder = here, n.id, module
                elif (isinstance(n, ast.Attribute) and isinstance(n.ctx, ast.Load) and isinstance(n.value, ast.Name)
                      and n.value.id in aliases and _is_global_name(n.attr)):
                    owner = aliases[n.value.id]
                    attr, holder = n.attr, _engine_module(owner)
                else:
                    continue
                if not hasattr(holder, attr):
                    continue   # a local or a parameter spelt in capitals
                value = getattr(holder, attr)
                if callable(value) or isinstance(value, types.ModuleType):
                    continue
                out.setdefault((owner, attr), set()).add(here or PACKAGE)
    return out


def _aliases(tree):
    """{local name: kinetic_core module} for every import of a src.kinetic_core module in a file (any scope)."""
    out = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module in (PACKAGE, "src.kinetic_core"):
            for alias in node.names:
                out[alias.asname or alias.name] = alias.name
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith(PACKAGE + ".") and alias.asname:
                    out[alias.asname] = alias.name[len(PACKAGE) + 1:]
    return out


def _assigned_switches():
    """Every (kinetic_core module, UPPERCASE attribute) that a script or a test assigns: ``alias.NAME = ...``,
    ``setattr(alias, "NAME", ...)`` and ``monkeypatch.setattr(alias, "NAME", ...)``."""
    modules = set(_module_names())
    found = {}
    for path in [*ROOT.joinpath("scripts").rglob("*.py"), *ROOT.joinpath("tests").rglob("*.py")]:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        aliases = {k: v for k, v in _aliases(tree).items() if v in modules}
        for node in ast.walk(tree):
            targets = []
            if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                    if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name):
                        targets.append((t.value.id, t.attr))
            elif isinstance(node, ast.Call) and len(node.args) >= 2:
                func = node.func
                is_setattr = (isinstance(func, ast.Name) and func.id == "setattr") or (
                    isinstance(func, ast.Attribute) and func.attr == "setattr")
                obj, attr = node.args[0], node.args[1]
                if is_setattr and isinstance(obj, ast.Name) and isinstance(attr, ast.Constant) \
                        and isinstance(attr.value, str):
                    targets.append((obj.id, attr.value))
            for alias, attr in targets:
                if alias in aliases and attr.isupper() and not attr.startswith("_"):
                    found.setdefault((aliases[alias], attr), str(path.relative_to(ROOT)))
    return found


def _module_level_reimports():
    """Every module-level ``from <a module binding it> import NAME`` of a switch in the engine: another binding."""
    by_attr = {}
    for key in config.KEYS:
        for module, attr in key.bindings:
            by_attr.setdefault(attr, set()).add(module)
    found = set()
    for path in ROOT.joinpath("src").rglob("*.py"):
        rel = path.relative_to(ROOT).with_suffix("")
        parts = rel.parts
        if parts[:2] != ("src", "kinetic_core"):
            continue   # only the engine's modules bind switches at module level today; others import in functions
        here = ".".join(parts[2:]).replace("__init__", "").rstrip(".")
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if not isinstance(node, ast.ImportFrom):
                continue
            source = node.module or ""
            if node.level == 0 and source.startswith(PACKAGE + "."):
                source = source[len(PACKAGE) + 1:]
            for alias in node.names:
                if alias.name in by_attr and source in by_attr[alias.name]:
                    found.add((here, alias.asname or alias.name))
    return found


def test_every_global_an_engine_function_reads_is_a_key_or_a_reviewed_constant():
    reads = _call_time_reads()
    covered, allowed = _covered(), engine_constants.allowlisted()
    # the scan sees the switches the survey saw (A2_MODULE_CONTRACT sec. 1) and the keyed solver numerics
    assert {("sulfur", "PH_COUPLING"), ("ph_state", "B51"), ("parameters_sulfur", "TDP_OWN_DECAY"),
            ("sulfur", "PH_FIXED_POINT_TOLERANCE"), ("lipid", "Q10_ASSUMPTION"), ("engine", "RAMP_STEP_K")} <= set(reads)
    unclassified = sorted(f"{m or PACKAGE}.{a} (read in {', '.join(sorted(where))})"
                          for (m, a), where in reads.items() if (m, a) not in covered and (m, a) not in allowed)
    assert not unclassified, ("globals the engine reads at call time that are neither a ConfigKey binding "
                              "(src/kinetic_core/config.py) nor a reviewed constant (tests/engine_constants.py): "
                              f"{unclassified}")
    both = sorted(set(allowed) & covered)
    assert not both, f"a key's binding is not a constant; drop it from tests/engine_constants.py: {both}"


def test_the_constants_allowlist_has_no_stale_entry():
    reads = _call_time_reads()
    stale = sorted(f"{m or PACKAGE}.{a}" for m, a in engine_constants.allowlisted() if (m, a) not in reads)
    assert not stale, f"no longer read at call time; drop from tests/engine_constants.py: {stale}"
    assert set(engine_constants.CONSTANTS) == set(engine_constants.WHY)   # every category says why


def test_every_switch_a_script_or_test_sets_has_a_config_key():
    covered, allowed = _covered(), engine_constants.allowlisted()
    assigned = {s: where for s, where in _assigned_switches().items() if s not in NOT_SWITCHES}
    unkeyed = sorted(f"{m}.{a} (set in {where})" for (m, a), where in assigned.items() if (m, a) not in covered)
    assert not unkeyed, f"switches a script or test sets, without a ConfigKey: {unkeyed}"
    assert not set(assigned) & set(allowed), "a 'constant' that a script or test sets is a switch: key it"
    assert {("sulfur", "PH_COUPLING"), ("ph_state", "B51")} <= set(assigned)


def test_every_binding_of_a_switch_is_listed_and_exists():
    covered = _covered()
    reimports = _module_level_reimports()
    assert {("lipid", "DECADIENAL_BREAKAGE"), ("", "Q10_ASSUMPTION")} <= reimports   # the scan sees them
    unlisted = sorted(reimports - covered)
    assert not unlisted, f"a switch re-imported at module level is a second binding; add it to the key: {unlisted}"
    for key in config.KEYS:
        for module, attr in key.bindings:
            mod = sys.modules[f"{PACKAGE}.{module}" if module else PACKAGE]
            assert hasattr(mod, attr), f"{key.name}: {module}.{attr} does not exist"


def test_shipped_is_the_engine_as_first_imported_and_from_globals_reads_it():
    pristine = {(module.__name__, name): value for (module, name), value in support._PRISTINE_SWITCHES.items()}
    shipped = config.Config.shipped()
    for key in config.KEYS:
        module, attr = key.bindings[0]
        assert key.freeze(pristine[(f"{PACKAGE}.{module}", attr)]) == shipped.values[key.name], key.name
    assert config.Config.from_globals() == shipped
    assert shipped.diff() == {}


def test_a_config_is_frozen_hashable_picklable_and_validated():
    shipped = config.Config.shipped()
    assert hash(shipped) == hash(config.Config.from_globals())
    flipped = _flipped()
    assert pickle.loads(pickle.dumps(flipped)) == flipped and hash(flipped) != hash(shipped)
    assert set(flipped.diff()) == set(config.KEY_BY_NAME) == set(FLIPPED)   # FLIPPED covers and flips every key
    assert flipped["ph.trunk_exponent_override"] == {"k_ama_odg": 0.2, "k_schiff": 0.5}
    assert shipped.replace(ph__coupling="inline") == shipped.replace({"ph.coupling": "inline"})
    with pytest.raises(dataclasses.FrozenInstanceError):
        shipped.items = ()
    with pytest.raises(KeyError):
        shipped.replace({"ph.no_such_switch": 1})
    with pytest.raises(ValueError, match="not one of"):
        shipped.replace({"ph.coupling": "inline_typo"})
    with pytest.raises(ValueError, match="every key"):
        config.Config(shipped.items[:-1])


def test_a_config_validates_types_bands_and_field_bands_on_construction():
    shipped = config.Config.shipped()
    # a number stands for a float, and is stored as one (40 and 40.0 are one value, one digest)
    as_int = shipped.replace({"ph.fixed_point.relaxation": 1})
    assert type(as_int.values["ph.fixed_point.relaxation"]) is float
    assert as_int == shipped.replace({"ph.fixed_point.relaxation": 1.0})
    assert as_int.digest() == shipped.replace({"ph.fixed_point.relaxation": 1.0}).digest()
    with pytest.raises(TypeError, match="ph.fixed_point.relaxation"):
        shipped.replace({"ph.fixed_point.relaxation": "0.7"})
    with pytest.raises(TypeError, match="ph.trunk_terms"):
        shipped.replace({"ph.trunk_terms": 1})   # a bool is not a number, nor a number a bool
    with pytest.raises(TypeError, match="ph.fixed_point.max_sweeps"):
        shipped.replace({"ph.fixed_point.max_sweeps": 40.5})
    with pytest.raises(TypeError, match="sulfur.tdp_own_decay"):
        shipped.replace({"sulfur.tdp_own_decay": ps.TdpCysteineTrap()})   # another holder
    with pytest.raises(TypeError, match="sulfur.ddp_formation_ph.slope"):
        shipped.replace({"sulfur.ddp_formation_ph": ps.DdpFormationPh(law="slope", slope="0.5")})
    with pytest.raises(ValueError, match="outside its band"):
        shipped.replace({"ph.fixed_point.relaxation": 0.2})
    # field bands: the dataclass fields that fits and feasibility draws sample
    with pytest.raises(ValueError, match=r"sulfur.nf_hydroxide_decay.mgo_share = 1.5"):
        shipped.replace({"sulfur.nf_hydroxide_decay": dataclasses.replace(ps.NfHydroxideDecay(), mgo_share=1.5)})
    with pytest.raises(ValueError, match=r"lipid.q10_assumption.default"):
        shipped.replace({"lipid.q10_assumption": dataclasses.replace(pl.Q10Assumption(), default=4.0)})
    with pytest.raises(ValueError, match=r"aa_pk"):
        shipped.replace({"sulfur.pentose_amine_entry": ps.PentoseAmineEntry(aa_law="saturating", aa_pk=8.0)})
    for value in (1.6, 3.3):   # the band is inclusive: CR-02's alternatives sit on its edges
        shipped.replace({"lipid.q10_assumption": dataclasses.replace(pl.Q10Assumption(), default=value)})
    banded = {name for name, key in config.KEY_BY_NAME.items() if key.bands}
    assert {"sulfur.nf_hydroxide_decay", "sulfur.pentose_amine_entry", "sulfur.cys_amadori_route",
            "sulfur.tdp_own_decay", "lipid.q10_assumption"} <= banded
    for name in banded:   # every band names a real field and holds the shipped value
        key = config.KEY_BY_NAME[name]
        for field_name, (lo, hi) in key.bands.items():
            assert lo <= getattr(key.default, field_name) <= hi, (name, field_name)


def test_digest_is_a_stable_sha256_unlike_hash():
    shipped, flipped = config.Config.shipped(), _flipped()
    assert len(shipped.digest()) == 64 and shipped.digest() == config.Config.from_globals().digest()
    assert flipped.digest() != shipped.digest()
    assert pickle.loads(pickle.dumps(flipped)).digest() == flipped.digest()
    for name in config.KEY_BY_NAME:   # every key moves the digest
        assert shipped.replace({name: FLIPPED[name]}).digest() != shipped.digest(), name
    # another process, another string-hash salt: the same digest
    code = ("from src.kinetic_core import config; s = config.Config.shipped(); "
            "print(s.digest(), s.replace({'ph.coupling': 'inline', 'ph.fixed_point.relaxation': 0.5}).digest())")
    env = dict(os.environ, PYTHONHASHSEED="12345")
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    assert out.stdout.split() == [shipped.digest(),
                                  shipped.replace({"ph.coupling": "inline", "ph.fixed_point.relaxation": 0.5}).digest()]


def test_build_phase_keys_are_read_only_facts(monkeypatch):
    shipped = config.Config.shipped()
    assert set(config.BUILD_KEYS) == {"dicarbonyl.shipped_b39", "methionine.shipped_b22", "methionine.shipped_b22b",
                                      "proline.shipped_b24", "proline.shipped_b24b",
                                      "uncertainty.b55_single_temperature_ea_prior",
                                      "uncertainty.b56_structural_discrepancy"}
    before = config._snapshot()
    for name in config.BUILD_KEYS:
        changed = shipped.replace({name: FLIPPED[name]})   # a Config may hold it ...
        with pytest.raises(config.BuildPhaseKeyError, match=name.replace(".", r"\.")):
            config.install(changed)                        # ... but never install it
        with pytest.raises(config.BuildPhaseKeyError, match=name.replace(".", r"\.")):
            with config.applied(changed):
                pass
        with pytest.raises(config.BuildPhaseKeyError, match=name.replace(".", r"\.")):
            config.update({name: FLIPPED[name], "ph.coupling": "inline"})
        after = config._snapshot()
        assert all(after[b] is before[b] for b in before), name   # refused before anything was set
    # a run-time change installs with the build keys at their shipped values, which install never touches
    from src.kinetic_core import parameters_methionine as pm
    monkeypatch.setattr(pm, "METHIONINE_SHIPPED", True)   # what a process that imported it so would hold
    held = config.Config.from_globals()
    assert held.diff() == {"methionine.shipped_b22": (False, True)}
    with pytest.raises(config.BuildPhaseKeyError, match="methionine.shipped_b22"):
        config.install(held)
    monkeypatch.undo()
    with config.applied(_flipped_runtime()):
        assert pm.METHIONINE_SHIPPED is False


def test_update_sets_every_binding_from_the_current_state_and_returns_the_config():
    from src.kinetic_core import lipid
    import src.kinetic_core as package

    with config.applied(config.Config.shipped()):
        first = config.update({"ph.coupling": "inline"})
        q10 = dataclasses.replace(pl.Q10Assumption(), default=2.0)
        second = config.update({"lipid.q10_assumption": q10})
        assert second == config.Config.from_globals() == first.replace({"lipid.q10_assumption": q10})
        assert pl.Q10_ASSUMPTION is lipid.Q10_ASSUMPTION is package.Q10_ASSUMPTION is q10   # every binding
        with pytest.raises(ValueError, match="outside its band"):
            config.update({"lipid.q10_assumption": dataclasses.replace(q10, default=5.0)})
        assert config.Config.from_globals() == second   # a refused update sets nothing


def test_applied_sets_every_binding_and_restores_the_very_objects_also_on_an_exception():
    before = config._snapshot()
    flipped = _flipped_runtime()
    with config.applied(flipped):
        assert config.Config.from_globals() == flipped
        from src.kinetic_core import lipid, sulfur
        assert sulfur.PH_COUPLING == "inline" and lipid.DECADIENAL_BREAKAGE.enabled is False
        assert pl.DECADIENAL_BREAKAGE.enabled is False   # both bindings
    after = config._snapshot()
    assert after.keys() == before.keys() and all(after[b] is before[b] for b in before)
    with pytest.raises(RuntimeError, match="inside"):
        with config.applied(flipped):
            raise RuntimeError("inside")
    after = config._snapshot()
    assert all(after[b] is before[b] for b in before)
    assert config.Config.from_globals() == config.Config.shipped()


def test_from_globals_refuses_bindings_that_disagree(monkeypatch):
    from src.kinetic_core import lipid

    monkeypatch.setattr(lipid, "DECADIENAL_BREAKAGE", dataclasses.replace(pl.DECADIENAL_BREAKAGE, enabled=False))
    with pytest.raises(ValueError, match="lipid.decadienal_breakage"):
        config.Config.from_globals()


def test_map_ordered_installs_a_jobs_config_serially_and_in_forked_workers(monkeypatch):
    flipped = _flipped_runtime()
    read = lambda _x: config.Config.from_globals()   # noqa: E731
    monkeypatch.setenv("MAILLARD_WORKERS", "1")
    assert parallel.map_ordered(read, [0, 1], config=flipped) == [flipped, flipped]
    assert parallel.map_ordered(read, [0, 1], n_workers=2, config=flipped) == [flipped, flipped]
    mixed = [parallel.Configured(0, config.Config.shipped()), 1]
    assert parallel.map_ordered(read, mixed, n_workers=2, config=flipped) == [config.Config.shipped(), flipped]
    assert config.Config.from_globals() == config.Config.shipped()   # the caller's state is untouched


@pytest.mark.skipif("fork" not in mp.get_all_start_methods(), reason="needs fork to compare against")
def test_a_spawned_worker_given_the_config_reproduces_the_forked_one_with_every_runtime_switch_flipped():
    flipped = _flipped_runtime()
    seeds = [0, 1]
    with config.applied(flipped):   # the forked workers inherit the state: today's mechanism
        forked = parallel.map_ordered(support.engine_config_probe, seeds, n_workers=2)
    serial_shipped = [support.engine_config_probe(s) for s in seeds]
    with mp.get_context("spawn").Pool(2) as pool:
        spawned = pool.map(parallel.run_pickled,
                           [(support.engine_config_probe, parallel.Configured(s, flipped)) for s in seeds])
        spawned_bare = pool.map(support.engine_config_probe, seeds)
    assert [r["config"] for r in forked] == [flipped, flipped]
    assert spawned == forked   # every number, bit for bit
    # the trap the Config closes: a spawned worker without it computes with the defaults, silently
    assert spawned_bare == serial_shipped and spawned_bare != forked
    assert [r["config"] for r in spawned_bare] == [config.Config.shipped()] * 2
    for numbers in ("rate_constants", "trunk_ph_exponents", "concentrations", "ph_series"):
        assert spawned_bare[0][numbers] != forked[0][numbers], numbers   # the probe reads the switches


@pytest.mark.skipif("fork" not in mp.get_all_start_methods(), reason="needs fork to compare against")
def test_the_persistent_pool_is_bit_identical_to_fork_and_never_stale(monkeypatch):
    """A3 phase 2: MAILLARD_POOL=persistent keeps one spawned pool for the process; every job carries the caller's
    Config read at the call, so the same workers follow the caller's state from call to call."""
    flipped, seeds = _flipped_runtime(), [0, 1, 2]
    monkeypatch.setenv("MAILLARD_POOL", "fork_per_call")
    with config.applied(flipped):
        forked = parallel.map_ordered(support.engine_config_probe, seeds, n_workers=2)
    forked_shipped = parallel.map_ordered(support.engine_config_probe, seeds, n_workers=2)
    monkeypatch.setenv("MAILLARD_POOL", "persistent")
    try:
        with config.applied(flipped):
            persistent = parallel.map_ordered(support.engine_config_probe, seeds, n_workers=2)
        persistent_shipped = parallel.map_ordered(support.engine_config_probe, seeds, n_workers=2)   # same pool
        pids = set(parallel.map_ordered(support.worker_pid, range(8), n_workers=2))
        again = set(parallel.map_ordered(support.worker_pid, range(8), n_workers=2))
        assert persistent == forked and persistent_shipped == forked_shipped   # every number, bit for bit
        assert persistent != persistent_shipped
        assert len(pids) <= 2 and pids == again and os.getpid() not in pids   # reused, not re-forked
        # a closure cannot be pickled: it runs on the per-call fork, which carries it
        offset = 3
        assert parallel.map_ordered(lambda x: x + offset, [1, 2], n_workers=2) == [4, 5]
    finally:
        parallel.shutdown_persistent_pool()
    monkeypatch.setenv("MAILLARD_POOL", "persistent_typo")
    with pytest.raises(ValueError, match="MAILLARD_POOL"):
        parallel.map_ordered(support.worker_pid, [0, 1], n_workers=2)


def _legacy_joint_install(joint, job):
    """fit_joint_deoxyosone_thiol.install as it was before A3 (origin/main 7b1bc7cd), for the equivalence below."""
    joint.sm3.install(job)
    _ps = ps
    _ps.STRECKER_SULFUR = dataclasses.replace(_ps.STRECKER_SULFUR,
                                              h2s_share=float(job.get("strecker_share", joint.STRECKER_SHARE)))
    _ps.MGO_CYSTEINE_STRECKER = dataclasses.replace(_ps.MGO_CYSTEINE_STRECKER,
                                                    k_l_per_mmol_min_100c_ph5=float(job.get("mgo_strecker_k", 0.0)))
    _ps.TDP_CYSTEINE_TRAP = dataclasses.replace(_ps.TDP_CYSTEINE_TRAP,
                                                k_l_per_mmol_min_100c_ph5=float(job.get("tdp_cys_trap_k", 0.0)))
    _ps.PENTOSE_AMINE_ENTRY = _ps.PentoseAmineEntry(**job.get("entry", {}))
    _ps.DDP_FORMATION_PH = _ps.DdpFormationPh(**job.get("ddp_law", {"law": "base"}))
    _default = _ps.CysAmadoriRoute()
    _split = job.get("cama_split", {})
    _ps.CYS_AMADORI_ROUTE = dataclasses.replace(
        _ps.CYS_AMADORI_ROUTE, split_law=str(_split.get("split_law", _default.split_law)),
        dpo_share=float(_split.get("dpo_share", _default.dpo_share)),
        split_ph_c=float(_split.get("split_ph_c", _default.split_ph_c)), dpo_share_ref_ph=_default.dpo_share_ref_ph)
    _ps.TDP_OWN_DECAY = _ps.TdpOwnDecay(**job.get("tdp_exit", {}))
    ph_state.TRUNK_SOLUTES_IN_BALANCE = bool(job.get("trunk_solutes", False))


def test_the_joint_fits_install_through_a_config_sets_what_it_set_before(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts" / "generators"))
    import fit_joint_deoxyosone_thiol as joint

    with config.applied(config.Config.shipped()):   # every binding back afterwards (the package's too)
        for arm in joint.ARMS:
            monkeypatch.setattr(joint, "TDP_EXIT_STATIC", joint.tdp_exit_static(joint.ARMS[arm]))
            keys, _kinds, _lo, _hi, start = joint.coordinates(False, joint.ARMS[arm])
            job = dict(joint.trial(list(start), keys), strecker_share=0.042)
            for variant in (job, dict(job, trunk_solutes=True, cama_split={"split_law": "base_saturating",
                                                                           "dpo_share": 0.4, "split_ph_c": 6.5})):
                for before in (config.Config.shipped(), _flipped_runtime()):
                    config.install(before)
                    _legacy_joint_install(joint, variant)
                    legacy = config.Config.from_globals()
                    config.install(before)
                    joint.install(variant)
                    assert config.Config.from_globals() == legacy, arm
                    assert legacy != before


def _captured_tasks(monkeypatch, module, jobs, rows, keys, point):
    """The tasks a fit's ``make_raw_many`` hands to ``parallel.map_ordered`` (none is run)."""
    captured = []

    def fake_map(fn, items, *args, **kwargs):
        captured.append((fn, list(items)))
        return [{"failed": "not run"} for _ in captured[-1][1]]

    monkeypatch.setattr(module.parallel, "map_ordered", fake_map)
    monkeypatch.setattr(module, "FAILURES", [])
    module.make_raw_many(jobs, keys)([point], rows)
    (fn, tasks), = captured
    return fn, tasks


def test_the_fits_tasks_carry_their_config_built_in_the_caller(monkeypatch):
    """A3 phase 2: the joint fit's and sub-model 3's tasks are ``parallel.Configured`` with the Config the legacy
    ``run_task`` installed in a worker forked from the caller (shipped + the inline coupling + sub-model 3's routes +
    the job's switches, the arm's statics included), and the worker function installs nothing itself. So the
    module's statics (STRECKER_SHARE, TDP_EXIT_STATIC ...) and the coupling are read in the caller, never inherited."""
    monkeypatch.syspath_prepend(str(ROOT / "scripts" / "generators"))
    import fit_deoxyosone_submodel as sm3
    import fit_joint_deoxyosone_thiol as joint

    rows = [types.SimpleNamespace(table="thiol_submodel", pot=p, ledger_id="x", species="MFT") for p in ("a", "b")]
    jobs = {f"thiol_submodel:{p}": {"pot": p, "strecker_share": 0.08} for p in ("a", "b")}
    caller = config.Config.shipped().replace({"ph.coupling": "inline"})
    with config.applied(caller):
        for arm in ("B7ASST", "R3En", "E"):
            monkeypatch.setattr(joint, "TDP_EXIT_STATIC", joint.tdp_exit_static(joint.ARMS[arm]))
            monkeypatch.setattr(joint, "DDP_LAW_STATIC", {"law": "none"} if "ddp_law_none" in joint.ARMS[arm] else {})
            keys, _kinds, _lo, _hi, start = joint.coordinates(False, joint.ARMS[arm])
            fn, tasks = _captured_tasks(monkeypatch, joint, jobs, rows, keys, start)
            assert fn is joint.run_pot and len(tasks) == 2
            for task in tasks:
                assert isinstance(task, parallel.Configured) and task.config["ph.coupling"] == "inline"
                config.install(caller)
                _legacy_joint_install(joint, task.item)
                assert task.config == config.Config.from_globals(), arm
                config.install(caller)
            assert tasks[0].config["sulfur.strecker_sulfur"].h2s_share == 0.08   # the job's, not the module's
        keys, *_rest, start = sm3.coordinates("A3")
        sm3_jobs = {p: {"pot": p} for p in ("a", "b")}
        fn, tasks = _captured_tasks(monkeypatch, sm3, sm3_jobs, rows, keys, start)
        assert fn is sm3.run_pot
        for task in tasks:
            for name in sm3.ROUTE_KEYS:
                setattr(ps, config.KEY_BY_NAME[name].bindings[0][1],
                        dataclasses.replace(config.KEY_BY_NAME[name].default, enabled=True))
            ps.NF_HYDROXIDE_DECAY = dataclasses.replace(ps.NF_HYDROXIDE_DECAY, mgo_share=task.item["mgo_share"])
            ps.PENTOSE_MGO_ROUTES = ps.PentoseMgoRoutes(**task.item["routes"])   # the legacy sm3.install
            assert task.config == config.Config.from_globals()
            config.install(caller)
    assert config.Config.from_globals() == config.Config.shipped()


def test_the_probe_arm_config_is_what_the_command_line_set(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts" / "generators"))
    import build_submodel_probe as probe

    args = types.SimpleNamespace(ph_coupling="inline", trunk_solutes_in_balance=True, trunk_ph_dynamic=True,
                                 switch_on=["TTCA_FORMATION", "NF_HYDROXIDE_DECAY"], nf_mgo_share=0.05)
    cfg = probe.arm_config(args, config.Config.shipped())
    assert cfg.diff() == {
        "ph.coupling": ("fixed_point_nodes", "inline"), "ph.trunk_solutes_in_balance": (False, True),
        "ph.trunk_terms": (False, True),
        "sulfur.ttca_formation": (ps.TtcaFormation(), dataclasses.replace(ps.TtcaFormation(), enabled=True)),
        "sulfur.nf_hydroxide_decay": (ps.NfHydroxideDecay(),
                                      dataclasses.replace(ps.NfHydroxideDecay(), enabled=True, mgo_share=0.05))}


def test_the_choice_robustness_q10_sweep_sets_every_binding(monkeypatch):
    """The bug the review found: build_choice_robustness set only parameters_lipid.Q10_ASSUMPTION, so the rate used
    the alternative while lipid.integrate_lipid reported ``q10_used`` (and the decades of extrapolation) from its own,
    re-imported binding, still the shipped 2.449. Now the sweep installs a Config: both agree, and nothing is left."""
    monkeypatch.syspath_prepend(str(ROOT / "scripts" / "generators"))
    import build_choice_robustness as robustness
    from src.kinetic_core import lipid
    from src.kinetic_core.engine import core_lipid_model

    branch, composition = core_lipid_model()
    charge = lipid.charge_from_carrier(pl.LIPID_CARRIERS["frankel_pure_hydroperoxide"], composition)

    def measure():
        segment = lipid.integrate_lipid(charge, [(1.0, 180.0)], branch).metadata["segments"][0]
        return segment["q10_used"], segment["k_LOOH_decomp_per_min"], segment["decades_of_extrapolation_from_25C"]

    alternatives = [{"name": "moist_matrix", "value": 1.6}, {"name": "bulk_oil", "value": 3.3}]
    out = robustness.q10_sweep(alternatives, measure)
    for a in alternatives:
        used, k, decades = out[a["name"]]
        q10 = dataclasses.replace(pl.Q10Assumption(), default=a["value"])
        assert used == a["value"]                                              # the report ...
        assert k == pl.k_looh_decomp_per_min(180.0, a["value"])                 # ... names the rate's Q10
        assert decades == q10.decades_of_extrapolation(180.0)
    assert out["moist_matrix"][1] != out["bulk_oil"][1]
    assert config.Config.from_globals() == config.Config.shipped()            # every binding restored
