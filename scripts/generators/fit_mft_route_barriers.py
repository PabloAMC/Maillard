#!/usr/bin/env python
"""
E6 + E11 -- MFT'S TWO INTACT-C5 ROUTES WITH PER-ROUTE FORMATION BARRIERS: stage 0 and the refit (2026-11-01).

Pre-registered, binding: ``results/validation/kinetic_core_mft_route_barriers_prereg.md`` (written before any
number here existed). Origin: ``docs/reference/pathway_atlas/ENGINE_VS_LITERATURE.md`` sec. 3 rows E6 / E11 and the
red-team amendment, sec. 5 item 2. Nothing here ships and nothing in ``src/`` is edited.

THE CHANGE (prereg sec. 1)
==========================
The per-STEP barriers go through the B62 hook: ``parameters_sulfur.FurfuralBarriers.barriers_kj_mol`` takes any
fitted key and ``with_fitted_sulfur`` reads it through ``measured_ea_overrides()``, holding k(145 C) at the fit's
log10 k, so a barrier moves only the off-145 C rates. Every member's engine state is ONE ``Config`` (A3,
``src/kinetic_core/config.py``) that differs from the shipped one in the ``sulfur.furfural_barriers`` key only, and
it is installed with ``config.applied`` around the work; no engine global is assigned by hand.

* B62's furfural barriers (k_tdp_fur 147, k_fur_decay 46.5 kJ/mol) in every member;
* k_dpo_nf at 83.8 kJ/mol (Lau 2003 / Pandit 2006 inverse-variance mean; 93.0 as the sensitivity);
* k_nf_mft PROFILED at {64.08 (the lump), 100, 133} kJ/mol, not freed;
* two coordinates freed beside B60's 23: log10 k_dpo_nf (145 C) and ONE log10 offset on k_ddp_mft and
  k_ddp_mft_hs together (the identifiable product; k_dpo_ddp stays frozen).

Members: N1 (B62 only, B60's free set), A64 / A100 / A133 (start 0, the B60 optimum), and, for the lowest-cost
converged A: A*-s1 (start 1, B8's perturbation protocol) and A*-93 (k_dpo_nf at 93.0).

THE OPTIMISER (prereg sec. 3)
=============================
B60's 52-row objective, its declared sigmas FIXED, its quick setting (``generate_kinetic_core_b60_fit``), its
bands. The fit runs through ``src/kinetic_core/fitting.fit`` with the converged-fit rules the prereg names
(EX-024): a budget of 1200 evaluations (Jacobian columns counted), ftol 1e-4, and only status 2 / 3 counts as
converged. ``fitting.fit`` uses B8's settings except the finite-difference step, which is B8's Laplace's absolute
step per coordinate kind: B8's relative ``diff_step`` collapses at a coordinate near zero (fitting.py's docstring),
and the new DDP offset STARTS at zero. Its Jacobian points run in parallel under ``MAILLARD_WORKERS``. A
budget-stopped member is re-launched warm from its last point once (``--warm``); it is not read before that.

THE TESTS (prereg sec. 4, frozen)
=================================
Route tests at a vector are the read-only probes of the prereg (sec. 10), ported unchanged: the same pots, time
grids and pH modes (clamped: the probes passed no drift, so a member's pH-drift coordinates do not enter them).
D3 is scored through the panel path (``scoring.score_panel`` on the two Hofmann bundles) with the engine reading
the member's vector from a temporary copy of the shipped report (the B62 ship rule's pattern).

Usage (every command through ./scripts/docker_maillard.sh run "..."):
    python scripts/generators/fit_mft_route_barriers.py --vcheck            # the prereg sec. 8 V-checks
    python scripts/generators/fit_mft_route_barriers.py --stage0            # incumbent constants, no refit
    python scripts/generators/fit_mft_route_barriers.py --member A100       # one member's fit
    python scripts/generators/fit_mft_route_barriers.py --member A100 --warm   # the one warm relaunch
    python scripts/generators/fit_mft_route_barriers.py --queue             # every member, two at a time, then --verdict
    python scripts/generators/fit_mft_route_barriers.py --verdict           # score the members, ledger, sec. 5 rule
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b60_fit as B60  # noqa: E402  (B60's rows, physics, numerics and quick setting)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import config as engine_config  # noqa: E402
from src.kinetic_core import engine, fitting, sulfur  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

B23, B8 = B60.B23, B60.B8
B24 = B8.B24
V = data_paths.VALIDATION_DIR
GENERATED_BY = "scripts/generators/fit_mft_route_barriers.py"
PREREG = V / "kinetic_core_mft_route_barriers_prereg.md"
STAGE0_OUT = V / "kinetic_core_mft_route_barriers_stage0.json"
REFIT_OUT = V / "kinetic_core_mft_route_barriers_refit.json"
MEMBER_DIR = V / "kinetic_core_mft_route_barriers_members"
LOG_DIR = ROOT / "scratch" / "mft_route_barriers"

#: The base is what ships (prereg sec. 1). If B64 (or anything else) has shipped, sec. 7 must be re-run first.
SHIPPED_REPORT = Path(engine._B2_FIT_REPORT)
BASE_REPORT = V / "kinetic_core_b60_fit_report.json"

# ---------------------------------------------------------------------------------------------------------------
# 1. The declared barriers (prereg sec. 1 table; NOT fitted)
# ---------------------------------------------------------------------------------------------------------------
B62_BARRIERS: Tuple[Tuple[str, float], ...] = tuple(ps.FurfuralBarriers().barriers_kj_mol)  # k_tdp_fur, k_fur_decay
NF_FORMATION_EA = 83.8          # k_dpo_nf: inverse-variance mean of Lau 2003 81.5 +- 3.2 and Pandit 2006 93.0 +- 6.4
NF_FORMATION_EA_SENSITIVITY = 93.0
#: k_nf_mft's profile points: the lump (B60's 64.08, read from the report), 100 and 133 kJ/mol
NF_MFT_POINTS: Dict[str, Optional[float]] = {"A64": None, "A100": 100.0, "A133": 133.0}
SOURCE = ("E6 + E11 prereg (kinetic_core_mft_route_barriers_prereg.md sec. 1): k_tdp_fur Schirle-Keller & Reineccius 1992 "
          "Table 2; k_fur_decay Yaghmur 2005; k_dpo_nf Lau 2003 / Pandit 2006 (inverse-variance mean 83.8, sensitivity "
          "93.0); k_nf_mft profiled in the thiol-assembly class [40, 133] (thiol prereg sec. 13)")

#: The converged-fit rules (prereg sec. 3, EX-024)
BUDGET = 1200
FTOL = 1e-4
CONVERGED_STATUS = (2, 3)
START1_COST_AGREEMENT = 2.0

#: The tied DDP-route coordinate (the joint fit's DDP_COORDINATE pair): one log10 offset on both keys
DDP_OFFSET = "offset_ddp_route"
DDP_KEYS = ("k_ddp_mft", "k_ddp_mft_hs")
NF_KEY = "k_dpo_nf"

#: The compensation ledger (prereg sec. 4.3): rows that may get worse by > 0.3 dex, in their shipped direction only
LEDGER: Dict[str, str] = {
    "zhou_pH7_MFT": "E5 (no buffer catalysis) and E2 (no pH on the sulfur lane's trunk steps)",
    "zhou_pH7_FFT": "E5 (no buffer catalysis) and E2 (no pH on the sulfur lane's trunk steps)",
    "zhou_pH7_MFT_over_FFT": "E5 (no buffer catalysis) and E2 (no pH on the sulfur lane's trunk steps)",
    "zhou_pH7_FUR_arp_alone": "E5 (no buffer catalysis) and E2 (no pH on the sulfur lane's trunk steps)",
}
GUARD_DEX = 0.3


def lumped_ea() -> float:
    return float(json.loads(BASE_REPORT.read_text(encoding="utf-8"))["frozen_parameters"]["lumped_formation_Ea_kJ_mol"])


# ---------------------------------------------------------------------------------------------------------------
# 2. Members and their configurations
# ---------------------------------------------------------------------------------------------------------------
_MEMBER_RE = re.compile(r"^(N1|A(64|100|133))(-s1|-93)?$")


@dataclasses.dataclass(frozen=True)
class Member:
    name: str
    barriers: Tuple[Tuple[str, float], ...]
    frees_route_coordinates: bool
    start: int
    nf_mft_ea: Optional[float]
    dpo_nf_ea: Optional[float]


def member(name: str) -> Member:
    """The member a name declares (prereg sec. 1): N1, A64 / A100 / A133, and A<point>-s1 / A<point>-93."""
    m = _MEMBER_RE.match(name)
    if not m or (m.group(1) == "N1" and m.group(3)):
        raise ValueError(f"{name!r} is not a pre-registered member")
    if name == "N1":
        return Member(name, B62_BARRIERS, False, 0, None, None)
    point = "A" + m.group(2)
    nf_mft = NF_MFT_POINTS[point]
    nf_mft = lumped_ea() if nf_mft is None else nf_mft
    dpo_nf = NF_FORMATION_EA_SENSITIVITY if m.group(3) == "-93" else NF_FORMATION_EA
    barriers = B62_BARRIERS + ((NF_KEY, dpo_nf), ("k_nf_mft", float(nf_mft)))
    return Member(name, barriers, True, 1 if m.group(3) == "-s1" else 0, float(nf_mft), dpo_nf)


def barrier_config(barriers: Sequence[Tuple[str, float]]) -> engine_config.Config:
    """The shipped Config with the B62 hook carrying ``barriers`` (enabled); every other key shipped."""
    hook = ps.FurfuralBarriers(enabled=True, barriers_kj_mol=tuple((k, float(v)) for k, v in barriers), source=SOURCE)
    return engine_config.Config.shipped().replace({"sulfur.furfural_barriers": hook})


def member_config(name: str) -> engine_config.Config:
    return barrier_config(member(name).barriers)


#: Stage 0 (prereg sec. 2): B60's frozen vector, no refit, under each configuration. ``shipped`` is the baseline.
def stage0_configs() -> Dict[str, engine_config.Config]:
    return {"shipped": engine_config.Config.shipped(), "N1_B62_alone": member_config("N1"),
            "A64": member_config("A64"), "A100": member_config("A100"), "A133": member_config("A133")}


# ---------------------------------------------------------------------------------------------------------------
# 3. The coordinates
# ---------------------------------------------------------------------------------------------------------------
def incumbent_vector() -> np.ndarray:
    """B60's frozen optimum, clipped into its bands (B60.incumbent_vector reads B51's; this is B60's own)."""
    lower, upper = B60.full_bounds()
    return np.clip(frozen_vector(json.loads(BASE_REPORT.read_text(encoding="utf-8"))), lower, upper)


def _kind(key: str) -> str:
    if key.startswith("Ea_"):
        return "Ea"
    if key == "ph_acid_yield_per_sink_event":
        return "acid_yield"
    if key == "ph_arp_secondary_ammonium_pKa":
        return "pka"
    return "log10k"


def coordinates(name: str) -> List[Tuple[str, str, float, float]]:
    """``[(label, kind, lower, upper), ...]``: B60's 23 free keys, then (A members) log10 k_dpo_nf and the DDP
    offset. The offset's band keeps BOTH keys inside their own bounds ("bounds those of each key")."""
    lower, upper = B60.full_bounds()
    base = incumbent_vector()
    out = [(k, _kind(k), float(lower[B60.ALL_KEYS.index(k)]), float(upper[B60.ALL_KEYS.index(k)]))
           for k in B60.FREE_KEYS]
    if member(name).frees_route_coordinates:
        i = B60.ALL_KEYS.index(NF_KEY)
        out.append((NF_KEY, "log10k", float(lower[i]), float(upper[i])))
        idx = [B60.ALL_KEYS.index(k) for k in DDP_KEYS]
        out.append((DDP_OFFSET, "log10k", float(max(lower[i] - base[i] for i in idx)),
                    float(min(upper[i] - base[i] for i in idx))))
    return out


def to_full(theta: Sequence[float], name: str, base: Optional[np.ndarray] = None) -> np.ndarray:
    base = incumbent_vector() if base is None else base
    x = base.copy()
    labels = [c[0] for c in coordinates(name)]
    for label, value in zip(labels, theta):
        if label == DDP_OFFSET:
            for k in DDP_KEYS:
                i = B60.ALL_KEYS.index(k)
                x[i] = base[i] + float(value)
        else:
            x[B60.ALL_KEYS.index(label)] = float(value)
    return x


def theta_of(x_full: np.ndarray, name: str, base: Optional[np.ndarray] = None) -> np.ndarray:
    base = incumbent_vector() if base is None else base
    out = []
    for label, *_ in coordinates(name):
        if label == DDP_OFFSET:
            i = B60.ALL_KEYS.index(DDP_KEYS[0])
            out.append(float(x_full[i] - base[i]))
        else:
            out.append(float(x_full[B60.ALL_KEYS.index(label)]))
    return np.array(out, dtype=float)


def start_theta(name: str) -> np.ndarray:
    """Start 0: B60's optimum (the two new coordinates at their shipped values). Start 1: B8's perturbation protocol
    (``generate_kinetic_core_b8_fit.start_vector``: a uniform local draw, 2 decades / 40 kJ / 0.15 / 1 pKa, seed
    SEED_BASE + SEED_STRIDE x start), extended in coordinate order to the two new log10 coordinates."""
    coords = coordinates(name)
    theta = theta_of(incumbent_vector(), name)
    start = member(name).start
    if start == 0:
        return theta
    rng = np.random.default_rng(B8.SEED_BASE + B8.SEED_STRIDE * start)
    for j, (label, _kind_, lo, hi) in enumerate(coords):
        half = B8._local_half_width(label)
        theta[j] = float(rng.uniform(max(lo, theta[j] - half), min(hi, theta[j] + half)))
    return theta


def problem(name: str) -> fitting.Problem:
    coords = coordinates(name)
    base = incumbent_vector()
    row_ids = tuple(row["id"] for row in B23.ACTIVE_FIT_ROWS)

    def residuals(theta, _rows):
        return B60.residual_vector(to_full(theta, name, base), True)

    return fitting.Problem(keys=tuple(c[0] for c in coords), kinds=tuple(c[1] for c in coords),
                           lower=np.array([c[2] for c in coords]), upper=np.array([c[3] for c in coords]),
                           rows=row_ids, residuals=residuals)


# ---------------------------------------------------------------------------------------------------------------
# 4. Scoring a vector: the 52 rows (raw and weighted) and the frozen route tests
# ---------------------------------------------------------------------------------------------------------------
def rows_at(x_full: np.ndarray) -> Dict[str, Any]:
    """B60's 52 rows at ``x_full`` under the INSTALLED configuration: prediction, raw residual (dex; pH units for a
    pH endpoint; one-sided rows clipped as the objective clips them) and the objective's weighted residual."""
    parameters = B23.build_parameters(x_full)
    drift = B23.unpack(x_full)[3]
    predicted, _ = B23.observables(parameters, True, drift)
    sigma_ph = B24.SIGMA_LOG_REFERENCE / float(B8.WEIGHT_E)
    out, weighted = {}, []
    for row in B23.ACTIVE_FIT_ROWS:
        p, t = float(predicted[row["id"]]), float(row["target"])
        if not np.isfinite(p):
            raw, w = float("nan"), 25.0
        elif row["kind"] == "ph_endpoint":
            raw = p - t
            w = raw / sigma_ph
        else:
            lr = math.log10((p + B23.FLOOR) / (t + B23.FLOOR))
            raw = max(0.0, lr) if row["kind"] in ("ceiling", "peak_fraction_ceiling") else (
                min(0.0, lr) if row["kind"] == "floor" else lr)
            w = raw / float(row["sigma_log"])
        w = float(np.clip(w, -25.0, 25.0))
        weighted.append(w)
        out[row["id"]] = {"kind": row["kind"], "predicted": p, "target": t, "raw": raw, "weighted": w,
                          "units": "pH" if row["kind"] == "ph_endpoint" else "dex"}
    r = np.array(weighted)
    return {"rows": out, "cost": float(0.5 * r @ r), "n_rows": len(weighted)}


MFT_STEPS = ("r_ddp_mft", "r_ddp_mft_hs", "r_nf_mft", "r_ha_mp_mft", "r_hmp_mft")
PHOS05 = BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="probe")
CLAMPED = BufferSpec(kind="clamped")
NO_BUFFER = BufferSpec(kind="none")


def _run(P, init, t_c, minutes, ph, buf):
    times = sorted({0.0, *np.geomspace(1e-5, minutes, 300).tolist(), *np.linspace(0.0, minutes, 121).tolist()})
    r = sulfur.integrate_sulfur(P, t_c + 273.15, init, times, ph=ph, buffer_spec=buf)
    return r, sulfur.reaction_flux_integrals(r)


def _mft(fl) -> float:
    return float(sum(fl.get(k, 0.0) for k in MFT_STEPS))


def _share(fl, key="r_nf_mft") -> float:
    m = _mft(fl)
    return float(fl.get(key, 0.0) / m) if m else float("nan")


def fed_norfuraneol_tracer(run) -> float:
    """``fit_joint_deoxyosone_thiol.fed_norfuraneol_tracer`` (the joint fit's exact unlabelled share of Cerny's pot
    D), ported so this module does not import the joint fit's sub-model chain."""
    node_times, k_nodes = run.rate_nodes
    times = np.asarray(run.times_min, dtype=float)
    keys = [x.key for x in sulfur.FULL_REACTIONS]
    nu_nf = np.array([x.products.get("NF", 0) - x.reactants.get("NF", 0) for x in sulfur.FULL_REACTIONS], float)
    nu_mft = np.array([x.products.get("MFT", 0) for x in sulfur.FULL_REACTIONS], float)
    i = keys.index("r_nf_mft")
    nf = run.series("NF")
    rates = np.array([sulfur.sulfur_reaction_rates(run.concentrations[j], sulfur.rate_vector_at(node_times, k_nodes, t))
                      for j, t in enumerate(times)])
    production = rates @ np.where(nu_nf > 0, nu_nf, 0.0)
    g = np.where(nf > 0, production / np.maximum(nf, 1e-300), 0.0)
    f = np.exp(-np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(times))]))
    total = float(np.trapezoid(rates @ nu_mft, times))
    return float(np.trapezoid(rates[:, i] * f, times)) / total if total > 0 else float("nan")


def route_tests(P: Mapping[str, Any]) -> Dict[str, Any]:
    """The prereg's route tests (sec. 4) at one parameter set: the probes of sec. 10, unchanged."""
    out: Dict[str, Any] = {}
    # D1 / D5: Cerny 2003 pots D and C (0.5 M K phosphate pH 5, 95 C, 240 min)
    rD, fD = _run(P, {"NF": 145.0, "PENT": 145.0, "Cys": 96.6}, 95.0, 240.0, 5.0, PHOS05)
    _rC, fC = _run(P, {"PENT": 290.0, "Cys": 96.6}, 95.0, 240.0, 5.0, PHOS05)
    out["D1_cerny_potD_fed_nf_share"] = fed_norfuraneol_tracer(rD)
    out["D1_cerny_potD_nf_route_flux_share"] = _share(fD)
    out["D5_cerny_potC_mp3p_over_mft"] = float(fC.get("r_nf_mp3p", 0.0) / _mft(fC))
    out["cerny_potC_nf_route_flux_share"] = _share(fC)
    # D2: Givaudan (Cys : carbonyl 1 : 1, 37 % water, pH 5 clamped, 100 C, 90 min); 2 M declared, 0.3 M sensitivity
    for c in (2000.0, 300.0):
        _, fR = _run(P, {"PENT": c, "Cys": c}, 100.0, 90.0, 5.0, CLAMPED)
        _, fN = _run(P, {"NF": c, "Cys": c}, 100.0, 90.0, 5.0, CLAMPED)
        out[f"D2_givaudan_ribose_over_nf_mft_{int(c)}mM"] = _mft(fR) / _mft(fN)
        out[f"givaudan_ribose_nf_route_share_{int(c)}mM"] = _share(fR)
    # D4: the in-situ NF-route share in Hofmann's ribose + Cys pot (145 C, 20 min, pH 5, 0.5 M phosphate); 100 / 120 C
    _, fRC = _run(P, {"PENT": 100.0, "Cys": 33.0}, 145.0, 20.0, 5.0, PHOS05)
    out["D4_hofmann145_ribose_nf_route_share"] = _share(fRC)
    out["hofmann145_ribose_ddp_route_share"] = float((fRC.get("r_ddp_mft", 0) + fRC.get("r_ddp_mft_hs", 0)) / _mft(fRC))
    for t_c, minutes in ((100.0, 60.0), (120.0, 40.0)):
        _, f = _run(P, {"PENT": 100.0, "Cys": 33.0}, t_c, minutes, 5.0, PHOS05)
        out[f"ribose_cys_{int(t_c)}C_nf_route_share"] = _share(f)
    # H1: Liu 2023b (CAAS), 168 C, 0.5 M K3PO4 pH 5.88; volume unstated: 10 mL declared, 5 / 20 swept
    for vol in (5.0, 10.0, 20.0):
        rib, cys = 0.13055 / (vol / 1000), 1.6342 / (vol / 1000)
        for minutes in (20.0, 60.0):
            _, f = _run(P, {"PENT": rib, "Cys": cys}, 168.0, minutes, 5.88, PHOS05)
            m = _mft(f)
            out[f"H1_liu168_{int(vol)}mL_{int(minutes)}min"] = {
                "fragment_C2C3_share": float(f.get("r_ha_mp_mft", 0.0) / m), "nf_share": _share(f),
                "ddp_share": float((f.get("r_ddp_mft", 0) + f.get("r_ddp_mft_hs", 0)) / m)}
    # H2: Wang 2012 (BTBU), unbuffered pH 4.8, 132 C, 90 min; pot C as pre-registered (Gly 12.1), the dossier's
    # Gly 9.80 as a sensitivity; pot A (GSH 7.42) reported. Pot B carries Glu, which is not an engine species.
    pots = {"potC": {"PENT": 26.0, "Cys": 12.1, "Gly": 12.1, "THI": 1.19},
            "potC_dossier_gly9.8": {"PENT": 26.2, "Cys": 12.15, "Gly": 9.80, "THI": 1.19},
            "potA_gsh": {"PENT": 26.2, "GSH": 7.42, "THI": 1.19}}
    for label, init in pots.items():
        _, f = _run(P, init, 132.0, 90.0, 4.8, NO_BUFFER)
        m = _mft(f)
        out[f"H2_wang2012_{label}"] = {
            "xylose_share_of_mft": float(1.0 - f.get("r_hmp_mft", 0.0) / m) if m else None,
            "nf_share": float(f.get("r_nf_mft", 0.0) / m) if m else None,
            "ddp_share": float((f.get("r_ddp_mft", 0) + f.get("r_ddp_mft_hs", 0)) / m) if m else None,
            "thiamine_mft_flux_mM": float(f.get("r_hmp_mft", 0.0)), "mft_total_flux_mM": m}
    return out


HOFMANN_BUNDLES = ("hofmann1998_norfuraneol_cysteine_145C_20min_pH5", "hofmann1998_ribose_cysteine_145C_20min_pH5")
MFT_NAME = "2-Methyl-3-furanthiol (MFT)"


@contextlib.contextmanager
def engine_reads(x_full: Optional[np.ndarray]) -> Iterator[None]:
    """The panel path reads the fitted vector from the engine's report. For a member's vector, a temporary copy of
    the shipped report with its ``frozen_parameters`` replaced (outside the repo); the pointer is restored on exit.
    (The report pointer is not an engine switch, so it is not a Config key; the B62 ship rule's pattern.)"""
    if x_full is None:
        yield
        return
    report = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    report["frozen_parameters"] = B24._b23_schema_frozen(np.asarray(x_full, dtype=float))
    previous = engine._B2_FIT_REPORT
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "kinetic_core_b60_fit_report.json"
        path.write_text(json.dumps(report))
        engine._B2_FIT_REPORT = path
        try:
            yield
        finally:
            engine._B2_FIT_REPORT = previous


def hofmann_panel_ratio() -> Dict[str, float]:
    """D3 on the panel path: MFT(NF + Cys) / MFT(ribose + Cys) in ppb, predicted and measured."""
    from src.kinetic_core.scoring import PASS_BAND_LEVEL, panel_bundles, score_panel

    bundles = [(p, t) for p, t in panel_bundles()[0] if Path(p).stem in HOFMANN_BUNDLES]
    assert len(bundles) == 2, bundles
    got = {}
    for b in score_panel(bundles=bundles, pass_band=PASS_BAND_LEVEL)["benchmarks"]:
        for r in b["compounds"]:
            if r.get("compound") == MFT_NAME:
                got[b["benchmark_id"]] = (float(r["predicted"]), float(r["measured"]))
    nf, rb = got[HOFMANN_BUNDLES[0]], got[HOFMANN_BUNDLES[1]]
    return {"predicted": nf[0] / rb[0], "measured": nf[1] / rb[1]}


def score_vector(x_full: np.ndarray, cfg: engine_config.Config, *, panel_from_report: bool) -> Dict[str, Any]:
    """Everything a stage or a member is judged on, at one vector under one Config."""
    with engine_config.applied(cfg):
        rows = rows_at(x_full)
        tests = route_tests(B23.build_parameters(x_full))
        with engine_reads(x_full if panel_from_report else None):
            tests["D3_hofmann145_panel"] = hofmann_panel_ratio()
        barriers = {k: getattr(p, "ea_kj_mol", None) for k, p in B23.build_parameters(x_full).items()
                    if k in ("k_dpo_nf", "k_dpo_ddp", "k_nf_mft", "k_ddp_mft", "k_ddp_mft_hs", "k_tdp_fur",
                             "k_fur_decay", "k_cys_h2s", "k_nf_decay", "k_pent_dpo")}
    return {"rows": rows, "tests": tests, "barriers_kj_mol": barriers}


# ---------------------------------------------------------------------------------------------------------------
# 5. Reading the tests (prereg secs. 4 - 5)
# ---------------------------------------------------------------------------------------------------------------
def read_tests(t: Mapping[str, Any], b60_thiamine_flux: float) -> Dict[str, Dict[str, Any]]:
    h1 = {k: v["fragment_C2C3_share"] for k, v in t.items() if k.startswith("H1_")}
    h1_decl = [h1["H1_liu168_10mL_20min"], h1["H1_liu168_10mL_60min"]]
    h2 = t["H2_wang2012_potC"]
    thi_x = h2["thiamine_mft_flux_mM"] / b60_thiamine_flux if b60_thiamine_flux > 0 else float("nan")
    d3 = t["D3_hofmann145_panel"]
    d1 = t["D1_cerny_potD_fed_nf_share"]
    return {
        "H1": {"value": h1_decl, "sweep": h1, "band": "<= 0.12 (corrected reading)",
               "pass": all(v <= 0.12 for v in h1_decl),
               "printed_reading_pass": all(0.10 <= v <= 0.38 for v in h1_decl), "binding": True,
               "note": "guard: ~0 power against E6 / E11 (no member touches r_ha_mp_mft)"},
        "H2": {"value": h2["xylose_share_of_mft"], "band": "<= 0.04", "pass": h2["xylose_share_of_mft"] <= 0.04,
               "nf_share": h2["nf_share"], "ddp_share": h2["ddp_share"],
               "thiamine_flux_x_vs_b60": thi_x, "through_thiamine_only": bool(
                   h2["xylose_share_of_mft"] <= 0.04 and thi_x > 3.0),
               "sensitivity_dossier_gly9.8": t["H2_wang2012_potC_dossier_gly9.8"]["xylose_share_of_mft"],
               "pot_A_gsh_reported": t["H2_wang2012_potA_gsh"]["xylose_share_of_mft"], "binding": True},
        "D1": {"value": d1, "band": "[0.03, 0.15]", "pass": 0.03 <= d1 <= 0.15,
               "route_flux_share": t["D1_cerny_potD_nf_route_flux_share"]},
        "D2": {"value": t["D2_givaudan_ribose_over_nf_mft_2000mM"], "band": ">= 10",
               "pass": t["D2_givaudan_ribose_over_nf_mft_2000mM"] >= 10.0,
               "sensitivity_0.3M": t["D2_givaudan_ribose_over_nf_mft_300mM"]},
        "D3": {"value": d3["predicted"], "measured": d3["measured"], "band": "within 3x of the measured",
               "pass": d3["measured"] / 3.0 <= d3["predicted"] <= 3.0 * d3["measured"]},
        "D4": {"value": t["D4_hofmann145_ribose_nf_route_share"], "band": "[0.30, 1.0]",
               "pass": 0.30 <= t["D4_hofmann145_ribose_nf_route_share"] <= 1.0,
               "at_100C": t["ribose_cys_100C_nf_route_share"], "at_120C": t["ribose_cys_120C_nf_route_share"]},
        "D5": {"value": t["D5_cerny_potC_mp3p_over_mft"], "band": "<= 0.1",
               "pass": t["D5_cerny_potC_mp3p_over_mft"] <= 0.1},
    }


def ledger(rows: Mapping[str, Any], b60_rows: Mapping[str, Any],
           stage0_rows: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    """The compensation-ledger guard (prereg sec. 4.3): a row worse than its B60 raw |residual| by > 0.3 is allowed
    only if it is in the ledger AND moved in its shipped residual's direction."""
    worse, allowed, violations = {}, {}, {}
    for rid, row in rows.items():
        raw, base = row["raw"], b60_rows[rid]["raw"]
        delta = (abs(raw) - abs(base)) if np.isfinite(raw) else float("inf")
        entry = {"raw": raw, "b60_raw": base, "delta_abs": delta, "units": row["units"]}
        if stage0_rows is not None:
            entry["stage0_raw"] = stage0_rows[rid]["raw"]
            entry["delta_abs_vs_stage0"] = abs(raw) - abs(stage0_rows[rid]["raw"])
        if delta > GUARD_DEX:
            worse[rid] = entry
            same_direction = np.sign(raw) == np.sign(base) and base != 0.0
            if rid in LEDGER and same_direction:
                allowed[rid] = {**entry, "named_error": LEDGER[rid]}
            else:
                violations[rid] = entry
    return {"holds": not violations, "worse_than_b60_by_0.3": worse, "allowed_by_ledger": allowed,
            "violations": violations}


# ---------------------------------------------------------------------------------------------------------------
# 6. Stage 0
# ---------------------------------------------------------------------------------------------------------------
def _inputs() -> List[Path]:
    return [ROOT / "src" / "kinetic_core" / "sulfur.py", ROOT / "src" / "kinetic_core" / "parameters_sulfur.py",
            ROOT / "src" / "kinetic_core" / "config.py", BASE_REPORT, PREREG,
            data_paths.BENCHMARKS_DIR / f"{HOFMANN_BUNDLES[0]}.json", data_paths.BENCHMARKS_DIR / f"{HOFMANN_BUNDLES[1]}.json"]


def _check_base() -> None:
    assert PREREG.exists(), "E6 + E11 is pre-registered; the prereg must exist before any score"
    assert SHIPPED_REPORT.resolve() == BASE_REPORT.resolve(), (
        f"the shipped report is {SHIPPED_REPORT.name}, not B60's: the prereg's sec. 7 must be re-run on the new base")
    assert engine_config.Config.from_globals() == engine_config.Config.shipped(), "the process does not start shipped"


def run_stage0() -> Dict[str, Any]:
    _check_base()
    x = incumbent_vector()
    t0 = time.time()
    scored = {}
    for label, cfg in stage0_configs().items():
        s = time.time()
        scored[label] = score_vector(x, cfg, panel_from_report=False)
        print(f"  stage 0 {label}: cost {scored[label]['rows']['cost']:.4f} ({time.time() - s:.0f} s)", flush=True)
    b60_rows = scored["shipped"]["rows"]["rows"]
    b60_thi = scored["shipped"]["tests"]["H2_wang2012_potC"]["thiamine_mft_flux_mM"]
    configs = {}
    for label, s in scored.items():
        configs[label] = {
            "config_diff_from_shipped": {k: repr(v[1]) for k, v in stage0_configs()[label].diff().items()},
            "barriers_kj_mol": s["barriers_kj_mol"], "cost": s["rows"]["cost"],
            "tests": read_tests(s["tests"], b60_thi), "raw_tests": s["tests"],
            "ledger": ledger(s["rows"]["rows"], b60_rows), "rows": s["rows"]["rows"]}
    return {
        "artifact": "kinetic_core_mft_route_barriers_stage0",
        "provenance": provenance.provenance_block("kinetic_core_mft_route_barriers_stage0", generated_by=GENERATED_BY,
                                                  wave="E6 + E11 stage 0", inputs=_inputs()),
        "prereg": data_paths.rel(PREREG), "base_report": data_paths.rel(BASE_REPORT),
        "what": ("B60's frozen vector, NOT refitted, under each pre-registered barrier configuration (prereg sec. 2); "
                 "every test of sec. 4 and all 52 B60 rows; the ledger guard against the shipped row residuals."),
        "vector": "B60's frozen optimum (clipped into its bands)", "quick_setting": B60.QUICK_SETTING,
        "workers": int(os.environ.get("MAILLARD_WORKERS", "1")), "wall_seconds": round(time.time() - t0, 1),
        "configs": configs}


def _fmt(v: Any) -> str:
    if isinstance(v, bool):
        return "PASS" if v else "FAIL"
    if isinstance(v, float):
        return f"{v:.3g}"
    if isinstance(v, (list, tuple)):
        return " / ".join(_fmt(x) for x in v)
    return str(v)


def _tests_table(per: Mapping[str, Mapping[str, Any]]) -> List[str]:
    labels = list(per)
    out = ["| test | band | " + " | ".join(labels) + " |", "|---|---|" + "---|" * len(labels)]
    for t in ("H1", "H2", "D1", "D2", "D3", "D4", "D5"):
        band = per[labels[0]][t]["band"]
        cells = [f"{_fmt(per[l][t]['value'])} ({_fmt(per[l][t]['pass'])})" for l in labels]
        out.append(f"| {t} | {band} | " + " | ".join(cells) + " |")
    return out


def render_stage0(p: Mapping[str, Any]) -> str:
    c = p["configs"]
    lines = [
        "# E6 + E11 stage 0: the barrier change at the INCUMBENT constants (no refit)", "",
        f"*Generated by `{GENERATED_BY} --stage0` on {p['provenance']['generated_on']}. Pre-registration: "
        f"[{Path(p['prereg']).name}]({Path(p['prereg']).name}) sec. 2. Base: `{p['base_report']}`. Nothing ships.*", "",
        "B60's frozen vector under each configuration: `shipped` (the baseline); `N1_B62_alone` (k_tdp_fur 147, "
        "k_fur_decay 46.5 kJ/mol); `A64` / `A100` / `A133` (B62 + k_dpo_nf 83.8 + k_nf_mft at the profile point). "
        "Each configuration is one A3 `Config`, differing from the shipped one in `sulfur.furfural_barriers` only.", "",
        "## Tests (value, pass)", ""] + _tests_table({k: v["tests"] for k, v in c.items()}) + [
        "", "H1 is the 10 mL declared volume at 20 / 60 min; H2 is Wang 2012 pot C as pre-registered (Gly 12.1 mM).", "",
        "## Cost and the compensation-ledger guard", "",
        "| config | cost (52 rows, quick) | rows worse by > 0.3 | allowed by the ledger | violations | guard |",
        "|---|---|---|---|---|---|"]
    for k, v in c.items():
        L = v["ledger"]
        lines.append(f"| {k} | {v['cost']:.3f} | {len(L['worse_than_b60_by_0.3'])} | "
                     f"{', '.join(L['allowed_by_ledger']) or '-'} | "
                     + (", ".join(f"{r} ({e['b60_raw']:+.2f} -> {e['raw']:+.2f})" for r, e in L["violations"].items())
                        or "-") + f" | {'holds' if L['holds'] else 'FAILS'} |")
    lines += ["", "## Supplementary reads", "", "| config | H2 NF share | H2 thiamine flux x B60 | H2 dossier Gly 9.8 | "
              "H2 pot A (GSH) | D1 route flux share | D2 at 0.3 M | D4 at 100 / 120 C |", "|---|---|---|---|---|---|---|---|"]
    for k, v in c.items():
        t = v["tests"]
        lines.append(f"| {k} | {_fmt(t['H2']['nf_share'])} | {_fmt(t['H2']['thiamine_flux_x_vs_b60'])} | "
                     f"{_fmt(t['H2']['sensitivity_dossier_gly9.8'])} | {_fmt(t['H2']['pot_A_gsh_reported'])} | "
                     f"{_fmt(t['D1']['route_flux_share'])} | {_fmt(t['D2']['sensitivity_0.3M'])} | "
                     f"{_fmt(t['D4']['at_100C'])} / {_fmt(t['D4']['at_120C'])} |")
    lines += ["", "Barriers in force (kJ/mol):", ""]
    for k, v in c.items():
        lines.append(f"- {k}: " + ", ".join(f"{b} {e:.1f}" for b, e in v["barriers_kj_mol"].items() if e is not None))
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------------------------
# 7. The refit: one member, the queue, the verdict
# ---------------------------------------------------------------------------------------------------------------
def member_path(name: str, warm: bool = False) -> Path:
    return MEMBER_DIR / f"{name}{'_warm' if warm else ''}.json"


def read_member(name: str) -> Optional[Dict[str, Any]]:
    """The member as it is read: its warm relaunch if it had one, else its first run."""
    for warm in (True, False):
        p = member_path(name, warm)
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))
    return None


def fit_one(name: str, warm: bool, budget: int = BUDGET) -> Dict[str, Any]:
    _check_base()
    m = member(name)
    cfg = member_config(name)
    prob = problem(name)
    if warm:
        first = json.loads(member_path(name).read_text(encoding="utf-8"))
        assert not first["converged"] and first["budget_exhausted"], f"{name}: only a budget-stopped member is relaunched"
        start = np.array([first["theta"][k] for k in prob.keys], dtype=float)
    else:
        start = start_theta(name)
    with engine_config.applied(cfg):
        print(f"[{name}] {'warm ' if warm else ''}start, {len(prob.keys)} free, budget {budget}, ftol {FTOL}, "
              f"workers {os.environ.get('MAILLARD_WORKERS', '1')}", flush=True)
        f_ = fitting.fit(prob, start, budget=budget, ftol=FTOL)
    x_full = to_full(f_.theta, name)
    out = {
        "artifact": f"kinetic_core_mft_route_barriers_member_{name}",
        "provenance": provenance.provenance_block(f"kinetic_core_mft_route_barriers_member_{name}",
                                                  generated_by=GENERATED_BY, wave="E6 + E11 refit", inputs=_inputs()),
        "member": name, "warm_relaunch": warm, "barriers_kj_mol": dict(m.barriers), "start": m.start,
        "config_diff_from_shipped": {k: repr(v[1]) for k, v in cfg.diff().items()},
        "coordinates": [{"key": k, "kind": kd, "lower": lo, "upper": hi} for k, kd, lo, hi in coordinates(name)],
        **f_.as_dict(prob.keys), "ftol": FTOL,
        "converged": f_.status in CONVERGED_STATUS, "x_full": [float(v) for v in x_full],
        "parameters": {k: float(x_full[i]) for i, k in enumerate(B60.ALL_KEYS)},
        "residual_by_row": {rid: float(r) for rid, r in zip(prob.rows, f_.residuals)},
        "workers": int(os.environ.get("MAILLARD_WORKERS", "1")),
    }
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    artifact_io.write_artifact(out, member_path(name, warm))
    print(f"[{name}] cost {f_.start_cost:.3f} -> {f_.cost:.3f}, status {f_.status}, evals {f_.evaluations}, "
          f"{f_.seconds:.0f} s, converged {out['converged']}", flush=True)
    return out


def best_a(names: Sequence[str] = ("A64", "A100", "A133")) -> Optional[str]:
    done = [(n, read_member(n)) for n in names]
    conv = [(n, m) for n, m in done if m is not None and m["converged"]]
    return min(conv, key=lambda nm: nm[1]["cost"])[0] if conv else None


def run_queue(workers: int, slots: int = 2) -> int:
    """Every member, ``slots`` at a time, each with ``workers`` workers: N1 + A64, A100 + A133 (a budget-stopped
    member relaunched warm once), then start 1 and the 93 kJ/mol sensitivity of the best converged A; then the
    verdict. Each member is its own process; its log goes to scratch/mft_route_barriers/ in the worktree."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "MAILLARD_WORKERS": str(workers), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    pending: List[Tuple[str, bool]] = [(n, False) for n in ("N1", "A64", "A100", "A133") if read_member(n) is None]
    running: Dict[Tuple[str, bool], subprocess.Popen] = {}
    #: members that need nothing more: converged, relaunched once, or exited without a record (reported, not retried)
    finished = {n for n in ("N1", "A64", "A100", "A133") if read_member(n) is not None
                and (read_member(n)["converged"] or read_member(n)["warm_relaunch"] or not read_member(n)["budget_exhausted"])}
    pending = [(n, True) for n in ("N1", "A64", "A100", "A133") if n not in finished and read_member(n) is not None] + pending
    derived_queued = False

    def launch(job):
        name, warm = job
        log = open(LOG_DIR / f"{name}{'_warm' if warm else ''}.log", "w")
        cmd = [sys.executable, __file__, "--member", name] + (["--warm"] if warm else [])
        running[job] = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
        print(f"{time.strftime('%H:%M:%S')} launched {name}{' (warm)' if warm else ''}", flush=True)

    while pending or running:
        while pending and len(running) < slots:
            launch(pending.pop(0))
        time.sleep(30)
        for job, proc in list(running.items()):
            if proc.poll() is None:
                continue
            del running[job]
            name, warm = job
            m = read_member(name)
            print(f"{time.strftime('%H:%M:%S')} {name}{' (warm)' if warm else ''} exited {proc.returncode}: "
                  + ("no record" if m is None else f"cost {m['cost']:.3f} status {m['status']} converged {m['converged']}"),
                  flush=True)
            if m is not None and not warm and m["budget_exhausted"]:
                pending.insert(0, (name, True))
            else:
                finished.add(name)
        if all(n in finished for n in ("A64", "A100", "A133")) and not derived_queued:
            derived_queued = True
            best = best_a()
            print(f"{time.strftime('%H:%M:%S')} best converged A: {best}", flush=True)
            if best is not None:
                pending += [(d, False) for d in (f"{best}-s1", f"{best}-93") if read_member(d) is None]
    print(f"{time.strftime('%H:%M:%S')} members done; scoring the verdict", flush=True)
    return run_verdict_and_write()


def verdict(members: Mapping[str, Mapping[str, Any]], best: Optional[str]) -> Dict[str, Any]:
    """Prereg sec. 5, read literally. Clauses are checked from the most conservative; every one that matches is
    listed; the primary verdict is the first that matches."""
    a_names = [n for n in ("A64", "A100", "A133") if n in members]
    conv_a = [n for n in a_names if members[n]["converged"]]
    clauses: List[str] = []
    notes: List[str] = []
    if not conv_a:
        clauses.append("DO NOT SHIP: no A member converges")
    elif all(not members[n]["ledger"]["holds"] for n in conv_a):
        clauses.append("DO NOT SHIP: the ledger guard fails at every (converged) A member")
    if conv_a and all(not members[n]["tests"]["D1"]["pass"] for n in conv_a):
        clauses.append("REFUSE (declared band): D1 fails at every k_nf_mft point after the refit")
    s1_ok = None
    if best is not None:
        s1 = members.get(f"{best}-s1")
        if s1 is None or not s1["converged"]:
            s1_ok = False
            notes.append(f"start 1 of {best} {'missing' if s1 is None else 'not converged'}")
        else:
            s1_ok = abs(s1["cost"] - members[best]["cost"]) <= START1_COST_AGREEMENT
            notes.append(f"start 1 cost {s1['cost']:.3f} vs {best} {members[best]['cost']:.3f}")
        if not s1_ok:
            clauses.append("INCONCLUSIVE: start 1 does not agree within 2 cost units (sec. 3)")
        t = members[best]["tests"]
        guard = members[best]["ledger"]["holds"]
        dev_ok = t["D1"]["pass"] and t["D2"]["pass"]
        if dev_ok and guard and (not t["H2"]["pass"] or t["H2"]["through_thiamine_only"]):
            clauses.append("INCONCLUSIVE: D1 and D2 pass and the guard holds, but H2 fails or passes only through "
                           "the thiamine route")
        d4_blocks = (not t["D4"]["pass"]) and t["D1"]["pass"] and t["D1"]["value"] < 0.03
        if (t["H2"]["pass"] and not t["H2"]["through_thiamine_only"] and t["H1"]["pass"] and guard and dev_ok
                and t["D3"]["pass"] and not d4_blocks and s1_ok):
            clauses.append("SHIP-eligible (owner approval required; this wave does not ship)")
    order = ("DO NOT SHIP", "REFUSE", "INCONCLUSIVE", "SHIP-eligible")
    primary = next((c for o in order for c in clauses if c.startswith(o)), None)
    if primary is None:
        primary = "NO sec. 5 CLAUSE MATCHES (not SHIP-eligible)"
        notes.append("the rule is silent on this case; the failing conditions at the best A are listed")
    return {"primary": primary, "clauses_matched": clauses, "best_converged_A": best, "start1_agrees": s1_ok,
            "notes": notes}


def run_verdict() -> Dict[str, Any]:
    _check_base()
    t0 = time.time()
    stage0 = json.loads(STAGE0_OUT.read_text(encoding="utf-8"))
    b60_rows = stage0["configs"]["shipped"]["rows"]
    b60_thi = stage0["configs"]["shipped"]["raw_tests"]["H2_wang2012_potC"]["thiamine_mft_flux_mM"]
    stage0_of = {"N1": "N1_B62_alone", "A64": "A64", "A100": "A100", "A133": "A133"}
    names = [p.stem.replace("_warm", "") for p in sorted(MEMBER_DIR.glob("*.json"))]
    members: Dict[str, Dict[str, Any]] = {}
    for name in dict.fromkeys(names):
        rec = read_member(name)
        x = np.array(rec["x_full"], dtype=float)
        s = score_vector(x, member_config(name), panel_from_report=True)
        base0 = stage0_of.get(re.sub(r"-s1$", "", name))
        s0_rows = stage0["configs"][base0]["rows"] if base0 else None
        members[name] = {
            "converged": rec["converged"], "status": rec["status"], "cost": rec["cost"],
            "rescored_cost": s["rows"]["cost"], "start_cost": rec["start_cost"], "evaluations": rec["evaluations"],
            "budget_exhausted": rec["budget_exhausted"], "warm_relaunch": rec["warm_relaunch"],
            "seconds": rec["seconds"], "barriers_kj_mol": rec["barriers_kj_mol"], "theta": rec["theta"],
            "route_constants": {k: rec["parameters"][k] for k in ("k_dpo_nf", "k_dpo_ddp", "k_nf_mft", "k_ddp_mft",
                                                                    "k_ddp_mft_hs", "k_tdp_fur", "k_fur_decay")},
            "ddp_offset_vs_b60": rec["theta"].get(DDP_OFFSET),
            "tests": read_tests(s["tests"], b60_thi), "raw_tests": s["tests"],
            "ledger": ledger(s["rows"]["rows"], b60_rows, s0_rows), "stage0_reference": base0,
            "rows": s["rows"]["rows"]}
        print(f"  scored {name}: cost {s['rows']['cost']:.3f}", flush=True)
    best = best_a()
    lap = None
    if best is not None:
        rec = read_member(best)
        prob = problem(best)
        with engine_config.applied(member_config(best)):
            lap = fitting.laplace(prob, np.array([rec["theta"][k] for k in prob.keys], dtype=float))
        lap = {k: v for k, v in lap.items() if k != "correlation"} | {
            "correlation_new_coordinates": {
                k: dict(zip(prob.keys, map(float, np.asarray(lap["correlation"])[prob.keys.index(k)])))
                for k in (NF_KEY, DDP_OFFSET, "k_nf_mft")}}
    v = verdict(members, best)
    return {
        "artifact": "kinetic_core_mft_route_barriers_refit",
        "provenance": provenance.provenance_block(
            "kinetic_core_mft_route_barriers_refit", generated_by=GENERATED_BY, wave="E6 + E11 refit",
            inputs=_inputs() + [STAGE0_OUT] + sorted(MEMBER_DIR.glob("*.json"))),
        "prereg": data_paths.rel(PREREG), "base_report": data_paths.rel(BASE_REPORT),
        "optimiser": {"harness": "src/kinetic_core/fitting.fit", "budget": BUDGET, "ftol": FTOL,
                      "converged_status": list(CONVERGED_STATUS),
                      "finite_difference": "forward, absolute step per coordinate kind (fitting.STEPS)",
                      "objective": "B60's 52 rows, declared sigmas fixed, quick setting"},
        "members": members, "laplace_best_A": lap, "verdict": v, "ships": False,
        "wall_seconds": round(time.time() - t0, 1)}


def run_verdict_and_write() -> int:
    payload = run_verdict()
    artifact_io.write_artifact(payload, REFIT_OUT, render=render_refit)
    print(f"wrote {REFIT_OUT}: {payload['verdict']['primary']}")
    return 0


def render_refit(p: Mapping[str, Any]) -> str:
    ms = p["members"]
    v = p["verdict"]
    lines = [
        "# E6 + E11 refit: per-route MFT formation barriers on B60's objective", "",
        f"*Generated by `{GENERATED_BY} --verdict` on {p['provenance']['generated_on']}. Pre-registration: "
        f"[{Path(p['prereg']).name}]({Path(p['prereg']).name}). Stage 0: "
        "[kinetic_core_mft_route_barriers_stage0.md](kinetic_core_mft_route_barriers_stage0.md). Nothing ships.*", "",
        f"**Verdict (sec. 5, read literally): {v['primary']}.**", ""]
    lines += [f"- clause matched: {c}" for c in v["clauses_matched"]] + [f"- note: {n}" for n in v["notes"]]
    lines += ["", "## Members", "", "| member | status | converged | evals | cost (start -> end) | k_dpo_nf | "
              "DDP offset | guard |", "|---|---|---|---|---|---|---|---|"]
    for n, m in ms.items():
        lines.append(f"| {n}{' (warm)' if m['warm_relaunch'] else ''} | {m['status']} | {m['converged']} | "
                     f"{m['evaluations']} | {m['start_cost']:.2f} -> {m['cost']:.2f} | "
                     f"{m['route_constants']['k_dpo_nf']:+.2f} | "
                     f"{_fmt(m['ddp_offset_vs_b60']) if m['ddp_offset_vs_b60'] is not None else '-'} | "
                     f"{'holds' if m['ledger']['holds'] else 'FAILS'} |")
    lines += ["", "## Tests (value, pass)", ""] + _tests_table({n: m["tests"] for n, m in ms.items()})
    lines += ["", "## Compensation ledger", ""]
    for n, m in ms.items():
        L = m["ledger"]
        lines.append(f"- **{n}**: allowed {', '.join(L['allowed_by_ledger']) or 'none'}; violations "
                     + (", ".join(f"{r} (B60 {e['b60_raw']:+.2f} -> {e['raw']:+.2f} {e['units']})"
                                  for r, e in L["violations"].items()) or "none"))
    lap = p.get("laplace_best_A")
    if lap:
        lines += ["", f"## Laplace at the best A ({v['best_converged_A']})", "",
                  f"chi2_red {lap['chi2_red']:.3f}. Route coordinates: "
                  + ", ".join(f"{k} sigma {lap['sigma'][k]:.3g} (identified {lap['identified'][k]}, on bound "
                              f"{lap['at_bound'][k]})" for k in (NF_KEY, DDP_OFFSET) if k in lap["sigma"])]
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------------------------
# 8. V-checks (prereg sec. 8)
# ---------------------------------------------------------------------------------------------------------------
def run_vcheck(probes: Optional[Path] = None) -> int:
    """The cost identity always; the route scorer against the prereg's probe outputs (sec. 10) when ``probes`` (the
    directory holding probe_b.json and probe_b_heldout.json) is given."""
    _check_base()
    x = incumbent_vector()
    shipped = engine_config.Config.shipped()
    with engine_config.applied(shipped):
        c_ship = 0.5 * float(np.sum(B60.residual_vector(x, True) ** 2))
        mine = rows_at(x)["cost"]
    lump = lumped_ea()
    ident = barrier_config(((NF_KEY, lump), ("k_nf_mft", lump)))
    with engine_config.applied(ident):
        theta = theta_of(x, "A64")
        assert np.array_equal(to_full(theta, "A64"), x), "the A coordinates do not round-trip at the incumbent"
        c_ident = 0.5 * float(np.sum(problem("A64").evaluate(theta) ** 2))
    print(f"shipped cost {c_ship:.9f} (rows_at {mine:.9f}); extra barriers at the lump, offsets at shipped: "
          f"{c_ident:.9f}; |diff| {abs(c_ident - c_ship):.2e}")
    assert abs(mine - c_ship) < 1e-9 and abs(c_ident - c_ship) < 1e-9
    if probes is not None:
        with engine_config.applied(shipped):
            t = route_tests(B23.build_parameters(x))
            d3 = hofmann_panel_ratio()
        pb = json.loads((probes / "probe_b.json").read_text())
        ph = json.loads((probes / "probe_b_heldout.json").read_text())
        pairs = [(t["D1_cerny_potD_fed_nf_share"], pb["cerny_potD_fed_nf_share"]),
                 (t["D5_cerny_potC_mp3p_over_mft"], pb["cerny_potC_mp3p_over_mft"]),
                 (t["D2_givaudan_ribose_over_nf_mft_2000mM"], pb["givaudan_ribose_over_nf_mft_2000mM"]),
                 (t["D2_givaudan_ribose_over_nf_mft_300mM"], pb["givaudan_ribose_over_nf_mft_300mM"]),
                 (t["D4_hofmann145_ribose_nf_route_share"], pb["hofmann145_ribose_nf_route_share"]),
                 (t["H2_wang2012_potC"]["xylose_share_of_mft"], ph["wang2012_potC_xylose_share_of_mft"]),
                 (t["H1_liu168_10mL_20min"]["nf_share"], ph["liu168_10mL_20min"]["nf_share"]),
                 (d3["predicted"], 11.471713632077048)]
        worst = max(abs(a - b) / max(abs(b), 1e-300) for a, b in pairs)
        print(f"route scorer vs the prereg probes: worst relative difference {worst:.2e}")
        assert worst < 1e-6
    print("V-checks pass")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcheck", action="store_true")
    parser.add_argument("--probes", type=Path, help="--vcheck: the prereg probes' output directory (sec. 10)")
    parser.add_argument("--stage0", action="store_true")
    parser.add_argument("--member")
    parser.add_argument("--warm", action="store_true")
    parser.add_argument("--budget", type=int, default=BUDGET)
    parser.add_argument("--queue", action="store_true")
    parser.add_argument("--queue-workers", dest="queue_workers", type=int, default=4)
    parser.add_argument("--verdict", action="store_true")
    parser.add_argument("--render", action="store_true", help="re-render the .md twins from the written JSON")
    args = parser.parse_args(argv)
    if args.render:
        for path, render in ((STAGE0_OUT, render_stage0), (REFIT_OUT, render_refit)):
            if path.exists():
                path.with_suffix(".md").write_text(render(json.loads(path.read_text(encoding="utf-8"))), encoding="utf-8")
        return 0
    if args.vcheck:
        return run_vcheck(args.probes)
    if args.stage0:
        payload = run_stage0()
        artifact_io.write_artifact(payload, STAGE0_OUT, render=render_stage0)
        print(f"wrote {STAGE0_OUT}")
        return 0
    if args.member:
        fit_one(args.member, args.warm, args.budget)
        return 0
    if args.queue:
        return run_queue(args.queue_workers)
    if args.verdict:
        return run_verdict_and_write()
    parser.error("one of --vcheck, --stage0, --member, --queue, --verdict")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
