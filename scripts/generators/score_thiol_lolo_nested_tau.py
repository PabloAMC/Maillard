#!/usr/bin/env python
"""
score_thiol_lolo_nested_tau.py -- H0: the thiol lane's (R3En) between-laboratory tau for a laboratory it has never
seen, estimated out of sample from EX-056's leave-one-laboratory-out records and H0's ten leave-two-laboratories-out
pair records; the shipped and the nested predictive arms scored against the two chemistry-free nulls.

Pre-registration (v3, frozen before any H0 number): results/validation/kinetic_core_between_lab_tau_h0_prereg.md;
its review: kinetic_core_between_lab_tau_h0_prereg_review.md. "sec." below always means the prereg's sections.

    python scripts/generators/score_thiol_lolo_nested_tau.py --guard      # sec. 4.1: pass / fail only
    python scripts/generators/score_thiol_lolo_nested_tau.py --design     # sec. 7.1: the design table, no residual
    python scripts/generators/score_thiol_lolo_nested_tau.py --controls   # secs. 8, 3.5 and the D2 grid, synthetic
    python scripts/generators/score_thiol_lolo_nested_tau.py --score      # secs. 2 - 6: the real scoring, run ONCE

--controls writes results/validation/kinetic_core_thiol_lolo_nested_tau_controls.{json,md}; --score writes
results/validation/kinetic_core_thiol_lolo_nested_tau.{json,md}.

Blinding (sec. 7; review N2). --guard rebuilds EX-056's published per-row values and prints pass / fail and the
largest relative difference only. --design reads row keys, classes, groups and sigmas. --controls draws synthetic
residuals on the real design (row sigma plus EX-056's saved class noise, a published aggregate). Only --score
computes an H0 quantity from real residuals; it refuses to run unless all ten pair records exist, each complete and
with a non-empty freeze commit, the controls record exists, and no hard control failed. Before any nested number it
checks that training_terms reproduces EX-056's saved terms and recomputes engine guard 4.2(c) from the guard record
(amendment A1, `kinetic_core_thiol_lolo_nested_tau_amendment_A1.md`). It runs once: with the primary record present it
refuses unless --amendment REASON is given, and then writes `..._amendment_<k>.{json,md}`. fit_thiol_lolo.py, the
frozen EX-056 generator, is imported and never modified.
"""
from __future__ import annotations

import argparse
import dataclasses
import functools
import itertools
import json
import math
import subprocess
import sys
import tempfile
import time
import zlib
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows, observations  # noqa: E402

import fit_deoxyosone_submodel as sm3  # noqa: E402
import fit_joint_deoxyosone_thiol as fj  # noqa: E402
import fit_thiol_lolo as lolo  # noqa: E402
import fit_thiol_lolo_pairs as pair_runner  # noqa: E402
from fit_thiol_lolo import (LEVELS, NULLS, engine_state, lab_of, level_of, load_fit_table, make_folds,  # noqa: E402
                            null_residuals, predictive_model, training_terms)

# ----------------------------------------------------------------------------------------------------------------
# Constants (every number here is the prereg's, unless marked "choice")
# ----------------------------------------------------------------------------------------------------------------

STEM = "kinetic_core_thiol_lolo_nested_tau"
CONTROLS_STEM = f"{STEM}_controls"
GENERATED_BY = "scripts/generators/score_thiol_lolo_nested_tau.py"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_between_lab_tau_h0_prereg.md"
REVIEW = data_paths.VALIDATION_DIR / "kinetic_core_between_lab_tau_h0_prereg_review.md"
#: where the EX-056 fold records and the pair records are read (module-level so tests and the T1b canary can point
#: them at a temporary directory)
FOLD_DIR = data_paths.VALIDATION_DIR
PAIR_DIR = data_paths.VALIDATION_DIR
OUT_DIR = data_paths.VALIDATION_DIR

SCORABLE = tuple(pair_runner.SCORABLE)          # sec. 1: the five scorable laboratories
NEVER_OUTER = tuple(pair_runner.NEVER_HELD)     # Nestle: in every training set, never outer, only in D-all
SETS = ("warm", "reset_where_leak_prone")       # EX-056's two sets (sec. 3.3)
MODELS = ("engine",) + tuple(NULLS)
N_CONSTANTS_ENGINE = 16                         # R3En's 16 coordinates (sec. 3.1)
ABS, FIG, FREE = "ABS", "FIG", "FREE"
FIG_METHOD = "semi_figure"
FIG_FLOOR = 0.30                                # tau_FIG,used = max(0.30, tau_FIG) (sec. 2)
RATIO_FLOOR = 0.10                              # zero-safe ratios (sec. 3.7)
D2_LINE = 0.8                                   # D2's 0.8 dex line
ONE_DEEP_LEAK = 0.3                             # |tau_onedeep - tau_strict| printed as the measured leak (sec. 3.4)
D2_GRID = tuple(round(0.10 + 0.05 * k, 2) for k in range(31))   # 0.10, 0.15, ..., 1.60 dex
D2_REPLICATES = 2000
D2_ALPHA = 0.05
R_CONTROLS = 1000
R_RERUN = 4000
SEED = 20261101                                 # choice: the controls' fixed seed (the prereg's date)
GUARD_TOL = 1e-9
T1_TOL = 1e-12
GATE_LEVEL = 0.90
D1_FLOOR = 0.70
ACCEPTANCE_BAND = tuple(lolo.GATE_BAND)
SIGN_FLIP_THRESHOLD = 2.0 / 32.0
SIGN_FLIP_THRESHOLD_FALLBACK = 1.0 / 32.0       # if D3a fails (sec. 8)
UNSTABLE_GATE = 0.3                             # P(UNSTABLE | homogeneous truth, tau = 0.8) <= 0.3 (sec. 6)
CONCENTRATION_TOP = 0.05
CONCENTRATION_SHARE = 0.5
OMEGA_NOTICE = 0.5
#: sec. 3.5: the leak toy's row counts (each laboratory's held-out rows in EX-056); asserted against the records
TOY_COUNTS = {"jiangnan_zhang_x": 454, "dalian_polytechnic_chen_d": 19, "dfa_hofmann_schieberle": 16,
              "nestle_blank": 19, "ochanomizu_murata": 5, "reading_mottram": 8}
ST_NODES = 64                                   # choice: quadrature nodes for the St row marginal (controls)
ST_NODES_REAL = 400                             # choice: the same quadrature, finer, for the one real run
LOG2PI = math.log(2.0 * math.pi)

#: positive / negative / toy / D3 configurations (sec. 8)
P_TAUS, P_OMEGAS, P_FS = (0.3, 0.8, 1.2), (0.0, 0.3), (0.0, 0.5, 1.0)
P1_HARD_TAUS = (0.8, 1.2)
TOY_KAPPAS, TOY_TAUS, TOY_DS = (0.0, 0.5, 1.0), (0.3, 0.8, 1.2), (-3.0, -1.0, 1.0, 3.0)

CHOICES = [
    "Class noise for the controls' design: each row's s_i^2 = sigma_row^2 + the class noise saved in EX-056's fold "
    "record for the row's own laboratory (a published aggregate, not a residual), through fit_thiol_lolo."
    "predictive_model. The warm set uses `models.engine_warm.terms`; the reset set uses `models.engine_reset.terms` "
    "where that fold has a reset arm (DFA, Jiangnan, Nestle) and the warm terms elsewhere, so the two sets differ in "
    "the design exactly as the real records do. The same s_i serve a laboratory's rows whether they are outer "
    "(fold L) or inner (a pair fit): the pair fits' training_terms do not exist yet.",
    "Inner masks in the controls: an inner laboratory L''s rows are EX-056's paired rows of L' AND the rows both "
    "strict nulls (trained without L and L') can answer, the latter computed from row keys and species presence "
    "only (`null_answers`, no residual). In --score the mask is: finite under EVERY engine variant of the pair "
    "record AND both strict nulls (mirroring EX-056's paired definition, so one mask serves both sets and all three "
    "models).",
    "Synthetic truth: one residual per row and replicate, the same in every record where the row appears (kappa = 0 "
    "outside the leak toy); the ABS and FIG level offsets are drawn independently per class with the same tau; FREE "
    "groups get an N(0, 1 dex) level. Every class shares the configuration's omega.",
    "A laboratory enters omega(S) only if it holds at least one group with >= 2 rows; with none in S, omega(S) = 0 "
    "(flagged).",
    "P1's six estimates are tau_ABS(all) and each of the five outer laboratories' nested tau_ABS(S). For Dalian and "
    "Reading (no ABS row of their own) tau_ABS(S) rests on 3 inner laboratories, for DFA, Jiangnan and Ochanomizu on "
    "2. An estimate that takes the minimum-data fallback has no untruncated value and is listed, not checked.",
    "The negative controls N1 - N3 are run at omega in {0, 0.3}, both sets, and must pass in every one of those "
    "configurations (the prereg names tau = 0 only).",
    "P2's test inversion uses, per replicate, the f = 1 grid at the configuration's TRUE omega in place of the "
    "replicate's omega_hat(all) (recomputing a 31 x 2000 grid per replicate is out of budget); P2 is run at omega in "
    "{0, 0.3}. In --score the grid is drawn at the real omega_hat(all) of each set.",
    "D2's bounds: U = the largest grid tau0 with P(tau_hat <= tau_obs | tau0) >= 0.05; L = the smallest tau0 with "
    "P(tau_hat >= tau_obs | tau0) >= 0.05 (the mirror); tau_hat is the TRUNCATED tau_ABS(all). U is computed per set "
    "at that set's tau_obs and omega_hat(all), and the ship rule reads the larger U, as it reads the larger tau. "
    "Grid edges (amendment A1): a U at the grid top is reported 'U >= 1.60' (and used as 1.60); a U that does not "
    "exist because tau_obs is below the 5 % quantile even at tau0 = 0.10 is reported 'U <= 0.10' and used as 0.10; an "
    "L at the grid bottom is reported 'L <= 0.10'; an L that does not exist (no grid tau0 qualifies) 'L: none'. Every "
    "two-set maximum is NaN-safe and order-free; an undefined tau_ABS(all) refuses D2 rather than reading as small.",
    "UNSTABLE (sec. 3.7) compares the three drop-one estimates only (max t~ / min t~ > 2), and 'the other side of "
    "0.8' is relative to the full tau_ABS(all). For D2 the flag counts if it fires in either set; it acts only if "
    "the controls' P(UNSTABLE) at tau = 0.8 is <= 0.3 in EVERY configuration (f in {0, 0.5, 1}, omega in {0, 0.3}, "
    "both sets): the largest is the gating number.",
    "St (sec. 5): df = the number of inner laboratories holding the class (2 for most nests, 3 for tau(all)); arm (a)'s "
    "St uses df = the number of training laboratories holding the class in that fold. A class on the minimum-data "
    "fallback keeps a Gaussian level. The row marginal's CDF is computed by quadrature (Gauss-Hermite over the row "
    "error when the level scale >= the row SD, Gauss-Legendre over the t quantile otherwise; 64 nodes in the "
    "controls, 400 in --score); FREE rows keep the Gaussian rule.",
    "S2 (the parameter term) is optional and is NOT implemented.",
    "The relative branch (sec. 5.6) centres each group with >= 2 rows on its 1/s^2-weighted mean with the variance "
    "s_i^2 - 1/W + c_i omega^2(S) (the model's own nested omega); its LPD is the restricted Gaussian density with "
    "covariance diag(s^2 + omega^2) (invariant to the centring weights). A laboratory with no such group is left out "
    "of that comparison's mean.",
    "Concentration check: the per-row contribution is (marginal log density of (b) engine - of (b) null) / n_L / K; "
    "the flag reads |sum of the top 5 % rows by |contribution|| / |sum of all| > 0.5 (the absolute-value share is "
    "reported beside it).",
    "D3a and D3b simulate the engine and the null as two independent synthetic worlds on the same design (the "
    "engine's s_i), each with its own strict nested calibration; D3a passes if P(LPD criterion passes) <= 0.10 in both "
    "sets. 'D3 power = x' is P(the LPD criterion passes) at engine tau 0.8 vs null tau 1.2, the smaller of the two "
    "sets.",
    "Exact bootstrap: the 90 % interval is (the smallest multiset mean whose cumulative multinomial weight reaches "
    "0.05, the smallest whose cumulative weight reaches 0.95). The sign-flip p counts patterns whose mean is >= the "
    "observed mean minus 1e-12 (ties count).",
    "D-all (sec. 5.8): every finite engine row of each fold. (a) is EX-056's predictive model on those rows; (b) uses "
    "the laboratory's strict nested omega / tau(S) for the five scorable laboratories, and tau(all) / omega(all) "
    "(S = the five scorable laboratories) for Nestle. Nestle is in EVERY EX-056 training set but its own fold, so "
    "each of those five fits was trained with Nestle: Nestle's (b) calibration carries the one-deep kind of leak "
    "(sec. 3.4), acceptable only because D-all is a diagnostic (corrected in amendment A1).",
    "The jackknife (FRAGILE) drops one outer laboratory from the SCORING (D1 and D3 recomputed on 4 laboratories, "
    "16 sign patterns, the same numeric threshold); the dropped laboratory stays an inner laboratory of the others' "
    "calibrations. Per dropped laboratory the two sets' jackknife verdicts are combined by the same precedence as "
    "the reported verdict; FRAGILE means that combined verdict differs from the reported combined verdict for some "
    "dropped laboratory (amendment A1). The per-set comparison is reported beside it.",
    "Guard 4.1: relative difference |a - b| / max(|a|, |b|) (0 when equal); integer coverage counts must be equal. "
    "The curated-input check compares the provenance inputs under data/ with their EX-056 sha256. 'Read them at "
    "7b1bc7cd' is implemented for data/keys/laboratories.yml (changed since EX-056 by appended laboratory groups "
    "only): its 7b1bc7cd version must hash as EX-056 recorded, and the fit table loaded with it must equal today's "
    "row for row; then today's file is used. Any other changed curated input fails the guard.",
    "The ship-rule quantities tau_ABS(all), omega(all), tau_FIG(all) are computed on EX-056's paired rows of the "
    "five scorable laboratories.",
    "The leak toy's FREE groups are perturbed in T1 too (D added to their level), with every class perturbed "
    "separately; T1b adds +3 dex to every residual of lab L in every synthetic record (fold rows, pair held rows, "
    "pair training rows) and divides L's synthetic data values by 10^3 (the nulls), then reruns the scoring step's "
    "calibration (training_terms, strict nulls, the estimator) for all three models and both sets.",
    "--score refuses unless every pair record carries the same freeze commit, that commit exists and is an ancestor "
    "of HEAD, and this script, the pair runner, the controls record, the prereg, fit_thiol_lolo.py, "
    "fit_joint_deoxyosone_thiol.py, fit_deoxyosone_submodel.py, build_submodel_probe.py and src/kinetic_core/"
    "fitting.py are unchanged since it; a later change runs only with --amendment REASON, which is recorded (prereg "
    "sec. 7.5). Run once (amendment A1): --score refuses if the primary record exists, unless --amendment REASON is "
    "given; then it writes `..._amendment_<k>.{json,md}` (the next free k). Without a primary record, --amendment "
    "writes the primary record and stores the reason (a pre-scoring amendment).",
    "Reproduction guards inside --score (amendment A1), before any nested number; a failure is INVALID (refusal): "
    "(i) training_terms(null_residuals(kind, train, train), train, n) reproduces each of the 12 saved EX-056 null "
    "terms (6 folds x 2 nulls) to a relative 1e-9; (ii) training_terms on the guard record's train_engine_raw['warm'] "
    "(row keys equal to split(data, ('reading_mottram',))'s training rows) reproduces Reading's saved "
    "models.engine_warm.terms to a relative 1e-9; (iii) engine guard 4.2(c) recomputed from the guard record against "
    "EX-056's Reading record: the warm cost to a relative 1e-3, the held-out raw residuals to 0.01 dex with the same "
    "finite mask. All three, and the freeze check, are recorded in the output.",
    "One-deep mask (amendment A1): the one-deep sensitivity's inner blocks use the STRICT inner mask of sec. 3.2 (the "
    "rows finite under every variant of the pair (L, L') AND answered by both strict nulls), intersected with fold "
    "L''s paired rows where its residuals live, so |tau_onedeep - tau_strict| measures the leak and not a mask "
    "difference. The controls' T2 one-deep path is the same code (Nesting.calibrate), so T2 uses it too.",
    "Missing inner laboratories (amendment A1): only a NON-CONVERGED pair fit (or impossible strict nesting) caps a "
    "set's verdict at INCONCLUSIVE. An inner laboratory whose rows are all removed by the strict masks leaves S and the "
    "minimum-data rule of sec. 2 applies, flagged; it does not cap the verdict.",
    "Widths (amendment A1): the coverage table's 'median of lab medians' is the per-laboratory median 90 % width, "
    "then the median over the five laboratories (labelled so); beside it the pooled row median (every scored row of "
    "the five laboratories) and the lab-weighted mean width (the mean over laboratories of each laboratory's mean "
    "width).",
    "Owner note on the relative branch (amendment A1): when D2 takes branch (i), the note fires unless the relative "
    "branch is shown better than BOTH nulls in EVERY set (the sets combined by precedence, so the worse set governs).",
    "Labels and flags (amendment A1): verdict 4 reads 'CALIBRATED, NOT SHOWN BETTER THAN A NULL (D3 power = x, engine "
    "and null errors simulated independent)'. The controls' N2 / N3 / N5 outcomes are carried into the report as "
    "flags (N2: the nested tau_ABS(S) flagged; N3: the real tau_ABS flagged 'may be inflated by outliers'; N5: D2's "
    "multi-row-group fallback). The pair runner's unit tests passed at the freeze commit; the scorer cannot re-run "
    "them and says so.",
    "Pair-record precondition (amendment A1): every one of the ten pair records must carry complete: true and the keys "
    "row_keys, held_engine_raw and train_engine_raw; otherwise --score refuses, naming the record and the key.",
]


def rng_for(name: str, seed: int = SEED) -> np.random.Generator:
    """A fixed, named random stream (one per control), so a control's draws do not depend on the order of the rest."""
    return np.random.default_rng(np.random.SeedSequence([seed, zlib.crc32(name.encode("utf-8"))]))


# ----------------------------------------------------------------------------------------------------------------
# 1. Classes and layouts (row keys, groups, classes; no residual)
# ----------------------------------------------------------------------------------------------------------------

def class_of_method(method: str) -> str:
    """ABS (scale tau 0: sida, calibrated), FIG (semi_figure, a finite figure scale), FREE (scale tau inf)."""
    t = float(sm3.SCALE_TAU[method])
    if t == 0.0:
        return ABS
    if math.isinf(t):
        return FREE
    if method == FIG_METHOD:
        return FIG
    raise ValueError(f"method class {method!r} (scale tau {t}) is none of H0's three classes")


class Layout:
    """The level groups (fit_thiol_lolo.level_of: paper x species for ABS, the series for FIG and FREE) of a row set,
    with each group's method and class, as an n x G indicator."""

    def __init__(self, group_keys: Sequence[str], methods: Sequence[str]):
        self.n = len(group_keys)
        self.keys = list(dict.fromkeys(group_keys))
        index = {k: j for j, k in enumerate(self.keys)}
        self.G = len(self.keys)
        self.gid = np.array([index[k] for k in group_keys], dtype=int)
        self.M = np.zeros((self.n, self.G))
        if self.n:
            self.M[np.arange(self.n), self.gid] = 1.0
        self.size = self.M.sum(axis=0)
        method_g: List[Optional[str]] = [None] * self.G
        for g, m in zip(self.gid, methods):
            if method_g[g] is None:
                method_g[g] = m
            elif method_g[g] != m:
                raise ValueError(f"level group {self.keys[g]!r} mixes method classes")
        self.method_g = np.array(method_g, dtype=object)
        self.class_g = np.array([class_of_method(m) for m in method_g], dtype=object)
        self.row_method = np.array(list(methods), dtype=object)
        self.row_class = self.class_g[self.gid] if self.n else np.array([], dtype=object)

    @classmethod
    def of_rows(cls, rows: Sequence[Any]) -> "Layout":
        return cls([level_of(r) for r in rows], [fj.method_of(r) for r in rows])


@dataclasses.dataclass
class Block:
    """One laboratory's rows in one calibration context: layout, base variance s_i^2 (row sigma + the fit's class
    noise, no omega) and residuals r (R x n, predicted minus measured, dex)."""
    lab: str
    layout: Layout
    s2: np.ndarray
    r: np.ndarray
    rows: tuple = ()
    note: str = ""


# ----------------------------------------------------------------------------------------------------------------
# 2. The estimator (sec. 2)
# ----------------------------------------------------------------------------------------------------------------

def within_moments(lay: Layout, s2: np.ndarray, r: np.ndarray):
    """Sec. 2 step 1, per row of a group with >= 2 rows: centred_i = r_i - sum(w r)/W, e_i = centred_i^2 - (s_i^2 -
    1/W), c_i = 1 - 2 w_i/W + sum_j w_j^2/W^2 (w = 1/s^2), so E[e_i] = c_i omega^2. Returns (centred, e, c, rows in
    multi-row groups)."""
    r = np.atleast_2d(np.asarray(r, dtype=float))
    w = 1.0 / s2
    W = w @ lay.M
    sel = lay.size[lay.gid] >= 2 if lay.n else np.zeros(0, dtype=bool)
    Wsafe = np.where(W > 0, W, 1.0)
    mean = ((r * w) @ lay.M) / Wsafe
    centred = r - mean[:, lay.gid]
    Wr = Wsafe[lay.gid]
    e = centred ** 2 - (s2 - 1.0 / Wr)
    c = 1.0 - 2.0 * w / Wr + ((w ** 2) @ lay.M)[lay.gid] / Wr ** 2
    return centred, e, c, sel


def lab_omega2(block: Block) -> Optional[np.ndarray]:
    """omega_l^2 = sum_i e_i / sum_i c_i over the laboratory's rows in groups with >= 2 rows; None if it has none."""
    if not block.layout.n:
        return None
    _, e, c, sel = within_moments(block.layout, block.s2, block.r)
    if not sel.any():
        return None
    return e[:, sel].sum(axis=1) / float(c[sel].sum())


def class_d(block: Block, omega2: np.ndarray, cls: str, methods: Optional[Sequence[str]] = None,
            multi_only: bool = False) -> Optional[np.ndarray]:
    """d_g = m_g^2 - v_g for the laboratory's class-c groups (R x G_c), with s~^2 = s^2 + omega^2(S); None if none."""
    lay = block.layout
    if not lay.n:
        return None
    gsel = lay.class_g == cls
    if methods is not None:
        gsel &= np.isin(lay.method_g, list(methods))
    if multi_only:
        gsel &= lay.size >= 2
    if not gsel.any():
        return None
    st2 = block.s2[None, :] + np.asarray(omega2, dtype=float).reshape(-1, 1)
    winv = 1.0 / st2
    Ms = lay.M[:, gsel]
    W = winv @ Ms
    m = ((block.r * winv) @ Ms) / W
    return m ** 2 - 1.0 / W


def estimate(blocks: Sequence[Block], R: int, fallback_abs: Optional[Mapping[str, float]] = None,
             omega2_fixed: Optional[float] = None, abs_methods: Optional[Sequence[str]] = None,
             abs_multi_only: bool = False) -> Dict[str, Any]:
    """Sec. 2: omega^2(S) = max(0, mean_l omega_l^2); tau_c^2(S) = max(0, mean_{l in S_c} mean_g d_g) for c in ABS,
    FIG; the minimum-data rule (>= 2 groups in S; one laboratory: 'single-lab', flagged; none or < 2 groups: EX-056's
    shipped in-sample tau for the fold, flagged); tau_FIG,used = max(0.30, tau_FIG). ``taus`` maps each method to
    the level width the nested predictive gives it (R,)."""
    blocks = [b for b in blocks if b.layout.n]
    flags: List[str] = []
    per_lab_omega = {b.lab: lab_omega2(b) for b in blocks}
    have = {lab: v for lab, v in per_lab_omega.items() if v is not None}
    if omega2_fixed is not None:
        omega2_raw = np.full(R, float(omega2_fixed))
    elif have:
        omega2_raw = np.mean(np.stack(list(have.values())), axis=0)
    else:
        omega2_raw = np.zeros(R)
        flags.append("no group with >= 2 rows in S: omega = 0")
    omega2 = np.maximum(0.0, omega2_raw)
    out: Dict[str, Any] = {"omega2_raw": omega2_raw, "omega2": omega2, "omega": np.sqrt(omega2),
                           "omega_labs": sorted(have), "labs": sorted(b.lab for b in blocks)}
    taus: Dict[str, np.ndarray] = {}
    for cls in (ABS, FIG):
        per_lab = {}
        for b in blocks:
            d = class_d(b, omega2, cls, abs_methods if cls == ABS else None, abs_multi_only if cls == ABS else False)
            if d is not None:
                per_lab[b.lab] = d
        n_groups = int(sum(d.shape[1] for d in per_lab.values()))
        rec: Dict[str, Any] = {"labs": sorted(per_lab), "n_labs": len(per_lab), "n_groups": n_groups,
                               "per_lab_dbar": {lab: d.mean(axis=1) for lab, d in per_lab.items()},
                               "per_lab_d": per_lab}
        if n_groups < 2 or not per_lab:
            rec.update({"status": "fallback", "tau2_raw": None, "tau": None})
            flags.append(f"{cls}: {len(per_lab)} lab(s), {n_groups} group(s) in S -> EX-056's shipped in-sample tau")
        else:
            raw = np.mean(np.stack([d.mean(axis=1) for d in per_lab.values()]), axis=0)
            rec.update({"status": "normal" if len(per_lab) >= 2 else "single-lab", "tau2_raw": raw,
                        "tau": np.sqrt(np.maximum(0.0, raw))})
            if len(per_lab) == 1:
                flags.append(f"{cls}: single-lab tau ({next(iter(per_lab))})")
        out[cls] = rec
    abs_methods_all = [m for m, t in sm3.SCALE_TAU.items() if float(t) == 0.0]
    if out[ABS]["status"] == "fallback":
        for m in abs_methods_all:
            v = (fallback_abs or {}).get(m)
            taus[m] = np.full(R, float(v) if v is not None else float("nan"))
    else:
        for m in abs_methods_all:
            taus[m] = out[ABS]["tau"]
    if out[FIG]["status"] == "fallback":
        out["tau_fig_used"] = np.full(R, float(sm3.SCALE_TAU[FIG_METHOD]))
    else:
        out["tau_fig_used"] = np.maximum(FIG_FLOOR, out[FIG]["tau"])
    taus[FIG_METHOD] = out["tau_fig_used"]
    out["taus"] = taus
    out["flags"] = flags
    return out


def paule_mandel_tau(est: Dict[str, Any], blocks: Sequence[Block], cls: str = ABS) -> Optional[float]:
    """Sensitivity (i) (R = 1): the lab-weighted Paule-Mandel root mean_l mean_g m_g^2 / (v_g + tau^2) = 1."""
    from scipy.optimize import brentq
    omega2 = est["omega2"]
    parts = []
    for b in blocks:
        lay = b.layout
        if not lay.n:
            continue
        gsel = lay.class_g == cls
        if not gsel.any():
            continue
        st2 = b.s2 + float(omega2[0])
        W = (1.0 / st2) @ lay.M[:, gsel]
        m = ((b.r[0] / st2) @ lay.M[:, gsel]) / W
        parts.append((m ** 2, 1.0 / W))
    if not parts:
        return None

    def f(t2):
        return float(np.mean([np.mean(m2 / (v + t2)) for m2, v in parts])) - 1.0

    if f(0.0) <= 0.0:
        return 0.0
    hi = 1.0
    while f(hi) > 0.0 and hi < 1e6:
        hi *= 4.0
    return math.sqrt(brentq(f, 0.0, hi, xtol=1e-14))


def row_weighted_tau_blocks(est: Dict[str, Any], blocks: Sequence[Block], cls: str = ABS) -> Optional[float]:
    """Sensitivity (ii) (R = 1): tau^2 = max(0, sum_g n_g d_g / sum_g n_g) over every class-c group of S."""
    num = den = 0.0
    for b in blocks:
        d = class_d(b, est["omega2"], cls)
        if d is None:
            continue
        sizes = b.layout.size[b.layout.class_g == cls]
        num += float(np.sum(sizes * d[0]))
        den += float(np.sum(sizes))
    return math.sqrt(max(0.0, num / den)) if den else None


# ----------------------------------------------------------------------------------------------------------------
# 3. Arms and scores (vectorised over replicates; the real run cross-checks against MeasurementModel)
# ----------------------------------------------------------------------------------------------------------------

def gaussian_arm(lay: Layout, s2: np.ndarray, r: np.ndarray, omega2, taus: Mapping[str, Any]) -> Dict[str, Any]:
    """The Gaussian predictive of sec. 5 for one laboratory's rows: s~^2 = s^2 + omega^2; ABS / FIG groups N(0, tau^2)
    levels (tau per method from ``taus``), FREE groups free. Per group the joint log density (MeasurementModel.
    log_density's formulas), per row the centred residual and its predictive variance (fit_thiol_lolo.row_intervals'
    formulas; NaN for a single-row FREE group), the level scale, and the per-row marginal log density."""
    r = np.atleast_2d(np.asarray(r, dtype=float))
    R = r.shape[0]
    om = np.broadcast_to(np.asarray(omega2, dtype=float).reshape(-1, 1), (R, 1))
    st2 = np.broadcast_to(s2[None, :] + om, r.shape)
    winv = 1.0 / st2
    M = lay.M
    A = (r * r * winv) @ M
    B = (r * winv) @ M
    W = winv @ M
    Lg = np.log(st2) @ M
    free = lay.class_g == FREE
    t2 = np.zeros((R, lay.G))
    for m in sorted({x for x in lay.method_g[~free]}):
        t = np.broadcast_to(np.asarray(taus[m], dtype=float).reshape(-1), (R,))
        t2[:, lay.method_g == m] = (t ** 2)[:, None]
    Wsafe = np.where(W > 0, W, 1.0)
    Q = np.where(free, A - B * B / Wsafe, A - t2 * B * B / (1.0 + t2 * W))
    logdet = np.where(free, Lg + np.log(Wsafe), Lg + np.log1p(t2 * W))
    neff = np.where(free, lay.size - 1.0, lay.size)
    lpd_g = -0.5 * (Q + logdet + neff * LOG2PI)
    fr = free[lay.gid] if lay.n else np.zeros(0, dtype=bool)
    centred = np.array(r, dtype=float)
    var = st2 + t2[:, lay.gid]
    if fr.any():
        mean = (B / Wsafe)[:, lay.gid]
        centred = np.where(fr, r - mean, r)
        var = np.where(fr, st2 - 1.0 / Wsafe[:, lay.gid], var)
        single = fr & (lay.size[lay.gid] < 2)
        centred[:, single] = np.nan
        var[:, single] = np.nan
    scale = np.where(fr, np.nan, np.sqrt(t2[:, lay.gid]))
    with np.errstate(invalid="ignore"):
        marg = -0.5 * (centred ** 2 / var + np.log(2.0 * math.pi * var))
    return {"lpd_g": lpd_g, "lpd": lpd_g.sum(axis=1), "centred": centred, "var": var, "st2": np.array(st2),
            "scale": scale, "marg": marg, "layout": lay}


def interval_summary(centred: np.ndarray, var: np.ndarray, levels=LEVELS, inside_of: Optional[Dict] = None,
                     half_width90: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Coverage counts per level, the 90 % interval score (Gneiting & Raftery 2007: 2h + 20 max(0, |c| - h)) and the
    mean / median full 90 % width, over the scored rows (finite variance). ``inside_of`` / ``half_width90`` replace
    the Gaussian rule (St, St4)."""
    from scipy.stats import norm
    ok = np.isfinite(var[0]) if var.size else np.zeros(0, dtype=bool)
    n = int(ok.sum())
    out: Dict[str, Any] = {"scored": n}
    for lv in levels:
        if inside_of is not None:
            inside = inside_of[lv][:, ok]
        else:
            q = float(norm.ppf(0.5 + lv / 2.0))
            inside = np.abs(centred[:, ok]) <= q * np.sqrt(var[:, ok])
        out[f"{lv:.2f}"] = inside.sum(axis=1)
    h = half_width90[:, ok] if half_width90 is not None else float(norm.ppf(0.95)) * np.sqrt(var[:, ok])
    ab = np.abs(centred[:, ok])
    out["is90"] = (2.0 * h + 20.0 * np.maximum(0.0, ab - h)).mean(axis=1) if n else np.full(centred.shape[0], np.nan)
    out["width90_mean"] = (2.0 * h).mean(axis=1) if n else np.full(centred.shape[0], np.nan)
    out["width90_median"] = np.median(2.0 * h, axis=1) if n else np.full(centred.shape[0], np.nan)
    return out


@functools.lru_cache(maxsize=None)
def _gh(n: int):
    """Gauss-Hermite nodes for E over N(0, 1) (numpy's weights overflow above ~150 nodes, so capped at 120)."""
    z, w = np.polynomial.hermite.hermgauss(min(n, 120))
    return z, w / math.sqrt(math.pi)


@functools.lru_cache(maxsize=None)
def _gl_prob(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


@functools.lru_cache(maxsize=None)
def _t_quantiles(df: float, n: int) -> np.ndarray:
    from scipy.special import stdtrit
    return stdtrit(df, _gl_prob(n)[0])


def st_cdf_abs(x, sd, scale, df: float, nodes: int = ST_NODES) -> np.ndarray:
    """P(|b + e| <= x) for b ~ scale * t_df, e ~ N(0, sd^2) (St's row marginal). Gauss-Hermite over e where scale >=
    sd (the t CDF is then smooth on the error's scale), Gauss-Legendre over the t quantile otherwise."""
    from scipy.special import ndtr, stdtr
    x = np.asarray(x, dtype=float).ravel()
    sd = np.asarray(sd, dtype=float).ravel()
    scale = np.asarray(scale, dtype=float).ravel()
    out = np.empty(x.size)
    if not np.isfinite(df):
        tot = np.sqrt(sd ** 2 + scale ** 2)
        return 2.0 * ndtr(x / tot) - 1.0
    z, wz = _gh(nodes)
    p, wp = _gl_prob(nodes)
    q = _t_quantiles(float(df), nodes)
    chunk = max(1, 400_000 // nodes)
    for s in range(0, x.size, chunk):
        sl = slice(s, s + chunk)
        xs, ss, cs = x[sl], sd[sl], scale[sl]
        res = np.empty(xs.size)
        a = cs >= ss
        if a.any():
            e = math.sqrt(2.0) * ss[a, None] * z[None, :]
            T = stdtr(df, (xs[a, None] - e) / cs[a, None]) - stdtr(df, (-xs[a, None] - e) / cs[a, None])
            res[a] = T @ wz
        b = ~a
        if b.any():
            bv = cs[b, None] * q[None, :]
            F = ndtr((xs[b, None] - bv) / ss[b, None]) - ndtr((-xs[b, None] - bv) / ss[b, None])
            res[b] = F @ wp
        out[sl] = res
    return out


def st_inside(arm: Dict[str, Any], df_row: np.ndarray, levels=LEVELS, nodes: int = ST_NODES) -> Dict[float, np.ndarray]:
    """St's coverage: a levelled row is inside its central p interval iff P(|X| <= |centred|) <= p; FREE rows keep
    the Gaussian rule."""
    from scipy.stats import norm
    c, st2, scale, var = arm["centred"], arm["st2"], arm["scale"], arm["var"]
    R, n = c.shape
    F = np.full((R, n), np.nan)
    levelled = np.isfinite(scale[0]) if n else np.zeros(0, dtype=bool)
    for df in sorted(set(df_row[levelled].tolist())):
        cols = levelled & (df_row == df)
        F[:, cols] = st_cdf_abs(np.abs(c[:, cols]), np.sqrt(st2[:, cols]), scale[:, cols], df, nodes
                                ).reshape(R, int(cols.sum()))
    out = {}
    for lv in levels:
        q = float(norm.ppf(0.5 + lv / 2.0))
        gauss = np.abs(c) <= q * np.sqrt(var)
        with np.errstate(invalid="ignore"):
            out[lv] = np.where(levelled[None, :], F <= lv, gauss)
    return out


def st_quantile(level: float, sd: float, scale: float, df: float, nodes: int = ST_NODES_REAL) -> float:
    """The half-width q with P(|b + e| <= q) = level (R = 1, the real run)."""
    from scipy.optimize import brentq
    from scipy.stats import t as tdist
    if not np.isfinite(df) or scale == 0.0:
        from scipy.stats import norm
        return float(norm.ppf(0.5 + level / 2.0)) * math.sqrt(sd * sd + scale * scale)
    hi = 10.0 * sd + 10.0 * scale * float(tdist.ppf(0.5 + level / 2.0, df))
    return brentq(lambda x: float(st_cdf_abs([x], [sd], [scale], df, nodes)[0]) - level, 0.0, hi, xtol=1e-12)


def st_group_logdensity(r: np.ndarray, st2: np.ndarray, scale: float, df: float, nodes: int = ST_NODES_REAL) -> float:
    """St's joint density of one ABS / FIG group: the rows' Gaussian likelihood integrated over a t_df(scale) level
    (1-D quadrature): -0.5(Q0 + sum log 2 pi s~^2) + 0.5 log(2 pi / W) + log E[N(b; m, 1/W)]."""
    from scipy.stats import t as tdist
    W = float(np.sum(1.0 / st2))
    m = float(np.sum(r / st2)) / W
    q0 = float(np.sum(r * r / st2)) - W * m * m
    base = -0.5 * (q0 + float(np.sum(np.log(2.0 * math.pi * st2)))) + 0.5 * math.log(2.0 * math.pi / W)
    sd_m = 1.0 / math.sqrt(W)
    if scale == 0.0 or not np.isfinite(df):
        v = sd_m ** 2 + scale ** 2
        return base - 0.5 * (math.log(2.0 * math.pi * v) + m * m / v)
    if sd_m <= scale:
        z, w = _gh(nodes)
        b = m + math.sqrt(2.0) * sd_m * z
        integral = float(np.sum(w * tdist.pdf(b / scale, df) / scale))
    else:
        p, w = _gl_prob(nodes)
        b = scale * _t_quantiles(float(df), nodes)
        integral = float(np.sum(w * np.exp(-0.5 * (b - m) ** 2 / sd_m ** 2) / (math.sqrt(2.0 * math.pi) * sd_m)))
    return base + math.log(integral)


# ----------------------------------------------------------------------------------------------------------------
# 4. Records (EX-056's folds, H0's pairs)
# ----------------------------------------------------------------------------------------------------------------

KEY_FIELDS = ("ledger_id", "species", "t_min", "method")


def row_keys(rows: Sequence[Any]) -> List[Dict[str, Any]]:
    return [pair_runner.row_key(r) for r in rows]


@dataclasses.dataclass
class Fold:
    """One EX-056 fold (laboratory L held out), columnar: residual arrays (R x n) keyed engine_<v> / null_<k>."""
    lab: str
    held: tuple
    paired: np.ndarray
    raw: Dict[str, np.ndarray]
    terms: Dict[str, Dict[str, Any]]
    leak_prone: List[str]
    record: Optional[Dict[str, Any]] = None


@dataclasses.dataclass
class Pair:
    """One H0 pair record (A, B held out), columnar."""
    labs: Tuple[str, str]
    held: Dict[str, tuple]
    train: tuple
    raw_held: Dict[str, Dict[str, np.ndarray]]
    raw_train: Dict[str, np.ndarray]
    converged: Dict[str, bool]
    leak_prone: List[str]
    freeze_commit: str = ""


def fold_path(lab: str, directory: Optional[Path] = None) -> Path:
    return Path(directory or FOLD_DIR) / f"{lolo.STEM}_{lab}.json"


def pair_path(a: str, b: str, directory: Optional[Path] = None) -> Path:
    return Path(directory or PAIR_DIR) / f"{pair_runner.STEM}_{a}__{b}.json"


def load_fold(lab: str, held: Sequence[Any], directory: Optional[Path] = None, residuals: bool = True) -> Fold:
    """An EX-056 fold record; its rows must equal make_folds' held-out rows of ``lab`` by keys AND order.
    ``residuals=False`` (the design, the controls) leaves every residual out of memory."""
    rec = json.loads(fold_path(lab, directory).read_text(encoding="utf-8"))
    keys = [{k: x[k] for k in KEY_FIELDS} for x in rec["rows"]]
    if keys != row_keys(held):
        raise SystemExit(f"{lab}: the fold record's rows differ from make_folds' held-out rows (keys or order)")
    leak = list(rec.get("leak_prone") or [])
    variants = ["warm"] + (["reset"] if leak else [])
    for x in rec["rows"]:
        assert sorted(x["engine_raw"]) == sorted(variants), f"{lab}: engine variants differ from leak_prone"
    paired = np.array([bool(x["paired"]) for x in rec["rows"]], dtype=bool)
    terms = {name: m["terms"] for name, m in rec["models"].items() if "terms" in m}
    raw: Dict[str, np.ndarray] = {}
    if residuals:
        for v in variants:
            raw[f"engine_{v}"] = np.array([[float(x["engine_raw"][v]) for x in rec["rows"]]], dtype=float)
        for k in NULLS:
            raw[f"null_{k}"] = np.array([[float(x["null_raw"][k]) for x in rec["rows"]]], dtype=float)
    return Fold(lab=lab, held=tuple(held), paired=paired, raw=raw, terms=terms, leak_prone=leak,
                record=rec if residuals else None)


#: the keys every pair record must carry before --score reads it (amendment A1)
PAIR_RECORD_KEYS = ("row_keys", "held_engine_raw", "train_engine_raw")


def check_pair_records(directory: Optional[Path] = None) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """--score's precondition (blinding): all ten pair records exist, are not guard records, are complete, carry the
    keys row_keys / held_engine_raw / train_engine_raw, and each carries a non-empty freeze commit. Otherwise
    SystemExit (non-zero), naming what is missing."""
    problems, recs = [], {}
    for a, b in pair_runner.pairs():
        p = pair_path(a, b, directory)
        if not p.exists():
            problems.append(f"missing {p.name}")
            continue
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            problems.append(f"{p.name}: unreadable ({exc})")
            continue
        if not isinstance(rec.get("freeze_commit"), str) or not rec["freeze_commit"].strip():
            problems.append(f"{p.name}: no freeze_commit")
        if rec.get("guard"):
            problems.append(f"{p.name}: a guard record")
        if list(rec.get("labs") or []) != [a, b]:
            problems.append(f"{p.name}: labs {rec.get('labs')} != {[a, b]}")
        if rec.get("complete") is not True:
            problems.append(f"{p.name}: not complete (complete: {rec.get('complete')!r}; the runner has not finished it)")
        if int(rec.get("failed_pot_runs") or 0) > 0:
            # amendment A2: a pot run that failed (in practice a wall-clock PotTimeout under a paused or
            # memory-starved machine) enters the fit as a penalty residual, so the record depends on machine state
            problems.append(f"{p.name}: {rec.get('failed_pot_runs')} failed pot runs (amendment A2: re-run it from "
                            "the freeze commit, same settings, unpaused)")
        absent = [k for k in PAIR_RECORD_KEYS if k not in rec]
        if absent:
            problems.append(f"{p.name}: missing key(s) {absent}")
        recs[(a, b)] = rec
    if problems:
        raise SystemExit("--score refuses to run (prereg sec. 7: the ten pair records, each with its freeze commit, "
                         "come first): " + "; ".join(problems))
    commits = {r["freeze_commit"] for r in recs.values()}
    if len(commits) != 1:
        raise SystemExit(f"--score refuses to run: the pair records carry different freeze commits {sorted(commits)}")
    return recs


def load_pair(a: str, b: str, data: Sequence[Any], rec: Optional[Dict[str, Any]] = None,
              directory: Optional[Path] = None) -> Pair:
    """A pair record; its held and training rows must equal fit_thiol_lolo_pairs.split's by keys AND order."""
    rec = rec if rec is not None else json.loads(pair_path(a, b, directory).read_text(encoding="utf-8"))
    train, held = pair_runner.split(data, (a, b))
    if rec["row_keys"]["train"] != row_keys(train):
        raise SystemExit(f"pair {a}+{b}: training rows differ from data minus both laboratories")
    for lab in (a, b):
        if rec["row_keys"]["held"][lab] != row_keys(held[lab]):
            raise SystemExit(f"pair {a}+{b}: held-out rows of {lab} differ (keys or order)")
    leak = list(rec.get("leak_prone") or [])
    variants = ["warm"] + (["reset"] if leak else [])
    if sorted(rec["held_engine_raw"]) != sorted(variants) or sorted(rec["train_engine_raw"]) != sorted(variants):
        raise SystemExit(f"pair {a}+{b}: variants {sorted(rec['held_engine_raw'])} do not match leak_prone")
    raw_held = {v: {lab: np.array([[float(x) for x in rec["held_engine_raw"][v][lab]]], dtype=float) for lab in (a, b)}
                for v in variants}
    for v in variants:
        for lab in (a, b):
            assert raw_held[v][lab].shape[1] == len(held[lab]), f"pair {a}+{b} {v} {lab}: misaligned residuals"
    raw_train = {v: np.array([float(x) for x in rec["train_engine_raw"][v]], dtype=float) for v in variants}
    for v in variants:
        assert raw_train[v].size == len(train), f"pair {a}+{b} {v}: misaligned training residuals"
    converged = {v: pair_runner.converged(rec[v]["fit"]) for v in variants}
    return Pair(labs=(a, b), held=held, train=train, raw_held=raw_held, raw_train=raw_train, converged=converged,
                leak_prone=leak, freeze_commit=str(rec.get("freeze_commit") or ""))


def null_answers(train: Sequence[Any], rows: Sequence[Any]) -> np.ndarray:
    """Which rows a null fitted on ``train`` can answer (fit_thiol_lolo.null_residuals' NaN condition: the species has
    a positive training row, and the row's own value is positive) -- from row keys and species presence only."""
    species = {x.species for x in train if x.value > 0}
    return np.array([(x.species in species and x.value > 0) for x in rows], dtype=bool)


def sigma2(terms: Mapping[str, Any], rows: Sequence[Any]) -> np.ndarray:
    """s_i^2 = sigma_row^2 + the fit's pooled class noise^2 (fit_thiol_lolo.predictive_model, its fallbacks kept)."""
    if not rows:
        return np.zeros(0)
    pm, _ = predictive_model(terms, rows)
    return np.array([pm.sigma(x) ** 2 for x in rows], dtype=float)


def shipped_taus(terms: Mapping[str, Any], rows: Sequence[Any]) -> Dict[str, float]:
    """EX-056's level widths for these rows: the training tau per absolute method (predictive_model's fallback for a
    method with none) and the declared figure width."""
    if rows:
        _, info = predictive_model(terms, rows)
        taus = dict(info["between_tau_used"])
    else:
        taus = dict(terms["between_tau"])
    for m, t in sm3.SCALE_TAU.items():
        if float(t) == 0.0:
            taus.setdefault(m, max(taus.values()) if taus else 0.0)
    taus[FIG_METHOD] = float(sm3.SCALE_TAU[FIG_METHOD])
    return taus


# ----------------------------------------------------------------------------------------------------------------
# 5. Nesting (sec. 3)
# ----------------------------------------------------------------------------------------------------------------

class Nesting:
    """The strict code path (sec. 3.1 - 3.3) and the one-deep sensitivity (3.4) over a set of fold and pair records,
    real or synthetic. ``pair_terms(pair, variant, lab)`` gives a pair fit's training terms (real: training_terms on
    its training residuals, n_constants 16); ``null_fn(kind, train, rows)`` the null residuals (real:
    fit_thiol_lolo.null_residuals); ``null_terms(kind, pair)`` the strict null's training terms."""

    def __init__(self, data: Sequence[Any], folds: Dict[str, Fold], pairs: Dict[Tuple[str, str], Pair], pots_all,
                 pair_terms: Optional[Callable] = None, null_fn: Optional[Callable] = None,
                 null_terms: Optional[Callable] = None, check_null_masks: bool = True,
                 design_cache: Optional[Dict[Any, Any]] = None):
        self.data, self.folds, self.pairs, self.pots_all = tuple(data), folds, pairs, pots_all
        self._pair_terms = pair_terms or self._real_pair_terms
        self._null_fn = null_fn or (lambda kind, train, rows: null_residuals(kind, train, rows, self.pots_all))
        self._null_terms = null_terms or self._real_null_terms
        self.check_null_masks = check_null_masks
        self._cache: Dict[Any, Any] = {}
        #: layouts and s_i^2 depend on row keys and saved terms only, never on residuals: synthetic worlds on one
        #: design share them
        self._design_cache: Dict[Any, Any] = design_cache if design_cache is not None else {}

    # -- the real providers
    def _real_pair_terms(self, pair: Pair, variant: str, lab: str) -> Dict[str, Any]:
        key = ("pair_terms", pair.labs, variant)
        if key not in self._cache:
            self._cache[key] = training_terms(pair.raw_train[variant], pair.train, N_CONSTANTS_ENGINE)
        return self._cache[key]

    def _real_null_terms(self, kind: str, pair: Pair) -> Dict[str, Any]:
        key = ("null_terms", kind, pair.labs)
        if key not in self._cache:
            raw, n = self._null_fn(kind, pair.train, pair.train)
            self._cache[key] = training_terms(raw, pair.train, n)
        return self._cache[key]

    def _null_held(self, kind: str, pair: Pair, lab: str) -> np.ndarray:
        key = ("null_held", kind, pair.labs, lab)
        if key not in self._cache:
            raw, _ = self._null_fn(kind, pair.train, pair.held[lab])
            self._cache[key] = np.atleast_2d(np.asarray(raw, dtype=float))
        return self._cache[key]

    # -- keys
    def fold_key(self, lab: str, model: str, set_name: str) -> str:
        if model == "engine":
            return "engine_reset" if set_name == SETS[1] and self.folds[lab].leak_prone else "engine_warm"
        return f"null_{model}"

    @staticmethod
    def pair_variant(pair: Pair, set_name: str) -> str:
        return "reset" if set_name == SETS[1] and "reset" in pair.raw_held else "warm"

    def pair_of(self, a: str, b: str) -> Pair:
        return self.pairs[tuple(sorted((a, b)))]

    def _layout(self, key, rows) -> Layout:
        k = ("layout", key)
        if k not in self._design_cache:
            self._design_cache[k] = Layout.of_rows(rows)
        return self._design_cache[k]

    def _s2(self, key, terms, rows) -> np.ndarray:
        k = ("s2", key, id(terms))
        if k not in self._design_cache:
            self._design_cache[k] = (terms, sigma2(terms, rows))   # terms kept alive so its id stays unique
        return self._design_cache[k][1]

    # -- blocks
    def outer(self, lab: str, model: str, set_name: str, every_finite: bool = False) -> Block:
        """L's scored rows (EX-056's paired mask; ``every_finite``: every finite engine row, D-all) with fold L's
        residuals and fold L's class noise."""
        fold = self.folds[lab]
        key = self.fold_key(lab, model, set_name)
        r_all = fold.raw[key]
        mask = np.all(np.isfinite(r_all), axis=0) if every_finite else fold.paired
        rows = tuple(x for x, m in zip(fold.held, mask) if m)
        mk = ("outer", lab, key, mask.tobytes())
        return Block(lab=lab, layout=self._layout(mk, rows), s2=self._s2(mk, fold.terms[key], rows),
                     r=r_all[:, mask], rows=rows)

    def inner_mask(self, outer: str, inner: str) -> np.ndarray:
        """Sec. 3.2: L''s rows that are finite under the strict engine fit (every variant of the pair record) AND
        that both strict nulls (trained without L and L') can answer."""
        key = ("mask", outer, inner)
        if key not in self._cache:
            pair = self.pair_of(outer, inner)
            finite = np.ones(len(pair.held[inner]), dtype=bool)
            for v in pair.raw_held:
                finite &= np.all(np.isfinite(pair.raw_held[v][inner]), axis=0)
            answers = null_answers(pair.train, pair.held[inner])
            if self.check_null_masks:
                for kind in NULLS:
                    got = np.all(np.isfinite(self._null_held(kind, pair, inner)), axis=0)
                    assert np.array_equal(got, answers), f"null_answers disagrees with null_residuals ({kind})"
            self._cache[key] = finite & answers
        return self._cache[key]

    def strict_inner(self, outer: str, inner: str, model: str, set_name: str
                     ) -> Tuple[Optional[Block], str, str]:
        """Sec. 3.1: inner L''s residuals from the fit trained without L AND L' (the pair record; for a null, the null
        retrained on data minus both). Returns (block, note, kind); kind is '' (a block), 'not_converged' (the
        engine's pair fit did not converge: caps the verdict, sec. 3.1) or 'mask_emptied' (every row removed by the
        strict masks: the minimum-data rule, flagged, sec. 2)."""
        pair = self.pair_of(outer, inner)
        mask = self.inner_mask(outer, inner)
        rows = tuple(x for x, m in zip(pair.held[inner], mask) if m)
        if model == "engine":
            v = self.pair_variant(pair, set_name)
            if not pair.converged[v]:
                return None, f"pair {'+'.join(pair.labs)} {v}: not converged (status <= 0)", "not_converged"
            r = pair.raw_held[v][inner][:, mask]
            terms = self._pair_terms(pair, v, inner)
            mk = ("strict", pair.labs, inner, "engine", v, mask.tobytes())
        else:
            r = self._null_held(model, pair, inner)[:, mask]
            terms = self._null_terms(model, pair)
            mk = ("strict", pair.labs, inner, model, mask.tobytes())
        if not rows:
            return None, f"{inner}: no row left by the strict masks", "mask_emptied"
        return Block(lab=inner, layout=self._layout(mk, rows), s2=self._s2(mk, terms, rows), r=r, rows=rows), "", ""

    def onedeep_inner(self, outer: str, inner: str, model: str, set_name: str) -> Tuple[Optional[Block], str, str]:
        """Sec. 3.4: inner L''s residuals from EX-056's fold L' (trained WITH L), on the STRICT inner mask of sec. 3.2
        (amendment A1) intersected with fold L''s paired rows, so that one-deep and strict differ by the leak only."""
        fold = self.folds[inner]
        key = self.fold_key(inner, model, set_name)
        mask = fold.paired & self.inner_mask(outer, inner)
        rows = tuple(x for x, m in zip(fold.held, mask) if m)
        if not rows:
            return None, f"{inner}: no row left by the strict masks (one-deep)", "mask_emptied"
        mk = ("onedeep", outer, inner, key, mask.tobytes())
        return Block(lab=inner, layout=self._layout(mk, rows), s2=self._s2(mk, fold.terms[key], rows),
                     r=fold.raw[key][:, mask], rows=rows), "", ""

    def calibrate(self, outer: str, model: str, set_name: str, nesting: str = "strict",
                  inner_labs: Optional[Sequence[str]] = None, **kw) -> Dict[str, Any]:
        """omega(S), tau_ABS(S), tau_FIG(S) for outer L from the inner laboratories (every scorable laboratory but L).
        A missing inner laboratory leaves S, and the minimum-data rule applies. Only a non-converged pair fit is
        listed in ``nonconverged_inner`` (it caps the verdict, sec. 3.1); a laboratory emptied by the strict masks is
        flagged (sec. 2) and caps nothing (amendment A1)."""
        labs = [x for x in (inner_labs or SCORABLE) if x != outer]
        blocks, notes, nonconverged, emptied = [], [], [], []
        for lp in labs:
            if nesting == "strict":
                blk, note, kind = self.strict_inner(outer, lp, model, set_name)
            elif nesting == "onedeep":
                blk, note, kind = self.onedeep_inner(outer, lp, model, set_name)
            else:
                raise ValueError(nesting)
            if blk is None:
                notes.append(note)
                (nonconverged if kind == "not_converged" else emptied).append(note)
            else:
                blocks.append(blk)
        fold = self.folds[outer]
        key = self.fold_key(outer, model, set_name)
        fallback = shipped_taus(fold.terms[key], tuple(x for x, m in zip(fold.held, fold.paired) if m))
        R = max([b.r.shape[0] for b in blocks] + [1])
        est = estimate(blocks, R, fallback_abs=fallback, **kw)
        est["flags"] += [f"{n} -> left out of S; the minimum-data rule applies (sec. 2)" for n in emptied]
        est["missing_inner"] = notes
        est["nonconverged_inner"] = nonconverged
        est["blocks"] = blocks
        return est

    def calibrate_all(self, model: str, set_name: str, labs: Sequence[str] = SCORABLE, **kw) -> Dict[str, Any]:
        """The ship-rule quantities (sec. 3.6): every laboratory's EX-056 single-out residuals."""
        blocks = [self.outer(lab, model, set_name) for lab in labs]
        R = max([b.r.shape[0] for b in blocks] + [1])
        est = estimate(blocks, R, **kw)
        est["blocks"] = blocks
        return est


def arm_b(block: Block, est: Dict[str, Any]) -> Dict[str, Any]:
    return gaussian_arm(block.layout, block.s2, block.r, est["omega2"], est["taus"])


def arm_a(block: Block, taus: Mapping[str, float]) -> Dict[str, Any]:
    return gaussian_arm(block.layout, block.s2, block.r, 0.0, taus)


def df_rows(block: Block, est: Dict[str, Any]) -> np.ndarray:
    """St's df per row: the number of laboratories informing the row's class (inf: a fallback, Gaussian)."""
    out = np.full(block.layout.n, np.inf)
    for cls in (ABS, FIG):
        rec = est[cls]
        if rec["status"] != "fallback":
            out[block.layout.row_class == cls] = float(rec["n_labs"])
    return out


def unstable_flag(full_tau: np.ndarray, drop_taus: Sequence[np.ndarray]) -> np.ndarray:
    """Sec. 3.7: t~ = max(tau, 0.10); UNSTABLE if max t~ / min t~ > 2 over the drop-one estimates, or a drop-one
    estimate lies on the other side of 0.8 dex from the full estimate."""
    t = np.maximum(np.stack(drop_taus), RATIO_FLOOR)
    ratio = t.max(axis=0) / t.min(axis=0)
    side = np.stack([(d > D2_LINE) != (full_tau > D2_LINE) for d in drop_taus]).any(axis=0)
    return (ratio > 2.0) | side


def abs_labs(nest: Nesting) -> List[str]:
    """The scorable laboratories holding ABS groups on their paired rows."""
    out = []
    for lab in SCORABLE:
        fold = nest.folds[lab]
        rows = [x for x, m in zip(fold.held, fold.paired) if m]
        if any(class_of_method(fj.method_of(x)) == ABS for x in rows):
            out.append(lab)
    return out


# ----------------------------------------------------------------------------------------------------------------
# 6. Statistics: sign flip, exact bootstrap, D3, D1
# ----------------------------------------------------------------------------------------------------------------

def sign_patterns(K: int) -> np.ndarray:
    return np.array(list(itertools.product((1.0, -1.0), repeat=K)))


def multisets(K: int) -> Tuple[np.ndarray, np.ndarray]:
    """Every multiset of K laboratories drawn from K (as count vectors) and its multinomial probability."""
    counts, weights = [], []
    for combo in itertools.combinations_with_replacement(range(K), K):
        k = np.bincount(combo, minlength=K)
        counts.append(k)
        weights.append(math.factorial(K) / float(np.prod([math.factorial(int(x)) for x in k])) / K ** K)
    return np.array(counts, dtype=float), np.array(weights)


def bootstrap_interval(d: np.ndarray, level: float = 0.90) -> Tuple[float, float]:
    counts, w = multisets(d.size)
    vals = counts @ d / d.size
    order = np.argsort(vals, kind="stable")
    cum = np.cumsum(w[order])
    lo = float(vals[order][np.searchsorted(cum, (1.0 - level) / 2.0 - 1e-15)])
    hi = float(vals[order][np.searchsorted(cum, 1.0 - (1.0 - level) / 2.0 - 1e-15)])
    return lo, hi


def d3_decision(d: np.ndarray, is_d: np.ndarray, threshold: float = SIGN_FLIP_THRESHOLD) -> Dict[str, Any]:
    """Sec. 6 D3 on per-laboratory paired differences (..., K): LPD (b) - (c_k) and IS (b) - (c_k). Vectorised over
    leading axes."""
    d = np.atleast_2d(d)
    is_d = np.atleast_2d(is_d)
    K = d.shape[-1]
    mean = d.mean(axis=-1)
    se = d.std(axis=-1, ddof=1) / math.sqrt(K)
    flipped = d @ sign_patterns(K).T / K
    eps = 1e-12 * np.maximum(1.0, np.abs(mean))
    p_better = (flipped >= (mean - eps)[:, None]).mean(axis=-1)
    p_worse = (flipped <= (mean + eps)[:, None]).mean(axis=-1)
    better = (mean > 2.0 * se) & (p_better <= threshold + 1e-15)
    null_better = (mean < -2.0 * se) & (p_worse <= threshold + 1e-15)
    is_mean = is_d.mean(axis=-1)
    trigger = better & (is_mean > 0.0)
    return {"mean": mean, "se": se, "p_better": p_better, "p_worse": p_worse, "lpd_pass": better,
            "null_better": null_better, "is_mean": is_mean, "sign_trigger": trigger}


def d3_label(res: Dict[str, Any], i: int = 0) -> str:
    if res["sign_trigger"][i]:
        return "INCONCLUSIVE (LPD passes, the interval score favours the null)"
    if res["lpd_pass"][i]:
        return "shown better"
    if res["null_better"][i]:
        return "the null is better"
    return "not shown better"


# ----------------------------------------------------------------------------------------------------------------
# 7. The reproduction guard (sec. 4.1)
# ----------------------------------------------------------------------------------------------------------------

def _rel(a, b) -> float:
    if a is None or b is None:
        return 0.0 if (a is None and b is None) else float("inf")
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    if not (math.isfinite(a) and math.isfinite(b)):
        return float("inf")
    return abs(a - b) / max(abs(a), abs(b))


def guard(data=None, folds_rows=None) -> Tuple[bool, float, List[str]]:
    """Rebuild, from the reloaded rows and every fold record's saved terms: per model x fold x set the per-row
    (centred, var), the coverage counts and the joint log_density on the paired rows; and the engine_*_all_rows
    log_density and coverage. Plus the curated inputs' sha256 against EX-056's provenance. Returns (pass, max relative
    difference, the names of failed checks)."""
    if data is None:
        data, _, _ = load_fit_table()
    folds_rows = folds_rows or make_folds(data)
    worst, failed = 0.0, []
    changed: Dict[str, str] = {}
    for lab, (_, held) in folds_rows.items():
        fold = load_fold(lab, held)
        rec = fold.record
        paired_rows = tuple(x for x, m in zip(held, fold.paired) if m)
        if int(rec["paired_rows"]) != len(paired_rows):
            failed.append(f"{lab}: paired_rows")
        for name, m in rec["models"].items():
            if name.endswith("_all_rows") or name.endswith("_sensitivity"):
                continue
            pm, _ = predictive_model(m["terms"], paired_rows)
            new = lolo.score(pm, fold.raw[name][0][fold.paired], paired_rows)
            for i, (p, q) in enumerate(zip(new["per_row"], m["per_row"])):
                worst = max(worst, _rel(p["centred"], q["centred"]), _rel(p["var"], q["var"]))
            if len(new["per_row"]) != len(m["per_row"]):
                failed.append(f"{lab}: {name} per_row length")
            for lv, c in m["coverage"].items():
                if (new["coverage"][lv]["inside"], new["coverage"][lv]["scored"]) != (c["inside"], c["scored"]):
                    failed.append(f"{lab}: {name} coverage {lv}")
            rel = _rel(new["log_density"], m["log_density"])
            worst = max(worst, rel)
            if rel > GUARD_TOL:
                failed.append(f"{lab}: {name} log_density")
        for v in ("warm", "reset"):
            name = f"engine_{v}_all_rows"
            if name not in rec["models"]:
                continue
            m = rec["models"][name]
            ok = np.isfinite(fold.raw[f"engine_{v}"][0])
            rows_ok = tuple(x for x, f in zip(held, ok) if f)
            pm, _ = predictive_model(rec["models"][f"engine_{v}"]["terms"], rows_ok)
            new = lolo.score(pm, fold.raw[f"engine_{v}"][0][ok], rows_ok)
            rel = _rel(new["log_density"], m["log_density"])
            worst = max(worst, rel)
            if rel > GUARD_TOL:
                failed.append(f"{lab}: {name} log_density")
            for lv, c in m["coverage"].items():
                if (new["coverage"][lv]["inside"], new["coverage"][lv]["scored"]) != (c["inside"], c["scored"]):
                    failed.append(f"{lab}: {name} coverage {lv}")
        for entry in rec["provenance"].get("inputs", []):
            path = (ROOT / entry["path"]).resolve()
            if path.is_relative_to(data_paths.DATA_ROOT) and provenance.sha256_of(path) != entry["sha256"]:
                changed[entry["path"]] = entry["sha256"]
    for rel_path, sha in sorted(changed.items()):
        if not _curated_input_at_ex056(rel_path, sha, data):
            failed.append(f"curated input {rel_path} changed since EX-056, and the fit table read with its "
                          f"{EX056_COMMIT} version differs")
    if worst > GUARD_TOL:
        failed.append("per-row or log_density relative difference above 1e-9")
    return (not failed), worst, sorted(set(failed))


#: the commit of EX-056's records (prereg secs. 0, 4.1)
EX056_COMMIT = "7b1bc7cd"


def _curated_input_at_ex056(rel_path: str, sha: str, data: Sequence[Any]) -> bool:
    """Sec. 4.1: a curated input whose sha256 differs from EX-056's is read at 7b1bc7cd. Implemented for the
    laboratory registry (the only curated input load_fit_table reads outside the fit-row tables): its 7b1bc7cd
    version must hash as EX-056 recorded, and the fit table rebuilt with it must equal today's row for row (every
    field, the laboratory included). Any other changed curated input fails the guard."""
    if (ROOT / rel_path).resolve() != Path(fit_rows.LABORATORIES).resolve():
        return False
    old = subprocess.run(["git", "show", f"{EX056_COMMIT}:{rel_path}"], cwd=str(ROOT), capture_output=True)
    if old.returncode != 0:
        return False
    import hashlib
    if hashlib.sha256(old.stdout).hexdigest() != sha:
        return False
    saved = fit_rows.LABORATORIES
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "laboratories.yml"
        p.write_bytes(old.stdout)
        try:
            fit_rows.LABORATORIES = p
            data_old, _, _ = load_fit_table()
        finally:
            fit_rows.LABORATORIES = saved
    return tuple(data_old) == tuple(data)


# ----------------------------------------------------------------------------------------------------------------
# 8. The design (sec. 7.1; no residual)
# ----------------------------------------------------------------------------------------------------------------

class Design:
    """The real design: the fit table's rows, EX-056's paired masks, its saved aggregate terms; no residual."""

    def __init__(self):
        data, _, pots_all = load_fit_table()
        self.data, self.pots_all = data, pots_all
        self.folds_rows = make_folds(data)
        self.labs = sorted(self.folds_rows)
        self.folds = {lab: load_fold(lab, held, residuals=False) for lab, (_, held) in self.folds_rows.items()}
        counts = {lab: len(held) for lab, (_, held) in self.folds_rows.items()}
        assert counts == TOY_COUNTS, f"the toy's declared row counts differ from the records: {counts}"
        self.counts = counts
        self.layouts = {lab: Layout.of_rows(f.held) for lab, f in self.folds.items()}
        self._s2: Dict[Tuple[str, str], np.ndarray] = {}
        self.nest_cache: Dict[Any, Any] = {}

    def terms_key(self, lab: str, set_name: str) -> str:
        return "engine_reset" if set_name == SETS[1] and self.folds[lab].leak_prone else "engine_warm"

    def terms(self, lab: str, set_name: str) -> Dict[str, Any]:
        return self.folds[lab].terms[self.terms_key(lab, set_name)]

    def s2(self, lab: str, set_name: str) -> np.ndarray:
        """Every held-out row's s_i^2 (row sigma + EX-056's saved class noise for that fold and set)."""
        key = (lab, set_name)
        if key not in self._s2:
            self._s2[key] = sigma2(self.terms(lab, set_name), self.folds[lab].held)
        return self._s2[key]

    def inner_mask(self, outer: str, inner: str) -> np.ndarray:
        train, held = pair_runner.split(self.data, (outer, inner))
        return self.folds[inner].paired & null_answers(train, held[inner])

    def table(self) -> Dict[str, Any]:
        rows_out = []
        for lab in self.labs:
            fold = self.folds[lab]
            lay = self.layouts[lab]
            for set_name in SETS:
                s = np.sqrt(self.s2(lab, set_name))
                for cls in (ABS, FIG, FREE):
                    sel_rows = (lay.row_class == cls) & fold.paired
                    if not sel_rows.any():
                        continue
                    pl = Layout.of_rows([x for x, m in zip(fold.held, sel_rows) if m])
                    rows_out.append({
                        "lab": lab, "set": set_name, "class": cls, "paired_rows": int(sel_rows.sum()),
                        "held_rows": int((lay.row_class == cls).sum()), "groups": pl.G,
                        "single_row_groups": int(np.sum(pl.size == 1)),
                        "methods": sorted(set(pl.method_g.tolist())),
                        "s_median": float(np.median(s[sel_rows])), "s_min": float(np.min(s[sel_rows])),
                        "s_max": float(np.max(s[sel_rows])), "terms": self.terms_key(lab, set_name)})
        inner = []
        for outer in SCORABLE:
            for lp in SCORABLE:
                if lp == outer:
                    continue
                m = self.inner_mask(outer, lp)
                lay = self.layouts[lp]
                inner.append({"outer": outer, "inner": lp,
                              **{cls: int(np.sum(m & (lay.row_class == cls))) for cls in (ABS, FIG, FREE)},
                              "lost_to_strict_nulls": int(np.sum(self.folds[lp].paired & ~m))})
        s_c = {}
        for outer in SCORABLE:
            s_c[outer] = {cls: sorted({x["inner"] for x in inner if x["outer"] == outer and x[cls] > 0})
                          for cls in (ABS, FIG, FREE)}
        return {"rows": rows_out, "inner_masks": inner, "inner_labs_per_class": s_c, "toy_counts": self.counts,
                "paired_rows_total": int(sum(int(f.paired.sum()) for f in self.folds.values()))}


def render_design(tab: Dict[str, Any]) -> str:
    lines = ["| lab | set | class | paired rows (held) | groups (single-row) | methods | s_i median [min, max] |",
             "|---|---|---|---|---|---|---|"]
    for r in tab["rows"]:
        lines.append(f"| {r['lab']} | {r['set']} | {r['class']} | {r['paired_rows']} ({r['held_rows']}) | {r['groups']} "
                     f"({r['single_row_groups']}) | {', '.join(r['methods'])} | {r['s_median']:.3f} "
                     f"[{r['s_min']:.3f}, {r['s_max']:.3f}] |")
    lines += ["", f"Paired rows in all: {tab['paired_rows_total']}.", "",
              "Strict inner masks (rows of L' kept for outer L: EX-056 paired AND both strict nulls answer):", "",
              "| outer L | inner L' | ABS | FIG | FREE | lost to strict nulls |", "|---|---|---|---|---|---|"]
    for x in tab["inner_masks"]:
        lines.append(f"| {x['outer']} | {x['inner']} | {x[ABS]} | {x[FIG]} | {x[FREE]} | {x['lost_to_strict_nulls']} |")
    lines += ["", "Inner laboratories per class (S_c) for each outer laboratory:", ""]
    for outer, d in tab["inner_labs_per_class"].items():
        lines.append(f"- {outer}: " + "; ".join(f"{c} {len(v)} ({', '.join(v) or '-'})" for c, v in d.items()))
    lines += ["", "Leak-toy row counts (EX-056 held-out rows): "
              + ", ".join(f"{k} {v}" for k, v in tab["toy_counts"].items())]
    return "\n".join(lines)


# ----------------------------------------------------------------------------------------------------------------
# 9. Synthetic worlds on the real design (sec. 8; the leak toy, sec. 3.5)
# ----------------------------------------------------------------------------------------------------------------

def draw(design: Design, set_name: str, R: int, rng: np.random.Generator, *, tau: float, omega: float, f: float,
         errors: str = "gauss", noise_scale: float = 1.0, rho: float = 0.0) -> Dict[str, Any]:
    """One synthetic truth for every laboratory's held-out rows: b_g = a_l + u_g, a_l ~ N(0, f tau^2) per class,
    u_g ~ N(0, (1 - f) tau^2); FREE groups an N(0, 1) level; row errors with SD sqrt(noise_scale^2 s_i^2 + omega^2)
    (Gaussian, a scaled t(3) of the same variance, or an AR(1) of correlation rho within a time course)."""
    out: Dict[str, Any] = {"a": {}, "u": {}, "free": {}, "e": {}, "R": R}
    for lab in design.labs:
        lay = design.layouts[lab]
        out["a"][lab] = {cls: math.sqrt(f) * tau * rng.standard_normal(R) for cls in (ABS, FIG)}
        out["u"][lab] = math.sqrt(1.0 - f) * tau * rng.standard_normal((R, lay.G))
        out["free"][lab] = rng.standard_normal((R, lay.G))
        sd = np.sqrt(noise_scale ** 2 * design.s2(lab, set_name) + omega ** 2)
        z = rng.standard_normal((R, lay.n))
        if errors == "t3":
            z = rng.standard_t(3, size=(R, lay.n)) / math.sqrt(3.0)
        elif errors == "ar1" and rho:
            held = design.folds[lab].held
            courses: Dict[Any, List[int]] = defaultdict(list)
            for i, x in enumerate(held):
                courses[(x.table, x.pot, x.species, level_of(x))].append(i)
            for idx in courses.values():
                idx = sorted(idx, key=lambda i: float(held[i].t_min))
                for j in range(1, len(idx)):
                    z[:, idx[j]] = rho * z[:, idx[j - 1]] + math.sqrt(1.0 - rho * rho) * z[:, idx[j]]
        out["e"][lab] = sd[None, :] * z
    return out


def residual_of(design: Design, draws: Dict[str, Any], lab: str, train_labs: Sequence[str], kappa: float) -> np.ndarray:
    """Predicted minus measured for every held-out row of ``lab`` under a fit trained on ``train_labs``: the fit
    absorbs kappa sum_{l in T} pi_l a_l per class (pi_l = l's share of T's rows), the row carries its own level and
    error. kappa = 0: the same residual in every record."""
    lay = design.layouts[lab]
    level = np.empty((draws["R"], lay.n))
    for cls in (ABS, FIG):
        cols = lay.row_class == cls
        if cols.any():
            shift = 0.0
            if kappa:
                tot = float(sum(design.counts[t] for t in train_labs))
                shift = kappa * sum((design.counts[t] / tot) * draws["a"][t][cls] for t in train_labs)
            lv = draws["a"][lab][cls][:, None] + draws["u"][lab][:, lay.gid[cols]]
            level[:, cols] = (np.asarray(shift).reshape(-1, 1) if kappa else 0.0) - lv
    cols = lay.row_class == FREE
    if cols.any():
        level[:, cols] = -draws["free"][lab][:, lay.gid[cols]]
    return level - draws["e"][lab]


def build_world(design: Design, set_name: str, draws: Dict[str, Any], kappa: float = 0.0,
                with_pairs: bool = True) -> Nesting:
    """Synthetic fold and pair records with the real records' structure, read by the real Nesting code path. A
    laboratory's non-paired EX-056 rows are NaN in every record (the engine-finite part of the strict mask)."""
    labs = design.labs
    cache: Dict[Tuple[str, Tuple[str, ...]], np.ndarray] = {}

    def res(lab, train):
        key = (lab, tuple(sorted(train)) if kappa else ())
        if key not in cache:
            r = residual_of(design, draws, lab, train, kappa)
            r[:, ~design.folds[lab].paired] = np.nan
            cache[key] = r
        return cache[key]

    folds = {}
    for lab in labs:
        meta = design.folds[lab]
        r = res(lab, [x for x in labs if x != lab])
        raw = {"engine_warm": r, **({"engine_reset": r} if meta.leak_prone else {})}
        folds[lab] = Fold(lab=lab, held=meta.held, paired=meta.paired, raw=raw, terms=meta.terms,
                          leak_prone=meta.leak_prone)
    pairs: Dict[Tuple[str, str], Pair] = {}
    if with_pairs:
        for a, b in pair_runner.pairs():
            train, held = pair_runner.split(design.data, (a, b))
            tl = [x for x in labs if x not in (a, b)]
            rh = {lab: res(lab, tl) for lab in (a, b)}
            pairs[(a, b)] = Pair(labs=(a, b), held=dict(held), train=train, raw_held={"warm": rh, "reset": rh},
                                 raw_train={}, converged={"warm": True, "reset": True}, leak_prone=["synthetic"])
    return Nesting(design.data, folds, pairs, design.pots_all,
                   pair_terms=lambda pair, v, lab: design.terms(lab, set_name), check_null_masks=False,
                   design_cache=design.nest_cache)


def outer_scores(nest: Nesting, set_name: str, ests: Dict[str, Dict[str, Any]], want_st: bool = False,
                 labs: Sequence[str] = SCORABLE) -> Dict[str, Any]:
    """Arm (b) on each outer laboratory: per-lab LPD per row, IS, 90 % coverage (and every level), width; St
    coverage if asked."""
    per = {}
    for lab in labs:
        blk = nest.outer(lab, "engine", set_name)
        arm = arm_b(blk, ests[lab])
        s = interval_summary(arm["centred"], arm["var"])
        d = {"lpd": arm["lpd"] / blk.layout.n, "is90": s["is90"], "width90": s["width90_median"],
             "cov": {lv: s[f"{lv:.2f}"] / s["scored"] for lv in LEVELS}, "scored": s["scored"],
             "inside90": s["0.90"]}
        if want_st:
            ins = st_inside(arm, df_rows(blk, ests[lab]), levels=(GATE_LEVEL,))
            d["st_cov90"] = ins[GATE_LEVEL][:, np.isfinite(arm["var"][0])].sum(axis=1) / s["scored"]
        per[lab] = d
    lw = np.mean(np.stack([per[lab]["cov"][GATE_LEVEL] for lab in labs]), axis=0)
    pooled = np.sum(np.stack([per[lab]["inside90"] for lab in labs]), axis=0) / sum(per[lab]["scored"] for lab in labs)
    out = {"per_lab": per, "cov90_lab_weighted": lw, "cov90_pooled": pooled}
    if want_st:
        out["st_cov90_lab_weighted"] = np.mean(np.stack([per[lab]["st_cov90"] for lab in labs]), axis=0)
    return out


# ----------------------------------------------------------------------------------------------------------------
# 10. The controls (sec. 8, sec. 3.5, the D2 grid)
# ----------------------------------------------------------------------------------------------------------------

def _mean_se(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    return float(np.mean(x)), float(np.std(x, ddof=1) / math.sqrt(x.size))


def positive_world(design: Design, set_name: str, tau: float, omega: float, f: float, R: int, rng,
                   want_st: bool = True, **draw_kw) -> Dict[str, Any]:
    """One configuration: tau_ABS(all), each outer laboratory's strict nested estimates, (b)'s scores and St's
    coverage, and the UNSTABLE flag."""
    dr = draw(design, set_name, R, rng, tau=tau, omega=omega, f=f, **draw_kw)
    nest = build_world(design, set_name, dr)
    est_all = nest.calibrate_all("engine", set_name)
    ests = {lab: nest.calibrate(lab, "engine", set_name) for lab in SCORABLE}
    sc = outer_scores(nest, set_name, ests, want_st=want_st)
    drops = [nest.calibrate_all("engine", set_name, labs=[x for x in SCORABLE if x != lab])[ABS]["tau"]
             for lab in abs_labs(nest)]
    unstable = unstable_flag(est_all[ABS]["tau"], drops)
    return {"est_all": est_all, "ests": ests, "scores": sc, "unstable": unstable}


def p1_ratios(world: Dict[str, Any], tau: float) -> Dict[str, Optional[np.ndarray]]:
    out = {"all": world["est_all"][ABS]["tau2_raw"] / tau ** 2 if world["est_all"][ABS]["tau2_raw"] is not None
           else None}
    for lab in SCORABLE:
        raw = world["ests"][lab][ABS]["tau2_raw"]
        out[lab] = raw / tau ** 2 if raw is not None else None
    return out


def run_positive(design: Design, R: int) -> Dict[str, Any]:
    configs, checks = [], []
    for set_name in SETS:
        for tau in P_TAUS:
            for omega in P_OMEGAS:
                for f in P_FS:
                    name = f"positive|{set_name}|{tau}|{omega}|{f}"
                    w = positive_world(design, set_name, tau, omega, f, R, rng_for(name))
                    ratios = p1_ratios(w, tau)
                    rec = {"set": set_name, "tau": tau, "omega": omega, "f": f, "R": R, "p1": {}}
                    for est_name, x in ratios.items():
                        if x is None:
                            rec["p1"][est_name] = {"status": "fallback (not checked)"}
                            continue
                        m, se = _mean_se(x)
                        tol = max(0.10, 3.0 * se)
                        rec["p1"][est_name] = {"mean_ratio": m, "mc_se": se, "tolerance": tol,
                                               "pass": abs(m - 1.0) <= tol, "se_units": abs(m - 1.0) / se}
                        if tau in P1_HARD_TAUS:
                            checks.append((set_name, tau, omega, f, est_name, rec["p1"][est_name]))
                    t_all = w["est_all"][ABS]["tau"] / tau
                    rec["p1b"] = {"median": float(np.median(t_all)), "q05": float(np.quantile(t_all, 0.05)),
                                  "q95": float(np.quantile(t_all, 0.95))}
                    rec["p1b_nested"] = {}
                    for lab in SCORABLE:
                        t = w["ests"][lab][ABS]["tau"]
                        if t is not None:
                            rec["p1b_nested"][lab] = {"n_labs": w["ests"][lab][ABS]["n_labs"],
                                                      "median": float(np.median(t / tau)) if tau else None,
                                                      "q05": float(np.quantile(t / tau, 0.05)) if tau else None,
                                                      "q95": float(np.quantile(t / tau, 0.95)) if tau else None}
                    fig = w["est_all"][FIG]["tau2_raw"]
                    rec["fig_all_mean_ratio"] = float(np.mean(fig / tau ** 2)) if fig is not None else None
                    sc = w["scores"]
                    rec["b_cov90_lab_weighted"] = float(np.mean(sc["cov90_lab_weighted"]))
                    rec["b_cov90_pooled"] = float(np.mean(sc["cov90_pooled"]))
                    rec["st_cov90_lab_weighted"] = float(np.mean(sc["st_cov90_lab_weighted"]))
                    rec["shortfall_b"] = GATE_LEVEL - rec["b_cov90_lab_weighted"]
                    rec["shortfall_st"] = GATE_LEVEL - rec["st_cov90_lab_weighted"]
                    rec["b_cov90_per_lab"] = {lab: float(np.mean(sc["per_lab"][lab]["cov"][GATE_LEVEL])) for lab in SCORABLE}
                    rec["p_unstable"] = float(np.mean(w["unstable"]))
                    rec["nested_flags"] = {lab: w["ests"][lab]["flags"] for lab in SCORABLE if w["ests"][lab]["flags"]}
                    configs.append(rec)
                    print(f"  positive {set_name} tau {tau} omega {omega} f {f}: done", flush=True)
    # P1 (hard at tau >= 0.8): one failure is re-run at R = 4000 and must then pass; >= 2, or one by > 5 MC SE, fail
    failures = [c for c in checks if not c[5]["pass"]]
    p1 = {"checks": len(checks), "failures": [{"set": c[0], "tau": c[1], "omega": c[2], "f": c[3], "estimate": c[4],
                                               **c[5]} for c in failures]}
    if not failures:
        p1["pass"] = True
    elif len(failures) >= 2 or failures[0][5]["se_units"] > 5.0:
        p1["pass"] = False
    else:
        s, tau, omega, f, est_name, _ = failures[0]
        w = positive_world(design, s, tau, omega, f, R_RERUN, rng_for(f"p1_rerun|{s}|{tau}|{omega}|{f}"), want_st=False)
        x = p1_ratios(w, tau)[est_name]
        m, se = _mean_se(x)
        p1["rerun"] = {"R": R_RERUN, "mean_ratio": m, "mc_se": se, "pass": abs(m - 1.0) <= max(0.10, 3.0 * se)}
        p1["pass"] = p1["rerun"]["pass"]
    p3 = {s: next(c for c in configs if c["set"] == s and c["tau"] == 0.8 and c["omega"] == 0.3 and c["f"] == 0.5)
          for s in SETS}
    p3_out = {s: {"b_cov90_lab_weighted": c["b_cov90_lab_weighted"], "st_cov90_lab_weighted": c["st_cov90_lab_weighted"],
                  "pass": c["b_cov90_lab_weighted"] >= D1_FLOOR} for s, c in p3.items()}
    p5 = [{"set": c["set"], "tau": c["tau"], "omega": c["omega"], "f": c["f"], "p_unstable": c["p_unstable"]}
          for c in configs if c["tau"] in (0.8, 1.2)]
    gate = max(x["p_unstable"] for x in p5 if x["tau"] == 0.8)
    return {"configs": configs, "P1": p1, "P3": p3_out, "P5": {"table": p5, "gate_p_unstable_tau08": gate,
                                                                "flag_informative": gate <= UNSTABLE_GATE}}


def run_negative(design: Design, R: int) -> Dict[str, Any]:
    out: Dict[str, Any] = {"N1": [], "N2": [], "N3": [], "N4": [], "N5": [], "N5b": []}
    for set_name in SETS:
        for omega in P_OMEGAS:
            for kind in ("gauss", "t3"):
                w = positive_world(design, set_name, 0.0, omega, 0.5, R, rng_for(f"neg|{set_name}|{omega}|{kind}"),
                                   want_st=False, errors=kind)
                frac_all = float(np.mean(w["est_all"][ABS]["tau"] > 0.3))
                if kind == "gauss":
                    out["N1"].append({"set": set_name, "omega": omega, "frac_tau_all_gt_0.3": frac_all,
                                      "pass": frac_all <= 0.10})
                    for lab in SCORABLE:
                        t = w["ests"][lab][ABS]["tau"]
                        frac = float(np.mean(t > 0.3)) if t is not None else None
                        out["N2"].append({"set": set_name, "omega": omega, "outer": lab, "frac_gt_0.3": frac,
                                          "pass": frac is not None and frac <= 0.15})
                else:
                    out["N3"].append({"set": set_name, "omega": omega, "frac_tau_all_gt_0.3": frac_all,
                                      "pass": frac_all <= 0.15})
            w = positive_world(design, set_name, 0.0, omega, 0.5, R, rng_for(f"neg_ar1|{set_name}|{omega}"),
                               want_st=False, errors="ar1", rho=0.5)
            out["N4"].append({"set": set_name, "omega": omega, "rho": 0.5,
                              "frac_tau_all_gt_0.3": float(np.mean(w["est_all"][ABS]["tau"] > 0.3)),
                              "mean_tau_all": float(np.mean(w["est_all"][ABS]["tau"]))})
        for scale in (1.5, 0.67):
            w = positive_world(design, set_name, 0.0, 0.0, 0.5, R, rng_for(f"neg_n5|{set_name}|{scale}"),
                               want_st=False, noise_scale=scale)
            frac = float(np.mean(w["est_all"][ABS]["tau"] > 0.3))
            if scale == 1.5:
                out["N5"].append({"set": set_name, "noise_scale": scale, "frac_tau_all_gt_0.3": frac,
                                  "pass": frac <= 0.20,
                                  "mean_omega_all": float(np.mean(w["est_all"]["omega"]))})
            else:
                sc = w["scores"]
                out["N5b"].append({"set": set_name, "noise_scale": scale, "frac_tau_all_gt_0.3": frac,
                                   "b_cov90_lab_weighted": float(np.mean(sc["cov90_lab_weighted"])),
                                   "b_cov90_pooled": float(np.mean(sc["cov90_pooled"])),
                                   "b_median_width90_per_lab": {lab: float(np.mean(sc["per_lab"][lab]["width90"]))
                                                                for lab in SCORABLE}})
        print(f"  negative {set_name}: done", flush=True)
    for k in ("N1", "N2", "N3", "N5"):
        out[f"{k}_pass"] = all(x["pass"] for x in out[k])
    return out


def run_d3(design: Design, R: int) -> Dict[str, Any]:
    out = {"D3a": [], "D3b": []}
    for set_name in SETS:
        for tau_e, tau_n, which in ((0.8, 0.8, "D3a"), (0.8, 1.2, "D3b"), (0.5, 1.2, "D3b")):
            res = {}
            for model, tau in (("engine", tau_e), ("null", tau_n)):
                w = positive_world(design, set_name, tau, 0.3, 0.5, R,
                                   rng_for(f"d3|{set_name}|{tau_e}|{tau_n}|{model}"), want_st=False)
                res[model] = w["scores"]["per_lab"]
            d = np.stack([res["engine"][lab]["lpd"] - res["null"][lab]["lpd"] for lab in SCORABLE], axis=1)
            isd = np.stack([res["engine"][lab]["is90"] - res["null"][lab]["is90"] for lab in SCORABLE], axis=1)
            dec = d3_decision(d, isd)
            out[which].append({"set": set_name, "engine_tau": tau_e, "null_tau": tau_n, "omega": 0.3, "f": 0.5,
                               "P_lpd_pass": float(np.mean(dec["lpd_pass"])),
                               "P_shown_better": float(np.mean(dec["lpd_pass"] & ~dec["sign_trigger"])),
                               "P_null_better": float(np.mean(dec["null_better"]))})
        print(f"  D3 {set_name}: done", flush=True)
    out["D3a_pass"] = all(x["P_lpd_pass"] <= 0.10 for x in out["D3a"])
    out["sign_flip_threshold"] = SIGN_FLIP_THRESHOLD if out["D3a_pass"] else SIGN_FLIP_THRESHOLD_FALLBACK
    out["power_x"] = min(x["P_lpd_pass"] for x in out["D3b"] if x["engine_tau"] == 0.8)
    return out


def d2_distribution(design: Design, set_name: str, omega: float, f: float, R: int, name: str,
                    grid: Sequence[float] = D2_GRID) -> np.ndarray:
    """The null distribution of the truncated tau_ABS(all) at each grid tau0 (len(grid) x R), on the real design."""
    out = np.empty((len(grid), R))
    for k, tau0 in enumerate(grid):
        dr = draw(design, set_name, R, rng_for(f"{name}|{set_name}|{omega}|{f}|{tau0}"), tau=tau0, omega=omega, f=f)
        nest = build_world(design, set_name, dr, with_pairs=False)
        out[k] = nest.calibrate_all("engine", set_name)[ABS]["tau"]
    return out


def invert(tau_obs: np.ndarray, dist: np.ndarray, grid: Sequence[float] = D2_GRID) -> Tuple[np.ndarray, np.ndarray]:
    """Test inversion (sec. 6): U = the largest tau0 with P(tau_hat <= tau_obs | tau0) >= 0.05 (NaN: none); L = the
    smallest tau0 with P(tau_hat >= tau_obs | tau0) >= 0.05 (NaN: none)."""
    tau_obs = np.atleast_1d(np.asarray(tau_obs, dtype=float))
    g = np.asarray(grid)
    srt = np.sort(dist, axis=1)
    n = srt.shape[1]
    p_le = np.stack([np.searchsorted(srt[k], tau_obs, side="right") / n for k in range(len(g))], axis=1)
    p_ge = np.stack([1.0 - np.searchsorted(srt[k], tau_obs, side="left") / n for k in range(len(g))], axis=1)
    ok_u = p_le >= D2_ALPHA - 1e-15
    ok_l = p_ge >= D2_ALPHA - 1e-15
    U = np.where(ok_u.any(axis=1), g[np.where(ok_u, np.arange(len(g)), -1).max(axis=1)], np.nan)
    L = np.where(ok_l.any(axis=1), g[np.where(ok_l, np.arange(len(g)), len(g)).min(axis=1).clip(max=len(g) - 1)], np.nan)
    return U, L


def run_p2(design: Design, R: int) -> Dict[str, Any]:
    rows = []
    grids = {}
    for set_name in SETS:
        for omega in P_OMEGAS:
            grid1 = d2_distribution(design, set_name, omega, 1.0, D2_REPLICATES, "d2grid")
            grid0 = d2_distribution(design, set_name, omega, 0.0, D2_REPLICATES, "d2grid")
            grids[(set_name, omega)] = {"f1_median": np.median(grid1, axis=1).tolist(),
                                        "f0_median": np.median(grid0, axis=1).tolist()}
            for f in (1.0, 0.0):
                dr = draw(design, set_name, R, rng_for(f"p2|{set_name}|{omega}|{f}"), tau=0.8, omega=omega, f=f)
                t_obs = build_world(design, set_name, dr, with_pairs=False).calibrate_all("engine", set_name)[ABS]["tau"]
                U, L = invert(t_obs, grid1)
                U0, _ = invert(t_obs, grid0)
                covered = np.where(np.isnan(U), False, U >= 0.8 - 1e-12)
                rows.append({"set": set_name, "omega": omega, "truth_f": f, "P_U_ge_tau": float(np.mean(covered)),
                             "P_U_ge_tau_with_f0_grid": float(np.mean(np.where(np.isnan(U0), False, U0 >= 0.8 - 1e-12))),
                             "P_L_le_tau": float(np.mean(np.where(np.isnan(L), True, L <= 0.8 + 1e-12))),
                             "pass": float(np.mean(covered)) >= 0.93})
            print(f"  P2 {set_name} omega {omega}: done", flush=True)
    return {"rows": rows, "pass": all(r["pass"] for r in rows),
            "grid_medians": {f"{k[0]}|omega={k[1]}": v for k, v in grids.items()}}


def _calib_values(est: Dict[str, Any]) -> Dict[str, np.ndarray]:
    out = {"omega2_raw": est["omega2_raw"], "omega": est["omega"]}
    for cls in (ABS, FIG):
        if est[cls]["tau2_raw"] is not None:
            out[f"tau2_raw_{cls}"] = est[cls]["tau2_raw"]
            out[f"tau_{cls}"] = est[cls]["tau"]
    out["tau_fig_used"] = est["tau_fig_used"]
    return out


def _max_rel(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> float:
    if set(a) != set(b):
        return float("inf")
    worst = 0.0
    for k in a:
        x, y = np.asarray(a[k], float), np.asarray(b[k], float)
        den = np.maximum(np.abs(x), np.abs(y))
        rel = np.where(x == y, 0.0, np.abs(x - y) / np.where(den > 0, den, 1.0))
        worst = max(worst, float(np.max(rel)) if rel.size else 0.0)
    return worst


def perturbed(draws: Dict[str, Any], lab: str, cls: str, D: float) -> Dict[str, Any]:
    """The same draw with ONLY lab's own offset of one class moved by D dex (T1)."""
    out = {"a": {k: dict(v) for k, v in draws["a"].items()}, "u": draws["u"], "free": dict(draws["free"]),
           "e": draws["e"], "R": draws["R"]}
    if cls == FREE:
        out["free"][lab] = draws["free"][lab] + D
    else:
        out["a"][lab][cls] = draws["a"][lab][cls] + D
    return out


def run_toy(design: Design, R: int) -> Dict[str, Any]:
    """The leak toy (sec. 3.5): T1 (strict invariance to the outer laboratory's own offset, hard), T2 (reported)."""
    t1_worst, t1_cases, t2 = 0.0, 0, []
    seventh = _seventh_lab(design)
    for set_name in SETS:
        for kappa in TOY_KAPPAS:
            for tau in TOY_TAUS:
                name = f"toy|{set_name}|{kappa}|{tau}"
                rng = rng_for(name)
                dr = draw(design, set_name, R, rng, tau=tau, omega=0.0, f=1.0)
                nest = build_world(design, set_name, dr, kappa=kappa)
                strict = {lab: nest.calibrate(lab, "engine", set_name) for lab in SCORABLE}
                onedeep = {lab: nest.calibrate(lab, "engine", set_name, nesting="onedeep") for lab in SCORABLE}
                rec = {"set": set_name, "kappa": kappa, "tau": tau, "per_lab": {}}
                for lab in SCORABLE:
                    base = _calib_values(strict[lab])
                    base1 = _calib_values(onedeep[lab])
                    slope = {}
                    for cls in (ABS, FIG, FREE):
                        for D in TOY_DS:
                            nest_p = build_world(design, set_name, perturbed(dr, lab, cls, D), kappa=kappa)
                            got = _calib_values(nest_p.calibrate(lab, "engine", set_name))
                            t1_worst = max(t1_worst, _max_rel(base, got))
                            t1_cases += 1
                            if cls != FREE and D == 1.0:
                                got1 = nest_p.calibrate(lab, "engine", set_name, nesting="onedeep")
                                for c2 in (ABS, FIG):
                                    if onedeep[lab][c2]["tau"] is not None and got1[c2]["tau"] is not None:
                                        slope[f"d_tau_{c2}_per_dex_of_D_{cls}"] = float(
                                            np.mean(got1[c2]["tau"] - onedeep[lab][c2]["tau"]))
                    lay = design.layouts[lab]
                    corr = {}
                    for cls in (ABS, FIG):
                        if not np.any((lay.row_class == cls) & design.folds[lab].paired):
                            continue
                        tot = float(sum(design.counts[t] for t in design.labs if t != lab))
                        err = kappa * sum((design.counts[t] / tot) * dr["a"][t][cls] for t in design.labs if t != lab) \
                            - dr["a"][lab][cls]
                        for nm, est in (("strict", strict[lab]), ("onedeep", onedeep[lab])):
                            t = est[cls]["tau"] if est[cls]["tau"] is not None else None
                            corr[f"{cls}_{nm}"] = _corr(np.abs(err), t) if t is not None else None
                    rec["per_lab"][lab] = {"corr_abs_level_error_vs_tau_S": corr, "onedeep_slope": slope}
                sc_s = outer_scores(nest, set_name, strict)
                sc_1 = outer_scores(nest, set_name, onedeep)
                rec["onedeep_minus_strict"] = {
                    "cov90_lab_weighted": float(np.mean(sc_1["cov90_lab_weighted"] - sc_s["cov90_lab_weighted"])),
                    "lpd_lab_weighted": float(np.mean(np.mean(np.stack([sc_1["per_lab"][l]["lpd"] - sc_s["per_lab"][l]["lpd"]
                                                                        for l in SCORABLE]), axis=0)))}
                rec["seventh_lab"] = seventh(nest, dr, kappa, tau, set_name, rng)
                rec["mean_abs_tau_strict_minus_onedeep_ABS"] = {
                    lab: (float(np.mean(np.abs(strict[lab][ABS]["tau"] - onedeep[lab][ABS]["tau"])))
                          if strict[lab][ABS]["tau"] is not None and onedeep[lab][ABS]["tau"] is not None else None)
                    for lab in SCORABLE}
                t2.append(rec)
                print(f"  toy {set_name} kappa {kappa} tau {tau}: done", flush=True)
    return {"T1": {"cases": t1_cases, "max_relative_difference": t1_worst, "pass": t1_worst <= T1_TOL}, "T2": t2}


def _corr(x: np.ndarray, y: np.ndarray) -> Optional[float]:
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _seventh_lab(design: Design):
    """Sec. 3.5's synthetic 7th laboratory: never trained on, ABS rows only, the median scorable laboratory's paired
    row and group counts, predicted by the all-six fit; calibrated by tau(all), omega(all)."""
    n_rows = int(np.median([int(design.folds[lab].paired.sum()) for lab in SCORABLE]))
    n_groups = int(np.median([Layout.of_rows([x for x, m in zip(design.folds[lab].held, design.folds[lab].paired) if m]).G
                              for lab in SCORABLE]))
    abs_s2 = np.concatenate([design.s2(lab, SETS[0])[(design.layouts[lab].row_class == ABS) & design.folds[lab].paired]
                             for lab in SCORABLE])
    s2_med = float(np.median(abs_s2))
    gkeys = [f"seventh|g{i % n_groups}" for i in range(n_rows)]
    lay = Layout(gkeys, ["calibrated"] * n_rows)

    def run(nest, dr, kappa, tau, set_name, rng):
        R = dr["R"]
        delta7 = tau * rng.standard_normal(R)
        tot = float(sum(design.counts.values()))
        err = kappa * sum((design.counts[t] / tot) * dr["a"][t][ABS] for t in design.labs) - delta7
        e = math.sqrt(s2_med) * rng.standard_normal((R, n_rows))
        r = err[:, None] - e
        est_all = nest.calibrate_all("engine", set_name)
        arm = gaussian_arm(lay, np.full(n_rows, s2_med), r, est_all["omega2"], est_all["taus"])
        s = interval_summary(arm["centred"], arm["var"])
        return {"rows": n_rows, "groups": n_groups, "s_median": math.sqrt(s2_med),
                "corr_abs_level_error_vs_tau_all": _corr(np.abs(err), est_all[ABS]["tau"]),
                "b_cov90": float(np.mean(s["0.90"] / s["scored"]))}
    return run


def run_t1b(design: Design, seeds: Sequence[int] = (1, 2)) -> Dict[str, Any]:
    """T1b (hard): synthetic residual FILES with the real record layout, read back by the real loaders and scored by
    the real calibration (training_terms on the pair fits' training residuals, the strict nulls on synthetic data
    values); +3 dex added to every row of lab L in every record (and L's data values / 10^3) must leave outer L's
    tau(S), omega(S), tau_FIG(S) unchanged to 1e-12, for every L, model and set."""
    worst, cases = 0.0, 0
    for seed in seeds:
        rng = rng_for(f"t1b|{seed}")
        kappa = 0.5 if seed == 1 else 1.0
        dr = draw(design, SETS[0], 1, rng, tau=0.8, omega=0.3, f=0.5)
        zvals = {id(x): rng.normal(-1.0, 0.7) for x in design.data}
        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp) / "base"
            base_dir.mkdir()
            data_syn = _write_synthetic_records(design, dr, kappa, zvals, base_dir, rng, perturb_lab=None)
            base = _t1b_calibrations(base_dir, data_syn, design)
            for lab in SCORABLE:
                pdir = Path(tmp) / f"perturbed_{lab}"
                pdir.mkdir()
                data_p = _write_synthetic_records(design, dr, kappa, zvals, pdir, rng, perturb_lab=lab)
                got = _t1b_calibrations(pdir, data_p, design, only_outer=lab)
                for key, vals in got.items():
                    worst = max(worst, _max_rel(base[key], vals))
                    cases += 1
    return {"cases": cases, "max_relative_difference": worst, "pass": worst <= T1_TOL and cases > 0}


def _synthetic_value(row, z: float) -> float:
    try:
        return float(observations.amount_in(row.unit, row.species, 1.0)) * 10.0 ** z
    except Exception:  # noqa: BLE001 - a unit with no molar conversion: any positive value will do
        return 10.0 ** z


def _write_synthetic_records(design: Design, dr: Dict[str, Any], kappa: float, zvals: Dict[int, float], directory: Path,
                             rng, perturb_lab: Optional[str]) -> tuple:
    """Six EX-056-layout fold records and ten pair records, synthetic residuals only; returns the synthetic data rows
    (values replaced). ``perturb_lab``: +3 dex on every residual of that laboratory everywhere, its values / 10^3."""
    D = 3.0
    data_syn = tuple(dataclasses.replace(x, value=_synthetic_value(x, zvals[id(x)] - (D if lab_of(x) == perturb_lab else 0.0)))
                     for x in design.data)
    labs = design.labs

    def res(lab, train):
        r = residual_of(design, dr, lab, train, kappa)[0].copy()
        if lab == perturb_lab:
            r = r + D
        return r

    for lab in labs:
        meta = design.folds[lab]
        r = res(lab, [x for x in labs if x != lab])
        r[~meta.paired] = np.nan
        variants = ["warm"] + (["reset"] if meta.leak_prone else [])
        rows = []
        for i, x in enumerate(meta.held):
            rows.append({**pair_runner.row_key(x), "engine_raw": {v: float(r[i]) for v in variants},
                         "null_raw": {k: float(r[i] + 0.1) for k in NULLS}, "paired": bool(meta.paired[i])})
        rec = {"artifact": f"synthetic_{lab}", "laboratory": lab, "leak_prone": meta.leak_prone, "rows": rows,
               "models": {k: {"terms": t} for k, t in meta.terms.items()}}
        fold_path(lab, directory).write_text(json.dumps(rec), encoding="utf-8")
    for a, b in pair_runner.pairs():
        train, held = pair_runner.split(design.data, (a, b))
        tl = [x for x in labs if x not in (a, b)]
        leak = ["synthetic"] if ("jiangnan_zhang_x" in (a, b) or "dfa_hofmann_schieberle" in (a, b)) else []
        variants = ["warm"] + (["reset"] if leak else [])
        hr = {}
        for lab in (a, b):
            r = res(lab, tl)
            r[~design.folds[lab].paired] = np.nan
            hr[lab] = r
        tr = {}
        for lab in tl:
            tr[lab] = res(lab, tl)
        pos = {lab: 0 for lab in tl}
        train_vals = []
        for x in train:
            lab = lab_of(x)
            train_vals.append(float(tr[lab][pos[lab]]))
            pos[lab] += 1
        rec = {"artifact": f"synthetic_{a}__{b}", "labs": [a, b], "freeze_commit": "synthetic", "leak_prone": leak,
               **{v: {"fit": {"status": 2}, "converged": True} for v in variants},
               "row_keys": {"held": {lab: row_keys(held[lab]) for lab in (a, b)}, "train": row_keys(train)},
               "held_engine_raw": {v: {lab: [float(t) for t in hr[lab] + (0.05 if v == "reset" else 0.0)]
                                       for lab in (a, b)} for v in variants},
               "train_engine_raw": {v: [t + (0.05 if v == "reset" else 0.0) for t in train_vals] for v in variants},
               "complete": True}
        pair_path(a, b, directory).write_text(json.dumps(rec), encoding="utf-8")
    return data_syn


def _t1b_calibrations(directory: Path, data_syn: tuple, design: Design, only_outer: Optional[str] = None
                      ) -> Dict[Tuple[str, str, str], Dict[str, np.ndarray]]:
    folds_rows = make_folds(data_syn)
    folds = {lab: load_fold(lab, held, directory=directory) for lab, (_, held) in folds_rows.items()}
    pairs = {(a, b): load_pair(a, b, data_syn, directory=directory) for a, b in pair_runner.pairs()}
    nest = Nesting(data_syn, folds, pairs, design.pots_all)
    out = {}
    for lab in ([only_outer] if only_outer else SCORABLE):
        for model in MODELS:
            for set_name in SETS:
                out[(lab, model, set_name)] = _calib_values(nest.calibrate(lab, model, set_name))
    return out


def run_controls(R: int = R_CONTROLS) -> Dict[str, Any]:
    t0 = time.time()
    design = Design()
    print("controls: design loaded", flush=True)
    result: Dict[str, Any] = {"artifact": CONTROLS_STEM, "prereg": PREREG.name, "R": R, "seed": SEED,
                              "design": design.table()}
    result["positive"] = run_positive(design, R)
    result["negative"] = run_negative(design, R)
    result["D3"] = run_d3(design, R)
    result["P2"] = run_p2(design, R)
    result["toy"] = run_toy(design, R)
    result["T1b"] = run_t1b(design)
    hard = []
    if not result["positive"]["P1"]["pass"]:
        hard.append("P1")
    if not result["negative"]["N1_pass"]:
        hard.append("N1")
    if not result["toy"]["T1"]["pass"]:
        hard.append("T1")
    if not result["T1b"]["pass"]:
        hard.append("T1b")
    result["summary"] = {
        "hard_failures": hard,
        "P1_pass": result["positive"]["P1"]["pass"], "N1_pass": result["negative"]["N1_pass"],
        "T1_pass": result["toy"]["T1"]["pass"], "T1b_pass": result["T1b"]["pass"],
        "P2_pass": result["P2"]["pass"],
        "P3_pass": {s: v["pass"] for s, v in result["positive"]["P3"].items()},
        "P5_gate_p_unstable_tau08": result["positive"]["P5"]["gate_p_unstable_tau08"],
        "unstable_flag_informative": result["positive"]["P5"]["flag_informative"],
        "N2_pass": result["negative"]["N2_pass"], "N3_pass": result["negative"]["N3_pass"],
        "N5_pass": result["negative"]["N5_pass"],
        "D3a_pass": result["D3"]["D3a_pass"], "sign_flip_threshold": result["D3"]["sign_flip_threshold"],
        "D3b_power_x": result["D3"]["power_x"]}
    result["choices"] = CHOICES
    result["wall_seconds"] = round(time.time() - t0, 1)
    result["provenance"] = _provenance(CONTROLS_STEM, score_inputs=False)
    return result


# ----------------------------------------------------------------------------------------------------------------
# 11. The real scoring (--score, run ONCE)
# ----------------------------------------------------------------------------------------------------------------

def _f(x, fmt="{:.3f}"):
    if x is None:
        return "-"
    try:
        if isinstance(x, float) and not math.isfinite(x):
            return "nan"
    except TypeError:
        pass
    return fmt.format(x)


def _scal(x) -> Optional[float]:
    if x is None:
        return None
    a = np.asarray(x, dtype=float).ravel()
    return float(a[0]) if a.size else None


def _est_summary(est: Dict[str, Any]) -> Dict[str, Any]:
    out = {"omega": _scal(est["omega"]), "omega2_raw": _scal(est["omega2_raw"]), "omega_labs": est["omega_labs"],
           "labs": est["labs"], "flags": est["flags"], "missing_inner": est.get("missing_inner", []),
           "nonconverged_inner": est.get("nonconverged_inner", []),
           "tau_fig_used": _scal(est["tau_fig_used"])}
    for cls in (ABS, FIG):
        rec = est[cls]
        out[cls] = {"tau": _scal(rec["tau"]), "tau2_raw": _scal(rec["tau2_raw"]), "status": rec["status"],
                    "labs": rec["labs"], "n_groups": rec["n_groups"],
                    "per_lab_dbar": {k: _scal(v) for k, v in rec["per_lab_dbar"].items()}}
    out["taus_used"] = {k: _scal(v) for k, v in est["taus"].items()}
    return out


def _cross_check(block: Block, arm: Dict[str, Any], omega2: float, taus: Mapping[str, float]) -> None:
    """The vectorised arm must equal fitting.MeasurementModel (log_density) and fit_thiol_lolo.row_intervals on the
    real rows (sec. 5 'Implementation'), to a relative 1e-9."""
    if not block.rows:
        return
    from fit_thiol_lolo import _model
    noise = {}
    for x, s2 in zip(block.rows, block.s2):
        e2 = s2 - float(x.sigma) ** 2 + omega2
        noise[fj.series_of(x)] = math.sqrt(max(e2, 0.0))
    abs_taus = {m: float(np.asarray(taus[m]).ravel()[0]) for m, t in sm3.SCALE_TAU.items() if float(t) == 0.0}
    pm = _model(noise, abs_taus, level_of)
    pm = dataclasses.replace(pm, scale_tau={**sm3.SCALE_TAU, FIG_METHOD: float(np.asarray(taus[FIG_METHOD]).ravel()[0])})
    r = block.r[0]
    ld = float(pm.log_density(r, block.rows))
    assert _rel(ld, float(arm["lpd"][0])) <= 1e-9 or abs(ld - float(arm["lpd"][0])) < 1e-9, \
        f"{block.lab}: vectorised log density disagrees with MeasurementModel"
    for p, c, v in zip(lolo.row_intervals(pm, r, block.rows), arm["centred"][0], arm["var"][0]):
        if p["var"] is None:
            assert not np.isfinite(v)
        else:
            assert _rel(p["var"], v) <= 1e-9 and abs(p["centred"] - c) <= 1e-9 * max(1.0, abs(c)), \
                f"{block.lab}: vectorised row intervals disagree with row_intervals"


def _lab_scores(block: Block, arm: Dict[str, Any], df_row: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """One laboratory under one arm (R = 1): the Gaussian scores, St and St4."""
    from scipy.stats import t as tdist
    n = block.layout.n
    s = interval_summary(arm["centred"], arm["var"])
    out = {"n": n, "log_density": float(arm["lpd"][0]), "lpd_per_row": float(arm["lpd"][0]) / n if n else None,
           "is90": _scal(s["is90"]), "width90_mean": _scal(s["width90_mean"]),
           "width90_median": _scal(s["width90_median"]),
           "coverage": {f"{lv:.2f}": {"inside": int(s[f"{lv:.2f}"][0]), "scored": s["scored"]} for lv in LEVELS},
           "abs_log_density": float(arm["lpd_g"][0][block.layout.class_g == ABS].sum()) if n else 0.0,
           "abs_n": int(np.sum(block.layout.row_class == ABS))}
    ok = np.isfinite(arm["var"][0])
    # St4: per-row t(4) marginal with the same variance (scale^2 = var / 2)
    c, v = arm["centred"][0][ok], arm["var"][0][ok]
    sc4 = np.sqrt(v / 2.0)
    st4 = {"log_density": float(np.sum(tdist.logpdf(c / sc4, 4) - np.log(sc4)))}
    for lv in LEVELS:
        st4[f"{lv:.2f}"] = int(np.sum(np.abs(c) <= tdist.ppf(0.5 + lv / 2.0, 4) * sc4))
    h4 = tdist.ppf(0.95, 4) * sc4
    st4["is90"] = float(np.mean(2 * h4 + 20 * np.maximum(0, np.abs(c) - h4))) if c.size else None
    st4["scored"] = int(ok.sum())
    out["St4"] = st4
    if df_row is not None:
        lay = block.layout
        st2 = arm["st2"][0]
        scale = arm["scale"][0]
        lpd = 0.0
        for g in range(lay.G):
            idx = np.where(lay.gid == g)[0]
            if lay.class_g[g] == FREE:
                lpd += float(arm["lpd_g"][0][g])
            else:
                lpd += st_group_logdensity(block.r[0][idx], st2[idx], float(scale[idx[0]]), float(df_row[idx[0]]))
        ins = st_inside(arm, df_row, nodes=ST_NODES_REAL)
        h = np.full(n, np.nan)
        for i in range(n):
            if not ok[i]:
                continue
            if np.isfinite(scale[i]):
                h[i] = st_quantile(GATE_LEVEL, math.sqrt(st2[i]), float(scale[i]), float(df_row[i]))
            else:
                h[i] = 1.6448536269514722 * math.sqrt(arm["var"][0][i])
        s_t = interval_summary(arm["centred"], arm["var"], inside_of=ins, half_width90=h[None, :])
        out["St"] = {"log_density": lpd, "lpd_per_row": lpd / n if n else None, "is90": _scal(s_t["is90"]),
                     "width90_median": _scal(s_t["width90_median"]),
                     "coverage": {f"{lv:.2f}": int(s_t[f"{lv:.2f}"][0]) for lv in LEVELS}, "scored": s_t["scored"]}
    return out


def _relative_branch(block: Block, omega2: float) -> Dict[str, Any]:
    """Sec. 5.6: every group with >= 2 rows centred on its 1/s^2-weighted mean, variance s^2 - 1/W + c omega^2; the
    restricted (shape-only) Gaussian density with covariance diag(s^2 + omega^2)."""
    from scipy.stats import norm
    lay = block.layout
    centred, _, c, sel = within_moments(lay, block.s2, block.r)
    if not sel.any():
        return {"n": 0}
    W = (1.0 / block.s2) @ lay.M
    var = block.s2 - 1.0 / W[lay.gid] + c * omega2
    sub_keys = [lay.keys[g] for g in lay.gid[sel]]
    sub = Layout(sub_keys, ["zhou_free"] * int(sel.sum()))   # every kept group free: the restricted density
    arm = gaussian_arm(sub, block.s2[sel], block.r[:, sel], omega2, {})
    out = {"n": int(sel.sum()), "log_density": float(arm["lpd"][0]), "lpd_per_row": float(arm["lpd"][0]) / int(sel.sum())}
    cc, vv = centred[0][sel], var[sel]
    for lv in LEVELS:
        q = float(norm.ppf(0.5 + lv / 2.0))
        out[f"{lv:.2f}"] = int(np.sum(np.abs(cc) <= q * np.sqrt(vv)))
    return out


def _paired_summary(d: np.ndarray, isd: np.ndarray, labs: Sequence[str], threshold: float) -> Dict[str, Any]:
    dec = d3_decision(d[None, :], isd[None, :], threshold)
    lo, hi = bootstrap_interval(d)
    return {"per_lab": dict(zip(labs, map(float, d))), "per_lab_is": dict(zip(labs, map(float, isd))),
            "mean": float(dec["mean"][0]), "se": float(dec["se"][0]), "sign_flip_p_better": float(dec["p_better"][0]),
            "sign_flip_p_worse": float(dec["p_worse"][0]), "bootstrap90": [lo, hi], "is_mean": float(dec["is_mean"][0]),
            "lpd_pass": bool(dec["lpd_pass"][0]), "sign_trigger": bool(dec["sign_trigger"][0]),
            "label": d3_label(dec), "n_patterns": 2 ** d.size}


def _set_verdict(d1_pass: bool, d3: Dict[str, Dict[str, Any]], incomplete: bool, p3_pass: bool, power: float) -> str:
    if incomplete or not p3_pass or any(r["sign_trigger"] for r in d3.values()):
        return "INCONCLUSIVE"
    if not d1_pass:
        return "NOT CALIBRATED"
    if not all(r["lpd_pass"] for r in d3.values()):
        return (f"CALIBRATED, NOT SHOWN BETTER THAN A NULL (D3 power = {power:.2f}, engine and null errors simulated "
                "independent)")
    return "ADEQUATE on EX-056's 5 labs, nested"


PRECEDENCE = ("INVALID", "INCONCLUSIVE", "NOT CALIBRATED", "CALIBRATED, NOT SHOWN BETTER THAN A NULL", "ADEQUATE")


def _rank(v: str) -> int:
    return next(i for i, p in enumerate(PRECEDENCE) if v.startswith(p))


def combine_by_precedence(verdicts: Sequence[str]) -> str:
    """Sec. 6: the verdict that comes FIRST in the precedence (the worse set governs)."""
    return min(verdicts, key=_rank)


def combined_fragility(jack_by_set: Mapping[str, Mapping[str, str]], reported: str) -> Dict[str, Any]:
    """FRAGILE (sec. 6; amendment A1): per dropped laboratory, the two sets' jackknife verdicts combined by precedence;
    FRAGILE iff that combined verdict differs from the reported combined verdict for some dropped laboratory."""
    drops = sorted({d for j in jack_by_set.values() for d in j})
    combined = {d: combine_by_precedence([jack_by_set[s][d] for s in jack_by_set if d in jack_by_set[s]])
                for d in drops}
    changed = sorted(d for d, v in combined.items() if _rank(v) != _rank(reported))
    return {"combined_per_drop": combined, "changes_when_dropped": changed, "fragile": bool(changed)}


def relative_note_fires(branch: str, relative_by_set: Mapping[str, Mapping[str, Mapping[str, Any]]]) -> bool:
    """Sec. 6 'reported, not ruled' (amendment A1): when D2 ships the relative branch (i), the owner is told that it has
    no evidence of beating a null unless it is shown better than BOTH nulls in EVERY set (the worse set governs)."""
    if branch != "(i)":
        return False
    shown = all(rel[null]["lpd_pass"] and not rel[null]["sign_trigger"]
                for rel in relative_by_set.values() for null in NULLS)
    return not shown


def _num(x) -> float:
    """None / NaN / array -> a float (NaN when undefined)."""
    if x is None:
        return float("nan")
    a = np.asarray(x, dtype=float).ravel()
    return float(a[0]) if a.size else float("nan")


def _nanmax(values: Sequence[Any]) -> float:
    """The maximum over the finite values, NaN if there is none: order-free, unlike max() with a NaN."""
    v = [x for x in (_num(y) for y in values) if math.isfinite(x)]
    return max(v) if v else float("nan")


def d2_bound_report(U, L, grid: Sequence[float] = D2_GRID) -> Dict[str, Any]:
    """D2's bounds as reported and used (amendment A1). U NaN means tau_obs lies below the 5 % quantile even at the
    smallest tau0, so U <= grid[0] (used as grid[0]); U at the grid top is 'U >= top' (used as top). L NaN: 'L:
    none'; L at the grid bottom: 'L <= bottom'."""
    lo, hi = float(grid[0]), float(grid[-1])
    u, l_ = _num(U), _num(L)
    if not math.isfinite(u):
        u_used, u_text = lo, f"U <= {lo:.2f}"
    elif u >= hi - 1e-12:
        u_used, u_text = u, f"U >= {hi:.2f}"
    else:
        u_used, u_text = u, f"U = {u:.2f}"
    if not math.isfinite(l_):
        l_text = "L: none"
    elif l_ <= lo + 1e-12:
        l_text = f"L <= {lo:.2f}"
    else:
        l_text = f"L = {l_:.2f}"
    return {"U": None if not math.isfinite(u) else u, "L": None if not math.isfinite(l_) else l_, "U_used": u_used,
            "U_text": u_text, "L_text": l_text}


def d2_decision(per_set: Mapping[str, Mapping[str, Any]], n5_pass: bool, unstable_informative: bool
                ) -> Dict[str, Any]:
    """D2's branch from the per-set tau_ABS(all), U (f = 1, as used), multi-row-group tau and UNSTABLE, every two-set
    maximum NaN-safe and order-free. An undefined tau_ABS(all) in any set refuses (it must not read as small)."""
    taus = [_num(per_set[s]["tau_ABS_all"]) for s in per_set]
    if not all(math.isfinite(t) for t in taus):
        raise SystemExit(f"D2 refuses: tau_ABS(all) is undefined in a set ({dict(zip(per_set, taus))})")
    tau_star = _nanmax(taus)
    U_star = _nanmax([per_set[s]["f1"]["U_used"] for s in per_set])
    multi_star = _nanmax([per_set[s]["multi_row_tau_ABS_all"] for s in per_set])
    unstable = any(bool(per_set[s]["unstable"]) for s in per_set)
    multi_ok = math.isfinite(multi_star) and multi_star > D2_LINE
    if tau_star > D2_LINE and (n5_pass or multi_ok):
        branch = "(i)"
    elif U_star > D2_LINE:
        branch = "(ii)"
    else:
        branch = "(iii)"
    moved = False
    if branch in ("(ii)", "(iii)") and unstable and unstable_informative:
        branch, moved = "(i)", True
    return {"tau_ABS_all_used": tau_star, "U_used": U_star, "multi_row_tau_ABS_all_used": multi_star,
            "unstable_any_set": unstable, "moved_by_unstable": moved, "branch": branch}


def score_set(nest: Nesting, set_name: str, gates: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"set": set_name}
    threshold = gates["sign_flip_threshold"]
    labs_abs = abs_labs(nest)
    models: Dict[str, Any] = {}
    lab_arms: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for model in MODELS:
        est_all = nest.calibrate_all(model, set_name)
        drops = {lab: nest.calibrate_all(model, set_name, labs=[x for x in SCORABLE if x != lab]) for lab in SCORABLE}
        un = bool(unstable_flag(est_all[ABS]["tau"], [drops[lab][ABS]["tau"] for lab in labs_abs])[0])
        multi = nest.calibrate_all(model, set_name, abs_multi_only=True)
        sens = {"paule_mandel_ABS": paule_mandel_tau(est_all, est_all["blocks"]),
                "row_weighted_ABS": row_weighted_tau_blocks(est_all, est_all["blocks"]),
                "per_method": {m: _scal(nest.calibrate_all(model, set_name, abs_methods=(m,))[ABS]["tau"])
                               for m in ("sida", "calibrated")}}
        m_rec: Dict[str, Any] = {"all": _est_summary(est_all), "drop_one_ABS": {k: _scal(v[ABS]["tau"]) for k, v in drops.items()},
                                 "drop_one_omega": {k: _scal(v["omega"]) for k, v in drops.items()},
                                 "drop_one_FIG": {k: _scal(v[FIG]["tau"]) for k, v in drops.items()},
                                 "unstable": un, "tau_ABS_all_multi_row_groups": _scal(multi[ABS]["tau"]),
                                 "sensitivities_all": sens, "nested": {}, "onedeep": {}, "per_lab": {}}
        lab_arms[model] = {}
        for lab in SCORABLE:
            strict = nest.calibrate(lab, model, set_name)
            one = nest.calibrate(lab, model, set_name, nesting="onedeep")
            blk = nest.outer(lab, model, set_name)
            fold = nest.folds[lab]
            key = nest.fold_key(lab, model, set_name)
            taus_a = shipped_taus(fold.terms[key], blk.rows)
            a = arm_a(blk, taus_a)
            b = arm_b(blk, strict)
            b1 = arm_b(blk, one)
            _cross_check(blk, a, 0.0, taus_a)
            _cross_check(blk, b, float(strict["omega2"][0]), {k: _scal(v) for k, v in strict["taus"].items()})
            # arm (a)'s St df: training laboratories holding the class in fold L
            train_labs = {cls: len({lab_of(x) for x in nest.data if lab_of(x) != lab
                                    and class_of_method(fj.method_of(x)) == cls}) for cls in (ABS, FIG)}
            df_a = np.full(blk.layout.n, np.inf)
            for cls in (ABS, FIG):
                df_a[blk.layout.row_class == cls] = float(train_labs[cls])
            lab_arms[model][lab] = {"a": a, "b": b, "b1": b1, "block": blk}
            pm_strict = {"nested": _est_summary(strict), "onedeep": _est_summary(one)}
            sens_l = {"paule_mandel_ABS": paule_mandel_tau(strict, strict["blocks"]),
                      "row_weighted_ABS": row_weighted_tau_blocks(strict, strict["blocks"])}
            m_rec["nested"][lab] = {**pm_strict["nested"], "sensitivities": sens_l}
            m_rec["onedeep"][lab] = pm_strict["onedeep"]
            m_rec["per_lab"][lab] = {"a": _lab_scores(blk, a, df_a), "b": _lab_scores(blk, b, df_rows(blk, strict)),
                                     "b_onedeep": _lab_scores(blk, b1),
                                     "relative": _relative_branch(blk, float(strict["omega2"][0]))}
            leak = {}
            for cls in (ABS, FIG):
                ts, to = _scal(strict[cls]["tau"]), _scal(one[cls]["tau"])
                if ts is not None and to is not None:
                    leak[cls] = {"strict": ts, "onedeep": to, "difference": to - ts,
                                 "measured_leak": abs(to - ts) > ONE_DEEP_LEAK}
            m_rec["per_lab"][lab]["onedeep_vs_strict"] = leak
        models[model] = m_rec
    out["models"] = models
    eng = models["engine"]
    # D1 on (b), engine
    cov_lab = [eng["per_lab"][lab]["b"]["coverage"]["0.90"] for lab in SCORABLE]
    lw = float(np.mean([c["inside"] / c["scored"] for c in cov_lab]))
    pooled = sum(c["inside"] for c in cov_lab) / sum(c["scored"] for c in cov_lab)
    out["D1"] = {"cov90_lab_weighted": lw, "cov90_pooled": pooled, "pass": lw >= D1_FLOOR and pooled >= D1_FLOOR,
                 "acceptance_band": ACCEPTANCE_BAND,
                 "in_acceptance_band": ACCEPTANCE_BAND[0] <= lw <= ACCEPTANCE_BAND[1]}
    # coverage at every level, every arm, every model: lab-weighted and pooled
    cov_tab = {}
    arm_key = {"a": "a", "b": "b", "b_onedeep": "b1"}
    for model in MODELS:
        for arm in ("a", "b", "b_onedeep"):
            row = {}
            for lv in LEVELS:
                cs = [models[model]["per_lab"][lab][arm]["coverage"][f"{lv:.2f}"] for lab in SCORABLE]
                row[f"{lv:.2f}"] = {"lab_weighted": float(np.mean([c["inside"] / c["scored"] for c in cs])),
                                    "pooled": sum(c["inside"] for c in cs) / sum(c["scored"] for c in cs)}
            # widths (amendment A1): the lab-weighted mean, the median of the per-lab medians, the pooled row median
            row["width90_mean_lab_weighted"] = float(np.mean([models[model]["per_lab"][lab][arm]["width90_mean"]
                                                              for lab in SCORABLE]))
            row["width90_median_of_lab_medians"] = float(np.median([models[model]["per_lab"][lab][arm]["width90_median"]
                                                                    for lab in SCORABLE]))
            pooled_w = np.concatenate([2.0 * 1.6448536269514722 * np.sqrt(v[np.isfinite(v)])
                                       for v in (lab_arms[model][lab][arm_key[arm]]["var"][0] for lab in SCORABLE)])
            row["width90_median_pooled_rows"] = float(np.median(pooled_w)) if pooled_w.size else float("nan")
            row["lpd_lab_weighted"] = float(np.mean([models[model]["per_lab"][lab][arm]["lpd_per_row"] for lab in SCORABLE]))
            row["is90_lab_weighted"] = float(np.mean([models[model]["per_lab"][lab][arm]["is90"] for lab in SCORABLE]))
            cov_tab[f"{model}|{arm}"] = row
    out["coverage_table"] = cov_tab
    # D3 and the reported comparisons
    d3, abs_only, relative, st_d3, st4_d3, conc = {}, {}, {}, {}, {}, {}
    for null in NULLS:
        nm = null
        d = np.array([eng["per_lab"][lab]["b"]["lpd_per_row"] - models[nm]["per_lab"][lab]["b"]["lpd_per_row"]
                      for lab in SCORABLE])
        isd = np.array([eng["per_lab"][lab]["b"]["is90"] - models[nm]["per_lab"][lab]["b"]["is90"] for lab in SCORABLE])
        d3[null] = _paired_summary(d, isd, SCORABLE, threshold)
        da = np.array([eng["per_lab"][lab]["b"]["abs_log_density"] / eng["per_lab"][lab]["b"]["abs_n"]
                       - models[nm]["per_lab"][lab]["b"]["abs_log_density"] / models[nm]["per_lab"][lab]["b"]["abs_n"]
                       for lab in labs_abs])
        isa = np.zeros_like(da)
        abs_only[null] = {**_paired_summary(da, isa, labs_abs, threshold),
                          "note": f"{len(labs_abs)} laboratories, {2 ** len(labs_abs)} sign patterns: the minimum p "
                                  f"is 1/{2 ** len(labs_abs)}, so this can never be significant"}
        rl = [lab for lab in SCORABLE if eng["per_lab"][lab]["relative"]["n"] and models[nm]["per_lab"][lab]["relative"]["n"]]
        dr = np.array([eng["per_lab"][lab]["relative"]["lpd_per_row"] - models[nm]["per_lab"][lab]["relative"]["lpd_per_row"]
                       for lab in rl])
        relative[null] = {**_paired_summary(dr, np.zeros_like(dr), rl, threshold),
                          "engine_cov90": float(np.mean([eng["per_lab"][lab]["relative"]["0.90"] / eng["per_lab"][lab]["relative"]["n"] for lab in rl])),
                          "null_cov90": float(np.mean([models[nm]["per_lab"][lab]["relative"]["0.90"] / models[nm]["per_lab"][lab]["relative"]["n"] for lab in rl]))}
        dst = np.array([eng["per_lab"][lab]["b"]["St"]["lpd_per_row"] - models[nm]["per_lab"][lab]["b"]["St"]["lpd_per_row"]
                        for lab in SCORABLE])
        st_d3[null] = _paired_summary(dst, np.zeros_like(dst), SCORABLE, threshold)
        d4 = np.array([(eng["per_lab"][lab]["b"]["St4"]["log_density"] - models[nm]["per_lab"][lab]["b"]["St4"]["log_density"])
                       / eng["per_lab"][lab]["b"]["St4"]["scored"] for lab in SCORABLE])
        st4_d3[null] = _paired_summary(d4, np.zeros_like(d4), SCORABLE, threshold)
        contrib = []
        for lab in SCORABLE:
            me = lab_arms["engine"][lab]["b"]["marg"][0]
            mn = lab_arms[nm][lab]["b"]["marg"][0]
            ok = np.isfinite(me) & np.isfinite(mn)
            contrib.append((me[ok] - mn[ok]) / lab_arms["engine"][lab]["block"].layout.n / len(SCORABLE))
        cvec = np.concatenate(contrib)
        k = max(1, int(math.ceil(CONCENTRATION_TOP * cvec.size)))
        top = np.argsort(-np.abs(cvec))[:k]
        total = float(np.sum(cvec))
        signed = float(np.sum(cvec[top]) / total) if total != 0 else float("inf")
        conc[null] = {"rows": int(cvec.size), "top_rows": k, "signed_share": signed,
                      "abs_share": float(np.sum(np.abs(cvec[top])) / np.sum(np.abs(cvec))) if np.any(cvec) else None,
                      "driven_by_few_rows": (not math.isfinite(signed)) or abs(signed) > CONCENTRATION_SHARE}
    out["D3"] = d3
    out["abs_only"] = abs_only
    out["relative_branch"] = relative
    out["St_paired"] = st_d3
    out["St4_paired"] = st4_d3
    out["concentration"] = conc
    # (b) vs (a), engine: does calibration help
    dba = np.array([eng["per_lab"][lab]["b"]["lpd_per_row"] - eng["per_lab"][lab]["a"]["lpd_per_row"] for lab in SCORABLE])
    isba = np.array([eng["per_lab"][lab]["b"]["is90"] - eng["per_lab"][lab]["a"]["is90"] for lab in SCORABLE])
    out["b_vs_a_engine"] = _paired_summary(dba, isba, SCORABLE, threshold)
    # incomplete strict nesting (sec. 3.1): a NON-CONVERGED pair fit only; a laboratory emptied by the strict masks
    # takes the minimum-data rule, flagged, and caps nothing (amendment A1)
    missing = sorted({n for m in models.values() for lab in SCORABLE for n in m["nested"][lab]["nonconverged_inner"]})
    out["strict_nesting_incomplete"] = missing
    out["strict_masks_emptied"] = sorted({f"outer {lab}, {model}: {n}" for model, m in models.items() for lab in SCORABLE
                                          for n in m["nested"][lab]["missing_inner"]
                                          if n not in m["nested"][lab]["nonconverged_inner"]})
    # the jackknife (FRAGILE)
    power = gates["D3b_power_x"]
    verdict = _set_verdict(out["D1"]["pass"], d3, bool(missing), gates["P3_pass"][set_name], power)
    out["verdict"] = verdict
    jack = {}
    for drop in SCORABLE:
        keep = [lab for lab in SCORABLE if lab != drop]
        cs = [eng["per_lab"][lab]["b"]["coverage"]["0.90"] for lab in keep]
        lwj = float(np.mean([c["inside"] / c["scored"] for c in cs]))
        pj = sum(c["inside"] for c in cs) / sum(c["scored"] for c in cs)
        d3j = {}
        for null in NULLS:
            d = np.array([d3[null]["per_lab"][lab] for lab in keep])
            isd = np.array([d3[null]["per_lab_is"][lab] for lab in keep])
            dec = d3_decision(d[None, :], isd[None, :], threshold)
            d3j[null] = {"lpd_pass": bool(dec["lpd_pass"][0]), "sign_trigger": bool(dec["sign_trigger"][0])}
        jack[drop] = _set_verdict(lwj >= D1_FLOOR and pj >= D1_FLOOR, d3j, bool(missing), gates["P3_pass"][set_name], power)
    out["jackknife"] = jack
    # reported only: FRAGILE itself is read on the two sets combined by precedence (run_score; amendment A1)
    out["fragile_within_set"] = any(_rank(v) != _rank(verdict) for v in jack.values())
    # D-all (diagnostic, engine only, every finite held-out row, six laboratories)
    dall = {}
    est_all = nest.calibrate_all("engine", set_name)
    for lab in sorted(nest.folds):
        blk = nest.outer(lab, "engine", set_name, every_finite=True)
        key = nest.fold_key(lab, "engine", set_name)
        a = arm_a(blk, shipped_taus(nest.folds[lab].terms[key], blk.rows))
        est = nest.calibrate(lab, "engine", set_name) if lab in SCORABLE else est_all
        b = arm_b(blk, est)
        dall[lab] = {}
        for nm_arm, arm in (("a", a), ("b", b)):
            s = interval_summary(arm["centred"], arm["var"])
            dall[lab][nm_arm] = {"rows": blk.layout.n, "scored": s["scored"],
                                 **{f"cov{lv:.2f}": (float(s[f"{lv:.2f}"][0]) / s["scored"]) if s["scored"] else None
                                    for lv in LEVELS},
                                 "width90_median": _scal(s["width90_median"])}
    out["D_all"] = dall
    return out


def _git(*args) -> Optional[str]:
    p = subprocess.run(["git", *args], cwd=str(ROOT), capture_output=True, text=True)
    return p.stdout.strip() if p.returncode == 0 else None


#: the files whose change since the freeze makes --score an amendment (sec. 7.5); the engine files were added in A1
FROZEN_FILES = (GENERATED_BY, "scripts/generators/fit_thiol_lolo_pairs.py", "scripts/generators/fit_thiol_lolo.py",
                "scripts/generators/fit_joint_deoxyosone_thiol.py", "scripts/generators/fit_deoxyosone_submodel.py",
                "scripts/generators/build_submodel_probe.py", "src/kinetic_core/fitting.py")


def _freeze_check(commit: str) -> Dict[str, Any]:
    full = _git("rev-parse", "--verify", f"{commit}^{{commit}}")
    files = list(FROZEN_FILES) + [data_paths.rel(OUT_DIR / f"{CONTROLS_STEM}.json"), data_paths.rel(PREREG)]
    out = {"freeze_commit": commit, "found": full is not None, "files_checked": files}
    if full:
        changed = (_git("diff", "--name-only", full, "--", *files) or "").split()
        out["changed_since_freeze"] = changed
        out["ancestor_of_head"] = subprocess.run(["git", "merge-base", "--is-ancestor", full, "HEAD"],
                                                 cwd=str(ROOT)).returncode == 0
    return out


# ---- the reproduction guards inside --score (amendment A1) ----------------------------------------------------

#: the single-out fold re-run through the pair runner for engine guard 4.2(c)
GUARD_RECORD_LAB = pair_runner.GUARD_FOLD


def guard_record_path(directory: Optional[Path] = None) -> Path:
    return Path(directory or PAIR_DIR) / f"{pair_runner.STEM}_guard_{GUARD_RECORD_LAB}.json"


def _terms_rel(new: Any, saved: Any, path: str = "") -> Tuple[float, List[str]]:
    """The largest relative difference between two terms structures (recursive; both-NaN and both-None count as
    equal); a structural mismatch is returned as a problem."""
    if isinstance(new, Mapping) or isinstance(saved, Mapping):
        if not (isinstance(new, Mapping) and isinstance(saved, Mapping)):
            return float("inf"), [f"{path or '/'}: one side is not a mapping"]
        if set(new) != set(saved):
            return float("inf"), [f"{path or '/'}: keys differ {sorted(set(new) ^ set(saved))}"]
        worst, probs = 0.0, []
        for k in sorted(saved):
            w, p = _terms_rel(new[k], saved[k], f"{path}/{k}")
            worst, probs = max(worst, w), probs + p
        return worst, probs
    if isinstance(new, (list, tuple)) or isinstance(saved, (list, tuple)):
        if not (isinstance(new, (list, tuple)) and isinstance(saved, (list, tuple))) or len(new) != len(saved):
            return float("inf"), [f"{path}: list lengths differ"]
        worst, probs = 0.0, []
        for i, (a, b) in enumerate(zip(new, saved)):
            w, p = _terms_rel(a, b, f"{path}[{i}]")
            worst, probs = max(worst, w), probs + p
        return worst, probs
    if new is None or saved is None:
        return (0.0, []) if new is None and saved is None else (float("inf"), [f"{path}: None on one side"])
    a, b = float(new), float(saved)
    if math.isnan(a) and math.isnan(b):
        return 0.0, []
    return _rel(a, b), []


def reproduction_guards(data: Sequence[Any], folds_rows: Mapping[str, Tuple[tuple, tuple]], folds: Mapping[str, Fold],
                        pots_all, guard_rec: Mapping[str, Any], null_fn: Optional[Callable] = None,
                        terms_fn: Optional[Callable] = None) -> Dict[str, Any]:
    """Amendment A1, before any nested number. (i) training_terms(null_residuals(kind, train, train), train, n)
    reproduces each saved EX-056 null `terms` (6 folds x 2 nulls) to a relative 1e-9; (ii) training_terms on the guard
    record's train_engine_raw['warm'] (rows aligned to split(data, (Reading,))'s training rows) reproduces Reading's
    saved models.engine_warm.terms to 1e-9; (iii) engine guard 4.2(c) recomputed from the guard record against EX-056's
    Reading record (cost to a relative 1e-3; held-out raw to 0.01 dex, same finite mask). Returns the record; its
    'pass' is False on any failure (the caller refuses: INVALID)."""
    null_fn = null_fn or (lambda kind, train, rows: null_residuals(kind, train, rows, pots_all))
    terms_fn = terms_fn or training_terms
    checks: List[Dict[str, Any]] = []
    for lab, (train, _) in sorted(folds_rows.items()):
        for kind in NULLS:
            raw, n = null_fn(kind, train, train)
            rel, probs = _terms_rel(terms_fn(np.asarray(raw, dtype=float), train, n), folds[lab].terms[f"null_{kind}"])
            checks.append({"fold": lab, "model": f"null_{kind}", "max_relative_difference": rel, "problems": probs,
                           "pass": not probs and rel <= GUARD_TOL})
    lab = GUARD_RECORD_LAB
    problems: List[str] = []
    train, held = pair_runner.split(data, (lab,))
    if list(guard_rec.get("labs") or []) != [lab]:
        problems.append(f"guard record labs {guard_rec.get('labs')} != {[lab]}")
    if guard_rec.get("complete") is not True:
        problems.append("guard record not complete")
    keys = guard_rec.get("row_keys") or {}
    if keys.get("train") != row_keys(train):
        problems.append("guard record training row keys differ from split(data, (Reading,))'s training rows")
    if (keys.get("held") or {}).get(lab) != row_keys(held[lab]):
        problems.append("guard record held-out row keys differ from split(data, (Reading,))'s held-out rows")
    engine_terms = {"fold": lab, "model": "engine_warm (guard record train_engine_raw['warm'])"}
    cost = {"tolerance": float(lolo.REPRODUCTION_TOL)}
    raw_c = {"tolerance_dex": float(pair_runner.GUARD_RAW_TOL_C)}
    if not problems:
        rel, probs = _terms_rel(terms_fn(np.asarray(guard_rec["train_engine_raw"]["warm"], dtype=float), train,
                                         N_CONSTANTS_ENGINE), folds[lab].terms["engine_warm"])
        engine_terms.update({"max_relative_difference": rel, "problems": probs, "pass": not probs and rel <= GUARD_TOL})
        c_new, c_old = float(guard_rec["warm"]["fit"]["cost"]), float(folds[lab].record["warm"]["fit"]["cost"])
        cost.update({"guard_record": c_new, "ex056": c_old, "relative_difference": abs(c_new - c_old) / abs(c_old)})
        cost["pass"] = cost["relative_difference"] <= cost["tolerance"]
        new = np.asarray(guard_rec["held_engine_raw"]["warm"][lab], dtype=float)
        old = np.asarray(folds[lab].raw["engine_warm"][0], dtype=float)
        same = new.shape == old.shape and bool(np.array_equal(np.isfinite(new), np.isfinite(old)))
        diff = float(np.nanmax(np.abs(new - old))) if same and np.isfinite(new).any() else (0.0 if same else float("inf"))
        raw_c.update({"same_finite_mask": same, "max_abs_difference_dex": diff,
                      "pass": same and diff <= raw_c["tolerance_dex"]})
    else:
        engine_terms["pass"] = cost["pass"] = raw_c["pass"] = False
    checks.append(engine_terms)
    terms_ok = all(c["pass"] for c in checks)
    out = {"terms": {"checks": checks, "n_checks": len(checks), "tolerance": GUARD_TOL, "pass": terms_ok,
                     "max_relative_difference": max(c.get("max_relative_difference", float("inf")) for c in checks)},
           "engine_guard_4_2c": {"fold": lab, "cost": cost, "held_out_raw": raw_c,
                                 "pass": bool(cost["pass"] and raw_c["pass"])},
           "guard_record_problems": problems}
    out["pass"] = bool(terms_ok and out["engine_guard_4_2c"]["pass"] and not problems)
    return out


def _failed_guard_checks(rg: Mapping[str, Any]) -> List[str]:
    failed = [f"{c['fold']} {c['model']}" for c in rg["terms"]["checks"] if not c["pass"]]
    failed += [f"4.2(c) {k}" for k in ("cost", "held_out_raw") if not rg["engine_guard_4_2c"][k]["pass"]]
    return failed + list(rg["guard_record_problems"])


# ---- run once (sec. 7.5; amendment A1) ---------------------------------------------------------------------------

def score_output(amendment: str = "", directory: Optional[Path] = None) -> Dict[str, Any]:
    """Where --score writes, and under what label. Refuses if the primary record exists and no --amendment REASON is
    given; with one, the next free `..._amendment_<k>`. Without a primary record, --amendment writes the primary
    record with its reason (a pre-scoring amendment)."""
    d = Path(directory or OUT_DIR)
    reason = (amendment or "").strip()
    primary = d / f"{STEM}.json"
    if not primary.exists():
        return {"stem": STEM, "kind": "pre-scoring amendment" if reason else "primary", "amendment": reason}
    if not reason:
        raise SystemExit(f"--score refuses to run: {primary.name} exists (prereg sec. 7.5: the scoring step runs "
                         "once); a later run is an amendment: --amendment REASON, written to its own record")
    k = 1
    while (d / f"{STEM}_amendment_{k}.json").exists() or (d / f"{STEM}_amendment_{k}.md").exists():
        k += 1
    return {"stem": f"{STEM}_amendment_{k}", "kind": f"post-scoring amendment {k}", "amendment": reason, "k": k}


def run_score(amendment: str = "", output: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    output = output or score_output(amendment)
    recs = check_pair_records()
    ctrl_path = OUT_DIR / f"{CONTROLS_STEM}.json"
    if not ctrl_path.exists():
        raise SystemExit(f"--score refuses to run: the controls record {ctrl_path.name} does not exist (sec. 7.2)")
    controls = json.loads(ctrl_path.read_text(encoding="utf-8"))
    gates = controls["summary"]
    if gates["hard_failures"]:
        raise SystemExit(f"--score refuses to run: INVALID, hard control failures {gates['hard_failures']} (sec. 6): "
                         "repair the estimator or the code by amendment before any real number is read")
    g_path = guard_record_path()
    if not g_path.exists():
        raise SystemExit(f"--score refuses to run: the engine-guard record {g_path.name} (sec. 4.2(c)) does not exist")
    guard_rec = json.loads(g_path.read_text(encoding="utf-8"))
    t0 = time.time()
    data, _, pots_all = load_fit_table()
    folds_rows = make_folds(data)
    ok, worst, failed = guard(data, folds_rows)
    if not ok:
        raise SystemExit(f"--score stops: INVALID, the reproduction guard (sec. 4.1) failed: {failed}")
    commits = {r["freeze_commit"] for r in recs.values()}
    if len(commits) != 1:
        raise SystemExit(f"--score refuses to run: the pair records carry different freeze commits {sorted(commits)}")
    freeze = _freeze_check(commits.pop())
    if not freeze["found"] or not freeze.get("ancestor_of_head"):
        raise SystemExit(f"--score refuses to run: the freeze commit is missing or not an ancestor of HEAD: {freeze}")
    if freeze.get("changed_since_freeze") and not output["amendment"]:
        raise SystemExit(f"--score refuses to run: frozen files changed since the freeze {freeze['changed_since_freeze']}; "
                         "rerun with --amendment REASON (prereg sec. 7.5: an amendment, scored separately)")
    freeze["amendment"] = output["amendment"]
    folds = {lab: load_fold(lab, held) for lab, (_, held) in folds_rows.items()}
    # amendment A1: the functions the scorer recomputes fresh must reproduce EX-056, before any nested number
    repro = reproduction_guards(data, folds_rows, folds, pots_all, guard_rec)
    if not repro["pass"]:
        raise SystemExit("--score stops: INVALID, a reproduction guard inside --score failed (training_terms or "
                         f"engine guard 4.2(c)): {_failed_guard_checks(repro)}")
    pairs = {(a, b): load_pair(a, b, data, rec=recs[(a, b)]) for a, b in pair_runner.pairs()}
    nest = Nesting(data, folds, pairs, pots_all)
    result: Dict[str, Any] = {"artifact": output["stem"], "prereg": PREREG.name, "run": output,
                              "guard": {"pass": ok, "max_relative_difference": worst},
                              "reproduction_guards": repro, "freeze": freeze, "controls_gates": gates,
                              "runner_unit_tests": (f"tests/unit/test_thiol_lolo_pairs.py passed at the freeze commit "
                                                    f"{freeze['freeze_commit']} (prereg secs. 4.3, 7.2); the scorer "
                                                    "cannot re-run them"),
                              "pair_fits": {f"{a}__{b}": {v: p.converged[v] for v in p.converged}
                                            for (a, b), p in pairs.items()},
                              "sets": {}}
    for set_name in SETS:
        result["sets"][set_name] = score_set(nest, set_name, gates)
        print(f"score: set {set_name} done", flush=True)
    # D2 (the ship rule; not nested), with its test-inversion bounds at each set's omega_hat(all)
    design = Design()
    d2: Dict[str, Any] = {"per_set": {}}
    for set_name in SETS:
        eng = result["sets"][set_name]["models"]["engine"]
        t_obs, om = _num(eng["all"][ABS]["tau"]), _num(eng["all"]["omega"])
        if not (math.isfinite(t_obs) and math.isfinite(om)):
            raise SystemExit(f"D2 refuses: tau_ABS(all) or omega(all) undefined in set {set_name} ({t_obs}, {om})")
        rec = {"tau_ABS_all": t_obs, "omega_all": om}
        for f in (1.0, 0.0):
            dist = d2_distribution(design, set_name, om, f, D2_REPLICATES, "d2score")
            U, L = invert(np.array([t_obs]), dist)
            rec[f"f{int(f)}"] = d2_bound_report(U, L)
        rec["unstable"] = eng["unstable"]
        rec["multi_row_tau_ABS_all"] = _num(eng["tau_ABS_all_multi_row_groups"])
        d2["per_set"][set_name] = rec
    informative = bool(gates["unstable_flag_informative"])
    dec = d2_decision(d2["per_set"], bool(gates["N5_pass"]), informative)
    branch = dec["branch"]
    # the reported U is the larger set's; on a tie a defined U is preferred (order-free)
    best = sorted(SETS, key=lambda s: (-d2["per_set"][s]["f1"]["U_used"], d2["per_set"][s]["f1"]["U"] is None, s))[0]
    d2.update({**dec, "U_text": d2["per_set"][best]["f1"]["U_text"], "unstable_informative": informative,
               "N5_fallback_applied": not gates["N5_pass"],
               "bounds_label": "" if gates["P2_pass"] else "anti-conservative (P2 failed); reported from f = 1",
               "unstable_note": ("" if informative else
                                 f"UNSTABLE fires under a homogeneous truth with P = {gates['P5_gate_p_unstable_tau08']:.2f} "
                                 "and changes nothing")})
    d2["text"] = {"(i)": "the thiol lane ships relative within-lab predictions only (within-study ratios, ranks, "
                         "time-course shapes, with omega in their intervals); absolute concentrations for a new "
                         "laboratory are refused",
                  "(ii)": "absolute predictions ship only with tau_ABS-widened intervals, labelled 'provisional; tau "
                          "biased low'",
                  "(iii)": "absolute predictions ship with tau_ABS-widened intervals, labelled provisional"}[branch]
    result["D2"] = d2
    sv = {s: result["sets"][s]["verdict"] for s in SETS}
    final = combine_by_precedence(list(sv.values()))
    frag = combined_fragility({s: result["sets"][s]["jackknife"] for s in SETS}, final)
    result["verdict"] = {"per_set": sv, "reported": final, "sets_disagree": len({_rank(v) for v in sv.values()}) > 1,
                         "fragile": frag["fragile"], "fragility": frag,
                         "fragile_within_set": {s: result["sets"][s]["fragile_within_set"] for s in SETS},
                         "provisional": True}
    notes = []
    om_all = _nanmax([result["sets"][s]["models"]["engine"]["all"]["omega"] for s in SETS])
    if math.isfinite(om_all) and om_all > OMEGA_NOTICE:
        notes.append(f"omega(all) = {om_all:.2f} dex > 0.5: relative predictions also carry a large shape error")
    if relative_note_fires(branch, {s: result["sets"][s]["relative_branch"] for s in SETS}):
        notes.append("D2 ships the relative branch, and the relative branch is not shown better than both nulls in "
                     "every set: it has no evidence of beating a null either")
    result["owner_notes"] = notes
    flags = []
    if not gates["N2_pass"]:
        flags.append("N2 failed: each outer laboratory's nested tau_ABS(S) is flagged (it exceeds 0.3 dex too often "
                     "under tau = 0)")
    if not gates["N3_pass"]:
        flags.append("N3 failed: the real tau_ABS is flagged 'may be inflated by outliers'")
    if not gates["N5_pass"]:
        flags.append("N5 failed: D2's multi-row-group fallback applies (branch (i) also needs the multi-row-group "
                     "tau_ABS(all) > 0.8)")
    result["control_flags"] = {"N2_pass": bool(gates["N2_pass"]), "N3_pass": bool(gates["N3_pass"]),
                               "N5_pass": bool(gates["N5_pass"]), "flags": flags}
    result["choices"] = CHOICES
    result["wall_seconds"] = round(time.time() - t0, 1)
    result["provenance"] = _provenance(output["stem"], score_inputs=True)
    return result


# ----------------------------------------------------------------------------------------------------------------
# 12. Rendering and provenance
# ----------------------------------------------------------------------------------------------------------------

def _provenance(name: str, score_inputs: bool) -> Dict[str, Any]:
    tables = fj.TABLES + (fj.ACTZ_TABLE, fj.ENTRY_TABLE, fj.ENTRY_PH_TABLE)
    inputs = ([fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in tables] + [fit_rows.LABORATORIES]
              + [fold_path(lab) for lab in sorted(TOY_COUNTS)]
              + [ROOT / GENERATED_BY, ROOT / "scripts/generators/fit_thiol_lolo.py",
                 ROOT / "scripts/generators/fit_thiol_lolo_pairs.py", ROOT / "src/kinetic_core/fitting.py", PREREG, REVIEW])
    if score_inputs:
        inputs += ([OUT_DIR / f"{CONTROLS_STEM}.json"] + [pair_path(a, b) for a, b in pair_runner.pairs()]
                   + [guard_record_path()])
    return provenance.provenance_block(name, generated_by=GENERATED_BY, inputs=inputs)


def _write(stem: str, payload: Dict[str, Any], md: str) -> None:
    (OUT_DIR / f"{stem}.json").write_text(provenance.dumps_compact(payload, default=float) + "\n", encoding="utf-8")
    (OUT_DIR / f"{stem}.md").write_text(md + "\n", encoding="utf-8")


def _pf(x) -> str:
    return "PASS" if x else "FAIL"


def render_controls(res: Dict[str, Any]) -> str:
    s = res["summary"]
    pos, neg, d3, p2, toy = res["positive"], res["negative"], res["D3"], res["P2"], res["toy"]
    L = ["# H0 controls: synthetic residuals on the real design (prereg secs. 8, 3.5; the D2 grid)", "",
         f"*Generated by `{GENERATED_BY} --controls`; pre-registration `{PREREG.name}` (v3). Every residual here is "
         f"synthetic; no real residual was read. R = {res['R']} per configuration (P1 re-run at {R_RERUN}), D2 grids "
         f"{D2_REPLICATES} per tau0, seed {SEED}. Design noise: row sigma plus EX-056's saved engine class noise per "
         "fold (warm terms; reset terms where a fold has a reset arm), a published aggregate.*", "",
         f"**Hard criteria: P1 {_pf(s['P1_pass'])}, N1 {_pf(s['N1_pass'])}, T1 {_pf(s['T1_pass'])}, "
         f"T1b {_pf(s['T1b_pass'])}** -> " + ("no hard failure" if not s["hard_failures"]
                                              else f"INVALID ({', '.join(s['hard_failures'])})"), "",
         f"- P3 (D1 interpretable): " + ", ".join(f"{k} {_pf(v)}" for k, v in s["P3_pass"].items()),
         f"- P2 (test-inversion U >= tau in >= 93 %): {_pf(s['P2_pass'])}",
         f"- P5: largest P(UNSTABLE | tau = 0.8) = {s['P5_gate_p_unstable_tau08']:.3f} -> the flag is "
         + ("informative (<= 0.3): it can move D2 to branch (i)" if s["unstable_flag_informative"]
            else "uninformative (> 0.3): reported, changes nothing"),
         f"- N2 {_pf(s['N2_pass'])} (flag only), N3 {_pf(s['N3_pass'])} (flag only), N5 {_pf(s['N5_pass'])} "
         "(failure: D2's multi-row-group fallback)",
         f"- D3a (size <= 0.10): {_pf(s['D3a_pass'])}; sign-flip threshold {s['sign_flip_threshold'] * 32:.0f}/32; "
         f"D3b power x = {s['D3b_power_x']:.3f}", "",
         "## Design table (sec. 7.1)", "", render_design(res["design"]), "",
         "## P1: mean untruncated tau_hat^2 / tau^2 (hard at tau >= 0.8; tolerance max(10 %, 3 MC SE))", "",
         f"{pos['P1']['checks']} hard checks; failures: {len(pos['P1']['failures'])}"
         + (f"; re-run: {pos['P1']['rerun']}" if "rerun" in pos["P1"] else "") + f" -> {_pf(pos['P1']['pass'])}", "",
         "| set | tau | omega | f | all | " + " | ".join(SCORABLE) + " | P1b median [5 %, 95 %] of tau_hat/tau | FIG(all) |",
         "|---|---|---|---|---|" + "---|" * len(SCORABLE) + "---|---|"]
    for c in pos["configs"]:
        cells = []
        for k in ("all",) + SCORABLE:
            x = c["p1"][k]
            cells.append(f"{x['mean_ratio']:.3f} +- {x['mc_se']:.3f}{'' if x['pass'] else ' (X)'}" if "mean_ratio" in x
                         else x["status"])
        L.append(f"| {c['set']} | {c['tau']} | {c['omega']} | {c['f']} | " + " | ".join(cells)
                 + f" | {c['p1b']['median']:.2f} [{c['p1b']['q05']:.2f}, {c['p1b']['q95']:.2f}] | "
                 f"{_f(c['fig_all_mean_ratio'])} |")
    L += ["", "chi^2_2/2 benchmark for tau_hat/tau at f = 1 (one lab-level truth, 2 labs): median 0.83, 5 - 95 % 0.23 - 1.73.",
          "", "P1b for the nested estimates at f = 1 (median [5 %, 95 %] of tau_hat(S)/tau; n inner ABS labs):", ""]
    for c in pos["configs"]:
        if c["f"] == 1.0 and c["tau"] in P1_HARD_TAUS:
            L.append(f"- {c['set']} tau {c['tau']} omega {c['omega']}: " + "; ".join(
                f"{lab.split('_')[0]} ({x['n_labs']}) {x['median']:.2f} [{x['q05']:.2f}, {x['q95']:.2f}]"
                for lab, x in c["p1b_nested"].items()))
    L += [
          "", "## P3 / P4: (b)'s lab-weighted 90 % coverage and St's (P3 at tau 0.8, omega 0.3, f 0.5; floor 0.70)", "",
          "| set | tau | omega | f | (b) lab-weighted | (b) pooled | St lab-weighted | shortfall (b) | shortfall St | "
          "P(UNSTABLE) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for c in pos["configs"]:
        L.append(f"| {c['set']} | {c['tau']} | {c['omega']} | {c['f']} | {c['b_cov90_lab_weighted']:.3f} | "
                 f"{c['b_cov90_pooled']:.3f} | {c['st_cov90_lab_weighted']:.3f} | {c['shortfall_b']:+.3f} | "
                 f"{c['shortfall_st']:+.3f} | {c['p_unstable']:.3f} |")
    L += ["", "P3: " + "; ".join(f"{k}: (b) {v['b_cov90_lab_weighted']:.3f}, St {v['st_cov90_lab_weighted']:.3f} -> "
                                 f"{_pf(v['pass'])}" for k, v in pos["P3"].items()), "",
          "## Negative controls (tau = 0)", "", "| id | set | omega / scale | P(tau > 0.3) | criterion | result |",
          "|---|---|---|---|---|---|"]
    for x in neg["N1"]:
        L.append(f"| N1 | {x['set']} | {x['omega']} | {x['frac_tau_all_gt_0.3']:.3f} | <= 0.10 | {_pf(x['pass'])} |")
    for x in neg["N2"]:
        L.append(f"| N2 {x['outer']} | {x['set']} | {x['omega']} | {_f(x['frac_gt_0.3'])} | <= 0.15 | {_pf(x['pass'])} |")
    for x in neg["N3"]:
        L.append(f"| N3 (t3) | {x['set']} | {x['omega']} | {x['frac_tau_all_gt_0.3']:.3f} | <= 0.15 | {_pf(x['pass'])} |")
    for x in neg["N4"]:
        L.append(f"| N4 (AR1 0.5) | {x['set']} | {x['omega']} | {x['frac_tau_all_gt_0.3']:.3f} (mean tau "
                 f"{x['mean_tau_all']:.3f}) | reported | - |")
    for x in neg["N5"]:
        L.append(f"| N5 (1.5 x s) | {x['set']} | {x['noise_scale']} | {x['frac_tau_all_gt_0.3']:.3f} (mean omega "
                 f"{x['mean_omega_all']:.3f}) | <= 0.20 | {_pf(x['pass'])} |")
    for x in neg["N5b"]:
        L.append(f"| N5b (0.67 x s) | {x['set']} | {x['noise_scale']} | {x['frac_tau_all_gt_0.3']:.3f}; (b) cov90 "
                 f"{x['b_cov90_lab_weighted']:.3f} lab-weighted / {x['b_cov90_pooled']:.3f} pooled; median widths "
                 + ", ".join(f"{k.split('_')[0]} {v:.2f}" for k, v in x["b_median_width90_per_lab"].items())
                 + " | reported | - |")
    L += ["", "## D3 size and power", "", "| id | set | engine tau | null tau | P(LPD criterion) | P(shown better) | "
          "P(null better) |", "|---|---|---|---|---|---|---|"]
    for k in ("D3a", "D3b"):
        for x in d3[k]:
            L.append(f"| {k} | {x['set']} | {x['engine_tau']} | {x['null_tau']} | {x['P_lpd_pass']:.3f} | "
                     f"{x['P_shown_better']:.3f} | {x['P_null_better']:.3f} |")
    L += ["", "## P2: test-inversion upper bound at tau = 0.8 (grid f = 1 at the true omega)", "",
          "| set | omega | truth f | P(U >= 0.8) | with the f = 0 grid | P(L <= 0.8) | >= 0.93 |", "|---|---|---|---|---|---|---|"]
    for x in p2["rows"]:
        L.append(f"| {x['set']} | {x['omega']} | {x['truth_f']} | {x['P_U_ge_tau']:.3f} | "
                 f"{x['P_U_ge_tau_with_f0_grid']:.3f} | {x['P_L_le_tau']:.3f} | {_pf(x['pass'])} |")
    L += ["", "## Leak toy (sec. 3.5)", "",
          f"- **T1** (strict invariance to the outer laboratory's own offset; {toy['T1']['cases']} perturbations): max "
          f"relative difference {toy['T1']['max_relative_difference']:.2e} -> {_pf(toy['T1']['pass'])}",
          f"- **T1b** (the real pipeline on synthetic record files; {res['T1b']['cases']} calibrations): max relative "
          f"difference {res['T1b']['max_relative_difference']:.2e} -> {_pf(res['T1b']['pass'])}", "",
          "T2 (reported): corr(|outer level error|, tau(S)) strict / one-deep per outer laboratory, against the synthetic "
          "7th laboratory's corr(|level error|, tau(all)); the one-deep slope d tau(S) / dD; one-deep minus strict "
          "(b) coverage and LPD.", "",
          "| set | kappa | tau | " + " | ".join(SCORABLE) + " | 7th lab corr (cov90) | 1-deep - strict cov90 / LPD |",
          "|---|---|---|" + "---|" * len(SCORABLE) + "---|---|"]
    for rec in toy["T2"]:
        cells = []
        for lab in SCORABLE:
            c = rec["per_lab"][lab]["corr_abs_level_error_vs_tau_S"]
            sl = rec["per_lab"][lab]["onedeep_slope"]
            parts = [f"{k.split('_')[0]} {_f(c.get(k + '_strict'), '{:+.2f}')}/{_f(c.get(k + '_onedeep'), '{:+.2f}')}"
                     for k in (ABS, FIG) if f"{k}_strict" in c]
            sk = [f"{_f(v, '{:+.2f}')}" for k, v in sl.items() if k.endswith(f"D_{ABS}") and k.startswith(f"d_tau_{ABS}")]
            cells.append("; ".join(parts) + (f"; slope {sk[0]}" if sk else ""))
        sv = rec["seventh_lab"]
        L.append(f"| {rec['set']} | {rec['kappa']} | {rec['tau']} | " + " | ".join(cells)
                 + f" | {_f(sv['corr_abs_level_error_vs_tau_all'], '{:+.2f}')} ({sv['b_cov90']:.2f}) | "
                 f"{rec['onedeep_minus_strict']['cov90_lab_weighted']:+.3f} / {rec['onedeep_minus_strict']['lpd_lab_weighted']:+.3f} |")
    L += ["", "## Choices (where the prereg is ambiguous; the reading taken)", ""]
    L += [f"{i + 1}. {c}" for i, c in enumerate(res["choices"])]
    L += ["", f"*Wall time {res['wall_seconds']} s.*"]
    return "\n".join(L)


def render_score(res: Dict[str, Any]) -> str:
    v = res["verdict"]
    run = res["run"]
    rg = res["reproduction_guards"]
    d2 = res["D2"]
    L = ["# H0: out-of-sample between-laboratory tau for the thiol lane (R3En), nested"
         + ("" if run["kind"] == "primary" else f" ({run['kind']})"), "",
         f"*Generated by `{GENERATED_BY} --score`; pre-registration `{PREREG.name}` (v3, frozen at "
         f"{res['freeze']['freeze_commit']}). Residual = predicted minus measured (dex). Every verdict is PROVISIONAL "
         "(sec. 9).*", ""]
    if run["amendment"]:
        L += [f"*{run['kind'].capitalize()}; reason: {run['amendment']}*", ""]
    L += [f"**Verdict: {v['reported']} -- PROVISIONAL**" + (" (FRAGILE)" if v["fragile"] else "")
          + (f"; the sets disagree: {v['per_set']}" if v["sets_disagree"] else ""), "",
          f"- Guard 4.1: {_pf(res['guard']['pass'])} (max rel {res['guard']['max_relative_difference']:.1e}); "
          f"controls: hard failures {res['controls_gates']['hard_failures'] or 'none'}",
          f"- Reproduction guards inside --score: training_terms {_pf(rg['terms']['pass'])} ({rg['terms']['n_checks']} "
          f"terms, max rel {rg['terms']['max_relative_difference']:.1e}); engine guard 4.2(c) "
          f"{_pf(rg['engine_guard_4_2c']['pass'])} (cost rel {rg['engine_guard_4_2c']['cost']['relative_difference']:.1e}, "
          f"held-out raw max {rg['engine_guard_4_2c']['held_out_raw']['max_abs_difference_dex']:.1e} dex); freeze files "
          f"changed since {res['freeze']['freeze_commit']}: {res['freeze'].get('changed_since_freeze') or 'none'}",
          f"- Pair runner: {res['runner_unit_tests']}",
          f"- D2: tau_ABS(all) = {d2['tau_ABS_all_used']:.3f} dex (larger set), {d2['U_text']} -> "
          f"branch {d2['branch']}: {d2['text']}" + (f" {d2['bounds_label']}" if d2['bounds_label'] else "")
          + (f"; {d2['unstable_note']}" if d2['unstable_note'] else ""),
          "- D2 bounds per set (f = 1 / f = 0): " + "; ".join(
              f"{s} {x['f1']['U_text']}, {x['f1']['L_text']} / {x['f0']['U_text']}, {x['f0']['L_text']}"
              for s, x in d2["per_set"].items()),
          f"- FRAGILE (two sets combined by precedence per dropped laboratory): {v['fragility']['combined_per_drop']}", ""]
    L += [f"- owner: {n}" for n in res["owner_notes"]]
    L += [f"- control flag: {n}" for n in res["control_flags"]["flags"]]
    for s, sr in res["sets"].items():
        eng = sr["models"]["engine"]
        L += ["", f"## Set '{s}': {sr['verdict']}"
              + (" (its own jackknife changes the verdict)" if sr["fragile_within_set"] else ""), "",
              f"- D1: (b) 90 % coverage {sr['D1']['cov90_lab_weighted']:.3f} lab-weighted, {sr['D1']['cov90_pooled']:.3f} "
              f"pooled -> {_pf(sr['D1']['pass'])} (Acceptance band {ACCEPTANCE_BAND}: {sr['D1']['in_acceptance_band']})",
              f"- strict nesting incomplete (non-converged pair fits): {sr['strict_nesting_incomplete'] or 'no'}",
              f"- inner laboratories emptied by the strict masks (minimum-data rule, flagged): "
              f"{sr['strict_masks_emptied'] or 'none'}",
              f"- engine tau_ABS(all) {_f(eng['all'][ABS]['tau'])}, omega(all) {_f(eng['all']['omega'])}, tau_FIG(all) "
              f"{_f(eng['all'][FIG]['tau'])}; drop-one ABS {eng['drop_one_ABS']}; UNSTABLE {eng['unstable']}; "
              f"multi-row-group tau_ABS(all) {_f(eng['tau_ABS_all_multi_row_groups'])}; sensitivities "
              f"{eng['sensitivities_all']}", "",
              "| model | outer lab | omega(S) | tau_ABS(S) [status] | tau_FIG(S) used | one-deep ABS / FIG | leak |",
              "|---|---|---|---|---|---|---|"]
        for model in MODELS:
            m = sr["models"][model]
            for lab in SCORABLE:
                n_ = m["nested"][lab]
                o = m["onedeep"][lab]
                leak = m["per_lab"][lab]["onedeep_vs_strict"]
                L.append(f"| {model} | {lab} | {_f(n_['omega'])} | {_f(n_[ABS]['tau'])} [{n_[ABS]['status']}] | "
                         f"{_f(n_['tau_fig_used'])} | {_f(o[ABS]['tau'])} / {_f(o[FIG]['tau'])} | "
                         f"{', '.join(k for k, x in leak.items() if x['measured_leak']) or '-'} |")
        L += ["", "| model / arm | cov80 lw / pooled | cov90 lw / pooled | cov95 lw / pooled | 90 % width: lab-weighted "
              "mean / median of lab medians / pooled row median | LPD lw | IS90 lw |", "|---|---|---|---|---|---|---|"]
        for k, row in sr["coverage_table"].items():
            L.append(f"| {k} | " + " | ".join(f"{row[f'{lv:.2f}']['lab_weighted']:.2f} / {row[f'{lv:.2f}']['pooled']:.2f}"
                                             for lv in LEVELS)
                     + f" | {row['width90_mean_lab_weighted']:.2f} / {row['width90_median_of_lab_medians']:.2f} / "
                     f"{row['width90_median_pooled_rows']:.2f} | {row['lpd_lab_weighted']:+.3f} | "
                     f"{row['is90_lab_weighted']:.3f} |")
        L += ["", "| comparison | mean | SE | sign-flip p | exact bootstrap 90 % | IS mean | result |", "|---|---|---|---|---|---|---|"]
        for null in NULLS:
            for name, d in (("D3 (b) engine - (b) " + null, sr["D3"][null]),
                            ("ABS-only " + null, sr["abs_only"][null]),
                            ("relative branch " + null, sr["relative_branch"][null]),
                            ("St " + null, sr["St_paired"][null]), ("St4 " + null, sr["St4_paired"][null])):
                L.append(f"| {name} | {d['mean']:+.3f} | {d['se']:.3f} | {d['sign_flip_p_better']:.3f} "
                         f"({d['n_patterns']}) | [{d['bootstrap90'][0]:+.3f}, {d['bootstrap90'][1]:+.3f}] | "
                         f"{d['is_mean']:+.3f} | {d['label']} |")
            c = sr["concentration"][null]
            L.append(f"| concentration {null} | top {c['top_rows']} of {c['rows']} rows carry {c['signed_share']:+.2f} "
                     f"(abs {_f(c['abs_share'], '{:.2f}')}) | | | | | {'driven by few rows' if c['driven_by_few_rows'] else '-'} |")
        d = sr["b_vs_a_engine"]
        L.append(f"| (b) - (a), engine | {d['mean']:+.3f} | {d['se']:.3f} | {d['sign_flip_p_better']:.3f} | "
                 f"[{d['bootstrap90'][0]:+.3f}, {d['bootstrap90'][1]:+.3f}] | {d['is_mean']:+.3f} | reported |")
        L += ["", f"Jackknife: {sr['jackknife']}", "", "D-all (engine, every finite row; coverage 90 % (a) / (b), median "
              "width (a) / (b)): " + "; ".join(f"{lab} {x['a']['rows']} rows {_f(x['a']['cov0.90'], '{:.2f}')} / "
                                              f"{_f(x['b']['cov0.90'], '{:.2f}')}, {_f(x['a']['width90_median'], '{:.2f}')} / "
                                              f"{_f(x['b']['width90_median'], '{:.2f}')}" for lab, x in sr["D_all"].items())]
    L += ["", "## Choices", ""] + [f"{i + 1}. {c}" for i, c in enumerate(res["choices"])]
    return "\n".join(L)


# ----------------------------------------------------------------------------------------------------------------
# 13. CLI
# ----------------------------------------------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--guard", action="store_true", help="sec. 4.1 only: pass / fail and the max relative difference")
    g.add_argument("--design", action="store_true", help="the design table (row keys, classes, groups, sigmas)")
    g.add_argument("--controls", action="store_true", help="secs. 8, 3.5 and the D2 grid on synthetic residuals")
    g.add_argument("--score", action="store_true", help="the real scoring, run ONCE after the ten pair records")
    parser.add_argument("--replicates", type=int, default=R_CONTROLS, help=argparse.SUPPRESS)
    parser.add_argument("--amendment", default="", help="--score after a frozen file changed: the recorded reason")
    args = parser.parse_args(argv)
    output = None
    if args.score:                    # refuse before anything is loaded
        output = score_output(args.amendment)   # run once (sec. 7.5; amendment A1)
        check_pair_records()
    with engine_state():
        if args.guard:
            ok, worst, failed = guard()
            print(f"guard 4.1: {'PASS' if ok else 'FAIL'}; max relative difference {worst:.2e}")
            if failed:
                print("failed checks: " + "; ".join(failed))
            return 0 if ok else 1
        if args.design:
            print(render_design(Design().table()))
            return 0
        if args.controls:
            res = run_controls(args.replicates)
            _write(CONTROLS_STEM, res, render_controls(res))
            s = res["summary"]
            print(f"controls: hard failures {s['hard_failures'] or 'none'}; wrote {CONTROLS_STEM}.{{json,md}}")
            return 1 if s["hard_failures"] else 0
        res = run_score(args.amendment, output=output)
        _write(output["stem"], res, render_score(res))
        print(f"score: {res['verdict']['reported']} (PROVISIONAL); wrote {output['stem']}.{{json,md}}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
