#!/usr/bin/env python
"""
Which mechanism steps make the engine's within-study contrast errors? (2026-11-07; a DIAGNOSTIC on spent data.)

WHY. Step 0 of the decision-relevance pre-registration (``results/validation/decision_relevance_step0.md``, PR #110)
measured the engine's out-of-sample within-study ratio error sigma_r ~ 1.27 dex on the sulfur lane, twice the iid
0.59 shipped, with a negative implied between-arm correlation: the engine over-responds to levers. The independent
check of step 0 found that two pH-lever studies scored under the FALSIFIED B0 charge balance (Xu 2007, Hofmann 1998
pH 3 / 7) carry most of it (sigma_r ~ 0.51 without them). This record asks which STRUCTURES carry the rest.

WHAT IT DOES. The units are step 0c's: every out-of-sample within-study pair on main (``contrast_scores.json``) and
every answered within-study row of the frozen contrast benchmark (``kinetic_core_contrast_benchmark_score_v1.json``),
grouped per (study, compound, lane) and connected component, each group's arms centred (the estimator's own
algebra; the grouping, centring, sigma_r and slope code is the step-0 checker's ``recompute.py``, vendored). For
every arm's pot it re-runs the SHIPPED NumPy engine (``engine.predict``, main's constants) and takes, by central
differences, the sensitivity of log10(prediction) to every coordinate the pot's lane integrates: log10 k_ref
(h = 0.1 dex) and Ea (h = 2 kJ/mol) of each step, the sulfur lane's two pH-drift constants, the acrylamide lane's pH
exponents and a_w multiplier. A second run per pot holds the pH at its label (``BufferSpec(kind="clamped")``): the
B0 balance switched off. Nothing is fitted; no constant moves.

    s_aj = d log10 pred_a / d theta_j        (arm a, coordinate j; centred within its group)
    e_a  = centred (pred - measured) residual of arm a (the sigma_r estimator's units)

* the share of the pooled residual sum of squares sum e^2 that ONE coordinate could remove, (s_j . e)^2 / (s_j . s_j),
  the shift that would do it, -(s_j . e) / (s_j . s_j), and its leave-one-study-out (LOSO) predictive version (the
  shift fitted on the other studies, applied to the held-out one); step groups by joint projection (chance share
  ~ rank / rows);
* the B0 switch: the residual with the balance off (the pH held), per pair and pooled;
* "no lever response": groups whose prediction does not move between arms and no coordinate moves it (a lever the
  structure lacks); their SS is attributable to no constant;
* the slope beta of observed on predicted (centred arms, through the origin, cluster-robust by study): per lane, per
  lever class, LOSO; beta << 1 = over-response;
* every coordinate cross-referenced against its S1 tier, S4 band, the frozen-fit drift record (PR #107) and the open
  CF entries that name it.

Every subset (B0-pH-dependent or not, a lever class) is chosen at PAIR level, then regrouped. A pair is
B0-pH-dependent when switching the balance off moves its predicted log ratio by more than the lane's within-laboratory
term omega = 0.419 dex or flips its sign (S3's refusal clause applied to the switch); PR #108's own refusal of each
source-1 pair is listed beside it.

    ./scripts/docker_maillard.sh run "nice -n 10 python scripts/generators/build_contrast_discrepancy_decomposition.py --time-one"
    ./scripts/docker_maillard.sh run "nice -n 10 python scripts/generators/build_contrast_discrepancy_decomposition.py --workers 3"
    ./scripts/docker_maillard.sh run "python scripts/generators/build_contrast_discrepancy_decomposition.py --check"
    # writes results/validation/contrast_discrepancy_decomposition.{json,md}; --check re-derives the analysis from the
    # json's stored sensitivities (no integration) and compares
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "scripts" / "generators"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "contrast_discrepancy_decomposition.json"
CONTRAST_SCORES = data_paths.VALIDATION_DIR / "contrast_scores.json"
PANEL_SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
CONTRAST_BENCHMARK = data_paths.VALIDATION_DIR / "kinetic_core_contrast_benchmark_score_v1.json"
TIERS = data_paths.VALIDATION_DIR / "incumbent_evidence_tiers.json"
DRIFT = data_paths.VALIDATION_DIR / "frozen_fit_reproducibility.json"
CF_REGISTER = ROOT / "docs" / "reference" / "conflicts_register.yml"
PR108_COMMIT = "5e61fed0"          # PR #108's head (branch s3-ph-balance-refusal) when this record was built
PR108_CONTRAST = "results/validation/contrast_scores.json"

OMEGA_DEX = 0.419          # H0's within-laboratory term (EX-059): S3's tolerance for a relative output
SIGN_DEX = 0.01            # ph_balance_refusal.SIGN_FLOOR_DEX: a log ratio this small has no sign
H_LOG10K = 0.1             # central-difference step on log10 k_ref (dex)
H_EA = 2.0                 # on Ea (kJ/mol); the drift record's Ea step
H_PH_DRIFT_LOG = 0.1       # on log10 acid_yield
H_PKA = 0.1                # on the ARP amine pKa
H_EXPONENT = 0.02          # on a pH exponent (decades per pH unit)
H_AW_SCALE = 0.05
ZERO_UG = 1e-9             # a prediction at or below this is a zero: no log, no sensitivity
FLAT_SENS = 1e-3           # a centred group whose every coordinate moves it less than this (per unit) is flat
ARM_TIMEOUT_S = 300
CV_SIDA, CV_OTHER = 0.10, 0.20
MAX_WORKERS = 3


# ============================================================================== the step-0 checker's algebra (vendored)
# From the independent step-0 check (rev/recompute.py, 2026-11-07): load / arms_from_pairs / groups / est / beta, kept
# in what they compute; the additions are the pot keys each arm carries and the lever class / B0 tag of each pair.

def _cv_var(cv: float, n: int = 2) -> float:
    return math.log10(1.0 + cv) ** 2 / n


def load_units(cs: Mapping[str, Any], panel: Mapping[str, Any], cb: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """The out-of-sample within-study pairs: source 1 (contrast_scores, out of sample) and source 3 (the contrast
    benchmark's answered, uncensored rows). p = engine log ratio, m = measured, both log10(a / b)."""
    qv = {}
    for b in panel["benchmarks"]:
        q = str(b.get("quantification_source") or "").lower()
        qv[b["benchmark_id"]] = _cv_var(CV_SIDA if ("isotope" in q or "sida" in q or "stable" in q) else CV_OTHER)
    pairs = []
    for p in cs["pairs"]:
        if not p["out_of_sample"]:
            continue
        pairs.append(dict(source=1, study=p["study"], comp=p["compound"], lane=p["lane"], a=p["a"], b=p["b"],
                          pot_a=("bundle", p["a"]), pot_b=("bundle", p["b"]), target=p["compound"],
                          p=float(p["engine_dex"]), m=float(p["measured_dex"]),
                          va=qv.get(p["a"], _cv_var(CV_OTHER)), vb=qv.get(p["b"], _cv_var(CV_OTHER))))
    for r in cb["primary"]["decision_set"]["rows"]:
        if r.get("m") is None or r.get("p") is None or r.get("kind") == "censored":
            continue
        va = r["sigma_m"] ** 2 / 2 if r.get("sigma_m") is not None else _cv_var(CV_OTHER)
        pairs.append(dict(source=3, study=r["table"], comp=r["species"], lane="sulfur", a=r["a"], b=r["b"],
                          pot_a=("table", r["table"], r["a"]), pot_b=("table", r["table"], r["b"]),
                          target=r["species"], p=float(r["p"]), m=float(r["m"]), va=va, vb=va))
    return pairs
