#!/usr/bin/env python
"""
fit_joint_deoxyosone_thiol.py -- the joint sub-model 3 + 4 fit (results/validation/kinetic_core_thiol_submodel_prereg.md
sec. 4, and its design review; 2026-10-30).

Sub-model 3's parent arm (k_dpo_ddp slowed, k_dpo_c2c3 raised) closes the calibrated methylglyoxal but starves the
H2S that MFT (from norfuraneol) and FFT (from furfural) need: the cysteine Strecker step releases H2S (B59a). The
constants are fitted together on the deoxyosone and thiol tables, with the H2S loss free. Engine state as sub-model 3
(route switches B58 / B59b / B59a / B61 on, sub-model 1's off, inline coupling).

After its design review (thiol prereg sec. 5): coordinates are log10 multipliers on the shipped k_ref of k_dpo_ddp
(down to 10^-8), k_dpo_c2c3, k_dpo_nf, k_nf_mft, BOTH FFT branches together (k_fur_fft, k_fur_fft_hs), k_h2s_loss and
k_mft_decay (``--mft-decay-cap``: the CF-013-capped rival, at most 1/5 of the shipped loss), and the norfuraneol ->
MGO share. The B59a Strecker H2S share is held at a registered CR-11 alternative (``--strecker-share``).
Rows: both tables, deduplicated; OUT as declared hold-outs (FIT_HOLDOUT_DECLARATION): every Zhai 140 C pot, Yiltirak,
Hofmann T1 and T2 pH 3 / 7; SCORE-ONLY (no coordinate moves them): Kumazawa survival, Hofmann T8 and the T10
thiamine pot; dropped: glyoxal and the hexose + cysteine pots. Thiols semi-quantitative per figure with tau = inf
(CR-15), Zhou's thiols and MGO a per-paper scale (CF-014, CF-060), the quinoxalines tau 0.3 (CR-16). Extra noise
per data set x observable from within-series scatter, with a declared 0.10 dex floor where it cannot be estimated.
2-Acetylthiazole in Zhou 2023's ARP + Cys pot (11.70 ug/L) is scored as a falsifier, not fitted.

    python scripts/generators/fit_joint_deoxyosone_thiol.py --strecker-share 0.042|0.080 [--mft-decay-cap] [--lopo]
Writes results/validation/kinetic_core_joint_deoxyosone_thiol_fit_s<share x 1000>[_capped].{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, fit_rows, fitting, parallel, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_deoxyosone_submodel as sm3  # noqa: E402
from fit_ph_acid_submodel import pooled_noise, raw_table  # noqa: E402

TABLES = ("deoxyosone_submodel", "thiol_submodel")
#: coordinate -> the parameter keys its multiplier moves
MULTIPLIED = {"k_dpo_ddp": ("k_dpo_ddp",), "k_dpo_c2c3": ("k_dpo_c2c3",), "k_dpo_nf": ("k_dpo_nf",),
              "k_nf_mft": ("k_nf_mft",), "fft_formation": ("k_fur_fft", "k_fur_fft_hs"), "k_h2s_loss": ("k_h2s_loss",),
              "k_mft_decay": ("k_mft_decay",)}
START = {"k_dpo_ddp": -1.0, "k_dpo_c2c3": 1.5, "k_dpo_nf": 0.0, "k_nf_mft": 0.0, "fft_formation": 0.0,
         "k_h2s_loss": 0.0, "k_mft_decay": -0.3}
LOWER = {"k_dpo_ddp": -8.0}
#: round-2 arms (thiol prereg sec. 6): R2a frees the methylglyoxal Strecker of cysteine (log10 of its k at 100 C /
#: pH 5, L/(mmol min)); R2b frees k_cys_actz on Zhou 2023's pH 7 row and scores pH 6 / 8
ARMS = {"none": frozenset(), "R2a": frozenset({"mgo_strecker"}), "R2b": frozenset({"actz"}),
        "R2ab": frozenset({"mgo_strecker", "actz"}),
        # the pentose + amine entry (kinetic_core_pentose_amine_entry_prereg.md secs. 7-8) on top of R2b
        "E": frozenset({"actz", "entry"}),
        # round 3 (thiol prereg sec. 9): MFT through the intact C5 skeleton, on R2b's state or on arm E's
        "R3": frozenset({"actz", "ddp"}), "R3E": frozenset({"actz", "entry", "ddp"}),
        # round 3 after its review (sec. 10): the DDP-forming step's pH law as the rival structures
        "R3n": frozenset({"actz", "ddp", "ddp_law_none"}), "R3s": frozenset({"actz", "ddp", "ddp_law_slope"}),
        "R3En": frozenset({"actz", "entry", "ddp", "ddp_law_none"}),
        "R3Es": frozenset({"actz", "entry", "ddp", "ddp_law_slope"})}
#: round 3 (review B2): the DDP route's PRODUCT (log10 multiplier of k_dpo_ddp x both DDP -> MFT branches) is the
#: coordinate; the DDP -> MFT multiplier is the product less k_dpo_ddp's own offset
DDP_COORDINATE = ("offset_ddp_route", ("k_ddp_mft", "k_ddp_mft_hs"), -4.0, 8.0, 5.5)
DDP_SLOPE_COORDINATE = ("ddp_ph_slope", "ph_slope", -1.0, 2.0, 0.0)
#: Cerny & Davidek 2003 (hold-out, sec. 9): pot D norfuraneol 72.5 + [13C5]ribose 72.5 + cysteine 48.3 umol in a
#: ~0.5 g charge (nominal mmol/L), pot C ribose 145 + cysteine 48.3 umol; 0.5 M K phosphate, pH 5.00, 95 C, 4 h
CERNY_POTS = {"cerny2003_potD_nf_ribose_cys": {"NF": 145.0, "PENT": 145.0, "Cys": 96.6},
              "cerny2003_potC_ribose_cys": {"PENT": 290.0, "Cys": 96.6}}
#: review of round 3 (sec. 10, item 3): the exact tracer share is scored two-sided; a share of 0 is no pass either
CERNY_NF_SHARE_BAND = (0.03, 0.15)
#: pot C's 2-mercapto-3-pentanone, a blind companion: Cerny 2003 finds none from ribose; in the lane it comes only
#: from norfuraneol (r_nf_mp3p). Declared: MP3P / MFT made in pot C <= 0.1 ("not detected" beside a detected MFT)
CERNY_MP3P_RATIO_MAX = 0.1
#: a pH law set by the arm rather than fitted (R3n): read by trial() in the main process, carried in every job
DDP_LAW_STATIC: Dict[str, Any] = {}
ARM_MULTIPLIED = {"k_cys_actz": ("k_cys_actz",), "k_aa_frag": ("k_aa_frag",), "ddp_mft": ("k_ddp_mft", "k_ddp_mft_hs")}
ARM_LOWER, ARM_UPPER, ARM_START = ({"k_cys_actz": -6.0, "k_aa_frag": -6.0}, {"k_cys_actz": 1.0, "k_aa_frag": 0.5},
                                   {"k_cys_actz": -3.0, "k_aa_frag": -1.5})
#: the entry's own coordinates: (key, kind, lower, upper, start). Rates as log10 of k at 100 C / pH 6 (entry
#: L/(mmol min), acetic 1/min); slopes in decades per in-situ pH unit; amine_return a share in [0, 1]
ENTRY_COORDINATES = (("log10_k_pent_amine", "log10_offset", -7.0, -2.0, -4.15),
                     ("entry_ph_slope", "ph_slope", -0.5, 1.5, 0.3),
                     ("log10_k_dpo_aa", "log10_offset", -4.0, 0.5, -0.5),
                     ("aa_ph_slope", "ph_slope", -0.5, 1.5, 0.3),
                     ("amine_return", "acid_yield", 0.0, 1.0, 0.5))
ENTRY_TABLE = "ph_acid_davidek_acids"
ENTRY_FIT_PREFIX = "davidek2008"
#: scored with the entry (prereg sec. 8): Cao 2017's cysteine + xylose + glycine pot (a hold-out) and Zhou 2023's
#: final pH (score-only: B60 fitted the ARP amine pKa and the acid yield on those endpoints)
ENTRY_PH_TABLE = "ph_acid_submodel"
ENTRY_PH_HOLDOUT_POTS = ("cao2017_gly_x1",)
ENTRY_PH_SCORE_PREFIX = "zhou2023"
MGO_STRECKER_BOUNDS, MGO_STRECKER_START = (-6.0, 0.0), -3.0
ACTZ_TABLE = "thiol_zhou2023_actz"
ACTZ_FIT_PH = 7.0
#: CF-014: 2-acetylthiazole in the same nominal pot, Zhou 2022 vs Zhou 2023, agrees to x1.4 (MFT's x7.8 is not
#: its spread; round-2 design review B3)
ACTZ_FLOOR_DEX = 0.15
#: declared hold-outs (FIT_HOLDOUT_DECLARATION.md) and score-only rows (design review B3, S3)
HOLDOUT_LEDGER_PREFIXES = ("yiltirak2026", "hofmann1998_t1_", )
HOLDOUT_POTS = ("hofmann1998_cys_ribose_ph3", "hofmann1998_cys_ribose_ph7")
SCORE_ONLY_LEDGER_PREFIXES = ("kumazawa2003", "hofmann1998_t8_")
SCORE_ONLY_POTS = ("hofmann1998_thiamin_ph5",)
NOISE_FLOOR_DEX = 0.10
#: the falsifier (design review B4): Zhou 2023 Table 1, ARP + Cys 20 + 20 mM, 120 C, 60 min, pH 7, unbuffered
ACTZ_FALSIFIER = {"pot": "zhou2023_arp_cys_ph7_actz", "ledger_id": "zhou2023 T1 (B60 fit row zhou_pH7_ACTZ)",
                  "t_c": 120.0, "label_ph": 7.0, "charge": {"ARP": 20.0, "Cys": 20.0}, "dropped": [], "times": [0.0, 60.0],
                  "buffer": "none (water)", "notes": [], "measured_ug_l": 11.70}
HEXOSES = ("Glc", "Fru")
OUT = data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit"
POT_SECONDS = 120
FAILURES: List[str] = []


def coordinates(mft_decay_cap: bool, arms=frozenset()):
    keys = ["mgo_share_nf"] + [f"offset_{k}" for k in MULTIPLIED]
    kinds = ["acid_yield"] + ["log10_offset"] * len(MULTIPLIED)
    lower = [0.0] + [LOWER.get(k, -4.0) for k in MULTIPLIED]
    upper = [1.0] + [math.log10(0.2) if (k == "k_mft_decay" and mft_decay_cap) else 3.0 for k in MULTIPLIED]
    start = [0.0142] + [min(START[k], u) for k, u in zip(MULTIPLIED, upper[1:])]
    if "actz" in arms:
        for k in ("k_cys_actz",):
            keys.append(f"offset_{k}"); kinds.append("log10_offset")
            lower.append(ARM_LOWER[k]); upper.append(ARM_UPPER[k]); start.append(ARM_START[k])
    if "entry" in arms:
        keys.append("offset_k_aa_frag"); kinds.append("log10_offset")
        lower.append(ARM_LOWER["k_aa_frag"]); upper.append(ARM_UPPER["k_aa_frag"]); start.append(ARM_START["k_aa_frag"])
        for key, kind, lo, hi, st in ENTRY_COORDINATES:
            keys.append(key); kinds.append(kind); lower.append(lo); upper.append(hi); start.append(st)
    if "ddp" in arms:
        key, _, lo, hi, st = DDP_COORDINATE
        keys.append(key); kinds.append("log10_offset"); lower.append(lo); upper.append(hi); start.append(st)
    if "ddp_law_slope" in arms:
        key, kind, lo, hi, st = DDP_SLOPE_COORDINATE
        keys.append(key); kinds.append(kind); lower.append(lo); upper.append(hi); start.append(st)
    if "mgo_strecker" in arms:
        keys.append("log10_k_mgo_cys_strecker"); kinds.append("log10_offset")
        lower.append(MGO_STRECKER_BOUNDS[0]); upper.append(MGO_STRECKER_BOUNDS[1]); start.append(MGO_STRECKER_START)
    return tuple(keys), tuple(kinds), np.array(lower), np.array(upper), np.array(start, dtype=float)


def trial(theta, keys) -> Dict[str, Any]:
    t = dict(zip(keys, map(float, theta)))
    base = engine.core_parameters(engine.SULFUR)
    multiplied = {**MULTIPLIED, **ARM_MULTIPLIED}
    k_ref = {key: base[key].k_ref * 10.0 ** v for coord, v in t.items()
             if coord.startswith("offset_") and coord != DDP_COORDINATE[0]
             for key in multiplied[coord[len("offset_"):]]}
    if DDP_COORDINATE[0] in t:   # the route product less the forming step's own offset
        for key in DDP_COORDINATE[1]:
            k_ref[key] = base[key].k_ref * 10.0 ** (t[DDP_COORDINATE[0]] - t["offset_k_dpo_ddp"])
    routes = {"k_tdp_mgo_100c_ph5": 0.0, "ea_tdp_mgo_kj_mol": 100.0, "k_pent_mgo_100c_ph5": 0.0, "ea_pent_mgo_kj_mol": 100.0}
    out = {"k_ref": k_ref, "routes": routes, "mgo_share": t["mgo_share_nf"]}
    if "log10_k_mgo_cys_strecker" in t:
        out["mgo_strecker_k"] = 10.0 ** t["log10_k_mgo_cys_strecker"]
    if "ddp_ph_slope" in t:
        out["ddp_law"] = {"law": "slope", "slope": t["ddp_ph_slope"]}
    elif DDP_LAW_STATIC:
        out["ddp_law"] = dict(DDP_LAW_STATIC)
    if "log10_k_pent_amine" in t:
        out["entry"] = {"k_entry_l_per_mmol_min_100c_ph6": 10.0 ** t["log10_k_pent_amine"],
                        "entry_ph_slope": t["entry_ph_slope"], "k_dpo_aa_per_min_100c_ph6": 10.0 ** t["log10_k_dpo_aa"],
                        "aa_ph_slope": t["aa_ph_slope"], "amine_return": t["amine_return"]}
        out["trunk_solutes"] = True     # design review B2: glycine (0.1 M, amine pKa ~8 at 120 C) in the balance
    return out


STRECKER_SHARE = 0.042


def run_task(job):
    install(job)   # every switch set on every task: a worker runs many
    try:
        return probe.run_pot_limited(job, POT_SECONDS)
    except Exception as exc:  # noqa: BLE001 - counted; its rows become failed residuals
        return {"failed": f"{type(exc).__name__}: {exc}"}


def key_of(row) -> str:
    return f"{row.table}:{row.pot}"


def method_of(row) -> str:
    """The quinoxalines (deoxyosone table) keep sub-model 3's per-figure tau 0.3 (CR-16); the thiol table's
    semi-quantitative rows are per-figure tau = inf (CR-15, CF-065); Zhou's thiols and MGO a per-paper scale."""
    if row.ledger_id.startswith("zhou") and row.species in ("MFT", "FFT", "MFTD", "FFTD", "FUR", "MGO"):
        return "zhou_free"
    if row.method in ("semi_quantitative", "relative"):
        return "semi_figure" if row.table == "deoxyosone_submodel" else "zhou_free"
    return row.method


def series_of(row) -> str:
    """A scaled row's calibration series (thiol prereg sec. 5, B6: one scale PER FIGURE): paper x figure x species;
    Zhou's thiols and MGO one per paper x species (CF-014, CF-060). Round 1 keyed it by ledger id, which splits a
    figure's arms (Zhai JAFC Fig. 3: TTCA vs TTCA + Cys) and Zhou 2023's pH 7 column from its pH 6 / 8 columns
    (round-2 design review B2). Absolute rows keep their own series."""
    method = method_of(row)
    if method not in ("semi_figure", "zhou_free"):
        return row.series
    paper = sm3.paper_of(row)
    if method == "zhou_free" and row.ledger_id.startswith("zhou"):
        return f"{paper}|{row.species}"
    return f"{paper}|{row.ledger_id.split('_')[1]}|{row.species}"


def noise_with_floor(model, raw_residuals, rows, n_constants=0):
    """pooled_noise, and NOISE_FLOOR_DEX (declared; CF-016 / CF-055 between-run scatter) for every series whose
    noise cannot be estimated (design review S1: Garching's 10 % sd would otherwise dominate the cost)."""
    out = pooled_noise(model, raw_residuals, rows, n_constants)
    for row in rows:
        out.setdefault(model.series_of(row), NOISE_FLOOR_DEX)
        if row.species == "ACTZ":   # R2b: one row, so its floor is the laboratory's reproducibility (CF-014)
            out[model.series_of(row)] = max(out[model.series_of(row)], ACTZ_FLOOR_DEX)
    return out


def dedupe(rows):
    """Rows printed in both tables (design review B7) counted once: by (ledger id, species, time, value)."""
    seen, out = set(), []
    for r in rows:
        k = (r.ledger_id, r.species, float(r.t_min), round(float(r.value), 6))
        if k not in seen:
            seen.add(k)
            out.append(r)
    return tuple(out)


def make_raw_many(jobs, keys):
    def raw_many(points, rows):
        pots = sorted({key_of(r) for r in rows})
        trials = [trial(x, keys) for x in points]
        results = parallel.map_ordered(run_task, [dict(jobs[p], **tr) for tr in trials for p in pots])
        out = []
        for i, tr in enumerate(trials):
            res = dict(zip(pots, results[i * len(pots):(i + 1) * len(pots)]))
            r_i = []
            for r in rows:
                got = res[key_of(r)]
                if "failed" in got:
                    FAILURES.append(f"{key_of(r)}: {got['failed']}")
                    r_i.append(float("nan"))
                    continue
                predicted = probe.predict(r, dict(jobs[key_of(r)], **tr), got)
                r_i.append(fit_rows.residual(r, predicted) if predicted > 0 else float("nan"))
            out.append(np.array(r_i, dtype=float))
        return out

    return raw_many


def after_scale(model, rows, raw_vals, offsets):
    """Raw residuals minus each series' fitted calibration offset: what a semi-quantitative series' done-test reads.
    A hold-out row takes the offset fitted on its series' fit rows; a scaled series with none is left out (its level
    is unknown, so a raw residual there says nothing)."""
    kept, vals = [], []
    for r, v in zip(rows, raw_vals):
        series = model.series_of(r)
        if model.scale_tau[model.method_of(r)] != 0.0 and series not in offsets:
            continue
        kept.append(r)
        vals.append(v - offsets.get(series, 0.0))
    return raw_table(kept, vals)


def scale_label(model, rows) -> Dict[str, str]:
    """data set | observable -> 'absolute', 'free' (tau inf) or 'tau 0.3'."""
    out = {}
    for r in rows:
        tau = float(model.scale_tau[model.method_of(r)])
        out[f"{r.ledger_id} | {r.species}"] = "absolute" if tau == 0 else ("free" if math.isinf(tau) else f"tau {tau:g}")
    return out


def cerny_route_checks(theta, keys) -> Dict[str, Any]:
    """Round 3's hold-out (thiol prereg sec. 9), scored once at the optimum: in Cerny 2003's pot D the share of MFT
    formed through free norfuraneol (r_nf_mft; measured unlabelled share 0.07, pass <= 0.15), with the fed share of
    the norfuraneol pool; in pot C the C2 + C3 route's share (Cerny 2004: "not relevant")."""
    import numpy as np
    from build_sulfur_flux_attribution import _ranked

    out = {}
    for pot, charge in CERNY_POTS.items():
        job = dict(trial(theta, keys), pot=pot, ledger_id="cerny2003 (hold-out)", t_c=95.0, label_ph=5.0,
                   charge=charge, dropped=[], notes=[], buffer="potassium phosphate 0.5 M",
                   strecker_share=STRECKER_SHARE,
                   times=sorted({0.0, *np.geomspace(1e-5, 240.0, 300).tolist(), *np.linspace(0.0, 240.0, 121).tolist()}))
        install(job)
        run = probe.integrate_pot(job)
        flux = sulfur.reaction_flux_integrals(run)
        made = _ranked(flux, "MFT")["made"]
        share = {t["reaction"]: t["share"] for t in made["top"]}
        tracer = fed_norfuraneol_tracer(run) if charge.get("NF") else None
        out[pot] = {"mft_made_mmol_L": made["total_mmol_L"], "share_by_step": share,
                    "nf_route_flux_share": share.get("r_nf_mft", 0.0),
                    "c2_c3_route_share": share.get("r_ha_mp_mft", 0.0),
                    "mp3p_over_mft": (_ranked(flux, "MP3P")["made"]["total_mmol_L"] / made["total_mmol_L"]
                                      if made["total_mmol_L"] > 0 else None),
                    **(tracer or {})}
    d, c = out["cerny2003_potD_nf_ribose_cys"], out["cerny2003_potC_ribose_cys"]
    lo, hi = CERNY_NF_SHARE_BAND
    out["verdict"] = (f"pot D unlabelled (fed-norfuraneol) MFT share {d['fed_nf_mft_share']:.3f} vs measured 0.07: "
                      + ("PASS" if lo <= d["fed_nf_mft_share"] <= hi else "FAIL")
                      + f"; pot C MP3P / MFT {c['mp3p_over_mft']:.3g} (blind; <= {CERNY_MP3P_RATIO_MAX}): "
                      + ("PASS" if (c["mp3p_over_mft"] or 0.0) <= CERNY_MP3P_RATIO_MAX else "FAIL"))
    return out


def fed_norfuraneol_tracer(run) -> Dict[str, float]:
    """The exact unlabelled share of Cerny's pot D (review of round 3, item 3): the fed norfuraneol is a tracer pool
    lost at the pool's per-molecule rate and never replenished, so its fraction f(t) of the norfuraneol pool obeys
    d ln f / dt = -(production of NF) / NF; MFT from fed norfuraneol is the integral of r_nf_mft x f, over all MFT
    formation."""
    import numpy as np

    node_times, k_nodes = run.rate_nodes
    times = np.asarray(run.times_min, dtype=float)
    keys = [r.key for r in sulfur.FULL_REACTIONS]
    nu_nf = np.array([r.products.get("NF", 0) - r.reactants.get("NF", 0) for r in sulfur.FULL_REACTIONS], float)
    nu_mft = np.array([r.products.get("MFT", 0) for r in sulfur.FULL_REACTIONS], float)
    i_nf_mft = keys.index("r_nf_mft")
    nf = run.series("NF")
    rates = np.array([sulfur.sulfur_reaction_rates(run.concentrations[i], sulfur.rate_vector_at(node_times, k_nodes, t))
                      for i, t in enumerate(times)])
    production = rates @ np.where(nu_nf > 0, nu_nf, 0.0)
    g = np.where(nf > 0, production / np.maximum(nf, 1e-300), 0.0)
    log_f = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(times))])
    f = np.exp(-log_f)
    from_fed = float(np.trapezoid(rates[:, i_nf_mft] * f, times))
    total = float(np.trapezoid(rates @ nu_mft, times))
    return {"fed_nf_mft_share": from_fed / total if total > 0 else float("nan"),
            "fed_fraction_of_nf_pool_at_end": float(f[-1])}


def install(job) -> None:
    """The engine state a job runs in, set in whichever process runs it (every switch set, none left over)."""
    sm3.install(job)
    import dataclasses as _dc
    from src.kinetic_core import parameters_sulfur as _ps
    from src.kinetic_core import ph_state as _ph
    _ps.STRECKER_SULFUR = _dc.replace(_ps.STRECKER_SULFUR, h2s_share=float(job.get("strecker_share", STRECKER_SHARE)))
    _ps.MGO_CYSTEINE_STRECKER = _dc.replace(_ps.MGO_CYSTEINE_STRECKER,
                                            k_l_per_mmol_min_100c_ph5=float(job.get("mgo_strecker_k", 0.0)))
    _ps.PENTOSE_AMINE_ENTRY = _ps.PentoseAmineEntry(**job.get("entry", {}))
    _ps.DDP_FORMATION_PH = _ps.DdpFormationPh(**job.get("ddp_law", {"law": "base"}))
    _ph.TRUNK_SOLUTES_IN_BALANCE = bool(job.get("trunk_solutes", False))


def entry_checks(theta, keys, jobs) -> Dict[str, Any]:
    """Prereg sec. 8 (design review S6, S7): Davidek 2008's norfuraneol as a sanity bound (<= 5 mol%; Desclaux's ~2.7
    %) and the censored xylose at pH 8 (60 - 240 min, on the axis: <= 0.5 mmol/L)."""
    out = {}
    for key, job in jobs.items():
        if not key.split(":")[1].startswith(ENTRY_FIT_PREFIX):
            continue
        j = dict(job, **trial(theta, keys), times=sorted({*job["times"], 60.0, 120.0, 240.0}))
        res = run_task(j)
        if "failed" in res:
            out[key] = {"failed": res["failed"]}
            continue
        states = dict(zip(j["times"], res["states"]))
        out[key] = {"nf_max_mol_percent": 100.0 * max(s["NF"] for s in res["states"]) / j["charge"]["PENT"],
                    "xylose_mmol_l": {f"{t:g}": states[t]["PENT"] for t in (60.0, 120.0, 240.0)},
                    "acetic_mmol_l": {f"{t:g}": states[t]["AA"] for t in (60.0, 120.0, 240.0)}}
    return out


def high_correlations(lap, keys, cutoff: float = 0.9) -> List[str]:
    """Pairs of coordinates whose Laplace correlation exceeds the declared cutoff (prereg sec. 8: not identified)."""
    corr = lap.get("correlation")
    if not corr:
        return []
    return [f"{keys[i]} x {keys[j]}: {corr[i][j]:+.2f}" for i in range(len(keys)) for j in range(i + 1, len(keys))
            if corr[i][j] is not None and abs(corr[i][j]) > cutoff]


def actz_falsifier(theta, start, keys) -> Dict[str, float]:
    """2-Acetylthiazole in Zhou 2023's ARP + Cys pot at the optimum and at the start (design review B4)."""
    from src.kinetic_core.observations import amount_in

    out = {}
    for label, x in (("start", start), ("optimum", theta)):
        job = dict(ACTZ_FALSIFIER, **trial(x, keys), strecker_share=STRECKER_SHARE)
        res = run_task(job)
        out[f"{label}_ug_l"] = float("nan") if "failed" in res else amount_in("ug_L", "ACTZ", res["states"][-1]["ACTZ"])
    out["log10_ratio"] = math.log10(max(out["optimum_ug_l"], 1e-30) / ACTZ_FALSIFIER["measured_ug_l"])
    return out


def _mr(cell) -> str:
    return "" if not cell else f"{cell['mean']:+.2f} / {cell['rms']:.2f}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--budget", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--lopo", action="store_true")
    parser.add_argument("--strecker-share", type=float, choices=(0.042, 0.080), default=0.042)
    parser.add_argument("--mft-decay-cap", action="store_true")
    parser.add_argument("--arm", choices=tuple(ARMS), default="none", help="round-2 arm (thiol prereg sec. 6)")
    parser.add_argument("--check-fixed-point", action="store_true",
                        help="re-score the optimum on the shipped fixed-point coupling (slow: ~1 h for these 0.4 M pots)")
    args = parser.parse_args(argv)
    global STRECKER_SHARE, OUT
    STRECKER_SHARE = float(args.strecker_share)
    # no dot in the stem (a ".042" would be taken for the suffix)
    OUT = data_paths.VALIDATION_DIR / (f"kinetic_core_joint_deoxyosone_thiol_fit_s{round(args.strecker_share * 1000):03d}"
                                       + ("_capped" if args.mft_decay_cap else "")
                                       + ("" if args.arm == "none" else f"_{args.arm.lower()}"))
    arms = ARMS[args.arm]
    global DDP_LAW_STATIC
    DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    t0 = time.time()
    sulfur.PH_COUPLING = "inline"
    predictable = set(sulfur.SULFUR_STATE_KEYS)
    jobs, rows_all, pots_all = {}, [], {}
    extra = ((ACTZ_TABLE,) if "actz" in arms else ()) + ((ENTRY_TABLE, ENTRY_PH_TABLE) if "entry" in arms else ())

    def wanted(table, pot):
        if table == ENTRY_TABLE:
            return pot.startswith(ENTRY_FIT_PREFIX)
        if table == ENTRY_PH_TABLE:
            return pot in ENTRY_PH_HOLDOUT_POTS or pot.startswith(ENTRY_PH_SCORE_PREFIX)
        return True

    for table in TABLES + extra:
        pots, rows = fit_rows.load(table)
        for j in probe.plan(pots, rows):
            if "skip" in j or any(pots[j["pot"]]["initial_mmol_L"].get(h) for h in HEXOSES) or not wanted(table, j["pot"]):
                continue
            jobs[f"{table}:{j['pot']}"] = j
            pots_all[f"{table}:{j['pot']}"] = pots[j["pot"]]
        rows_all += [r for r in rows if r.species in predictable or (table == ENTRY_PH_TABLE and r.species == "pH")]
    rows_all = list(dedupe(r for r in rows_all if key_of(r) in jobs))

    def is_holdout(r):
        if r.table == ENTRY_PH_TABLE:
            return r.pot in ENTRY_PH_HOLDOUT_POTS
        return ((r.ledger_id.startswith("zhai2023") and float(pots_all[key_of(r)]["t_c"]) >= 140.0)
                or r.ledger_id.startswith(HOLDOUT_LEDGER_PREFIXES) or r.pot in HOLDOUT_POTS)

    def is_score_only(r):
        if r.table == ENTRY_PH_TABLE:
            return True
        if r.table == ACTZ_TABLE:   # R2b: pH 7 fitted, pH 6 / 8 scored (a test of the step's pH structure)
            return float(pots_all[key_of(r)]["label_ph"]) != ACTZ_FIT_PH
        return r.ledger_id.startswith(SCORE_ONLY_LEDGER_PREFIXES) or r.pot in SCORE_ONLY_POTS

    holdout = tuple(r for r in rows_all if is_holdout(r))
    score_only = tuple(r for r in rows_all if not is_holdout(r) and is_score_only(r))
    data = tuple(r for r in rows_all if not is_holdout(r) and not is_score_only(r))
    keys, kinds, lower, upper, start = coordinates(args.mft_decay_cap, arms)
    raw_many = make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    model = fitting.MeasurementModel(
        series_of=series_of,
        sigma_of=lambda r: r.sigma, method_of=method_of, scale_tau=sm3.SCALE_TAU)

    def make_problem(m, rows=data):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: m.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [m.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=sm3.paper_of, nuisance_dof=m.nuisance_dof)

    raw_start = raw(start, data)
    best, final_model = fitting.fit_with_series_noise(make_problem, model, raw, start, budget=args.budget,
                                                      rounds=args.rounds, estimator=noise_with_floor)
    problem = make_problem(final_model)
    lap = fitting.laplace(problem, best.theta)
    raw_best = raw(best.theta, data)
    raw_holdout = raw(best.theta, holdout) if holdout else []
    offsets = final_model.offsets(raw_best, data)
    report: Dict[str, Any] = {
        "artifact": OUT.name,
        "provenance": provenance.provenance_block(
            OUT.name, generated_by="scripts/generators/fit_joint_deoxyosone_thiol.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in TABLES]
                   + [ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                      ROOT / "src/kinetic_core/fitting.py", ROOT / "scripts/generators/build_submodel_probe.py"]),
        "engine_state": "B60 + B58/B59b/B59a/B61 on; sub-model 1 switches off; inline coupling"
                        + ("" if args.arm == "none" else f"; round-2 arm {args.arm}"),
        "rows": len(data), "pots": len({key_of(r) for r in data}), "start": dict(zip(keys, map(float, start))),
        "fit": best.as_dict(keys), "laplace": lap, "series_extra_noise": dict(final_model.extra_noise),
        "failed_pot_runs": len(FAILURES), "failures_sample": FAILURES[:10],
        "raw_at_start": raw_table(data, raw_start), "raw_at_optimum": raw_table(data, raw_best),
        "strecker_share": STRECKER_SHARE, "mft_decay_cap": bool(args.mft_decay_cap),
        "holdout_raw_at_optimum": raw_table(holdout, raw_holdout) if holdout else {},
        "scale_of": scale_label(final_model, data + holdout + score_only),
        "after_scale_at_start": after_scale(final_model, data, raw_start, final_model.offsets(raw_start, data)),
        "after_scale_at_optimum": after_scale(final_model, data, raw_best, offsets),
        "holdout_after_scale_at_optimum": after_scale(final_model, holdout, raw_holdout, offsets) if holdout else {},
        "holdout_raw_at_start": raw_table(holdout, raw(start, holdout)) if holdout else {},
        "score_only_raw_at_optimum": raw_table(score_only, raw(best.theta, score_only)) if score_only else {},
        "actz_falsifier": actz_falsifier(best.theta, start, keys),
        "entry_checks": entry_checks(best.theta, keys, jobs) if "entry" in arms else None,
        "cerny_route_checks": cerny_route_checks(best.theta, keys),
        "high_correlations": high_correlations(lap, keys),
        "diagnostics_whitened": fitting.diagnostics(problem, best.theta, {
            "paper": sm3.paper_of, "data_set": lambda r: r.ledger_id, "observable": lambda r: r.species,
            "temperature_C": lambda r: float(pots_all[key_of(r)]["t_c"]),
            "label_pH": lambda r: float(pots_all[key_of(r)]["label_ph"]), "time_min": lambda r: r.t_min}),
    }
    if args.check_fixed_point:
        sulfur.PH_COUPLING = "fixed_point_nodes"
        report["shipped_coupling_max_abs_difference"] = float(np.nanmax(np.abs(raw(best.theta, data) - raw_best)))
        sulfur.PH_COUPLING = "inline"
    if args.lopo:
        out, sig = {}, np.array([lap["sigma"][k] for k in keys])
        for paper in sorted({sm3.paper_of(r) for r in data}):
            kept = tuple(r for r in data if sm3.paper_of(r) != paper)
            held = tuple(r for r in data if sm3.paper_of(r) == paper)
            f_, _m = fitting.fit_with_series_noise(lambda mm: make_problem(mm, kept), model, raw, best.theta,
                                                   budget=args.budget // 2, rounds=1, estimator=noise_with_floor)
            out[paper] = {"rows_held": len(held), "theta": dict(zip(keys, map(float, f_.theta))),
                          "shift_sigma": dict(zip(keys, map(float, (f_.theta - best.theta) / sig))),
                          "held_out_raw": raw_table(held, raw(f_.theta, held))}
        report["leave_one_paper_out"] = out
    report["wall_seconds"] = round(time.time() - t0, 1)
    OUT.with_suffix(".json").write_text(provenance.dumps_compact(report, default=float) + "\n", encoding="utf-8")
    lines = ["# Joint sub-model 3 + 4 fit (deoxyosones + thiols)", "",
             "*Generated by `scripts/generators/fit_joint_deoxyosone_thiol.py`; thiol prereg sec. 4. "
             f"{report['rows']} rows on {report['pots']} pots; engine state: {report['engine_state']}. Raw residual = "
             "predicted minus measured (dex), before any scale; 'after scale' subtracts each series' fitted calibration offset "
             "(absolute rows: none; 'free': the series' weighted mean; tau: its shrunk mean) and is what a done-test reads. "
             "A hold-out row of a scaled series takes the offset fitted on that series' fit rows; with none it is not "
             "scored.*", "", "## Result", ""]
    lines += [f"- {line}" for line in fitting.summary_lines(problem, best, lap)]
    lines += [f"- failed pot runs: {len(FAILURES)}",
              ("- largest raw-residual difference, inline vs the shipped fixed point, at the optimum: "
               f"{report['shipped_coupling_max_abs_difference']:.2g}" if "shipped_coupling_max_abs_difference" in report
               else "- inline vs fixed point not re-checked in this run (--check-fixed-point; 0.004 on the smoke run)"), "",
              "## Residuals per data set x observable", "",
              "| data set \\| observable | n | scale | start, after scale | optimum, raw | optimum, after scale |",
              "|---|---|---|---|---|---|"]
    lines += [f"| {k} | {v['n']} | {report['scale_of'][k]} | {_mr(report['after_scale_at_start'].get(k))} | {_mr(v)} | "
              f"{_mr(report['after_scale_at_optimum'].get(k))} |" for k, v in report["raw_at_optimum"].items()]
    fa = report["actz_falsifier"]
    lines += ["", "## 2-acetylthiazole, Zhou 2023 ARP + Cys pH 7 (measured 11.70 ug/L; "
              + ("FITTED in this arm: the pH 6 / 8 rows are under score-only)" if "actz" in arms else "a falsifier, not fitted)"), "",
              f"- predicted at the start {fa['start_ug_l']:.4g} ug/L, at the optimum {fa['optimum_ug_l']:.4g} ug/L "
              f"(log10 ratio at the optimum {fa['log10_ratio']:+.2f})"]
    if report["holdout_raw_at_optimum"]:
        lines += ["", "## Declared hold-outs (every Zhai 140 C pot, Yiltirak, Hofmann T1 / T2 pH 3 / 7), scored once", "",
                  "| data set \\| observable | n | scale | start, raw | optimum, raw | optimum, after the fit's scale |",
                  "|---|---|---|---|---|---|"]
        lines += [f"| {k} | {v['n']} | {report['scale_of'][k]} | {_mr(report['holdout_raw_at_start'][k])} | {_mr(v)} | "
                  f"{_mr(report['holdout_after_scale_at_optimum'].get(k)) or 'not scored (no fitted scale)'} |"
                  for k, v in report["holdout_raw_at_optimum"].items()]
    if report.get("entry_checks"):
        lines += ["", "## Entry checks (Davidek 2008 pots at the optimum)", ""]
        lines += [f"- {k}: {v}" for k, v in report["entry_checks"].items()]
    cr = report.get("cerny_route_checks") or {}
    if cr:
        lines += ["", "## Hold-out: Cerny 2003 MFT routes (95 C, pH 5; measured unlabelled MFT 0.07)", "",
                  f"- {cr.get('verdict')}"]
        lines += [f"- {k}: norfuraneol-route flux share {v['nf_route_flux_share']:.3f}, C2 + C3 route "
                  f"{v['c2_c3_route_share']:.3f}, MFT made {v['mft_made_mmol_L']:.3g} mmol/L, steps {v['share_by_step']}"
                  + (f", fed-norfuraneol share {v['fed_nf_mft_share']:.3f}" if "fed_nf_mft_share" in v else "")
                  for k, v in cr.items() if isinstance(v, dict)]
    if report.get("high_correlations"):
        lines += ["", "## Laplace correlations above 0.9 (declared: not identified)", ""]
        lines += [f"- {c}" for c in report["high_correlations"]]
    if report["score_only_raw_at_optimum"]:
        lines += ["", "## Score-only rows at the optimum", ""]
        lines += [f"- {k}: n {v['n']}, mean {v['mean']:+.2f}, RMS {v['rms']:.2f}" for k, v in report["score_only_raw_at_optimum"].items()]
    if "leave_one_paper_out" in report:
        lines += ["", "## Leave one paper out", ""]
        for paper, res in report["leave_one_paper_out"].items():
            lines.append(f"- **{paper}** ({res['rows_held']} rows): " + "; ".join(
                f"{k} {v['mean']:+.2f} / {v['rms']:.2f}" for k, v in res["held_out_raw"].items()))
    lines += ["", fitting.diagnostics_markdown(report["diagnostics_whitened"], "Residual diagnostics at the optimum (whitened)")]
    OUT.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"joint: cost {best.start_cost:.4g} -> {best.cost:.4g}; {len(FAILURES)} failed; {report['wall_seconds']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
