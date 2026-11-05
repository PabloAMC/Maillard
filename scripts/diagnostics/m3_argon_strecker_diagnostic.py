#!/usr/bin/env python
"""
m3_argon_strecker_diagnostic.py -- why every M3 rival over-predicts hofmann2000b Table 2's ARGON phenylacetaldehyde +
phenylacetic acid (2026-11-06; the open diagnostic of the M3 prereg sec. 23, "CF candidate").

A DIAGNOSTIC. Nothing is fitted, no band, prereg, src/ module or shipped output changes. It re-runs the M3 phase-1
layer's Table 2 pots (ARP-Phe 0.1 M, phosphate 0.5 M pH 7, 100 C, 120 min, closed vial; argon, air, air + Cu 5 mM)
with one shared piece swapped at a time for constants taken from cited literature on disk, and asks whether the
argon median comes within the layer's own noise of the measured 0.09 mM.

IMPORT. The layer is ``src/m3_strecker_layer.py`` on branch ``m3-phase1-build`` (not yet merged). Until the M3 PR
merges this script loads it BY PATH from ``$MAILLARD_M3_LAYER``; once merged it should import
``from src import m3_strecker_layer as L`` (tried first). The layer's sha256 is written into the record.

DATA (seen development data only). hofmann2000 Tables 3 - 4 and hofmann2000b Tables 2 - 3 and Fig. 1, through the
layer's own ``load_data`` (EX-049's calibration inputs), and the fed-dicarbonyl / fed-Amadori decay constants of the
dossiers cited per constant below. hofmann2000b Figs. 2 - 3 and their keyed rows are NOT read (the layer refuses the
keyed table; this script never names it).

PASS RULE (declared before any hypothesis was computed; committed in this file before the run)
=============================================================================================
A hypothesis EXPLAINS the miss if, with constants taken only from cited measurements (no fitting to Table 2's argon
value), for EVERY rival R1 - R4 whose feasible share is >= 5 %:
  (a) the median predicted argon PA + PAA lies within +-2 s_Ar (dex) of the measured 0.09 mM, where
      s_Ar^2 = s_read^2 + CV_hi^2 + sigma_run^2 is the layer's own noise for that total: s_read the Table 2 print
      error with s_rep (``arp_t2_err`` (+) ``S_REP_DEX``, propagated to the sum on the linear scale), CV_hi the top of
      the layer's analyser CV band (``CV_ARP_BAND[1]`` = 10 %, as dex), sigma_run the layer's ``SIGMA_RUN_DEX``
      (omega, the model-discrepancy term, is NOT in the band: it would make a 5x miss look tolerable);
  (b) the in-sample air arm is not moved outside its own band: the median air PA + PAA within +-2 s_air of the
      measured air total, s_air built the same way.
A partial reduction is reported as the fraction of the log gap closed: (log m_base - log m_hyp) / (log m_base -
log 0.09) per rival, m the median over that rival's feasible draws (the same Sobol draws in every scenario).

Run (wrapper, threads pinned, niced):
    MAILLARD_M3_LAYER=<path to m3_strecker_layer.py> OMP_NUM_THREADS=1 nice -n 10 \
        ./scripts/docker_maillard.sh run "python scripts/diagnostics/m3_argon_strecker_diagnostic.py"
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

OUT_STEM = data_paths.VALIDATION_DIR / "m3_argon_strecker_diagnostic"
RIVALS = ("R1", "R2", "R3", "R4")
N_DRAWS = 1024
SEED = 3                      # the checks record's v_time seed
MIN_FEASIBLE_SHARE = 0.05
BAND_Z = 2.0

# ============================================================================== constants, each with its source
#: first-order losses (/min, 100 C) of the ARP-made donors. "fed" = a fed-donor decay; "fit" = a multiresponse
#: literature fit (a measured system, not a direct decay); the layer's own k0_DG1 band is Davidek 2006a.
K0 = {
    "TDG_ghiron_0_60": (0.0067, "Ghiron 1988 Fig. 2: fed 3-DG alone, 0.4 M, water (unbuffered), 100 C, 0-1 h; "
                                "first-order [Z] (data/lit/extraction_dossiers/ghiron1988_extraction.md sec. 3)", "fed"),
    "TDG_ghiron_60_120": (0.0036, "Ghiron 1988 Fig. 2, 1-2 h (the dossier's 'safer' window)", "fed"),
    "TDG_martins2005": (0.0345, "Martins & van Boekel 2005 Table 2 step 5, 3-DG -> formic acid, 100 C, pH 6.8 "
                                "(glucose / glycine; martins2005_extraction.md sec. 1)", "fit"),
    "TDG_martins2003": (0.0907 + 0.0274, "Martins & van Boekel 2003 Table 3 col. C (fed DFG 10 mM, 0.1 M phosphate, "
                                         "100 C, pH 6.8): k5 3-DG -> C_n 0.0907 + k6 3-DG -> FA 0.0274 "
                                         "(martins2003_extraction.md sec. 5.2)", "fit"),
    "MGO_kanzler_0_120": (0.0098, "Kanzler 2026: fed MGO alone 100 mM, pH 5, 100 C, first-order reading 0-120 min "
                                  "(kanzler2026_extraction.md sec. 4.7 item 1; the data are second order)", "fed"),
    "MGO_kanzler_0_30": (0.0144, "Kanzler 2026, first-order reading 0-30 min (fastest window)", "fed"),
    "MGO_second_order": (0.0, "Kanzler 2026 second-order reading k2 1.85e-4 /mM/min: at the ARP pot's sub-mM to "
                              "few-mM MGO the first-order equivalent is <= ~1e-3 /min; taken as 0", "fed"),
    "DG1_martins2005": (1.45, "Martins 2005 Table 2 step 8, 1-DG -> acetic acid, 100 C, pH 6.8", "fit"),
    "DG1_martins2003": (1.9085, "Martins 2003 Table 3 col. C, k9 1-DG -> AA (k8 = 0)", "fit"),
}
#: Martins 2003 Table 3 col. C (fed DFG, 100 C, pH 6.8): the Amadori compound's anoxic exit split (S1, context)
M03 = {"k1": 0.0057, "k2": 0.0156, "k3": 0.0155, "k4": 0.0794, "k10": 0.0707, "k11": 0.1131, "k7": 0.2125,
       "k16": 0.0134}


def load_layer():
    try:
        from src import m3_strecker_layer as L  # noqa: F401  (after the M3 PR merges)
        return L, Path(L.__file__)
    except ImportError:
        pass
    path = os.environ.get("MAILLARD_M3_LAYER")
    if not path or not Path(path).is_file():
        raise SystemExit("src.m3_strecker_layer not importable and $MAILLARD_M3_LAYER unset: point it at the layer")
    spec = importlib.util.spec_from_file_location("m3_strecker_layer_by_path", path)
    L = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = L
    spec.loader.exec_module(L)
    return L, Path(path)


# ================================================================================================= noise bands

def noise_band(L, data, arm: str) -> Dict[str, float]:
    pa, paa = data.arp_t2[arm]
    epa, epaa = (math.hypot(e, L.S_REP_DEX) for e in data.arp_t2_err[arm])
    ln10 = math.log(10.0)
    sd_lin = math.hypot(pa * epa * ln10, paa * epaa * ln10)
    s_read = sd_lin / (pa + paa) / ln10
    cv = L.CV_ARP_BAND[1] / ln10
    s = math.sqrt(s_read ** 2 + cv ** 2 + L.SIGMA_RUN_DEX ** 2)
    tot = pa + paa
    return {"measured_mM": tot, "PA_mM": pa, "PAA_mM": paa, "s_read_dex": s_read, "cv_hi_dex": cv,
            "sigma_run_dex": L.SIGMA_RUN_DEX, "s_dex": s, "z": BAND_Z,
            "lo_mM": tot * 10 ** (-BAND_Z * s), "hi_mM": tot * 10 ** (BAND_Z * s)}


# ============================================================================ the fed calibration with a k0

def fed_courses_k0(L, kappa, nu, k0):
    """Fed donor with dD/dt = -kappa D^2 - k0 D (D0 = 100 mM); amine bookkeeping as EX-049 (nu per donor lost)."""
    g = L.GRID[None, :]
    if k0 <= 0:
        return L.fed_donor_courses(kappa, nu)
    e = np.exp(-k0 * g)
    d = L.D0 * k0 * e / (k0 + kappa[:, None] * L.D0 * (1.0 - e))
    aa = np.maximum(L.AA0 - nu[:, None] * (L.D0 - d), 0.0)
    return d, aa


def make_calibrate_k0(L, k0_fed: Dict[str, float]):
    """EX-049's route-A inversion (L.calibrate_route_a) with a first-order donor loss in the fed pots."""
    def calibrate_route_a(th, data):
        k_loss = 10.0 ** th["log10_k_loss"]
        i60 = L.idx(60.0)
        kap = L.kappas(th)
        cal = {"k": {}, "phi": {}, "omega": {}, "kappa": kap, "k_loss": k_loss}
        courses = {}
        for d in L.DONORS:
            dd, aa = fed_courses_k0(L, kap[d], th["nu"], k0_fed.get(d, 0.0))
            f = dd * aa
            integ, lossy = L.cum_trapz(f), L.cum_with_loss(f, k_loss)
            courses[d] = (integ, lossy)
            pa, paa = data.fed[d]["argon"]
            a, b = pa / lossy[:, i60], paa / integ[:, i60]
            cal["k"][d] = a + b
            cal["phi"][d] = b / (a + b)
        for cls, d in L.OMEGA_CALIBRATOR.items():
            cal["omega"][cls] = data.fed[d]["air"][1] / (cal["k"][d] * courses[d][0][:, i60]) - cal["phi"][d]
        integ, lossy = courses["MGO"]
        k, phi, om = cal["k"]["MGO"][:, None], cal["phi"]["MGO"][:, None], cal["omega"]["C23"][:, None]
        cal["t3_pa"] = (1 - phi) * k * lossy
        cal["t3_total"] = cal["t3_pa"] + (phi + om) * k * integ
        return cal
    return calibrate_route_a


@contextlib.contextmanager
def patched_calibration(L, k0_fed: Optional[Dict[str, float]]):
    if not k0_fed:
        yield
        return
    orig = L.calibrate_route_a
    L.calibrate_route_a = make_calibrate_k0(L, k0_fed)
    try:
        yield
    finally:
        L.calibrate_route_a = orig


# ================================================================================================= scenarios

def k(name: str) -> float:
    return K0[name][0]


def martins2003_split() -> Dict[str, float]:
    tot = M03["k1"] + M03["k2"] + M03["k3"]
    e1, e2, cn = M03["k1"] / tot, M03["k2"] / tot, M03["k3"] / tot
    to_tdg = M03["k4"] / (M03["k4"] + M03["k10"] + M03["k11"])
    to_dg1 = M03["k7"] / (M03["k7"] + M03["k16"])
    a3, a1, af = e1 * to_tdg, e2 * to_dg1, cn
    inert = 1.0 - a1 - a3 - af
    return {"a1_DG1": a1, "a3_TDG": a3, "aF_fragments": af, "inert_hexose_reversion": inert,
            "k_d_total_per_min": tot}


SCENARIOS: Dict[str, Dict[str, Any]] = {
    "base": {"what": "the layer as built (primary form)"},
    # H1: a first-order loss of the ARP-made donors
    "H1_fed": {"what": "H1, fed-donor decays: 3-DG Ghiron 0-1 h, MGO Kanzler first-order 0-120 min; 1-DG the layer's "
                       "own k0 band", "k0": {"TDG": k("TDG_ghiron_0_60"), "MGO": k("MGO_kanzler_0_120")}},
    "H1_fed_low": {"what": "H1 low bracket: 3-DG Ghiron 1-2 h, MGO second-order reading (0)",
                   "k0": {"TDG": k("TDG_ghiron_60_120"), "MGO": k("MGO_second_order")}},
    "H1_martins2005": {"what": "H1 with Martins 2005's 3-DG -> FA and 1-DG -> AA (fit), MGO Kanzler 0-120",
                       "k0": {"TDG": k("TDG_martins2005"), "MGO": k("MGO_kanzler_0_120"),
                              "DG1": k("DG1_martins2005")}},
    "H1_max": {"what": "H1 upper bracket: every largest sourced value (3-DG Martins 2003 k5 + k6; MGO Kanzler 0-30; "
                       "1-DG Martins 2003 k9), ARP-made donors only",
               "k0": {"TDG": k("TDG_martins2003"), "MGO": k("MGO_kanzler_0_30"), "DG1": k("DG1_martins2003")}},
    "H1_max_consistent": {"what": "H1_max with the same 3-DG / MGO losses ALSO in the fed-donor calibration "
                                  "(Table 4 re-inverted)",
                          "k0": {"TDG": k("TDG_martins2003"), "MGO": k("MGO_kanzler_0_30"),
                                 "DG1": k("DG1_martins2003")},
                          "k0_fed": {"TDG": k("TDG_martins2003"), "MGO": k("MGO_kanzler_0_30")}},
    # H2: donor efficiency
    "H2a_dg1_as_3dg": {"what": "H2: the 1-DG donor carries 3-DG's calibrated Strecker share and phi (1-DG has no "
                               "Table 4 arm)", "eff": {"DG1": "TDG"}},
    "H2b_anoxic_as_3dg": {"what": "H2: every anoxic-exit donor (1-DG and fragment MGO) carries 3-DG's share and phi",
                          "eff": {"DG1": "TDG", "MGO": "TDG"}},
    "H2_low_as_glucosone": {"what": "H2 low bracket: 1-DG and fragment MGO carry glucosone's share and phi (the "
                                    "lowest measured)", "eff": {"DG1": "G", "MGO": "G"}},
    "H12_max": {"what": "H1_max + H2b: the most favourable combination of sourced constants",
                "k0": {"TDG": k("TDG_martins2003"), "MGO": k("MGO_kanzler_0_30"), "DG1": k("DG1_martins2003")},
                "eff": {"DG1": "TDG", "MGO": "TDG"}},
    # H3: an oxygen-independent PA sink (no measured PA loss on disk: the layer's own k_loss band top)
    "H3_kloss_top": {"what": "H3: PA first-order loss at the top of EX-049's k_loss band (1e-2 /min; in-sample band, "
                             "not a measurement), re-inverted consistently", "log10_k_loss": -2.0},
    # S1 (outside H1 - H3; context for the probability statement)
    "S1_martins2003_split": {"what": "S1: the ARP's anoxic exit split set to Martins 2003's fed-DFG split (100 C, "
                                     "pH 6.8), renormalised over the three exits the layer has", "split": True},
    "S1_plus_H1_max": {"what": "S1 + H1_max", "split": True,
                       "k0": {"TDG": k("TDG_martins2003"), "MGO": k("MGO_kanzler_0_30"),
                              "DG1": k("DG1_martins2003")}},
}


def run(L, data, rival: str, sc: Dict[str, Any], u: np.ndarray) -> Dict[str, np.ndarray]:
    th = L.from_unit(u, rival, L.PRIMARY)
    if "log10_k_loss" in sc:
        th["log10_k_loss"] = np.full(u.shape[0], sc["log10_k_loss"])
    with patched_calibration(L, sc.get("k0_fed")):
        p = L.derive(th, rival, L.PRIMARY, data)
        for x, v in (sc.get("k0") or {}).items():
            p.donor[x]["k0"] = np.full(p.n, float(v))
        for x, src in (sc.get("eff") or {}).items():
            s = p.donor[src]["k"] / p.donor[src]["kappa"]
            p.donor[x]["k"] = s * p.donor[x]["kappa"]
            p.donor[x]["phi"] = p.donor[src]["phi"].copy()
        if sc.get("split"):
            sp = martins2003_split()
            tot = sp["a1_DG1"] + sp["a3_TDG"] + sp["aF_fragments"]
            p.a = tuple(np.full(p.n, sp[key] / tot) for key in ("a1_DG1", "a3_TDG", "aF_fragments"))
        p = L.calibrate(p, data)
    sig = np.where(np.isfinite(p.sigma_ox), p.sigma_ox, 0.0)
    cu = np.where(np.isfinite(p.c_cu5), p.c_cu5, 1.0)
    out = {"feasible": p.feasible, "sigma_ox": p.sigma_ox, "c_cu5": p.c_cu5}
    for atm in ("argon", "air", "air_cu5"):
        r = L.solve_arp_pot(p, L.POT_T2[atm], sig, cu)
        out[f"{atm}_PA"], out[f"{atm}_PAA"] = r["PA"][:, 0], r["PAA"][:, 0]
        out[f"{atm}_tot"] = r["PA"][:, 0] + r["PAA"][:, 0]
        if atm == "argon":
            out["argon_decayed"] = 1.0 - r["ARP"][:, 0] / p.th["arp0"]
    if sc is SCENARIOS["base"]:   # per-donor attribution of the argon total (each donor's k zeroed in turn)
        for x in ("DG1", "TDG", "MGO"):
            p2 = copy.copy(p)
            p2.donor = {kk: dict(v) for kk, v in p.donor.items()}
            p2.donor[x]["k"] = np.zeros(p.n)
            r = L.solve_arp_pot(p2, L.POT_T2["argon"], sig, cu)
            out[f"attr_{x}"] = out["argon_tot"] - (r["PA"][:, 0] + r["PAA"][:, 0])
    return out


def q(x) -> list:
    x = np.asarray(x)[np.isfinite(x)]
    return [round(float(v), 4) for v in np.percentile(x, [10, 50, 90])] if x.size else [None] * 3


def check_bit_exact(L, data) -> bool:
    th = L.from_unit(L.sobol_unit(64, len(L.coords("R1")), 7), "R1")
    a, b = L.calibrate_route_a(th, data), make_calibrate_k0(L, {})(th, data)
    return all(np.array_equal(a[key][d], b[key][d]) for key in ("k", "phi") for d in L.DONORS) and \
        all(np.array_equal(a["omega"][c], b["omega"][c]) for c in ("C23", "C6"))


def memory_free_pct() -> Optional[int]:
    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, timeout=60).stdout
        return int(out.strip().splitlines()[-1].split(":")[-1].strip().rstrip("%"))
    except Exception:
        return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--draws", type=int, default=N_DRAWS)
    ap.add_argument("--only", nargs="*", default=None, help="scenario names (default: all)")
    args = ap.parse_args(argv)
    free = memory_free_pct()
    if free is not None and free < 30:
        raise SystemExit(f"memory free {free} % < 30 %: not running")
    L, layer_path = load_layer()
    data = L.load_data()
    t0 = time.time()
    bands = {arm: noise_band(L, data, arm) for arm in ("argon", "air", "air_cu5")}
    names = args.only or list(SCENARIOS)
    res: Dict[str, Dict[str, Any]] = {}
    base_med: Dict[str, float] = {}
    for name in ["base"] + [n for n in names if n != "base"]:
        sc = SCENARIOS[name]
        res[name] = {"what": sc["what"], "rivals": {}}
        for rival in RIVALS:
            u = L.sobol_unit(args.draws, len(L.coords(rival)), SEED)
            o = run(L, data, rival, sc, u)
            f = o["feasible"]
            r: Dict[str, Any] = {"feasible_share": round(float(f.mean()), 4), "n_feasible": int(f.sum())}
            for atm in ("argon", "air", "air_cu5"):
                r[f"{atm}_tot_q10_50_90_mM"] = q(o[f"{atm}_tot"][f])
            r["argon_PA_q10_50_90_mM"] = q(o["argon_PA"][f])
            r["argon_PAA_q10_50_90_mM"] = q(o["argon_PAA"][f])
            r["argon_decayed_q10_50_90"] = q(o["argon_decayed"][f])
            r["sigma_ox_q10_50_90"] = q(o["sigma_ox"][f])
            med = float(np.median(o["argon_tot"][f])) if f.any() else float("nan")
            if name == "base":
                base_med[rival] = med
                for x in ("DG1", "TDG", "MGO"):
                    r[f"argon_attribution_{x}_share_q10_50_90"] = q((o[f"attr_{x}"] / o["argon_tot"])[f])
                r["argon_within_2x_share"] = round(float(np.mean(np.abs(np.log10(o["argon_tot"][f] / 0.09))
                                                               <= math.log10(2))), 4)
            gap = math.log10(base_med[rival]) - math.log10(bands["argon"]["measured_mM"])
            r["log_gap_closed"] = round((math.log10(base_med[rival]) - math.log10(med)) / gap, 3)
            dev = abs(math.log10(med / bands["argon"]["measured_mM"]))
            r["argon_dev_dex"] = round(dev, 3)
            air_med = float(np.median(o["air_tot"][f])) if f.any() else float("nan")
            r["air_dev_dex"] = round(abs(math.log10(air_med / bands["air"]["measured_mM"])), 3)
            r["pass_argon"] = bool(dev <= BAND_Z * bands["argon"]["s_dex"])
            r["pass_air"] = bool(r["air_dev_dex"] <= BAND_Z * bands["air"]["s_dex"])
            res[name]["rivals"][rival] = r
        counted = [v for v in res[name]["rivals"].values() if v["feasible_share"] >= MIN_FEASIBLE_SHARE]
        res[name]["passes"] = bool(counted) and all(v["pass_argon"] and v["pass_air"] for v in counted)
        res[name]["log_gap_closed_range"] = [min(v["log_gap_closed"] for v in counted),
                                             max(v["log_gap_closed"] for v in counted)] if counted else None
        print(name, res[name]["passes"], res[name]["log_gap_closed_range"],
              {rv: v["argon_tot_q10_50_90_mM"] for rv, v in res[name]["rivals"].items()}, flush=True)
    record = {
        "title": "M3 diagnostic: why every rival over-predicts the argon Strecker aldehyde (seen data only)",
        "generated_by": "scripts/diagnostics/m3_argon_strecker_diagnostic.py",
        "status": "DIAGNOSTIC: no fit, no band, prereg or shipped output changed",
        "layer": {"path_basename": layer_path.name,
                  "sha256": hashlib.sha256(layer_path.read_bytes()).hexdigest(),
                  "branch": "m3-phase1-build (unmerged; loaded by path)"},
        "draws_per_rival": args.draws, "sobol_seed": SEED, "rivals": list(RIVALS), "form": "primary",
        "pass_rule": __doc__.split("PASS RULE")[1].split("Run (wrapper")[0].strip(),
        "noise_bands": bands, "constants": {kk: {"value": v[0], "source": v[1], "kind": v[2]} for kk, v in K0.items()},
        "martins2003_exit_split": martins2003_split(),
        "calibration_with_k0_bit_exact_at_zero": check_bit_exact(L, data),
        "scenarios": res, "wall_s": round(time.time() - t0, 1), "memory_free_pct_at_launch": free,
    }
    OUT_STEM.with_suffix(".json").write_text(json.dumps(record, indent=1) + "\n")
    print("wrote", OUT_STEM.with_suffix(".json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
