"""`maillard calibrate` on the real engine (calibration_prereg.md T2, T4, T5).

T2: measurements generated from the engine itself along a time series at 145 C with the FFT decay
rate shifted and a threefold response factor on FFT; the fit must recover the factor within 0.1 dex
and the shift within twice its posterior sigma, and the hold-out must improve. (A first version
shifted the norfuraneol-to-MFT rate: at 145 C that rate scales MFT almost uniformly, so the
contrasts could not separate it from the factor and the fit rightly left it there. A sink rate
changes the SHAPE of a time series, which is what a contrast can see.)
T5: Yiltirak 2026's four-temperature ladder as a pretend laboratory: 100 and 120 C fit, 110 and
130 C validate; the hold-out median fold error must improve. The temperature trend is not expected
to flip (a response factor cannot fix a structural miss), and the card says which it did.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import calibration as C
from src.kinetic_core import user_fit as F

ROOT = Path(__file__).resolve().parents[2]
MFT = "2-methyl-3-furanthiol"
FFT = "2-furfurylthiol"
COORD = "b8.k_fft_decay.log10_k_ref_145C"
TRUE_SHIFT = -0.2          # about two shipped sigmas; the prior pulls back, the posterior must still cover the truth
TRUE_FACTOR = 3.0
SHIFTED = FFT               # the compound the shifted coordinate acts on, and the one that gets the factor

pytestmark = pytest.mark.slow


def _truth_calibration() -> C.Calibration:
    cands = {c.name: (c, v, s, b) for c, v, s, b in C.candidate_coordinates("sulfur")}
    coord, value, sigma, band = cands[COORD]
    return C.Calibration("truth", "b9", "2026-09-08", "water", {},
                         (C.Override(coord, value, sigma, value + TRUE_SHIFT, sigma, band),), (), ())


def _synthetic_document():
    """Six Hofmann-like pots along a time series at 145 C (a sink rate shows in the shape)."""
    truth = _truth_calibration()
    systems = []
    grid = [(145.0, 5.0), (145.0, 10.0), (145.0, 20.0), (145.0, 40.0), (145.0, 60.0), (145.0, 90.0)]
    for i, (T, t) in enumerate(grid):
        spec = {"name": f"pot_{i}", "precursors": {"L-Cysteine": 33.0, "D-Ribose": 100.0}, "temp_C": T, "time_min": t,
                "ph": 5.0, "aw": 0.98, "matrix": "water", "buffer": {"kind": "phosphate", "phosphate_mol_l": 0.5}}
        p = F.engine_predict(spec, [MFT, FFT], truth)
        assert p.values.get(MFT) and p.values.get(FFT), (i, p)
        systems.append({**spec, "measured": {MFT: {"value": p.values[MFT], "uncertainty_pct": 10},
                                             FFT: {"value": p.values[FFT] * TRUE_FACTOR, "uncertainty_pct": 10}},
                        "quantification_class": "stable_isotope_dilution_gcms", "source": {"lab": "synthetic"}})
    return {"systems": systems}


@pytest.fixture(scope="module")
def synthetic():
    """The unrestricted run: the contrasts choose among every calibratable sulfur coordinate."""
    return F.calibrate(_synthetic_document(), "synthetic", max_coordinates=2, max_nfev=25)


@pytest.fixture(scope="module")
def synthetic_named():
    """The run a laboratory makes when it knows which step its pot differs in."""
    return F.calibrate(_synthetic_document(), "synthetic", max_coordinates=2, max_nfev=25, coordinates=[COORD])


def test_t2_the_factor_is_recovered_on_the_engine(synthetic):
    cal, card = synthetic
    assert abs(cal.response_factors[FFT].log10 - math.log10(TRUE_FACTOR)) < 0.1, (card["response_factors"], card["diagnostics"], card["overrides"])
    assert abs(cal.response_factors[MFT].log10) < 0.1, (card["response_factors"], card["diagnostics"], card["overrides"])


def test_t2_the_shifted_coordinate_is_recovered_when_named(synthetic_named):
    cal, card = synthetic_named
    moved = {o.coordinate.name: o for o in cal.overrides}
    assert COORD in moved, card["diagnostics"]
    o = moved[COORD]
    assert abs(o.shift - TRUE_SHIFT) < max(2 * o.sigma, 0.1), (o.shift, o.sigma, card["diagnostics"])
    assert abs(cal.response_factors[FFT].log10 - math.log10(TRUE_FACTOR)) < 0.1, card["response_factors"]


def test_t2_the_unrestricted_run_reports_its_aliasing(synthetic):
    """One time series at one temperature cannot tell a sink rate from the formation and osone rates:
    the unrestricted run picks the best-explaining set and the card lists what it did not choose."""
    cal, card = synthetic
    d = card["diagnostics"]
    assert d["identified"] and len(d["identified"]) <= 2
    assert COORD in d["identified"] + d["not_identified"]
    assert cal.overrides


def test_t2_the_holdout_improves_on_the_engine(synthetic):
    _, card = synthetic
    before, after = card["holdout"]["before"], card["holdout"]["after"]
    assert after["median_fold"] < before["median_fold"], (before["median_fold"], after["median_fold"])


def test_t3_no_validate_record_was_read_during_the_fit(synthetic):
    cal, card = synthetic
    assert cal.validate_records
    assert not set(card["reads_during_fit"]) & set(cal.validate_records)


def _yiltirak_document():
    """The four buffer pots of Yiltirak 2026 as one laboratory's records, from the tracked bundles."""
    base = ROOT / "data" / "benchmarks" / "external_validation" / "maillard_path"
    files = sorted(base.glob("mp_holdout_ribose_cysteine_buffer_*_Yiltirak2026.json"))
    assert len(files) == 4, files
    systems = []
    for f in files:
        b = json.loads(f.read_text(encoding="utf-8"))
        c = b["conditions"]
        T = float(c["temp_C"])
        measured = {}
        for name, rec in b["reference_volatiles"].items():
            if "furanthiol" in name.lower() or "furfurylthiol" in name.lower():
                key = MFT if "furanthiol" in name.lower() else FFT
                measured[key] = {"value": float(rec["conc_ppb"]), "uncertainty_pct": float(rec.get("uncertainty_pct") or 15)}
        assert measured, f.name
        systems.append({
            "name": f"yiltirak_{int(T)}C", "precursors": {k: float(v["concentration_mM"]) for k, v in b["precursors"].items()},
            "temp_C": T, "time_min": float(c["time_min"]), "ph": float(c["ph"]), "aw": float(c.get("water_activity") or 0.98),
            "matrix": "water", "buffer": {"kind": "phosphate", "phosphate_mol_l": 0.5},
            "role": "fit" if int(T) in (100, 120) else "validate",
            "measured": measured, "quantification_class": "stable_isotope_dilution_gcms", "source": {"lab": "Reading (pretend)"},
        })
    return {"systems": systems}


@pytest.fixture(scope="module")
def yiltirak():
    return F.calibrate(_yiltirak_document(), "Reading (pretend)", max_coordinates=2, max_nfev=25)


def test_t5_a_real_ladder_as_a_pretend_laboratory(yiltirak):
    cal, card = yiltirak
    assert set(cal.validate_records) == {"yiltirak_110C", "yiltirak_130C"}
    before, after = card["holdout"]["before"], card["holdout"]["after"]
    assert before and after and after["n"] == before["n"] > 0
    assert after["median_fold"] < before["median_fold"], (before["median_fold"], after["median_fold"])
    # the card must state a factor for both thiols: the pretend laboratory's levels sit far below the model's
    assert cal.response_factors[MFT].log10 < 0 and cal.response_factors[FFT].log10 < 0


def test_the_cli_template_calibrates_unvalidated(tmp_path):
    """One record: the verb runs, writes a calibration, and says the calibration is unvalidated."""
    import subprocess
    import sys

    template = subprocess.run([sys.executable, "scripts/maillard.py", "score", "--template"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    doc = tmp_path / "one.yml"
    doc.write_text(template, encoding="utf-8")
    out = subprocess.run([sys.executable, "scripts/maillard.py", "calibrate", str(doc), "--lab", "smoke", "--out", str(tmp_path)],
                         cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "UNVALIDATED" in out.stdout
    written = list((tmp_path / "smoke").glob("calibration_*.json"))
    assert written, out.stderr
    cal = C.Calibration.load([p for p in written if "card" not in p.name][0])
    assert cal.lab == "smoke" and cal.response_factors
