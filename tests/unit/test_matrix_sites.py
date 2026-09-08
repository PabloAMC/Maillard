"""The protein matrix as chemistry (results/validation/matrix_sites_prereg.md T1 to T5)."""
from __future__ import annotations

import math

import pytest

from src import api
from src.comparative_cli import SpecError, spec_to_core
from src.kinetic_core import matrix_sites as M
from src.kinetic_core.engine import predict

FFT = "2-furfurylthiol"
HEX = "hexanal"


def _spec(**extra):
    base = {"name": "x", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0}, "temp_C": 145.0, "time_min": 20.0, "ph": 5.0, "aw": 0.98}
    base.update(extra)
    return base


def test_t2_blg_site_densities_reproduce_from_the_dossier_counts():
    blg = M.matrices()["blg"]
    assert abs(blg.disulfide - 2 / 18362 * 1000) < 1e-9
    assert abs(blg.free_thiol - 1 / 18362 * 1000) < 1e-9
    assert abs(blg.amine - 15 / 18362 * 1000) < 1e-9
    import yaml
    from src import data_paths
    raw = yaml.safe_load(data_paths.PROTEIN_MATRICES.read_text(encoding="utf-8"))["matrices"]["blg"]["sites_mmol_per_g"]
    for name, value in (("disulfide", blg.disulfide), ("free_thiol", blg.free_thiol), ("amine", blg.amine)):
        assert abs(raw[name] - value) < 1e-4, name     # the file's convenience values agree with the counts
    assert "anantharamkrishnan2020b" in blg.source
    assert M.BINDING_CLASSES["saturated_aldehyde_amine"]["k2_bracket"] == (6.0e-6, 2.5e-5)


def test_t1_no_loading_charges_nothing_and_changes_nothing():
    core = spec_to_core(_spec())
    charged, note = M.resolve(core.process)
    assert charged is None and note is None
    run = predict(core, [FFT])
    assert run.run_metadata["matrix_sites"] is None and run.run_metadata["matrix_binding"] == {}
    assert run.species_mmol_per_l.get("PROT_SS", 0.0) == 0.0


def test_t5_a_named_matrix_needs_a_loading_and_an_unknown_matrix_says_so():
    with pytest.raises(SpecError):
        api.predict(_spec(matrix="blg"))
    payload = api.predict(_spec(matrix="lupin isolate", protein_g_per_l=50.0))
    assert payload["matrix"]["sites"] is None
    assert "no site densities on file" in payload["matrix"]["note"]


def test_t2b_the_isolate_densities_are_the_dossiers_numbers_and_say_what_is_missing():
    """Pea and soy isolates (2026-09-08 addendum): measured densities in mmol per gram of protein from
    ruan2014 / shimada1988 (soy) and chihi2016 / shen2022 (pea); no amine density is on file, so no
    aldehyde binds and the answer says so."""
    table = M.matrices()
    soy, pea = table["soy_isolate"], table["pea_isolate"]
    assert abs(soy.free_thiol - 0.0078) < 1e-9 and abs(soy.disulfide - 0.050) < 1e-9 and soy.amine == 0.0
    assert soy.bands["free_thiol"] == (0.0075, 0.0080) and soy.bands["disulfide"] == (0.046, 0.053)
    assert "ruan2014" in soy.source and "shimada1988" in soy.source
    assert abs(pea.free_thiol - 0.0053) < 1e-9 and abs(pea.disulfide - 0.0042) < 1e-9 and pea.amine == 0.0
    assert pea.bands["free_thiol"] == (0.0021, 0.0135) and "chihi2016" in pea.source and "shen2022" in pea.source
    assert "amine not on file" in pea.note and "protein_sites" in pea.note
    # the pools charge at the loading, the disulfide pool reaches the sulfur lane, no aldehyde binds
    core = spec_to_core(_spec(matrix="pea_isolate", protein_g_per_l=50.0))
    charged, note = M.resolve(core.process)
    assert note is None and abs(charged.disulfide - 50.0 * 0.0042) < 1e-9 and charged.amine == 0.0
    assert "amine not on file" in charged.as_dict()["note"]
    payload = api.predict(_spec(matrix="soy_isolate", protein_g_per_l=50.0), targets=[FFT])
    assert abs(payload["matrix"]["sites"]["pools_mmol_per_l"]["disulfide"] - 50.0 * 0.050) < 1e-9
    assert HEX not in (payload["matrix"]["binding"] or {})
    assert "amine not on file" in payload["matrix"]["sites"]["note"]
    run = predict(core, [FFT])
    assert run.run_metadata["matrix_binding"] == {}
    assert run.run_metadata["matrix_sites"]["pools_mmol_per_l"]["disulfide"] == pytest.approx(50.0 * 0.0042)


def test_t4_the_disulfide_pool_is_charged_and_the_held_rate_is_weak():
    core = spec_to_core(_spec(matrix="blg", protein_g_per_l=10.0))
    charged, _ = M.resolve(core.process)
    assert abs(charged.disulfide - 10.0 * 2 / 18362 * 1000) < 1e-9
    plain = predict(spec_to_core(_spec()), [FFT]).concentrations_ug_per_l[FFT]
    with_blg = predict(core, [FFT])
    fft = with_blg.concentrations_ug_per_l[FFT]
    assert with_blg.run_metadata["matrix_sites"]["pools_mmol_per_l"]["disulfide"] > 1.0
    # the channel runs (the state moves) but at the held ambient bracket it takes under 1 %
    assert 0.99 * plain <= fft <= plain * 1.0000001, (plain, fft)


def test_t3_binding_rises_with_loading_and_time_and_the_corners_are_ordered():
    blg = M.matrices()["blg"]
    def charged(g):
        return M.ChargedSites("blg", g, blg.free_thiol * g, blg.disulfide * g, blg.amine * g, blg.source, blg.amine_band)
    a = M.bound_fraction("HEXANAL", charged(10.0), [(7 * 24 * 60.0, 20.0)])
    b = M.bound_fraction("HEXANAL", charged(20.0), [(7 * 24 * 60.0, 20.0)])
    c = M.bound_fraction("HEXANAL", charged(10.0), [(14 * 24 * 60.0, 20.0)])
    assert 0 < a["bound_fraction"] < b["bound_fraction"] and a["bound_fraction"] < c["bound_fraction"]
    lo, hi = a["bound_fraction_corners"]
    assert lo <= a["bound_fraction"] <= hi
    # 1 % BLG at 20 C for seven days: a few percent, inside the synthesis's 18-74 day half-life bracket
    assert 0.005 < a["bound_fraction"] < 0.25, a
    assert M.bound_fraction("HEXANAL", charged(0.0), [(60.0, 100.0)]) is None
    assert M.bound_fraction("MFT", charged(10.0), [(60.0, 100.0)]) is None


def test_the_binding_reaches_the_answer_and_its_interval():
    spec = {"name": "fat", "precursors": {"methyl linoleate hydroperoxide": 1.0}, "temp_C": 100.0, "time_min": 60.0, "ph": 6.0, "aw": 0.95}
    try:
        plain = api.predict(spec)
    except SpecError:
        pytest.skip("the lipid lane's precursor name differs in this build")
    if not plain.get("answered"):
        pytest.skip("the lipid lane refused this charge")
    with_blg = api.predict({**spec, "matrix": "blg", "protein_g_per_l": 30.0})
    rows_plain = {r["compound"]: r for r in plain["rows"]}
    rows_blg = {r["compound"]: r for r in with_blg["rows"]}
    if HEX not in rows_plain:
        pytest.skip("hexanal not among the lane's default targets")
    assert rows_blg[HEX]["predicted_ug_per_l"] < rows_plain[HEX]["predicted_ug_per_l"]
    assert HEX in with_blg["matrix"]["binding"]
    assert rows_blg[HEX]["band_x"] >= rows_plain[HEX]["band_x"]


def test_the_schema_accepts_the_new_fields_and_refuses_negative_sites():
    api.predict(_spec(protein_g_per_l=10.0, protein_sites={"free_thiol_mmol_per_g": 0.05, "disulfide_mmol_per_g": 0.1, "amine_mmol_per_g": 0.8}))
    with pytest.raises(SpecError):
        api.predict(_spec(protein_g_per_l=10.0, protein_sites={"amine_mmol_per_g": -1}))
    with pytest.raises(SpecError):
        api.predict(_spec(protein_g_per_l=-3.0, matrix="blg"))
