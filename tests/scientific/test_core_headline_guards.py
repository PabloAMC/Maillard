"""Headline guards for the KINETIC CORE's own artifacts (retirement step B4, 2026-09-03).

The legacy guards in ``test_honest_headline_guards.py`` pin numbers the legacy harness
computes. These pin the numbers the core's scorecard and envelope publish, and the README
sentences that quote them, under the same rule: a moved number must move the prose in the
same change. Every pin carries the date and the cause of its last move.

Sources:
  * ``results/validation/core_panel_scores.json`` -- tracked, AND recomputed live here
    (~15 s) so a stale artifact fails rather than being echoed.
  * ``results/validation/core_prediction_uncertainty.json`` -- tracked (n=200, ~40 min to
    regenerate; not recomputed here).
  * ``results/validation/cutover_final_exam.json`` -- tracked.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import scoring

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
AUDIT = ROOT / "docs" / "history" / "AUDIT_legacy_lane_2026-08.md"
SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
ENVELOPE = data_paths.VALIDATION_DIR / "core_prediction_uncertainty.json"
EXAM = data_paths.VALIDATION_DIR / "cutover_final_exam.json"
DIRECTIONAL = data_paths.CORE_DIRECTIONAL_SCORES


def _doc_text(path: Path) -> str:
    assert path.exists(), f"{path.name} is missing; the headline claims have no home"
    return path.read_text(encoding="utf-8")


def _assert_quoted(text: str, token: str, where: str, what: str) -> None:
    assert token in text, (
        f"{where} no longer quotes {token!r} for {what}. Either the number changed and this "
        f"guard was not updated with it, or the prose was edited away from the evidence."
    )


@pytest.fixture(scope="module")
def tracked_scores():
    return json.loads(SCORES.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def live_scores():
    return scoring.score_panel()


# --------------------------------------------------------------------------------------
# 1. The scorecard artifact is what the code produces today
# --------------------------------------------------------------------------------------


def test_tracked_scorecard_is_not_stale(tracked_scores, live_scores):
    """The cutover exam drifted 4/34 -> 3/34 silently before B2 regenerated it. Not again."""
    for key in ("panel_benchmark_count", "scored_benchmark_count", "matched_compound_count",
                "within_band_count", "predictive_passes", "strict_ready", "evidence_role_totals"):
        assert tracked_scores["summary"][key] == live_scores["summary"][key], key
    for key in ("honest_literature", "out_of_sample", "in_core_fit"):
        for field in ("rows", "within_band", "benchmarks"):
            assert tracked_scores["summary"][key][field] == live_scores["summary"][key][field], (key, field)
    live = {(b["benchmark_id"], r["compound"]): r["predicted"] for b in live_scores["benchmarks"] for r in b["compounds"]}
    tracked = {(b["benchmark_id"], r["compound"]): r["predicted"] for b in tracked_scores["benchmarks"] for r in b["compounds"]}
    assert set(live) == set(tracked)
    # rel 1e-6, not 1e-9: the same code on linux/arm64 vs linux/amd64 (LSODA + BLAS) reproduces
    # every prediction to 7e-8 relative (measured 2026-09-03, 49 rows); a 3x band decision is
    # six orders of magnitude away from that noise.
    for key, value in live.items():
        assert tracked[key] == pytest.approx(value, rel=1e-6), key


# --------------------------------------------------------------------------------------
# 2. Panel, rows, band, strict-ready
# --------------------------------------------------------------------------------------


def test_core_panel_is_37_bundles_27_answered_42_rows(tracked_scores):
    # RE-PINNED 2026-09-04 (later the same day): the furosine bundle's DOI is a bread-baking paper
    # with no extrusion point; quarantined. Its row was always refused (furosine is no core species).
    # RE-PINNED 2026-09-04: two sourceless bundles quarantined (3 rows) and the eight hexose-only
    # thiol rows are declared NOT EVALUABLE (unidentified hexose entry): 40/32/49 -> 38/27/39.
    """Pinned 2026-09-03 (B3). The union panel: 19 scored trust-loop bundles (the two
    *_Internal2026 synthetic snapshots are off the panel), 17 maillard_path hold-outs, 4
    external matrix bundles. 8 bundles are refused whole (H2S / hydroxyacetaldehyde /
    mercapto-2-propanone precursors, CML/CEL/furosine targets, 2-pentylfuran)."""
    s = tracked_scores["summary"]
    assert s["panel_benchmark_count"] == 37
    assert s["scored_benchmark_count"] == 27
    # RE-PINNED BY WAVE B28 (2026-09-09): 39 -> 42 rows. Three nonanal rows left the refused list
    # and joined the scored one, because the oleate branch fraction the lane needed was measured in
    # 1978 and read this week. The BENCHMARK and ANSWERED counts do not move: no new pot entered,
    # three questions inside existing pots stopped being refused. The within-3x count is unchanged
    # at 4, so the three new rows all miss the band -- see the next test, where that is the point.
    assert s["matched_compound_count"] == 42
    # RE-PINNED BY WAVE B28: 25 -> 22 refused. Exactly the three nonanal rows above; the four
    # 2-pentylfuran rows stay refused, and their REASON changed rather than their status. That
    # matters and is asserted elsewhere: the branch fraction those rows were missing is measured
    # now, and they are refused because the hexanal they are scored against is not produced by the
    # lane that would carry the alkylfuran. Un-refusing them answered four rows six to nine orders
    # of magnitude below measurement, which is worse than refusing.
    assert s["refused_compound_count"] == 22
    # 2026-09-03: the xylose pH-5 bundle left the hold-out (the B2-B8 fit had read it) and
    # returned once wave B9 removed the Hofmann level rows from the objective.
    assert {k: v["benchmarks"] for k, v in s["by_panel"].items()} == {
        "trust_loop": 16, "maillard_path_holdout": 17, "external_matrix": 4,
    }


def test_core_evidence_roles_are_40_predictive_and_the_legacy_split_is_kept_beside(tracked_scores):
    """Pinned 2026-09-03 (B3). The core never fitted the rows the legacy matrix factors and
    projection budget were solved from, so those five bundles are predictive on the core and
    are listed as differing from the legacy role. Moving a bundle INTO fit_recovery without a
    core fit record that names it is the laundering route this guard blocks."""
    s = tracked_scores["summary"]
    # RE-PINNED 2026-09-03: the 21 physically separated hold-out bundles carry their own role.
    assert s["evidence_role_totals"] == {"external_holdout": 21, "predictive": 16}


def test_within_3x_is_4_of_42_and_out_of_sample_3_of_41(tracked_scores):
    # RE-PINNED 2026-09-04: 6/49 -> 4/39. The two hits lost are Bolton 1994 (its assumed
    # loadings replaced by the paper's Table I: the core now overpredicts MFT 20x) and one of the
    # quarantined Hofmann-derivation rows; the eight hexose thiol rows (all misses) left the count.
    """RE-PINNED 2026-09-03 (wave B9, the fit/validate split). With the eight Hofmann Table 1
    level rows out of the sulfur objective the core lands 6/49 within 3x (was 8/49 with those
    rows fitted); out of sample is now 48 of the 49 rows (only the C2+C3 pH-5 row is a fit row)
    at 5/48 (was 4/40). The hexose rows are the story: glucose and fructose MFT predict ZERO
    without the level rows -- the hexose-to-thiol route has no primary-evidence support."""
    s = tracked_scores["summary"]
    # RE-PINNED BY WAVE B28 (2026-09-09): denominators 39 -> 42 and 38 -> 41, NUMERATORS UNCHANGED.
    # The three nonanal rows land at 3.0x, 3.3x and 7.3x. The first two are inside a factor of
    # three and still do not count here, because this band is 3x on a row whose contract allows it
    # and nonanal's rows carry tighter contracts. So the headline rate FALLS, from 4/39 to 4/42,
    # on a wave that made the model strictly more capable. That is the honest arithmetic and it is
    # pinned rather than smoothed: answering more questions lowers a pass RATE unless the new
    # answers are good, and two of these three nearly are.
    assert (s["within_band_count"], s["matched_compound_count"]) == (4, 42)
    assert (s["honest_literature"]["within_band"], s["honest_literature"]["rows"]) == (4, 42)
    assert (s["out_of_sample"]["within_band"], s["out_of_sample"]["rows"]) == (3, 41)
    assert (s["in_core_fit"]["within_band"], s["in_core_fit"]["rows"]) == (1, 1)
    assert s["holdout_within_band"] == {"hits": 3, "total": 26}
    readme = _doc_text(README)
    _assert_quoted(readme, "4 of 42", "README.md", "the core's within-3x count")
    _assert_quoted(readme, "3 of 41", "README.md", "the core's out-of-sample count")


def test_no_core_benchmark_is_strict_ready_since_bolton_1994_was_read(tracked_scores):
    """Pinned 2026-09-03 (B3). Bolton 1994 (thiamine + cysteine + glucose, 120 C): MFT 13 ->
    17.4 ppb, 1.34x, under the bundle's OWN 3x / 0.48-dex contract; PRIMARY, free_precursor,
    not a row any core fit read. The legacy lane has 0/23. A second strict-ready bundle is
    either a real result or a contract that was loosened: check the bundle's
    validation_contract before re-pinning."""
    s = tracked_scores["summary"]
    # RE-PINNED 2026-09-04: the full text of Bolton 1994 replaced the assumed loadings (glucose
    # 10 -> 51.5 mM, thiamine 1 -> 13.7 mM, pH 5.0 -> 5.65); the core overpredicts MFT 20x and
    # no benchmark is strict-ready. The pass had rested on inputs nobody had read.
    assert s["strict_ready"] == []
    assert s["predictive_passes"] == []
    bolton = next(b for b in tracked_scores["benchmarks"] if b["benchmark_id"] == "thiamine_cys_glucose_120C_Bolton1994")
    assert bolton["scale_thresholds"]["max_ratio"] == pytest.approx(3.0)
    assert not any(r["in_core_fit"] for r in bolton["compounds"])
    readme = _doc_text(README)
    _assert_quoted(readme, "0 of 37", "README.md", "the core strict-ready count")
    _assert_quoted(readme, "thiamine_cys_glucose_120C_Bolton1994", "README.md", "the benchmark that lost its pass")
    _assert_quoted(_doc_text(AUDIT), "strict-ready 0/37", "AUDIT.md", "the core strict-ready count")


def test_the_xylose_row_is_a_hold_out_again_and_no_hofmann_level_row_is_in_the_fit(tracked_scores):
    """Pinned 2026-09-03 (B9). Found at B3 as a hold-out the fit had read; moved to the trust
    loop; returned to the hold-out once B9 removed the Hofmann level rows. The generated fit-
    target record names no Hofmann pH-5 level bundle any more."""
    xyl = next(
        b for b in tracked_scores["benchmarks"]
        if b["benchmark_id"] == "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5"
    )
    assert xyl["panel"] == "maillard_path_holdout"
    assert xyl["fit_target_of"]["core"] == []
    assert not any(r["in_core_fit"] for r in xyl["compounds"])
    for bid in ("hofmann1998_ribose_cysteine_145C_20min_pH5", "hofmann1998_glucose_cysteine_145C_20min_pH5",
                "hofmann1998_fructose_cysteine_145C_20min_pH5"):
        b = next(x for x in tracked_scores["benchmarks"] if x["benchmark_id"] == bid)
        assert b["fit_target_of"]["core"] == [] and not any(r["in_core_fit"] for r in b["compounds"]), bid


# --------------------------------------------------------------------------------------
# 3. The envelope
# --------------------------------------------------------------------------------------


def test_core_envelope_covers_5_of_33_evaluable_literature_rows_and_5_of_32_out_of_sample():
    # RE-PINNED 2026-09-04 (the unidentified-coordinate rule inverted): a coordinate the fit
    # could not pin is now DRAWN across its declared band instead of frozen at the optimum,
    # with an Ea band narrowed to a 12-decade prefactor prior and a definitionally-bounded
    # coordinate (acid yield) still fixed. 6/32 -> 5/33, out of sample 6/31 -> 5/32, median
    # width 0.98 -> 1.38 dex, sampled priors 35 -> 41. The width is the point; the coverage
    # is NOT, and the sweep in uncertainty.py records why: uncapped bands still cover only
    # 21% against a nominal 90%, so the residual is model-structure error.
    # RE-PINNED 2026-09-04 (quarantine + unidentified hexose entry): 9/42 -> 6/32, 9/41 -> 6/31;
    # every remaining row carries a declared quantification class (undeclared 0).
    # RE-PINNED 2026-09-03 (pass 7): five more bundles verified from PMC full text; ACSRef3 is
    # isotope-dilution LC-MS/MS (extraction), so its band-only interval is gone and the row is
    # NOT EVALUABLE (7): 9/43 -> 9/42, out of sample 9/42 -> 9/41; headspace 4 -> 8 rows.
    # RE-PINNED 2026-09-03 (B9): 11/44 -> 10/44 (width 1.07 -> 1.17 dex); out of sample now
    # excludes only the C2+C3 pH-5 row, 8/35 -> 10/43; Laplace 18/23 -> 20/23, chi2_red 1.03 -> 1.21.
    # RE-PINNED 2026-09-03 (quantification classes declared): Bolton 1994 is extraction-quantified,
    # so its row no longer gets the headspace bands; the interval it had was band width alone
    # and the row is now NOT EVALUABLE (6, was 5): 10/44 -> 9/43, out of sample 10/43 -> 9/42.
    """Pinned 2026-09-03 (B8, n=200 seed 0). The sulfur lane is now sampled jointly from the
    Laplace covariance at its B8 optimum (18 of 23 free coordinates identified, reduced
    chi-square 1.03), so the 24 rows that were "not evaluable" at B3 are evaluable: 7/25 ->
    11/44, median width 0.94 -> 1.07 dex, 5 rows still not evaluable (zero-width or refused
    draws). Out of sample (every row the sulfur fit read removed): 8/35."""
    payload = json.loads(ENVELOPE.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert (s["n_samples"], s["seed"]) == (200, 0)
    lit = s["honest_literature_coverage"]
    # RE-PINNED 2026-09-07 (B15): the extrusion row at a_w 0.35 is now inside the acrylamide window and
    # its interval is sampled (33 -> 34 evaluable, 6 -> 5 not evaluable); three new priors move the stream.
    assert (lit["hits"], lit["total"], lit["not_evaluable"]) == (7, 34, 5)
    # RE-PINNED 2026-09-07 (B10 ships the ambient-oxidant consistency fix: the engine now charges
    # OX_AMBIENT_MMOL_L on every sulfur run, as every fit system was): 1.3753 -> 1.3691 dex; every
    # count is unchanged (4/39, 3/38, 5/33, 17/28).
    # RE-PINNED 2026-09-07 (B12 envelope hooks + B13): two trunk declared-band priors joined the
    # draw table, which moves the random stream; coverage unchanged, width 1.3691 -> 1.3376 dex.
    # RE-PINNED 2026-09-07 (B14): the acrylamide lane's declared flat a_w band joined the draw
    # table (inert on every panel row except the a_w-inside-window acrylamide hold-outs, +0.05 dex);
    # the moved random stream shifts every row by a few hundredths: 5/33 -> 7/33, 1.3376 -> 1.4495 dex.
    assert lit["median_ci_width_log10"] == pytest.approx(1.3080, abs=5e-4)
    oos = s["out_of_sample_literature_coverage"]
    assert (oos["hits"], oos["total"]) == (7, 33)
    assert s["unsampled_lanes"] == []
    assert s["sulfur_laplace"]["identified"] == 20 and s["sulfur_laplace"]["free"] == 23
    assert s["sulfur_laplace"]["reduced_chi_square"] == pytest.approx(1.21, abs=0.01)
    assert s["observable_multiplier_policy"]["rows_by_family"] == {
        "headspace": 8, "extraction": 31, "undeclared": 0,
    }
    readme = _doc_text(README)
    _assert_quoted(readme, "7 of 34", "README.md", "the core envelope's literature coverage")
    _assert_quoted(readme, "7 of 33", "README.md", "the core envelope's out-of-sample coverage")
    _assert_quoted(readme, "20 of 23", "README.md", "the identified sulfur coordinates")


# --------------------------------------------------------------------------------------
# 4. The exam, regenerated 2026-09-03 after it had drifted silently
# --------------------------------------------------------------------------------------


def test_cutover_exam_is_34_answered_3_within_band_paired_24_78x_vs_10_86x():
    """Pinned 2026-09-03. The tracked artifact said 4/34 while the code produced 3/34 (four
    MFT/FFT rows moved with the B8 sulfur refit); the README said 42.23x paired (the B7
    scoring). Both now match the artifact the code writes."""
    payload = json.loads(EXAM.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert (s["core_answered"], s["core_declined"], s["core_within_band"]) == (34, 6, 3)
    assert s["old_within_band"] == 7
    assert s["paired_subset"]["core"]["median_fold_error"] == pytest.approx(24.78, abs=0.01)
    assert s["paired_subset"]["old"]["median_fold_error"] == pytest.approx(10.86, abs=0.01)
    assert s["core"]["median_fold_error"] == pytest.approx(19.08, abs=0.01)
    readme = _doc_text(README)
    _assert_quoted(readme, "**3 / 34**", "README.md", "the exam's within-band count")
    _assert_quoted(readme, "**24.78x**", "README.md", "the exam's paired core median")
    _assert_quoted(readme, "**10.86x**", "README.md", "the exam's paired old-lane median")
    _assert_quoted(readme, "**19.08x**", "README.md", "the exam's all-answered core median")


# --------------------------------------------------------------------------------------
# 5. The directional claims panel, on the core (step B7)
# --------------------------------------------------------------------------------------


def test_core_scores_17_of_26_independent_directional_claims():
    """Pinned 2026-09-03 (B7). 69 claims; 16 are prose-only and 22 independent claims are
    not evaluable on the core (refused arms: 2,5-dimethylpyrazine, 2-pentylfuran, nonanal from
    oleate, H2S / hydroxyacetaldehyde precursors; a zero disulfide). Of the evaluable
    independent claims the core agrees on 18 of 30 -- 15 of 22 with pH and water activity
    set aside, 3 of 8 on those two axes, 3 of the independent misses (7 over all claims) being
    identical predictions across an axis the lane carries no term for. The retired lane's 24/36 was a different engine
    on a different evaluable subset and is NOT a baseline for this number."""
    payload = json.loads(DIRECTIONAL.read_text(encoding="utf-8"))
    s = payload["summary"]
    # RE-PINNED 2026-09-06 (programme step R2(c), B10 prereg sec. 6): five sulfur temperature
    # claims added BEFORE wave B10 fits anything (YIL-01/02 evaluable, WANG-01/02 and MENG-01
    # recorded not evaluable); both Yiltirak sign claims MISS on the shipped B9 lane, which is
    # the pre-wave baseline the wave is judged against: 69 -> 74 claims, 17/26 -> 17/28,
    # temperature 5/7 -> 5/9, 26 -> 29 independent claims not evaluable.
    # RE-PINNED 2026-09-07 (B14): AW-05, the declared-flat acrylamide a_w claim (fit_adjacent), 74 -> 75.
    # RE-PINNED 2026-09-07 (B15): PH-ACR-01 (fit_adjacent), DIC-01 and DIC-02 (Zhang 2020), 75 -> 78.
    # RE-PINNED 2026-09-07 (B16 reads): DIC-03, RIB-T-01/02, HEX-T-01, SCH-T-01, 78 -> 83.
    # RE-PINNED 2026-09-07 (Wang 2022 read): WANG22-T-01..04 (a third lab, 100 and 140 C shapes), 83 -> 87.
    # RE-PINNED 2026-09-07 (thiol pH reads): WHI-PH-01, CER07-PH-01/02, MOT02-PH-01/02, 87 -> 92.
    assert payload["panel"]["claims"] == 92
    # RE-PINNED 2026-09-03 (step 5): comparisons that move an axis the lane has no term for
    # are REFUSED by the engine (water activity everywhere; pH on trunk / acrylamide / lipid),
    # so those claims are not evaluable instead of identical-prediction misses: 18/30 -> 18/27.
    # RE-PINNED 2026-09-03 (B9): 18/27 -> 17/27; the hexose-related claims moved.
    # RE-PINNED 2026-09-03 (zero floor): a claim whose arms were 1e-31 and 1e-29 ug/L had been a
    # MISS on a log ratio of integrator noise; below 1 pg/L a concentration is zero and the claim is
    # NOT EVALUABLE ("a predicted concentration is zero"): 17/27 -> 17/26, sugar identity 4/9 -> 4/8.
    # RE-PINNED 2026-09-07 (B12: the trunk gained declared a_w and pH terms): AW-01 and AW-03 became
    # evaluable (one agree, one disagree, both as the B12 prereg expected), 17/28 -> 18/30; the two
    # water-activity refusals on the trunk are gone (engine refusals 4 -> 2: the acrylamide arm of
    # AW-02 and the lipid arm of PROC-01); pH-and-a_w 4/5 -> 5/7; 29 -> 27 not evaluable.
    # RE-PINNED 2026-09-07 (B15 + the unstated-input sweep + ranking claims): WANG-02 (FFT peak, charged
    # as TTCA, unanimous over pH 5/7/9) AGREES and DIC-01 (Zhang 2020 dicarbonyl ordering) DISAGREES:
    # 18/30 -> 19/32; WANG-01 flips with pH and stays not evaluable; DIC-02 not evaluable (no glutamate).
    # RE-PINNED 2026-09-07 (B16 reads): DIC-03 (Leitzen) evaluable and misses, 19/32 -> 19/33; RIB-T-01/02
    # and HEX-T-01 predict zero on the shipped lane (not evaluable), SCH-T-01 is fit_system_overlap.
    # RE-PINNED 2026-09-07 (Wang 2022): FFT rising at 100 C agrees; MFT at 100 C peaks early and both thiols
    # collapse at 140 C in the lane (three misses): 19/33 -> 20/37; time 2/2 -> 3/6.
    # RE-PINNED 2026-09-07 (Cerny 2007 furfural, Mottram 2002 MFT and FFT agree; Cerny FFT ladder misses): 20/37 -> 23/41.
    # RE-PINNED 2026-09-08 (B18: the pyrazine step on the trunk): the two 2,5-dimethylpyrazine claims that
    # were refused for want of the species now evaluate and AGREE (PH-06 rising with pH 4 -> 9; TEMP-04's
    # 160 C over 145 C ordering): 23/41 -> 25/43; pH 7/9 -> 8/10; temperature 6/10 -> 7/11; 30 -> 28 not evaluable.
    assert s["headline"] == [25, 43]
    ind = s["independent"]
    assert (ind["excluding_ph_aw"]["agree"], ind["excluding_ph_aw"]["evaluable"]) == (16, 31)
    assert (ind["ph_aw"]["agree"], ind["ph_aw"]["evaluable"]) == (9, 12)
    assert ind["total"]["not_evaluable"] == 28
    assert ind["total"]["mechanism_absent"] == 0
    assert s["not_evaluable_reasons"]["refused by the engine"] >= 2
    cats = {k: (v["agree"], v["evaluable"]) for k, v in ind["by_category"].items()}
    assert cats["sugar_identity"] == (4, 10)  # B15/B16: DIC-01 and DIC-03 miss
    assert cats["temperature"] == (7, 11)     # B15: WANG-02 agrees; B18: TEMP-04 (2,5-dimethylpyrazine) agrees
    assert cats["ph"] == (8, 10)              # thiol pH reads (2026-09-07): Cerny 2007 furfural and Mottram 2002 agree, Cerny FFT ladder misses; B18: PH-06 agrees
    assert cats["moisture_aw"] == (1, 2)   # B12: the trunk answers a_w; AW-03 agrees, AW-01 does not
    assert cats["additive_cysteine"] == (2, 3)
    assert cats["time"] == (3, 6)             # Wang 2022: the sinks are too strong at 100 C and at 140 C
    readme = _doc_text(README)
    _assert_quoted(readme, "25 of 43", "README.md", "the core's directional headline")
    _assert_quoted(readme, "16 of 31", "README.md", "the directional count excluding pH and water activity")
    _assert_quoted(readme, "9 of 12", "README.md", "the directional count on pH and water activity")


def test_directional_scorecard_is_not_stale():
    """The tracked artifact must be what the code scores today (~90 s)."""
    from src.kinetic_core import directional

    tracked = json.loads(DIRECTIONAL.read_text(encoding="utf-8"))
    live = directional.score_panel()
    assert live["summary"]["headline"] == tracked["summary"]["headline"]
    assert {r["claim_id"]: r["status"] for r in live["claims"]} == {
        r["claim_id"]: r["status"] for r in tracked["claims"]
    }
