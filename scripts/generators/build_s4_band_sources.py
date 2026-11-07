#!/usr/bin/env python
"""
build_s4_band_sources.py -- S4 band sourcing for the 43 sensitive fitted-unidentified coordinates (plan sec. 7, S4).

Reads the incumbent evidence audit (``results/validation/incumbent_evidence_tiers.json``, PR #98: its
``summary.s4_sensitive_unidentified`` list and its rows) and writes ``results/validation/s4_band_sources.{json,md}``.

For each coordinate (log10 k_ref and Ea separately) it records the candidate band sources found on disk (dossiers,
ledger, registers, the B55 measured-barrier classes), each with its paper key, table / figure, conditions, laboratory,
relation to the engine step and reader status, and proposes a band under S4: "a band comes from measured bounds, a
same-class analogue, or the reaction-class distribution of measured constants; a search bound is never a band; with no
class, the output is refused". The five direct measurements the audit records as S4 sources (DM-01 / 02 / 03 / 04 / 08)
are included; the two outside the 43 (k_aa_frag, k_mgo_mel) as extra entries.

The arithmetic is done here, from the numbers recorded below: a measured bound at temperature T is carried to the
coordinate's reference temperature over the WHOLE Ea band of its class (the loose edge, since k and Ea are sampled
jointly), and its 95 % edge combines the reported noise, the within-lab between-paper spread, the assay spread and the
between-lab term in quadrature (S2 (b)'s terms; the between-lab term is the step's measured spread with >= 2
laboratories, else the lane's: tau_ABS 0.964 dex for the sulfur lane, H0; undefined for a lane without H0 or a pooled
class spread, so there the 95 % edge is undefined and the bound cannot make a band).

The band sources below were first read from the dossiers for #100 and re-read row by row by the S1 second reader
(2026-11-06); the second reader's corrections are applied (SECOND_READING_CORRECTIONS, the class ruling, the acrylamide
between-lab term, held-out rows never band sources). The shipped (k, Ea) pair is checked jointly against every bound.
Moves no constant, runs no engine integration (the engine's pH-factor function is called as a pure formula).

    python scripts/generators/build_s4_band_sources.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402

AUDIT = data_paths.VALIDATION_DIR / "incumbent_evidence_tiers.json"
OUT = data_paths.VALIDATION_DIR / "s4_band_sources.json"
R_KJ = 8.314462618e-3
LN10 = math.log(10.0)
Z95 = 1.96

CARB, SULF = "carbohydrate", "sulfur"

# ---------------------------------------------------------------------------------------------------------------------
# The candidate band sources. value: in the coordinate's own unit (linear), at temperature_c. relation: two_sided (a
# same-channel measurement) | upper (a measured superset loss, or a ceiling) | lower (a measured subset channel) |
# none (recorded, not a bound on this coordinate). terms: one-sigma dex terms (reported noise, within-lab
# between-paper spread, assay spread). usable: whether it enters the band, with the reason when not.
# ---------------------------------------------------------------------------------------------------------------------
SOURCES: Dict[str, Dict[str, Any]] = {
    # --- norfuraneol loss -------------------------------------------------------------------------------------------
    "bork_nf_alone": dict(
        paper_keys=["bork2022", "bork2023"], dm="DM-01",
        where="Bork 2022 Fig. 1D (bork2022_extraction / bork2022b_extraction sec. 4.3); Bork 2023 Fig. 4B "
              "(bork2023_extraction sec. 4.3); CF-112, CF-113",
        lab="tu_berlin_kanzler", conditions="NF 20 mM alone in water, pH 5 initial, unbuffered, 130 C, 0-300 min; "
        "pH drift and O2 / headspace unstated", relation="lower",
        value=4.8e-3, temperature_c=130.0, conversion="first-order k, weighted fit 0-300 min (4.9e-3 refit, 5.1e-3 "
        "unweighted); Bork 2023 two-point 0-120 min 6.1e-3",
        # CF-113's x1.4 within-lab spread enters once (it covers Bork 2022 vs 2023, 4.8 vs 6.1e-3); S1 second reading
        terms=dict(noise=math.log10(4.98 / 4.8), within_lab_between_paper=math.log10(1.4)),
        record="no ledger row (proposed in the bork2022 / bork2023 dossiers, CF-112)",
        reader="second", reader_note="CF-112: second reader agrees on Fig. 1D; Bork 2023 Fig. 4B first reader",
        usable=True, why="the subset rule applied consistently: k_nf_decay lumps the carbonyl-sink family, and Bork "
        "2022's NF + MGO pot shows a co-reactant channel the engine lacks, so the NF-alone loss is a SUBSET of the "
        "lumped step: a lower bound (S1 second reading; the first reading took it as two-sided)"),
    "mikami_nf_ph8": dict(
        paper_keys=["mikami2017"], where="mikami2017_extraction; ledger trunk.yml (Mikami NF loss); "
        "source_verification/deoxyosone_submodel.yml", lab="ochanomizu",
        conditions="NF at pH 8, 100 C: 0.032 /min", relation="none", reader="second",
        usable=False, why="a different pH: carrying it to the engine's pH-free step needs a pH law (B61's hydroxide "
        "term is off): a rival structure (S3), not a band; CF-117 records an O2 / metal rival"),
    "lau_pandit_nf_matrix": dict(
        paper_keys=["lau2003", "pandit2006"], where="lau2003_extraction, pandit2006_extraction",
        lab="wsu_tang", conditions="NF in a gel / mashed potato, 116-121 C, pH not reported", relation="none",
        reader="first", usable=False,
        why="a matrix, pH unreported, a loss inferred from a formation plateau (<~6e-3 /min): recorded, not used"),
    # --- deoxypentosone exits (1-DX) --------------------------------------------------------------------------------
    "hou2026_1dx_exit": dict(
        paper_keys=["hou2026"], where="hou2026_extraction (lifetime bounds, 40-60 min window)", lab="jiangnan",
        conditions="glycine-xylose ARP 5 mM fed, initial pH 7.0, 120 C; pH MEASURED (~4.5 over the 40-60 min "
        "window); 1-DX by ex-situ OPD, semi-quantitative (CF-062)",
        relation="upper", value=0.29, temperature_c=120.0,
        conversion="upper end of the 1-DX lifetime bracket (all ARP loss routed to 1-DX): total 1-DX exit "
        "<= 0.29 /min, so each exit channel (r_dpo_nf) <= 0.29 /min", terms=dict(),
        record="dossier only", reader="first", usable=True, why="a superset (total exit) bounds the channel above"),
    "zhou2024b_1dx_exit": dict(
        paper_keys=["zhou2024b"], where="zhou2024b_extraction (derived lifetimes, 90-120 min window)", lab="jiangnan",
        conditions="ARP fed with Gly / Ala, 120 C, in-situ pH < 7 unknown; 1-DX by OPD, semi-quantitative (CF-062)",
        relation="upper", value=0.14, temperature_c=120.0,
        conversion="upper end of k1 bracket (all ARP loss -> 1-DX): total 1-DX exit <= 0.14 /min", terms=dict(),
        record="dossier only", reader="first", usable=True, why="a superset (total exit) bounds the channel above"),
    "nf_appearance_ea_wsu": dict(
        paper_keys=["lau2003", "pandit2006", "bornhorst2017b"], where="lau2003_extraction (81.5 +/- 3.2 kJ/mol), "
        "pandit2006_extraction (93.0 +/- 6.4), bornhorst2017b_extraction (104.9-122.3 at pH 8.4-9.5); ledger "
        "trunk.yml NF rows (status_uncertain)", lab="wsu_tang", conditions="ribose + amine -> NF, overall appearance",
        relation="none", reader="first", usable=False,
        why="the barrier of OVERALL NF appearance (supply + formation + loss): a rival structure, not the step's Ea"),
    # --- 3-DX -> furfural -------------------------------------------------------------------------------------------
    "furfural_formation_ea": dict(
        paper_keys=["schirlekeller1992", "garrett1969", "kang2026"], where="ledger sulfur.yml (Schirle-Keller 147, "
        "refit 156; B62); garrett1969_extraction (128.9 kJ/mol, specific-acid route, pH ~0); kang2026_SI_extraction "
        "(37.9 from a net 120-min level)", lab="minnesota / florida / jiangnan",
        conditions="net furfural formation from pentose", relation="none", reader="first", usable=False,
        why="net-formation barriers lump supply and loss (rival structures); they lie inside the class band"),
    # --- furfural loss ----------------------------------------------------------------------------------------------
    "almhofer_fur_water": dict(
        paper_keys=["almhofer2022"], where="almhofer2022_extraction Fig. 2B (printed Arrhenius fit; 145 C "
        "interpolated); ledger trunk.yml almhofer2022_fig2_furfural_loss_arrhenius_water", lab="jku_linz",
        conditions="furfural 52-520 mM in water only, 125-200 C, order 0.8-1.3", relation="lower",
        value=1.87e-4, temperature_c=145.0, conversion="the printed fit at 145 C (inside the measured 125-200 C)",
        terms=dict(), record="ledger", reader="first", usable=True,
        why="the water-only channel is a SUBSET of the lumped pot-independent furfural loss: a lower bound. Caveat: "
            "stainless-steel batch reactors; wall catalysis is a vessel term the bound carries unseparated"),
    "balagiannis_fur_k9": dict(
        paper_keys=["balagiannis2010"], where="balagiannis2010_extraction Table I (k9: furfural t1/2 25.4 min)",
        lab="reading", conditions="1:1 beef-sirloin / water extract, 130 C, pH not reported", relation="upper",
        value=2.73e-2, temperature_c=130.0,
        conversion="first-order k9 +/- 31 % as a 95 % HPD interval (one sigma = log10(1.31) / 1.96)",
        terms=dict(noise=math.log10(1.31) / 1.96), record="dossier (k9); ledger trunk.yml has the ribose row only",
        reader="first", usable=True,
        why="total furfural loss in a matrix (amines, thiols present) is a superset of the pot-independent channel. "
            "Caveat: k9 is a fitted constant of a nine-step lumped scheme, so the bound is inferred, not a direct "
            "measurement"),
    "yaghmur_fur_cys": dict(
        paper_keys=["yaghmur2005"], dm="DM-08", where="yaghmur2005_extraction (k_obs 3.2 / 6.7 / 12.6 / 16.0 e-2 /h "
        "at 40-70 C; Ea 46.5 kJ/mol); ledger sulfur.yml yaghmur rows", lab="hebrew_university",
        conditions="furfural + ~600 mM cysteine, pH 5, 40-70 C", relation="none", reader="first", usable=False,
        why="k_obs is proportional to [Cys] (a different rate law): a rival (B62's furfural barriers), not a band; "
            "note that B62 as built (k(145 C) held, Ea 46.5) gives ~0.021 /min at 70 C, ~0.9 dex above Yaghmur's own "
            "total-loss value there"),
    "furfural_loss_ea_other": dict(
        paper_keys=["almhofer2022", "garrett1969"], where="almhofer2022_extraction (49.8 +/- 3.6, water-only); "
        "garrett1969_extraction (89.1, acid solvolysis); ledger trunk.yml Danon 2014 review (48-135: a quote of "
        "quotes)", lab="jku_linz / florida", conditions="furfural loss channels", relation="none", reader="first",
        usable=False, why="barriers of SUBSET channels bound nothing on the lumped barrier; carried with DM-08 as the "
        "B62 rival (S3)"),
    # --- osone / pentose caramel ------------------------------------------------------------------------------------
    "hexose_osone_ea": dict(
        paper_keys=["martins2005", "knol2010", "deng2021", "tang2019"], where="martins2005_extraction (3-DG -> FA "
        "29.6, 1-DG -> AA 75.7); knol2010 (75 +/- 10); deng2021_extraction (pentose ARP ~75 +/- 10); tang2019_"
        "extraction (GSH-xylose ARP -> 3-DX 110.5, -> 1-DX 124.3, zero order)", lab="wageningen / jiangnan",
        conditions="deoxyosone / Amadori decay barriers", relation="none", reader="first", usable=False,
        why="no pentose-deoxyosone decay barrier exists; the hexose values are already members of the class"),
    "larralde_ribose": dict(
        paper_keys=["larralde1995"], where="larralde1995_extraction (Ea 107.4 kJ/mol; k0 2.86e-3 /min at 100 C); "
        "data/lit/fit_rows/heldout_larralde1995.yml", lab="ucsd_miller",
        conditions="amine-free ribose total loss, pD 7.4 phosphate, 60-118 C", relation="none", reader="second",
        usable=False, heldout=True, why="a SUPERSET barrier (all ribose losses) bounds nothing on one channel's barrier; a "
        "declared hold-out (TEST) row"),
    # --- pentose + cysteine steps -----------------------------------------------------------------------------------
    "zhang2025b_xylose_loss": dict(
        paper_keys=["zhang2025b"], where="zhang2025b_extraction (printed k_Xyl 0.073 /min; verified points ~0.061); "
        "ledger sulfur.yml; source_verification/route_submodel.yml (units corrected)", lab="btbu_xie",
        conditions="Cys 40 + Xyl 40 mM, 0.2 M phosphate pH 6.5, 120 C", relation="none", reader="second",
        usable=False, why="a superset upper bound on k x f_pH x [Cys]free (~20 mM at 5 min) whose free-cysteine "
        "trajectory needs an engine evaluation; Cao 2017's blank pot gives the closed-form ceiling instead"),
    "pentose_amine_analogues": dict(
        paper_keys=["balagiannis2010", "blank2006", "us8920862", "desclaux2006"], where="balagiannis2010 (ribose "
        "5.2e-2 /min, 130 C, beef extract); blank2006 (Xyl + Gly 1.1e-3 /min, 70 C); us8920862 (Xyl + Gly 0.50 /h, "
        "120 C, pH 6); desclaux2006 (ARP -> 1-DX 2.79e-2, -> 3-DX 5.5e-4 /min, 100 C)", lab="several",
        conditions="pentose with amines other than cysteine", relation="none", reader="first", usable=False,
        why="the amine is not cysteine and no measured law transfers between amines"),
    "amadori_step_ea": dict(
        paper_keys=["tang2019", "zhang2026", "martins2005"], where="tang2019 (110.5 / 124.3), zhang2026 k16 (85.7), "
        "martins2005 (DFG -> 3-DG 97.1, -> 1-DG 107.3)", lab="jiangnan / wageningen",
        conditions="Amadori-step barriers", relation="none", reader="first", usable=False,
        why="Amadori decay, not free pentose + Cys: a rival structure (the route sub-model)"),
    "cao2017_cys_ceiling": dict(
        paper_keys=["cao2017"], where="cao2017_extraction (every cysteine-destroying channel <= ~4e-4 /min at 90 C)",
        lab="jiangnan", conditions="pH 4.2-4.6, 90 C", relation="none", reader="second", usable=False,
        why="bounds k x [DPO], and DPO was not measured: no bound on k alone"),
    "cao2017_xylose_pool": dict(
        paper_keys=["cao2017"], where="cao2017_extraction Figs. 4-5 (blank pot: Xyl + TTCA + Cys-Amadori 0.478 M at "
        "1 h -> 0.438 M at 5 h; free Xyl 0.077 -> 0.069, free Cys 0.068 -> 0.056 M; pH 4.56 -> 4.15)",
        lab="jiangnan", conditions="xylose + cysteine 0.40 M each, blank (no glycine), 90 C, pH 5.49 -> 4.15",
        relation="upper",
        value=((0.478 - 0.438) + math.sqrt(6.0) * 0.004) * 1000.0 / 240.0 / (69.0 * 56.0), temperature_c=90.0,
        conversion="closed-form superset ceiling: every xylose-destroying channel (the pool decline plus its read "
        "error, 0.040 + 0.010 M over 240 min) over the lowest free [Xyl][Cys] (69 x 56 mM^2) bounds k x f_pH; the "
        "engine's pH factor is taken at its smallest over the window's pH 4.15-4.56 (the loose edge)",
        terms=dict(), record="dossier (digitised Figs. 4-5)", reader="second",
        reader_note="Cao 2017 Figs. 4-5 read by two readers (dossier); the ceiling is the S1 second reader's, "
        "re-derived here", usable=True,
        why="a closed xylose pool: its decline bounds every pentose-consuming channel above, so each enolisation "
            "step (r_pent_tdp, r_pent_dpo) is bounded"),
    "hexosulose_cys_analogue": dict(
        paper_keys=["wedzicha1991", "edwards1992"], where="wedzicha1991_extraction, edwards1992_extraction; ledger "
        "sulfur.yml", lab="leeds", conditions="3-deoxyhexosulose + Cys, thiol liberation ~0.07 /h, 55 C, pH 5.5",
        relation="none", reader="first", usable=False,
        why="a hexose deoxyosone at 55 C, no measured transfer law to DPO + Cys at 145 C"),
    # --- TTCA -------------------------------------------------------------------------------------------------------
    "zhai2021_ttca": dict(
        paper_keys=["zhai2021"], dm="DM-07", where="zhai2021_extraction (zero-order fits; first-order equivalents "
        "0.0027 / 0.0065 / 0.0081 /min at 100 / 120 / 140 C); ledger sulfur.yml zhai2021_fig3a_ttca_tladder; "
        "source_verification/route_submodel.yml (second reader)", lab="jiangnan",
        conditions="fed TTCA ~10 mM, unbuffered water, natural pH 5.5, 100-140 C (x9.6 faster at pH 8, CF-110)",
        relation="upper", value=math.log(10.0 / 7.59) / 20.0, temperature_c=140.0,
        conversion="the EARLY first-order rate, Fig. 3a 140 C: 10 -> 7.59 mM in 20 min (0.0138 /min; local 20-40 min "
        "0.0128), not the zero-order line's 0.0081 (S1 second reading)",
        terms=dict(), record="ledger", reader="second", usable=True,
        why="every TTCA loss (a superset of r_ttca_deg) bounds it above"),
    "cao2017_ttca_amadori": dict(
        paper_keys=["cao2017"], where="cao2017_extraction (TTCA -> Cys-Amadori >= ~1.0e-3 /min at 90 C)",
        lab="jiangnan", conditions="90 C, pH 4.2-4.6", relation="none", reader="second", usable=False,
        why="counts toward k_ttca_deg only if the Cys-Amadori route is non-releasing; the engine's r_cama_* release "
        "Cys: model-conditional"),
    "zhai2021_ttca_ea": dict(
        paper_keys=["zhai2021"], where="zhai2021_extraction; CF-108 (printed 80.99 not reproducible; slope ratio "
        "~35; first-order whole-window ~70)", lab="jiangnan", conditions="fed TTCA 100-140 C", relation="none",
        reader="second", usable=False, why="a SUPERSET barrier that depends on the assumed rate law"),
    # --- glucose C2/C3 ----------------------------------------------------------------------------------------------
    "li2011_glucose_loss": dict(
        paper_keys=["li2011"], where="li2011_extraction sec. 4 (K_glc 0.0013 / 0.0032 / 0.0073 / 0.0151 /min at "
        "100-130 C); ledger sulfur.yml li2011 rows", lab="jiangnan_chen",
        conditions="glucose + cysteine, 100-130 C", relation="upper", value=0.0151, temperature_c=130.0,
        conversion="total first-order glucose loss", terms=dict(), record="ledger", reader="first", usable=True,
        why="total glucose loss is a superset of the C2/C3 cleavage"),
    "thornalley1999_mgo": dict(
        paper_keys=["thornalley1999"], where="thornalley1999_extraction (glucose 50 mM, MGO initial rate 0.098 "
        "+/- 0.003 uM/day); data/lit/fit_rows/heldout_thornalley1999_figs1b_2e_2f_3c.yml", lab="essex",
        conditions="glucose alone, 100 mM phosphate pH 7.4, 37 C", relation="none", reader="second", usable=False,
        heldout=True,
        why="a lower bound (~1.4e-9 /min) only if the engine's other glucose -> MGO routes are minor (model-"
        "conditional), and a declared hold-out TEST row: using it as a band source spends it (a plan-level call)"),
    "glucose_ea_rivals": dict(
        paper_keys=["li2011", "martins2005"], where="li2011 (refit 102.4; printed 87.85 does not reproduce); "
        "martins2005 (Glc -> FA + AA 236.7 +/- 63.4; DFG -> MG 124.5)", lab="jiangnan / wageningen",
        conditions="glucose loss barriers", relation="none", reader="first", usable=False,
        why="superset or different-step barriers; inside the class band"),
    # --- thiols -----------------------------------------------------------------------------------------------------
    "grosch1992_mft_inpot": dict(
        paper_keys=["grosch1992"], where="grosch1992_extraction Table VIII (inferred in-pot lambda 0.044, "
        "0.035-0.053 /min at 100 C); ledger sulfur.yml grosch1992 row", lab="dfa_garching",
        conditions="cysteine / ribose pot, 100 C; loss inferred from formation and level", relation="upper",
        record="ledger", reader="first", usable=False,
        why="removed as an upper source (S1 second reading): the lambdas are inferred under constant formation and "
            "read closer to LOWER bounds on the in-pot loss; the 0.044 /min is the thiamine-only pot, not a cysteine / "
            "ribose pot; the Method I rows are a declared hold-out (heldout_grosch1992_tab2.yml)"),
    "patent_mft_standard": dict(
        paper_keys=["us10076126"], dm="DM-02", where="us10076126 Methods (co-heated [2H3]-MFT survives 240 min); "
        "results/validation/m0_checks_batch2.md CK-23; ledger us10076126_fig1_mft_time_course", lab="patent_assignee",
        conditions="xylose + cysteine 0.25 M, 0.2 M phosphate pH 5.5, 110 C", relation="upper", value=0.058,
        temperature_c=110.0, conversion="lambda <= ln(1/f_min)/240 min at the loose corner f_min 1e-6 (declared, "
        "not printed)", terms=dict(), record="ledger", reader="second", usable=True,
        why="a patent counts as a bound only (S2 (d)); with Grosch 1992 it nominates a cysteine-slowed loss (case "
            "M1 C2): a rival route source as well"),
    "vanseeventer_mft": dict(
        paper_keys=["vanseeventer2001"], where="vanseeventer2001_extraction Table 1 (MFT 59 %/day, FFT 28 %/day, "
        "50 C)", lab="nestle", conditions="process flavour, pH 5.0, 0.5 M phosphate", relation="none",
        reader="first", usable=False, why="MFT loss by C-5 oligomerisation (ch_oligomer, carried at zero rate): a "
        "rival structure; a declared hold-out"),
    "thiol_net_loss": dict(
        paper_keys=["zhai2023jafc", "zhang2024"], where="k6a_sulfur_ladders_synthesis (Zhai 140 C k_net MFT 0.0065, "
        "FFT 0.0141 /min); Zhang2024_extraction (21-55 % in 60 min, 115 C)", lab="jiangnan",
        conditions="net thiol turnover in pots", relation="none", reader="first", usable=False,
        why="a net loss is a LOWER bound on the total loss, which bounds no single channel from below"),
    "gigl2021_capture": dict(
        paper_keys=["gigl2021"], where="k6a_sulfur_ladders_synthesis (k 1.17e-5 / 2.57e-4 / 7.83e-4 /min at "
        "279-333 K; Ea 60.2, range 7-102); ledger sulfur.yml; source_verification/thiol_submodel.yml (conditions "
        "FAIL; Fig. 4 not read by the second reader)", lab="tu_munich",
        conditions="FFT covalent capture by coffee melanoidins, 6-60 C", relation="none", reader="first",
        usable=False, why="a same-class analogue with no measured transfer law (melanoidin capture at 6-60 C to a "
        "pot-independent sink at 100-145 C, 3-5 decades): not a legitimate analogue band; its central value is "
        "already a member of the sulfur class"),
    "kumazawa_fft_ph3": dict(
        paper_keys=["kumazawa2003"], where="kumazawa2003_extraction Fig. 3 (pH 3.0: 99.5 % residual, read error "
        "+/- 1.5 pp); ledger sulfur.yml (in_model_fit); source_verification/thiol_submodel.yml (adjudicator "
        "re-digitised)", lab="ogawa", conditions="FFT 1 ppm in McIlvaine buffer pH 3.0, 121 C, 10 min, canned in air",
        relation="none", reader="second", usable=False,
        why="bounds k_fft_decay.log10_k_ref (not one of the 43) jointly with the family barrier: see the joint note"),
    # --- H2S --------------------------------------------------------------------------------------------------------
    "li2022_h2s": dict(
        paper_keys=["li2022"], where="li2022_extraction (sink 0.023-0.031 /min if all release deficit is loss); "
        "ledger sulfur.yml", lab="copenhagen_arla", conditions="Cys 5.7 mM, pH 7, 10 mM phosphate, sealed vial, "
        "air headspace, 90 C, n = 1", relation="none", reader="first", usable=False,
        why="inferred through the Zheng & Ho release law (model-conditional), not a measured bound"),
    "cheng2020_h2s": dict(
        paper_keys=["cheng2020"], where="cheng2020_extraction (H2S falls 3-36 %, 100 C, 60-90 s); ledger sulfur.yml",
        lab="cheng_lab", conditions="mandarin juice, sealed vial, semi-quantitative (EMS)", relation="none",
        reader="first", usable=False, why="a net loss in a matrix (metals, ascorbate): a lower bound on a total, "
        "not the channel"),
    # --- furfural + sulfide -> FFT ----------------------------------------------------------------------------------
    "zhai2023_fur_sulfide": dict(
        paper_keys=["zhai2023jafc"], where="zhai2023jafc_extraction sec. 11, Fig. 6a; ledger sulfur.yml; "
        "source_verification/thiol_submodel.yml (verified)", lab="jiangnan",
        conditions="furfural + (NH4)2S 10 mM each, label pH 7.0 (in situ 6.33), 120 C; dissolved sulfide not "
        "measured (34-70 %, CF-065)", relation="upper", value=math.log(1.00 / 0.22) / 10.0 / (10.0 * 0.34),
        temperature_c=120.0, conversion="total furfural loss ln(1.00/0.22)/10 min = 0.151 /min over dissolved "
        "sulfide 10 mM x 0.34 (the loose end of CF-065's range); a ceiling on the SUM of the furfural + sulfide "
        "channels, taken per channel at speciation weight 1 (the tightest case)", terms=dict(),
        record="ledger", reader="second", usable=True,
        why="a total loss bounds each furfural + sulfide channel above"),
    "munch1997_order": dict(
        paper_keys=["munch1997"], where="ledger sulfur.yml (Munch 1997 T5: order in H2S 1.04 [0.90, 1.18]); "
        "CF-065", lab="dfa_garching", conditions="furfural + H2S", relation="none", reader="first", usable=False,
        why="a measured LAW (first order in H2S) that supports the rate law; not a value"),
    "formation_ea_rivals": dict(
        paper_keys=["meng", "feng2022", "zhai2023jafc", "ames"], where="k6a_sulfur_ladders_synthesis (net-yield "
        "barriers: Meng 57.8 / 72.9, Feng 30.3 / 44.4, Zhai 14.5-74.8, Ames 103.5 / 107.7)", lab="several",
        conditions="apparent thiol net-yield barriers", relation="none", reader="first", usable=False,
        why="net-yield barriers lump supply, formation and loss: rival structures"),
    "zhai2023_nf_mft": dict(
        paper_keys=["zhai2023jafc"], where="zhai2023jafc_extraction Fig. 6b (NF 10.6 mM gone by 10 min; MFT 3.2 mM "
        "at 10 min)", lab="jiangnan", conditions="NF + (NH4)2S, 120 C", relation="none", reader="second",
        usable=False, why="bounds k_nf_mft.log10_k_ref (not one of the 43) from below; product rising after the "
        "substrate is gone points to an intermediate (a rival to the one-step route)"),
    "zhai2023_mgo_sulfide": dict(
        paper_keys=["zhai2023jafc"], where="zhai2023jafc_extraction Fig. 6e (MGO + GO 0.67 -> 0.27 g/L in 30 min)",
        lab="jiangnan", conditions="dicarbonyls + (NH4)2S, 120 C; MP not measured", relation="none",
        reader="second", usable=False, why="bounds k_mgo_mp.log10_k_ref (not one of the 43); non-binding"),
    "dimer_subset": dict(
        paper_keys=["zhang2024", "zhou2023", "kumazawa2003", "hofmann2001"], where="Zhang2024_extraction Fig. 1 "
        "(dimer shares, 115 C; verification conditions FAIL); ledger sulfur.yml zhou2023 (6.5-9.6 %, pH-invariant); "
        "kumazawa2003 Table 3; hofmann2001 (< 1.5 % at 30 C)", lab="several",
        conditions="disulfide shares", relation="none", reader="first", usable=False,
        why="the oxidant pool OX is never measured, so k (L^2/mmol^2/min) cannot be separated from [OX]; a share "
        "bounds the cumulative flux, not k; the encounter ceiling (log10 9.47) is physical and non-binding"),
    # --- acrylamide scavenging --------------------------------------------------------------------------------------
    "claeys2005_gln": dict(
        paper_keys=["claeys2005"], where="claeys2005_extraction Table 2 (glutamine row: k_Eref 274.1 +/- 81.6 e-3 "
        "/min at T_ref 160 C); ledger acrylamide.yml claeys2005 rows", lab="ku_leuven",
        conditions="Asn + Glc 0.01 M, + Gln 0.01 M, 0.05 M citrate pH 6, 140-200 C", relation="upper",
        value=0.2741 / 10.0, temperature_c=160.0, conversion="lumped first-order elimination k_E over [Gln] 10 mM "
        "(the bound holds while [Gln] stays near 10 mM; glucose consumes it)",
        terms=dict(noise=math.log10((274.1 + 81.6) / 274.1)), record="ledger", reader="first", usable=True,
        why="the lumped elimination k_E includes every acrylamide loss: a ceiling on the Gln channel"),
    "claeys2005_lys": dict(
        paper_keys=["claeys2005"], where="claeys2005_extraction Table 2 (lysine row: k_Eref 280.2 +/- 43.4 e-3 /min "
        "at T_ref 160 C)", lab="ku_leuven", conditions="Asn + Glc 0.01 M, + Lys 0.01 M, 0.05 M citrate pH 6",
        relation="upper", value=0.2802 / 10.0, temperature_c=160.0, conversion="k_E over [Lys] 10 mM",
        terms=dict(noise=math.log10((280.2 + 43.4) / 280.2)), record="ledger", reader="first", usable=True,
        why="a ceiling on the Lys channel"),
    "devleeschouwer2009_gln": dict(
        paper_keys=["devleeschouwer2009"], where="devleeschouwer2009_extraction Table 3 (Gln column; k_E2 not in "
        "the model)", lab="ku_leuven", conditions="freeze-dried powder, a_w 0.92, 160 C", relation="none",
        reader="first", usable=False, why="glutamine deamidates (half-life ~1 min), so no floor on [Gln]; a declared "
        "hold-out; a different a_w"),
    "acr_cys_analogue": dict(
        paper_keys=["devleeschouwer2009", "claeys2005", "zamora2010"], where="devleeschouwer2009 Table 3 Cys column "
        "(49.36 +/- 1.18 /M/min, Ea 51.3 +/- 1.5); zamora2010_extraction (Ea Gly 52, butylamine 44; rates figure-"
        "only)", lab="ku_leuven / seville", conditions="acrylamide + Cys / amines", relation="none", reader="first",
        usable=False, why="a same-class analogue with no measured thiol -> amine law (Friedman's 100-300x is quoted "
        "inside Claeys, not measured in the corpus); Zamora's rates are figure-only"),
    # --- Strecker ---------------------------------------------------------------------------------------------------
    "zhou2024_ala": dict(
        paper_keys=["zhou2024"], where="zhou2024_extraction Table 2 (pyrazine 0.0279, 2,5-DMP 0.0035 umol/L/min at "
        "100 C); ledger trunk.yml (in_model_fit)", lab="jiangnan_zhang_ho",
        conditions="GO / MGO + alanine 20 + 20 mM, water, initial pH 8, 100-120 C; no SE printed", relation="none",
        reader="first", usable=False, why="a lower bound for ALANINE at pH 8 (>= -6.85 for GO); the alanine -> "
        "glycine transfer is declared (+/- 0.5 dex), not measured: no analogue band"),
    "hofmann2000_phe": dict(
        paper_keys=["hofmann2000"], where="hofmann2000 Tables 3-4; source_verification/strecker_hofmann2000.yml "
        "(second reader + adjudicator)", lab="dfa_garching", conditions="GO / MGO + Phe 100 + 100 mM, phosphate "
        "pH 7, 98 C", relation="none", reader="second", usable=False,
        why="Strecker products are a subset: lower bounds for PHENYLALANINE (>= -5.72 GO, argon); no Phe -> Gly law"),
    "sedehi2013_gly": dict(
        paper_keys=["sedehi2013"], where="sedehi2013_extraction (transport to 100 C with the class Ea 50 +/- 11); "
        "data/lit/fit_rows/heldout_sedehi2013_tab1.yml (second reader)", lab="usd_dehaan",
        conditions="GO + Gly total consumption, 4-25 C, 0.2-1.5 M, D2O, pH 2-7", relation="upper",
        heldout=True, record="hold-out fit row", reader="second", usable=False,
        why="a declared hold-out row (heldout_sedehi2013_tab1.yml): held-out rows are never band sources (coordinator "
            "ruling, as Thornalley 1999 and Larralde 1995). As a ceiling it would be ~-3.7 at 100 C (central -4.4), "
            "above the search box anyway"),
    "strecker_cascade_ea": dict(
        paper_keys=["zhou2024", "cremer2000", "yu2018", "leahy1989"], where="zhou2024 printed Ea (no SE; printed vs "
        "refit gives the B18 bands); cremer2000 (115-124); yu2018 (99.8 +/- 6.7, 104.1 +/- 5.1, pH 10); "
        "leahy1989 (150-177)", lab="several", conditions="pyrazine-cascade barriers", relation="none",
        reader="first", usable=False, why="the alanine barrier has no measured transfer to glycine; the rest are "
        "cascade barriers (rival structures); all lie inside the class band"),
    # --- the two DM sources outside the 43 --------------------------------------------------------------------------
    "davidek2008_aa": dict(
        paper_keys=["davidek2008"], dm="DM-03", where="davidek2008 (acetic acid plateau 26.9 -> 27.6 mol %, "
        "60-240 min); CF-069", lab="nestle_davidek", conditions="xylose + glycine, 0.2 M phosphate pH 8, 120 C",
        relation="upper", value=4e-4, temperature_c=120.0, conversion="loss <= ~4e-4 /min if formation has stopped",
        terms=dict(), record="no ledger row", reader="first", usable=True,
        why="an inferred ceiling (formation over the plateau unstated); pentose and pH 8 vs the fit's glucose, pH 6.8"),
    "kanzler_mgo_ea": dict(
        paper_keys=["kanzler2026", "kanzler2021"], dm="DM-04", where="kanzler2026_extraction secs. 4.7-6; CF-122",
        lab="tu_berlin_kanzler", conditions="MGO + Ala / Phe / Pro 100 mM, pH 5 unbuffered, 60-140 C",
        relation="none", reader="first", usable=False,
        why="second order in MGO, ~[amine]^0.75: a rival MGO-sink structure (CF-122); inside the class band"),
}

# ---------------------------------------------------------------------------------------------------------------------
# Per coordinate: its sources and, for an Ea, the B55 class whose measured distribution is the band (None = no class);
# for a log10 k_ref, the class whose Ea band carries a bound across temperature.
# ---------------------------------------------------------------------------------------------------------------------
NONE_ON_DISK = "nothing on disk beyond the search bound"

COORDINATES: Dict[str, Dict[str, Any]] = {
    # log10 k_ref
    "k_dpo_nf.log10_k_ref": dict(sources=["hou2026_1dx_exit", "zhou2024b_1dx_exit"], transport=CARB),
    "k_nf_decay.log10_k_ref": dict(sources=["bork_nf_alone", "mikami_nf_ph8", "lau_pandit_nf_matrix"],
                                   transport=CARB),
    "k_fur_fft_hs.log10_k_ref": dict(sources=["zhai2023_fur_sulfide", "munch1997_order"], transport=SULF,
                                     factors={"zhai2023_fur_sulfide": ("hs_anion", (6.33,), 393.15)},
                                     note="k_ref is normalised at pH 5: the engine's HS- factor at the in-situ pH "
                                          "6.33, 120 C carries the bound to the coordinate (S1 second reading; with "
                                          "Suleimenov's pK1 6.45 the second reader gets x13.8 and raw -1.44; at the "
                                          "label pH 7.0, raw -1.67, 95 % 0.22)"),
    "k_mft_decay.log10_k_ref": dict(sources=["grosch1992_mft_inpot", "patent_mft_standard", "vanseeventer_mft",
                                             "thiol_net_loss"], transport=SULF),
    "k_h2s_loss.log10_k_ref": dict(sources=["li2022_h2s", "cheng2020_h2s"], transport=SULF),
    "k_pent_tdp.log10_k_ref": dict(sources=["cao2017_xylose_pool", "zhang2025b_xylose_loss", "pentose_amine_analogues"],
                                   transport=SULF, factors={"cao2017_xylose_pool": ("acid", (4.15, 4.56), 363.15)}),
    "k_pent_dpo.log10_k_ref": dict(sources=["cao2017_xylose_pool", "zhang2025b_xylose_loss", "pentose_amine_analogues"],
                                   transport=SULF, factors={"cao2017_xylose_pool": ("base", (4.15, 4.56), 363.15)}),
    "k_dpo_ddp.log10_k_ref": dict(sources=["cao2017_cys_ceiling", "hexosulose_cys_analogue"], transport=SULF),
    "k_acr_gln.log10_k_ref": dict(sources=["claeys2005_gln", "devleeschouwer2009_gln", "acr_cys_analogue"],
                                  transport=None),
    "k_acr_lys.log10_k_ref": dict(sources=["claeys2005_lys", "acr_cys_analogue"], transport=None),
    "k_ddp_mft_hs.log10_k_ref": dict(sources=[], transport=SULF, note=NONE_ON_DISK + " (no fed 1,4-dideoxyosone "
                                     "anywhere; Cerny 2003 route shares are ratios)"),
    "k_ddp_mft.log10_k_ref": dict(sources=[], transport=SULF, note=NONE_ON_DISK + " (as k_ddp_mft_hs)"),
    "k_fur_decay.log10_k_ref": dict(sources=["almhofer_fur_water", "balagiannis_fur_k9", "yaghmur_fur_cys"],
                                    transport=CARB),
    "k_ttca_deg.log10_k_ref": dict(sources=["zhai2021_ttca", "cao2017_ttca_amadori"], transport=SULF),
    "k_glc_ha.log10_k_ref": dict(sources=["li2011_glucose_loss", "thornalley1999_mgo"], transport=CARB),
    "k_fur_fft.log10_k_ref": dict(sources=["zhai2023_fur_sulfide", "munch1997_order"], transport=SULF,
                                  factors={"zhai2023_fur_sulfide": ("neutral_h2s", (6.33,), 393.15)},
                                  note="the engine's neutral-H2S factor at pH 6.33, 120 C (the second reader's 0.55)"),
    "k_dimer_mft.log10_k_ref": dict(sources=["dimer_subset"], transport=SULF),
    "k_go_ak.log10_k_ref": dict(sources=["sedehi2013_gly", "zhou2024_ala", "hofmann2000_phe"], transport=CARB),
    "k_dimer_fft.log10_k_ref": dict(sources=["dimer_subset"], transport=SULF),
    # Ea: the thiol-sink family (one shared barrier)
    "k_fft_decay.ea_kj_mol": dict(sources=["gigl2021_capture", "kumazawa_fft_ph3"], ea_class=SULF,
                                  family="thiol_sink"),
    "k_mft_decay.ea_kj_mol": dict(sources=["gigl2021_capture"], ea_class=SULF, family="thiol_sink"),
    "k_thiolate_loss.ea_kj_mol": dict(sources=["gigl2021_capture"], ea_class=SULF, family="thiol_sink"),
    # Ea: Strecker (trunk)
    "k_mgo_ak.ea_kj_mol": dict(sources=["strecker_cascade_ea"], ea_class=CARB),
    "k_go_ak.ea_kj_mol": dict(sources=["strecker_cascade_ea"], ea_class=CARB),
    # Ea: the carbonyl-sink family (one shared barrier)
    "k_nf_decay.ea_kj_mol": dict(sources=[], ea_class=CARB, family="carbonyl_sink",
                                 note="no measured NF-decay barrier (Mikami 2017, Bork print none)"),
    "k_fur_decay.ea_kj_mol": dict(sources=["yaghmur_fur_cys", "furfural_loss_ea_other"], ea_class=CARB,
                                  family="carbonyl_sink"),
    "k_osone_decay.ea_kj_mol": dict(sources=["hexose_osone_ea"], ea_class=CARB, family="carbonyl_sink"),
    # Ea: the B55 single-temperature constants (already sampled from their class)
    "k_dpo_nf.ea_kj_mol": dict(sources=["nf_appearance_ea_wsu"], ea_class=CARB),
    "k_tdp_fur.ea_kj_mol": dict(sources=["furfural_formation_ea"], ea_class=CARB),
    "k_ttca_deg.ea_kj_mol": dict(sources=["zhai2021_ttca_ea"], ea_class=SULF),
    "k_glc_ha.ea_kj_mol": dict(sources=["glucose_ea_rivals"], ea_class=CARB),
    "k_pent_tdp.ea_kj_mol": dict(sources=["amadori_step_ea"], ea_class=SULF),
    "k_pent_dpo.ea_kj_mol": dict(sources=["amadori_step_ea"], ea_class=SULF),
    "k_dpo_ddp.ea_kj_mol": dict(sources=[], ea_class=SULF, note="no measured barrier"),
    "k_pent_caramel.ea_kj_mol": dict(sources=["larralde_ribose"], ea_class=CARB),
    "k_nf_mft.ea_kj_mol": dict(sources=["formation_ea_rivals", "zhai2023_nf_mft"], ea_class=SULF),
    "k_fur_fft_hs.ea_kj_mol": dict(sources=["formation_ea_rivals"], ea_class=SULF),
    "k_h2s_loss.ea_kj_mol": dict(sources=[], ea_class=SULF, note="no measured H2S-loss barrier; no class member "
                                 "is an H2S-loss step"),
    "k_fur_fft.ea_kj_mol": dict(sources=["formation_ea_rivals"], ea_class=SULF),
    "k_mgo_mp.ea_kj_mol": dict(sources=["zhai2023_mgo_sulfide"], ea_class=SULF),
    "k_nf_mp3p.ea_kj_mol": dict(sources=[], ea_class=SULF, note="no measured barrier"),
    "k_ddp_mft_hs.ea_kj_mol": dict(sources=[], ea_class=SULF, note="no measured barrier"),
    "k_ddp_mft.ea_kj_mol": dict(sources=[], ea_class=SULF, note="no measured barrier"),
}

#: The two DM sources outside the 43 (the audit's contradicting measurements).
EXTRA: Dict[str, Dict[str, Any]] = {
    "k_aa_frag.log10_k_ref": dict(sources=["davidek2008_aa"], transport=CARB),
    "k_mgo_mel.ea_kj_mol": dict(sources=["kanzler_mgo_ea"], ea_class=CARB),
}

CLASS_REASON = {
    CARB: "a sugar-derived carbonyl or dicarbonyl reactant and no sulfur-bearing reactant (B55's carbohydrate class)",
    SULF: "a sulfur-bearing reactant (B55's sulfur class)",
}

#: The S1 class ruling (both readers, 2026-11-06), B55's own definition: a step with a sulfur-bearing reactant is class
#: sulfur. B55's CLASS_OF_ROUTE keys the class by the ROUTE a constant serves (FORMATION_ROUTE_OF), which conflates two
#: principles (the route served vs the step's reactant chemistry): it put k_ttca_deg, k_pent_tdp, k_pent_dpo, k_dpo_ddp
#: and the class member Zhang 2026 k16 (Cys-Amadori -> alpha-dicarbonyl) in the carbohydrate class. Here the member
#: moves to the sulfur class (ranges unchanged; n 19 / 16); the four coordinates are ruled sulfur in COORDINATES.
#: ea_class_prior.py itself is NOT changed (that would move shipped outputs; an S4 case applies it).
CLASS_RULING_MEMBER_MOVES: Mapping[str, str] = {"parameters_sulfur.ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL": SULF}
CLASS_RULING_COORDINATES = ("k_ttca_deg.ea_kj_mol", "k_pent_tdp.ea_kj_mol", "k_pent_dpo.ea_kj_mol",
                            "k_dpo_ddp.ea_kj_mol")

#: A between-lab term for the acrylamide class (coordinator ruling under S2 (b): "the measured spread of that step or its
#: reaction class where >= 2 laboratories measured it"): first-order acrylamide elimination at 160 C, Knol 2005
#: (Wageningen; 88.1e-3 /min, knol2005_extraction) vs Claeys 2005 (Leuven; control 111.1e-3 /min, Table 2).
ACRYLAMIDE_ELIMINATION_BETWEEN_LAB_DEX = math.log10(111.1 / 88.1)
EXTRA_LANE_BETWEEN_LAB = {"acrylamide": ACRYLAMIDE_ELIMINATION_BETWEEN_LAB_DEX}
EXTRA_LANE_BETWEEN_LAB_SOURCE = {
    "acrylamide": "acrylamide elimination at 160 C, Knol 2005 (Wageningen) vs Claeys 2005 (Leuven), 0.10 dex apart "
                  "(coordinator ruling under S2 (b): the reaction class's measured spread with >= 2 laboratories)"}

#: The S1 second reading (2026-11-06) re-read every row of the first version of this record (#100, c70038d0). Rows not
#: listed were agreed as read; a listed row was corrected and the correction is applied here.
SECOND_READING_DATE = "2026-11-06"
SECOND_READING_CORRECTIONS: Mapping[str, str] = {
    "k_nf_decay.log10_k_ref": "A1: one-sided LOWER bound by the subset rule; CF-113's x1.4 counted once",
    "k_fur_fft_hs.log10_k_ref": "A2: the HS- speciation factor (k_ref normalised at pH 5) carries the bound",
    "k_mft_decay.log10_k_ref": "A3: Grosch 1992 removed as an upper source; the edge stands on the patent alone",
    "k_fur_decay.log10_k_ref": "A4: Balagiannis's +/- 31 % is a 95 % HPD; Almhofer's steel-wall term and k9's "
                               "inferred status noted",
    "k_ttca_deg.log10_k_ref": "A5: the early first-order rate (0.0138 /min), carried over the sulfur class",
    "k_fur_fft.log10_k_ref": "A6: the neutral-H2S factor at pH 6.33 carries the bound",
    "k_dpo_nf.log10_k_ref": "A7: CF-062 semi-quantitative 1-DX flag; Hou's pH was measured (~4.5)",
    "k_pent_tdp.log10_k_ref": "D: Cao 2017's blank-pot xylose-pool ceiling (not refused)",
    "k_pent_dpo.log10_k_ref": "D: Cao 2017's blank-pot xylose-pool ceiling (not refused)",
    "k_acr_gln.log10_k_ref": "D: the acrylamide-class between-lab term (Knol 2005 vs Claeys 2005) defines the edge",
    "k_acr_lys.log10_k_ref": "D: the acrylamide-class between-lab term (Knol 2005 vs Claeys 2005) defines the edge",
    "k_go_ak.log10_k_ref": "D: Sedehi 2013 dropped (a held-out row is never a band source); stays refused",
    "k_ttca_deg.ea_kj_mol": "C: class sulfur (B55's own definition)",
    "k_pent_tdp.ea_kj_mol": "C: class sulfur (B55's own definition)",
    "k_pent_dpo.ea_kj_mol": "C: class sulfur (B55's own definition)",
    "k_dpo_ddp.ea_kj_mol": "C: class sulfur (B55's own definition)",
    "k_dpo_ddp.log10_k_ref": "C: transported over the sulfur class",
}


# ---------------------------------------------------------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------------------------------------------------------
def arrhenius_shift_dex(ea_kj_mol: float, from_c: float, to_c: float) -> float:
    """log10 k(to) - log10 k(from) for barrier ``ea_kj_mol``."""
    return ea_kj_mol / (R_KJ * LN10) * (1.0 / (from_c + 273.15) - 1.0 / (to_c + 273.15))


def quadrature(terms: Iterable[Optional[float]]) -> float:
    return math.sqrt(sum(float(t) ** 2 for t in terms if t))


def engine_factor(spec: Mapping[str, Any], source_id: str) -> Tuple[float, Optional[str]]:
    """The engine's own pH factor for this coordinate in that source's pot (k_ref is normalised at pH 5), at its
    smallest over the pot's pH range (the loose edge for a ceiling on k); 1 when the step has none there."""
    f = (spec.get("factors") or {}).get(source_id)
    if not f:
        return 1.0, None
    from src.kinetic_core.parameters_sulfur import ph_factor

    kind, phs, t_k = f
    vals = [ph_factor(kind, float(ph), float(t_k)) for ph in phs]
    return min(vals), (f"the engine's {kind} factor at pH {'-'.join(f'{p:g}' for p in phs)}, {t_k - 273.15:.0f} C: "
                       + ", ".join(f"{v:.3g}" for v in vals) + " (the smallest taken)")


def joint_check(src: Mapping[str, Any], *, factor: float, shipped_k: float, shipped_ea: Optional[float],
                t_ref_c: float, between_lab: Optional[float]) -> Optional[Dict[str, Any]]:
    """The SHIPPED (k, Ea) pair against one bound at the bound's own temperature (no transport over a band)."""
    if shipped_ea is None and abs(float(src["temperature_c"]) - t_ref_c) > 1e-9:
        return None
    t = float(src["temperature_c"])
    at_t = shipped_k + (0.0 if shipped_ea is None else arrhenius_shift_dex(shipped_ea, t_ref_c, t))
    bound = math.log10(float(src["value"]) / factor)
    sig_raw = quadrature(src.get("terms", {}).values())
    sig_95 = None if between_lab is None else quadrature([sig_raw, between_lab])
    out = {"shipped_log10_k_at_T": at_t, "bound_log10_k_at_T": bound, "temperature_c": t, "shipped_ea": shipped_ea}
    sides = ("upper",) if src["relation"] == "upper" else ("lower",) if src["relation"] == "lower" else ("upper", "lower")
    for side in sides:
        excess = (at_t - bound) if side == "upper" else (bound - at_t)  # > 0: the shipped pair is beyond the bound
        out[side] = {"excess_dex": excess, "margin_raw_dex": Z95 * sig_raw,
                     "margin_95_dex": None if sig_95 is None else Z95 * sig_95,
                     "violates_raw": excess > Z95 * sig_raw,
                     "violates_95": None if sig_95 is None else excess > Z95 * sig_95}
    return out


def evaluate_bound(src: Mapping[str, Any], *, t_ref_c: float, ea_band: Optional[Tuple[float, float]],
                   between_lab: Optional[float], search: Tuple[float, float], shipped: float,
                   factor: float = 1.0) -> Dict[str, Any]:
    """One measured bound carried to the coordinate (log10 at t_ref_c): raw and S2 (b) 95 % edges, binding, violation.
    ``factor``: the engine's pH factor in the source's pot (the measured rate is k x factor)."""
    log_meas = math.log10(float(src["value"]) / factor)
    t = float(src["temperature_c"])
    if abs(t - t_ref_c) < 1e-9:
        shifts = (0.0, 0.0)
    elif ea_band is None:
        return {"transport": "undefined: the bound is at another temperature and the coordinate's barrier has no band"}
    else:
        s = [arrhenius_shift_dex(e, t, t_ref_c) for e in ea_band]
        shifts = (min(s), max(s))
    sig_raw = quadrature(src.get("terms", {}).values())
    sig_95 = None if between_lab is None else quadrature([sig_raw, between_lab])
    out: Dict[str, Any] = {"log10_at_measured_T": log_meas, "transport_dex": list(shifts), "sigma_raw_dex": sig_raw,
                           "sigma_95_dex": sig_95}
    lo_s, hi_s = search
    for side in (("upper",) if src["relation"] == "upper" else ("lower",) if src["relation"] == "lower"
                 else ("upper", "lower")):
        sign = 1.0 if side == "upper" else -1.0
        centre = log_meas + (shifts[1] if side == "upper" else shifts[0])
        raw = centre + sign * Z95 * sig_raw
        e95 = None if sig_95 is None else centre + sign * Z95 * sig_95
        binding = None if e95 is None else (e95 < hi_s if side == "upper" else e95 > lo_s)
        out[side] = {"raw_edge": raw, "edge_95": e95, "binding": binding,
                     "shipped_outside_raw": shipped > raw if side == "upper" else shipped < raw,
                     "shipped_outside_95": None if e95 is None else (shipped > e95 if side == "upper" else shipped < e95)}
    return out


def class_stats(values: Mapping[str, Tuple[float, ...]]) -> Dict[str, Dict[str, Any]]:
    out = {}
    for k, v in values.items():
        s = sorted(v)
        n = len(s)
        med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
        out[k] = {"n": n, "median": med, "lo": s[0], "hi": s[-1]}
    return out


def _current_band(row: Mapping[str, Any]) -> Dict[str, Any]:
    b = row.get("band") or {}
    legit = b.get("distribution") == "empirical_measured_class"
    return {"source": b.get("source"), "distribution": b.get("distribution"), "lo": b.get("lo"), "hi": b.get("hi"),
            "sigma": b.get("sigma"), "sampled": bool(b.get("sampled")), "width_dex": b.get("width_dex"),
            "legitimate_under_s4": legit,
            "why": ("a B55 measured-class distribution" if legit else
                    "fixed: no band" if b.get("distribution") == "fixed" or not b else
                    "a fit interval, a search bound or a search bound clipped by the declared prefactor span: "
                    "never a band under S4")}


def _weight(row: Mapping[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    s = row.get("sensitivity") or {}
    return s.get("max_odorant_path_share"), s.get("max_pot_flux_share"), s.get("b38_sensitivity_norm")


def propose(cid: str, spec: Mapping[str, Any], row: Mapping[str, Any], *, classes: Mapping[str, Dict[str, Any]],
            lane_between: Mapping[str, float], dex_per_kj: float, shipped_ea: Optional[float] = None,
            b55_route: Optional[Tuple[str, str]] = None) -> Dict[str, Any]:
    """The S4 proposal for one coordinate. ``shipped_ea``: the constant's shipped barrier (for the joint check of a
    log10 k_ref); ``b55_route``: (route, class) from B55's route map, when the constant is in it."""
    lane = row["lane"]
    kind = row["coordinate"]
    shipped = float(row["value"])
    ident = row.get("identification") or {}
    search = tuple(ident.get("bounds") or (None, None))
    current = _current_band(row)
    srcs = [dict(SOURCES[s], id=s) for s in spec.get("sources", [])]
    between = lane_between.get(lane)
    entry: Dict[str, Any] = {
        "coordinate": cid, "lane": lane, "coordinate_kind": kind, "shipped_value": shipped, "unit": row["unit"],
        "reference_temperature_c": row["reference_temperature_c"], "search_bounds": list(search),
        "current_band": current, "note": spec.get("note"),
        "lane_between_lab_dex": between,
        "outputs": {"species": (row.get("sensitivity") or {}).get("species_of_its_steps", []),
                    "max_odorant_path_share": _weight(row)[0], "max_pot_flux_share": _weight(row)[1],
                    "b38_sensitivity_norm": _weight(row)[2]},
    }
    flags: List[str] = []
    if kind == "ea_kj_mol":
        klass = spec.get("ea_class")
        c = classes.get(klass) if klass else None
        for s in srcs:
            s["evaluation"] = None
        if c is None:
            prop = {"route": "refused", "basis": "no reaction class: the dependent outputs are refused"}
        else:
            prop = {"route": "class_distribution", "class": klass, "lo": c["lo"], "hi": c["hi"], "n": c["n"],
                    "median": c["median"], "width_dex": (c["hi"] - c["lo"]) * dex_per_kj,
                    "basis": f"the empirical distribution of the {c['n']} measured barriers of the {klass} class "
                             f"(src/kinetic_core/ea_class_prior.py, B55): {CLASS_REASON[klass]}; k at the reference "
                             f"temperature held",
                    "shipped_inside": c["lo"] <= shipped <= c["hi"]}
            if spec.get("family"):
                prop["basis"] += f"; one draw for the whole {spec['family']} family (one shared barrier)"
        same_as_current = (current["legitimate_under_s4"] and c is not None and current["lo"] is not None
                           and abs(current["lo"] - c["lo"]) < 1e-6 and abs(current["hi"] - c["hi"]) < 1e-6)
        if prop["route"] == "refused":
            consequence = "refuse"
        elif same_as_current:
            consequence = "unchanged"
        else:
            cur_lo, cur_hi = current["lo"], current["hi"]
            if cur_lo is None or cur_hi is None:
                consequence = "widen"
            elif c["lo"] <= cur_lo and c["hi"] >= cur_hi:
                consequence = "widen"
            else:
                consequence = "shift"
                flags.append(f"narrower than today's band on one side ({cur_lo:.1f}-{cur_hi:.1f} -> "
                             f"{c['lo']:.1f}-{c['hi']:.1f} kJ/mol): an upgrade on that side")
        if cid in CLASS_RULING_COORDINATES:
            flags.append(f"class ruling (S1, both readers): a sulfur-bearing reactant puts it in class {klass} by "
                         f"B55's own definition; B55's route map ({b55_route[0] if b55_route else '?'} -> "
                         f"{b55_route[1] if b55_route else '?'}) conflated the route served with the step's chemistry")
        elif prop["route"] == "class_distribution" and not current["legitimate_under_s4"]:
            if b55_route:
                flags.append(f"B55's route map already assigns it (route {b55_route[0]} -> class {b55_route[1]}), "
                             f"but B55 does not sample decay-family barriers; class agreed by both readers")
            else:
                flags.append("class assignment new in this record (not in B55's route map); agreed by both readers")
        if prop["route"] == "class_distribution" and not prop["shipped_inside"]:
            flags.append(f"the SHIPPED barrier ({shipped:.1f}) lies outside the class's measured range: an S4 case "
                         f"moves the point as well as the band")
    else:
        klass = spec.get("transport")
        ea_band = (classes[klass]["lo"], classes[klass]["hi"]) if klass else None
        uppers, lowers, undefined_95 = [], [], False
        for s in srcs:
            if s.get("heldout") and s.get("usable"):
                raise ValueError(f"{s['id']}: a held-out row is never a band source")
            if not s.get("usable") or s.get("value") is None:
                s["evaluation"] = None
                continue
            factor, factor_note = engine_factor(spec, s["id"])
            ev = evaluate_bound(s, t_ref_c=float(row["reference_temperature_c"]), ea_band=ea_band,
                                between_lab=between, search=search, shipped=shipped, factor=factor)
            ev["engine_factor"] = factor
            ev["engine_factor_note"] = factor_note
            ev["joint_shipped"] = joint_check(s, factor=factor, shipped_k=shipped, shipped_ea=shipped_ea,
                                              t_ref_c=float(row["reference_temperature_c"]), between_lab=between)
            s["evaluation"] = ev
            if "transport" in ev:
                continue
            if between is None:
                undefined_95 = True
            for side, bucket in (("upper", uppers), ("lower", lowers)):
                if side in ev and ev[side]["edge_95"] is not None:
                    bucket.append((ev[side]["edge_95"], s["id"]))
        hi = min(uppers) if uppers else None
        lo = max(lowers) if lowers else None
        if hi and lo:
            route = "measurement_two_sided"
        elif hi:
            route = "measurement_one_sided_upper"
        elif lo:
            route = "measurement_one_sided_lower"
        else:
            route = "refused"
        prop = {"route": route, "lo": None if lo is None else lo[0], "hi": None if hi is None else hi[0],
                "lo_from": None if lo is None else lo[1], "hi_from": None if hi is None else hi[1],
                "transport_class": klass,
                # whether a measured edge lies inside the fit's search box (a band edge outside it is still a band
                # edge: the search bound is never a band, and neither is it a limit on one)
                "lo_inside_search_box": None if lo is None else lo[0] > search[0],
                "hi_inside_search_box": None if hi is None else hi[0] < search[1]}
        if route == "refused":
            usable = [s for s in srcs if s.get("usable")]
            if undefined_95:
                prop["basis"] = ("only a raw measured bound: its 95 % edge is undefined because the lane has no "
                                 "between-lab term (no H0, no pooled class spread in the S1 record); no band, the "
                                 "dependent outputs are refused")
            elif usable:
                prop["basis"] = ("the usable bounds on disk cannot be carried to the reference temperature; no band, "
                                 "the dependent outputs are refused")
            else:
                prop["basis"] = "no measured bound, analogue with a measured law, or class distribution of k on disk; " \
                                "the dependent outputs are refused"
        elif route == "measurement_two_sided":
            prop["width_dex"] = hi[0] - lo[0]
            prop["basis"] = "a measured two-sided band at S2 (b)'s 95 % edges, carried over the class's Ea band"
            prop["shipped_inside"] = lo[0] <= shipped <= hi[0]
        else:
            inside = prop["hi_inside_search_box"] if hi else prop["lo_inside_search_box"]
            prop["basis"] = ("a one-sided measured bound at S2 (b)'s 95 % edge"
                             + ("" if inside else " (outside the fit's search box: it trims nothing of today's "
                                                  "search range)")
                             + "; the other side has no source and is refused")
            prop["shipped_inside"] = (shipped <= hi[0]) if hi else (shipped >= lo[0])
        consequence = {"refused": "refuse", "measurement_two_sided": "widen" if not current["sampled"] or
                       (current["width_dex"] or 0) < prop.get("width_dex", 0) else "narrow",
                       "measurement_one_sided_upper": "one_sided", "measurement_one_sided_lower": "one_sided"}[route]
        joint = []
        for s in srcs:
            j = (s.get("evaluation") or {}).get("joint_shipped")
            if not j:
                continue
            for sd in ("upper", "lower"):
                if sd in j and j[sd]["violates_raw"]:
                    joint.append({"source": s["id"], "side": sd, "temperature_c": j["temperature_c"],
                                  "excess_dex": j[sd]["excess_dex"], "margin_95_dex": j[sd]["margin_95_dex"],
                                  "beyond_95": j[sd]["violates_95"]})
        entry["joint_violations"] = joint
        for v in joint:
            flags.append(f"the SHIPPED (k, Ea) pair is {v['excess_dex']:.2f} dex beyond {v['source']}'s {v['side']} "
                         f"bound at {v['temperature_c']:.0f} C ("
                         + ("beyond the 95 % margin" if v["beyond_95"] else
                            "inside the 95 % margin" if v["beyond_95"] is not None else "95 % margin undefined")
                         + (f", {v['margin_95_dex']:.2f} dex)" if v["margin_95_dex"] is not None else ")"))
    used = [s for s in srcs if s.get("usable") and s.get("evaluation") and "transport" not in s["evaluation"]]
    deciding = [s for s in used if s["id"] in (prop.get("lo_from"), prop.get("hi_from"))] if kind != "ea_kj_mol" else []
    labs = sorted({s["lab"] for s in deciding})
    if deciding and len(labs) == 1:
        flags.append(f"the band rests on a single laboratory ({labs[0]})")
    if deciding and all(s["reader"] == "first" for s in deciding):
        flags.append("the band rests on a first-reader reading")
    correction = SECOND_READING_CORRECTIONS.get(cid)
    entry["second_reading"] = {"date": SECOND_READING_DATE,
                               "verdict": "corrected and applied" if correction else "agreed",
                               "correction": correction}
    is_upgrade = (kind != "ea_kj_mol" and prop["route"] != "refused") or consequence in ("shift", "narrow")
    entry["upgrade"] = is_upgrade
    entry.setdefault("joint_violations", [])
    # S1's two readers: the record row was re-read (all rows, SECOND_READING_DATE); a measured band also needs every
    # deciding source's number read from the paper by two readers (verification file or register); a class band's
    # assignment is agreed by both readers (its members are B55's registry values).
    if prop["route"] == "refused":
        entry["meets_s1_two_readers"] = None
    elif kind == "ea_kj_mol":
        entry["meets_s1_two_readers"] = True
    else:
        entry["meets_s1_two_readers"] = bool(deciding) and all(s["reader"] == "second" for s in deciding)
    if is_upgrade:
        flags.append("an upgrade from no legitimate band: S1's second reader "
                     + ("corrected it (applied)" if correction else "agreed"))
    entry.update({"sources": srcs, "proposal": prop, "consequence": consequence, "flags": flags,
                  "dm": sorted({s["dm"] for s in srcs if s.get("dm")})})
    return entry


def kumazawa_joint_note(rows: Mapping[str, Mapping[str, Any]], between: Optional[float]) -> Dict[str, Any]:
    """Kumazawa 2003 Fig. 3 pH 3: k_fft_decay(121 C) <= -ln(0.98)/10 min, with k(145 C) held at the shipped value."""
    k145 = float(rows["k_fft_decay.log10_k_ref"]["value"])
    bound = math.log10(-math.log(0.98) / 10.0)
    per_kj = arrhenius_shift_dex(1.0, 121.0, 145.0)
    need_raw = (k145 - bound) / per_kj
    need_95 = None if between is None else max(0.0, (k145 - (bound + Z95 * between)) / per_kj)
    return {"source": "kumazawa_fft_ph3", "k_fft_decay_log10_k_145C": k145, "bound_log10_k_121C": bound,
            "ea_needed_raw_kj_mol": need_raw, "ea_needed_95_kj_mol": need_95,
            "reading": "with k(145 C) held, the pH 3 survival ceiling at 121 C needs the thiol-sink barrier >= "
                       f"{need_raw:.0f} kJ/mol at its raw edge (outside every class value; B60 scores the row in "
                       "dex of survival, where the 1.1 dex rate gap shows as -0.11 dex); once the lane's between-lab "
                       "term is added it binds nothing. A joint (k, Ea) constraint for the S4 case to sample within, "
                       "and an S2 (b) record for k_fft_decay.log10_k_ref (not one of the 43)"}


def build() -> dict:
    from src.kinetic_core import ea_class_prior

    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    rows = {r["id"]: r for r in audit["rows"]}
    listed = [e["coordinate"] for e in audit["summary"]["s4_sensitive_unidentified"]]
    if sorted(listed) != sorted(COORDINATES):
        missing = sorted(set(listed) - set(COORDINATES))
        extra = sorted(set(COORDINATES) - set(listed))
        raise ValueError(f"the audit's S4 list and this record disagree: missing {missing}, extra {extra}")
    from src.kinetic_core import parameters_sulfur as ps

    b55_classes = class_stats(ea_class_prior.class_values())
    members: Dict[str, List[float]] = {CARB: [], SULF: []}
    for m in ea_class_prior.measured_activation_energies():
        members[CLASS_RULING_MEMBER_MOVES.get(m.label, m.klass)].append(m.ea_kj_mol)
    classes = class_stats({k: tuple(v) for k, v in members.items()})
    lane_between = {**audit["rules"]["lane_between_lab_dex"], **EXTRA_LANE_BETWEEN_LAB}
    dex_per_kj = float(audit["rules"]["dex_per_kj_mol_100_145C"])

    def one(c: str, spec: Mapping[str, Any]) -> dict:
        constant = c.split(".")[0]
        ea_row = rows.get(f"{constant}.ea_kj_mol")
        route = ps.FORMATION_ROUTE_OF.get(constant)
        b55_route = (route, ea_class_prior.CLASS_OF_ROUTE[route]) if constant in ps.FORMATION_ROUTE_OF else None
        return propose(c, spec, rows[c], classes=classes, lane_between=lane_between, dex_per_kj=dex_per_kj,
                       shipped_ea=None if ea_row is None else float(ea_row["value"]), b55_route=b55_route)

    entries = [one(c, COORDINATES[c]) for c in listed]
    extras = [one(c, EXTRA[c]) for c in EXTRA]
    return {"artifact": "s4_band_sources",
            "provenance": provenance.provenance_block(
                "s4_band_sources", generated_by="scripts/generators/build_s4_band_sources.py",
                inputs=[AUDIT, ROOT / "src" / "kinetic_core" / "ea_class_prior.py",
                        ROOT / "scripts" / "generators" / "build_s4_band_sources.py"]),
            "rules": {"s4": "a band comes from measured bounds, a same-class analogue, or the reaction-class "
                            "distribution of measured constants; a search bound is never a band; with no class, the "
                            "output is refused (sec. 7 S4, 066684b7)",
                      "edge_95": "a bound's 95 % edge: 1.96 x the quadrature of the reported noise, the within-lab "
                                 "between-paper and assay spreads and the between-lab term (the lane's, S2 (b)); "
                                 "undefined for a lane without one",
                      "transport": "a bound at temperature T is carried to the reference temperature over the whole "
                                   "Ea band of the coordinate's class, at its loose edge (k and Ea are sampled jointly)",
                      "reader": "a source's `reader` is the reading of its NUMBER from the paper (second = a "
                                "verification file or register records two readers); every row of this record was "
                                "re-read by the S1 second reader on " + SECOND_READING_DATE + " and its corrections "
                                "are applied (`second_reading` per row)",
                      "joint_check": "the SHIPPED (k, Ea) pair is checked against every bound at the bound's own "
                                     "temperature and pot (the engine's pH factor included); a violation is reported "
                                     "against the raw margin and the S2 (b) 95 % margin",
                      "heldout": "a declared hold-out row is never a band source (coordinator ruling)",
                      "lane_between_lab_dex": lane_between,
                      "lane_between_lab_source": {"sulfur": audit["rules"].get("lane_between_lab_source"),
                                                  **EXTRA_LANE_BETWEEN_LAB_SOURCE},
                      "dex_per_kj_mol_100_145C": dex_per_kj},
            "classes": classes, "b55_classes": b55_classes,
            "class_ruling": "S1 (both readers, " + SECOND_READING_DATE + "), B55's own definition: a step with a "
                            "sulfur-bearing reactant is class sulfur. B55's route map (FORMATION_ROUTE_OF) and "
                            "CLASS_OF_ROUTE conflated two principles, the route a constant serves and its step's "
                            "reactant chemistry: k_ttca_deg, k_pent_tdp, k_pent_dpo and k_dpo_ddp go to sulfur, and the "
                            "class member Zhang 2026 k16 (Cys-Amadori -> alpha-dicarbonyl) moves to sulfur (ranges "
                            "unchanged; n " + f"{classes[CARB]['n']} / {classes[SULF]['n']}" + "). "
                            "src/kinetic_core/ea_class_prior.py is unchanged here: an S4 case applies the ruling",
            "coordinates": entries, "extra_dm_coordinates": extras,
            "joint_constraints": [kumazawa_joint_note(rows, lane_between.get("sulfur"))],
            "audit_record_corrections": [
                {"record": "incumbent_evidence_tiers DM-07 (k_ttca_deg.log10_k_ref)",
                 "issue": "Zhai 2021's 0.0271 is a ZERO-order slope in mmol/L/min (c = 10.331 - 0.0271 t; probably "
                          "0.0274, zhai2021_extraction), not a first-order constant: its first-order equivalent at "
                          "10 mM is ~0.0027 /min. The record reads it as 'first order, k 0.0271 /min'. This record "
                          "uses the early first-order rate at 140 C instead (0.0138 /min)",
                 "effect": "none on the verdict (the bound stays satisfied); the S1 record's DM-07 text should read "
                           "'zero order, 0.0271 (probably 0.0274) mmol/L/min'"}],
            "summary": _summary(entries, extras)}


def _summary(entries: List[dict], extras: List[dict]) -> dict:
    def count(es):
        out: Dict[str, int] = {}
        for e in es:
            out[e["proposal"]["route"]] = out.get(e["proposal"]["route"], 0) + 1
        return out

    def weight(e):
        o = e["outputs"]
        share = max(x for x in (o["max_odorant_path_share"], o["max_pot_flux_share"], 0.0) if x is not None)
        return share if e["lane"] == "sulfur" else min(1.0, (o["b38_sensitivity_norm"] or 0.0) / 10.0)

    order = {"refuse": 0, "one_sided": 1, "widen": 2, "shift": 3, "narrow": 4, "unchanged": 5}
    moved = sorted((e for e in entries if e["consequence"] != "unchanged"),
                   key=lambda e: (order[e["consequence"]], -weight(e)))
    routes = count(entries)
    return {"coordinates": len(entries), "by_route": routes,
            "banded_from_measurement": sum(v for k, v in routes.items() if k.startswith("measurement")),
            "banded_from_class_distribution": routes.get("class_distribution", 0),
            "refused": routes.get("refused", 0),
            "measured_edge_inside_search_box": [
                e["coordinate"] for e in entries if e["proposal"]["route"].startswith("measurement")
                and (e["proposal"]["lo_inside_search_box"] or e["proposal"]["hi_inside_search_box"])],
            "by_consequence": {k: sum(1 for e in entries if e["consequence"] == k) for k in order},
            "extra_dm_by_route": count(extras),
            "no_legitimate_band": [e["coordinate"] for e in entries if e["proposal"]["route"] == "refused"],
            "most_outputs_moved": [{"coordinate": e["coordinate"], "consequence": e["consequence"],
                                    "weight": weight(e), "species": e["outputs"]["species"]} for e in moved],
            "single_lab_or_first_reader": [e["coordinate"] for e in entries
                                           if any("single laboratory" in f or "first-reader" in f for f in e["flags"])],
            "class_ruling_coordinates": list(CLASS_RULING_COORDINATES),
            "joint_violations": [{"coordinate": e["coordinate"], **v} for e in entries + extras
                                 for v in e.get("joint_violations", [])],
            "second_reading_corrected": sorted(c for c in SECOND_READING_CORRECTIONS),
            "meets_s1_two_readers": [e["coordinate"] for e in entries if e.get("meets_s1_two_readers")],
            "band_without_two_source_readers": [e["coordinate"] for e in entries
                                                if e.get("meets_s1_two_readers") is False]}


# ---------------------------------------------------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------------------------------------------------
def _f(x: Any, nd: int = 2) -> str:
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if abs(x) >= 10 ** (-nd) or x == 0 else f"{x:.2e}"
    return str(x)


def _prop(e: dict) -> str:
    p = e["proposal"]
    r = p["route"]
    if r == "class_distribution":
        return f"class `{p['class']}` (n {p['n']}): {p['lo']:.1f}-{p['hi']:.1f} kJ/mol ({p['width_dex']:.2f} dex)"
    if r == "refused":
        return "**no legitimate band: dependent outputs refused**"
    lo = "refused" if p["lo"] is None else f"{p['lo']:.2f}" + ("" if p["lo_inside_search_box"] else "°")
    hi = "refused" if p["hi"] is None else f"{p['hi']:.2f}" + ("" if p["hi_inside_search_box"] else "°")
    return f"{r.replace('measurement_', 'measured ').replace('_', ' ')}: [{lo}, {hi}]"


def _cur(e: dict) -> str:
    c = e["current_band"]
    if not c["source"] or c["distribution"] == "fixed":
        return "fixed (none)"
    rng = f" {_f(c['lo'], 1)}-{_f(c['hi'], 1)}" if c["lo"] is not None else ""
    sig = f" sigma {_f(c['sigma'])}" if c["sigma"] and c["lo"] is None else ""
    return f"{c['distribution']}{rng}{sig}{'' if c['legitimate_under_s4'] else ' (not S4)'}"


def _outputs(e: dict) -> str:
    o = e["outputs"]
    parts = [", ".join(o["species"])]
    if o["max_odorant_path_share"] is not None:
        parts.append(f"odorant {_f(o['max_odorant_path_share'])}")
    if o["max_pot_flux_share"] is not None:
        parts.append(f"flux {_f(o['max_pot_flux_share'])}")
    if o["b38_sensitivity_norm"] is not None:
        parts.append(f"B38 {_f(o['b38_sensitivity_norm'], 1)}")
    return "; ".join(parts)


def _src_line(s: dict, e: dict) -> str:
    ev = s.get("evaluation")
    txt = (f"- `{s['id']}`{' (' + s['dm'] + ')' if s.get('dm') else ''}: {', '.join(s['paper_keys'])}; {s['where']}; "
           f"lab {s['lab']}; {s['conditions']}; reader: {s['reader']}"
           + (f" ({s['reader_note']})" if s.get("reader_note") else "") + ". ")
    if s.get("usable"):
        txt += f"**Used** ({s['relation']}): {s['why']}."
        if s.get("conversion"):
            txt += f" Value {s['value']:.3g} at {s['temperature_c']:.0f} C ({s['conversion']})."
        if ev and "transport" in ev:
            txt += f" {ev['transport']}."
        elif ev:
            if ev.get("engine_factor_note"):
                txt += f" Divided by {ev['engine_factor_note']}."
            j = ev.get("joint_shipped")
            if j:
                txt += (f" Shipped pair at {j['temperature_c']:.0f} C: {j['shipped_log10_k_at_T']:.2f} against the "
                        f"bound {j['bound_log10_k_at_T']:.2f}.")
            for side in ("upper", "lower"):
                if side in ev:
                    d = ev[side]
                    txt += (f" {side} edge at {e['reference_temperature_c']:.0f} C: raw {d['raw_edge']:.2f}, 95 % "
                            f"{_f(d['edge_95'])} ({'binding' if d['binding'] else 'non-binding' if d['binding'] is not None else 'undefined'}).")
    else:
        txt += f"Not used: {s['why']}."
    return txt


def markdown(p: dict) -> str:
    s = p["summary"]
    L = ["# S4 band sources: the 43 sensitive fitted-unidentified coordinates (plan sec. 7, S4)", "",
         "*Generated by `scripts/generators/build_s4_band_sources.py` from the incumbent evidence audit "
         "(`results/validation/incumbent_evidence_tiers.json`, PR #98) and the B55 measured-barrier classes "
         "(`src/kinetic_core/ea_class_prior.py`). Moves no constant; the evidence for a later S4 case. Do not edit "
         "by hand.*", "",
         "**Rule (S4).** " + p["rules"]["s4"] + ". **95 % edge.** " + p["rules"]["edge_95"] + ". **Transport.** "
         + p["rules"]["transport"] + ". **Readers.** " + p["rules"]["reader"] + ".", "",
         f"**{s['coordinates']} coordinates: {s['banded_from_measurement']} banded from a measurement "
         f"(" + ", ".join(f"{k.replace('measurement_', '')} {v}" for k, v in s["by_route"].items()
                          if k.startswith("measurement")) + f"), {s['banded_from_class_distribution']} from a class "
         f"distribution, {s['refused']} refused.** Consequences: "
         + ", ".join(f"{k} {v}" for k, v in s["by_consequence"].items() if v) + ". Of the measured bands, only "
         + ", ".join(f"`{c}`" for c in s["measured_edge_inside_search_box"]) + " have an edge inside the fit's search "
         "box; the other measured edges are legitimate but trim nothing. Every log10 k_ref band would be new (none is "
         "legitimate today: no class distribution of rate constants exists on disk, so a k with no measured bound is "
         "refused); every Ea coordinate has a class (the B55 classes), so none is refused.", "",
         "Classes (after the S1 class ruling): " + "; ".join(
             f"`{k}` n={v['n']}, median {v['median']:.1f}, range {v['lo']:.1f}-{v['hi']:.1f} kJ/mol"
             for k, v in p["classes"].items()) + " (B55 as shipped: " + "; ".join(
             f"`{k}` n={v['n']}" for k, v in p["b55_classes"].items()) + "). **Class ruling.** " + p["class_ruling"]
         + ".", "",
         "**Second reading (S1).** Every row was re-read on " + SECOND_READING_DATE + "; "
         + str(len(s["second_reading_corrected"])) + " rows were corrected and the corrections are applied (column "
         "`flags`, and `second_reading` in the JSON). Rows meeting S1's two-reader requirement (record row re-read, "
         "and every deciding number read from the paper by two readers): "
         + ", ".join(f"`{c}`" for c in s["meets_s1_two_readers"]) + ". Bands whose deciding number still has one "
         "reader: " + (", ".join(f"`{c}`" for c in s["band_without_two_source_readers"]) or "none") + ".", "",
         "**Rule: held-out rows.** " + p["rules"]["heldout"] + ". **Joint check.** " + p["rules"]["joint_check"]
         + ". Between-lab terms: " + "; ".join(f"{k} {v:.3f} dex ({p['rules']['lane_between_lab_source'].get(k)})"
                                              for k, v in p["rules"]["lane_between_lab_dex"].items()) + ".", "",
         "## 1. The proposals", "",
         "| coordinate | shipped | today's band | proposed S4 band | consequence | outputs it moves | flags |",
         "|---|---|---|---|---|---|---|"]
    for e in p["coordinates"]:
        L.append(f"| `{e['coordinate']}` | {_f(e['shipped_value'])} | {_cur(e)} | {_prop(e)} | {e['consequence']} | "
                 f"{_outputs(e)} | {'; '.join(e['flags']) or '—'} |")
    L += ["", "° a measured edge outside the fit's search box (still a band edge; it trims nothing of today's search "
          "range). A measured log10 k_ref edge is at the coordinate's reference temperature, in its unit.", "",
          "**Consequence.** `refuse`: no legitimate band; the outputs whose sensitivity to it exceeds the lane's "
          "within-lab term are refused (today the envelope holds it fixed or samples a non-S4 band, so the shipped "
          "intervals omit it). `one_sided`: the bounded side widens to the measured edge, the other side is refused. "
          "`widen`: a wider band than today's, away from the reference temperature for an Ea (k there is held). "
          "`shift`: the class band is narrower than today's on one side (an upgrade there). `unchanged`: today's B55 "
          "class band is already the S4 band.", "",
          "## 2. The bands that would change the most outputs", "",
          "Ordered refuse, one-sided, widen, shift; within each by the largest odorant-path / pot-flux share (sulfur) "
          "or B38 norm / 10 (other lanes).", "", "| coordinate | consequence | weight | species |", "|---|---|---|---|"]
    for m in s["most_outputs_moved"]:
        L.append(f"| `{m['coordinate']}` | {m['consequence']} | {_f(m['weight'])} | {', '.join(m['species'])} |")
    L += ["", "## 3. Coordinates with no legitimate band", "",
          ", ".join(f"`{c}`" for c in s["no_legitimate_band"]) + ".", "",
          "## 4. The shipped (k, Ea) pairs against the bounds, joint constraints, corrections to the audit record", "",
          "| coordinate | bound | side | at | shipped beyond it by (dex) | 95 % margin (dex) | verdict |",
          "|---|---|---|---|---|---|---|"]
    for v in s["joint_violations"]:
        verdict = ("beyond the 95 % margin" if v["beyond_95"] else "inside the 95 % margin" if v["beyond_95"]
                   is not None else "95 % margin undefined")
        L.append(f"| `{v['coordinate']}` | `{v['source']}` | {v['side']} | {v['temperature_c']:.0f} C | "
                 f"{v['excess_dex']:.2f} | {_f(v['margin_95_dex'])} | {verdict} |")
    L.append("")
    for j in p["joint_constraints"]:
        L.append(f"- **{j['source']}**: {j['reading']}.")
    for c in p["audit_record_corrections"]:
        L.append(f"- **{c['record']}**: {c['issue']}. Effect: {c['effect']}.")
    L += ["", "## 5. The contradicting measurements outside the 43 (DM-03, DM-04)", "",
          "| coordinate | shipped | today's band | proposed S4 band | consequence | flags |", "|---|---|---|---|---|---|"]
    for e in p["extra_dm_coordinates"]:
        L.append(f"| `{e['coordinate']}` | {_f(e['shipped_value'])} | {_cur(e)} | {_prop(e)} | {e['consequence']} | "
                 f"{'; '.join(e['flags']) or '—'} |")
    L += ["", "## 6. Sources, per coordinate", ""]
    for e in p["coordinates"] + p["extra_dm_coordinates"]:
        L.append(f"### `{e['coordinate']}` ({e['lane']}; shipped {_f(e['shipped_value'], 3)} {e['unit']} at "
                 f"{e['reference_temperature_c']:.0f} C; search {e['search_bounds']})")
        L.append("")
        L.append(f"Proposal: {_prop(e)}. {e['proposal'].get('basis', '')}." + (f" Note: {e['note']}." if e["note"] else ""))
        L.append("")
        if not e["sources"]:
            L.append("- no candidate source on disk.")
        L += [_src_line(src, e) for src in e["sources"]]
        L.append("")
    return "\n".join(L) + "\n"


def main() -> int:
    payload = build()
    OUT.write_text(json.dumps(payload, indent=1, default=str) + "\n", encoding="utf-8")
    OUT.with_suffix(".md").write_text(markdown(payload), encoding="utf-8")
    s = payload["summary"]
    print(f"wrote {data_paths.rel(OUT)}: {s['coordinates']} coordinates; {s['by_route']}; {s['by_consequence']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
