# Cysteine fate in sugar + cysteine pots (2026-10-31; read-only, no engine runs)

**Core finding.** The engine runs fed-TTCA pots at Cao's K = 60 L/mol equilibrium (at 10 mM: ~23 % closed, ~7 mM free Cys). Its k_ttca_cys (~0.028 /min) and r_cama_form carry no pH or catalyst term. The 10 mM water pots are **kinetically far from that equilibrium**, and the decisive test is already on disk. In Zhai 2021 Fig. 6 (route_submodel.yml), 10 mM extra xylose RAISES free Cys at 120 min: 0.41 vs 0.33 mM at pH 5.5, and 0.55/0.43, 1.47/1.23, 2.15/1.59 at pH 6/7/8. Any equilibrium model predicts a fall (-26 % at K 60, ~-90 % at K 10^4). So at 100 C, k_f ≲ 0.17 L/mol/min (mine).

Zhai 2019 (water, 83 mM Cys + 165 mM Xyl, pH 7.4) relaxes 14.5 -> 41 -> ~54 mM TTCA in 20 min at 70 C. That gives k_f 0.34-0.44 L/mol/min (mine). Cao gives k_f ≥ 0.44 at 90 C. Same lab, same method: k_f differs ≥ 8-30x (T-corrected) between 83 and 10 mM.

**Data-file conflict:** prereg sec. 13 calls Zhai 2021's pots 0.2 M phosphate. The dossier and fit rows say unbuffered water.

## 1. Hypotheses (P = main cause of Cys-high + TTCA-low)

| # | hypothesis | P | for / against |
|---|---|---|---|
| H1 | Ring hydrolysis, recombination and the Amadori step are **catalysed**, so their pseudo-constants scale with solute/buffer concentration. Anchored on concentrated pots, they are 10-100x too fast at 10 mM. Split: amino-acid self-catalysis 0.30, phosphate 0.10, acid limb of the pH-rate profile 0.15 | **0.55** | For: xylose test; Zhai 2019 vs 2021; de Roos 2005 (pre-equilibrium ≤ 1 h at 50 C at ~2 M; "general-acid catalysed by other amino acids"); Cao's glycine arms; linear free-Cys release. Against: C-S ring-chain opening is fast at room temperature (Schubert 1939, Nagasawa 1981, Roberts 1987), but that is tautomerism, not hydrolysis |
| H2 | K(pH, T) wrong (flat in the engine) | 0.10 | For: Kallen 1971 K_app x20 from pH 5.6 to 7.5. Against: pH 5.5 would need K ≥ 10^4; fails the xylose test |
| H3 | Non-volatile sink: Cys + osone/GO/MGO/furfural adducts, melanoidin | 0.12 (secondary 0.7) | For: Cys-moiety gap (TTCA lost - free Cys - 3-DX - 1-DX) is ~2.3 mM at pH 5.5 and ~8 mM at pH 8 (Zhai 2021 Figs. 4b/5/6); added Cys cuts osones/GO 55-80 % (Zhai 2019 Fig. 3B). Against: it lowers TTCA further, the wrong sign for the joint misfit |
| H4 | Headspace O2 -> cystine | 0.08 (secondary 0.3) | For: the gap grows with pH. Against: the excess is largest at pH 5.5; headspace unstated |
| H5 | Analytical (ELSD low end; Cao dissociation during analysis) | 0.10 | Against: Cao's K is flat (±4 %) while TTCA falls 34 %; an artefact would drift it ~+50 %. (Other causes, e.g. a TTCA ammonium salt or heat-up: 0.05) |

H2S is not the sink: thermal release ~2.7 % per 120 min at 100 C (Zheng & Ho 1994); volatiles < 0.5 % of S (Tai & Ho).

## 2. Discriminating observations

- **H1 vs H2:** the ±xylose contrast. On disk and decisive. The 13C5-xylose label fraction in TTCA (Zhai 2021) would give k_f at 10 mM directly; only fragments are printed, so ask the authors.
- **H1 sub-splits:** Zhai 2019 shows phosphate is not needed. Settling the rest needs a design: TTCA 10 vs 100 mM in water, and 10 mM ± 0.2 M phosphate, at 100 C.
- **H3 vs H4:** Cys-Amadori in the same Zhai 2021 pots (measured, plotted nowhere), plus cystine and total thiol after TCEP. CAMA/TTCA ≤ 0.06 at pH 7.4 (Zhai 2019, 2023c) argues against CAMA. Vessel headspace is unstated (Zhai, Kang); pH drifts 7.0 -> 4.9 (Zhai 2023 FC, Kang).

## 3. Anchors

| quantity | value | conditions |
|---|---|---|
| K_app | 60 (45-75) L/mol | Cao: 90 C, pH 4.2-4.7, 0.2 M Pi, 0.4 M |
| k_r / k_f floors | ≥ 7e-3 /min / ≥ 0.44 L/mol/min | Cao (derived) |
| k_f; K | 0.34-0.44 L/mol/min; ≥ 20-41 L/mol | Zhai 2019: 70 C, water, 83 mM (mine) |
| free-Cys release | 2.5e-4 / 3.5e-4 / ~1e-3 / ~2.5e-3 /min | Zhai 2021: pH 5.5/6/7/8, 100 C, 10 mM |
| TTCA loss (0-40 min) | 0.0017 / 0.0021 / 0.0049 / 0.0164 /min | same |
| CAMA decay | 4.5e-3 - 1.4e-2 /min | de Roos: pH 3-6, ~100 C, 0.5 M Pi |
| Cys consumption | Ea 55.1 kJ/mol, 16 % / 2 h | Kang: pH 7, 10 mM |
| H2S release | Ea 133 kJ/mol | Zheng & Ho: pH 3-9, 80-110 C |

## 4. Papers not on disk

DOIs checked against Crossref unless noted. "Abstract not seen" means the title is the only content verified.

1. Kallen RG 1971, JACS 93:6236-6248, **10.1021/ja00752a040**. Title: Schiff-base mechanism of Cys + formaldehyde thiazolidine formation. Abstract not seen.
2. Fife TH, Natarajan R, Shen CC, Bembi R 1991, JACS 113:3071-3079, **10.1021/ja00008a041**. Title: thiazolidine ring opening AND hydrolysis. Abstract not seen.
3. Pesek JJ, Frost JH 1975, Tetrahedron 31:907-913, **10.1016/0040-4020(75)80099-8**. Title: thiazolidine decomposition in acid and base. Abstract not seen.
4. Schubert MP 1936, J Biol Chem 114:341-350, **10.1016/S0021-9258(20)65226-7**. Title: thiol acids + aldehydes (Schubert 1939's companion). Abstract not seen.
5. de Roos KB 1992, ACS Symp Ser 490:203-216, **10.1021/bk-1992-0490.ch016**. Title: meat flavour from cysteine + sugars. Abstract not seen.
6. Lo TW, Westwood ME, McLellan AC, Selwood T, Thornalley PJ 1994, J Biol Chem 269:32299-32305, **10.1016/S0021-9258(18)31635-1** (web-search record). Reported to contain: MGO reacts rapidly and reversibly with N-acetylcysteine to a hemithioacetal (pH 7.4, 37 C).
7. Wei S, Cui H, Hayat K, Zhang X, Ho C-T 2022, JAFC 70:12164-12171, **10.1021/acs.jafc.2c04961**. Title: TTCA forms from a deoxyosone + exogenous cysteine (the trap). Abstract not seen.
8. Guan Y-G et al. 2012, J Dairy Sci 95:590-601, **10.3168/jds.2011-4860** (PubMed). Reported to contain: phosphate, unlike water, accelerates the reversible early stages by bifunctional proton transfer (glucose + proline).

**B7b proposal.** k_ttca_cys and k_cama_form become k0 + k_aa[amino acid] + k_P[phosphate]. Fit Cao + Zhai 2019 (concentrated) jointly with Zhai 2021 / Kang (dilute), with the ±xylose contrast as a must-pass signature. Add Cys + osone (second order) as the sink. P(Cys excess at pH 5.5 drops ≥ 1 dex) ~0.5.
