# Fit / Hold-out Split Declaration — Kinetic-Core Rebuild (Phase 2)

**Declared: 2026-08-28, BEFORE any parameterization of the newly extracted corpus.**
This document pre-registers the role of every dataset the rebuild may use. It was
drafted by the Wave K3 extraction synthesis and ratified with the four orchestrator
decisions recorded in §5. After this commit, no dataset may change columns without a
dated amendment appended here explaining why — and a dataset that has been fit may
NEVER be promoted to hold-out.

Rules (all verified at declaration time):
1. No dataset appears in both columns; the five within-paper cuts each follow one
   declared axis (§4) and no arm/row appears twice.
2. The 21 frozen hold-out bundles (`evidence_class: external_validation_only`) are
   untouched: 4 matrix-path + 17 maillard-path files under
   `data/benchmarks/external_validation/`. None is moved, re-scoped or promoted.
3. Every module holds out at least one dataset (§3 coverage table).
4. The existing 23-file fit panel in `data/benchmarks/` is unchanged.

## D.3 — MODULE 1: SULFUR FORMATION

| dataset | **role** | one-line reasoning |
|---|---|---|
| Hofmann 1998 T1/T2 **pH 5.0** rows (ribose, glucose, fructose, xylose already split) | **FIT** *(already)* | the corpus's only SIDA absolute anchors at the repo's reference pH; unchanged |
| Hofmann 1998 T1/T2 **pH 3 / pH 7** rows | **HOLD-OUT** *(already frozen)* | the pH axis is the model's weakest (6/10); keeping it out is the only honest test of the pH term |
| Hofmann 1998 T3/T4/T6/T7/T8/T10 (step-level fed-precursor) | **FIT** | step-level constraints are what the network is *for*; they identify individual barriers rather than lumps |
| Hofmann 1998 T2 **dry 180 °C / 6 min** rows | **HOLD-OUT (new)** | a genuinely different regime (silica, 6 min, 180 °C); the only a_w extrapolation test the sulfur branch has |
| **Zhou 2023 Table 1, pH 7 column (ARP and ARP+Cys)** | **FIT (new)** | the fed-Amadori initial condition the corpus has never had; pH 7 is its midpoint and its MFT maximum |
| **★ Zhou 2023 Table 1, pH 6 and pH 8 columns** | **★ HOLD-OUT (new)** | **the strongest new hold-out available.** A model fitted at the pH-7 maximum must *predict* the fall on both sides — and it can only do so if the two-factor law (§B2.4 / §A.3.3) is structurally present rather than fitted |
| **Zhou 2023 Figure 3, Cys+MGO time × pH grid** | **HOLD-OUT (new)** | a completely different feedstock (α-dicarbonyl, not Amadori) and the only time axis; nothing in the fit panel constrains it |
| **Zhou 2023 §A.3.3 cross-system pair (582 vs 665 µg/L)** | **HOLD-OUT (new)** | a *derived consistency* test, not a level; it scores the ARP → MGO flux with no free downstream parameter |
| Whitfield 1999 (fed NF, pH 4.5, 140 °C) | **FIT** | the only independent replication of the Hofmann NF chemistry; needed to pin the NF channel |
| Whitfield 2001 (fed NF, **pH 6.5**) | **HOLD-OUT** | the ≥150× pH collapse is the sharpest single pH prediction the model can be asked to make ⚠️ its H₂S column should NOT carry a standalone row (Methods omits HMF — likely a typo) |
| Cerny 2003 T2/T3 (isotope, 95 °C / 4 h) | **HOLD-OUT** | 95 °C is 50 °C below the fit panel and the route mix is measured to change there; the NF ≤ 7 % ceiling must be **evaluated at Cerny's conditions** |
| Cerny 2004 (in-situ branching 54:46 etc.) | **FIT** | pure ratios, free of the fed-precursor artefact; they identify topology, not magnitude |
| van Seeventer 2001 precursor conversion (55 % / 75 %) | **FIT** | a *reactant*-side constraint; it cannot be traded against any product-side residual, so it adds information without competing |
| Meynier 1995 (~20 directional rows) | **HOLD-OUT** | directional only, and the pH axis again; excellent as a shape test, useless as a level |
| Shu 1988 | **neither** | GC area %, no internal standard; pH-shape source only. Record as a shape constraint (§B2.2), not a scored dataset |

## D.4 — MODULE 2: THIOL CONSUMPTION

| dataset | **role** | one-line reasoning |
|---|---|---|
| Charles-Bernard 2005 T2 ladder, **25 °C** | **FIT** | ten first-order constants spanning 755× at one well-specified condition — the module's calibration backbone |
| Charles-Bernard 2005 T3 attenuation matrix | **FIT** | branch weights, not levels; they identify the channel decomposition the ladder alone cannot |
| Hofmann 2002 Fig. 6 + T2, **30 °C** | **FIT** | the cross-validating partner (agreement within 25 %); fitting one without the other wastes the corroboration |
| **Hofmann 2002 Fig. 1, real coffee brew, 80 °C** | **★ HOLD-OUT** | the only temperature extrapolation available — **and it is the one the model will get WRONG in the informative direction** (the brew is *slower* than the 30 °C models). A pass here would be strong evidence; a fail localises the electrophile-pool depletion term |
| **van Seeventer 2001 Table 1, 50 °C zero-order MFT/FFT** | **★ HOLD-OUT** | a third temperature, a third mechanism, and **zero order** — it tests the *functional form*, not just the magnitude |
| **Zhang 2024 Figs. 2d/e/f dimer & MMFT fractions, 115 °C** | **HOLD-OUT** | the highest temperature available and a channel the 30 °C data explicitly exclude; the sharpest test of the "named channels, not one Arrhenius" architecture |
| Zhang 2024 Fig. 1 (four sulfur additives, redox-state axis) | **FIT** | one condition, four levels — usable to identify the redox term that nothing else in the corpus constrains |
| **Zhou 2023 dimer/monomer fractions, 120 °C, pH 6–8** | **HOLD-OUT** | independent lab; **its agreement with Zhang to 1.3× is the score to beat**, and holding both out keeps that agreement as an out-of-sample fact |
| Hofmann 1996 (6 °C, organic solvent) | **neither** | an artefact-control study in ether/DCM/pentane; ingest only the MFT ≫ FFT **ordering** (§A.4), with Z2's magnitude warning |
| Zhu 2023 kafirin binding constants | **neither** | 53 % ethanol; not a Maillard matrix. Use to bound the reversible-partition share of an observed headspace loss, not as a scored row |

## D.5 — MODULE 3: ACRYLAMIDE / SAFETY

| dataset | **role** | one-line reasoning |
|---|---|---|
| Claeys 2005 T2 **control** system | **FIT** | the cleanest complete two-step set with SEs and a non-isothermal correction, and it self-validates |
| Claeys 2005 T2 **competitor** systems (Gln, Cys, Lys, Ala) | **FIT** | they identify the competition term; without them the control row alone is unidentifiable in a mixture |
| De Vleeschouwer I, glucose non-italic subset | **FIT** | the only bimolecular initiation constant; needed to give the trunk a real second-order form |
| **De Vleeschouwer I, FRUCTOSE and SUCROSE systems** | **★ HOLD-OUT** | a different sugar at the same a_w — the sugar-transfer test — and fructose's HPDs span zero, which the model should reproduce as *low confidence*, not as a fitted value |
| De Vleeschouwer II, **cysteine** system (k_E2, Ea_E2) | **FIT** | the tightest parameters in the corpus; they anchor the thiol/Michael-acceptor family that the sulfur module also needs |
| **De Vleeschouwer II, GLUTAMINE system** | **HOLD-OUT** | it carries the B5.5 sign-crossing (promotion grows with T in liquid, shrinks at a_w 0.92) — a shape no fitted term should be allowed to see first |
| Knol 2005 T1 (30 k + 6 Ea) | **FIT** | the largest self-consistent Arrhenius set for the Asn/Glc trunk |
| **Knol 2010 T2 (7 steps × 5 T)** | **★ HOLD-OUT** | a *third* lab on the same trunk; the only genuine cross-lab extrapolation test the acrylamide module can have |
| Knol 2009 real-food band (9.3 × 10³–2.6 × 10⁴ µg/kg dm) | **HOLD-OUT** | real food, not a model system, and the author himself refused to transfer (B8.5) — exactly what a hold-out is for |
| Knol 2009 **degradation** parameters | **neither** | SD ≥ estimate on every one; unidentifiable (C.6) |
| Quan 2020 | **neither** | no units, no orders (C.5) |
| existing `acrylamide_spi_extrusion_130C_ACSRef3` | **FIT** *(already)* | unchanged |

## D.6 — MODULES 4–8

### Module 4 — trunk / browning
| dataset | role | reasoning |
|---|---|---|
| Martins 2005 T2, all 10 steps + HPDs | **FIT** | the repo already fits this corpus; now with uncertainties |
| **Martins 2005 browning (step 9, ε 0.64)** | **★ HOLD-OUT** | the browning lane has **never** had a parameter; holding it out is the only way to learn whether the trunk predicts colour or merely accommodates it |
| Knol 2005 ε = 282 L mol⁻¹ cm⁻¹ (Glc/Asn) | **HOLD-OUT** | tests the amine-specificity of ε (§A.1) rather than fitting around it |
| Zheng 1994 Tables I/III/V (36 k + 8 Ea) | **FIT** | kinetic reference for cysteine thermolysis and β-elimination; **not a benchmark source** (no absolute concentrations) |
| Hemmler 2018 | **neither** | ESI intensity ≠ concentration; pH uncontrolled and drifting 2–3 units. Ordinal only |

### Module 5 — lipid oxidation
| dataset | role | reasoning |
|---|---|---|
| Frankel 1989, the 26-column distribution at zero additive | **FIT** | the only measured branch distribution; it must replace the shipped `hexanal 0.37` |
| **Frankel 1989, the α-tocopherol arms** | **★ HOLD-OUT** | tocopherol moves total and hexanal-share in **opposite** directions — a two-sided test a fitted split cannot fake |
| **Frankel 1989, the nonanal ABSENCE** | **HOLD-OUT (negative test)** | the model must predict ~zero nonanal from fed linoleate hydroperoxide; a fitted split would never be asked this otherwise |
| existing matrix-path hexanal hold-outs (Bi 2020 ×2) | **HOLD-OUT** *(already frozen)* | unchanged |

### Module 6 — protein binding
| dataset | role | reasoning |
|---|---|---|
| Damodaran 1981 soy K's (stated 100 kDa basis) | **FIT** | the strongest provenance in the batch; the basis is printed, not inferred |
| Andriot 2000 β-lactoglobulin K's | **FIT** | needed for the chain-length slope; label `recovered_by_arithmetic` |
| **Barallat-Pérez 2024 lupin + mucin constants** | **★ HOLD-OUT** | a third protein, a third method, and it carries the **method-dependence** finding (B6.7); fitting it would hide the very effect it proves |
| **Andriot 2000 sensory-intensity arm** | **HOLD-OUT** | the model predicts headspace; whether that maps to perceived intensity is precisely the open question (B6.4/B6.5) |
| Starkenmann 2008 saliva binding | **neither** | stranded, no basis, mechanism unresolved (§A.6) |

### Module 7 — thresholds / matrix correction
| dataset | role | reasoning |
|---|---|---|
| **Zhou 2023 SI Table S2 thresholds (15 compounds, water)** | **reference table, NOT a scored dataset** | thresholds are model *inputs*, not predictions; scoring them would be a category error |
| **Zhou 2023 SI Table S2 OAVs** | **HOLD-OUT (arithmetic check)** | already verified exact on 4 of 15 spot-checks; useful as a regression test on the OAV code path |
| Vega 1994 gelatin ladder (6 compounds × 4 T) | **FIT** *(as lookup-table entries)* | the cleanest matrix threshold set — no thermal step after dosing |
| **Brewer 1995 beef set** | **HOLD-OUT, and RECLASSIFIED `dose_added_pre_cook`** | it is not a threshold in the repo's sense (C.2); holding it out avoids fitting thermal loss into a perception term |
| Tian 2020 milk set | **neither, until the unit is settled** | a literal `?` in the units cell; factor-of-1000 basis risk |
| Guadagni 1963/72 aqueous references | **reference only** | held second-hand through Vega; every ratio in §A.7.5 divides by them |

### Module 8 — extrusion / matrix path
| dataset | role | reasoning |
|---|---|---|
| **Xin 2026 (Food Hydrocolloids) 9-formulation carbohydrate factorial** | **FIT** | the only dosed-precursor factorial in a true HME; it is what the extrusion lane exists to reproduce |
| **Xin 2026 xylose and ribose arms specifically** | **★ HOLD-OUT (carved out of the above)** | these are the two arms that **invert** the pentose ≫ hexose ordering (B1.2). Fitting them would let the model absorb the inversion as a parameter instead of explaining it. ⚠️ **This is the one place the split cuts *within* a paper — the cut is by ARM, and no arm appears twice** |
| Conti 2025 I + II (thiamine, 3 severities) | **FIT** | the only dosed-thiamine extrusion data; needed for the formation term |
| **Conti 2025 hexanal absolute pair (1.95 vs 4.70 µg/g)** | **HOLD-OUT** | the only absolute values in either Conti paper (the IS is hexanal-D12, so hexanal alone is properly quantified) — the strongest single number in the lane deserves to be out-of-sample |
| Guo 2020 retention series | **FIT** | pure retention; it identifies the loss term that Conti's net cannot separate |
| existing matrix-path hold-outs (Li 2026, Liu 2023) | **HOLD-OUT** *(already frozen)* | unchanged |
| **Xin 2026b six-formulation PROTEIN-substitution series** (Food Res. Int. 233, 119010) | **★ HOLD-OUT (new)** | it varies **protein composition at zero added carbohydrate** — precisely the axis the fitted Xin 2026 companion holds constant. Fit the carbohydrate axis, **predict** the protein axis. Same extruder, same settings, same group ⇒ a genuinely controlled cross |
| **Xin 2026b total free amino acids (HILIC-MS/MS, 1613.74 → 3347.27 µg/g)** | **FIT** | the only *validated absolute* measurement in either Xin paper; it constrains the amine precursor pool that every Maillard module divides by, and nothing else in the corpus measures it in an extrudate |
| **Xin 2026 vs Xin 2026b same-sample 10–23× discrepancy** | **⚠️ NOT A DATASET — a CALIBRATION FACT** (§E.2.2, §B9.1) | it must be applied as an `absolute_concentration: false` flag on other rows. **Scoring it would be a category error** |

### Module 9 — ★ NEW: the Z3 block (norfuraneol, thiamine route, furfural sink, browning Ea)

| dataset | **role** | one-line reasoning |
|---|---|---|
| **Bornhorst 2017b Ea(norfuraneol) = 104.9 / 121.1 / 122.3 kJ/mol + z** | **FIT** | the corpus's only norfuraneol Ea; it must go in as a *prior*, not a target, and it must carry the label "approach-to-plateau accumulation in an alkaline mashed-potato gel" |
| **Bornhorst 2017 90 °C matrix pair (egg white vs mashed potato, M-2∞ and k)** | **★ HOLD-OUT** | a *matrix*-transfer test at fixed T and fixed formulation — the model should predict a 3.9× k difference between two food gels it was not fitted to |
| Bornhorst 2017 structural zero (M-2 = 0 with no precursors, all three matrices) | **FIT** | a free, unambiguous zero; fitting it costs nothing and catches sign errors |
| **Bornhorst 2017b `2_R,2_L` a\* row at 100 °C, and the Ea 245.6 it drives** | **neither** | **D = 0.5 min < the 1.75 min come-up time.** Not a kinetic measurement (§F #24) |
| **Cerny 2007b full ternary MFT split (54 : 46)** | **FIT** | the branching ratio the thiamine lane has never had a number for |
| **★ Cerny 2007b the two single-route controls (no-cysteine > 99 : < 1; no-thiamine < 5 : > 95)** | **★ HOLD-OUT** | **the sharpest structural test in the whole split.** A model fitted on the ternary must *predict* both limiting cases; getting 54 : 46 right while missing either control means the routes are wrong and the ratio was fitted |
| **★ Cerny 2007 Table 5 concentration pair (85 : 15 at 1× vs 54 : 46 at 2×)** | **★ HOLD-OUT** | **the single highest-value hold-out row Z3 adds** (B10.1). It scores whether the model's branch fractions *respond to concentration at all* — the exact defect Z1 §E diagnosed. A model with fixed branch fractions fails it by construction, which is the point |
| Cerny 2007 Table 2 pH ladder (65 peak-area cells, 5 pH) | **HOLD-OUT (directional only)** | peak areas, never concentrations; but six sign-level switches and a third MFT-vs-pH shape (B10.10, B10.11) — a shape test the pH term should not see first |
| Cerny 2007 Table 4 isotope splits across pH | **FIT** | fractions, response-factor-immune; they identify which route supplies MFT at each pH, which is a prerequisite for any pH term at all |
| **Yaghmur 2005 water arm (k_obs ×4, Ea 46.50 ± 1.0)** | **neither — `audit_flag` only** | 40–70 °C, and a lump over ≳98.8 % non-FFT flux. **Not comparable to an Eyring barrier** (§C.14, §F #29). Record it; do not score against it |
| Yaghmur FFT-share bound (≲ 1.2 % of the furfural flux) | **FIT (as a ceiling)** | a one-sided constraint on the furfural → FFT branch, and the corpus has no other |
| **Ajandouz 2008 own-measurement Ea set (24 values)** | **FIT (as priors on Ea only)** | glucose-loss and amino-loss Ea transfer to pH 5–7 at ±10 % by the paper's own licence; **no rate transfers, and the browning Ea explicitly do not transfer** (§C.13) |
| **★ Ajandouz caramelisation partition (25–80 % of A₂₉₄; 7–55 % of A₄₂₀)** | **★ HOLD-OUT** | it sizes the amine-independent lane — a quantity the model computes but has never been scored on. Holding it out is the only way to learn whether `structural_zero` is right |
| Hofmann 1996b Sotolon anchors (13.5 / 764.7 / 273.1 µg, pH 3/5/7) | **FIT** | the corpus's first Sotolon numbers; nothing else constrains that node |
| **Hofmann 1996b oxidant series (argon 1.7 / air 9.8 / Cu²⁺ 20.1, 38.4 µg)** | **★ HOLD-OUT** | the only required-oxidant measurement in the corpus, and a 5.8× effect the model currently cannot express at all. ⚠️ ingest the second Cu²⁺ dose as "a higher level" (§F #19) |
| Hofmann 1996b AT time × temperature surface (Fig. 2, 30 points) | **HOLD-OUT** | `digitised_from_figure`, ±0.3 µg — but it is the corpus's **only** measured T × t surface with a **sign reversal**, and B10.4's Ea inequality is derived from it, so scoring it would be circular if it were also fitted |
| **Hofmann 1996b Tables 1 and 2** | **neither — DUPLICATES** | verbatim re-publications of Hofmann 1998 Tables 2, 8, 10 (§F #16). **Double-counting hazard** |
| Amrani-Hemaimi 1995 Table 2 isotope fractions (40 cells) | **FIT** | fractions, response-factor-immune, and the only pyrazine carbon-origin bookkeeping the corpus has |
| **Amrani-Hemaimi 1995 the alanine-vs-glycine ethyl-pyrazine ON/OFF switch** | **★ HOLD-OUT** | 20 / 19 % with alanine, **0 / 0 %** with glycine, and 100 % ¹³C-labelled ⇒ *"one single reaction route exists"*. An on/off amino-acid switch is unfakeable by a fitted continuous term |
| Amrani-Hemaimi Table 1 row 6+7 (co-elution) | **neither** | not a single-species value, and its labelling figure is an average of two compounds (§F #28) |
| **van Boekel 2005** | **neither** | one-page poster abstract, zero parameters (§C.16) |
| van Seeventer Z3 addendum §2/§3 tables | **neither** | arithmetic exercises run to prove an extrapolation cannot be made (§C.17, §F #30) |

## D.7 — COVERAGE CHECK

| module | ≥1 hold-out? | which |
|---|---|---|
| 1. sulfur formation | ✅ | Zhou pH 6/8 · Zhou Fig. 3 · Hofmann dry-180 · Whitfield 2001 · Cerny 2003 · Meynier (+ 5 frozen Hofmann rows) |
| 2. thiol consumption | ✅ | Hofmann 80 °C brew · van Seeventer 50 °C · Zhang 115 °C · Zhou 120 °C |
| 3. acrylamide / safety | ✅ | De Vleeschouwer fructose+sucrose · De Vleeschouwer glutamine · Knol 2010 · Knol 2009 real food (+ 12 frozen mp_holdouts) |
| 4. trunk / browning | ✅ | Martins step 9 + ε · Knol ε |
| 5. lipid oxidation | ✅ | Frankel tocopherol arms · Frankel nonanal absence (+ 2 frozen Bi 2020) |
| 6. protein binding | ✅ | Barallat-Pérez lupin + mucin · Andriot sensory arm |
| 7. thresholds / matrix | ✅ | Brewer beef (reclassified) · Zhou OAV arithmetic |
| 8. extrusion / matrix path | ✅ | Xin xylose + ribose arms · **Xin 2026b six-formulation protein series** · Conti hexanal pair (+ 2 frozen) |
| **9. Z3 block (norfuraneol / thiamine route / browning Ea)** | ✅ | **Cerny 2007 concentration pair · Cerny 2007b two single-route controls · Bornhorst matrix pair · Ajandouz caramelisation partition · Hofmann 1996b oxidant series + Fig. 2 surface · Amrani-Hemaimi on/off switch** |

**Disjointness check: no dataset appears in both columns.** The five places the split cuts inside
a single paper are stated explicitly and each cuts along a declared axis:

| # | paper | axis of the cut | fit side | hold-out side |
|---:|---|---|---|---|
| i | **Hofmann 1998** | **pH** | pH 5.0 rows | pH 3 / pH 7 / xylose (already frozen) |
| ii | **Zhou 2023** | **pH column** | pH 7 | pH 6 and pH 8 |
| iii | **Xin 2026** | **carbohydrate arm** | 7 sugars + control | **xylose, ribose** |
| iv | **Cerny 2007 / 2007b** | **system composition** | full ternary + isotope splits | **the two single-route controls, and the 1×/2× concentration pair** |
| v | **Bornhorst 2017 / 2017b** | **what is scored** | the Ea (as a prior) and the structural zero | **the 90 °C matrix pair (egg white vs mashed potato)** |

Convention (i) is the repo's existing one; (ii) mirrors it deliberately; (iii)–(v) each hold out
the *limiting case* or the *transfer axis* while fitting the interior, which is the only
arrangement under which a pass is evidence of structure rather than of fitting.

## §5 — Orchestrator decisions (ratified 2026-08-28)

The four questions flagged by the drafting wave, resolved:

1. **Zhou 2023 pH cut stands, with a diagnostic guard.** The pH labels are initial
   (unbuffered; pH 6 and 7 runs converge to within 0.2 units by end of heating), so
   the pH-6 hold-out column is scored DIAGNOSTIC, not gating, until the model carries
   a pH-trajectory state. The pH-8 column gates normally.
2. **Knol 2010 is split by step, not held out whole.** The acrylamide steps are
   HOLD-OUT (the cross-lab trunk test); the organic-acid and isomerisation steps
   (isomerisation 61±8, acetic 75±10, formic 84±14 kJ/mol) are FIT as Module 4
   (trunk) constants. The two roles live in different modules, so rule 1 holds.
3. **The alkaline Z3 block enters as priors only.** Every Bornhorst/alkaline-source
   Ea is ingested with a mandatory `pH_of_measurement` field and
   `rate_transfer: not_licensed`. The Bornhorst 90 °C matrix pair and the Ajandouz
   caramelisation partition are scored DIAGNOSTIC until a pH-dependent barrier
   exists. Only the parameters with an out-of-sample referee may be described as
   validated; the rest must be labelled unvalidated priors.
4. **The two `diagnostic_only` internal pilots keep their exact current role.**
   No quiet promotion.

Source of record for the underlying extraction evidence:
`data/lit/extraction_dossiers/k3_final_parameter_inventory.md`, whose §A–§C carry
per-row units, conditions, anchors and provenance verdicts, with the K1/K2 and
per-paper dossiers alongside it in the same directory.

## Amendment 1 — 2026-08-28 (research round 2, before any use)

New sources discovered after the original declaration; roles assigned BEFORE any
wave parameterizes them. Direction of change is additive only — no existing row
moves columns.

| dataset | role | reasoning |
|---|---|---|
| Martins 2003 Wageningen thesis (edepot.wur.nl/121418), Table 4.2.3 reverse rates (Amadori -> parent sugar, 100/120 C x pH 5.5/6.8) + Table 5.1 model discrimination (SB <-> DFG reversible, dAICc 276 over irreversible) | **FIT** | parent document of the already-FIT Martins 2005 dataset; same system, same lab; closes the declared reverse-rate structural gap with measured values |
| Martins 2003 thesis Table 4.1.1 glycine-release yields (65-95%) | **FIT** | model-free reactant-side constraint, same system |
| Martins 2003 thesis Table 3.3.1 melanoidin epsilon (1.00 +/- 0.03 at 420 nm, 0.65 +/- 0.02 at 470 nm; T- and DP-invariant) | **FIT (measured input, not a scored target)** | epsilon is an observable-model constant; the Martins step-9 browning HOLD-OUT is unaffected and stays held out |
| Martins 2003 thesis Ch. 6 pH ladder (k at pH 4.8/5.5/6.0/6.8/7.5 + fitted pH exponents) | **HOLD-OUT except pH 5.5/6.8 columns (which are the already-FIT conditions)** | the model is pH-fixed at 6.8; the off-pH columns are exactly the extrapolation test the trunk lacks. Cut axis: pH, mirroring the Hofmann convention |
| Cai 2024 (10.1016/j.foodres.2024.114591) 2-AP multiresponse 100/120/140 C | **role deferred until retrieved** | abstract-verified only; structural facts usable now: glyoxal rate-determining, non-monotone optimum ~100 C |
| Weykamp & Penders 1982 Schiff k-1 = 0.435 /h at 37 C (biomedical) | **FIT (as a prior, rate_transfer: not_licensed to food T without the Ea)** | only measured Schiff reverse rate in any literature; internal k1/k-1/K consistency verified to 4 digits |
| Ge & Lee 1997 reverse-rate Ea | **quarantined until tables retrieved** | method disputed in print by van Boekel; abstract has a unit defect (kcal/J confusion) |

## Amendment 2 — 2026-08-28 (orchestrator ruling on a discovered disjointness violation)

Wave B2 found that Cerny 2007 Table 4's MFT column (declared FIT) and Table 5's
1x arm (declared HOLD-OUT) are the SAME measurement (85:15) published twice —
violating rule 1. RULING: the FIT column shrinks. Cerny 2007 Table 4's MFT
column is excluded from fitting; topology identification is carried by Cerny
2007b's ternary split and Table 4's non-MFT columns. Table 5's concentration
pair keeps its full hold-out role. This is the conservative direction (nothing
held out was ever fit).

Disclosure, for the record: the K3 inventory document itself prints the Zhou
pH-6/pH-8 columns and Zhang Fig. 2d/e values, so Wave B2's builder saw those
hold-out values during directed reads. They entered no parameter, bound, or
initialisation (enforced by a literal-grep firewall test), but "seen" is
recorded honestly here. The frozen mp_holdout_* bundles were never opened.

## Amendment 3 — 2026-08-28 (correction to Amendment 1, from the Martins-cluster extraction)

Direct extraction of the four Martins papers found Amendment 1 partly unsupported.
Corrections, effective immediately and BEFORE any wave uses the affected rows:

1. **FROZEN: the "Martins thesis Table 4.2.3 reverse rates (Amadori -> parent
   sugar), FIT" row.** No numeric reverse rate constant exists in Part I, Part II,
   the epsilon paper, the pH-ladder paper, or Food Chem. 90. The authors DELETED
   DFG -> parent sugars as quantitatively unimportant; parent-sugar re-formation is
   modelled as FORWARD steps out of the Schiff intermediate (Part II Table 3
   k10/k11/k16). The row stays frozen until the thesis PDF is read directly.
2. **Corrected attribution of the model discrimination:** the reversibility result
   is Food Chem. 90 (martins2005.pdf) Table 1, not "thesis Table 5.1": winning
   hypothesis A (Glu+Gly -> E1 <-> DFG) beats hypothesis B by dAICc 276 and the
   FULLY IRREVERSIBLE variant by 287.46. The structural claim (Schiff <-> Amadori
   reversible) STANDS; the citation and the "Amadori -> parent sugar" framing were
   wrong.
3. **Same-experiment de-duplication (binding):** thesis chapters and the four
   journal papers are the SAME experiments published twice. The journal-paper
   versions are canonical; thesis-based rows in Amendment 1 that duplicate them
   (epsilon Table 3.3.1; pH ladder Ch. 6; glycine yields Table 4.1.1) are hereby
   MERGED into the paper rows, declared once. Also binding: Martins Part I and
   Part II are one experiment (fit either the time courses or the derived rate
   constants, never both); martins2005b Tables 1/2/4(1:1) print one experiment
   three times — it counts once.
4. **Roles adopted for the four papers as proposed by the extraction dossiers**
   (Part II Table 3: pH 6.8 columns FIT / pH 5.5 columns HOLD-OUT, quarantining
   k5, k8, k12, k13, k14, k15; Part I: pH 6.8 FIT / pH 5.5 HOLD-OUT, Fig. 5
   diagnostic_only; epsilon: measured input + induction times as HOLD-OUT;
   pH ladder: Table 2 pH 5.5+6.8 FIT once, pH 4.8/6.0/7.5 HOLD-OUT, Table 4
   ratio arms HOLD-OUT as a ratio-invariance test, Table 1 pH-stat-vs-drifting
   as a methodological HOLD-OUT).
5. **Recorded hazard:** Martins 2003 a/b/c letter suffixes are permuted
   differently between the authors' own reference lists — resolve Martins 2003
   citations only by journal + volume + page, never by letter.

## Amendment 4 — 2026-08-28 (roles for the research-round-2 corpus, waves K4a/K4b/K4c)

Roles assigned BEFORE any wave parameterizes these sources; per-dataset detail
lives in each dossier's final section (data/lit/extraction_dossiers/), adopted
as proposed except where amended below. Flagship assignments, binding:

FIT (priority order): Pereyra Gonzales milk-powder k set (15 rate constants,
real food, twice-verified); Meynier partition RATIOS + 17/20 enthalpies (never
its absolute K_aw — static-headspace absolute-scale suspect, spread 9.5x);
Miao's 4 Ea; Leksrisompong K ratios; Lievonen pH-drift values; Stack 2018
NAC arm forward+reverse constants (CORRECTED set — the paper's printed
activation parameters carry a spurious ln(10); dossier's refit is canonical);
li2016 k2 ladder; Cai 2024 100 C course + Ea pair; kumazawa pH grid at 121 C;
Bell 1996 17 k with the 2.4x deconfounded Tg effect (not the paper's "7-fold");
weykamp Schiff set as prior (pH never stated — flag carried); Kang 2026 Table 1
Tier A columns ONLY (mu-g/L basis per the 900-dpi verification, NOT mg/L)
— Tier B semi-quant is never scorable as levels.

HOLD-OUT (gating unless noted): Hong 2020's 10 paired soy/water threshold
ratios — pass criterion: >=7/10 within 5x AND correct sign on all 10 including
the ethyl-4-methylpentanoate inversion; Cai 2024's 120/140 C 2-AP courses (the
model must predict the sign reversal); Stack's GSH arm (nucleophile-identity
transfer); sun2019 pH-9 column (temperature-ordering inversion); Zhou 2025
50 C row (interpolation test; its printed FFT Arrhenius is REFUTED — fit to
methanethiol data — the dossier refit Ea 7.92 kJ/mol is the reference and
means ambient fade is transport-limited).

NEITHER: anantharamkrishnan2020 (mechanism_reference: no rates, saturating
dose; its pH-3 adduct gate is a structural constraint only); bagiyan2004
(initial rates the authors disclaim as qualitative); Zhou 2025 OAV table
(irreconcilable with its own concentrations); Pereyra Gonzales Ea=498 kJ/mol
(viscosity transition, must never enter an Ea prior); Kang Fig 1a sulfur bars
at 120/140 C (corrupted — near-identical to the pH series).

RULING — figure-digitised data: may hold GATING hold-out status only when the
digitisation is cross-validated against text-quoted values from the same paper
(Nakamura qualifies: 14/14 reproduced); otherwise DIAGNOSTIC only.

DEFERRED: Kang 2026 per-compound temperature ladder (Tables S4/S5 are SI-only,
not yet retrieved) — the sulfur T-ladder gap stays OPEN until the SI arrives.
REOPENED: the Hofmann 1998 temperature/water-content series does NOT exist
(two fully-confounded process conditions only) — Amendment 1's hope withdrawn.

PROCESS HAZARD (binding on all future ingestion): pdftotext renders mu as "m"
in RSC/Arbortext PDFs — a silent 1000x unit corruption. Any value from such a
PDF must be verified against a raster before ingestion.

NEW PROVENANCE FIELDS required on matrix records: cross_study_cross_method,
dispersion_scale, absolute_scale_suspect, pH on aldehyde-binding records,
response_variable, dose_saturating.

## Amendment 5 — 2026-08-28 (Kang 2026 SI + adduct survey, wave K4d)

**The sulfur temperature ladder EXISTS.** Kang 2026 Table S4 gives Tier-A
absolute MFT and FFT at 100/120/140 C (pH 7, 120 min, calibration R^2 > 0.998,
subtotals closing to +/-0.003). Roles, binding:
- FIT: MFT/FFT/furfural at 100 and 120 C; Ea(free-cysteine depletion) = 55.1
  kJ/mol from digitised Fig S4 (R^2 0.994); TTCA->Cys yield ceiling 16.3 mol%
  as a one-sided bound.
- ★ HOLD-OUT (gating): MFT and FFT at 140 C — a true extrapolation test where
  a single-Ea model fitted on 100/120 under-predicts 3.8x/2.5x. The measured
  behaviour is NON-ARRHENIUS (apparent Ea rises 7->98 kJ/mol across the legs
  while the sulfur class as a whole decelerates): free-thiol formation switches
  on between 120 and 140 C. A single-Ea sulfur branch is expected to fail this
  hold-out; passing requires the switch-on to emerge structurally.
- QUARANTINED: Table S5's pH-7 column (duplicates S4's 100 C column with
  different SDs — one of the two is mis-pasted); all Kang SDs (use means with
  the dossier's replacement uncertainty); Tier-B polysulfide rows as levels;
  the main paper's Fig 1a sulfur bars at 120/140 C (confirmed corrupted, true
  values recovered in the dossier); the main paper's H2S claim (nothing in the
  SI measures H2S).
- Unit authority: the SI (true-mu Aspose typesetting, raster-verified) outranks
  the mu->m-corrupted main PDF.

**anantharamkrishnan2020b (the parent adduct survey):** hexanal DOES form
covalent lysine adducts on beta-lactoglobulin — the mechanism behind the
matrix hexanal over-prediction is demonstrated but not sizable here
(saturating dose, no rates, no y-axes). NO FIT ROWS. Seven ordinal/binary
gates adopted as proposed, including the 32-of-47 no-adduct negative gate
(falsifies any generic protein-binding loss term). New required channel for
the sulfur branch: FFT is consumed by protein disulfides (1:1 and 2:1 adducts)
— a protein-matrix thiol sink the model does not yet have. The derived
ambient-rate bound (t1/2 37-760 days at 3 wt%) says the covalent channel
CANNOT explain ambient losses and bites only at process temperature.

## Amendment 6 — 2026-08-29 (Meynier 2004, wave K4e — the covalent-sink verdict)

meynier2004 (Int. Dairy J. 14:681, figure-only paper, Fig. 6 digitised at
+/-0.1 mol%): NO FIT ROWS — third consecutive covalent-adduct paper to yield
none. Ten hold-out gates and twelve do-not-use entries adopted as proposed in
the dossier §13. Binding rulings:
1. The covalent aldehyde-protein channel is now BRACKETED two-sided: K4d's
   lower bound and this paper's upper bound overlap at t1/2 ~ 37-74 days
   ambient (two labs, two decades, two protein families, two techniques).
   Carry k2(hexanal-Lys) <= 2.5e-5 M^-1 s^-1 at 20 C as a ceiling with the
   t-2-hexenal pair (5.3-7.9e-5 Lys / 1.5-2.4e-4 His) as measured analogues.
2. VERDICT, computed not asserted: covalent + reversible binding does NOT
   close the hexanal matrix gap — covalent supplies ~0.06% of the 1304x
   log-shift on a threshold-panel timescale; reversible ~25%. The matrix
   module's ratio-reporting design stands. The one surviving corner: covalent
   extent is uncapped (Lys in vast excess), so it becomes relevant at process
   temperature ONLY if Ea >~ 70 kJ/mol — and that Ea is UNMEASURED in every
   corpus source. NEW NAMED WET-LAB GAP: Ea of aldehyde-lysine adduct
   formation on food proteins.
3. K4d's binary gate G-2 is REPLACED by the ordinal form (hexanal DOES
   cross-link at high dose; enals do so at >=4x lower dose).
4. Derived covalent/reversible split adopted as a constraint: hexanal binding
   >=98% reversible on headspace timescales, ~2/3 reversible for t-2-hexenal;
   ~22-33% of Meynier2002's "partition coefficient" is irreversible chemistry
   (sizes that dossier's M-4 flag).
5. Do-not-use: the paper's abstract partition claim (imported from Meynier
   2003 — double-count hazard); its Fig. 2 Trp-fluorescence dataset (authors'
   own disclaimer); fluorescence as a covalent-extent proxy anywhere
   (measured 7x anti-correlated in one experiment here).

## Amendment 7 — 2026-08-29 (pH-trajectory state, pre-registered before Wave B2.2)

For the pH-trajectory state (owner-approved option (a)): Zhou 2023's three
measured final-pH endpoints (initial 6/7/8 -> final 3.22/3.42/5.07, unbuffered,
120 C/2 h) are declared FIT as calibration anchors for the drift model (acid
production from the tracked FA/AA pools + textbook buffer capacity). They are
process-state observables, not volatile levels — no volatile hold-out is
affected. The Chang water-vs-buffer arm pair REMAINS HOLD-OUT (it is the
drift model's exam). Kang decay-barrier identification uses ONLY the declared
Kang 100/120 C FIT rungs; the 140 C column and the Yiltirak/van Seeventer/
Zhang/Zhou-dimer rows remain untouchable hold-outs.

## Amendment 8 — 2026-08-29 (Blank & Fay 1996, DMHF route)

blank1996 (JAFC 44:531, one condition, no yield tables): route-topology source
only. Roles as the dossier proposes: isotopomer branch fractions (4/12/84 Gly;
15/47/38 Ala; ~70/30 Strecker-vs-fragment) = PRIOR on the C1-donor split;
alanine/glycine HEMF truth table + norfuraneol >> DMHF ordering = structural
HOLD-OUTS; the "low mg/kg" prose estimate = DO-NOT-INGEST (no method/basis —
the cys_ribose_140C failure class); NO barrier may ever be fit to this paper
(prohibited derivation). The DMHF channel remains topologically enabled but
uncalibrated until the quantitative companion (Blank, Fay, Lakner & Schlosser
1997, JAFC 45:2642, isotope dilution assays) is retrieved.

Recorded defects for the code-quality pass: the repo's "the amino acid is the
reductant (Blank & Fay 1996)" claim (barrier_constants.py:325 and 4 other
sites) misquotes the paper (it names dismutation OR enoloxo compounds); the
ledger's "stoichiometry the model uses" claim is false (the model lacks the
step); the repo's Strecker-aldehyde-leaves co-product contradicts the paper's
consumed-into-ring route.

## Amendment 9 — 2026-08-29 (pre-declaration of Wave B2.3, before it runs)

1. CHARGE-CONSERVATION FIX (invariant, not a parameter): B2.2 found that
   solute charge is not conserved — sink steps delete carboxylate-bearing
   species without depositing the acid equivalent (diagnosis §3 names every
   site). B2.3 will fix ALL such sites, refit on the unchanged FIT rows, and
   blind re-sit the full panel + exam. Trigger disclosure: the defect was
   noticed via hold-out scoring, so the Kang 140 C MFT row is demoted to
   seen-diagnostic PERMANENTLY (its B2.2 nominal pass was ruled not-counted);
   the fix itself is a conservation law, defensible with no reference to any
   hold-out value.
2. BUFFER-FIELD COMPLETION: the 21 frozen bundles record no buffer species/
   concentration; the exam therefore scores buffered experiments as water
   (measured 11x predicted-FFT sensitivity). B2.3 will complete the CONDITION
   records (buffer identity, concentration, initial pH scale) from each
   bundle's SOURCE PAPER with per-field provenance notes — measured values,
   compound lists, and roles remain byte-identical. The exam is then reported
   BOTH WAYS (buffer-completed and as-was) in the same artifact, permanently.

## Amendment 10 — 2026-08-29 (Wave B6, the lipid-oxidation module: disclosure + two findings)

Written by the build wave itself, additive only. **No dataset changes columns.**

1. **EXPOSURE DISCLOSURE — the α-tocopherol arms are now `seen_diagnostic`,
   permanently.** Frankel 1989 prints the zero-additive column (FIT) and the
   α-tocopherol / 1,4-cyclohexadiene columns (HOLD-OUT) **in the same table
   rows** of Tables 1 and 2, and states the hold-out result in prose in its
   ABSTRACT. Extracting the FIT column from `data/articles/frankel1989.pdf`
   without seeing the hold-out column is not possible. B6's builder saw them.
   Under the Amendment 9 clause 1 precedent (the Kang 140 °C row) the arms are
   demoted to **seen-diagnostic and may never gate**. The mitigation adopted is
   STRUCTURAL rather than procedural: B6's hydrogen-donor term has **no fitted
   parameter and no stored value**, so its two-sided prediction is a
   monotonicity theorem over the whole donor range and there was nothing to tune
   toward the seen numbers even in principle. A literal-grep firewall test
   (`tests/unit/test_kinetic_core_b6.py`) asserts that no hold-out-only literal
   appears in the lipid package or in the frozen fit report.
   **Unaffected and still GATING:** the Frankel **nonanal absence** — a
   structural zero fixed by molecular topology, not by any number in the paper.
   **Never opened:** the two frozen Bi 2020 bundles, reached only through
   `scripts/generators/generate_cutover_final_exam.py`.
   The 1,4-cyclohexadiene arms, unassigned in D.6, are hereby recorded as
   **HOLD-OUT (seen_diagnostic)** — the conservative direction; nothing held out
   was ever fit.

2. **DEFECT FIXED: the cutover exam never scored the four matrix_path bundles.**
   `generate_cutover_final_exam.py` read a measured value from
   `holdout_targets[c].target_value` and then from
   `reference_volatiles[c].conc_ppb`. The four matrix_path bundles carry neither:
   their value is in `measured_volatiles[c].conc_ppb`. Every matrix_path point
   therefore returned `measured = None` and was reported **with no fold error in
   BOTH lanes** — the two scorers shared the bug. It went unnoticed because the
   core REFUSED all four bundles before B6, so a missing referee cost nothing
   visible. This is a schema-reading fix with the same standing as Amendment 9
   clause 2's buffer-field completion: defensible with no reference to any
   hold-out value. **Consequence: the OLD lane's reported numbers move too**, and
   its matrix-path performance is out of sample for the first time.

3. **REPORTED FOR ORCHESTRATOR RULING — two benchmark-record problems B6 found
   and did NOT act on unilaterally:**
   * `bi_2020_raw_pea_hexanal` and `liu_2023_ppi_offnote_baseline` record
     **40 °C for 10 min**. That is an ambient headspace measurement of an
     as-received isolate, so the measured hexanal is the isolate's **accumulated
     storage oxidation**, not a process output. Scoring a formation model
     against an inventory measurement is a category error of the same class as
     Brewer 1995's `dose_added_pre_cook` reclassification (D.6 Module 7).
     Candidate remedy: reclassify both as `inventory_not_process` and score them
     only against a model that carries an ambient-storage formation term.
   * Those two bundles declare **identical** conditions (40 °C, 10 min, pH 6.0,
     a_w 0.95) and **identical** precursors (`Pea Protein Isolate`), yet their
     measured hexanal differs **9.0×**. Whatever separates the two samples is not
     written down, so no model reading only the recorded fields can distinguish
     them. This is a record defect, not a model defect.

4. **New declared gap, named:** *the rate of lipid-hydroperoxide decomposition at
   cooking temperature.* B6 carries it as a bounded, labelled INPUT (Schroën
   2022's 25 °C, hand-fitted, no-SD, lumped `k4`, with the authors' own Q10 =
   2–3 as an explicit assumption with a band). No Q10 number is baked into any
   stated constant, and every lipid prediction returns
   `in_envelope_extrapolated`. Re-affirms `k3` §C.9 and research-round-3 §F.3.

## Amendment 11 — 2026-08-29 (orchestrator rulings on Wave B6's flags)

1. The two 40 C / 10 min matrix-path bundles (bi_2020_raw_pea,
   liu_2023_ppi_offnote) measure the AMBIENT HEADSPACE INVENTORY of an
   as-received isolate, not the output of a thermal process — scoring a
   formation model against them is the Brewer-1995 category error class.
   RULING: they stay in the exam table (frozen, untouched) but carry
   category_mismatch: inventory_measurement and are excluded from gating
   aggregates, reported separately. The two 160 C process points (3.7x, 8.7x,
   both inside interval) are the lane's real exam and gate normally.
2. RECORD DEFECT, no edit: bi_2020_raw_pea and liu_2023_ppi_offnote declare
   identical conditions and precursors yet differ 9.0x in measured hexanal —
   consistent with the corpus's same-sample 10-23x dispersion finding; the
   pair is itself evidence for interval-only absolute reporting.
3. Ratified: Wave B6's self-recorded Amendment 10 (Frankel tocopherol arms
   demoted to seen_diagnostic — the paper prints FIT and hold-out columns in
   the same rows); the two-sided signature HOLDING for every donor split in
   (0,1) is recorded as structural corroboration, not a scored pass.
4. Adopted finding: the historic 36x/1304x hexanal OVER-prediction was a
   FAST-lane artifact (unbounded extrapolation from the invented 0.37 share),
   not lipid chemistry — the measured-slate lane UNDER-predicts. The matrix
   gap's sign was the old lane's, not the food's.

## Amendment 12 — 2026-08-29 (K5a/K5b roles: HMF + DMHF clusters; corrections to Amendment 8)

Per-dataset detail in k5a_hmf_synthesis.md / k5b_dmhf_synthesis.md, adopted as
proposed except as amended. Binding highlights:

DMHF (K5b): blank1997's 39 SIDA cells = FIT (pentose formation edge — the
channel's calibration). CORRECTIONS to Amendment 8's blank1996 roles (a GC-O
dash is a non-detection, not a zero): the HEMF alanine switch is demoted to a
~10-25x preference test; "DMHF from pentose alone" is SIGN-REVERSED into a
ceiling (<0.01 ug/mmol). Corroborated: 70/30 Strecker split (73/27 by the
non-isotopic null control). Wang & Ho: MGO edge structurally separated
(no shared intermediate with the intact-skeleton edge; resolves Poisson's
open b-vs-c); its level is digitised — prior only. NAMED HOLD-OUTS:
apriyantono1993 held-vs-drifting pH pair scored as ONE paired log-ratio test
(the corpus's only pH-trajectory pair — overturns round-3's "ratio-only"
verdict: Table 1 is absolute nmol/mol xylose, and it CLOSES the declared
no-pentose-furfural-yields gap, 274x pH effect); shu1988-vs-wang2008 paired
sink/net pH test. PROHIBITED: any thiol_addition_dmhf fit to Shu's 6.0% GC
area. Edge C (cysteine sink) has no magnitude — carried structural.

HMF (K5a): Kocadagli glucose system (Kocada2016.pdf = jafc.6b01862 — NOTE THE
FILE SWAP) = primary FIT; its NaCl arm + Gursul's 27 C zero-accumulation row =
sharpest hold-outs. Hamzalioglu 2018 (DOI 10.1016/j.foodchem.2017.07.131) =
the HMF SINK: first second-order HMF+Cys constants (3.95/5.15/23.3 M^-1 day^-1
at 5/25/50 C, pH 3.5) + HMF self-degradation + a SAME-METHOD matrix-vs-water
rate-ratio pair (Cys selectivity 11.4x water -> 1.2x coffee) — the class
previously declared nonexistent. Use the dossier's REFIT prefactors, not the
published ones (2/6 reproduce; third confirmed real-Ea-bad-prefactor case).
REFUSED ENTIRELY: Gursul Aktag Table 2 Ea (0/43 reproduce; six mathematically
underivable). Agcam 2022 = only internally reproducible Ea table (9/9) but
uniresponse synthetic juice — priors with scope flags. The fructose-limb vs
3-DG question is RESOLVED AS TWO QUANTITIES: fructose wins on flux (model
discrimination, unmeasured intermediate), 3-DG wins on published constants
(measured species) — the node design carries both parallel sources + one
sink, branch fraction matrix-dependent, never hard-coded.

PROCESS: Kocadagli2016.pdf (wheat flour) has a cipher-garbled text layer —
greps CANNOT see its content (round-3's negative sweep was blind to it); the
dossier publishes the solved glyph map. QUEUED for the quality pass: audit
every A_value in the corpus against refits from source k-tables (three
confirmed cases of correct Ea bolted to defective prefactors).

## Amendment 13 — 2026-08-29 (Wave B2.4: a declared weighting, two permanent scorer conditions, and the shared rows)

Written by the build wave itself, additive only. **No stoichiometry, no
species, no reaction, no benchmark, no hold-out row, no pass band, and no role
changes.** Pre-registered in `results/validation/kinetic_core_b2_4_prereg.md`
before any fit ran; mandated by `d1_exam_panel_reconciliation.md` §8 items 1–5.

1. **THE pH-UNIT-vs-LOG-FOLD EXCHANGE RATE IS NOW DECLARED.** Through B2.3 the
   objective scored three `zhou_final_pH_*` rows in **pH units** at σ = 0.25
   inside the same sum of squares as 55 rows in **log₁₀ folds** at σ_log =
   0.20–0.60. Nobody chose the resulting exchange rate; D1 §3 measured it at ~9×
   (one pH unit priced at nine times a 3× level miss) and showed one of those
   rows carrying **44% of the entire B2.2→B2.3 refit**. From B2.4 the rate is a
   declared scalar `PH_ENDPOINT_WEIGHT`, `E` = decades of level error per pH
   unit, implemented as `σ_ph = 0.35 / E`, and the fit is run and **reported at
   three declared values**: **1.40** (B2.3's accidental rate made explicit —
   reproduces its objective exactly, the control), **0.70** (a 4×-in-cost
   down-weight), **0.28** (the fit corpus's own measured rate: Kumazawa 2003's
   six declared-FIT retention rungs move 2-furfurylthiol 0.956 decades across
   3.4 pH units). The shipping value is chosen by a criterion written before any
   fit ran and **blind to both scorecards** (prereg §3).

   **CORRECTION, DISCLOSED.** `E = 0.28` replaces a first draft of **0.32** that
   was derived from D1 §4's Hofmann pH-3-to-pH-7 slope — **a hold-out-derived
   quantity**. It was caught by this wave's own firewall test as the literal
   `1.28` on the fit side, before any member at that weighting ran, and the
   basis was re-derived from declared FIT rows only. Recorded as prereg
   Correction C1.

2. **TWO PERMANENT SCORER CONDITIONS, on both scorecards.** The cutover exam
   now reports the **geometric-mean fold** beside the paired median (D1 §2: four
   rows crossing the middle of a lumpy 23-point pool swung the median 4.7×
   between B2.2 and B2.3 while the geometric mean moved 1.78×, and the reported
   median row did not move in either wave). The kinetic-core hold-out panel now
   reports **median |log₁₀ fold| and geometric-mean fold beside the gating
   count** (D1 §5: 22 of 34 gating rows were already failing, and 1.42 net
   decades of degradation landed entirely on them, invisible to a censored
   count). Both are additive; no existing number changes.

3. **THE FOUR SHARED ROWS ARE DECLARED IN BOTH ARTIFACTS.**
   `hofmann_ribose_pH{3,7}_{FFT,MFT}` on the panel and the four
   `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH{3,7}` MFT/FFT points on
   the exam are **the same measurements**, not analogues. The exam and the panel
   are therefore **not independent evidence** on that axis, and agreement
   between them there is one measurement counted twice. Declared in both
   artifacts from this wave onward.

4. **`ph_state.AMINE_FATE_BASIS`'s CITATION IS REPAIRED.** D1 §7 found the cited
   pre-freeze probe (`scratch/b23_encoding_probe.py`) cannot run — `KeyError:
   'AMN'`, a species not present in `src/kinetic_core/` — so the wave's
   self-declared strongest assumption cited evidence nobody could reproduce, and
   two of that probe's three "encodings" were the same code path. B2.4 rebuilds
   the axis against the **current** species set
   (`scripts/generators/probe_amine_fate_b2_4.py`), deriving the released amino
   nitrogen as `(Cys + ARP + TTCA at t=0) − (Cys + ARP + TTCA at t)` with no new
   species and no new parameter. **The rebuild reproduces the B2.3
   pre-registration's published three-row table in every cell to two decimals**,
   so the defect was in the script and not in the evidence: the declaration
   keeps its basis, its citation now points at a probe that runs, and its
   self-referential last sentence is removed. Frozen result:
   `results/validation/kinetic_core_b2_4_amine_fate_probe.json`.

5. **A CORPUS CONTRADICTION, NAMED (D1 §8 item 5): KUMAZAWA vs HOFMANN AT
   pH ≥ 6.5.** Kumazawa 2003 Fig. 3 (a declared FIT row) says 89% of spiked
   2-furfurylthiol is destroyed in **10 min at 121 °C at pH 6.4**. Hofmann &
   Schieberle 1998 (a declared HOLD-OUT) reports **12 ppb of 2-furfurylthiol
   still present after 20 min at 145 °C at pH 7**. No parameter choice satisfies
   both: a destruction law fitted to Kumazawa's top rung and extrapolated in pH,
   temperature and time annihilates Hofmann's pH-7 pot. **This is a different
   conflict from the Hofmann-vs-Zhou MFT sign conflict already recorded** (that
   one is about the DIRECTION of the MFT-vs-pH response between two systems;
   this one is about a thiol's survival at pH ≥ 6.5) and it is recorded here as
   a corpus contradiction in its own right. **Reported, not resolved.** Neither
   dataset is demoted and neither role changes.

6. **TWO PROCESS DEFECTS FOUND WHILE DOING THE ABOVE, both reported and both
   repaired only in the narrowest way:**
   * **The exam's pool moved under a concurrent wave.** D1's geometric means are
     over **23** answered points; the exam now answers **27**, because Wave B6's
     lipid lane answers four `matrix_path` points the core previously declined,
     at a family median near 1900×. Any B2.4-vs-D1 comparison of a pool-wide
     statistic is therefore **not like-for-like**, and every cross-wave number in
     the B2.4 reports is computed on **D1's own 23-point pool** with the
     27-point number printed beside it.
   * **The firewall's literal list cannot express a wave named B2.4.** The
     shipped list contains `"2.4"` (a Yiltirak hold-out fold), which occurs
     inside the string `B2.4`, so as first written the guard made it impossible
     for any fit-side file to name the fourth wave of its own module. The regex
     gained one lookbehind — a match glued to a letter is exempt — which costs
     the guard no reach, because a measured value is never immediately preceded
     by a letter. An explicit, reason-carrying `FIREWALL-OK` marker is the only
     other escape and a test asserts every such marker states its reason.

## Amendment 14 — 2026-08-29 (Wave B7, the furanic channel: two disclosures, two findings, one request for a ruling)

Written by the build wave itself, additive only. **No dataset changes columns,
no hold-out role changes, no benchmark record is edited, and no B1/B2/B3/B6
fitted constant is moved.** Pre-registered in
`results/validation/kinetic_core_b7_prereg.md` before any fit ran.

1. **EXPOSURE DISCLOSURE — the seven exam rows are `seen_diagnostic`,
   permanently.** `results/validation/cutover_final_exam.md` prints the MEASURED
   value of all 40 exam points in its per-point table, and B7's builder opened
   it while locating the five HMF and two DMHF bundles the wave was asked to
   convert. No file under `data/benchmarks/external_validation/` was opened.
   Under the Amendment 9 clause 1 / Amendment 10 clause 1 precedent the seven
   rows are demoted to **seen-diagnostic and may never gate**. The mitigation
   adopted is STRUCTURAL, as B6's was: **the HMF node has no fitted parameter at
   all** — its seven formation/sink constants are ingested whole from Kocadagli
   & Gokmen 2016's glucose Table 2 and Hamzalioglu & Gokmen 2018's Table 1 refit
   — and **the DMHF node has exactly one**, `k_dpo_af`, calibrated on three
   pentose cells of Blank 1997 at 90 °C, a system with no glucose, no alanine
   and no relation to any exam bundle. A literal-grep firewall
   (`tests/unit/test_kinetic_core_b7.py`) asserts that no hold-out-only literal
   appears in the furanic package or in the frozen fit report.
   **Also seen, and disclosed for the same reason:** every B7 hold-out value is
   printed in K5a §9.1 / K5b §8.3, which this wave was instructed to read. Each
   hold-out therefore carries a written, quantitative prediction made before its
   scorer existed.

2. **THE FURANIC CHANNEL HANGS ON THE TRUNK, AND THE TRUNK MOVED.** All four of
   its parents — fructose, 3-deoxyglucosone, 1-deoxyglucosone, methylglyoxal —
   are B1 trunk species, so the eleven new steps run inside the trunk, sulfur
   AND acrylamide lanes and there is no new lane conflict. The cost is that four
   of them put a new drain on a B1-fitted species. **B1 is NOT refit.** Measured
   consequences, all reported rather than absorbed: Martins' 100 °C regression
   moves ≤ 5.8 % on any species (fructose most, which is the right signature);
   18 previously-answered exam rows move, the largest by **1.26×**, inside the
   pre-registered 1.5× ceiling.

3. **★ REPORTED FOR ORCHESTRATOR RULING — B2.3's objective at its own frozen
   vector has moved 2.6 % (8.1754 → 8.3862) without one of its 48 constants
   being touched.** The sulfur network runs the trunk, so a trunk change moves
   B2.3's residuals. **This wave did NOT refit anything to absorb it**: refitting
   the sulfur lane to swallow a trunk change is precisely the undeclared
   exchange-rate move that `d1_exam_panel_reconciliation.md` found accounted for
   96 % of the B2.2→B2.3 exam regression. The B2.4 vector is therefore now
   fitted against a slightly different trunk from the one that ships. The
   discrepancy is pinned, with its size and its reason, in
   `tests/unit/test_kinetic_core_b2_4.py`. **A ruling on whether to re-run
   B2.3/B2.4 against the B7 trunk is requested; nothing presumes one.**

4. **TWO PRE-REGISTERED PREDICTIONS WERE WRONG, AND BOTH ARE DISCLOSED IN THE
   ARTEFACTS RATHER THAN QUIETLY RESTATED:**
   * **HMF direction.** The pre-registration argued that a node with no
     validated sink at cooking temperature (K5a G2: the 50–150 °C window is
     empty) must OVER-predict. **All five exam rows UNDER-predict**, by
     2.0–11.9×. The missing sink is therefore not the binding constraint — the
     SOURCE terms are, and the melt→matrix transfer loses more flux than the
     absent sink adds back. A directional result about a transfer nobody had
     tested.
   * **Gursul Aktag's 27 °C row.** Pre-registered PASS on the declared 100 µg/L
     floor; **it FAILS at 1 171 µg/L.** The cause is diagnosed and is not a
     defect in the ingested activation energies: at 27 °C the model's HMF is
     99.9 % 3-DG-limb, and that limb's terminal step carries **Ea = 0 by
     Kocadagli's own choice** (footnote: "does not follow Arrhenius equation").
     A zero-barrier terminal step cannot switch off as temperature falls.
     Amendment 12 called this row "the cheapest and most informative single test
     in the module"; it earned the description. **G1 measured for the first
     time.**

5. **A CONTRADICTION FOUND WHILE DOING THE ABOVE, REPORTED NOT IMPROVISED.**
   Amendment 12 records the HEMF alanine preference as "~10-25x". Computed from
   Blank 1997 Table 1 itself, the three sugars give **5.2× (arabinose), 25×
   (xylose) and 14.3× (ribose)**. The 10–25× band is the xylose/ribose pair;
   **arabinose sits at half its lower edge.** The correction Amendment 12 makes
   to Amendment 8 stands entirely — it is a preference and not a switch — but
   the range should read **≈5–25×**. Nothing was changed on this account.

6. **DECLARED, AND CARRIED ON EVERY ANSWER:** the module has NO pH term on any
   furanic edge (K5a G8 — six pH values across seven papers and no single paper
   varies pH, which `k3` §B.2 already forbids fitting across); NO activation
   energy on any furanone edge (all five K5b papers are single-temperature —
   the partition barrier is a declared assumption banded at ±50 kJ/mol and
   priced by re-integration); NO magnitude on the DMHF cysteine sink (Edge C
   ships at exactly zero, and fitting Shu & Ho's 6.0 % GC area is a named
   prohibited derivation); and Hamzalioglu's sink constant is **CLAMPED at
   50 °C**, never extrapolated. The two Schibilsky bundles differ only in pH and
   the model returns identical HMF and DMHF for both — a pre-registered
   structural miss whose size (1.76× on HMF, 5.1× on DMHF) is now the measured
   size of G8.

## Amendment 15 — 2026-08-30 (orchestrator rulings on Wave B7's flags)

1. The 2.6% movement of B2.3's frozen objective at its own frozen vector
   (the furanic channel adds trunk consumption edges; the sulfur lane runs
   the trunk) is ACCEPTED as physics propagation — no refit occurs in
   response (the D1 sigma-exchange failure class); the drift is pinned with
   its reason and will be absorbed by the next pre-registered sulfur wave.
2. Amendment 12's HEMF alanine-preference range is corrected to ~5-25x
   (computed per-sugar from Blank 1997 Table 1); the Amendment 8 correction
   stands.
3. B7's seven exam rows are seen_diagnostic per its disclosed exposure of
   cutover_final_exam.md; ratified with the structural mitigation noted
   (zero fitted parameters in the HMF node; the one fitted DMHF constant is
   identified on a pentose system appearing in no exam bundle).
