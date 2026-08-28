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
