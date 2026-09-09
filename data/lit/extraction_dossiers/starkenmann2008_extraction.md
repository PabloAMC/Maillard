# Starkenmann et al. 2008 — EXTRACTION (cysteine-S-conjugates released to thiols by mouth microflora; NO Maillard kinetics)

**Source on disk:** `data/articles/starkenmann2008.pdf` (owner's download, 2026-08-28). Read-only extraction,
2026-09-07, from `pdftotext -layout`; the paper has no numeric tables (Figs 2–4 only) and every number below is
from the text or a caption. The repo already carries its two usable items: the 22 ng/L retronasal threshold
(`k2_matrix_and_thresholds.md`, `k3` — USE-Q, "already an in-mouth number") and the saliva thiol quench
(STRANDED; `docs/reference/FIT_HOLDOUT_DECLARATION.md`: "neither"). This dossier confirms both from the full text
and records that the paper contains nothing thermal.

## 0. Identity

| field | value |
|---|---|
| Title | "Olfactory Perception of Cysteine-S-Conjugates from Fruits and Vegetables" |
| Authors | Christian Starkenmann, Bénédicte Le Calvé, Yvan Niclass, Isabelle Cayeux, Sabine Beccucci, Myriam Troccaz (Firmenich SA, Geneva) |
| Venue | J. Agric. Food Chem. 2008, 56 (20), 9575–9580 |
| DOI | 10.1021/jf801873h (received 18 Jun 2008, accepted 24 Aug 2008, web 24 Sep 2008) |

## 1. Why it matters — and why it holds no kinetics for this model

The paper is about the mouth as a reactor: odourless S-alkyl-L-cysteines in wine, onion and bell pepper are cleaved
to volatile thiols by anaerobic oral bacteria (C–S β-lyase activity) over 20–30 s to minutes, and saliva proteins
bind free thiols. **No sugar, no heating, no Maillard step, no rate constant.** The cysteine-S-conjugate is the
free amino acid S-alkylated on the thiol — a stable, non-thermal precursor class — and its cleavage is enzymatic
(bacterial), not chemical. Nothing here bears on the sulfur lane's cysteine chemistry at 100–145 °C. What survives
for a model is perception-side: a retronasal threshold and evidence that saliva binds thiols strongly and
compound-class-specifically.

## 2. Methods as they matter to a model

| item | value | where |
|---|---|---|
| Compounds | S-(R/S)-3-(1-hexanol)-L-cysteine (1) → 3-sulfanylhexan-1-ol (4); S-propyl-L-cysteine (2) → 1-propanethiol (5); S-(2-heptyl)-L-cysteine (3) → 2-heptanethiol (6) | Fig. 1 |
| Thiol assay | Acrylodan derivatisation, HPLC-fluorescence (ex 390 / em 500 nm); LOD in water **0.001 mg/L**; conjugates by UPLC-ESI(+) SIM, calibration 0.02–10 mg/L (R² 0.9998) | p. 9575–9576 |
| Saliva | pooled from 4 adults; crude (centrifuged) or sterile (pasteurised 60 °C 1 h); **pH 7.6 ± 0.1; protein 0.59 ± 0.2 mg/mL; 8 × 10⁷ cfu/mL anaerobes in crude saliva**, none in sterile | p. 9576 |
| Binding experiment | thiol 4 at 1 mg/L (or 0.1 mg/L) in water + Na₂HPO₄ 0.2 M, saliva 1–10 % v/v, **45 min at 22 °C**, then Acrylodan | p. 9576 |
| Bacterial conversion | conjugate 5 mg/L (2.5 mg/L at t = 0 after dilution = 100 %) in crude saliva, sterile saliva or 0.9 % NaCl, ± Fusobacterium nucleatum DSM 20482 (5 × 10⁶ cfu), 37 °C, 2 h – 4 days | p. 9577 |
| Sensory | 30 trained Firmenich panellists; 30 mL sample held in mouth 5 s then spat (retronasal); 3-AFC threshold at 4, 8, 16, 32, 64, 120, 250, 500 ng/L, two sessions, geometric mean; dose–response at 0, 2.2, 4.4, 8.8, 17, 35, 70 ng/L; time–intensity (Tbegin, Tmax, Imax, Tend) | p. 9577, Fig. 2 |
| Replicates | sensory: two sessions / duplicate tastings (t-test p < 0.05); chemistry: ± values printed, n not stated | p. 9577–9579 |

## 3. Every number the text holds

| quantity | value | where |
|---|---|---|
| Retronasal detection threshold, 3-sulfanylhexan-1-ol, water | **22 ng/L** (lit. 1–60 ng/L) | p. 9577 |
| Time–intensity, 3-sulfanylhexan-1-ol 0.01 mg/L vs conjugate 1 at 1 mg/L | onset immediate vs 12 ± 2 s; Tmax 16 ± 5 s vs 39 ± 8 s; Tend 100 ± 12 s vs 163 ± 9 s | p. 9577 |
| Time–intensity, 1-propanethiol 0.01 mg/L vs S-propyl-cysteine 20 mg/L | Tmax 20 ± 8 s vs 37 ± 5 s; Tend 76 ± 10 s vs 156 ± 16 s | p. 9577 |
| Time–intensity, 2-heptanethiol 0.01 mg/L vs conjugate 3 at 1 mg/L | Tend 94 ± 17 s vs > 180 s | p. 9577 |
| Analytical LOD of thiol 4 in crude saliva | **60 ± 10 mg/L** (vs 0.001 mg/L in water — a 6 × 10⁴ quench; the paper's "3 × 10⁶ above its odor threshold" divides an analytical LOD by a sensory threshold) | p. 9578 |
| Saliva quench | 10 % crude saliva quenches 1 mg/L thiol; 2 % quenches 0.1 mg/L; sterile saliva quenches less | p. 9578, Fig. 3 |
| Conjugate 1 conversion, crude saliva, 37 °C | 20 % consumed at 2 h, 80 % at 24 h (same with added F. nucleatum) | p. 9579 |
| Conjugate 1, sterile saliva | < 15 ± 10 % consumed after 4 days | p. 9579 |
| Conjugate 1, sterile saliva + F. nucleatum | 20 ± 10 % at 2 h; + anaerobes 60 ± 10 % at 4 days; aerobes inactive | p. 9579, Fig. 4A |
| Conjugate 1 (5 mg/L) + F. nucleatum in 0.9 % NaCl | thiol 0.087 ± 0.005 mg/L at 2 h = **1.7 % yield**; 3.9 % at 1 day | p. 9579, Fig. 4B |
| Same in minimal medium, 3 days | 75 % consumed, **40 % thiol yield** | p. 9579 |
| S-propyl-L-cysteine in commercial onion powder | 1.9 ± 0.2 mg/kg (p. 9576; "1.8 mg/kg (± 0.2)" on p. 9577 — ⚠ internal inconsistency) | |
| S-(2-heptyl)-L-cysteine in green bell pepper | 0.057 ± 0.002 mg/kg (one lot); 0.001–0.09 mg/kg across five lots/seasons | p. 9576 |
| Real-time APCI-MS breath monitoring | no thiol signal even at 10 mg/L conjugate (instrument LOD too high) | p. 9577 |

## 4. What the repo could take

No FIT rows; no thermal or Maillard number exists in the paper. Items a model could use, all already recorded:

1. **22 ng/L retronasal threshold for 3-sulfanylhexan-1-ol** (USE-Q in k2/k3). It is measured with the sample in the
   mouth, so it already contains the salivary interaction — applying a saliva correction on top double-counts.
2. **Saliva binds free thiols at ≈ 6 × 10⁴ (analytical) — but no basis, stoichiometry or mechanism** (authors:
   "Whether this absorption results from physicochemical interactions, chemical transformations, or covalent linkage
   to glycoproteins remains unclear"). STRANDED, as `k2` §B.6 says; Baek 1999 shows the effect is compound-class-
   specific, so no single saliva factor.
3. Pathway fact for the sulfur lane's vocabulary: **S-alkyl-cysteines are stable, odourless conjugates whose thiol is
   released only by (bacterial) C–S lyase**; their thermal stability was not tested here, so they are not a thermal
   thiol reservoir the model can charge.
4. Perception-side timing: conjugate-derived thiols appear at 12–20 s and persist to 2.5–3 min — a delayed-release
   channel outside the model's scope.

## 5. Caveats

- All conversion figures are at 37 °C in saliva or saline with live bacteria; they are microbiological, not chemical.
- The saliva-binding series (Fig. 3) is undigitised; the two quench points in the text are the only numbers.
- Pooled saliva from four donors; bacterial counts were standardised but donor variability is not reported.
- The 3-AFC threshold protocol is retronasal (mouth), not orthonasal; do not compare it directly with orthonasal
  water thresholds in the corpus without noting the modality.
- Two conflicting onion-powder values (1.8 vs 1.9 mg/kg) — trivial, but a sign the text was not proof-read for
  numbers; read any single figure with that in mind.
