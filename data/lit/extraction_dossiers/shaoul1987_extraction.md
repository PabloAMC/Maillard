# Shaoul & Sporns 1987 — EXTRACTION (file stem `shaoul1987`)
Orly Shaoul & Peter Sporns, "Hydrolytic stability at intermediate pHs of the common purine nucleotides in food,
inosine-5'-monophosphate, guanosine-5'-monophosphate and adenosine-5'-monophosphate", *J. Food Sci.* **52**(3):810-812
(1987). DOI: 10.1111/j.1365-2621.1987.tb06733.x (no DOI printed on the scan). **Crossref 2026-11-03**: title, authors
(Orly Shaoul, Peter Sporns), J. Food Sci. 52(3) 810-812, issued 1987-05: **matches the PDF.** Crossref also lists a
conference abstract by the same authors, "Stability of 5'-purine nucleotides under canning conditions", Can. Inst.
Food Sci. Technol. J. 19(4):xxxvi (1986), 10.1016/s0315-5463(86)71510-1: presumably the same data, not on disk.
No naming caveat. Page 3 of the PDF carries the tail of an unrelated NMR/sucrose paper (references only).

### IMP, GMP, AMP (30 mg/mL) in 0.02 M phosphate, sealed Pyrex tubes, silicone-oil bath. (i) 121 °C, pH 3.0-8.0 (9 pH), sampled every 15 min for 3 h; (ii) pH 5, 121 / 106 / 91 / 76 / 53 °C, run until ≥ half hydrolysed. Pseudo-first-order k (HPLC loss of parent) with half-lives and r; Arrhenius and "shelf-life" regressions at pH 5. **Products named, not quantified; no ribose, no R5P, no Pi, no pH after heating, no replicates stated.**

**Read 2026-11-03 from `data/articles/shaoul1987.pdf`** (3 pp. scan with OCR layer). Prose from `pdftotext -layout`.
**Table 1 (p. 1) read by eye from 300 dpi crops** (scratchpad `ingest/shaoul1987_work/t1-1.png`, `t1b-1.png`);
**Tables 2-4 and Fig. 1 (p. 2) from 300 dpi crops** (`t234-2.png`, `f1-2.png`); methods paragraph (`meth-1.png`).
[!] The text layer misreads Table 1 IMP pH 4.0 as "0.07 11"; the print is **0.0111** (and the half-life 62.4 min
agrees only with 0.0111). Arrhenius re-fit of Table 2 by a short numpy script (scratchpad
`scripts/shaoul_fit.py`) [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Science, **University of Alberta**, Edmonton; senior author **Peter Sporns** (no
corresponding author marked; Sporns is the PI, Shaoul the student). **Not in `data/keys/laboratories.yml`: NEW**,
proposed id **`alberta_sporns`**, lab string "Univ. Alberta (Sporns)". Same group: Nguyen & Sporns 1984 (HPLC
method) and 1985 (J. Food Sci. 50:812, nucleotide/MSG decomposition during canning) [C], neither on disk; Yaylayan &
Sporns 1989 is cited in `apriyantono1993_extraction.md` (a citation, not a dossier). No same-lab dossier exists.

## 1. What this paper answers for the model

| model need | what Shaoul 1987 gives | tag |
|---|---|---|
| IMP loss rate vs pH (aqueous) | Table 1: k (min⁻¹), t½ (min), r at 121 °C, pH 3.0 / 4.0 / 4.5 / 5.0 / 5.5 / 6.0 / 6.5 / 7.0 / 8.0; also GMP, AMP | [M] |
| IMP loss rate vs T | Table 2: k (h⁻¹), t½ (h), r at pH 5, 121 / 106 / 91 / 76 / 53 °C | [M] |
| Ea / A | Table 3: slope (K) and ln-intercept (h⁻¹) at pH 5 → Ea 125.5 kJ/mol (IMP) | [M] → [Z] |
| which bond breaks | pH ≥ 4: "virtually the only product" is the nucleoside (inosine) [Q, V]; pH 3: nucleoside major plus "significant" base (hypoxanthine) [Q] | [Q] |
| product amounts (inosine, hypoxanthine, ribose, R5P, Pi) | **none printed**; checked "after about one half-life" only | [!] |
| pH during / after heating | **not measured**; buffer 0.02 M phosphate at 57-86 mM nucleotide | [!] |
| pH × T interaction (Ea at other pH) | **none**: T series at pH 5 only | [!] |
| replicates / ± on k | **none**; r only | [!] |

**Answers.**
1. **Rate constants vs pH and T: yes, a clean, internally consistent set** (confidence ~95 % that the printed numbers
   are what was measured). Every k ↔ t½ pair satisfies t½ = ln 2 / k to the printed digits [Z, all 27 + 15 pairs
   checked]. IMP at 121 °C: k = 0.0214 min⁻¹ (pH 3) → a plateau 0.0097-0.0111 (pH 4-6) → 0.0067 / 0.0055 / 0.0019
   (pH 6.5 / 7 / 8). At pH 5: k(IMP) = 0.662 h⁻¹ (121 °C) to 2.32 × 10⁻⁴ h⁻¹ (53 °C). My Arrhenius fit of Table 2
   reproduces Table 3: **Ea(IMP) = 125.5 ± 1.6 kJ/mol, ln A = 37.95 (A in h⁻¹)**; GMP 122.9 ± 1.5, AMP 126.4 ± 4.2
   (± = regression s.e., 5 points, 3 dof) [Z]. Curvature is small (IMP max |residual| 0.11 in ln k).
2. **Products and buffer** (confidence ~85 % on the mechanism reading). Products identified by HPLC against bought
   standards (nucleotides, nucleosides, bases), but **reported only qualitatively**: the reaction at pH 4-8 is
   phosphate-monoester hydrolysis IMP → inosine + Pi; at pH 3 glycosidic cleavage also runs (base seen), which the
   authors use to explain the pH-3 rise. Ribose, ribose-5-phosphate and Pi were not measured. Buffer: 0.02 M
   phosphate at every pH, including pH 3-5.5 where phosphate barely buffers (pKa 2.1 / 7.2) [Z]; the nucleotide
   itself (57-86 mM, pKa ~6) out-buffers the buffer [Z].
3. **Usable as an IMP → ribose(-5-P) source term? Not for R5P; yes for the first step of an IMP → inosine → ribose
   chain** (confidence ~80 %). The paper's own product call says the C-O-P bond breaks first and the N-glycosidic
   bond survives at pH ≥ 4, so IMP releases **inosine + Pi, not R5P**; inosine is not a reducing sugar. A ribose
   source needs a second step, inosine → hypoxanthine + ribose, which this paper does not measure
   (`davidek1972_extraction.md` gives it at 90 °C vs pH, without Ea). P(direct IMP → hypoxanthine + R5P is a ≥ 10 %
   branch at pH 4-7) ≈ 25 %: the sibling `matoba1988_extraction.md` estimates ~18 % (15-22 %) of IMP going straight
   to base + R5P at pH 4, 100 °C from its Figs. 3-4 [Z there], which Shaoul's "virtually only the nucleoside" call
   at pH ≥ 4 does not rule out at that level; at pH 3 unknown (base seen, origin not resolved). If R5P forms it
   dephosphorylates to ribose fast (Matoba Table 3: 8.3 × 10⁻² h⁻¹, pH 4, 100 °C).
   A defensible first-step term: k1(T, pH) = k121(pH) · exp[−(Ea/R)(1/T − 1/394.15 K)], k121 from Table 1, Ea =
   125.5 kJ/mol **measured at pH 5 only** (extending it to other pH is an assumption; Davidek 1972's IMP Ea runs
   142 → 118 kJ/mol from pH 3 → 5, so a ±15 kJ/mol band is honest) [Z, P]. Illustration [Z]: at pH 5.5, 140 °C,
   k1 ≈ 0.061 min⁻¹ (t½ ≈ 11 min); at 100 °C, t½ ≈ 8.4 h.
4. **Cross-check against Davidek 1972 (90 °C, independent lab, ion exchange + UV)** [Z]. Shaoul extrapolated to 90 °C
   (pH-5 Ea at all pH) vs Davidek Fig. 2 read, log10(k/s⁻¹) + 7: pH 4: 1.88 vs 1.83; pH 5: 1.88 vs 1.54; pH 6: 1.82 vs
   1.31; pH 7: 1.58 vs 1.15; pH 8: 1.11 vs 1.02. **Agreement within 0.1 log at pH 4 and 8; Davidek 2-3× slower at
   pH 5-7.** Shaoul's pH 4-6 plateau is absent in Davidek (monotone fall). Candidate causes: phosphate general
   acid/base catalysis (Davidek reports buffer catalysis) or pH drift in Shaoul's weakly buffered tubes. Between-
   paper spread ≤ 0.5 log units: the right width for a prior on k1.
4b. **Against Matoba 1988 (sibling dossier, 0.01 M acetate / MOPS, 0.4 mM, 80-100 °C)** [Z]: Shaoul extrapolated to
   100 °C gives IMP t½ 8.3 h (pH 4) and 16.7 h (pH 7) vs Matoba's printed 8.7 and 13.1 h; Matoba's IMP Ea
   1.1 / 1.3 × 10⁵ J/mol (pH 4 / 7) brackets Shaoul's 125.5 kJ/mol. **Three labs agree within ×1.3 at pH 4** (Davidek
   9.2 h at 100 °C, pH 4 [Z]); at pH 7 Davidek is the outlier (~48 h [Z], ×3-4 slower).
5. **What it does NOT give**: no product yields, no ribose/R5P, no Ea at any pH but 5, no pH trajectory, no matrix
   (buffer only), no replicates or ± on k, no temperatures above 121 °C (cooking at 140-180 °C is an extrapolation of
   ~20-60 K beyond the data).

## 2. Conditions

- **Reactants** [M]: IMP disodium salt and hypoxanthine (Sigma); GMP disodium salt, AMP, adenosine, adenine, guanine
  (Aldrich); inosine, guanosine (Terochem). Purity not stated ("analytical grade or better" for other reagents).
  Each nucleotide alone, **30 mg/mL** [M, read by eye] → ~57-77 mM for Na2IMP (hydrate state not stated) [Z].
- **Buffer / pH** [M]: 0.02 M phosphate "at the appropriate pH"; pH 3.0-8.0. Temperature of the pH reading not stated,
  presumably room temperature; **no pH after heating** [!].
- **Vessel / headspace** [M]: 0.5 mL in 7 cm × 6 mm o.d. Pyrex tubes, flame-sealed (air headspace, small).
- **Heating** [M]: silicone-oil bath at 121 °C, stirred, five thermocouples logged (Fluke 2240B); heat-up of a 6 mm
  tube is seconds to ~1 min, negligible against t½ ≥ 31 min [Z]. Quench: ice bath, then freezer.
- **Times** [M]: 121 °C: every 15 min for 3 h (13 points presumably); lower T: "as convenient", until ≥ half
  hydrolysed (53 °C: t½ 2980 h ≈ 124 d).
- **Replicates**: not stated [!].
- Water activity: dilute aqueous (aw ≈ 1).

## 3. Analytical methods

HPLC per Nguyen & Sporns 1984 (J. Assoc. Off. Anal. Chem. 67(4), page garbled in the scan [C]; column and detector not given here). Rate from the decrease of the
parent nucleotide; first-order fit, linear correlation coefficient r printed (no ±). Quantification basis: external
standards presumably (standards purchased); LOD not stated. The pH 8 runs reach only ~30-40 % conversion in the 3 h
window (t½ 239-365 min) [Z], so those k rest on shallow decays (AMP pH 8 r = 0.894).

## 4. Data

### 4.1 Table 1 — first-order k at 121 °C vs pH [M, read by eye]; k in min⁻¹ (t½ in min), r

| pH | IMP | r | GMP | r | AMP | r |
|---|---|---|---|---|---|---|
| 3.0 | 0.0214 (32.4) | 0.988 | 0.0220 (31.5) | 0.990 | 0.0153 (45.3) | 0.987 |
| 4.0 | 0.0111 (62.4) | 0.990 | 0.0151 (45.9) | 0.982 | 0.0120 (57.8) | 0.991 |
| 4.5 | 0.0109 (63.6) | 0.992 | 0.0167 (41.5) | 0.995 | 0.0121 (57.3) | 0.981 |
| 5.0 | 0.0110 (63.0) | 0.994 | 0.0169 (41.0) | 0.996 | 0.0136 (51.0) | 0.977 |
| 5.5 | 0.0104 (66.6) | 0.995 | 0.0158 (43.9) | 0.992 | 0.0106 (65.4) | 0.986 |
| 6.0 | 0.0097 (71.4) | 0.993 | 0.0153 (45.3) | 0.990 | 0.0100 (69.3) | 0.987 |
| 6.5 | 0.0067 (103) | 0.996 | 0.0111 (62.4) | 0.992 | 0.0079 (87.7) | 0.985 |
| 7.0 | 0.0055 (126) | 0.938 | 0.0067 (103) | 0.992 | 0.0049 (141) | 0.983 |
| 8.0 | 0.0019 (365) | 0.974 | 0.0029 (239) | 0.981 | 0.0027 (257) | 0.894 |

Fig. 1 (p. 2) is a smoothed curve of these same points (k vs pH, pH ~2.7-8.3); no extra data.

### 4.2 Table 2 — first-order k at pH 5 vs T [M]; k in h⁻¹ (t½ in h), r

| T (°C) | IMP | r | GMP | r | AMP | r |
|---|---|---|---|---|---|---|
| 121 | 0.662 (1.05) | 0.994 | 1.015 (0.68) | 0.996 | 0.816 (0.85) | 0.977 |
| 106 | 0.158 (4.39) | 0.988 | 0.228 (3.04) | 0.998 | 0.116 (5.97) | 0.975 |
| 91 | 0.0339 (20.4) | 0.995 | 0.0516 (13.4) | 0.996 | 0.0348 (19.9) | 0.981 |
| 76 | 0.00516 (134) | 0.959 | 0.00895 (77.4) | 0.969 | 0.00469 (148) | 0.968 |
| 53 | 0.000232 (2,980) | 0.924 | 0.000398 (1,740) | 0.963 | 0.000231 (3,000) | 0.958 |

### 4.3 Table 3 — Arrhenius (ln k vs 1/T, k in h⁻¹) at pH 5 [M]; Table 4 — "shelf-life plot" [M]

| nucleotide | slope (K) | intercept ln(h⁻¹) | r | Ea [Z] kJ/mol | shelf-life slope (°C⁻¹) | intercept (log10 h) | r |
|---|---|---|---|---|---|---|---|
| IMP | −1.51 × 10⁴ | 38.0 | 0.999 | 125.5 | −0.0509 | 6.07 | 0.997 |
| GMP | −1.48 × 10⁴ | 37.6 | 0.999 | 123.0 | −0.0499 | 5.77 | 0.997 |
| AMP | −1.52 × 10⁴ | 38.3 | 0.998 | 126.4 | −0.0514 | 6.12 | 0.997 |

Text extrapolations to 23 °C, pH 5 [M]: Arrhenius t½ 35.7 / 18.7 / 39.7 yr (abstract: 36 / 19 / 40); shelf-life plot
3.99 / 4.83 / 10.0 yr (IMP / GMP / AMP). Fig. 2 (AMP, Arrhenius vs shelf-life extrapolations at 43 / 33 / 23 °C) adds
no data.

## 5. Flags and internal consistency checks

- **k ↔ t½: all consistent** to the printed digits (Tables 1, 2) [Z]. Table 2 121 °C row = Table 1 pH 5 row × 60
  (0.0110 → 0.660 vs 0.662; 0.0169 → 1.014 vs 1.015; 0.0136 → 0.816) [Z]: probably the same run reused although the
  text calls the T series a "second series" (~70 %) [!].
- **Table 3 reproduces from Table 2** (my fit: −15091 K, 37.95 for IMP) [Z]. Arrhenius 23 °C half-lives reproduce
  within rounding (my 35.3 / 18.4 / 39.0 yr) [Z].
- [!] **IMP shelf-life value 3.99 yr does not follow from Table 4**: 10^(6.07 − 0.0509 × 23) h = 9.05 yr; GMP (4.78
  vs 4.83) and AMP (9.89 vs 10.0) do reproduce [Z]. A text or table typo for IMP; irrelevant to the model.
- [!] **"Maximum rate at about pH 5"** [Q] holds for GMP and AMP; IMP is flat at pH 4-5.5 (0.0104-0.0111, within
  their scatter) [Z].
- [!] **Weak buffering**: 0.02 M phosphate under ~60-80 mM nucleotide; releasing Pi (pKa2 7.2) from a monoester
  (pKa2 ~6.2) takes up protons, so pH likely rose in the pH 6-8 tubes during reaction [Z, P]; pH not re-measured.
- pH 8 and the 53-76 °C runs have the lowest r (0.89-0.97); shallow decays (pH 8: ~30-40 % conversion in 3 h) [Z].
- Air in sealed tubes; purines are not O2-sensitive in this regime (~90 % no effect) [P].

## 6. Existing records in the repo that cite this paper

`grep -rIn -i -E "shaoul|tb06733"` over `data/ docs/ src/ scripts/ tasks/` (excluding the alias line):

| record | claim | check vs paper |
|---|---|---|
| `data/lit/extraction_dossiers/madruga1998_extraction.md:97-99` (read from the text layer) | IMP phosphate hydrolysis to inosine at 121 °C, t½ ~32 / ~64 / ~71 min at pH 3.0 / 4.5 / 6.0 | **match** by eye (32.4 / 63.6 / 71.4). Wording nit: the pH 3.0 value is total IMP loss (phosphate + glycosidic), not phosphate hydrolysis alone. |
| same, l. 99 | glycosidic cleavage significant only at pH 3 | **match** [Q]. |
| same, l. 99-100 | "ribose supply rises steeply toward pH 3 on its own" | **inference, consistent** with the paper's product call; Shaoul measured no ribose. |
| `data/lit/extraction_dossiers/zhang1991_extraction.md:157` | Zhang 1991 cites Matoba 1988 and Shaoul & Sporns 1987 for hydrolysis kinetics | no Shaoul number to check; not checked against the zhang1991 PDF. |
| `docs/guides/LITERATURE_SEARCH_2026-10_ROUND3.md:37` | "121 C, pH 3-8, Arrhenius" | **match**, with a caveat: the Arrhenius series is **pH 5 only**, 53-121 °C. |
| same, l. 355 | t½ 63 min at pH 5 (GMP 41, AMP 51); glycosidic cleavage too at pH 3; Arrhenius over several T | **match** (63.0 / 41.0 / 51.0; five temperatures). |
| same, l. 355 | "citers: rate at pH 7 double that at pH 8" | **partial**: k(pH 7)/k(pH 8) = 2.9 (IMP), 2.3 (GMP), 1.8 (AMP). |
| same, l. 416 | half-lives at 121 C, pH 3-8, plus Arrhenius data | **match**. |
| `data/keys/pdf_dossier_aliases.yml:93` | IMP / GMP / AMP phosphate hydrolysis at 121 C, pH 3-8, half-lives | **match**. |

**Related prior, not citing this paper** (`data/lit/computational_priors.json:956-972`,
`matoba_1988_nucleotide_hydrolysis_v1`, IMP t½ at 100 °C by pH and at pH 7 by T; Matoba 1988 is the sibling
dossier's job). Shaoul-derived values for comparison [Z]: 100 °C t½ = 8.3 h (pH 4), 8.4 h (pH 5), 16.7 h (pH 7,
assuming the pH-5 Ea) vs the prior's 8.7 / 12.0 / 13.1 h; pH 7, 121 °C: 2.1 h vs the prior's log-interpolated 1.3 h;
the prior's 100 → 150 °C span implies Ea ≈ 142 kJ/mol at pH 7 vs Shaoul's 125.5 at pH 5. Agreement within ×1.5.

## 7. Proposed ledger entries (trunk.yml; not written)

```yaml
  - id: shaoul1987_tab1_nucleotide_k_vs_ph_121c
    dossier: shaoul1987_extraction.md
    paper: "Shaoul & Sporns 1987 J Food Sci 52:810"
    lab: "Univ. Alberta (Sporns)"
    system: "IMP, GMP or AMP alone, 30 mg/mL (~57-86 mM)"
    medium: "0.02 M phosphate; pH set at room temperature presumably, not re-measured; sealed Pyrex tubes, oil bath"
    temperature_c: [121]
    ph: {initial: [3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0], measured_during: false, final: null}
    time_min: [0, 180]
    species: [IMP, GMP, AMP]
    kind: rate_constant
    location: "Table 1 (k min-1, t1/2 min, r; 27 cells); Fig. 1 smooths the same points"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "buffered at the appropriate pH with 0.02M phosphate buffer"
    flags: [ph_start_only, ph_not_measured, no_replicates, no_error_bars]
    note: "Parent loss only (HPLC). Product at pH >= 4 is the nucleoside (inosine + Pi), not R5P; base also at pH 3. No IMP species in the engine yet; first step of a future IMP -> inosine -> ribose chain. Text layer misreads IMP pH 4.0 k; dossier value read by eye (0.0111)."
  - id: shaoul1987_tab2_3_nucleotide_k_vs_t_ph5
    dossier: shaoul1987_extraction.md
    paper: "Shaoul & Sporns 1987 J Food Sci 52:810"
    lab: "Univ. Alberta (Sporns)"
    system: "IMP, GMP or AMP alone, 30 mg/mL (~57-86 mM)"
    medium: "0.02 M phosphate, pH 5; sealed Pyrex tubes"
    temperature_c: [53, 76, 91, 106, 121]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: []
    species: [IMP, GMP, AMP]
    kind: rate_constant
    location: "Table 2 (k h-1, t1/2 h, r); Table 3 Arrhenius slope/intercept (Ea 125.5 kJ/mol IMP); Table 4 shelf-life regression (same data)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "The pH for these hydrolyses was set at 5"
    flags: [ph_start_only, condition_not_stated:time_min, no_replicates, internal_inconsistency]
    note: "Run until >= half hydrolysed (53 C ~4 months). 121 C row equals Table 1 pH 5 (likely same run). Ea at pH 5 only. internal_inconsistency = IMP shelf-life 3.99 yr vs 9.05 yr from Table 4 (no bearing on k)."
```

## 8. Supplementary material

None referenced (1987 research paper, 3 pp.). Nothing model-relevant lies elsewhere, except the method (Nguyen &
Sporns 1984) and the group's earlier canning study (Nguyen & Sporns 1985, J. Food Sci. 50:812), neither on disk.
