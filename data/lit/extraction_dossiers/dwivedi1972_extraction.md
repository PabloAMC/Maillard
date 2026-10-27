# Dwivedi & Arnold 1972 — EXTRACTION
### *"Chemistry of Thiamine Degradation. Mechanisms of Thiamine Degradation in a Model System"*
### The ³⁵S sulfur balance behind the review's pH split. At 121 °C, pH 3.5–6.0, 96–100 % of the degraded thiamine's sulfur is the intact thiazole (so H2S plus every other sulfur product is ≤ 4 %); at pH 7.0 only 9–16 % is, and at pH 8.0 none. One temperature. No H2S number of its own. Read 2026-10-27.

**Full read of `data/articles/dwivedi1972.pdf`** (3 pages, Acrobat Capture OCR). Table I and the
figure captions were checked on a 200-dpi raster, and Fig. 5 on a 400-dpi crop, because the text
layer misreads two Table I cells: pH 7.0 / 30 min "% destroyed" is **46.5** (text layer "16.5") and
pH 5.0 / 30 min thiazole is **99.0** (text layer "e9.0"). This is the ³⁵S study that
`dwivedi1973_extraction.md` lists as "not in the corpus" (cited there as J. Food Sci., in press).
Companions: `arnold1971_extraction.md` (the H2S data B50 uses) and `dwivedi1972b_extraction.md`
(the identification of the thiazole).

**Provenance:** **[M]** measured, printed as a number · **[Q]** qualitative · **[C]** cited from
another paper · **[Z]** derived by me from the paper's own numbers, never printed · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Basant K. Dwivedi, Roy G. Arnold |
| Title | "Chemistry of Thiamine Degradation" / subtitle "Mechanisms of Thiamine Degradation in a Model System" |
| Venue | *J. Food Sci.* **1972**, 37 (6), **886–888** (volume and pages from the running heads; issue from CrossRef) |
| DOI | **10.1111/j.1365-2621.1972.tb03694.x**. **Not printed in the PDF.** Found via the CrossRef REST API (`api.crossref.org/works?query.bibliographic=...`) on 2026-10-27; CrossRef metadata for this DOI: title "CHEMISTRY OF THIAMINE DEGRADATION", J. Food Sci. 37(6):886–888, Nov 1972, authors Dwivedi, Arnold. Matches. |
| Affiliation | Dept. of Food Science & Technology, University of Nebraska-Lincoln, Lincoln NE 68503 |
| Dates | Ms received 5/31/72, revised 7/15/72, accepted 7/16/72; presented at the 32nd IFT Annual Meeting, Minneapolis; Nebraska AES Journal Series Paper No. 3394, Project 16-19; Nutrition Foundation grant 427 |
| PDF character | 3 pages: Figs. 1–3 (radiochromatogram scans), Table I (12 rows), Fig. 4 (acid–base equilibria scheme), Fig. 5 (volatile structures). Duplicates stated for 15 and 30 min; no SD, no scatter |

**[!] Citation-key trap.** This paper's own reference list calls the J. Food Sci. 37:689 thiazole
paper "Dwivedi et al. 1972a" (in the corpus as `dwivedi1972b.pdf`) and a minor-volatiles paper
"Dwivedi et al. 1972b, submitted". The latter appeared as Dwivedi, Arnold & Libbey 1973, *J. Food
Sci.* 38:450–452, "Some minor volatile components from thermally degraded thiamine",
doi 10.1111/j.1365-2621.1973.tb01451.x (CrossRef). **That paper is not in the corpus.** The repo
file names do not follow the paper's a/b labels.

---

## 1. THE ONE-PARAGRAPH ANSWER

Thiamine-³⁵S hydrochloride (0.5 mg/mL, [Z] 1.48 mM) in 0.1 M phosphate at **pH 3.5, 5.0, 6.0, 7.0,
8.0** (25 °C labels), sealed vials, **autoclaved at 121 °C** for 15 and 30 min (pH 6.0 and 8.0 also
2 h). Products separated by TLC and quantified by radio-scan peak areas. Thiamine loss at 121 °C:
**4–8 % at pH 3.5, 9–18 % at pH 5–6 (41 % after 2 h at pH 6), 28–47 % at pH 7, 33–59 % at pH 8
(93 % after 2 h)**. **At pH 3.5–6.0 the intact 4-methyl-5-(β-hydroxyethyl)thiazole carries 96.0–100 %
of the destroyed thiamine's ³⁵S**, i.e. methylene-bridge cleavage is essentially the only route
and **all other sulfur products together, H2S included, are ≤ 4 % (≤ 1.5 % at 30 min)**. **At pH 7.0
the thiazole is only 16.4 % (15 min) and 8.9 % (30 min)** of the destroyed thiamine; the rest is
missing radioactivity the authors attribute to H2S and other volatile ³⁵S compounds lost in TLC,
plus two small peaks (C, D) and a broad peak near thiamine. **At pH 8.0 no thiazole peak at all.**
H2S itself was **not measured** in this paper. The authors interpret the pH ≥ 7 switch as
degradation through the minor pseudobase/thiol forms, with an Ea for H2S removal **lower** than
for bridge cleavage, and an Ea for the 5-carbon-fragment volatiles "much higher" than for H2S
elimination. **No numbers for either Ea.** 2-Methyl-3-furanthiol is not mentioned.

---

## 2. METHODS [M] (p. 886)

| item | as printed |
|---|---|
| Substrate | thiamine-³⁵S hydrochloride (Amersham/Searle), **0.5 mg/mL**, specific activity **0.05 mc/mL** [Z: 0.5/337.27 = 1.48 mM, ~12× the 1.2 × 10⁻⁴ M of Dwivedi & Arnold 1971] |
| Buffer / pH | "buffered to pH **3.5, 5.0, 6.0, 7.0 and 8.0** with **0.1 M phosphate** buffer". **Preparation-temperature (room-temperature) labels; no measurement after heating, no hot-pH correction** [!] |
| Vessel | **0.5 mL** per sealed screw-cap vial, Teflon-lined caps (closed system) |
| Heating | autoclave, **121 °C**, **15 and 30 min** (duplicates); one sample each of pH 6.0 and 8.0 for **2 hr**. Time counted from the solution **reaching 121 °C**; come-up time not counted [!] |
| Cooling | refrigerated **2 hr** immediately after autoclaving |
| Controls | unheated thiamine-³⁵S |
| Separation | TLC (Waring et al. 1968), silica gel 6060 sheets with fluorescent indicator; acetonitrile/H2O/formic acid 40/10/to pH **2.54**; run in a fume hood because H2S and volatile S compounds are produced. Spots assigned by Rf against non-radioactive 10 mg/mL thiamine run in parallel |
| Quantification | Vanguard 880 autoscanner; "% of total radioactivity ... contributed by individual spots was estimated by measuring the areas under each peak"; thiamine loss = heated vs unheated thiamine peak area. Thiamine peak (A) re-scanned at higher attenuation to keep it on scale |
| H2S | **not measured**; volatile ³⁵S escapes during TLC |

---

## 3. RESULTS

### 3.1 Table I, as printed [M] (p. 887), checked on the raster

"Production of 4-methyl-5-(β-hydroxyethyl) thiazole-S³⁵ in relation to thiamine-S³⁵ destruction at
various pH's and autoclaving times" (121 °C).

| pH | time of autoclaving (min) | % thiamine-S³⁵ destroyed | thiazole-S³⁵ derivative as % of thiamine-S³⁵ destroyed |
|---|---|---|---|
| 3.5 | 15 | 4.1 | 100.0 |
| 3.5 | 30 | 7.8 | 98.5 |
| 5.0 | 15 | 8.5 | 99.6 |
| 5.0 | 30 | 15.1 | 99.0 |
| 6.0 | 15 | 9.6 | 96.0 |
| 6.0 | 30 | 18.2 | 98.6 |
| 6.0 | 120 | 41.0 | 97.9 |
| 7.0 | 15 | 28.0 | 16.4 |
| 7.0 | 30 | 46.5 | 8.9 |
| 8.0 | 15 | 32.7 | – |
| 8.0 | 30 | 58.5 | – |
| 8.0 | 120 | 92.6 | – |

"–" = no thiazole peak (text: "absence of 4-methyl-5-(β-hydroxyethyl) thiazole-S³⁵ peak").

### 3.2 Text and figures [Q]

- pH 6.0, 15 min (Fig. 1): **only two radioactive peaks**, A (thiamine) and B (thiazole). pH 6.0 at
  30 min and 2 h: "similar". The thiazole's radioactivity "was nearly equal to thiamine-S³⁵ destroyed".
- pH 3.5 and 5.0: also only the two peaks; thiazole peaks at pH 3.5 smaller than at 5.0 and 6.0
  "but equal to thiamine-S³⁵ destroyed".
- pH 7.0, 15 min (Fig. 2): total radioactivity of all peaks **less than** that of unheated thiamine;
  difference attributed to H2S (citing Dwivedi & Arnold 1971) and other volatile ³⁵S compounds that
  escape during TLC. **The deficit is not printed** [!]. Small peaks C and D (Rf above the thiazole),
  possibly elemental sulfur and other non-polar compounds (Morfee & Liska 1972 [C]); a broad peak
  near thiamine, possibly condensation complexes of aldehydes from thiazole-ring opening, matching
  no available sulfur standard. pH 7.0 at 30 min: "similar to Figure 2".
- pH 8.0 (15 and 30 min): like pH 7.0 but **no thiazole peak** (Fig. 3).
- Chemistry (p. 887–888): pK₁ ≈ **4.8** (protonated aminopyrimidine, I → II); a second titration
  step with midpoint **pH 9.2**, two equivalents of base, "slow", forming the pseudobase (III), which
  ring-opens to the ionised thiol form (IV) [C: Williams & Ruehle 1935; Clarke & Gurin 1935;
  Watanabe & Asahi 1955]. Less breakdown at pH 3.5 than at 5.0/6.0 → the protonated form I is
  less labile than II; the mechanism at pH ≤ 6 is "the same for both forms": C–N bond cleavage of
  the methylene bridge → pyrimidine ("probably 2-methyl-4-amino-5-hydroxymethyl pyrimidine") +
  the thiazole.
- At pH 7.0 and 8.0 breakdown is "governed by the minor forms (III and IV)"; **"Hydrogen sulfide is
  probably the major degradation product of the thiazole moiety of thiamine under these
  conditions"** (citing 1971, not measured here). Authors' belief: **Ea(H2S removal) < Ea(bridge
  cleavage) at pH 7.0 and 8.0**, which would explain the missing thiazole at pH 8 and less of it at
  pH 7. H2S release pulls the II ⇌ III ⇌ IV equilibrium toward III/IV.
- Volatiles (Fig. 5, structures only, no amounts): from fragment B of the thiazole (the dashed line
  cuts the S–C2 and N–C4 bonds, leaving S, C4(CH3), C5 and the CH2CH2OH side chain, a 5-carbon unit).
  Read off the 400-dpi crop: **(I) 2-methyltetrahydrothiophen-3-one, (II) 2-methylthiophene,
  (III) 2-methyl-4,5-dihydrothiophene, (IV) 2-methylfuran, (V) 2-methyl-5-(methylthio)thiophene**
  (structure names mine, from the drawings [!]). "The reaction sequence ... is not known."
  Formation by cleavage of N–C or C–S, rearrangement and cyclisation of fragment B, with
  elimination of **one H2S or one water**. "Comparatively few degradation products were obtained
  from pH 8.0 ... indicating that energies of activation for the formation of these compounds ...
  is much higher than that for the main reaction of hydrogen sulfide elimination."

### 3.3 Inconsistencies found [!]

1. **Fig. 3 caption** reads "pH 6.0 autoclaved for 15 min", the same as Fig. 1, but the text cites
   Fig. 3 for the **pH 8.0** sample, and the trace has no thiazole (B) peak, which a pH 6.0 trace
   must have. The caption is wrong; Fig. 3 is pH 8.0.
2. **pH 3.5 / 15 min = 100.0 % thiazole** is the full mass balance to the stated precision; with 4.1 %
   destroyed, the peak-area method cannot resolve a ±1–2 % share. Treat the pH ≤ 6 "non-thiazole"
   shares (0–4 %) as **upper bounds of order the method's resolution**, not as measurements.
3. **pH 6.0 thiamine loss disagrees with Dwivedi & Arnold 1971** at the same temperature and
   buffer: 18.2 % at 30 min here vs **6.3 %** (1971 Fig. 2, 120 °C autoclave, 30 min). pH 7.0 agrees
   within ~10 points (46.5 % vs 51.1 %; 28.0 % vs 40.1 % at 15 min). Concentration differs 12×.
4. **pH 6.0 sulfur fate disagrees with Dwivedi & Arnold 1971 at 120 °C.** 1971 Fig. 4 reads
   22.5 mmol H2S/mol thiamine at pH 6.0, 120 °C, 30 min, against 6.3 % loss: [Z] a share of **~0.36**
   of degraded thiamine. Here, at pH 6.0, 121 °C, 30 min, **≤ 1.4 %** of the degraded thiamine's sulfur
   is anything other than the thiazole ([Z] ≤ 2.5 mmol/mol of initial thiamine). The ³⁵S balance is
   the more direct measurement of where the sulfur went; the 1971 autoclave H2S at pH ≤ 6.5 is the
   suspect number. Not resolvable from either paper.

---

## 4. DERIVED QUANTITIES [Z] (mine, not printed)

**Apparent first-order constants at 121 °C** (k = −ln(1 − destroyed)/t), split into the thiazole
(bridge-cleavage) channel and everything else (H2S + volatile S + peaks C/D + broad peak):

| pH | t (min) | k_total (min⁻¹) | k_thiazole (min⁻¹) | k_non-thiazole (min⁻¹) | non-thiazole S, % of initial thiamine |
|---|---|---|---|---|---|
| 3.5 | 15 / 30 | 2.79 / 2.71 × 10⁻³ | 2.79 / 2.67 × 10⁻³ | 0 / 4 × 10⁻⁵ | 0 / 0.12 |
| 5.0 | 15 / 30 | 5.92 / 5.46 × 10⁻³ | 5.90 / 5.40 × 10⁻³ | 2 / 5 × 10⁻⁵ | 0.03 / 0.15 |
| 6.0 | 15 / 30 / 120 | 6.73 / 6.70 / 4.40 × 10⁻³ | 6.46 / 6.60 / 4.30 × 10⁻³ | 2.7 / 0.94 / 0.92 × 10⁻⁴ | 0.38 / 0.25 / 0.86 |
| 7.0 | 15 / 30 | 2.19 / 2.09 × 10⁻² | 3.6 / 1.9 × 10⁻³ | **1.83 / 1.90 × 10⁻²** | 23.4 / 42.4 |
| 8.0 | 15 / 30 / 120 | 2.64 / 2.93 / 2.17 × 10⁻² | 0 | 2.64 / 2.93 / 2.17 × 10⁻² | 32.7 / 58.5 / 92.6 |

- **Thiamine loss is close to first order over 15–30 min at every pH** (k at 15 and 30 min within
  10 %); at 2 h k falls by ~35 % (pH 6) and ~25 % (pH 8). No pH drift was measured, so the cause is open.
- **pH slope of the non-thiazole channel, 6.0 → 7.0, 121 °C:** ×68 (15 min) to ×200 (30 min), i.e.
  **d log k / d pH ≈ 1.8–2.3**. Dwivedi & Arnold 1971 found ≈ 1.9–2.4 for the H2S channel between
  pH 6.5 and 7.0 at ~100 °C (`arnold1971_extraction.md` §4). Two independent methods and two
  temperatures, same slope of about two decades per pH unit.
- The thiazole channel at pH 7 (2–4 × 10⁻³ min⁻¹) is **lower** than at pH 6 (6.5 × 10⁻³), and the
  absolute thiazole stays at ~4–5 % of initial thiamine from 15 to 30 min. Either the bridge channel
  is suppressed at pH 7 or the thiazole is itself consumed there (and fully at pH 8). The paper
  cannot tell these apart; the authors favour competition.
- **Temperature dependence: none inside the paper** (one temperature). A cross-paper two-point
  estimate at pH 7 (1971 sealed tube ~96–100 °C, k = 5.0–5.8 × 10⁻³ min⁻¹; here 121 °C,
  2.1–2.2 × 10⁻² min⁻¹; ratio 3.6–4.4) gives **Ea(total loss, pH 7) ≈ 62–86 kJ/mol**, consistent with
  the 1971 internal estimate (78–95). The same exercise at pH 6.0 (1971 loss 0.2–1.9 % in 30 min)
  gives a ratio of 10–100 and Ea 110–270 kJ/mol, which only shows that the 1971 pH 6 loss is at the
  floor of its method. **Neither is defensible as a prior:** different concentrations (12×),
  uncounted autoclave come-up, 25 °C pH labels at two temperatures.
- **Upper bound on the H2S share at 121 °C:** pH 7.0 ≤ 0.84 (15 min) – 0.91 (30 min) of degraded
  thiamine; pH 6.0 ≤ 0.04 (15 min), ≤ 0.014 (30 min), ≤ 0.021 (2 h); pH 5.0 ≤ 0.004–0.010; pH 3.5
  ≤ 0.015. Compare 1971 at 120 °C, pH 7: H2S / degraded = 0.31 (5 min), 0.22 (15 min), 0.17 (30 min)
  ([Z] from 1971 Figs. 3 and 5), so at 120 °C free H2S is about a fifth to a third of the non-thiazole
  sulfur, and the rest is other S compounds (C, D, the broad peak, elemental S, volatiles).

---

## 5. WHAT THE MODEL CAN TAKE FROM IT

**Context.** B50 (`THIAMINE_H2S_SHARE_POINTS`, `thiamine_h2s_share` in
`src/kinetic_core/parameters_sulfur.py`) puts a share of thiamine's thermal loss into H2S: 0.39 at
pH 7.0, 0.04 at pH 6.5, log-linear between, **zero below pH 6.25**, constant with temperature, from
Dwivedi & Arnold 1971 at ~100 °C. `r_thi_h2s = share/(1 − share) × (r_thi_hmp + r_thi_mesh)`.

**Verdict: (b), a sulfur-fate fraction at one temperature (121 °C) with pH dependence; it bounds
rather than measures the H2S share. Not (a).**

1. **B50's floor is confirmed, and extended to 121 °C.** At pH ≤ 6.0, ≥ 96 % (≥ 98.6 % at 30 min) of
   degraded thiamine's sulfur is the intact thiazole at 121 °C, so the H2S share there is
   ≤ 0.01–0.04, consistent with zero below 6.25. This **contradicts the 1971 autoclave point at
   pH 6.0 (share ≈ 0.36)**, which B50 already ignores.
2. **The pH 7 share is neither confirmed nor refined numerically.** Here only the bound ≤ 0.84–0.91
   at 121 °C. Combined with 1971's 120 °C data (0.17–0.31, falling with time), **holding 0.39 at
   120 °C overstates the free-H2S share by ~1.3–2.3×** at 30 min. The direction of the temperature
   effect on the *share* is down, not flat, but two papers at two concentrations cannot fix a slope.
3. **The pH shape is confirmed:** about two decades per pH unit in the non-thiazole channel (6 → 7,
   121 °C), matching the 1971 slope (6.5 → 7, 100 °C). The 6.25 floor plus log-linear rise is a fair
   summary; no data between 6.0 and 7.0 here.
4. **pH scale:** 25 °C preparation labels of 0.1 M phosphate. No in-situ value. Same convention as
   1971, so B51's `label_scale_tables` (read the tables at the cooled pH) applies to both.
5. **Structural gap in the thiamine steps (main contradiction):** the model has **no
   bridge-cleavage channel** (THI → thiazole + pyrimidine; no thiazole species). At pH 3.5–6 that
   channel is ≥ 96 % of thiamine loss at 121 °C, with k ≈ 3–7 × 10⁻³ min⁻¹. So (i) the model's
   thiamine depletion at pH ≤ 6 is whatever `k_thi_hmp + k_thi_mesh` were fitted to, not the measured
   loss; (ii) **at pH ≤ 6 essentially none of the degraded thiamine's sulfur becomes HMP/MFT, MeSH
   or H2S**: ≤ 4 % in total at 121 °C. If `k_thi_hmp` was fitted at pH ≤ 6 to MFT yields, the MFT
   must come from a sub-percent channel. That is not a contradiction, but the fitted constant
   cannot be read as "the" thiamine loss rate. (iii) B50 multiplies the share by
   `r_thi_hmp + r_thi_mesh`, which at pH 7 approximates the non-thiazole channel (84–91 % of loss);
   the 1971 share is defined against **total** loss, so the mismatch is ≤ 10–20 % at pH 7. That is
   inside B50's band.
6. **Thiamine route to MFT (d):** nothing direct. 2-Methyl-3-furanthiol is not mentioned. The
   volatiles all come from **fragment B, the 5-carbon S-containing unit** (C4-CH3, C5, CH2CH2OH) that
   is the skeleton of the model's HMP (5-hydroxy-3-mercapto-2-pentanone). 2-Methylfuran is fragment
   B *minus H2S* (a desulfurised competitor of MFT); the thiophenes keep the sulfur. Authors say
   these form with a **higher Ea than H2S elimination** and are scarce at pH 8 (qualitative only).
7. **Usable kinetics:** the k_total table above (121 °C, five pH values, 25 °C labels, 1.5 mM). It can
   check a thiamine loss rate if the model ever adds the bridge-cleavage channel. It is not a fit
   target as it stands (no matching species).

**Not usable:** any Ea (none printed; the cross-paper estimates are not priors); any H2S yield
(not measured).
