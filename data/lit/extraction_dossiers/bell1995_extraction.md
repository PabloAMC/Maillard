# Bell 1995/1996 — COMPLETE TRANSCRIPTION (Wave K4c)

**Full extraction of every number in `data/articles/bell1995.pdf`.**
Extraction date 2026-08-28. Tables 1–3 and Figs 1–4 re-read off 300 dpi page rasters; the
Acrobat-Capture OCR text layer was used only for prose and was **not** trusted for digits
(it mangled the Tg column — see §0b).

---

## 0. PAPER IDENTITY — MATCHES THE WAVE BRIEF, WITH ONE YEAR CAVEAT

| field | value |
|---|---|
| Author | **Leonard N. Bell** (sole author) |
| Title | **"Kinetics of non-enzymatic browning in amorphous solid systems: distinguishing the effects of water activity and the glass transition"** |
| Journal | **Food Research International, Vol. 28, No. 6, pp. 591–597** |
| **Year printed on the page** | **1996** — *"Food Research International, Vol. 28, No. 6, pp. 591–597, **1996**"* and *"Copyright © **1996** Canadian Institute of Food Science and Technology"* (p. 591, verified on raster) |
| Received / accepted | **Received 18 May 1995; accepted 24 June 1995** (p. 597) |
| DOI | **10.1016/0963-9969(95)00052-6** ✔ matches the brief |
| PII / ISSN line | `0963-9969(95)00052-6`; `0963-9969/96 $15.00 + 0.00` |
| Affiliation | Nutrition and Food Science, **Auburn University, AL 36849, USA**; e-mail lbell@humsci.auburn.edu |
| Funding | Alabama Agricultural Experiment Station (**AAES Journal No. 10-955008**) + Auburn University Grant-in-Aid |
| PDF | 7 pages, Acrobat 3.0 Capture (OCR over a clean scan); text layer usable for prose, unreliable for tabular digits |

**⚠ Year caveat for the citation graph.** The repo's filename and any `bell1995` key imply
1995. **The article itself is stamped 1996** (volume 28 issue 6, 1996 per the running head and
copyright line), while the DOI/PII prefix and the submission dates are 1995. Both years appear
in the literature. **Cite as Bell, L. N. (1996) *Food Res. Int.* **28**, 591–597, DOI
10.1016/0963-9969(95)00052-6** and keep `bell1995` only as an internal key, with an alias.

### 0a. What the brief asked for vs what the paper is

The brief asked for *"all rate constants vs a_w/Tg and the deconfounding design"*. **Both are
here in full and nothing is missing.** The paper contains **17 pseudo-zero-order rate
constants** (13 in Table 2 + 4 in Table 3), each with a **95 % confidence limit** and an **R²**,
plus the complete composition table (Table 1) that makes the deconfounding design auditable.
There are **no** rate constants at any temperature other than **25 °C**, and therefore **no
activation energy anywhere in this paper.**

### 0b. Text-layer damage — why every digit was re-read off the raster

The OCR text layer corrupts the Table 2 Tg column systematically: it renders `8` as `08`,
`15` as `15` but `−7` as `-07`, `−30` as `-30`, and drops several minus signs into the
adjoining column. It also produces `PIT-K30`, `PVP-Kl5`, `PVP-LM W`, `Table 3 … 0.00 I3`, and
`0.0016` variants. **Every value in §3–§5 below was read off the 300 dpi page raster of PDF
pages 4 and 5 (printed pp. 594–595). Not one cell was unreadable.**

### 0c. Provenance codes

- **[M]** measured / reported directly by the author in a table
- **[F]** fitted by the author (the pseudo-zero-order regression parameters he printed)
- **[Z]** derived by me from the paper's own numbers — never printed by the author
- **±** here always means the author's printed **95 % confidence limit** on the regression
  slope, *not* a replicate SD. **Replicate structure (n) is never stated anywhere.**

---

## 1. THE ONE-PARAGRAPH ANSWER

All **17 rate constants** are legible and transcribed, with their 95 % CLs and R². The units
are **OD/g/d** — *optical density at 420 nm per gram of dry solids per day* — which makes them
**pseudo-zero-order pigment-accumulation rates in an extract, not molar rates**, and therefore
**not directly convertible to a molar Maillard flux without the extraction geometry** (200 mg
solid → 4 mL water, §2c). The deconfounding design is genuinely clean and is the paper's whole
point: **three PVP molecular weights with near-identical moisture sorption isotherms give three
different Tg at the same a_w**, so a_w and Tg can be varied independently — the first time this
had been done for browning. The headline results: at **a_w 0.33 all three systems are glassy
and their rate constants are statistically indistinguishable (p > 0.05)**; at **a_w 0.54 they
straddle the glass transition and differ significantly (p ≤ 0.05)** with the glassy system
slowest; the glassy→rubbery jump is **"7-fold"** as stated by the author (§6a shows this is a
max/min comparison — the mean ratio is 3.0×); and the **largest single effect in the whole
paper is neither a_w nor Tg but reactant concentration**, which spans **0.0034 → 0.068 OD/g/d,
a 20× range** (Table 3).

---

## 2. SYSTEM COMPOSITION AND METHOD — applies to every number below

### 2a. The matrix and the deconfounding logic (p. 592, "Description of model system") **[M]**

| variable | value as printed |
|---|---|
| Inert matrix | **Polyvinylpyrrolidone (PVP)** — *"a polar, water soluble polymer"* |
| Molecular weights | **PVP-K30 (MW ≈ 40 000)**, **PVP-K15 (MW = 10 000)**, **PVP-LMW (MW < 3500)** |
| **Deconfounding claim** | Bell & Hageman (1995) showed the moisture sorption isotherms of all three *"were virtually identical (i.e. at a given water activity all three molecular weight polymers had similar moisture contents)"*. **At a given moisture content each has a different Tg.** Therefore a_w and moisture can be held constant while Tg is moved by MW alone. |
| Reactants | **D-glucose + glycine** |
| Buffer | **0.1 M phosphate, pH 7** |
| **Design intent** | *"such that the final internal buffer concentration after dehydration and moisture sorption would be close to **0.1 molal**"* and *"the final internal reactant concentration … would be close to **1 molal** for both glucose and glycine"* — i.e. **the reactant molality in the sorbed aqueous phase, not the mass ratio, is what is held fixed.** This is the explicit correction to Karmas & Karel (1994) and Buera & Karel (1995), *"both of these studies used a constant mass ratio of reactants which would result in different soluble reactant concentrations after moisture sorption at the different water activities"* (p. 592). |

### 2b. Sample preparation (p. 592) **[M]**
PVP-K15 and PVP-K30 dissolved in water, neutralised to **pH 7 with NaOH**; **dialysed through
3500 MWCO bags** into a single **3 L** volume of water, itself neutralised to pH 7 — the
dialysate is what becomes **PVP-LMW**. Both retentate and dialysate **lyophilised 4 d at room
temperature under < 20 mtorr**, stored over anhydrous CaSO₄. Lyophilised powders redissolved
in purified water; **glucose and glycine each dissolved separately at 1.0 M**; a measured volume
of 0.1 M pH 7 phosphate plus the two reactant solutions added; the solution **flash-frozen by
dripping into liquid nitrogen through a 22 gauge needle**; frozen pellets **lyophilised ≥ 2 d at
room temperature under < 20 mtorr**; dried pellets stored in a vacuum desiccator over anhydrous
CaSO₄.

### 2c. Kinetic assay (p. 593, "Analysis and kinetic study") **[M]**
| variable | value as printed |
|---|---|
| Sample mass | **200 mg** per vial, multiple vials per condition |
| Equilibration | **closed desiccators at 25 °C for one week** over saturated salt solutions |
| Salts named | **0.33 [MgCl₂], 0.54 [Mg(NO₃)₂], 0.76 [NaCl]** — *only these three salts are named; the salts for a_w 0.44, 0.65/0.66 and 0.85 are never given* ⚠ |
| Storage temperature | **25 °C throughout. There is no second temperature anywhere in this paper.** |
| t = 0 definition | *"the first sample being designated as time zero"* — i.e. **t = 0 is after the one-week equilibration, not at lyophilisation.** Any browning during equilibration is baked into the intercept, not the slope. |
| Moisture | determined **gravimetrically** as moisture sorbed into the dry system |
| Extraction | **4 mL purified water added to the 200 mg sample**, dissolved completely, filtered through a **0.20 μm syringe filter** |
| Read-out | **absorbance at 420 nm**, Spectronic 20 (Bausch & Lomb, Rochester NY) |
| Kinetic model | **pseudo zero order**, justified by citation (Mizrahi et al. 1970; Warmbier et al. 1976a,b; Labuza & Saltmarch 1981; Labuza & Baisier 1992) — **not tested against alternatives in this paper** |
| Regression | **computerised least squares with 95 % confidence limits**, per Labuza & Kamman (1983) |
| Longest time point | **≈ 67 days** (from the Fig. 1 x-axis, which runs 0–80 d with the last datum at ≈ 67 d) |

### 2d. DSC (p. 593) **[M]**
TA Instruments **2910** DSC; **4–11 mg** hermetically sealed in aluminium pans;
**scan rate 5 °C/min**; duplicates and rescans used to verify the endothermic baseline shift.
**Tg reported as the ONSET temperature.** For the a_w 0.96 solution, conventional DSC failed
and Tg was obtained in **modulation mode, 5 °C/min, amplitude ± 1 °C over a 60 s period**.
Author's note: results *"were similar to those obtained previously"* (Karmas et al. 1992;
Buera et al. 1992; Bell & Hageman 1995).
**No uncertainty is given on any Tg. All Tg values are printed as whole degrees Celsius.**

### 2e. The "other experiments" — how the concentration arm was built (p. 593) **[M]**
This paragraph is the key to Table 3 and is easy to misread:
- PVP-LMW was **also** prepared for and stored at **a_w 0.44, 0.65 and 0.85**.
- **PVP-K30 samples intended for a_w 0.33 and 0.76 were instead equilibrated at a_w 0.76 and
  0.54 respectively**; **PVP-K15 samples intended for a_w 0.76 were equilibrated at a_w 0.33.**
  Formulating for one a_w and equilibrating at another is exactly how the internal molality is
  moved off 1 molal while keeping the same solid.
- A PVP-LMW sample at **a_w 0.85 with a reactant concentration of 0.77 molal**.
- *"These experiments were to evaluate the importance of the reactant concentration in the
  internal aqueous microenvironment."*
- Separately, a **solution** control: **1 molal glucose + 1 molal glycine + 0.1 molal phosphate
  pH 7 + 1 % PVP-LMW**, measured a_w **0.96** on an **Aqualab** instrument (Decagon, Pullman WA).

---

## 3. TABLE 1 — COMPOSITION OF THE MODEL SYSTEMS **[M]**

**Anchor: Table 1, p. 593 (PDF page 3).** Printed title: *"Composition of model systems
(PVP with 0.1 M phosphate buffer pH 7)"*. Column headers exactly as printed:
`Sample | aw | Moisture (% db) | PVP (g) | Water^a (g) | Added solutions^b (μl) | Buffer conc.^c (molal) | Reactant conc.^c (molal)`.
Footnotes as printed: **ᵃ Sorbed water. ᵇ Volume of 0.1 M buffer, 1.0 M glycine, and 1.0 M
glucose. ᶜ Internal concentration after moisture sorption assuming complete dissolution in
aqueous phase.**

| Sample | a_w | Moisture (% db) | PVP (g) | Water (g) | Added solutions (μl) | Buffer conc. (molal) | **Reactant conc. (molal)** |
|---|---|---|---|---|---|---|---|
| PVP-LMW | 0.33 | 10.0 | 3.02 | 0.302 | 320 | 0.11 | **1.1** |
| PVP-K15 | 0.33 | 10.1 | 3.02 | 0.304 | 304 | 0.10 | **1.0** |
| PVP-K30 | 0.33 | 11.4 | 2.99 | 0.341 | 333 | 0.098 | **0.98** |
| PVP-LMW | 0.44 | 14.4 | 2.34 | 0.337 | 339 | 0.10 | **1.0** |
| PVP-LMW | 0.54 | 18.5 | 3.02 | 0.558 | 612 | 0.11 | **1.1** |
| PVP-K15 | 0.54 | 18.3 | 2.98 | 0.543 | 579 | 0.11 | **1.1** |
| PVP-K30 | 0.54 | 18.5 | 3.02 | 0.558 | 604 | 0.11 | **1.1** |
| PVP-LMW | 0.66 | 26.1 | 3.39 | 0.884 | 949 | 0.11 | **1.1** |
| PVP-LMW | 0.76 | 38.4 | 3.02 | 1.161 | 1130 | 0.097 | **0.97** |
| PVP-K15 | 0.76 | 36.6 | 3.06 | 1.121 | 1136 | 0.10 | **1.0** |
| PVP-K30 | 0.76 | 36.1 | 3.04 | 1.099 | 1104 | 0.10 | **1.0** |
| PVP-LMW | 0.85 | 57.7 | 2.27 | 1.310 | 1272 | 0.097 | **0.97** |

**Twelve rows. The a_w 0.96 solution has no Table 1 row.**

**Verification of the deconfounding claim [Z].** At each a_w the moisture contents of the three
MWs are: 0.33 → 10.0 / 10.1 / 11.4 % db (spread **14 %** relative); 0.54 → 18.5 / 18.3 / 18.5
(spread **1 %**); 0.76 → 38.4 / 36.6 / 36.1 (spread **6 %**). **The isotherms really are
near-identical, and the reactant molalities really are held at 0.97–1.1** (spread ± 6 %) across
every row. **The design does what it claims.** The one soft spot is a_w 0.33, where PVP-K30
holds 14 % more water than PVP-LMW — and a_w 0.33 is precisely the condition where the author
reports no significant rate difference, so the residual moisture mismatch is not doing work.

**Sorbed-water consistency check [Z]:** water (g) ÷ PVP (g) reproduces the printed % db to
≤ 0.3 percentage points in every row (e.g. LMW 0.33: 0.302/3.02 = 10.0 % ✔; LMW 0.85:
1.310/2.27 = 57.7 % ✔). **Table 1 is arithmetically self-consistent.**

**⚠ a_w label inconsistency [Z].** The Methods paragraph (p. 593) says PVP-LMW was also stored
at *"water activities 0.44, **0.65** and 0.85"*, but Tables 1 and 2 both print **0.66**. A
0.01 difference is immaterial physically but matters for machine-matching. **Use 0.66 (the
tables).**

---

## 4. TABLE 2 — 13 RATE CONSTANTS AT 1 MOLAL REACTANT **[M]/[F]** — the deconfounding core

**Anchor: Table 2, p. 594 (PDF page 4).** Printed title: *"Rate of brown pigment formation in
PVP systems containing 1 molal reactant concentration as influenced by water activity (a_w) and
glass transition temperature (T_g)"*. Column headers as printed:
`System | a_W | T_g (°C) | State of the system at 25°C | Rate constant w/ 95% C.L. (OD/g/d) | R²`.
**Units quoted exactly as printed: T_g in `°C`, k in `OD/g/d`.** All 13 rows legible.
`T − T_g` is **[Z]**, computed as 25 °C − T_g (the abscissa of Fig. 2).

| System | a_w | **T_g (°C)** | **T − T_g (°C) [Z]** | State at 25 °C | **k ± 95 % C.L. (OD/g/d)** | R² | rel. CL **[Z]** |
|---|---|---|---|---|---|---|---|
| PVP-LMW | 0.33 | **38** | −13 | **Glassy** | **0.0053 ± 0.0016** | 0.938 | 30 % |
| PVP-K15 | 0.33 | **43** | −18 | **Glassy** | **0.0073 ± 0.0023** | 0.932 | 32 % |
| PVP-K30 | 0.33 | **55** | −30 | **Glassy** | **0.0054 ± 0.0008** | 0.986 | 15 % |
| PVP-LMW | 0.44 | **20** | +5 | Rubbery | **0.027 ± 0.004** | 0.982 | 15 % |
| PVP-LMW | 0.54 | **8** | +17 | Rubbery | **0.034 ± 0.001** | 0.998 | 3 % |
| PVP-K15 | 0.54 | **15** | +10 | Rubbery | **0.024 ± 0.004** | 0.982 | 17 % |
| PVP-K30 | 0.54 | **28** | −3 | **Glassy** | **0.014 ± 0.003** | 0.961 | 21 % |
| PVP-LMW | 0.66 | **−7** | +32 | Rubbery | **0.023 ± 0.003** | 0.991 | 13 % |
| PVP-LMW | 0.76 | **−30** | +55 | Rubbery | **0.023 ± 0.002** | 0.996 | 9 % |
| PVP-K15 | 0.76 | **−26** | +51 | Rubbery | **0.024 ± 0.003** | 0.990 | 12 % |
| PVP-K30 | 0.76 | **−14** | +39 | Rubbery | **0.013 ± 0.002** | 0.984 | 15 % |
| PVP-LMW | 0.85 | **−62** | +87 | **Viscous solution** | **0.035 ± 0.004** | 0.991 | 11 % |
| PVP-LMW | 0.96 | **−76** | +101 | **Solution** | **0.032 ± 0.008** | 0.972 | 25 % |

**The three-state vocabulary matters:** the author uses **four** state labels — Glassy, Rubbery,
Viscous solution, Solution — not two. Only the a_w 0.85 and 0.96 rows are liquid.

### 4a. The deconfounding result, stated as the author states it **[M]**

- **At a_w 0.33 (all three glassy):** *"the rates of brown pigment formation in the systems were
  not significantly different (**P > 0.05**) from each other. However, the browning reaction did
  occur, indicating that even in an 'immobilized' glassy state, sufficient mobility exists for
  two molecules to diffuse and react."* (p. 594)
  → **k = 0.0053 / 0.0073 / 0.0054 across Tg = 38 / 43 / 55 °C. A 17 °C Tg span produces no
  resolvable rate difference *inside* the glass.**
- **At a_w 0.54 (straddling the transition):** *"the PVP-K30 and PVP-K15 systems were on either
  side of the T_g (glassy and rubbery, respectively) while PVP-LMW was more rubbery. The rate
  constants at a water activity of 0.54 were significantly different (**P ≤ 0.05**) from each
  other with the slowest reaction occurring in the glassy system and the fastest reaction
  occurring in the most rubbery system."* (p. 594)
  → **k = 0.014 (glassy, T−Tg = −3) < 0.024 (rubbery, +10) < 0.034 (rubbery, +17). Same a_w,
  same moisture (18.3–18.5 % db), same molality (1.1) — a 2.4× rate span attributable to Tg
  alone. This is the paper's single cleanest measurement.**
- **At a_w 0.76 (all three rubbery):** *"all systems were in a rubbery state and the rate of
  browning was about the same as at a water activity of 0.54 for a given type of PVP. However,
  the slowest browning rate again occurred in the system with the highest T_g."*
  → **0.023 / 0.024 / 0.013; the K30 system (highest Tg, −14 °C) is again ≈ 1.8× slower.**
- **Overall:** *"the rate of brown pigment formation is influenced more significantly by the
  state of the system than its water activity."* (p. 594)

### 4b. The "7-fold" claim — what it actually is **[M] + [Z]**

Printed (p. 594): *"As the system changed from a glassy state to a rubbery state, the reaction
rate increased **7-fold**, which is less than the **100-fold** increase in mobility of the ESR
probe (Roozen et al., 1991)."*
**[Z] The 7-fold is a max/min comparison, not a mean ratio.** Across all 13 rows,
**k_max / k_min = 0.035 / 0.0053 = 6.6×** — that is the "7-fold". The **mean glassy k =
0.0080 OD/g/d (4 rows) and mean rubbery k = 0.0240 OD/g/d (7 rows), a ratio of 3.0×**. The
**within-a_w-0.54 glassy→most-rubbery ratio is 2.4×**. **All three numbers are defensible
readings of "the glass transition effect"; they differ by 2.7×. Use 2.4× if you want the
deconfounded value, 6.6× only if you are reproducing the author's sentence.**

### 4c. Figure 2 — the k vs T−T_g curve **[M]**
**Anchor: Fig. 2, p. 594.** *"Rate constants for brown pigment formation in polyvinylpyrrolidone
containing 1 molal reactant concentration and 0.1 molal phosphate buffer at pH 7 and 25 °C as a
function of the distance from the glass transition temperature."* y-axis **Rate Constant
(OD/g/d), 0.00–0.05**; x-axis **T−T_g, −40 to +120**; a vertical divider at T−T_g = 0 labelled
**Glassy** (left) / **Rubbery** (right). Error bars are drawn (the 95 % CLs of Table 2).
**Fig. 2 plots exactly the 13 Table 2 rows against my [Z] T−T_g column and adds no new data**
— which independently confirms the T−T_g arithmetic and hence the Tg transcription.

**[M] The author's own caveats on Fig. 2 (pp. 594–595):**
- *"The minima in the rate constant around **T − T_g = 40** are not well understood, but could
  result from the large quantity of water present; the mass action of water could reduce the
  rate of those reaction pathways from which water is released, such as the formation of the
  glycosylamine, the dehydration to yield the Schiff base, and the dehydration to yield
  reductones (Hodge, 1953; Labuza & Baisier, 1992). If this was the case, however, one would
  expect the rate to continue to decrease with larger values of T − T_g … It seems more likely
  that in the rubbery state, the rates of browning are roughly constant."*
  → **[Z] The "minimum at T−T_g ≈ +40" is the PVP-K30 a_w 0.76 point (k = 0.013, T−Tg = +39).
  It is a single system, and it is the same PVP-K30 that is also slowest at a_w 0.54 and 0.33.
  The dip is better read as a persistent PVP-K30 offset than as a T−T_g feature.**
- *"Figure 2 suggests that the glass transition influences the reaction rate, but not to the
  extent expected for a bimolecular reaction such as browning … The less-than-dramatic effect
  of the glass transition on the rate of non-enzymatic browning suggests that other factors,
  perhaps water activity, are also influencing the reaction."*

### 4d. Figure 3 — k vs a_w at fixed 1 molal **[M]**
**Anchor: Fig. 3, p. 595.** Same caption but *"as influenced by water activity"*; x-axis
**Water Activity 0.2–1.0**, y-axis **0.00–0.05 OD/g/d**. Plots the **seven PVP-LMW rows** of
Table 2 (a_w 0.33, 0.44, 0.54, 0.66, 0.76, 0.85, 0.96) with their CLs, plus a hand-drawn
sigmoid that rises to a plateau at a_w ≈ 0.5.
**[M]** *"This figure displays an increase in the browning rate from a water activity of 0.33 to
0.54 … Above a water activity of 0.54, the reaction rate appears to form a rough plateau
(Fig. 3), which is appropriate because no reactant dilution would be occurring due to the
constant reactant concentration of 1 molal."*
**[Z] Quantified: 0.0053 → 0.027 → 0.034 over a_w 0.33 → 0.44 → 0.54 (a 6.4× rise), then
0.034 / 0.023 / 0.023 / 0.035 / 0.032 over 0.54 → 0.96 (flat to ± 25 %). The classic
bell-shaped a_w curve with a maximum at 0.55–0.7 is ABSENT — deliberately so, because the
descending limb of that classic curve is reactant dilution, which this design removes.**
**This is the paper's second headline and the more important one for any model that carries an
a_w-dependence term:** the textbook browning-vs-a_w maximum is **not intrinsic**; it is an
artefact of holding mass ratio rather than molality constant.
**[M]** *"A maximum in the browning rate between a water activity of 0.55 and 0.7 has been
frequently demonstrated (Loncin et al., 1968; Warmbier et al., 1976a; Karel, 1984; Karmas &
Karel, 1994). This maximum has been attributed to mobility limitations with decreasing water
activity below the maximum and increasing reactant dilution at water activities above the
maximum (Eichner & Karel, 1972; Labuza, 1980)."*

---

## 5. TABLE 3 — 4 RATE CONSTANTS AT NON-1-MOLAL REACTANT **[M]/[F]**

**Anchor: Table 3, p. 595 (PDF page 5).** Printed title: *"Rate of brown pigment formation in
PVP systems containing various reactant concentrations"*. Headers as printed:
`System | a_W | Reactant concentration (molal) | State of the system at 25°C | Rate constant w/ 95% C.L. (OD/g/d) | R²`.
All 4 rows legible.

| System | a_w | **Reactant conc. (molal)** | State at 25 °C | **k ± 95 % C.L. (OD/g/d)** | R² |
|---|---|---|---|---|---|
| PVP-K15 | **0.33** | **4.9** | Glassy | **0.068 ± 0.011** | 0.982 |
| PVP-K30 | **0.56** | **2.31** | Glassy | **0.045 ± 0.007** | 0.981 |
| PVP-K30 | **0.76** | **0.28** | Rubbery | **0.0034 ± 0.0013** | 0.932 |
| PVP-LMW | **0.85** | **0.77** | Viscous solution | **0.013 ± 0.001** | 0.993 |

**⚠ Three inconsistencies between Table 3 and the surrounding text — all flagged [Z]:**
1. **a_w = 0.56 appears only here.** Every other a_w in the paper is 0.33 / 0.44 / 0.54 / 0.66 /
   0.76 / 0.85 / 0.96, and the Methods say the PVP-K30 concentration sample was *"equilibrated
   at water activity 0.54"*. **0.56 is almost certainly a typo for 0.54**, but it is what is
   printed. Read it as **0.54 (typo)** and flag.
2. **The running text says 0.29 molal; Table 3 says 0.28 molal.** p. 595: *"the concentration of
   reactants in the 0.76 water activity system was **0.29 molal**"* and *"at **0.29 molal** the
   rate constant decreased to 0.0034"*. **Table 3 prints 0.28.** Use **0.28** (the table).
3. **The PVP-K30 a_w 0.76 / 0.28 molal row is stated to be glassy nowhere and rubbery here** —
   consistent with Table 2's PVP-K30 a_w 0.76 Tg = −14 °C. No conflict; noted for completeness.

**[M] The author's own paired comparisons (p. 595) — the concentration effect isolated:**
- *"In PVP-K15 at water activity 0.33, the browning rate constant was **0.0073 OD/g/d at a
  reactant concentration of 1 molal** whereas at a reactant concentration of **4.9 molal**, the
  rate constant increased to **0.068 OD/g/d**."* → **9.3× for a 4.9× concentration change.**
- *"In PVP-K30 at water activity 0.76, the rate constant was **0.013 OD/g/d at 1 molal** while
  at **0.29 [table: 0.28] molal** the rate constant decreased to **0.0034**."* → **3.8× for a
  3.6× concentration change.**
- *"An interesting observation from the current study was that the **fastest browning rate
  occurred in the glassy state at the lowest water activity (0.33)** while the **slowest reaction
  occurred in a rubbery system at water activity 0.76**, which is contrary to that expected for a
  diffusionally-controlled reaction. However, these results can be explained by the reactant
  concentrations in each case (Table 3)."*

### 5a. Reaction order in reactant concentration **[Z]** — the paper stops short of computing this

The two matched pairs above are the only true like-for-like concentration comparisons (same
polymer, same a_w, same state):

| pair | conc. ratio | k ratio | **apparent order n = ln(k₂/k₁)/ln(c₂/c₁)** |
|---|---|---|---|
| PVP-K15, a_w 0.33, 1.0 → 4.9 molal | 4.90 | 9.32 | **1.40** |
| PVP-K30, a_w 0.76, 0.28 → 1.0 molal | 3.57 | 3.82 | **1.05** |

**Apparent order 1.05–1.40 in *each* reactant (glucose and glycine move together, so this is
the order in the pair).** For a genuinely bimolecular initiation step, holding [Glc] = [Gly] = c
should give k ∝ c², i.e. n = 2. **Neither pair reaches 2.** Whatever the repo assumes about
browning order in a reduced-moisture matrix, this is the only published deconfounded
measurement of it and it says **sub-bimolecular, n ≈ 1.0–1.4**.

### 5b. Figure 4 — k vs [Glucose]×[Glycine] **[M]**
**Anchor: Fig. 4, p. 596.** *"Rate constants for brown pigment formation in polyvinylpyrrolidone
at pH 7 and 25 °C as a function of reactant concentration."* x-axis **[Glucose]×[Glycine]
(molal²)**, with a **broken axis: 0, 1, 2, 3, 4, 5, then a break, then 20, 25**; y-axis
**0.00–0.08 OD/g/d**. It plots all 17 rate constants (Tables 2 + 3) against the product of the
two molalities.
**[M]** *"The spread of the data around 1 molal² can be attributed to the different water
activities and glass transition temperatures leading to different rate constants. Despite the
dispersed data, the trend is clear that reactant concentration has a significant impact on the
rate of browning, which is consistent with that found in solution (Baisier & Labuza, 1992)."*
**[Z] Quantified: over all 17 points, Pearson r(k, [Glc]×[Gly]) = 0.783; a log–log regression
gives slope 0.53 with R² = 0.51** — i.e. **the product-of-concentrations axis explains only about
half the variance, and the vertical scatter at [Glc][Gly] ≈ 1 molal² spans 0.0053 → 0.035
(6.6×), which is precisely the Tg/a_w effect.** Fig. 4 therefore does **not** collapse the data;
concentration and state are **both** required.

---

## 6. THE COMPLETE PICTURE OF EFFECT SIZES **[Z]** — what dominates what

Ranked by the rate-constant span each factor produces, at 25 °C and pH 7:

| factor | held constant | span in k (OD/g/d) | **fold** |
|---|---|---|---|
| **Reactant molality, 0.28 → 4.9** | (spans a_w and state) | 0.0034 → 0.068 | **20.0×** |
| Reactant molality, matched pair (K15, a_w 0.33) | polymer, a_w, state | 0.0073 → 0.068 | **9.3×** |
| Glass transition, max/min across all rows | — | 0.0053 → 0.035 | **6.6×** (the author's "7-fold") |
| Water activity, 0.33 → 0.54 at 1 molal (LMW) | polymer, molality | 0.0053 → 0.034 | **6.4×** |
| Reactant molality, matched pair (K30, a_w 0.76) | polymer, a_w, state | 0.0034 → 0.013 | **3.8×** |
| Glass transition, glassy vs rubbery means | — | 0.0080 → 0.0240 | **3.0×** |
| **Glass transition, fully deconfounded (a_w 0.54)** | **a_w, moisture, molality all fixed** | 0.014 → 0.034 | **2.4×** |
| Water activity, 0.54 → 0.96 at 1 molal (LMW) | polymer, molality | 0.023 → 0.035 | **1.5×** (i.e. flat) |
| Tg within the glass (a_w 0.33, Tg 38→55 °C) | a_w, moisture, molality | 0.0053 → 0.0073 | **1.4×, n.s. (p > 0.05)** |

**The single most important line for a Maillard model is the last-but-two: with a_w, moisture
and reactant molality genuinely held fixed, the entire glass-transition effect on browning rate
is 2.4×.** Every larger number quoted in the glass-transition literature — including this
paper's own "7-fold" — mixes in a_w or concentration.

**The author's own conclusion (p. 596), transcribed in full:**
> *"Several important conclusions can be derived from this study. Chemical reactions, even
> where molecular diffusion and collision are required, can occur in a glassy matrix. While both
> water activity and glass transition influence the rate of brown pigment formation, the state
> of the system as determined by the glass transition temperature is the rate-limiting factor
> with water activity playing a lesser role. However, perhaps more critical than the state of
> the system, with respect to brown pigment formation as well as other bimolecular reactions, is
> the reactant concentration within the aqueous phase of reduced-moisture solid system. Thus,
> when dealing with a complex reaction, such as non-enzymatic browning, broad generalizations
> should be avoided as to what factors control the reaction rate unless the properties of the
> experimental system are well defined."*

---

## 7. FIGURE 1 — the only raw kinetic trace shown **[M]**

**Anchor: Fig. 1, p. 594.** *"Pseudo-zero order plot of brown pigment formation in
polyvinylpyrrolidone containing 1 molal reactant concentration and 0.1 molal phosphate buffer
at pH 7, a_w 0.54 and 25 °C."*
y-axis **Absorbance/g dry solids, 0.0–3.0**; x-axis **Time (days), 0–80**.
Legend is by **T_g (°C): open circle 8, filled circle 15, open triangle 28** — i.e. the three
a_w 0.54 systems (PVP-LMW, PVP-K15, PVP-K30 respectively).

| series (T_g) | intercept at t = 0 **[D]** | last point (≈ 67 d) **[D]** | slope **[Z]** | Table 2 k **[M]** |
|---|---|---|---|---|
| **8 °C** (PVP-LMW) | ≈ 0.50 | ≈ 2.75 | ≈ **0.0336** | **0.034** ✔ |
| **15 °C** (PVP-K15) | ≈ 0.50 | ≈ 2.12 | ≈ **0.0242** | **0.024** ✔ |
| **28 °C** (PVP-K30) | ≈ 0.31 | ≈ 1.21 | ≈ **0.0134** | **0.014** ✔ |

**All three digitised slopes reproduce Table 2 to ≤ 4 %. Fig. 1 corroborates Table 2 and adds
no new parameters.** Six or seven time points per series; data run to ≈ 67 days.

**⚠ Two things Fig. 1 reveals that the tables hide [Z]:**
1. **The t = 0 intercept is NOT zero — it is 0.31–0.50 OD/g.** That is **9–19 % of the final
   absorbance already present at "time zero"**, i.e. browning that happened during lyophilisation
   and the one-week 25 °C equilibration. **The rate constants are therefore slopes on an already
   partly-browned system, and any model asked to reproduce them must not start from zero
   pigment.**
2. **The intercepts differ between systems (0.50 / 0.50 / 0.31)** — the PVP-K30 system starts
   *less* browned as well as browning *slower*. Some of the "glass transition effect" may have
   been established during the equilibration week, before t = 0.

---

## 8. THINGS THAT ARE CITED, NOT MEASURED — do not import numbers from these

- **"100-fold increase in mobility"** through Tg — **Roozen, Hemminga & Walstra (1991)**,
  *Carbohydr. Res.* **215**, 229–237, ESR of TEMPOL in glassy maltodextrin. **Cited only;
  no ESR was done in this paper.**
- The pseudo-zero-order form itself — Mizrahi et al. 1970; Warmbier et al. 1976a,b;
  Labuza & Saltmarch 1981; Labuza & Baisier 1992. **Assumed, not validated here.**
- Aspartame rearrangement is governed by a_w, **not** Tg — **Bell & Hageman (1994)**,
  *J. Agric. Food Chem.* **42**, 2398–2401. Used to argue that the *Amadori rearrangement*
  step of browning should be a_w-sensitive while the *diffusion/collision* step is Tg-sensitive
  (p. 595). **A hypothesis in this paper, with no measurement behind it.**
- Sucrose hydrolysis governed by internal pH, not state — **Buera, Chirife & Karel (1995)**,
  *Food Res. Int.* **28**, 359–365.
- The moisture-sorption-isotherm equivalence of the three PVPs — **Bell & Hageman (1995)**,
  *J. Food Qual.* **18**, 141–147. **The isotherms themselves are not reproduced in this paper**;
  only Table 1's moisture contents stand in for them.
- Tg values *"similar to those obtained previously"* — Karmas et al. 1992; Buera et al. 1992;
  Bell & Hageman 1995. **The Tg values in Table 2 are this paper's own DSC measurements.**
- 28 references total; the full list is on pp. 596–597 and contains no further numeric data.

---

## 9. USABILITY CAVEATS THAT APPLY TO EVERY NUMBER ABOVE

1. **One temperature only: 25 °C. There is no Eₐ, no Arrhenius pair, and no route to one from
   this paper.** Anything in the repo that wants a temperature dependence of browning must get
   it elsewhere; Bell gives only the 25 °C level and its a_w/Tg/concentration dependence.
2. **The units are OD/g/d, not mol L⁻¹ s⁻¹.** Conversion needs the extraction geometry
   (200 mg solid → 4 mL water → A₄₂₀ of the filtrate), a melanoidin extinction coefficient, and
   an assumption that all pigment is water-extractable and passes a 0.20 μm filter. **None of
   those is provided. Treat k as an instrument-referenced pigment rate.**
3. **Pseudo-zero order is assumed, not tested.** R² is 0.932–0.998 over ≤ 67 d, but zero-order
   fits look good over any short window of a sigmoidal curve.
4. **t = 0 is post-equilibration, and the t = 0 absorbance is already 0.31–0.50 OD/g** (§7).
5. **No replicate count.** *"Multiple 200 mg samples"* is the only statement. The ± values are
   **regression 95 % CLs on the slope**, which conflate measurement scatter with model error.
   **Relative CL ranges from 3 % to 32 %** (§4) — the a_w 0.33 glassy rows are the noisiest,
   which is exactly why "no significant difference" there is a weak result.
6. **The salts for a_w 0.44, 0.66 and 0.85 are never named**, so those three conditions cannot
   be reproduced exactly.
7. **"Room-temperature" pH 7 phosphate, in a solid.** The internal pH of a reduced-moisture
   phosphate-buffered solid is not 7; the author himself published on this (Bell & Labuza 1994,
   *J. Food Eng.* **22**, 291–312, cited on p. 592) and does not correct for it here.
8. **PVP is an inert *polar* matrix, not a food.** No protein, no lipid, no competing amine, no
   metal ions. The absolute k values are not transferable to a real matrix; the **ratios** are
   what this paper contributes.
9. **Tg is the DSC ONSET at 5 °C/min.** Midpoint Tg would be several degrees higher and
   scan-rate dependent. Any repo Tg model must match the onset convention.
10. **The a_w 0.96 "solution" row is a different physical system** (1 % PVP, not a solid) and its
    Tg was measured by modulated DSC, a different method from the other twelve. Keep it separate.
11. **PVP-K30 is systematically the slowest system at every a_w** (0.0054 at 0.33, 0.014 at 0.54,
    0.013 at 0.76). The paper attributes this to its higher Tg. **[Z] An alternative reading —
    that MW itself (viscosity, chain entanglement, or residual dialysis differences) matters
    independently of Tg — is not excluded by this design**, because MW and Tg are perfectly
    collinear in it. The deconfounding separates a_w from Tg; **it does not separate Tg from MW.**
    ⚠ This is the design's one real limitation and the paper does not acknowledge it.

---

## 10. VERDICT — what is usable

### NOW USABLE

| block | count | status |
|---|---|---|
| **Table 2 — 13 k with 95 % CL and R², at 25 °C** | 13 × 3 | **FULLY USABLE.** Every cell legible; three of them independently reproduced from Fig. 1 to ≤ 4 % (§7). |
| **Table 2 — 13 T_g values (DSC onset, 5 °C/min)** | 13 | **FULLY USABLE**; confirmed against Fig. 2's T−T_g abscissa. |
| **Table 3 — 4 k at non-1-molal, with CL and R²** | 4 × 3 | **USABLE**, with the a_w 0.56 → 0.54 typo and the 0.28-vs-0.29-molal discrepancy flagged. |
| **Table 1 — full composition, 12 systems × 6 quantities** | 72 | **FULLY USABLE and arithmetically self-consistent** (§3). This is what makes the design auditable. |
| **T − T_g for all 13 systems** | 13 | **[Z]**, trivially derived, matches Fig. 2. |
| **Apparent reaction order in reactant molality, n = 1.05 and 1.40** | 2 | **NEW — [Z], never computed by the author.** §5a. |
| **Deconfounded glass-transition effect size, 2.4×** | 1 | **NEW — [Z].** §6. The number the repo actually needs. |
| **The a_w plateau above 0.54 at fixed molality** | 7 points | **[M]**, §4d. Directly contradicts the textbook bell-shaped a_w curve and explains why. |

**Score against the brief: 17 of 17 rate constants, 13 of 13 T_g, and the complete
deconfounding design transcribed. Not one cell was unreadable. Nothing is missing from this
paper — it is short, complete, and fully extracted.**

### STILL MISSING / NOT OBTAINABLE FROM THIS PAPER

1. **Any temperature other than 25 °C** → no Eₐ, no Arrhenius pair, ever.
2. **Any molar rate.** OD/g/d cannot be converted without an extinction coefficient.
3. **Replicate n and any SD** — only regression CLs.
4. **The moisture sorption isotherms** — in Bell & Hageman (1995), not here.
5. **Any measurement of an intermediate** — no Amadori, no dicarbonyl, no fluorescence, no
   A₂₉₄. **A₄₂₀ pigment only.** The paper's central hypothesis (Tg governs diffusion/collision
   while a_w governs the Amadori rearrangement, p. 595) is therefore **untested by its own data**
   and rests entirely on the aspartame analogy.
6. **Separation of Tg from polymer MW** — structurally impossible in this design (§9.11).
7. **The salt identities for a_w 0.44 / 0.66 / 0.85.**

---

*Extraction performed 2026-08-28 (Wave K4c). Tables 1–3 and Figures 1–4 re-read off 300 dpi
rasters of PDF pages 3–5 (printed pp. 593–595); page 1 header re-read to settle the 1995/1996
year question. The OCR text layer was used for prose only. No cell in any table was
unreadable.*
