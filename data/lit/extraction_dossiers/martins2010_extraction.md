# Martins, Leussink, Rosing, Desclaux & Boucon 2010 — EXTRACTION
### *"Meat Flavor Generation in Complex Maillard Model Systems: Reactants Synergetic and/or Competition Effect?"*
### A 32-run D-optimal fractional factorial at 100 °C, pH 6 (0.2 M pyrophosphate): xylose, glucose, 5-ketogluconic acid, cysteine, glycine, glutamic acid at two levels, time at three. Six responses (MFT, FFT, norfuraneol, furaneol, HMF, furfural) measured, but **only the response ranges and the signs/strengths of significant regression effects are printed** — no per-run data. One bar chart (furaneol) and one prediction profiler (norfuraneol) carry readable numbers. Read 2026-10-27.

**Full read of `data/articles/martins2010.pdf`** (13 pages, book pp. 71–83; publisher PDF with a clean
text layer). Prose from `pdftotext`. **Tables I–VII were re-read from 120-dpi page renders** because the
text layer separates each table's header row from its sign row (e.g. Table II prints the six factor
names, then the six sign strings, with no alignment). **Figures 2 and 3 were read from 300/150-dpi
renders**; Fig. 1 (scheme), Fig. 4 (PCA loading plot) inspected, no numbers taken. No supplementary
material.

**Provenance:** **[M]** printed number · **[S]** printed sign/strength code · **[F]** read off a figure
by me · **[Q]** qualitative text · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | **S. I. F. S. Martins**\* (sara.martins@unilever.com), A. Leussink, E. A. E. Rosing, G. A. Desclaux, C. Boucon |
| Affiliation | Unilever R&D Vlaardingen, Olivier van Noortlaan 120, P.O. Box 114, 3130 AC Vlaardingen, The Netherlands |
| Title | "Meat Flavor Generation in Complex Maillard Model Systems: Reactants Synergetic and/or Competition Effect?" |
| Venue | Chapter 8 (pp. **71–83**) in *Controlling Maillard Pathways To Generate Flavors*; Mottram, D., et al., Eds.; **ACS Symposium Series 1042**; American Chemical Society: Washington, DC, **2010** |
| DOI | **10.1021/bk-2010-1042.ch008** (printed in the page margin); web publication date August 10, 2010 |
| PDF character | 7 tables (one ranges table, six effect-sign tables), 4 figures. **No raw or per-run data, no replicate count, no SDs.** |
| Relation to the repo | Same first author as the trunk's kinetic source (Martins & van Boekel 2003/2005, glucose–glycine; `species.py` "Glycine in the Martins system"), here at Unilever rather than Wageningen; **a different study** — a screening design, not a kinetic fit |

---

## 1. THE ONE-PARAGRAPH ANSWER

Six reactants at two levels and time at three (0–3 h), 100 °C, pH 6, 0.2 M pyrophosphate, 50 mL,
oil bath; 32 runs in 4 blocks. MFT and FFT by HS-SPME–GC–MS with deuterated internal standards and DTT
derivatisation; norfuraneol, furaneol, HMF and furfural by ethyl-acetate extraction and external
calibration. Responses were log10-transformed and fitted by OLS with main effects, time², and two-way
interactions; only significant effects are printed, as signs (+/−, one to three). **What the model can
use:** (i) response **ranges** across the design (MFT 0–348.4, FFT 0–352.4, norfuraneol 0–10 003.1,
HMF 0–1.6, furfural 0–8.9, furaneol 0–42.3 "ppm"; the thiol unit is probably ng/g, see Flag 1);
(ii) **directional claims**: cysteine lowers HMF (−−−) and furfural (−−); glycine raises norfuraneol
(++), furaneol (+), MFT (++) and FFT (++); cysteine × glycine is **synergistic** for both thiols (+);
xylose × KGA compete for norfuraneol (−); glucose has **no significant effect** on MFT, FFT,
norfuraneol or furfural; (iii) **Fig. 2 furaneol at 3 h**: glucose ≈ 3.8, xylose ≈ 1.2, glucose +
glycine ≈ 32, xylose + glycine ≈ 19.7 ppm [F] (reactants 0.2 M). No rate constant, no time course per
run, no level for a specific composition other than those four bars.

---

## 2. METHODS [M] (pp. 72–76)

| item | as printed |
|---|---|
| Medium | **pyrophosphate buffer solution, 50 mL, 0.2 M, pH 6**; "constant pH6"; how pH set/measured and pH after heating **not stated** [!] |
| Heating | **oil bath at 100 °C**; vessel not described (open/closed not stated) [!]; come-up time not stated |
| a_w | aqueous, not stated (≈ 0.99 [Z]) |
| Reactants | D-glucose, D-xylose, L-cysteine, L-glutamic acid, glycine (Sigma-Aldrich); 5-ketogluconic acid (KGA) "donated by a supplier" |
| Design | **D-optimal customised fractional factorial, 32 combinations, 4 blocks of 8**; six ingredients at 2 levels, time at 3 levels; supports main effects, 2-factor interactions, time² |
| Regression | OLS on **log10** responses (heteroscedasticity); **null values replaced by the lowest measured value**; mean of each experiment used; non-significant terms dropped and model refitted (JMP v8); α = 0.05 |
| Thiols (MFT, FFT) | 5 g sample + **2 % DTT** (1,4-dithioerythritol; the text also writes "DDT"/"DTT") + [²H3]-MFT (**12 ng/g**) and [²H2]-FFT (**21 ng/g**) (aromaLAB); vortex 1 min; **linearity 0.5–800 ng/g** in pyrophosphate + 2 % DTT; HS-SPME 100 µm PDMS, 2 min pre-heat + 30 min at **60 °C**, 250 rpm; Agilent 6890/5973, Rxi-5MS 30 m × 0.25 mm × 0.25 µm, 0 °C 2 min → 3 °C/min → 75 °C → 30 °C/min → 250 °C (15 min); SIM on molecular ion + qualifiers; autosampler tray 5 °C |
| Non-sulfur (norfuraneol, furaneol, HMF, furfural) | 2 mL sample + 1 mL saturated NaCl, extracted with 1 mL ethyl acetate (vortex 1 min, 2500 rpm); 1 µL into a CIS at 80 °C → 250 °C; splitless except norfuraneol (split 1:200); Wax 58 25 m × 0.15 mm; 60 °C 2 min → 10 °C/min → 240 °C (9 min); **external calibration** |
| Sensory | expert panel, 0–5 categories, PCA (Varimax) |
| Replicates | **not stated** ("the mean value for each experiment was used") [!] |

### Table I (p. 75), verbatim [M]

| | Name | Lower Level | Higher Level |
|---|---|---|---|
| Factors | D-Xylose (%) | 0.05 | 3 |
| | D-Glucose (%) | 0 | 3.6 |
| | 5-ketogluconic acid (%) | 0 | 3.0 |
| | L-Cysteine (%) | 0.05 | 0.75 |
| | L-Glycine (%) | 0 | 1.5 |
| | L-Glutamic acid (%) | 0 | 1.5 |
| | Time (h) | 0 | 3.0 |
| Responses | MFT (ppm) | 0 | 348.4 |
| | FFT (ppm) | 0 | 352.4 |
| | Norfuraneol (ppm) | 0 | 10003.1 |
| | 5-HMF (ppm) | 0 | 1.6 |
| | 2-Furfural (ppm) | 0 | 8.9 |
| | Furaneol (ppm) | 0 | 42.3 |

Footnote: MFT = 2-methyl-3-furanthiol; FFT = 2-furfurylthiol; norfuraneol = 4-hydroxy-5-methyl-3(2H)-furanone;
5-HMF = 5-hydroxymethylfurfural; furaneol = 4-hydroxy-2,5-dimethyl-3(2H)-furanone.

Molar equivalents [Z] (assuming % = g/100 mL in the 50 mL buffer):

| factor | low → high (mM) |
|---|---|
| xylose (150.13) | 3.3 → **200** |
| glucose (180.16) | 0 → **200** |
| KGA (194.14) | 0 → 155 |
| cysteine (121.16) | 4.1 → **61.9** |
| glycine (75.07) | 0 → **200** |
| glutamic acid (147.13) | 0 → 102 |

The high levels of xylose, glucose and glycine are each **20 mmol/100 mL**, matching the Fig. 2 caption.
The **middle time level is not printed**; Fig. 3's profiler centre line sits at **1.75 h**, which is not
the midpoint of 0–3 [!].

---

## 3. RESULTS

### 3.1 Effect tables, verbatim (from the page images)

Common footnote: "Main factors as significant & relative importance. (+) positive effect or interaction
(synergy between two components); (−) negative effect or interaction (competition between components).
The number of signs indicates the strength of the effect." R² and RMSE refer to the **log10** response.

**Table II. Main factors in 5-HMF formation** (p. 77)

| Xyl | Glc | Cys | Xyl\*Glu | Glc\*Gly | Cys\*Glu | R² | RMSE |
|---|---|---|---|---|---|---|---|
| +++ | ++ | −−− | ++ | − | −− | 0.86 | 0.11 |

**Table III. Main factors in Furaneol formation** (p. 77)

| Gly | Time | Xyl\*Glc | R² | RMSE |
|---|---|---|---|---|
| + | +++ | − | 0.82 | 0.78 |

**Table IV. Main factors in Norfuraneol formation** (p. 79)

| Xyl | KGA | Gly | Time | Xyl\*KGA | Cys\*Glu | R² | RMSE |
|---|---|---|---|---|---|---|---|
| ++ | ++ | ++ | ++ | − | − | 0.85 | 0.42 |

**Table V. Main factors in 2-Furfural formation** (p. 79)

| Xyl | KGA | Gly | Cys | Time | Glc\*Gly | R² | RMSE |
|---|---|---|---|---|---|---|---|
| +++ | ++ | − | −− | ++ | + | 0.79 | 0.43 |

**Table VI. Main factors in MFT formation** (p. 81)

| Cys | Gly | Time | Xyl\*KGA | Xyl\*Gly | Cys\*Gly | R² | RMSE |
|---|---|---|---|---|---|---|---|
| ++ | ++ | ++ | − | − | + | 0.82 | 0.52 |

**Table VII. Main factors in FFT formation** (p. 81)

| Xyl | Cys | Gly | Time | Xyl\*Gly | Cys\*Gly | R² | RMSE |
|---|---|---|---|---|---|---|---|
| ++ | ++ | ++ | − | + | + | 0.82 | 0.52 |

RMSE in log10 units [Z]: 0.11 → typical prediction error ×1.3 (HMF); 0.42–0.52 → ×2.6–3.3 (norfuraneol,
furfural, MFT, FFT); 0.78 → ×6 (furaneol).

### 3.2 Figure 2 (p. 78) — furaneol after 3 h at 100 °C, pH 6; "reactants concentration 20 mmol/100 mL" [F]

| arm | furaneol (ppm) [F, ±0.5] |
|---|---|
| Glc + NO Gly | 3.8 |
| Xyl + NO Gly | 1.2 |
| Glc + Gly | 32.1 |
| Xyl + Gly | 19.7 |

What the "NO Gly" arms contain besides the sugar (cysteine? glutamic acid?) is **not stated** [!]. The
32.1 ppm bar is below the Table I maximum of 42.3, so these are design runs or a subset.
[Z] Glc + Gly: 32.1 mg/L / 128.13 = **0.25 mM furaneol = 0.13 mol % of 200 mM glucose**;
Xyl + Gly: 0.154 mM = 0.077 mol % of xylose.

### 3.3 Figure 3 (p. 80) — prediction profiler for log10(norfuraneol) [F, ±0.1 log unit]

Upper row, xylose held at **0.05 %**; lower row, xylose held at **3 %**; other factors at the dashed
centre lines (KGA 1.5 %, glycine 0.75 %, time 1.75 h, cysteine 0.4 %, glutamic acid 0.75 %).

| setting | predicted log10 norfuraneol (ppm) at the centre | KGA 0 → 3 % |
|---|---|---|
| xylose 0.05 % | ≈ 2.3 (≈ 200 ppm) | ≈ 1.9 → 2.8 (rising) |
| xylose 3 % | ≈ 3.0 (≈ 1000 ppm) | ≈ 3.0 → 2.95 (flat) |

Other profiles, xylose 3 % row [F]: glycine 0 → 1.5 % ≈ 2.6 → 3.3; time 0 → 3 h ≈ 2.7 → 3.2; cysteine
and glutamic acid ≈ flat to slightly rising. These are **model predictions**, not measurements.

### 3.4 Text claims [Q] (pp. 77–83)

- Cysteine "inhibits sugar degradation via the 1,2-enolisation step" (negative on HMF and furfural);
  glycine promotes 2,3-enolisation (positive on furaneol and norfuraneol). Labelled "preliminary".
- In glucose mixtures with no KGA and xylose at 0.05 %, norfuraneol and furfural "were formed in very
  low amounts".
- "Regarding FFT formation, xylose had a clear positive effect whereas the impact of glucose was minimal.
  **The same was observed with MFT formation**" — but **Table VI does not list xylose** among MFT's
  significant main effects [!].
- Glutamic acid: no significant effect on MFT.
- Sensory PCA: first two components 20.4 % and 19.8 % of variance; "meaty" correlates with xylose and
  cysteine; without cysteine "caramellic notes prevailed". Authors call it qualitative.
- Conclusion: R² 0.8–0.9; "intended as a screening experiment … findings would need further
  validation … if we were to predict compound formation with statistical exactitude".

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Cannot anchor:** any rate constant, activation energy, pH term or absolute yield for a defined pot.
There are no per-run data; the printed ranges are the extremes of 32 different compositions and times
(including t = 0, hence the zeros). The regression signs are on log10 responses with strengths
expressed only as 1–3 signs; they are not elasticities. One temperature (100 °C), one pH (6).

**Can test (directional, within-study, hold-out only):** the conditions — 100 °C, pH 6, 0.2 M
pyrophosphate, aqueous, 0–3 h, cysteine 4–62 mM with pentose up to 200 mM — sit inside the sulfur
lane's natural envelope and close to the trunk's. Candidate directional claims:

| # | claim | lane / species | sign printed | strength |
|---|---|---|---|---|
| D1 | adding cysteine **lowers furfural** | sulfur (`FUR`, via `FUR → FFT`, `HMFAD`-type adducts) | Cys −− (Table V) | robust in sign; magnitude unknown |
| D2 | adding cysteine **lowers HMF** | trunk/furanic (`HMF`) in a xylose/glucose pot | Cys −−− (Table II) | strongest single effect printed |
| D3 | adding glycine **raises MFT and FFT** in a xylose + cysteine pot, and Cys × Gly is synergistic | sulfur (`MFT`, `FFT`; amine feeds ARP/DPO → `NF`) | Gly ++, Cys\*Gly + (VI, VII) | the paper's headline |
| D4 | adding glycine **raises norfuraneol** | sulfur (`NF`) | Gly ++ (IV) | — |
| D5 | glucose **does not** measurably change MFT, FFT, norfuraneol or furfural | sulfur | Glc absent from IV–VII | an absence of significance, weak |
| D6 | xylose and KGA **compete** for norfuraneol | — (KGA not in the model) | Xyl\*KGA − (IV; Fig. 3) | not testable |
| D7 | FFT **falls** with time over 0–3 h at 100 °C, MFT rises | sulfur (FFT sinks vs formation) | FFT Time −, MFT Time ++ | interesting: the lane must show net FFT loss within 3 h at 100 °C at some design points |

And one level test: **Fig. 2 furaneol** (0.2 M glucose ± 0.2 M glycine, 3 h, 100 °C, pH 6 pyrophosphate):
≈ 3.8 → 32 ppm with glycine (×8.4). The trunk's furanic channel (`parameters_furanic.py`, DMHF) could be
compared against the glucose + glycine bar, **but** the "NO Gly" arm's composition is unknown, the
amine in the trunk is glycine-generic, and the bars carry no SD — so only the order of magnitude
(≈ 0.1 mol % of glucose at 3 h) and the glycine ratio are defensible.

**Unit caveat that decides any level use (Flag 1):** the thiol calibration was linear **0.5–800 ng/g**
against internal standards spiked at 12 and 21 ng/g, yet Table I prints MFT/FFT maxima of ~350 "ppm".
348.4 ppm MFT would be 3.05 mM, ~5 mol % of the highest cysteine level at 100 °C — far above any
literature yield at 100 °C; 348.4 ng/g (ppb) sits inside the calibrated range. **I read the thiol maxima
as ng/g (ppb), not ppm.** The non-sulfur units (external calibration, ethyl-acetate extract) are
plausibly ppm; norfuraneol 10 003 ppm (≈ 88 mM, ≈ 25 mol % of xylose + KGA at their highs) is large but
the split 1:200 injection says the analysts expected it to be large.

**Under the observable-first rule:** nothing here can promote a prior. D1–D4 and D7 are honest
within-study directional hold-outs at one temperature/pH; they need the design's actual run compositions
(not printed) to be posed precisely, so as printed they are **qualitative sign checks**, not scored
claims.

---

## 5. FLAGS [!]

1. **Thiol units**: "ppm" in Table I versus a 0.5–800 ng/g calibration; probably ng/g.
2. **No per-run data, no replicate count, no SD.** Zeros replaced by the minimum before log-transform
   (biases effects of absent reactants).
3. **Time levels**: Table I says 0–3 h, three levels, middle not printed; Fig. 3 centre at 1.75 h.
4. **Table VI omits xylose** for MFT although the text says xylose's positive effect "was observed with
   MFT formation". Either the table dropped it or the text is loose; since the minimum xylose is 0.05 %
   (never zero), the pentose was present in every run.
5. Table VI/VII footnotes omit glutamic acid from the abbreviation list. Throughout, "Glu" is glutamic
   acid and "Glc" glucose — the opposite of the Cha 2019 convention (`cha2019_extraction.md`, where
   "Glu" is glucose).
6. Vessel (open/closed), come-up time, pH after heating: not stated.
7. Fig. 2's "NO Gly" arm composition not stated.
8. "DDT" in the text is a typo for DTT; the paper also names the reagent "1,4-dithioerythritol (DTT)",
   but 1,4-dithioerythritol is DTE and DTT is its threo isomer — either way a dithiol reductant. 2 % DTT in the sample also reduces any
   MFTD/FFTD back to the thiols, so **MFT and FFT here are total thiol (free + disulfide-bound)**, not free
   thiol — relevant if ever compared with a model's free `MFT` alone.
