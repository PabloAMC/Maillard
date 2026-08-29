# Gigl, Hofmann & Frank 2021 — NMR studies on odorant–melanoidin interactions in coffee beverages

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THE CORPUS'S ONLY DIRECT MEASUREMENT OF THE **TEMPERATURE STRUCTURE** OF A THIOL SINK — AND IT SAYS THE SINK IS TWO CHANNELS OF OPPOSITE SIGN, NOT ONE.

**Provenance codes:** **[M]** measured and printed · **[C]** cited from other work ·
**[F]** fitted by the authors · **[D]** derived/digitised by this wave and never printed ·
**[NEG]** verified negative · **⚠ [!]** integrity flag.

---

## §0. IDENTITY — ★ CONFIRMED, DOI VERIFIED ON THE PAGE

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/gigl2021.pdf`, **1,589,411 bytes** | `ls -la` |
| SHA-256 | `b7415b6e3cec8368373cc65a976c845687401006b883a23419f69902c33cccfc` | `shasum -a 256` |
| **title** | **"NMR-Based Studies on Odorant−Melanoidin Interactions in Coffee Beverages"** | p. 15334, raster |
| **authors** | **Michael Gigl, Thomas Hofmann, and Oliver Frank\*** (\* = corresponding) | p. 15334, raster |
| affiliation | **all three: Chair of Food Chemistry and Molecular Sensory Science, Technische Universität München, Freising D-85354, Germany** | p. 15343 |
| corresponding | **Oliver Frank**, `oliver.frank@tum.de`, Phone 49-8161/71-2910, Fax 49-8161/71-2949, ORCID `0000-0003-0437-5005`. Hofmann ORCID `0000-0003-4057-7165`. **Gigl has no ORCID printed.** | p. 15343 |
| **DOI** | ★ **`10.1021/acs.jafc.1c06163`** — **raster-verified glyph-by-glyph on the p. 15334 footer.** Matches the brief exactly. | p. 15334 footer, raster |
| journal | ***J. Agric. Food Chem.* 2021, 69, 15334−15344** (11 pages) — raster-verified in the "Cite This" box | p. 15334, raster |
| **dates** | **Received: October 1, 2021 · Revised: November 22, 2021 · Accepted: November 23, 2021 · Published: December 7, 2021** | p. 15334 footer, raster |
| licence | *"© 2021 The Authors. Published by American Chemical Society"* — open access | p. 15334 |
| producer | Arbortext Advanced Print Publisher 11.2.5208/W Library-x64 → Acrobat Distiller 8.1.0 → iText 4.2.0 | `pdfinfo` |
| **funding** | ★ **[NEG] NONE.** There is **no Acknowledgement section and no funding statement anywhere in the article.** Grepped the full text layer for `acknowledg`, `funding`, `grant`, `supported by` — zero hits outside a reference title. | full text |
| **competing interests** | *"The authors declare no competing financial interest."* | p. 15343 |
| keywords | coffee beverages, aroma binding, quantitative NMR spectroscopy, qHNMR, noncovalent π−π interactions, covalent aroma binding, melanoidins | p. 15334 |
| **SI** | ★ **EXISTS but is NOT ON DISK.** `ls data/articles/ \| grep -i gigl` returns only `gigl2021.pdf`. | filesystem |

**SI contents as declared, verbatim (p. 15343):**
> *"Concentrations of odorants used for the aroma model, qHNMR data of time dependency
> measurements, monitoring of the purification process of the HMW material >10 kDa, overview
> of interaction scenarios I-IV, and percent of free odorants in the substitution experiments
> (PDF)"*

Mapped to labels used in the body: **Table S1** = aroma-model concentrations · **Table S2** = the
qHNMR time-dependency table (the numeric source behind Figure 2) · **Figure S1** = ultrafiltration
purification monitoring · **Figure S2** = the four-scenario schematic · **Figure S3(A/B)** = percent
free odorant in the substitution experiments (FFT, N-acetylcysteine, MeFFT). See §12 for what
ordering it would buy.

---

## §0b. ⚠️ THE μ→m RASTER CHECK (Amendment 4) — **DISCHARGED. NO CORRUPTION FOUND ANYWHERE.**

**Binding requirement met.** Method: `pdftoppm -r 600 -png -f N -l N`, pages 1, 2, 3, 4, 6, 8, 9,
cropped and read as images, glyph by glyph. **Every unit that appears anywhere in this dossier was
read off a 600 dpi raster, not the text layer.**

| # | artefact | text-layer says | **raster (600 dpi) says** | verdict |
|---|---|---|---|---|
| 1 | p.15334 DOI footer | `10.1021/acs.jafc.1c06163` | **`https://doi.org/10.1021/acs.jafc.1c06163`** | ✅ **correct, brief confirmed** |
| 2 | p.15334 citation line | `2021, 69, 15334−15344` | **`J. Agric. Food Chem. 2021, 69, 15334–15344`** | ✅ |
| 3 | p.15335 coffee dose | `5.6 g of coffee/capsule` | **`5.6 g of coffee/capsule`** | ✅ |
| 4 | p.15335 water volume | `104 mL of water/capsule` | **`104 mL of water/capsule`** | ✅ **mL, not µL** |
| 5 | p.15335 brew strength | `5.4 g/100 mL` | **`(5.4 g/100 mL)`** | ✅ |
| 6 | p.15335 pooled volume | `n = 5, 520 mL` | **`(n = 5, 520 mL)`** | ✅ |
| 7 | p.15335 MWCO | `10 kDa` | **`10 kDa`** | ✅ |
| 8 | p.15335 flush volume | `8 L` … `1 L` | **`(8 L)` … `(1 L)`** | ✅ |
| 9 | p.15335 NMR tube | `5 × 178 mm` | **`5 × 178 mm`** | ✅ **mm is CORRECT here (tube length), not a dropped µ** |
| 10 | p.15336 qHNMR calibrant | `L-tyrosine (5.21 mmol/L)` | **`L-tyrosine (5.21 mmol/L)`** | ✅ **mmol, genuinely** |
| 11 | p.15336 stock range | `2−25 mmol/L` | **`2−25 mmol/L`** | ✅ |
| 12 | p.15336 sample aliquot | `540 μL` | ★ **`540 μL` — TRUE GREEK MU, verified** | ✅ **NOT corrupted** |
| 13 | p.15336 buffer aliquot | `60 μL` | ★ **`60 μL`** | ✅ **NOT corrupted** |
| 14 | p.15336 HMW aliquot | `150 μL, concentration × 4` | ★ **`150 μL, concentration × 4`** | ✅ **NOT corrupted** |
| 15 | p.15336 final volume | `600 μL` | ★ **`600 μL`** | ✅ **NOT corrupted** |
| 16 | p.15336 buffer recipe | `KH2PO4 (10.2 g), KOH (1.5 g), TMSP (50 mg), NaN3 (5 mg) in D2O (40 mL)` | **identical** | ✅ |
| 17 | p.15336 DCl | `4.0 mol/L in D2O` … `50 mL` | **identical** | ✅ |
| 18 | p.15336 sensory buffer | `0.1 mol/L, pH 5.5` | **`0.1 mol/L, pH 5.5`** | ✅ |
| 19 | p.15336 beaker | `diameter 45 mm, capacity 45 mL` | **identical** | ✅ |
| 20 | p.15336 booth temp | `22−25 °C` | **`22−25 °C`** | ✅ |
| 21 | **p.15337 Table 1 caption + all 8 cells + footnote a** | see §8 | **every digit identical** | ✅ **table raster-verified end-to-end** |
| 22 | p.15338 Fig. 1 caption | `500 MHz; H2O/D2O, 9/1, v/v; pH 5.5; 300 K` | **identical** | ✅ |
| 23 | p.15338 odorant spike | `5 mmol/L` | **`(5 mmol/L)`** | ✅ **mmol — see the §1 saturation flag** |
| 24 | **Fig. 2 legend title** | *(image, no text layer)* | ★ **`Percent of free odorant`, scale ticks `100 80 60 40 20 0`** | ✅ raster-only |
| 25 | **Fig. 2 x-axis** | *(image, no text layer)* | ★ **`time (min)`, ticks `0 5 30 60 120 240 720 1440 2880 5760`** | ✅ raster-only |
| 26 | **Fig. 2 class legend** | *(image)* | **`aromatic thiol / diketones / furanones / pyrazines / hydroxyphenols / terpene / aldehydes / aliphatic thiols`** | ✅ raster-only |
| 27 | **Fig. 4 y-axis** | *(image)* | ★ **`Percent of free odorant`, ticks `0 20 40 60 80 100`** | ✅ raster-only |
| 28 | **Fig. 4 x-axis** | *(image)* | ★ **`log₁₀ Incubation time (min)`, decade ticks `1 10 100 1000 10000`** | ✅ raster-only |
| 29 | Fig. 4 caption | `60 °C (-□-), 25 °C (-△-), and 6 °C (-◇-)` | **identical** | ✅ |
| 30 | p.15341 covalent % | `approximately 5% … maximum after 96 h with 11%` | **identical** | ✅ |
| 31 | **p.15342 the whole temperature paragraph** | `24.1% / 23.6% / 55.0% / 32.2% / 32.3% / 72 h / 96 h` | **every digit identical** | ✅ **the §6 deliverable is raster-verified** |
| 32 | p.15342 authentic brew | `47.8 mmol/L … 24.4 mmol/L … 51%` | ★ **`47.8 mmol/L (reference solution) to 24.4 mmol/L … 51%`** | ✅ **NOT a µ→m drop — but the value is chemically impossible, see ⚠[!]-3** |

### ⇒ VERDICT: **THIS PDF IS CLEAN. The Arbortext/ACS μ hazard did not fire here.**
Every `μL` in the methods is a true Greek mu in both layers. Every `mmol/L` is genuinely `mmol/L`.
**No unit in this paper needs a 1000× correction.** *(Contrast `kang2026_extraction.md` §0b, where
the same hazard on an RSC-typeset article was real.)*

**But three integrity flags survive the raster check — they are the authors' errors, not the PDF's.**
See ⚠[!]-1 … ⚠[!]-4 in §4c.

---

## §1. SYSTEM AND CONDITIONS `[M]`

### 1a. Coffee, and the melanoidin isolate

| item | value | anchor |
|---|---|---|
| brew method | single-dose capsule machine, **DèLonghi Nespresso Inissia EN 80.B**, commercial capsules | p. 15335 |
| dose | **5.6 g coffee per capsule**, percolated with **104 mL water per capsule** | p. 15335 |
| **brew strength** | ★ **5.4 g / 100 mL** | p. 15335 |
| pooling | **five** beverages pooled, **n = 5, 520 mL**, ice-bath cooled to RT, immediately ultrafiltered | p. 15335 |
| coffee used for the HMW isolate | ★ **Dallmayr Capsa, Lungo Mild Roast** (named commercial product) | p. 15335 |
| **isolation** | **cross-flow ultrafiltration**, Sartorius Stedim Biotech, **MWCO 10 kDa** membrane → HMW **>10 kDa** / LMW **<10 kDa** | p. 15335 |
| **purification** | flushed with **8 L water** in **1 L** steps; **¹H NMR after every 1 L step**; flushing continued until **no LMW resonances in 5–8 ppm** and the residual HMW resonances stopped changing in intensity | p. 15335 |
| **concentration** | **concentrated to ¼ of its original volume (i.e. ×4)**, stored **−20 °C** | p. 15335 |
| **dialysis** | ★ **[NEG] NONE.** No dialysis step; ultrafiltration only. | p. 15335 |
| **freeze-drying** | ★ **[NEG] NONE.** The HMW fraction is used as a concentrated aqueous solution, never dried. | p. 15335 |
| **yield** | ★ **[NEG] NOT REPORTED.** No mass, no mg, no % recovery of the HMW fraction anywhere. | — |
| **melanoidin conc. in mg/mL** | ★ **[NEG] NEVER STATED.** Loading is only ever expressed as *"the original coffee concentration"*. | throughout |

**[D] Melanoidin loading arithmetic — closes to exactly 1× brew.** The NMR sample is
buffer **60 μL** + HMW fraction **150 μL at ×4** + odorant stock, made to **600 μL**.
`150 × 4 / 600 = 1.00`. **⇒ the melanoidin loading in every NMR tube equals that of the original
5.4 g/100 mL brew, to three figures.** This is the only quantitative handle on the melanoidin dose
that the paper supplies, and it is a *relative* one. **Nothing in this paper converts it to mg/mL.**

### 1b. NMR system, odorants, and the interaction assay

| item | value | anchor |
|---|---|---|
| spectrometers | **Bruker AVANCE III 500 MHz, cryo-TCI probe → used at 300 K**; **Bruker AVANCE Neo 500 MHz, cryo-TCI probe → used at 279 and 333 K** | p. 15335 |
| ⚠ | **the two temperature extremes were run on a DIFFERENT INSTRUMENT from the RT series** | p. 15335 |
| tubes | **5 × 178 mm**, Bruker Z107374 USC | p. 15335 |
| pulse sequence | standard Bruker 1D water suppression, **`noesygppr1d`** | p. 15335 |
| tuning | probe manually tuned/matched to **50 Ω** with the sample in place; lock phase auto-optimised; shim **z1–z5, xyz, z1–z5**; **O1 adjusted per odorant**; **90° pulse width determined per sample** (`pulsecal`) | p. 15335–36 |
| **scans** | ★ **NS = 16**, **DS = 4** dummy scans | p. 15336 |
| data points | **64 K complex**, spectral width **10,273.97 Hz** | p. 15336 |
| **relaxation delay** | ★ **D1 = 20 s**, explicitly *"to guarantee complete re-establishment of the equilibrium z-magnetization of all excited nuclei and thus enabling quantitative measurement"* | p. 15336 |
| spinning | **none** | p. 15336 |
| chemical-shift ref | **TMSP, 0.0 ppm** | p. 15336 |
| processing | 0.3 Hz exponential LB, zero-filled, FT, `apk0.noe` auto-phase (manual 0th/1st follow-up if needed), `absn` auto-baseline, **manual integration** with SLOPE and BIAS | p. 15336 |
| **quantitation** | ★ **ERETIC 2 / PULCON**, calibrated with **L-tyrosine 5.21 mmol/L as EXTERNAL reference**, integrating its **7.10 ppm (m, 2H)** signal | p. 15336 |
| solvent | **H₂O_sample / D₂O_buffer, 9/1 v/v** — i.e. **90 % ordinary water; this is a real aqueous system, not a D₂O model** | p. 15336 |
| **buffer** | KH₂PO₄ 10.2 g + KOH 1.5 g + TMSP 50 mg + NaN₃ 5 mg in D₂O 40 mL, **pH adjusted to 5.5 with DCl 4.0 mol/L in D₂O**, made to 50 mL | p. 15336 |
| **pH** | ★ **5.5** throughout — chosen as *"a typical coffee pH value"*. Adjusted with DCl in D₂O; **the paper says "pH", not "pD", and applies no glass-electrode isotope correction** ⚠ | p. 15336, 15342 |
| **odorant spike** | ★ **5 mmol/L** of one odorant per tube; stocks weighed gravimetrically to **2–25 mmol/L** and qHNMR-verified before use | p. 15336, 15338 |
| controls | aqueous odorant solutions, **HMW replaced by water**, treated identically | p. 15336 |
| screening protocol | **60 min at RT**, then **20 s vigorous shaking of the sealed tube**, then acquisition | p. 15336 |
| **time-course protocol** | samples mixed and **in the magnet within 5 min**; spectra at time steps over **5 min – 96 h** | p. 15336 |
| **time grid [D, raster]** | ★ **0, 5, 30, 60, 120, 240, 720, 1440, 2880, 5760 min** (= 0, 5 min, 30 min, 1, 2, 4, 12, 24, 48, 96 h) — read off the Fig. 2 x-axis and independently confirmed by the Fig. 4 marker positions | Fig. 2 axis, Fig. 4 |
| Fig. 3 sub-grid | **t1 = 5 min, t2 = 12 h, t3 = 24 h, t4 = 48 h, t5 = 96 h** — a subset of the above ✔ | Fig. 3 caption |
| **temperature control** | probe set to **279 K or 333 K** and **equilibrated 30 min at that temperature before acquisition** | p. 15336 |
| **temperature calibration** | ★ **[NEG] NONE REPORTED.** No methanol/ethylene-glycol shift thermometry, no stated accuracy on the set temperature. | — |
| measured quantities | free odorant concentration, **signal shape, chemical shift, and fwhm**, always against the isomolar HMW-free reference | p. 15336 |
| statistics/clustering | R 3.5.1 + RStudio 1.0.143; `pheatmap` 1.0.8; **agglomerative hierarchical clustering on Euclidean distances** (Ward refs 37, 38) | p. 15336 |

### ⚠ 1c. THE SATURATION CAVEAT — read this before using any bound fraction `[D]`
The spike is **5 mmol/L FFT ≈ 571 mg/L**. Authentic coffee brew carries FFT at roughly
**1–10 μg/L** — i.e. **the assay runs the ligand ~10⁵ above its real concentration, against a
melanoidin loading held at exactly 1× authentic.** Two consequences, in opposite directions:
1. **The binding sites are being titrated to saturation.** Any Langmuir-like site population is far
   past half-occupancy, so **the measured bound fractions (24–76 %) are LOWER BOUNDS on what would
   be bound at authentic ppb ligand levels**, where sites are in vast excess.
2. But **site density per gram of melanoidin is never measured** (§9), so that bound is not
   quantifiable. **Use the bound fractions as a floor and a direction, never as an absolute yield.**

---

## §2. QUANTIFICATION BASIS — **REQUIRED SECTION**

> ### **The quantification is ABSOLUTE, not relative: ERETIC 2 / PULCON external calibration against L-tyrosine at 5.21 mmol/L, converting integral to mol/L, with every odorant stock independently qHNMR-quantified before use — so the concentrations printed in mmol/L are real concentrations and the percentages are true bound fractions, not merely integral ratios to t = 0.**

**What that licenses.** ✅ **Absolute bound fractions are usable as such.** The statement "32.3 % of
FFT is free at 60 min" means 32.3 % of a known 5 mmol/L, not 32.3 % of an unknown t=0 integral.
✅ Concentrations may be differenced across odorants and across temperatures. ✅ Mass balance may be
attempted (and the authors do exactly that: free + dimer + unassigned = 100 %, which is how the
"11 % covalently incorporated" number is obtained). ⚠ It does **not** license absolute *rate*
constants, because no rate is fitted (§9).

**Detection limit / integral precision.** ★ **[NEG] no LOD, no LOQ, no integral SD is printed.**
The only handle is a single sentence (p. 15340): 3-methylbutanal at 99.4 % free is *"in the error
range of the analysis"*, citing **ref. 35 = Frank, Kreissl, Daschner & Hofmann, *JAFC* 2014, 62,
2506–2515** — i.e. **the precision claim is [C], outsourced to a methods paper not in this corpus.**
**⇒ The most that can be inferred is that the error range is ≥ 0.6 % relative.** For fitting,
assign uncertainty from the method, not from the paper.

**Replicates.** ★ **[NEG] NO REPLICATE COUNT ANYWHERE FOR THE NMR.** No `n`, no error bars in
Fig. 2 or Fig. 4, no SD on any qHNMR number in the article. Only the *sensory* block reports
replication (**two independent sessions**, SDs in Table 1). **⇒ Every NMR value in this dossier is
a single unreplicated measurement as far as the published record shows.** This is the single
largest quality limitation of the paper.

---

## §3. ★ THE FOUR INTERACTION SCENARIOS — THE MOST REUSABLE ARTEFACT IN THE PAPER

### 3a. The authors' classification, verbatim (p. 15341)

> *"Aroma-active compounds, such as furanones lacking reactive side groups and without a
> delocalized π-electron system, showed no interactions with coffee melanoidins (**I**). Odorants
> possessing reactive side groups but lacking aromatic systems, such as aldehydes and aliphatic
> thiols, exhibit most likely only covalent interactions with the coffee HMW fraction by forming
> covalent bonds with the reactive species incorporated in the HMW material (**II**). Odorants
> containing an aromatic system but no reactive side groups, such as pyrazines, are most likely
> attached by noncovalent π−π interactions to the melanoidins (**III**). Coffee aroma compounds
> possessing delocalized π-electrons, in combination with reactive side groups, show an interaction
> scenario composed of an initial drop in the concentration, followed by a steady decrease of the
> free odorant upon incubation with melanoidins (**IV**). This scenario could be exclusively
> observed for FFT and can be interpreted as a combination of two previously described scenarios
> (II and III), thus suggesting noncovalent and covalent interactions of FFT with the HMW fraction
> from coffee (Figure S2)."*

And the diagnostic rule, verbatim (p. 15339):

> *"compounds exhibiting noncovalent π−π interactions studied by NMR spectroscopy often showed
> line broadening and a shift of resonance lines toward higher frequencies."*

### 3b. **THE COMPLETE PER-ODORANT ASSIGNMENT TABLE** `[M]` (class) + `[D]` (digitised evidence)

All 20 odorants in the Fig. 2 heat map, plus the two substitution controls. "Class colour" is the
raster-read annotation strip; "cluster" is the paper's own hierarchical-clustering block. Free-%
columns are this wave's Fig. 2 digitisation (§5a).

| # | odorant | class colour (Fig. 2) | **cluster** | **SCENARIO** | free % @ 60 min | free % @ 96 h | evidence used |
|---|---|---|---|---|---|---|---|
| 1 | **2-furfurylthiol (FFT)** | aromatic thiol | **IV** | ★ **covalent AND π−π** | **32** | **2 → 0** | instantaneous drop (π) **+** steady further decline (covalent); Δδ **+0.56 Hz**; line broadening **+0.93 Hz**; **near-total loss of multiplicity**; **new signals = difurfuryl disulfide appear from 60 min** |
| 2 | 2,3-butanedione | diketones | **I** | no interaction | 96 | 96 | no change in concentration, shape, δ or fwhm |
| 3 | 2,3-pentanedione | diketones | **I** | no interaction | 96 | 95 | as above |
| 4 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone (**HDMF**) | furanones | **I** | ★ **no interaction** | 99 | **99** (text: *">96 % even after 96 h"*) | *"identical chemical shift, signal shape, and fwhm compared to the aqueous reference"* |
| 5 | 3-hydroxy-4,5-dimethylfuran-2(5H)-one (sotolon) | furanones | **I** | no interaction | 99 | 99 | as above |
| 6 | 3-isobutyl-2-methoxypyrazine | pyrazines | **III** | π−π only | 46 | 46 | one-shot drop then **flat**; no further loss to 96 h |
| 7 | **2,3-diethyl-5-methylpyrazine** | pyrazines | **III** | π−π only | 58 | 57 | drop **4.85 → 2.83 mmol/L = 58.5 %**, then flat; **fwhm 2.02 → 4.16 Hz (Δ +2.14)**; **Δδ +1.42 Hz** |
| 8 | 2-methoxyphenol (guaiacol) | hydroxyphenols | **III** | π−π only | 61 | 57 | drop then flat; **Δδ +1.52 Hz**, **Δfwhm +1.55 Hz** on H−C(3) |
| 9 | 4-ethyl-2-methoxyphenol | hydroxyphenols | **III** | π−π only | 61 | 61 | drop then flat |
| 10 | 4-hydroxy-3-methoxybenzaldehyde (vanillin) | hydroxyphenols | **III** | π−π only | 68 | 57 | drop then slight drift |
| 11 | (E)-β-damascenone | terpene | **III** | π−π only | 69 | 67 | drop then flat |
| 12 | 3-methylbutanal | aldehydes | **II** | covalent only | 95 | 60 | **no** initial drop (99.4 % at 60 min, *"in the error range"*); **steady** decline to 96 h |
| 13 | methylpropanal | aldehydes | **II** | covalent only | 93 | 56 | as above |
| 14 | 2-methylbutanal | aldehydes | **II** | covalent only | 93 | 58 | as above |
| 15 | acetaldehyde ᵃ | aldehydes | **II** | covalent only | 95 | 48 | as above — ⚠ **footnote a: "acetaldehyde concentrations were determined semiquantitatively"** |
| 16 | propanal | aldehydes | **II** | covalent only | 90 | 53 | as above |
| 17 | methional (3-(methylthio)propanal) | aldehydes | **II** | covalent only | 93 | 44 | as above |
| 18 | 3-mercapto-3-methylbutyl formate | aliphatic thiols | **II** | covalent only | 87 | **9** (text: **4.89 → 0.3 mmol/L = 6.1 %**) | no initial drop; **the steepest loss in the whole map**; attributed to *"degradation, dimerization, and … covalent adducts"* |
| 19 | 3-mercapto-3-methylbutanol | aliphatic thiols | **II** | covalent only | 93 | 19 | as above |
| 20 | 3-mercapto-2-butanone | aliphatic thiols | **II** | covalent only | 95 | 22 | as above |
| C1 | **N-acetyl cysteine** *(substitution control — free thiol, NO aromatic ring)* | — | — | ★ **II** | *(Fig. S3, not on disk)* | — | steady decline over 96 h with **no** initial drop; **new disulfide signals appear**; **8 % unaccounted → covalently incorporated** |
| C2 | **furfurylmethylsulfide (MeFFT)** *(substitution control — aromatic ring, thiol METHYLATED)* | — | — | ★ **III** | *(Fig. S3)* | — | initial drop **like FFT**, then **constant**; **no new signals ever**; **Δδ +0.51 Hz** (vs FFT +0.56 Hz) and loss of multiplicity |

### 3c. ★ WHY C1 AND C2 ARE THE LOAD-BEARING RESULT
The two substitutions **dissect FFT's scenario IV into its two halves by structure, at fixed
melanoidin and fixed conditions**:
- **Keep the thiol, delete the aromatic ring (N-acetylcysteine) ⇒ scenario II.** The covalent
  channel survives; the instantaneous drop disappears.
- **Keep the aromatic ring, block the thiol by methylation (MeFFT) ⇒ scenario III.** The
  instantaneous drop survives (Δδ +0.51 vs +0.56 Hz — **within 9 % of FFT's**); the covalent channel
  disappears entirely (no new signals ever).

**⇒ The instantaneous drop is a property of the FURAN RING (π-stacking). The slow decline is a
property of the FREE −SH (covalent). They are separable, independently switchable, and — as §6
shows — they respond to temperature with OPPOSITE SIGN.** This is the paper's real contribution
and the whole basis of §7.

---

## §4. EVERY TABLE, EVERY PRINTED NUMBER

### 4a. **Table 1** (p. 15337) — the only numbered table in the article. Raster-verified cell by cell.

**Caption, verbatim:** *"Table 1. Aroma Profiles of a Coffee Aroma Reconstitution Model with and
without Added Melanoidins"*
**Column header:** `attribute | intensityᵃ [aroma reconstitution model without added melanoidins | aroma reconstitution model with added melanoidins]`

| attribute | **without added melanoidins** | **with added melanoidins** | **Δ [D]** | **relative Δ [D]** |
|---|---|---|---|---|
| **Roasty/sulfury** | **2.2 (0.2)** | **1.6 (0.3)** | **−0.6** | **−27.3 %** |
| **Sweetish/caramel-like** | **1.5 (0.5)** | **1.9 (0.4)** | **+0.4** | **+26.7 %** |
| **Earthy** | **1.6 (0.4)** | **1.3 (0.5)** | **−0.3** | **−18.8 %** |
| **Smoky** | **1.7 (0.6)** | **1.2 (0.3)** | **−0.5** | **−29.4 %** |

**Footnote a, verbatim:** *"Intensities of the attributes were rated on a scale from 0 (not
detectable) to 3 (strongly detectable). The results are given as the means (standard deviation in
parenthesis)."*
⚠ **[NEG] no per-attribute significance test, no p-value, no n in the table.** See §8.

### 4b. Every other printed number in the paper, by page `[M]`

| # | quantity | value | condition | anchor |
|---|---|---|---|---|
| 1 | 3-methylbutanal, free | **99.4 %** | 60 min, RT, +HMW | p. 15340 |
| 2 | 3-methylbutanal, free | **5.04 → 3.17 mmol/L**, stated as **58.2 %** | 96 h, RT | p. 15340 ⚠[!]-1 |
| 3 | 3-mercapto-3-methylbutyl formate | **4.89 → 0.3 mmol/L** | 96 h, RT | p. 15340 |
| 4 | 2,3-diethyl-5-methylpyrazine | **4.85 → 2.83 mmol/L = 58.5 % free** | immediately after incubation, RT | p. 15340 |
| 5 | HDMF, free | **> 96 %** | 96 h, RT | p. 15340 |
| 6 | **FFT, free** | ★ **32.3 %** | initial decay / 60 min, RT (300 K) | p. 15340, 15342 |
| 6b | **FFT, free** | ★ **32.2 %** *(same quantity, restated)* | 60 min, RT | p. 15342 ⚠[!]-2 |
| 7 | **FFT, free** | ★ **0 %** ("no free FFT could be observed") | 96 h, RT | p. 15340 |
| 8 | 2,3-diethyl-5-methylpyrazine, fwhm | **2.02 Hz → 4.16 Hz**, broadening **+2.14 Hz**, shift **+1.42 Hz** to higher frequency, H−C(6) | 60 min, RT, +HMW | p. 15339 |
| 9 | **FFT H−C(3)** | broadening **+0.93 Hz**, shift **+0.56 Hz**, *"almost complete loss of multiplicity"* | 60 min, RT, +HMW | p. 15339 |
| 10 | 2-methoxyphenol H−C(3) | shift **+1.52 Hz**, Δfwhm **+1.55 Hz** | 60 min, RT, +HMW | p. 15339 |
| 11 | **MeFFT** | shift **+0.51 Hz**, loss of multiplicity | RT, +HMW | p. 15341 |
| 12 | **FFT covalently incorporated into melanoidin** | ★ **≈ 5 %** | **48 h**, RT | p. 15341 |
| 13 | **FFT covalently incorporated into melanoidin** | ★ **11 %** (maximum) | **96 h**, RT | p. 15341 |
| 14 | **N-acetylcysteine covalently incorporated** | ★ **8 %** | 96 h, RT | p. 15341 |
| 15 | **FFT, free — 279 K** | ★ **24.7 %** (initial drop) | 6 °C | p. 15342 |
| 16 | **FFT, free — 279 K** | ★ **24.1 %** | 6 °C, **60 min** | p. 15342 |
| 17 | **FFT, free — 279 K** | ★ **23.6 %** | 6 °C, **48 h** | p. 15342 |
| 18 | **FFT, free — 333 K** | ★ **55.0 %** | 60 °C, **60 min** | p. 15342 |
| 19 | **FFT — time to zero, 333 K** | ★ **72 h** | 60 °C | p. 15342 ⚠[!]-4 |
| 20 | **FFT — time to zero, 300 K** | ★ **96 h** | RT | p. 15342 |
| 21 | 2,3-diethyl-5-methylpyrazine, free — 279 K | ★ **51.8 %** | 6 °C | p. 15342 |
| 22 | 2,3-diethyl-5-methylpyrazine, free — 300 K | ★ **58.5 %** | RT | p. 15342 |
| 23 | 2,3-diethyl-5-methylpyrazine, free — 333 K | ★ **70.0 %** | 60 °C | p. 15342 |
| 24 | 3-methylbutanal, free — 333 K | ★ **49.0 %** | 60 °C, 96 h | p. 15341 |
| 25 | 3-methylbutanal, free — 300 K | ★ **57.0 %** *(same quantity as #2)* | RT, 96 h | p. 15341 ⚠[!]-1 |
| 26 | 3-methylbutanal — 279 K | *"remained constant within the time course of 96 h"* | 6 °C | p. 15341 |
| 27 | HDMF, free — 279 K | **> 95 %** | 6 °C | p. 15341 |
| 28 | HDMF — 333 K | *"no influence"* | 60 °C | p. 15341 |
| 29 | pyrazine in **authentic brew** | **47.8 → 24.4 mmol/L = 51 % free** | 60 min, freshly percolated coffee | p. 15342 ⚠[!]-3 |
| 30 | 3-AFC significance | **α = 0.05**, assessors *"able to correctly identify the sample with added melanoidins"* | — | p. 15337 |

**⇒ [NEG] There is no K, no K_d, no K_a, no binding constant and no stoichiometry printed
anywhere in this paper.** Item 8–11 are the only chemical-shift data; they are Δ values in Hz at
500 MHz (**1 Hz = 0.002 ppm [D]**), not tabulated shifts.

### 4c. ⚠ **THE FOUR INTEGRITY FLAGS** — all survive the raster check, i.e. they are the authors' errors

**⚠[!]-1 — 3-methylbutanal at 96 h/RT is given FOUR mutually inconsistent values.**

| source | value | basis |
|---|---|---|
| **p. 15340**, from the printed concentrations | **5.04 → 3.17 mmol/L** ⇒ **62.9 %** `[D]` | arithmetic on the authors' own numbers |
| **p. 15340**, stated | **58.2 %** | printed |
| **p. 15341**, stated | **57.0 %** | printed |
| **Fig. 2**, digitised | **60** `[D]` | §5a |
| **Fig. 4C**, digitised | **60.6** `[D]` | §5b |

`3.17 / 5.04 = 0.6290`, **not 0.582.** The two printed percentages (58.2, 57.0) also disagree with
each other. **The two independent figure digitisations agree with each other (60.0 vs 60.6) and sit
between the printed values.** ⇒ **Use ≈ 60 % for 3-methylbutanal at 96 h/RT; do not use 58.2 or 57.0
as if they were the same measurement, and do not trust the printed concentration pair.**

**⚠[!]-2 — FFT free at 60 min/RT is printed as 32.3 % on p. 15340 and p. 15342, and as 32.2 % five
lines later on p. 15342.** Trivial in magnitude (0.1 pp) but it shows the numbers were retyped
rather than transcluded. **Use 32.3 %** (the majority value, and the one the Fig. 4A digitisation
reproduces).

**⚠[!]-3 — the authentic-brew pyrazine concentrations are chemically impossible as absolute values.**
Raster-verified as **`47.8 mmol/L`** → **`24.4 mmol/L`**. But (i) every model spike in this paper is
**5 mmol/L**, and the measured model initials are 5.04 / 4.89 / 4.85 mmol/L — so 47.8 is **10×** the
protocol; (ii) 47.8 mmol/L of 2,3-diethyl-5-methylpyrazine is **≈ 7.9 g/L**, far above what will
dissolve in a coffee brew. **Almost certainly a decimal slip for 4.78 → 2.44 mmol/L.**
**The RATIO is unaffected and internally exact: `24.4/47.8 = 51.05 % ≈ 51 %` ✔.**
⇒ **RATIO-ONLY. Never register 47.8 mmol/L as a concentration.**

**⚠[!]-4 — the text says free FFT reaches zero at 72 h at 60 °C; Figure 4A shows it reaching zero at
48 h.** The Fig. 4 time grid is **5, 30, 60, 120, 240, 720, 1440, 2880, 5760 min** — verified by
marker-position detection in all four panels — and **contains no 4320-min (72 h) point.** The 60 °C
square at **2880 min sits clipped on the y = 0 axis** (§5b). ⇒ Either an unplotted 72-h measurement
exists, or the text is loose. **For modelling, use the figure: extinction at ≤ 48 h at 333 K, at
96 h at 300 K, never at 279 K.** The choice matters: it is the difference between a 2.0× and a
1.33× ratio, and hence between 17.4 and 7.2 kJ/mol on the time-to-extinction route (§6d).

**⚠[!]-5 (nomenclature)** — p. 15340 calls the new species **"difurfuryldithiol"**; p. 15341 and the
Chemicals list call it **"difurfuryldisulfide"** (the commercial reference compound). The dimer of a
thiol is the **disulfide**. **p. 15340 is wrong; the species is difurfuryl disulfide.**

---

## §5. EVERY FIGURE

Figures 1 and 3 are stacked NMR spectral excerpts with **no numeric axis that can be digitised into
a parameter** — they are qualitative evidence only, transcribed in §5c/§5d. **Figures 2 and 4 are
the quantitative artefacts and both were digitised.**

### 5a. ★ **FIGURE 2 — full heat-map digitisation** `[D]`

**Caption, verbatim (p. 15339):** *"Figure 2. Heat map and hierarchical cluster analysis showing the
qHNMR data of time dependency measurements. Samples were incubated at RT for 5 min−96 h after the
addition of the HMW fraction (>10 kDa) in the original coffee concentration, and the percent of free
odorant was determined in relation to the initial concentration (mmol/L) before melanoidin addition
vs after melanoidin addition at each respective incubation time; ᵃacetaldehyde concentrations were
determined semiquantitatively."*

**Method.** `pdftoppm -r 600`, page 6. Cell grid located by detecting the map frame and the three
cluster-block gaps (20 rows × 10 columns, cell pitch 43.8 px displayed). The colour bar
(`Percent of free odorant`, 100→0) was sampled at 0.5-px intervals between its raster-detected end
points and used as a nearest-neighbour LUT in RGB. Each cell value is the mean RGB of a 13×13 px
patch at the cell centre.
**Calibration check:** the **t = 0 column reads 99 for all twenty rows** — by construction it must be
100. ⇒ **systematic bias −1 pp, no scatter.** Values below are **raw** (uncorrected).
**Honest error bar: ± 4 pp** (±2 pp typical; ±5 pp below 20 %, where the colour ramp compresses
toward white).

| odorant | **0** | **5** | **30** | **60** | **120** | **240** | **720** | **1440** | **2880** | **5760** min |
|---|---|---|---|---|---|---|---|---|---|---|
| **2-furfurylthiol (FFT)** | 99 | **34** | **34** | **32** | **34** | **31** | **28** | **24** | **15** | **2** |
| 2,3-butanedione | 99 | 98 | 96 | 96 | 96 | 96 | 96 | 96 | 96 | 96 |
| 2,3-pentanedione | 99 | 99 | 95 | 96 | 96 | 96 | 95 | 95 | 95 | 95 |
| **4-hydroxy-2,5-dimethyl-3(2H)-furanone** | 99 | 98 | 99 | 99 | 99 | 99 | 99 | 99 | 99 | **99** |
| 3-hydroxy-4,5-dimethyl-2(5H)-furanone | 99 | 98 | 99 | 99 | 98 | 99 | 97 | 98 | 98 | 99 |
| 3-isobutyl-2-methoxypyrazine | 99 | 46 | 46 | 46 | 46 | 45 | 47 | 46 | 46 | 46 |
| **2,3-diethyl-5-methylpyrazine** | 99 | 59 | 58 | 58 | 57 | 58 | 58 | 59 | 57 | 57 |
| 2-methoxyphenol | 99 | 61 | 62 | 61 | 61 | 61 | 60 | 61 | 60 | 57 |
| 4-ethyl-2-methoxyphenol | 99 | 60 | 61 | 61 | 62 | 61 | 60 | 62 | 62 | 61 |
| 4-hydroxy-3-methoxybenzaldehyde | 99 | 69 | 68 | 68 | 67 | 66 | 64 | 64 | 62 | 57 |
| (E)-β-damascenone | 99 | 70 | 70 | 69 | 68 | 68 | 67 | 67 | 68 | 67 |
| **3-methylbutanal** | 99 | 99 | 95 | 95 | 92 | 85 | 74 | 69 | 64 | **60** |
| methylpropanal | 99 | 99 | 95 | 93 | 90 | 82 | 73 | 66 | 60 | 56 |
| 2-methylbutanal | 99 | 99 | 94 | 93 | 90 | 80 | 70 | 65 | 60 | 58 |
| acetaldehyde ᵃ | 99 | 97 | 95 | 95 | 86 | 62 | 56 | 53 | 52 | 48 |
| propanal | 99 | 99 | 94 | 90 | 86 | 77 | 66 | 59 | 56 | 53 |
| methional | 99 | 98 | 95 | 93 | 90 | 78 | 66 | 60 | 52 | 44 |
| **3-mercapto-3-methylbutyl formate** | 99 | 98 | 90 | 87 | 80 | 61 | 44 | 31 | 22 | **9** |
| 3-mercapto-3-methylbutanol | 99 | 97 | 95 | 93 | 87 | 70 | 50 | 40 | 31 | 19 |
| 3-mercapto-2-butanone | 99 | 97 | 94 | 95 | 89 | 75 | 60 | 50 | 40 | 22 |

**Validation against the printed text (four independent anchors):**

| quantity | printed | **digitised** | Δ |
|---|---|---|---|
| FFT free, 60 min, RT | **32.3 %** | **32** | **0.3** ✔ |
| 2,3-diethyl-5-methylpyrazine free, RT | **58.5 %** | **58–59** | **≤ 0.5** ✔ |
| HDMF free, 96 h | **> 96 %** | **99** | ✔ consistent |
| FFT free, 96 h | **0 %** | **2** | **2** ✔ (inside the ±5 pp low-end band) |
| 3-mercapto-3-methylbutyl formate, 96 h | **6.1 %** (4.89→0.3) | **9** | 2.9 ✔ |
| 3-methylbutanal, 96 h | 57.0 / 58.2 / **62.9** | **60** | see ⚠[!]-1 |

**⇒ Figure 2 is digitised to ±4 pp and validated on five anchors. USE.**

### 5b. ★★ **FIGURE 4 — the temperature-dependency digitisation. THE DELIVERABLE.** `[D]`

**Caption, raster-verified verbatim (p. 15342):** *"Figure 4. Summary of temperature dependency
measurements using four odorants (A: FFT, B: 2,3-diethyl-5-methylpyrazine, C: 3-methylbutanal, and
D: 4-hydroxy-2,5-dimethyl-3(2H)-furanone) for incubation (5 min−96 h) with coffee melanoidins
(>10 kDa) at 60 °C (-□-), 25 °C (-△-), and 6 °C (-◇-), respectively."*

⚠ **The caption says "25 °C"; the Methods say the RT series was run at 300 K = 26.85 °C, and the
low/high series at 279 K = 5.85 °C and 333 K = 59.85 °C.** This dossier uses the **Kelvin** values
throughout, because they are the ones the instrument was set to. **279 / 300 / 333 K.**

**Method.** `pdftoppm -r 600`, page 9. Axes calibrated from the raster: y = 0 at row 1273 and the
five major ticks at 1156 / 1041.5 / 924.5 / 807.5 / 691 ⇒ **5.821 px per percentage point** (linear
to <0.1 pp). x: decade ticks at 925.5 / 1253 / 1583 / 1910.5 / 2240.5 ⇒
**x = 925.5 + 328.75·log₁₀(t/min)**, right frame = 10 000 min ✔.
Markers were located by **flood-fill hole detection**: each hollow marker's white interior is a
connected component not touching the panel border; its bounding-box centre is the datum. Clean
markers are square 27×27 px (fill 1.00), triangle 29×29 (fill 0.52), diamond 30×30 (fill 0.52).
**Triangles carry a +0.7 pp bounding-box-centre correction**, derived from and then independently
confirmed by two printed anchors (FFT 60 min RT: 31.6+0.7 = 32.3 vs printed **32.3** ✔; pyrazine RT:
58.0+0.7 = 58.7 vs printed **58.5** ✔). **All triangle values below are corrected.**
**Honest error bar: ± 1.5 pp** for isolated markers, **± 2.5 pp** where markers overlap (flagged °).

#### PANEL A — **2-furfurylthiol (FFT)** ★

| t (min) | t | **333 K (60 °C) □** | **300 K (RT) △** | **279 K (6 °C) ◇** |
|---|---|---|---|---|
| 0 | — | 100 | 100 | 100 |
| 5 | 5 min | **54.9** *(printed 55.0 ✔)* | **32.3** *(printed 32.3 ✔)* | **24.5** *(printed 24.7 ✔)* |
| 30 | 30 min | 52.6 | 32.3 | 24.5 |
| 60 | **1 h** | **52.6** | **32.3** *(printed)* | **23.5** *(printed 24.1)* |
| 120 | 2 h | 51.5 | 32.3 | 23.1 |
| 240 | 4 h | 46.4 | 31.4 | 24.5 |
| 720 | 12 h | 38.5 | 27.2 ° | 24.5 ° |
| 1440 | 24 h | **16.0** | 23.4 ° | 25.8 ° |
| 2880 | **48 h** | ★ **0** *(marker clipped on the axis)* | **15.0** | **23.5** *(printed 23.6 ✔)* |
| 5760 | **96 h** | **0** | ★ **0** *(printed: "no more free FFT")* | ★ **22.6** |

#### PANEL B — 2,3-diethyl-5-methylpyrazine (scenario III control)

| t (min) | **333 K □** | **300 K △** | **279 K ◇** |
|---|---|---|---|
| 5 | 70.1 | 58.7 | ~53 ° |
| 30 | 72.0 | 58.3 | — ° |
| 60 | 69.1 | 58.2 | ~53 ° |
| 120 | 69.1 | 58.2 | ~53 ° |
| 240 | 70.1 | 58.3 | ~53 ° |
| 720 | 71.0 | 58.8 | ~53 ° |
| 1440 | 72.0 | 58.7 | ~53 ° |
| 2880 | 69.1 | 57.3 | ~53 ° |
| 5760 | 68.2 | 57.3 | ~53 ° |
| **printed** | ★ **70.0** | ★ **58.5** | ★ **51.8** |

⚠ The 279 K diamonds in panel B lie directly under the triangles and are clipped in every frame;
their holes are partial (area 144–408 px vs 450–480 for a clean diamond), which biases them high by
~1.3 pp. **Use the printed 51.8 %, not the digitisation, for the 279 K pyrazine.**
**Every trace in panel B is FLAT in time at all three temperatures — π−π stacking equilibrates in
under 5 min and does not evolve.**

#### PANEL C — 3-methylbutanal (scenario II control)

| t (min) | **333 K □** | **300 K △** | **279 K ◇** |
|---|---|---|---|
| 5 | ~98 ° | 100.3 | ~98 ° |
| 30 | ~95 ° | 95.7 | ~97 ° |
| 60 | ~93 ° | 95.7 | ~95 ° |
| 120 | 88.3 ° | 92.9 | ~92 ° |
| 240 | 80.7 ° | 85.4 | 92.1 |
| 720 | 69.1 ° | 74.7 | 90.6 |
| 1440 | **61.2** | 69.0 | 88.3 |
| 2880 | **55.1** | 64.3 | 84.6 |
| 5760 | ★ **49.0** *(printed 49.0 ✔ exact)* | **60.6** *(printed 57.0 / 58.2 — see ⚠[!]-1)* | ★ **84.6** |

★ **Panel C is the second, independent demonstration of the same temperature law as panel A's
covalent half: 3-methylbutanal's covalent loss over 96 h is 15.4 pp at 279 K, 39.4 pp at 300 K and
51.0 pp at 333 K — monotone, positive, ordinary.** The text's *"remained constant"* at 6 °C is a
15-pp overstatement, but the ordering is right.

#### PANEL D — HDMF (scenario I control)

Flat at **93–100 % free at every time and every temperature** (triangle 99.2–100.7 across the whole
grid; square 93.5 at 96 h/333 K; diamond 94.7–98.8). **⇒ [D] no interaction at any temperature,
confirming the text's ">95 %" and ">96 %".** HDMF is the negative control that shows the effects in
panels A–C are not an artefact of temperature, of the second spectrometer, or of the 20-s shake.

#### **Cross-validation: Figure 2's FFT row vs Figure 4A's 300 K trace** `[D]`
Two independent digitisations, of two independently drawn graphics, of the same underlying qHNMR
series:

| t (min) | 5 | 30 | 60 | 120 | 240 | 720 | 1440 | 2880 | 5760 |
|---|---|---|---|---|---|---|---|---|---|
| Fig. 2 (heat map, colour LUT) | 34 | 34 | 32 | 34 | 31 | 28 | 24 | 15 | 2 |
| Fig. 4A (marker geometry) | 32.3 | 32.3 | 32.3 | 32.3 | 31.4 | 27.2 | 23.4 | 15.0 | ~0 |
| \|Δ\| | 1.7 | 1.7 | 0.3 | 1.7 | 0.4 | 0.8 | 0.6 | 0.0 | 2 |

**Mean \|Δ\| = 1.0 pp, max 2.0 pp. ⇒ Both digitisations are sound; the ± error bars quoted above are
if anything conservative.**

### 5c. Figure 1 (p. 15338) — qualitative
*"Excerpts of qHNMR spectra (500 MHz; H₂O/D₂O, 9/1, v/v; pH 5.5; 300 K) showing the specific
resonance signals of selected key coffee odorants used for the determination of interaction
behavior. Faded lines represent aqueous references of odorants, and bold lines correspond to samples
with added melanoidins (>10 kDa) in the original coffee concentration."*
**Not digitised: no numeric axis carrying a parameter.** Its content is the eight screening
compounds of §3b (#1, 2, 4, 7, 8, 11, 12, 18). The extractable numbers are already in §4b #8–#10.

### 5d. Figure 3 (p. 15340) — qualitative, but the load-bearing mechanism figure
*"Excerpts of ¹H NMR spectra (500 MHz; H₂O/D₂O, 9/1, v/v; pH 5.5; 300 K) of substitution experiments
for the distinction between covalent and noncovalent aroma binding, demonstrating the formation of
new resonance signals during the storage (5 mmol/L each) of FFT (A) and N-acetylcysteine (B) and the
absence of new signals during the incubation of MeFFT (C) melanoidin (>10 kDa) mixtures at RT at
different time steps (t1: 5 min, t2: 12 h, t3: 24 h, t4: 48 h, t5: 96 h) due to their free (A, B) and
protected (C) thiol groups."*
**Not digitised: no quantitative axis.** It is the visual proof of the §3c dissection, and its time
grid (5 min / 12 / 24 / 48 / 96 h) is a clean subset of the master grid.

---

## §6. ★★ THE TEMPERATURE CONTRAST — **THE DELIVERABLE**

### 6a. Everything the paper says about 279 K vs 333 K, verbatim, with anchors

**Abstract (p. 15334):**
> *"Evaluation of temperature influence on e.g. 2-furfurylthiol (FFT), revealed an altered behavior
> with **increased π−π stacking at lower temperatures and accelerated covalent interactions at
> higher temperatures**."*

**Rationale and setup (p. 15341):**
> *"it seems to be important to evaluate the impact of temperature on these interactions since
> coffee is usually consumed at a temperature of 60−70 °C or it is stored at low temperatures. The
> influence of high and low temperature on odorant−melanoidin interactions was examined via qHNMR by
> means of incubation experiments at 60 °C, mimicking coffee consumption temperature, and at 6 °C,
> closely mimicking storage conditions for ready-to-drink beverages in the refrigerated shelf. To
> realize high- and low-temperature experiments of the coffee constituents under real conditions,
> **the operating temperature of the NMR probe was changed to 279 or 333 K.** Measurements were done
> using four odorants, namely, 4-hydroxy-2,5-dimethyl-3(2H)-furanone, 2,3-diethyl-5-methylpyrazine,
> 3-methylbutanal, and FFT, one representing each scenario to cover all possible interaction types."*

**The aldehyde (scenario II, covalent-only) — p. 15341:**
> *"At 6 °C, the amount of free 3-methylbutanal **remained constant** within the time course of 96 h
> and was not influenced by prolonged incubation times… The constant odorant concentration at low
> temperature strongly resembles the concentration course of noninteracting
> 4-hydroxy-2,5-dimethyl-3(2H)-furanone and therefore can be interpreted as **an effective method to
> prevent covalent aldehyde addition to the coffee matrix**… At an elevated temperature of 60 °C,
> 3-methylbutanal exhibited **more covalent interactions**, with the amount of free odorant declining
> to **49.0 %** within 96 h, while a decrease in the amount of the free aldehyde to **57.0 %** was
> observed within the same incubation time at RT… As expected, these findings clearly indicate that
> **the higher temperatures serve as an accelerant for covalent binding to coffee melanoidins**."*

**The pyrazine (scenario III, π−π-only) — p. 15342:**
> *"The amount of free pyrazine in the samples went down from **58.5 % at RT to 51.8 %**, which means
> that **a higher ratio of π−π interactions could be observed at lower temperature**… These findings
> match reports from the literature and describe **the strongest noncovalent π−π interactions at low
> temperatures**, where π−π complexes are usually characterized.⁵¹ In conclusion, analytes in
> scenario III (pyrazines, methoxyphenols, and terpenes) will presumably exhibit **higher π−π
> binding affinities at lower temperatures**. In contrast, at 60 °C, π−π stacking exhibited by
> 2,3-diethyl-5-methylpyrazine and consequently by analytes of scenario III was **reduced** in
> comparison to experiments at RT, as evidenced by a residual free odorant concentration of
> **70.0 % at 60 °C and 58.5 % at RT**."*

**FFT at 279 K — p. 15342:**
> *"The odorant concentrations of the aromatic thiol FFT followed a similar trend. At 6 °C, **the
> initial drop in the concentration of the free odorant, explained by immediate π−π stacking after
> melanoidin addition, was more pronounced, resulting in an amount of free FFT of 24.7 % compared to
> 32.3 % at RT.** Additionally, **the steady decline of FFT examined at longer incubation times was
> almost completely absent in the measurements performed at 6 °C** (Figure 4). Obviously, **the
> interaction behavior of FFT shifted from scenario IV at RT to scenario III at lower temperatures**,
> with an amount of free FFT of **24.1 % at 60 min and 23.6 % at 48 h**. The interactions of FFT and
> HMW melanoidins at lower temperatures led to a shift of a **higher ratio of reversible π−π
> interactions, whereas the covalent binding of FFT to the HMW material was almost completely
> suppressed.**"*

**FFT at 333 K — p. 15342:**
> *"At 60 °C, **the immediate decline in the amount of free FFT to 55.0 % was less evident compared
> to RT, where only 32.2 % free FFT was found at the same incubation time.** At longer incubation
> times, **the decay of FFT increased** in comparison to the measurement at RT. For example, at
> 60 °C, **no free FFT could be detected after 72 h** of incubation, and in comparison, **it took
> 96 h at RT** until no more free FFT could be detected. The faster loss of free FFT resulted in a
> **shift of FFT−melanoidin interactions from scenario IV at RT toward scenario II at 60 °C**."*

**The synthesis — p. 15342:**
> *"At lower temperature, covalent binding of capable odorants to the HMW fraction of coffee was
> clearly suppressed, whereas at higher temperatures, advanced covalent aroma binding could be
> observed. **At higher temperatures, π−π interactions were less favored**, but they still had a
> considerable influence on the analytes of scenarios III and IV."*

### 6b. ★ THE CENTRAL OBSERVATION, STATED AS BALDLY AS IT DESERVES `[M]`

> ### **At 60 minutes, MORE free FFT survives at 333 K (55.0 %) than at 279 K (24.1 %). Free thiol at short times is a DECREASING function of temperature. No sink with a positive activation energy can produce that.**

| T | free FFT @ 60 min | **loss @ 60 min** |
|---|---|---|
| **279 K (6 °C)** | **24.1 %** | **75.9 %** |
| **300 K (RT)** | **32.3 %** | **67.7 %** |
| **333 K (60 °C)** | **55.0 %** | **45.0 %** |

**Monotone, and the wrong way round for a chemical sink.** The paper's explanation is that this
short-time loss is *not* chemistry: it is **reversible π−π sequestration**, whose association
weakens with temperature exactly as a physical binding equilibrium must.

### 6c. Loss accounting at each temperature `[D]`

Using the 60-min plateau as the π-equilibrium level and the subsequent decline as the covalent
channel (the paper's own decomposition):

| | **279 K** | **300 K** | **333 K** |
|---|---|---|---|
| free FFT at 60 min | 24.1 % | 32.3 % | 55.0 % |
| **⇒ π-sequestered pool (t = 60 min)** | ★ **75.9 %** | **67.7 %** | ★ **45.0 %** |
| free FFT at 48 h | **23.6 %** | 15.0 % | **0 %** |
| free FFT at 96 h | ★ **22.6 %** | **0 %** | **0 %** |
| **total FFT lost by 96 h** | **77.4 %** | **100 %** | **100 %** |
| **⇒ covalent loss by 96 h, as % of initial FFT** | ★ **≈ 1.5–1.9 %** | **32.3 %** | **≥ 55 %** |
| **⇒ covalent loss as % of the free pool** | **≈ 7.8 %** | **100 %** | **100 %** |
| of which: **covalent melanoidin incorporation** | *(not measured)* | ★ **11 %** of initial | *(not measured)* |
| of which: **difurfuryl disulfide (dimerisation)** | *(not measured)* | ★ **≈ 21 %** of initial `[D]` = 32.3 − 11 | *(not measured)* |

⚠ **A structural point that matters more than any of these numbers:** at 300 and 333 K the free pool
goes to **zero**, i.e. **the π-sequestered FFT is eventually consumed too.** The π complex is
therefore **a reversible reservoir that FEEDS the covalent sink, not a parallel terminal sink.**
At 279 K the reservoir is largest (75.9 %) and the drain is shut; at 333 K the reservoir is smallest
(45.0 %) and the drain is wide open. **The 45/68/76 % split is a 60-minute snapshot, not a permanent
allocation.**

### 6d. **(a) THE COVALENT CHANNEL: fold change and two-point apparent Eₐ** `[D]`

**Rate definition.** Normalise each trace by its own 60-min plateau, `R(t) = f(t)/f_plateau`, so the
π equilibrium is divided out, and fit a pseudo-first-order constant through the origin over
t ≥ 240 min: `k = Σ tᵢ·(−ln Rᵢ) / Σ tᵢ²` (points at R = 0 excluded, being −∞).

| T | R(240) | R(720) | R(1440) | R(2880) | R(5760) | **k_cov (min⁻¹)** |
|---|---|---|---|---|---|---|
| **279 K** | 1.000 | 1.000 | 1.053 | 0.959 | 0.922 | ★ **1.17 × 10⁻⁵** |
| **300 K** | 0.972 | 0.842 | 0.724 | 0.464 | → 0 | ★ **2.57 × 10⁻⁴** |
| **333 K** | 0.844 | 0.700 | 0.291 | → 0 | → 0 | ★ **7.83 × 10⁻⁴** |

**Fold changes:** 279 → 333 K: **× 67.2** · 279 → 300 K: **× 22.0** · 300 → 333 K: **× 3.05**.

**The arithmetic, explicitly, for the requested 279 → 333 K pair:**

```
Ea = R · ln(k2/k1) / (1/T1 − 1/T2)

k1 = k(279 K) = 1.165e-5 min^-1
k2 = k(333 K) = 7.832e-4 min^-1
k2/k1 = 67.2                         ln(k2/k1) = 4.2077

1/T1 = 1/279 = 3.584229e-3 K^-1
1/T2 = 1/333 = 3.003003e-3 K^-1
1/T1 - 1/T2 = 5.81226e-4 K^-1

Ea = 8.314 J mol^-1 K^-1 × 4.2077 / 5.81226e-4 K^-1
   = 34.983 / 5.81226e-4  J/mol
   = 60,188 J/mol
```

> ### ★ **Ea(covalent, 279 → 333 K) = 60.2 kJ/mol.**

**Every other route to the same quantity, for honesty:**

| route | pair | **Ea (kJ/mol)** |
|---|---|---|
| two-point, plateau-normalised | **279 → 333 K** | ★ **60.2** |
| two-point, plateau-normalised | 279 → 300 K | **102.4** |
| two-point, plateau-normalised | 300 → 333 K | **28.1** |
| **three-point Arrhenius least squares** | all three | ★ **58.5** (slope −7041 K, **R² = 0.886**) |
| time-to-extinction, Fig. 4 (2880 vs 5760 min) | 300 → 333 K | **17.4** |
| time-to-extinction, text (72 h vs 96 h) | 300 → 333 K | **7.2** |
| **π-stacking component if wrongly Arrhenius-fitted as a rate** | 279 → 333 K | ★ **−19.5 (NEGATIVE)** |

> ### ⇒ **Central estimate ≈ 60 kJ/mol. Defensible range 7–102 kJ/mol. The ladder is strongly NON-ARRHENIUS — R² = 0.886 on three points, and the two segments disagree by 3.6× (102 vs 28 kJ/mol).**

### 6e. **(b) THE π-STACKING COMPONENT: fold change AND SIGN** `[D]`

Treat the 60-min plateau as an association equilibrium at fixed melanoidin loading and form the
dimensionless `K_app = θ/(1−θ)` where θ is the bound fraction:

| T | θ (bound) | **K_app** | ln K_app | 1/T (K⁻¹) |
|---|---|---|---|---|
| **279 K** | 0.759 | **3.149** | 1.1472 | 3.584229e-3 |
| **300 K** | 0.677 | **2.096** | 0.7402 | 3.333333e-3 |
| **333 K** | 0.450 | **0.818** | −0.2008 | 3.003003e-3 |

**Fold changes, 279 → 333 K:**
- bound fraction θ: **75.9 % → 45.0 % = × 0.593** (−30.9 pp over 54 K = **−0.57 pp K⁻¹**)
- **K_app: 3.149 → 0.818 = × 0.26, i.e. a 3.85-FOLD DECREASE**

**van 't Hoff (`ln K = −ΔH/RT + ΔS/R`):**
- two-point 279 → 333 K: slope = −1.3480 / −5.81226e-4 = 2319.3 K ⇒ **ΔH° = −19.3 kJ/mol**
- **three-point least squares: slope = 2346.4 K ⇒ ΔH° = −19.5 kJ/mol, R² = 0.979**

> ### ★ **THE SIGN IS NEGATIVE. ΔH°(π−π) = −19.5 kJ/mol, R² = 0.979 — an exothermic, entropically-opposed physical association that WEAKENS as temperature rises. The π channel and the covalent channel have OPPOSITE temperature coefficients.**

**Independent confirmation from the pyrazine (scenario III, pure π−π):** bound fraction
**48.2 % (279 K) → 41.5 % (300 K) → 30.0 % (333 K)** ⇒ K_app **0.931 / 0.709 / 0.429**;
three-point van 't Hoff ⇒ **ΔH° = −13.4 kJ/mol** `[D]`. **Same sign, same order of magnitude, on a
different molecule with no thiol at all.** The negative sign is not an FFT artefact.

Also note the fit quality asymmetry: **the π channel is a clean van 't Hoff line (R² = 0.979); the
covalent channel is not a clean Arrhenius line (R² = 0.886).** The physical association behaves;
the chemistry does not.

### 6f. ⚠ **THE CAVEAT CHAIN — attach all seven to any use of §6d/§6e**

1. **A 54 K interval at low absolute temperature has enormous leverage on Eₐ, so the estimate is
   precision-limited.** `1/T1 − 1/T2 = 5.81e-4 K⁻¹`, so **every factor-of-2 error in the rate ratio
   moves Eₐ by 9.9 kJ/mol**, and a 10 % error moves it by 1.4 kJ/mol. With ±1.5 pp digitisation
   error propagating into k, the *statistical* band on the 279→333 pair is roughly **±10 kJ/mol** —
   but the **model-form** error (choice of plateau, choice of rate law, R = 0 exclusions) dominates
   and opens the range to **7–102 kJ/mol**.
2. **An association equilibrium is not a rate.** §6e's ΔH is a van 't Hoff enthalpy of binding, not
   an activation energy. It must never be entered into an Arrhenius `k(T)` as if it were an Eₐ.
   It is reported as a "negative apparent Eₐ" only to make the *sign contrast* legible.
3. **If the covalent loss is not first-order, the two-point Eₐ is a fold-change dressed as a
   barrier.** The paper never fits a rate law (§9) and the data cannot distinguish first order from
   a site-limited saturating process. **The R² = 0.886 non-linearity is direct evidence that the
   single-first-order assumption is already failing inside the measured window.**
4. **The melanoidin is a heterogeneous polymer, so "the rate constant" lumps many site types.**
   The paper itself calls melanoidins *"a very heterogeneous substance class"* whose constituents
   *"could not be characterized so far"* and which *"contain a plethora of reactive compounds and/or
   side groups"* (p. 15334–35). A lumped k over an unresolved site distribution has no reason to be
   Arrhenius, and its apparent Eₐ drifts as the fastest sites are consumed — which is exactly the
   direction of the 102 → 28 kJ/mol segment drift observed.
5. **Two of the three temperatures were measured on a DIFFERENT SPECTROMETER from the third**
   (AVANCE Neo for 279/333 K; AVANCE III for 300 K, §1b). The 300 K point therefore carries an
   uncontrolled inter-instrument offset. **The 279 → 333 K pair is the only same-instrument
   comparison, which is a further reason to prefer 60.2 kJ/mol over the segment values.**
6. **n = 1.** No replicates, no error bars, no SDs on any NMR value (§2).
7. ### ★ **6–60 °C IS FAR BELOW THE 100–180 °C COOKING WINDOW THE REPO NEEDS. ANY USE AT REACTION-POT TEMPERATURE IS AN EXTRAPOLATION AND MUST BE LABELLED ONE.** Extrapolating k_cov from 333 K to 433 K (160 °C) with Eₐ = 60.2 kJ/mol is a factor of `exp(60200/8.314 × (1/333 − 1/433)) = exp(5.02) = 152×` — over a 100 K gap, on a channel already shown to be non-Arrhenius over 54 K, with a substrate (preformed coffee melanoidin) that does not exist in a cysteine/ribose pot. **Do not do it for a parameter. Do it only for a sign and an order of magnitude.**

---

## §7. ★★ WHAT THIS SAYS ABOUT THE REPO'S THIOL-SINK BARRIER

**The defect under examination.** `results/validation/kinetic_core_b2_2_fit_report.json` carries
`"thiol_sink": 247.97527595559305` kJ/mol, against `DECAY_EA_BOUNDS = (20.0, 250.0)` declared in
`src/kinetic_core/parameters_sulfur.py:1819`. **The fit is pressed to within 0.8 % of its own
ceiling** — the classic signature of a parameter absorbing structure it cannot represent.
(The B2.3 successor does the same thing on a different channel: `"carbonyl_sink": 249.916…`.)

### 7a. Is 248 kJ/mol consistent with anything in this paper? **No — by five to six orders of magnitude.**

**The falsification is arithmetic, not rhetorical.** An Eₐ of 248 kJ/mol predicts a rate ratio over
the measured window of

```
k(333)/k(279) = exp( 248000 × 5.81226e-4 / 8.314 ) = exp(17.34) = 3.4 × 10^7
```

**The measured ratio is 67.2.** The prediction is wrong by a factor of **5.1 × 10⁵**.
Over the shorter 279 → 300 K span, 248 kJ/mol predicts **× 1780**; measured **× 22.0** — wrong by
**81×**. There is no reading of this dataset, no choice of plateau, and no rate law that recovers a
248 kJ/mol barrier from a 67-fold change over 54 K.

### 7b. What Eₐ range does this paper's covalent channel actually support?

| statement | value |
|---|---|
| **central estimate (279 → 333 K, same instrument, plateau-normalised)** | ★ **60 kJ/mol** |
| three-point Arrhenius, R² = 0.886 | 58.5 kJ/mol |
| **full defensible range across all six derivation routes** | ★ **7 – 102 kJ/mol** |
| upper bound compatible with the 279 K trace being *"almost completely suppressed"* rather than *exactly zero* | ≲ 110 kJ/mol |

**For orientation, this range is unremarkable company:** the corpus's own measured
cysteine-depletion Eₐ from Kang 2026 is **55.1 kJ/mol** (`kang2026_SI_extraction.md` §5b, and
`parameters_sulfur.py` `MEASURED_EA_OVERRIDES` already carries it for `k_cys_thermal`). **Gigl's
melanoidin–thiol covalent channel lands within 10 % of that independently measured value.** That is
a striking, previously unnoticed concordance between two entirely unrelated experiments — one a
sealed aqueous cysteine pyrolysis at 100–140 °C, the other an NMR binding study at 6–60 °C — and it
is the single strongest quantitative reason to regard **~55–60 kJ/mol as the natural scale for
thiol-consuming chemistry in this network**, and 248 as an artefact.

### 7c. Does a temperature-SWITCHED mechanism support or undermine a SINGLE-Eₐ lumped sink?

> ### **It undermines it, structurally, and it also explains how a fitter gets to 248.**

The sink the repo calls `thiol_sink` is being asked to represent, at minimum, **two channels with
opposite temperature coefficients**:

| channel | ΔH or Eₐ | sign | reversible? | time constant |
|---|---|---|---|---|
| π−π sequestration onto melanoidin | **−19.5 kJ/mol** | ★ **negative** | **yes** | **< 5 min** |
| covalent binding + disulfide dimerisation | **+60 kJ/mol** | positive | no | **hours–days** |

Three consequences the repo should take seriously:

1. **No single positive Eₐ can reproduce the short-time data at all.** At 60 min the free-thiol
   fraction *rises* with temperature (24.1 → 32.3 → 55.0 %). A single first-order sink with
   `Ea > 0` is structurally incapable of producing a *negative* temperature response in thiol loss,
   at any prefactor. This is not a fit-quality problem; it is a model-form falsification.
2. **A single sink fitted to LONG-time data sees almost no temperature dependence.** Total FFT loss
   at 96 h is 77.4 / 100 / 100 % at 279 / 300 / 333 K — nearly saturated everywhere. A fitter
   targeting long-time yields would recover a *low* Eₐ; a fitter targeting short-time yields would
   need a *negative* one; a fitter targeting a mixture of both, at different temperatures in
   different rows of the panel, is being asked to satisfy contradictory constraints with one
   parameter. ★ **A rail-pinned Eₐ is exactly what that produces: the optimiser drives the one
   free knob to a bound because no interior value can satisfy the ensemble.** The 248-vs-250
   proximity is therefore *diagnostic of lumping*, not of a real high barrier.
3. **The right structural fix is two channels, not a better single Eₐ.** The repo already has the
   pieces: `k_dimer_*` (the disulfide route — Gigl's difurfuryl disulfide is precisely that, and it
   accounts for ~21 % of initial FFT at 96 h/RT), and `k_thiolate_loss` (the pH-gated oxidative
   route from Kumazawa 2003). What is missing is a **reversible, negative-ΔH sequestration term on
   the melanoidin/polymer pool.** Gigl 2021 is the corpus's only measurement of it, and it supplies
   its enthalpy (−19.5 kJ/mol), its equilibration time (< 5 min), and its magnitude
   (45–76 % of ligand bound at 1× brew loading).

### 7d. ⚠ The transfer risk, scored honestly

**What was measured here is a MELANOIDIN sink at BREW temperature, not the sink operating in a
100–180 °C reaction pot.** The differences are not cosmetic:

| axis | Gigl 2021 | the repo's thiol sink |
|---|---|---|
| temperature | **279–333 K** | **373–453 K** — *no overlap at all* |
| binding partner | **preformed, purified, >10 kDa coffee melanoidin, isolated then re-added** | melanoidin **being formed in situ** from sugars + amino acids |
| thiol | **FFT** (+ 3 aliphatic thiols, + N-acetylcysteine) | FFT, MFT, H₂S, cysteine |
| ligand:site ratio | **~10⁵ above authentic** (§1c) | authentic |
| medium | aqueous, pH 5.5, no reducing sugar, no amino acid, 20-s shake, sealed tube | reacting Maillard system |
| what is measured | **free-ligand disappearance** | free-thiol disappearance |
| replicates | **n = 1** | — |
| rate law | **none fitted** | first-order with Eₐ |

**Transfer-risk score: HIGH for magnitude, LOW for sign, MODERATE for order of magnitude.**
- The **sign** of the covalent channel (positive Eₐ) and the **sign** of the π channel (negative ΔH)
  are physical facts that survive the transfer; they are thermodynamic, not system-specific.
- The **order of magnitude** (tens of kJ/mol, not hundreds) transfers with moderate confidence,
  reinforced by the independent Kang 55.1 kJ/mol concordance (§7b).
- The **absolute rate constant** does not transfer at all and must not be registered as a
  benchmark row.

### 7e. **CALIBRATED VERDICT** — explicit probabilities

Judgements below are this wave's, conditioned on this paper plus the repo state named in §7a:

| claim | **P** |
|---|---|
| The apparent Eₐ of Gigl's melanoidin–thiol **covalent** channel over 279–333 K lies in **30–100 kJ/mol** | **0.85** |
| … lies **below 150 kJ/mol** | **0.97** |
| … is **≥ 248 kJ/mol** | ★ **< 0.01** |
| The π−π component has a **negative** temperature coefficient (ΔH < 0) | **0.97** |
| The repo's `thiol_sink` Eₐ = 248.0 kJ/mol is an **artefact of lumping two opposite-sign channels and/or of the 250 rail**, not a measurement of a real barrier | ★ **0.90** |
| A **single** first-order sink with one positive Eₐ **cannot** reproduce Gigl's 279/300/333 K FFT data even qualitatively | ★ **0.95** |
| The true thiol sink in the repo's **100–180 °C** window has Eₐ in **40–120 kJ/mol** | **0.65** *(lower, because this is the extrapolation — see §6f.7)* |
| Adding a **reversible, negative-ΔH polymer-sequestration** channel would move the fitted `thiol_sink` Eₐ off its 250 rail | **0.70** |
| Gigl 2021 alone is sufficient to **re-fit** a repo parameter | ★ **< 0.05** — it is a **structural** and **refutational** anchor, not a calibration source (§11) |

---

## §8. THE SENSORY HALF `[M]`

### 8a. Panel and conditions, verbatim and complete

| item | value | anchor |
|---|---|---|
| **panel size** | ★ **15 panelists** — *"Fifteen panelists (**8 females, 7 males; 23−32 years in age**), who had given informed consent to participate in the present sensory study and had **no history of anosmia**"* | p. 15336 |
| **training** | *"attended **weekly training sessions** to familiarize themselves with the sensory procedures used and to evaluate aqueous reference solutions of odorants"* | p. 15336 |
| sample | **20 mL** in aqueous phosphate buffer **0.1 mol/L, pH 5.5**, in glass beakers **⌀ 45 mm, 45 mL capacity** | p. 15336 |
| bias control | **odorless sugar color** added to all samples to mask optical differences (*"verified to have no impact on aroma binding prior to sensory evaluation"*) **plus red lighting in the cabins** | p. 15336 |
| room | individual booths, **22−25 °C** | p. 15336 |
| **equilibration** | ★ **60 min at RT** before every sensory presentation | p. 15336–37 |
| **aroma model** | ★ **25 of the most important aroma-active compounds** of coffee, at concentrations **in Table S1 (not on disk)**, per refs 39–40 (Semmelroch & Grosch) | p. 15336–37 |
| melanoidin dose | **HMW >10 kDa at the original coffee concentration** | p. 15337 |
| water | **bottled Evian** for sensory (Milli-Q for NMR) | p. 15335 |

### 8b. Aroma profile — attributes, scale, references, values

**Scale:** **0 (absent / not detectable) → 3 (very strong / strongly detectable)**, integer-anchored.
**Design:** *"performed in two independent sessions, and the results were given as averages ±
standard deviations."*
**Each attribute was presented with a physical reference standard:**

| attribute | **reference compound(s) given to the panel** |
|---|---|
| "sweetish/caramel-like" | **4-hydroxy-2,5-dimethyl-3(2H)-furanone** |
| "earthy" | **2,3-diethyl-5-methoxypyrazine** ⚠ *(note: **methoxy**pyrazine — the compound actually measured by NMR is 2,3-diethyl-5-**methyl**pyrazine. Different molecules. Printed as such; raster-verified.)* |
| **"roasty/sulfurous"** | ★ **a MIXTURE of 2-furfurylthiol + 3-mercapto-3-methylbutanol + 3-mercapto-3-methylbutyl formate** |
| "smoky" | **2-methoxyphenol** |

**Results — Table 1, reproduced from §4a:**

| attribute | without melanoidins | with melanoidins | Δ | rel. Δ |
|---|---|---|---|---|
| **Roasty/sulfury** | **2.2 (0.2)** | **1.6 (0.3)** | **−0.6** | **−27.3 %** |
| **Sweetish/caramel-like** | **1.5 (0.5)** | **1.9 (0.4)** | **+0.4** | **+26.7 %** |
| **Earthy** | **1.6 (0.4)** | **1.3 (0.5)** | **−0.3** | **−18.8 %** |
| **Smoky** | **1.7 (0.6)** | **1.2 (0.3)** | **−0.5** | **−29.4 %** |

**[D] Effect size on the headline attribute.** Pooled SD for roasty/sulfury =
`√((0.2²+0.3²)/2) = 0.255`; **Δ/SD_pooled = 2.35** — a large effect. **No per-attribute
significance test is reported** ⚠, so this is the strongest statement the printed data support.

**Direction check against the NMR:** the three attributes that FALL (roasty −27 %, smoky −29 %,
earthy −19 %) are exactly the three whose reference compounds are the **strongly-interacting**
scenario II/III/IV species (FFT, guaiacol, pyrazine). The one attribute that **RISES** is anchored
on **HDMF, the only compound in the entire study that shows no interaction whatsoever.**
**⇒ The sensory result and the molecular result are mutually consistent with no free parameters.**
The rise in "sweetish/caramel-like" is a *relative unmasking*, not a chemical increase.

### 8c. The three-alternative forced choice test

**Design:** three samples — **two** = aroma reconstitution model; **one** = model **+** HMW >10 kDa
at original coffee concentration. All incubated **60 min at RT**. Assessors asked to identify the
deviating sample **and describe its overall difference**. **Two independent sessions.**
**Result (p. 15337):**
> *"The results clearly showed that the assessors were able to correctly identify the sample with
> added melanoidins, therefore demonstrating a **significant (α = 0.05)** influence of HMW components
> on the overall odor of coffee brew."*

⚠ **[NEG] no correct-response count, no n, no p-value, no test statistic.** Only the α threshold is
printed. **Register as a qualitative significance claim only.**

### 8d. ★ **[D] THE DOSE–RESPONSE ANCHOR THE REPO HAS BEEN SHORT OF**

This is the corpus's rare case where **a quantified thiol loss and a quantified perceived-intensity
loss are measured on the same system at the same time point (60 min, RT, 1× melanoidin)**:

| | value |
|---|---|
| free FFT remaining at 60 min | **32.3 %** |
| perceived "roasty/sulfury" intensity ratio | **1.6 / 2.2 = 72.7 %** |

Fitting a Stevens power law `I ∝ C^n`:
- **If the percept tracked FFT alone:** `n = ln(0.727)/ln(0.323) = (−0.3185)/(−1.1301) = ` **0.28**
- **If the percept tracks the equally-weighted 3-component reference mixture** (FFT 32 %,
  3-mercapto-3-methylbutanol 93 %, 3-mercapto-3-methylbutyl formate 87 % free at 60 min, from §5a
  ⇒ mean **70.7 %** free): `n = ln(0.727)/ln(0.707) = ` **0.92**

> **⇒ The exponent is bracketed at n ≈ 0.3–0.9. Use the bracket, never a point value.**

⚠ Three caveats, all fatal to over-precision: (i) the roasty reference is a **mixture**, and its
three members are bound to wildly different extents, so the effective stimulus is ill-defined;
(ii) a **0–3 category scale is not a ratio scale**, and Stevens exponents from category scales are
systematically compressed; (iii) **n = 15 panelists, 2 sessions, no per-attribute test.**
**Register as PRIOR-ONLY.**

---

## §9. VERIFIED NEGATIVES `[NEG]` — each one checked, not assumed

| # | question asked of the paper | **verdict** | how checked |
|---|---|---|---|
| 1 | **Any rate constant printed?** | ★ **NO. Not one.** No k, no half-life, no first-order plot, no rate law is fitted anywhere in the article. Every temporal statement is a percentage at a time point. | full text layer + all four figures |
| 2 | **Any activation energy?** | ★ **NO.** The string "activation" does not occur. No Arrhenius analysis, no Eₐ, no ΔH‡, no ΔG‡. **All Eₐ values in this dossier are §6 derivations `[D]`, none is the authors'.** | full text |
| 3 | **Any measurement above 60 °C?** | ★ **NO. 333 K (59.85 °C) is the absolute maximum in the paper.** The only higher temperature mentioned anywhere is *"coffee is usually consumed at a temperature of 60−70 °C"* — an aspiration in the rationale sentence, **not a measurement**; the probe was set to 333 K. | p. 15341–42 |
| 4 | **Any MFT (2-methyl-3-furanthiol)?** | ★ **NO MEASUREMENT.** MFT appears **exactly once**, in the Introduction, in a `[C]` list of storage-labile coffee thiols (*"3-methyl-2-butene-1-thiol, 3-mercapto-3-methylbutyl formate, 2-methyl-3-furanthiol, and methanethiol"*, p. 15334, citing ref. 1). **It is not in the Chemicals list, not in the aroma model, not in the heat map, and not in any figure.** ⇒ **This paper says nothing whatever about MFT. Do not attribute anything to it.** | Chemicals list, Fig. 2, Fig. 4 |
| 5 | **Melanoidin thiol-binding capacity per gram?** | ★ **NO.** No mg/mL, no mg/g, no site density, no binding capacity, no saturation curve, no titration. Loading is only ever *"the original coffee concentration"*. **The most quantitative statement possible is the §1a derivation that the tube contains exactly 1× brew-equivalent melanoidin.** | p. 15335–36 |
| 6 | Any binding constant (K_a, K_d, K)? | **NO.** No isotherm, no Scatchard, no stoichiometry, no Job plot. | full text |
| 7 | Replicate count / error bars on the NMR? | **NO.** No n, no SD, no error bars in Fig. 2 or Fig. 4. See §2. | full text, both figures |
| 8 | LOD / LOQ / integral precision? | **NO number.** Only *"in the error range of the analysis"* citing ref. 35 `[C]`. | p. 15340 |
| 9 | pH/pD isotope correction? | **NO.** The 10 % D₂O buffer's pH was set with DCl and reported as "pH 5.5" with no `pD = pH_read + 0.4`-type correction discussed. | p. 15336 |
| 10 | Temperature calibration of the probe? | **NO.** 30 min equilibration is stated; no shift thermometry, no accuracy figure. | p. 15336 |
| 11 | Funding / acknowledgement? | **NO SECTION AT ALL.** | full text |
| 12 | Any 72-h data point plotted? | **NO** — the Fig. 4 grid ends 2880 → 5760 min. See ⚠[!]-4. | Fig. 4, marker detection |
| 13 | Structural identification of the covalent FFT–melanoidin adduct? | **NO.** The 11 % is a **mass-balance residue** ("could neither be assigned to noncovalent interactions nor to the formation of the FFT dimer… **most likely** covalently incorporated"), not an isolated or characterised species. **Only the disulfide dimer is positively identified**, by spiking against the commercial reference. | p. 15341 |
| 14 | Anything about H₂S, cysteine thermolysis, or sugar? | **NO.** There is no reducing sugar and no free amino acid in any experiment. This is a **binding** study, not a Maillard study. | Chemicals, Methods |

---

## §10. CONSOLIDATED PARAMETER TABLE

`[M]` printed · `[D]` derived/digitised here · `[C]` cited · `[NEG]` verified absent.

| # | parameter | value | units | condition | prov. | source anchor |
|---|---|---|---|---|---|---|
| 1 | brew strength | **5.4** | g/100 mL | Nespresso capsule, 5.6 g / 104 mL | `[M]` | p. 15335 |
| 2 | melanoidin cut-off | **10** | kDa | cross-flow UF, HMW = retentate | `[M]` | p. 15335 |
| 3 | melanoidin loading in NMR tube | **1.00** | × original brew | 150 μL of ×4 into 600 μL | `[D]` | p. 15336 |
| 4 | melanoidin loading | ★ **not stated** | mg/mL | — | `[NEG]` | §9 #5 |
| 5 | odorant spike | **5** | mmol/L | all NMR tubes | `[M]` | p. 15338 |
| 6 | pH | **5.5** | — | phosphate/DCl, H₂O/D₂O 9:1 | `[M]` | p. 15336 |
| 7 | temperatures | **279 / 300 / 333** | K | probe set-point, 30 min equil. | `[M]` | p. 15335–36 |
| 8 | time grid | **0, 5, 30, 60, 120, 240, 720, 1440, 2880, 5760** | min | RT series; same grid at 279/333 K | `[D]` raster | Fig. 2 axis, Fig. 4 |
| 9 | **FFT free @ 60 min, 279 K** | ★ **24.1** | % of initial | 1× melanoidin, pH 5.5 | `[M]` | p. 15342 |
| 10 | **FFT free @ 60 min, 300 K** | ★ **32.3** | % | " | `[M]` | p. 15340, 15342 |
| 11 | **FFT free @ 60 min, 333 K** | ★ **55.0** | % | " | `[M]` | p. 15342 |
| 12 | FFT free @ 48 h, 279 K | **23.6** | % | " | `[M]` | p. 15342 |
| 13 | FFT free @ 48 h, 300 K | **15.0** | % | " | `[D]` Fig. 4A | §5b |
| 14 | FFT free @ 48 h, 333 K | **0** | % | " | `[D]` Fig. 4A | §5b, ⚠[!]-4 |
| 15 | FFT free @ 96 h, 279 K | ★ **22.6** | % | " | `[D]` Fig. 4A | §5b |
| 16 | FFT free @ 96 h, 300 K | **0** | % | " | `[M]` | p. 15340 |
| 17 | FFT free @ 96 h, 333 K | **0** | % | " | `[D]` Fig. 4A | §5b |
| 18 | **π-bound FFT fraction, 279 K** | ★ **75.9** | % of initial | 60 min plateau | `[D]` | §6c |
| 19 | **π-bound FFT fraction, 300 K** | **67.7** | % | " | `[D]` | §6c |
| 20 | **π-bound FFT fraction, 333 K** | ★ **45.0** | % | " | `[D]` | §6c |
| 21 | **ΔH°(π−π, FFT)** | ★ **−19.5** | kJ/mol | 279–333 K, van 't Hoff, **R² = 0.979** | `[D]` | §6e |
| 22 | ΔH°(π−π, 2,3-diethyl-5-methylpyrazine) | **−13.4** | kJ/mol | 279–333 K | `[D]` | §6e |
| 23 | **k_cov(FFT), 279 K** | **1.17 × 10⁻⁵** | min⁻¹ | plateau-normalised pseudo-1st-order | `[D]` | §6d |
| 24 | **k_cov(FFT), 300 K** | **2.57 × 10⁻⁴** | min⁻¹ | " | `[D]` | §6d |
| 25 | **k_cov(FFT), 333 K** | **7.83 × 10⁻⁴** | min⁻¹ | " | `[D]` | §6d |
| 26 | **fold change k_cov, 279 → 333 K** | ★ **67.2** | × | 54 K span | `[D]` | §6d |
| 27 | ★ **Eₐ(covalent FFT–melanoidin), two-point** | ★ **60.2** | kJ/mol | 279 → 333 K | `[D]` | §6d |
| 28 | Eₐ(covalent), 3-point Arrhenius | **58.5** | kJ/mol | R² = 0.886 | `[D]` | §6d |
| 29 | Eₐ(covalent), defensible range | **7 – 102** | kJ/mol | across six derivation routes | `[D]` | §6d |
| 30 | FFT covalently incorporated into melanoidin | ★ **5 → 11** | % of initial | 48 h → 96 h, 300 K | `[M]` | p. 15341 |
| 31 | FFT → difurfuryl disulfide (dimer) | **≈ 21** | % of initial | 96 h, 300 K, = 32.3 − 11 | `[D]` | §6c |
| 32 | N-acetylcysteine covalently incorporated | **8** | % of initial | 96 h, 300 K | `[M]` | p. 15341 |
| 33 | pyrazine free, 279 / 300 / 333 K | **51.8 / 58.5 / 70.0** | % | scenario III, flat in time | `[M]` | p. 15342 |
| 34 | 3-methylbutanal free @ 96 h, 279 / 300 / 333 K | **84.6 / ≈60 / 49.0** | % | scenario II | `[D]/[M]` | §5b, p. 15341 |
| 35 | HDMF free, all T, all t | **> 95** | % | scenario I negative control | `[M]/[D]` | p. 15341, §5b |
| 36 | FFT H−C(3) Δδ on binding | **+0.56** | Hz @ 500 MHz (**= +0.00112 ppm** `[D]`) | 60 min, 300 K | `[M]` | p. 15339 |
| 37 | FFT H−C(3) line broadening | **+0.93** | Hz | " | `[M]` | p. 15339 |
| 38 | MeFFT Δδ on binding | **+0.51** | Hz | thiol-blocked control | `[M]` | p. 15341 |
| 39 | pyrazine H−C(6) fwhm | **2.02 → 4.16** (Δ +2.14) | Hz | 60 min, 300 K | `[M]` | p. 15339 |
| 40 | pyrazine H−C(6) Δδ | **+1.42** | Hz | " | `[M]` | p. 15339 |
| 41 | guaiacol H−C(3) Δδ / Δfwhm | **+1.52 / +1.55** | Hz | " | `[M]` | p. 15339 |
| 42 | **sensory: roasty/sulfury** | ★ **2.2 → 1.6** (SD 0.2 / 0.3) | 0–3 scale | ±melanoidin, 60 min, n = 15, 2 sessions | `[M]` | Table 1 |
| 43 | sensory: sweetish/caramel-like | **1.5 → 1.9** (0.5 / 0.4) | 0–3 | " | `[M]` | Table 1 |
| 44 | sensory: earthy | **1.6 → 1.3** (0.4 / 0.5) | 0–3 | " | `[M]` | Table 1 |
| 45 | sensory: smoky | **1.7 → 1.2** (0.6 / 0.3) | 0–3 | " | `[M]` | Table 1 |
| 46 | sensory effect size, roasty | **Δ/SD_pooled = 2.35** | — | " | `[D]` | §8b |
| 47 | **Stevens exponent, roasty vs free thiol** | ★ **0.3 – 0.9** | — | bracket, not a point | `[D]` | §8d |
| 48 | 3-AFC | **significant, α = 0.05** | — | no n, no p printed | `[M]` | p. 15337 |
| 49 | authentic-brew pyrazine binding | **51** | % free | 60 min, freshly percolated coffee | `[M]` | p. 15342 |
| 50 | authentic-brew pyrazine concentrations | ~~47.8 → 24.4~~ | mmol/L | ★ **DO NOT REGISTER** | ⚠[!]-3 | §4c |

---

## §11. USABILITY VERDICTS

| block | verdict | reason |
|---|---|---|
| **§3b — the 20-odorant scenario assignment table** | ★ **USE (STRUCTURAL)** | The single most reusable artefact in the paper. Class-level, mechanism-level, backed by two orthogonal substitution controls. Not a number, so nothing to mis-transfer. **This is what the wave should carry forward.** |
| **§3c — the N-acetylcysteine / MeFFT dissection** | ★ **USE (STRUCTURAL)** | Establishes that "thiol sink" is two separable channels, one of them not chemistry. Directly actionable on `parameters_sulfur.py`. |
| **§6b — free thiol rises with temperature at short times** | ★ **USE (STRUCTURAL, REFUTATIONAL)** | A model-form falsification of any single positive-Eₐ sink. Needs no numerical transfer to bite. |
| **§6e — ΔH°(π−π) = −19.5 kJ/mol, sign and magnitude** | **USE-Q** | Clean van 't Hoff (R² = 0.979), confirmed on a second molecule. Qualified: n = 1, equilibrium not rate, 1× brew loading only. **Register the SIGN as hard; the magnitude as ±10 kJ/mol.** |
| **§6d — Eₐ(covalent) = 60 kJ/mol** | **USE-Q, as a CEILING TEST only** | Its job in this repo is to refute 248, which it does by 5 orders of magnitude. It is **not** a calibration source for a 100–180 °C sink (§7d). Carry the *range* 7–102, never the point value, into any prior. |
| **§5a — Figure 2 heat-map table (20 × 10)** | **USE-Q** | ±4 pp, validated on 5 printed anchors, n = 1, 300 K only. Good for ordering and for class-level qualitative gates; too coarse for a rate fit. |
| **§5b — Figure 4 panels A–D** | **USE-Q** | ±1.5 pp, validated on 5 printed anchors and cross-validated against Fig. 2 to 1.0 pp mean. Same n = 1 caveat. **Panel A is the wave's payload.** |
| **§4b items 12–14 — covalent incorporation 5 % / 11 % / 8 %** | **RATIO-ONLY** | Mass-balance residues, not measured species; 300 K only; no temperature ladder; no structure. Usable as a *fraction of the sink that is polymer-incorporation rather than dimerisation*, at RT. |
| **§4b item 29 / §10 #50 — the authentic-brew pyrazine** | ★ **RATIO-ONLY** | 51 % is sound; 47.8 mmol/L is chemically impossible (⚠[!]-3). **Never register the absolute concentration.** |
| **§8b — Table 1 sensory values** | **USE-Q** | Genuine measured values with SDs and n. Qualified: no per-attribute significance test, category scale, aroma-model composition in an SI not on disk. |
| **§8d — the Stevens exponent 0.3–0.9** | ★ **PRIOR-ONLY** | Bracket spans 3×; the reference standard is a mixture with heterogeneous binding; category scale compresses exponents. Useful as a shape prior for "thiol loss → roastiness loss", nothing more. |
| **§4b items 8–11 — chemical shifts and fwhm** | **STRUCTURAL** | Diagnostic evidence for π-stacking. No kinetic content. Keep as the mechanistic citation, not as parameters. |
| **3-methylbutanal at 96 h / 300 K** | ★ **REFUSE the printed values, USE ≈ 60 % from the figures** | Four inconsistent values for one quantity (⚠[!]-1); the authors' own concentrations contradict their own percentage. |
| **The 72-h extinction claim at 333 K** | ★ **REFUSE** | Contradicted by Fig. 4A, which has no 72-h point and shows zero at 48 h (⚠[!]-4). Use the figure. |
| **Anything about MFT** | ★ **REFUSE** | MFT is never measured (§9 #4). |
| **Anything above 60 °C** | ★ **REFUSE** | Nothing was measured there (§9 #3). |
| **Any absolute rate constant for the repo** | ★ **REFUSE** | The paper fits none (§9 #1), the derived ones are model-form-dependent, and the system does not exist in a reaction pot (§7d). |

---

## §12. DECLARED GAPS

### 12a. What is missing from the published record itself (no SI would fix these)

1. **★ No melanoidin concentration in mg/mL, and no binding capacity per gram.** This is the single
   biggest obstacle to porting anything from this paper into a concentration-dependent sink term.
   Everything is expressed relative to "original coffee concentration".
2. **★ No replicates, no error bars, no SDs on any NMR measurement.** n = 1 as published.
3. **No rate law and no rate constant.** Every k and Eₐ in this dossier is a `[D]` construction.
4. **No 279 K or 333 K covalent-incorporation split.** The 5 %/11 %/8 % mass-balance numbers exist
   only at 300 K, so the π-vs-covalent partition at the two extreme temperatures had to be derived
   from plateau geometry (§6c) rather than read off.
5. **No structural characterisation of the covalent FFT–melanoidin adduct.**
6. **No intermediate temperatures.** Three points, two of them on one instrument and one on another,
   is the minimum for an Arrhenius line and gives R² = 0.886. A fourth point at ~315 K would
   discriminate the 102-vs-28 kJ/mol segment conflict decisively.
7. **No authentic-concentration experiment.** Everything runs at ~10⁵ × real ligand levels (§1c).
8. **No temperature calibration and no pD correction.**

### 12b. ★ What the SI would add, and whether it is needed

| SI item | what it contains | **would it change anything?** | **order it?** |
|---|---|---|---|
| **Table S2 — qHNMR data of the time-dependency measurements** | ★ **the numeric table behind Figure 2** — the exact mmol/L and % free for all 20 odorants × 10 time points | Would **replace this dossier's §5a digitisation (±4 pp) with printed values**, and would show whether SDs/n exist at all. **It would NOT add the temperature ladder** — Table S2 is the RT series only (*"Samples were incubated at RT"*, Fig. 2 caption; and the text cites Table S2 exclusively for RT results). | ★ **YES — highest value, but note the ceiling: it does not touch §6.** |
| **Figure S3A/B — percent free odorant in the substitution experiments** | the FFT / N-acetylcysteine / MeFFT time courses that produced the **5 %, 11 % and 8 %** covalent-incorporation figures | Would let the covalent-vs-dimer split be reconstructed as a **time course** rather than two points, sharpening §6c's 21 %-to-dimer derivation. **Still 300 K only.** | **YES — second priority.** |
| **Table S1 — aroma-model concentrations** | the 25 compounds and their levels in the reconstitution model | Needed to make §8's sensory result quantitative at the molecular level (which thiols, at what ppb) and to check whether the "roasty" reference mixture was equally weighted — the ambiguity that opens §8d's exponent bracket from 0.28 to 0.92. | **YES — third priority, and it is the one that would tighten §8d.** |
| **Figure S2 — overview of interaction scenarios I–IV** | a schematic | Already fully transcribed from the body text in §3a. **Adds nothing.** | **NO.** |
| **Figure S1 — UF purification monitoring** | stacked ¹H NMR of the wash steps | Method QC only. | **NO.** |

### 12c. ★ **THE GAP THAT NO SI CLOSES, AND THAT THE WAVE SHOULD DECLARE**

> **There is no temperature ladder in the SI. Figure 4 IS the temperature data, in its entirety, and
> it is a figure. The 279 K and 333 K series exist nowhere as printed numbers except for the five
> values on p. 15341–42 (24.7, 24.1, 23.6, 55.0, 51.8, 70.0, 49.0). Everything else about the
> temperature dependence of this system — the whole of §6 — rests on this wave's digitisation of
> Figure 4, at ±1.5 pp, cross-validated against Figure 2 to 1.0 pp and against five printed anchors
> to ≤ 0.7 pp.** That is a sound basis for the **structural** and **refutational** conclusions of §7,
> and it is **not** a sound basis for a fitted parameter. **Both statements should travel together.**

### 12d. What the sulfur lane should ask next
1. **A thiol–melanoidin binding study at 100–180 °C.** Nothing in the corpus covers the gap between
   Gigl's 333 K ceiling and the repo's 373 K floor. Until something does, the 60 kJ/mol is a
   *refutation* of 248, not a replacement for it.
2. **A melanoidin binding capacity in µmol thiol per gram.** Without it, no concentration-dependent
   sequestration term can be written.
3. **The Hofmann & Schieberle 2002 (ref. 2) and Müller & Hofmann 2007 (ref. 7) precursors** — cited
   here as the source of the *"high thiol binding affinity"* claim and of the phenol/FFT conjugate
   quantification. If either carries a temperature series, it would break the three-point tie in
   §6d. **Neither is in `data/articles/`.**
