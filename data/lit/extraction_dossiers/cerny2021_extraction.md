# Cerny, Schlichtherle-Cerny, Gibe & Yuan 2021 — EXTRACTION (in-bean: washed Arabica green beans soaked with furfuryl alcohol (FFA) or d2-FFA 1-3 mg/g, freeze-dried, roasted 150 → 200 °C in 5 min (IKAWA sample roaster); espresso brew; FFA and 2-furfurylthiol (FFT) by HS-SPME-GC-MS SIDA; isotopologue ratios of FFT and furfural)

### A DRY-ROAST, IN-BEAN LABELLING STUDY, NOT AN AQUEOUS ONE. d2-FFA → d2-FFT is shown directly (label on the CH2); furfural formed in the same beans is **entirely unlabelled**, so FFA is not oxidised to furfural. **Furfural → FFT was NOT tested here** (the "furfural is not a precursor" conclusion rests on Poisson et al. 2009's d4-furfural in-bean result, cited). Spiked FFA accounts for **at most ~1/3-1/2 of the FFT** even when it is 2/3 of the FFA pool, and the molar conversion is tiny: **≈ 0.1 mol % of the d2-FFA lost during roasting, 0.025 % of the d2-FFA present before roasting, 0.004 % of what was spiked** (my arithmetic). Mechanism: proposed, not shown — FFA → furfuryl cation → protein-bound S-furfuryl-cysteine → β-elimination to FFT + protein dehydroalanine.

**Source on disk:** `data/articles/cerny2021.pdf` (5 pp., Elsevier publisher PDF, clean text layer). Prose
and Tables 1-3 from `pdftotext -layout`; Fig. 2 read off a 170 dpi render (scratchpad `lit5/cer-3.png`)
and tagged **[D]**; Figs. 1 (mass spectra) and 3 (scheme) qualitative. No SI.

## 0. Identity

| field | value |
|---|---|
| Title | "Furfuryl alcohol is a precursor for furfurylthiol in coffee" |
| Authors | **Christoph Cerny** (corresponding, christoph.cerny@firmenich.com), Hedwig Schlichtherle-Cerny, Romelo Gibe, Yuan Yuan — Corporate Research, Firmenich Aromatics China, Shanghai |
| Venue | Food Chemistry **337** (2021) 128008. Received 8 May 2020, revised 1 Sep 2020, accepted 1 Sep 2020, online 6 Sep 2020 |
| DOI as printed | **10.1016/j.foodchem.2020.128008** ✔ matches the brief. **Identity verified.** |
| Same author, in corpus | `cerny2013_extraction.md` (S-furfuryl-L-cysteine and S-(2-methyl-3-furyl)-L-cysteine from Cys + Xyl and from FFA + Cys, 100 °C reflux, identification only); `cerny2015_extraction.md` (review chapter); `cerny2003`, `cerny2004`, `cerny2007`, `cerny2007b`, `cerny2008` |

## 1. Why it matters

The engine makes FFT 95-99 % from furfural + H2S (`docs/guides/SULFUR_FLUX_MAP.md` sec. 2, 5) and has no
furfuryl-alcohol species. `docs/guides/LITERATURE_SEARCH_2026-10_PATHS.md` (Poisson 2009 row) flagged
that, in coffee, free furfural is not the FFT precursor. This paper is the positive half of that claim:
FFA is a precursor in the bean. The question for the model is whether anything transfers to buffered
aqueous pots.

## 2. Methods as they matter to a model

- **Matrix:** washed Arabica green beans, Yunnan. Green-coffee protein: 0.2 g/kg cysteine, 4.7 g/kg cystine
  (Meichelbeck & Zahn 1968, cited).
- **Spiking:** 50 g beans soaked in 20 g water + FFA or d2-FFA (**50, 100, 150 mg → 1, 2, 3 mg/g**) at 50 °C,
  4 h (rotavapor, atmospheric), + 8 g water, 1 h more; freeze-dried; −20 °C.
- **Retention before roasting:** beans spiked at 2 mg/g held **342 µg/g FFA** (Table 1, 0 min) — **17 %
  of the spike**; "the majority of added FFA was lost during drying, possibly reacting with coffee solids
  rather than by evaporation".
- **Roast:** 17-18 g, IKAWA sample roaster, **150 → 200 °C within 5 min** (airflow 70-60 %), rapid cool to
  50 °C. Time series stopped at 1, 2, 3, 4, 5 min = **final air temperatures 160, 170, 180, 190, 200 °C**.
  **Bean temperature, moisture, pH not measured.** Roast-profile temperatures are the roaster set-points.
- **Brew:** 7 g ground (Baratza Encore, setting 10), DeLonghi EC 680 pump espresso, 30-38 mL, iced.
- **Analysis:** HS-SPME (75 µm DVB/CAR/PDMS, 4 mL brew, 10 + 30 min at 40 °C), DB-WAX GC-MS (EI, scan),
  SIDA with d2-FFA (ions m/z 81/83) and d2-FFT (m/z 114/116); FFA brew diluted 12.5×; **200 mg L-cysteine
  added to the brew before FFT analysis "to prevent FFT binding"**. Concentrations **per g green coffee**.
  Duplicate analyses; one or two roasting trials (Table 1 footnotes).

## 3. Tables (verbatim) and figure reads

### Table 1 — FFA and FFT in brew, per g green coffee (mean ± deviation)

| added FFA (mg/g) | roast time (min) | FFA (µg/g) | FFT (ng/g) |
|---|---|---|---|
| untreated ᵇ | 5 | 223 ± 14 | 246 ± 9 |
| 1 ᵇ | 5 | 321 ± 42 | 291 ± 20 |
| 3 ᵇ | 5 | 468 ± 79 | 390 ± 23 |
| 2 | 0 | 342 ± 1 | n.d. |
| 2 | 1 | 322 ± 8 | n.d. |
| 2 | 2 | 303 ± 4 | n.d. |
| 2 | 3 | 323 ± 24 | 21 ± 0 |
| 2 | 4 | 311 ± 1 | 70 ± 2 |
| 2 | 5 | 380 ± 5 | 296 ± 15 |

ᵇ two roasting trials × two brews; others one trial × two brews.

### Table 2 — FFT isotopologues (m/z), 5 min / 200 °C, d2-FFA spiked

| m/z | untreated | 1 mg/g | 2 mg/g | 3 mg/g |
|---|---|---|---|---|
| 114 | 89 % | 70 % | 60 % | 51 % |
| 115 | 6 % | 6 % | 6 % | 6 % |
| 116 | 5 % | 21 % | 31 % | 38 % |
| 117 | — | 2 % | 2 % | 3 % |
| 118 | — | 1 % | 1 % | 2 % |

### Table 3 — d2 share vs roast time (2 mg/g d2-FFA)

| min (final °C) | 0 | 1 (160) | 2 (170) | 3 (180) | 4 (190) | 5 (200) |
|---|---|---|---|---|---|---|
| d2-FFA, % of FFA | > 98 | > 98 | > 98 | 93 | 82 | 67 |
| d2-FFT, % of FFT | n.d. | n.d. | n.d. | **46** | **38** | **34** |

### Fig. 2 [D] (read error ± 5 ng/g or ± 5 µg/g; the authors computed these from Tables 1-3)

- 2A, d2-FFT / unlabelled FFT (ng/g) at 1 / 2 / 3 mg/g: **60 / 230, 95 / 200, 163 / 227**.
- 2B, FFA (µg/g), 0-5 min: total 342, 322, 303, 323, 311, 380; d2-FFA ≈ 342, 322, 300, 300, 255, 255;
  unlabelled ≈ 0, 0, 0, 22, 56, **125**.
- 2C, FFT (ng/g), 0-5 min: total 0, 0, 0, 21, 70, 296; unlabelled ≈ 0, 0, 0, 11, 43, **195**; d2 ≈ 0, 0, 0,
  10, 27, **101**.

Fig. 1B (furfural mass spectrum, 3 mg/g d2-FFA, 5 min): molecular ions m/z 95/96 only — **no m/z 97/98,
"exclusively unlabelled furfural"**.

## 4. Arithmetic (mine)

1. **Molar conversion d2-FFA → d2-FFT (2 mg/g arm, 5 min).** d2-FFT 101 ng/g ÷ 116.2 = **0.87 nmol/g**.
   d2-FFA spiked 2 mg/g ÷ 100.1 = 20 µmol/g → **0.004 mol %**; present before roasting 342 µg/g = 3.4 µmol/g
   → **0.025 mol %**; lost during roasting (342 → 255 µg/g, 24 %, partly evaporation) 0.87 µmol/g → **≈ 0.1
   mol %**. At 3 mg/g: 163 ng/g = 1.40 nmol/g vs 30 µmol/g spiked → 0.005 %.
2. **Route share ceiling.** At 5 min the FFA pool is 67 % d2, time-averaged over 2-5 min ≈ 80-90 %, yet
   d2-FFT is 34-46 %. If all FFT came from free FFA with one rate law, the d2 share of FFT would track the
   d2 share of FFA. It is about half that, so **≥ ~50 % of FFT comes from something other than the free FFA
   pool** (native FFA in a different compartment, or a non-FFA route). The labelled experiment gives a
   *lower bound* on the FFA route only for spiked beans; **the FFA-route share in unspiked coffee is not
   determined.**
3. **Competition for the sulfur donor.** Untreated beans give 246 ng/g FFT; 2 mg/g spiked beans give 296
   total but only ~195 unlabelled. Spiking **displaced ~50 ng/g of native FFT** while adding ~100 ng/g d2.
   Consistent with a shared limiting sulfur donor (protein Cys), which is my inference, not the authors'.
4. **Natural M+2 not subtracted.** Unlabelled beans show 5 % at m/z 116 (³⁴S, ¹³C₂), which the authors
   count as d2 in Fig. 2A → d2-FFT overstated by ~5 points (≈ 15 ng/g; ~25 % relative at 1 mg/g).
5. **Apparent temperature law (for the record, NOT usable):** FFT increments 21 / 49 / 226 ng/g over
   min 2-3 / 3-4 / 4-5 (mean air T 175 / 185 / 195 °C) give "Ea" of 145 and 273 kJ/mol between successive
   intervals — confounded by precursor build-up, bean-vs-air temperature lag and drying.

## 5. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| FFA → FFT (in bean) | molar yield | 0.004 (of spike) / 0.025 (of pre-roast pool) / ~0.1 (of FFA lost) | mol % | green coffee, dry roast 150 → 200 °C, 5 min | Tables 1, 3; arithmetic 1 | derived; matrix-specific |
| FFT from spiked FFA | share of FFT | 34-46 | % | 2 mg/g d2-FFA, 3-5 min | Table 3 | measured (isotopologue) |
| FFA → furfural | — | none detected | — | same, 3 mg/g | Fig. 1B | measured, qualitative |
| FFT onset | first detected | > 2 min (> 170 °C air) | — | same | Table 1, 3 | measured |
| furfural → FFT | — | **NOT TESTED in this paper** | — | — | cites Poisson 2009 | — |
| FFA + Cys → FFT in water | relative | FFT from FFA ≈ 1/10 of FFT from furfural | — | aqueous Cys model (conditions not given here) | Introduction, citing Parliment & Stahl 1995 | secondary, unverified |

## 6. What is NOT reported

- No aqueous experiment, no buffer, no pH, no water activity, no bean temperature.
- No rate constant, no Ea, no S-furfuryl-cysteine measurement (the intermediate is proposed by analogy to
  Cerny 2013 and Kreppenhofer 2011; protein-bound conjugate not detected).
- No furfural quantities; no furfural spiking; no H2S; no cysteine/cystine balance.
- No labelled-FFT result for native (unspiked) FFA → the FFA-route share in normal coffee is unknown.

## 7. Flags

1. **"Results unambiguously show that FFA serves as a precursor"** — yes for spiked FFA; the abstract's
   "furfural ... plays no major role as precursor" rests on another paper (Poisson 2009), not on data here.
2. The 200 mg cysteine added to the brew before FFT analysis is a sample-preparation step, not a reactant;
   it could in principle liberate FFT from S-conjugates/disulfides in the brew (it is meant to).
3. Per-gram values refer to green coffee before roasting; roast loss not stated.
4. Different numbers of roasting trials across rows of Table 1 (footnote b).

## 8. Use in the model

- **Transfer to aqueous pots: low.** Dry roasting, protein-bound cysteine, a bean matrix that retains FFA
  on solids, no water phase. The authors themselves argue the route works *because* H2S would volatilise
  in a dry roast and free Cys is scarce — conditions opposite to a sealed buffered pot.
- **Furfural → FFT in water is not challenged by this paper.** The only aqueous statement it carries
  (citing Parliment & Stahl 1995) is that cysteine makes ~10× MORE FFT from furfural than from FFA in
  water — consistent with the engine's furfural + H2S FFT route for buffer pots. **No contradiction with
  `SULFUR_FLUX_MAP.md` for aqueous pots**; the in-bean result says only that the engine's FFT route
  should not be carried over to a coffee/dry-roast matrix.
- **Mechanistic option for a later wave:** an FFA species (from furfural reduction or from pentose) + Cys
  → S-furfuryl-Cys → FFT (β-elimination) would give a thiol route that does not draw on free H2S — relevant
  to the FFT/MFT mis-weighting (flux map sec. 5) only if an aqueous rate is found. Cerny 2013 shows the
  conjugates do form in a 100 °C Cys + Xyl syrup; Liu 2021 (`liu2021`, in the literature search) is the
  candidate aqueous source. Nothing here is a number the engine can take.
