# Poisson, Schmalzried, Davidek, Blank & Kerler 2009 — EXTRACTION (green Arabica coffee: water-extracted beans reconstituted with the natural extract or a biomimetic recombinate, with omission of free sugars / free amino acids, spiking of sugars, alanine, cysteine, and isotope-labelled [U-¹³C5]-arabinose, [U-¹³C6-fructose]-sucrose, L-[3-¹³C]-alanine and [²H4]-furfural; fluidised-bed roast 236 °C, 380 s; SIDA/SPME-GC-MS of furfural, 2-furfurylthiol and ten other odorants, mg/kg dry matter)

### THE IN-BEAN TEST OF "FFT = FURFURAL + H2S": in roasted coffee, furfural and 2-furfurylthiol move in OPPOSITE directions under every manipulation (sugar omission: furfural −66 %, FFT +62 %; sucrose spike: furfural +61 %, FFT −32 %), and labelled furfural gave **0 % labelled FFT**. The authors conclude the furfural route is "only a minor pathway" in coffee. **No route share is quantified**, and the decisive [²H4]-furfural test is weaker than it looks: the labelled furfural itself was not found in the roasted beans either.

**Source on disk:** `data/articles/poisson2009.pdf` (9 pp., J. Agric. Food Chem. 2009, 57, 9923-9931).
Read from the `pdftotext -layout` text layer; **Tables 5, 6 and 7 were checked cell by cell against a
200 dpi render of p. 9927** (scratchpad `fork_nf/poi-5.png`). Pyrazine/diketone sections (Table 8,
pp. 9928-9929) were read for context only. No SI.

## 0. Identity

| field | value |
|---|---|
| Title | "Study on the Role of Precursors in Coffee Flavor Formation Using In-Bean Experiments" |
| Authors | Luigi Poisson, Frank Schmalzried (Univ. Hohenheim), Tomas Davidek, Imre Blank, **Josef Kerler** (corresponding) — Nestlé Product Technology Centre Orbe |
| Venue | J. Agric. Food Chem. 2009, **57** (21), 9923-9931. Received 20 May 2009, revised 23 July, accepted 10 Sept; published on web 9 Oct 2009 |
| DOI as printed | **10.1021/jf901683v** ✔ matches the brief. **Identity verified.** |
| Sample codes | RGB reference green beans; WEB water-exhausted beans; NREB WEB + natural extract; BREB WEB + biomimetic recombinate; FFT 2-furfurylthiol; 2-FA 2-furaldehyde (furfural) |
| Companion | `poisson2019_extraction.md` (same group, later in-bean work) |

## 1. Why it matters

Gap (e). In the engine FFT is **95-99 % furfural + H2S** (`r_fur_fft`, `r_fur_fft_hs`;
`SULFUR_FLUX_MAP.md` §5), and FFT/MFT is 6-60× too high in Hofmann pH 5 and Yiltirak, 370× too high
in Zhang 2025b. Poisson et al. is the main literature statement that, **in a real roasting matrix**,
FFT does **not** track furfural: the sulfur and the furan skeleton come from the water-insoluble bean
fraction (bound cysteine; polysaccharide), and an alternative intermediate, furfuryl alcohol, is
"> 10×" more abundant than furfural in roasted coffee (own data, not shown).

## 2. Methods as they matter to a model

- **Beans.** Colombian Arabica, wet-processed. Hot-water extraction 4 × 280 L at 95 °C for 4 h total
  (70 kg beans); soluble solids 19.4 % (21.7 % dm). Free amino acids 96.8 % extracted, sucrose 82.9 %,
  total sugars 17.7 %, total amino acids 4.7 % (Table 3) — i.e. **WEB keeps 95 % of total amino acids
  (protein-bound, incl. bound cysteine) and 82 % of total sugars (polysaccharides)**.
- **Reconstitution.** 33.95 g dried extract (or the Table 1 recombinate, pH adjusted to **5.5** with
  KOH) into 125 g WEB, soaked ≥ 5 h at 50 °C; spikes dissolved in 60 g water, 150 g RGB soaked 2 h at
  50 °C + 2 h RT. Spikes: sugars 0.025 mol (~4 g) or sucrose 4 g; **L-cysteine 0.45 g / 150 g beans**
  (3.0 g/kg = 24.8 mmol/kg); L-alanine 0.48 g / 150 g; [U-¹³C5]-arabinose 0.9 g; [U-¹³C6-fructose]-sucrose
  2 g; **[²H4]-furfural 30 mg / 150 g** (200 mg/kg). All beans then **frozen at −80 °C and freeze-dried**
  to 10 ± 0.5 % moisture.
- **Roast.** Neuhaus "Signum" rotating fluidised bed, **380 s at 236 °C**, same for all samples;
  colour CTN 60-144 (Table 2) — so roast degree varies between samples at fixed time/temperature.
- **Analysis.** R&G 5 g (FFT etc.) in 100 mL boiling water; **cysteine 25 mg/g added before
  quantifying FFT** "to suppress reactions of FFT with coffee matrix"; labelled internal standards
  ([²H2]-FFT, [²H4]-furfural etc.); SPME DVB/CAR/PDMS 40 °C, 10 min; ZB-Wax; SIM + full scan;
  duplicates, RSD < 20 % (exceptions noted). Labelling results "corrected for the ¹³C content of the
  natural isotope"; corrected values < 0.5 % set to 0.

## 3. Tables re-typed (furan / thiol rows; mg/kg dry matter)

### Table 5 (verbatim rows). Roasted RGB vs WEB, NREB, BREB and omission recombinates

| no. | compound | (1) RGB | (2) WEB | (3) NREB | (4) BREB | (5) BREB no amino acids | (6) BREB no free sugars |
|---|---|---|---|---|---|---|---|
| 1 | 2-furaldehyde | 165.8 | 133.8 | 146.8 | 112.7 | 145.9 | **38.3** |
| 2 | 2-furfurylthiol (FFT) | 1.15 | 1.85 | 1.05 | 1.07 | 1.12 | **1.73** |
| 3 | methanethiol | 1.90 | 3.54 | 2.50 | 1.23 | 1.77 | 3.41 |
| 4 | dimethyl sulfide | 1.56 | 0.36 | 0.67 | 0.12 | 0.24 | 0.16 |
| 6 | 2,3-butanedione | 79.0 | 47.2 | 82.9 | 69.5 | 65.6 | 30.6 |
| 7 | 2,3-pentanedione | 46.2 | 7.88 | 48.0 | 35.0 | 36.1 | 3.26 |
| 8 | HDMF (furaneol) | 230.9 | 149.3 | 237.0 | 232.3 | 375.6 | 105.0 |

### Table 6 (verbatim rows). RGB spiked with precursors

| compound | (7) +SUC | (8) +GLU | (9) +ARA | (10) +RHA | (11) +ALA | (12) +CYS |
|---|---|---|---|---|---|---|
| 2-furaldehyde | **266.9** | 196.4 | 156.8 | 173.6 | 152.2 | 163.1 |
| 2-furfurylthiol | **0.78** | 1.08 | 1.18 | 1.12 | 1.14 | **1.71** |
| methanethiol | 1.84 | 2.54 | 2.57 | 2.66 | 2.25 | 2.89 |
| 2,3-pentanedione | 74.7 | 56.5 | 44.7 | 54.0 | 54.2 | 86.7 |

### Table 7 (verbatim). Labelling distribution (%) of furfural and FFT

| ion | furfural m/z | BREB^b + [¹³C6]-SUC | BREB^b + [¹³C5]-ARA | RGB + [²H4]-2-FA | FFT m/z | BREB^b + [¹³C6]-SUC | BREB^b + [¹³C5]-ARA | RGB + [²H4]-2-FA |
|---|---|---|---|---|---|---|---|---|
| M | 96 | 95 | 98 | 97 | 114 | 95 | 89 | 98 |
| M+1 | 97 | 1 | 1 | 3 | 115 | 0 | 3 | 0 |
| M+2 | 98 | 0 | 1 | 0 | 116 | 5 | 4 | 2 |
| M+3 | 99 | 0 | 0 | 0 | 117 | 0 | 4 | 0 |
| M+4 | 100 | 0 | 0 | **0** | 118 | 0 | 0 | **0** |
| M+5 | 101 | 4 | 0 | 0 | 119 | 0 | **0** | 0 |

^b recombinate omitted in free sugars.

### Arithmetic (mine)

1. **Furfural vs FFT, every manipulation** (relative to the paper's own reference):
   | manipulation | furfural | FFT |
   |---|---|---|
   | WEB vs RGB | −19 % | **+61 %** |
   | no free sugars vs BREB | **−66 %** | **+62 %** |
   | no amino acids vs BREB | +29 % | +5 % |
   | +sucrose vs RGB | **+61 %** | **−32 %** |
   | +cysteine vs RGB | −2 % | **+49 %** |
   Correlation is negative or nil in all five.
2. **FFT/furfural molar ratio in RGB:** 1.15/114.2 ÷ 165.8/96.1 = 10.1 µmol/kg ÷ 1.73 mmol/kg =
   **0.58 mol %**.
3. **Cysteine spike efficiency:** +0.56 mg/kg dm FFT (4.9 µmol/kg) from 24.8 mmol/kg cysteine added →
   **≈ 0.02 mol %** of the added cysteine ends as extra FFT (bean wet/dry basis difference ~10 % ignored).
4. **Was the [²H4]-furfural there at roast?** 200 mg/kg of [²H4]-furfural against 165.8 mg/kg of
   roast-generated furfural would show as tens of % M+4 in furfural if it survived; Table 7 shows
   **0 %** at m/z 100. So the spike was **lost before or during roasting** (furfural is volatile and the
   beans were freeze-dried after soaking), or exchanged its ring D. Therefore "0 % [²H4]-FFT" means
   only that the spike did not reach FFT; it does not by itself show that furfural formed *in situ*
   cannot. **The paper does not discuss this.**
5. **Free [¹³C5]-arabinose** (0.9 g vs 0.011 g natural free arabinose in the recombinate, i.e. the free
   pool almost fully labelled) gave **0 % M+5 in furfural and 0 % M+5 in FFT**: free arabinose is not
   the precursor of *either* compound in the bean. This argues that both come from the insoluble
   fraction — not, by itself, that FFT bypasses furfural. FFT's partial labels (M+1 3 %, M+2 4 %,
   M+3 4 %) indicate fragment recombination from the arabinose.
6. **[¹³C6-fructose]-sucrose:** furfural 4 % M+5 (intact C5 from fructose); FFT 5 % M+2, 0 % M+5.

## 4. Numbers the repository can use

| observable | value | conditions | source | evidence class |
|---|---|---|---|---|
| FFT in roasted coffee | 1.15 mg/kg dm (RGB) | Arabica, 236 °C, 380 s, fluidised bed | Table 5 | measured (SIDA, duplicate) |
| furfural in roasted coffee | 165.8 mg/kg dm (RGB) | same | Table 5 | measured |
| FFT/furfural | 0.58 mol % | RGB | arithmetic 2 | derived |
| furfural vs FFT response to sugars | anti-correlated (table in arithmetic 1) | in-bean | Tables 5, 6 | within-study ratios |
| cysteine → extra FFT | ~0.02 mol % of added Cys | +3 g/kg Cys | Table 6 | derived |
| intact C5 from free arabinose into FFT | 0 % (< 0.5 %) | [U-¹³C5]-Ara in recombinate | Table 7 | measured (labelling) |
| [²H4]-furfural → [²H4]-FFT | 0 % — **spike itself not recovered** | 200 mg/kg spike | Table 7 | measured but **inconclusive** (arithmetic 4) |
| **share of FFT via furfural** | **NOT QUANTIFIED** ("only a minor pathway") | — | p. 9926 | qualitative |

## 5. What this can and cannot anchor

**Can:** a qualitative, matrix-level statement that, in a roasting bean, FFT formation is not
rate-limited by the free-furfural pool (anti-correlation under five manipulations), and that bound
(insoluble) cysteine and non-extractable carbohydrate are the main precursors.

**Cannot:**
- Give a route share, a rate, or anything in aqueous solution. Dry roast at 236 °C with water
  ~10 % → steam; pressure and pH unknown; colour varies (CTN 60-144) at fixed time/T.
- Be scored by the engine: no engine pot represents a bean.
- Settle the furfural question cleanly (arithmetic 4: the labelled furfural was not recovered).

## 6. Flags

1. **[²H4]-furfural spike not recovered** (Table 7, m/z 100 = 0 %) — the key "ruled out" experiment is
   inconclusive as reported.
2. **Roast degree differs** between samples (Table 2 CTN 60-144); sample 2 (WEB) is much lighter (119)
   and sample 6 (no sugars) the lightest (144). Colour partly reflects the sugar charge itself, so the
   omission arms are not at equal reaction extent; that confounds the omission comparisons (direction
   of the bias on FFT vs furfural not determinable from the paper).
3. **Cysteine added before FFT quantitation** (25 mg/g R&G) to protect FFT — standard for coffee SIDA,
   but means FFT values are "recoverable FFT".
4. Furfuryl alcohol is proposed as the alternative intermediate with "> 10×" furfural levels — from
   unpublished data; FFT yield from FF-OH ">10 lower" than from furfural (ref. 37, Parliment & Stahl 1995).
5. Duplicates only; RSD < 20 %.

## 7. Use in the model

- **Gap (e), FFT route shares.** The engine's assumption (FFT ≈ 95-99 % furfural + H2S) is **contradicted
  in direction for coffee roasting**, but only qualitatively and in a matrix the engine does not model.
  It is consistent with — and the only literature support for — the flux map's diagnosis that the
  FFT/MFT mis-weighting sits in the furfural branch: if part of FFT in real systems comes from a
  route other than free furfural + H2S (furfuryl alcohol + H2S; a bound/S-glycoside precursor), then a
  furfural-only FFT route that is fitted to FFT endpoints will over-weight furfural's FFT branch
  wherever furfural is abundant (low pH, low T, long times — exactly where FFT is most over-predicted
  in Yiltirak).
- **Not usable as a rate or hold-out.** Use as a declared caveat in any pre-registration that splits
  furfural's decay barrier (§5 of the flux map): "in a roasting matrix the furfural route is reported
  minor (Poisson 2009), unquantified; this wave does not test it". A furfuryl-alcohol + H2S → FFT
  branch would need its own measured constant (Parliment & Stahl 1995, ref. 37, not on disk).
- **Caveats:** coffee bean, 236 °C dry roast, roast degree confounded, labelled-furfural spike lost.
