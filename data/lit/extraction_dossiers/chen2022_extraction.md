# Chen, Wang, Zhang, Yang, Du, Yu & Xie 2022 — EXTRACTION (lab pea protein isolate, 87.4 % protein, heated at 95 °C conventionally or ohmically: SH in native buffer and in 8 M urea, H0, gels)
### A pea isolate's free SH after 95 °C / 10 min is printed (2.95 native-buffer, 7.8 urea µmol/g protein); the native values are figure-only; no S-S, no amine.

**Source on disk:** `data/articles/Chen2022.pdf` (owner's download, 2026-09-08). Read from the pdftotext
layer in the scratchpad (`articles/Chen2022.txt`, 1485 lines) and re-checked page 8 with pypdf: the
SH data are a bar chart (Fig. 6a) whose text layer carries only the axis tick labels and the
significance letters, so the native and most treated values are **FIGURE-ONLY**; the four SH values
quoted in §3.5 of the text are the only numbers. Tables 1 (FTIR) and 2 (gel texture) have a clean
text layer and are re-typed. Repo status before this dossier: `pea_isolate` in `protein_matrices.yml`
carries free thiol 0.0159 (band 0.0021-0.0174) and disulfide 0.0257 (0.0042-0.0297) mmol per g
protein, native; no amine; no after-heat factor.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of moderate electric fields on the structural and gelation properties of pea protein isolate" |
| Authors | Yan Chen, Tong Wang, Yinfeng Zhang, Xinrun Yang, Jing Du, Dianyu Yu*, Fengying Xie* (School of Food Science, Northeast Agricultural University, Harbin) |
| Venue | Innovative Food Science and Emerging Technologies 77 (2022) 102959; received 11 Oct 2021, revised 20 Feb 2022, accepted 21 Feb 2022, online 23 Feb 2022 |
| DOI | 10.1016/j.ifset.2022.102959 |
| Naming | "Native" = untreated PPI dispersion; CH = conventional water-bath heating; OH = ohmic heating (the "moderate electric field", MEF) at 20 kHz or 50 Hz and 5 / 10 / 20 V/cm, same temperature programme as CH. **"Total SH" in this paper is the SH titrated in 8 M urea WITHOUT any reducing agent** (§2.9) — it is the unfolded-protein free thiol, NOT half-cystine; "free SH" is the SH titrated in Tris-Gly without denaturant (surface thiol). No disulfide was measured; every "S-S" statement in the paper is an inference from SH decrease. |
| Method parent | Beveridge, Toma & Nakai 1974 (J Food Sci 39, 49) "with some modifications" |

## 1. Why it matters

The repo's pea free-thiol density is a native, 8 M-urea Ellman's number from Gao 2020 (whole
isolate), bracketed below by Chihi 2016's globulin fraction; the pre-registration (§5) records that
heating moves these densities and that this is not modelled. This paper gives, for a laboratory
whole pea isolate of stated protein content, the urea-unfolded SH (the quantity comparable to Gao's
"free SH") and the native-buffer SH, after a defined 95 °C / ~15 min treatment of a 1 % dispersion
at pH 7 — the closest thing on disk to a "pea isolate after cooking" thiol density — plus the same
after ohmic heating. Its usable content is small: two printed SH pairs (CH and the lowest OH
condition), a native baseline that is figure-only but bounded by the text, and a within-study
OH/CH ratio. It measures no amine and no disulfide, so it leaves the amine pool unfilled.

## 2. Methods as they matter to a model

- **Peas and flour (§2.1):** yellow peas (Yantai Oriental); "48.7% total starch, 7.6% moisture,
  23.4% crude protein, 1.5% fat, and 2.3% ash on a dry basis"; Kjeldahl, **nitrogen factor 6.25**.
- **PPI (§2.2, Cui 2020 method):** flour sieved 0.180 mm; defatted with n-hexane 1:5 (w/v) three
  times; 100 g flour in water 1:15 (w/v), pH 8.5 (1 M NaOH), 2 h at 25 °C, 5000 rpm 20 min 4 °C;
  supernatant to pH 4.5 (1 M HCl), 5000 rpm 15 min; precipitate washed three times with water,
  brought to pH 7.0, freeze-dried 24 h; "a PPI with a protein content of 87.43 g/100 g" (protein
  assay for the isolate not restated; by §2.1, Kjeldahl N × 6.25). Alkaline-extraction / isoelectric
  whole isolate — the same class as Gao 2020's (pH 8.5-9.5, pI 4.5, 83-85 % protein) and Shen 2022's.
- **CH (§2.3.1):** 1 % (w/v) PPI in 0.01 M phosphate pH 7.0, stirred 2 h at 25 °C; 80 mL heated in a
  95 °C water bath, "approximately 5 min of come-up time and 10 min of holding time"; then 4 °C
  ice-water bath. Thermocouple-logged.
- **OH (§2.3.2):** same dispersion, 80 mL in a quartz cell (5.1 × 4.5 × 5.5 cm), two platinized
  titanium electrodes 4 cm apart, magnetic stirring; sinusoidal field from a signal generator +
  power amplifier; voltage adjusted to hold the same programme (~5 min come-up + 10 min at 95 °C);
  frequency 50 Hz or 20 kHz; field 5, 10, 20 V/cm; "No electrode corrosion or electrolyte gas
  occurred". Six OH conditions + CH + native = 8 samples.
- **Gels (§2.4):** 12 % (w/v) PPI + 0.2 mol/L CaCl2, 95 °C for 1 h by CH or OH, 4 °C overnight. Gel
  data (WHC, TPA, SEM, rheology) are on these gels, not on the 1 % dispersions.
- **SH (§2.9), quoted:** "The samples were diluted to 4 mg/mL with phosphate buffer (0.01 M, pH
  7.0). Then, a 1 mL sample solution was dissolved in 5 mL Tris-glycine buffer containing ...
  (4 mmol/L EDTA, 0.086 mol/L Tris, 0.09 mol/L glycine, pH 8.0), and then 50 µL Ellman's reagent was
  added. The Tris-glycine buffer with (total SH) or without 8 mol/L urea (free SH). The solution was
  oscillated for 30 min at room temperature and centrifuged at 6000 rpm for 10 min. The solution
  without Ellman reagent was used as a reagent blank. The absorbance ... at 412 nm". Equation:
  SH (µmol/g) = 73.53 × A412 × D / C, "C is the final sample concentration (g/mL)". ⚠ With C in
  g/mL the equation returns µmol/mg; the printed magnitudes (1.8-7.8 µmol/g) require C in mg/mL as
  in Beveridge 1974. D is not stated (nominally 6.05). Ellman's reagent concentration not stated.
  So: **"free SH" = surface thiol in native buffer at pH 8; "total SH" = thiol after 30 min in
  ~6.6 M urea** (5 mL of 8 M into 6.05 mL), no reduction step. The "4 mg/mL" is the diluted 1 %
  dispersion, i.e. 4 mg of PPI POWDER per mL; the paper labels the result "µmol/g protein" — either
  a division by 0.8743 was applied silently or the label is loose (both readings given in §4).
- **H0 (§2.10):** ANS (Voutsinas 1983), 0.05-0.40 mg/mL in 0.01 M phosphate pH 7.0, ex 370 / em 470
  nm, slope of intensity vs concentration. Arbitrary instrument scale.
- **Other:** DLS particle size (Mastersizer 2000, 1 mg/mL); SDS-PAGE ± β-ME; FTIR on powder (second
  derivative, amide I); intrinsic fluorescence 0.2 mg/mL, ex 290 nm.
- **Replicates (§2.15):** triplicate, mean ± SD, Duncan P < 0.05.

## 3. Tables re-typed

### Table 1. "Secondary structure content of untreated, conventional heat (CH)-treated and ohmic heat (OH)-treated pea protein isolate (PPI)" (%; band assignments α-helix 1646-1662, β-sheet 1615-1636 and 1682-1700, β-turn 1663-1681, random coil 1637-1645 cm⁻¹)

| sample | α-helix | β-sheet | β-turn | random coil |
|---|---:|---:|---:|---:|
| Native | 20.45 ± 0.15 a | 35.28 ± 0.18 g | 21.17 ± 0.14 a | 23.10 ± 0.18 f |
| CH | 19.13 ± 0.20 b | 37.11 ± 0.14 f | 19.55 ± 0.15 b | 24.21 ± 0.17 e |
| OH 20 kHz 5 V/cm | 18.76 ± 0.14 c | 37.52 ± 0.19 e | 19.24 ± 0.14 c | 24.48 ± 0.12 d |
| OH 20 kHz 10 V/cm | 18.58 ± 0.15 d | 37.86 ± 0.23 d | 19.02 ± 0.17 d | 24.54 ± 0.21 d |
| OH 20 kHz 20 V/cm | 18.24 ± 0.23 e | 38.15 ± 0.18 c | 18.71 ± 0.14 e | 24.90 ± 0.13 c |
| OH 50 Hz 5 V/cm | 18.17 ± 0.21 e | 38.23 ± 0.11 c | 18.68 ± 0.18 ef | 24.92 ± 0.22 c |
| OH 50 Hz 10 V/cm | 17.87 ± 0.18 f | 38.48 ± 0.10 b | 18.51 ± 0.16 f | 25.14 ± 0.19 b |
| OH 50 Hz 20 V/cm | 17.53 ± 0.18 g | 38.83 ± 0.14 a | 18.25 ± 0.18 g | 25.39 ± 0.15 a |

### Table 2. "Textural properties of pea protein isolate (PPI) gels induced by conventional heating (CH) and ohmic heating (OH)" (12 % PPI + 0.2 M CaCl2, 95 °C 1 h)

| sample | hardness (g) | springiness | cohesiveness | gumminess (g) | chewiness (g) |
|---|---:|---:|---:|---:|---:|
| CH | 286.59 ± 9.89 a | 0.90 ± 0.01 ab | 0.28 ± 0.01 c | 81.09 ± 0.79 a | 73.22 ± 0.78 a |
| OH 20 kHz 5 V/cm | 218.64 ± 7.40 g | 0.90 ± 0.01 ab | 0.32 ± 0.02 a | 70.55 ± 0.87 d | 63.28 ± 0.52 d |
| OH 20 kHz 10 V/cm | 224.48 ± 8.46 f | 0.88 ± 0.01 c | 0.32 ± 0.01 a | 72.60 ± 0.70 c | 63.83 ± 0.61 d |
| OH 20 kHz 20 V/cm | 243.62 ± 7.36 d | 0.88 ± 0.02 c | 0.31 ± 0.01 ab | 74.73 ± 0.74 b | 65.56 ± 0.47 c |
| OH 50 Hz 5 V/cm | 238.87 ± 8.59 e | 0.89 ± 0.00 bc | 0.31 ± 0.01 ab | 74.49 ± 1.01 b | 66.44 ± 0.53 c |
| OH 50 Hz 10 V/cm | 257.39 ± 8.94 c | 0.91 ± 0.01 a | 0.29 ± 0.02 bc | 75.21 ± 0.84 b | 68.20 ± 0.60 bc |
| OH 50 Hz 20 V/cm | 268.66 ± 9.45 b | 0.90 ± 0.01 ab | 0.29 ± 0.02 bc | 77.37 ± 0.81 b | 69.34 ± 0.65 b |

### Fig. 6a (SH) and 6b (H0): FIGURE-ONLY except the text values

Bar charts, eight bars each. Text layer holds only the tick labels (free-SH axis 1-5 µmol/g;
"total"-SH axis 2-10 µmol/g; H0 axis 0-1000) and the Duncan letters. Values printed in §3.5-3.6:

| sample | free SH, native buffer (µmol/g protein) | "total" SH, 8 M urea, no reductant (µmol/g protein) | H0 (ANS slope) |
|---|---:|---:|---:|
| Native | FIGURE-ONLY (text: significantly higher than all heated samples) | FIGURE-ONLY (same) | 377.15 |
| CH, 95 °C ~5 + 10 min | 2.95 ("highest" of the heated) | 7.8 ("highest" of the heated) | 633.65 ("PPI after heat treatment"; which treatment is not said) |
| OH 20 kHz 5 V/cm | 1.82 ("lowest") | 3.63 ("lowest") | FIGURE-ONLY |
| other five OH | FIGURE-ONLY; text: free SH rises with field strength (10 and 20 V/cm not different), 50 Hz > 20 kHz | FIGURE-ONLY | FIGURE-ONLY; 50 Hz > 20 kHz > CH; rises with field |

Other printed numbers: particle size, native bimodal peaks at 106 and 4800 nm; OH 50 Hz 20 V/cm
D43 211.67 nm, PDI 0.54 (lowest); OH 20 kHz 5 V/cm 233.51 nm. Fluorescence λmax 334 (native) →
338 (CH) → 340 nm (OH 50 Hz); intensity 2862 → 3072 a.u. from 5 to 20 V/cm at 50 Hz. Gel WHC:
CH 86.73 g/100 g; OH 50 Hz 90.12 (5 V/cm) → 89.07 (20 V/cm); all OH > CH, OH conditions not
different from each other. Fig. 7 (WHC), Fig. 8 (SEM), Fig. 9 (G', G") figure-only otherwise.

## 4. Site densities the repository can use

Basis assumption: "per g protein" as printed; if the 4 mg/mL in Eq. 2 was powder (as the method
reads), divide by 0.8743 for the per-protein value (second number in parentheses). The urea column
is the quantity comparable to the repo's `free_thiol` (Gao 2020: Ellman's in 8 M urea, no
reductant); the native-buffer column is surface thiol, a quantity the repo table does not hold.
Nothing in this paper is a disulfide or a half-cystine.

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN (arithmetic; basis) | conditions | source | evidence |
|---|---|---:|---|---|---|---|---|
| pea protein isolate, lab AE-IEP, 87.43 % protein (N × 6.25) | free SH in 8 M urea, no reductant ("total SH") | FIGURE-ONLY; **> 7.8** (text: native significantly above CH) and below the axis top of 10 | µmol/g | > 0.0078 and < 0.010 as printed (> 0.0089 and < 0.0114 if per powder) | native, 1 % dispersion pH 7, room temperature | Fig. 6a; §3.5 | figure_only, text-bounded |
| same | free SH, native buffer (surface) | FIGURE-ONLY; **> 2.95** and below the axis top of 5 | µmol/g | > 0.0030 and < 0.005 (> 0.0034 / < 0.0057) | native | Fig. 6a; §3.5 | figure_only, text-bounded |
| same, **after CH 95 °C (~5 min come-up + 10 min hold), 1 % w/v, 0.01 M phosphate pH 7.0** | free SH in 8 M urea | 7.8 (sd not printed) | µmol/g protein | **0.0078** (0.0089 per powder) | heated, then ice bath | §3.5 text | measured, level printed in text |
| same, after CH | free SH, native buffer | 2.95 (sd not printed) | µmol/g protein | **0.00295** (0.00337) | heated | §3.5 text | measured |
| same, after OH 20 kHz 5 V/cm, same temperature programme | free SH in 8 M urea | 3.63 | µmol/g protein | 0.00363 (0.00415) | ohmic, lowest of all | §3.5 text | measured |
| same, after OH 20 kHz 5 V/cm | free SH, native buffer | 1.82 | µmol/g protein | 0.00182 (0.00208) | ohmic | §3.5 text | measured |
| ratio OH(20 kHz, 5 V/cm) / CH | urea SH; native-buffer SH | 3.63/7.8 = **0.47**; 1.82/2.95 = **0.62** | — | basis-independent | same heating profile, field added | derived from §3.5 | within-study ratio |
| ratio CH / native | either SH | not recoverable (native is figure-only); < 1 by the text | — | — | | §3.5 | directional |
| same | disulfide, half-cystine | NOT MEASURED (no reductant used; "S-S" statements are inferences) | — | — | — | — | absent |
| same | amine (free NH2 / lysine) | NOT MEASURED | — | — | — | — | absent |
| same | H0 (ANS slope) | 377.15 native; 633.65 heated | arbitrary | n/a | 0.05-0.4 mg/mL, pH 7 | §3.6 text | measured, instrument scale |
| same | protein content | 87.43 | g/100 g | f_p = 0.8743 (N × 6.25) | freeze-dried | §2.2 | stated |
| same | pea flour composition | 23.4 % protein, 48.7 % starch, 1.5 % fat, 2.3 % ash, 7.6 % moisture (db) | % | | | §2.1 | stated |

Comparison with the repo table (native pea isolate, mmol per g protein): repo `free_thiol` centre
0.0159 (band 0.0021-0.0174). This paper's native urea-SH lies between 0.0078 and 0.010 as printed
(0.009-0.011 per powder): **inside the band, 0.5-0.7x the centre**, closer to the book chapter's
PPI (0.0063, `xiao2024_extraction.md`) than to Gao 2020's 0.0147-0.0174, and far above Chihi's
globulin fraction (0.0021). After a 95 °C / 10 min hold at 1 % the urea-SH is 0.0078, i.e. still at
or above the repo band's lower half — heating this whole isolate does not collapse its thiol to the
globulin value, it trims it. The native-buffer (surface) thiol is roughly a third of the urea
value both before and after heating.

## 5. Flags

1. **"Total SH" is not total.** No reductant (β-ME, DTT) and no TCA step were used; the "total"
   column is Ellman's after 30 min in ~6.6 M urea. Consequently the paper contains **no disulfide
   measurement**, and its repeated statements that SH was "oxidized to S–S bonds" are inferences.
   Do not enter any S-S number from this paper.
2. **The native baseline is figure-only.** Only CH (7.8 / 2.95) and the lowest OH condition
   (3.63 / 1.82) are printed. The native value is bounded below by the CH value (text says every
   heated sample is significantly lower) and above by the axis range in the figure's text layer
   (5 and 10 µmol/g); that is a bracket, not a value.
3. **Basis of "per g protein" unstated**; the method dilutes a 1 % w/v POWDER dispersion to 4 mg/mL,
   and Eq. 2 has a unit typo ("C ... (g/mL)"). Both readings given; they differ by 14 %.
4. **No amine, no lysine, no TNBS / OPA.** Nothing here for the amine pool; the Introduction cites
   Li 2018 for OH raising "free amino groups" of a legume protein, not measured here.
5. **Only one heating condition, and it is mild for a Maillard model** (95 °C, ~10 min hold, 1 %
   protein, pH 7, dilute phosphate, no sugar). The SH loss on heating a whole pea isolate here
   runs opposite in sign to Chihi 2016's pea GLOBULINS at 85 °C / 60 min (which gained free thiol
   and lost S-S); different fraction, buffer and assay — not a contradiction the repo needs to
   resolve, but the after-heat factor for pea is not a single number.
6. **Ohmic heating is a confound the repo does not model.** The 0.47-0.62 OH/CH ratio at 20 kHz,
   5 V/cm is attributed to field-accelerated SH oxidation; the text says free SH is higher at 50 Hz
   and at higher field, so the six OH conditions bracket the CH value from below only (all six are
   below CH by the text). No OH number should stand in for "heated pea isolate".
7. **SDs of the four text values are not printed** (triplicates were run; the bars in Fig. 6a carry
   error bars that are figure-only).
8. **H0 "after heat treatment (633.65)"** does not say which treatment; likely CH (the sentence
   contrasts heated vs untreated before introducing OH). Arbitrary ANS slope; only the 1.68x ratio
   is transferable, and only within this instrument.
9. Gel measurements (12 % PPI + 0.2 M CaCl2, 95 °C 1 h) are a different system from the 1 %
   dispersions on which SH and H0 were measured; do not pair them.
10. The paper's 87.43 % protein was obtained with N × 6.25 (by the flour sentence); the book's
    Chapter 9 lists pea factors of 5.40-5.44, which would put the true protein nearer 76 % — the
    same caveat the repo already carries for Gao 2020.
