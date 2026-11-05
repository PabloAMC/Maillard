# Gupta, Mikhaylenko, Balasubramaniam & Tang 2011 (LWT 44:2141) — EXTRACTION (file stem `gupta2011`)
### *"Combined pressure–temperature effects on the chemical marker (4-hydroxy-5-methyl-3(2H)-furanone) formation in whey protein gels"*
### Whey protein gel (WPC 24 + WPI 11 + **D-ribose 1** + NaCl 0.8 + water 63.2 g/100 g; ~29 % protein, ~2 % lactose), **paste pH 6.1 (natural) or 8.3 (acetate + bicarbonate solutions)**, gel set 80 °C / 40 min; vacuum-packed ~3.5 g slices treated **0.1 MPa / 105 °C (thermal)**, 350 or 700 MPa at 105 °C (PATP) or 30 °C (HPP), **0 / 5 / 10 / 20 min**; M-2 = norfuraneol by anion-exclusion HPLC-DAD 285 nm in **mg / 100 g gel, calibrated on 4-hydroxy-2,5-DIMETHYL-3(2H)-furanone (furaneol), not on norfuraneol**. Formation side: k(105 °C, 0.1 MPa) = **0.0772 min⁻¹ at pH 6.1, 0.056 at pH 8.3**; at 20 min NF is **4.5× lower at pH 8.3** (14.9 vs 67.0 mg/100 g).

**Read 2026-11-05 (first reader) from `data/articles/gupta2011.pdf`** (6 pp., Elsevier publisher PDF, clean text layer).
Prose and Tables 1-3 from `pdftotext -raw`, Table 2 checked on a 120 dpi page render (the minus signs are lost in the text
layer). **Figs. 2 and 3 digitised from 300 dpi crops of p. 2144** (scratchpad `img/gupta_f2-4.png`, `img/gupta_f3-4.png`):
Fig. 2 6.75 px per mg/100 g, read error ±0.5; Fig. 3 28.4 px per mg/100 g, read error ±0.15.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Rockendra Gupta (OSU), Galina Mikhaylenko (WSU), **V.M. Balasubramaniam*** (OSU, corresponding), Juming Tang (WSU) |
| Venue | *LWT – Food Science and Technology* **44** (2011) 2141-2146 |
| DOI | **10.1016/j.lwt.2011.05.007** (printed) — matches the brief |
| Dates | received 16 Sep 2010; revised 10 May 2011; accepted 11 May 2011 |

**Laboratory.** Ohio State (Balasubramaniam, high-pressure processing) with Washington State (Tang, Mikhaylenko — the
M-2 marker group of `lau2003` and `pandit2006`; Mikhaylenko is a co-author of Pandit 2006). The gel recipe is Lau's,
modified; the HPLC method is the Tang lab's. **For leave-one-lab-out, place with the WSU Tang marker papers**
(`aroma_other.yml` uses the lab string "Washington State (Tang)"; `laboratories.yml` has no group for it yet). P ≈ 0.6 that
the M-2 analyses were run at WSU [Z, from the author list only].

## 1. What this paper answers for the model

| model need | what Gupta 2011 gives | tag |
|---|---|---|
| NF formation vs time, ribose + protein, 105 °C | **pH 6.1: 9.6 → 35.6 → 47.3 → 67.0 mg/100 g at 0 / 5 / 10 / 20 min**; pH 8.3: 10.1 → 11.4 → 12.7 → 14.9 | [F] |
| first-order formation k | **0.0772 min⁻¹ (pH 6.1, R² 0.99)**, **0.056 min⁻¹ (pH 8.3, R² 0.99)** at 105 °C, 0.1 MPa | [M] |
| NF level vs pH | 4.5× lower at pH 8.3 at 20 min; the authors attribute it to other by-products at alkaline pH [Q] | [F] / [Q] |
| NF decay | **not separable**: no fed-NF arm. Pressure arms at pH 8.3 rise then fall (PATP 350 MPa 105 °C: 11.5 → 15.2 → 13.0 → 8.0) | [F] |
| pH of a Lau-type whey gel | **paste pH 6.1, measured** (natural; before gelation and heating). The first measured pH for a gel of this family | [M] |
| O2 | slices vacuum-packed in polypropylene pouches (low O2 [Z]) | [M] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Gel | per 100 g: WPC 24 g (78.4 % protein, 8.41 % lactose), WPI 11 g (93.2 % protein), **D-ribose 1 g**, NaCl 0.8 g, water 63.2 g; protein ~29 g, lactose ~2 g |
| Ribose | 1 g / 100 g = **66.6 mmol/kg gel ≈ 105 mmol/kg water** [Z] |
| pH | "natural pH of the paste was 6.1" [V]; pH 8.3 by adding 10 g/100 mL sodium acetate ("pH 7.0") and 5 g/100 mL potassium bicarbonate ("pH 10.0") solutions; **measured on the paste before gelation; not after heating** |
| Gel set | 200 g paste in a foil-covered Nalgene bottle, 80 °C / 40 min water bath, ice quench, 4 °C |
| Samples | 7 × 2 × 0.1 cm slices (~3.5 g), vacuum-packed (polypropylene pouches) |
| Thermal (0.1 MPa) | pouches immersed in propylene glycol at **105 °C**, 0 / 5 / 10 / 20 min, ice-water quench; come-up not stated (a 1 mm slice heats in seconds [Z]) |
| PATP / HPP | 54 mL pressure vessel; target 105 °C or 30 °C at 350 / 700 MPa; Table 1 temperature history (pre-heat 84 or 57 °C, 326 s; come-up 35 s; holding 104.7-105.5 °C) |
| Replicates | n = 3, mean ± SD |

## 3. Assay [M]

1.2 g gel + 10 mL 10 mM H2SO4, ground, centrifuged, 0.45 µm; Bio-Rad fast acid analysis column, 10 mM H2SO4 1 mL/min,
DAD 285 nm; M-2 elutes at 5-8 min. **Standard: "4-hydroxy-2,5-dimethyl-3(2H)-furanone" (Sigma, purum 99 %)** [V,
abridged] — i.e. **furaneol (HDMF), not norfuraneol**, although the text names it as the M-2 standard. Concentrations are
"mg M-2 / 100 g gel". **No LOD / LOQ.** Unprocessed gel: **9.17 mg/100 g (pH 6.1), 6.1 mg/100 g (pH 8.3)** [M].

## 4. Data

### 4.1 Fig. 2, pH 6.1 (mg M-2 / 100 g gel) [F ±0.5; error bars SD, n = 3]

| treatment | 0 min | 5 min | 10 min | 20 min |
|---|---|---|---|---|
| **Thermal 0.1 MPa, 105 °C** | 9.6 | 35.6 (±0.8) | 47.3 (±2.6) | **67.0** (±1.8) |
| PATP 350 MPa, 105 °C | 8.4 | 20.0 | 25.0 | 36.6 |
| PATP 700 MPa, 105 °C | 5.6 | 10.0 | 10.4 | 11.3 |
| HPP 350 MPa, 30 °C | 11.1 | 10.3 | 9.6 | 9.3 |
| HPP 700 MPa, 30 °C | 7.8 | 7.3 | 7.0 | 5.6 |

### 4.2 Fig. 3, pH 8.3 [F ±0.15]

| treatment | 0 min | 5 min | 10 min | 20 min |
|---|---|---|---|---|
| **Thermal 0.1 MPa, 105 °C** | 10.1 | 11.4 | 12.7 | **14.9** |
| PATP 350 MPa, 105 °C | 11.5 | 15.2 | 13.0 | 8.0 |
| PATP 700 MPa, 105 °C | 4.8 | 4.9 | 4.1 | 3.4 |
| HPP 350 MPa, 30 °C | 7.3 | 7.8 | 5.7 | 2.5 |
| HPP 700 MPa, 30 °C | 2.7 | 4.7 | 4.45 | 4.4 |

### 4.3 Table 2 [M] — k (min⁻¹), R²

pH 6.1: PATP 350 0.0232 (0.99); PATP 700 0.0033 (0.69); **thermal 0.0772 (0.99)**; HPP 350 −0.0013 (0.87); HPP 700
−0.0014 (0.97). pH 8.3: PATP 350 −0.04 (0.3); PATP 700 −0.022 (0.89); **thermal 0.056 (0.99)**; HPP 350 −0.006 (0.92);
HPP 700 0.004 (0.34). Negative = net loss. Model: dC/dt = k(C∞ − C)ⁿ, first order (Lau 2003 [C]); C∞ not printed.

### 4.4 Arithmetic [Z]

1. **Molar yield, pH 6.1, 20 min:** 67.0 mg/100 g = 5.87 mmol/kg if the mass is NF (114.1 g/mol) → **8.8 mol % of the
   ribose**; if the number is HDMF-equivalent mass (128.1), 5.2 mmol/kg, 7.9 mol %. Still rising at 20 min. Order of
   Pandit 2006's plateau (2.4 mol %, potato) and above it; Lau 2003 printed no molar yield.
2. **Check of k = 0.0772 (pH 6.1):** first-order approach from C0 = 9.6 through the 5 min point gives C∞ ≈ 91 and predicts
   53 / 74 at 10 / 20 min vs 47 / 67 measured — the fit is approximate (four points, C∞ free).
3. **The 0-min values scatter by treatment group** although every slice comes from one gel: 5.6-11.1 (pH 6.1) and
   2.7-11.5 (pH 8.3), against the printed unprocessed 9.17 / 6.1. **That spread (±~40 %) is the effective "time zero"
   noise** of the pressure arms (pre-heating and come-up are inside it), larger than the error bars.
4. **NF at pH 8.3 vs B61** (hydroxide-catalysed decay, OFF; `parameters_sulfur.py` `NF_HYDROXIDE_DECAY`): at 105 °C B61
   gives k = 5.2 × 10⁻⁴ min⁻¹ at pH 6.1 and **0.082 min⁻¹ at pH 8.3** (t½ 8 min). With formation k ≈ 0.06-0.08 min⁻¹ and
   such a sink, NF at pH 8.3 would sit near a quasi-steady level well below the pH 6.1 curve and creep up with the
   formation flux — qualitatively what Fig. 3 shows. **Consistent with B61 in direction; not a test** (the pH 8.3 gel's
   pH after heating is unknown, the formation flux at pH 8.3 is unknown, and the 0-min level already differs).

## 5. What it can and cannot test

**Can:**
- **The formation side, at one temperature:** NF appearance from ribose + whey protein at 105 °C, pH 6.1, k ≈ 0.08 min⁻¹ —
  a fifth temperature to set against Lau 2003's 116-131 °C series (0.110-0.281 min⁻¹, Ea 81.5 kJ/mol). Lau's Arrhenius
  line extrapolated to 105 °C gives 0.110 × exp(−81.5/R (1/378.15 − 1/389.15)) = **0.053 min⁻¹** [Z]; Gupta's 0.077 is
  1.5× higher — **consistent within the method spread** given a different gel recipe (more protein, salt, lactose) and a
  four-point fit with C∞ free.
- **Lau 2003's missing pH.** Gupta's paste pH 6.1 for a near-Lau gel (same lab family, higher protein, 0.8 % NaCl) is the
  best available estimate of Lau's gel pH; with it, B61 at 116 °C gives ~2.0 × 10⁻³ min⁻¹ (scaling `lau2003` sec. 7's
  1.6 × 10⁻³ at pH 6.0 by 10^0.1), **inside Lau's k_d ≲ 6 × 10⁻³ bound**. Caveat: Lau's recipe differs; paste pH ≠ gel pH after
  80 °C setting.
- **A pH contrast of NF level in a protein matrix** (6.1 vs 8.3: 4.5× lower at 20 min) — direction as Mikami 2017
  (Xyl + Lys, 100 °C: lower NF at pH 8 than 6.5-7) and Nakamura 2020.

**Cannot:**
- **No NF decay constant**: no fed-NF arm, and the pH 8.3 thermal curve still rises.
- **Absolute NF uncertain**: calibrated on furaneol, not norfuraneol; molar absorptivities at 285 nm differ [P]; unit
  basis (NF mass or HDMF-equivalent mass) not stated.
- **pH 8.3 is a paste value** with acetate + bicarbonate added; the gel pH during heating (CO2 loss from bicarbonate, acids
  from ribose) is unknown.
- **Pressure arms** have no engine counterpart (no pressure term) and their 0-min scatter is large.
- Free amine pool not characterised; ~2 % lactose co-present.

## 6. Conflicts and agreements with the repo's records

1. **[!] Standard identity:** the paper calibrates "M-2" on furaneol (HDMF). `pandit2006_extraction.md` records only "a
   commercial standard" for the same lab family; whether that was also HDMF is not checked here (pandit2006.pdf not
   re-opened). If it was, Pandit's
   mg/g plateau (and its 2.4 mol %) carries the same basis question. Flag for the next reader of pandit2006.
2. **Agrees** with Lau 2003's first-order, plateau-type formation and its Ea band (sec. 5).
3. **SUPPORTS hub direction "NF level falls at high pH"** (Mikami, Nakamura, Meynier 140 °C in the F13 atlas), in a protein
   matrix at 105 °C; does not distinguish lower formation from faster decay.
4. Bristow & Isaacs 1999 [C] (pH 7.0 NF rise-then-fall under pressure) — not on disk.

## 7. Proposed ledger entries (aroma_other.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: gupta2011_fig2_fig3_norfuraneol_whey_gel_105c_ph61_ph83
    dossier: gupta2011_extraction.md
    paper: "Gupta, Mikhaylenko, Balasubramaniam & Tang 2011 LWT 44:2141"
    lab: "Washington State (Tang)"
    system: "whey protein gel (WPC 24 + WPI 11 %), D-ribose 1 g/100 g, NaCl 0.8 %, ~2 % lactose"
    medium: "gel, paste pH 6.1 (natural) or 8.3 (sodium acetate + potassium bicarbonate); vacuum-packed slices"
    temperature_c: [105]
    ph: {initial: [6.1, 8.3], measured_during: false, final: null}
    time_min: [0, 5, 10, 20]
    species: [norfuraneol]
    kind: time_series
    location: "Figs. 2-3 p. 2144 (digitised, dossier sec. 4) and Table 2 p. 2145"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "The natural pH of the paste was 6.1"
    flags: [protein_matrix, calibrated_on_furaneol_not_norfuraneol, ph_paste_before_gelation, amine_pool_unknown, lactose_present]
    note: "Thermal 0.1 MPa arm only is engine-relevant: k_form 0.0772 (pH 6.1) / 0.056 (pH 8.3) 1/min; NF 67.0 vs 14.9 mg/100 g at 20 min. Pressure arms: no engine term."
```

## 8. Supplementary material

None referenced.
