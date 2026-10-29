# Garrett & Dvorchik 1969 (J. Pharm. Sci.) — EXTRACTION (four aldopentoses 0.1-1 mM in 0.125-1.00 M HCl, 60-80 °C; furfural by A277 for up to 1400 h; first-order in pentose and in a_H+, Ea 30.3-30.9 kcal/mol, furfural yield 70-92 % via a parallel-loss constant; furfural's own acid solvolysis, Ea 21.3 kcal/mol)

### THE SPECIFIC-ACID PENTOSE → FURFURAL LAW, WITH AN ACTIVATION ENERGY AND A SELECTIVITY — but measured at pH ≈ 0. Extrapolated to food pH it is tiny: at 145 °C and pH 5 the xylose → furfural acid term is ~3 × 10⁻⁷ min⁻¹, four decades below the model's amine-free pentose rate. Specific-acid catalysis is real and first-order in a_H+, but it cannot be what carries pentose to furfural at pH 3-7.

**Source on disk:** `data/articles/garrett1969.pdf` (8 pp., scanned; text layer usable for prose).
**Table I re-typed from a 300 dpi render of p. 815** (scratchpad `pent/gar_tab1.png`,
`gar_tab1foot.png`); Fig. 3 read from the 100 dpi page image (**[D]**, coarse). No supplementary
material.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics and Mechanisms of the Acid Degradation of the Aldopentoses to Furfural" ✔ |
| Authors | **Edward R. Garrett** and Barry H. Dvorchik — College of Pharmacy, University of Florida, Gainesville |
| Venue | Journal of Pharmaceutical Sciences **58**(7), 813-820, July 1969 ✔. Received 30 Jan 1969, accepted 19 Mar 1969; presented APhA Academy, Miami Beach, May 1968 |
| DOI as printed | **not printed**; volume/issue/pages match 10.1002/jps.2600580703 |
| Funding | institutional grant (American Cancer Society); PHS grant CA-10738-06 (NCI) |
| Verdict | **Identity verified.** |

## 1. Why it matters

The model routes pentose to furfural as PENT →(1,2-enolisation)→ TDP →(`r_tdp_fur`, k = 10^-3.49 =
3.2 × 10⁻⁴ min⁻¹ at 145 °C, no pH factor)→ FUR. The 1,2-enolisation is either amine-catalysed
(`r_pent_tdp`, pH factor `"acid"` = 10^(5−pH)) or amine-free (`r_pent_caramel`, pH factor `""`, 5.6 ×
10⁻³ min⁻¹ at 145 °C). Everything sits on the lumped 64 kJ/mol. This paper gives, for the
**amine-free** case only:

1. The acid law is **first order in pentose and first order in HCl / a_H+** (0.125-1.00 M HCl; Eq. 1,
   Eq. 3, Fig. 6).
2. **ΔH‡ = 30.3-30.9 kcal/mol (127-129 kJ/mol), the same for all four pentoses** (Table I).
3. Reactivity **ribose > xylose ≈ lyxose > arabinose, 5.4 : 2.2 : 2.2 : 1**.
4. **Furfural selectivity 70-92 %**: "largely transformed to furfural with small parallel reactions
   (7-25 %)" — the rest goes to "nonchromophoric products" by a parallel step (Eq. 8), not via
   furfural.
5. Furfural itself decays by "specific hydrogen ion-catalyzed solvolysis", Ea 21.3 kcal/mol.
6. The intermediate between pentose and furfural is not seen; the first dehydration (in the pyranose
   form) is rate-determining and the rest is "fast" (Eq. 5, Scheme III).

## 2. Methods as they matter to a model

- **Charge.** Pentose master solution 1 × 10⁻² M; pentose 1-10 × 10⁻⁴ M in the concentration series,
  **1 × 10⁻³ M** in the acid and temperature series; HCl **0.25-1.00 M** (0.125 M for xylose only).
  Medium: water + HCl, **no buffer**, no amine.
- **Temperature.** 60.0, 70.0, 75.0, 80.0 °C (baths ± 0.1 °C). Sealing/atmosphere **not stated**.
- **Observable.** Absorbance at 277 mμ (Beckman DU-2, 1 cm quartz). Furfural ε measured **15,200**
  (lit. 14,800); **15,000** used in all calculations. Spectra of degraded pentose = furfural's.
- **Rates.** "In the time interval 0 to 24 hr., the rate of increase in absorbance at 277 mμ is linear
  with time"; dA/dt from 5-20 h (Table I footnote a). Apparent first-order constants from dA/dt / ε /
  [P]₀. Full time courses to completion at 80 °C, 1 M HCl (Fig. 3, to ~1400 h) fitted on an analog
  computer to Eq. 8 for k₁ (→ furfural) and k₂ (→ other), with k₃ (furfural decay) fixed from separate
  furfural runs.
- **Furfural decay runs.** Furfural 0.021-0.944 × 10⁻⁴ M in acid (Fig. 7), first order.
- **Furfural-addition test.** Ribose 5 × 10⁻⁴ M, 1.00 M HCl, 80 °C, with furfural 0, 1.32, 2.64, 6.6 ×
  10⁻⁴ M added: "no change whatsoever" in furfural production rate or maximum yield (Fig. 8) — the
  intermediate does not react with furfural (Eq. 7 rejected).
- Replicates/SD: **not reported** (n = 3 per Fig. 6 line, n = 6-7 per Fig. 5 line). Activities a_H+
  of HCl at 80 °C taken from Harned & Owen; **the a_H+ values themselves are not printed**.

## 3. Table I re-typed (verbatim; units per footnotes)

dA/dt in absorbance units s⁻¹ (×10⁶ printed); k₁, k₂ in s⁻¹ (×10⁶); k_HCl, k_aH+ in L mol⁻¹ s⁻¹ (×10⁶);
ΔH_a kcal/mol; ΔS_a cal mol⁻¹ K⁻¹.

**(a) 10⁶ dA/dt in 1.00 M HCl at 80.0 °C vs 10⁴[pentose] (M)**

| sugar | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 10 | 10⁶k₁ expt | 10⁶k₁ compt | 10⁶k₂ compt |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ribose | 2.59 | 4.39 | 6.64 | 8.96 | 10.6 | 13.63 | — | 22.69 | 1.54 | 1.67 | 0.545 |
| xylose | 1.28 | 2.08 | 2.97 | 4.12 | 5.20 | 6.11 | — | 8.89 | 0.678 | 0.699 | 0.112 |
| lyxose | 1.17 | 2.10 | 3.12 | 4.38 | 5.55 | 6.25 | — | 9.21 | 0.711 | 0.677 | 0.292 |
| arabinose | 0.48 | 0.94 | 1.42 | 1.85 | — | 2.61 | 3.03 | 4.81 | 0.310 | 0.312 | 0.028 |

**(b) 10⁶ dA/dt at 1 × 10⁻³ M pentose, 80.0 °C vs [HCl] (M)**

| sugar | 0.125 | 0.25 | 0.375 | 0.50 | 0.75 | 1.00 | 10⁶k_HCl | 10⁶k_aH+ |
|---|---|---|---|---|---|---|---|---|
| ribose | — | 4.41 | — | 11.83 | 16.77 | 22.69 | 1.94 | 1.49 |
| xylose | 0.916 | 1.25 | 3.47 | 3.89 | 5.69 | 8.89 | 0.736 | 0.539 |
| lyxose | — | 1.81 | — | 4.23 | 6.73 | 9.21 | 0.808 | 0.597 |
| arabinose | — | 0.90 | — | 2.33 | 2.70 | 4.81 | 0.354 | 0.314 |

**(c) 10⁶ dA/dt at 1 × 10⁻³ M pentose, 1.00 M HCl vs T (°C)**

| sugar | 80.0 | 75.0 | 70.0 | 60.0 | ΔH_a | log P | ΔS_a |
|---|---|---|---|---|---|---|---|
| ribose | 22.69 | 11.75 | 9.54 | 1.96 | 30.9 | 17.5 | 21.2 |
| xylose | 8.89 | 4.48 | — | 0.60 | 30.3 | 17.0 | 18.9 |
| lyxose | 9.21 | 4.83 | — | 0.63 | 30.3 | 17.0 | 18.9 |
| arabinose | 4.27 | 2.24 | — | 0.31 | 30.8 | 16.7 | 17.5 |

Footnote g: log k = −(ΔH_a/2.303R)(1/T) + log P, where k is "the second-order rate constant for the
degradation of the sugars" from (dA/dt)/[pentose][HCl] at 1 mM, 1 M HCl; "all ΔH_a's can be considered
to be the same"; log P values computed with a **common ΔH_a of 30.8 kcal/mol**.

**Furfural degradation (text p. 814):** k₃ = **0.582 × 10⁻⁶ s⁻¹** (80 °C, 1 M HCl); k_aH+ = **0.79 ×
10⁻⁶** and k_HCl = **0.58 × 10⁻⁶ L mol⁻¹ s⁻¹**; **ΔH_a = 21.3 kcal/mol**, log P = 7.9.

**Fig. 3 [D], peak furfural** (1 mM pentose, 1 M HCl, 80 °C; read ± 0.15 × 10⁻⁴ M, ± 50 h): ribose
≈ 4.9 × 10⁻⁴ M at ≈ 300 h; xylose ≈ 3.8 × 10⁻⁴ M at ≈ 480 h; lyxose ≈ 3.7 × 10⁻⁴ M at ≈ 420 h;
arabinose ≈ 2.6 × 10⁻⁴ M at ≈ 700 h. (Peak ≠ yield: furfural decays concurrently.)

## 4. Arithmetic (mine)

1. **Furfural selectivity k₁/(k₁+k₂)** (compt values, 80 °C, 1 M HCl): ribose **75 %**, xylose
   **86 %**, lyxose **70 %**, arabinose **92 %**. The abstract's "7-25 %" side reaction matches all but
   lyxose (30 %).
2. **Ea in kJ/mol:** ribose 129.3, xylose 126.8, lyxose 126.8, arabinose 128.9; common **128.9**.
   Furfural solvolysis **89.1 kJ/mol**. Two-point checks from (c): ribose 60/80 °C → 120 kJ/mol;
   xylose → 132 kJ/mol. Ribose 70 → 75 °C gives 41 kJ/mol: **the 70 °C point is anomalous** (Flag 2).
3. **Extrapolated specific-acid furfural formation** (k_aH+ × a_H+, a_H+ = 10^−pH, own Ea; min⁻¹):

| sugar | T | k_aH+ (L mol⁻¹ s⁻¹) | pH 3 | pH 5 |
|---|---|---|---|---|
| xylose | 100 °C | 5.5 × 10⁻⁶ | 3.3 × 10⁻⁷ | 3.3 × 10⁻⁹ |
| xylose | 120 °C | 4.4 × 10⁻⁵ | 2.6 × 10⁻⁶ | 2.6 × 10⁻⁸ |
| xylose | 145 °C | 4.4 × 10⁻⁴ | **2.7 × 10⁻⁵** | **2.7 × 10⁻⁷** |
| ribose | 145 °C | 1.4 × 10⁻³ | 8.4 × 10⁻⁵ | 8.4 × 10⁻⁷ |

   Extrapolation spans 3-5 decades of a_H+ and 20-65 K; order-of-magnitude only.
4. **Crossover with Larralde 1995** (ribose, 80 °C; Larralde's slope-1 base line log k = pD − 11.53,
   pD = pH + 0.4): the specific-acid term is **95 % of total loss at pH 2, 67 % at pH 2.5, 17 % at pH 3,
   2 % at pH 3.5, 0.2 % at pH 4**. Free-pentose loss has its minimum near **pH 2.5-3** at 80 °C.
   At 145 °C (ribose), acid share at pH 3 ≈ 45 % (8.4 × 10⁻⁵ vs Larralde's 1.0 × 10⁻⁴ min⁻¹).
5. **Furfural acid solvolysis at food pH:** k_aH+ at 145 °C = 8.8 × 10⁻⁵ L mol⁻¹ s⁻¹ → **5 × 10⁻⁸ min⁻¹
   at pH 5** — negligible against any furfural sink in the corpus (Yaghmur 2005's water-arm k_obs at
   65 °C is 2.1 × 10⁻³ min⁻¹).

## 5. Kinetic numbers the repository can use

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| xylose → furfural (specific acid) | k_aH+ | **0.539 × 10⁻⁶** | L mol⁻¹ s⁻¹ | 80 °C, 0.125-1 M HCl, 1 mM | Table I(b) | printed |
| ribose / lyxose / arabinose → furfural | k_aH+ | 1.49 / 0.597 / 0.314 × 10⁻⁶ | L mol⁻¹ s⁻¹ | same | Table I(b) | printed |
| aldopentose → furfural | ΔH‡ | **30.3-30.9** (common 30.8) | kcal/mol | 60-80 °C, 1 M HCl | Table I(c) | printed |
| xylose → furfural (1 M HCl) | k₁ | 0.678 × 10⁻⁶ (expt), 0.699 × 10⁻⁶ (fit) | s⁻¹ | 80 °C | Table I(a) | printed |
| xylose → non-furfural (parallel) | k₂ | 0.112 × 10⁻⁶ | s⁻¹ | 80 °C, 1 M HCl | Table I(a) | printed (analog-computer fit) |
| furfural selectivity | k₁/(k₁+k₂) | 75 / 86 / 70 / 92 % (rib/xyl/lyx/ara) | — | 80 °C, 1 M HCl | arithmetic 1 | derived |
| furfural → non-chromophoric (acid) | k₃; k_aH+; ΔH‡ | 0.582 × 10⁻⁶ s⁻¹; 0.79 × 10⁻⁶ L mol⁻¹ s⁻¹; **21.3 kcal/mol** | — | 80 °C, 1 M HCl | p. 814 | printed |
| relative reactivity | ratio | 5.4 : 2.2 : 2.2 : 1 (rib : xyl : lyx : ara) | — | 80 °C | abstract | within_study_ratio |

## 6. Can / cannot

**Can:** fix the **form and magnitude of the specific-acid term** of amine-free pentose → furfural
(and its Ea ≈ 129 kJ/mol), the **furfural selectivity of that route** (70-92 %), and show that
furfural acid solvolysis (Ea 89 kJ/mol) is irrelevant above pH ~2.

**Cannot:** anything at pH > ~0.5 directly; any amine-catalysed step; the 3-deoxypentosone (never
observed, so no `k_tdp_fur`); FFT; route shares.

## 7. Flags

1. **pH ≈ −0.1 to 0.9 only.** Every number is in 0.125-1 M HCl. The model's pots are pH 3-7.5.
2. **Ribose 70 °C point (9.54) is inconsistent** with the 75 °C point and the common ΔH‡.
3. **k_HCl ≠ k₁ expt / 1 M** (ribose 1.94 vs 1.54): k_HCl is the slope over the acid series, which
   is not proportional at the low end (xylose 0.25 → 0.375 M jumps 2.8x). Prefer k_aH+.
4. **a_H+ values not printed**; reproduced only via Harned & Owen (ref. 9).
5. Analog-computer k₁, k₂ come with no uncertainties.
6. "Furfural" = A277; any 277-nm-absorbing co-product would be counted (spectra matched furfural's).

## Use in the model

- **Gap (c), the 1,2- vs 2,3-enolisation pH law at low pH — contradicts the model's framing.** The
  model's only acid factor is a full decade per pH unit on the amine-catalysed `r_pent_tdp`,
  normalised at pH 5, so going pH 5 → 3 multiplies the 1,2-branch 100x. Garrett confirms a
  slope-1 specific-acid law for the **amine-free** pentose → furfural reaction but sizes it: at
  145 °C xylose, **2.7 × 10⁻⁷ min⁻¹ at pH 5 and 2.7 × 10⁻⁵ min⁻¹ at pH 3**, i.e. 4 and 2 decades below
  the model's amine-free `k_pent_caramel` (5.6 × 10⁻³ min⁻¹, flat in pH). With Larralde 1995 it puts
  the minimum of free-pentose loss near pH 2.5-3 and makes the acid route a minor share above pH 3.
  So for the free sugar, **the low-pH partition is not acid-dominated at pH 3**; a correct pH law
  needs the acid term at Garrett's magnitude plus a base/water term at Larralde's. Whether the
  *amine-catalysed* 1,2-enolisation (via the Schiff base / Amadori) is as steep is not addressed;
  that is where the model's pH-3 error (MFT 22x low, FFT 6.6x high) lives, and this paper is
  consistent with the flux map's "too steep" reading but does not prove it.
- **Gap (e), furfural branch temperature dependence.** Pentose → furfural Ea ≈ **129 kJ/mol**
  (specific-acid route) is twice the model's lumped 64 kJ/mol on the formation lane; a steeper
  formation Ea makes furfural (and FFT) fall faster as T drops — the direction the Yiltirak FFT
  ladder (234x high at 100 °C) asks for. Usable as a **declared band** for an un-lumped
  pentose → furfural Ea, caveat: different catalysis regime from the pots.
- **Furfural decay barrier.** Furfural's acid solvolysis Ea 89 kJ/mol is below the model's
  143 kJ/mol family barrier, but at pH 5 this route is ~10⁻⁸ min⁻¹ and cannot be the furfural sink;
  it does **not** support splitting `r_fur_decay`'s barrier toward 89 kJ/mol (Yaghmur's 46.5 kJ/mol
  cysteine-mediated sink remains the relevant number).
- **Selectivity.** 70-92 % of acid-degraded pentose becomes furfural (Eq. 8 parallel loss
  7-30 %), so if a free-sugar acid route is carried it should carry a non-furfural branch of
  ~10-25 % (xylose 14 %).
- No FFT route shares; nothing for gaps (a), (b), (d), (f).
