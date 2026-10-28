# Mikami, Nakamura, Yamada & Murata 2017 — EXTRACTION (xylose / ribose / arabinose 13.3 mmol/L + L-lysine·HCl 34 mmol/L in 0.2 M Na/K phosphate pH 5.5-8.0, 100 °C, 30-120 min: norfuraneol (HMFO) by HPLC-DAD 280 nm; isolated norfuraneol 1.7-1.8 mmol/L heated alone at pH 8, 100 °C, 0-6 h: its decay and its methylglyoxal / diacetyl fragments by OPD trapping)

### THE ONE PAPER IN THE CORPUS THAT HEATS PURE NORFURANEOL IN BUFFER AND FOLLOWS ITS LOSS: at pH 8, 100 °C it decays first-order with **k ≈ 0.032 min⁻¹ (t½ ≈ 22 min)** — about **700× faster than the engine's `k_nf_decay` at 100 °C** — and fragments to diacetyl (≤ 5.5 mol %) and methylglyoxal (≤ 0.6 mol %). In the xylose-lysine pot, norfuraneol peaks at pH 6.5 at ~13-16 mol % of the pentose, and it is still rising at 2 h at pH 6 but already falling after 1 h at pH 7-8.

**Source on disk:** `data/articles/mikami2017.pdf` (7 pp., Food Sci. Technol. Res. 23 (2), 283-289,
2017; publisher PDF, clean text layer). Prose and methods from `pdftotext -layout`. **Every number in
Figs. 3 and 6 is figure-read** from 200 dpi renders (scratchpad `fork_nf/mik-4.png`, `mik-6.png`)
and tagged **[D]** with its read error. Fig. 4 (photos + UV spectra) and Fig. 5 (GPC) are
qualitative. No SI.

## 0. Identity

| field | value |
|---|---|
| Title | "4-Hydroxy-5-Methyl-3(2H)-Furanone (HMFO) Contributes to Browning in the Xylose-Lysine Maillard Reaction System" |
| Authors | Yoko Mikami, Miki Nakamura, Shinji Yamada, **Masatsune Murata** (corresponding, murata.masatsune@ocha.ac.jp) — Ochanomizu University, Tokyo |
| Venue | Food Science and Technology Research **23** (2), 283-289, 2017. Received 17 Oct 2016, accepted 17 Nov 2016 |
| DOI as printed | **10.3136/fstr.23.283** ✔ matches the brief. **Identity verified.** |
| Naming | HMFO = 4-hydroxy-5-methyl-3(2H)-furanone = norfuraneol = the engine's `NF` (C5H6O3, M 114.1; structure by X-ray, CCDC 1500343) |

## 1. Why it matters

Gap (d), norfuraneol formation/decay vs pH. In the engine, `NF` is made by `r_dpo_nf` (DPO → NF,
first order, no pH tag; log10 k(145 °C) = +0.43 min⁻¹ in B60) and lost by `r_nf_decay` (NF →
FRAG_C 5, `k_nf_decay`, **no pH tag**, log10 k(145 °C) = −2.18 min⁻¹, in the shared carbonyl-sink
barrier family at **143.5 kJ/mol**), plus the H2S-consuming `r_nf_mft`, `r_nf_mp3p`, `r_nf_shet`.
That puts the engine's NF thermal decay at 100 °C at **4.5 × 10⁻⁵ min⁻¹** (t½ ≈ 250 h). This paper
measures it at pH 8 as **3.2 × 10⁻² min⁻¹**, and in the lysine pot at pH 7 a net loss that bounds it
from below at ~2 × 10⁻³ min⁻¹. It also identifies what NF fragments to — **diacetyl and
methylglyoxal** — whereas the engine sends it to inert FRAG_C.

It also gives a clean amine-sugar NF yield–pH curve (pH 5.5-8, phosphate, 100 °C), a time series
at three pH values, and the pentose ranking ribose > xylose ≈ arabinose.

## 2. Methods as they matter to a model

- **Pots.** "Buffer solutions ... containing 13.3 mM sugar (xylose, ribose, arabinose, or glucose) and
  34 mM l-Lys monohydrochloride" — 0.2 M phosphate-Na/K (pH 5.5-8.0), 0.2 M acetate-Na (pH 3.5-6.0)
  or Tris-HCl (pH 7.0-9.0); "heated for 30-120 min in boiling water" in a capped test tube (volume
  not stated). **Only phosphate data are shown** (Figs. 3A-C). 13.3 mM xylose = 200 mg/100 mL.
- **pH measured?** Initial buffer pH only; no post-heating pH.
- **NF analysis.** RP-HPLC-DAD (YMC ODS-A, 0.1 % TFA water/MeOH gradient), NF at ~7.5 min, λmax 280 nm.
  **Calibration not described**; the authors had ~330 mg of pure crystalline NF (isolated from a 0.3 M
  xylose + 0.3 M Lys, pH 7.0 phosphate, 60 min pot), so an external standard is likely but unstated.
- **NF stability.** "A solution of 0.2 M phosphate buffer (pH 8.0) containing 1.8 mM HMFO was heated at
  100 °C for 0.5, 1, 3, and 6 h" (Fig. 6 caption says 1.7 mM), then reacted with an equal volume of 1.8 mM
  o-phenylenediamine (OPD), 5 h at 20 °C; quinoxalines by HPLC-DAD and LC-MS, identified against
  synthesised MGO- and diacetyl-quinoxaline standards. Browning series (Fig. 4): 1.8 mM NF ± 3.4 mM Lys,
  pH 6/7/8, 0-6 h.
- **Replicates.** Fig. 3: n = 3, error bars shown. Fig. 6: error bars on some points.

## 3. Figures digitised

### Fig. 3A — NF vs pH, xylose 13.3 mM + Lys 34 mM, phosphate 0.2 M, 100 °C, 60 min [D]

100 % = **24 mg/100 mL = 2.10 mmol/L** (caption). Read error ± 3 % of scale.

| pH | 5.5 | 6.0 | 6.5 | 7.0 | 7.5 | 8.0 |
|---|---|---|---|---|---|---|
| NF, % of max | 27 | 61 | 100 | 83 | 60 | 21 |
| NF, mmol/L (mine) | 0.57 | 1.28 | 2.10 | 1.74 | 1.26 | 0.44 |
| mol % of xylose (mine) | 4.3 | 9.6 | **15.8** | 13.1 | 9.5 | 3.3 |

### Fig. 3B — NF vs time at pH 6 / 7 / 8, same pot [D]

100 % = **20 mg/100 mL = 1.75 mmol/L**. Read error ± 3 %.

| pH | 30 min | 1 h | 2 h |
|---|---|---|---|
| 6 | 30 % (0.53 mM) | 67 % (1.17 mM) | 102 % (1.79 mM) |
| 7 | 52 % (0.91 mM) | 80 % (1.40 mM) | 72 % (1.26 mM) |
| 8 | 23 % (0.40 mM) | 30 % (0.53 mM) | 19 % (0.33 mM) |

Text: maximum "with 1 h heating at pH 7 or 8 ... compared to 4 h heating at pH 6 (data not shown)".

### Fig. 3C — pentose comparison, 13.3 mM sugar + 34 mM Lys, 100 °C, 1 h [D]

100 % = **33 mg/100 mL = 2.89 mmol/L**. Read error ± 3 %.

| sugar | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| arabinose | 33 % | 54 % | 25 % |
| ribose | **100 %** (2.89 mM, 21.7 mol %) | 83 % | 31 % |
| xylose | 43 % | 64 % | 23 % |

### Fig. 6B — pure NF (1.7-1.8 mM), 0.2 M phosphate pH 8.0, 100 °C [D]

NF on the right axis (0-2.0 mM, read ± 0.03 mM); quinoxalines on the left (0-0.15 mM, read ± 0.003 mM).

| t (h) | 0 | 0.5 | 1 | 3 | 6 |
|---|---|---|---|---|---|
| NF (mM) | 1.72 | 0.63 | 0.25 | ~0.02 | ~0 |
| diacetyl-quinoxaline II (mM) | 0 | 0.030 | 0.063 | **0.094** | 0.046 |
| MGO-quinoxaline I (mM) | 0 | **0.011** | 0.004 | 0.003 | ~0.001 |

### Arithmetic (mine)

1. **NF first-order decay, pH 8, 100 °C:** ln-linear fit of 0 / 0.5 / 1 h → **k = 0.032 min⁻¹
   (t½ 22 min)**; pointwise 0.034 (0-0.5 h) and 0.032 (0-1 h) — first order holds within read error.
2. **Engine comparison:** `k_nf_decay` = 10^−2.18 = 6.6 × 10⁻³ min⁻¹ at 145 °C, Ea 143.5 kJ/mol →
   **4.5 × 10⁻⁵ min⁻¹ at 100 °C. Mikami / engine ≈ 710×** (pH 8). Even at 145 °C the engine's value is
   5× below this paper's 100 °C value.
3. **Lower bound at pH 7 (pot):** NF falls 1.40 → 1.26 mM from 1 to 2 h while still being formed →
   net loss k ≥ ln(1.40/1.26)/60 = **≥ 1.8 × 10⁻³ min⁻¹** at pH 7, 100 °C (≥ 40× the engine; in the
   presence of 34 mM Lys, which may react with NF).
4. **pH 6:** no net loss to 2 h; maximum at ~4 h (not shown). Consistent with decay ≪ 10⁻³ min⁻¹ at
   pH 6 or with slower formation; the paper says both are slower ("formed and decomposed more slowly at
   pH 6"). **Decay is pH-dependent**; only the pH 8 value is a clean number.
5. **Fragment yields (pH 8):** diacetyl peak 0.094 / 1.72 = **5.5 mol %** of NF at 3 h; MGO peak
   0.011 / 1.72 = **0.6 mol %** at 0.5 h. Both are transients (they themselves decay/polymerise), so
   these are **lower bounds** on the fragment fluxes; the OPD derivatisation efficiency is not reported.
   Most NF carbon goes to "colorless polymers" and brown pigments (Fig. 5, GPC) — unquantified.
6. **Yield conversions.** The text's "about 10 %" NF yield from xylose is **w/w** (20 mg from 200 mg);
   molar it is 13.2 mol % (× 150.1/114.1). At the pH 6.5 optimum (Fig. 3A) it is 15.8 mol %.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| **NF thermal decay** | first-order k | **0.032 (± ~10 %)** | min⁻¹ | 1.7-1.8 mM NF alone, 0.2 M Na/K phosphate pH 8.0, 100 °C | Fig. 6B [D], arithmetic 1 | **derived_from_digitised** (clean, fed-precursor) |
| NF net loss in Lys pot | k floor | ≥ 1.8 × 10⁻³ | min⁻¹ | Xyl 13.3 + Lys 34 mM, pH 7 phosphate, 100 °C, 1-2 h | Fig. 3B [D] | derived bound |
| NF → diacetyl | transient max | ≥ 5.5 | mol % of NF | pH 8, 100 °C, 3 h | Fig. 6B [D] | lower bound (transient) |
| NF → MGO | transient max | ≥ 0.6 | mol % of NF | pH 8, 100 °C, 0.5 h | Fig. 6B [D] | lower bound (transient) |
| NF yield vs pH | 60-min level | Fig. 3A table | mmol/L, mol % | Xyl + Lys, phosphate 5.5-8, 100 °C | Fig. 3A [D] | digitised (n = 3) |
| NF vs time | levels | Fig. 3B table | mmol/L | pH 6 / 7 / 8 | Fig. 3B [D] | digitised time series (3 points) |
| NF by pentose | 60-min level | Fig. 3C table | % of 2.89 mM | Rib / Xyl / Ara, pH 6-8 | Fig. 3C [D] | within-study ratio |
| NF optimum | pH of max at 1 h | 6.5 | — | Xyl + Lys | Fig. 3A | measured |
| NF Ea | — | **NOT REPORTED** (one temperature) | — | — | — | — |

## 5. What this can and cannot anchor

**Can:**
1. **A measured `k_nf_decay` at one (T, pH)** — the first fed-norfuraneol-alone decay in the corpus
   (the Bornhorst constants are formation constants and `parameters_sulfur.py` forbids reading them as
   decay). It is far outside the engine's current value (arithmetic 2).
2. **A pH dependence for NF decay** (qualitative below pH 8; floor at pH 7).
3. **NF → diacetyl / MGO** as products of the NF sink, with lower-bound yields.
4. **A pre-registered hold-out pot** (Xyl + Lys, phosphate, pH 5.5-8, 100 °C, 60 min) for the
   amine-sugar → DPO → NF branch and its pH optimum — **if** the engine's lysine-generic amine route
   can run it; the sulfur lane's DPO entry is cysteine-catalysed (`r_pent_dpo`), so a lysine pot
   tests the ARP (`r_arp_dpo`) route instead.

**Cannot:**
- No temperature dependence (100 °C only) — no Ea for formation or decay.
- No sulfur: nothing on NF + H2S.
- Phosphate-only shown; acetate/Tris runs described but not reported.
- Calibration of NF unstated; pH not measured after heating; NF in the Lys pot may be consumed by Lys
  (Lys "seemed to repress the colorization").

## 6. Flags

1. **Three different 100 % scales** in Fig. 3 (24, 20, 33 mg/100 mL). Cross-check: xylose pH 7, 1 h
   reads 1.74 (3A), 1.40 (3B) and 1.85 mM (3C) — **±15 % between panels of the same condition**;
   use that as the replicate-level uncertainty, not the error bars.
2. **"30 mg/100 mL of HMFO ... from 100 mg/100 mL of ribose"** (text) contradicts the methods
   (13.3 mM = 200 mg/100 mL sugar). Taking the methods, ribose gives 21.7 mol %; taking the text it
   would be 43 mol %, implausible. Use the methods.
3. **Fig. 4B's 280 nm absorbance** falls only to ~50 % at 0.5-1 h, versus Fig. 6B's NF falling to
   37 % and 15 %. The UV trace includes product absorbance; do not read NF from Fig. 4B.
4. **Initial NF 1.7 vs 1.8 mM** (Fig. 6 caption vs methods). Immaterial to k.
5. **Lys·HCl at 34 mM in 0.2 M phosphate** — moderate buffer capacity; no measured pH drift.
6. Registry: `NF` exists in the engine; "HMFO" is a synonym to add to `data/keys/compounds.yml` if the
   paper is keyed. Paper not in `data/keys/papers.yml` (not checked further).

## 7. Use in the model

- **Gap (d) — direct contradiction.** The engine's `r_nf_decay` is pH-independent with a 143.5 kJ/mol
  barrier shared with furfural and the osones, giving t½ ≈ 250 h at 100 °C. Mikami measures t½ ≈ 22 min
  at pH 8 / 100 °C and a net loss at pH 7 within 1 h. Either NF decay needs a **base-catalysed term**
  (pH-tagged, like the 2,3-enolisation steps), or its barrier must leave the carbonyl-sink family. A
  pre-registered test: give `k_nf_decay` a hydroxide/`base` factor anchored on this pH 8 point (declared
  band ±0.3 dex), and score Hofmann's pH 3/5/7 sweep and Zhou 2023c (pH 7.5, MFT high 14-161×) as
  hold-outs — a faster NF sink at pH ≥ 7 pushes MFT down exactly where the model is too high in Zhou.
- **Gap (f) — a cross-link.** NF → diacetyl / MGO means the NF sink feeds the MGO pool that
  `r_cys_actz` draws on (2-acetylthiazole). The engine's NF decay goes to FRAG_C, so this source of
  MGO is missing; it is small (≥ 0.6 mol % MGO, transient) but not zero, and diacetyl (≥ 5.5 %) is the
  thiophene/thiazole precursor in other routes.
- **Gap (c) — formation pH law.** The NF optimum at pH 6.5 with Lys (falling 4× to pH 5.5 and 5× to
  pH 8 at 60 min) is a net (formation − decay) curve; with the pH 8 decay known, the formation side
  at pH 8 can be backed out (rough steady state: formation ≈ 0.032 × 0.44 mM ≈ 0.014 mM/min).
- **Caveats:** lysine, not cysteine (no thiazolidine sequestration, different Amadori basicity);
  one temperature; figure-read values (± 3 % of scale, ± 15 % between panels); NF quantitation method
  not described.
