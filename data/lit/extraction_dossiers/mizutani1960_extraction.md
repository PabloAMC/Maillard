# Mizutani, Obata & Ishikawa 1960 — EXTRACTION (crystalline cysteine–aldehyde thiazolidines of formaldehyde, acetaldehyde, n-butyraldehyde, benzaldehyde and furfural, plus free cysteine and cystine, 1 mmol S or N in 10 mL water, heated in a 120 °C oil bath under a CO2 (or N2) sweep for 0.5 / 1 / 3 / 6 h; cumulative H2S (methylene blue) and NH3 (Nessler))

### THE ONLY SIDE-BY-SIDE MEASUREMENT IN THE CORPUS OF H2S RELEASE FROM FREE CYSTEINE AND FROM ITS THIAZOLIDINES UNDER ONE PROTOCOL — and two things fall out. Free cysteine's H2S (0.95 % of S in 6 h) matches the engine's Zheng & Ho thermolysis constant evaluated at 100 °C, pH 5, to within 1.5× at all four times. And binding cysteine in a thiazolidine lowers H2S release by only 1.2-3.5× (simple aldehydes) — while the furfural thiazolidine releases 2× MORE than free cysteine, and 6× more after a month's storage.

**Source on disk:** `data/articles/mizutani1960.pdf` (5 pp.: a Taylor & Francis cover + journal pp. 382-385).
Read from 150-dpi renders via the Read tool; **Tables I and II re-checked on 250-dpi crops** (scratchpad
`img/miz_t1d.png`, `img/miz_t2.png`). The scan carries handwritten annotations (underlining, Japanese marginalia)
that do not touch any number. Figs. 1 and 3 are calibration curves; Fig. 2 is the apparatus.

## 0. Identity

| field | value |
|---|---|
| Title (cover) | "Studies on the Degradation Mechanisms of Proteins and their Derivatives in the Foods" |
| Title (article) | same, **"Part III. The Decompositions of Several Cysteine-aldehyde Compounds in Aqueous Solutions"** ✔ |
| Authors | Junya Mizutani, Yataro Obata, Yoshinori Ishikawa — Dept. of Agricultural Chemistry, Faculty of Agriculture, Hokkaido University, Sapporo ✔ |
| Venue | Bulletin of the Agricultural Chemical Society of Japan **24**(4), Suppl., 382-385 (1960); received 3 Oct 1959 ✔ |
| DOI as printed | **10.1080/03758397.1960.10857681** ✔ matches the brief |
| Context | Part of a series on the browning of canned "Yamatoni" beef seasoned with soy-bean meal hydrolysate; the authors were chasing heat-labile sulfur compounds |

## 1. Why it matters

Gap (b), the de Roos path (`SULFUR_FLUX_MAP.md` §4, §8): B59 plans cysteine + pentose ⇌ TTCA as a **reservoir**
from which cysteine is released slowly; the engine's TTCA (`r_ttca_cys`, `r_ttca_deg`) releases **no H2S of its
own** — its sulfur goes to `FRAG_S` on degradation, and H2S comes only from free cysteine (`r_cys_h2s`). This paper
measures the H2S that thiazolidines release in hot water, against free cysteine, under one protocol:

| compound (1 mmol S) | H2S at 6 h, % of S | relative to free cysteine (6 h / 3 h) |
|---|---|---|
| L-cystine | 1.40 | 1.46 / 1.42 |
| **L-cysteine** | **0.95** | 1 / 1 |
| thiazolidine of formaldehyde | 0.78 | 0.81 / 0.74 |
| thiazolidine of acetaldehyde | 0.68 | 0.72 / 0.69 |
| thiazolidine of n-butyraldehyde | 0.29 | 0.30 / 0.24 |
| thiazolidine of benzaldehyde | 0.57 | 0.60 / 0.50 |
| **thiazolidine of furfural (fresh)** | **2.06** | **2.15 / 1.96** |
| thiazolidine of furfural (stored > 1 month, coloured) | n.d. (2.13 at 3 h) | — / **5.98** |

(all mine from Table I; % = µg H2S / 34.08 / 998 µmol S)

Two readings matter for the model: (i) the thiazolidine of a simple aldehyde is **not** an H2S-tight reservoir — it
releases 30-80 % as much H2S as free cysteine; (ii) the **furfural** thiazolidine is a **net H2S source** relative to
cysteine, which bears on the furfural/FFT branch (gap e): furfural binding cysteine does not sequester its sulfur.

## 2. Methods as they matter to a model

- **Thiazolidine synthesis (p. 383).** L-cysteine 8.3 mmol in 20 mL 50 % ethanol + aldehyde 8.3 mmol, overnight
  in an ice-box; crystals recrystallised. Yields: formaldehyde 0.9 g (82 %), acetaldehyde 0.9 g (73 %),
  n-butyraldehyde 1.0 g (69 %), benzaldehyde 1.3 g (77 %), furfural 1.1 g (66 %). N and S analyses within
  ~0.5 % of theory. Nitroprusside-negative (no free SH), ninhydrin yellow. The furfural compound "was unstable
  very much and gradually changed into grayish brown". **The 5-HMF compound gave no crystals** and was not tested.
  L-cysteine was made from hair cystine via Sn/HCl reduction.
- **H2S (p. 383-384).** Sample "containing 32 mg of sulfur" (= 0.998 mmol S) + 10 mL water in a flask;
  CO2 bubbled through; heated **0.5, 1, 3, 6 h (separate runs) in an oil bath at 120 °C**; then 10 mL
  hydrochloric acid (concentration not stated) added and CO2 passed 30 min more; H2S trapped in 2 % zinc acetate,
  measured as methylene blue at 745 mµ. "Recovery of hydrogen sulfide of this apparatus was 95±5%".
  N2 instead of CO2 "the same results were obtained".
- **NH3 (p. 384-385).** Sample "containing 14 mg N" (= 0.999 mmol N), 10 mL water, oil bath 120 °C, **3 and 6 h**;
  bath cooled below 80 °C, 10 mL 10 % NaOH, aerated 30 min into 0.1 N HCl, Nessler (Folin-Wu), 420 mµ.
- **Concentration:** **0.1 mol/L** cysteine-equivalent in 10 mL (S basis; cystine 0.05 mol/L).
- **pH: NOT reported, NOT buffered.** Free cysteine in water sits near its pI (~5.1); the thiazolidine-4-carboxylic
  acids are weaker bases (ring N pKa ~6.2) so their solutions are more acidic (~pH 4, my estimate); the CO2 sweep
  pushes all towards pH 4-5. Treat every number as **unbuffered, pH ~4-5.5, drifting**.
- **Temperature, the key caveat:** 120 °C is the **oil bath**. The flask is swept by gas at ~1 atm through an
  outlet to the traps (Fig. 2), so the liquid **cannot exceed its boiling point, ~100 °C**. I read the solution as
  at **~100 °C (gently boiling, ~85 % confidence)**; the vertical tube over the flask in Fig. 2 suggests reflux.
- **Statistics:** none — single determinations, no replicates.

## 3. Tables re-typed (verbatim)

### Table I. Hydrogen sulfide (µg) liberated from L-cysteine-aldehyde compounds in heated aqueous solutions (oil bath 120 °C); compound containing 32 mg S

| compound | 0.5 h | 1 h | 3 h | 6 h |
|---|---|---|---|---|
| L-cystine | 8.5 | 29 | 172 | 475 |
| L-cysteine | 19 | 47 | 121 | 325 |
| L-cysteine-formaldehyde | 10 | 20 | 90 | 264 |
| L-cysteine-acetaldehyde | 8 | 19 | 84 | 232.5 |
| L-cysteine-n-butyraldehyde | — | 7 | 29 | 97 |
| L-cysteine-benzaldehyde | 7 | 15 | 60 | 195 |
| L-cysteine-furfural (I) | 32 | 68 | 237 | 700 |
| L-cysteine-furfural (II) | 106 | 195 | 723 | — |

Footnote: furfural (II) "was storaged for more than a month and considerably colored".

### Table II. Ammonia (µg) liberated from L-cysteine-aldehyde compounds (oil bath 120 °C); compound containing 14 mg N

| compound | 3 h | 6 h |
|---|---|---|
| L-cystine | 99 | 160 |
| L-cysteine | 127 | 228 |
| L-cysteine-formaldehyde | 69 | 147 |
| L-cysteine-acetaldehyde | 52 | 132 |
| L-cysteine-n-butyraldehyde | 24 | 65 |
| L-cysteine-benzaldehyde | 56 | 119.5 |
| L-cysteine-furfural | 146 | 311 |

### Converted (mine): % of S as H2S, % of N as NH3

| compound | H2S 0.5 h | 1 h | 3 h | 6 h | NH3 3 h | NH3 6 h |
|---|---|---|---|---|---|---|
| cystine | 0.025 | 0.085 | 0.51 | 1.40 | 0.58 | 0.94 |
| **cysteine** | **0.056** | **0.138** | **0.356** | **0.955** | **0.75** | **1.34** |
| formaldehyde adduct | 0.029 | 0.059 | 0.265 | 0.776 | 0.41 | 0.86 |
| acetaldehyde adduct | 0.024 | 0.056 | 0.247 | 0.683 | 0.31 | 0.78 |
| n-butyraldehyde adduct | — | 0.021 | 0.085 | 0.285 | 0.14 | 0.38 |
| benzaldehyde adduct | 0.021 | 0.044 | 0.176 | 0.573 | 0.33 | 0.70 |
| furfural adduct (I) | 0.094 | 0.200 | 0.697 | 2.058 | 0.86 | 1.83 |
| furfural adduct (II, aged) | 0.312 | 0.573 | 2.125 | — | — | — |

### Arithmetic (mine)

**1. Free cysteine vs the engine's thermolysis constant.** `cysteine_thermolysis_k` (Zheng & Ho 1994 matched
pairs, 1/min) gives cumulative H2S at 0.5 / 1 / 3 / 6 h of:

| solution T, pH (assumed) | 0.5 h | 1 h | 3 h | 6 h | measured / model at 6 h |
|---|---|---|---|---|---|
| 100 °C, pH 5 | 0.084 % | 0.168 % | 0.50 % | 1.00 % | **0.95** |
| 100 °C, pH 3 pair | 0.076 % | 0.152 % | 0.46 % | 0.91 % | 1.05 |
| 120 °C, pH 5 | 0.74 % | 1.47 % | 4.36 % | 8.52 % | 0.11 |
| 100 °C, pH 7 | 0.46 % | 0.91 % | 2.71 % | 5.35 % | 0.18 |

Measured: 0.056 / 0.138 / 0.356 / 0.955 %. **At 100 °C and pH 3-5 the model is within 0.67-1.05× at every time.**
At a true 120 °C the model would be 9× high. Since the solution cannot exceed ~100 °C in this apparatus, this is
an **independent confirmation of the Zheng & Ho constant at 100 °C, pH ~5, 0.1 M**, from a different lab, 34 years
earlier, with a different sulfide assay. Interval rates (fraction of S per min): 1.9e-5 (0-0.5 h), 2.7e-5 (0.5-1 h),
1.8e-5 (1-3 h), 3.3e-5 (3-6 h) — roughly constant to mildly accelerating; the model's 2.8e-5 min⁻¹ sits in the middle.

**2. Thiazolidine release, relative.** Averaged over 3 and 6 h: formaldehyde 0.78×, acetaldehyde 0.70×,
benzaldehyde 0.55×, n-butyraldehyde 0.27×, furfural 2.05× free cysteine. Apparent first-order H2S constants
at ~100 °C (6 h, mine): cysteine 2.7e-5, formaldehyde adduct 2.2e-5, acetaldehyde 1.9e-5, benzaldehyde 1.6e-5,
n-butyraldehyde 7.9e-6, furfural adduct 5.7e-5 min⁻¹.

**3. NH3 / H2S.** Cysteine 6 h: 13.4 µmol NH3 vs 9.5 µmol H2S → **1.4**; furfural adduct 18.3 / 20.5 → 0.89;
formaldehyde adduct 8.6 / 7.7 → 1.1. Cysteine's NH3 exceeds the 1:1 of β-elimination (H2S + NH3 + pyruvate) by 40 %,
so ~30 % of its NH3 comes from non-sulfide channels — consistent in sign with the engine's `r_cys_thermal` lump
(non-sulfide cysteine loss), though the magnitude here is small (1.3 % N in 6 h).

**4. Can the thiazolidine H2S come only through free cysteine?** For the formaldehyde adduct (thioproline),
ring-opening equilibrium lies far to the ring (my qualitative knowledge, K not in this paper). Releasing 0.8× the
H2S of free cysteine through a tiny free-cysteine fraction would need that fraction to be ~80 % — implausible unless
the gas sweep strips formaldehyde/acetaldehyde and pulls the ring open. n-Butyraldehyde (also volatile, b.p. 75 °C)
gives the **lowest** release, which argues against stripping as the whole story. **Either thiazolidines release H2S
directly, or the sweep opens them; this paper cannot tell which** (~50/50).

## 4. Kinetic numbers the repository can use

Shared conditions: 0.1 mol/L (S basis) in 10 mL water, unbuffered (pH ~4-5.5, not measured), oil bath 120 °C under a
CO2 or N2 sweep (solution ~100 °C, inferred), single determinations, sulfide recovery 95 ± 5 %.

| step / observable | quantity | value | unit | evidence class |
|---|---|---|---|---|
| Cys → H2S (neat) | cumulative H2S at 0.5/1/3/6 h | 0.056 / 0.138 / 0.356 / 0.955 | % of S | **measured time series** (Table I) |
| Cys → H2S | apparent first-order k (6 h) | **2.7e-5** | min⁻¹ | derived; T assumed 100 °C |
| thiazolidine (aldehyde) → H2S | relative to free Cys | **0.27-0.81** | — | within-study ratio, 4 aldehydes |
| furfural thiazolidine → H2S | relative to free Cys | **2.0-2.2** (fresh), **6.0** (aged, 3 h) | — | within-study ratio |
| Cys → NH3 | cumulative NH3 at 3/6 h | 0.75 / 1.34 | % of N | measured (Table II) |
| thiazolidine → NH3 | relative to free Cys (6 h) | 0.29-0.65; furfural 1.36 | — | within-study ratio |

No Ea (one temperature), no pH series.

## 5. What this can and cannot anchor

**Can:**
1. **A hold-out on `r_cys_h2s` at 100 °C** (time series, four points), independent of Zheng & Ho — if the solution
   temperature is accepted as ~100 °C. It passes already (0.67-1.05×); it is a check, not a fit row.
2. **A constraint on any thiazolidine H2S channel:** under the engine's present structure (TTCA releases H2S only via
   free cysteine), a thiazolidine with a large formation constant would release ≪ free cysteine. Measured: 0.3-0.8×.
   If B59 makes TTCA a tight reservoir, this paper says the reservoir still leaks sulfide at a rate comparable to
   free cysteine at ~100 °C (with the stripping caveat).
3. **Sign for the furfural branch:** cysteine bound to furfural releases H2S faster than free cysteine.

**Cannot:**
- No pentose or deoxyosone thiazolidine (TTCA itself was not made; the HMF adduct failed to crystallise).
- No pH, no buffer, no measured solution temperature.
- No product identification beyond H2S and NH3 (no FFT, no thiazoline, no MAA).
- Single determinations; acid step at the end (10 mL HCl, 30 min sweep) may release acid-labile sulfide formed
  during heating — included in every number.

## 6. Flags

1. **Oil bath 120 °C ≠ solution temperature** (open, gas-swept, aqueous → ≤ ~100 °C). The citation in
   `LITERATURE_SEARCH_2026-10_ROUND3.md` ("at 120 C") should be read as "~100 °C solution". Under a 120 °C reading the
   paper contradicts Zheng & Ho by 9×; under 100 °C it confirms it.
2. **Gas sweep strips volatile aldehydes** (formaldehyde, acetaldehyde), which can pull thiazolidines open; the
   relative numbers are protocol-specific.
3. **Aged furfural adduct (II) releases 3× more than fresh (I)** — the furfural thiazolidine decomposes on the shelf;
   fresh-sample values are themselves possibly affected.
4. Cystine releases more H2S than cysteine per S at 3-6 h but less at 0.5-1 h (lag) — not discussed by the authors.
5. The NH3 assay ends with NaOH aeration: anything base-labile that yields NH3 is included.
6. Registry: paper absent from `data/keys/papers.yml`; the thiazolidines are not engine species.

## Use in the model

- **Gap (b), de Roos path.** The usable facts are relative: thiazolidines of simple aldehydes release H2S at
  **0.27-0.81×** and NH3 at **0.29-0.65×** the free-cysteine rate in hot water (~100 °C, pH ~4-5.5). For B59's TTCA
  reservoir, this argues for either (i) a direct TTCA → H2S channel at a declared fraction of `k_cys_h2s`, or (ii) a fast
  ring-opening so that release through free cysteine is not suppressed by the binding constant — and it argues
  against a TTCA that withholds sulfide entirely. It does **not** measure TTCA (pentose adduct), so any number carried
  over is a declared assumption with a band of ~0.3-1× free cysteine.
- **Gap (a).** Not addressed: no dicarbonyl, no Strecker, no MAA.
- **Gap (e), furfural.** The cysteine-furfural thiazolidine releases **~2×** the H2S of free cysteine (fresh), more
  when aged. The engine has no Cys + furfural adduct; Yaghmur 2005's thiazolidine products (B, D) are the same
  family. This is a sign-level hint that furfural + cysteine can **feed** the H2S pool rather than only consume it.
- **Validation check that already passes:** free-cysteine H2S at ~100 °C matches `cysteine_thermolysis_k` within
  0.67-1.05× at four times (pH 3-5 pairs). Candidate **pre-registered hold-out** for the cysteine thermolysis channel, with
  the solution temperature declared as 100 °C (band 95-100 °C) and pH 4.5 ± 0.5.
- **Caveats:** unbuffered, pH unmeasured; solution temperature inferred; single determinations; volatile-aldehyde
  stripping; 1960 colorimetry (methylene blue, Nessler).
