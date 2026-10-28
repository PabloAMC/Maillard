# Danon, Marcotullio & de Jong 2014 — EXTRACTION (critical review: pentose → furfural in aqueous media with homogeneous catalysts — mechanism (Figs. 2-7), furfural-loss Arrhenius parameters (Table 1, 8 studies, 140-240 °C, strong acid or plain water), pentose-dehydration rate constants (Table 2, ~35 rows, 80-400 °C), biomass studies (Table 3))

### A BIOREFINERY REVIEW: EVERYTHING IS AT pH ≲ 2 (50 mM-5 M mineral acid) OR IN UNBUFFERED HOT WATER (180-400 °C). No buffered datum at pH > 3. Its answer to the model's question is nevertheless clear: **in the literature, furfural's DECAY barrier is BELOW pentose → furfural FORMATION's**: decay Ea 48-125 kJ/mol (median ~84; HCl and plain water 48-68), acid-catalysed xylose dehydration 120-133 kJ/mol (and Garrett & Dvorchik 1969, in the corpus: 127-129 vs 89). The engine has it the other way round (decay 143.5 on the carbonyl-sink family vs formation 64 lumped). The engine's `k_fur_decay` at 145 °C (0.38 min⁻¹) is also **70-1 700× faster** than any pure-furfural thermal decay in Table 1 extrapolated to 145 °C, so it is not, and cannot be read as, furfural's own thermal decay.

**Source on disk:** `data/articles/Danon2014.pdf` (16 pp., RSC publisher PDF, CC-BY 3.0, clean text
layer). Tables 1-3 and prose from `pdftotext -layout` (Table 2's exponents checked against the layout);
Fig. 6 read from a 150 dpi render (scratchpad `lit5/dan-05.png`). No SI. **Every kinetic number here is
SECONDARY** (the review's transcription of the primary paper); none of the primaries is on disk except
Garrett & Dvorchik 1969 (`garrett1969_extraction.md`).

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanistic and kinetic aspects of pentose dehydration towards furfural in aqueous media employing homogeneous catalysis" |
| Authors | Bart Danon, Gianluca Marcotullio, **Wiebren de Jong** (corresponding, wiebren.dejong@tudelft.nl) — TU Delft Energy Technology; Marcotullio at SEA Servizi Energia Ambiente, L'Aquila |
| Venue | Green Chemistry **16**, 39-54 (2014), Critical Review. Received 9 Jul 2013, accepted 15 Oct 2013, published online 15 Oct 2013 |
| DOI as printed | **10.1039/c3gc41351a** ✔ matches the brief. **Identity verified.** |
| Related primaries named in the brief | Danon, van der Aa & de Jong 2013, Carbohydr. Res. 375:145 (ref. 56: furfural degradation ± **glucose**, 50 mM HCl + 500 mM NaCl, 160-200 °C); Marcotullio & de Jong 2010 Green Chem. 12:1739 (ref. 39); Marcotullio 2011 PhD thesis (ref. 47) |

## 1. Why it matters

`docs/guides/SULFUR_FLUX_MAP.md` sec. 5 names, at ~55 % confidence, **furfural's decay barrier** as the
suspect for the FFT/MFT mis-weighting: `r_fur_decay` shares the carbonyl-sink barrier (143.5 kJ/mol,
fitted, unidentified) while furfural formation `r_tdp_fur` runs on the lumped 64 kJ/mol, so furfural
survives too long as T falls. The brief asks: is decay's barrier above or below formation's in the
literature, are there data at pH > 3, and what about the second-order furfural + sugar loss.

## 2. What the review covers (scope and conditions)

- Aqueous, homogeneous catalysis only (no heterogeneous, no ionic liquids). Catalysts: H2SO4, HCl,
  HCOOH, fumaric/maleic acid, NaOH (one row), halide salts, FeCl3, CrCl3.
- Furfural-loss studies: **140-240 °C**, initial furfural 10-160 mM.
- Pentose → furfural: **80-400 °C**; only Garrett & Dvorchik 1969 (ref. 30) goes below 100 °C (60-80 °C,
  0.125-1 M HCl).
- **pH > 3:** none buffered. The uncatalysed hot-water rows (Jing 2007, Kim 2011, Aida 2010, Gairola 2012,
  Oefner 1992, Kootstra 2009) start near neutral but are unbuffered and acidify by organic acids (Oefner
  identifies pyruvic, formic, glycolic, lactic, acetic acids); Lamminpää 2012 estimates [H+] ≈ 50 mM at
  200 °C in its formic-acid runs. One alkaline row (Oefner, 100 mM NaOH, 220 °C).

## 3. Table 1 (verbatim) — furfural loss, first-order unless stated; k = A·exp(−E/RT), A in s⁻¹

| ref. | primary | [F]₀ (mM) | catalyst (mM) | T (°C) | ln A (s⁻¹) | E (kJ/mol) |
|---|---|---|---|---|---|---|
| 53 | Williams & Dunlop 1948 | 10-20 | H2SO4 (50) | 150-210 | 12.36 | 83.6 |
| 78 | Root et al. 1959 | 47.8 | H2SO4 (3.1-400) | 160-240 | [H+]·17.15 | 92.3 |
| 54 | Rose et al. 2000 | 10 | HCl (100) | 150-169 | 3.07 | 48.1 |
| 52 | Marcotullio et al. 2009 | 60.4-72.5 | H2SO4 (36.4-145.5) | 150-200 | a_H+·26.64 | 125.1 |
| 79 | Lamminpää et al. 2012 | 50 | HCOOH (50-80) | 160-200 | uncatalysed 13.82 | 75.5 |
| | | | | | acid-catalysed [H+]·31.44 | 135.0 |
| 80 | Weingarten et al. 2010 | 160 | HCl (100) | 140-160 | 8.43 | 67.6 |
| 81 | Jing & Lu 2007 | 34 | **none** (water, 10 MPa) | 180-220 | 7.59 | 58.8 |
| 56 | Danon et al. 2013 | 50 | HCl (50) + NaCl (500) | 160-200 | 16.84 | 102.1 |

The review's own reading (p. 46): reaction order often reported ≠ 1 but first order fits; "the highest
[E] values are obtained in studies where H2SO4 was used"; uncatalysed and HCl values "in the same range,
and ... lower than those with H2SO4". Root and Williams & Dunlop both find the rate rises with initial
furfural concentration (higher order; Root's mechanism adds a second-order furfural + intermediate
term). Products: formic acid + resinous tar (Williams & Dunlop); Marcotullio 2009 saw no formic acid.

### Arithmetic (mine): Table 1 evaluated at the engine's temperatures (min⁻¹; scratchpad `lit5/danon_tab1.py`, run in the maillard env)

| source (catalyst as evaluated) | E | k(100 °C) | k(145 °C) | k(180 °C) |
|---|---|---|---|---|
| Williams & Dunlop (50 mM H2SO4) | 83.6 | 2.8e-5 | 5.1e-4 | 3.2e-3 |
| Root ([H+] = 0.1 M) | 92.3 | 2.0e-5 | 5.0e-4 | 3.9e-3 |
| Rose (100 mM HCl) | 48.1 | 2.4e-4 | 1.3e-3 | 3.7e-3 |
| Marcotullio 2009 (a_H+ = 0.1) | 125.1 | 6.9e-6 | 5.3e-4 | 8.5e-3 |
| Lamminpää, uncatalysed term | 75.5 | 1.6e-3 | 2.2e-2 | 1.2e-1 |
| Lamminpää, acid term ([H+] = 0.05 M) | 135.0 | 1.7e-5 | 1.9e-3 | 3.7e-2 |
| Weingarten (100 mM HCl) | 67.6 | 9.5e-5 | 9.9e-4 | 4.4e-3 |
| Jing (water, no catalyst) | 58.8 | 7.0e-4 | 5.4e-3 | 2.0e-2 |
| Danon 2013 (HCl + NaCl) | 102.1 | 6.3e-6 | 2.2e-4 | 2.1e-3 |
| **engine `k_fur_decay`** (B60, pH-free) | **143.5** | **2.6e-3** | **0.38** | **9.3** |
| engine `k_tdp_fur` (3-deoxypentosone → furfural) | 64.1 | 3.5e-5 | 3.2e-4 | 1.3e-3 |

Most of these are extrapolations below the measured range (100 and 145 °C are below every study's window
except Weingarten's 140-160). Lamminpää's uncatalysed term exceeds its own acid term at every T, which
looks like a units/definition issue in the transcription — **treat that row as suspect**.

## 4. Table 2 (key rows, verbatim) — pentose disappearance, R1 (→ furfural) + R2 (→ loss), k in s⁻¹

| ref. | pentose (mM) | catalyst (mM) | T (°C) | k1 | k2 | max yield % |
|---|---|---|---|---|---|---|
| 81 Jing | Xyl 72 | none | 200 | 1.98e-4 | 2.05e-4 | ~50 |
| 84 Kim | Xyl 200 | none | 200 | 4.30e-6 | [X]·7.73e-9 (k2a), 6.39e-6 (k2b) | ~52 |
| 86 Gairola | Ara 14-350 | none | 200 | 1.56e-4 | 1.74e-4 | 8 |
| 87 Kootstra | Ara 33 | none / FA / MA / H2SO4 (50) | 170 | k1+k2 = 5.81e-3 / 4.52e-3 / 5.81e-3 / 1.56e-2 | | n.d. |
| 30 Garrett | Rib / Xyl / Lyx / Ara 0.1-1 | HCl (1000) | 80 | 1.67e-6 / 6.99e-7 / 6.77e-7 / 3.12e-7 | 5.45e-7 / 1.12e-7 / 2.92e-7 / 0.28e-7 | 48 / 38 / 35 / 26 |
| 50 Oefner | Xyl 67 | none / H2SO4 (50) / NaOH (100) | 190 / 190 / 220 | k1+k2 = 2.75e-4 / 4.05e-4 / 1.03 | | ~13 / ~52 / ~49 |
| 39 Marcotullio 2010 | Xyl 35 | HCl 50 + NaCl 340 | 200 | 5.28e-3 | 1.17e-3 | 73 |
| 79 Lamminpää | Xyl 20-80 | HCOOH (50) | 200 | 6.85e-4 | 2.28e-4 | ~65 |
| 82 Schoenemann | Xyl 84-350 | HCl (277) | 150 | 3.8e-4 | 3.0e-5 + [X]^1.16·[F]^0.5 term | ~46 |
| 80 Weingarten | Xyl 740 | HCl (100) | 170 | 6.11e-4 | [F]·2.04e-3 | 30 |

(Full table: 35 rows, also high-salt, FeCl3, KBr/KI, supercritical 350-400 °C, biphasic.) **Table 2 gives
rate constants at one temperature per row, not Arrhenius parameters.**

### Activation energies for FORMATION quoted in the text (p. 43 and p. 49)

| step | E (kJ/mol) | conditions | primary |
|---|---|---|---|
| acid-catalysed xylose dehydration | **133** | dilute HCl/H2SO4 | refs. 36, 39, 47 (Marcotullio thesis) |
| xylose degradation (k1 + k2), water | **119.4** | 180-220 °C | Oefner 1992 |
| same, H2SO4 0.5-50 mM | **120.7** | 180-220 °C | Oefner 1992 |
| same, 100 mM NaOH | **63.7** | 180-220 °C | Oefner 1992 |
| xylose → xylulose isomerisation | 96.7 | aqueous CrCl3 | Choudhary 2012 |
| xylulose dehydration | 64.9 | aqueous HCl | Choudhary 2012 |
| (corpus cross-check) pentose → furfural | 127-129 (ΔH‡) | 0.125-1 M HCl, 60-80 °C | Garrett & Dvorchik 1969 (`garrett1969_extraction.md`) |
| (corpus cross-check) furfural acid solvolysis | 89 (21.3 kcal) | same | Garrett & Dvorchik 1969 |

## 5. Second-order furfural loss terms (as summarised)

- **Danon 2013 (ref. 56)**: furfural 50 mM, 50 mM HCl + 500 mM NaCl, 160-200 °C, **± glucose** (a hexose,
  not a pentose). Pure furfural fits first or second order equally; **with glucose, "clear second-order
  reaction kinetics"**, attributed to a **Diels-Alder** reaction. **Only the first-order pure-furfural
  parameters are tabulated (ln A 16.84, E 102.1); the second-order furfural + glucose constants are NOT
  given in the review** — the primary (Carbohydr. Res. 375:145) is needed.
- Schoenemann & Hofmann 1957: "Kondensation" of furfural with a xylose-derived intermediate, orders
  1.16 in xylose and 0.5 in furfural (150 °C, 277 mM HCl).
- Root 1959: same two-route structure (self-degradation + furfural × intermediate); Garrett 1969 (in
  corpus) found **no** effect of added furfural on furfural yield from ribose at 80 °C, 1 M HCl.
- Weingarten 2010: xylose loss ∝ [F]·2.04 × 10⁻³ (170 °C, 100 mM HCl).
- Lamminpää 2012: compared mechanisms; the xylose-alone loss route beats the intermediate + furfural
  route.

## 6. Mechanism (Figs. 2-7, the proposed synthesis in Fig. 6)

All routes start from **acyclic** xylose (1); the 2,5-anhydro-furanose route (Fig. 4, Nimlos/Qian) is
judged weaker. Three entries, coexisting, weighted by conditions:
- **1,2-enolisation** (H⁺, halide-assisted) → 1,2-enediol (2) → dehydration → **2,3-unsaturated aldehyde
  (4)** → −2 H2O → furfural. (4) is in tautomeric equilibrium with **(5) = 3-deoxypentosulose**, i.e. the
  engine's TDP (3-deoxypentosone).
- **β-elimination** directly from (1) to (4) (Hurd 1932), halide-assisted.
- **1,2-hydride shift** isomerisation to xylulose (3) at pH < 1 or with M³⁺ (Cr, Al, Fe), then dehydration.
- From the ketose, **2,3-enolisation → 2,3-enediol (3e)** is drawn as an **irreversible exit to side
  products** (in Maillard terms, the 1-deoxyosone direction the engine calls DPO), not to furfural.
- Rate-limiting step: enolisation/isomerisation (no solvent H in furfural's C-H in strong acid; partial
  D-exchange at pD 4.5 → enolisation reversible at higher pD). Enolisation is base-favoured, dehydration
  acid-favoured, so "furfural formation from aldoses cannot be considered a simple ... specific acid
  catalysis". Side products: organic acids, acetaldehyde, formaldehyde, crotonaldehyde, pyruvaldehyde
  (Antal 1991), reductic acid (from furfural), humins.
- Pentose reactivity: ribose > xylose ≈ lyxose > arabinose (Garrett); xylose ~30 % faster than arabinose
  at 200 °C (Jing vs Gairola).

## 7. What is NOT reported

- No buffered data, no data at pH 3-7, no amine, no thiol, no H2S. No Maillard conditions at all.
- No Arrhenius parameters for formation in a table (text values only), no uncertainty on any E.
- The Danon 2013 second-order (furfural + glucose) rate constants.
- No furfural decay with amino acids (the dominant sink in Maillard pots).

## 8. Flags

1. Inter-study agreement is "merely in the order of magnitude" (review's own words); e.g. Jing vs Kim in
   plain water at 200 °C differ 46× in k1.
2. Lamminpää's uncatalysed ln A looks inconsistent with its acid term (sec. 3 arithmetic).
3. Oefner's NaOH row: k1 + k2 = 1.03 s⁻¹ with "~49 %" max furfural yield in 100 mM NaOH is chemically
   odd (alkaline xylose gives acids, not furfural); possibly a transcription issue — do not use.
4. Table 1 and 2 values are secondary; cite the primary when any of them enters a parameter file.

## 9. Use in the model

- **Barrier ordering — the answer to the question:** across the review, furfural decay E (48-125,
  median ~84 kJ/mol; HCl/water 48-68) is **lower** than acid-catalysed pentose → furfural formation
  (120-133 kJ/mol), by ~40-80 kJ/mol; Garrett 1969 in the corpus agrees (89 vs 127-129). Only
  Lamminpää's split acid term (135) exceeds the formation values. **The engine inverts this** (decay 143.5
  on the shared carbonyl-sink barrier, formation `r_tdp_fur` 64). This **supports, does not contradict,**
  `SULFUR_FLUX_MAP.md` sec. 5's suspicion that a too-high furfural decay barrier lets furfural survive
  as T falls. Transfer caveats: the literature steps are pentose → furfural overall (the engine's
  `r_tdp_fur` starts from TDP, after the enolisation), at pH 0-2, without amines; the engine's lumped
  decay includes amine/thiol consumption, whose barrier (Yaghmur 2005: 46.5 kJ/mol with cysteine, pH 5,
  40-70 °C) is even lower. Every independent line puts furfural's loss barrier below ~100 kJ/mol.
- **Magnitude:** the engine's `k_fur_decay` (0.38 min⁻¹ at 145 °C) is 70× the plain-water value (Jing,
  extrapolated) and 200-1 700× the strong-acid values. The engine's furfural sink is therefore
  dominated by non-thermal chemistry (amines, thiols, H2S) if it is right at all. A split of the furfural
  decay barrier from the carbonyl-sink family with a declared band of ~50-100 kJ/mol (Yaghmur's 46.5 at
  the bottom; Table 1's water/HCl 48-68 and H2SO4 84-125 above) would be grounded in this table. **No
  datum here can anchor its pH dependence** (pH > 3 absent).
- **Mechanism mapping:** Fig. 6's (5) = TDP ⇌ (4) → furfural matches the engine's PENT → TDP → FUR, and
  its 3e exit = the 2,3-enolisation (DPO) branch. The review supports treating the enolisation as rate
  limiting and base-favoured with an acid-favoured dehydration — i.e. the engine's single `"acid"` pH
  factor on `r_pent_tdp` is a simplification that this review explicitly says is not "simple specific acid
  catalysis".
- **Second-order furfural + sugar loss:** exists (Danon 2013 with glucose; Schoenemann/Root with a
  xylose intermediate) but no usable constant is printed here; Garrett 1969 found none at 80 °C. Low
  priority for aqueous pots at ≤ 0.1 M sugar.
