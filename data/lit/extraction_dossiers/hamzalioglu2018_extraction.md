# Hamzalıoğlu & Gökmen 2018 — reactions of HMF with amino and thiol groups of amino acids

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THIS IS THE CORPUS'S ONLY HMF-SINK SOURCE, AND IT ALSO FEEDS THE SULFUR LANE'S HMF–THIOL CROSSLINK.

---

## §0. IDENTITY — ★ THE ROUND-3 LEAD IS CONFIRMED AND ITS AUTHORS ARE NOW KNOWN

`research_round3_channels.md` §A.3.5 recorded an unattributed lead:

> *"Investigation and kinetic evaluation of the reactions of hydroxymethylfurfural with amino and
> thiol groups of amino acids"*, Food Chem., ScienceDirect PII `S0308814617312852`.
> ⚠️ **I did not verify its authors — do not attribute it.** … its identity must be Crossref-verified
> before it is ordered."

**It is on disk, and the identity is now established from the document itself:**

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/Hamzalioglu2018.pdf` (749,637 bytes) | `ls` |
| title | *"Investigation and kinetic evaluation of the reactions of hydroxymethylfurfural with amino and thiol groups of amino acids"* ✔ **exact match to the round-3 lead** | p.1 |
| **authors** | ★ **Aytül Hamzalıoğlu, Vural Gökmen** — Food Quality and Safety (FoQuS) Research Group, Department of Food Engineering, Hacettepe University, 06800 Beytepe, Ankara, Turkey | p.1 |
| **DOI** | ★ **`10.1016/j.foodchem.2017.07.131`** | printed on p.1, footer |
| journal | ***Food Chemistry* 240 (2018) 354–360** | running head on all 7 pages |
| dates | Received 29 May 2017 · Revised 13 Jul 2017 · Accepted 25 Jul 2017 · online 26 Jul 2017 | p.1 footer |
| version | published version of record, `0308-8146/© 2017 Elsevier Ltd.` | p.1 |
| funding | ★ **none** — *"This research did not receive any specific grant from funding agencies"* | Acknowledgement |
| competing interests | *"There are no conflicts of interest to declare."* | p.6 |

**⇒ The round-3 lead resolves to a Gökmen-school paper, so it belongs to the same methodological
family as the five multiresponse papers of this cluster. The PII `S0308814617312852` and the DOI
`10.1016/j.foodchem.2017.07.131` are the same object; the round-3 caution ("do not attribute") can
be lifted.**

---

## §1. ★ WHY THIS PAPER MATTERS TWICE OVER

**(a) THE HMF SINK.** Round-3 §A.3.5 states the problem exactly: *"The repo needs an HMF sink and
currently has none — Han 2025's k15 is the only published one and it is unidentifiable."* Since
then, wave K5a has examined the five Gökmen multiresponse papers and found:

| paper | HMF sink in its model | usable? |
|---|---|---|
| Kocadağlı JAFC 2016 (glucose ± NaCl) | `HMF → P6`, first order | **NO** — Ea HPD (±154.4) > estimate; k = 0 at 160 °C in the NaCl system |
| Kocadağlı Food Chem 2016 (wheat flour) | `HMF → P8`, first order | **NO** — k = 0 at 180 °C, HPD 71 % at 160 °C; and the bimolecular amine sink was **rejected** |
| Göncüoğlu Taş 2017 (hazelnut) | `HMF → P7`, first order | **partly** — usable at 150 and 160 °C only |
| Gürsul Aktağ 2020 (fruit juice) | ★ **none at all** | — |
| Şen 2022 (nuts/seeds) | `HMF + AA → P3`, **bimolecular** | **NO** — zero or HPD-dominated in 6 of 10 cells |
| Han 2025 (round-3 §A.3.1) | `HMF → melanoidin`, first order | **NO** — unidentifiable (HPD 50–95 % of estimate; non-monotonic in T) |

**⇒ Of seven published HMF sinks, none survives. This paper is the eighth, it is a dedicated study
rather than a by-product of a network fit, and its constants DO survive — at 5–50 °C.**

**(b) THE SULFUR LANE.** The paper measures **HMF + cysteine** explicitly and finds it the fastest
of the three amino acids **because of its thiol group**. That is a direct, quantified
HMF–thiol crosslink for the sulfur lane, and it is the only one in the corpus.

---

## §2. SYSTEM AND CONDITIONS `[M]`

| item | value | anchor |
|---|---|---|
| **HIGH-MOISTURE model** | **8 µmol HMF + 30 µmol amino acid in 1.5 mL of 0.05 % benzoic acid solution** | §2.2 |
| **pH (high moisture)** | ★ **3.5** — *"The pH of these model systems was measured as 3.5."* | §2.2 |
| amine excess | **30 : 8 = 3.75-fold molar excess of amino acid over HMF**, deliberately: *"Since HMF content of foods is minor compared to amino acids and not to limit the possible reactions of HMF with amino acids, excess amount of amino acids was added"* | §2.2 |
| benzoic acid | **antimicrobial only** — *"used in order to prevent microbial growth during the reaction period"* | §2.2 |
| **LOW-MOISTURE model** | **1 g ground roasted coffee + 8 µmol arginine, cysteine or lysine** in a glass flask | §2.2 |
| coffee preparation | green coffee beans roasted in an oven at **220 °C for 10 min**, chosen because *"it cause to the formation of HMF at highest amounts in roasted coffee (Kocadağlı, Göncüoğlu, Hamzalıoğlu, & Gökmen, 2012)"* | §2.2 |
| **initial HMF in coffee** | ★ **1.706 ± 0.061 µmol** (per 1 g roasted coffee) | §3.1 |
| controls | **"HMF"** (no amino acid) and **"Coffee"** (no amino acid) — *"to monitor the self-degradation of HMF during storage"* | §2.2 |
| **temperatures** | ★ **5, 25 and 50 °C**, water bath — *"selected to reveal the reactions of HMF with amino acids at cold and high storage conditions **besides the body conditions**"* | §2.2 |
| **times** | **1, 2, 3, 5 and 7 days** (5 points) | §2.2 |
| amino acids | **L-arginine, L-cysteine (>99 %), L-lysine** (Sigma-Aldrich) | §2.1 |
| HMF analysis | Shimadzu HPLC, Atlantis dC18 250×4.6 mm 5 µm, isocratic **10 mM aqueous formic acid / acetonitrile 90:10 v/v**, 1.0 mL/min, 25 °C, DAD **285 nm**; calibration **0.5–50 µg/mL** (0.5, 1, 2, 5, 10, 50) | §2.3 |
| adduct confirmation | **Thermo Q Exactive Orbitrap HRMS**, positive ESI, **R = 100 000**, m/z 50–600, mass accuracy **Δ < 1 ppm** | §2.4 |
| statistics | one-way ANOVA, Duncan test, SPSS 17.0, **p < 0.05** | §2.5 |
| fitting | **MATLAB Curve Fitting Toolbox**, non-linear regression; **R² between 0.80 and 0.99** | §3.3 |

⚠️ **Note the pH match.** The high-moisture system is at **pH 3.5** — the **same pH as Ağçam 2022's
juice mediums (§0 of that dossier) and within 0.1 of Gürsul Aktağ 2020's juices (pH 3.4).**
★ **Three independent HMF papers in this cluster sit at pH 3.4–3.5. That is the corpus's densest
pH agreement anywhere in the HMF lane.**

---

## §3. THE KINETIC MODEL — verbatim from §3.3

> "In order to determine the rate constants, the reaction of HMF with amino acids was considered to
> proceed as;
>
> `A + B --k--> C`   (1)
>
> A represents HMF, B represents amino acid with amino and/or thiol compounds, and C represents
> HMF-amino acid adducts. … The differential rate equation … could be described as follows;
>
> `d[A]/dt = −k[A][B]`   (2)
>
> In both model systems, amino acid concentration is excess in reaction medium compared to HMF.
> Pseudo-first order, although the actual reaction mechanism is bimolecular, can happen when one of
> the reactants is excess (van Boekel, 2008). For the reaction given in Eq. (1), B is excess and
> **the rate constant k′ = k[B] is constant** as long as [B] does not change significantly and
> thereby reaction is considered as pseudo-first order. The Eq. (2) given becomes as follows;
>
> `d[A]/dt = −k[A][B] = −k′[A]`   (3)"

**★ `k′ = k·[B]`. The printed constants are PSEUDO-first-order and carry the amino-acid
concentration inside them.** To recover a second-order `k`, divide by `[B]`. In the high-moisture
system `[B] = 30 µmol / 1.5 mL = 20 mM`, so `k = k′/0.020 M`. **In the low-moisture coffee system
`[B]` is 8 µmol in 1 g of ground coffee with no defined volume, so NO second-order constant can be
recovered there.**

**⚠️ AND THE AUTHORS DECLARE THE PSEUDO-FIRST-ORDER ASSUMPTION BROKEN AT 50 °C**, verbatim:

> "The reaction rate constants were calculated as pseudo-first order reaction rate constants due to
> the assumption of constant amino acid concentration compared to HMF. **However, since some part of
> amino acids might be degraded at 50 °C, measured rate constants become apparent pseudo-first order
> rate constants that also bear the change in amino acid concentration in time in the case of
> reactions at 50 °C.**"

★ **This matters enormously for the Arrhenius fit**, because — as §6 shows — **the 50 °C point
carries essentially all of the temperature signal.**

---

## §4. TABLE 1 — the pseudo-first-order rate constants

**Transcribed verbatim from p.5.** Caption verbatim: *"Kinetic rate constants for the reaction of
5-hydroxymethylfurfural with amino acids in both high moisture ("HMF-Arg", "HMF-Cys" and "HMF-Lys")
and low moisture ("Coffee-Arg", "Coffee-Cys", "Coffee-Lys") model systems."*
Column header verbatim: **`k′ (×10⁻¹ day⁻ı)`** (the `ı` is a dotless-i typo for `1`).
Provenance **`[F]`**.

### HIGH MOISTURE (pH 3.5, 8 µmol HMF + 30 µmol amino acid in 1.5 mL)

| T (°C) | HMF-Arg | HMF-Cys | HMF-Lys |
|---|---|---|---|
| 50 | 0.55 ± 0.03 | ★ **4.65 ± 0.23** | 0.160 ± 0.10 |
| 25 | 0.15 ± 0.01 | 1.03 ± 0.05 | 0.090 ± 0.05 |
| 5 | 0.12 ± 0.01 | 0.79 ± 0.04 | 0.088 ± 0.04 |

### LOW MOISTURE (roasted coffee, 1 g + 8 µmol amino acid)

| T (°C) | Coffee-Arg | Coffee-Cys | Coffee-Lys |
|---|---|---|---|
| 50 | 1.18 ± 0.06 | 1.23 ± 0.06 | 1.22 ± 0.06 |
| 25 | 1.01 ± 0.05 | 1.03 ± 0.05 | 0.84 ± 0.04 |
| 5 | 0.97 ± 0.05 | 0.61 ± 0.03 | 0.58 ± 0.03 |

**In `day⁻¹` (multiply by 0.1):**

| system | 5 °C | 25 °C | 50 °C |
|---|---|---|---|
| HMF-Arg | 0.012 | 0.015 | 0.055 |
| **HMF-Cys** | **0.079** | **0.103** | ★ **0.465** |
| HMF-Lys | 0.0088 | 0.0090 | 0.0160 |
| Coffee-Arg | 0.097 | 0.101 | 0.118 |
| Coffee-Cys | 0.061 | 0.103 | 0.123 |
| Coffee-Lys | 0.058 | 0.084 | 0.122 |

**★ SECOND-ORDER `k = k′/[B]` for the HIGH-MOISTURE system only, `[B] = 0.020 M` `[D]`:**

| system | k(5 °C) | k(25 °C) | k(50 °C) | units |
|---|---|---|---|---|
| HMF + Arg | 0.60 | 0.75 | 2.75 | **M⁻¹ day⁻¹** |
| **HMF + Cys** | **3.95** | **5.15** | **23.3** | **M⁻¹ day⁻¹** |
| HMF + Lys | 0.44 | 0.45 | 0.80 | **M⁻¹ day⁻¹** |

(In SI: `k(HMF+Cys, 25 °C) = 5.15 M⁻¹ day⁻¹ = 5.96 × 10⁻⁵ M⁻¹ s⁻¹`.)
★ **These nine numbers are the first genuine second-order HMF-sink constants the corpus has ever
held.** Derived `[D]` — the authors never print them.

### 4.1 ⚠️ THE HMF-Lys SD COLUMN IS MIS-SCALED BY ~10×

Every other row in Table 1 has an SD of **exactly 5 %** of the estimate (0.55±0.03, 4.65±0.23,
0.15±0.01, 1.03±0.05, 0.12±0.01, 0.79±0.04, and all six coffee rows). The HMF-Lys row instead reads
**0.160 ± 0.10 (63 %), 0.090 ± 0.05 (56 %), 0.088 ± 0.04 (45 %)** — i.e. the SDs are ~10× larger
than the paper's own uniform pattern, and would make all three estimates statistically
indistinguishable from zero. **Almost certainly the ×10⁻¹ scaling was not applied to this column:
the intended values are ±0.010 / ±0.005 / ±0.004, restoring ~5 %.** Transcribed verbatim above.
⚠️ **As printed, the HMF-Lys row must be REFUSED; under the 10× correction it is usable. Flag the
ambiguity rather than choosing.**

---

## §5. ACTIVATION ENERGIES AND PRE-EXPONENTIALS `[F]`

The Arrhenius section, verbatim §3.4:

> "The plots of ln k versus 1/T were shown in Fig. 4 for the reaction of HMF in both high and low
> moisture model systems. From the slopes of linear fitting lines, **pre-exponential factor A and
> activation energies were calculated. Pre-exponential factor A was 612.59, 23980.59 and 1.63
> day⁻¹ for "HMF-Arg", "HMF-Cys" and "HMF-Lys" model systems whereas, this value was 2.55, 12.08
> and 9.67 day⁻¹ for "Coffee-Arg", "Coffee-Cys" and "Coffee-Lys" model systems, respectively. The
> activation energies were calculated as 25.457, 29.639 and 10.018 kJ mol⁻¹ for the reactions of HMF
> with arginine, cysteine and lysine in high moisture model systems ("HMF-Arg", "HMF-Cys" and
> "HMF-Lys"), respectively. Similarly, the activation energies were calculated as 3.275, 11.55 and
> 12.33 kJ mol⁻¹ for the reactions of HMF with arginine, cysteine and lysine in low moisture model
> systems ("Coffee-Arg", "Coffee-Cys", "Coffee-Lys"), respectively.**"

| system | **Ea (kJ/mol)** | **A (day⁻¹)** |
|---|---|---|
| HMF-Arg | 25.457 | 612.59 |
| HMF-Cys | 29.639 | 23980.59 |
| HMF-Lys | 10.018 | 1.63 |
| Coffee-Arg | 3.275 | 2.55 |
| Coffee-Cys | 11.55 | 12.08 |
| Coffee-Lys | 12.33 | 9.67 |

The authors' own commentary, verbatim:

> "However, since calculated activation energies for the reactions of HMF with amino or thiol groups
> **ranged between 3 and 30 kJ mol⁻¹, temperature dependence of these reactions was found as
> relatively low** … Particularly, compared to reactions of HMF with arginine and cysteine,
> reactions of HMF with lysine are less affected by the temperature increase. Considering the
> reaction conditions, **comparably higher activation energies were found in high moisture
> conditions than low moisture conditions. This means that limited water in reaction system
> restricted the temperature dependence of the reaction.** … **The activation energies of chemical
> reactions in foods are generally reported ranging between 50 and 150 kJ mol⁻¹ (van Boekel, 2008).**"

**⇒ THE AUTHORS THEMSELVES FLAG THAT THEIR Ea ARE 2–50× BELOW THE FOOD-CHEMISTRY NORM.**

---

## §6. ★★ WAVE K5a CROSS-CHECK — I refitted every Arrhenius line from Table 1. THE Ea ALL REPRODUCE; **FOUR OF THE SIX PRE-EXPONENTIALS DO NOT.**

**Method `[D]`.** OLS `ln k′` vs `1/T` over 278.15 / 298.15 / 323.15 K, with `k′` in `day⁻¹`.

| system | **Ea refit** | Ea published | **A refit (day⁻¹)** | A published | **R²** | verdict |
|---|---|---|---|---|---|---|
| HMF-Arg | **25.489** | 25.457 | **615.82** | 612.59 | 0.872 | ✔ **both reproduce** |
| HMF-Cys | **29.675** | 29.639 | **24115.1** | 23980.59 | 0.874 | ✔ **both reproduce** |
| **HMF-Lys** | **10.036** | 10.018 | **0.6156** | **1.63** | 0.795 | ⚠️ **A is WRONG — and `1/0.6156 = 1.624 ≈ 1.63`** |
| **Coffee-Arg** | **3.276** | 3.275 | **0.3926** | **2.55** | 0.909 | ⚠️ **A is WRONG — and `1/0.3926 = 2.547 ≈ 2.55`** |
| **Coffee-Cys** | **11.571** | 11.550 | **9.689** | **12.08** | 0.913 | ⚠️ **A is WRONG — it equals the published Coffee-LYS value** |
| **Coffee-Lys** | **12.343** | 12.330 | **12.115** | **9.67** | 0.9996 | ⚠️ **A is WRONG — it equals the published Coffee-CYS value** |

### 6.1 ★ TWO INDEPENDENT, DIAGNOSABLE ERRORS IN THE PRE-EXPONENTIAL ROW

**ERROR 1 — A SIGN FLIP ON EVERY NEGATIVE INTERCEPT (2 of 6 systems).**
For HMF-Lys the true intercept is `ln A = −0.4852`, giving `A = 0.6156`. The paper prints **1.63 =
e^{+0.4852}**. For Coffee-Arg the true intercept is `ln A = −0.9349`, giving `A = 0.3926`; the paper
prints **2.55 = e^{+0.9349}**. **These are exactly the two systems whose fitted `A` falls below 1
day⁻¹ — i.e. the two whose `ln A` is negative. The published value is the reciprocal of the correct
one, in both cases, to three significant figures.** A spreadsheet that took `exp(|intercept|)`
reproduces this precisely.

**ERROR 2 — THE COFFEE-CYS AND COFFEE-LYS PRE-EXPONENTIALS ARE SWAPPED.**
My refit gives Coffee-Cys `A = 9.689` and Coffee-Lys `A = 12.115`. The paper prints Coffee-Cys
**12.08** and Coffee-Lys **9.67** — i.e. **each system is printed with the other's value, matching my
refits to 0.1 % once exchanged.** The corresponding `Ea` (11.55 and 12.33) are **not** swapped and
are correct as printed.

**⇒ ONLY 2 OF THE 6 PUBLISHED PRE-EXPONENTIAL FACTORS ARE CORRECT AS PRINTED.** All six activation
energies are correct. **This defect could not have been found without reproducing the fit from the
paper's own k table; it is invisible on inspection.**
★ **This is the same forensic signature the repo already carries twice** — `k3` §0.2:
*"Two constants in `arrhenius_params.yml` follow the same forensic signature — a confirmed Ea bolted
to an invented prefactor."* **Here it is the SOURCE LITERATURE doing it: a confirmed Ea bolted to a
mis-derived prefactor. Three of three audited cases now show the same pattern.**

### 6.2 ⚠️ THE THIRD DEFECT — THE TEMPERATURE SIGNAL COMES ALMOST ENTIRELY FROM ONE POINT

Look at the 5 → 25 °C step, a 20 °C rise:

| system | k(5) | k(25) | ratio over 20 °C | k(50) | ratio 25→50 |
|---|---|---|---|---|---|
| HMF-Arg | 0.012 | 0.015 | **1.25×** | 0.055 | 3.67× |
| HMF-Cys | 0.079 | 0.103 | **1.30×** | 0.465 | 4.51× |
| HMF-Lys | 0.0088 | 0.0090 | ★ **1.02×** | 0.0160 | 1.78× |
| Coffee-Arg | 0.097 | 0.101 | ★ **1.04×** | 0.118 | 1.17× |
| Coffee-Cys | 0.061 | 0.103 | 1.69× | 0.123 | 1.19× |
| Coffee-Lys | 0.058 | 0.084 | 1.45× | 0.122 | 1.45× |

**In four of six systems the rate barely moves between 5 and 25 °C (1.02–1.30×), and essentially all
the Arrhenius slope is set by the single 50 °C point** — **the very point at which the authors
declare the pseudo-first-order assumption invalid** (§3, verbatim: *"since some part of amino acids
might be degraded at 50 °C, measured rate constants become apparent pseudo-first order rate
constants that also bear the change in amino acid concentration in time"*).

**⇒ This is the Ma-2022 defect class of round-3 §A.3.2, arriving from the other direction.** There,
a flat k over 40 °C produced an artefactually tiny Ea (1.84 kJ/mol for MGO). Here, three
temperatures with two of them nearly coincident produce Ea of **3.275** (Coffee-Arg) and **10.018**
(HMF-Lys) — and the refit R² is **0.795** (HMF-Lys), **0.872** (HMF-Arg), **0.874** (HMF-Cys).
**An R² of 0.795 on a three-point line is not a fit.**

**⚠️ Ea(Coffee-Arg) = 3.275 kJ/mol is a chemically impossible barrier** — it says the HMF+arginine
reaction in roasted coffee is essentially temperature-independent from 5 to 50 °C — and it comes
from k = 0.097 / 0.101 / 0.118, a **1.22× move over 45 °C**. **REFUSE it for exactly the reason
round-3 §A.3.2 refuses Ma's Ea(MGO) = 1.84.**

### 6.3 What survives

| system | Ea | k spread over 45 °C | refit R² | verdict |
|---|---|---|---|---|
| **HMF-Cys** | 29.639 | **5.9×** | 0.874 | ★ **the best row — largest signal, thiol-specific, and the Ea reproduces** |
| HMF-Arg | 25.457 | 4.6× | 0.872 | **USE-Q** |
| Coffee-Lys | 12.330 | 2.1× | 0.9996 | **USE-Q** — highest R², but small signal |
| Coffee-Cys | 11.550 | 2.0× | 0.913 | **USE-Q** |
| HMF-Lys | 10.018 | 1.8× | **0.795** | **PRIOR-ONLY** — worst fit; and its SD column is broken (§4.1) |
| **Coffee-Arg** | 3.275 | **1.22×** | 0.909 | ★ **REFUSE — plateau artefact** |

---

## §7. ABSOLUTE DEPLETION — the measured outcome, independent of any model `[M]`

**LOW MOISTURE (roasted coffee), % of HMF depleted after 7 days:**

| system | 5 °C | 50 °C |
|---|---|---|
| Coffee (no amino acid) | **24.5 %** | **39.7 %** |
| Coffee-Arg | 52.8 % | 57.6 % |
| Coffee-Cys | 42.7 % | 58.0 % |
| Coffee-Lys | 36.1 % | 62.1 % |

*(25 °C values are in Fig. 1b, a raster; the text states only that elimination "gradually increased".)*

**HIGH MOISTURE (pH 3.5), % of HMF depleted after 7 days:**

| system | 5 °C | 25 °C | 50 °C |
|---|---|---|---|
| **HMF alone (self-degradation)** | ★ **0.9 %** | — | — |
| HMF-Arg | — | 9.7 % | 31.9 % |
| **HMF-Cys** | ★ **41 %** | **52.8 %** | ★ **97.2 %** |
| HMF-Lys | — | 6.5 % | 10.2 % |

★★ **THE HEADLINE NUMBERS FOR THE REPO:**
- **HMF self-degradation in water at pH 3.5, 5 °C, 7 days = 0.9 %.** The amine sink is therefore
  **essentially the whole sink** in this regime. Verbatim: *"elimination of HMF through its
  self-degradation might be considered as minor compared to its reactions with amino or thiol groups
  in the presence of amino acids."*
- **Cysteine alone removes 41 % of the HMF in 7 days at 5 °C, and 97.2 % at 50 °C.**
- **Cysteine at 5 °C (41 %) removes more HMF than lysine does at 50 °C (10.2 %) — a 4× inversion
  across a 45 °C temperature gap.** ⇒ **Amine identity outranks temperature by a wide margin over
  this window. Any model carrying one lumped "HMF + AA" constant is wrong by ≥6× depending on which
  amino acid dominates.**
- **In roasted coffee the no-amine control still loses 24.5–39.7 %**, i.e. the low-moisture matrix
  has its own large HMF sink (matrix binding, melanoidin incorporation, or reaction with the
  coffee's endogenous amines).

### 7.1 Ordering `[F]`, verbatim

> "Calculated pseudo-first order reaction rate constants were in the following order;
> **k_Cysteine > k_Arginine > k_Lysine for high moisture model systems.** Comparing to these rate
> constants, **the k_Cysteine decreased whereas, k_Arginine and k_Lysine increased under the low
> moisture conditions** of Coffee-amino acid model systems."

Check `[D]`: high moisture 25 °C → 1.03 (Cys) > 0.15 (Arg) > 0.090 (Lys) ✔, **ratios 11.4 : 1.7 : 1**.
Low moisture 25 °C → 1.03 (Cys) ≈ 1.01 (Arg) > 0.84 (Lys), **ratios 1.23 : 1.20 : 1**.
★ **The 11.4× cysteine advantage in water COLLAPSES TO 1.2× in the low-moisture coffee matrix.**
**⇒ Amino-acid selectivity toward HMF is nearly abolished by the matrix.** That is a first-class
matrix constraint and it is measured, one lab, one method, one pair of experiments.

The mechanistic reading, verbatim: *"Since cysteine bears reactive -SH group, it has **strong avidity
for the double bond of conjugated vinyl compounds** (Friedman, 2003)."*

---

## §8. THE ADDUCTS — HRMS structural evidence `[M]`, Δ < 1 ppm

All confirmed in the **high-moisture** systems at **50 °C for 7 days**, `[M+H]⁺`, positive ESI.
(*"Presence of HMF-amino acid adducts were confirmed in high and low moisture model systems with very
high mass accuracy (Δ < 1 ppm) but only data recorded for high moisture model systems was
illustrated in Fig. 3."*)

| system | adduct | formula | m/z `[M+H]⁺` | Δ (ppm) | confirmed? |
|---|---|---|---|---|---|
| **HMF-Arg** | Michael adduct 1 (1 HMF + 1 Arg) | C₁₂H₂₀N₄O₅ | **301.1506** | 0.0633 | ✔ |
| | **Schiff base** | C₁₂H₁₈N₄O₄ | **283.1401** | 0.0208 | ✔ |
| | Michael adduct 2 (2 HMF + 1 Arg) | C₁₈H₂₆N₄O₈ | 427.1823 | — | ✘ **not confirmed** |
| | Adduct 3 (Arg + Michael adduct 1) | C₁₈H₃₄N₈O₇ | **475.2622** | 0.4608 | ✔ |
| **HMF-Cys** | ★ **Michael adducts 1 and 1′** (via **–NH₂** and via **–SH**; same formula) | C₉H₁₃NO₅S | **248.0587** | 0.1688 | ✔ |
| | Michael adduct 2 (2 HMF + 1 Cys) | C₁₅H₁₉NO₈S | **374.0906** | 0.9233 | ✔ |
| | Adducts 3 / 3′ (Cys + Michael adduct 1/1′) | C₁₂H₂₀N₂O₇S₂ | 369.0785 | — | ✘ **not found** |
| | **Schiff base** | C₉H₁₁NO₄S | **230.0482** | 0.0736 | ✔ |
| **HMF-Lys** | Michael adduct 1 | C₁₂H₂₀N₂O₅ | **273.1445** | 0.2326 | ✔ |
| | Michael adduct 2 | C₁₈H₂₆N₂O₈ | **399.1762** | 0.8053 | ✔ |
| | **Schiff base** | C₁₂H₁₈N₂O₄ | **255.1339** | 0.5101 | ✔ |
| | Adduct (Lys + Michael adduct 1) | C₁₈H₃₄N₄O₇ | 419.2500 | — | ✘ **not observed** |

**★ STOICHIOMETRY, for the sulfur lane:** cysteine forms **both** a **thiol-Michael** adduct and an
**amine-Michael** adduct of identical formula (C₉H₁₃NO₅S, m/z 248.0587), **plus** a Schiff base
(C₉H₁₁NO₄S, m/z 230.0482), **plus** a 2:1 HMF:Cys Michael adduct (C₁₅H₁₉NO₈S, m/z 374.0906).
⇒ **The HMF + cysteine stoichiometry is not 1:1.** A model carrying a single 1:1 HMF+Cys edge
under-counts the HMF consumed per cysteine.

**Mechanism, verbatim §3.1:** *"nucleophilic groups (–SH, –NH₂) of amino acids side chains are easily
added to **β-carbon of HMF through Michael type addition**. Besides, if the addition occurs between
amino group of amino acids and carbonyl group of HMF, reaction proceeds through **elimination of one
molecule of water resulting the formation of Schiff base**."*

---

## §9. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL`

| # | constraint | anchor |
|---|---|---|
| S1 | ★ **HMF self-degradation is negligible (0.9 % / 7 d / 5 °C / pH 3.5) — the amine reaction IS the sink** | §3.1 |
| S2 | ★ **Cysteine is the dominant HMF sink among Arg/Cys/Lys, by 11.4× over Lys in water, because of the thiol** | §3.3, §7.1 |
| S3 | ★ **The 11.4× cysteine selectivity collapses to 1.2× in a low-moisture roasted-coffee matrix** | §7.1 ⇒ **matrix abolishes amine selectivity** |
| S4 | **The reaction runs at 5 °C.** *"The reactions of HMF with cysteine were found to take place even at temperatures as low as 5 °C."* ⇒ **there is no temperature threshold to switch this sink on** | Conclusion |
| S5 | **Both Michael addition and Schiff-base formation occur, confirmed by HRMS for all three amino acids** | §3.2 |
| S6 | **Multi-HMF adducts form** (2 HMF : 1 Cys; 2 HMF : 1 Lys; 2 Arg : 1 HMF) ⇒ **stoichiometry ≠ 1:1** | §3.2 |
| S7 | **Amine blocking suppresses HMF self-dimerisation** — *"presence of amino compound in reaction medium limited the decomposition of HMF. Carbonyl group of HMF was blocked due to Schiff base formation with amino compounds thus dimerization of HMF was inhibited"* `[C]` Nikolov & Yaylayan 2011b | §3.1 ⇒ **the amine sink and the self-degradation sink are NOT additive; they compete for the same carbonyl** |
| S8 | **Low moisture reduces the temperature dependence of the reaction** — *"limited water in reaction system restricted the temperature dependence"* | §3.4 |
| S9 | **The HMF sink is chemically the same reaction class as the acrylamide sink** — literature `[C]`: Ea(acrylamide + cysteine) = **24.4 kJ/mol** (Hidalgo 2011); Ea(acrylamide + butylamine) = **44 kJ/mol** (Zamora 2010) — bracketing this paper's HMF+Cys 29.6 and HMF+Arg 25.5 | §3.4 ★ **an independent external plausibility check on the Ea magnitude, from a different lab and a different Michael acceptor** |
| S10 | `[C]` **Zou et al. 2015**: cysteine + glycine simultaneously reduced acrylamide and HMF; **HMF depleted by 94.9 %** after 0.1 mmol HMF + 1 mmol cysteine at **160 °C for 15 min** — *"cysteine depleted most of the HMF dominantly by reacting with HMF through Michael adduction"* | §3.3 ★ **the only high-temperature datum for this reaction anywhere in the cluster; cited, not measured here** |
| S11 | `[C]` **Zou et al. 2016**: 25 µmol/mL cysteine depleted **91 %** of 315.3 µg/mL HMF at **40 °C** | §3.1 |
| S12 | `[C]` **Gökmen et al. 2012**: *"carbonyl group of HMF rapidly reacted with amino acid asparagine leading to acrylamide formation through Schiff base formation"* | §1 ⇒ **the HMF sink feeds the acrylamide lane** |

---

## §10. VERIFIED NEGATIVES `[NEG]`

- **No temperature above 50 °C.** The entire dataset is 5–50 °C. **Extrapolating these Ea to
  120–200 °C cooking spans a 150–250 °C gap with a 3-point, R² = 0.80–0.91 line whose slope rests on
  a single point the authors declare compromised.** ⇒ see §11.
- **No pH ladder.** High-moisture system fixed at pH 3.5; **the coffee system's pH is never
  measured or reported.**
- **No amino acid other than arginine, cysteine and lysine.** No glycine, no asparagine, no
  histidine, no proline. (Asparagine is the one the acrylamide lane needs.)
- **No adduct quantification.** The adducts are confirmed by HRMS but **never quantified** — no
  yields, no mass balance, no molar recovery.
- **No amino-acid consumption measurement.** Only HMF is followed by HPLC.
- **No second-order constant is published.** The §4 second-order values are `[D]`.
- **No `[B]` for the coffee system** — 8 µmol in 1 g of solid, no volume ⇒ **no second-order
  constant can be recovered for any low-moisture row.**
- **No furfural, no α-dicarbonyls, no reaction network, no multiresponse fit.**
- **No supplementary material referenced or held.**

---

## §11. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| **Table 1, HIGH-MOISTURE k′ (Arg, Cys rows)** | `[F]` | ★ **USE** — `×10⁻¹ day⁻¹`, pseudo-first order at `[AA] = 20 mM`, pH 3.5, 5/25/50 °C |
| **Table 1, HMF-Lys row** | `[F]` | **USE-Q with the SD flag of §4.1** — the printed SDs (45–63 %) are ~10× the paper's uniform 5 % pattern; as printed the row is REFUSE, under the ×10 correction it is usable |
| **Table 1, COFFEE k′ (all six)** | `[F]` | **USE-Q — pseudo-first-order only; no `[B]` volume exists, so no second-order constant can be derived. Coffee pH unknown.** |
| ★ **`k(HMF + Cys) = 3.95 / 5.15 / 23.3 M⁻¹ day⁻¹` at 5 / 25 / 50 °C, pH 3.5** | `[D]` | ★★ **USE — the corpus's FIRST genuine second-order HMF-sink constant, and the sulfur lane's HMF–thiol crosslink.** Label: *"second-order rate constant for HMF consumption by cysteine in 0.05 % benzoic acid at pH 3.5, derived from the published pseudo-first-order k′ and the stated 20 mM amino-acid loading; author-declared apparent at 50 °C."* |
| `k(HMF + Arg) = 0.60 / 0.75 / 2.75` and `k(HMF + Lys) = 0.44 / 0.45 / 0.80 M⁻¹ day⁻¹` | `[D]` | **USE-Q** (Lys with the §4.1 flag) |
| **Ea = 25.457 (Arg) and 29.639 (Cys) kJ/mol, high moisture** | `[F]`, **reproduced `[D]`** | ★ **USE-Q.** Label: *"apparent activation energy for HMF consumption by arginine/cysteine, aqueous pH 3.5, 5–50 °C, 3 points, refit R² = 0.87, slope dominated by the 50 °C point at which the authors declare the pseudo-first-order assumption compromised."* **Never "the HMF-sink barrier."** |
| Ea = 10.018 (HMF-Lys), 11.55 (Coffee-Cys), 12.33 (Coffee-Lys) | `[F]`, reproduced `[D]` | **PRIOR-ONLY** — small signal (1.8–2.1× over 45 °C) |
| **Ea = 3.275 kJ/mol (Coffee-Arg)** | `[F]` | ★ **REFUSE — chemically impossible barrier from a 1.22× k spread over 45 °C. Same defect class as Ma 2022's Ea(MGO) = 1.84 (round-3 §A.3.2). Record the refusal.** |
| **A = 1.63 (HMF-Lys) and 2.55 (Coffee-Arg) day⁻¹** | `[F]` | ★ **REFUSE — sign-flipped intercepts. Correct values: 0.6156 and 0.3926 day⁻¹ `[D]`.** |
| **A = 12.08 (Coffee-Cys) and 9.67 (Coffee-Lys) day⁻¹** | `[F]` | ★ **REFUSE AS LABELLED — the two are swapped. Correct: Coffee-Cys 9.689, Coffee-Lys 12.115 day⁻¹ `[D]`.** |
| A = 612.59 (HMF-Arg) and 23980.59 (HMF-Cys) day⁻¹ | `[F]`, reproduced `[D]` | **USE** — the only two correct as printed |
| **any extrapolation of these Ea above 50 °C** | — | ★ **REFUSE — RECORD AS A PROHIBITED DERIVATION**, in the register of `k3` §C.1 and round-3 §A.3.4. The data span 5–50 °C; cooking is 120–200 °C. |
| **Depletion percentages (§7)**, incl. **HMF self-degradation 0.9 % / 7 d / 5 °C** | `[M]` | ★ **USE — model-free, and the 0.9 % self-degradation figure is the load-bearing one** |
| **Selectivity ratios 11.4 : 1.7 : 1 (water) vs 1.23 : 1.20 : 1 (coffee)** | `[D]` | ★ **USE — the cluster's cleanest same-method matrix-vs-water pair on a rate ratio.** Directly relevant to `k3` §C.2's complaint that no same-method matrix/water pair exists. |
| HRMS adduct table (§8) | `[M]` | ★ **STRUCTURAL — USE.** Δ < 1 ppm, R = 100 000. **Non-1:1 stoichiometry is the key row.** |
| S1–S12 | `[M]`/`[F]`/`[C]` | **STRUCTURAL** — S1, S2, S3, S4 and S6 are the highest-value; S10 and S11 are `[C]` and must stay labelled as such |
