# K6a — THE SULFUR TEMPERATURE-LADDER CLUSTER

### Consolidated synthesis for wave K6a. 2026-08-29 / 30.
### **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or `FIT_HOLDOUT_DECLARATION.md` was written, staged or modified. Nothing was committed.**

**Scope.** The nine PDFs named in the K6a brief, read end to end from their text layers and from
600–900 dpi rasters of every figure that carries data.

**Per-paper dossiers** (all new, all in this directory):
`zhai2023jafc_extraction.md` · `zhai2023foodchem_extraction.md` · `wang2026_extraction.md` ·
`zhang2026_extraction.md` · `chan1994_extraction.md` · `feng2022_extraction.md` ·
`meng2017_extraction.md` · `ames2001_extraction.md` · `gigl2021_extraction.md`

**Provenance codes** are those of `k3` §"HOW TO READ THIS FILE": `[M]` measured · `[C]` cited ·
`[F]` fitted by the authors · `[D]` derived by this wave and never printed · `[NEG]` verified
negative · `[!]` defect. **Verdicts:** `USE` · `USE-Q` · `RATIO-ONLY` · `STRUCTURAL` ·
`PRIOR-ONLY` · `REFUSE`.

---

# §1. ★ WRONG-FILE IDENTITIES AND NAMING TRAPS — REPORT THESE FIRST

| # | trap | correction |
|---|---|---|
| **1** | ★★ **THE TWO ZHAI FILES ARE SWAPPED RELATIVE TO THE BRIEF.** | **`Zhai2023b.pdf`** (3.9 MB) = **JAFC `10.1021/acs.jafc.3c04166`** — the TTCA + **cysteine** paper the brief called "Zhai2023, THE PRIORITY". **`Zhai2023.pdf`** (1.6 MB) = ***Food Chem.* `10.1016/j.foodchem.2022.134420`** — the TTCA + **xylose** paper the brief called "Zhai2023b". Any ingestion keyed on the filename stem has them the wrong way round. Dossiers here are named by **paper identity** (`zhai2023jafc` / `zhai2023foodchem`), the same convention K5a used for the two Kocadağlı files. |
| **2** | ★★ **THE BRIEF'S "SEVEN TEMPERATURES 90–150 °C" IS A BROWNING LADDER, NOT A VOLATILE LADDER.** | Zhai JAFC runs seven temperatures **only for A₄₂₀** (Fig. 1b). **Volatiles are measured at 100, 120 and 140 °C only.** There is no seven-rung MFT/FFT ladder in that paper, or anywhere in the corpus. What the paper *does* have is better: a **39 compound × 8 condition × 3 temperature numeric heat map — 936 values** (§6 of that dossier), never before ingested. |
| **3** | ⚠️ **`Zhai2023.pdf` is a 2022 paper carrying a 2023 volume number.** | Online **28 Sep 2022**, assigned to *Food Chem.* **404 (2023)**. Cite as Zhai et al. 2023; **date the data to Sept 2022** when reasoning about precedence — which matters a great deal (§2). |
| **4** | ⚠️ **`meng2017.pdf` self-identifies as 2016 throughout.** | Published online **3 Oct 2016**, DOI `10.1080/09168451.2016.1238295`, **no volume, no issue, no page numbers** — it is the online-first version. Keep `meng2017` as the repo key; cite as "published online 3 October 2016"; **verify the volume via Crossref before freezing a bibliography.** |
| **5** | ⚠️ **`meng2017` has THREE temperatures, not two.** | The brief expected 95/120 °C. The paper runs **80, 95 and 120 °C × two hold times (5 and 20 min)** — a genuine 3-rung ladder, and its curvature is a live result (§4.4). |
| **6** | ✔ **`Wang2026.pdf`, `Zhang2026.pdf`, `feng2022.pdf`, `ames2001.pdf`, `gigl2021.pdf`, `chan1994.pdf` are all exactly what the brief expected**, DOIs raster-verified on the page. |
| **7** | ⚠️ **`chan1994.pdf` is an ABBYY FineReader OCR scan.** | Digits were re-verified against a 600 dpi raster cell by cell; **mangling is confined to Greek letters and spaces, never digits.** The OCR layer is safe for numbers **after** that check, not before. |
| **8** | ★ **μ → m raster check: the hazard did NOT fire in any of the nine papers.** | Verified on 600–900 dpi rasters for the three Arbortext/ACS files (`Zhai2023b`, `feng2022`, `gigl2021`), the two Elsevier files, `ames2001` (Xyvision/Distiller) and `meng2017`. **`Kang2026.pdf` remains the only corrupted file in the corpus.** One unrelated typo found: `1.8 m M` on meng2017 p. 2, self-corrected three lines above. |

---

# §2. ★★★ THE HEADLINE: THE CORPUS'S SULFUR LADDER IS ONE EXPERIMENT, PUBLISHED FOUR TIMES

The repo's only sulfur temperature ladder — Kang 2026's Table S4, on which the `kang_140C_*` and
`kang_switch_on_*` hold-out rows and the 55.1 kJ mol⁻¹ cysteine-depletion override all rest — is a
**re-publication of a 2022 data set**.

| # | publication | date | where the ladder appears | what it adds |
|---|---|---|---|---|
| **1** | ★ **Zhai et al., *Food Chem.* 404:134420** | **Sep 2022** | **Table 1, TTCA columns**, 34 compounds × 3 T, **with SDs and Duncan letters** | **the primary source** |
| 2 | Zhai et al., *JAFC* 71:14300 | Sep 2023 | **Fig. 3, column E** | +5 compounds, +3 hold times, +Cys arm; SDs dropped |
| 3 | Kang et al., *Sust. Food Technol.* 4:3239, **Table S4** | Mar 2026 | the SI table the repo ingested | corrects one cell |
| 4 | Kang, **Table S5, pH-7 column** | Mar 2026 | the same 100 °C column again | (already flagged in `kang2026_SI_extraction.md` §4c) |

**Verification.** Over the 102 shared cells (34 rows × 3 temperatures), sources 1, 2 and 3 agree
**cell for cell to the last printed decimal in 101 of 102**; all six class subtotals and all three
grand totals agree exactly, **including an arithmetic error at 120 °C that is reproduced
identically in sources 1 and 3**; the **printed SDs are identical**; and the FFT/MFT ratios
(3.02 / 2.96 / 1.94) are identical. Kang's corresponding authors are **Yuying Fu and Yun Zhai** —
the same Yun Zhai who is first author of both Zhai papers.

**The single discrepant cell is diagnosed.** 2-Methylthiophene at 120 °C reads `1.128 ± 0.126ᵇ` in
Zhai Table 1 for **both** the TTCA and the TTCA+Xyl column — identical mean, SD *and* Duncan letter,
i.e. a copy-paste — while the printed subtotal excludes it. **Kang's `0.000` is the correction; the
two Zhai papers carry the artefact. The true value is `not detected`.**

### 2.1 ★★ THE SEVEN CONSEQUENCES — all binding

1. **The ladder is n = 1.** `research_round4_nulls.md` §C.2's null ("no second sulfur temperature
   ladder for MFT/FFT") is not merely unfalsified — **the apparent existence of a second source was
   an artefact of duplication.** *(But see §3: wave K6a has now found two genuinely independent
   ladders.)*
2. ★ **Re-point the provenance** of every `kang_*` sulfur row and of the 55.1 kJ mol⁻¹ cysteine
   override to `zhai2023foodchem` (Table 1 / p. 6), recording Kang as a re-publication.
3. ★ **Any weighting that counts Kang S4 and S5 as two observations is wrong by 2×**; adding either
   Zhai paper unwittingly makes it 3× or 4×.
4. ★ **A fit/hold-out collision is now possible.** If a later wave ingests Zhai JAFC Fig. 3 **column
   E** as a new anchor while the `kang_140C_*` rows remain hold-outs, **the model is fitted on its
   own hold-out.** Column E must be marked ALREADY-SEEN. **The other 819 of the 936 values in that
   grid are genuinely new.**
5. ★★ **The Tier A / Tier B uncertainty split for the ladder is wrong.** `kang2026_SI_extraction.md`
   §4d assigned Tier A (±15 %) on the strength of Kang's Table S3 calibration curves. **The primary
   source states the opposite basis, in a printed equation: `Wi = f′ × (Ai·ms)/(As·V)` with
   "the correction factor was 1".** Identical numbers cannot come from two quantification routes.
   **⇒ Treat the whole ladder as Tier B semi-quant: × ÷ 2 to × ÷ 3 on absolute magnitude.
   Eₐ and shape are unaffected** (a constant response factor cancels in a slope).
6. ★ **The "measured pH 4.9" comparator** used by `kinetic_core_b2_2_diagnosis.md` §3 originates in
   Zhai Food Chem p. 5, **attached to no temperature and no time**. Fine as a diagnostic; never a
   scored row.
7. ★ **The K5 SD objection should be restated.** Kang's SDs are **not** fabricated — they are Zhai's
   original triplicate SDs, carried through. They remain implausibly tight for HS-SPME, and Table S5
   disagrees with S4 because S5's column is a re-typeset duplicate, not a re-measurement.

---

# §3. ★★ THE QUANTIFICATION-BASIS AUDIT — the wave's second structural finding

Every dossier was required to name the internal standard and state whether quantification is
absolute or single-IS semi-quant, **and to say explicitly what a semi-quant ladder does and does not
license.** The answers vary more than expected:

| paper | internal standard | basis | Eₐ-usable? | absolute yields? |
|---|---|---|---|---|
| **Zhai Food Chem 2022** | 1,2-DCB, 0.054 μg | ★ **SEMI-QUANT, `f′ = 1`, equation printed** | ✅ yes (shape) | ❌ no |
| **Zhai JAFC 2023** | 1,2-DCB, 0.054 μg | ⚠ claims **calibration curves (Table S2)** — **contradicts the 2022 paper for the same numbers** | ✅ | ❌ (§2.1 pt 5) |
| **Kang 2026** (ingested) | 1,2-DCB | ⚠ Table S3 curves exist but **were not used for S4** | ✅ | ❌ |
| **Wang 2026** | 1,2-DCB, 0.1 μg | ⚠ **NEVER STATED**; and 0.5 g of lyophilisate is reported in μg/L with no volume | ✅ | ❌ **unreconstructible** |
| **Zhang 2026** | 1,2-DCB (volatiles) | ★ **non-volatiles ABSOLUTE** (14 MRM curves, R² ≥ 0.993); **volatiles semi-quant** (ng/g) | ✅ | ✅ for non-volatiles only |
| **Feng 2022** | 1,2-DCB, 0.018 μg/μL | ★ **ABSOLUTE** — authentic external standards, per-compound curves (in the missing Table S1); **1-DO and 3-DO semi-quant** | ✅ | ✅ (unauditable) |
| **Meng 2016/17** | — | ★ **SPLIT**: MFT absolute (y = 0.5002x, r = 0.999, 0–20 μg/L); 2FM/ET2MP/BM **no printed curve** | ✅ | ✅ for MFT only |
| **Ames 2001** | 1,2-DCB, 65 ng, ★ **spiked onto the trap AFTER trapping** | ★ **SEMI-QUANT, stated verbatim**: *"Semiquantitative data were obtained from the mass spectral integration report, with reference to the internal standard"* | ✅ | ❌ — the ng/10 g unit is an area ratio in disguise |
| **Chan 1994** | 4-methylthiazole, 500 ppm, GC/AED | ★ **SEMI-QUANT**; ⚠ **the ppm basis is never defined** (aqueous vs 0.5–1.0 mL extract = a 50–100× ambiguity); ⚠ **DMDS has 2 S atoms vs 1** so even in-paper cross-compound k are on different scales | ✅ | ❌ |
| **Gigl 2021** | ERETIC 2 / PULCON, L-tyrosine 5.21 mmol/L | ★ **ABSOLUTE qHNMR** — true bound fractions | ✅ | ✅ ⚠ **n = 1, no SDs, no LOD** |

### 3.1 ★ THE RULE, STATED ONCE FOR THE WHOLE WAVE
> **A single-IS semi-quant ladder IS activation-energy-usable as within-study shape, because a
> constant unknown response factor cancels in a ratio and therefore in an Arrhenius slope.
> It does NOT license absolute yields, cross-compound magnitude comparison, cross-paper magnitude
> comparison, or class subtotals read as physical totals.**

**⇒ Six of the nine papers are semi-quant in whole or in part, and all six are still usable for the
thing this wave was sent to get.** The constraint bites only on magnitudes — which is exactly where
the corpus was already weakest.

---

# §4. ★★★ DOES THE KANG SWITCH-ON REPLICATE? — **NO. It fails on four independent axes, three of them inside its own source laboratory.**

**The claim under test** (`kang2026_SI_extraction.md` §7a): MFT and FFT are flat from 100→120 °C
(×1.12, ×1.10) then **jump** to 140 °C (×4.26, ×2.78); apparent Eₐ climbs from ~6–7 to
97.8 / 69.2 kJ mol⁻¹, **against** a sulfur class that decelerates (57.5 → 35.2). The dossier argued
a saturation artefact could not produce this, because depletion *depresses* top-leg Eₐ.

## 4.1 ★ AXIS 1 — HOLD TIME. **The feature exists at 120 min and at no other hold time in its own experiment.**
From the Zhai JAFC 936-value grid, TTCA arm, same pot, same three temperatures:

| hold | **MFT** ×(100→120) / ×(120→140) | **FFT** ×(100→120) / ×(120→140) |
|---|---|---|
| 40 min | **5.54× / 2.61×** ← *decelerating* | 1.59× / **4.50×** |
| 80 min | 1.84× / 3.25× | 2.50× / 4.12× |
| **120 min** ← *Kang's column* | ★ **1.12× / 4.26×** | ★ **1.10× / 2.79×** |
| 180 min | 2.69× / **1.62×** ← *decelerating* | 1.37× / **1.15×** ← *flat* |

**The mechanism is visible in the raw time courses:** at 140 °C **FFT peaks at 80 min (12.439 μg/L)
and falls to 4.906 by 180 min**, while at 120 °C it is still climbing. **A fixed-time slice across a
family of curves whose maxima move left as temperature rises measures where each curve sits relative
to its own peak, not an activation energy.** The K5 dossier's falsification argument is correct for
a monotone-rising system and **fails when the top-temperature curve has a maximum inside the
window** — which is what the grid shows.

## 4.2 ★ AXIS 2 — CO-SUBSTRATE. **Add an equimolar sugar and it vanishes.**
Zhai Food Chem Table 1, same lab, same three temperatures, same 120-min hold:

| | TTCA arm ×(100→120) / ×(120→140) | **TTCA + Xyl arm** |
|---|---|---|
| **MFT** | 1.12× / **4.26×** | **1.44× / 1.30×** |
| **FFT** | 1.10× / **2.79×** | **1.35× / 0.96×** ← *FFT falls* |

## 4.3 ★★ AXIS 3 — AN INDEPENDENT LABORATORY, FIVE RUNGS. **The thiols peak in mid-ladder and fall.**
Wang 2026 (Ningxia + SJTU; Cys-Amadori + Glu-Amadori in **0.2 M PBS**), 85/95/105/115/125 °C:

| species | ×(85→95) | ×(95→105) | ×(105→115) | ×(115→125) |
|---|---:|---:|---:|---:|
| **MFT** | — | ≈2.5× | **2.6×** | ★ **0.30× ↓** |
| **FFT** | 2.7× | ★ **0.56× ↓** | ★ **0.15× ↓** | ★ **→ 0** |
| Furfural (FFT's precursor) | 2.44× | 2.03× | **1.10×** | **1.03×** ← *saturating* |

**⇒ FFT collapses while its own precursor saturates. That is a SINK observation, in a buffered pot,
in a lab with no connection to the Zhai/Kang chain.**

## 4.4 ⚠️ AXIS 4 — meng2017, three rungs, and the curvature runs the **other way**
MFT and FFT in soy sauce model heating at **80 / 95 / 120 °C**: apparent Eₐ **falls** from 122.8 to
25.0 kJ mol⁻¹ (MFT, 5 min) — **concave, the opposite sign to Kang's convex 6.9 → 97.8** — with an
identical `− + −` residual pattern in all four series. ⚠️ **But the 95 °C rung is time-saturated
(5 → 20 min gains only ×1.03) and the vessel changes at exactly the leg boundary (open glass
cylinder → autoclave).** Curvature magnitude co-varies with the confound.
**Registered as UNRESOLVED, not as a refutation.**

## 4.5 ★ WHAT DOES REPLICATE — three things, and they are worth more than the switch-on
1. ★★ **The CLASS CONTRAST is real and replicates in a completely different medium.**
   In **Ames 2001's extrusion ladder** (xylose/pH 7.5, the one clean ladder of six), the switch-on
   question asked properly — as **excess Eₐ over the bulk**, which cancels response factor,
   residence time, shear, pH drift and common die loss — gives:

   | leg | aliphatic S | **FFT** | S-substituted furans | **MFT** | thiophenes | thiazoles |
   |---|---:|---:|---:|---:|---:|---:|
   | **120→150 °C** | **+108.5** | **+46.7** | **+36.7** | **+33.5** | −6.5 | **−42.7** |
   | 150→180 °C | +5 | ±5 | ±5 | ±5 | ±5 | ±5 |

   **The direct H₂S + carbonyl products (mercaptoketones, furanthiols) switch ON; the ring-closure
   products (thiazoles, thiophenes) do not — the same split Kang found, reached by a completely
   different route.** ★ **And the die-volatilisation bias runs the favourable way** (loss rises with
   die temperature and the thiols are among the more volatile species, so the bias pushes their
   excess Eₐ *more negative*): **the positive low-leg excess is a lower bound.**
2. ★ **Kang's flat LOW leg is substantively vindicated once precursor depletion is divided out.**
   Feng 2022 (GSH-ribose ARP, 100 → 120 °C) measures **93.5 % ARP consumed at 120 °C vs ~48 % at
   100 °C (×1.95 ± 0.16)**. Normalising the thiol fold-changes by that extent gives
   **MFT ×0.84 and FFT ×1.06**, bracketing Kang's ×1.12 / ×1.10. **The furan class does *not*
   flatten under the same correction (×3.28), so the normalisation is discriminating, not a
   rescaling artefact.** Meng's 5-min 95→120 leg independently gives **MFT ×1.68 / FFT ×1.94**,
   within 4–6 % of Feng's raw ARP-alone leg. **Three systems triangulate the low leg at ×1.1 – ×2.1.**
3. ★★ **A DIFFERENT, BETTER-MEASURED THRESHOLD DOES EXIST — and it is at 100→120 °C, not 120→140.**
   Zhai Food Chem's **¹³C₅-xylose isotope table** measures the fraction of each product's carbon
   skeleton coming from *exogenous* sugar:

   | compound | 100 °C | 120 °C | 140 °C |
   |---|---:|---:|---:|
   | 3-Thiophenethiol | **0 %** | 19 % | 34 % |
   | Thieno[3,2-b]thiophene | **0 %** | 21 % | 45 % |
   | 2-Methylthieno[2,3-b]thiophene | **<1 %** | 21 % | 45 % |

   **An isotope ratio is immune to response factor, hold time and headspace partitioning — every
   artefact that contaminates the yield ladders.** ⚠️ MFT and FFT were **not** traced, which is
   precisely the measurement that would have settled the question.

## 4.6 ★★ THE VERDICT, AND WHAT TO DO WITH IT
> **The MFT/FFT "switch-on between 120 and 140 °C" is a HOLD-TIME ARTEFACT of a single 120-minute
> slice through curves that peak and turn over. It does not survive changing the hold time, adding
> a co-substrate, or moving to another laboratory. `kang_switch_on_*` should be RETIRED as a scored
> row, not merely failed — fitting a two-regime or threshold Eₐ form to it, as
> `kang2026_SI_extraction.md` §7a recommended, would fit a sampling artefact into the physics.**
>
> **What survives, and is better supported than the original claim: (i) free thiols behave
> differently from ring-closure products, replicated in water and in an extruder; (ii) something
> turns on between ~100 and ~150 °C on the direct H₂S-addition channel, located by isotope tracing
> at 100→120 °C and by extrusion excess-Eₐ at 120→150 °C; (iii) above ~120–140 °C thiol consumption
> outruns thiol formation.**

---

# §5. ★★ THE Eₐ STRUCTURE THE EVIDENCE NOW SUPPORTS — formation vs decay

## 5.1 FORMATION channels

| channel | Eₐ (kJ mol⁻¹) | source | conditions | verdict |
|---|---:|---|---|---|
| ★ **Cys + sugar → Cys-Amadori** | **100.9** (R² **1.000**, legs 98.9/103.0) | Zhang k15 `[F/D]` | 130–150 °C, pH 6.5 buffered | ★ `USE-Q` — **the best-conditioned constant of the wave** |
| ★ **Cys-Amadori → α-dicarbonyl** | **85.7** (R² **1.000**, legs 86.7/84.7) | Zhang k16 | " | ★ `USE-Q` |
| **Amine + sugar → Amadori (GSH, Gly)** | **107.0**, **73.9** | Zhang k1, k13 | " | `USE-Q` |
| ⚠ **Amadori → α-DC (generic)** | published **53–68**; ★ **low-leg 92–106** | Zhang k4/k8/k11/k14 | " | ★ `PRIOR-ONLY` **at the low-leg value only** — the published values are plateau artefacts (§6.2) |
| **Free cysteine degradation** | **54.7** (legs 62.3 / 46.2, R² 0.996) | Zhai Food Chem p. 6 `[D]` | 100–140 °C, unbuffered | ★ `USE-Q` — **this is the true origin of the repo's 55.1** |
| **GSH → free Cys release** | **66.4** | Feng 2022 `[D]` | 100–120 °C | `USE-Q` |
| **Peptide hydrolysis (GSH, Cys-Gly)** | **121.1**, **138.6** | Zhang k3, k9 | 130–150 °C | ★ `STRUCTURAL` — ⇒ **in a peptide system H₂S supply is hydrolysis-limited, not Cys-degradation-limited** |
| **MFT / FFT in water, whole ladder** | ⚠ **14.5 – 74.8**, hold-time dependent | Zhai JAFC `[D]` | 100–140 °C | ⚠ **no single value exists** (§4.1) |
| **MFT / FFT, 20-min hold** | **57.8**, **72.9** (R² 0.983, 0.984) | Meng `[D]` | 80–120 °C | `USE-Q` — the cleanest aqueous thiol fits in the corpus |
| **MFT / FFT, low leg only** | **30.3**, **44.4** | Feng `[D]` | 100–120 °C | `USE-Q` |
| ★ **MFT / FFT in a LOW-MOISTURE MATRIX** | ★ **103.5**, **107.7** (R² 0.941, 0.966) | Ames `[D]` | 120–180 °C die, xylose pH 7.5 | ★ `PRIOR-ONLY` — see the warning below |
| **aliphatic sulfur, matrix** | **145.6** (R² **0.999**) | Ames `[D]` | " | `PRIOR-ONLY` |
| **Methional, DMDS (Strecker route)** | **74.0 – 103.4** whole-ladder; ⚠ per-leg span **6–15×** | Chan `[D]` | 75–115 °C, pH 6/7/8, **0.1 M phosphate** | `PRIOR-ONLY` |
| **Browning (A₄₂₀)** | ★ **28.3** (R² 0.969) and ★ **30.6** (R² 0.971) | Zhai JAFC + Wang `[D]` | 90–150 °C unbuffered; 85–125 °C buffered | ★★ `USE-Q` — **the best-replicated Eₐ of the wave: two labs, two systems, two buffer states, 8 % apart** |

★★ **THE MOST CONSEQUENTIAL SINGLE COMPARISON IN §5.1:** **the same two thiols carry an apparent Eₐ
of ~30–75 kJ mol⁻¹ in water and ~104–108 kJ mol⁻¹ in a starch extrudate.** Even allowing that the
extrusion value is the weakest class in the wave, the gap is a factor of ~2 and it is in the
direction of a **matrix-slowed** reaction. **⇒ An aqueous-fitted thiol formation barrier should not
be transported into the matrix lane without an explicit transfer term.**

⚠️ **And Ames' pH dependence has the WRONG SIGN relative to aqueous work.** MFT/FFT go from
**literally zero at feed pH 5.5 to 3038 / 5346 (semi-quant units) at pH 7.5**, while the aqueous
literature (and Wang 2026 §5.2) has them *falling* with rising pH. **An aqueous-fitted pH term
applied in a low-moisture matrix would have the wrong sign.** That may be the single most
consequential transfer warning this wave produced.

## 5.2 ★★ DECAY channels — **and the 248 kJ mol⁻¹ thiol sink is refuted**

`kinetic_core_b2_2_diagnosis.md` F-5: the fit put `Ea_decay_thiol_sink` at **248.0 kJ mol⁻¹**,
pressed against its 250 ceiling, against a formation barrier of 76.2 — buying "a sink whose absolute
rate is tiny across the whole 100–145 °C window". **Wave K6a produced four independent measurements
of thiol decay, and none of them is compatible with that.**

| evidence | number | source | what it says |
|---|---|---|---|
| ★★ **Thiol → melanoidin COVALENT binding** | **Eₐ = 60.2 kJ mol⁻¹** (`R ln 67.2 / (1/279 − 1/333)`); 3-point fit 58.5, R² 0.886; ★ **defensible range 7–102** | Gigl `[D]` from k = 1.17e−5 / 2.57e−4 / 7.83e−4 min⁻¹ at 279 / 300 / 333 K | ★ **248 kJ mol⁻¹ predicts k(333)/k(279) = 3.4 × 10⁷; the measured ratio is 67.2 — wrong by 5.1 × 10⁵** |
| ★★ **π-stacking reservoir** | **ΔH° = −19.5 kJ mol⁻¹** (R² 0.979) — ★ **NEGATIVE** | Gigl `[D]` | at 60 min **free FFT is 24.1 % at 279 K but 55.0 % at 333 K — more thiol survives HOT than COLD.** ★ **No single sink with any positive Eₐ can produce that at any prefactor: it is a model-form falsification, not a fit-quality problem** |
| ★★ **Net FFT turnover at 140 °C** | **k_net = 0.0141 min⁻¹, t½ ≈ 49 min** (11.439 → 4.906 μg/L over 120→180 min) — a **lower bound**, since formation continues | Zhai JAFC `[D]` | the sink runs **fast** at 140 °C, in a pot where furfural is still rising (11.039 → 12.386) |
| ★ **FFT collapse over 105→125 °C** | ×0.15 then → 0, while furfural saturates | Wang `[D]` | an independent lab, buffered, sign-consistent |
| ★ **Thiol → disulfide oxidation** | **Eₐ = 122.2 kJ mol⁻¹** (R² 0.971) | Zhang k17 `[F/D]` | ★ **the corpus's only measured barrier for the dimerisation channel** |
| ⚠ **Lumped multi-order sink, as a cautionary example** | **k18 → melanoidin, Eₐ = 118.0**, from a ★ **FIFTH-ORDER rate law** | Zhang `[F]` | **4× the directly measured browning Eₐ (28.3, 30.6).** A lumped sink with unmodelled co-reactants inflates its own fitted "barrier" — ★ **exactly the failure mode that produces a rail-pinned 248** |

### 5.2.1 ★★ THE MECHANISM GIGL SUPPLIES, AND WHY IT MATTERS
Gigl's FFT is the **only** one of 20 odorants in **Scenario IV — covalent *and* π–π**. Two controls
dissect it: **N-acetylcysteine (thiol, no ring) → covalent only; methylated FFT (ring, no free −SH)
→ π–π only.** The instantaneous drop is a **ring** property; the slow decline is a **free −SH**
property. And at 300 and 333 K the free pool goes to **zero**, while at 279 K the π-bound reservoir
is largest and the covalent drain is shut.
**⇒ The π-stacked pool is a reversible reservoir that FEEDS the covalent sink, not a parallel
terminal sink.**

> ★★ **THE STRUCTURAL CONCLUSION: the repo is asking one lumped sink to represent two channels of
> opposite temperature sign (−19.5 and +60 kJ mol⁻¹), one of them reversible and equilibrating in
> under five minutes. A rail-pinned 248-against-250 is the diagnostic signature of that lumping,
> not a physical barrier.**

### 5.2.2 The honest transfer caveat
Gigl's measurements are at **6–60 °C on coffee melanoidins**, three to five decades below the
100–180 °C window. **Any use there is an extrapolation and must be labelled one.** Calibrated
judgement carried over from that dossier: **P(true Eₐ of the covalent channel in 30–100 kJ mol⁻¹)
= 0.85; P(≥ 248) < 0.01; P(the repo's 248 is a lumping/rail artefact) = 0.90; P(a single
positive-Eₐ sink cannot reproduce these data even qualitatively) = 0.95; P(this paper alone suffices
to re-fit a repo parameter) < 0.05.** ★ **It is a refutational and structural anchor, not a
calibration source.** Note also that **60.2 kJ mol⁻¹ lands within 10 % of the 54.7 kJ mol⁻¹
cysteine-degradation barrier already in `MEASURED_EA_OVERRIDES`.**

## 5.3 ★ THE RESULTING PICTURE — formation vs decay, side by side

| | barrier (kJ mol⁻¹) | evidence |
|---|---|---|
| **thiol FORMATION, aqueous** | **~30–75** | Zhai, Feng, Meng |
| **thiol DECAY, covalent capture** | **~60** (range 7–102) | Gigl |
| **thiol DECAY, oxidative dimerisation** | **~122** | Zhang k17 |
| **precursor supply (Amadori → α-DC)** | **~86–106** | Zhang |
| **precursor supply saturates above** | **~140 °C** | Zhang legs; Wang furfural plateau |

★★ **This is a coherent, and quite different, picture from the one the B2.2 fit landed on. Decay
barriers are of the SAME ORDER as formation barriers — decay is somewhat higher, but not by
170 kJ mol⁻¹. That is precisely the configuration that produces a peak-and-decline: at low
temperature formation wins, at high temperature decay catches up, and above ~140 °C precursor
supply saturates while consumption keeps accelerating.** ⇒ The corpus can now explain the turnover
it measures, without a 250 kJ mol⁻¹ rail.

---

# §6. ★ CROSS-PAPER AUDIT RESULTS — what the refits found

## 6.1 Reproducibility scoreboard
| paper | published Eₐ audited | reproduce from the authors' own data? |
|---|---|---|
| **Zhang 2026** | **19 of 19** | ★ **ALL. Mean \|dev\| 2.1 %, median 1.7 %, worst +6.4 % (k2, rounding-limited).** No arithmetic error |
| **Chan 1994** | **6 of 8** (2-acetylthiophene publishes no k, so its Eₐ is unauditable and was refused) | ✔ within **4–12 %**; ★ two refits reproduce the authors' `r²` to the third decimal, which also **proves Figs. 4–5 plot the 2nd replicate** (never stated) and identifies Table I's unlabelled `r²` column as the **Arrhenius** one |
| every other paper | — | **publishes no Eₐ at all** |

## 6.2 ★★ THE PLATEAU-ARTEFACT CLASS (Ma-2022, K5a §6.3) — the cleanest example yet
**Zhang's four Amadori→α-DC steps saturate on the 140→150 °C leg:**

| step | leg 130→140 | leg 140→150 | published Eₐ | 3-pt R² |
|---|---:|---:|---:|---:|
| k8 | 99.2 | ⚠ **4.7** | 53.5 | **0.794** |
| k14 | 96.9 | ⚠ **7.3** | 53.4 | **0.811** |
| k4 | 105.6 | **22.0** | 68.1 | **0.879** |
| k11 | 92.3 | 35.5 | 66.25 | 0.940 |

★ **The authors noticed and wrote it down** — *"k values increased rapidly from 130 to 140 °C but
slowly from 140 to 150 °C, including k4, k8, k11, k13, k14, and k16"* — **and published a single
Arrhenius Eₐ for each step anyway.** The same table contains steps that do *not* show it (k15, k16:
legs within 4 %), so this is a property of the Amadori-degradation family, not of the fit.
**⇒ Ingest the low-leg values (92–106), never the published 53–68.**

## 6.3 ★ A DEFECT CLASS THIS WAVE ADDS — **the arbitrary-scale node**
**8 of Zhang's 19 constants (k4, k8, k11, k13, k14, k16, k18, k19) touch a node measured as a UV
absorbance or a GC total, not as a concentration.** Their **Eₐ survive** (a constant scale cancels in
`ln k₂/k₁`) but their **magnitudes do not**, and neither does any comparison between a k on a
measured node and a k on an absorbance node.
★ **Consequence: the paper's central magnitude claim — that the largest rate constants are the
Amadori→α-DC degradations — is a scale artefact plus a plateau artefact and should not be carried
forward as chemistry.** This is the same distinction K5a §2.3 drew for the fructose limb, now
generalised: **in a multiresponse fit, the Eₐ ordering is transferable and the k ordering is not.**

## 6.4 Other defects worth recording
| paper | defect |
|---|---|
| **Zhai JAFC** | the text's `63.035 μg/L` (140 °C, TTCA, 180 min) is a typo for **65.035**, which the figure and Fig. 4c's bar both give. 7 of 8 other text-quoted totals reproduce exactly |
| **Zhai Food Chem** | 2-methylthiophene copy-paste (§2); Conclusion's pyrazine total `10.020` vs Table 1's `10.296`; one inverted Duncan letter |
| **Zhang 2026** | ★ **the α-DC mass balance is inconsistent with equations 14–15** (first order in one place, fifth and third in another); k20–k22 declared and never reported; **no units on any k**; **no confidence intervals anywhere**, so identifiability is untestable; Cys-Amadori's calibration intercept crosses zero at 0.475 mmol/L, above its lowest standard |
| **Chan 1994** | ★ **the Conclusions claim methional's Eₐ is "17–19 kcal/mol"; its own Table I says 19.39–27.07** (both raster-verified) — the paper's most quotable sentence is wrong. Methods give the 95 °C series as 20/40/60/100/120 min; all three 95 °C figures label 20/40/60/80/100. Between-replicate spread reaches **21.8 %** — *that*, not r², is the right error bar. And a mass balance the paper never does shows **the three measured products account for < 0.5 % of the methionine**, so the unexplained plateaus are **not** bulk substrate exhaustion |
| **Ames 2001** | ★★ **three of the six ladders have their ENTIRE 150 °C column collapse ~10× across every chemical class at once** (glucose pH 6.5; xylose pH 6.5 and 5.5) — one bad run, unfalsifiable because there is **one extrusion per condition** (the "triplicate" is three isolates from one extrudate). **Those three ladders are REFUSED; only xylose pH 7.5 is clean.** Table 5 passes all 99 subtotal checks; **Table 4 has 5 genuine errors**, one a diagnosable column shift. **A published erratum found:** the 66 % maximum xylose sulfur RA is at pH 7.5, not pH 6.5 |
| **Feng 2022** | one arithmetic slip in the text ("1.98×" should be 1.612 from its own Table 2) |
| **Gigl 2021** | four author-side integrity flags: 3-methylbutanal at 96 h has **four inconsistent values**; FFT 32.3 vs 32.2 %; the authentic-brew pyrazine "47.8 → 24.4 mmol/L" is **chemically impossible** (10× the protocol spike ≈ 7.9 g/L) though its 51 % ratio is exact → **RATIO-ONLY**; and the text's "72 h to zero at 60 °C" is contradicted by Fig. 4A, which shows zero at **48 h** on a grid with no 72-h point |

---

# §7. ★ CONSOLIDATED PARAMETER TABLE — the wave's ingestible output

`[M]` measured · `[F]` fitted by authors · `[D]` derived here. **All Eₐ in kJ mol⁻¹.**

## 7.1 `USE` / `USE-Q` — ingestible with the stated qualification

| # | parameter | value | conditions | prov | source | qualification |
|---:|---|---|---|---|---|---|
| 1 | ★ **Eₐ, Cys + glucose → Cys-Amadori** | **100.9** (R² 1.000) | 130–150 °C, pH 6.5, 0.2 M phosphate | F/D | Zhang k15 | no CI published; extrapolation below 130 °C |
| 2 | ★ **Eₐ, Cys-Amadori → α-dicarbonyl** | **85.7** (R² 1.000) | " | F/D | Zhang k16 | product node is an absorbance — Eₐ survives, k does not |
| 3 | **Eₐ, amine + glucose → Amadori (GSH)** | **107.0** | " | F/D | Zhang k1 | |
| 4 | ★ **Eₐ, thiol → disulfide (Cys-Gly dimer)** | **122.2** (R² 0.971) | " | F/D | Zhang k17 | ★ the corpus's only such barrier |
| 5 | **Eₐ, GSH hydrolysis** | **121.1** (R² 0.991) | " | F/D | Zhang k3 | |
| 6 | ★ **Eₐ, free-Cys degradation** | **54.7** (R² 0.996; legs 62.3/46.2) | 100–140 °C, 10 mM, pH 7, unbuffered, 120 min | D | Zhai FoodChem p. 6 | ★ **not independent of the repo's existing 55.1 — it is the same experiment** |
| 7 | **Eₐ, GSH → Cys release** | **66.4** | 100–120 °C | D | Feng | 2-point |
| 8 | ★ **Eₐ, browning (A₄₂₀)** | ★ **28.3** (R² 0.969) and **30.6** (R² 0.971) | 90–150 °C unbuffered / 85–125 °C buffered | D | Zhai JAFC + Wang | ★ **two labs, 8 % apart — the wave's best replication** |
| 9 | **Eₐ, MFT / FFT at a 20-min hold** | **57.8 / 72.9** (R² 0.983 / 0.984) | 80–120 °C | D | Meng | ⚠ hold-time dependent (see #10) |
| 10 | ⚠ **Eₐ, MFT, hold-time sensitivity** | **25.0 → 47.7** on the SAME leg, purely by changing the hold from 5 to 20 min | 95→120 °C | D | Meng | ★★ **any repo parameter named "Eₐ of MFT formation" is undefined without a hold time** |
| 11 | **Eₐ, MFT / FFT low leg** | **30.3 / 44.4** | 100–120 °C, GSH-ribose ARP | D | Feng | |
| 12 | ★ **k_net, FFT decay at 140 °C** | ★ **0.0141 min⁻¹** (t½ 49 min) | 120→180 min, TTCA/Cys pot | D | Zhai JAFC | ★ **lower bound** — net of ongoing formation; headspace basis; 3 post-peak points are visibly NOT single-exponential |
| 13 | **k_net, MFT decay at 140 °C** | 0.0065 min⁻¹ | " | D | Zhai JAFC | " |
| 14 | ★★ **Eₐ, thiol → melanoidin covalent capture** | ★ **60.2** (range **7–102**) | **279–333 K**, coffee HMW melanoidins, qHNMR | D | Gigl | ★ **a 3–5 decade extrapolation to the cooking window — label it one** |
| 15 | ★★ **ΔH°, thiol–melanoidin π-stacking** | ★ **−19.5** (R² 0.979) — **negative** | " | D | Gigl | confirmed independently on a pyrazine (−13.4) |
| 16 | ★ **furan + H₂S → thiol molar yield** | ★ **22 %** (furfural→FFT), **35 %** (4-OH-5-Me-furanone→MFT), 29 % (furan→thiophene), 15 % (2-Me-furan→2-Me-thiophene) | (NH₄)₂S 10 mM, **120 °C, 20 min** | D | Zhai JAFC Fig. 6 | ★ **the corpus's first O→S branch fraction**; neat 10 mM model, 3 decades above pot levels; factor-2 interval |
| 17 | **furan-pool consumption rate** | **> 78 % in 10 min; complete by 20 min** | " | M | Zhai JAFC | pseudo-first-order k ≳ 0.15 min⁻¹ |
| 18 | ★ **α-dicarbonyl hierarchy at peak** | **3-DX 0.81–1.13 ≫ 1-DX 0.11–0.18 > MGO 0.072–0.082 ≫ GO 0.006–0.010 ≫ DA 0.003–0.004 mmol/L** | 100–140 °C, TTCA ± Xyl | D | Zhai ×2 | ★ **replicated across both Zhai papers and both sugar arms — a ~200× span** |
| 19 | **3-DX peak time** | **60 min (+Xyl) / 80–100 min (TTCA)** | 100–140 °C | D | Zhai ×2 | |
| 20 | **exogenous-carbon incorporation** | **0 / 19–21 / 34–45 %** at 100/120/140 °C | ¹³C₅-Xyl, 120 min | M | Zhai FoodChem Table 2 | ★ response-factor-free |
| 21 | **Cys-Amadori vs Glu-Amadori thermal optimum** | **105 °C vs 115 °C** | buffered, peak-area basis | M | Wang Fig. 5B | ratio-only scale |
| 22 | **Ames extrudate pH offset** | ★ **1.4–2.6 units BELOW feed pH**, measured | extrusion, 120–180 °C die | M | Ames Table 1 | ★ `USE` — a measured matrix pH shift the repo has nowhere else |
| 23 | **Ames measured vs target die temperature** | ★ **−7 to +14 °C**; the xylose pH 7.5 upper leg is **153→174 °C (21 K)**, not 30 K | " | M | Ames Table 1 | ★ **the nominal ladder is not the real one; Eₐ on the true axis are ~40 % lower on that leg** |
| 24 | **Meng in-matrix MFT sensory JND** | **0.2 μg/L** — ★ **50× the cited 4 ng/L wine-model threshold** | soy sauce | M | Meng | ★ a matrix-vs-water threshold pair |

## 7.2 `RATIO-ONLY`
Every Wang 2026 volatile value; every Ames ng/10 g; every Chan ppm; Zhang's volatile ng/g;
Gigl's authentic-brew pyrazine concentrations; the FFT:MFT ratio at any single (T, t) cell.

## 7.3 ★ `STRUCTURAL` — shape constraints, no number to fit
1. ★★ **Free thiols in a Cys/pentose pot PEAK and then DECLINE; the peak moves earlier as
   temperature rises** (≈ 80 min at 140 °C for FFT). *Zhai JAFC, Wang, Meng.*
2. ★★ **Above ~120–140 °C precursor supply saturates while consumption keeps accelerating.**
   *Zhang's leg structure; Wang's furfural plateau against FFT collapse.*
3. ★★ **The π-stacked odorant pool is a reversible reservoir feeding the covalent sink, not a
   parallel terminal sink.** *Gigl.*
4. ★ **Direct H₂S-addition products (mercaptoketones, furanthiols) and ring-closure products
   (thiazoles, thiophenes) have different temperature responses** — replicated in water and in an
   extruder. *Kang/Zhai + Ames.*
5. ★ **Excess H₂S abolishes the furan pool**: furfural = 0 in the +Cys arm at 120 and 140 °C.
   *Zhai JAFC.*
6. ★ **N-heterocycles are exactly zero in a TTCA-only pot at 100/120/140 °C and appear as soon as
   free xylose is present** — pyrazine formation is **dicarbonyl-limited, not amine-limited**.
   *Zhai Food Chem.* Corroborated by **Wang: pyrazine not detected at pH 5 or 6 at any temperature.**
7. ★ **Added cysteine DECREASES 2-acetylthiazole** (16× at 120 °C) while increasing every other
   S-heterocycle. *Zhai JAFC.* Corroborated by Wang's 115 °C dip and the authors' own
   "suppressed by high temperature" annotation.
8. ★ **In a peptide system, H₂S supply is hydrolysis-limited** (121–139 kJ mol⁻¹) rather than
   Cys-degradation-limited (55). *Zhang.*
9. **Scheme discrimination:** both the GSH-hydrolysis route and the GSH-Amadori route are needed;
   neither alone fits. *Zhang, 17 responses, rRMSE 20–30× lower.*
10. ⚠ **Fructose is undetectable in a buffered pH-6.5 glucose pot at 130–150 °C.** *Zhang.*
    ★ **This conflicts with K5a's structural fructose-limb constraint** (established at 160–200 °C in
    melts and low-moisture matrices). **Open contradiction — not adjudicated here.**

## 7.4 ★ `REFUSE` — with the reason recorded so a later wave does not re-ingest

| item | reason |
|---|---|
| ★ **The MFT/FFT "switch-on" Eₐ (6.9→97.8; 5.8→69.2) as a physical two-regime barrier** | hold-time artefact; fails on four axes (§4) |
| ★ **Kang Table S4 / S5 as independent observations** | duplicates of Zhai Food Chem Table 1 (§2) |
| ★ **The Tier A ±15 % uncertainty on the Kang/Zhai ladder** | the primary source used `f′ = 1` (§3) |
| **`63.035 μg/L`** (Zhai JAFC text) | arithmetic error; the figure gives 65.035 |
| **`1.128` for 2-methylthiophene, TTCA, 120 °C** | copy-paste artefact; true value `not detected` |
| **Zhang's 19 k MAGNITUDES** | ★ **no units printed**, mixed reaction orders, 8 of 19 on arbitrary-scale nodes |
| **Zhang k18 (118.0) as a barrier** | fifth-order rate law; 4× the measured browning Eₐ |
| **Zhang k2 (138.9), k5 (184.3)** | rounding-limited to ±40 and ±18 kJ mol⁻¹ |
| **Zhang's published 53–68 for Amadori→α-DC** | plateau artefact; use the low-leg 92–106 |
| **Chan's 2-acetylthiophene Eₐ (30.8, 32.7 kcal)** | no k published anywhere — unauditable |
| **Chan's pre-exponentials** | `ln A` = `ln(f·A)` on an undefined ppm basis |
| **Ames' three collapsed ladders** (glucose pH 6.5; xylose pH 6.5 and 5.5) | whole 150 °C column drops ~10× across every class at once — one bad run, unfalsifiable at n = 1 extrusion |
| **Ames' glucose-vs-xylose magnitude comparisons** | no batch / run-order / QC statement |
| **Any absolute μg/L from Wang 2026** | basis undescribed; 0.5 g solid reported per litre |
| **Gigl's authentic-brew pyrazine concentrations** | chemically impossible (10× the spike) |
| **Any Eₐ from Wang's MFT/FFT ladders** | ★ **non-monotone in T — no Arrhenius form exists.** The refusal *is* the finding |
| **Every +Cys Eₐ in Zhai JAFC that spans a zero cell** | non-monotone; unfittable |

---

# §8. ★ CROSS-CHECKS AND CONTRADICTIONS — reported, not adjudicated

| # | contradiction | the two sides |
|---|---|---|
| 1 | ★★ **Quantification basis of the corpus's ladder** | Zhai 2022: `f′ = 1` semi-quant, equation printed · Zhai 2023 JAFC: "calibration curves (Table S2)" — **for identical numbers** |
| 2 | ★★ **FFT:MFT branching vs pH** | Kang: **invariant** (2.73–3.02 over pH 5.5–8) · Wang: *"at pH < 7 … mainly 2-furanethiol; at pH > 7 … 2-methyl-3-furanthiol"* |
| 3 | ★★ **pH sign in a matrix** | aqueous literature + Wang: thiols **fall** with rising pH · Ames extrusion: MFT/FFT go from **zero at pH 5.5 to maximum at pH 7.5** |
| 4 | ★ **Fructose limb** | K5a: a hard structural constraint, three matrices, 160–200 °C · Zhang: fructose **undetectable** at pH 6.5, 130–150 °C, and omitted with R² 0.90–0.96 |
| 5 | ★ **Direction of Arrhenius curvature for MFT/FFT** | Kang/Zhai 120-min slice: **convex** (Eₐ rises) · Meng: **concave** (Eₐ falls 122.8 → 25.0) · Zhai at other holds: **both signs** |
| 6 | **Browning Eₐ, direct vs lumped** | direct A₄₂₀ ladders: **28.3 and 30.6** · Zhang's fitted melanoidin lump k18: **118.0** — resolved by the fifth-order rate law (§6.3) |
| 7 | ⚠ **Nonanal in a lipid-free system** | Wang reports nonanal rising **2.6 → 9.4 μg/L** with temperature in a pot containing **no lipid**. Either an artefact (septum/fibre bleed) or a non-lipid route. ★ **Flag to the lipid lane** |

---

# §9. ★★ DRAFT FIT / HOLD-OUT ROLES — **FOR THE ORCHESTRATOR. NO DECLARATION FILE WAS OPENED OR EDITED.**

## 9.1 Corrections to EXISTING rows (do these first)
| action | row(s) | reason |
|---|---|---|
| ★ **RETIRE** | `kang_switch_on_*` (both) | the target is not a property of the chemistry (§4.6) |
| ★ **RE-LABEL provenance** | `kang_140C_MFT`, `kang_140C_FFT`, all `kang_*` sulfur rows, and the 55.1 kJ mol⁻¹ cysteine override | the source is `zhai2023foodchem`; Kang is a re-publication (§2) |
| ★ **RE-BAND uncertainty** | the whole Kang/Zhai ladder | Tier B, × ÷ 2–3 on magnitude, not Tier A ±15 % (§3) |
| ★ **MARK ALREADY-SEEN** | Zhai JAFC Fig. 3 **column E**, all three temperatures | it *is* Kang Table S4 (§2.1 pt 4) |
| **CORRECT** | 2-methylthiophene, TTCA, 120 °C → `not detected` | §2 |

## 9.2 Recommended **FIT** candidates (all genuinely unpublished or newly derived)
| candidate | value | why |
|---|---|---|
| ★ **Eₐ(Cys + sugar → Cys-Amadori)** | **100.9**, R² 1.000 | best-conditioned constant of the wave; measured node both sides |
| ★ **Eₐ(Cys-Amadori → α-DC)** | **85.7**, R² 1.000 | flat legs; scale-free |
| ★ **Eₐ(thiol → disulfide)** | **122.2** | the dimer channel currently has no measured barrier |
| **furan → thiol molar yield** | **15–35 %**, four independent reactions at 120 °C | the O→S branch fraction is currently unconstrained |
| **MFT/FFT time courses at 100 and 120 °C, TTCA arm** (8 points) | Zhai JAFC §6 | constrains formation on the two rungs whose pH state is least wrong; contains no 140 °C information |
| **3-DX peak height and peak time** | 0.81–1.13 mmol/L, 60–100 min | trunk dicarbonyl |

## 9.3 Recommended **HOLD-OUT** — ranked by how much they would teach
| rank | candidate | what it tests |
|---|---|---|
| **1** | ★★ **FFT at 140 °C over four hold times: 5.079 / 12.439 / 11.439 / 4.906 μg/L** | ★ **the turnover the model has never been asked to reproduce.** A formation-only core predicts monotone growth and must fail |
| **2** | ★★ **k_net(FFT, 140 °C) ≥ 0.01 min⁻¹**, as a one-sided test | directly contradicts the B2.2 fit's near-inert 248 kJ mol⁻¹ sink |
| **3** | ★★ **The ¹³C exogenous-carbon fractions (0 / 19–21 / 34–45 %)** | response-factor-free; tests **carbon routing**, not yield — the cleanest hold-out in the wave |
| **4** | ★ **Gigl's sign test: free FFT is HIGHER at 333 K than at 279 K at 60 min** | ★ **no single positive-Eₐ sink can pass it** — a model-form test, not a magnitude test |
| **5** | ★ **Wang: FFT declines 105 → 125 °C while furfural saturates** (sign test) | an independent lab, buffered |
| **6** | ★ **Wang: FFT peaks at 80–100 min and returns to its 60-min level by 120 min** | ★ **the Yiltirak failure in miniature, in a clean two-substrate pot** |
| **7** | **The +Xyl arm of Zhai Food Chem Table 1** (34 compounds × 3 T) | a co-substrate perturbation never seen; sulfur class must go 20.19 → 91.81 → 163.43 |
| **8** | **N-heterocycles = 0 in TTCA-only; ≠ 0 with Xyl; = 0 at pH 5–6 (Wang)** | a hard on/off structural gate |
| **9** | **Ames: MFT/FFT excess Eₐ over bulk is positive on 120→150 and ~zero on 150→180** | the matrix lane's only ladder; a **lower bound** because die volatilisation biases it downward |
| **10** | **Zhang: Amadori→α-DC steps saturate on 140→150 while formation steps accelerate** | the supply/consumption crossover |

## 9.4 ★ Why the Yiltirak failure now has a mechanism
`cutover_final_exam.md` reports the Yiltirak family 0/8, **over**-predicting at 4/4 rungs, worst at
long holds, and names the time axis. **Wave K6a supplies the missing measurement three times over**
— thiols peak and decline, with the peak moving from ≈ 80 min at 140 °C to ≥ 180 min at 120 °C
(Zhai), from 80–100 min in a buffered ARP pot (Wang), and with precursor supply saturating above
140 °C while consumption accelerates (Zhang). **A core with a near-inert sink cannot produce a peak
at all, so it must over-predict every long hold — which is what the exam measures.** The fix
indicated by the evidence is **not** a bigger sink barrier but a **smaller** one (≈ 60–120 rather
than 248) **plus** a co-reactant-bearing rate law, because §6.3 shows a lumped multi-order sink
inflates its own fitted Eₐ.

---

# §10. ★ WHICH SI FILES ARE ACTUALLY NEEDED — ranked, with what each would settle

| rank | SI | what it contains | what it settles |
|---|---|---|---|
| **1** | ★★ **Wang 2026, `10.1016/j.foodres.2026.118535` — Table S2** | the **31-compound numeric volatile data across five pH levels** | ★ **the only independent MFT/FFT ladder in the corpus is currently a digitisation.** Would also give the pH-dependent FFT↔MFT branching that contradicts Kang |
| **2** | ★★ **Zhang 2026, `10.1016/j.foodchem.2026.148681` — Tables S6–S8** | volatiles + **GC-O FD values at 130 / 140 / 150 °C** | ★ **a fourth independent MFT ladder, in a peptide system, at temperatures above every other** |
| **3** | ★ **Zhai JAFC, `10.1021/acs.jafc.3c04166` — Table S2** | **the calibration curves for the 39 volatiles** | ★ **settles the §3 quantification contradiction outright**, and gives the Tier A/B split for the 936-value grid |
| **4** | **Zhai Food Chem, `10.1016/j.foodchem.2022.134420` — Table S2 + Fig. S2 + Table S1** | MGO−GO−Cys model flavour table; **Cys-alone degradation CURVES**; ROAV/thresholds | the dicarbonyl→pyrazine link; shape (not just endpoint) for the 54.7 kJ mol⁻¹ cysteine barrier |
| **5** | **Feng 2022, `10.1021/acs.jafc.2c02949` — Table S1** | the calibration curves, LODs and ranges | makes Feng's "absolute" claim auditable |
| **6** | **Gigl 2021, `10.1021/acs.jafc.1c06163` — Tables S1–S2, Figs S1–S3** | NMR assignment + the underlying time courses | would let the 7–102 kJ mol⁻¹ range be narrowed |
| — | **Ames 2001** | ⚠ **no SI**; screw configuration, die geometry, product moisture and expansion are all deferred to **Bredie et al. 1997/1998** — ref 19, *JAFC* 1998, **46**, 1479–1487, extruded maize flour. ★ **A K6-follow-up candidate in its own right** |
| — | **Chan 1994, Meng 2016** | no SI exists | — |

---

# §11. DECLARED GAPS FROM WAVE K6a

| # | gap | note |
|---|---|---|
| G1 | ★ **MFT and FFT were never isotope-traced.** The one measurement that would settle the switch-on directly | would require new experiments |
| G2 | ★ **No measured H₂S anywhere in nine papers.** Every mechanism in this cluster routes through H₂S and none of them measures it | (NH₄)₂S surrogate models (Zhai JAFC Fig. 6) are the closest available |
| G3 | ★ **The corpus still has no thiol-sink measurement in the 100–180 °C window.** Gigl is 6–60 °C; Zhai's k_net is a lower bound net of formation | the highest-value new experiment the repo could request |
| G4 | **Zhang's rate-constant units are unresolved**, so 19 magnitudes are unusable | ask the authors, or reconstruct the ODE from digitised Fig. 3 |
| G5 | **No confidence intervals on any of Zhang's 19 constants** — identifiability untestable | the Han-2025 defect class cannot even be checked |
| G6 | ⚠ **Wang never states the fixed levels of the untested factors** in each one-factor series, so the ladder's other coordinates are unknown | not recoverable from the paper |
| G7 | ⚠ **The fructose contradiction with K5a** (§8 row 4) needs adjudication at orchestrator level | |
| G8 | ⚠ **Nonanal in a lipid-free system** (§8 row 7) — flag to the lipid lane | |
| G9 | **Whether Kang 2026 cites Zhai as the source of Table S4** — the independence bookkeeping depends on it | read Kang's reference list |
| G10 | **Zhai's TTCA/Cys/Xyl depletion data are deferred to Zhai et al. 2021**, which is not on disk and which also holds the LODs and HPLC methods both papers omit | order it |
