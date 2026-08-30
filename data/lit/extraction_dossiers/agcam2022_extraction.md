# Ağçam 2022 — HMF and furfural formation kinetics in fruit-juice-based model mediums

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### Wave brief: *"content was NEVER read by anyone — approach with full skepticism, verify identity first."* Identity verified. **The round-3 expectation of what this paper contains is substantially wrong.**

---

## §0. IDENTITY — VERIFIED, DOI CORRECT, BUT THE DESCRIPTION IS NOT

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/agcam2022.pdf` (808,521 bytes) | `ls` |
| **DOI** | **`10.1007/s12161-021-02214-x`** ✔ matches the round-3 expectation | printed on p.1 of the PDF |
| title | *"A Kinetic Approach to Explain Hydroxymethylfurfural and Furfural Formations Induced by Maillard, Caramelization, and Ascorbic Acid Degradation Reactions in Fruit Juice-Based Mediums"* | p.1 |
| **author** | ★ **Erdal Ağçam — SINGLE AUTHOR** | p.1 (`Erdal Agcam¹`, *"Extended author information available on the last page"*) |
| journal | ***Food Analytical Methods*** **15** (2022) 1286–1299 | running head, all pages |
| dates | Received 19 Aug 2021 · Accepted 15 Dec 2021 · online 13 Jan 2022 | p.1 |
| version | published version of record, Springer | p.1 |

**⚠️ THIS IS NOT A GÖKMEN-SCHOOL PAPER.** It is by a different, unaffiliated single author in a
different journal (*Food Analytical Methods*, not *Food Chemistry*/*JAFC*). Every other paper in the
K5a cluster comes from Hacettepe FOQUS. This one does not, and its methodology is entirely
different.

---

## §1. ★★ FOUR WAYS THE ROUND-3 EXPECTATION IS WRONG

`research_round3_channels.md` §E TIER 2 listed this as *"fruit-juice HMF/furfural kinetics, three
source reactions"*. Having now read it:

**(1) IT IS NOT A FRUIT JUICE.** "Fruit juice–based medium" means a synthetic aqueous solution:
**glucose (12 % w/v) OR fructose (12 % w/v), ± asparagine 1500 mg/L, OR ascorbic acid 500 mg/L,
titrated to pH 3.5 with 0.1 N citric acid.** No juice, no pulp, no phenolics, no minerals.
Verbatim, §Handling: *"sugar, asparagine, and ascorbic acid contents in fruit juice-based mediums
were fixed as 12 %, 1500 mg/L, and 500 mg/L, respectively, **to mimic fruit juices**."*
The author names the limitation himself in the Conclusion: *"**natural fruit juices are complex
mediums due to including numerous bioactive ingredients that may affect formation of HMF/furfural.
Therefore, more studies are necessary to reveal effect of medium pH, aw, and presence of other
ingredients**"*.

**(2) IT IS NOT A MULTIRESPONSE KINETIC MODEL AND HAS NO REACTION NETWORK.** There is no ODE system,
no Athena Visual Studio, no HPD interval, no model discrimination. It is a set of **independent
uniresponse curve fits** in **Excel 2010** — verbatim, §Data Analyses: *"The collected data from the
experimental studies were evaluated by **Excel 2010 software (Microsoft Corp., USA)** for
determination of the best kinetic models."*

**(3) ★ IT CONTAINS NO 3-DEOXYGLUCOSONE, NO FRUCTOFURANOSYL CATION, AND NO α-DICARBONYL OF ANY KIND
AS A MEASURED SPECIES.** Only HMF, furfural, ascorbic acid and asparagine are quantified (§Sample
Analyzing). **⇒ THIS PAPER CANNOT SPEAK TO THE FRUCTOSE-LIMB-vs-3-DG QUESTION AT ALL.** It discusses
3-DG and FFC only in prose, citing Perez Locas & Yaylayan 2008 and Göncüoğlu Taş & Gökmen 2017. Its
"three source reactions" are **caramelization / Maillard / ascorbic-acid degradation**, i.e. three
*precursor classes*, not three branches of the HMF network.

**(4) The comparison it DOES make is fructose-vs-glucose, at food pH, which is genuinely useful —
just not for the branch question.** See §4.

---

## §2. SYSTEM AND CONDITIONS `[M]`

| item | value | anchor |
|---|---|---|
| mediums (5) | **Glu** = glucose 12 % w/v · **Fru** = fructose 12 % w/v · **Glu+Asn** = glucose 12 % + asparagine 1500 mg/L · **Fru+Asn** = fructose 12 % + asparagine 1500 mg/L · **AA** = ascorbic acid 500 mg/L | §Handling |
| **pH** | ★ **3.5**, set with **0.1 N citric acid** for the sugar mediums and **0.1 N NaOH** for the ascorbic-acid medium | §Handling |
| temperatures | **90, 105 and 120 °C** | §Heat Treatment |
| times | **0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 min** (11 sampling points) | §Heat Treatment |
| vessel | **10 mL sample in a 12 mL glass tube, hermetically sealed with Teflon screw caps.** ⚠️ *"The air in head space of the tubes **was not modified**."* — i.e. **aerobic, with a 2 mL headspace** | §Heat Treatment |
| heating | oil bath (IKA); **come-up time excluded from the treatment time**; cooled under cold water to ~25 °C | §Heat Treatment |
| replication | **triplicate for all mediums** | §Heat Treatment |
| analysis | HPLC-PDA/RID (Shimadzu), C18 ACE 4.6×250 mm 5 µm, isocratic **methanol/acetic acid/water 20/1/79 v/v/v**, 0.5 mL/min, 30 °C; PDA at **250 nm (AA), 275 nm (furfural), 285 nm (HMF)** | §Sample Analyzing |
| standards | authentic HMF (Merck), furfural (Sigma-Aldrich), ascorbic acid; **5-point calibration, R² ≥ 0.998** | §Sample Analyzing |
| LOD / LOQ | **HMF 1.36 / 1.51 µg/L**; **furfural 1.06 / 1.21 µg/L**; **AA 30.8 / 39.5 µg/L** | §Sample Analyzing |
| asparagine | followed by the **dansyl chloride** method (Sherovski 2018) | §Sample Analyzing |

⚠️ **PRINTED INCONSISTENCY IN THE HEATING TIME.** Methods say **0–90 min** and list 11 sampling
points ending at 90 min. The Results say **"0–120 min"** at least four times (e.g. *"thermal
treatment (90–120 °C, 0–120 min)"*, *"exposed to 90, 105, and 120 °C heating for 0–120 min"*).
**The Methods are authoritative: 0–90 min.** The "0–120 min" appears to be the temperature range
mis-copied into the time slot. Also, line 720: *"for heating at 90, 105, and **105 °C**"* — should
read 120 °C. **Three transcription errors in the Results section; treat the Results prose as
lower-confidence than the tables.**

---

## §3. THE KINETIC MODELS — six candidates, transcribed verbatim

The paper fits **two families × three orders** to each accumulation curve.

**Mass-formation-based** (Eq. 1: `d[C]/dt = k[C]ⁿ`):
- Eq. 2, n = 0: `[C_HMF,F] = [C_HMF,F]₀ + k·t`
- **Eq. 3, n = 1: `[C_HMF,F] = [C_HMF,F]₀ · e^{k·t}`** ← **the model actually used for the sugar mediums**
- Eq. 4, n = 2: `1/[C_HMF,F] = 1/[C_HMF,F]₀ − k·t`

**Precursor-based** (integrating over the precursor's own decay):
- Eq. 5, n = 0: `[C] = [C]₀ + k_i·t`
- **Eq. 6, n = 1: `[C] = [C]₀ + [C_i]₀·(1 − e^{−k_i·t})`** ← **the model used for the ascorbic-acid medium**
- Eq. 7, n = 2: `[C] = [C]₀ + [C_i]₀·(1 − 1/(1 + k_i·t·[C_i]₀))`

Maillard step (Eq. 8): `−d[C_Asn]/dt = −d[C_Sug]/dt = k₂·[C_Asn]·[C_Sug]`, simplified to
pseudo-first order (Eqs. 10–13) under `[C_Sug]₀ ≫ [C_Asn]₀`, giving `k′ = k₂·[C_Sug]`.
Combined (Eq. 15): `d[C_HMF,F]/dt = k₁·[C_Sug]ⁿ + k′·[C_Asn]`.

Derived outputs: **Eq. 19** `k = A·e^{−Ea/RT}` · **Eq. 20** `Q₁₀ = (k_{T2}/k_{T1})^{10/(T2−T1)}`
· **Eq. 21** `[C] = [C]₀·2^{t/t₂}` (doubling time).

**Model selection, verbatim:** *"For all the samples exposed to heating at various temperatures,
**the best fitting performances were determined for first-order mass formation–based kinetic
model** (Table 1)."* For furfural from ascorbic acid: *"the kinetic model … **was decided to be the
first order of the precursor-based kinetic model**, since the highest fitting performances were
determined for ascorbic acid degradation at the studied temperatures (the average R²-value was
determined as **0.991** for first-order degradation kinetic of ascorbic acid)."*

### 3.1 ⚠️ THE STRUCTURAL DEFECT THAT GOVERNS EVERY CONSTANT IN THIS PAPER

**Eq. 3 is a pure exponential GROWTH law with no sink term.** But the paper's own headline finding is
that **HMF is being destroyed**:

> "HMF decomposed remarkably in the mediums including asparagine and heated > 90 °C." (Highlights)
> "The possible explanation is that **HMF is quickly decomposing in fruit juice-based medium
> containing amino acid for heating conditions higher than 90 °C.** In other words, HMF
> concentration formed by caramelization and Maillard reactions is **transforming to another
> compound**, if processing temperature is higher than 90 °C in the presence of amino acid."
> "it can be said that HMF is **unstable compound under acidic and high temperature conditions
> especially in the presence of amino acids**, and therefore, it is interacting with other compounds
> in the medium, **which results in high rate of decomposition**."

**⇒ Every `k` in Table 1 for a `+Asn` medium is a NET (formation − destruction) lump, not a
formation constant, and the author says the destruction term is large.** The sugar-only mediums are
cleaner but not clean (self-degradation and dimerisation are documented by the author's own cited
literature). **This must travel with any ingestion.**

### 3.2 ⚠️ THE "PRECURSOR-BASED" R² COLUMNS ARE DEGENERATE AND CARRY NO INFORMATION

In Tables 1 and 2, the two precursor-based columns (n = 1 and n = 2) are **numerically identical to
each other and to the mass-formation n = 0 column** in essentially every row — e.g. glucose 90 °C:
`0.926 | 0.926 | 0.926`; glucose 105 °C: `0.815 | 0.815 | 0.815`; fructose 120 °C:
`0.892 | 0.892 | 0.892`. The paper's footnote explains only the n = 0 case: *"The
zero-order-precursor-based kinetic model simplifies to the zero-order-mass formation-based kinetic
model."* The n = 1 and n = 2 collapse is separately explained in the Conclusion, verbatim:

> "The modeling study showed that the choosing of hexose sugars or amino acid as precursors to
> explain HMF and furfural formation kinetics **were not resulted with high performance outputs
> because of excess sugar concentration in the fruit juice–based mediums.**"

i.e. with `[Sug]₀ = 12 %` and HMF at µmol/L, `(1 − e^{−k t})` never leaves its linear regime, so
Eqs. 6 and 7 both degenerate to the straight line of Eq. 5. **⇒ 60 % of the R² table (24 of 40
cells across both tables) is a repeated value, not four independent model comparisons. The
"first-order mass formation is best" conclusion is a comparison between essentially two distinct
models, not six.**

---

## §4. TABLE 1 — HMF kinetic outputs

**Transcribed verbatim from p.1293.** Caption verbatim: *"The kinetic outputs for HMF formation in
different medium-based fruit juice precursor compounds\*"*. Footnote `*`: *"The kinetic outputs were
calculated according to first-order-mass formation-based model. **± standard deviation**"*.
Provenance **`[F]`**.

| precursor | T (°C) | R² mass n=0 | n=1 | n=2 | R² prec. n=1 | n=2 | **k (1/min ×10⁻⁴)** | **t₂ (min)** | **Ea (kJ/mol)** | **Q₁₀** |
|---|---|---|---|---|---|---|---|---|---|---|
| **Glucose** | 90 | 0.926 | 0.956 | 0.974 | 0.926 | 0.919 | **52.3 ± 0.5** | 116.3 ± 2.7 | **107.7 ± 0.4** | **2.48 ± 0.01** |
| | 105 | 0.815 | 0.982 | 0.896 | 0.815 | 0.815 | **455.2 ± 1.1** | 16.6 ± 0.0 | | |
| | 120 | 0.871 | 0.908 | 0.488 | 0.870 | 0.870 | **781.2 ± 0.7** | 8.9 ± 0.1 | | |
| **Fructose** | 90 | 0.904 | 0.976 | 0.979 | 0.904 | 0.904 | **140.1 ± 1.3** | 46.8 ± 2.2 | **71.0 ± 1.3** | **1.82 ± 0.02** |
| | 105 | 0.824 | 0.978 | 0.724 | 0.824 | 0.824 | **575.1 ± 0.9** | 12.7 ± 0.3 | | |
| | 120 | 0.892 | 0.904 | 0.517 | 0.892 | 0.892 | **833.7 ± 34.6** | 8.7 ± 0.6 | | |
| **Glucose + Asn** | 90 | 0.942 | 0.967 | 0.982 | 0.942 | 0.942 | **61.5 ± 1.4** | 111.8 ± 2.8 | **98.5 ± 1.0** | **2.29 ± 0.02** |
| | 105 | 0.891 | 0.987 | 0.848 | 0.891 | 0.891 | **414.4 ± 5.2** | 18.2 ± 0.4 | | |
| | 120 | 0.893 | 0.913 | 0.506 | 0.893 | 0.893 | **729.5 ± 2.0** | 9.6 ± 0.0 | | |
| **Fructose + Asn** | 90 | 0.927 | 0.995 | 0.966 | 0.927 | 0.927 | **213.0 ± 8.5** | 34.0 ± 1.8 | **49.0 ± 2.6** | **1.51 ± 0.03** |
| | 105 | 0.880 | 0.969 | 0.701 | 0.880 | 0.880 | **547.9 ± 4.7** | 13.3 ± 0.1 | | |
| | 120 | 0.816 | 0.904 | 0.344 | 0.816 | 0.816 | **735.5 ± 18.3** | 8.5 ± 0.0 | | |
| **Ascorbic acid** | — | — | — | — | — | — | ★ **NO HMF FORMED AT ALL** | — | — | — |

## §5. TABLE 2 — furfural kinetic outputs

**Transcribed verbatim from p.1296.** Footnote `*`: *"The kinetic outputs were calculated according
to the precursor-based kinetic model for ascorbic acid, while the first-order-mass formation–based
kinetic model was used in calculation of kinetic outputs of the other precursors. ± standard
deviation"*. Provenance **`[F]`**.

| precursor | T (°C) | R² mass n=0 | n=1 | n=2 | R² prec. n=1 | n=2 | **k (1/min ×10⁻⁴)** | **t₂ (min)** | **Ea (kJ/mol)** | **Q₁₀** |
|---|---|---|---|---|---|---|---|---|---|---|
| **Glucose** | 90 | 0.966 | 0.973 | 0.979 | 0.968 | 0.965 | **28.8 ± 8.1** | 237.4 ± 68.7 | **109.7 ± 9.5** | **2.53 ± 0.20** |
| | 105 | 0.871 | 0.958 | 0.962 | 0.871 | 0.871 | **162.2 ± 0.9** | 39.2 ± 1.3 | | |
| | 120 | 0.913 | 0.976 | 0.784 | 0.913 | 0.913 | **447.2 ± 14.8** | 16.7 ± 0.7 | | |
| **Fructose** | 90 | 0.959 | 0.973 | 0.976 | 0.959 | 0.959 | **52.7 ± 4.6** | 122.0 ± 10.7 | **101.3 ± 2.2** | **2.35 ± 0.04** |
| | 105 | 0.878 | 0.995 | 0.901 | 0.878 | 0.878 | **343.0 ± 7.8** | 21.6 ± 0.3 | | |
| | 120 | 0.894 | 0.929 | 0.571 | 0.894 | 0.894 | **672.6 ± 23.8** | 10.7 ± 0.6 | | |
| **Glucose + Asn** | 90 | 0.965 | 0.970 | 0.972 | 0.965 | 0.965 | **23.8 ± 1.7** | 254.2 ± 41.3 | **120.9 ± 2.2** | **2.77 ± 0.05** |
| | 105 | 0.887 | 0.970 | 0.995 | 0.887 | 0.887 | **143.8 ± 1.5** | 47.4 ± 0.2 | | |
| | 120 | 0.914 | 0.966 | 0.764 | 0.914 | 0.914 | **501.8 ± 8.0** | 15.4 ± 0.1 | | |
| **Fructose + Asn** | 90 | 0.965 | 0.975 | 0.943 | 0.965 | 0.965 | **96.4 ± 2.9** | 65.3 ± 1.0 | **76.7 ± 0.6** | **1.91 ± 0.01** |
| | 105 | 0.875 | 0.989 | 0.878 | 0.875 | 0.875 | **381.0 ± 2.6** | 19.9 ± 0.1 | | |
| | 120 | 0.819 | 0.915 | 0.474 | 0.818 | 0.816 | **665.0 ± 9.7** | 9.8 ± 0.2 | | |
| **Ascorbic acid** | 90 | 0.950 | 0.945 | 0.613 | 0.950 | 0.950 | **0.033 ± 0.001** | **15.2 ± 0.2** ⚠️ | **130.3 ± 0.1** | **3.00 ± 0.00** |
| | 105 | 0.975 | 0.856 | 0.406 | 0.975 | 0.975 | **0.244 ± 0.005** | **11.3 ± 0.2** ⚠️ | | |
| | 120 | 0.982 | 0.792 | 0.272 | 0.982 | 0.982 | **0.886 ± 0.018** | **9.3 ± 0.2** ⚠️ | | |

---

## §6. ★ WAVE K5a CROSS-CHECK — this paper's internal arithmetic, tested `[D]`

Unlike every other paper in the K5a cluster, **this one publishes both the per-temperature k and
the derived Ea/Q₁₀/t₂, so all three derivations are independently testable.**

### 6.1 Ea and Q₁₀ — ★ **ALL EIGHTEEN VALUES REPRODUCE**

OLS `ln k` vs `1/T` over 363.15 / 378.15 / 393.15 K, and `Q₁₀ = (k₁₂₀/k₉₀)^{10/30}`:

| row | **Ea recomputed** | Ea published | **Q₁₀ recomputed** | Q₁₀ published | R² of the 3-pt fit |
|---|---|---|---|---|---|
| HMF Glucose | **107.8** | 107.7 | **2.46** | 2.48 | 0.906 |
| HMF Fructose | **71.1** | 71.0 | **1.81** | 1.82 | 0.911 |
| HMF Glc+Asn | **98.5** | 98.5 | **2.28** | 2.29 | 0.923 |
| HMF Fru+Asn | **49.3** | 49.0 | **1.51** | 1.51 | 0.928 |
| FUR Glucose | **108.8** | 109.7 | **2.49** | 2.53 | 0.984 |
| FUR Fructose | **101.3** | 101.3 | **2.34** | 2.35 | 0.942 |
| FUR Glc+Asn | **120.8** | 120.9 | **2.76** | 2.77 | 0.994 |
| FUR Fru+Asn | **76.8** | 76.7 | **1.90** | 1.91 | 0.954 |
| FUR Ascorbic acid | **130.5** | 130.3 | **2.99** | 3.00 | 0.990 |

**Every Ea reproduces to within 0.9 kJ/mol and every Q₁₀ to within 0.04.** ★ **This is the only
activation-energy table in the entire K5a cluster that is internally reproducible from its own
published rate constants.** (Compare: Gürsul Aktağ 2020 — 0 of 43 reproduce, six are mathematically
underivable; Kocadağlı JAFC — 11 of 17 mismatch by >15 %; Göncüoğlu Taş and Şen and Kocadağlı-flour
publish no Ea at all.)

### 6.2 Is it in the Ma-2022 plateau-artefact regime? — ★ **NO**

Round-3 §A.3.2 refuses Ma's `Ea(5-HMF) = 14.85 kJ/mol` because its k moved only **1.4×** over 40 °C
(a saturated response fitted with an exponential). Here:

| row | k(90) → k(120) | fold change over 30 °C |
|---|---|---|
| HMF Glucose | 52.3 → 781.2 | **14.9×** |
| HMF Fructose | 140.1 → 833.7 | **6.0×** |
| HMF Glc+Asn | 61.5 → 729.5 | **11.9×** |
| HMF Fru+Asn | 213.0 → 735.5 | **3.5×** |
| FUR Ascorbic | 0.033 → 0.886 | **26.8×** |

**Every row moves ≥3.5×, and most move 6–27×. These are genuine temperature responses, not
plateaus.** ⇒ **The Ma-2022 defect does NOT apply to this paper.** Note also that the Ea values sit
in the physically ordinary 49–130 kJ/mol band, bracketing van Boekel's ~120 kJ/mol norm.

⚠️ **The one row to watch is `HMF Fru+Asn`, Ea = 49.0.** It has the smallest fold-change (3.5×), it
is the medium where the author asserts HMF destruction is strongest, and it is the lowest Ea in the
set. **Its low value is most plausibly the net-of-destruction lump of §3.1, not a lower barrier.**
The systematic pattern supports this: **adding asparagine LOWERS Ea in both sugars**
(glucose 107.7 → 98.5; fructose 71.0 → 49.0) — exactly what a growing destruction term does to an
apparent formation Ea.

### 6.3 ★ THE DOUBLING TIMES — eight rows check out, the ascorbic-acid row is broken by 10⁴×

Eq. 21 with first-order growth gives `t₂ = ln2 / k`. Testing every row `[D]`:

| row | published t₂ (90/105/120 °C) | **ln2/k recomputed** | verdict |
|---|---|---|---|
| HMF Glucose | 116.3 / 16.6 / 8.9 | 132.5 / 15.2 / 8.9 | ✔ ~1–14 % |
| HMF Fructose | 46.8 / 12.7 / 8.7 | 49.5 / 12.1 / 8.3 | ✔ ~5 % |
| HMF Glc+Asn | 111.8 / 18.2 / 9.6 | 112.7 / 16.7 / 9.5 | ✔ ~1–9 % |
| HMF Fru+Asn | 34.0 / 13.3 / 8.5 | 32.5 / 12.7 / 9.4 | ✔ ~5–10 % |
| FUR Glucose | 237.4 / 39.2 / 16.7 | 240.7 / 42.7 / 15.5 | ✔ ~1–9 % |
| FUR Fructose | 122.0 / 21.6 / 10.7 | 131.5 / 20.2 / 10.3 | ✔ ~4–8 % |
| FUR Glc+Asn | 254.2 / 47.4 / 15.4 | 291.2 / 48.2 / 13.8 | ✔ ~2–15 % |
| FUR Fru+Asn | 65.3 / 19.9 / 9.8 | 71.9 / 18.2 / 10.4 | ✔ ~6–10 % |
| **FUR Ascorbic acid** | **15.2 / 11.3 / 9.3** | **210 045 / 28 408 / 7 823** | ★ **BROKEN BY ~10⁴×** |

**And it cannot be rescued by any unit reinterpretation.** If the AA row's k is read as plain
`1/min` rather than `1/min × 10⁻⁴`, `ln2/k` = **21.0 / 2.8 / 0.8 min** — still not 15.2 / 11.3 / 9.3,
and now non-monotone against the published values. **⇒ The ascorbic-acid `t₂` triplet is
internally inconsistent with its own `k` under every reading. REFUSE those three numbers.**

Note the AA row uses the *precursor-based* model (Eq. 6), not Eq. 3, so `t₂` as defined by Eq. 21
is arguably inapplicable to it — but the paper applies it anyway and quotes the values in the text
(*"The lowest t₂ value was determined for the ascorbic medium heated at 120 °C"*). **The AA row's
`k` and `Ea` are mutually consistent (§6.1) and survive; only its `t₂` is refused.**

---

## §7. ABSOLUTE CONCENTRATIONS — the paper's most transferable content `[M]`

All ranges are min–max over 0–90 min at that temperature, in **µmol/L**.

### 7.1 HMF

| medium | 90 °C | 105 °C | 120 °C |
|---|---|---|---|
| Glucose | 0.25–0.44 | 0.27–**11.45** | 0.28–**318.33** |
| **Fructose** | 0.46–1.74 | 0.50–**66.65** | 1.54–**1541.27** |
| Glucose + Asn | 0.28–0.48 | 0.28–**8.52** | 0.30–**197.58** |
| Fructose + Asn | 0.43–2.69 | 0.48–**52.10** | 0.74–**1148.66** |
| **Ascorbic acid** | ★ **none** | ★ **none** | ★ **none** |

### 7.2 Furfural

| medium | 90 °C | 105 °C | 120 °C |
|---|---|---|---|
| Glucose | 0.27–0.35 | 0.23–1.15 | 0.28–**11.77** |
| Fructose | 0.05–0.09 | 0.06–1.04 | 0.06–**21.29** |
| Glucose + Asn | 0.26–0.33 | 0.27–0.99 | 0.28–**15.88** |
| Fructose + Asn | 0.05–0.12 | 0.05–1.26 | 0.05–**30.72** |
| **Ascorbic acid** | 0.01–0.84 | 0.02–5.06 | 0.02–**16.87** |

### 7.3 ★ THE RATIOS — the form that survives every caveat `[D]`

**Fructose : glucose HMF yield** (same tube, same method, same day):

| condition | ratio |
|---|---|
| no Asn, 90 °C | 1.74/0.44 = **3.95×** |
| no Asn, 105 °C | 66.65/11.45 = **5.82×** |
| no Asn, 120 °C | 1541.27/318.33 = **4.84×** |
| +Asn, 120 °C | 1148.66/197.58 = **5.81×** |

★ **Fructose out-yields glucose on HMF by 4–6× at pH 3.5, and the ratio is stable across 30 °C and
across the amine on/off condition.** This is a clean, unit-free, single-lab number.

⚠️ **AND IT CONTRADICTS THE LITERATURE VALUE THE AUTHOR QUOTES APPROVINGLY.** Verbatim: *"Similar
results were determined by **Lee and Nagy (1990)** in model systems. They reported **36 times higher
HMF formation for fructose-based medium acidified by citric acid to the pH of 3.5** than glucose."*
**The author's own experiment — same acid, same pH, same comparison — gives 4–6×, i.e. 6–9× lower
than the number he cites as corroborating.** `[M]` vs `[C]` discipline: **use 4–6× (measured here);
refuse 36× (cited, and contradicted by the citing paper's own data).**

**HMF : furfural ratio at 120 °C** `[D]`: glucose **27.0**, fructose **72.4**,
glucose+Asn **12.4**, fructose+Asn **37.4**.
★ **Compare Lee et al. 2022 (round-3 §A.3.3), a solid sponge cake at 200 °C: HMF:furfural = 10:1
(glucose alone) and 6.25:1 (glucose + leucine).** Same ordering (amine lowers the ratio: 27→12.4
here, 10→6.25 there) but the **absolute ratio is 1.2–7× higher in the aqueous acid medium.**
⇒ **The HMF:furfural ratio is matrix- and pH-dependent by up to ~7×, from two independent labs.**

**Amine effect on HMF `[D]`:** +Asn/no-Asn = 1.09 (Glu, 90 °C), 0.74 (105 °C), 0.62 (120 °C);
1.55 (Fru, 90 °C), 0.78 (105 °C), 0.75 (120 °C).
★ **Asparagine INCREASES HMF at 90 °C and DECREASES it at 105 and 120 °C — a sign crossing in
temperature, in one pot.** The author's explanation is HMF consumption by the amine (§3.1).
**⇒ A directly relevant constraint for the HMF-sink question: the amine sink turns on somewhere
between 90 and 105 °C at pH 3.5.** ★ This dovetails with **Hamzalıoğlu & Gökmen 2018**, which
measures the HMF + amino-acid reaction at 5–50 °C and finds it already running at 5 °C.

---

## §8. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL`

| # | constraint | anchor |
|---|---|---|
| S1 | ★ **ASCORBIC ACID PRODUCES NO HMF AT ALL** — measured at three temperatures over 90 min, LOD 1.36 µg/L | Results, §HMF: *"no HMF formation was determined after AA (500 mg/L) degradation under 90–120 °C heating and 0–90 min holding time conditions"*; corroborated `[C]` by Yuan & Chen 1998 |
| S2 | ★ **ASCORBIC ACID IS A MAJOR FURFURAL SOURCE** — up to 16.87 µmol/L at 120 °C, the highest furfural of any medium except Fru+Asn | Table 2, Fig. 4 ⇒ **AA is a furfural-only precursor. A model that routes ascorbate to HMF is falsified.** |
| S3 | **Fructose > glucose for HMF by 4–6× at pH 3.5, stable across T and across amine on/off** | §7.3 |
| S4 | ★ **Asparagine's effect on HMF crosses sign between 90 and 105 °C** | §7.3 — a real sign crossing for `k3` §B.5 |
| S5 | **Furfural IS formed from hexoses alone**, at pH 3.5, with no pentose and no cation catalyst — 0.28–11.77 µmol/L (glucose, 120 °C) | Table 2 ⚠️ **This contradicts Gökmen & Şenyuva 2006, which the author cites**: *"Interestingly, in the absence of cations, they reported **no furfural formation** in the medium containing asparagine as amino acid and glucose as reducing sugar."* ⇒ **a direct measured-vs-cited conflict on furfural's hexose origin; record both.** |
| S6 | **For furfural, the amine effect is also temperature-crossing**: *"formation rates are lower in the medium including glucose + asparagine than the medium containing just glucose. However, interestingly, its concentration was higher at 120 °C"* | Results, §Furfural |
| S7 | **First-order mass formation beats zero- and second-order for HMF and furfural from sugars**; precursor-based models are degenerate under excess sugar | Tables 1–2 + Conclusion, §3.2 |
| S8 | **Q₁₀ ordering: ascorbic-acid furfural (3.00) > glucose (2.48–2.53) > glucose+Asn (2.29/2.77) > fructose (1.82/2.35) > fructose+Asn (1.51/1.91)** | Tables 1–2 ⇒ **the more reactive the precursor, the flatter the temperature response** — the generic signature of an approach to a diffusion/saturation limit |
| S9 | Mechanistic prose `[C]`, not measured here: *"the critical precursor compound for HMF formation from glucose is 3-deoxyglucosone, while it is fructofuranosyl cation from fructose … **the fructose moiety and glucose in the medium contributed to 90 % and 10 % in HMF formation** at high temperature"* (Perez Locas & Yaylayan 2008) | Results, §HMF — ⚠️ **`[C]` ONLY. This paper measured no 3-DG and no FFC. Do not attribute the 90/10 split to this DOI.** |

---

## §9. VERIFIED NEGATIVES `[NEG]`

- **No 3-deoxyglucosone, no fructofuranosyl cation, no glucosone, no 1-DG, no glyoxal, no
  methylglyoxal, no α-dicarbonyl of any kind measured.**
- **No reaction network, no ODE system, no multiresponse fit, no HPD interval, no model
  discrimination.** Excel curve fits only.
- **No sucrose** — despite the paper's own argument that sucrose is the FFC source in juices.
- **No real fruit juice.** Synthetic mediums only.
- **No browning / colour / melanoidin.**
- **No acrylamide**, despite asparagine being present (the paper discusses it only from citations).
- **No water activity, no a_w variation.** pH fixed at 3.5; **no pH ladder.**
- **Supplementary files 1–4** (medium composition, dansyl method, standards, reaction-pathway
  scheme) are **not in the PDF on disk**. Supplementary file 4 is the paper's only reaction scheme.

---

## §10. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| **Table 1 k (12 HMF values), Table 2 k (15 furfural values)** | `[F]` | **USE-Q** — `1/min ×10⁻⁴`, first-order apparent **accumulation** constants at pH 3.5 in a 12 % hexose solution. ⚠️ **The four `+Asn` HMF rows are net-of-destruction lumps (§3.1).** |
| **Ea: HMF 107.7 / 71.0 / 98.5 / 49.0; furfural 109.7 / 101.3 / 120.9 / 76.7 / 130.3 kJ/mol** | `[F]`, **reproduced `[D]`** | ★ **USE-Q — the only internally reproducible Ea set in the K5a cluster.** Label discipline, mandatory: *"apparent lumped accumulation Ea for HMF (or furfural) in a 12 % w/v aqueous hexose solution at pH 3.5 (citric acid), 90–120 °C, 3 points, triplicate, fitted as first-order growth with no sink term."* **Never "the HMF barrier."** |
| **Ea(HMF, Fru+Asn) = 49.0** | `[F]` | **PRIOR-ONLY** — lowest fold-change (3.5×), strongest destruction term; most likely a net lump |
| **Q₁₀: 2.48 / 1.82 / 2.29 / 1.51 (HMF) and 2.53 / 2.35 / 2.77 / 1.91 / 3.00 (furfural)** | `[F]`, reproduced `[D]` | ★ **USE — dimensionless, reproduced exactly, and immune to the k-unit question** |
| **t₂ for the 8 sugar rows** | `[F]` | **USE-Q** — reproduce to 1–15 % |
| **t₂ for the ascorbic-acid row (15.2 / 11.3 / 9.3 min)** | `[F]` | ★ **REFUSE — inconsistent with its own k by ~10⁴× under every unit reading** |
| **precursor-based R² columns (n = 1, n = 2) in both tables** | `[F]` | **REFUSE — degenerate; identical to the n = 0 mass-formation column (§3.2)** |
| absolute HMF and furfural ranges (§7.1–7.2) | `[M]` | **USE — authentic standards, R² ≥ 0.998, LOD 1.06–1.36 µg/L** |
| **fructose:glucose HMF ratio 4–6×** | `[D]` | ★ **USE** |
| **"36 times higher" (Lee & Nagy 1990)** | `[C]` | **REFUSE — contradicted 6–9× by the citing paper's own measurement under identical conditions** |
| **HMF:furfural = 27.0 / 72.4 / 12.4 / 37.4 at 120 °C** | `[D]` | **USE-Q** — pairs against Lee 2022's 10:1 / 6.25:1 in cake |
| S1 (ascorbate → no HMF) and S2 (ascorbate → furfural) | `[M]` | ★ **STRUCTURAL — USE. The cleanest ascorbate-route constraint the corpus has.** |
| S4, S5, S6 (sign crossings; furfural from hexose) | `[M]` | **STRUCTURAL — USE**, S5 with the cited-conflict flag |
| S9 (the 90 %/10 % fructose/glucose split) | `[C]` | **DO NOT ATTRIBUTE TO THIS DOI** — cited from Perez Locas & Yaylayan 2008; nothing in this paper measures it |
| Results-section prose | — | **lower confidence** — three transcription errors ("0–120 min" ×4; "105 °C" for 120 °C). Prefer the tables and the Methods. |
