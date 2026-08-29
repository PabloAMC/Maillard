# K5a — THE HMF / CARAMELIZATION KINETICS CLUSTER

### Consolidated synthesis for wave K5a. 2026-08-29.
### **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was written, staged or modified. Nothing was committed.**

**Scope.** The seven PDFs named in the K5a brief, read end to end from their text layers:
the two Kocadağlı & Gökmen 2016 papers, Göncüoğlu Taş & Gökmen 2017, Gürsul Aktağ & Gökmen 2020,
Şen & Gökmen 2022, Ağçam 2022, and Hamzalıoğlu & Gökmen 2018.

**Per-paper dossiers** (all new, all in this directory):
`kocadagli2016jafc_extraction.md` · `kocadagli2016foodchem_extraction.md` ·
`goncuoglutas2017_extraction.md` · `gursulaktag2020_extraction.md` · `sen2022_extraction.md` ·
`agcam2022_extraction.md` · `hamzalioglu2018_extraction.md`

**Provenance codes** are those of `k3` §"HOW TO READ THIS FILE": `[M]` measured · `[C]` cited ·
`[F]` fitted by the authors · `[D]` derived by this analyst and never printed by the authors.
**Verdicts:** `USE` · `USE-Q` · `RATIO-ONLY` · `STRUCTURAL` · `PRIOR-ONLY` · `REFUSE`.

---

# §1. ★ WRONG-FILE IDENTITIES AND NAMING TRAPS — REPORT THESE FIRST

| # | trap | correction |
|---|---|---|
| **1** | ★ **THE TWO KOCADAĞLI FILES ARE SWAPPED RELATIVE TO THE NAÏVE GUESS.** | **`Kocada2016.pdf`** (shorter stem, 810 kB) = **JAFC `10.1021/acs.jafc.6b01862`**, glucose ± NaCl caramelization — **the HMF pick**. **`Kocadagli2016.pdf`** (longer stem, 613 kB) = **Food Chem. `10.1016/j.foodchem.2016.05.150`**, glucose/wheat flour. Any ingestion keyed on the filename stem has them the wrong way round. |
| **2** | ★★ **`Kocadagli2016.pdf` HAS A CIPHER-GARBLED TEXT LAYER.** 876 of 2 035 lines are mojibake; **`grep` for `HMF`, `Ea`, `Arrhenius` or any English word returns ZERO HITS REGARDLESS OF CONTENT.** | The full glyph cipher is solved and published in `kocadagli2016foodchem_extraction.md` §0.1. ⚠️ **`research_round3_channels.md` §A.1's corpus-wide HMF sweep — "I grepped all 101 PDFs … Only six files pass" — could not have seen this file. It contains an HMF node and eleven HMF-bearing rate constants.** Any other grep-derived "verified negative" on this file is void. |
| **3** | ⚠️ **`Goncouglu2016.pdf` is a 2017 paper.** | Accepted 30 Nov 2016, published *Food Chemistry* **221**, 1911–1922 (**2017**). Cite as **Göncüoğlu Taş & Gökmen 2017**. Dossier named accordingly. |
| **4** | ⚠️ **A near-identical stem exists and is a DIFFERENT PAPER.** | **`Goncouglu2026.pdf`** = *"Effect of extrusion parameters and feed composition on … pea protein-enriched corn extrudates"*, **Food Chemistry: X 36 (2026) 104019**, same first author. It is an **extrusion/acrylamide** paper and belongs to the matrix lane, **not** the HMF cluster. `2016` vs `2026` differ by one digit. |
| **5** | ⚠️ **`agcam2022.pdf` is not what round-3 §E expected.** | Identity ✔ (`10.1007/s12161-021-02214-x`, *Food Analytical Methods* **15**:1286–1299) but it is a **single-author, non-Gökmen, Excel-fitted uniresponse study of synthetic 12 % hexose solutions at pH 3.5** — **no reaction network, no 3-DG, no FFC, no α-dicarbonyl of any kind.** **It cannot address the fructose-vs-3-DG partition.** |
| **6** | ★ **The round-3 §A.3.5 "unattributed HMF-sink lead" is RESOLVED.** | PII `S0308814617312852` = **Hamzalıoğlu & Gökmen 2018, *Food Chemistry* 240:354–360, DOI `10.1016/j.foodchem.2017.07.131`.** It is a Gökmen-school paper. The round-3 "do not attribute" caution can be lifted. |
| **7** | ⚠️ Two of the five multiresponse papers publish their activation energies **only in supplementary material that is not on disk**. | Göncüoğlu Taş 2017's Ea live in **Table S4**, not held. Şen 2022 and Kocadağlı-flour 2016 publish **no Ea at all**. |

---

# §2. ★★ THE HEADLINE CORRECTION: THE ROUND-3 PREMISE FOR THE FRUCTOSE LIMB NEEDS AMENDING

`research_round3_channels.md` §0a is **right in its conclusion and wrong in one of its supports.**

### 2.1 What §0a got right
The 3-DG → HMF edge is **not** sufficient on its own. **Three independent matrices, one lab, one
method, three model-discrimination tests, one answer:** deleting the fructose/cation → HMF edge
makes the model under-predict HMF badly.

| paper | test | verbatim outcome |
|---|---|---|
| Kocadağlı JAFC 2016 (glucose melt) | omit `Fru → Int → HMF` (Fig. S5) | *"did not fit to the experimentally observed data … **by no means**"* |
| Kocadağlı Food Chem 2016 (wheat flour) | omit `Fru → Int → HMF` (Fig. 5b) | *"the amount of HMF was **remarkably underestimated** through the formation only via 3-deoxyglucosone"* |
| Göncüoğlu Taş 2017 (hazelnut) | omit `FFC → HMF` (Fig. 4) | *"The predicted values … were found to be **far below the experimental values**"* |

**⇒ `USE` as a hard structural constraint. A 3-DG-only HMF node is falsified three times over.**

### 2.2 ★ What §0a got wrong — the Gürsul Aktağ significance claim

Round-3 §A.2 row 5 calls Gürsul Aktağ 2020 *"the only one where the fructose-vs-3-DG partition was
tested for statistical significance."* **It was not.** Three findings, all in
`gursulaktag2020_extraction.md` §1:

1. **The Methods declare ANOVA/Tukey for the FREE AMINO ACID data only.** No test of any kind is
   described for rate constants. The abstract's "(p < 0.05)" has no described procedure behind it.
2. **At the terminal HMF-forming step the table says the opposite in 2 of 3 juices.**
   `k₈ (3-DG→HMF)` vs `k₁₁ (FFC→HMF)` at 37 °C: orange **2.0 ± 0.68 vs 0.0 ± 0.02**;
   peach **11.8 ± 2.79 vs 0.4 ± 0.22**; apple 0.7 ± 1.98 vs 5.3 `ind` — and the only cell favouring
   the fructose limb is the one the authors marked **indeterminate**.
3. **The supporting "7, 4 and 2 times higher" ratio set mixes two different reaction steps.**
   Reproduced `[D]`: apple uses `k₁₁/k₈` = 7.6; orange and peach use `k₁₀/k₈` = 3.6 and 1.9, where
   `k₁₀` is `FRU → FFC` — **the entry to the limb, not the step that makes HMF.** Using `k₁₁/k₈`
   consistently gives **7.6 / 0 / 0.034**.

**★ PROPOSED EDIT TO `research_round3_channels.md` §A.2 row 5:** delete "the only one where the
fructose-vs-3-DG partition was tested for statistical significance"; replace with *"the abstract
asserts a significant fructose-limb advantage, but the paper describes no test on rate constants,
its terminal HMF-forming constant favours the 3-DG limb in 2 of 3 juices, and the supporting ratio
set mixes two different reaction steps (K5a §2.2)."*

### 2.3 ★★ AND THE REASON THE APPARENT CONTRADICTION IS NOT A CONTRADICTION — the governing finding of this wave

**In every single one of these five networks, the fructose/cation intermediate is UNMEASURED.**

| paper | fructose-limb intermediate | measured? | 3-DG limb species | measured? |
|---|---|---|---|---|
| Kocadağlı JAFC 2016 | `Int` ("undetermined intermediate") | ✘ | 3-DG ✔, 3,4-DG semi-quant | ✔ / ⚠️ |
| Kocadağlı Food Chem 2016 | `Int` | ✘ | 3-DG ✔, 3,4-DG semi-quant | ✔ / ⚠️ |
| Göncüoğlu Taş 2017 | `FFC` | ✘ | 3-DG ✔, 3,4-DG ✔ | ✔ |
| Gürsul Aktağ 2020 | `FFC` | ✘ | 3-DG ✔ (no 3,4-DG node) | ✔ |
| Şen 2022 | `FFC` | ✘ | 3-DG ✔ (no 3,4-DG node) | ✔ |

In a fitted ODE, an unmeasured node's concentration scale is **not identified**; only the product
`k·[node]` is. **A rate constant on an unmeasured node is therefore not commensurable in magnitude
with one on a measured node.** ★ **Göncüoğlu Taş & Gökmen state this themselves**, verbatim
(§3.5 of that paper, transcribed in full in `goncuoglutas2017_extraction.md` §5.3):

> "The rate constants of 3,4-dideoxyglucosone formation from 3-deoxyglucosone (k₁₈) and
> 5-hydroxymethylfurfural formation from 3,4-dideoxyglucosone (k₁₉) **were higher than** the rate
> constants of 5-hydroxymethylfurfural formation through fructofuranosyl cation (k₁₁). **The lower
> reaction rate constants of 5-hydroxymethylfurfural formation through fructofuranosyl cation could
> be attributed to the fact that concentration of fructofuranosyl cation could not be measured
> because of experimental restrictions and the reaction steps … were reduced to one dehydration
> step, in order to simplify the proposed model.**"

**⇒ THE CORRECT FRAMING FOR THE BUILD WAVE:**
> **"The fructose limb wins on FLUX, established by model discrimination in three matrices. The
> 3-DG limb wins on published RATE CONSTANT in four of six comparisons, because its species are the
> measured ones. These are not in conflict; they are two different quantities, and only the flux
> statement is a chemistry claim."**

---

# §3. ★ THE FRUCTOSE-LIMB-vs-3-DG VERDICT TABLE — every paper, every condition

**This is the deliverable the brief asked for: the branch fraction's matrix dependence, bounded by
evidence.** Ratios are of the **terminal HMF-forming rate constants**, `[D]` from each paper's own
table. ⚠️ **Read the whole table with §2.3 in force: these ratios compare a measured node against an
unmeasured one, so their magnitude is not transferable — only their PATTERN is.**

| paper | matrix | amine? | pH | T (°C) | verdict as STATED | k(3-DG limb) / k(fructose limb) **as computed `[D]`** | what actually decides it |
|---|---|---|---|---|---|---|---|
| **Kocadağlı JAFC 2016** | freeze-dried glucose melt | ★ **NONE** | not measured | 160/180/200 | *"HMF is **primarily formed from fructose** dehydration"* | k₅/k₇ = **516 / 59 / 15** → the 3-DG constant is larger, but it is the **only step in the paper the authors set Ea = 0 on** ("does not follow Arrhenius") | model discrimination (Fig. S5) + the 3-DG pool being starved: k₃ (Glc→3-DG) is the **smallest** constant in the table (0.91 ×10⁻³ min⁻¹) |
| **Kocadağlı JAFC 2016 + NaCl** | same + 0.1 M NaCl | NONE | — | 160/180/200 | *"rate constants increase 4 fold"* on the fructose limb | k₅/k₇ = **40 / 16 / 10** | ★ **k₆(NaCl)/k₆(glucose) = 3.91 / 3.88 / 4.06 — flat across 40 °C.** NaCl **switches flux to the cyclic route** and simultaneously **halves k₃** |
| **Kocadağlı Food Chem 2016** | glucose + wheat flour, 7 % moisture | ✔ 37.7 µmol/g (free + bound lump) | not measured | 160/180/200 | *"Formation of HMF from fructose was found to be a **key step**"* | k₁₃/k₁₅ = **0.13 / 41 / 4.2** — **non-monotonic and sign-inverting** | model discrimination (Fig. 5b). ⚠️ k₁₅ falls **77× from 160→180 °C**: this limb is unidentified here |
| **Göncüoğlu Taş 2017** | roasted hazelnut | ✔ 7513 mg/kg | measured, **not reported** | 150/160/170 | *"mainly proceeded via **fructofuranosyl cation** dehydration rather than 3-deoxglucosone"* | k₁₉/k₁₁ = **0 / 235 / 193** | ★ **model discrimination (Fig. 4) — and the authors explicitly say the k comparison is confounded by FFC being unmeasured** |
| **Gürsul Aktağ 2020** | apple / orange / peach juice | ★ **ANOVA-flat: none** | ★ **3.4** | 27 / 37 | *"significantly higher (p < 0.05)"* for the fructose limb | k₈/k₁₁ at 37 °C = **0.13 (apple, and k₁₁ is `ind`) / ∞ (orange) / 29.5 (peach)** | ★ **pool size**: *"fructose … was the dominant form of sugars"*, 2.3× and 9.3× more fructose in apple than orange/peach. **The significance claim itself fails (§2.2).** |
| **Şen 2022** | 5 nuts/seeds, sucrose-rich | ✔ 7–12 g/kg bound Lys | ★ **6.62–7.14** | 160 / 180 | ★ **the dissent**: *"3-DG pathway was **predominant** in HMF formation for all samples"* | k₁₀/k₂ = **∞ ×6, then 1.2 / 21.2 / 217 / 950** | ★ **FFC is drained before it can dehydrate**: k₃ (FFC+bLYS→HP) and k₂₀ (FFC+AA→P4) are 10–300× k₂. **There is no free-fructose node at all in this network.** |
| **Ağçam 2022** | 12 % glucose or fructose solution | ✔/✘ (asparagine on/off) | ★ **3.5** | 90/105/120 | — **no branch measured** | **not applicable** — no 3-DG, no FFC in the paper | but **fructose out-yields glucose 4–6× on HMF** `[D]`, stable across T and across amine on/off |

### 3.1 ★★ THE THREE STRUCTURAL RULES THAT SURVIVE ALL SEVEN PAPERS

**RULE 1 — THE BRANCH IS SET BY PRECURSOR SUPPLY, NOT BY THE TERMINAL RATE CONSTANTS.**
Every paper that names *why* its limb wins names a supply reason: pool size (Gürsul Aktağ), a
starved 3-DG source (Kocadağlı JAFC k₃), or a drained FFC pool (Şen k₃, k₂₀). **No paper argues from
the terminal constant, and in four of six comparisons the terminal constant points the other way.**

**RULE 2 — SUCROSE AND HEXOSE MATRICES ARE NOT TESTING THE SAME PARTITION.**
In sucrose-rich matrices (hazelnut, nuts/seeds, juice) the cation comes from **sucrose cleavage**:
Göncüoğlu Taş measured `FRU → FFC` and found it **zero**, removing the step —
*"Contribution of free fructose to the fructofuranosyl cation formation was not included in the
model … as the rate constants (k₂₇) were found to equal to zero."* Şen 2022's network has **no
fructose node at all**. In hexose matrices (Kocadağlı ×2) the intermediate comes from
isomerisation-derived **free fructose**. ⇒ **"Fructose limb" and "FFC limb" are two different edges
with two different parents. Do not merge them into one branch fraction.**

**RULE 3 — THE 3-DG LIMB'S INTERNAL RATE-DETERMINING STEP IS `3-DG → 3,4-DG`, IN TWO INDEPENDENT
MATRICES, AND IT IS THE STRONGEST TRANSFERABLE STRUCTURAL FACT IN THE CLUSTER.**

| paper | k(3,4-DG→HMF) / k(3-DG→3,4-DG) | verbatim |
|---|---|---|
| Göncüoğlu Taş 2017, hazelnut | 0 (150 °C) / **4.56** (160) / **4.43** (170) `[D]` | *"almost 5 times higher … 3,4-dideoxyglucosone formation from 3-deoxyglucosone was the **rate determining step**"* |
| Kocadağlı Food Chem 2016, flour | **2.17 / 3.50 / 2.18** at 160/180/200 °C `[D]` | *"dehydration of 3-deoxyglucosone to 3,4-dideoxyglucosone is a **rate-limiting step**"* |

★ **And the fructose limb's rate-determining step is in the OPPOSITE position:** Gürsul Aktağ,
verbatim — *"formation of FFC from fructose was found to be the **fast step** and the rate
determining step was the **HMF formation from FFC**."* ⇒ **The two limbs have mirror-image internal
structure. A model that lumps either limb into one edge loses this.**

---

# §4. ★★ THE HMF SINK — the corpus's gap, now partly closed

`research_round3_channels.md` §A.3.5: *"The repo needs an HMF sink and currently has none."*

## 4.1 The eight published sinks, and why seven fail

| source | sink form | status |
|---|---|---|
| Han 2025 (round-3 §A.3.1) | `HMF → melanoidin`, 1st order | **REFUSE** — HPD 50–95 % of estimate; non-monotonic in T; ~10⁶× every other constant |
| Kocadağlı JAFC 2016, glucose | `HMF → P6`, 1st order | **REFUSE** — Ea = 152.8 **± 154.4** (HPD > estimate) |
| Kocadağlı JAFC 2016, NaCl | `HMF → P6`, 1st order | **REFUSE** — **k = 0.0 exactly at 160 °C** ⇒ published Ea 137.8 **cannot be derived from the paper's own table** |
| Kocadağlı Food Chem 2016 | `HMF → P8`, 1st order | **REFUSE** — **k = 0 at 180 °C**; 71 % HPD at 160 °C |
| Kocadağlı Food Chem 2016 | `HMF + AA → P`, bimolecular | ★ **TESTED AND REJECTED BY THE AUTHORS** — *"HMF did not show a good fit when its degradation included amino acids in a bimolecular reaction. Therefore, this step was excluded"* |
| **Göncüoğlu Taş 2017** | `HMF → P7`, 1st order | ⚠️ **PARTLY USABLE** — 12 ± 3.7 (150 °C), 21 ± 11 (160 °C) usable; 103 ± 63.7 (170 °C) REFUSE. **The only first-order HMF sink in the cluster with a usable HPD at more than one temperature.** |
| Gürsul Aktağ 2020 | ★ **none** | HMF is a strictly accumulating terminal species |
| Şen 2022 | `HMF + AA → P3`, bimolecular | **REFUSE as a set** — zero or HPD-dominated in 6 of 10 cells |

## 4.2 ★★ THE SINK THAT WORKS — Hamzalıoğlu & Gökmen 2018

A dedicated study rather than a by-product of a network fit. **Second-order constants derived `[D]`
from the published pseudo-first-order `k′` and the stated 20 mM amino-acid loading, pH 3.5:**

| reaction | k(5 °C) | k(25 °C) | k(50 °C) | Ea (kJ/mol) | verdict |
|---|---|---|---|---|---|
| **HMF + cysteine** | **3.95** | **5.15** | **23.3** M⁻¹ day⁻¹ | 29.639 (reproduced 29.675) | ★ **USE-Q — the corpus's first genuine second-order HMF-sink constant** |
| HMF + arginine | 0.60 | 0.75 | 2.75 M⁻¹ day⁻¹ | 25.457 (reproduced 25.489) | **USE-Q** |
| HMF + lysine | 0.44 | 0.45 | 0.80 M⁻¹ day⁻¹ | 10.018 (reproduced 10.036) | **PRIOR-ONLY** — refit R² = 0.795 |

**Model-free anchors `[M]`, 7 days:**
- ★ **HMF self-degradation at pH 3.5, 5 °C = 0.9 %.** The amine reaction is essentially the whole
  sink in aqueous acid.
- **Cysteine removes 41 % (5 °C), 52.8 % (25 °C), 97.2 % (50 °C).**
- **Cysteine at 5 °C removes more HMF than lysine at 50 °C — amine identity outranks a 45 °C
  temperature gap by 4×.**
- **Selectivity Cys : Arg : Lys = 11.4 : 1.7 : 1 in water, collapsing to 1.23 : 1.20 : 1 in roasted
  coffee.** ★ **A same-method matrix-vs-water pair on a rate ratio — the class `k3` §C.2 complains
  does not exist.**

## 4.3 ★ THE THREE-WAY CONFLICT, AND ITS RESOLUTION

Three Gökmen papers give three different answers about the HMF + amine edge:

| | says |
|---|---|
| Kocadağlı flour 2016 (160–200 °C, low moisture) | the edge **does not improve the fit; excluded** |
| Şen 2022 (160–180 °C, low moisture) | the edge is **in the model** (and is its only HMF sink) |
| Hamzalıoğlu 2018 (5–50 °C, aqueous + coffee) | the edge is **large, measured, and dominant** |

**RESOLUTION — competitive kinetics, not absence.** Kocadağlı's own explanation, verbatim:
*"amino acids were primarily consumed by **other more reactive carbonyl compounds** under the
investigated conditions … It could be proposed that the reaction of amino compounds with HMF is of
**minor importance in a quantitative way**."* Its own Table 1 confirms the competition: at 160 °C
`GO + AA` = 616, `1-DG + AA` = 373, `MGO + AA` = 115, versus a rejected `HMF + AA`. ★ **And Ağçam
2022 independently locates the crossover in temperature**: asparagine **raises** HMF at 90 °C and
**lowers** it at 105 and 120 °C, at pH 3.5 `[D]`.

**⇒ THE SINK IS REAL AND ITS RATE IS SET BY WHAT ELSE IS COMPETING FOR THE AMINE POOL. Model it as
a competitive bimolecular edge, not as a fixed first-order decay.**

---

# §5. ★ CROSS-PAPER CONSISTENCY CHECKS — same lab, same method, different answers

The cluster is unusually well suited to same-lab comparison: **six of the seven papers are from
Hacettepe FOQUS, using the same o-phenylenediamine/quinoxaline LC-MS SIM method, the same Athena
Visual Studio determinant-criterion fitting, and the same HMF HPLC-DAD method at 285 nm.** That
makes the following disagreements *method-controlled*, and therefore attributable to matrix.

| # | quantity | values | spread |
|---|---|---|---|
| 1 | **3-DG : 1-DG ratio** | glucose melt **4.6 / 4.0 / 3.9**; glucose-NaCl **5.6 / 5.2 / 5.3**; wheat flour **10.6 / 7.8 / 6.9**; hazelnut *"3,4-DG ≈ 3-DG/5"* | **2.3× between the amine-free melt and the flour matrix** |
| 2 | **Is 1,2-enediol a required node?** | JAFC melt: **NO** (*"the interconversion was obviously fast"*); wheat flour: **YES** (*"a rate-limiting step"*); hazelnut: **YES**; juice: **YES** | ★ **A topology switch driven by amine presence and physical state — Gürsul Aktağ names it: *"the physical form of the reactants in food systems is a strong determinant in whether the enolisation step was important or not."*** |
| 3 | **Where does MGO come from?** | JAFC melt: **3-DG** (1-DG route → 0); wheat flour: **1-DG** (3-DG route → 0); hazelnut: **1-DG**; Şen: **1-DG**; juice: **3-DG** | ★ **A parent switch. Kocadağlı reconciles it himself: *"the main source of methylglyoxal may quantitatively depend on the amount of precursor α-dicarbonyl compound formed."*** |
| 4 | **Where does 1-DG come from?** | JAFC melt: **fructose only** (2,3-enolisation); wheat flour: **Amadori only** (`Fru→1-DG` = 0); hazelnut: **AP and HP**; juice: **not detected at all** (acidic) | ★ **Four matrices, four answers, one lab.** |
| 5 | **Does glucosone exist?** | juice: **the dominant dicarbonyl** in apple and orange; Şen: **the least abundant**; hazelnut: ★ **not detected at any time-temperature**; melt/flour: present but semi-quantitated | **from "dominant" to "absent"** |
| 6 | **HMF : furfural ratio** | Ağçam pH 3.5, 120 °C: **27.0** (Glu), **72.4** (Fru), **12.4** (Glu+Asn), **37.4** (Fru+Asn) `[D]`; Lee 2022 cake 200 °C (round-3 §A.3.3): **10 : 1** and **6.25 : 1** | **1.2–7× between an aqueous acid medium and a solid cake, two independent labs** — but **the same direction: adding amine lowers the ratio in both** |
| 7 | **Fructose : glucose HMF advantage** | Ağçam measured `[D]`: **3.95 / 5.82 / 4.84 / 5.81×** | ⚠️ **vs the 36× he approvingly cites from Lee & Nagy 1990 under identical conditions (pH 3.5, citric acid). Use 4–6× `[M]`; refuse 36× `[C]`.** |

★ **All seven rows are new instances of `k3` §0's finding that branch fractions are not constants.
Rows 2, 3 and 4 are stronger than that: they show that network TOPOLOGY, not just parameter values,
is matrix-determined.** Şen 2022 §3.3 makes it explicit by fitting **two different networks** to
five nuts, on PCA grounds.

---

# §6. ★★ THE Ea AUDIT — the Ma-2022 and Han-2025 defect classes, tested across all seven papers

The brief asked: *"test each paper's Ea sets for the same defects (reproduce Arrhenius fits from
their own k tables where possible)."* **I did this for every paper that publishes both a k table
and an Ea, and refitted the two that publish k but no Ea.**

## 6.1 Reproducibility scoreboard

| paper | publishes Ea? | Ea reproduce from its own k? | verdict |
|---|---|---|---|
| **Ağçam 2022** | ✔ 9 values + 9 Q₁₀ | ★ **9/9 Ea to ≤0.9 kJ/mol; 9/9 Q₁₀ to ≤0.04** | ★ **the ONLY internally reproducible Ea table in the cluster** |
| **Hamzalıoğlu 2018** | ✔ 6 Ea + 6 A | ★ **6/6 Ea reproduce to 4 decimal places** — but ★ **only 2/6 pre-exponentials do** | Ea `USE-Q`; **4 of 6 A REFUSE** |
| **Kocadağlı JAFC 2016** | ✔ 34 values (2 systems × 17 steps) | ⚠️ **only 6/34 agree within 15 % of a naïve 3-pt refit** — but the published values come from a *global simultaneous* fit, so disagreement is expected; **the size of the disagreement is the diagnostic** | mixed — see §6.2 |
| **Gürsul Aktağ 2020** | ✔ 43 values | ★★ **0 of 43. Zero.** Six are **mathematically underivable** (k fixed at 0.0 at 27 °C ⇒ ln 0); the R² column is impossible for a 2-point fit; **7 outright sign flips**; the best-behaved step is off by 1.8× | ★ **REFUSE THE ENTIRE TABLE** |
| **Göncüoğlu Taş 2017** | ⚠️ in **Table S4, not held** | n/a — but the authors report the range as **0–1174 kJ/mol with six zeros**, and refitting Table 1 shows **six steps contain an exact zero and admit no Arrhenius line**, matching their "six zero" count | ★ **REFUSE — author-disclaimed twice, and 1174 kJ/mol is absurd** |
| **Şen 2022** | ✘ none | n/a — implied 2-pt Ea for `3-DG→HMF` spans **−134 to +220 kJ/mol** across five matrices, **flipping sign three times** | ★ **no Ea exists; and the implied ones falsify any fixed value** |
| **Kocadağlı Food Chem 2016** | ✘ none | n/a — refit gives `Int→HMF` **Ea = −97.6, R² = 0.322** (k falls **77×** from 160→180 °C) | ★ **no Ea exists; do not attribute one** |

## 6.2 ★ The single most informative result of the whole wave

In the **amine-free glucose melt** (Kocadağlı JAFC 2016), the fructose limb is the **only** part of
the paper whose temperature dependence reproduces:

| step | naïve 3-pt refit `[D]` | published global fit | refit R² |
|---|---|---|---|
| `Fru → Int` (glucose) | **100.5** | **100.4** | **1.000** |
| `Fru → Int` (NaCl) | **102.1** | **94.9** | **1.000** |
| `Int → HMF` (glucose) | **145.0** | **151.4** | **1.000** |
| `Int → HMF` (NaCl) | **152.7** | **149.8** | **0.991** |
| **`3,4-DG → HMF` (glucose)** | **−7.0** | ★ **0, author-fixed: *"does not follow Arrhenius equation"*** | **0.189** |

**But in all three REAL-MATRIX systems, that same fructose limb collapses:**

| paper | fructose-limb terminal step | refit `[D]` |
|---|---|---|
| Kocadağlı wheat flour | `Int → HMF` | **Ea = −97.6 kJ/mol, R² = 0.322** (775 → 10.1 → 71.4) |
| Göncüoğlu Taş hazelnut | `FFC → HMF` | **non-monotonic** (0.58 → 0.57 → 2.02), R² = 0.728 |
| Şen nuts/seeds | `FFC → HMF` | **k = 0 in 6 of 10 cells** — no fit possible |

**⇒ CONCLUSION FOR THE BUILD WAVE: there is exactly ONE well-behaved HMF activation energy in the
Gökmen corpus — `Int → HMF` ≈ 145–152 kJ/mol in an AMINE-FREE amorphous glucose melt, over an
UNMEASURED intermediate pool. It reproduces four ways. And every real food matrix destroys it.**

## 6.3 The Ma-2022 plateau-artefact class — where it does and does not apply

Round-3 §A.3.2 refuses Ma's `Ea(5-HMF) = 14.85` because k moved only **1.4×** over 40 °C.
Applying the same test:

| candidate | k fold-change | verdict |
|---|---|---|
| **Kocadağlı flour `1-DG→MGO`, Ea 16.2 `[D]`** | **1.5× over 40 °C** | ★ **PLATEAU ARTEFACT — refuse** (R² = 0.991 on three nearly flat points is not evidence of a small barrier) |
| **Hamzalıoğlu `Coffee-Arg`, Ea 3.275** | **1.22× over 45 °C** | ★ **PLATEAU ARTEFACT — refuse.** Same class as Ma's Ea(MGO) = 1.84 |
| Hamzalıoğlu `HMF-Lys`, Ea 10.018 | 1.8× over 45 °C | **PRIOR-ONLY**, refit R² = 0.795 |
| **Ağçam, all 9 rows** | **3.5–27× over 30 °C** | ★ **NOT plateau-limited. The Ma defect does not apply.** |
| Kocadağlı JAFC `3-DG→3,4-DG` glucose, Ea 36.9 | 2.1× over 40 °C | **USE-Q with the caveat** |

⚠️ **A second defect class, specific to Hamzalıoğlu 2018:** in 4 of 6 systems the rate barely moves
between 5 and 25 °C (**1.02–1.30×**), so **essentially the whole Arrhenius slope rests on the single
50 °C point — the very point at which the authors declare the pseudo-first-order assumption
compromised** (*"since some part of amino acids might be degraded at 50 °C, measured rate constants
become apparent pseudo-first order rate constants"*).

## 6.4 The Han-2025 unidentifiability class — HPD ≥ estimate

| paper | cells with HPD ≥ estimate, or `ind`, or `fixed` | fraction |
|---|---|---|
| Kocadağlı JAFC 2016 | 13 of 102 | 12.7 % |
| Kocadağlı Food Chem 2016 | 10 of 78 | 12.8 % |
| Şen 2022 | ~52 of ~230 | ~22 % |
| **Göncüoğlu Taş 2017** | **~28 of 78** | ★ **36 %** |
| **Gürsul Aktağ 2020** | **42 of 90** (21 `fixed`, 2 `ind`, 19 HPD-dominated) | ★★ **47 %** |

⚠️ **And the stiff-unidentified-pair signature (one constant 10³–10⁶× its neighbours) appears three
more times:** Gürsul Aktağ step 15 peach 27 °C = **94 365.5 ± 58 130**; Şen 2022 k₂₃ sunflower
160 °C = **18 591 ± 19 730**; Göncüoğlu Taş k₁₆ at 160 °C = **33 016, `ind`**. **Same pathology as
Han 2025's k8/k15.**

## 6.5 ★ A NEW DEFECT CLASS THIS WAVE FOUND — a correct Ea bolted to a mis-derived prefactor

`k3` §0.2 records: *"Two constants in `arrhenius_params.yml` follow the same forensic signature — a
confirmed Ea bolted to an invented prefactor. **Two of two audited.**"*

**Hamzalıoğlu & Gökmen 2018 is the third case, and here it is the SOURCE LITERATURE doing it.**
All six Ea reproduce exactly from Table 1; **only two of six pre-exponentials do.** Two diagnosable
errors, both invisible without reproducing the fit:
1. **A sign flip on every negative intercept** (2 of 6): published `A` = **1/(correct A)** to three
   significant figures, for exactly the two systems whose `ln A` is negative
   (HMF-Lys 1.63 vs 0.6156; Coffee-Arg 2.55 vs 0.3926).
2. **The Coffee-Cys and Coffee-Lys pre-exponentials are SWAPPED** (published 12.08 / 9.67; correct
   9.689 / 12.115) — while their Ea are correct and not swapped.

**⇒ Three of three audited cases now show the same signature. RECOMMEND: audit every remaining
`A_value` in `arrhenius_params.yml` against a refit of its cited source's own k table, not just
against the cited Ea.** (`k3` §0.2 already asks for this; this wave supplies the third motivating
case.)

---

# §7. CONSOLIDATED PARAMETER TABLE

## 7.1 — `USE` / `USE-Q` — ingestible with the stated qualification

| # | parameter | value | units | step | conditions | source | prov. | verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | **k(HMF + cysteine)** | **3.95 / 5.15 / 23.3** | **M⁻¹ day⁻¹** | HMF sink, 2nd order | aq., pH 3.5, 5/25/50 °C, [Cys] 20 mM | Hamzalıoğlu 2018 T1 | `[D]` | ★ **USE — the corpus's first 2nd-order HMF-sink constant; also the sulfur lane's HMF–thiol crosslink** |
| 2 | k(HMF + arginine) | 0.60 / 0.75 / 2.75 | M⁻¹ day⁻¹ | " | " | " | `[D]` | **USE-Q** |
| 3 | k(HMF + lysine) | 0.44 / 0.45 / 0.80 | M⁻¹ day⁻¹ | " | " | " | `[D]` | **PRIOR-ONLY** (refit R² 0.795; SD column broken, §4.1 of that dossier) |
| 4 | **Ea(HMF + Cys)** | **29.639** (refit 29.675) | kJ/mol | " | 5–50 °C, 3 pts, refit R² 0.874 | " | `[F]`+`[D]` | ★ **USE-Q** — *"apparent Ea for HMF consumption by cysteine, aq. pH 3.5, 5–50 °C, slope dominated by the 50 °C point"* |
| 5 | Ea(HMF + Arg) | 25.457 (refit 25.489) | kJ/mol | " | refit R² 0.872 | " | `[F]`+`[D]` | **USE-Q** |
| 6 | **A(HMF-Arg), A(HMF-Cys)** | **612.59, 23 980.59** | day⁻¹ | " | " | " | `[F]`+`[D]` | **USE — the only 2 of 6 correct as printed** |
| 7 | **A(HMF-Lys) corrected** | **0.6156** (not 1.63) | day⁻¹ | " | " | " | `[D]` | **USE-Q — the published value is sign-flipped** |
| 8 | **A(Coffee-Arg) corrected** | **0.3926** (not 2.55) | day⁻¹ | " | " | " | `[D]` | **USE-Q — sign-flipped** |
| 9 | **A(Coffee-Cys / Coffee-Lys) corrected** | **9.689 / 12.115** (published pair is swapped) | day⁻¹ | " | " | " | `[D]` | **USE-Q** |
| 10 | **Ea(Fru → Int)** | **100.4 ± 6.6 (glucose); 94.9 ± 8.1 (NaCl)**; refit 100.5 / 102.1, R² 1.000 | kJ/mol | fructose-limb entry | amine-free amorphous glucose melt, 160–200 °C, T_b 180 °C | Kocadağlı JAFC T2 | `[F]`+`[D]` | ★ **USE-Q** — *"apparent Ea for lumped fructose dehydration over an UNMEASURED cyclic intermediate, amine-free, 3 points"* |
| 11 | **Ea(Int → HMF)** | **151.4 ± 34.3 (glucose); 149.8 ± 19.0 (NaCl)**; refit 145.0 / 152.7, R² 1.000 / 0.991 | kJ/mol | fructose-limb terminal | " | Kocadağlı JAFC T2 | `[F]`+`[D]` | ★ **USE-Q — the single best-corroborated HMF Ea in the cluster. Never "the HMF barrier."** |
| 12 | k_b(Fru→Int), k_b(Int→HMF) | 330 ± 22.8, 1.84 ± 0.70 (glucose); 1402 ± 122, 8.79 ± 1.67 (NaCl) | min⁻¹×10³ at T_b 180 °C | " | " | " | `[F]` | **RATIO-ONLY** — `[Int]` unmeasured |
| 13 | ★ **NaCl catalysis ratio, k₆** | **3.91 / 3.88 / 4.06** at 160/180/200 °C | dimensionless | Fru→Int | amine-free glucose melt ± 0.1 M NaCl | Kocadağlı JAFC T1 | `[D]` | ★ **USE — flat across 40 °C; scale-free; the paper's best number** |
| 14 | NaCl catalysis ratio, k₇ | 3.71 / 5.35 / 4.41 | dimensionless | Int→HMF | " | " | `[D]` | **USE-Q** (1.4× spread) |
| 15 | NaCl effect on 3-DG formation | **0.43 / 0.24 / 1.13** | dimensionless | Glc→3-DG | " | " | `[D]` | **USE-Q — ⚠️ SIGN FLIPS at 200 °C. Quote with temperature or not at all.** |
| 16 | ★ **glucose → HMF mole conversion** | **0.4 / 1.6 / 3.5 %** (glucose); **1.4 / 3.1 / 3.7 %** (+NaCl) | mol % | overall HMF yield | 160/180/200 °C, maxima | Kocadağlı JAFC §3 | `[M]` | ★ **USE — dimensionless; immune to every response-factor caveat. NaCl enhancement collapses from 3.5× to 1.06× as T rises.** |
| 17 | **RDS ratio, 3-DG limb** | **4.56 / 4.43** (hazelnut, 160/170 °C); **2.17 / 3.50 / 2.18** (flour, 160/180/200 °C) | dimensionless | k(3,4-DG→HMF)/k(3-DG→3,4-DG) | two matrices | Göncüoğlu Taş T1; Kocadağlı FC T1 | `[D]` | ★ **USE — two independent matrices agree on the ordering; magnitude 2–5×** |
| 18 | **Ea(HMF formation), glucose / fructose / +Asn** | **107.7 ± 0.4 / 71.0 ± 1.3 / 98.5 ± 1.0 / 49.0 ± 2.6** | kJ/mol | apparent HMF accumulation | 12 % w/v hexose, **pH 3.5**, 90–120 °C, 3 pts, triplicate | Ağçam T1 | `[F]`+`[D]` reproduced | ★ **USE-Q** — *"apparent lumped accumulation Ea, fitted as 1st-order growth with NO sink term"*; ⚠️ the two `+Asn` rows are net-of-destruction lumps |
| 19 | **Ea(furfural formation)** | **109.7 ± 9.5 / 101.3 ± 2.2 / 120.9 ± 2.2 / 76.7 ± 0.6 / 130.3 ± 0.1 (ascorbate)** | kJ/mol | " | " | Ağçam T2 | `[F]`+`[D]` reproduced | **USE-Q**, same label |
| 20 | **Q₁₀ set** | HMF **2.48 / 1.82 / 2.29 / 1.51**; furfural **2.53 / 2.35 / 2.77 / 1.91 / 3.00** | dimensionless | " | " | Ağçam T1/T2 | `[F]`+`[D]` reproduced exactly | ★ **USE — dimensionless and immune to the k-unit question** |
| 21 | k(HMF), k(furfural) | 12 + 15 values, 52.3–833.7 and 0.033–672.6 | min⁻¹×10⁻⁴ | apparent formation | " | Ağçam T1/T2 | `[F]` | **USE-Q** |
| 22 | **fructose : glucose HMF advantage** | **3.95 / 5.82 / 4.84 / 5.81×** | dimensionless | HMF yield | pH 3.5, 90/105/120 °C, ± Asn | Ağçam §Results | `[D]` | ★ **USE** — ⚠️ **refuse the 36× he cites from Lee & Nagy 1990** |
| 23 | k(HMF sink), hazelnut | **12 ± 3.7 (150 °C), 21 ± 11 (160 °C)** | min⁻¹×10³ | HMF → P7, 1st order | roasted hazelnut | Göncüoğlu Taş T1 | `[F]` | **USE-Q — the only usable 1st-order HMF sink at cooking temperature.** REFUSE the 170 °C cell (103 ± 63.7) |
| 24 | k(3,4-DG→HMF), hazelnut | 390 ± 111 (170 °C only) | min⁻¹×10³ | 3-DG-limb terminal | " | " | `[F]` | **USE-Q** — REFUSE 150 °C (=0) and 160 °C (134 ± 127) |
| 25 | k(3-DG→HMF), Şen | 6 of 10 cells well determined; range **37.53–2528.66** | min⁻¹×10³ | lumped 3-DG→HMF | 5 nuts/seeds, pH 6.6–7.1, a_w 0.43–0.62, 160/180 °C | Şen T2 | `[F]` | **USE-Q — ⚠️ LUMPS 3-DG→3,4-DG→HMF; not commensurable with #17, #24** |
| 26 | k(3-DG→HMF), juice | **2.0 ± 0.68 (orange), 11.8 ± 2.79 (peach)**, 37 °C | week⁻¹×10³ | lumped 3-DG→HMF | pH 3.4, 37 °C | Gürsul Aktağ T1 | `[F]` | **USE-Q** — the only 2 well-determined HMF-forming cells in that paper; also lumped |
| 27 | k₁ sucrose cleavage, hazelnut | **6.9 ± 0.8 / 15 ± 3.1 / 22 ± 2.0** | min⁻¹×10³ | SUC → GLC + FFC | 150/160/170 °C | Göncüoğlu Taş T1 | `[F]` | **USE** — tight HPD at all 3 T; cross-verified against the running text |
| 28 | k₁ sucrose hydrolysis, juice | 36.3–147.4 | week⁻¹×10³ | SUC → FRU + GLU | pH 3.4, 27/37 °C | Gürsul Aktağ T1 | `[F]` | **USE**; ★ **corrected Ea = 82.6–94.9 kJ/mol `[D]`, NOT the published 48–56** |
| 29 | k₆ (Glc→3-DG), wheat flour | 4.50 ± 0.9 / 12.2 ± 4.1 / 27.6 ± 7.8 | min⁻¹×10³ | 3-DG source | 7 % moisture, 160/180/200 °C | Kocadağlı FC T1 | `[F]` | **USE — the best-behaved row in that paper (refit R² 0.999, Ea 71.8 `[D]`)** |
| 30 | 3-DG : 1-DG ratios | **4.6±0.5 / 4.0±1.2 / 3.9±1.0** (melt); **5.6±0.3 / 5.2±1.5 / 5.3±1.2** (+NaCl); **10.6±2.8 / 7.8±2.7 / 6.9±1.5** (flour) | dimensionless | — | 160/180/200 °C | Kocadağlı ×2 | `[M]` | **USE — 2.3× matrix spread, same lab, same method** |

## 7.2 — `PRIOR-ONLY` — weak, directional, or single-condition

| item | value | why demoted |
|---|---|---|
| Ea(HMF-Lys) 10.018; Ea(Coffee-Cys) 11.55; Ea(Coffee-Lys) 12.33 | kJ/mol | small signal (1.8–2.1× over 45 °C); Lys refit R² 0.795 |
| Ea(HMF, Fru+Asn) = 49.0 | kJ/mol | smallest fold-change (3.5×), strongest destruction term ⇒ most likely a net lump |
| 3-pt Ea refits of Kocadağlı wheat flour, steps 1, 3, 6, 8, 12, 18 | 43.7–117.5 kJ/mol `[D]` | derived, never published by the authors; use only where monotone and R² ≥ 0.88 |
| 3-pt Ea refits of Göncüoğlu Taş, steps 1, 2, 6, 25 | 90.6–216.4 kJ/mol `[D]` | same; ⚠️ 212–216 is far above the food norm |
| Kocadağlı JAFC Table 2 Ea, steps 3, 4, 8–12, 14 | as printed | not HMF-lane; carry the §6.2 refit column as a companion uncertainty |

## 7.3 — ★ `REFUSE` — with the reason recorded so a later wave does not re-ingest

| item | value | ★ reason for refusal |
|---|---|---|
| ★★ **Gürsul Aktağ 2020 Table 2 — ALL 43 activation energies** | −230.31 to +121.55 kJ/mol | **NOT ONE reproduces from the paper's own Table 1.** The six HMF-step Ea are **mathematically underivable** (k = 0.0 `fixed` at 27 °C ⇒ ln 0). The R² column is impossible for a 2-point fit (must be 1.000; prints 0.00–1.00, with **R² = 0.51 repeated on exactly the six underivable rows**). **Seven sign flips** on rows that are computable. |
| ★ **Göncüoğlu Taş 2017 — any Ea for this DOI** | range **0–1174 kJ/mol, six zeros** | Values live in **Table S4, not in the file on disk.** Authors disclaim Arrhenius **twice** (abstract + conclusion). **1174 kJ/mol is chemically absurd**; the six zeros correspond to steps whose k contains an exact zero (ln 0 undefined). **PROHIBITED DERIVATION.** |
| ★ **Kocadağlı Food Chem 2016 — any Ea** | — | **The paper publishes none.** Only a disclaimer: *"Reaction rate constants found in certain steps do not necessarily follow the Arrhenius equation."* |
| ★ **Şen 2022 — any Ea** | — | **The paper publishes none.** Implied 2-pt Ea for `3-DG→HMF` spans **−134 to +220 kJ/mol** across five matrices with three sign flips `[D]`. |
| ★ **Kocadağlı JAFC Ea(3,4-DG→HMF)** | 0 (glucose) / 169.8 (NaCl) | Glucose: **authors set Ea = 0 with the footnote "does not follow Arrhenius equation"**; k is non-monotonic (160→110→137). NaCl: published 169.8 vs refit 94.2 = **1.80× discrepancy** on a trace with R² = 0.997. |
| ★ **Kocadağlı JAFC Ea(HMF→P6)** | 152.8 ± 154.4 / 137.8 ± 24.6 | Glucose: **HPD > estimate**. NaCl: **k(160 °C) = 0.0 exactly ⇒ not derivable from the paper's own table.** |
| ★ **Kocadağlı JAFC Ea(Glc↔Fru)** | 151.5 ± 99.8 / 146.4 ± 104.2 / 263.6 / 280.8 | ★ **AUTHOR-DECLARED NOT A BARRIER**: *"the activation energy estimated here should not be considered as measure of an energy barrier for the reaction."* HPD 66–71 % of estimate; both k_b are `ind`. |
| ★ **Hamzalıoğlu Ea(Coffee-Arg)** | 3.275 kJ/mol | Chemically impossible barrier from a **1.22× k spread over 45 °C**. Same defect class as Ma 2022's Ea(MGO) = 1.84 (round-3 §A.3.2). |
| ★ **Hamzalıoğlu A(HMF-Lys) = 1.63 and A(Coffee-Arg) = 2.55 day⁻¹** | — | **Sign-flipped intercepts.** Published = 1/(correct) to 3 s.f. Correct: **0.6156** and **0.3926**. |
| ★ **Hamzalıoğlu A(Coffee-Cys) = 12.08 and A(Coffee-Lys) = 9.67 day⁻¹** | — | **The two are SWAPPED.** Correct: 9.689 and 12.115. |
| ★ **any Hamzalıoğlu Ea extrapolated above 50 °C** | — | Data span **5–50 °C**; cooking is 120–200 °C. 3-point line, R² 0.80–0.91, slope resting on the single point the authors declare compromised. **PROHIBITED DERIVATION**, register of `k3` §C.1. |
| ★ **Ağçam ascorbic-acid t₂ = 15.2 / 11.3 / 9.3 min** | — | **Inconsistent with its own k by ~10⁴×** (ln2/k = 210 045 / 28 408 / 7 823 min), and **not rescuable by any unit reinterpretation** (as 1/min: 21.0 / 2.8 / 0.8). |
| ★ **Ağçam precursor-based R² columns (n = 1, n = 2), both tables** | — | **Degenerate** — numerically identical to the n = 0 mass-formation column in essentially every row, because `[Sug]₀ = 12 %` keeps `(1 − e^{−kt})` in its linear regime. **24 of 40 cells are a repeated value.** |
| ★ **"36× fructose over glucose" (Lee & Nagy 1990)** | 36× | `[C]`, and **contradicted 6–9×** by the measurement of the paper that cites it, under identical stated conditions (pH 3.5, citric acid): **4–6×**. |
| ★ **"fructose 90 % / glucose 10 % of HMF"** | — | `[C]` from Perez Locas & Yaylayan 2008, quoted in Ağçam. **Ağçam measured no 3-DG and no FFC. DO NOT ATTRIBUTE THIS SPLIT TO `10.1007/s12161-021-02214-x`.** |
| ★ **Gürsul Aktağ abstract "significantly higher (p < 0.05)" for the fructose limb** | — | **No test on rate constants is described; the table contradicts it in 2 of 3 juices; the ratio set mixes two reaction steps (§2.2).** |
| ★ **Gürsul Aktağ HMF LOD/LOQ "10 mg/L and 30 µg/L"** | — | **Internally impossible** — an LOD three orders above its own LOQ and above the calibration range. |
| ★ **Kocadağlı flour Ea(1-DG→MGO) = 16.2 `[D]`** | — | **Plateau artefact**: k moves 1.5× over 40 °C. Do not promote the refit. |
| **Table cells with HPD ≥ estimate, `ind`, or `fixed`** | 13 % (Kocadağlı ×2) → **47 %** (Gürsul Aktağ) | interval spans zero; same pathology the repo refused in Knol (`k3` §C.6) and Han 2025 |

## 7.4 — `STRUCTURAL` — the shape constraints, no number to fit

| # | constraint | evidence |
|---|---|---|
| **C1** | ★ **A 3-DG-only HMF node is falsified.** Deleting the fructose/cation edge under-predicts HMF badly. | **3 independent matrices**, Kocadağlı ×2 + Göncüoğlu Taş |
| **C2** | ★ **The fructose-limb intermediate (`Int`/`FFC`) is UNMEASURED in all five networks; its rate constants are identified only up to the pool scale.** | author-declared in Göncüoğlu Taş §3.5 and Gürsul Aktağ §3.2(i) and Şen §3.4 |
| **C3** | ★ **`3-DG → 3,4-DG` is the RDS of the 3-DG limb; `3,4-DG → HMF` is fast (2–5×).** | hazelnut + wheat flour |
| **C4** | ★ **`FFC → HMF` is the RDS of the fructose limb; `FRU → FFC` is fast.** — **the mirror image of C3** | Gürsul Aktağ §3.2(ii) |
| **C5** | ★ **In sucrose matrices, FFC comes from SUCROSE CLEAVAGE, not from free fructose** (`FRU→FFC` fitted to zero and removed). | Göncüoğlu Taş §3.5; Şen 2022 has no fructose node at all |
| **C6** | ★ **NaCl switches glucose degradation from the open-chain α-dicarbonyl route to the cyclic route**: k(Fru→Int) ×3.9–4.1, k(Glc→3-DG) ×0.24–0.43. | Kocadağlı JAFC |
| **C7** | ★ **HMF self-degradation is negligible in aqueous acid (0.9 % / 7 d / 5 °C); the amine reaction IS the sink.** | Hamzalıoğlu §3.1 |
| **C8** | ★ **The HMF + amine sink runs at 5 °C — no temperature threshold.** But at 160–200 °C in flour it is out-competed by GO, 1-DG and MGO for the amine pool. | Hamzalıoğlu; Kocadağlı FC §3.7 |
| **C9** | ★ **Amine identity outranks temperature: Cys : Arg : Lys = 11.4 : 1.7 : 1 in water. Cysteine at 5 °C beats lysine at 50 °C.** | Hamzalıoğlu §7.1 |
| **C10** | ★ **That 11.4× selectivity collapses to 1.23 : 1.20 : 1 in a low-moisture roasted-coffee matrix** — a same-method matrix/water pair on a rate ratio (`k3` §C.2's missing class). | Hamzalıoğlu |
| **C11** | ★ **HMF + cysteine stoichiometry is NOT 1:1** — a thiol-Michael adduct, an amine-Michael adduct, a Schiff base and a 2 HMF : 1 Cys adduct all confirmed at Δ < 1 ppm. | Hamzalıoğlu §3.2 |
| **C12** | ★ **Ascorbic acid produces NO HMF at all** (3 temperatures, 90 min, LOD 1.36 µg/L) **but IS a major furfural source.** | Ağçam S1/S2 |
| **C13** | ★ **Furfural IS formed from hexoses alone at pH 3.5, with no pentose and no cation catalyst** — contradicting Gökmen & Şenyuva 2006, which the same paper cites. | Ağçam S5 |
| **C14** | ★ **Asparagine's effect on HMF CROSSES SIGN between 90 and 105 °C** (raises at 90, lowers above). | Ağçam `[D]` |
| **C15** | ★ **At pH 3.4 and 4–37 °C the Maillard reaction does not run (amino acids ANOVA-flat); sugar dehydration alone accounts for the α-dicarbonyls and HMF.** | Gürsul Aktağ §3.1 |
| **C16** | ★ **1-DG is absent at acidic pH** (2,3-enolisation needs alkali) — a clean pH sign constraint for `k3` §B.2. | Gürsul Aktağ §3.2(iii) |
| **C17** | ★ **Whether 1,2-enolisation is a required node is matrix-determined** — omitted in the amine-free melt, required in flour, hazelnut and juice. | 4 papers; Gürsul Aktağ names the mechanism |
| **C18** | ★ **The parent of MGO, and the parent of 1-DG, each switch between matrices in the same lab.** | §5 rows 3–4 |
| **C19** | ★ **Network TOPOLOGY, not just parameters, is matrix-determined** — Şen fits two different networks to five nuts on PCA grounds. | Şen §3.3 |
| **C20** | ⚠️ **Volatile-loss artefact**: GO/MGO/DA degradation constants *decrease* with temperature in both Kocadağlı papers, author-attributed to headspace escape. Three Ea are set to 0 for this reason. **Transport, not chemistry** — the Lee-2024 defect class of round-3 §A.3.4. | Kocadağlı ×2 |
| **C21** | ⚠️ **3-DG and 3-deoxygalactosone were NOT chromatographically resolved** (author-declared). Applies to **every 3-DG number in all five Gökmen multiresponse papers** — same method throughout. | Kocadağlı FC §3.4 |
| **C22** | ⚠️ **3,4-DG, 1-DG and glucosone are SEMI-QUANTITATED** against the 3-DG (or glucosone) response factor. **3,4-DG is the immediate precursor of HMF on the 3-DG limb**, so k(3-DG→3,4-DG) and k(3,4-DG→HMF) both inherit an unknown scale. Flag `absolute_concentration: false`. | Kocadağlı ×2 (author-declared), Gürsul Aktağ (threosone) |

---

# §8. ★ THE HMF NODE DESIGN THE EVIDENCE SUPPORTS

## 8.1 The four published architectures

| architecture | who | HMF sink |
|---|---|---|
| **A.** `d[HMF]/dt = k_a[3,4-DG] + k_b[Int] − k_c[HMF]` | Kocadağlı JAFC 2016; Kocadağlı Food Chem 2016; Göncüoğlu Taş 2017 | 1st order, unnamed products |
| **B.** `d[HMF]/dt = k_a[3-DG] + k_b[FFC]` | Gürsul Aktağ 2020 | ★ **none** |
| **C.** `d[HMF]/dt = k_a[FFC] + k_b[3-DG] − k_c[HMF][AA]` | Şen 2022 | **bimolecular in amine** |
| **D.** `d[HMF]/dt = k₇[Fru] + k₈[3-DG] − k₁₅[HMF]` | Han 2025 (round-3 §A.3.1) | 1st order, to melanoidin |

**Every one of the four agrees on the topology of the SOURCES: exactly two parallel first-order
inputs, one from the 3-DG/3,4-DG chain and one from the fructose/cation chain.** ★ **Four
independent groups, four matrices, one source topology. That is the strongest architectural
agreement in the K5a corpus and it should be adopted without modification.**

## 8.2 ★ RECOMMENDED NODE

```
                     k_dg1              k_dg2
  glucose  ──────►  3-DG  ──────────►  3,4-DG  ──────────►  HMF
                     ▲     (RDS, C3)              (fast, 2-5x k_dg1)
                     │                                       ▲
              (Amadori/Heyns,                                │  k_ffc  (RDS of this limb, C4)
               matrix-dependent, C18)                        │
                                                             │
  fructose / sucrose  ──────────►  FFC or Int  ──────────────┘
                        k_ent      (UNMEASURED — C2)
                       (fast, C4)
                                                             │
                                                             ▼
                                        ┌────────────────────┴──────────────────┐
                                        │  k_amine·[HMF]·[amine]   (dominant)   │
                                        │  k_self·[HMF]            (≤1 % / 7 d  │
                                        │                          at pH 3.5)   │
                                        └───────────────────────────────────────┘
```

### The edges the evidence supports, and what to attach to each

| edge | ship it? | measured constant available | caveats that MUST travel |
|---|---|---|---|
| **`3-DG → 3,4-DG`** | ★ **YES — both endpoints measured** | Kocadağlı FC k₁₂ (46.5/117/139), JAFC k₄ (23.1/30.5/49.3 and 10.6/43.3/101), hazelnut k₁₈ (4.27/29.4/88.1) | **RDS of the limb (C3).** 3,4-DG is **semi-quantitated (C22)**; 3-DG may contain 3-deoxygalactosone (C21) |
| **`3,4-DG → HMF`** | ★ **YES** | JAFC k₅, flour k₁₃, hazelnut k₁₉ | **fast, 2–5× k(3-DG→3,4-DG).** ⚠️ **NO usable Ea exists for this edge in any paper** (§7.3) |
| **`FFC/Int → HMF`** | ★ **YES — required by C1** | JAFC k₇ **only** (the one place it is well determined) | ★★ **RATIO-ONLY. `[FFC]` is unmeasured (C2), so ship this edge with the pool concentration folded into the constant, or as an explicit lumped `fructose → HMF` pseudo-step. NEVER as a transferable elementary constant.** |
| **`fructose/sucrose → FFC/Int`** | ★ **YES** | JAFC k₆ (the only reproducible one, R² = 1.000 both systems); hazelnut k₁ for the sucrose route | **fast (C4).** In sucrose matrices the parent is SUCROSE, not free fructose (C5) |
| **`3-DG → HMF` as a single lumped edge** | ⚠️ **ONLY if 3,4-DG is not carried** | Gürsul Aktağ k₈, Şen k₁₀ | **These lump two steps and are NOT commensurable with the two-step constants above.** Pick one representation and state it. |
| **`HMF + amine → adduct`** | ★ **YES — this should be the primary sink** | ★ **Hamzalıoğlu, 2nd order, 3 amino acids, 3 temperatures, `[D]`** | **Amine identity outranks temperature by 4× (C9); matrix collapses selectivity (C10); stoichiometry ≠ 1:1 (C11); it is OUT-COMPETED by GO/1-DG/MGO for the amine pool at 160–200 °C (C8); Ea must not be extrapolated above 50 °C** |
| **`HMF → products` (1st order)** | **YES, small** | hazelnut k₂₆ at 150/160 °C only | **≤1 %/7 d at pH 3.5, 5 °C (C7); and it is not additive with the amine sink — amine blocking SUPPRESSES self-dimerisation** |
| **`Fru → 3-DG`** | ★ **NO** | — | fitted to **zero in every case** (Kocadağlı JAFC) |
| **`ascorbate → HMF`** | ★ **NO** | — | **measured zero at three temperatures** (C12). Ascorbate is a **furfural-only** precursor |
| **`melanoidin` as the named HMF product** | ⚠️ **only if browning is modelled** | Han 2025 only, and its k₁₅ is unidentifiable | Every Gökmen paper routes HMF to an unnamed `P`; **none of the five models melanoidins at all** |

### 8.3 ★ THE FIVE THINGS THE HMF NODE MUST NOT DO

1. ★ **MUST NOT carry a fixed branch fraction between the two limbs.** Six papers, six matrices,
   verdicts spanning "fructose limb dominant" to "3-DG limb dominant by ∞" — and the deciding
   factor is precursor supply, not the terminal constants (§3.1 Rule 1).
2. ★ **MUST NOT carry an activation energy on either HMF-forming edge except the amine-free
   `Int→HMF` value (145–152 kJ/mol), and that one must be labelled as an amine-free glucose-melt
   lump over an unmeasured pool.** Every real matrix destroys it (§6.2). **Zero of 43 Gürsul Aktağ
   Ea reproduce; Göncüoğlu Taş's are not in the file and are author-disclaimed; two papers publish
   none at all.**
3. ★ **MUST NOT compare a rate constant on `FFC`/`Int` with one on `3-DG` as if they were
   commensurable** (C2).
4. ★ **MUST NOT model the HMF sink as a fixed first-order decay.** It is competitive, bimolecular,
   amine-identity-dependent by ≥11×, and matrix-suppressed (C7–C11).
5. ★ **MUST NOT merge "fructose limb" and "FFC limb" into one branch** — different parents in
   hexose vs sucrose matrices (C5, §3.1 Rule 2).

---

# §9. PROPOSED FIT vs HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR**

> **STATUS: DRAFT. NOT APPLIED. No benchmark file was created, renamed or modified. No repo file
> was touched. This section is a proposal for the orchestrator to accept, amend or reject.**

**Rules inherited from `k3` §D.0, and obeyed here:** (1) no dataset in both columns; (2) the 21
frozen hold-outs of `k3` §D.1 are untouched; (3) the module holds out at least one dataset.
**A fourth rule specific to this module:** *no dataset whose Ea table failed the §6 reproducibility
audit may enter the FIT panel on the strength of that Ea; it may enter only on its concentrations
or its per-temperature k.*

## 9.1 — MODULE: HMF / CARAMELIZATION (new)

| dataset | **proposed role** | one-line reasoning |
|---|---|---|
| **Kocadağlı JAFC 2016, glucose system, Table 1 + Table 2 (steps 4, 6, 7, 8, 11, 12, 14)** | ★ **FIT** | **the only amine-free system in the corpus.** Pure sugar chemistry, un-confoundable by the amino-acid lane; `Fru→Int` and `Int→HMF` are the only HMF constants that reproduce four independent ways (§6.2) |
| **Kocadağlı JAFC 2016, glucose-NaCl system** | ★ **HOLD-OUT** | ★ **the sharpest available single-variable perturbation test.** A model fitted on the salt-free system must *predict* a **3.9–4.1× rise in k(Fru→Int)** and a **2.3–4.2× fall in k(Glc→3-DG)** with no free parameter. It also predicts the **mole-conversion collapse from 3.5× to 1.06× as T rises** (§7.1 #16) — a two-sided test the formation term alone cannot pass |
| **Kocadağlı Food Chem 2016 (wheat flour), Table 1** | **FIT** | the amine-on companion at the same three temperatures from the same lab; supplies the `Glc+AA→AP` and Amadori-degradation block the melt lacks. ⚠️ **fit on the k only — the paper has no Ea, and its `Int→HMF` column must not be fitted (77× non-monotonic)** |
| **★ Göncüoğlu Taş 2017 hazelnut — `HMF = 104 ± 0.5 / 238 ± 1.9 / 278 ± 0.7 mg/kg dw @ 120 min, 150/160/170 °C`** | ★ **HOLD-OUT** | ★ **the strongest new hold-out in this module.** Absolute, authentic-standard, three temperatures, a real solid matrix, **and independently re-cited by Şen 2022** — so a pass is corroborated out-of-sample. **Its Table 1 rate constants must NOT enter the fit** (36 % defective cells, no Ea in the file) |
| **★ Gürsul Aktağ 2020 — `max HMF = 16.2 ± 0.7 / 3.8 ± 0.2 / 12.2 ± 0.5 mg/L` @ 37 °C, 24 wk, and ZERO accumulation at 27 °C** | ★ **HOLD-OUT** | ★ **the only food-pH (3.4), storage-temperature, absolutely-quantified HMF anchor in existence, with a measured amine-off control.** It sits **120–170 °C below every fit-panel temperature** — the module's only genuine temperature extrapolation. **The "no accumulation at 27 °C" row is the cheapest and most informative single test in the module: a model with an over-large HMF Ea will fail it immediately.** ⚠️ **Table 2 must never enter either column (§7.3)** |
| **★ Şen 2022 — `max HMF = 247.0 ± 8.3 mg/kg (sunflower, 180 °C, 30 min)` and `42.7 ± 0.4 (pumpkin)`, plus Table 1 (pH 6.62–7.14, a_w 0.43–0.62)** | ★ **HOLD-OUT** | ★ **the matrix-envelope test.** Five matrices, one method, one pH band, one a_w band — a model must reproduce the **5.8× spread between sunflower and pumpkin** from composition alone. It is also **the 3-DG-dominant end of the branch range** (§3), so it scores the branch logic rather than the branch value |
| **Şen 2022 Table 2 (rate constants)** | **NEITHER — reference only** | ~22 % defective cells, no Ea, and `k₁₀` lumps two steps. Use it for the §3 ordering, not as a target |
| **Ağçam 2022 Table 1 (HMF k, Ea, Q₁₀) + the glucose/fructose absolute series** | ★ **FIT** | ★ **the only internally reproducible Ea table in the cluster (§6.1), at food pH 3.5, with a clean amine on/off × 2-sugar × 3-temperature design and triplicate replication.** It is also the module's only 90–120 °C bridge between the 27–37 °C storage data and the 150–200 °C roasting data |
| **★ Ağçam 2022 — the ASCORBIC ACID medium (no HMF; furfural 0.01–16.87 µmol/L)** | ★ **HOLD-OUT** | ★ **a structural falsification test that costs nothing to score.** Any model that routes ascorbate to HMF fails outright (C12); any model that gives ascorbate no furfural route also fails. **Binary, unambiguous, and orthogonal to every fitted parameter** |
| **Ağçam 2022 Table 2 (furfural)** | **HOLD-OUT** | the furfural lane is out of this module's scope but shares the 3-DG node; holding it out keeps the HMF:furfural ratio (12.4–72.4) an out-of-sample quantity |
| **★ Hamzalıoğlu 2018 — high-moisture k′ (Arg, Cys) + Ea + the 0.9 % self-degradation anchor** | ★ **FIT** | ★ **the only measured HMF sink in existence that survives audit.** The node needs a sink and there is no alternative source. ⚠️ **fit the 5–50 °C constants only; the Ea must be flagged `no_extrapolation_above: 323.15 K`** |
| **★ Hamzalıoğlu 2018 — the LOW-MOISTURE (roasted coffee) series, all six k′ + the 24.5/39.7 % no-amine control** | ★ **HOLD-OUT** | ★ **the module's matrix test for the sink, and it is a hard one.** A model fitted on the aqueous system must *predict* that **cysteine's 11.4× selectivity collapses to 1.2×** in a low-moisture matrix (C10) — the same-method matrix/water pair `k3` §C.2 says the corpus lacks. ⚠️ **its Ea are separately REFUSED (Coffee-Arg = 3.275, plateau artefact) and its A values are swapped/sign-flipped (§7.3), so hold out the k′ and the depletion percentages, not the derived parameters** |

## 9.2 — Check against the `k3` §D.0 rules

- **Rule 1 (no dataset in both columns).** ✔ By construction. The two Kocadağlı JAFC *systems* are
  split (salt-free → fit, NaCl → hold-out) along a **single named variable**, exactly as the
  existing Hofmann-1998 pH partition of `k3` §D.1 splits on pH. Hamzalıoğlu is split along the
  **moisture** axis on the same principle.
- **Rule 2 (frozen hold-outs untouched).** ✔ None of the 21 bundles in `k3` §D.1 is a Gökmen-school
  HMF paper; nothing is moved, renamed, re-scoped or promoted.
- **Rule 3 (every module holds out at least one dataset).** ✔ **Seven** hold-outs proposed.
- **Rule 4 (this module).** ✔ No dataset enters FIT on a failed Ea table. Gürsul Aktağ and
  Göncüoğlu Taş enter only on concentrations; Kocadağlı-flour and Şen enter only on k (or not at
  all); **only Ağçam 2022, Hamzalıoğlu 2018 and the Kocadağlı-JAFC glucose system contribute an Ea
  to the fit, and all three passed the §6 audit.**

## 9.3 — ⚠️ Two open questions for the orchestrator

**(i) The Kocadağlı-JAFC glucose/NaCl split assumes the two systems are iso-loaded. They are not:
initial glucose is 47.1 ± 0.66 µmol (salt-free) vs 56.4 ± 0.77 µmol (NaCl), a 19.7 % difference**,
despite the authors' stated intent. Per `k3` §0, a 2× precursor change moves a branch fraction by
3×. **The mole-conversion targets (§7.1 #16) are normalised and are safe; the raw µmol are not.
Recommend scoring the NaCl hold-out on the mole conversions and the k-ratios, not on absolute µmol.**

**(ii) The `3-DG → HMF` edge has two incompatible representations in this literature** — a two-step
chain via 3,4-DG (Kocadağlı ×2, Göncüoğlu Taş) and a single lumped edge (Gürsul Aktağ, Şen).
**The fit panel as proposed uses the two-step form; two of the four hold-outs were fitted by their
authors with the lumped form.** That is not a problem for scoring *concentrations* (the hold-out
targets above are all concentrations, deliberately), but **it does mean no hold-out in this module
can score the 3-DG-limb rate constants directly.** Flag as a declared gap.

---

# §10. DECLARED GAPS FROM WAVE K5a

| # | gap |
|---|---|
| **G1** | ★ **No usable activation energy exists for either HMF-forming edge in any real food matrix.** The only reproducible one (`Int→HMF`, 145–152 kJ/mol) is from an amine-free freeze-dried glucose melt over an unmeasured intermediate pool, and all three real-matrix systems destroy it (§6.2). |
| **G2** | ★ **No HMF-sink constant exists above 50 °C that survives audit.** Hamzalıoğlu covers 5–50 °C; the hazelnut `HMF→P7` covers 150–160 °C but is a lumped first-order decay with no amine dependence. **The 50–150 °C window is empty.** |
| **G3** | ★ **The fructose-limb intermediate has never been measured in any of the five networks.** Until `[FFC]` or `[Int]` is quantified, no absolute rate constant on that limb is recoverable, from any of these papers or their successors. |
| **G4** | **Göncüoğlu Taş 2017's Table S4 (the hazelnut Ea) and Gürsul Aktağ 2020's Table S4 (the comprehensive-model constants) are not on disk.** Both are the only place those numbers exist. ⚠️ Göncüoğlu Taş's are disqualified in advance by the authors' disclaimer and the 1174 kJ/mol bound; **Gürsul Aktağ's comprehensive-model Table S4 has not been audited and should be treated as unknown, not as usable.** |
| **G5** | **No paper in the cluster models melanoidins, browning or colour.** Every HMF sink terminates in an unnamed `P`. The one paper that names melanoidins (Han 2025) has an unidentifiable sink constant. |
| **G6** | **`3-DG` and `3-deoxygalactosone` were not chromatographically resolved in ANY of the five Gökmen multiresponse papers** (author-declared once, same method throughout). Every 3-DG number in this cluster carries that ambiguity. |
| **G7** | **No paper measures both HMF and an amino-acid consumption curve in the same experiment above 50 °C.** Hamzalıoğlu measures HMF only; the multiresponse papers measure a lumped "total AA". **The competitive-sink hypothesis of §4.3 is therefore inferred, not directly measured.** |
| **G8** | **No pH ladder anywhere in the cluster.** Six distinct pH values exist across papers (3.4, 3.5, 3.5, 6.62–7.14, and two unmeasured), but **no single paper varies pH**. Any pH term in the HMF node would be fitted across labs and matrices simultaneously — which `k3` §B.2 already forbids at family level. |
