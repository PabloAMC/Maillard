# Zheng 2023 — EXTRACTION (GSH 0.25 mM and kaempferol 0.25 mM vs methylglyoxal at 1:2, 1:5, 1:10 in **100 mM sodium phosphate pH 7.4, 1.5 % DMSO, 37 °C, 0–48 h**; GSH, **GSSG** and the GSH-MGO adduct all quantified by LC-TQ-MS; a purified-adduct reversibility test; a GSH-vs-kaempferol competition; and MTT EC50s in SH-SY5Y cells)

### THE HEADLINE: this paper prints the corpus's largest thiol × α-dicarbonyl second-order rate constant — **4.1 × 10⁴ M⁻¹ s⁻¹ for N-acetylcysteine + methylglyoxal, with a reverse first-order constant of 7.5 × 10⁻³ s⁻¹** — but it prints them as **CITATIONS to Lo et al. 1994, not as its own measurements**, and both belong to the **ADDUCT branch (reversible hemithioacetal), not to a redox branch.** Zheng 2023 measures **no rate constant of its own at all.** Its own contribution to the B27 question is different and two-edged: unlike Zheng 2022, it **quantified GSSG on a calibration curve and watched it rise steadily over 48 h at 37 °C** — while attributing that rise to autoxidation of GSH by O₂ rather than to methylglyoxal, **without showing an MGO-free control**.

**Source on disk:** `data/articles/Zheng2023.pdf` (8 pp., Food Chemistry: X 20 (2023) 100920; open access CC BY).
Read from the `pdftotext -layout` text layer, which came through clean (single layer, unlike `Zheng2022.pdf`).
**This paper prints NO TABLES in the main text.** Its two tables — **Table S1** (the MRM/SIM acquisition parameters) and **Table S2** (the identification information for kaempferol and its MGO adducts) — are **supplementary and are NOT on disk** (`data/articles/` holds no Zheng2023 supplementary file). **Figures S1 (extracted ion chromatograms), S2 (GSH / GSSG / kaempferol in the competition incubation), S3 (scavenger-alone viability) and S4 (proposed adduct structures) are also off disk.** Figures 1–6 are in the PDF and are plots; per house rule their axis values are recorded as figure-only.
**Repo status before this dossier:** Zheng 2023 is cited nowhere in `src/kinetic_core/`, nowhere in `data/lit/reaction_rules.yml`, and has no extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Comparison of the methylglyoxal scavenging effects of kaempferol and glutathione and the consequences for the toxicity of methylglyoxal in SH-SY5Y cells" |
| Authors | Liang Zheng (corresponding, liang.zheng@wur.nl), Wouter Bakker, Ignacio Miro Estruch, Frances Widjaja, Ivonne M.C.M. Rietjens — Division of Toxicology, Wageningen University and Research, Stippeneng 4, 6708 WE Wageningen, The Netherlands |
| Venue | **Food Chemistry: X 20 (2023) 100920.** Received 10 May 2023; revised 27 September 2023; accepted 2 October 2023; online 6 October 2023. 2590-1575 / © 2023 The Author(s), Elsevier Ltd, CC BY |
| DOI | **10.1016/j.fochx.2023.100920** — printed twice, in the front-matter block on p. 1 (`https://doi.org/10.1016/j.fochx.2023.100920`) and in the Appendix A supplementary-data line on p. 7 |
| Funding / interest | China Scholarship Council, Grant Number 202008510115 (L. Zheng), declared as a competing-interest disclosure |
| Computational content | **None.** No DFT, no molecular modelling, no docking; the "highest electron density at C8" statement in the Discussion is a **citation to Zhu et al. 2019** and is not a number this paper computes or prints. Nothing was excluded under the no-DFT policy. |
| The α-dicarbonyl | **methylglyoxal (MGO)**, 40 % in water, Merck. Only one; no glyoxal, no 3-DG |
| The thiol | **L-glutathione reduced (GSH)**, ≥98 %, Merck. **L-glutathione oxidized (GSSG)**, ≥98 %, was also purchased **and used as a quantitative reference compound** (§3.2) |
| The competitor | **kaempferol**, >99 %, MedChemExpress — a flavonol, not a thiol |
| Companion on disk | `zheng2022_extraction.md` — the same first author, same group, the direct predecessor; this paper explicitly reuses that paper's GSH-MGO adduct *m/z* for its MRM method (§2.6) |

## 1. Why it matters

**Against `results/validation/kinetic_core_b27_prereg.md`, this paper does three things, and only the third is unambiguously helpful.**

**(1) It supplies the batch's only thiol × α-dicarbonyl rate constants in the fast regime — and they are
CITED, not measured, and they are ADDUCT branch.** The Discussion (p. 7) prints, in one paragraph:

> the apparent second order rate constant for the reaction of MGO with kaempferol **6.3 × 10⁻² M⁻¹ s⁻¹**
> (Zhu et al. 2020a); the second order rate constant for the reaction of MGO with **N-acetylcysteine**
> as a thiol scavenger **4.1 × 10⁴ M⁻¹ s⁻¹** (Lo et al. 1994); and the **reverse first order rate
> constant** of the thiol adduct **7.5 × 10⁻³ s⁻¹** (Lo et al. 1994).

All three are attributed to other papers. **None is measured in this study.** All three describe
**hemithioacetal formation and its reverse** — the Introduction says so in terms ("MGO can also react
reversibly with the thiol group of N-acetylcysteine and glutathione to form **hemithioacetals**",
p. 1, citing the same Lo 1994). **There is no redox constant, no disulfide-forming rate, no order for
one, and no barrier for one, anywhere in this paper.**

**(2) It complicates Zheng 2022's negative result — the pre-registration should know this.**
Zheng 2022 reported that over **6 hours** at 37 °C "the formation of GSSG from GSH was limited"
(data not shown). Zheng 2023 ran the same chemistry for **48 hours** with GSSG on its **own
calibration curve against a commercial reference compound** (§3.2, p. 3) and reports the opposite
shape: "A further continuous decrease in GSH levels **accompanied by an increase in GSSG content**
during the subsequent 48 h was observed" (§3.2, p. 3), and, in the competition incubation, "The
reductions in GSH levels during the subsequent 48 h of the incubations were **mainly due to the
formation of GSSG from GSH**" (§3.4, p. 5). **Over 48 hours, the disulfide is not a minor product —
it is the main sink for the thiol.**
The authors attribute this to autoxidation ("likely mainly due to the autooxidation of GSH to GSSG"),
and that attribution may well be right. **But no MGO-free GSH control incubation is reported
anywhere in this paper**, and the α-dicarbonyl is present at up to 2.5 mM. So the experiment that
would separate "O₂ oxidised the thiol slowly" from "MGO oxidised the thiol slowly" was not run, or
was not shown, in either Zheng paper. **The two papers together therefore establish that a thiol +
α-dicarbonyl + O₂ pot at 37 °C makes disulfide slowly and steadily, and neither establishes what
made it.** That is materially weaker counter-evidence than Zheng 2022 alone reads as, and the
pre-registration's record should say so.

**(3) It is the corpus's cleanest demonstration that "the thiol wins the race and loses the war."**
The paper's central finding is a *kinetic-versus-thermodynamic* split that maps directly onto the
repository's sulfur lane. The thiol adduct forms **instantaneously** and is **reversible**; the
flavonoid adduct forms **slowly** and is **irreversible**; over 48 h the flavonoid strips the MGO out
of the thiol adduct. Measured: the GSH-MGO adduct in the three-way competition falls to **16.2 % of
its starting peak area by 48 h** while the kaempferol-monoMGO rises monotonically (§3.4, p. 5); the
purified kaempferol-monoMGO releases **<0.8 %** free kaempferol in 24 h (§3.3, p. 5). **A reversible
thiol adduct is not a sink.** The repository's `k_thioether` channel is an irreversible matrix
thioether; B17 variant (b) tried a reversible disulfide and was refused; B25 tried irreversible
addition to deoxypentosones and was refused. **This paper says the reversible-adduct structure is the
physically correct one for a thiol + α-oxoaldehyde at 37 °C and that it is, by construction, a poor
sink** — which is consistent with all three refusals and is an argument for looking elsewhere, i.e.
for the pre-registration's instinct even though not for its specific mechanism.

**What it does NOT do for the wave.** Nothing here makes an α-dicarbonyl from norfuraneol
(pre-reg §3(a) is untouched). Nothing here is above 37 °C. Nothing here is at pH 4.5. Nothing here
involves a heteroaromatic thiol. The α-dicarbonyl is again an α-**oxoaldehyde** with an aldehyde
carbon, not an alkyl diketone. And the two headline constants are somebody else's.

## 2. Methods as they matter to a model

**Pot A — the GSH/kaempferol × MGO kinetic series (§2.2, p. 2). The one a kinetic model can read.**

- **Reactants.** GSH **or** kaempferol at a **final 0.25 mM** (constant); MGO at a **final 0.5, 1.25 or
  2.5 mM**, giving **scavenger : MGO molar ratios of 1:2, 1:5 and 1:10**. One scavenger per
  incubation.
- **Buffer, co-solvent, pH.** **100 mM sodium phosphate, pH 7.4**, with **1.5 % DMSO in every
  incubation** (kaempferol's 5 mM stock is in DMSO; the DMSO is matched into the GSH incubations so
  the two are comparable). Stocks of MGO, GSH and GSSG were made at 5 mM in the phosphate buffer.
- **Temperature.** **37 °C in a water bath.**
- **Time.** **0, 1, 2, 4, 8, 24 and 48 h.**
- **Atmosphere.** **Not controlled and not stated.** A water bath, no degassing, no inert gas. The
  authors' own explanation for the progressive GSH → GSSG conversion is autoxidation, which is a
  direct admission that O₂ was present and active for the full 48 h.
- **Quench.** **100 µL of stop solution (ethanol containing 2 % acetic acid) added to an equal volume
  of sample**, i.e. a **1:1 dilution into ~50 % ethanol / 1 % acetic acid (mine)**; immediately
  stored at **−80 °C**.
- **How the thiol and the DISULFIDE were quantified — the load-bearing method statement.**
  **LC-TQ-MS**: Shimadzu Nexera XR LC-20AD XR UHPLC + **Shimadzu 8050** triple quadrupole, ESI, with
  **MRM and SIM in ionisation-polarity switching**, for the simultaneous determination of GSH,
  kaempferol and the reaction products. Column **Supelco Discovery HS F5-3, 15 cm × 2.1 mm, 3 µm**;
  A = water + 0.1 % formic acid, B = acetonitrile + 0.1 % formic acid; **0.25 mL/min**; gradient 0–2
  min 100 % A, 2–5 min 100→40 % A, 5–11 min 40→5 % A, 11–14 min 5 % A, 14–14.1 min 5→100 % A,
  14.1–24 min 100 % A. Nebulising gas 3.0 L/min; drying and heating gas 10.0 L/min; interface 300 °C;
  heat block 400 °C. **The acquisition parameters are in Table S1, which is off disk.**
  **§3.2, p. 3, states the calibration status explicitly: "The concentrations of GSH and GSSG were
  determined via calibration curves using their commercially available reference compounds, while
  the amount of GSH-MGO adduct was expressed as peak area due to a lack of reference compounds."**
  So: **GSH and GSSG are absolute molar concentrations; the GSH-MGO adduct is a peak area only.**
- **How the kaempferol adducts were handled.** **LC-TOF-MS** (Agilent 1200 + Bruker micro-TOF),
  Acquity UPLC BEH C18 50 × 2.1 mm 1.7 µm, 0.18 mL/min, **negative** ESI, *m/z* 100–1500, capillary
  +3200 V, nebuliser 2 bar, drying gas 8 L/min at 200 °C. Full-scan acquisition was used for
  quantification; kaempferol was quantified against its own reference and **its products were
  quantified using kaempferol as the reference in extracted-ion mode**, giving total mass recoveries
  of **91.2 % to 107.6 %** (§3.2, p. 4). **That mass-balance closure exists for kaempferol and does
  not exist for GSH** (Flags 5).
- **What was NOT measured.** **MGO itself is never quantified** in any incubation. No hydroxyacetone,
  no lactaldehyde, no α-hydroxyketone, no D-lactate. No dissolved oxygen. No second temperature. No
  second pH. No MGO-free GSH blank.
- **Replication.** Mean ± SEM of **three replications** (Figs. 1, 2, 4 captions).

**Pot B — the purified-adduct reversibility test (§2.3, p. 2).** The **kaempferol-monoMGO adduct,
isolated and purified by preparative LC**, incubated at a **final 5 µM** in **100 mM phosphate buffer
pH 7.4, 37 °C, 24 h**, same quench, analysed for free kaempferol by LC-TQ-MS in SIM. **Note what is
NOT in this paper: the equivalent test was never run on the GSH-MGO adduct.** The reversibility of
the thiol adduct is asserted on a citation (Andreeva et al. 2019, Lo et al. 1994) and inferred from
the competition experiment, not measured directly here.

**Pot C — the competition (§2.4, p. 2).** GSH **0.25 mM** *and* kaempferol **0.25 mM** *and* MGO
**0.25 mM** (a true 1:1:1) in **100 mM phosphate pH 7.4 with 1.5 % DMSO, 37 °C, 0–48 h**, quenched
and analysed as above for kaempferol, GSH, **GSSG**, and both adducts.

**Pot D — the cells (§2.7–§2.8, p. 3).** **SH-SY5Y human neuroblastoma (ATCC CRL-2266)**, DMEM/F12
GlutaMAX + 10 % FCS + 1 % NEAA + 1 % pen/strep, 5 % CO₂, 37 °C; 1.5 × 10⁴ cells/well in 96-well
plates, 24 h to attach, then **24 h exposure** to the pre-incubated or non-pre-incubated mixture
supplemented with 10 % FCS. **Pre-incubations were 48 h at 37 °C in cell-free HBSS.** Viability by
**MTT**: 3 h in HBSS + 10 % FCS + 0.5 mg/mL MTT, formazan dissolved in 100 µL DMSO, absorbance
**562 nm minus a 620 nm reference**, expressed as % of solvent control. Two designs: (i) an MGO
concentration–response ± a fixed 0.25 mM scavenger; (ii) a fixed **2 mM MGO** against scavenger
0.05–0.25 mM, **with and without** the 48 h pre-incubation.

**Statistics.** Mean ± SEM, ≥3 independent experiments (≥4 for Fig. 5); GraphPad Prism 9;
**EC50 by non-linear regression**, compared by Student's *t*; two-way ANOVA with Dunnett's or Tukey's
post hoc; p < 0.05.

## 3. Tables re-typed

**THE MAIN TEXT PRINTS NO TABLES.** The only two tables in the paper are **Table S1** and **Table S2**,
both supplementary and **both off disk**. Table S1 holds every MRM/SIM transition and collision
energy — so **not one *m/z* transition used for GSH, GSSG or the GSH-MGO adduct quantification in
this paper is recoverable from what is on disk** (the *m/z* values quoted below for kaempferol and
its adducts come from the running text of §3.1). Table S2 holds the retention times, formulae and
mass errors for kaempferol's nine adducts; only the fragments quoted in §3.1 survive.

### Every number printed in the running text

| quantity | value | class | where |
|---|---|---|---|
| **GSH lost immediately at t = 0, GSH:MGO = 1:2** | **19.4 % (0.05 mM)** | `[M]` | §3.2, p. 3 |
| **GSH lost immediately at t = 0, GSH:MGO = 1:5** | **26.0 % (0.07 mM)** | `[M]` | §3.2, p. 3 |
| **GSH lost immediately at t = 0, GSH:MGO = 1:10** | **34.7 % (0.09 mM)** | `[M]` | §3.2, p. 3 |
| **GSSG over the subsequent 48 h** | **"a further continuous decrease in GSH levels accompanied by an increase in GSSG content"** — the curves are Fig. 1B, **no numeric value is printed** | `[M]`, figure-only | §3.2, p. 3 |
| the cause the authors assign to that GSSG | "**likely mainly due to the autooxidation of GSH to GSSG**" | `[M]` (attribution, not measurement) | §3.2, p. 3 |
| GSH-MGO adduct trajectory | "a slight tendency to increase in the first few hours followed by a decrease", ascribed to reversibility plus GSH autoxidation; **peak area only, no concentration** | `[M]`, figure-only | §3.2, p. 3; Fig. 1C |
| kaempferol remaining after 48 h, 1:2 / 1:5 / 1:10 | **47.7 % / 24.1 % / 6.9 %** | `[M]` | §3.2, p. 4 |
| total mass recovery, kaempferol + its products | **91.2 % to 107.6 %** | `[M]` | §3.2, p. 4 |
| kaempferol, parent ion | *m/z* **285.1** [M−H]⁻ | `[M]` | §3.1, p. 3 |
| kaem-monoMGO | *m/z* **357.1** [M−H]⁻ (+72), retention **12.0 min** | `[M]` | §3.1, p. 3 |
| kaem-diMGO a / b / c | *m/z* **429.1** [M−H]⁻ (+72 on monoMGO), retentions **11.1, 11.5, 12.1 min**; **a dominant** | `[M]` | §3.1, p. 3 |
| oxidized kaem-monoMGO a / b / c | *m/z* **355.1** [M−H]⁻ (−2 from monoMGO), retentions **12.0, 13.0, 14.0 min** | `[M]` | §3.1, p. 3 |
| oxidized kaem-diMGO a / b | *m/z* **427.1** [M−H]⁻ (−2 from diMGO), retentions **11.5, 12.0 min** | `[M]` | §3.1, p. 3 |
| number of kaempferol × MGO products | **nine major reaction products** after 0.25 mM kaempferol + 2.5 mM MGO, 48 h | `[M]` | §3.1, p. 3 |
| **kaempferol-monoMGO, decrease over 24 h at 5 µM** | **22.3 %** | `[M]` | §3.3, p. 5 |
| **free kaempferol released from the purified adduct in 24 h** | **< 0.8 %, no significant change** → "the monoMGO adduct formation of kaempferol was **not reversible**" | `[M]` | §3.3, p. 5 |
| **GSH-MGO adduct remaining at 48 h in the three-way competition** | **16.2 %** of its starting peak area; the decrease starts after **8 h** | `[M]` | §3.4, p. 5 |
| GSH lost at t = 0 in the competition (1:1:1, 0.25 mM each) | **17.9 %**, "followed by a rapid decrease" | `[M]` | §3.4, p. 5 (Fig. S2A) |
| GSH loss over 48 h in the competition, cause | "**mainly due to the formation of GSSG from GSH**" | `[M]` | §3.4, p. 5 (Fig. S2B) |
| kaempferol lost over 48 h in the competition | **16.3 %** | `[M]` | §3.4, p. 5 (Fig. S2C) |
| **EC50, MGO pre-incubated alone, SH-SY5Y, MTT, 24 h** | **1.37 mM** | `[F]` (non-linear regression) | §3.5, p. 6 |
| **EC50, MGO pre-incubated 48 h with 0.25 mM GSH** | **1.68 mM** (p < 0.01 vs MGO alone) | `[F]` | §3.5, p. 6 |
| **EC50, MGO pre-incubated 48 h with 0.25 mM kaempferol** | **2.05 mM** (p < 0.01 vs MGO alone; also p < 0.05 vs the GSH value) | `[F]` | §3.5, p. 6 |
| GSH and kaempferol alone, up to 0.25 mM | no effect on SH-SY5Y viability | `[M]` | §3.5, p. 6 (Fig. S3) |
| kaempferol protection threshold, 2 mM MGO, after 48 h pre-incubation | **≥ 0.15 mM** kaempferol | `[M]` | §3.5, p. 6 |
| GSH protection, 2 mM MGO | significant **with or without** pre-incubation, no difference between them | `[M]` | §3.5, p. 6 |
| kaempferol protection, 2 mM MGO, **without** pre-incubation | **not significant** | `[M]` | §3.5, p. 6 |
| **k₂, MGO + kaempferol (apparent second order)** | **6.3 × 10⁻² M⁻¹ s⁻¹** | **`[C]`** — Zhu et al. 2020a | §4, p. 7 |
| **k₂, MGO + N-acetylcysteine (second order)** | **4.1 × 10⁴ M⁻¹ s⁻¹** | **`[C]`** — Lo et al. 1994 | §4, p. 7 |
| **k₋₁, reverse first-order constant of the thiol adduct** | **7.5 × 10⁻³ s⁻¹** | **`[C]`** — Lo et al. 1994 | §4, p. 7 |
| flavonoid MGO-trapping sites | **C6 and C8 of the A ring**, C8 highest reactivity | `[C]` — Lv et al. 2011, Zhu et al. 2019 | §4, p. 6 |
| intracellular GSH levels | "generally high and in the **mM range**" | `[C]` — Meister & Anderson 1983 | §4, p. 7 |
| physiological flavonoid levels | "generally in the **low µM range**" | `[C]` — Erlund et al. 2000, 2001 | §4, p. 7 |

**Figure-only** (per house rule, not typed): every point of Fig. 1A–C (GSH, GSSG and adduct vs time
at three ratios — **including the entire GSSG trajectory, which is the quantity this dossier most
wants**), Fig. 2A–C, Fig. 3A–B, Fig. 4, Fig. 5, Fig. 6A–B, and all of Figs. S1–S4.

### Arithmetic on the printed numbers (all mine)

**1. Internal check on the three printed GSH losses — it passes.** 19.4 % × 0.25 mM = 0.0485 ≈ the
printed 0.05 mM; 26.0 % × 0.25 = 0.065 ≈ 0.07; 34.7 % × 0.25 = 0.0868 ≈ 0.09. Consistent to the
printed rounding.

**2. The instantaneous plateau read as a 1:1 equilibrium — and it does NOT close (mine).** Assuming
`MGO + GSH ⇌ adduct`, 1:1, free MGO = total − adduct:

| ratio | [adduct] mM | [GSH]free mM | [MGO]free mM | K_app = [adduct]/([MGO][GSH]) (mine) |
|---|---:|---:|---:|---:|
| 1:2 | 0.0485 | 0.2015 | 0.4515 | **≈ 5.3e2 M⁻¹** |
| 1:5 | 0.065 | 0.185 | 1.185 | **≈ 3.0e2 M⁻¹** |
| 1:10 | 0.0868 | 0.1632 | 2.413 | **≈ 2.2e2 M⁻¹** |

A single 1:1 equilibrium constant would give the **same** K at all three ratios. It **falls by 2.4×
from 1:2 to 1:10**. So either the stoichiometry is not 1:1, or the free-MGO approximation fails
because most MGO is hydrated and unavailable, or the plateau is not a true equilibrium. **Recorded as
a check that fails; none of these three numbers should be shipped.**

**3. Cross-study consistency with Zheng 2022 (mine, and it is striking).** The same arithmetic
applied to Zheng 2022's equimolar MGO run (18.1 % of 0.5 mM GSH consumed) gives **≈ 5.4e2 M⁻¹**.
This paper's lowest-MGO ratio gives **≈ 5.3e2 M⁻¹**. **Two papers, two buffers (25 mM potassium vs
100 mM sodium phosphate), two DMSO states (0 vs 1.5 %), two quenches (aqueous acetic acid vs
ethanolic acetic acid) — and the low-MGO apparent constant agrees to 2 %.** That is a real
reproducibility result for the adduct branch and it is worth recording, with the caveat that both
numbers rest on the same four unstated assumptions (§3 item 2).

**4. The two cited constants, taken as a pair, do not reproduce either paper's own data (mine).**
Lo 1994's forward 4.1e4 M⁻¹ s⁻¹ over its reverse 7.5e-3 s⁻¹ gives **K = 5.5 × 10⁶ M⁻¹**. At the
1:10 ratio (2.5 mM MGO) that predicts the fraction of thiol bound as
`K[MGO]/(1+K[MGO]) ≈ 0.99997`, i.e. **essentially all** of the GSH as adduct. The paper measures
**34.7 %**. Even granting that only ~1 % of MGO is unhydrated (Zheng 2022 §4, cited from refs 37–38),
the effective K would still be ~5.5e4 M⁻¹ and predict ~99 % bound. **The gap is at least two orders
of magnitude.** Possible reasons, none of them tested here: N-acetylcysteine is not GSH; Lo's two
constants may not be a matched forward/reverse pair for the same species; the ethanolic acid quench
may reverse the hemithioacetal before analysis (Flags 4). **This is the strongest single reason not
to import the Lo pair into the repository on Zheng 2023's authority.**

**5. Thiol vs flavonoid forward-rate ratio, as the paper's own citations state it (mine).**
4.1e4 / 6.3e-2 = **6.5 × 10⁵×** in favour of the thiol. The paper's qualitative reading ("the initial
scavenging of MGO by a thiol reagent like N-acetylcysteine or GSH will be **highly favored** over
that by kaempferol") is consistent with it.

**6. Adduct half-life from the cited reverse constant (mine).** ln2 / 7.5e-3 s⁻¹ = **92 s**. That is
the timescale on which a thiol–MGO hemithioacetal falls apart at 37 °C, pH 7.4 — **conditional on
the cited constant and on N-acetylcysteine standing in for GSH**. It is the quantitative content of
"a reversible adduct is not a sink", and it is why the kaempferol adduct eventually wins.

**7. The EC50 shifts as within-study ratios (mine).** GSH: 1.68 / 1.37 = **1.23×**. Kaempferol:
2.05 / 1.37 = **1.50×**. Kaempferol over GSH: 2.05 / 1.68 = **1.22×**. Both scavengers were at the
same 0.25 mM against a 1.37 mM EC50, so **neither is stoichiometrically able to remove the MGO** —
these are partial protections, and the ranking, not the magnitude, is the result.

## 4. Numbers the repository can use

Every measured row below shares: **cell-free, 0.25 mM scavenger, 100 mM sodium phosphate pH 7.4,
1.5 % DMSO, 37 °C water bath, atmosphere uncontrolled, ethanolic 2 % acetic acid quench, LC-TQ-MS or
LC-TOF-MS, mean ± SEM of three replications.** The cell rows are SH-SY5Y, MTT, 24 h exposure after a
48 h cell-free HBSS pre-incubation at 37 °C.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **k₂, N-acetylcysteine + MGO** | **4.1 × 10⁴** | M⁻¹ s⁻¹, order 2 | "under physiological conditions" per the cited title; **conditions not restated by Zheng** | §4, p. 7 — **cited to Lo et al. 1994, *J. Biol. Chem.* 269:32299** | **`level_only` / CITED — NOT a Zheng measurement.** If wanted, extract Lo 1994 directly; do not ship on this anchor |
| **k₋₁, thiol–MGO adduct reverse** | **7.5 × 10⁻³** | s⁻¹, order 1 | as above | §4, p. 7 — cited to Lo et al. 1994 | **`level_only` / CITED** — same caution |
| k₂, kaempferol + MGO (apparent) | **6.3 × 10⁻²** | M⁻¹ s⁻¹, order 2 | as above | §4, p. 7 — cited to Zhu et al. 2020a | **`level_only` / CITED**; a flavonoid, not a thiol |
| **any rate, order or barrier MEASURED in this paper** | **— none** | — | — | — | **absent** |
| **any rate, order or barrier for thiol → disulfide** | **— none, anywhere, cited or measured** | — | — | — | **absent** |
| GSH consumed instantaneously, GSH:MGO 1:2 / 1:5 / 1:10 | **19.4 / 26.0 / 34.7** | % of 0.25 mM GSH | pH 7.4, 37 °C, t = 0 | §3.2, p. 3 | **`measured_ratio`** |
| dose response of the instantaneous capture | **19.4 → 26.0 → 34.7 %** across a **5× MGO span** (0.5 → 2.5 mM) | — | as above | §3.2, p. 3 | **`within_study_ratio`** — a **1.79× (mine)** capture rise for a 5× dicarbonyl rise; strongly sub-proportional, i.e. saturating |
| **GSSG formed over 48 h** | rises continuously; **quantified on a GSSG calibration curve but printed only as Fig. 1B** | mol/L | 0.25 mM GSH, 0.5–2.5 mM MGO, pH 7.4, 37 °C, **aerobic** | §3.2, p. 3; Fig. 1B | **`figure_only`** — the single most useful number in the paper for B27 and it is not printed |
| GSSG as the dominant 48 h GSH sink (competition pot) | "mainly due to the formation of GSSG from GSH" — **no number** | — | 0.25 mM GSH + 0.25 mM kaempferol + 0.25 mM MGO, 48 h | §3.4, p. 5; Fig. S2B (off disk) | **`level_only`** — direction only |
| GSH-MGO adduct remaining at 48 h under kaempferol competition | **16.2** | % of starting peak area | 1:1:1 at 0.25 mM, pH 7.4, 37 °C | §3.4, p. 5 | **`measured_ratio`** — the reversibility demonstration |
| onset of the adduct's decline | **> 8 h** | h | as above | §3.4, p. 5 | `measured_bound` |
| kaempferol-monoMGO reversibility | **< 0.8 %** free kaempferol released in 24 h; adduct itself down **22.3 %** (to isomerisation and oxidation, not hydrolysis) | % | 5 µM purified adduct, 100 mM phosphate pH 7.4, 37 °C, 24 h | §3.3, p. 5 | **`measured_bound`** — an upper bound on reversibility |
| kaempferol remaining at 48 h, 1:2 / 1:5 / 1:10 | **47.7 / 24.1 / 6.9** | % | pH 7.4, 37 °C | §3.2, p. 4 | `measured_ratio` |
| mass balance closure, kaempferol lane | **91.2–107.6** | % | as above | §3.2, p. 4 | **`measured_bound`** — and note **no equivalent closure exists for the GSH lane** |
| EC50, MGO alone / +GSH / +kaempferol | **1.37 / 1.68 / 2.05** | mM | SH-SY5Y, MTT, 24 h, after 48 h pre-incubation with 0.25 mM scavenger in HBSS | §3.5, p. 6 | **`level_only`** (fitted EC50s; a cytotoxicity endpoint, not a chemical rate) |
| EC50 shift ratios | **1.23× (GSH) / 1.50× (kaempferol) / 1.22× (kaempferol over GSH)** | — | as above | (mine) | `within_study_ratio` |
| apparent 1:1 adduct equilibrium constant at 1:2 / 1:5 / 1:10 | **≈ 5.3e2 / 3.0e2 / 2.2e2** | M⁻¹ | as above; assumes 1:1, free MGO = total, all GSH loss is adduct | §3 item 2 (mine) | **`derived_assumption` — DO NOT SHIP**; the three values should be equal and are not |
| agreement of the low-MGO apparent K with Zheng 2022's | **5.3e2 vs 5.4e2 M⁻¹ (2 %)** | — | two papers, two buffers, two quenches, same 37 °C / pH 7.4 | §3 item 3 (mine) | **`within_study_ratio`** (cross-study; both legs are my arithmetic) |
| adduct half-life implied by the cited reverse constant | **≈ 92 s** | s | 37 °C, pH 7.4, **N-acetylcysteine, per Lo 1994** | ln2/7.5e-3 (mine) | **`derived_assumption`** |

### Adduct branch or redox branch?

**Every constant this paper prints is ADDUCT branch. Its only redox observation is a figure.**

| quantity | branch | how it was decided |
|---|---|---|
| **4.1 × 10⁴ M⁻¹ s⁻¹ (NAC + MGO)** | **ADDUCT** | The paper's own Introduction names the product: "MGO can also react reversibly with the thiol group of N-acetylcysteine and glutathione to form **hemithioacetals** (Lo … Thornalley, 1994)" (p. 1). The Discussion re-uses the same citation for the constant. A hemithioacetal is a C–S bond at the aldehyde carbon. **No disulfide.** |
| **7.5 × 10⁻³ s⁻¹ (reverse)** | **ADDUCT** | It is explicitly "the reversible nature of the **thiol adduct** formation with a reverse first order rate constant" (§4, p. 7). It is the hydrolysis of the same hemithioacetal. |
| 6.3 × 10⁻² M⁻¹ s⁻¹ (kaempferol + MGO) | **ADDUCT**, and not a thiol at all | C-nucleophilic attack at the flavonol A-ring C8/C6 (§4, p. 6, cited). |
| the 19.4 / 26.0 / 34.7 % instantaneous GSH losses | **ADDUCT** | They are accompanied by the appearance of the GSH-MGO adduct peak (Fig. 1C), whose identity was established in Zheng 2022 by LC-TOF-MS as the *m/z* 380 conjugate, and which this paper reuses by *m/z* for its MRM method (§2.6). |
| the 16.2 % remaining adduct at 48 h | **ADDUCT**, reversing | The MGO leaves the thiol and reappears on kaempferol; the thiol is released, not oxidised, in this step. |
| **the rising GSSG in Fig. 1B and Fig. S2B** | **REDOX — the real thing, measured, and unprinted** | GSSG was quantified against a commercial reference on its own calibration curve (§3.2, p. 3). Over 48 h it becomes the dominant fate of the thiol (§3.4, p. 5). **The paper attributes it to O₂ autoxidation, not to MGO** — an attribution with no control experiment behind it in this paper (Flags 2). This is the only redox measurement in the batch besides Zheng 2022's unquantified "limited", and **it points the opposite way from that phrase once the window is 48 h instead of 6 h**. |
| the "oxidized kaem-monoMGO" and "oxidized kaem-diMGO" species | **REDOX, but of the FLAVONOID, not the thiol** | *m/z* −2 from their parent adducts, i.e. a two-electron oxidation of the kaempferol–MGO conjugate (§3.1, p. 3; proposed structures in the off-disk Fig. S4). **Important context, not a thiol result:** something in this aerobic pot at 37 °C was oxidising organic substrates over 48 h. It does not identify what. |

**The honest summary of the branch question for this paper.** Zheng 2023 supplies **no** redox rate,
order or barrier. It supplies a **cited** adduct-branch forward and reverse constant pair which
cannot be reconciled with its own measurements (§3 item 4). And it supplies an **unprinted**
measurement of the disulfide's growth over 48 h that the pre-registration would very much like to
have as a number.

## 5. Flags

1. **TEMPERATURE TRANSFER.** Everything here is at **37 °C** against the wave's **140 °C**, a 103 K
   gap, with **no second temperature anywhere in the paper**. What that costs is concrete and
   arithmetic: to carry the cited 4.1e4 M⁻¹ s⁻¹ to 413 K one must import a barrier, and the resulting
   multiplier is **≈ 1.1e2 at Ea = 50 kJ/mol, ≈ 1.5e3 at 80, ≈ 1.1e5 at 122.2 kJ/mol** (all mine,
   Arrhenius). **The modeller chooses three decades of the answer.** Worse, the *reverse* constant
   would also move, by a different and unknown factor, so **even the adduct's equilibrium position at
   140 °C is undetermined by anything here.** And a hemithioacetal with a 92 s half-life at 37 °C
   (§3 item 6) has essentially no lifetime at 140 °C: at 122 kJ/mol its half-life would fall to
   **≈ 1 ms (mine)**. **Any use of these constants at Maillard temperature is a `derived_assumption`,
   and the adduct branch may simply not exist as a persisting species there at all** — which, if
   true, removes the very competitor that made Zheng 2022's negative result look decisive.
2. **The GSSG attribution has no control behind it in this paper.** "Likely mainly due to the
   autooxidation of GSH to GSSG" (§3.2) is an inference. **No MGO-free GSH incubation over 0–48 h is
   reported**, in the main text or (as far as the figure list shows) in the supplementary. Without
   that blank, the paper cannot distinguish O₂-driven from MGO-driven thiol oxidation, and the
   MGO was present at up to 2.5 mM. **This is exactly the control the pre-registration's question
   turns on, in both Zheng papers, and neither ran it.**
3. **THE THREE HEADLINE CONSTANTS ARE NOT THIS PAPER'S.** 4.1e4, 7.5e-3 and 6.3e-2 are citations to
   Lo et al. 1994 (*J. Biol. Chem.* 269:32299–32305) and Zhu et al. 2020a. Zheng 2023 restates them
   in one Discussion paragraph and does not restate their conditions — **no temperature, no pH, no
   buffer, no ionic strength, no method is given for any of the three in this paper.** If the
   repository wants them, **extract Lo 1994 directly**; citing them to Zheng 2023 would put an
   unanchored number into the registry behind a wrong anchor.
4. **The cited constants and the measured data are inconsistent by ≥2 orders of magnitude** (§3
   item 4). The Lo pair implies ~100 % of the thiol bound at 2.5 mM MGO; the paper measures 34.7 %.
   The paper does not notice this. Candidate explanations — the ethanolic acid quench reversing the
   adduct on a 92 s timescale, NAC ≠ GSH, or the two Lo constants not being a matched pair — are all
   untested. **Until it is resolved, neither the constants nor the measured percentages should be
   used to calibrate the other.**
5. **There is no mass balance for the GSH lane.** Kaempferol's lane closes at 91.2–107.6 %; the GSH
   lane has **GSH and GSSG in molar units and the adduct in peak-area units only**, with no reference
   compound. So the paper cannot state, and does not state, what fraction of the lost GSH is adduct
   versus disulfide at any time point. **The one arithmetic that would answer B27's question is
   arithmetically unavailable from the published data.**
6. **The quench is 50 % ethanol with ~1 % acetic acid.** That is a much harsher perturbation than
   Zheng 2022's 2 % aqueous acetic acid, applied to an adduct the paper itself says is reversible
   with a ~92 s half-life. **The measured adduct peak areas are post-quench quantities**; no recovery
   or stability control is reported for the quenched adduct.
7. **1.5 % DMSO is present in every incubation**, including the GSH-only ones (added to match the
   kaempferol arm). DMSO is mildly reducing, is a hydroxyl-radical scavenger, and is not present in
   Zheng 2022's pots — which is one uncontrolled difference between the two papers' GSSG results.
8. **The α-dicarbonyl is again an α-oxoaldehyde, not an alkyl diketone.** MGO has an aldehyde carbon
   and that carbon is where the hemithioacetal forms. The pre-registration's proposed species
   (2,3-pentanedione, 2,4-pentanedione, 3,4-hexanedione) have none. **The dominant branch measured
   here is structurally unavailable to the pot's own diketones.**
9. **The thiol is glutathione and N-acetylcysteine, both aliphatic peptidic thiols.** The disulfide
   share the wave must reproduce belongs to **2-methyl-3-furanthiol**, a heteroaromatic thiol. See
   `zheng2022_extraction.md` Flags 9 for the same point at greater length; it applies identically
   here.
10. **pH 7.4 only, and 100 mM phosphate.** The Whitfield pot is pH 4.5 in 0.5 M phosphate. Both the
    hemithioacetal and any thiol oxidation are thiolate-gated; a ~3 unit pH drop is ~3 orders of
    magnitude on the thiolate fraction, and the two branches need not scale with it identically.
11. **All tables and four of the figures are supplementary and off disk.** Table S1 (**every MRM/SIM
    transition**), Table S2 (adduct identifications), Fig. S1 (extracted ion chromatograms), Fig. S2
    (**the GSH / GSSG / kaempferol trajectories in the competition pot — the figure that would show
    how much GSSG**), Fig. S3, Fig. S4. Enough of S2's content is quoted in §3.4 to recover the three
    percentages typed above; **nothing else is recoverable.**
12. **The EC50s are fitted, not measured, and the scavengers are sub-stoichiometric.** 0.25 mM
    scavenger against a 1.37 mM EC50 cannot remove the MGO; the shifts (1.23× and 1.50×) are partial
    protections and the paper's conclusion rests on their ranking, which is statistically supported
    (p < 0.05 between the two) but small.
13. **What to request from the authors:** (i) **the numeric GSSG trajectories behind Fig. 1B and
    Fig. S2B** — the µM of GSSG at 0, 1, 2, 4, 8, 24 and 48 h at each of the three MGO ratios;
    (ii) **an MGO-free GSH blank** over the same 48 h in the same buffer with the same DMSO, which
    would settle the autoxidation attribution and, with (i), would give the first real bound on
    dicarbonyl-driven thiol oxidation in the corpus; (iii) a GSH-MGO adduct reference compound or a
    response factor, so the GSH lane can be mass-balanced; (iv) an adduct-stability control for the
    ethanolic quench; (v) whether any incubation was ever run above 37 °C.
14. **What this paper does NOT contain:** any main-text table; any rate constant of its own; any
    activation energy; any second temperature; any second pH; any measurement of MGO itself; any
    mass balance for the thiol lane; any printed GSSG number; any MGO-free control; any α-diketone;
    any heteroaromatic thiol; any food matrix; any volatile or headspace measurement; any DFT or
    other computation.
