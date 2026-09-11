# Shi et al. 2022 — EXTRACTION (free versus protein-bound hexanal in soymilk and soy isolate: a reversible non-covalent channel the model does not have)

**Source on disk:** `data/articles/shi2022.pdf` (10 pp.; downloaded 2026-09-11 at this repository's
request). Read 2026-09-11 via `pdftotext -layout`; **Table 1's text layer is scrambled and must not
be used** — it was transcribed from a 400 dpi render instead. Wave B45.

| field | value |
|---|---|
| Title | "Changes of hexanal content in fermented soymilk: Induced by lactic acid bacterial fermentation and thermal treatment" |
| Authors | X. Shi, Z. Hao, R. Wang, Z. Chen, F. Zuo, Y. Wan, S. Guo (China Agricultural University) |
| Venue | Journal of Food Processing and Preservation **46**(5) (2022) e16555 |
| DOI | 10.1111/jfpp.16555 |
| Units | **µg/L, absolute** — HS-SPME/GC-FID against an external standard curve with 2-methyl-3-heptanone as internal standard |
| Thermal treatments | single fixed holds only: **95 °C for 5 min** or **95 °C for 10 min**. No time series, no temperature series. |

## 1. Why this paper was fetched

For a **free-versus-bound hexanal split**. The corpus scores hexanal in protein matrices measured by
headspace SPME, and the engine predicts a total; if a large share of hexanal is protein-bound and
invisible to the headspace, every such row compares two different quantities.

## 2. What is refused

- **No binding constant, no isotherm, no K_d, no B_max, no Scatchard or Klotz analysis.** No hexanal
  was added exogenously; it is endogenous throughout, so no binding experiment was performed. The
  only thermodynamic number in the paper is **−2 kcal/mol**, and it is a **citation to Damodaran &
  Arora 2013**, is a binding *energy* not a constant, and comes with the statement that "their
  binding constants were low" and no value.
- **No usable bound fraction.** The conclusion's *"nearly half of the hexanal was bound to soymilk
  protein"* carries **no number, no table reference and no error**, and cannot be traced to any table.
  It also sits in a sentence that contradicts itself: the same paragraph says *"the bound hexanal was
  not dissociated after heating"* and then *"Further heat treatment is required to dissociate hexanal
  from protein"*. Quoted as printed; not reconciled.
- **Structurally, a bound fraction is unobtainable from this design.** Total hexanal is never
  measured — only the headspace-accessible share, before and after a treatment that itself changes
  the partitioning. Every "bound" figure in the paper is a **difference**, i.e. "hexanal that became
  headspace-accessible because of this heat treatment", not "hexanal that was bound".
- **No rate.** Every heat treatment is one fixed point. There is no intermediate time and no second
  temperature, so no release kinetics and no Arrhenius handle.
- **The system cannot discriminate the candidate mechanisms.** Lipoxygenase was deliberately
  inactivated at the very first step (*"blanched with 300 ml of NaHCO₃ solution (0.04 mol/L) at 80 °C
  for 3 min to deactivate lipoxygenases"*), so enzyme destruction is **pre-empted, not tested**; and
  headspace loss during the 95 °C hold is never quantified.

## 3. What is accepted — Table 2, the one internally consistent dataset

Acidic soy protein isolate hydrolysates, 30 mg/mL protein, pH 4.5, n = 3, µg/L. This table satisfies
its own stated arithmetic in every cell.

| protease (enzyme/protein, m/m) | before heat | **after 95 °C** | released |
|---|---:|---:|---:|
| 0 % | 23 ± 0.59 | **66 ± 1.7** | 43 ± 1.1 ᵃ |
| 0.17 % | 23 ± 0.73 | **80 ± 5.2** | 57 ± 4.4 ᵇ |
| 0.5 % | 29 ± 3.6 | **85 ± 2.9** | 56 ± 0.68 ᵇ |
| 1 % | 28 ± 3.5 | **95 ± 1.8** | 67 ± 1.7 ᶜ |
| 1.7 % | 35 ± 0.44 | **104 ± 8.2** | 70 ± 7.7 ᵈ |

The fermented counterparts are below detection before heat and 14 ± 0.80, 6.7 ± 0.52, 7.0 ± 0.17,
6.4 ± 0.51, 7.0 ± 1.6 µg/L after it.

**Table 1 is partly corrupt and is used only in part.** Its row 2 prints "9 ± 2.9" for acidic SPI and
"2 ± 0.24" for the *L. delbrueckii* arm, both of which violate the table's own stated arithmetic
(39 + 56 = 95 and 0 + 28 = 28 respectively) and both of which look like dropped digits. Verified at
400 dpi: the page genuinely prints those strings. **No substitute value is invented here**; rows 1, 3
and 4 are used and row 2 is not. Its significance letters are also defective — the sequence skips `e`
and assigns `f` to both 56 ± 1.9 and 21 ± 2.4, which cannot both hold.

**The 7S / 11S comparison is figure only** (Fig. 4, no table, no standard deviations), and the y-axis
is broken (0–60, then 150–300). Values as quoted in the running text, 30 mg/mL, pH 4.5:

| | before heat | after 95 °C / 5 min |
|---|---:|---:|
| acidic 7S | 9 µg/L | 37 µg/L |
| acidic 11S | 40 µg/L | **220 µg/L** |

The abstract's claim that 11S binds "approximately five times" what 7S binds appears **only in the
abstract**; the Results give no such ratio and no calculation.

Other endpoints, all figure-only with the text values quoted: soymilk 24 µg/L; fermented soymilk not
detected; heated fermented soymilk ~7 µg/L; whey fraction after heat 23.2 µg/L against 1.6 µg/L for
the non-whey fraction, i.e. **the hexanal off-note sits in the whey**. Hexanal's sensory threshold is
quoted as 4.5 µg/L.

## 4. The finding this repository takes, and what it costs

**The direction is opposite to the model's.** `matrix_sites.py` binds hexanal to lysine amines as a
**one-way covalent adduct**, so more thermal load always means *less* measured hexanal. Shi measures
the reverse in every single comparison — soymilk, isolate, hydrolysates at all five protease levels,
7S and 11S alike: **heating an acidified plant-protein system raises headspace hexanal**, by 2.9× in
acidic SPI and 5.5× in 11S.

**And the sizes are not comparable.** Wave B45's probe P1 ran the engine's own binding block on Shi's
pot: at 30 g/L and 95 °C for 5 min the model binds **0.0108 %** of the hexanal (bracket corners
0.0027–0.0431 %), and even under a 140 °C / 60 min cook it never reaches 0.25 %. Against a measured
+187 % release that is about **four orders of magnitude, with the sign reversed**.

**The resolution is that these are two different channels, and the model has only the smaller one.**
The declared bracket (Meynier's and Anantharamkrishnan's `k₂ ≤ 2.5 × 10⁻⁵ M⁻¹s⁻¹` at 20 °C) measures
*covalent* Schiff-base adduction with lysine, and nothing here disputes it. Shi is measuring
**reversible non-covalent hydrophobic sequestration** — the authors' own mechanism, released by acid
plus heat, priced by their cited source at about −2 kcal/mol with "low" binding constants. That
channel is **absent from the model entirely.**

Shi's own four-state taxonomy, verbatim, is the clearest statement of what a complete model would
need: *"I, free hexanal; II, bound hexanal, which can be removed by lactic acid bacteria during
fermentation but in a small amount; III, bound hexanal, which were not removed by lactic acid
bacteria, but its binding ability to protein was decreased by way of acidification and enzymatic
hydrolysis, and could be released after heating; and IV, bound hexanal, which was tightly bound to
globulin and could not be removed by lactic acid bacteria nor be released under acidic or
high-temperature conditions."* The model has none of II, III or IV.

## 5. What was refused on purpose

Shi licenses an argument this wave declined to use: that HS-SPME under-reports total hexanal by a
history-dependent factor, and that every hexanal row in the corpus could therefore be judged against
a wider tolerance. Applying it would raise this model's headline score **without a single new fact
about its chemistry**, so B45 registered the refusal in advance and holds to it. The bundles are not
edited and no tolerance moves.

What is recorded instead is a **declared measurement-channel mismatch**: in a protein matrix, an
HS-SPME hexanal row measures the headspace-accessible share, the engine predicts the total, and the
gap between them depends on the sample's pH and thermal history in a way this model cannot currently
compute. That is a completeness debt, and naming it is worth more than a looser threshold.

## 6. What would close it

An equilibrium binding experiment: exogenous hexanal added at several concentrations to a fixed
protein concentration, at two or more temperatures and at both pH 4.5 and pH 7, with **total** hexanal
measured on the same aliquot (exhaustive solvent extraction or purge-and-trap to completion) beside
the headspace value — yielding an isotherm, a constant, and its temperature dependence. Entered in
`EXPERIMENTS.md`. Nothing in this paper substitutes for it.
