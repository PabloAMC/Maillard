# Jansson et al. 2020 — EXTRACTION (hexanal at seven temperatures in a wet protein solution: the right experiment, figure-only data — REFUSED as a barrier source, kept as a direction)

**Source on disk:** `data/articles/jansson2020.pdf` (7 pp., typeset; downloaded 2026-09-11 at this
repository's request). Read 2026-09-11 via `pdftotext -layout` plus 300 dpi renders of the four
figure pages. Wave B45.

| field | value |
|---|---|
| Title | "Temperature-dependency of unwanted aroma formation in reconstituted whey protein isolate solutions" |
| Authors | Therese Jansson, Søren B. Nielsen, Mikael A. Petersen, Marianne N. Lund (U. Copenhagen; Arla Foods Ingredients) |
| Venue | International Dairy Journal **104** (2020) 104653 |
| DOI | 10.1016/j.idairyj.2020.104653 |
| System | **3 % (w/v) whey protein isolate** (Lacprodan DI-9224) in Milli-Q water, **pH 7**, 20 mL in a sealed 20 mL microwave vial |
| Temperatures | **50, 60, 70, 75, 80, 85, 90 °C for 10 min**, plus a 5-min series; microwave heating, passive cooling at ≈0.3 °C/s to 50 °C |
| Hexanal units | **µg/L, absolute** — one of six compounds quantified against authentic external standards |

## 1. Why this paper was fetched

Same reason as `li2020_extraction.md`: the lipid lane needs a **hexanal barrier in a hot, wet,
protein-containing system**, between bulk oil (114–122 kJ/mol) and moist nut paste (61–65 kJ/mol).

**Structurally this is the right experiment** — and the only one found so far that is. Hexanal, in
real concentration units, at **seven temperatures spanning 40 K**, in a wet protein solution at pH 7,
with two durations at four of them.

## 2. Why it still cannot yield a barrier — six reasons, in order of severity

1. **Every hexanal value is figure-only.** Fig. 4 carries no data labels, there is no data table in
   the article, and no hexanal number appears in the running text. Supplementary Tables S1 (the 48
   identified volatiles) and S2 (the standard curves) and Fig. S1 (the heating/cooling profile) are
   all cited and **none is in the 7-page PDF or anywhere in `data/articles/`**. There is nothing to
   fit. This is the one blocker that is *fixable* — see §6.
2. **No unheated control in Fig. 4.** The lowest bar is 50 °C/10 min. The authors state volatiles are
   already formed at 50 °C, so a rate needs the increment above an unmeasured baseline. Since the
   whole 50→80 °C block is statistically flat, a baseline of comparable size would swing any apparent
   barrier by an unbounded factor.
3. **There is no measurable lipid substrate.** Table 1, verbatim: "Fat — Below detection limit."
   No oil was added. The hexanal comes from **trace, unquantified residual lipid in the powder**.
   There is no hydroperoxide pool, no known substrate concentration, and no way to distinguish a
   genuinely slow rate from exhaustion of a trace contaminant. A rate law needs a substrate.
4. **The thermal history is not isothermal.** Microwave ramp (profile in the missing Fig. S1) plus
   "passive cooling (approximately 0.3 °C s⁻¹ to 50 °C)". At that rate the 90 °C sample spends
   ~133 s above 50 °C on the way down alone, against a 300 s nominal hold — and that excess **grows
   with set-point**, so the nominal (T, t) pairs understate the hotter samples' load and would bias a
   naive fit toward a *lower* barrier.
5. **The measurand is headspace release from a protein solution, and protein state is the variable.**
   The authors themselves cite Kühn et al. 2006 on volatile–protein binding altering release.
   Irreversible denaturation runs across exactly the fitted range, so measured ≠ formed with a
   temperature-correlated efficiency error.
6. **The companion analytes are non-monotonic.** Benzaldehyde, dimethyl disulphide and dimethyl
   trisulphide all peak at 85 °C/10 min and **fall** at 90 °C/10 min. Whatever governs the top of the
   range is not one Arrhenius channel.

## 3. Fig. 4 — the treatment set and the significance letters, transcribed

Caption verbatim: "Concentration of hexanal in 3% reconstituted whey protein isolate heated at
different temperatures and times. Mean values of independent triplicates are presented and error bars
represent standard deviations; different letters denote statistical difference (p < 0.05)."
y-axis "Concentration of hexanal (µg L⁻¹)", gridlines 0 to 4.5 in steps of 0.5. No unheated bar.

| bar | T (°C) | t (min) | hexanal | letter |
|---|---:|---:|---|---|
| 70/5 | 70 | 5 | **FIGURE ONLY** | bcd |
| 75/5 | 75 | 5 | **FIGURE ONLY** | cd |
| 80/5 | 80 | 5 | **FIGURE ONLY** | cd |
| 85/5 | 85 | 5 | **FIGURE ONLY** | ab |
| 50/10 | 50 | 10 | **FIGURE ONLY** | d |
| 60/10 | 60 | 10 | **FIGURE ONLY** | cd |
| 70/10 | 70 | 10 | **FIGURE ONLY** | cd |
| 75/10 | 75 | 10 | **FIGURE ONLY** | cd |
| 80/10 | 80 | 10 | **FIGURE ONLY** | bc |
| 85/10 | 85 | 10 | **FIGURE ONLY** | a |
| 90/10 | 90 | 10 | **FIGURE ONLY** | a |

Verbatim: "The concentration of hexanal in the heated samples significantly increased with increased
heat load (Fig. 4), and a significantly higher concentration was observed for samples heated at 85 or
90 °C for 10 min compared with 80 °C min and below. A similar pattern was found for other lipid
oxidation products such as heptanal and nonanal (data not shown)."

**Heptanal and nonanal: "data not shown" — no values, no figure.** Pentanal, octanal, 2-pentylfuran
and hexanol are not mentioned at all (grep-verified); they may sit in the missing Table S1.

**Discrepancy flagged, not resolved.** §2.3 says the 5-min series ran "at 80, 85 or 90 °C for 5 min",
but Figs. 2–4 each plot four 5-min bars — 70/5, 75/5, 80/5, 85/5 — and **no 90/5**, while the sensory
Table 2 compares a "90 °C, 5 min" sample that appears in no figure. Methods text, figures and Table 2
disagree about which 5-min cells were run. Not resolvable from the PDF.

## 4. The Strecker and sulfur analytes — the part that does bear on the trunk

- **2- and 3-methylbutanal (Fig. 2A, figure only): flat.** Verbatim: "No significant differences were
  observed for the concentrations of 2- and 3-methylbutanal in WPI samples heated in the range
  between 50 and 90 °C for 5 or 10 min". **All eleven bars carry the letter "a".** y-axis 0–0.7 µg/L.
- **Benzaldehyde (Fig. 2B, figure only)**, y-axis 0–6 µg/L, letters in bar order (70/5, 75/5, 80/5,
  85/5, 50/10, 60/10, 70/10, 75/10, 80/10, 85/10, 90/10): bc, b, bc, bc, d, d, d, cd, bcd, **a**, bc.
- **Dimethyl disulphide (Fig. 3A, figure only)**, y-axis 0–1.8 µg/L, same bar order: c, c, c, bc, bc,
  c, bc, c, bc, **a**, b.
- **Dimethyl trisulphide (Fig. 3B, figure only)**, y-axis 0–0.35 µg/L: b, nd, b, b, b, b, b, nd, b,
  **a**, a ("nd, not detected" at 75/5 and 75/10).

**The only volatile numbers printed anywhere in this paper**, both in the Discussion and both without
a standard deviation:

> "the concentration of DMTS significantly increased from **0.012 µg L⁻¹** in the sample heated at
> **50 °C for 10 min** to **0.20 µg L⁻¹** in the sample heated at **90 °C for 10 min**"

Flagged: that in-text 0.20 µg/L for 90/10 sits **above** what Fig. 3B draws for that bar (85/10 is
drawn taller). Text and figure are not obviously consistent; both are recorded as printed and neither
is reconciled here. Two points, no dispersion, and an internal inconsistency do not make a barrier.

**No hydroperoxide measurement of any kind** — no peroxide value, no conjugated dienes, no TBARS, no
lipid extraction (grep-verified). The paper constrains the aldehyde output with the pool feeding it
unquantified and no measurable substrate at all.

## 5. The direction the repository can bank, without reading a bar height

Across **50 → 90 °C, a 40 K span**, in a wet protein system at pH 7, hexanal moved **by less than one
significance grouping across eight of the eleven treatments** (50/10 "d" through 80/10 "bc"), with a
detectable rise only at 85–90 °C for 10 min; and doubling the hold from 5 to 10 min at fixed
temperature moved nothing clearly (85/5 "ab" vs 85/10 "a"; 80/5 "cd" vs 80/10 "bc"). That is the
fingerprint of a **weak apparent temperature dependence**, directionally consistent with the
**61–65 kJ/mol moist-system figure rather than the 114–122 kJ/mol bulk-oil figure.**

Given §2.2 (no baseline) and §2.3 (no substrate) this is **suggestive, not evidential, and must not be
converted into a number.** It is recorded as a direction only.

## 6. Kinetics, and what would rescue this paper

**No kinetics**, despite the title. No rate constant, no activation energy, no half-life, no Arrhenius
or Eyring treatment (grep-verified). The quantitative apparatus is a two-way ANOVA with batch, time
and temperature as factors, least-square means at p < 0.05, and sensR for the triangle tests. The
"temperature-dependency" of the title is a categorical threshold claim — "70 °C is the critical
temperature" — not a parameterisation.

**Precisely what is needed:** the numeric hexanal means ± SD for the eleven (T, t) cells plus an
unheated control.

**CORRECTED 2026-09-11, same day, before anyone acted on it.** An earlier draft of this dossier said
those numbers "should be in Supplementary Table S1". **That was an inference and it is very probably
wrong.** The paper says what each supplementary item contains, and none of them is a concentration
table: *"A total of 48 volatile compounds were identified in the heated WPI samples using DHS GC-MS
(Supplementary material Table S1)"* — an **identification** list; *"Standard curves were prepared for
selected volatile compounds (Supplementary material Table S2)"* — **calibration curves**; and
Fig. S1 is the **heating-and-cooling profile**. Nothing in the article states that any supplementary
item carries the per-treatment hexanal values.

So the retrieval is worth doing but is **not sufficient on its own**: S2 and Fig. S1 are both needed
for a fit and neither is on disk, while the eleven numbers behind Fig. 4 most likely exist only in the
authors' own records. **The reliable route is a request to the corresponding author**
(Marianne N. Lund, mnl@food.ku.dk), asking for the numeric means and standard deviations behind
Fig. 4 together with an unheated control. Anyone pursuing this should fetch the supplementary file
too, but should not expect it to close the gap by itself. Even with them, reasons 2–6 above mean any
resulting barrier must carry a wide band and be treated as a weak-dependence data point, not as an
equal-weight competitor to the two measurements already on disk. This is entered in `EXPERIMENTS.md`
as a **retrieval** item, not a measurement item — the one item on that list someone could close by
downloading a file.
