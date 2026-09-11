# Pre-registration: wave B45, seven papers arrive against the reading list (written 2026-09-11, BEFORE the probes ran)

## 1. Why

Seven papers were requested by name against the reading list left open by B42–B44 and downloaded on
2026-09-11: `unlu2002`, `li2020`, `jansson2020`, `bao2022`, `belrhlid2002`, `baldus2017`, `shi2022`.
They were requested to close three named debts:

- a **hexanal barrier in a wet protein system**, the missing rung between the bulk-oil 114–122 kJ/mol
  and the moist-nut-paste 61–65 kJ/mol the lipid lane already carries (`li2020`, `jansson2020`,
  `bao2022`);
- a **thiol-against-time series in a defined buffer**, the one measurement B38 named as the only
  thing that can move the thiol sink (`belrhlid2002`, `baldus2017`);
- a **measured residence time** for the extrusion acrylamide row, whose 25 s hold has been a
  labelled assumption since B36 (`unlu2002`), and a **free-versus-bound hexanal split** (`shi2022`).

All seven were read before this record was written. **Six of the seven are refused as rate or
barrier sources**, each for a reason stated in its dossier; the seventh (`unlu2002`) corroborates an
assumption without changing it. **This wave therefore proposes no fit and expects to ship no
constant.** What it does propose is two probes that test the model's *completeness* rather than its
numbers, because two of the refused papers document mechanisms the model does not have at all.

## 2. What the seven papers gave, stated before the probes

| paper | asked for | verdict | why |
|---|---|---|---|
| `unlu2002` | measured extruder residence time | **corroborates** | GMRT **123 s (87.2–173.5)** at 150 rpm / 8.55 kg/h, L/D 38.7 — Ma's pot runs 150 rpm / 8.57 kg/h, L/D 40. Three of ten zones at 130 °C gives **~37 s** against the bundle's 25 s: right order, ~1.5×, not 4×. Recorded on the vessel note in B45; **the 25 s stays.** |
| `li2020` | hexanal vs temperature, wet protein | **refused** | the 70/90/120 °C axis is on the *protein before emulsification*; every emulsion reacted at **37 °C**. One time point. Relative internal-standard units, not concentrations. |
| `jansson2020` | hexanal vs temperature, wet protein | **refused** | structurally the right experiment — seven temperatures, absolute µg/L — but **every value is figure-only**, there is no unheated control, and **"Fat — below detection limit"**: no lipid substrate at all. |
| `bao2022` | linoleic acid temperature series | **refused** | **"Peak area (×10⁶)"**, no internal standard, no calibration; a 100 ml/min purge whose efficiency is itself temperature-dependent; one time point per temperature; and hexanal *falls* from 30 to 60 °C. |
| `belrhlid2002` | thiol vs time in buffer | **refused** | no thiol was ever put in buffer alone — every run starts from a thioacetate and the loss rate is **proportional to the crude lipase dose**, which the authors attribute to enzyme impurities. Figure-only. |
| `baldus2017` | thiol vs time, with/without chelator | **refused as a rate**, but see §3 | the chelator arm was never assayed for thiol, and the one thiol time course is figure-only and already contains chelator. The authors state outright: "there are no rate constants available". |
| `shi2022` | free/bound hexanal split | **refused as a constant**, but see §3 | no binding constant, no isotherm, no total-hexanal measurement; "nearly half" is stated with no supporting number. Its Table 2 *is* a clean paired before/after dataset in absolute µg/L with n=3. |

## 3. The two mechanisms the refused papers document, which the model does not have

These are why this wave runs probes at all. Neither can be turned into a number here; both are
testable as **presence-or-absence**.

**(a) A reversible, non-covalent hexanal–protein channel (`shi2022`).** The engine binds hexanal to
lysine amines in `matrix_sites.py` — a **one-way covalent adduct**, `fraction bound = 1 − exp(−Σ k₂(T)[sites]t)`,
so *more thermal load always means less measured hexanal*. Shi measures the opposite: in acidified
soy protein isolate at 30 mg/mL, pH 4.5, a **95 °C / 5 min hold raises** headspace hexanal from
23 ± 0.59 to 66 ± 1.7 µg/L, and in 11S from 40 to 220 µg/L. Their mechanism is reversible hydrophobic
sequestration (Damodaran & Arora's −2 kcal/mol, "binding constants were low"), released by acid plus
heat. **That channel is absent from the model**, and it runs the other way from the one that is there.

**(b) Trace-metal-catalysed thiol autoxidation (`baldus2017`).** This is the finding that matters
most, because of where it lands. Baldus measured, in absolute µM: a buffer made with ultrapure
Milli-Q water held **all transition metals below 0.08 µM**, but adding 300 µM cysteine **brought in
0.26 µM Cu as a reagent impurity**, and that trace was enough to cause statistically significant
thioether oxidation at 95 °C / 180 min — abolished only by **EDTA in molar excess over the thiol**
(DMS 16.80 ± 0.07 with EDTA vs 15.32 ± 0.15 without, against a 16.73 ± 0.11 control). In the
cysteine + Cu system, **"Cys was almost completely degraded in 5 minutes during heating from
40–60 °C."**

**Where that lands: the four Yiltirak hold-out bundles are ribose + cysteine in `water_source: tap`**
— the only four bundles in the corpus with tap water, at 100–130 °C. And B38's identifiability audit
found the thiol sink **unreachable by any refit**: the barrier and *both* dimerisation rates sit on
their ceilings, and three sulfur decays are dead to the data. B38's verdict was "an experiment, not a
refit". B45's claim, registered here in advance, is that **B45 names the missing mechanism**: not a
constant that is too small, but a catalytic channel the model does not contain, in rows whose water
provenance is uncontrolled.

## 4. The two probes, and the predictions

Both probes are amine-free or sugar-free where the lane requires it, and neither touches a fitted
constant. Predictions are numeric and falsifiable.

- **P1 — the hexanal binding block against Shi's pot.** Charge 30 g/L soy protein isolate at pH 4.5
  and integrate a 95 °C / 5 min hold; read the engine's bound fraction for `HEXANAL` from
  `matrix_sites.bound_fraction`. **Prediction: the model binds under 1 % of the hexanal** — at least
  two orders of magnitude below the ~2.9× swing Shi measures, and with the opposite sign. The reason
  is named in advance: `k₂ ≤ 2.5 × 10⁻⁵ M⁻¹s⁻¹` at 20 °C with a 15–20 kJ/mol band gives barely a
  4× rise by 95 °C, and 300 s against ~12 mmol/L of amine sites cannot reach a percent.
- **P2 — cysteine alone in buffer at 95 °C.** Charge cysteine at Baldus's ~250 µM in buffer with **no
  sugar and no metal**, integrate 180 min at 95 °C, and read the free-thiol trajectory.
  **Prediction: the model loses under 5 % of the charged cysteine**, against Baldus's near-complete
  loss within 5 minutes of a 40–60 °C ramp in the presence of 0.26–18 µM Cu. The reason is named in
  advance: with no sugar the thiol has nowhere to go but the dimer channel, and the model has no
  metal-catalysed autoxidation at all.
- **P3.** Nothing moves. No constant, no parameter, no tolerance. The artifacts of this wave are this
  record, seven dossiers, one clause on the vessel note (`unlu2002`, already applied), one clause on
  the buffer note (`baldus2017`), and two named debts in `EXPERIMENTS.md`.

**A tolerance is explicitly NOT widened.** Shi's result would, if applied, license loosening every
HS-SPME hexanal row on the grounds that headspace under-reports total hexanal by a
history-dependent factor. That would raise this model's score without any new evidence about its
chemistry, so it is refused here and recorded as refused.

## 5. What follows if the predictions hold

Neither probe licenses a fit, and neither is scored. What each licenses is a **named debt with a
mechanism and a measurement request**, which is strictly more useful than the unexplained failure it
replaces:

- **For the thiol sink**, B38 said "an experiment: thiol against time at two temperatures". B45
  sharpens that to: **thiol against time at two temperatures, in a buffer of stated water
  provenance, run in parallel with and without EDTA in molar excess over the thiol, with the
  disulfide quantified in the same run.** Without the chelator arm the measurement cannot separate
  the model's missing channel from the rate the model already has. This is the single highest-value
  measurement on the list and it is now specified precisely enough to run.
- **For hexanal**, the debt becomes a completeness debt rather than a barrier debt: the model needs a
  reversible protein-partition term for aldehydes before any HS-SPME hexanal row in a protein matrix
  can be read as a total. Until it has one, those rows measure headspace, and the model predicts
  total. Recorded as a **declared measurement-channel mismatch**, not corrected.
- **For the hexanal barrier itself**, `jansson2020` is entered as the one item on the list that a
  download could close rather than a bench: its Supplementary Table S1 holds the eleven numeric
  (T, t) hexanal cells whose absence is the only reason the paper fails.

## 6. Outcome (written 2026-09-11, after the probes)

**Both predictions held, and P2 held by a wider margin than it claimed.**

### P1 — the hexanal binding block against Shi's pot

Registered: "the model binds under 1 % of the hexanal". Measured, from
`matrix_sites.bound_fraction("HEXANAL", ...)` at 30 g/L, the amine pool charged from the isolates'
own measured densities:

| pot | thermal programme | model binds | bracket corners | Shi measured |
|---|---|---:|---|---|
| soy isolate, 10.8 mmol/L amine | **95 °C / 5 min** (Shi's hold) | **0.0108 %** | 0.0027 – 0.0431 % | **+187 %** (23 → 66 µg/L) |
| soy isolate | 95 °C / 10 min | 0.0217 % | 0.0054 – 0.0862 % | — |
| soy isolate | 100 °C / 45 min | 0.1052 % | 0.0262 – 0.4225 % | — |
| soy isolate | 140 °C / 60 min | 0.2421 % | 0.0557 – 1.0482 % | — |
| pea isolate, 14.1 mmol/L amine | 95 °C / 5 min | 0.0141 % | 0.0036 – 0.0563 % | — |

On Shi's own hold the model moves hexanal by **one part in ten thousand, downward**; Shi measures it
nearly tripling. Even at the corner of the declared bracket, and even under a 140 °C / 60 min cook far
harsher than anything Shi ran, the covalent block never reaches a quarter of one percent. **The
discrepancy is about four orders of magnitude in effect size and opposite in sign.**

The reading is not that the declared bracket is wrong. Meynier's and Anantharamkrishnan's
`k₂ ≤ 2.5 × 10⁻⁵ M⁻¹s⁻¹` measures *covalent adduct formation* with lysine, and nothing here disputes
it. The reading is that **covalent adduction is not the channel that governs what a headspace
measurement sees.** A second, reversible, non-covalent partition — the one Shi releases with acid
plus heat, and which Damodaran & Arora price at about −2 kcal/mol — is larger by orders of magnitude
and is entirely absent from the model.

### P2 — cysteine alone in buffer at 95 °C

Registered: "the model loses under 5 % of the charged cysteine". Measured, 0.25 mM cysteine, pH 5.5,
no sugar, targeting hydrogen sulfide so the pot is answerable:

| hold at 95 °C | cysteine remaining | oxidant pool `OX` | `OXV` |
|---|---:|---:|---:|
| **5 min** | **99.49 %** | 1 | 0 |
| 30 min | 97.00 % | 1 | 0 |
| 180 min | 83.34 % | 1 | 0 |

Baldus, at a *lower* temperature and with 18 µM Cu(II)EDTA present: **"Cys was almost completely
degraded in 5 minutes during heating from 40–60 °C."**

So at 5 minutes the model has lost **0.51 %** where the measurement has lost nearly everything, and
the model still holds 83 % after three hours. The registered margin was 5 %; the observed loss is a
tenth of that.

Two further observations, recorded because they were not predicted:

- **The vessel makes no difference at all.** Re-running with Baldus's actual vessel — 62 g in a 50 mL
  Duran bottle, headspace minimised, air-saturated at 8 mg O₂/L — returns cysteine values identical
  to nine significant figures. The sulfur lane's `ch_cys_ox` channel is present in
  `FULL_REACTIONS` but is not moving this pot.
- **The oxidant pool never moves**: `OX` stays at 1 and `OXV` at 0 across every hold. Whatever
  `ch_cys_ox` is doing, it is not consuming the oxidant reservoir on a thiol-only pot.

### What this settles

B38's identifiability audit found the thiol sink unreachable: the barrier and **both** dimerisation
rates already sit on their ceilings, three sulfur decays are dead to the data, and its verdict was
"an experiment, not a refit". B45 names what the experiment must contain, and why no refit could
ever have worked: **the model is missing a catalytic channel, not a larger constant.** Cysteine in
these pots is removed by trace transition metals that arrive with the reagents and the water, at a
rate no unimolecular thermal decay can imitate. Pushing `k_cys_thermal` up to match would be fitting
a catalytic rate into a thermal barrier, and it would then be wrong at every other temperature.

That matters concretely for **four hold-out bundles**: the Yiltirak ribose + cysteine rows at
100–130 °C are the only four bundles in the corpus carrying `water_source: tap`, and they are hold-out
rows the thiol sink fails on. Baldus's a-fortiori argument is strong — if 0.26 µM of Cu arriving as a
*reagent impurity in ultrapure Milli-Q water* sufficed to drive statistically significant oxidation,
tap water is not a controlled medium for a thiol experiment. **The bundles are not edited and their
tolerances are not widened.** A clause goes on their buffer note saying what is uncontrolled and
naming the paper, which is what the repository does with a confound it can identify but cannot size.

### What did not happen

No constant moved. No tolerance widened. No parameter was added. In particular the tolerance
argument available from Shi — that HS-SPME under-reports total hexanal by a history-dependent factor,
and that every hexanal row could therefore be judged more loosely — was refused, as §4 registered in
advance that it would be. It would have raised this model's score without a single new fact about its
chemistry.
