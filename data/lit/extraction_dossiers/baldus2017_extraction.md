# Baldus et al. 2017 — EXTRACTION (trace transition metals from reagent impurities drive thiol-mediated oxidation in a nominally metal-free buffer — REFUSED as a rate, ADOPTED as a named confound on four hold-out rows)

**Source on disk:** `data/articles/baldus2017.pdf` (ACS "Just Accepted" unedited manuscript, 42 pp.;
downloaded 2026-09-11 at this repository's request). Read 2026-09-11 via `pdftotext -layout` plus
rendered figure pages. **Supporting Information Table S1 (the full response-surface dataset) is a
separate file on the ACS site and is not on disk.** Wave B45.

| field | value |
|---|---|
| Title | "Effect of L-Cysteine and Transition Metal Ions on Dimethyl Sulfide Oxidation" |
| Authors | M. Baldus, R. Klie, X. De, F.-J. Methner (TU Berlin, Chair of Brewing Science) |
| Venue | J. Agric. Food Chem., Just Accepted, web 18 Feb 2017 |
| DOI | 10.1021/acs.jafc.6b05472 |
| Medium | **0.1 M sodium acetate / acetic acid buffer, pH 5.5**, ultrapure (Milli-Q) water |
| Atmosphere | purged with synthetic air (20.5 ± 0.5 % O₂) to **8 mg/L O₂**, measured optically; 62 g in 50 mL Duran bottles, headspace minimised, PTFE-lined caps, foil-wrapped |
| Temperature | **95 °C**, 0–180 min (the authors note it took 20 min to reach 95 °C) |
| Quantification | **every value is an absolute µM concentration** — no peak areas anywhere |

## 1. Why this paper was fetched, and what it is not

Fetched for a **thiol-versus-time series with and without a chelator** — the measurement B38 named as
the only thing that can move the thiol sink. **It is not that**, for four reasons:

1. **The thiol is L-cysteine only.** No 2-methyl-3-furanthiol, no 2-furfurylthiol, no dimethyl
   disulfide. The model's two target thiols are absent from the paper.
2. **The only thiol-versus-time data is figure only** (Fig. 6, rendered and checked: axis ticks, no
   data labels). The only numbers in text for that run are the t = 0 value **220 µM**, the H₂O₂
   maximum **~71 µM**, and the qualitative "almost completely degraded in 5 minutes".
3. **That one quasi-rate statement is not a chelator comparison** — it is measured in a system that
   *already contains* 18 µM Cu(II)EDTA, over a **40 → 60 °C ramp**, not an isothermal hold.
4. **The one clean with/without-chelator pair does not measure the thiol.** Table 3 reports the
   thioether (DMS) and its sulfoxide only; cysteine was not assayed in those arms.

The authors say so themselves, verbatim: **"there are no rate constants available, either for DMSO
reduction by Cys as well as for the Cys-Cu(I)-Cys formation and autoxidation."** They present the
work as pre-kinetic groundwork.

## 2. What it *does* establish, in absolute concentrations — and why it matters here

This is the finding wave B45 adopts. Three measurements, all in µM against calibration curves:

**(a) The buffer was clean and the reagent was not.** Verbatim, by ICP-OES with a stated 0.08 µM
detection limit: *"In the buffer solution, all tested transition metal ions were below the detection
limit of 0.08 µM. However, the addition of **300 µM Cys led to an increase of Cu to 0.26 µM**,
whereas all other metal ions were still below the detection limit."* The copper arrived **as an
impurity in the cysteine**, at highest available purity grade, into ultrapure Milli-Q water.

**(b) That trace was enough to drive significant oxidation, and a chelator abolished it.** Table 3 —
16.95 µM initial DMS, 0.1 M acetate pH 5.5, air-saturated, **95 °C for 180 min**, n = 3, letters from
Tukey–Kramer HSD:

| treatment | DMS (µM) | DMSO (µM) |
|---|---:|---:|
| no addition | 16.73 ᵃ ± 0.11 | n.d. |
| **13 µM cysteine** | **15.32 ᵇ ± 0.15** | **0.72 ± 0.02** |
| **13 µM cysteine + 40 µM EDTA** | **16.80 ᵃ ± 0.07** | **n.d.** |

The EDTA arm is statistically indistinguishable from the no-cysteine control. The authors' reading,
verbatim: *"This prooxidative effect was eliminated by molar excess of EDTA. This indicates, that
certain amounts of transition metal ions were already present in the model solution, which were then
chelatively inactivated by the molar excess of EDTA."*

**The chelation logic is stated precisely and is not "EDTA suppresses oxidation".** Verbatim: *"DMS
oxidation could only be eliminated when EDTA was added in **relative molar excess to potentially
abundant transition metal ions and Cys**… The close proximity of the stability constants indicates
that the substance, which is present in relative molar excess to Cu, dominates its complexation."*
At 1:1 with copper and **sub-stoichiometric to the thiol, EDTA made things worse, not better** —
DMS loss went from ~40 % (Cu alone + cysteine) to **~72 %** (Cu(II)EDTA + cysteine), because
chelation merely delays the cysteine–copper exchange to a higher temperature. Any experiment that
adds "a chelator" without checking it against the thiol concentration can get the opposite of the
intended control.

**(c) Cysteine itself vanishes fast in the presence of that copper.** Verbatim: **"Cys was almost
completely degraded in 5 minutes during heating from 40–60 °C."** Starting point 220 µM; H₂O₂ peaked
at ~71 µM on reaching 95 °C, a cysteine-to-peroxide ratio of 3.22:1.

## 3. The other printed numbers, for the record

DMS consumption over 180 min at 95 °C, as percentages in the running text (the underlying bars in
Figs. 1 and 3 carry no data labels and are figure only):

| system | DMS consumed | DMSO formed |
|---|---:|---:|
| Fe(II) alone, or Cu(II) alone, no cysteine | "below 2 %" | — |
| Fe(II) + 250 µM cysteine | ~13 % (the abstract says ~12 %; both printed, neither reconciled) | ~8.3 % |
| Cu(II) + 250 µM cysteine | ~40 % | ~20 % |
| Fe(II)EDTA + cysteine | ~7.7 % | ~7 % |
| **Cu(II)EDTA + cysteine** | **~72 %** | **44 %** |

Response-surface models over 180 min at 95 °C (endpoint regressions, **not rate laws**):
`DMS = 14.99295 − 0.043752·Cys` and `DMSO = 4.21 − 0.57·Cu(II)EDTA + 2.60·Cys − 0.66·Cu(II)EDTA·Cys`.
ANOVA: cysteine dominates both responses (F = 197.1 and 172.2, p < 0.0001); Cu(II)EDTA is
insignificant for DMS (p = 0.39) and significant for DMSO (p = 0.019).

DMSO reduction by cysteine, deoxygenated (< 0.08 mg O₂/L), 1 mM EDTA, 250 µM cysteine, **24 h at
95 °C**: DMSO fell from **13 ± 0.02 to 11.8 ± 0.54 µM** with **1.47 ± 0.11 µM** DMS formed — "the
overall DMSO reduction was below 10 %".

Literature rate constants quoted (none measured here): DMS + H₂O₂ 1.4–8.1 × 10⁻² M⁻¹s⁻¹ at 22–25 °C;
DMSO + H₂O₂ 0.5–4.5 × 10⁻⁵; DMS + •OH 1.9 × 10¹⁰; stability constants cysteine–Cu K = 19.2,
cysteine–Fe K = 6.2. **All are citations and none is adopted here.**

## 4. Where this lands in this repository

The four **Yiltirak ribose + cysteine hold-out bundles** at 100–130 °C are the **only four bundles in
the corpus carrying `water_source: tap`**. They are hold-out rows, and they are rows the thiol sink
fails on.

Baldus's argument transfers a fortiori and the direction is unambiguous: if **0.26 µM** of copper,
arriving as a reagent impurity into **ultrapure water**, sufficed to cause statistically resolved
oxidation at 95 °C, then **tap water is not a controlled medium for a thiol experiment.** The size of
the effect in those pots is unknown and this dossier does not estimate it.

**What B45 does with it:** a clause on those bundles' buffer note naming the paper and the confound.
**What B45 does not do:** edit the bundles, widen their tolerances, or add a metal-catalysed channel
to the model. There is no rate constant to add — the authors state there is none — and inventing one
to close a gap would be exactly the failure this repository exists to avoid.

## 5. The consequence for the thiol sink, stated plainly

B38 found the thiol sink unreachable by refit: the barrier and **both** dimerisation rates already sit
on their ceilings. B45's probe P2 confirms the shape of the gap — the model retains **99.49 %** of
charged cysteine after 5 minutes at 95 °C, and **83 %** after three hours, where Baldus's system has
lost essentially all of it within 5 minutes of a 40–60 °C ramp.

So the missing quantity is **not a larger constant; it is a catalytic channel.** A metal-catalysed
thiol autoxidation depends on a catalyst concentration the model does not carry, and forcing
`k_cys_thermal` upward to match would fit a catalytic rate into a thermal barrier and then be wrong
at every other temperature. This is recorded as a **named missing mechanism**, and it sharpens B38's
experiment request: the thiol-against-time run must be **paired with and without a chelator in molar
excess over the thiol**, in water of stated provenance, with the disulfide quantified in the same run.
