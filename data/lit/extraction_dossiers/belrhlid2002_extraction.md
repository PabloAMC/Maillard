# Bel Rhlid et al. 2002 — EXTRACTION (2-methyl-3-furanthiol and 2-furfurylthiol stability in phosphate buffer — REFUSED as a thiol removal rate, enzyme-confounded and figure-only)

**Source on disk:** `data/articles/belrhlid2002.pdf` (downloaded 2026-09-11 at this repository's
request). Read 2026-09-11 via `pdftotext -layout` plus rendered figure pages. Wave B45.

| field | value |
|---|---|
| Title | "Lipase-Assisted Generation of 2-Methyl-3-furanthiol and 2-Furfurylthiol from Thioacetates" |
| Authors | R. Bel Rhlid, W. Matthey-Doret, I. Blank, L. B. Fay, M. A. Juillerat (Nestlé Research Center, Lausanne) |
| Venue | J. Agric. Food Chem. **50** (14) (2002) 4087–4090 |
| DOI | 10.1021/jf0202335 |
| System | thioacetate 0.064 mmol in 10 mL (**6.4 mM**) of water or **0.2 M phosphate buffer**; *Candida rugosa* lipase 6.5–262 units |
| pH | 5.8, 6.0, 7.0, 8.0 | 
| Temperatures | **4, 23, 37 °C** |
| Times | 1 min to 72 h |
| Atmosphere | **NOT STATED anywhere** — only "gentle magnetic stirring". For a paper about air-oxidation of thiols to disulfides this is a material omission. |

## 1. Why this paper was fetched

B38's identifiability audit named one measurement as the only thing that could move the thiol sink:
**thiol against time in a defined buffer, at two temperatures.** This paper measures both of the
model's target thiols — 2-methyl-3-furanthiol and 2-furfurylthiol — against time, in 0.2 M phosphate
buffer, at three temperatures. On the face of it, exactly the missing experiment.

## 2. Why it is refused — four reasons, any one sufficient

1. **No thiol was ever put in buffer alone.** Every experiment starts from a **thioacetate** and
   generates the thiol enzymatically in situ. What the curves show is the *net* of simultaneous
   formation and loss. There is no experiment in this paper in which a thiol is dissolved in buffer
   and watched.
2. **The loss rate is proportional to the enzyme dose, and the authors say it is an impurity
   artefact.** Verbatim: "the degradation rate of aroma compound 5 was proportional to the quantity
   of enzyme"; "The degradation of this volatile molecule was much higher when the amount of
   biocatalyst was increased"; and the attribution, verbatim: "This phenomenon could be explained by
   the presence of **side activities due to the presence of impurities in the commercial crude
   enzyme preparation**." Any rate read off these curves is a rate in buffer **plus a crude
   commercial enzyme of unknown side-activity**, not a rate in a defined buffer.
3. **All kinetic data are figure only.** There is no data table (Table 1 is a literature odour-
   threshold table). Figures 2–6 plot "Yield (%)" against a **categorical, non-linear** time axis.
   No point is printed and none has been estimated here.
4. **The disulfide was never plotted in the same run as any thiol curve.** So the mass balance that
   would attribute thiol loss to dimerisation rather than to volatilisation, extraction loss or
   enzyme side-activity cannot be closed.

## 3. The endpoints that are printed

**2-Furfurylthiol** (0.064 mmol, 65 units, pH 5.8, 0.2 M phosphate):
- "a yield of 80% after 1 h of reaction time" (abstract, room temperature, pH 5.8)
- "A maximum yield of 74% was obtained at pH 5.8, whereas the maximum yield at pH 8.0 was only 50%."
- **"A maximum yield of ∼70% of 2-furfurylthiol was obtained after 24 h at 4 °C and after 1 h and
  45 min at 23 and 37 °C, respectively."** — a *time-to-peak*, not a decay rate.
- "Moreover, 2-furfurylthiol was much more stable at 4 °C as compared to 23 and 37 °C, leading mainly
  to the dimer 3."

**2-Methyl-3-furanthiol** (0.064 mmol, 65 units, pH 5.8, 23 °C, 0.2 M phosphate):
- "a yield of 88% was obtained after 15 min of reaction time when 65 units of enzyme was used…
  whereas the yields were only 50 and 10% when 26 and 6.5 units of enzyme were used, respectively."
- **"odorant 5 is easily consumed as only 50% was left after an incubation period of 1 h"** — the
  single most quantitative stability statement in the paper, and see §2.1 and §2.2 for why it is not
  a removal rate.
- "thiol 5 was completely transformed into the corresponding disulfide 6 after 2 h of reaction time"
  (n-hexane / n-pentane arm).

**Bis(2-furfuryl) disulfide:** one number in the whole paper — "up to ∼10% of compound 3 was obtained
at all pH values studied after 24 h of reaction time". **Bis(2-methyl-3-furyl) disulfide: no
quantitative value anywhere.**

**Medium comparison, all asserted without a number:** "the stability of the generated 2-methyl-3-
furanthiol and 2-furfurylthiol was better in n-hexane, n-pentane, and the water/propylene glycol
mixture as compared to that in water or phosphate buffer." No half-life, no percentage lost, no rate
is given for any medium.

**Figure 4** is the only multi-temperature thiol-vs-time run — 2-furfurylthiol at 4, 23 and 37 °C,
pH 5.8, 0.2 M phosphate, 65 units. It has the right shape and is **figure only**, is the wrong thiol
for the sink, plots no disulfide, and carries the enzyme confound in full.

## 4. Kinetics

**None.** No rate constant, no half-life, no reaction order, no Arrhenius treatment, no activation
energy. The phrase "Kinetic studies were carried out" appears once, in the propylene-glycol section,
and no kinetic parameter is reported from it.

## 5. What is still missing, stated precisely

The measurement B38 asked for, restated so that it cannot be mistaken for this paper again:

> a control incubation of **pure 2-methyl-3-furanthiol (or 2-furfurylthiol) in a stated buffer at a
> stated pH, with no enzyme and no thioacetate present**, sampled against time at **≥2 temperatures**
> under a **stated atmosphere**, with the free thiol **and its disulfide** both quantified in the
> same run against a calibrated internal standard, and the numeric values printed.

`belrhlid2002` supplies none of those six things. Wave B45 adds a seventh requirement from
`baldus2017_extraction.md`: the run must be **paired with and without a chelator in molar excess over
the thiol**, or it cannot separate the model's missing metal-catalysed channel from the thermal rate
the model already has.

## 6. Verdict

**REFUSED as a thiol removal rate.** Kept as qualitative corroboration of an ordering the sulfur lane
already assumes: 2-methyl-3-furanthiol is markedly less stable than 2-furfurylthiol in aqueous
phosphate buffer at room temperature, both are more stable cold, and both are more stable in
non-aqueous media. All of it carries the crude-lipase confound the authors themselves identify.
