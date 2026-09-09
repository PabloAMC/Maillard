# Desclaux et al. 2006 — EXTRACTION (xylose 6 mmol + glycine 6 mmol in 25 g of 0.5 M phosphate pH 6, sealed stirred steel tubes, 80 / 100 / 120 C, sampled at 30 / 60 / 90 / 120 min; a multiresponse ODE network of deoxyosones, dicarbonyls, hydroxycarbonyls and furanoids fitted in Athena Visual Studio; four conference pages)
### Two first-order rate constants are printed — ARP -> 1-deoxyosone 2.79e-2 /min and ARP -> 3-deoxyosone 5.50e-4 /min at 100 C, pH 6 — and nothing else: no barrier, no dicarbonyl number, no table; the four time courses shown are figure-only, and glyoxal / methylglyoxal / diacetyl / pentanedione appear only as boxes in the network diagram.

**Source on disk:** `data/articles/desclaux2006.pdf` (4 pp., pp. 367-370 of the proceedings; a
scanned bitmap with an OCR text layer, not a born-digital PDF; owner's download, 2026-09-09). Read
from the OCR text layer (`scratchpad/articles/desclaux2006.txt`, 210 lines) and, because the OCR
garbles superscripts and the figures, from page renders of pp. 369-370 (Figure 1 network, Figure 2
fits, the two rate constants). No table exists in the paper. No supplement; the "detailed
mechanisms to be published separately" and the "subsequent paper [that] will report the kinetic
interpretation" are the Reading PhD thesis (Desclaux 2006), NOT on disk (flag 1). Repo status
before this dossier: `parker2013_extraction.md` (row "ref 55" and its section 5) names this paper
as "the only tabulated glyoxal / methylglyoxal / 2,3-butanedione / 2,3-pentanedione time courses at
80/100/120 C in water, pH 6"; `balagiannis2015_extraction.md` row 8 quotes the thesis as covering
xylose/glycine, glucose/glycine and isoleucine/xylose at pH 4 / 6 / 8 with "rate constants and
activation energies for the whole network" and HMF, furfural, furaneol, nor-furaneol, maltol and
2-methylbutanal; `results/validation/kinetic_core_b19_prereg_draft.md` row 33 lists it under "still
to fetch ... (glyoxal, methylglyoxal, butanedione time courses 80 to 120 C)". None of that is in
these four pages (flags 1-3).

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Modelling the formation of Maillard reaction intermediates for the generation of flavour" |
| Authors | Guillaume Desclaux (a), Tahir I. Malik (b), Chris Winkel (c), D. Leo Pyle (a), Donald S. Mottram (a) — (a) University of Reading, School of Food Biosciences; (b) ICI Strategic Technology Group, Wilton Centre, Redcar; (c) Quest Foods, Bussum, NL |
| Venue | W. L. P. Bredie and M. A. Petersen (eds), *Flavour Science: Recent Advances and Trends*, Developments in Food Science 43, Elsevier, 2006, pp. 367-370 (conference paper, 11th Weurman Flavour Research Symposium) |
| DOI | not printed on the pages; the PDF metadata carries PII S0167-4501(06)80087-7 |
| PDF on disk | `desclaux2006.pdf`, 4 scanned pages (CCITT bitmaps + OCR layer) |
| Tables / figures | 0 tables; Figure 1 (network diagram), Figure 2 (three panels of model-vs-data at 100 C, pH 6: xylose; nor-HDF; 3-deoxyosone + hydroxyacetone) |
| Systems in THIS paper | xylose + glycine only, pH 6 only, 80 / 100 / 120 C (results shown for 100 C only) |
| Systems NOT in this paper | glucose + glycine, isoleucine + xylose, pH 4 and pH 8 — thesis material quoted second-hand by Balagiannis 2015 (flag 2) |
| Naming | ARP = Amadori rearrangement product; Done-1 / DONE-1 / 1-done = 1-deoxyosone (1-deoxypentosone from xylose); Done-3 / 3D / 3-done = 3-deoxyosone (3-deoxypentosone); nor-HDF = 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol); Hyd / Hydroxyac = hydroxyacetone; Glycolald = glycolaldehyde; Methylgly = methylglyoxal; Pendione = 2,3-pentanedione; Diacetyl = 2,3-butanedione; Sugardegrad = amine-free sugar degradation; "Maillard Products" = the common sink |
| Companions on disk | `parker2013_extraction.md` (the review that cites this paper; ref 55), `balagiannis2015_extraction.md` (row 8, thesis summary), `hofmann1998b_extraction.md` / `hofmann1998_reconciliation.md` (the derivatisation and the dicarbonyl routes it cites as refs 3, 4, 6), `mottram2002_extraction.md` (ref 7), `zhou2025b_extraction.md` (the other xylose + amino acid Amadori rate on disk) |

## 1. Why it matters

The pyrazine step shipped in `results/validation/kinetic_core_b18_prereg.md` section 6 carries a
supply caveat: its two Strecker constants are measured on fed dicarbonyls, and from a sugar + amine
pot the trunk makes glyoxal only through the amine-free glass entry of
`src/kinetic_core/parameters_dicarbonyl.py` (Kocadagli & Gokmen 2016, 160-200 C), so a glucose +
glycine pot at 95 C in water gives almost no pyrazine (T3 and T4 there fail by three decades). What
that caveat asks for is the formation of the small dicarbonyls from a sugar + amine pot in water at
70-120 C. Parker 2013 pointed at this paper as the tabulated source. It is not one. The four pages
describe the network (Figure 1: ARP -> 1- and 3-deoxyosone -> glyoxal, methylglyoxal,
glycolaldehyde, hydroxyacetone, diacetyl, 2,3-pentanedione, furfural, nor-HDF -> "Maillard
Products") and the fitting procedure, show four fitted time courses at 100 C (xylose, nor-HDF,
3-deoxyosone, hydroxyacetone; Figure 2, figure-only), and print exactly two numbers: the
first-order constants for ARP -> 1-deoxyosone (2.79e-2 /min) and ARP -> 3-deoxyosone (5.50e-4
/min) at 100 C, pH 6. Those two are trunk quantities on the pentose side (the sulfur lane's pentose
Amadori step; `parameters.py` carries Martins' hexose equivalents k_ama_odg 1.57e-2 and k_ama_tdg
1.11e-2 /min at 100 C), and their ratio — 1-deoxyosone favoured fifty-fold over 3-deoxyosone at
pH 6 in phosphate — is a within-study statement the trunk's hexose ratio of 1.4 does not share
(section 4). For the dicarbonyl supply itself the paper gives no number; the wishlist item stands
and the thesis is what would fill it (flag 1).

## 2. Methods as they matter to a model

- **Pot (verbatim core).** "Samples (25 g) were prepared using xylose (6 mmol), glycine (6 mmol)
  and buffer (0.5 M phosphate at pH 6)." Conversions (mine): xylose 6 mmol = 0.901 g, glycine 6
  mmol = 0.450 g; 6 mmol in 25 g of solution is 0.24 mol/kg of each, about 0.24-0.25 mol/L if the
  0.5 M phosphate solution has density 1.03-1.05 kg/L (density not printed). Sugar : amine = 1 : 1.
  Phosphate 0.5 M is a strong catalyst of the Amadori rearrangement and of 2,3-enolisation; the
  constants are phosphate-pH-6 constants (flag 5).
- **Heating.** "Reaction was carried out in stirred sealed steel tubes placed in thermostatically
  controlled heating block. Batch experiments were performed at pH 6 and temperatures of 80, 100,
  and 120 C. Measurements were taken at 30, 60, 90 and 120 min. The time to rise to the final
  temperature was typically 5 min." Tube volume, headspace and pressure not stated; no time-zero
  sample; four time points per temperature.
- **Analytes and derivatisation (verbatim core).** "Dicarbonyls were analysed after derivatisation
  with 1,2-diaminobenzene; hydroxycarbonyls were analysed after derivatisation with ethoxyamine-HCl;
  deoxyosones were derivatised with 1,2-diaminobenzene prior to silylation with BSTFA ... and
  heterocyclic compounds were directly analysed [3]. All analyses were carried out by solvent
  extraction followed by GC or GC-MS analysis." Ref 3 is Hofmann 1999 (Eur Food Res Technol 209:113),
  the quinoxaline / OPD method. No calibration, recovery, internal standard, LOD or replicate
  statement. Which compounds were actually quantified is not listed; the text names, as network
  members, glyoxal, methylglyoxal, diacetyl, pentanedione, hydroxyacetone, glycolaldehyde,
  2-furfural and nor-HDF, plus the two deoxyosones, xylose, glycine and ARP.
- **Concentration unit.** Figure 2's y-axes read "µmol per mmol initial sugar" (xylose panel:
  presumably the same unit, i.e. xylose remaining per mmol charged). No absolute concentration is
  printed anywhere; with 6 mmol xylose in 25 g, 1 µmol per mmol initial sugar = 0.24 µmol per g of
  pot = about 0.24 mmol/L (mine, density assumption as above).
- **Model (verbatim core).** "The model consists of an ODE ... system ... Apart from the first stage
  reaction of sugar with amino acid to ARP, all reactions were assumed to follow first order
  kinetics with Arrhenius temperature dependence." Software Athena Visual Studio (Stewart &
  Associates, Madison); non-linear least squares; parameters first estimated for sugar + amino acid
  -> ARP on a reduced model, then all parameters estimated simultaneously on the full network with
  those as starting values, iterated to convergence. The network (Figure 1, redrawn as edges):
  Xylose -> Sugardegrad; Xylose + Amino Acid -> ARP; ARP -> DONE-3; ARP -> DONE-1; DONE-3 ->
  Furfural; DONE-3 -> Glyoxal; DONE-3 -> Glycolald; DONE-3 -> Methylgly; DONE-1 -> Nor-HDF; DONE-1 ->
  Diacetyl; DONE-1 -> Hydroxyac; DONE-1 -> Pendione; DONE-1 -> Glycolald; DONE-1 -> Methylgly;
  Glyoxal <- (a second arrow from the DONE-1 side); Amino Acid regenerated from ARP (an arrow back to
  Amino Acid); every product -> Maillard Products. Whether "Arrhenius temperature dependence" was
  fitted across 80 / 100 / 120 C in the reported run, or the 100 C set alone, is not stated; only
  100 C results are shown (flag 3).
- **Replicates, uncertainty.** None stated. No confidence interval on the two printed constants.

## 3. Tables re-typed

**The paper has no table.** Every quantitative statement is either the two rate constants in the
running text (p. 370) or a point in Figure 2. Re-typed here is everything the text prints.

### Numbers printed in the running text (p. 370)

| statement as printed | value | unit | conditions |
|---|---|---|---|
| "The model gave rate constant of 2.79x10^-02/min for the formation of 1-deoxyosone from ARP" | 2.79e-2 | /min | xylose + glycine, 0.5 M phosphate pH 6, 100 C, fit to 30-120 min data |
| "and 5.50x10^-04/min for the formation of 3-deoxyosone" | 5.50e-4 | /min | same |
| "meaning that 1-deoxyosone is formed faster at pH 6 (i.e. 1-done is favoured compared to 3-done)" | ratio 50.7 (mine) | — | same |

(The OCR layer renders these as "2.79x10-~" and "5.50xl 0~"; the exponents -02 and -04 and the
"/min" were read from the page image at 130 dpi, where they are unambiguous.)

### Figure 2 (FIGURE-ONLY): "Fit of the model (lines) to the experimental data (points) for Xylose-Glycine pH6 at 100 C showing xylose (sugar), nor-HDF, and 3-deoxyosone (3D) and hydroxyacetone (Hyd)"

Three panels; x-axis "Time (min)", four data points per species at the 30 / 60 / 90 / 120 min
sampling times; y-axis "µmol per mmol initial sugar" on every panel. Qualitative content the text
itself states: "The model fitted well the formation of hydroxyacetone and nor-HDF. 3-deoxyosone was
formed rapidly in the initial reaction and then suffered a net loss over the main course of the
reaction. This was also predicted well but ... the prediction for the loss of xylose could be
improved." No point value is typed here; the four series are figure_only.

### What the network diagram (Figure 1) contains and the text does not quantify

Glyoxal, methylglyoxal, glycolaldehyde, hydroxyacetone, diacetyl, 2,3-pentanedione, furfural,
nor-HDF, the two deoxyosones, ARP, xylose, glycine, an amine-free sugar-degradation sink and a
common "Maillard Products" sink. Of these, only xylose, nor-HDF, 3-deoxyosone and hydroxyacetone are
shown against time, and only the two ARP -> deoxyosone constants are printed. **There is no
glyoxal, methylglyoxal, diacetyl or 2,3-pentanedione number, curve or table in this paper.**

## 4. Numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): nor-HDF -> `norfuraneol`; furfural -> `furfural`;
diacetyl -> `2_3_butanedione`; 2-methylbutanal (thesis only) -> `2_methylbutanal`; HMF (thesis only)
-> `hmf`; furaneol (thesis only) -> `furaneol`. Glyoxal, methylglyoxal, 2,3-pentanedione,
hydroxyacetone, glycolaldehyde, maltol, the deoxyosones, the Amadori compound, xylose and glycine are
not registry keys (flag 8).

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| k, ARP -> 1-deoxyosone (1-deoxypentosone), first order in ARP | 2.79e-2 | /min | xylose + glycine 1:1 at ~0.24 mol/kg each, 0.5 M phosphate pH 6, 100 C, sealed stirred steel tube, fitted on 30-120 min data with the full network | text p. 370 | measured_rate (multiresponse fit; no CI printed) |
| k, ARP -> 3-deoxyosone (3-deoxypentosone), first order in ARP | 5.50e-4 | /min | same | text p. 370 | measured_rate (multiresponse fit; no CI printed) |
| k(ARP -> 1-DH) / k(ARP -> 3-DH) | 50.7 (mine) | — | same; pH 6 phosphate | derived from the two printed constants | within_study_ratio |
| activation energies, pre-exponential factors, any constant at 80 or 120 C | NOT PRINTED | — | the model has them ("Arrhenius temperature dependence"); none is given | — | — |
| xylose remaining, nor-HDF, 3-deoxyosone, hydroxyacetone vs time | — | µmol per mmol initial sugar | 100 C, pH 6, 30 / 60 / 90 / 120 min | Figure 2 | figure_only |
| glyoxal, methylglyoxal, diacetyl, 2,3-pentanedione, glycolaldehyde, furfural vs time | NOT SHOWN | — | measured per Methods 2.1, not shown or tabulated | Figure 1 names them only | — |
| 3-deoxyosone "formed rapidly in the initial reaction and then suffered a net loss" | — | — | 100 C, pH 6 | text p. 369-370 | level_only (qualitative; matches Parker 2013's "3DH peaks within 40 min" reading of the thesis) |
| glucose + glycine; isoleucine + xylose; pH 4 and 8; HMF, furaneol, maltol, 2-methylbutanal | NOT IN THIS PAPER | — | thesis only (Balagiannis 2015 row 8) | — | — |
| pot composition | xylose 6 mmol + glycine 6 mmol in 25 g, 0.5 M phosphate pH 6 (0.24 mol/kg each, mine) | — | — | Methods 2.1 | measured (as charged) |
| heat-up | "typically 5 min" to set point | min | steel tubes in a heating block | Methods 2.1 | measured (stated) |

**Reading against the trunk (arithmetic, not a fit).** The trunk's hexose Amadori compound
(Martins, glucose + glycine, `parameters.py` k_ama_odg / k_ama_tdg) splits 1-deoxyglucosone :
3-deoxyglucosone = 1.57e-2 : 1.11e-2 /min = 1.4 : 1 at 100 C. Desclaux's pentose Amadori compound
in 0.5 M phosphate at pH 6 splits 1-deoxypentosone : 3-deoxypentosone = 2.79e-2 : 5.50e-4 = 51 : 1,
with the 1-deoxy constant 1.8x the hexose one and the 3-deoxy constant 20x smaller. Two things
differ at once (pentose vs hexose; 0.5 M phosphate at pH 6 vs Martins' buffer), and a fitted
3-deoxyosone formation constant in a network where 3-deoxyosone is also being consumed is only as
good as the sink it was fitted against, so the ratio is recorded as a within-study observation and
not as a correction to the trunk. If the pentose Amadori step on the sulfur lane is ever checked
against it, the check is at 100 C, pH 6, phosphate, and the 1-deoxy branch is the one that
carries the nor-HDF, hydroxyacetone, diacetyl and pentanedione flux in this network.

**What this paper does NOT give the pyrazine supply question.** No glyoxal or methylglyoxal
formation constant, no level, no time course, no barrier. The supply caveat in the B18 report
section 6 and the wishlist entries in `parameters_dicarbonyl.py` (`DICARBONYL_WISHLIST`) are
unchanged by this paper.

## 5. Flags

1. **The numbers Parker 2013 and Balagiannis 2015 attribute to "Desclaux" are in the Reading PhD
   thesis (Desclaux, G., 2006, University of Reading), not in this conference paper.** The four
   pages print two rate constants and four fitted curves; "rate constants and activation energies
   for the whole network", the glucose + glycine and isoleucine + xylose systems, pH 4 and 8, and
   the HMF / furaneol / maltol / 2-methylbutanal outputs are thesis content. The paper itself says
   "detailed mechanisms to be published separately" and "A subsequent paper will report the kinetic
   interpretation"; no such journal paper is on disk and none is cited by Parker 2013 or Balagiannis
   2015 beyond the thesis. The item to fetch for the dicarbonyl supply is the thesis (British Library
   EThOS, or the Reading repository), not another copy of these pages.
2. **`parker2013_extraction.md` row "ref 55" and its section 5 overstate this paper**: "the only
   tabulated glyoxal / methylglyoxal / 2,3-butanedione / 2,3-pentanedione time courses at 80/100/120
   C in water, pH 6" — there is no table, no dicarbonyl curve, and results are shown for 100 C only.
   Parker's own summary row ("3DH peaks within 40 min; hydroxyacetone plateaus after ~3 h") also
   reads beyond these pages (the paper's data stop at 120 min), so Parker was reading the thesis or
   the authors' later work. Those files are outside this dossier's remit and were not edited; the
   `b19_prereg_draft.md` row 33 entry "Desclaux 2006 (glyoxal, methylglyoxal, butanedione time
   courses 80 to 120 C)" should be re-pointed at the thesis.
3. **Only 100 C is reported.** The 80 and 120 C runs are described in Methods and appear nowhere in
   Results; whether the two printed constants come from a 100 C-only fit or from an Arrhenius fit
   across three temperatures evaluated at 100 C is not stated. No barrier can be derived.
4. **No uncertainty on either constant**, no replicate statement, no time-zero point, no
   calibration or recovery for the OPD / ethoxyamine derivatisations. The 3-deoxyosone constant
   (5.50e-4 /min) is a small number fitted to a species that peaks before the first sample and
   declines through all four points; it is strongly correlated with the (unprinted) 3-deoxyosone
   sink constants and should be treated as order-of-magnitude.
5. **0.5 M phosphate.** Phosphate at this strength accelerates the Amadori rearrangement and the
   deoxyosone-forming enolisations by large factors; the constants are not transferable to an
   unbuffered or citrate pot without a phosphate term the trunk does not carry.
6. **Concentration basis.** All Figure 2 quantities are "µmol per mmol initial sugar"; converting to
   mol/L needs the pot density (not printed) and the assumption that 25 g is the whole reacting
   mass. The pot charge (6 mmol in 25 g) is the only absolute number.
7. **"Xylose (sugar)" panel** — the text says the xylose-loss prediction "could be improved"; the
   sugar + amino acid -> ARP step is the one step not first order ("apart from the first stage
   reaction"), and its order and constant are not printed.
8. **Registry gaps against `data/keys/compounds.yml`**: glyoxal, methylglyoxal, 2,3-pentanedione,
   hydroxyacetone, glycolaldehyde, maltol, 1-deoxypentosone, 3-deoxypentosone, the xylose-glycine
   Amadori compound, xylose and glycine have no key. Norfuraneol, furfural, 2,3-butanedione,
   2-methylbutanal, HMF and furaneol do.
9. **What to request from the authors / the thesis**: the full parameter table (rate constants
   with confidence intervals at 80 / 100 / 120 C, activation energies), the measured concentrations
   of glyoxal, methylglyoxal, diacetyl and 2,3-pentanedione against time at the three temperatures
   and three pHs, the calibration and recovery of the OPD derivatisation, and the pot density or
   volume.
10. **OCR layer**: the on-disk PDF is a scan; superscripts, the degree sign ("~" for "C") and all
    figure text are garbled in the text layer. Every number in this dossier was checked against the
    page image.
