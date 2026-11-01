# Designed time-course experiment — protocol draft (north star C8)

*Draft, 2026-10-28, for Pablo to review before it goes to a partner laboratory. Nothing has been sent. Why it
exists: on B60, 94 % of the sulfur lane's flux runs through six FITTED trunk steps
(`SULFUR_FLUX_MAP.md` sec. 5c), and six literature rounds found that the rate data that would fix them do not exist.
An independent review of the repository reached the same conclusion: "nothing else changes identifiability".*

## What it must identify

What the literature cannot give, ranked by the flux it decides:

1. **Pentose and cysteine consumption over time**, and its split into TTCA, the Cys-Amadori compound and the
   deoxyosones. This is the trunk's timescale: the model runs it in seconds, three pots say hours.
2. **The pH law of that split**, in the amine-catalysed regime. Hofmann's pH sweep says MFT is highest at pH 3;
   the model says the opposite.
3. **The temperature dependence** (barriers) of the same steps. 31 sulfur constants are known at 145 °C only.
4. **Norfuraneol and furfural** formation and loss, the MFT and FFT precursors.
5. **H2S supply over time**, which caps both thiols.
6. **MFT and FFT** themselves, by isotope dilution.

## Design (one factorial core, two small fed arms)

**Core:** D-ribose 30 mM + L-cysteine 30 mM in 0.2 M sodium phosphate. Factors:

- **Temperature:** 100, 120, 140 °C.
- **Label pH:** 5.0, 6.0, 7.0, set at 25 °C. Record the pH after cooling at every time point.
- **Vessel:** sealed glass pressure vials, a stated fill and headspace volume, air headspace (record volumes),
  heated in an oil bath or block with a logged sample temperature.
- **Sampling:** 0, 5, 10, 20, 40, 80 and 160 min (at 140 °C: 0, 2, 5, 10, 20, 40, 80). Quench on ice. Separate
  vials per time point.
- **Replicates:** duplicates; triplicates at the centre point (120 °C, pH 6).
- **Size:** 9 conditions × 7 times × 2 ≈ 130 vials, plus the centre replicates.
- **If the budget is tight**, keep the "plus" of the design: 120 °C at pH 5 / 6 / 7, plus pH 6 at 100 / 140 °C.
  That is 5 conditions, ~75 vials, and it still identifies the pH law and the barriers.

**Analytes per vial:**

| analyte | method (suggested) | why |
| --- | --- | --- |
| free cysteine | HPLC-FLD (OPA/thiol derivatisation) or LC-MS/MS | trunk timescale |
| ribose | HPAEC-PAD | trunk timescale |
| TTCA, Cys-Amadori compound | LC-MS/MS with synthesised standards (Cao 2017; Zhai 2019 routes) | the de Roos route |
| norfuraneol, furfural | HPLC-DAD (280 nm) or LC-MS | the MFT and FFT precursors |
| H2S | headspace, ³⁴S-sulfide SIDA (Li, Schieberle & Steinhaus 2016) or methylene blue | sulfur supply |
| MFT, FFT, bis(2-methyl-3-furyl) disulfide | SIDA GC-MS with labelled standards | the observables, quantified |
| pH | after cooling, at 25 °C | the pH trajectory |

**Fed arm A (decisive for the FFT branch):** furfural 1 mM + Na2S 1 mM, 0.2 M phosphate pH 6, at 100 / 120 /
140 °C, 20 min. Furfural and FFT by SIDA. Hofmann's furan-2-aldehyde + H2S protocol at three temperatures: the
experiment no paper has done, and the one that settles B62.

**Fed arm B (checks B61):** norfuraneol 1.75 mM alone, 0.2 M phosphate pH 6 and 7, 120 °C, 0 / 20 / 60 min.
Mikami 2017 measured pH 8 at 100 °C only.

## What we would do with it

The core's time courses become FIT rows for the trunk (B59b): the first rows in the corpus that constrain
timescales rather than endpoints. The fed arms enter as measured constants. Everything is pre-registered before
the data arrive, and the model's blind predictions for every vial are frozen and sent first, as in the earlier
blind test. Success criterion (the usefulness bar, `tasks/data_restructure_plan.md` sec. 7): the sulfur lane beats
the typical-yield baseline on these new-laboratory data, out of sample.

## Open points for Pablo

- The partner laboratory's analytical capability: SIDA standards for MFT and FFT, and LC-MS/MS for TTCA and
  Cys-Amadori.
- Air vs nitrogen headspace: air is food-relevant; nitrogen isolates oxidation. One arm of each at the centre
  point if affordable.
- Whether a xylose twin of the centre point is worth adding. The model treats ribose and xylose as one pentose;
  the data say ribose is faster (Mikami 2017, Furukawa 2013).
