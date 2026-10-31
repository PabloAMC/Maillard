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
| H2S | ³⁴S-sulfide SIDA after in-vial derivatisation (Li, Schieberle & Steinhaus 2017, Eur Food Res Technol 243:69; NOT headspace: cool the sealed vial, add the ³⁴S standard and 3-buten-2-one through the septum before opening, ~30 % of H2S sits in the headspace at pH 5; 1-10 mL aliquots cover ~0.01-0.9 mmol/L; validated at pH 7.4 only; one day per batch; `li2016b_extraction.md`) or methylene blue | sulfur supply |
| MFT, FFT, bis(2-methyl-3-furyl) disulfide | SIDA GC-MS with labelled standards | the observables, quantified |
| pH | after cooling, at 25 °C | the pH trajectory |

**Fed arm A (decisive for the FFT branch):** furfural 1 mM + Na2S 1 mM, 0.2 M phosphate pH 6, at 100 / 120 /
140 °C, 20 min. Furfural and FFT by SIDA. Hofmann's furan-2-aldehyde + H2S protocol at three temperatures: the
experiment no paper has done, and the one that settles B62.

**Fed arm B (checks B61):** norfuraneol 1.75 mM alone, 0.2 M phosphate pH 6 and 7, 120 °C, 0 / 20 / 60 min.
Mikami 2017 measured pH 8 at 100 °C only.

**Buffer arm C (added 2026-10-29; separates phosphate catalysis from pH).** The core is in 0.2 M phosphate, and
phosphate may itself catalyse the Amadori steps: Davidek 2002 measures 1.3-3.5× faster Amadori degradation in
0.1 M phosphate at 90 °C (same laboratory, pH-stat), and the two laboratories that measured this route's
intermediates disagree ~2-5× on its formation rate in exactly the way phosphate would explain (Xie group, 0.2 M
phosphate; Jiangnan group, unbuffered water). No paper separates the two (search record
`LITERATURE_SEARCH_2026-10-29_phosphate_thiazolidine.md`). At the centre point (120 °C) and at pH 6 and 7.4: the
core recipe in **unbuffered water** (pH adjusted with NaOH, pH logged at every time point), in **0.05 M and 0.2 M
phosphate**, and in a **non-carboxylate, non-phosphate buffer at equal pH** (a Good's buffer such as MES or HEPES;
its own reactivity checked in a sugar-free blank). Analytes as the core, TTCA and Cys-Amadori first. 4 media × 2
pH × 7 times × 2 ≈ 110 vials; if the budget is tight, pH 6 only (~55 vials). It also gives the pH trajectory of an
unbuffered pot, which no paper on disk reports.

**Cells added by the conflict analyses (2026-10-29, `docs/reference/conflict_analyses/`):** each settles a recorded
conflict the literature cannot.
- **pH 4 cell** at 140 C in the core (norfuraneol and MFT vs pH; CF-012: the model's norfuraneol at pH 3 is 55x below
  pH 5).
- **Standards before vs after heating**, in sister vials at the centre point: labelled MFT / FFT added before heating
  in one vial, after it in the other (CF-013, CF-015: the thiols' loss during heating, which bounds the sink).
- **A 3:1 xylose : cysteine cell** at the centre point (CF-015: the recipe that favours FFT).
- **Hofmann's exact pot in triplicate** (ribose 100 / cysteine 33 mmol/L, 0.5 M phosphate pH 5, 145 C ramp; CF-016:
  furfural 13x apart between two papers on the same run).

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

## Additions from the 2026-10-31 evidence map (mechanisms the literature cannot discriminate)

The within-study factor matrix and mechanism map ([EVIDENCE_MAP_2026-10-31.md](EVIDENCE_MAP_2026-10-31.md)) tell, per open mechanism, whether published labelling, fed-intermediate, trapping or measured-intermediate experiments can decide it.
- **Decidable from the literature (≥ 2 labs):** the hexose hub; the existence of a cysteine trap of 3-DX; MFT routes at 95 - 140 °C; disulfide reversibility; the Strecker donor ranking; the pyrazines' carbon origin.
- **Single-lab or no discriminating experiment:** the pentose split and its pH and phosphate laws; the MFT lifetime at process temperature; amine acceleration of β-cleavage; dilute cysteine fate; the Strecker O2 acid route; the pyrrole sink; Wolff vs Amadori oxidation; lipid carbonyls as Strecker donors.

Ranked by (mechanisms discriminated) × P(it changes a structural choice) ÷ cost. Each cell is either already in the core above or added to it:

| # | experiment | decides | relation to the core | P (changes a structure) |
|---|---|---|---|---|
| E1 | **Pentose hub ladder, in-situ vs ex-situ OPD.** Xylose + Gly 50 mM, pH-stat 5.5 / 6.5 / 7.5, phosphate 0 / 0.2 M, Cys 0 / 5 / 25 mM, 100 °C (+ 120 °C corner), 9 times to 240 min | the pentose split as a FLUX (not a standing level); the 3-DX cysteine trap and its pH law; phosphate catalysis; the OPD operator | new (the core has no OPD arm); its phosphate factor extends arm C | 0.6 |
| E3 | **Labelled MFT spiked into a running ribose + Cys pot (CAMOLA, [13C5]-ribose 1:1).** 100 vs 140 °C, air vs argon, [2H3]-MFT spike at t = 0 vs at sampling | the MFT loss law (lifetime, O2 and carbonyl dependence); the time-resolved route mix; the thiol work-up operator | runs ON the core recipe; adds the spike and label arms | 0.65 |
| E5 | **Forced-input Strecker / pyrazine pot, second lab.** MGO or butanedione 10 mM + Phe / Leu / Ala / Gly; pH 3.2 / 7.0; argon vs air; closed vial vs Likens-Nickerson; [13C2]-glyoxal and 15N-Ala arms | Strecker flux vs yield; the O2 acid route (cross-lab for Hofmann 2000); the aminoketone partition to pyrazine vs pyrrole; the SDE operator behind the 20x Weenen / Hofmann gap (CF-073) | new; feeds the forced-input Strecker layer (EX-049) directly | 0.6 |
| E4 | **Cysteine / sulfur mass balance with an H2S time course,** dilute (10 / 10 mM) and concentrated (83 / 165 mM) | dilute cysteine fate (catalysed hydrolysis vs K(pH) vs adducts); H2S supply | the core plus arm C cover most of it; add the dilute level and argon | 0.55 |
| E2 | **Fed 1-DG / 3-DG (3-DX if available) decay, ± Gly 0 / 20 / 100 mM, ± phosphate, ± Cys** | amine acceleration of β-cleavage; the deoxyosone lifetimes; the 3-DX sink's identity | new, small (published syntheses) | 0.55 |
| E6 | **Free sugar vs Amadori compound under argon / air / O2, α- vs β-amino acid, ± DTPA** | Wolff vs Amadori oxidation; Namiki C2; 1,4-dideoxyosones; oxidative α-cleavage acids (second lab) | new | 0.45 |
| E7 | **Lipid vs sugar dicarbonyl as Strecker donor in one pot, ± MFT spike** | lipid carbonyls as Strecker donors (cross-lab); as thiol traps | new | 0.4 |
| E8 | **Water activity × thermal history on one recipe** (ribose / Cys / Gly on cellulose, a_w 0.3 / 0.6 / 0.9, ramp vs isothermal) | water activity; thermal history; the low-water C2 + C3 MFT route | new; product-facing | 0.35 |

**If one experiment only:** E3 on the core recipe. It decides the likeliest missing thiol element (the MFT loss law) and the route mix over time in one run. It needs a [2H3]-MFT standard and [13C5]-ribose, both standard reagents for the groups that run SIDA. **If two:** add E1. The pentose hub is the most upstream open question, and its in-situ OPD arm tests a method artefact that otherwise confounds every published deoxyosone level. P(E1 + E3 together resolve at least two of the hub / thiol questions) ~0.7; P(the literature alone does) ~0.3.

Candidate partners have the methods: the Wageningen group (Fogliano), whose multiresponse kinetics come from van Boekel's tradition; Granvogl / Hofmann (Munich: SIDA, CAMOLA, Strecker); Reading (Parker / Mottram: thiols, pyrazine isotopomers); Nestlé / Davidek (deoxyosone feeding).
