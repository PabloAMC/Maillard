# Balagiannis 2015 — EXTRACTION (review chapter; a citation map of the primary pyrazine-kinetics sources, not a data source)
### Chapter 10 of a Woodhead volume: methodology of multi-response modelling plus a survey of flavour-kinetics studies. It prints no rate constant and only second-hand Ea ranges.

**Source on disk:** `data/articles/balagiannis2015.pdf` (23 pp., owner's download, 2026-09-08). Read
from the text layer (`scratchpad/articles/balagiannis2015.txt`), which is clean. The chapter contains
no numeric tables (Table 10.1 is a classification of volatile classes; Figures 10.1-10.8 are schemes
and a flowchart). Nothing here is primary; every number below is the chapter quoting someone else and
is recorded so the repo knows which primary to fetch, not to be used as a value.

## 0. Identity

| field | value |
|---|---|
| Title | "Predicting aroma formation with kinetic models" |
| Author | D. P. Balagiannis (University of Reading) |
| Venue | Chapter 10 in *Flavour Development, Analysis and Perception in Food and Beverages* (J. K. Parker, J. S. Elmore, L. Methven eds.), Woodhead Publishing / Elsevier, 2015, pp. 211-233 |
| DOI | 10.1016/B978-1-78242-103-0.00010-2 |
| Naming | "APR" (sic) for Amadori rearrangement product in one passage; "Leahy and Reineccius (1989a, b)" = ch. 7 of ACS Symp. Ser. 388 and ch. 18 of ACS Symp. Ser. 409 respectively (both now in the repo as `leahy1989_extraction.md`, `leahy1989a_extraction.md`); the reference list misprints the editors ("Battery" for Buttery, "Parlinent, R.J., McCorrin" for Parliment, T.H., McGorrin) and gives 1989a's pages as 78-91 (actual 76-91) |
| Already in repo | Jousse 2002 (`jousse2002_extraction.md`), Martins & van Boekel 2003/2005 (`martins2003*`, `martins2005*`), Chan & Reineccius 1994 (`chan1994_extraction.md`), Leahy 1989 both chapters, Yu 2018 is NOT cited (published after this chapter) |

## 1. Why it matters

The chapter is useful for exactly one thing here: it names the handful of primary studies that
report pyrazine kinetics (rate law, order, Ea) and says what kind of number each contains. The
engine's pyrazine lane needs measured rates with units; this chapter has none, and its one numeric
pyrazine statement (the Leahy Ea range) is a unit conversion of a sentence in the primary, not of
the primary's table. Its methodological content (multi-response modelling, determinant criterion,
Bayesian priors, "Maillard Ea about 100 kJ/mol") is background the repo already follows.

## 2. What the chapter says about pyrazine (and adjacent) kinetics, with the primary for each

Verbatim where the wording matters; each row is a fetch pointer.

| # | claim in the chapter (short quote) | number as printed | primary source named | what the primary should contain | repo status |
|---:|---|---|---|---|---|
| 1 | Reineccius 1999 "notes that pyrazines are the most studied compounds" and "indicates the need for studies in real food matrices" | — | Reineccius, G. A. 1999, in *Flavor Chemistry: Thirty Years of Progress*, Plenum, pp. 345-352 | a review of flavour-formation kinetics up to 1999 (likely tabulates Leahy, Huang, Chan, Schirle-Keller Ea) | not on disk |
| 2 | "Leahy and Reineccius (1989a, b) studied the formation of pyrazines in model systems where the type of sugar, type of amino acid, pH and aw were varied. ... activation energies were reported assuming zero-order kinetics and were ranged from 113 to 188 kJ/mol." | 113-188 kJ/mol, zero order | Leahy & Reineccius 1989a (ACS 388 ch. 7), 1989b (ACS 409 ch. 18) | k in ppm/h at 75/85/95 C, Ea in kcal/mol per compound and system; pH 5/7/9; aw 0.32-0.84 in NFDM | **on disk, extracted**. 113-188 = 27 x 4.184 and 45 x 4.184: the chapter converted the primary's prose range "27 to 45 kcal/mole"; the table range is 27.3-44.8 kcal/mol = 114.2-187.4 kJ/mol |
| 3 | "Huang et al. (1989) reported zero-order kinetics on the formation of pyrazines in aqueous model systems composed of glucose and various amino acids, heated from 120 to 140 C at pH 10." | zero order; 120-140 C; pH 10; no numbers | Huang, T. C.; Bruechert, L. J.; Ho, C. T. 1989, J. Food Sci. 54 (6), 1611-1614 | zero-order pyrazine rates and Ea for glucose + several amino acids at 120-140 C — the closest published window to the engine's 100-145 C; also cited by Zhou 2024 (ref 28: arginine-glucose, zero order, 120-140 C) and by Jousse 2002 as a Fig. 7 rate source | **not on disk — highest-priority fetch** |
| 4 | "they used three different model systems (aqueous, 80% propylene glycol and ethanol) and high pressure to report zero-order kinetics for the formation of tetramethylpyrazine (Huang et al., 1995). In all their studies, they calculated the relative activation energies." | zero order; no numbers here (Jousse 2002 quotes 79 kJ/mol, 25-55 C, from 3-hydroxy-2-butanone — see `jousse2002_extraction.md` §4.2) | Huang, T. C.; Fu, H. Y.; Ho, C. T. 1995, ACS Symp. Ser. 610, pp. 49-62 | tetramethylpyrazine from acetoin + ammonium acetate, pressure and solvent dependence; a single-product R28-type system | not on disk |
| 5 | "Jusino et al. (1997) studied the pyrazine formation on a dry model system consisting of amioca starch, lysine and glucose ... from 80 to 120 C and concluded that 2,5-dimethylpyrazine and 2-methylpyrazine share the same rate-limiting reaction step. Also, they reported that the kinetics of those compounds was best fitted to a first-order reaction with Ea 56.5 kJ/mol." | first order; Ea 56.5 kJ/mol; 80-120 C; solid | Jusino, M. G.; Ho, C. T.; Tong, C. H. 1997, JAFC 45 (8), 3164-3170 | first-order rate constants for 2-methylpyrazine and 2,5-dimethylpyrazine in a low-moisture starch matrix; Jousse 2002 borrowed this Ea for its R6 carbonyl step | not on disk; relevant to a low-aw / extrusion lane rather than the aqueous engine |
| 6 | Jousse et al. 2002 "were the first who attempted to model the formation of aroma using the multiresponse approach ... the authors manipulated, rather ambiguously 'by hand', the values of the rate constants ... This suggestion might be risky because the model was not statistically validated" | none repeated | Jousse, Jongen, Agterof, Russell, Braat 2002, J. Food Sci. 67, 2534-2542 | lumped PZ class; see `jousse2002_extraction.md` (verdict there agrees with this chapter) | on disk, extracted |
| 7 | "Low (2006) performed studies on a potato model system and observed similar rates in the formation of acrylamide and alkylpyrazines; she concluded that these compounds share a common rate determining intermediate ... extended the model proposed by Wedzicha et al. (2005) ... to include the formation of alkylpyrazines ... strained to include a range of temperatures (160-190 C)" | 160-190 C; no numbers | Low, M. Y. 2006, PhD thesis, University of Reading; the journal version is Low, Parker & Mottram 2007, JAFC 55 (10), 4087-4094 (Zhou 2024 ref 10) | multi-response rate constants for alkylpyrazines in a potato/glycine dry model at 160-190 C | not on disk; the 2007 JAFC paper is the fetchable version |
| 8 | "Desclaux modelled ... xylose/glycine, glucose/glycine and isoleucine/xylose ... 80-120 C ... pH (4, 6 and 8) ... estimated rate constants and activation energies for the whole network ... HMF, furfural, furaneol, nor-furaneol, maltol and 2-methylbutanal" | no numbers | Desclaux, G. 2006 PhD thesis, Reading; Desclaux et al. 2006, in *Flavour Science: Recent Advances and Trends*, pp. 367-370 | no pyrazines named; deoxyosone / dicarbonyl / hydroxycarbonyl kinetics at pH 4-8 — trunk-relevant, not pyrazine-relevant | not on disk |
| 9 | "Cremer and Eichner (2000) ... Strecker reaction followed pseudo-zero-order kinetics ... activation energies ... ranged from 115 to 124 kJ/mol" (70-110 C, low-moisture) | 115-124 kJ/mol | Cremer & Eichner 2000, Food Chem. 71, 37-43 | Strecker aldehyde (R07) rates — the step upstream of R28 | not on disk (also quoted by Jousse for R9) |
| 10 | Chan & Reineccius 1994a,b: Strecker compounds "best described by pseudo-zero-order kinetics, and the activation energies ranged from 60 to 129 kJ/mol" (pH 6-8, 75-115 C) | 60-129 kJ/mol | Chan & Reineccius 1994 (ACS) | see `chan1994_extraction.md` (peak-area caveat there) | on disk |
| 11 | Schirle-Keller & Reineccius 1992: furfural 147, 2-acetylfuran 151, 5-methylfurfural 155, DDMP 129, HMF 118 kJ/mol, zero order, glucose + cysteine 80-150 C | as listed | Schirle-Keller & Reineccius 1992, ACS Symp. Ser. 490, pp. 244-258 | furan channel, not pyrazine | not on disk |
| 12 | "in the Maillard reaction, the activation energy has a value of about 100 kJ/mol. So, the estimates for the activation energy should not deviate much from this value" | ~100 kJ/mol | chapter's own rule of thumb (van Boekel lineage) | a prior, not a measurement | — |
| 13 | Parker 2013, "The Kinetics of Thermal Generation of Flavour", J. Sci. Food Agric. 93 (2), 197-208: "an excellent and thorough review of the most significant kinetic studies" | — | Parker, J. K. 2013 | the most recent consolidated table of flavour-kinetics Ea/rates before 2015; one fetch may cover items 3-5, 7 | not on disk — second-priority fetch |
| 14 | Mundt & Wedzicha 2003 / Wedzicha 1984: the Maillard reaction "can be compressed to a series of reactions with only three rate determining steps" (Figure 10.5: aldose + amino acid -> Int 1 -> Int 2 -> melanoidins, with ketose in parallel) | — | Mundt & Wedzicha 2003, JAFC 51, 3651-3655 | browning kinetics; the skeleton Low 2006 extended to pyrazines | not on disk |

Nothing else in the chapter concerns pyrazines. The remaining citations (van Boekel's Bayesian /
determinant-criterion methodology; Brands & van Boekel sugar-casein; Martins & van Boekel
glucose/glycine; Davies / Leong / Wedzicha sulfite-trapped browning; Knol and De Vleeschouwer
acrylamide; Quintas caramelisation; De Wit & Nieuwenhuijse milk sulfur; Balagiannis 2009/2010
meat-extract Strecker aldehydes and furfural) are trunk or other-lane context.

## 3. Numbers the chapter prints (all second-hand)

| quantity | value as printed | unit | primary | usable as |
|---|---|---|---|---|
| pyrazine formation Ea, sugar + amino acid, zero order, 75-95 C | 113-188 | kJ/mol | Leahy & Reineccius 1989a,b | superseded by the primaries' tables (114.2-187.4 kJ/mol; see the two Leahy dossiers) |
| 2-methylpyrazine / 2,5-dimethylpyrazine Ea, solid starch-lysine-glucose, first order, 80-120 C | 56.5 | kJ/mol | Jusino 1997 | pointer only; fetch the primary |
| Strecker aldehyde Ea, pseudo-zero order, 70-110 C | 115-124 | kJ/mol | Cremer & Eichner 2000 | pointer only |
| Strecker compound Ea, pseudo-zero order, 75-115 C, pH 6-8 | 60-129 | kJ/mol | Chan & Reineccius 1994 | see `chan1994_extraction.md` |
| furan-class Ea, zero order, 80-150 C | 118-155 | kJ/mol | Schirle-Keller & Reineccius 1992 | pointer only |
| generic Maillard Ea prior | ~100 | kJ/mol | chapter | sanity bound only |

No rate constant, no pre-exponential factor, no reaction-order determination for pyrazines is
printed in the chapter.

## 4. What the repository should do with this

1. **Do not cite this chapter for any value.** Cite the primaries; where the chapter's number
   differs from the primary's table (Leahy range), the table wins.
2. **Fetch order for a pyrazine lane**, by closeness to the engine's aqueous 100-145 C window:
   (a) Huang, Bruechert & Ho 1989, J. Food Sci. 54, 1611 (aqueous, 120-140 C, pH 10, zero order,
   several amino acids); (b) Parker 2013, J. Sci. Food Agric. 93, 197 (consolidated review; may
   tabulate a-c); (c) Low, Parker & Mottram 2007, JAFC 55, 4087 (multi-response alkylpyrazines,
   160-190 C, dry potato matrix — for the extrusion/roast end); (d) Jusino, Ho & Tong 1997, JAFC 45,
   3164 (first order, solid, 80-120 C — low-aw lane); (e) Huang, Fu & Ho 1995, ACS 610 (a single
   R28-type condensation, tetramethylpyrazine from acetoin); (f) Cremer & Eichner 2000, Food Chem.
   71, 37 (R07 Strecker rates); (g) Reineccius 1999 (older consolidated review).
3. The chapter's methodological position — that pre-2002 pseudo-zero/first-order fits "do not
   provide any significant mechanistic insight" and that multi-response modelling is "the method of
   choice" — matches the repo's own rule that whole-cascade rates are stored as rates at fixed
   reactant concentrations, not as step constants.

## 5. Flags

1. The Leahy Ea range "113 to 188 kJ/mol" is the primary's prose range in kcal/mol converted at
   4.184, not the table range; harmless but a sign that the chapter did not re-read the tables.
2. Reference-list typos: editors' names of both Leahy chapters; 1989a page range 78-91 vs actual
   76-91; "Aillard Reactions" for Maillard (van Boekel & Berg 1994); "APR" for ARP in the text.
3. "Low (2006)" is a thesis; the repo should target the 2007 JAFC article (already cited by Zhou
   2024 as its ref 10).
4. Item 4's tetramethylpyrazine Ea (79 kJ/mol) is NOT in this chapter; it is quoted by Jousse 2002
   — recorded here only so the two dossiers are not read as independent confirmations.
5. The chapter predates Yu 2018 and Zhou 2024; it has no fed-dicarbonyl pyrazine kinetics at all,
   so the corpus's only measured R07 + R28 rates (Zhou 2024) are not among its pointers.
