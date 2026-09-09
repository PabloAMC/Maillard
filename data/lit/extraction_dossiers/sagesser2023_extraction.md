# Sägesser et al. 2023 — EXTRACTION (true protein of Roquette Nutralys F85M pea isolate, lot W007T, and Solae Alpha 8 soy CONCENTRATE, lot R880002119, alongside eight Chlorellaceae powders: total nitrogen, non-protein nitrogen, sample-specific nitrogen-to-protein factors kA, amino-acid profiles from two laboratories — the profiles printed as a bar chart, the factors as a bar chart, and the numbers "available on request")
### The amino-acid table the repository wanted for Nutralys F85M is not in this PDF: the paper prints one table (true protein, g/100 g: Nutralys 68.9 +/- 2.6, Alpha 8 54.0 +/- 2.3) and puts the amino-acid profiles (Fig. 5), the nitrogen factors (Fig. 3) and the nitrogen balance (Figs 2, 4) in raster figures, with the relative compositions in the online supplement — so no lysine, methionine or arginine per g protein can be typed from it; the only per-protein amino-acid number in the text is "cysteine of approx. 1.2 %" for every powder but one.

**Source on disk:** `data/articles/Sagesser2023.pdf` (10 pp., open access CC BY, owner's download,
2026-09-09). Read from the text layer (`scratchpad/articles/Sagesser2023.txt`, 622 lines, clean)
and, for the figures, from `pdfimages -list` and a page-8 render: Fig. 1 (analytical scheme), Fig. 2
(TN and PN by four routes), Fig. 3 (kp and kA, two laboratories), Fig. 4 (NPN composition; PN by
difference), Fig. 5 (amino-acid profiles, g/100 g protein, stacked bars for powders A-K) are each
one embedded bitmap with no text layer — FIGURE-ONLY throughout; nothing was read off them.
Table 1 (true vs crude protein) came through clean and is re-typed below. The **supplementary
material** (powder colours, primer sequences, microscopy, phylogenetic trees, and "the relative amino
acid compositions of both analyses") is NOT on disk; "Data will be made available on request."
Repo status before this dossier: `sagesser2024_extraction.md` flag 2 says this paper "holds, for
the same lots, the total amino-acid profile, the sample-specific nitrogen factor, the NPN and the
NSI ... it would give Nutralys F85M's lysine per g true protein, which no other dossier has" — it
holds them as figures, not as numbers (flag 1). `data/species/protein_matrices.yml` charges
`pea_isolate` with amine 0.47 [0.40, 0.70] and `soy_isolate` with 0.36 [0.31, 0.58] mmol lysine
per g protein from Jaeger 2023 and Gorissen 2018, both on an N x 6.25 basis.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "A novel approach for the protein determination in food-relevant microalgae" |
| Authors | Corina Sägesser, Johanna M. Kallfelz, Samy Boulos*, Laila Hammer, Lukas Böcker, Reto Portmann, Laura Nyström, Alexander Mathys* — ETH Zurich (Sustainable Food Processing; Food Biochemistry), Agroscope Bern, FFSH Zürich, Wageningen |
| Venue | Bioresource Technology 390 (2023) 129849; received 8 Sep 2023, revised 6 Oct, accepted 6 Oct, online 7 Oct 2023 |
| DOI | 10.1016/j.biortech.2023.129849 |
| Products | **K** = "pea protein isolate Nutralys F85M (Roquette, Lestem [sic; Lestrem], France; lot code: W007T)"; **J** = "soy protein concentrate Alpha 8 (Solae LLC, Missoursi [sic], USA; lot code: R880002119)" — a CONCENTRATE, not an isolate; **A-H** = eight commercial Chlorellaceae powders (A photoautotrophic *C. sorokiniana*, Roquette Klötze; B, C, D heterotrophic *C. sorokiniana*, Langyatai / Duplaco / Allmicroalgae; E *C. vulgaris*, Algenuity; F, G, H *Auxenochlorella protothecoides*, Daesang / Sophie's Bionutrients / Alver), species confirmed by rbcL / ITS sequencing. Same lots as `sagesser2024_extraction.md` (K = PPI, J = SPC, F = MAB1, H = MAB2) |
| Nitrogen factors | kp = sum of anhydrous amino acids / total N; kA = sum of anhydrous amino acids / protein N; **true protein = kA x (TN - NPN)** with a sample-specific kA (the average of the Eurofins-based and the in-house-based kA); crude protein for microalgae = TN x 4.78; 6.25 is discussed and not used |
| Naming | TN total nitrogen; PN protein nitrogen; NPN non-protein nitrogen (nitrate, ammonia, urea, nucleic acids, chlorophyll); NSI nitrogen solubility index; AAA anhydrous amino acids; Asx = Asn + Asp, Glx = Gln + Glu; IAA indispensable amino acids; IAAS indispensable amino acid score (WHO 2007) |
| Tables / figures | 1 table (Table 1, p. 8); Figs 1-5 (all raster) |
| Compound registry | none of the analytes is a volatile; lysine maps to the engine's amine pool (`protein_matrices.yml`), cysteine to the thiol / disulfide pools; `data/keys/compounds.yml` has `reactive_lysine` and no free amino-acid keys |

## 1. Why it matters

The repository wanted this paper for one thing: the total amino-acid composition of a NAMED
commercial pea isolate (Nutralys F85M, the product on which the Snel 2023 ketone-binding constants
in `data/lit/binding_constants.yml` were measured) with lysine, arginine, cysteine, methionine and
the amide pool per g protein, on a nitrogen factor better than 6.25 — to set beside Jaeger 2023's
anonymous Naturz isolates (pea 0.539, soy 0.410 mmol lysine per g N x 6.25 protein) and Gorissen
2018's pooled market means (0.40 / 0.31). The paper measured all of it, twice (Eurofins IC-UV with
oxidative hydrolysis for Cys and Met; in-house AccQTag UHPLC with internal standards), and printed
none of it as a number: Fig. 5 is the composition, in g/100 g protein, as a stacked bar chart, and
the supplement holds the "relative amino acid compositions". What the text does print for the two
plant powders is the true protein content (Table 1: Nutralys 68.9 +/- 2.6, Alpha 8 54.0 +/- 2.3
g/100 g), the NPN class ("below 2 %" of TN), the NSI band (22-34 %), the IAAS (Nutralys 60-65 %,
Alpha 8 about 80 %, both limited by the sulfur amino acids), and one amino acid: "Apart from powder
H, all samples contained comparable amounts of cysteine of approx. 1.2 %" — read as g per 100 g
protein, 0.099 mmol cysteine per g true protein (mine), by an oxidative hydrolysis that, unlike
Gorissen's, recovers cysteine. The one structural fact the paper adds to the isolate question is the
protein basis: Nutralys F85M at 68.9 g true protein per 100 g against a nominal 85 % on N x 6.25
means any per-g-protein density taken from this paper's data will read about 15-20 % higher than the
same powder's density on the 6.25 basis the matrix table uses (section 4, flag 3). The microalgae
work — the paper's actual subject — is out of scope here: it shows that neither TN x kp nor the
sum of amino acids gives the true protein of Chlorellaceae because NPN (3.4-15.4 % of TN, mostly
ammonia, nitrate and nucleic acids) and hydrolysis losses pull in opposite directions, and
proposes kA x (TN - NPN) with kA 4.8-5.7 (average 5.3), which happens to reproduce TN x 4.78
within +/- 10 % on average.

## 2. Methods as they matter to a model

- **Materials.** Powders as received; lots as in section 0. Microscopy (supplement) showed no intact
  cells in the plant powders and "powder K was very homogeneous". Nothing on the isolates' process.
- **Total nitrogen.** Chemiluminescence after combustion of suspensions (Shimadzu TNM-L, at least
  triplicate) and Kjeldahl (duplicate); "Elemental analysis and Kjeldahl lead to comparable results
  for total nitrogen." Values FIGURE-ONLY (Fig. 2).
- **Nitrogen solubility index.** 1 % (w/w) suspension, 90 min hydration at room temperature, 45 mL
  centrifuged at 10,000 rcf for 15 min; NSI = supernatant N / total N of the suspension. Plant
  proteins 22-34 % (text); per-powder values not printed.
- **Amino acids, laboratory 1 (Eurofins Scientific AG, Schönenwerd; verbatim core).** "Eurofins
  executed an alkaline hydrolysis for tryptophan based on the ISO method 13904:2016 and an oxidative
  hydrolysis for cysteine and methionine as well as an acidic hydrolysis for all remaining amino
  acids based on the ISO method 13903:2005. The hydrolysed amino acids were quantified via ion
  exchange chromatography (IC-UV) resp. liquid chromatography (LC-FLD) for tryptophan. Eurofins did
  single determinations without use of internal standards. Three samples were chosen to be analysed
  in duplicates." Asx and Glx as sums. For powders **K and J the Eurofins protein nitrogen came out
  HIGHER than the total nitrogen** ("notably higher than TN for two of the samples (powder K and
  J)"), which is what triggered the second analysis (flag 4).
- **Amino acids, laboratory 2 (in-house / Agroscope; verbatim core).** "The tryptophan analysis was
  performed as described in Walther et al. (2022), which involved alkaline hydrolysis at 110 C for
  20 h and subsequent detection by UHPLC-UV. All other amino acids were measured according to ISO
  4214 | IDF 254:2022 ... after acidic hydrolysis at 110 C for 24 h, hydrolysates were neutralized,
  derivatized with AccQTag Ultra reagent (Waters), and amino acids measured by UHPLC. The internal
  standards methyl-tryptophan (for tryptophan) and L-Norvaline (for all other amino acids) were used
  to correct for losses." No oxidative hydrolysis step is described for this second run, so its
  cysteine and methionine are lower bounds in the Gorissen sense (flag 5). Four samples in
  triplicate, the rest in duplicate. "PN of Eurofins was in average 12 % higher than PN determined
  internally"; "the relative amino acid compositions of both analyses were very similar (see
  supplementary material)".
- **Amide nitrogen.** Mossé et al. 1985: 1 g powder + 20 mL 2 M HCl, 115 C, 3 h, sealed Pyrex;
  steam distillation into boric acid, titration with 0.025 M HCl; free ammonia by the same
  distillation without hydrolysis at pH 10; Namide = difference; **degree of amidation = mol Namide /
  (mol Glx + mol Asx)**; triplicate; controls showed "high recovery of Asn and Gln amide nitrogen and
  no ammonia release of other amino acids". The Asn / Asp and Gln / Glu split in Fig. 5 assumes the
  same degree of amidation in both pairs. Per-powder degrees of amidation NOT printed (flag 2).
- **Non-protein nitrogen.** Urea (OPA / primaquine colorimetry at 430 nm after 100 C / 10 min in 1 M
  carbonate pH 10); nitrate and ammonia by MQuant strips (nitrate on unheated suspensions); nucleic
  acids by A260 after bead-beating and Zymo column isolation, N content taken as 16.48 % (mean of
  DNA 16.84 % and RNA 16.12 %); chlorophyll assumed (1 % of N dark green, 0.4 % bright green, 0
  yellow / white). Triplicate (nitrate duplicate). For the plant powders NPN "was below 2 %" of TN
  (Fig. 4a, text); nitrate "was only detected in powder C and K" (i.e. the pea isolate carries some
  nitrate; amount FIGURE-ONLY).
- **Protein content.** kp = sum AAA / TN (eq. 1); kA = sum AAA / PN (eq. 2); true protein = kA x (TN
  - NPN) (eq. 3), kA per powder = average of the two laboratories' kA. "kA of different samples
  ranged from 4.8 to 5.7"; microalgal average 5.3. The kA of K and J individually are FIGURE-ONLY
  (Fig. 3).
- **Nutritional score.** IAAS = share of the limiting IAA in the test protein / share in the WHO
  2007 requirement pattern (eq. 4); digestibility not included.
- **Statistics.** "Due to small sample size of n <= 3 no statistical analyses were executed";
  means +/- SD.

## 3. Tables re-typed

### Table 1. "Protein content (g/100 g) based on protein-nitrogen (PN = TN - NPN) multiplied with sample-specific kA, here called 'true protein', vs. total nitrogen multiplied by generic kp of 4.78 suggeted [sic] for microalgae by Lourenço et al. (2004) here called 'crude protein'."

| group | powder | True protein (PN x kA), g/100 g | Crude protein (TN x 4.78), g/100 g | Crude / true protein |
|---|---|---|---|---|
| Plant proteins | **K** (Nutralys F85M pea isolate) | **68.9 +/- 2.6** | — (not applied) | — |
| Plant proteins | **J** (Alpha 8 soy concentrate) | **54.0 +/- 2.3** | — (not applied) | — |
| *A. protothecoides* | F | 44.3 +/- 1.0 | 46.9 +/- 0.3 | 106 % |
| *A. protothecoides* | G | 41.7 +/- 1.6 | 44.7 +/- 0.6 | 107 % |
| *A. protothecoides* | H | 41.7 +/- 2.0 | 43.4 +/- 0.7 | 104 % |
| *C. sorokiniana* | B | 52.5 +/- 1.8 | 47.1 +/- 0.8 | 90 % |
| *C. sorokiniana* | A | 43.4 +/- 2.3 | 39.9 +/- 1.5 | 92 % |
| *C. sorokiniana* | D | 26.6 +/- 1.8 | 24.2 +/- 1.2 | 91 % |
| *C. sorokiniana* | C | 26.0 +/- 1.1 | 24.3 +/- 0.6 | 93 % |
| *C. vulgaris* | E | 25.1 +/- 1.4 | 25.2 +/- 0.7 | 100 % |
| Microalgal average | — | 37.7 | 36.9 | 98 % |

The crude-protein column is left blank for K and J in the print (4.78 is a microalgal factor). The
basis of "g/100 g" (as-received powder or dry matter) is not stated in the caption or the Methods
(flag 6). The +/- is the SD of n = 2-3 determinations of the nitrogen terms; the kA uncertainty is
not propagated in any stated way.

### Numbers printed in the running text (all figures otherwise FIGURE-ONLY)

| item | as printed | where | applies to |
|---|---|---|---|
| NSI, plant proteins | 22 % to 34 % | 3 (Results, first paragraph) | J and K as a band; which is which not stated |
| NSI, microalgae | 12 % to 58 % overall; *C. sorokiniana* 12-19 %; *A. protothecoides* 40-58 %; *C. vulgaris* 34 % | same | A-H |
| Eurofins PN vs TN | "notably higher than TN for two of the samples (powder K and J)" | 3.1 | K, J |
| Eurofins PN vs in-house PN | Eurofins "in average 12 % higher" | 3.1 | all |
| PN by difference (TN - NPN) vs the two amino-acid PN | "on average 4 % lower and 7 % higher than the PN based on Eurofins' amino acid profile and our internal amino acid quantification, respectively" | 3.3 (Fig. 4b) | all |
| kA, microalgae | average 5.3; range over all samples 4.8 to 5.7 | 3.2; Conclusion | A-H (range possibly including J, K; not stated) |
| kp for crude protein | 4.78 | 3.4 | microalgae |
| NPN, microalgae / plant proteins | 3.4 % to 15.4 % / "below 2 %" of TN | 3.3 | A-H / J, K |
| nucleic-acid N, microalgae | 1 % to 5.8 % of TN, "generally lower for plants" | 3.3 | — |
| nucleic acids, w/w | powder A 2.5 %, powder G 2.1 % (the two above the 2 % guideline) | 3.3 | A, G |
| nitrate detected | "only detected in powder C and K" | 3.3 | C, K |
| crude vs true protein | within +/- 10 %; 36.9 vs 37.7 g/100 g on average | 3.4 | microalgae |
| IAA content | *C. sorokiniana* "above 39.8 g / 100 g protein"; egg 40.5-42.7 (Attia 2020); *C. vulgaris* "comparable levels to the two plant proteins"; *A. protothecoides* lowest | 3.5 | — |
| IAAS | powders F to K limited by the sulfur amino acids; "For powder F, G and J the score is around 80 % while it is between 60 and 65 % for powders H and K" | 3.5 | J about 80 %; K 60-65 % |
| **cysteine** | "Apart from powder H, all samples contained comparable amounts of cysteine of approx. 1.2 %. Powder H contained 0.8 %" | 3.5 | all incl. J and K; basis presumably g/100 g protein (Fig. 5's axis) |
| Arg and Glx, powder H vs A | 3.8x and 3.2x | 3.5 | microalgae |

## 4. Numbers the repository can use

Molar masses: Lys 146.19, Arg 174.20, Met 149.21, Cys 121.16 g/mol. "True protein" = kA x (TN -
NPN) — a lower denominator than N x 6.25 for the same powder.

| product | quantity | value +/- sd | unit as printed | mmol per g protein (arithmetic) | method | source | evidence class |
|---|---|---|---|---|---|---|---|
| **Nutralys F85M (K), lot W007T** | **true protein** | **68.9 +/- 2.6** | g/100 g (basis unstated) | — | kA x (TN - NPN); TN by TNM-L and Kjeldahl; kA from two amino-acid analyses | Table 1 | measured |
| **Alpha 8 soy concentrate (J), lot R880002119** | **true protein** | **54.0 +/- 2.3** | g/100 g (basis unstated) | — | as above | Table 1 | measured |
| K, J | total lysine, arginine, methionine, Asx, Glx, every other amino acid | NOT PRINTED (Fig. 5 bars, g/100 g protein; supplement) | — | — | Eurofins IC-UV + in-house AccQTag UHPLC, averaged | Fig. 5 | figure_only |
| K, J | nitrogen-to-protein factor kA (and kp) | NOT PRINTED individually; all-sample range 4.8-5.7 | — | — | eq. 1-2 | Fig. 3; Conclusion | figure_only (range: level_only) |
| K, J | total nitrogen | NOT PRINTED | — | — | TNM-L; Kjeldahl | Fig. 2 | figure_only |
| K, J | NPN | "below 2 %" of TN | % of TN | — | nitrate, ammonia, urea, nucleic acids, chlorophyll | text 3.3; Fig. 4a | level_only |
| K, J | NSI | 22-34 % (band for the two, unassigned) | % of N | — | 1 % w/w, 10,000 rcf | text | level_only |
| K, J (and A-G) | **cysteine** | **approx. 1.2** | % (read as g/100 g true protein) | 12 mg/g / 121.16 = **0.099 mmol/g protein** (mine); if the Eurofins oxidative hydrolysis dominates the average this is total half-cystine as cysteic acid, i.e. free thiol + 2 x disulfide | oxidative hydrolysis (Eurofins) averaged with plain acid hydrolysis (in-house) | text 3.5 | level_only (one rounded value stated for nine powders; not powder-specific) |
| K | IAAS | 60-65 %, limiting = Met + Cys | — | — | eq. 4, WHO 2007 | text 3.5 | level_only |
| J | IAAS | around 80 %, limiting = Met + Cys | — | — | as above | text 3.5 | level_only |
| K | implied kA | about 5.1-5.2 | — | — | 68.9 g true protein / (13.6 g N per 100 g x 0.98), where 13.6 = 85 / 6.25 from the "F85" grade name (nominal, not measured here) and 0.98 = 1 - NPN | derived from Table 1 + the product name | derived_assumption (mine; falls with the nominal protein: at a true N x 6.25 of 80 % the implied kA is 5.5) |
| K, J | free amino acids, free sugars, moisture, ash | NOT MEASURED in this paper (free amino acids and ash are in `sagesser2024_extraction.md`, figure-only there) | — | — | — | — | — |

**Comparison the matrix layer can print (arithmetic, not a fit).** Nothing here changes the amine
pool numbers: this paper prints no lysine. What it prints is the denominator. For the same kind of
powder the three isolate sources sit on three protein bases — Jaeger 2023 PPI 81.22 g N x 6.25
protein per 100 g DM; Gorissen 2018 pea 80 % N x 6.25; Sägesser 2023 Nutralys F85M 68.9 g true
protein per 100 g. If Nutralys is near its nominal 85 % on N x 6.25, the true-protein basis is
0.81x the 6.25 basis, so a lysine density read from this paper's Fig. 5 or supplement would be
1.23x the same powder's density on the matrix table's basis (1.18x at a true 6.25-basis of 80 %).
Concretely: the table's pea amine 0.47 mmol per g N x 6.25 protein would appear as about 0.55-0.58
mmol per g true protein in this paper's units; the band [0.40, 0.70] would read [0.47, 0.86]. Any
number eventually taken from the supplement must be divided back by kA/6.25 x (1 - NPN/TN) before
it enters `protein_matrices.yml`, or the table's basis must be changed for every entry at once
(flag 3). For soy the comparison is between product classes: Alpha 8 is a concentrate at 54.0 g
true protein per 100 g against Jaeger's isolate at 89.22 and Gorissen's seven "isolates" at a
figure-only mean near 74 %; per g protein the lysine of a concentrate and an isolate from the same
process are close, per g powder they are not.

**Cysteine.** The text's "approx. 1.2 %" (0.099 mmol/g protein, mine) is the only cysteine value
on disk for a NAMED pea isolate by a method that recovers cysteine (oxidative hydrolysis, at least
in one of the two averaged runs). It sits at the top of the soy half-cystine values the matrix
table carries (Ruan 2014 / Shimada 1988: 0.100-0.114 mmol/g protein, Ellman's) and above the pea
values (Gao 2020: free SH + 2 x S-S = 0.061-0.074 mmol/g protein) — but it is one rounded number
the authors give for nine different powders at once, so it is a class level, not a Nutralys
measurement, and it is not entered as one. Gorissen's 0.021-0.022 mmol/g (no oxidation) is
confirmed as a method floor, not a composition.

## 5. Flags

1. **The amino-acid table is not in the PDF.** Fig. 5 carries the profiles of all eleven powders as
   stacked bars in g/100 g protein (raster, no text layer); the supplement carries "the relative
   amino acid compositions of both analyses"; the numbers are "available on request". No lysine,
   arginine, methionine, Asx or Glx value for Nutralys F85M or Alpha 8 can be typed from this paper.
   `sagesser2024_extraction.md` flag 2 ("It is not on disk; it is the item to fetch ... it would give
   Nutralys F85M's lysine per g true protein") should be amended: the paper is on disk and does not
   give it; the supplement (biortech.2023.129849 Appendix A) or the authors' data are what would.
   That file was not edited.
2. **The amide split is a measured quantity the paper does not print.** The degree of amidation
   (mol amide N per mol Asx + Glx) was determined in triplicate for every powder and is used to
   split Asn / Asp and Gln / Glu in Fig. 5, but no value appears in text or table — so the
   asparagine pool (acrylamide route) of the pea isolate, which this paper alone among the isolate
   sources resolves, is figure-only.
3. **Protein basis.** True protein = kA x (TN - NPN) with kA about 5.1-5.5 for the legume powders
   (inferred; the printed range over all samples is 4.8-5.7) against N x 6.25 in Jaeger 2023,
   Gorissen 2018, Gao 2020, Tang 2024. The same lysine content reads 15-25 % higher per g protein
   on this basis. Every per-g-protein number that ever comes from this paper's data must carry
   "per g true protein (kA basis)" and be converted before it meets the matrix table.
4. **The Eurofins amino-acid sums exceeded the total nitrogen for K and J.** The paper says so and
   used it as the reason to re-run the analysis in-house (which came out 12 % lower on average,
   without oxidative hydrolysis). The averaged profile in Fig. 5 therefore mixes a run that
   over-recovers (or over-calibrates: single determinations, no internal standard) with one that
   under-recovers. For the two plant powders specifically, the Eurofins absolute values are
   suspect; the relative profile is what the authors trust ("very similar" between runs).
5. **Cysteine and methionine come from two incompatible hydrolyses averaged.** Eurofins used an
   oxidative hydrolysis (cysteic acid / methionine sulfone, the right method); the in-house ISO
   4214 run describes a plain 24-h acid hydrolysis with norvaline correction, which destroys most
   cysteine. If the "approx. 1.2 %" cysteine is the average of the two, the oxidative value is
   higher; if it is Eurofins' alone the paper does not say. Treat 0.099 mmol/g protein as a
   class-level figure with a factor-of-1.5 uncertainty in either direction.
6. **Basis of "g/100 g" in Table 1 not stated** (as received or dry matter). Moisture was not
   measured in this paper (Karl Fischer water is in the 2024 companion, figure-only). For an
   isolate at 5-8 % moisture this is a 5-8 % ambiguity on 68.9 and 54.0.
7. **Alpha 8 is a soy protein CONCENTRATE**, not an isolate, at 54.0 g true protein per 100 g; the
   repository's `soy_isolate` entry should not take any per-g-powder number from it. Per g protein
   its lysine would be comparable to an isolate's, but that number is not printed either.
8. **No statistics, n = 2-3.** "No statistical analyses were executed"; the +/- in Table 1 is an
   SD of duplicates or triplicates of the nitrogen terms, and the kA averaging of two laboratories
   that disagree by 12 % is not propagated into it.
9. **Nitrate was detected in powder K** (the pea isolate) — a process residue that a nitrogen
   balance must subtract and that an N x 6.25 protein figure counts as protein; the amount is
   figure-only (Fig. 4a) and inside the "< 2 % NPN" statement.
10. **What to request from the authors**: the per-powder amino-acid table for K and J from both
    laboratories (g/100 g powder and g/100 g protein), the per-powder kA, kp, TN and NPN, the
    degree of amidation, and the moisture basis of Table 1. With those, Nutralys F85M becomes the
    one isolate on disk with named lot, resolved Asn / Asp, cysteine by oxidative hydrolysis and a
    measured nitrogen factor.
11. **Registry**: `data/keys/compounds.yml` is a volatile registry; lysine appears only as
    `reactive_lysine`; cysteine, methionine, arginine, asparagine and the nitrogen species have no
    key, and none is needed for this paper's use (the matrix table, not the registry, carries the
    pools).
12. **Two misprints**: "Lestem" for Lestrem (Roquette's site) and "Missoursi" for Missouri; lot
    codes as printed. The abstract's "Neither crude protein (kp x TN) nor sum of amino acids
    accurate" is the paper's finding for microalgae, not for the plant references.
