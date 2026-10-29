# Citation contexts of the remaining literature candidates (2026-10)

Run 2026-10-27. Purpose: before downloading more papers, check **how and why** each remaining candidate from the
citation-chaining search is cited, so that only papers cited for usable quantitative content are fetched. Search record
only; nothing was added to `data/`.

## Method

- **Candidates (124 unique DOIs or titles)**: every row of the round-1 class A and B tables and the round-2 class A and B
  tables in [LITERATURE_SEARCH_2026-10.md](LITERATURE_SEARCH_2026-10.md), plus the candidate primary papers named in
  [tasks/data_restructure_plan.md](../../tasks/data_restructure_plan.md) from the Sohail 2022 review (7) and from the
  Sun 2022 review's candidate list (section 6 of `sun2022_extraction.md`, 15). Duplicates across lists are merged; the
  "source list" column shows every list a paper was on (`R1-A` = round 1 class A, `Sun2022-lead`, and so on).
- **Exclusions (43 on disk)**: DOIs in `data/keys/papers.yml`; the DOIs printed on the first two pages of every
  `data/articles` PDF (`pdftotext -l 2`); and, for old papers that print no DOI, an author-year file-name match confirmed by
  title words on page 1 (13 papers, e.g. `zhang1988.pdf`, `farmer1992.pdf`, `mottram1995c.pdf`). The list is at the end.
  **81 candidates remain**, and are classified below.
- **Source**: Semantic Scholar Graph API, anonymous, no key or e-mail, about 1 request/s with back-off,
  `GET /graph/v1/paper/DOI:<doi>/citations?fields=title,year,externalIds,contexts,intents,isInfluential&limit=100`
  (a second page at offset 100 where there were more citers; `externalIds` was added to the requested fields only so that
  citers in our corpus could be recognised). This returns, per citing paper, the sentences around the citation
  (contexts) and S2's intent labels. Cache in the session scratchpad, not in the repo.
- **Reading**: per candidate, contexts were sorted with citers in our corpus first (DOI in the registry or on a PDF),
  then contexts with numbers or quantitative words, then most recent; the top 5-16 were read in full and every context was
  also keyword-scanned for thiol/H2S mentions and for contradiction language ("not confirmed", "contradict",
  "in contrast", ...). The "contexts / citers" column gives what S2 holds, not only what was read.
- **Where S2 has no contexts** (class ?): OpenAlex title + abstract only, and said so in the row.
- **Classes**: **Q** = citers use a quantitative result a kinetic model could use (yield, amount change of named compounds,
  rate, activation energy, T/pH/a_w series, labelling fraction); **M** = a mechanism or pathway claim only (including a
  direction stated without numbers); **B** = background, identification, method, or a context that is about something
  else; **?** = no contexts.
- **Recommendation**: DOWNLOAD = Q and inside the model's lanes (trunk, sulfur, lipid, methionine); MAYBE = Q outside
  the aqueous model or with weak evidence, or few contexts but a quantitative abstract; SKIP = everything else. A Q
  paper with no sulfur and only ARP-preparation yields (vacuum dehydration, spray drying, microwave) is SKIP: those are
  process-engineering numbers, not kinetics of a step the engine has.

Counts over the 81 remaining: Q 13, M 26, B 30, ? 12.

## Caveats

- **Context coverage is thin.** Only 663 of the 3429 citing papers S2 lists for the 80 queried candidates carry any
  context (19 %; median 0.29 contexts per citer), and the most informative citers (our own corpus) are sparse: 19 of 80
  candidates have a corpus citer with a context, often the same 2024 review ("Application of Animal Resources into the
  Maillard Reaction Model ...").
- **S2 mislinks contexts.** For Zhou 2024 (10.1021/acs.jafc.4c00021), Wang 2008 (both), Ze 2019, Yu 2011 and Liu 2015 a
  large share of the contexts appear to be about a different paper by the same first author and year; these were discounted.
- **A citer's use is not the paper's content.** Citers of Chinese-lab ARP papers quote yields; citers of the 1970s-1990s
  JAFC papers quote mechanisms, even where the paper tabulates amounts. Several M rows (Boelens 1974, Elmore 2000,
  Smuda 2013) are known to hold numbers; they are SKIP here only because nobody cites them for numbers.
- **Contradiction or non-reproduction**: only one, second-hand. A 2019 German dissertation (Bildung flüchtiger
  Verbindungen und Carbonsäuren in Modellsystemen) states that the formation route it discusses, cited to work
  including Yaylayan 1999 and Pfeifer 2010, "could not be confirmed in model reactions at 50 °C by Smuda and Glomb".
  One citer of Pripis-Nicolau 2000 notes that the Strecker-aldehyde concentrations "were not reported in model systems".
  No context disputes any Q candidate.
- **Data hygiene found on the way**: `data/articles/qi2022.pdf` is not Qi et al. 2022 (Food Chem 10.1016/j.foodchem.2021.131111)
  but a JSTOR literary article (Cacciari, Revista de Letras 2005); `gao2020.pdf` is a different Gao 2020 (pea protein),
  so the tilapia paper is not on disk.

## Download list (ranked)

1. **Xu 2011** JSFA, lard x cysteine-xylose, [10.1002/jsfa.4445](https://doi.org/10.1002/jsfa.4445), `xu2011.pdf`: citers quote
   measured drops of MFT, FFT and the methylthiophenes with lard; a direct lipid x sulfur-lane test in the core pot.
2. **de Roos 2005** ACS Symp. Ser. 905, cysteine Amadori compounds, [10.1021/bk-2005-0905.ch009](https://doi.org/10.1021/bk-2005-0905.ch009),
   `deroos2005.pdf`: cited for quantified thiols from cysteine + carbonyls and for dilute- vs concentrated-medium routes
   (the cysteine-ARP step the sulfur lane lacks). Book chapter: may need library access.
3. **Yaghmur 2002** JAFC, furfural-cysteine in microemulsions, [10.1021/jf011158l](https://doi.org/10.1021/jf011158l),
   `yaghmur2002.pdf`: FFT yield vs a water reference; the earlier paper of the group whose 2005 paper sets
   `YAGHMUR_FFT_SHARE_CEILING`.
4. **Mandin 1999** JAFC, potato-like models, [10.1021/jf981277+](https://doi.org/10.1021/jf981277+), `mandin1999.pdf`:
   linoleic acid raising DMDS/DMTS 13-19x and methional formation; a methionine-lane check.
5. **Gao 2020** Food Chem, tilapia hydrolysate + xylose + cysteine RSM, [10.1016/j.foodchem.2020.127353](https://doi.org/10.1016/j.foodchem.2020.127353),
   `gao2020b.pdf`: FFT rising with temperature 150-190 °C; a directional matrix point, weakest of the five.

MAYBE, in order: Yang 2015 (`yang2015.pdf`, hexose twin of Xu 2011), Cannon 2023 (`cannon2023.pdf`, no contexts but an
abstract that says quantified), Zou 2019 (`zou2019.pdf`, free vs peptide-bound cysteine as H2S donor), Fanun 2001
(`fanun2001.pdf`, microemulsion rate), Gasior 2021 (`gasior2021.pdf`, FFT by SIDA in roasted goose), Du 2024 LWT
(`du2024b.pdf`, 3 contexts only), Qi 2022 (replace the wrong `qi2022.pdf`).

## Q: cited for a quantitative result (13)

| first author, year | title | DOI | source list | contexts / citers | what citers use it for | representative context (trimmed, ... = cut; <=25 words) | recommendation |
|---|---|---|---|---|---|---|---|
| Xu 2011 | Effects of lard on the formation of volatiles from the Maillard reaction of cysteine with xylose | [10.1002/jsfa.4445](https://doi.org/10.1002/jsfa.4445) | R1-A | 23 / 51 | lard lowers MFT, FFT, 2-/3-methylthiophene amounts in cysteine-xylose and raises lipid-oxidation volatiles (lipid x sulfur-lane test) | "presence of lard in cysteine and xylose model system reduced the formation of sulfur-containing compounds such as 2-methyl-3-furanthiol, 2-furanmethanethiol" | DOWNLOAD #1 (`xu2011.pdf`) |
| de Roos 2005 | Amadori Compounds of Cysteine and Their Role in the Development of Meat Flavor | [10.1021/bk-2005-0905.ch009](https://doi.org/10.1021/bk-2005-0905.ch009) | R1-B | 5 / 14 | quantified high-impact thiols from cysteine + carbonyls; cysteine-ARP vs cysteine + other-ARP routes in dilute vs concentrated media | "de Roos et al. quantified several high impact thiols from the reaction of cysteine with different carbonyl compounds" | DOWNLOAD #2 (`deroos2005.pdf`) |
| Yaghmur 2002 | Furfural−Cysteine Model Reaction in Food Grade Nonionic Oil/Water Microemulsions for Selective Flavor Formation | [10.1021/jf011158l](https://doi.org/10.1021/jf011158l) | R1-A | 8 / 52 | FFT yield from furfural + cysteine, microemulsion vs water reference; lineage of the corpus Yaghmur 2005 FFT-share ceiling | "Furfurylthiol ... was produced at much higher concentrations in w/o microemulsions than in a water reference system" | DOWNLOAD #3 (`yaghmur2002.pdf`) |
| Mandin 1999 | Volatile Compounds from Potato-like Model Systems | [10.1021/jf981277+](https://doi.org/10.1021/jf981277+) | R2-A | 11 / 23 | linoleic acid raises DMDS/DMTS 13-19x and accelerates methionine Strecker (methional) in potato models; methionine lane | "the addition of linoleic acid resulted in unexpected higher levels (13-19 folds) of dimethyl sulfides" | DOWNLOAD #4 (`mandin1999.pdf`) |
| Gao 2020 | Optimization of the Maillard reaction of xylose with cysteine ... fermented tilapia fish head hydrolysate (RSM) | [10.1016/j.foodchem.2020.127353](https://doi.org/10.1016/j.foodchem.2020.127353) | Sun2022-lead | 41 / 66 | FFT rises with temperature 150-190 C in a xylose-cysteine-hydrolysate pot; FFT and MFT-MFT as main meaty markers (directional matrix data) | "tilapia fish head hydrolysate-based MRPs contained more 2-furfurylthiol as the heating temperature increased from 150 to 190 C" | DOWNLOAD #5 (`gao2020b.pdf`) |
| Fanun 2001 | Sucrose ester microemulsions as microreactors for model Maillard reaction | [10.1016/s0927-7757(01)00786-5](https://doi.org/10.1016/s0927-7757(01)00786-5) | R2-B | 7 / 39 | rate enhancement and yields of furfural-cysteine products in sugar-ester W/O microemulsions; pseudo-first-order kinetics (microemulsion, not aqueous) | "The interface enhanced the reaction rate, increased the yields of key odorants, and formed new volatile compounds" | MAYBE (`fanun2001.pdf`) |
| Yang 2015 | Aromatic effect of fat and oxidized fat on a meat‐like model reaction system of cysteine and glucose | [10.1002/ffj.3248](https://doi.org/10.1002/ffj.3248) | R1-A | 14 / 49 | chicken fat suppresses S-volatiles in cysteine-glucose, raises alkyl heterocycles and aldehydes (hexose twin of Xu 2011) | "the inclusion of chicken fat in a cysteine-glucose reaction system suppressed the formation of sulfur-containing VOCs" | MAYBE (`yang2015.pdf`) |
| Cui 2017 | Controlled formation of flavor compounds by preparation and application of Maillard reaction intermediate (MRI) derived from xylose and phenylalanine | [10.1039/c7ra09355a](https://doi.org/10.1039/c7ra09355a) | R1-B | 30 / 65 | Phe-xylose ARP gives 8.77-fold more volatiles than free precursors (no sulfur) | "ARPs produced 8.77-fold higher volatile compound concentrations than their precursor mixture" | SKIP |
| Cui 2018 | Improved controlled flavor formation during heat-treatment with a stable Maillard reaction intermediate derived from xylose-phenylalanine | [10.1016/j.foodchem.2018.07.161](https://doi.org/10.1016/j.foodchem.2018.07.161) | R1-B | 32 / 99 | Phe-xylose ARP preparation yield and storage stability (no sulfur) | "Phenylalanine-xylose-ARP only reduced 6.49% after 2 months storage at room temperature" | SKIP |
| Cui 2018 | Synergistic Effect of a Thermal Reaction and Vacuum Dehydration on Improving Xylose–Phenylalanine Conversion to N-(1-Deoxy- d -xylulos-1-yl)-phenylalanine during an Aqueous Maillard Reaction | [10.1021/acs.jafc.8b04448](https://doi.org/10.1021/acs.jafc.8b04448) | R1-B | 8 / 48 | ARP yield raised from 13.62 % by vacuum dehydration (Phe-xylose; no sulfur) | "the yield of Amadori compounds increased from 13.62% t[o] ..." | SKIP |
| Wang 2024 | Rapid preparation of the Amadori rearrangement product of glutamic acid - xylose through intermittent microwave heating and its browning formation potential in microwave thermal processing | [10.1016/j.foodres.2024.114075](https://doi.org/10.1016/j.foodres.2024.114075) | R2-B | 7 / 11 | Glu-xylose ARP 58 % by intermittent microwave heating (no sulfur) | "Wang et al. (2024) utilized intermittent microwave heating to yield Xylose-Glu with the productivity of 58.09%" | SKIP |
| Xu 2019 | Preparation of N ‐(1‐Deoxy‐Α‐D‐Xylulos‐1‐Yl)‐Glutamic Acid via Aqueous Maillard Reaction Coupled with Vacuum Dehydration and Its Flavor Formation Through Thermal Treatment of Baking Process | [10.1111/1750-3841.14733](https://doi.org/10.1111/1750-3841.14733) | R2-B | 13 / 45 | Glu-xylose ARP yield 1.85 % aqueous vs 75 % with vacuum dehydration (no sulfur) | "Xu et al. (2019) produced Xylose-Glu in the aqueous system, with the productivity of merely 1.85%" | SKIP |
| Zhang 2021 | Accelerated Dissipation of Free and Immobilized Water Facilitating the Intramolecular Dehydration of N-Xylosamine and Conversion Improvement of the Amadori Rearrangement Product of Aspartic Acid–Xylose Reaction | [10.1021/acs.jafc.1c05827](https://doi.org/10.1021/acs.jafc.1c05827) | R1-B | 17 / 20 | Asp-xylose ARP yields 21.5-30.5 % vs spray-dryer inlet temperature, 69-78 % with vacuum dehydration (process engineering, no sulfur) | "the output of Xylose-Asp was 21.5% at 150 C, while the output of Xylose-Asp was 30.5% at 180 C" | SKIP |

## M: cited for a mechanism or pathway claim only (26)

| first author, year | title | DOI | source list | contexts / citers | what citers use it for | representative context (trimmed, ... = cut; <=25 words) | recommendation |
|---|---|---|---|---|---|---|---|
| Zou 2019 | Flavour precursor peptide from an enzymatic beef hydrolysate Maillard reaction-II: Mechanism of the synthesis of flavour compounds from a sulphur-containing peptide through a Maillard reaction | [10.1016/j.lwt.2019.04.022](https://doi.org/10.1016/j.lwt.2019.04.022) | R2-B,Sun2022-lead | 17 / 46 | free cysteine releases H2S more easily than peptide-bound; Cys-Gly-Val as MFT precursor | "cysteine releases H2S more easily when present as a free amino acid rather than as a peptide form" | MAYBE (`zou2019.pdf`) |
| Boelens 1974 | Organic sulfur compounds from fatty aldehydes, hydrogen sulfide, thiols, and ammonia as flavor constituents | [10.1021/jf60196a036](https://doi.org/10.1021/jf60196a036) | R2-B | 6 / 73 | alkyltrithiolanes and S-heterocycles from aldehydes + H2S/NH3; cysteine -> H2S, NH3, ethanal | "the degradation of cysteine with water which can easily lead to H2S, NH3, CO2 and ethanal" | SKIP |
| Chen 2023 | Exogenous Threonine-Induced Conversion of Threonine-Xylose Amadori Compound to Heyns Compound for Efficiently Promoting the Formation of Pyrazines | [10.1021/acs.jafc.3c02311](https://doi.org/10.1021/acs.jafc.3c02311) | R2-B | 3 / 15 | exogenous Thr converts Thr-xylose ARP to Heyns, promotes pyrazines (no sulfur) | "exogenous Thr could be engaged in the decomposition pathway of Xylose-Thr Amadori compound, and facilitate the production of pyrazines" | SKIP |
| Elmore 2000 | Formation of 2-Alkyl-(2H)-thiapyrans and 2-Alkylthiophenes in Cooked Beef and Lamb | [10.1021/jf000124v](https://doi.org/10.1021/jf000124v) | R2-B | 10 / 34 | 2,4-dienals + H2S -> 2-alkylthiophenes/thiapyrans in cooked meat | "(E,E)-2,4-heptadienal or (E,E)-2,4-nonadienal can react with H2S to produce 2-pentylthiophene" | SKIP |
| Eric 2014 | Temperature and cysteine addition effect on formation of sunflower hydrolysate Maillard reaction products | [10.1016/j.foodres.2014.01.030](https://doi.org/10.1016/j.foodres.2014.01.030) | Sun2022-lead | 35 / 86 | more cysteine -> more sulfur odour, less meaty; pH falls during heating (sunflower hydrolysate, sensory) | "lower cysteine addition leads to less pronounced meaty taste; while, its higher levels increase the sulfur odor" | SKIP |
| Guentert 1990 | Identification and formation of some selected sulfur-containing flavor compounds in various meat model systems | [10.1021/jf00101a007](https://doi.org/10.1021/jf00101a007) | R1-B | 21 / 71 | thiamine routes to thiazoles; tetrahydrothiophen-3-one identification | "Breakdown of thiamine to produce substituted thiazoles (after Guentert et al., 1990)" | SKIP |
| Huang 2012 | Inhibiting the color formation by gradient temperature‐elevating Maillard reaction of soybean peptide‐xylose system based on interaction of l ‐cysteine and Amadori compounds | [10.1002/psc.2406](https://doi.org/10.1002/psc.2406) | R1-B | 12 / 43 | cysteine inhibits browning only if added early (interaction with intermediate) | "the effective inhibition only occurred when cysteine was added to the system at the early stage at a specific time" | SKIP |
| Liu 2015 | Effect of thermal treatment on the flavor generation from Maillard reaction of xylose and chicken peptide | [10.1016/j.lwt.2015.05.061](https://doi.org/10.1016/j.lwt.2015.05.061) | Sun2022-lead | 87 / 193 | temperature windows: ketones below 80 C, N-heterocycles need high T (peptide-xylose); some contexts belong to other Liu papers | "Ketones, which mainly derive from sugar degradation, can be formed at temperatures below 80 C" | SKIP |
| Pfeifer 2010 | Investigation of Reactive α-Dicarbonyl Compounds Generated from the Maillard Reactions of l -Methionine with Reducing Sugars via Their Stable Quinoxaline Derivatives | [10.1021/jf1008988](https://doi.org/10.1021/jf1008988) | R2-B | 15 / 43 | methional from methionine + alpha-dicarbonyls; same non-confirmation note as Yaylayan 1999 | "The established pathway for the formation of methional involves the reaction of methionine with a-dicarbonyl compounds" | SKIP |
| Pripis‐Nicolau 2000 | Formation of Flavor Components by the Reaction of Amino Acid and Carbonyl Compounds in Mild Conditions | [10.1021/jf991024w](https://doi.org/10.1021/jf991024w) | R2-B | 77 / 200 | Strecker aldehydes and S-heterocycles from dicarbonyls + amino acids at wine conditions; one citer notes the aldehyde amounts were not reported | "the concentrations of these Strecker aldehydes were not reported in model systems" | SKIP |
| Shibamoto 1977 | Formation of sulfur- and nitrogen-containing compounds from the reaction of furfural with hydrogen sulfide and ammonia | [10.1021/jf60209a055](https://doi.org/10.1021/jf60209a055) | R1-B | 6 / 50 | furfural + H2S -> FFT and thiophenes | "2-furfurylthiol contributing to the chicken aroma was probably formed by the reaction of H2S with ... 2-furfural" | SKIP |
| Smuda 2013 | Fragmentation Pathways during Maillard-Induced Carbohydrate Degradation | [10.1021/jf305117s](https://doi.org/10.1021/jf305117s) | R2-B | 31 / 115 | carbohydrate fragmentation to GO/MGO/DA (pathway) | "the carbohydrate cleavage yields a spectrum of highly reactive carbonyls, including glyoxal and methylglyoxal" | SKIP |
| Totlani 2005 | Reactivity of Epicatechin in Aqueous Glycine and Glucose Maillard Reaction Models: Quenching of C2, C3, and C4 Sugar Fragments | [10.1021/jf050044x](https://doi.org/10.1021/jf050044x) | R2-B | 31 / 135 | epicatechin traps C2-C4 sugar fragments | "epi-catechin ... was in 2005 demonstrated by Totlani and Peterson to inhibit Maillard's reaction by trapping a-dicarbonyls" | SKIP |
| Van Lancker 2010 | Mechanistic Insights into Furan Formation in Maillard Model Systems | [10.1021/jf102929u](https://doi.org/10.1021/jf102929u) | R2-B | 37 / 81 | furan formation mechanisms (labelling) | "the presence of lysine inhibited the formation of furan compounds in model Maillard systems when heated to 180 C" | SKIP |
| Wang 2020 | Reinvestigation of 2-acetylthiazole formation pathways in the Maillard reaction | [10.1016/j.foodchem.2020.128761](https://doi.org/10.1016/j.foodchem.2020.128761) | R1-B | 11 / 39 | 2-acetylthiazole from GO/MGO + H2S + NH3 (labelling) | "2-acetylthiazole can derive from GO and MGO from glucose reacting with H2S and NH3 via Strecker degradation of cysteine" | SKIP |
| Wang 2022 | Controlled Formation of Pyrazines: Inhibition by Ellagic Acid Interaction with N-(1-Deoxy- d -xylulos-1-yl)-glycine and Promotion through Ellagic Acid Oxidation | [10.1021/acs.jafc.1c07391](https://doi.org/10.1021/acs.jafc.1c07391) | R2-B | 4 / 22 | pyrazine inhibition depends on ellagic-acid addition time (no sulfur) | "the inhibition effect on pyrazine relayed on the addition time of ellagic acid" | SKIP |
| Weenen 1998 | Reactive intermediates and carbohydrate fragmentation in Maillard chemistry | [10.1016/s0308-8146(98)00074-0](https://doi.org/10.1016/s0308-8146(98)00074-0) | R2-B | 103 / 156 | MGO from retro-aldol of 3-DG; GO not from ARP/deoxyosones | "Weenen suggested that glyoxal formation could not be explained by only degradation of Amadori product or deoxyosones" | SKIP |
| Xia 2022 | Structural diversity and concentration dependence of pyrazine formation: Exogenous amino substrates and reaction parameters during thermal processing of l-alanyl-l-glutamine Amadori compound | [10.1016/j.foodchem.2022.133144](https://doi.org/10.1016/j.foodchem.2022.133144) | R2-B | 10 / 31 | pyrazine formation depends on T, time, concentration in Ala-Gln ARP (no sulfur) | "the formation of pyrazine compounds was affected by the reaction temperature, reaction time, and reactant concentration" | SKIP |
| Xia 2023 | Simultaneously Enhanced Formation of Pyrazines and Furans during Thermal Degradation of the Glycyl- l -glutamine Amadori Compound by Selected Exogenous Amino Acids and Appropriate Elevated Temperatures | [10.1021/acs.jafc.3c00085](https://doi.org/10.1021/acs.jafc.3c00085) | R2-B | 6 / 21 | ARP degradation promotes pyrazines and furans (no sulfur) | "the thermal degradation of Amadori compounds promotes the formation of pyrazines and furans" | SKIP |
| Yaylayan 1999 | Origin of 2,3-Pentanedione and 2,3-Butanedione in d -Glucose/ l -Alanine Maillard Model Systems | [10.1021/jf9902292](https://doi.org/10.1021/jf9902292) | R2-B | 21 / 100 | 13C origin of 2,3-butanedione/pentanedione; a 2019 thesis says the route was not confirmed at 50 C by Smuda & Glomb | "2,3-butanedione can be generated ... via a single pathway involving only glucose carbon atoms" | SKIP |
| Yaylayan 2000 | Origin of Carbohydrate Degradation Products in l -Alanine/ d -[ 13 C]Glucose Model Systems | [10.1021/jf000004n](https://doi.org/10.1021/jf000004n) | R2-B | 27 / 95 | 13C fragmentation of glucose/alanine (hydroxyacetone C4-C6, GO) | "dehydration from C3-C4 and C5-C6 of glucose molecules followed by retro-aldol cleavage can also generate GO" | SKIP |
| Yaylayan 2001 | Carbohydrate and Amino Acid Degradation Pathways in l -Methionine/ d -[13C] Glucose Model Systems | [10.1021/jf000986w](https://doi.org/10.1021/jf000986w) | R2-B | 14 / 44 | methionine/glucose degradation pathways; methional | "SD of methionine leads to the generation of methional, which smells like cooked potato" | SKIP |
| Yu 2011 | Mechanism of formation of sulphur aroma compounds from l-ascorbic acid and l-cysteine during the Maillard reaction | [10.1016/j.foodchem.2011.11.111](https://doi.org/10.1016/j.foodchem.2011.11.111) | R1-B | 33 / 85 | S-S compounds and trithiolanes from ascorbic acid + cysteine; many contexts actually belong to Yu 2012 browning papers | "Maillard-type reactions with a-dicarbonyls ... between ascorbic acid degradation products and Cys, have been reported to yield S-S compounds" | SKIP |
| Yu 2019 | Effective Mechanism of (−)-Epigallocatechin Gallate Indicating the Critical Formation Conditions of Amadori Compound during an Aqueous Maillard Reaction | [10.1021/acs.jafc.9b00034](https://doi.org/10.1021/acs.jafc.9b00034) | R2-B | 5 / 40 | EGCG traps deoxypentosones, shifts Ala-xylose ARP formation (no sulfur) | "the activation energy of ARP formation was decreased by EGCG, and the reaction equilibrium was pushed to the direction of ARP generation" | SKIP |
| Zhang 2023 | Controlled Selective Formation of Amadori Compounds from α/ε Mono- or Di-glycation of Lysine with Xylose | [10.1021/acs.jafc.3c00388](https://doi.org/10.1021/acs.jafc.3c00388) | R2-B | 3 / 8 | mono- vs di-glycated lysine ARPs (no sulfur) | "Prolonged heating time and high temperature facilitated the production of diXylose-a,e-Lys" | SKIP |
| Zou 2018 | Discovery of Amadori‐Type Conjugates in a Peptide Maillard Reaction and Their Corresponding Influence on the Formation of Pyrazines | [10.1111/1750-3841.14156](https://doi.org/10.1111/1750-3841.14156) | R1-B | 5 / 42 | dipeptide ARPs give more pyrazines than amino-acid ARPs (no sulfur) | "dipeptide-derived Amadori compound was degraded to yield more pyrazines, compared to amino acid-derived" | SKIP |

## B: background, identification or method only (30)

| first author, year | title | DOI | source list | contexts / citers | what citers use it for | representative context (trimmed, ... = cut; <=25 words) | recommendation |
|---|---|---|---|---|---|---|---|
| Du 2024 | Effect of 2-alkenals on the Maillard reaction of cysteine-glucose: Initial stage intermediate formation and flavor compound generation pathways | [10.1016/j.lwt.2024.115741](https://doi.org/10.1016/j.lwt.2024.115741) | R1-A | 3 / 20 | only 3 contexts: general temperature windows and aldehyde origins | "thiophenes and thiazoles ... are favored at high temperatures (>140 C)" | MAYBE (`du2024b.pdf`) |
| Gasior 2021 | Key Aroma Compounds in Roasted White Kołuda Goose | [10.1021/acs.jafc.1c01475](https://doi.org/10.1021/acs.jafc.1c01475) | Sohail-lead | 24 / 49 | SIDA/SAFE method and FFT occurrence in roasted goose; key-odorant identification | "the flavor contribution of FFT to foodstuffs has been demonstrated in ... roasted goose/duck" | MAYBE (`gasior2021.pdf`) |
| Qi 2022 | Key aroma compounds in stewed mutton (goat meat) with thyme | [10.1016/j.foodchem.2021.131111](https://doi.org/10.1016/j.foodchem.2021.131111) | Sohail-lead | 38 / 80 | lipid-derived odorants in stewed goat meat; sensory method. NB data/articles/qi2022.pdf is a wrong file (a JSTOR literary article) | "(E)-2-Octenal is recognized as an important contributor to meaty aroma and can be generated from the oxidation of linoleic acid" | MAYBE (`qi2022b.pdf (or replace the wrong qi2022.pdf)`) |
| Alim 2019 | Research of beef-meaty aroma compounds from yeast extract using carbon module labeling (CAMOLA) technique | [10.1016/j.lwt.2019.06.006](https://doi.org/10.1016/j.lwt.2019.06.006) | R1-B | 15 / 29 | yeast-extract meat-flavour background; retention-index method | "the characteristics of meat flavor are mainly derived from sulfur-containing nitrogen compounds" | SKIP |
| Brehm 2019 | Novel Taste-Enhancing 4-Amino-2-methyl-5-heteroalkypyrimidines Formed from Thiamine by Maillard-Type Reactions | [10.1021/acs.jafc.9b05896](https://doi.org/10.1021/acs.jafc.9b05896) | R2-B | 7 / 20 | thiamine-derived taste compounds, kokumi | "S-((4-amino-2methylpyrimidin-5-yl)methyl)-L-cysteine ... functions as a kokumi substance" | SKIP |
| Buttery 1984 | Thiamin odor and bis(2-methyl-3-furyl) disulfide | [10.1021/jf00123a061](https://doi.org/10.1021/jf00123a061) | R2-B | 7 / 54 | odour threshold of bis(2-methyl-3-furyl) disulfide (0.02 ng/kg); identification | "The odour threshold value of this disulfide has been reported as 0.02 ng/kg" | SKIP |
| Chen 2018 | Effect of temperature on flavor compounds ... Maillard reaction products derived from mushroom hydrolysate | [10.3390/molecules23020247](https://doi.org/10.3390/molecules23020247) | Sun2022-lead | 46 / 76 | umami amino acids and 5'-GMP in mushroom MRPs | "The concentration of free amino acids and 5'-GMP ... under 125 C were 3 to 4 times higher than those in the mushroom hydrolysate" | SKIP |
| Deng 2024 | Promoted formation of pyrazines by targeted precursor addition to improve aroma of thermally processed methionine-glucose Amadori compound | [10.1016/j.foodchem.2024.142033](https://doi.org/10.1016/j.foodchem.2024.142033) | R2-B | 8 / 16 | pyrazines/meat-note background | "Most pyrazines are derived from the Strecker reaction between amino acids and a-carboxylic compounds" | SKIP |
| Feng 2024 | Directed preparation of N-(1-deoxy-d-ribulos-1-yl)-glutathione during short-term high temperature dehydration: Suppressing the side reactions of glutathione | [10.1016/j.fbio.2024.104468](https://doi.org/10.1016/j.fbio.2024.104468) | R1-B | 1 / 2 | one uninformative context | "This was consistent with the results" | SKIP |
| Jerić 2009 | Screening for glucose‐triggered modifications of glutathione | [10.1002/psc.1159](https://doi.org/10.1002/psc.1159) | R2-B | 1 / 15 | GSH-glucose glycation products reported | "a direct reaction with GSH to form pyroglutamic acid related Amadori compounds ... have been reported" | SKIP |
| Lee 2011 | Determination of Volatile Sulfur Compounds Formed by the Maillard Reaction of Glutathione with Glucose | [10.1021/bk-2011-1068.ch011](https://doi.org/10.1021/bk-2011-1068.ch011) | R1-A | 1 / 3 | one context: thiophene pathway mention | "the formation pathways of 2-ethylthiophene, 2,5-dimethylthiophene ... were explained in the glutathione and [13C6]glucose system" | SKIP |
| Lei 2024 | Analysis of volatiles from the thermal decomposition of Amadori rearrangement products in the cysteine-glucose Maillard reaction and density functional theory study | [10.1016/j.foodres.2024.114454](https://doi.org/10.1016/j.foodres.2024.114454) | R1-B | 2 / 17 | ARP preparation method | "The general method for the preparation of ARPs is carried out in organic solvents such as methanol and glycerol" | SKIP |
| Li 2021 | The Maillard Reaction as Source of Meat Flavor Compounds in Dry Cured Meat Model Systems under Mild Temperature Conditions | [10.3390/molecules26010223](https://doi.org/10.3390/molecules26010223) | R1-A | 24 / 54 | dry-cure processing conditions, pH/temperature background (cited mostly by meat-quality papers) | "Previous studies found that the volatile generation from Maillard reaction is not only dependent on the concentration of the substrates but also ... pH" | SKIP |
| Liu 2023 | Effect of Aliphatic Aldehydes on Flavor Formation in Glutathione–Ribose Maillard Reactions | [10.3390/foods12010217](https://doi.org/10.3390/foods12010217) | R1-A | 8 / 32 | aldehydes from lipid oxidation/Strecker; Maillard-lipid interaction in general | "Aldehydes are derived from lipid oxidation during thermal treatment or Strecker degradation" | SKIP |
| Ma 2020 | Effects of different enzymatic hydrolyses of mutton tallow on the aroma of the xylose-cysteine Maillard reaction | [10.1002/ejlt.201900212](https://doi.org/10.1002/ejlt.201900212) | Sun2022-lead | 11 / 27 | lipid-oxidation aldehydes, enzyme choice for tallow | "excessive lipid oxidation should be avoided as it produces high content of lipid oxidation-derived VOCs" | SKIP |
| Nikolov 2012 | Role of the Ribose-Specific Marker Furfuryl-amine in the Formation of Aroma Active 1-(Furan-2-ylmethyl)-1H-pyrrole (or Furfuryl-pyrrole) Derivatives | [10.1021/jf3030002](https://doi.org/10.1021/jf3030002) | R1-B | 1 / 17 | one context: furfurylpyrrole occurrence | "N-furfurylpyrrole is found in various foods, such as coffee, chocolate, popcorn, and roasted chicken" | SKIP |
| Ota 2006 | Characterization of a New Maillard Type Reaction Product Generated by Heating 1-Deoxymaltulosyl-glycine in the Presence of Cysteine | [10.1021/jf060296d](https://doi.org/10.1021/jf060296d) | R1-B | 2 / 16 | NMR/MS identification of a 1-deoxyosone-cysteine product | "the production of 1-deoxyosone intermediate was confirmed ... in the presence of cysteine" | SKIP |
| Schieberle 2005 | The Carbon Module Labeling (CAMOLA) Technique: A Useful Tool for Identifying Transient Intermediates in the Formation of Maillard‐Type Target Molecules | [10.1196/annals.1333.029](https://doi.org/10.1196/annals.1333.029) | R1-B | 19 / 71 | CAMOLA method | "an advanced tool to quantify the relative contribution of carbohydrate fragments" | SKIP |
| Silwar 1989 | Gaschromatographisch-massenspektrometrische Untersuchung der im Modellsystem Cystein/Methionin/Furfural unter Röstbedingungen gebildeten Aromastoffe | [10.1007/bf01028064](https://doi.org/10.1007/bf01028064) | R2-B | 2 / 17 | furfural as precursor; occurrence | "Furfural was probably important as a precursor in the formation of other aroma compounds" | SKIP |
| Tuong 2003 | Generation of Thiols by Biotransformation of Cysteine−Aldehyde Conjugates with Baker's Yeast | [10.1021/jf026198j](https://doi.org/10.1021/jf026198j) | R2-B | 8 / 30 | yeast biotransformation of cysteine conjugates to thiols | "biotransformation products of the cysteine-aldehyde conjugate by baker's yeast" | SKIP |
| Wang 2008 | Amino acid-dependent formation pathways of 2-acetylfuran and 2,5-dimethyl-4-hydroxy-3[2H]-furanone in the Maillard reaction | [10.1016/j.foodchem.2008.12.014](https://doi.org/10.1016/j.foodchem.2008.12.014) | R2-B | 28 / 80 | contexts mostly mislinked to other Wang 2009 papers | "high temperature and long heating times were beneficial to the Maillard reaction of DMTS" | SKIP |
| Wang 2008 | Comparison of 2-Acetylfuran Formation between Ribose and Glucose in the Maillard Reaction | [10.1021/jf802683a](https://doi.org/10.1021/jf802683a) | R2-B | 6 / 19 | contexts all from one ghee thesis citing it as a method reference | "When the concentration of reactants was kept constant across aqueous, emulsion and cubic systems, the rates of product formation were found to be different" | SKIP |
| Wei 2020 | Aromatic effects of immobilized enzymatic oxidation of chicken fat on flaxseed-derived Maillard reaction products | [10.1016/j.foodchem.2019.125560](https://doi.org/10.1016/j.foodchem.2019.125560) | Sun2022-lead | 21 / 76 | oxidised-fat aldehydes as reactive carbonyls; model-system background | "lipid-derived aldehydes ... should be considered not only as odorants, but also as reactive carbonyl precursors" | SKIP |
| Wei 2022 | Glycine-Xylose Amadori Compound Formation Tracing through Maillard Browning Inhibition by 2-Threityl-thiazolidine-4-carboxylic Acid Formation from Deoxyosone and Exogenous Cysteine | [10.1021/acs.jafc.2c04961](https://doi.org/10.1021/acs.jafc.2c04961) | R1-A | 1 / 10 | NMR isomers of ARPs | "Similar isomers were found in many ACs derived from amino acids and monosaccharides by NMR" | SKIP |
| Yang 2012 | Response Surface Methodology for Meat‐like Odorants from Maillard Reaction with Glutathione I: the Optimization Analysis and the General Pathway Exploration | [10.1111/j.1750-3841.2012.02863.x](https://doi.org/10.1111/j.1750-3841.2012.02863.x) | R1-B | 4 / 25 | MFT disulfide origin unclear; RSM method | "The formation of 2-methyl-3-(methyldisulfanyl) furan is probably caused by the Maillard reaction with glutathione and xylose" | SKIP |
| Yang 2012 | The mechanism of peptide bonds cleavage and volatile compounds generated from pentapeptide to heptapeptide via Maillard reaction | [10.1016/j.foodchem.2012.01.044](https://doi.org/10.1016/j.foodchem.2012.01.044) | R1-B | 16 / 52 | peptides in Maillard: bond cleavage background | "polypeptides could be degraded into dipeptides and tripeptides and further form volatiles" | SKIP |
| Yang 2017 | Preparation of 1-Amino-1-deoxyfructose Derivatives by Stepwise Increase of Temperature in Aqueous Medium and Their Flavor Formation Compared with Maillard Reaction Products | [10.1007/s11947-017-2039-4](https://doi.org/10.1007/s11947-017-2039-4) | R2-B | 4 / 22 | ARP flavour profile resembles free precursors | "the flavor profile obtained from ARP is similar to that obtained from the corresponding amino acids and reducing sugars" | SKIP |
| Yuan 2024 | Mechanism of pyrazines and thioethers formation promoted by high oxygen concentration in the methionine‐glucose Maillard reaction system | [10.1002/jsfa.14094](https://doi.org/10.1002/jsfa.14094) | R2-B | 2 / 9 | 2 contexts, both off-topic (lysine as nucleophile) | "Free amino acids, especially lysine and arginine, act as the main nucleophiles initiating the MR" | SKIP |
| Ze 2019 | Studies on the Origin of Carbons in Aroma Compounds from [13C6]Glucose -Cysteine-(E)-2-Nonenal Model Reaction Systems | [10.3390/polym11030521](https://doi.org/10.3390/polym11030521) | R1-B | 2 / 12 | 2 contexts, both about ketones in meat (likely mislinked) | "Ketones ... serve a coordinating function in aroma development" | SKIP |
| Zhou 2024 | Exogenous Alanine Promoting the Reaction between Amadori Compound and Deoxyxylosone and Inhibiting the Formation of 2-Furfural during Thermal Treatment | [10.1021/acs.jafc.4c00021](https://doi.org/10.1021/acs.jafc.4c00021) | R1-A | 12 / 8 | contexts are about ARP bioactivity and belong to a different Zhou 2024 (S2 mislink) | "Fru-Gly ... suppress the generation of radicals induced by metal ions" | SKIP |

## ?: no citation contexts available (12)

| first author, year | title | DOI | source list | contexts / citers | what citers use it for | representative context (trimmed, ... = cut; <=25 words) | recommendation |
|---|---|---|---|---|---|---|---|
| Cannon 2023 | Quantitation of Key Odorants in Arginine/Cysteine-Glucose Maillard Reactions | [10.31665/jfb.2023.18338](https://doi.org/10.31665/jfb.2023.18338) | R1-A | 0 / 0 | no S2 citers; OpenAlex abstract: key odorants of Arg/Cys-glucose identified and quantified (GC-MS, GC-O) under several conditions | "(abstract) key odorants ... were identified and quantified by gas chromatography-mass spectrometry and gas chromatography-olfactometry" | MAYBE (`cannon2023.pdf`) |
| Chen 2026 | Exogenous xylose trapped threonine to regenerate Amadori compounds and accelerated pH decline impeding pyrazines formation upon thermal treatment of Amadori compounds revealed by sotope labeling | [10.1016/j.crfs.2026.101433](https://doi.org/10.1016/j.crfs.2026.101433) | R2-B | 0 / 0 | no citers (2026); no sulfur | - | SKIP |
| Deng 2022 | Mechanism of Pyrazine Formation Intervened by Oxidized Methionines during Thermal Degradation of the Methionine–Glucose Amadori Compound | [10.1021/acs.jafc.2c06458](https://doi.org/10.1021/acs.jafc.2c06458) | R2-B | 0 / 12 | 12 citers, none with contexts; Met-glucose ARP pyrazines | - | SKIP |
| Liu 2026 | Underlying substituted pyrazines formation in the glucose-alanine Maillard reaction: Carbon module labeling and carbonyl analysis strategy | [10.1016/j.foodres.2026.119409](https://doi.org/10.1016/j.foodres.2026.119409) | R2-B | 0 / 1 | 1 citer, no contexts; no sulfur | - | SKIP |
| Lu 2026 | Formation of characteristic aroma and inorganic fugacious products during pyrolysis of Maillard intermediates with different reaction degrees | [10.1016/j.jaap.2026.107660](https://doi.org/10.1016/j.jaap.2026.107660) | R2-B | 0 / 0 | no citers (2026); pyrolysis | - | SKIP |
| Ma 2025 | Preparation, characterization, and thermal pyrolysis of three Maillard intermediates | [10.1016/j.jaap.2025.107274](https://doi.org/10.1016/j.jaap.2025.107274) | R2-B | 0 / 5 | 5 citers, no contexts; pyrolysis | - | SKIP |
| Parker 2013 | Comparison of Ribose and Ascorbic Acid in Model Process Reactions Containing Cysteine | [10.1016/b978-0-12-398549-1.00040-4](https://doi.org/10.1016/b978-0-12-398549-1.00040-4) | R2-A | 0 / 1 | 1 citer, no contexts; OpenAlex has no abstract | - | SKIP |
| Schieberle 2003 | The carbohydrate modul labelling technique - a useful tool to clarify formation pathways of aroma compounds formed in Maillard-type reactions. | - | R2-B | - | no DOI, not queried (Weurman 2002 proceedings, CAMOLA) | - | SKIP |
| Xu 2014 | Thermal reaction of anisaldehyde in the presence of l-cysteine, a model reaction of Chinese stew meat flavor generation | [10.7282/t3mp51kh](https://doi.org/10.7282/t3mp51kh) | R2-B | 0 / 0 | thesis, no citers; abstract: anisaldehyde + cysteine at 60-150 C, 2 h, products identified | - | SKIP |
| Zhang 2024 | Control Formation of Furans and Pyrazines Resulting from Dual Glycation Sites in Nα,Nε-Di(1-deoxy- d -xylulos-1-yl)lysine via Elevating Thermal Degradation Temperatures | [10.1021/acs.jafc.4c07733](https://doi.org/10.1021/acs.jafc.4c07733) | R2-B | 0 / 2 | 2 citers, no contexts; no sulfur | - | SKIP |
| Zhao 2025 | Unraveling the formation pathways of key meat-flavor Spirocyclic sulfur compound from thiamin-cysteine model reaction: A combined experimental and computational study | [10.1016/j.foodchem.2025.147221](https://doi.org/10.1016/j.foodchem.2025.147221) | R1-B | 0 / 0 | no citers; OpenAlex has no abstract; title says combined experimental and computational (DFT part unusable under policy) | - | SKIP |
| Zhao 2026 | Unraveling the formation mechanism of pyrazines and their derivatives in the glucose-glycine Maillard reaction using carbon module labeling coupled with carbonyl analysis | [10.1016/j.crfs.2026.101552](https://doi.org/10.1016/j.crfs.2026.101552) | R2-B | 0 / 0 | no citers (2026); no sulfur | - | SKIP |

## Excluded: already on disk or in the registry (43)

`registry` = DOI in `data/keys/papers.yml`; `pdf:` = DOI printed on the PDF; `(title)` / `(name/title)` = matched by
file name and page-1 title words because the PDF prints no DOI.

| first author, year | DOI | source list | on disk as |
|---|---|---|---|
| Zhao 2018 | 10.1016/j.foodchem.2018.07.106 | R1-A,Sun2022-lead | registry |
| Zhang 2020 | 10.1016/j.foodres.2020.109918 | R1-A,Sun2022-lead | registry |
| Wang 2019 | 10.1016/j.foodchem.2019.125482 | R1-A,Sun2022-lead | registry |
| Zhou 2023 | 10.1021/acs.jafc.3c04424 | R1-A | registry |
| Zhou 2022 | 10.1021/acs.jafc.2c07026 | R1-A | registry |
| Du 2024 | 10.1016/j.crfs.2024.100719 | R1-A | registry |
| Elmore 2002 | 10.1021/jf0108718 | R1-A | registry |
| Madruga 1995 | 10.1002/jsfa.2740680308 | R1-A | registry |
| Shu 1985 | 10.1021/jf00063a029 | R1-A | registry |
| Xu 2008 | 10.1016/j.foodres.2008.05.005 | R1-A,Sun2022-lead | registry |
| Sun 2019 | 10.1111/1750-3841.14911 | R1-A,Sun2022-lead | registry |
| Liu 2022 | 10.1021/acs.jafc.2c04874 | R1-A | registry |
| Xu 2007 | 10.1002/jsfa.3093 | R1-A,Sun2022-lead | pdf:xu2007.pdf |
| Hartman 1984 | 10.1021/jf00125a016 | R1-A | pdf:hartman1984.pdf (name/title) |
| Mottram 1995 | 10.1021/jf00052a027 | R1-A | pdf:mottram1995b.pdf (title) |
| Hincelin 1992 | 10.1016/0308-8146(92)90272-4 | R1-A | pdf:hincelin1992.pdf (title) |
| Zhang 2023 | 10.1016/j.foodres.2023.113200 | R1-A | pdf:zhang2023b.pdf |
| Tang 2019 | 10.1021/acs.jafc.9b04694 | R1-A | pdf:tang2019.pdf |
| Hartman 1984 | 10.1111/j.1365-2621.1984.tb12480.x | R1-A | pdf:hartman1984b.pdf (name/title) |
| Hofmann 1997 | 10.1021/jf960456t | R1-B | registry |
| Zhang 1989 | 10.1021/jf00088a044 | R2-A | pdf:zhang1989.pdf (title) |
| Chen 2002 | 10.1021/jf011244l | R2-A | pdf:chen2002.pdf |
| Yu 2010 | 10.1016/j.foodchem.2010.01.049 | R2-A | pdf:yu2010.pdf |
| Zhang 1988 | 10.1021/jf00083a022 | R2-A | pdf:zhang1988.pdf (title) |
| Chen 2000 | 10.1021/jf991076l | R2-A | pdf:chen2000.pdf |
| Hofmann 1998 | 10.1007/s002170050324 | R2-A | pdf:hofmann1998c.pdf (title) |
| Lu 2019 | 10.1021/acs.jafc.9b02053 | R2-A | pdf:lu2019.pdf |
| Zhang 2022 | 10.1016/j.foodres.2022.111081 | R2-A | pdf:zhang2022.pdf |
| Dreher 2003 | 10.1021/jf034023j | R2-A | pdf:dreher2003.pdf |
| El-Massry 2003 | 10.1007/s00726-002-0318-4 | R2-A | pdf:el-massry2003.pdf |
| Zhang 1991 | 10.1021/jf00006a031 | R2-A | pdf:zhang1991.pdf (title) |
| Farmer 1992 | 10.1002/jsfa.2740600414 | R2-A | pdf:farmer1992.pdf (title) |
| Mottram 1995 | 10.1021/jf00053a033 | R2-A | pdf:mottram1995c.pdf (title) |
| Zhang 1991 | 10.1021/jf00004a029 | R2-A,Sun2022-lead | pdf:zhang1991b.pdf (title) |
| Schieberle 1996 | 10.1533/9781845698232.4.175 | R2-A | pdf:schieberle1996.pdf (title) |
| Guth 1995 | 10.1021/jf00056a044 | R2-A | pdf:guth1995.pdf (title) |
| Fan 2018 | 10.1021/acs.jafc.8b03297 | Sohail-lead | registry |
| Guth 1994 | 10.1021/jf00048a039 | Sohail-lead | registry |
| Liu 2019 | 10.1021/acs.jafc.9b01564 | Sohail-lead | pdf:liu2019.pdf |
| Du 2020 | 10.1002/ffj.3621 | Sohail-lead | pdf:du2020.pdf |
| Zhao 2017 | 10.1016/j.foodchem.2017.01.011 | Sohail-lead | pdf:zhao2017.pdf |
| Cao 2017 | 10.1016/j.foodres.2017.06.012 | Sun2022-lead | registry |
| Raza 2020 | 10.1039/d0fo01946a | Sun2022-lead | pdf:raza2020.pdf |
