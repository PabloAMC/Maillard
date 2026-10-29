# Literature search 2026-10, round 3: corpus co-citation and sulfur-lane gaps

Run 2026-10-27/28. Third round after [LITERATURE_SEARCH_2026-10.md](LITERATURE_SEARCH_2026-10.md) (rounds 1-2, citation
chaining from seeds) and [CITATION_CONTEXT_2026-10.md](CITATION_CONTEXT_2026-10.md) (how the remaining candidates are cited).
Three parts: (1) which papers the **whole corpus** cites that we do not have; (2) targeted searches for ten mechanisms the
sulfur lane lacks or holds weakly (G1-G10); (3) the citation-context filter applied to both. Search record only; nothing
was added to `data/`.

## Bottom line

- **Corpus co-citation (Part 1) finds no sulfur-lane kinetics we lack.** The papers our corpus cites most often without
  having them are reviews, books, acrylamide classics (Mottram 2002, Stadler 2002, Tareke 2002, Zyzak 2003) and
  flavour-binding papers. Of the top 60, 16 are Q (cited for a number), but none of the 17 fills G1-G10. The quantitative ones are
  trunk or acrylamide kinetics, marked MAYBE below. The most cited sulfur-lane paper we lack is Lee 2010 (GSH + sugars,
  14 corpus citers, all sulfur-lane). It is a labelling/AEDA paper, with no rates.
- **Gap-by-gap (Part 2)**, honestly:
  - G1 thiol sink: no aqueous-buffer stability kinetics of MFT/FFT at two or more temperatures or pH values exists in what
    we found. The quantitative thiol-loss data are all in coffee brews at one temperature: Hofmann 2001, Müller 2007, Sun 2024,
    Kraehenbuehl 2005. The model-meat exception is Grosch 1992, which gives MFT, its disulfide and FFT amounts.
  - G2 TTCA / cysteine Amadori: no new formation/decay rate constants at several temperatures. The useful items are de Roos 1992
    (pH dependence via neutral thiazolidine), Mizutani 1960 (H2S/NH3 release from cysteine thiazolidines at 120 C) and
    Edwards 1992 (cysteine + 3-DH: fast thiazolidine, slow release). Fan 2011 (Chinese, no DOI) claims multi-temperature
    first-order cysteine/glucose loss.
  - G3 cysteine + dicarbonyl: real rate constants exist, but at 37 C. Lo 1994 gives NAC + MGO hemithioacetal forward/reverse
    constants. Edwards & Wedzicha 1991/1992 give 3-deoxyhexosulose + thiol kinetics, which is the closest to a food-model
    step.
  - G4 GSH thermal degradation: nothing at two or more food-relevant temperatures. Only a pharmaceutical forced-degradation
    study (Liu 2020) and a 1994 stability note.
  - G5 thiamine: loss kinetics are already covered (corpus has Voelker 2018/2021, Dwivedi 1972/73). New multi-T/pH loss sets are Lee 1988
    (thesis) and Fox 1982. **H2S yield from thiamine vs T/pH: nothing beyond Dwivedi & Arnold 1971, which is already on disk as
    `arnold1971.pdf`.** Grosch 1992 gives MFT yield; Morfee 1972 says elemental S is a major product at basic pH.
  - G6 cysteine Strecker H2S: **no H2S yield or rate vs dicarbonyl, pH or T found.** Closest: Mizutani 1960 (H2S from
    cysteine-aldehyde adducts, 120 C), Edwards 1992 (slow Strecker-type release). Everything else is product identification.
  - G7 cysteamine/2-acetyl-2-thiazoline: one strong miss. **Hofmann & Schieberle 1995 (10.1021/jf00059a032) is not in the
    corpus.** The corpus `hofmann1995.pdf` is the AEDA paper (jf00056a042). Engel & Schieberle 2002 gives SIDA yields from
    fructose + cysteamine at 145 C. Cysteine decarboxylation to cysteamine has no kinetic paper.
  - G8 R5P/IMP: dephosphorylation kinetics exist. Matoba 1988 is already in the corpus; new are Shaoul & Sporns 1987
    (121 C, pH 3-8, Arrhenius) and Lee 1979 (100 C pH series). **R5P/IMP vs ribose flavour yields with cysteine: nothing beyond
    Mottram & Nobrega 2002 (in corpus).** Farmer 1996 is a title-only lead.
  - G9 thiol-protein binding: **no binding isotherms for MFT/FFT with proteins, and nothing for plant proteins.** The
    only protein study is Mottram, Szauman-Szumski & Dodson 1996 (JAFC, egg albumin/casein reduce flavour disulfides to
    thiols), already on disk as `mottram1996.pdf`. The search surfaced only its Weurman conference version. The rest is the
    coffee-melanoidin literature (mostly in corpus).
  - G10 metal/heme: Cu-catalysed thiol/H2S oxidation kinetics exist, from wine and biochemistry: Kreitman 2017, Kachur 1999,
    with Kreitman 2016a/b already in corpus. Fe/heme-catalysed thiol oxidation kinetics: nothing. Heme effects on
    plant-meat volatiles: one paper, Meng 2025, and its citers mention only colour.

## Method

- **Corpus (575 DOIs)**: every DOI in `data/keys/papers.yml` (535); DOIs printed on the first two pages of every
  `data/articles` PDF (`pdftotext -l 2`, 399 PDFs; 25 DOIs not in the registry); the 13 old PDFs that print no DOI,
  matched by name and title in CITATION_CONTEXT_2026-10.md; the 7 registry-gap dossiers of round 1. That makes 575 unique.
  **Exclusion set (649)**: the corpus plus every DOI written in the two earlier documents.
- **Sulfur-lane flag**: a corpus paper counts as sulfur-lane if its extraction dossier mentions cysteine, thiamine,
  glutathione, MFT, FFT, furanthiol, furfurylthiol or thiol (215 of 575).
- **Part 1**: OpenAlex `filter=doi:a|b|...` in batches of 50, `select=...,referenced_works`. 565 of 575 resolved (518 with a
  reference list). Ten failed: four `jf...+` DOIs whose `+` the batch filter does not accept, two truncated DOIs read from
  dossier headers, four others.
  For every referenced work not in the corpus, count the **distinct corpus papers citing it** (12 106 works; 964 cited by
  three or more). Metadata was fetched for those 964, and the exclusion set was removed. Rank by count, then sulfur-lane citers,
  then OpenAlex citation count. Name/year matches to `data/articles` were checked by hand. Shu 1985 jf00063a030 is the
  Furaneol paper and is not on disk (the disk copy is the cysteine paper), so it stays.
- **Part 2**: per gap, OpenAlex `filter=title_and_abstract.search:` with 6-10 phrasings (50 results each, then 25 for a
  second set). Plus one-hop chaining from the named corpus papers: all citers (cap 400) and all references. Then a per-gap
  topic regex on title and abstract, and the exclusion set. Then manual triage by title and abstract. Seeds that OpenAlex knows
  without a reference list (de Roos 2005, Kang 2026, Zheng 1994, Dwivedi 1973, Hartman 1984 JAFC, Zhang & Ho 1991,
  both Hofmann & Schieberle 1995 papers) only contributed citers.
  - G7's near papers resolved as Cerny & Grosch 1993 (10.1007/bf01190805) and Hofmann & Schieberle 1995
    (10.1021/jf00059a032). Neither is in the corpus; both were used as seeds **and** kept as candidates.
  - Pool sizes, as all hits / after the topic filter / second keyword set: G1 344/50/53, G2 343/21/61, G3 215/82/59,
    G4 234/68/83, G5 361/109/13, G6 270/19/-, G7 358/39/78, G8 229/66/79, G9 271/26/111, G10 137/76/106.
- **Budget failure, and the fallback**: the anonymous OpenAlex daily budget for this IP ran out twice. The first time was
  after G1-G5; it was waited out to the 00:00 UTC reset. The second time was after the G6-G10 chaining run. So the **second keyword set for G7-G10 was run on Crossref**
  (`query.bibliographic`, 20 rows, anonymous, no mailto). Crossref matches mostly titles and rarely has abstracts, so it is noisier
  ("found by" says `kw-Crossref`). Two papers were found by direct Crossref lookups: Shaoul & Sporns 1987 and the Matoba
  1988 check, which turned out to be in the corpus already.
- **Part 3**: Semantic Scholar `GET /graph/v1/paper/DOI:<doi>/citations?fields=title,year,externalIds,contexts,intents,isInfluential`
  (limit 100, a second page where needed, about 1.5 s/request with back-off), anonymous. The round-2 cache was reused.
  For each paper the contexts were ranked: citers in our corpus first, then contexts with numbers or kinetic words, then
  sulfur words. The top 2 (Part 1) or 4-5 (Part 2) were read, plus the OpenAlex/Crossref abstract.
  - Classes as in CITATION_CONTEXT_2026-10.md: **Q** = cited for a usable number. **M** = mechanism or direction only.
    **B** = background, identification or method. **?** = no contexts.
  - One class is added: **Qa** = the contexts are absent or silent, but the **abstract itself** reports measured rates or
    amounts. Without it, most pre-2000 primary papers would be "?" only because nobody quotes them with context.
- **DOWNLOAD** = Q or Qa, and the content plausibly fills G1-G10 or gives multi-temperature (or multi-pH) kinetics of a step
  the engine has or needs. MAYBE = quantitative but off-matrix, single-condition, or hard to obtain; or title-only leads on
  an exact gap. SKIP = the rest.

## Caveats

- **Co-citation favours the trunk.** Most corpus papers are trunk, lipid, acrylamide and binding papers. So raw co-citation counts
  rank their shared reviews first, and sulfur-lane misses sit lower (the supplement table below is ranked by sulfur-lane
  citers instead). A paper cited by one or two corpus papers never enters Part 1. The Part 1 pool is ">= 3 corpus citers".
- **OpenAlex reference lists are incomplete** for old JAFC papers and ACS chapters (0 references for several seeds). So both
  co-citation counts and backward chaining undercount the 1970s-1990s literature, which is where most sulfur kinetics would be.
- **"Qa" rests on abstracts**, which are written by the authors. "Qa" says that numbers are reported, not that they are good.
- **S2 contexts are thin and sometimes mislinked** (Mottram 1996 Weurman, Farmer 1996 Weurman: contexts from other papers).
- **Old PDFs without a printed DOI hide corpus papers.** Two Part 2 hits turned out to be on disk: `hofmann1996b.pdf` (Weurman
  1996) and, through its JAFC version, `mottram1996.pdf`. `arnold1971.pdf` is Dwivedi & Arnold 1971 (thiamine H2S). DOWNLOAD
  and MAYBE names were checked against all 108 no-DOI PDFs by name and by page-1 title where the name matched.
- **Book chapters** (ACS Symp. Ser. 490 chapters by Grosch, de Roos, Ho; Weurman 1996 chapters) may need library access.

## Part 1: top 60 co-cited works not in the corpus

"corpus citers" = distinct corpus papers whose OpenAlex reference list contains the work. "S2 contexts / citers" is what
Semantic Scholar holds; 200 citers is the fetch cap (two pages), not the true count.

| # | corpus citers | of which sulfur-lane | first author, year | title | venue | DOI | type | S2 contexts / citers | class | rec. | note |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 48 | 17 | Hodge 1953 | Dehydrated Foods, Chemistry of Browning Reactions in Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf60015a004](https://doi.org/10.1021/jf60015a004) | article | 72 / 200 | B | SKIP | early review/scheme; cited for the reaction scheme |
| 2 | 24 | 7 | Mottram 2002 | Acrylamide is formed in the Maillard reaction | Nature | [10.1038/419448a](https://doi.org/10.1038/419448a) | article | 83 / 200 | Q | MAYBE | acrylamide from Asn + sugars vs temperature; acrylamide lane, not a G gap |
| 3 | 24 | 4 | Boekel 1996 | Statistical Aspects of Kinetic Modeling for Food Science Problems | Journal of Food Science | [10.1111/j.1365-2621.1996.tb13138.x](https://doi.org/10.1111/j.1365-2621.1996.tb13138.x) | article | 82 / 200 | B | SKIP | statistics of kinetic modelling (method) |
| 4 | 22 | 12 | Parliment 1989 | Thermal Generation of Aromas | ACS symposium series | [10.1021/bk-1989-0409](https://doi.org/10.1021/bk-1989-0409) | conference-paper | 8 / 105 | B | SKIP | whole ACS Symp. Ser. 409 volume; corpus cites its chapters |
| 5 | 22 | 6 | Ledl 1990 | New Aspects of the Maillard Reaction in Foods and in the Human Body | Angewandte Chemie International Edition in English | [10.1002/anie.199005653](https://doi.org/10.1002/anie.199005653) | article | 65 / 200 | B | SKIP | review |
| 6 | 21 | 6 | Stadler 2002 | Acrylamide from Maillard reaction products | Nature | [10.1038/419449a](https://doi.org/10.1038/419449a) | article | 64 / 200 | Q | SKIP | acrylamide yields at one condition |
| 7 | 21 | 5 | Yaylayan 2003 | Recent Advances in the Chemistry of Strecker Degradation and Amadori Rearrangement: Implications to Aroma and Color Formation | Food Science and Technology Research | [10.3136/fstr.9.1](https://doi.org/10.3136/fstr.9.1) | article | 82 / 198 | B | SKIP | review (Strecker / Amadori) |
| 8 | 20 | 4 | Kroh 1994 | Caramelisation in food and beverages | Food Chemistry | [10.1016/0308-8146(94)90188-0](https://doi.org/10.1016/0308-8146(94)90188-0) | article | 115 / 200 | B | SKIP | caramelisation review |
| 9 | 18 | 11 | Werkhoff 1990 | Isolation and characterization of volatile sulfur-containing meat flavor components in model systems | Journal of Agricultural and Food Chemistry | [10.1021/jf00093a041](https://doi.org/10.1021/jf00093a041) | article | 17 / 64 | B | SKIP | identification of meat S-volatiles |
| 10 | 18 | 7 | Ajandouz 2001 | Effects of pH on Caramelization and Maillard Reaction Kinetics in Fructose‐Lysine Model Systems | Journal of Food Science | [10.1111/j.1365-2621.2001.tb08213.x](https://doi.org/10.1111/j.1365-2621.2001.tb08213.x) | article | 81 / 200 | Q | MAYBE | pH 4-12 kinetics of fructose-lysine browning and caramelisation (trunk; corpus has Ajandouz 2008) |
| 11 | 17 | 5 | Yaylayan 1994 | Chemistry of Amadori rearrangement products: Analysis, synthesis, kinetics, reactions, and spectroscopic properties | Critical Reviews in Food Science and Nutrition | [10.1080/10408399409527667](https://doi.org/10.1080/10408399409527667) | article | 43 / 200 | B | SKIP | ARP review |
| 12 | 16 | 5 | Wang 2021 | Insights into formation, detection and removal of the beany flavor in soybean protein | Trends in Food Science & Technology | [10.1016/j.tifs.2021.04.018](https://doi.org/10.1016/j.tifs.2021.04.018) | article | 92 / 197 | B | SKIP | beany-flavour review |
| 13 | 15 | 10 | Charalambous1 1976 | Phenolic, Sulfur, and Nitrogen Compounds in Food Flavors | ACS symposium series | [10.1021/bk-1976-0026](https://doi.org/10.1021/bk-1976-0026) | conference-paper | 7 / 55 | B | SKIP | ACS Symp. Ser. 26 volume |
| 14 | 15 | 8 | Boekel 2006 | Formation of flavour compounds in the Maillard reaction | Biotechnology Advances | [10.1016/j.biotechadv.2005.11.004](https://doi.org/10.1016/j.biotechadv.2005.11.004) | article | 68 / 200 | B | SKIP | review; Arrhenius generalities |
| 15 | 15 | 8 | Whitfield 1992 | Volatiles from interactions of Maillard reactions and lipids | Critical Reviews in Food Science and Nutrition | [10.1080/10408399209527560](https://doi.org/10.1080/10408399209527560) | review | 106 / 200 | B | SKIP | Maillard-lipid review |
| 16 | 15 | 5 | Tareke 2002 | Analysis of Acrylamide, a Carcinogen Formed in Heated Foodstuffs | Journal of Agricultural and Food Chemistry | [10.1021/jf020302f](https://doi.org/10.1021/jf020302f) | article | 101 / 200 | Q | SKIP | acrylamide occurrence in foods |
| 17 | 15 | 4 | Zyzak 2003 | Acrylamide Formation Mechanism in Heated Foods | Journal of Agricultural and Food Chemistry | [10.1021/jf034180i](https://doi.org/10.1021/jf034180i) | article | 83 / 200 | M | SKIP | acrylamide mechanism |
| 18 | 14 | 14 | Lee 2010 | Investigation of the Aroma-Active Compounds Formed in the Maillard Reaction between Glutathione and Reducing Sugars | Journal of Agricultural and Food Chemistry | [10.1021/jf9043327](https://doi.org/10.1021/jf9043327) | article | 29 / 68 | B | MAYBE | GSH + glucose/fructose AEDA and CAMOLA; all 14 corpus citers are sulfur-lane (G4) |
| 19 | 14 | 7 | Lund 2017 | Control of Maillard Reactions in Foods: Strategies and Chemical Mechanisms | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.7b00882](https://doi.org/10.1021/acs.jafc.7b00882) | review | 84 / 200 | B | SKIP | review |
| 20 | 14 | 7 | Gremli 1974 | Interaction of flavor compounds with soy protein | Journal of the American Oil Chemists Society | [10.1007/bf02542100](https://doi.org/10.1007/bf02542100) | article | 33 / 93 | Q | SKIP | soy protein binds aldehydes 2-4x ketones (carbonyl binding, not thiols) |
| 21 | 14 | 5 | Leong 2000 | A critical appraisal of the kinetic model for the Maillard browning of glucose with glycine | Food Chemistry | [10.1016/s0308-8146(99)00146-6](https://doi.org/10.1016/s0308-8146(99)00146-6) | article | 11 / 50 | Q | MAYBE | glucose-glycine browning kinetic model (trunk) |
| 22 | 13 | 9 | Shakoor 2022 | Maillard reaction chemistry in formation of critical intermediates and flavour compounds and their antioxidant properties | Food Chemistry | [10.1016/j.foodchem.2022.133416](https://doi.org/10.1016/j.foodchem.2022.133416) | review | 62 / 200 | B | SKIP | review |
| 23 | 13 | 9 | Blank 1992 | Potent odorants of the roasted powder and brew of Arabica coffee | European Food Research and Technology | [10.1007/bf01202802](https://doi.org/10.1007/bf01202802) | article | 103 / 200 | B | SKIP | coffee AEDA (identification) |
| 24 | 13 | 9 | Cui 2021 | Formation and fate of Amadori rearrangement products in Maillard reaction | Trends in Food Science & Technology | [10.1016/j.tifs.2021.06.055](https://doi.org/10.1016/j.tifs.2021.06.055) | article | 67 / 200 | B | SKIP | ARP review |
| 25 | 13 | 8 | Scalone 2015 | Influence of Free Amino Acids, Oligopeptides, and Polypeptides on the Formation of Pyrazines in Maillard Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.5b01129](https://doi.org/10.1021/acs.jafc.5b01129) | article | 33 / 137 | M | SKIP | peptides vs amino acids in volatile formation |
| 26 | 13 | 5 | Wang 2022 | Flavor challenges in extruded plant‐based meat alternatives: A review | Comprehensive Reviews in Food Science and Food Safety | [10.1111/1541-4337.12964](https://doi.org/10.1111/1541-4337.12964) | review | 128 / 200 | B | SKIP | review |
| 27 | 13 | 4 | Roland 2016 | Flavor Aspects of Pulse Ingredients | Cereal Chemistry | [10.1094/cchem-06-16-0161-fi](https://doi.org/10.1094/cchem-06-16-0161-fi) | article | 160 / 200 | B | SKIP | pulse-flavour review |
| 28 | 13 | 2 | Hollnagel 1998 | Formation of α-dicarbonyl fragments from mono- and disaccharides under caramelization and Maillard reaction conditions | Zeitschrift für Lebensmittel-Untersuchung und -Forschung. B, Referate und Lebensmittelrecht | [10.1007/s002170050294](https://doi.org/10.1007/s002170050294) | article | 50 / 121 | M | MAYBE | dicarbonyl fragments from sugars (corpus has Hollnagel 2000) |
| 29 | 12 | 9 | Verhue 1994 | Recent developments in flavor and fragrance chemistry | Trends in Food Science & Technology | [10.1016/0924-2244(94)90205-4](https://doi.org/10.1016/0924-2244(94)90205-4) | book-review | 0 / 34 | ? | SKIP | book review; no contexts |
| 30 | 12 | 7 | Teranishi 1992 | FLAVOR PRECURSORS: THERMAL AND ENZYMATIC CONVERSIONS | Journal of Food Biochemistry | [10.1111/j.1745-4514.1992.tb00440.x](https://doi.org/10.1111/j.1745-4514.1992.tb00440.x) | article | 14 / 69 | B | SKIP | precursor review |
| 31 | 12 | 4 | Thornalley 1999 | Formation of glyoxal, methylglyoxal and 3-deoxyglucosone in the glycation of proteins by glucose | Biochemical Journal | [10.1042/0264-6021:3440109](https://doi.org/10.1042/0264-6021:3440109) | article | not in S2 | ? | MAYBE | not in S2; GO/MGO/3-DG formation in glycation (quantitative, 37 C) |
| 32 | 12 | 4 | Heng 2004 | Protein–flavour interactions in relation to development of novel protein foods | Trends in Food Science & Technology | [10.1016/j.tifs.2003.09.018](https://doi.org/10.1016/j.tifs.2003.09.018) | article | 57 / 113 | Q | SKIP | pea vicilin binds aldehydes/ketones (not thiols) |
| 33 | 12 | 3 | Hofmann 1999 | Quantitative studies on the role of browning precursors in the Maillard reaction of pentoses and hexoses with l -alanine | European Food Research and Technology | [10.1007/s002170050468](https://doi.org/10.1007/s002170050468) | article | 43 / 73 | Q | MAYBE | glyoxal time course in glucose/alanine (dicarbonyl supply for G3/G6) |
| 34 | 12 | 1 | Locas 2008 | Isotope Labeling Studies on the Formation of 5-(Hydroxymethyl)-2-furaldehyde (HMF) from Sucrose by Pyrolysis-GC/MS | Journal of Agricultural and Food Chemistry | [10.1021/jf8010245](https://doi.org/10.1021/jf8010245) | article | 86 / 200 | Q | SKIP | HMF 13C labelling fractions |
| 35 | 11 | 10 | Mulders 1973 | Volatile components from the non-enzymic browning reaction of the cysteine/cystine-ribose system | European Food Research and Technology | [10.1007/bf01106493](https://doi.org/10.1007/bf01106493) | article | 12 / 71 | M | SKIP | 2-acetylthiazole from cysteine + ribose via MGO (G7-adjacent, mechanism) |
| 36 | 11 | 7 | Dool 1963 | A generalization of the retention index system including linear temperature programmed gas—liquid partition chromatography | Journal of Chromatography A | [10.1016/s0021-9673(01)80947-x](https://doi.org/10.1016/s0021-9673(01)80947-x) | article | 83 / 200 | B | SKIP | retention-index method |
| 37 | 11 | 7 | Gasser 1988 | Identification of volatile flavour compounds with high aroma values from cooked beef | European Food Research and Technology | - | article | - | ? | SKIP | no DOI; cooked-beef AEDA |
| 38 | 11 | 4 | Guichard 2002 | Interactions between flavor compounds and food ingredients and their influence on flavor perception | Food Reviews International | [10.1081/fri-120003417](https://doi.org/10.1081/fri-120003417) | article | 85 / 200 | B | SKIP | flavour-matrix interaction review |
| 39 | 11 | 4 | Yaylayan 2003 | Why Asparagine Needs Carbohydrates To Generate Acrylamide | Journal of Agricultural and Food Chemistry | [10.1021/jf0261506](https://doi.org/10.1021/jf0261506) | article | 94 / 200 | M | SKIP | acrylamide mechanism |
| 40 | 11 | 3 | Lam 2016 | Pea protein isolates: Structure, extraction, and functionality | Food Reviews International | [10.1080/87559129.2016.1242135](https://doi.org/10.1080/87559129.2016.1242135) | article | 81 / 200 | B | SKIP | pea protein review |
| 41 | 11 | 3 | Taeymans 2004 | A Review of Acrylamide: An Industry Perspective on Research, Analysis, Formation, and Control | Critical Reviews in Food Science and Nutrition | [10.1080/10408690490478082](https://doi.org/10.1080/10408690490478082) | article | 126 / 200 | B | SKIP | acrylamide review |
| 42 | 11 | 3 | Brands 2003 | Kinetic modelling of reactions in heated disaccharide–casein systems | Food Chemistry | [10.1016/s0308-8146(03)00031-1](https://doi.org/10.1016/s0308-8146(03)00031-1) | article | 21 / 64 | Q | MAYBE | kinetic model of disaccharide-casein at 120 C (trunk; corpus has Brands 2001/2002) |
| 43 | 11 | 2 | Waller 1983 | The Maillard Reaction in Foods and Nutrition | ACS symposium series | [10.1021/bk-1983-0215](https://doi.org/10.1021/bk-1983-0215) | conference-paper | 34 / 200 | B | SKIP | ACS Symp. Ser. 215 volume |
| 44 | 11 | 1 | Boekel 1998 | Effect of heating on Maillard reactions in milk | Food Chemistry | [10.1016/s0308-8146(98)00075-2](https://doi.org/10.1016/s0308-8146(98)00075-2) | article | 98 / 200 | Q | SKIP | milk Maillard, Ea about 120 kJ/mol |
| 45 | 11 | 0 | Hull 2011 | Nε-(carboxymethyl)lysine content of foods commonly consumed in a Western style diet | Food Chemistry | [10.1016/j.foodchem.2011.08.055](https://doi.org/10.1016/j.foodchem.2011.08.055) | article | 141 / 200 | Q | SKIP | CML in foods |
| 46 | 10 | 8 | McHale 1991 | Flavour science and technology | Food Chemistry | [10.1016/0308-8146(91)90139-f](https://doi.org/10.1016/0308-8146(91)90139-f) | book-review | 0 / 0 | ? | SKIP | book review |
| 47 | 10 | 6 | Boekel 2008 | Kinetic Modeling of Reactions In Foods | CRC Press eBooks | [10.1201/9781420017410](https://doi.org/10.1201/9781420017410) | reference-entry | 144 / 200 | B | SKIP | kinetics textbook (method) |
| 48 | 10 | 6 | Davídek 2017 | Preparation of activated flavor precursor DFG, N-(1-deoxy-1-fructosylglycine) by combination of vacuum evaporation and closed system heating steps | Food Chemistry | [10.1016/j.foodchem.2017.10.013](https://doi.org/10.1016/j.foodchem.2017.10.013) | article | 20 / 42 | Q | SKIP | ARP preparation yields (process engineering) |
| 49 | 10 | 5 | Nooshkam 2018 | The Maillard reaction products as food-born antioxidant and antibrowning agents in model and real food systems | Food Chemistry | [10.1016/j.foodchem.2018.09.083](https://doi.org/10.1016/j.foodchem.2018.09.083) | article | 114 / 200 | B | SKIP | review |
| 50 | 10 | 5 | Schieberle 1991 | Primary odorants in popcorn | Journal of Agricultural and Food Chemistry | [10.1021/jf00006a030](https://doi.org/10.1021/jf00006a030) | article | 62 / 178 | B | SKIP | popcorn odorants (identification) |
| 51 | 10 | 4 | Laroque 2008 | Kinetic study on the Maillard reaction. Consideration of sugar reactivity | Food Chemistry | [10.1016/j.foodchem.2008.05.033](https://doi.org/10.1016/j.foodchem.2008.05.033) | article | 90 / 200 | M | MAYBE | sugar reactivity in Maillard kinetics; cited for ribose > xylose |
| 52 | 10 | 4 | Murat 2013 | Characterisation of odour active compounds along extraction process from pea flour to pea protein extract | Food Research International | [10.1016/j.foodres.2013.03.049](https://doi.org/10.1016/j.foodres.2013.03.049) | article | 57 / 118 | Q | SKIP | pea-flour odorant amounts |
| 53 | 10 | 4 | Kühn 2008 | Binding of Flavor Compounds and Whey Protein Isolate as Affected by Heat and High Pressure Treatments | Journal of Agricultural and Food Chemistry | [10.1021/jf801810b](https://doi.org/10.1021/jf801810b) | article | 39 / 89 | Q | SKIP | whey protein-carbonyl binding |
| 54 | 10 | 2 | Zamora 2005 | Coordinate Contribution of Lipid Oxidation and Maillard Reaction to the Nonenzymatic Food Browning | Critical Reviews in Food Science and Nutrition | [10.1080/10408690590900117](https://doi.org/10.1080/10408690590900117) | review | 92 / 200 | B | SKIP | lipid oxidation x Maillard review |
| 55 | 10 | 1 | Capuano 2010 | Acrylamide and 5-hydroxymethylfurfural (HMF): A review on metabolism, toxicity, occurrence in food and mitigation strategies | LWT | [10.1016/j.lwt.2010.11.002](https://doi.org/10.1016/j.lwt.2010.11.002) | article | 91 / 200 | B | SKIP | acrylamide/HMF review |
| 56 | 10 | 1 | Fu 1996 | The Advanced Glycation End Product, N∊-(Carboxymethyl)lysine, Is a Product of both Lipid Peroxidation and Glycoxidation Reactions | Journal of Biological Chemistry | [10.1074/jbc.271.17.9982](https://doi.org/10.1074/jbc.271.17.9982) | article | 126 / 200 | M | SKIP | CML formation routes |
| 57 | 10 | 1 | O'Brien 1998 | The Maillard Reaction in Foods and Medicine | Woodhead Publishing Limited eBooks | [10.1533/9781845698447](https://doi.org/10.1533/9781845698447) | book | 36 / 153 | B | SKIP | book |
| 58 | 9 | 7 | Evers 1976 | Furans Substituted at the Three Position with Sulfur | ACS symposium series | [10.1021/bk-1976-0026.ch011](https://doi.org/10.1021/bk-1976-0026.ch011) | conference-paper | 3 / 27 | M | SKIP | meaty odour of 3-thio furans (MFT family) |
| 59 | 9 | 5 | Land 1979 | Progress in flavour research |  | - | book | - | ? | SKIP | book, no DOI |
| 60 | 9 | 5 | Ashoor 1984 | Maillard Browning of Common Amino Acids and Sugars | Journal of Food Science | [10.1111/j.1365-2621.1984.tb10432.x](https://doi.org/10.1111/j.1365-2621.1984.tb10432.x) | article | 78 / 200 | Q | SKIP | browning of amino acids + sugars at one condition |

### Supplement: works outside the top 60 with five or more sulfur-lane corpus citers

| corpus citers | sulfur-lane | first author, year | title | DOI |
|---|---|---|---|---|
| 8 | 6 | Shu 1985 | Volatile components of the thermal degradation of 2,5-dimethyl-4-hydroxy-3(2H)-furanone | [10.1021/jf00063a030](https://doi.org/10.1021/jf00063a030) |
| 7 | 6 | Engel 1999 | Solvent assisted flavour evaporation - a new and versatile technique for the careful and direct isolation of aroma compounds from complex food matrices | [10.1007/s002170050486](https://doi.org/10.1007/s002170050486) |
| 7 | 6 | Adams 2011 | Model Studies on the Pattern of Volatiles Generated in Mixtures of Amino Acids, Lipid-Oxidation-Derived Aldehydes, and Glucose | [10.1021/jf104091p](https://doi.org/10.1021/jf104091p) |
| 7 | 6 | Anantharamkrishnan 2020 | Method To Characterize and Monitor Covalent Interactions of Flavor Compounds with β‑Lactoglobulin Using Mass Spectrometry and Proteomics | [10.1021/acs.jafc.9b07978](https://doi.org/10.1021/acs.jafc.9b07978) |
| 7 | 6 | Zhan 2021 | Characteristic flavor formation of thermally processed N-(1-deoxy-α-d-ribulos-1-yl)-glycine: Decisive role of additional amino acids and promotional effect of glyoxal | [10.1016/j.foodchem.2021.131137](https://doi.org/10.1016/j.foodchem.2021.131137) |
| 6 | 6 | Semmelroch 1996 | Studies on Character Impact Odorants of Coffee Brews | [10.1021/jf9505988](https://doi.org/10.1021/jf9505988) |
| 9 | 5 | Tressl 1981 | Investigation of sulfur-containing components in roasted coffee | [10.1021/jf00107a045](https://doi.org/10.1021/jf00107a045) |
| 8 | 5 | Grosch 1993 | Detection of potent odorants in foods by aroma extract dilution analysis | [10.1016/0924-2244(93)90187-f](https://doi.org/10.1016/0924-2244(93)90187-f) |
| 8 | 5 | Davídek 2006 | Sugar Fragmentation in the Maillard Reaction Cascade: Formation of Short-Chain Carboxylic Acids by a New Oxidative α-Dicarbonyl Cleavage Pathway | [10.1021/jf060668i](https://doi.org/10.1021/jf060668i) |
| 8 | 5 | Nikolov 2011 | Thermal Decomposition of 5-(Hydroxymethyl)-2-furaldehyde (HMF) and Its Further Transformations in the Presence of Glycine | [10.1021/jf202470u](https://doi.org/10.1021/jf202470u) |
| 8 | 5 | Vernin 1982 | Chemistry of Heterocyclic Compounds in Flavours and Aromas | - |
| 7 | 5 | Mottram 1983 | The role of triglycerides and phospholipids in the aroma of cooked beef | [10.1002/jsfa.2740340513](https://doi.org/10.1002/jsfa.2740340513) |
| 7 | 5 | Adams 2011 | Amino Acid Catalysis of 2-Alkylfuran Formation from Lipid Oxidation-Derived α,β-Unsaturated Aldehydes | [10.1021/jf202448v](https://doi.org/10.1021/jf202448v) |
| 7 | 5 | Salter 1989 | Volatile compounds produced in maiilard reactions involving glycine, ribose and phospholipid | [10.1002/jsfa.2740460211](https://doi.org/10.1002/jsfa.2740460211) |
| 7 | 5 | Mundt 2003 | A Kinetic Model for the Glucose−Fructose−Glycine Browning Reaction | [10.1021/jf026027e](https://doi.org/10.1021/jf026027e) |
| 6 | 5 | Schieberle 1987 | Quantitative analysis of aroma compounds in wheat and rye bread crusts using a stable isotope dilution assay | [10.1021/jf00074a021](https://doi.org/10.1021/jf00074a021) |
| 6 | 5 | Schutte 1974 | Precursors of sulfur‐containing flavor compounds | [10.1080/10408397409527166](https://doi.org/10.1080/10408397409527166) |
| 6 | 5 | Davídek 2004 | Analysis of Amadori Compounds by High-Performance Cation Exchange Chromatography Coupled to Tandem Mass Spectrometry | [10.1021/ac048925a](https://doi.org/10.1021/ac048925a) |
| 6 | 5 | Troise 2016 | Evolution of protein bound Maillard reaction end-products and free Amadori compounds in low lactose milk in presence of fructosamine oxidase I | [10.1016/j.foodchem.2016.06.037](https://doi.org/10.1016/j.foodchem.2016.06.037) |
| 6 | 5 | Mottram 2017 | Meat | [10.1201/9780203734285-5](https://doi.org/10.1201/9780203734285-5) |
| 6 | 5 | Parisi 2018 | Maillard Reaction in Processed Foods—Reaction Mechanisms | [10.1007/978-3-319-95463-9_2](https://doi.org/10.1007/978-3-319-95463-9_2) |

Sulfur-relevant entries in the supplement: Shu 1985b (Furaneol thermal degradation), Tressl 1981 (coffee S-volatiles), Schutte
1974 (S-precursors review). All are identification or review papers, so none is a DOWNLOAD.

## Part 2: gap-targeted search

Each table lists the triaged, topic-relevant candidates, not the full pool. "found by": `kw` = OpenAlex keyword search,
`kw-Crossref` = Crossref fallback, `cites X` / `cited by X` = one-hop link to the named near paper, `P1 rank` = also in the
Part 1 co-citation list.

### G1 thiol sink (MFT/FFT/thiol loss kinetics, disulfide formation)

Near papers: Hofmann 1996 (jf9500703), Blank 2002 (jf011329m), Kumazawa 2003 (jf021025f), Guth & Grosch 1995 (jf00056a044).
**Verdict: no buffer kinetics at two or more temperatures or pH values. The quantitative loss data are coffee-brew, one temperature.**
Grosch 1992 is the one model-meat quantitation.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Grosch 1992 (Grosch & Zeiler-Hilgart) | Formation of meatlike flavor compounds (ACS Symp. Ser. 490) | [10.1021/bk-1992-0490.ch014](https://doi.org/10.1021/bk-1992-0490.ch014) | kw; P1 rank 398 (4 corpus citers) | 2 / 12 | Qa | amounts of MFT, its disulfide and FFT from thiamine vs ribose + cysteine at meat-level precursor concentrations; cysteine and H2S raise MFT from thiamine; disulfide and FFT very low (also G5) | DOWNLOAD #4 |
| Hofmann 2001 (Hofmann, Czerny, Calligaris, Schieberle) | Model studies on the influence of coffee melanoidins on flavor volatiles | [10.1021/jf0012042](https://doi.org/10.1021/jf0012042) | cited by Gigl2021, cited by Hofmann2001, cited by Kumazawa2003; P1 rank 150 (6 corpus citers) | 34 / 95 | Q | citers quote FFT down 50 % after 20 min and about 100 % after 30 min with melanoidins (headspace, model brew); 6 corpus citers (also G9) | DOWNLOAD #5 |
| Müller 2007 (Müller & Hofmann) | Phenol/2-furfurylthiol conjugates in coffee beverages (SIDA) | [10.1021/jf070095p](https://doi.org/10.1021/jf070095p) | kw; P1 rank 344 (4 corpus citers) | 19 / 39 | Qa | abstract: about 60 % of FFT reacts with hydroxyhydroquinone within 10 min of brewing; conjugate concentrations by SIDA over storage time | DOWNLOAD #15 |
| Sun 2024 | Aroma binding factors of coffee matrix in a model system | [10.1016/j.lwt.2024.116930](https://doi.org/10.1016/j.lwt.2024.116930) | cites Blank2002, cites Charles-Bernard2005, cites Hofmann2001, cites Kumazawa2003, kw | 4 / 7 | Qa | abstract: in 30 min of brew incubation FFT falls to 3.3 %, MFT to 34.4 %, methanethiol to 25.5 %; loss attributed among lipids, hydroxyhydroquinone, melanoidin size fractions, carbohydrates (also G9) | DOWNLOAD #16 |
| Kraehenbuehl 2005 | Coffee aroma quality: reaction pathways in the degradation of thiols (ASIC 2004 proceedings) | no DOI | kw (OpenAlex record without DOI) | - | Qa | abstract: k_obs for a set of thiols in coffee brew with pathway-specific inhibitors; aliphatic thiols by nucleophilic addition, FFT also by radical cascade. Grey literature (Nestlé, ASIC) | MAYBE |
| Müller 2005 (Müller & Hofmann) | Raw coffee screened for thiol binding-site precursors (in-bean roasting) | [10.1021/jf048027f](https://doi.org/10.1021/jf048027f) | cited by Gigl2021, cites Hofmann2001, kw | 5 / 26 | M | free FFT after 20 min at 30 C falls with roast degree; citers use the CGA-derived binding-site claim (also G9) | MAYBE |
| Kerler 1997 (Kerler & Grosch) | Character impact odorants of boiled chicken: refrigerated storage and reheating | [10.1007/s002170050157](https://doi.org/10.1007/s002170050157) | kw; P1 rank 164 (6 corpus citers) | 19 / 69 | M | citers: FFT decreases on storage/reheating (direction; SIDA numbers likely in paper); matrix, not buffer | MAYBE |
| Dulsat-Serra 2016 | Volatile thiols in coffee: formation, degradation (review) | [10.1016/j.foodres.2016.02.008](https://doi.org/10.1016/j.foodres.2016.02.008) | cites Blank2002, cites Hofmann1996, cites Kumazawa2003, kw | 19 / 77 | B | review; cited for thresholds | SKIP |
| Gigl 2022 | Reaction products of chlorogenic/caffeic acid and their degradation products with odor-active thiols | [10.1021/acs.jafc.2c01378](https://doi.org/10.1021/acs.jafc.2c01378) | kw | 2 / 19 | B | LC-MS markers of thiol-phenol adducts; corpus already has Gigl 2021 | SKIP |
| Hofmann 2003 | Influence of melanoidins on the aroma staling of coffee beverage | [10.1021/bk-2004-0871.ch015](https://doi.org/10.1021/bk-2004-0871.ch015) | kw | 2 / 6 | M | chapter restating Hofmann & Schieberle 2002 (in corpus) | SKIP |
| Müller 2006 | Synthesis of FFT/di- and trihydroxybenzene conjugates | [10.1021/jf062728q](https://doi.org/10.1021/jf062728q) | kw; P1 rank 653 (3 corpus citers) | 8 / 22 | M | structures only | SKIP |
| Buettner 2002 | Human salivary enzymes and odorant concentration (esters, thiols) | [10.1021/jf011586r](https://doi.org/10.1021/jf011586r) | kw | 71 / 118 | Q | FFT 70 to 20 % in 10 min in saliva: enzymatic, out of scope | SKIP |
| Leclercq 2007 | Oil type and delta-tocopherol on model flavor stability | [10.1021/jf071686i](https://doi.org/10.1021/jf071686i) | cites Blank2002 | 7 / 15 | B | oil-matrix storage; contexts about oil oxidation | SKIP |
| Leclercq 2008 | Antioxidants vs sulphur flavour oxidation in sunflower oil | [10.1002/ffj.1889](https://doi.org/10.1002/ffj.1889) | cites Blank2002, kw | 0 / 4 | ? | no contexts; oil matrix | SKIP |
| Zhang 2022 | Refrigerated-reheated xylose-cysteine-lecithin model (GCxGC) | [10.1111/1750-3841.16059](https://doi.org/10.1111/1750-3841.16059) | kw | 5 / 11 | B | FFT rises after refrigeration (direction only) | SKIP |
| Majcenovic 2002 | Ethanethiol / diethyl disulfide SIDA in wine | [10.1021/jf020478h](https://doi.org/10.1021/jf020478h) | cites Guth1995, cites Hofmann1996 | 8 / 30 | Q | fraction of disulfide formed from thiol during analysis; wine | SKIP |
| Nikolantonaki 2012 | 3-sulfanylhexanol + o-quinone adduct kinetics (wine) | [10.1021/jf204295s](https://doi.org/10.1021/jf204295s) | cites Blank2002 | 31 / 66 | Q | thiol + quinone adduct kinetics in wine model; mechanism analogue of the coffee thiol sink | SKIP |
| Hou 2025 | Stability of isolated GSH and cysteine polysulfides (wine) | [10.1016/j.foodchem.2025.144603](https://doi.org/10.1016/j.foodchem.2025.144603) | kw | 19 / 4 | Qa | polysulfide decay at 25/35/40 C; half-life of GS-S-SG under 9 h at room temperature (also G4) | MAYBE |

### G2 TTCA / cysteine Amadori formation and decay

Near papers: Cao 2017, Zhai 2019, Hou 2017, de Roos 2005, Kang 2026 (no OpenAlex references).
**Verdict: no new multi-temperature rate constants for TTCA or Cys-ARP.** The Chinese-lab TTCA papers were already
covered by rounds 1-2. de Roos 1992 and Mizutani 1960 add the pH dependence and the thiazolidine-to-H2S step.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| de Roos 1992 | Meat flavor generation from cysteine and sugars (ACS Symp. Ser. 490) | [10.1021/bk-1992-0490.ch016](https://doi.org/10.1021/bk-1992-0490.ch016) | kw; P1 rank 315 (4 corpus citers) | 4 / 17 | Qa | reactivity of cysteine-sugar systems scales with the neutral (not anionic) thiazolidinecarboxylic-acid concentration, i.e. a pH-dependent TCA step; cysteine inhibits browning via stable thiazolidines; Amadori + cysteine more reactive. Predecessor of de Roos 2005 (on disk) | DOWNLOAD #7 |
| Mizutani 1960 | Degradation mechanisms of proteins and derivatives in foods (cysteine-aldehyde compounds) | [10.1080/03758397.1960.10857681](https://doi.org/10.1080/03758397.1960.10857681) | kw | 0 / 0 | Qa | H2S and NH3 released from cysteine-aldehyde thiazolidines (formaldehyde ... furfural) heated in water at 120 C, measured; the furfural adduct is unstable and releases most H2S (also G6) | DOWNLOAD #6 |
| Edwards 1992 | see G3 | [10.1080/02652039209374098](https://doi.org/10.1080/02652039209374098) | kw | 2 / 9 | Qa | aminothiols (cysteine) react fast with 3-deoxyhexosulose, probably to a thiazolidine, then release thiol slowly | DOWNLOAD #2 |
| Canle 1996 | Rate and equilibrium constants for oxazolidine and thiazolidine ring opening | [10.1351/pac199668040813](https://doi.org/10.1351/pac199668040813) | kw | 8 / 13 | Q | ring-opening rate/equilibrium constants for simple thiazolidines (thiazolidine favoured over Schiff base by about 10^5); not sugar-derived | MAYBE |
| Radomski 1989 | Thiazolidine-4(R)-carboxylic acids from sugars: C-2 epimerisation in aqueous solution | [10.1016/0008-6215(89)80005-9](https://doi.org/10.1016/0008-6215(89)80005-9) | kw | 2 / 13 | B | epimerisation of sugar TCAs in water (NMR); contexts are synthesis citations | MAYBE |
| Roberts 1987 | Prodrugs of L-cysteine: 2-(polyhydroxyalkyl)thiazolidine-4-carboxylic acids | [10.1021/jm00393a034](https://doi.org/10.1021/jm00393a034) | kw | 29 / 99 | B | sugar-cysteine TCAs (incl. ribose-cysteine) release cysteine by non-enzymatic ring opening/hydrolysis; release rates possibly reported (37 C) | MAYBE |
| Belikov 1967 | Stability of the products of condensation of L-cysteine with carbohydrates | [10.1007/bf01171720](https://doi.org/10.1007/bf01171720) | kw | 0 / 1 | ? | title fits exactly; Russian, 1 citer, no abstract | MAYBE |
| Fan 2011 | Kinetics of Maillard reaction in glucose/L-cysteine model systems (Chinese) | no DOI | kw (OpenAlex record without DOI) | - | Qa | abstract: glucose and cysteine loss first order at several temperatures; browning Ea 164.8 kJ/mol, pH-drop Ea 149.2 kJ/mol. No DOI; access uncertain | MAYBE |
| Hou 2008 | Kinetics of the L-arabinose/cysteine Maillard reaction (Chinese thesis) | no DOI | kw (OpenAlex record without DOI) | - | Qa | arabinose loss by HPLC-ELSD under several conditions; arabinose not a modelled sugar | SKIP |
| Vauthey 2000 | Structured fluids as microreactors (cysteine/furfural, cysteine/ribose) | [10.1021/jf991254a](https://doi.org/10.1021/jf991254a) | kw; P1 rank 382 (4 corpus citers) | 12 / 65 | Q | FFT 5-7x higher in L2/cubic phases than water; microemulsion (corpus has Yaghmur) | SKIP |

### G3 cysteine + dicarbonyl kinetics

Near papers: Liu 2022 (acs.jafc.2c04874), Zhou 2022 (acs.jafc.2c07026).
**Verdict: rate constants exist, at 37 C** (Lo 1994; Rae 1990 for GSH). The 3-deoxyhexosulose + thiol kinetics of Wedzicha's
group (1991/1992) are the closest food-chemistry data. No cysteine + glyoxal/MGO/butanedione rate at cooking temperatures.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Edwards 1992 (Edwards & Wedzicha) | Kinetics and mechanism of the reaction between 3-deoxyhexosulose and thiols | [10.1080/02652039209374098](https://doi.org/10.1080/02652039209374098) | kw | 2 / 9 | Qa | rate-determining step is spontaneous conversion of 3-DH to a reactive intermediate, then fast thiol capture (1:2 stoichiometry for mercaptoethanol, GSH); N-acetylcysteine less reactive; cysteine reacts fast then slowly releases thiol, attributed to Strecker degradation (also G2, G6) | DOWNLOAD #2 |
| Wedzicha 1991 (Wedzicha & Edwards) | Kinetics of the reaction of 3-deoxyhexosulose with thiols | [10.1016/0308-8146(91)90020-o](https://doi.org/10.1016/0308-8146(91)90020-o) | kw | 0 / 8 | ? | companion paper; 8 citers, no contexts, no abstract in OpenAlex | MAYBE (fetch with #2) |
| Lo 1994 (Thornalley group) | Binding and modification of proteins by methylglyoxal: kinetics with N-acetyl amino acids | [10.1016/s0021-9258(18)31635-1](https://doi.org/10.1016/s0021-9258(18)31635-1) | kw; P1 rank 106 (7 corpus citers) | 112 / 200 | Q | citers quote k(MGO + N-acetylcysteine hemithioacetal) and the reverse first-order constant 7.5e-3 s-1 (pH 7.4, 37 C); 7 corpus citers (P1 rank 106) | DOWNLOAD #3 |
| Rae 1990 | Kinetic analysis of the erythrocyte glyoxalase system by 1H NMR | [10.1111/j.1432-1033.1990.tb19307.x](https://doi.org/10.1111/j.1432-1033.1990.tb19307.x) | kw | 31 / 76 | Q | rate constants of MGO hydration and GSH-MGO hemithioacetal formation; 37 C, physiological | MAYBE |
| Marchand 2000 | Release of aroma compounds from cysteine + carbonyls under wine conditions | [10.1021/jf000149u](https://doi.org/10.1021/jf000149u) | cited by Zhou2022, kw; P1 rank 605 (3 corpus citers) | 45 / 83 | M | cysteine + glyoxal/MGO/diacetyl/2,3-pentanedione at low T and pH: thiazoles, FFT, thiophene-2-thiol; contexts cite the pathway only (also G6) | MAYBE |
| Marchand 2002 | Cysteine in aroma production in wine (N-(2-sulfanylethyl)-2-oxopropanamide) | [10.1021/jf025604w](https://doi.org/10.1021/jf025604w) | kw | 5 / 30 | M | intermediate of 2-acetylthiazole from MGO + cysteine (structure) | SKIP |
| Marchand 2011 | Cysteine + diacetyl under wine-like conditions | [10.1111/j.1750-3841.2011.02261.x](https://doi.org/10.1111/j.1750-3841.2011.02261.x) | kw | 10 / 23 | M | product list; origin of the methyl-C2 fragment | SKIP |
| Zeng 2005 | Adducts and S-carboxymethylcysteine from alpha-dicarbonyls + thiols | [10.1021/tx050074u](https://doi.org/10.1021/tx050074u) | kw | 87 / 157 | M | biomedical; adduct structures | SKIP |
| Hou 2025 | Exogenous cysteine oxidation induced browning via glyoxal (glycine-xylose MRPs) | [10.1016/j.foodchem.2025.147363](https://doi.org/10.1016/j.foodchem.2025.147363) | cites Liu2022, kw | 0 / 2 | ? | 2 citers, no contexts; cysteine oxidation as a glyoxal source (also G6) | MAYBE |
| Hamzalıoğlu 2019 | Broccoli and wine sulphur compounds on glyoxal scavenging (digestion) | [10.1016/j.jff.2019.02.033](https://doi.org/10.1016/j.jff.2019.02.033) | kw | 0 / 5 | ? | simulated digestion | SKIP |
| Cömert 2020 | MGO scavenging kinetics of food matrices (intestinal conditions) | [10.1007/s00217-020-03589-7](https://doi.org/10.1007/s00217-020-03589-7) | kw | 4 / 5 | B | 37 C digestion | SKIP |
| Huang 2024 | Thiols inhibit browning precursors in sorbitol/glycine | [10.1002/jsfa.13387](https://doi.org/10.1002/jsfa.13387) | cites Zhou2023, kw | 0 / 2 | ? | no contexts | SKIP |
| Noda 2016 | Pyrrolothiazolate from cysteine + glucose | [10.1039/c5fo01625h](https://doi.org/10.1039/c5fo01625h) | kw | 12 / 17 | M | pigment structure | SKIP |
| Revel 2003 | Maillard-type aroma in wine-like cysteine-carbonyl models | [10.1021/bk-2004-0871.ch026](https://doi.org/10.1021/bk-2004-0871.ch026) | kw | 2 / 10 | M | chapter version of Marchand 2000 | SKIP |

### G4 glutathione thermal degradation

Near papers: Zhou 2023 (acs.jafc.3c04424), Zheng & Ho 1994 (no OpenAlex references), Wang 2020 (foodchem.2019.125482).
**Verdict: nothing at two or more food-relevant temperatures.** The keyword pool is dominated by biomedical "glutathione" hits.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Liu 2020 | Preformulation of L-glutathione: degradation kinetics | [10.1080/03639045.2020.1752708](https://doi.org/10.1080/03639045.2020.1752708) | kw | 11 / 15 | Qa | forced-degradation kinetics; main route first-order oxidation to GSSG; seven degradation products by LC-MS; pharmaceutical conditions, temperatures not stated in abstract | MAYBE |
| Yamamoto 1994 (Yamamoto & Ishihara) | Stability of glutathione in solution | [10.1007/978-1-4615-2674-2_62](https://doi.org/10.1007/978-1-4615-2674-2_62) | kw | 4 / 9 | M | cited only for 'GSH is not stable in aqueous solution' (oxidation) | MAYBE |
| Lee 2010 | Aroma-active compounds from GSH + reducing sugars (AEDA, CAMOLA) | [10.1021/jf9043327](https://doi.org/10.1021/jf9043327) | cited by Wang2019, cited by Zhou2022, cited by Zhou2023; P1 rank 18 (14 corpus citers) | 29 / 68 | B | P1 rank 18: 14 corpus citers, all sulfur-lane; pathway labelling, no rates | MAYBE |
| Ho 1992 | Peptides as flavor precursors in model Maillard reactions | [10.1021/bk-1992-0490.ch015](https://doi.org/10.1021/bk-1992-0490.ch015) | kw | 6 / 21 | M | a citer: GSH-glucose releases H2S and NH3 more slowly than Glu + Cys + Gly; 180 C, pH 4-5 | MAYBE |
| Feng 2023 | Ribose-GSH Amadori product via limited oxidation and degradation of GSH | [10.1021/acs.jafc.3c05593](https://doi.org/10.1021/acs.jafc.3c05593) | kw | 3 / 6 | B | GSH oxidation/degradation grows at low ribose:GSH and high temperature (LC-MS, direction) | MAYBE |
| Lee 2011 | Metabolomics of key volatiles in GSH-Maillard reaction products | [10.1016/j.aca.2011.07.028](https://doi.org/10.1016/j.aca.2011.07.028) | kw | 16 / 41 | B | volatile-sensory correlation | SKIP |
| Hou 2025 | see G1 | [10.1016/j.foodchem.2025.144603](https://doi.org/10.1016/j.foodchem.2025.144603) | kw | 19 / 4 | Qa | GSH polysulfide decay at 25/35/40 C | MAYBE |
| Ueda 1997 | Flavor characteristics of glutathione in raw and cooked foodstuffs | [10.1271/bbb.61.1977](https://doi.org/10.1271/bbb.61.1977) | cited by Zhou2023, kw; P1 rank 368 (4 corpus citers) | 99 / 197 | B | kokumi taste; GSH contents | SKIP |

### G5 thiamine thermal degradation with H2S / volatile yields

Near papers: Dwivedi & Arnold 1972 (two papers), 1973; Hartman 1984 (two papers).
**Verdict: loss kinetics are well covered. For H2S yield, the only source is Dwivedi & Arnold 1971 ("Hydrogen sulfide from heat
degradation of thiamine"), already on disk as `arnold1971.pdf` with no DOI printed. Nothing newer was found.**

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Lee 1988 | pH, water activity, reducing sugars and kinetics of thermal thiamin breakdown (Iowa State thesis) | [10.31274/rtd-180813-12233](https://doi.org/10.31274/rtd-180813-12233) | cites Hartman1984, kw | 2 / 2 | Qa | first-order thiamin loss at 75, 85, 95 C in phosphate-buffered models and ground pork vs pH, a_w and reducing sugars | DOWNLOAD #13 |
| Fox 1982 | a_w, pH and heat on thiamine breakdown in foods | [10.1111/j.1745-4557.1982.tb00742.x](https://doi.org/10.1111/j.1745-4557.1982.tb00742.x) | cites Dwivedi1973 | 1 / 12 | Qa | rate constants at 85 and 95 C for pH 3, 5, 7 and a_w 0.90-1.0 | DOWNLOAD #14 |
| Grosch 1992 | see G1 | [10.1021/bk-1992-0490.ch014](https://doi.org/10.1021/bk-1992-0490.ch014) | kw; P1 rank 398 (4 corpus citers) | 2 / 12 | Qa | MFT yield from thiamine vs ribose + cysteine | DOWNLOAD #4 |
| Pachapurkar 2005 (Pachapurkar & Bell) | Kinetics of thiamin degradation in solutions at ambient storage | [10.1111/j.1365-2621.2005.tb11463.x](https://doi.org/10.1111/j.1365-2621.2005.tb11463.x) | cites Dwivedi1972, kw | 18 / 40 | Q | pH 4-7, phosphate vs citrate, 0.02/0.1 M at 25 C only; corpus has Voelker 2018/2021 (multi-T, multi-pH) | MAYBE |
| Voelker 2020 | Thiamine mononitrate kinetics in the aqueous phase of bread dough | [10.1016/j.foodres.2020.110084](https://doi.org/10.1016/j.foodres.2020.110084) | cites Dwivedi1973, kw | 1 / 10 | Q | dough matrix; the solution papers of the same group are in the corpus | SKIP |
| Morfee 1972 (Morfee & Liska) | Elemental sulfur: degradation product of thiamine in model systems | [10.3168/jds.s0022-0302(72)85441-9](https://doi.org/10.3168/jds.s0022-0302(72)85441-9) | cited by Dwivedi1972; P1 rank 553 (3 corpus citers) | 5 / 7 | M | 121 C / 40 min: elemental S a major product at basic pH (sulfur mass balance: S0 vs H2S) | MAYBE |
| Güntert 1992 | Thermally degraded thiamin | [10.1021/bk-1992-0490.ch011](https://doi.org/10.1021/bk-1992-0490.ch011) | kw; P1 rank 396 (4 corpus citers) | not in S2 | ? | not in S2; thiamin solutions at several concentrations and pH heated for several times; volatiles identified (amounts uncertain) | MAYBE |
| Tai 1999 | Thiamin oxidation and thermal formation of meaty aroma compounds | [10.1007/978-1-4615-4731-0_12](https://doi.org/10.1007/978-1-4615-4731-0_12) | cites Dwivedi1972, cites Dwivedi1973, kw | 1 / 3 | ? | 3 citers, 1 context | MAYBE |
| Arnold 1969 | Volatile flavor compounds from heat degradation of thiamine | [10.1021/jf60162a007](https://doi.org/10.1021/jf60162a007) | cited by Dwivedi1972, kw; P1 rank 112 (7 corpus citers) | 1 / 35 | B | identification | SKIP |
| Dwivedi 1973 | Minor volatile components from thermally degraded thiamine | [10.1111/j.1365-2621.1973.tb01451.x](https://doi.org/10.1111/j.1365-2621.1973.tb01451.x) | cited by Hartman1984, cites Dwivedi1972; P1 rank 455 (4 corpus citers) | 1 / 14 | B | identification (Dwivedi 1972/1973 JAFC in corpus) | SKIP |
| Ames 1992 | Novel volatile thermal degradation products of thiamine | [10.1002/jsfa.2740580220](https://doi.org/10.1002/jsfa.2740580220) | cites Hartman1984, kw | 0 / 6 | ? | identification | SKIP |
| Bendix 1951 | Factors influencing thiamine stability during heat sterilization | [10.1111/j.1365-2621.1951.tb17411.x](https://doi.org/10.1111/j.1365-2621.1951.tb17411.x) | cited by Dwivedi1972; P1 rank 226 (5 corpus citers) | 3 / 26 | Q | pH effect in foods; superseded | SKIP |
| Goulette 2019 | Thiamine degradation parameters in spaceflight foods (endpoints method) | [10.1016/j.foodchem.2019.125365](https://doi.org/10.1016/j.foodchem.2019.125365) | cites Dwivedi1973, kw | 2 / 14 | Q | food matrices | SKIP |

### G6 Strecker degradation of cysteine by dicarbonyls (H2S, NH3, mercaptoacetaldehyde)

Near papers: Liu 2022, Zhou 2022, Hofmann & Schieberle 1995 AEDA (jf00056a042).
**Verdict: no H2S yield or rate vs dicarbonyl, pH or temperature.** The keyword pool is mostly plant and cell biology of
H2S/MGO signalling.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Edwards 1992 | see G3 | [10.1080/02652039209374098](https://doi.org/10.1080/02652039209374098) | kw | 2 / 9 | Qa | slow thiol release after fast cysteine-3-DH reaction attributed to Strecker degradation; kinetic data | DOWNLOAD #2 |
| Mizutani 1960 | see G2 | [10.1080/03758397.1960.10857681](https://doi.org/10.1080/03758397.1960.10857681) | kw | 0 / 0 | Qa | measured H2S and NH3 from cysteine-aldehyde adducts at 120 C | DOWNLOAD #6 |
| Ho 1982 (Ho & Hartman) | Oxazolines/oxazoles in Strecker degradation of alanine and cysteine with 2,3-butanedione | [10.1021/jf00112a046](https://doi.org/10.1021/jf00112a046) | kw | 1 / 17 | M | 1-page communication; product identification | SKIP |
| Marchand 2000 | see G3 | [10.1021/jf000149u](https://doi.org/10.1021/jf000149u) | cited by Zhou2022, kw; P1 rank 605 (3 corpus citers) | 45 / 83 | M | cysteine + dicarbonyls at wine conditions | MAYBE |
| Majcher 2007 | Cysteine and cystine addition to extruded potato snacks | [10.1021/jf0703147](https://doi.org/10.1021/jf0703147) | cites Hofmann1995 | 14 / 41 | B | odorant changes with 0.25-1 % cysteine; process | SKIP |
| Mall 2018 | 2-alkylthiazolidine-4-carboxylic acids in stored orange juice | [10.1021/acs.jafc.8b04333](https://doi.org/10.1021/acs.jafc.8b04333) | cites Hofmann1995 | 7 / 6 | B | cysteine traps saturated aldehydes at 0 C (quantified, but juice) | SKIP |
| Hou 2025 | see G3 | [10.1016/j.foodchem.2025.147363](https://doi.org/10.1016/j.foodchem.2025.147363) | cites Liu2022, kw | 0 / 2 | ? | cysteine oxidation as glyoxal source | MAYBE |

### G7 cysteamine route to 2-acetyl-2-thiazoline

Near papers: Hofmann & Schieberle 1995 (jf00059a032, **not in corpus**), Cerny & Grosch 1993 (bf01190805, not in corpus).
**Verdict: one strong primary paper missing (Hofmann 1995), plus one yield paper (Engel 2002).** No kinetics of cysteine
decarboxylation to cysteamine was found. The Hofmann & Schieberle 1996 Weurman chapter on 2-oxopropanal/cysteamine intermediates
(10.1533/9781845698232.4.182) is already on disk as `hofmann1996b.pdf` (no DOI printed), and is removed from the table.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Hofmann 1995 (Hofmann & Schieberle) | Studies on the formation and stability of the roast-flavor compound 2-acetyl-2-thiazoline | [10.1021/jf00059a032](https://doi.org/10.1021/jf00059a032) | kw-Crossref | 23 / 32 | Q | NOT in the corpus (corpus hofmann1995.pdf is the AEDA paper, jf00056a042). Cysteamine + 2-oxopropanal and cysteine/sugar routes to 2-acetyl-2-thiazoline; citers quote 96 % loss of 2-acetyl-2-thiazoline after 60 min refluxing in water (100 C) | DOWNLOAD #1 |
| Engel 2002 (Engel & Schieberle) | Key aroma compounds from fructose + cysteamine or isothiaproline | [10.1021/jf0203186](https://doi.org/10.1021/jf0203186) | kw | 2 / 15 | Qa | fructose + cysteamine (model A) at 145 C, pH 7.0 in buffer: AEDA then quantitation (SIDA) of 5-acetyl-3,4-dihydro-2H-1,4-thiazine, 2-acetyl-2-thiazoline etc.; a direct cysteamine-route yield | DOWNLOAD #8 |
| Cerny 1993 (Cerny & Grosch) | Quantification of character-impact odour compounds of roasted beef | [10.1007/bf01190805](https://doi.org/10.1007/bf01190805) | manual lookup; P1 rank 587 (3 corpus citers) | 15 / 89 | B | not in corpus by DOI; SIDA amounts of 2-acetyl-2-thiazoline and others in roasted beef (food occurrence, not kinetics); citers use it for labelled standards | MAYBE |
| Huang 1998 | Thiazolidine and related volatiles in cysteamine/carbonyl systems | [10.1016/s0167-4501(98)80073-3](https://doi.org/10.1016/s0167-4501(98)80073-3) | cites Hofmann1995, kw | 0 / 0 | ? | 0 citers in S2; Weurman-type chapter on cysteamine + carbonyl mechanisms | MAYBE |
| Martin 1959 | Hydrolysis of 2-methyl-thiazoline; formation from N-acetyl-cysteamine | [10.1021/ja01528a018](https://doi.org/10.1021/ja01528a018) | kw, kw-Crossref | 0 / 56 | ? | 56 citers, no contexts, no abstract; by title, kinetics of 2-methyl-thiazoline hydrolysis in water (a stability analogue for 2-acetyl-2-thiazoline) | MAYBE |
| Rhlid 2002 | Roasted notes from 2-acetyl-2-thiazoline and its precursor 2-(1-hydroxyethyl)-4,5-dihydrothiazole | [10.1021/jf011170d](https://doi.org/10.1021/jf011170d) | cites Hofmann1995 | 15 / 12 | M | bioconversion route (yeast) | SKIP |
| Kerler 1997 | alpha-Acetyl-N-heterocycles in the Maillard reaction (review) | [10.1080/87559129709541140](https://doi.org/10.1080/87559129709541140) | cites Cerny1993, cites Hofmann1995, kw | 7 / 27 | B | review; thresholds | SKIP |
| Adams 2006 (Adams & De Kimpe) | Chemistry of 2-acetyl-1-pyrroline, 2-acetyl-2-thiazoline ... (review) | [10.1021/cr040097y](https://doi.org/10.1021/cr040097y) | cites Cerny1993, cites Hofmann1995, kw-Crossref; P1 rank 185 (6 corpus citers) | 17 / 167 | B | review; may point to primary yields | SKIP |

### G8 ribose-5-phosphate and nucleotides as precursors

Near papers: Mottram & Nobrega 2002 (jf0200826), Zhang & Ho 1991 IMP (jf00006a031, no OpenAlex references). Matoba et al.
1988 (IMP/GMP degradation vs pH and T, 10.1111/j.1365-2621.1988.tb13551.x) is already in the corpus.
**Verdict: dephosphorylation kinetics available; flavour-yield comparisons beyond the corpus are title-only leads.**

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Shaoul 1987 (Shaoul & Sporns) | Hydrolytic stability at intermediate pHs of IMP, GMP, AMP | [10.1111/j.1365-2621.1987.tb06733.x](https://doi.org/10.1111/j.1365-2621.1987.tb06733.x) | manual lookup | 2 / 12 | Q | 121 C, pH 3-8: half-life of IMP phosphate-bond hydrolysis 63 min at pH 5 (GMP 41, AMP 51); glycosidic cleavage too at pH 3; Arrhenius plot over several temperatures; citers: rate at pH 7 double that at pH 8 | DOWNLOAD #9 |
| Lee 1979 | Kinetics of inosine 5'-monophosphate instability | [10.1111/j.1365-2621.1979.tb08550.x](https://doi.org/10.1111/j.1365-2621.1979.tb08550.x) | kw-Crossref | 0 / 4 | Qa | 100 C, pH 5.52-6.87, ionic strength 0.4: pseudo-first-order IMP decomposition, log k-pH profile with a maximum near a pKa | DOWNLOAD #12 |
| Kuchiba 1990 (Kuchiba, Kaizaki, Matoba, Hasegawa) | Salts depress thermal degradation of IMP and GMP in aqueous solution | [10.1021/jf00093a001](https://doi.org/10.1021/jf00093a001) | kw, kw-Crossref | 1 / 5 | ? | same group as Matoba 1988 (in corpus); salt effect on the rate; 1 context | MAYBE |
| Sandwick 2005 | Maillard reactions of ribose 5-phosphate and amino acids | [10.1196/annals.1333.011](https://doi.org/10.1196/annals.1333.011) | kw | 26 / 35 | Qa | R5P reacts much faster than common sugars and sugar phosphates at pH 8, 37 C (UV/browning rates, not flavour) | MAYBE |
| Farmer 1996 | Three sugars and IMP as precursors of meat aroma (Weurman 1996) | [10.1533/9781845698232.4.225](https://doi.org/10.1533/9781845698232.4.225) | kw-Crossref | 3 / 9 | ? | exact gap topic by title; contexts mislinked; no abstract | MAYBE |
| Meinert 2009 | Glucose, glucose 6-phosphate, ribose, mannose as flavour precursors in pork | [10.1016/j.meatsci.2008.08.018](https://doi.org/10.1016/j.meatsci.2008.08.018) | kw-Crossref | 38 / 102 | B | pork matrix; one corpus citer | MAYBE |
| Vani 2005 | IMP degradation in aqueous and chicken muscle fibre systems | [10.1016/j.lwt.2005.05.003](https://doi.org/10.1016/j.lwt.2005.05.003) | kw | 50 / 51 | B | contexts are IMP contents of meat | SKIP |
| Kuchiba-Manabe 1991 | Umami of heated IMP solution (sensory) | [10.1111/j.1365-2621.1991.tb04790.x](https://doi.org/10.1111/j.1365-2621.1991.tb04790.x) | kw | 14 / 44 | B | sensory difference thresholds | SKIP |

### G9 thiol binding to melanoidins and proteins

Near papers: Hofmann & Schieberle 2002 (jf010823n), Charles-Bernard 2005 part 2 (jf048020y), Gigl 2021 (acs.jafc.1c06163).
**Verdict: coffee-melanoidin data only (mostly in corpus). No binding isotherm of MFT/FFT with any protein, and none with plant proteins.**
The one protein-thiol-disulfide exchange study is Mottram et al. 1996 JAFC (10.1021/jf960170o), already on disk with a dossier.

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Mottram 1996 (Mottram, Szauman-Szumski, Dodson) | Changes in thiol and disulfide flavour compounds resulting from the interaction with proteins (Weurman 1996) | [10.1533/9781845698232.7.413](https://doi.org/10.1533/9781845698232.7.413) | kw-Crossref | 1 / 7 | ? | conference version of Mottram, Szauman-Szumski & Dodson 1996 JAFC (10.1021/jf960170o), which is on disk (mottram1996.pdf, with dossier); no new content expected | SKIP |
| Hofmann 2001 | see G1 | [10.1021/jf0012042](https://doi.org/10.1021/jf0012042) | cited by Gigl2021, cited by Hofmann2001, cited by Kumazawa2003; P1 rank 150 (6 corpus citers) | 34 / 95 | Q | melanoidin-driven FFT loss | DOWNLOAD #5 |
| Sun 2024 | see G1 | [10.1016/j.lwt.2024.116930](https://doi.org/10.1016/j.lwt.2024.116930) | cites Blank2002, cites Charles-Bernard2005, cites Hofmann2001, cites Kumazawa2003, kw | 4 / 7 | Qa | melanoidin size fractions vs MFT/FFT/MeSH loss | DOWNLOAD #16 |
| Charles-Bernard 2005a | Interactions between volatile and nonvolatile coffee components. 1. Screening | [10.1021/jf048021q](https://doi.org/10.1021/jf048021q) | cited by Charles-Bernard2005, cites Charles-Bernard2005, cites Hofmann2001 | 25 / 60 | M | part 2 (jf048020y) is in the corpus; part 1 is the screening; citers quote tertiary > secondary > primary thiol stability (direction) | SKIP |
| Müller 2005 | see G1 | [10.1021/jf048027f](https://doi.org/10.1021/jf048027f) | cited by Gigl2021, cites Hofmann2001, kw | 5 / 26 | M | binding-site precursors | MAYBE |
| Jin 2026 | Interactions of food proteins with sulfur-containing flavor compounds (review) | [10.1016/j.tifs.2026.105960](https://doi.org/10.1016/j.tifs.2026.105960) | kw | 1 / 2 | B | 2026 review; may list protein-thiol binding constants and plant-protein data to chase | MAYBE |
| Negishi 2002 | Food materials removing Allium volatile sulfur compounds | [10.1021/jf020038q](https://doi.org/10.1021/jf020038q) | cited by Charles-Bernard2005 | 22 / 51 | Q | milk captures about 90 % of diallyl (di)sulfides; sulfides, not thiols | SKIP |

### G10 metal- and heme-catalysed thiol oxidation; heme in plant-based meat

Near papers: Baldus 2017 (acs.jafc.6b05472); the Impossible Foods patents have no DOI.
**Verdict: Cu kinetics yes (wine/biochemistry); Fe/heme-catalysed thiol oxidation kinetics no; heme x Maillard volatiles:
Meng 2025 only.**

| first author, year | title (short) | DOI | found by | S2 contexts / citers | class | what it offers (abstract / contexts) | rec. |
|---|---|---|---|---|---|---|---|
| Kreitman 2017 | Cu(II)-mediated H2S and thiol oxidation to disulfides and polysulfanes, and their reductive cleavage (model wine) | [10.1021/acs.jafc.6b05418](https://doi.org/10.1021/acs.jafc.6b05418) | kw-Crossref | 71 / 76 | Q | companion to Kreitman 2016a/b (in corpus): Cu(II)-mediated oxidation of H2S and thiols to disulfides/polysulfanes and recovery of H2S/MeSH from Cu-sulfhydryl complexes and disulfides by reductants; 71 contexts | DOWNLOAD #10 |
| Kachur 1999 | Mechanism of copper-catalyzed autoxidation of cysteine | [10.1080/10715769900300571](https://doi.org/10.1080/10715769900300571) | cited by Baldus2017, kw-Crossref | 43 / 107 | Q | two-phase O2-consumption kinetics of Cu-catalysed cysteine/NAC/GSH autoxidation; H2O2 and OH radical yields | DOWNLOAD #11 |
| Pecci 1997 | Novel findings on the copper-catalysed oxidation of cysteine | [10.1007/bf01372599](https://doi.org/10.1007/bf01372599) | cited by Baldus2017 | 66 / 117 | Q | Cu(I)-cysteine complex stabilisation; kinetics context | MAYBE |
| Rigo 2004 | Interaction of copper with cysteine: cuprous complex stability and catalytic role | [10.1016/j.jinorgbio.2004.06.008](https://doi.org/10.1016/j.jinorgbio.2004.06.008) | cited by Baldus2017 | 101 / 200 | M | kinetic scheme for Cu-catalysed cysteine oxidation (excess cysteine) | MAYBE |
| Hanaki 1995 | Kinetics of dioxygen in Cu-catalysed cysteine autoxidation | [10.1246/bcsj.68.831](https://doi.org/10.1246/bcsj.68.831) | kw-Crossref | 2 / 12 | B | 2 contexts (general) | SKIP |
| Meng 2025 | Leghemoglobin and myoglobin addition and the aroma of plant-based meat alternatives | [10.1016/j.foodchem.2025.145299](https://doi.org/10.1016/j.foodchem.2025.145299) | kw-Crossref | 3 / 10 | B | only heme-in-plant-meat volatile study found; 3 contexts, about colour (1 % LegH or Mb); amounts not visible | MAYBE |
| Álvarez-Armenta 2022 | Sulfmyoglobin production by free cysteine during thermal treatment (heme iron) | [10.1016/j.foodchem.2022.135165](https://doi.org/10.1016/j.foodchem.2022.135165) | kw | 1 / 14 | ? | heme iron + cysteine on heating; 1 off-topic context | MAYBE |
| Yancey 2006 | Iron, myoglobin, hemoglobin, lipid oxidation and livery off-flavor | [10.1016/j.meatsci.2006.03.013](https://doi.org/10.1016/j.meatsci.2006.03.013) | kw-Crossref | 57 / 90 | B | sensory correlations in beef | SKIP |
| Bekker 2016 | H2S formation in wine: copper and SO2 | [10.3390/molecules21091214](https://doi.org/10.3390/molecules21091214) | kw-Crossref | 20 / 41 | Q | wine storage | SKIP |

## Part 3: citation-context classes (union of Part 1 top 60 and Part 2 rows)

151 unique papers (60 from Part 1, 91 from Part 2 after merging papers listed under more than one gap): Q 35, Qa 16, M 25, B 55, ? 20. DOWNLOAD 16; the class of every paper is in its table row.

## DOWNLOAD list (ranked across gaps)

File names: lowercase `<firstauthor><year>.pdf`, with b/c/d added where that name is already in `data/articles`.

| # | paper | DOI | gap | class | why | file |
|---|---|---|---|---|---|---|
| 1 | Hofmann & Schieberle 1995, JAFC 43:2946, formation and stability of 2-acetyl-2-thiazoline | [10.1021/jf00059a032](https://doi.org/10.1021/jf00059a032) | G7 (G1) | Q | cysteamine + 2-oxopropanal route and 2-acetyl-2-thiazoline stability (96 % lost after 60 min at 100 C, per citers); not the corpus `hofmann1995.pdf` | `hofmann1995b.pdf` |
| 2 | Edwards & Wedzicha 1992, Food Addit. Contam., 3-deoxyhexosulose + thiols | [10.1080/02652039209374098](https://doi.org/10.1080/02652039209374098) | G3, G6, G2 | Qa | kinetics and stoichiometry of thiols with 3-DH; cysteine: fast thiazolidine then slow (Strecker) release | `edwards1992.pdf` |
| 3 | Lo et al. 1994, J. Biol. Chem., MGO with N-acetyl amino acids | [10.1016/s0021-9258(18)31635-1](https://doi.org/10.1016/s0021-9258(18)31635-1) | G3 | Q | forward and reverse rate constants of thiol + MGO hemithioacetal (pH 7.4, 37 C); 7 corpus citers | `lo1994.pdf` |
| 4 | Grosch & Zeiler-Hilgart 1992, ACS Symp. Ser. 490 ch. 14 | [10.1021/bk-1992-0490.ch014](https://doi.org/10.1021/bk-1992-0490.ch014) | G1, G5 | Qa | MFT, MFT disulfide and FFT amounts from thiamine vs ribose + cysteine at meat-level concentrations; cysteine/H2S effects | `grosch1992.pdf` |
| 5 | Hofmann, Czerny, Calligaris & Schieberle 2001, JAFC, coffee melanoidins and volatiles | [10.1021/jf0012042](https://doi.org/10.1021/jf0012042) | G1, G9 | Q | FFT loss time course with melanoidins (50 % at 20 min, about 100 % at 30 min) | `hofmann2001b.pdf` |
| 6 | Mizutani 1960, Bull. Agric. Chem. Soc. Japan, cysteine-aldehyde compounds | [10.1080/03758397.1960.10857681](https://doi.org/10.1080/03758397.1960.10857681) | G6, G2 | Qa | measured H2S and NH3 from cysteine-aldehyde thiazolidines heated at 120 C (furfural adduct highest) | `mizutani1960.pdf` |
| 7 | de Roos 1992, ACS Symp. Ser. 490 ch. 16 | [10.1021/bk-1992-0490.ch016](https://doi.org/10.1021/bk-1992-0490.ch016) | G2 | Qa | reactivity vs neutral thiazolidinecarboxylic-acid concentration (pH dependence of the TCA step); precursor of de Roos 2005 | `deroos1992.pdf` |
| 8 | Engel & Schieberle 2002, JAFC, fructose + cysteamine / isothiaproline | [10.1021/jf0203186](https://doi.org/10.1021/jf0203186) | G7 | Qa | SIDA yields of 2-acetyl-2-thiazoline and acetyl-thiazines from cysteamine + fructose, 145 C, pH 7 | `engel2002.pdf` |
| 9 | Shaoul & Sporns 1987, J. Food Sci., hydrolytic stability of IMP/GMP/AMP | [10.1111/j.1365-2621.1987.tb06733.x](https://doi.org/10.1111/j.1365-2621.1987.tb06733.x) | G8 | Q | half-lives at 121 C, pH 3-8, plus Arrhenius data; complements Matoba 1988 (in corpus) | `shaoul1987.pdf` |
| 10 | Kreitman et al. 2017, JAFC, Cu(II)-mediated H2S and thiol oxidation | [10.1021/acs.jafc.6b05418](https://doi.org/10.1021/acs.jafc.6b05418) | G10 | Q | third paper of the Kreitman series (2016a/b in corpus): disulfide/polysulfane formation and H2S recovery | `kreitman2017.pdf` |
| 11 | Kachur et al. 1999, Free Radic. Res., Cu-catalysed cysteine autoxidation | [10.1080/10715769900300571](https://doi.org/10.1080/10715769900300571) | G10 | Q | two-phase kinetics of Cu-catalysed cysteine/NAC/GSH oxidation | `kachur1999.pdf` |
| 12 | Lee 1979, J. Food Sci., kinetics of IMP instability | [10.1111/j.1365-2621.1979.tb08550.x](https://doi.org/10.1111/j.1365-2621.1979.tb08550.x) | G8 | Qa | pseudo-first-order IMP loss at 100 C, pH 5.5-6.9 | `lee1979.pdf` |
| 13 | Lee 1988, Iowa State thesis, thermal thiamin breakdown | [10.31274/rtd-180813-12233](https://doi.org/10.31274/rtd-180813-12233) | G5 | Qa | first-order thiamin loss at 75/85/95 C vs pH, a_w, reducing sugars, buffer and pork (open repository) | `lee1988.pdf` |
| 14 | Fox 1982, J. Food Qual., a_w/pH/heat and thiamine breakdown | [10.1111/j.1745-4557.1982.tb00742.x](https://doi.org/10.1111/j.1745-4557.1982.tb00742.x) | G5 | Qa | k at 85 and 95 C x pH 3/5/7 x a_w 0.90-1.0 | `fox1982.pdf` |
| 15 | Müller & Hofmann 2007, JAFC, phenol/FFT conjugates in coffee | [10.1021/jf070095p](https://doi.org/10.1021/jf070095p) | G1 | Qa | 60 % of FFT conjugated within 10 min; SIDA time course (coffee) | `muller2007.pdf` |
| 16 | Sun et al. 2024, LWT, coffee-matrix aroma binding factors | [10.1016/j.lwt.2024.116930](https://doi.org/10.1016/j.lwt.2024.116930) | G1, G9 | Qa | MFT to 34 %, FFT to 3 % in 30 min; loss apportioned to lipids, HHQ, melanoidin fractions | `sun2024.pdf` |

Ranking logic: first the two items the sulfur lane can use directly as a mechanism step with a rate or yield (1, 2), then the
G3/G1 rate data (3-5), then G6/G2/G7 yield data (6-8), then G8/G10 kinetics (9-12), then G5 loss kinetics the corpus partly
covers (13-14), then coffee-only thiol-sink numbers (15-16).

**MAYBE, in order:** Wedzicha & Edwards 1991 ([10.1016/0308-8146(91)90020-o](https://doi.org/10.1016/0308-8146(91)90020-o),
`wedzicha1991.pdf`, fetch together with #2). Fan 2011, glucose/L-cysteine kinetics (Chinese, no DOI). Rae 1990 (`rae1990.pdf`).
Kraehenbuehl 2005 (ASIC, no DOI). Meng 2025, leghemoglobin/myoglobin in plant-based meat (`meng2025.pdf`). Farmer 1996, Weurman,
sugars vs IMP (`farmer1996.pdf`). Liu 2020, GSH degradation
(`liu2020.pdf`). Kuchiba 1990 (`kuchiba1990.pdf`). Martin 1959, thiazoline hydrolysis (`martin1959.pdf`). Lee 2010, GSH + sugars
(`lee2010.pdf`). Hou 2025, polysulfides (`hou2025.pdf`). Morfee 1972 (`morfee1972.pdf`). Cerny & Grosch 1993 (`cerny1993.pdf`).
Sandwick 2005 (`sandwick2005.pdf`). The trunk/acrylamide MAYBEs from Part 1 are Ajandouz 2001, Leong 2000, Brands 2003,
Hofmann 1999, Laroque 2008, Thornalley 1999 and Mottram 2002 Nature.
