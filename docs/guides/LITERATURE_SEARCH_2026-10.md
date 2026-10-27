# Literature search 2026-10: citation chaining for sulfur / meat-flavour formation

Run 2026-10-27. Purpose: measure and close the coverage gap exposed when a Wageningen expert's 24-paper list on meat-flavour
volatiles from amino acids, sugars, thiamine and glutathione overlapped our corpus in only 8 of 24. Our original search came
from LLM deep-research digests and keyword queries, which favour well-cited kinetics papers; the misses were mostly mechanism,
intermediate and labelling papers. This note is the search record only; nothing was added to `data/`.

## Method

- **Source**: OpenAlex REST API (`api.openalex.org`, anonymous, no key; list queries are rate-limited to 1000 credits/day).
  Lookups by DOI (`/works/doi:...`), backward links from each seed's `referenced_works`, forward links from
  `/works?filter=cites:<id>` (all citing works; per-seed cap 600, never reached, maximum 412), metadata for references via
  `filter=openalex_id:` batches. About 100 API calls in total.
- **Seeds A (our core, 9 DOIs)**: Hofmann & Schieberle 1998 (10.1021/jf9705983) and 1998b (10.1021/jf970990g), Whitfield & Mottram
  1999 (10.1021/jf980980v), Kang 2026 (10.1039/d5fb00932d), Cerny & Briffod 2007 (10.1021/jf062874w), Zhou 2023
  (10.1021/acs.jafc.2c08360), Yiltirak 2026 (10.1016/j.foodres.2026.118600), Meynier & Mottram 1995
  (10.1016/0308-8146(95)93282-V), Mottram & Nobrega 2002 (10.1021/jf0200826). OpenAlex has **no reference list and no citers for
  Kang 2026**, and one citer for Yiltirak 2026, so in practice the A set has 8 working seeds.
- **Seeds B (expert list, 24 DOIs)**: the references in `scratch/wur_blind_test/paper for MeatCODE.docx` (Wang 2012 appears twice in
  that list; 24 unique DOIs). All 24 resolved in OpenAlex.
- **Score** = number of distinct seeds a candidate is directly linked to (cites the seed, or is cited by it).
  Ties are broken by **coupling** (number of seed references the candidate also cites; only works with a reference list)
  and then by citation count. Table columns show `score/coupling`.
- **Topic filter**: keep a candidate only if its title or OpenAlex abstract matches at least one of: cystein(e), cystin(e),
  thiamin(e), glutathione, ribose, xylose, pentose, furanthiol, furfurylthiol, meat flavo(u)r / meaty, sulfur / sulphur,
  thiazolidine, Amadori, Maillard, H2S / hydrogen sulfide, thiophene, disulfide. (Many Elsevier abstracts are absent from
  OpenAlex, so for those the filter runs on the title alone.)
- **Corpus exclusion**: every DOI in `data/keys/papers.yml` (case-insensitive), plus the 24 expert papers. See the registry
  gap below for papers the DOI test misses.

## Recall test (seeds A only)

The question: had we run this method from our own core before the expert wrote, how many of the 16 expert papers we did
**not** have would it have surfaced? The topic-filtered, corpus-excluded candidate pool from seeds A has 261 works.

| top N | expert papers found (of 16) | recall |
|---|---|---|
| 50 | 6 | 38% |
| 100 | 7 | 44% |
| 200 | 8 | 50% |
| whole pool (261) | 9 | 56% |

Not excluding the corpus (pool 296) gives the same counts (6/7/8); only the ranks move down.

| rank in A-pool | score/coupling | expert paper | DOI |
|---|---|---|---|
| 3 | 4/8 | Martins 2010: Meat Flavor Generation in Complex Maillard Model Systems | [10.1021/bk-2010-1042.ch008](https://doi.org/10.1021/bk-2010-1042.ch008) |
| 8 | 3/9 | Farmer 1989: Volatile compounds produced in maillard reactions involving cysteine, ribose and phospholi | [10.1002/jsfa.2740490311](https://doi.org/10.1002/jsfa.2740490311) |
| 10 | 3/8 | Zhao 2018: Meat flavor generation from different composition patterns of initial Maillard stage inter | [10.1016/j.foodchem.2018.08.096](https://doi.org/10.1016/j.foodchem.2018.08.096) |
| 25 | 2/8 | Modi 2008 (OpenAlex lists Linforth first): Effect of pH and Water Activity in Generation of Selected Meaty Aroma Compounds in a Meat  | [10.3923/ajft.2008.68.78](https://doi.org/10.3923/ajft.2008.68.78) |
| 29 | 2/6 | Cao 2017: Effect of glycine on reaction of cysteine-xylose: Insights on initial Maillard stage inter | [10.1016/j.foodres.2017.06.012](https://doi.org/10.1016/j.foodres.2017.06.012) |
| 31 | 2/5 | Xu 2013: Formation mechanism of volatile and non-volatile compounds in peptide–xylose Maillard reac | [10.1016/j.foodres.2013.07.066](https://doi.org/10.1016/j.foodres.2013.07.066) |
| 77 | 1/6 | Zhang 2024: Effect of lysine on the cysteine-xylose Maillard reaction to form flavor compounds | [10.1016/j.foodchem.2024.142529](https://doi.org/10.1016/j.foodchem.2024.142529) |
| 107 | 1/4 | Cerny 2013: Formation of cysteine-S-conjugates in the Maillard reaction of cysteine and xylose | [10.1016/j.foodchem.2013.04.043](https://doi.org/10.1016/j.foodchem.2013.04.043) |
| 210 | 1/1 | Garbusov 1976: Volatile sulfur compounds contributing to meat flavour. Part. I. Components identified in  | [10.1002/food.19760200302](https://doi.org/10.1002/food.19760200302) |
| not reached | - | Sun 2021: Maillard reaction of food-derived peptides as a potential route to generate meat flavor co | [10.1016/j.foodres.2021.110823](https://doi.org/10.1016/j.foodres.2021.110823) |
| not reached | - | Zheng 1997: Formation of Sulfur-Containing Flavor Compounds from Reactions of Furaneol and Cysteine, G | [10.1021/jf960624h](https://doi.org/10.1021/jf960624h) |
| not reached | - | Sohail 2022: Aroma compounds identified in cooked meat: A review | [10.1016/j.foodres.2022.111385](https://doi.org/10.1016/j.foodres.2022.111385) |
| not reached | - | Zhai 2023: Intrinsic Molecular Mechanisms of Transformation between Isomeric Intermediates Formed at  | [10.1021/acs.jafc.3c06506](https://doi.org/10.1021/acs.jafc.3c06506) |
| not reached | - | Zhai 2019: Interaction of Added l -Cysteine with 2-Threityl-thiazolidine-4-carboxylic Acid Derived fr | [10.1021/acs.jafc.9b04374](https://doi.org/10.1021/acs.jafc.9b04374) |
| not reached | - | Du 2024: Effects of unsaturated C18 fatty acids on “glucose-glutathione” Maillard reaction: Compari | [10.1016/j.foodres.2024.115645](https://doi.org/10.1016/j.foodres.2024.115645) |
| not reached | - | Cha 2019: Analysis of α-dicarbonyl compounds and volatiles formed in Maillard reaction model systems | [10.1038/s41598-019-41824-8](https://doi.org/10.1038/s41598-019-41824-8) |

Reading the result:

- **One hop finds 9 of 16; the top 50 already holds 6.** The seven misses have no direct citation link to any working A seed.
  A second forward hop (works citing the top-50 A candidates) would reach four more (Du 2024/25, Zhai 2019, Sun 2021/22, Sohail
  2022: each cites 2 to 9 of the top 50), for a ceiling of about 13/16. Zhai 2023 cites one of the top 200, Cha 2019 and
  Zheng 1997 none; they need other seeds (Zheng 1997 is a furaneol + cysteine/H2S paper, reachable from the furaneol side).
- **Precision is low in absolute terms**: 6 expert hits in a top 50 means most of the top 50 are other papers. The full-search
  triage below suggests roughly one in five top-ranked papers is class A, so a reader screening the top 150 by title finds
  about 30 quantitative candidates.
- Both seed sets are small and partly overlapping (Cerny & Briffod 2007 and Meynier 1995 are in A and in the expert list), so
  these are rough numbers: one-hop recall is about 40-50% at N = 50-200, with a ceiling near 80% if a second hop is added.

## Registry gap found on the way

`data/keys/papers.yml` misses DOIs for papers we do have, so a DOI-only dedupe flags them as new:

- 3 of the 8 expert papers we count as known are not in the registry by DOI: Meynier & Mottram 1995, Wang 2012
  (dossier `wang2012_extraction.md` is untracked), Cerny 2015 (book chapter).
- 7 of the full-search top 165 are corpus papers with a dossier but no registry DOI (they are in the
  registry's `dossiers_without_printed_doi` list). Checked by hand against dossier headers and removed from the tables:
  Hofmann 1995 ([10.1021/jf00056a042](https://doi.org/10.1021/jf00056a042), `hofmann1995 dossier`);
  Farmer 1990 ([10.1002/jsfa.2740530409](https://doi.org/10.1002/jsfa.2740530409), `farmer1990 dossier`);
  Xu 2010 ([10.17221/189/2008-cjfs](https://doi.org/10.17221/189/2008-cjfs), `xu2010 dossier`);
  Whitfield 1988 ([10.1002/jsfa.2740420309](https://doi.org/10.1002/jsfa.2740420309), `whitfield1988 dossier`);
  Shu 1988 ([10.1021/jf00082a032](https://doi.org/10.1021/jf00082a032), `shu1988 dossier`);
  Hofmann 1996 ([10.1021/jf9500703](https://doi.org/10.1021/jf9500703), `hofmann1996 dossier`);
  Desclaux 2006 ([10.1016/s0167-4501(06)80087-7](https://doi.org/10.1016/s0167-4501(06)80087-7), `desclaux2006 dossier`).
- Name matches that turned out to be **different** papers (kept as candidates): Zhang 2023 DMDS/DMTS (10.1016/j.foodres.2023.113200
  vs dossier 10.1021/acs.jafc.3c02450), Zhou 2024 (10.1021/acs.jafc.4c00021 vs dossier 10.1021/acs.jafc.4c03706), Sun 2019,
  Schieberle 2005, Mottram 1995 (the cysteine-ribose-phospholipid paper vs the dione + H2S dossier).

## Full search (seeds A and B)

Topic-filtered pool after dropping the corpus and the 24 expert papers: **924 works**. Score distribution:
score 1: 625, score 2: 160, score 3: 54, score 4: 45, score 5: 17, score 6: 14, score 7: 4, score 8: 2, score 9: 2, score 11: 1.

The top 165 were triaged by title and abstract; 7 were corpus papers (above), leaving **158 triaged**:

| class | meaning | count |
|---|---|---|
| A | likely quantitative data on S-volatiles, thiols, H2S, ARP, TTCA, thiamine or GSH degradation with stated conditions | 29 |
| B | mechanism or isotope labelling | 22 |
| C | review, book chapter, identification-only, or food/hydrolysate application study with no reusable model step | 91 |
| X | off-topic | 16 |

Triage is from titles and abstracts only; "likely quantitative" has to be confirmed in the full text before anything is extracted.
Many abstracts were missing, so some class C papers may hold usable numbers. Ranks below are positions in the corpus-excluded
full ranking.

### Class A: likely quantitative, could anchor or test a model step

| rank | score/coupling | first author, year | title | journal | DOI | why |
|---|---|---|---|---|---|---|
| 1 | 9/24 | Zhao 2018 | Formation mechanism of aroma compounds in a glutathione-glucose reaction with fat or oxidized fat | Food Chemistry | [10.1016/j.foodchem.2018.07.106](https://doi.org/10.1016/j.foodchem.2018.07.106) | GSH-glucose +/- fat/oxidised fat: S-volatiles and intermediates in a GSH pot |
| 2 | 9/15 | Zhang 2020 | Effects of phospholipids and reheating treatment on volatile compounds in phospholipid-xylose-cysteine reaction systems | Food Research International | [10.1016/j.foodres.2020.109918](https://doi.org/10.1016/j.foodres.2020.109918) | xylose-cysteine-phospholipid pot incl. reheating: S-volatile amounts vs lipid |
| 7 | 7/16 | Wang 2019 | Characterization of initial reaction intermediates in heated model systems of glucose, glutathione, and aliphatic aldehydes | Food Chemistry | [10.1016/j.foodchem.2019.125482](https://doi.org/10.1016/j.foodchem.2019.125482) | glucose-GSH(+aldehyde) initial intermediates (ARP/thiazolidines) quantified |
| 21 | 5/20 | Zhou 2023 | Promotion or Inhibition Effects of Exogenous Glutathione-Degraded Amino Acids on the Formation of 2,3-Butanedione and Pyrazines via Varied Pathways of Interaction with α-Dicarbonyl Compounds Derived from N-(1-Deoxy- d -xylulos-1-yl)-alanine | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.3c04424](https://doi.org/10.1021/acs.jafc.3c04424) | Ala-xylose ARP + GSH-derived Glu/Cys/Gly: butanedione and pyrazine yields (seed-Zhou 2023 lab) |
| 24 | 5/18 | Zhou 2022 | Competitive Formation of 2,3-Butanedione and Pyrazines through Intervention of Added Cysteine during Thermal Processing of Alanine-Xylose Amadori Compounds | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.2c07026](https://doi.org/10.1021/acs.jafc.2c07026) | added cysteine vs Ala-xylose ARP: competitive butanedione/pyrazine formation, dose series |
| 25 | 5/18 | Du 2024 | The number and position of unsaturated bonds in aliphatic aldehydes affect meat flavorings system: Insights on initial Maillard reaction stage and meat flavor formation from thiazolidine derivatives | Current Research in Food Science | [10.1016/j.crfs.2024.100719](https://doi.org/10.1016/j.crfs.2024.100719) | cysteine-Amadori and thiazolidine derivatives quantified by UPLC-MS/MS vs aldehyde type |
| 27 | 5/15 | Elmore 2002 | Effect of Lipid Composition on Meat-like Model Systems Containing Cysteine, Ribose, and Polyunsaturated Fatty Acids | Journal of Agricultural and Food Chemistry | [10.1021/jf0108718](https://doi.org/10.1021/jf0108718) | cysteine-ribose + methyl linoleate/linolenate: quantified S-heterocycles (Mottram lab) |
| 28 | 5/14 | Li 2021 | The Maillard Reaction as Source of Meat Flavor Compounds in Dry Cured Meat Model Systems under Mild Temperature Conditions | Molecules | [10.3390/molecules26010223](https://doi.org/10.3390/molecules26010223) | meat model systems at mild (dry-cure) temperatures: Maillard volatiles quantified |
| 30 | 5/12 | Madruga 1995 | The effect of pH on the formation of maillard‐derived aroma volatiles using a cooked meat system | Journal of the Science of Food and Agriculture | [10.1002/jsfa.2740680308](https://doi.org/10.1002/jsfa.2740680308) | Madruga & Mottram 1995: pH 4.0-5.6 series in cooked minced beef, S-volatiles; matrix twin of Meynier 1995 |
| 36 | 5/0 | Shu 1985 | pH Effect on the volatile components in the thermal degradation of cysteine | Journal of Agricultural and Food Chemistry | [10.1021/jf00063a029](https://doi.org/10.1021/jf00063a029) | Shu 1985: thermal degradation of cysteine at pH 2.2/5.1/7.1, volatile yields - an H2S-source step |
| 46 | 4/21 | Xu 2008 | Effects of ribose to cysteine ratios on the formation of volatile compounds from the Maillard reaction in supercritical carbon dioxide | Food Research International | [10.1016/j.foodres.2008.05.005](https://doi.org/10.1016/j.foodres.2008.05.005) | ribose:cysteine ratio series (SC-CO2 vs control), quantified volatiles; sibling of xu2010 |
| 47 | 4/19 | Xu 2011 | Effects of lard on the formation of volatiles from the Maillard reaction of cysteine with xylose | Journal of the Science of Food and Agriculture | [10.1002/jsfa.4445](https://doi.org/10.1002/jsfa.4445) | cysteine-xylose +/- lard/oxidised lard: headspace S-volatiles |
| 54 | 4/17 | Sun 2019 | Aqueous Preparation of Maillard Reaction Intermediate from Glutathione and Xylose and its Volatile Formation During Thermal Treatment | Journal of Food Science | [10.1111/1750-3841.14911](https://doi.org/10.1111/1750-3841.14911) | GSH-xylose ARP (MGX) prepared, then heated: volatile formation from the ARP |
| 55 | 4/17 | Du 2024 | Effect of 2-alkenals on the Maillard reaction of cysteine-glucose: Initial stage intermediate formation and flavor compound generation pathways | LWT | [10.1016/j.lwt.2024.115741](https://doi.org/10.1016/j.lwt.2024.115741) | 2-alkenals in cysteine-glucose: cystine-Amadori/thiazolidine quantified (UPLC-MS/MS) + volatiles |
| 56 | 4/16 | Liu 2022 | Formation Priority of Pyrazines and 2-Acetylthiazole Dependent on the Added Cysteine and Fragments of Deoxyosones during the Thermal Process of the Glycine–Ribose Amadori Compound | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.2c04874](https://doi.org/10.1021/acs.jafc.2c04874) | Gly-Ribose ARP + cysteine dose: 2-acetylthiazole vs pyrazines, deoxyosone fragments |
| 60 | 4/14 | Xu 2007 | Effects of supercritical carbon dioxide on volatile formation from Maillard reaction between ribose and cysteine | Journal of the Science of Food and Agriculture | [10.1002/jsfa.3093](https://doi.org/10.1002/jsfa.3093) | ribose-cysteine 140 C under SC-CO2 vs control: volatiles quantified (xu2010 sibling) |
| 66 | 4/11 | Yang 2015 | Aromatic effect of fat and oxidized fat on a meat‐like model reaction system of cysteine and glucose | Flavour and Fragrance Journal | [10.1002/ffj.3248](https://doi.org/10.1002/ffj.3248) | cysteine-glucose +/- fat/oxidised fat: 76 volatiles with amounts |
| 78 | 4/0 | Hartman 1984 | Volatile products formed from the thermal degradation of thiamin at high and low moisture levels | Journal of Agricultural and Food Chemistry | [10.1021/jf00125a016](https://doi.org/10.1021/jf00125a016) | Hartman 1984: thiamine thermal degradation volatiles at high vs low moisture |
| 79 | 4/0 | Mottram 1995 | Volatile Compounds from the Reaction of Cysteine, Ribose, and Phospholipid in Low-Moisture Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf00052a027](https://doi.org/10.1021/jf00052a027) | Mottram & Whitfield 1995: cysteine-ribose-phospholipid at low moisture, volatiles quantified |
| 81 | 3/32 | Hincelin 1992 | The effect of xylose on the generation of volatiles from heated thiamine | Food Chemistry | [10.1016/0308-8146(92)90272-4](https://doi.org/10.1016/0308-8146(92)90272-4) | Hincelin 1992: xylose effect on heated thiamine volatiles |
| 87 | 3/21 | Cannon 2023 | Quantitation of Key Odorants in Arginine/Cysteine-Glucose Maillard Reactions | Journal of Food Bioactives | [10.31665/jfb.2023.18338](https://doi.org/10.31665/jfb.2023.18338) | Cannon 2023: quantitation of key odorants in arginine/cysteine-glucose reactions |
| 95 | 3/16 | Zhou 2024 | Exogenous Alanine Promoting the Reaction between Amadori Compound and Deoxyxylosone and Inhibiting the Formation of 2-Furfural during Thermal Treatment | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.4c00021](https://doi.org/10.1021/acs.jafc.4c00021) | Ala-xylose ARP + exogenous Ala: deoxyxylosone and 2-furfural yields (no S; ARP step) |
| 102 | 3/13 | Wei 2022 | Glycine-Xylose Amadori Compound Formation Tracing through Maillard Browning Inhibition by 2-Threityl-thiazolidine-4-carboxylic Acid Formation from Deoxyosone and Exogenous Cysteine | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.2c04961](https://doi.org/10.1021/acs.jafc.2c04961) | Gly-xylose ARP traced via TTCA formation from deoxyosone + exogenous cysteine |
| 104 | 3/12 | Yaghmur 2002 | Furfural−Cysteine Model Reaction in Food Grade Nonionic Oil/Water Microemulsions for Selective Flavor Formation | Journal of Agricultural and Food Chemistry | [10.1021/jf011158l](https://doi.org/10.1021/jf011158l) | furfural + cysteine at 65 C in microemulsions: FFT formation rates |
| 109 | 3/11 | Zhang 2023 | Factors influencing on the formation of dimethyl disulfide and dimethyl trisulfide in model systems | Food Research International | [10.1016/j.foodres.2023.113200](https://doi.org/10.1016/j.foodres.2023.113200) | factors controlling DMDS/DMTS formation in model systems (sibling of zhang2023 dossier, different paper) |
| 112 | 3/10 | Lee 2011 | Determination of Volatile Sulfur Compounds Formed by the Maillard Reaction of Glutathione with Glucose | ACS symposium series | [10.1021/bk-2011-1068.ch011](https://doi.org/10.1021/bk-2011-1068.ch011) | volatile S-compounds from GSH + glucose (GC-MS; check for amounts) |
| 149 | 2/17 | Tang 2019 | N -(1-Deoxy- d -xylulos-1-yl)-glutathione: A Maillard Reaction Intermediate Predominating in Aqueous Glutathione-Xylose Systems by Simultaneous Dehydration-Reaction | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.9b04694](https://doi.org/10.1021/acs.jafc.9b04694) | GSH-xylose ARP yield and 3-DX/1-DX kinetics (reaction kinetics applied) |
| 156 | 2/16 | Liu 2023 | Effect of Aliphatic Aldehydes on Flavor Formation in Glutathione–Ribose Maillard Reactions | Foods | [10.3390/foods12010217](https://doi.org/10.3390/foods12010217) | nine aldehyde-GSH-ribose pots: browning, intermediates, volatiles |
| 157 | 2/16 | Hartman 1984 | Effect of Water Activity on the Major Volatileś Produced in a Model System Approximating Cooked Meat | Journal of Food Science | [10.1111/j.1365-2621.1984.tb12480.x](https://doi.org/10.1111/j.1365-2621.1984.tb12480.x) | Hartman 1984: aw effect on volatiles of glutamate/ascorbate/thiamine/cystine meat model at 135 C |

### Class B: mechanism / labelling

| rank | score/coupling | first author, year | title | journal | DOI | why |
|---|---|---|---|---|---|---|
| 14 | 6/15 | Alim 2019 | Research of beef-meaty aroma compounds from yeast extract using carbon module labeling (CAMOLA) technique | LWT | [10.1016/j.lwt.2019.06.006](https://doi.org/10.1016/j.lwt.2019.06.006) | CAMOLA labelling of beef-meaty odorants from yeast extract |
| 15 | 6/14 | Ze 2019 | Studies on the Origin of Carbons in Aroma Compounds from [13C6]Glucose -Cysteine-(E)-2-Nonenal Model Reaction Systems | Polymers | [10.3390/polym11030521](https://doi.org/10.3390/polym11030521) | [13C6]glucose-cysteine-(E)-2-nonenal carbon-origin labelling |
| 17 | 6/6 | Hofmann 1997 | Identification of Potent Aroma Compounds in Thermally Treated Mixtures of Glucose/Cysteine and Rhamnose/Cysteine Using Aroma Extract Dilution Techniques | Journal of Agricultural and Food Chemistry | [10.1021/jf960456t](https://doi.org/10.1021/jf960456t) | Hofmann & Schieberle 1997: glucose/rhamnose-cysteine 145 C/20 min AEDA, FFT/thiazine routes; protocol twin of our Hofmann 1998 holdouts (no dossier; cited only in the quarantined-benchmarks README and an audit note) |
| 18 | 6/2 | Schieberle 2005 | The Carbon Module Labeling (CAMOLA) Technique: A Useful Tool for Identifying Transient Intermediates in the Formation of Maillard‐Type Target Molecules | Annals of the New York Academy of Sciences | [10.1196/annals.1333.029](https://doi.org/10.1196/annals.1333.029) | CAMOLA method paper (Schieberle 2005) |
| 23 | 5/19 | Lei 2024 | Analysis of volatiles from the thermal decomposition of Amadori rearrangement products in the cysteine-glucose Maillard reaction and density functional theory study | Food Research International | [10.1016/j.foodres.2024.114454](https://doi.org/10.1016/j.foodres.2024.114454) | volatiles from cysteine-glucose ARP thermal decomposition (+ DFT; DFT part unusable under policy) |
| 26 | 5/17 | Cui 2017 | Controlled formation of flavor compounds by preparation and application of Maillard reaction intermediate (MRI) derived from xylose and phenylalanine | RSC Advances | [10.1039/c7ra09355a](https://doi.org/10.1039/c7ra09355a) | xylose-phenylalanine MRI preparation (no S) |
| 31 | 5/12 | Yu 2011 | Mechanism of formation of sulphur aroma compounds from l-ascorbic acid and l-cysteine during the Maillard reaction | Food Chemistry | [10.1016/j.foodchem.2011.11.111](https://doi.org/10.1016/j.foodchem.2011.11.111) | S-aroma formation mechanism from ascorbic acid + cysteine |
| 52 | 4/17 | Zou 2018 | Discovery of Amadori‐Type Conjugates in a Peptide Maillard Reaction and Their Corresponding Influence on the Formation of Pyrazines | Journal of Food Science | [10.1111/1750-3841.14156](https://doi.org/10.1111/1750-3841.14156) | Amadori-type peptide conjugates, CAMOLA, pyrazines |
| 59 | 4/14 | Wang 2020 | Reinvestigation of 2-acetylthiazole formation pathways in the Maillard reaction | Food Chemistry | [10.1016/j.foodchem.2020.128761](https://doi.org/10.1016/j.foodchem.2020.128761) | 2-acetylthiazole formation pathways re-investigated (labelling) |
| 61 | 4/13 | Cui 2018 | Improved controlled flavor formation during heat-treatment with a stable Maillard reaction intermediate derived from xylose-phenylalanine | Food Chemistry | [10.1016/j.foodchem.2018.07.161](https://doi.org/10.1016/j.foodchem.2018.07.161) | xylose-Phe MRI (no S) |
| 64 | 4/13 | Yang 2012 | Response Surface Methodology for Meat‐like Odorants from Maillard Reaction with Glutathione I: the Optimization Analysis and the General Pathway Exploration | Journal of Food Science | [10.1111/j.1750-3841.2012.02863.x](https://doi.org/10.1111/j.1750-3841.2012.02863.x) | GSH Maillard RSM optimisation + general pathway |
| 68 | 4/11 | Zhao 2025 | Unraveling the formation pathways of key meat-flavor Spirocyclic sulfur compound from thiamin-cysteine model reaction: A combined experimental and computational study | Food Chemistry | [10.1016/j.foodchem.2025.147221](https://doi.org/10.1016/j.foodchem.2025.147221) | thiamine-cysteine spirocyclic S-compound formation (expt + computation) |
| 72 | 4/6 | Huang 2012 | Inhibiting the color formation by gradient temperature‐elevating Maillard reaction of soybean peptide‐xylose system based on interaction of l ‐cysteine and Amadori compounds | Journal of Peptide Science | [10.1002/psc.2406](https://doi.org/10.1002/psc.2406) | cysteine-Amadori interaction inhibiting colour (soy peptide-xylose) |
| 73 | 4/6 | Ota 2006 | Characterization of a New Maillard Type Reaction Product Generated by Heating 1-Deoxymaltulosyl-glycine in the Presence of Cysteine | Journal of Agricultural and Food Chemistry | [10.1021/jf060296d](https://doi.org/10.1021/jf060296d) | new 1-deoxyosone-cysteine intermediate from an Amadori compound |
| 77 | 4/0 | Guentert 1990 | Identification and formation of some selected sulfur-containing flavor compounds in various meat model systems | Journal of Agricultural and Food Chemistry | [10.1021/jf00101a007](https://doi.org/10.1021/jf00101a007) | Guentert 1990: formation of selected S-compounds in meat model systems |
| 80 | 4/0 | de Roos 2005 | Amadori Compounds of Cysteine and Their Role in the Development of Meat Flavor | ACS symposium series | [10.1021/bk-2005-0905.ch009](https://doi.org/10.1021/bk-2005-0905.ch009) | de Roos 2005: cysteine Amadori/thiazine intermediates in meat flavour |
| 90 | 3/18 | Yang 2012 | The mechanism of peptide bonds cleavage and volatile compounds generated from pentapeptide to heptapeptide via Maillard reaction | Food Chemistry | [10.1016/j.foodchem.2012.01.044](https://doi.org/10.1016/j.foodchem.2012.01.044) | peptide-bond cleavage and volatiles, penta-/heptapeptide MR |
| 91 | 3/18 | Feng 2024 | Directed preparation of N-(1-deoxy-d-ribulos-1-yl)-glutathione during short-term high temperature dehydration: Suppressing the side reactions of glutathione | Food Bioscience | [10.1016/j.fbio.2024.104468](https://doi.org/10.1016/j.fbio.2024.104468) | directed preparation of ribulosyl-GSH ARP |
| 93 | 3/17 | Zhang 2021 | Accelerated Dissipation of Free and Immobilized Water Facilitating the Intramolecular Dehydration of N-Xylosamine and Conversion Improvement of the Amadori Rearrangement Product of Aspartic Acid–Xylose Reaction | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.1c05827](https://doi.org/10.1021/acs.jafc.1c05827) | Asp-xylose ARP conversion by dehydration (no S) |
| 126 | 3/4 | Nikolov 2012 | Role of the Ribose-Specific Marker Furfuryl-amine in the Formation of Aroma Active 1-(Furan-2-ylmethyl)-1H-pyrrole (or Furfuryl-pyrrole) Derivatives | Journal of Agricultural and Food Chemistry | [10.1021/jf3030002](https://doi.org/10.1021/jf3030002) | ribose-specific furfurylamine marker, furfurylpyrrole formation |
| 129 | 3/0 | Shibamoto 1977 | Formation of sulfur- and nitrogen-containing compounds from the reaction of furfural with hydrogen sulfide and ammonia | Journal of Agricultural and Food Chemistry | [10.1021/jf60209a055](https://doi.org/10.1021/jf60209a055) | Shibamoto 1977: furfural + H2S + NH3 products |
| 158 | 2/15 | Cui 2018 | Synergistic Effect of a Thermal Reaction and Vacuum Dehydration on Improving Xylose–Phenylalanine Conversion to N-(1-Deoxy- d -xylulos-1-yl)-phenylalanine during an Aqueous Maillard Reaction | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.8b04448](https://doi.org/10.1021/acs.jafc.8b04448) | xylose-Phe ARP conversion (no S) |

## Reading the candidate lists

- The class A list clusters in three families. (1) **Mottram/Whitfield-era pots** that are siblings of seeds we already
  use: Madruga & Mottram 1995 (pH series in cooked beef), Mottram & Whitfield 1995 (low-moisture cysteine-ribose-phospholipid),
  Elmore 2002 (PUFA), Shu 1985 (cysteine degradation vs pH), Hartman 1984 x2 and Hincelin 1992 (thiamine). These are the most
  likely to give tabulated amounts. (2) **ARP and thiazolidine papers from the Jiangnan (Zhang/Cui/Ho) and BTBU (Xie) groups** in the same series as the
  Zhou 2023 seed and the expert's Zhai papers: quantified ARP/TTCA/deoxyosone kinetics. (3) **GSH-sugar systems**
  (Zhao 2018, Wang 2019, Sun 2019, Tang 2019, Lee 2011, Liu 2023), the glutathione lane the expert list emphasised.
- Class B is where the expert list's character lives (labelling and intermediates): CAMOLA papers, de Roos 2005 on cysteine
  Amadori/thiazine intermediates, Ota 2006, Hofmann & Schieberle 1997 (same 145 C / 20 min protocol as our Hofmann 1998
  holdouts).
- Reproduce: the scripts (`fetch.py`, `rank.py`, `triage.py`, `write_doc.py`) and the API cache were kept in the session
  scratchpad, not in the repo.

## Round 2 (2026-10-27)

Same day, after the expert papers and the round-1 class-A downloads were read. Question: does a second citation hop from a
larger seed set recover what one hop missed, and what new quantitative papers does it surface? Search record only; nothing
was added to `data/`.

### Method (changes from round 1)

- **Seeds (46 DOIs)**: the 9 core seeds A, the 24 expert papers, and the 15 round-1 class-A downloads read today
  (Zhao 2018/2019 10.1016/j.foodchem.2018.07.106, Zhang 2020/2021 10.1016/j.foodres.2020.109918, Wang 2019/2020
  10.1016/j.foodchem.2019.125482, Zhou 2023 10.1021/acs.jafc.3c04424, Zhou 2022 10.1021/acs.jafc.2c07026, Du 2024
  10.1016/j.crfs.2024.100719, Elmore 2002 10.1021/jf0108718, Madruga 1995 10.1002/jsfa.2740680308, Shu 1985 10.1021/jf00063a029,
  Xu 2008 10.1016/j.foodres.2008.05.005, Sun 2019 10.1111/1750-3841.14911, Liu 2022 10.1021/acs.jafc.2c04874, Hofmann 1997
  10.1021/jf960456t, Fan 2018 10.1021/acs.jafc.8b03297, Guth 1994 10.1021/jf00048a039). Cerny & Briffod 2007 and Meynier 1995
  are in both A and the expert list, so 48 listed seeds are 46 unique; all 46 resolved in OpenAlex. The dossiers with today's
  timestamp that are not in this list (the 21:48 and 14:39 batches: arnold1971 ... suleimenov1997, hidalgo2004, zamora2007/2008)
  are unmodified in git (timestamp only) and are not sulfur-lane reads, so they were not used as seeds.
- **Hop 1** (backward + forward) from every seed, as in round 1 (forward cap 600 citers per seed).
- **Hop 2**: the top 60 hop-1 candidates (topic-filtered, seeds excluded, corpus papers kept because they are valid bridges)
  were expanded both ways: their reference lists and their citers (forward cap 400 per intermediary, most-cited first).
- **Score** = number of distinct seeds linked directly + 0.5 x number of further distinct seeds reached through a hop-2
  intermediary (a seed counted directly is not counted again at 0.5). Ties: coupling (shared seed references), then citations.
  Tables show `score (direct + 2nd-hop seeds)`.
- Topic filter unchanged from round 1, for comparability.
- **Exclusions**: DOIs in `data/keys/papers.yml` (case-insensitive); the first DOI printed on page 1 of each `data/articles` PDF
  (`pdftotext -l 1`; 271 of 366 PDFs yield one, mostly post-2000 papers); every DOI in this document above; every DOI triaged
  in round 1 including the unlisted class C/X ones (83 such fell in the top 300); and two corpus papers caught by an
  author-year clash check against dossier stems (marked IN below).
- About 147 new OpenAlex calls on top of the round-1 cache; anonymous, no email or key sent.

### Recall check (core seeds A only, two hops)

Same 16 target expert papers as round 1; corpus excluded except the targets and today's class-A downloads (so the pool is
what the method would have shown before today). Hop-2 intermediaries = top 60 of the core-only hop-1 pool.

| procedure | pool | top 50 | top 100 | top 200 | anywhere in pool |
|---|---|---|---|---|---|
| one hop (round 1 recomputed; pool was 261 then) | 259 | 6 | 7 | 8 | 9 |
| two hops | 1503 | 6 | 8 | 11 | 13 |

| two-hop rank | one-hop rank | score | expert paper | DOI |
|---|---|---|---|---|
| 1 | 2 | 5.5 (4+3) | Meat Flavor Generation in Complex Maillard Model Systems (2010) | [10.1021/bk-2010-1042.ch008](https://doi.org/10.1021/bk-2010-1042.ch008) |
| 5 | 7 | 5 (3+4) | Volatile compounds produced in maillard reactions involving cysteine, ribose and phospholi (1989) | [10.1002/jsfa.2740490311](https://doi.org/10.1002/jsfa.2740490311) |
| 6 | 9 | 5 (3+4) | Meat flavor generation from different composition patterns of initial Maillard stage inter (2018) | [10.1016/j.foodchem.2018.08.096](https://doi.org/10.1016/j.foodchem.2018.08.096) |
| 17 | 28 | 4.5 (2+5) | Effect of glycine on reaction of cysteine-xylose: Insights on initial Maillard stage inter (2017) | [10.1016/j.foodres.2017.06.012](https://doi.org/10.1016/j.foodres.2017.06.012) |
| 20 | 30 | 4.5 (2+5) | Formation mechanism of volatile and non-volatile compounds in peptide–xylose Maillard reac (2013) | [10.1016/j.foodres.2013.07.066](https://doi.org/10.1016/j.foodres.2013.07.066) |
| 42 | 106 | 4 (1+6) | Formation of cysteine-S-conjugates in the Maillard reaction of cysteine and xylose (2013) | [10.1016/j.foodchem.2013.04.043](https://doi.org/10.1016/j.foodchem.2013.04.043) |
| 61 | not reached | 3.5 (0+7) | Interaction of Added l -Cysteine with 2-Threityl-thiazolidine-4-carboxylic Acid Derived fr (2019) | [10.1021/acs.jafc.9b04374](https://doi.org/10.1021/acs.jafc.9b04374) |
| 89 | not reached | 3.5 (0+7) | Formation of Sulfur-Containing Flavor Compounds from Reactions of Furaneol and Cysteine, G (1997) | [10.1021/jf960624h](https://doi.org/10.1021/jf960624h) |
| 108 | 24 | 3 (2+2) | Effect of pH and Water Activity in Generation of Selected Meaty Aroma Compounds in a Meat  (2008) | [10.3923/ajft.2008.68.78](https://doi.org/10.3923/ajft.2008.68.78) |
| 115 | 76 | 3 (1+4) | Effect of lysine on the cysteine-xylose Maillard reaction to form flavor compounds (2024) | [10.1016/j.foodchem.2024.142529](https://doi.org/10.1016/j.foodchem.2024.142529) |
| 116 | not reached | 3 (0+6) | Maillard reaction of food-derived peptides as a potential route to generate meat flavor co (2021) | [10.1016/j.foodres.2021.110823](https://doi.org/10.1016/j.foodres.2021.110823) |
| 229 | not reached | 2.5 (0+5) | Effects of unsaturated C18 fatty acids on “glucose-glutathione” Maillard reaction: Compari (2024) | [10.1016/j.foodres.2024.115645](https://doi.org/10.1016/j.foodres.2024.115645) |
| 949 | 209 | 1 (1+0) | Volatile sulfur compounds contributing to meat flavour. Part. I. Components identified in  (1976) | [10.1002/food.19760200302](https://doi.org/10.1002/food.19760200302) |
| not reached | not reached | - | Analysis of α-dicarbonyl compounds and volatiles formed in Maillard reaction model systems (2019) | [10.1038/s41598-019-41824-8](https://doi.org/10.1038/s41598-019-41824-8) |
| not reached | not reached | - | Aroma compounds identified in cooked meat: A review (2022) | [10.1016/j.foodres.2022.111385](https://doi.org/10.1016/j.foodres.2022.111385) |
| not reached | not reached | - | Intrinsic Molecular Mechanisms of Transformation between Isomeric Intermediates Formed at  (2023) | [10.1021/acs.jafc.3c06506](https://doi.org/10.1021/acs.jafc.3c06506) |

Reading the result:

- **The second hop reaches 13 of 16 (81%), as round 1 predicted, but only inside a pool six times larger.** The top-50 count
  does not move (6); the gain is at depth: 8 in the top 100, 11 in the top 200. The four newly reached papers (Zhai 2019,
  Zheng 1997, Sun 2021/22, Du 2024/25) enter at ranks 61-229. Cha 2019, Sohail 2022 and Zhai 2023 stay unreachable from the
  core: they neither cite nor are cited by anything within two hops of it.
- **The 0.5-weighted hop-2 term dominates the score.** A candidate linked to a hub intermediary inherits every seed of that
  hub, so reviews and food-application papers that cite many hubs rise. Direct links still order the targets that were already
  found (they move up a little: Martins 2010 rank 2 -> 1, Cerny 2013 106 -> 42), but the added depth is mostly noise.
- Practical rule: screen the one-hop list to about 100, then use the second hop only to rescue specific known-missing papers
  or families; a second hop is a recall tool at N = 200, not a precision tool at N = 50.

### Full search (46 seeds, two hops)

Topic-filtered two-hop pool: 2172 works; minus 105 registry DOIs, 2 PDF-page-1 DOIs and
43 DOIs already in this document leaves 2022; the round-1 triaged DOIs were then removed. The top 217 of what
remains were triaged by title and abstract (scores 18.5 down to 9.0); 2 turned out to be corpus papers, leaving 215 triaged:

| class | meaning | count |
|---|---|---|
| A | likely quantitative S-volatiles, thiols, H2S, ARP/TTCA, thiamine, GSH or cysteine kinetics with stated conditions | 18 |
| B | mechanism or isotope labelling | 40 |
| C | review, book chapter, identification-only, or food/hydrolysate application | 135 |
| X | off-topic | 22 |

Corpus papers found by the author-year check (dossier without a registry DOI): rank 51 [10.1007/bf01198411](https://doi.org/10.1007/bf01198411) (nedvidek1992 dossier); rank 97 [10.1021/jf00049a035](https://doi.org/10.1021/jf00049a035) (mottram1995 dossier (alkanediones + H2S)).

Class A density is about 8% (vs about 18% in round 1's top 158): the second hop mostly adds papers that cite the hubs.
Triage is from titles and abstracts only; many old JAFC abstracts are absent from OpenAlex, so "likely quantitative" for
the 1988-1998 papers rests on the known style of those labs and has to be confirmed in the full text. Suggested file names
follow `<firstauthor><year>.pdf` (OpenAlex first author and year), with b/c added where the name is taken in `data/articles`
or by an existing dossier stem.

#### Round 2 Class A: likely quantitative, could anchor or test a model step

| rank | score | first author, year | title | journal | DOI | why | file |
|---|---|---|---|---|---|---|---|
| 5 | 16.5 (4+25) | Zhang 1989 | Volatile compounds formed from thermal interaction of 2,4-decadienal with cysteine and glutathione | Journal of Agricultural and Food Chemistry | [10.1021/jf00088a044](https://doi.org/10.1021/jf00088a044) | Ho lab: 2,4-decadienal + cysteine or GSH heated in water; S-heterocycle yields (lipid-aldehyde x thiol step) | `zhang1989.pdf` |
| 6 | 16 (3+26) | Chen 2002 | Effects of Carnosine on Volatile Generation from Maillard Reaction of Ribose and Cysteine | Journal of Agricultural and Food Chemistry | [10.1021/jf011244l](https://doi.org/10.1021/jf011244l) | Ho lab: ribose-cysteine pot +/- carnosine; quantified S-volatiles (antioxidant perturbation of the trunk) | `chen2002.pdf` |
| 13 | 15 (4+22) | Yu 2010 | Aroma compounds generated from thermal reaction of l-ascorbic acid with l-cysteine | Food Chemistry | [10.1016/j.foodchem.2010.01.049](https://doi.org/10.1016/j.foodchem.2010.01.049) | ascorbic acid + cysteine pot: S-aroma compounds vs conditions (Yu lab; check for amounts) | `yu2010.pdf` |
| 20 | 15 (2+26) | Zhang 1988 | Comparison of the volatile compounds obtained from thermal degradation of cysteine and glutathione in water | Journal of Agricultural and Food Chemistry | [10.1021/jf00083a022](https://doi.org/10.1021/jf00083a022) | Ho lab: thermal degradation of cysteine vs GSH in water; volatile (H2S-derived) yields - an H2S-source step | `zhang1988.pdf` |
| 59 | 12.5 (0+25) | Chen 2000 | Effect of Urea on Volatile Generation from Maillard Reaction of Cysteine and Ribose | Journal of Agricultural and Food Chemistry | [10.1021/jf991076l](https://doi.org/10.1021/jf991076l) | Ho lab: cysteine-ribose pot +/- urea; quantified S-volatiles | `chen2000.pdf` |
| 61 | 12.5 (2+21) | Hofmann 1998 | Identification of key aroma compounds generated from cysteine and carbohydrates under roasting conditions | Zeitschrift für Lebensmittel-Untersuchung und -Forschung. B, Referate und Lebensmittelrecht | [10.1007/s002170050324](https://doi.org/10.1007/s002170050324) | Hofmann & Schieberle 1998 Z Lebensm Unters Forsch: cysteine + carbohydrates under roasting (dry) - key odorants, FFT/MFT; dry-heat twin of our Hofmann 1998 holdouts | `hofmann1998c.pdf` |
| 62 | 12 (2+20) | Lu 2019 | Timely Addition of Glutathione for Its Interaction with Deoxypentosone To Inhibit the Aqueous Maillard Reaction and Browning of Glycylglycine–Arabinose System | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.9b02053](https://doi.org/10.1021/acs.jafc.9b02053) | GSH (vs GSSG) trapping deoxypentosone in Gly-Gly-arabinose; timed addition, ARP/browning kinetics (Zhang/Cui lab) | `lu2019.pdf` |
| 65 | 12 (1+22) | Zhang 2022 | The potential application of vegetable oils in the D-xylose and L-cysteine Maillard reaction system for meaty aroma production | Food Research International | [10.1016/j.foodres.2022.111081](https://doi.org/10.1016/j.foodres.2022.111081) | xylose-cysteine pot + vegetable oils: meaty S-volatiles vs oil (sibling of round-1 Xu 2011 lard paper; check amounts) | `zhang2022.pdf` |
| 68 | 12 (1+22) | Dreher 2003 | GC−Olfactometric Characterization of Aroma Volatiles from the Thermal Degradation of Thiamin in Model Orange Juice | Journal of Agricultural and Food Chemistry | [10.1021/jf034023j](https://doi.org/10.1021/jf034023j) | thiamine degradation in model orange juice, 35 C up to 8 weeks: aroma volatiles over time (low-T thiamine step) | `dreher2003.pdf` |
| 69 | 12 (7+10) | El-Massry 2003 | Volatile constituents of glutathione – ribose model system and its antioxidant activity | Amino Acids | [10.1007/s00726-002-0318-4](https://doi.org/10.1007/s00726-002-0318-4) | El-Massry 2003: GSH-ribose pot volatiles (+ antioxidant activity); GSH lane, check for amounts | `elmassry2003.pdf` |
| 78 | 12 (1+22) | Zhang 1991 | Formation of meatlike aroma compounds from thermal reaction of inosine 5'-monophosphate with cysteine and glutathione | Journal of Agricultural and Food Chemistry | [10.1021/jf00006a031](https://doi.org/10.1021/jf00006a031) | Ho lab: IMP (ribose-phosphate source) + cysteine or GSH; meat-like volatile yields | `zhang1991.pdf` |
| 90 | 11.5 (2+19) | Farmer 1992 | Effect of cysteine and ribose on the volatile thermal degradation products of a triglyceride and three phospholipids | Journal of the Science of Food and Agriculture | [10.1002/jsfa.2740600414](https://doi.org/10.1002/jsfa.2740600414) | Farmer & Mottram 1992: cysteine + ribose effect on triglyceride/phospholipid degradation volatiles, quantified (sibling of farmer1989/1990) | `farmer1992.pdf` |
| 96 | 11.5 (2+19) | Mottram 1995 | Maillard-Lipid Interactions in Nonaqueous Systems: Volatiles from the Reaction of Cysteine and Ribose with Phosphatidylcholine | Journal of Agricultural and Food Chemistry | [10.1021/jf00053a033](https://doi.org/10.1021/jf00053a033) | Mottram & Whitfield 1995: cysteine-ribose-phosphatidylcholine in nonaqueous systems, quantified volatiles (not the round-1 low-moisture paper) | `mottram1995b.pdf` |
| 109 | 11 (1+20) | Mandin 1999 | Volatile Compounds from Potato-like Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf981277+](https://doi.org/10.1021/jf981277+) | Ames lab: methionine-glucose-linoleic acid-starch potato models; methional/S-volatile amounts (methionine, not cysteine) | `mandin1999.pdf` |
| 139 | 10.5 (3+15) | Zhang 1991 | Comparison of the volatile compounds formed from the thermal reaction of glucose with cysteine and glutathione | Journal of Agricultural and Food Chemistry | [10.1021/jf00004a029](https://doi.org/10.1021/jf00004a029) | Ho lab: glucose + cysteine vs glucose + GSH heated; S-volatile yields (hexose x thiol) | `zhang1991b.pdf` |
| 158 | 10 (1+18) | Parker 2013 | Comparison of Ribose and Ascorbic Acid in Model Process Reactions Containing Cysteine | Elsevier eBooks | [10.1016/b978-0-12-398549-1.00040-4](https://doi.org/10.1016/b978-0-12-398549-1.00040-4) | Parker 2013 Weurman: ribose vs ascorbic acid with cysteine, process-reaction volatiles (Reading lab) | `parker2013b.pdf` |
| 161 | 10 (2+16) | Schieberle 1996 | IDENTIFICATION OF THE KEY ODORANTS IN PROCESSED RIBOSE-CYSTEINE MAILLARD MIXTURES BY INSTRUMENTAL ANALYSIS AND SENSORY STUDIES | Elsevier eBooks | [10.1533/9781845698232.4.175](https://doi.org/10.1533/9781845698232.4.175) | Schieberle & Hofmann 1996 Weurman: key odorants in processed ribose-cysteine, instrumental quantitation (twin of hofmann1995) | `schieberle1996.pdf` |
| 188 | 9.5 (3+13) | Guth 1995 | Model Reactions on the Stability of Disulfides in Heated Foods | Journal of Agricultural and Food Chemistry | [10.1021/jf00056a044](https://doi.org/10.1021/jf00056a044) | Guth et al. 1995 (Grosch group): stability of flavour disulfides/thiols (FFT, MFT) in heated model foods - a thiol-loss step | `guth1995.pdf` |

#### Round 2 Class B: mechanism / labelling

| rank | score | first author, year | title | journal | DOI | why | file |
|---|---|---|---|---|---|---|---|
| 2 | 17.5 (1+33) | Wang 2008 | Amino acid-dependent formation pathways of 2-acetylfuran and 2,5-dimethyl-4-hydroxy-3[2H]-furanone in the Maillard reaction | Food Chemistry | [10.1016/j.foodchem.2008.12.014](https://doi.org/10.1016/j.foodchem.2008.12.014) | Wang & Ho: amino-acid-dependent pathways to 2-acetylfuran and furaneol (labelling; no S) | `wang2008b.pdf` |
| 12 | 15 (2+26) | Zhao 2026 | Unraveling the formation mechanism of pyrazines and their derivatives in the glucose-glycine Maillard reaction using carbon module labeling coupled with carbonyl analysis | Current Research in Food Science | [10.1016/j.crfs.2026.101552](https://doi.org/10.1016/j.crfs.2026.101552) | CAMOLA + carbonyl trapping: pyrazines in glucose-glycine (no S) | `zhao2026.pdf` |
| 15 | 15 (1+28) | Smuda 2013 | Fragmentation Pathways during Maillard-Induced Carbohydrate Degradation | Journal of Agricultural and Food Chemistry | [10.1021/jf305117s](https://doi.org/10.1021/jf305117s) | Glomb lab: carbohydrate fragmentation pathways, dicarbonyl yields (no S) | `smuda2013.pdf` |
| 17 | 15 (1+28) | Brehm 2019 | Novel Taste-Enhancing 4-Amino-2-methyl-5-heteroalkypyrimidines Formed from Thiamine by Maillard-Type Reactions | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.9b05896](https://doi.org/10.1021/acs.jafc.9b05896) | thiamine-derived taste-active pyrimidines by Maillard-type reactions (sibling of brehm2020 dossier) | `brehm2019.pdf` |
| 27 | 14.5 (3+23) | Schieberle 2003 | The carbohydrate modul labelling technique - a useful tool to clarify formation pathways of aroma compounds formed in Maillard-type reactions. | Flavour Research at the Dawn of the Twenty-first Century - Proceedings of the 10th Weurman Flavour Research Symposium, Beaune, France, 25-28 June, 2002. | - | Schieberle 2003 Weurman proceedings: CAMOLA technique | `schieberle2003.pdf` |
| 28 | 14 (2+24) | Zou 2019 | Flavour precursor peptide from an enzymatic beef hydrolysate Maillard reaction-II: Mechanism of the synthesis of flavour compounds from a sulphur-containing peptide through a Maillard reaction | LWT | [10.1016/j.lwt.2019.04.022](https://doi.org/10.1016/j.lwt.2019.04.022) | S-containing beef peptide in Maillard: mechanism of flavour formation | `zou2019.pdf` |
| 34 | 13.5 (4+19) | Xu 2019 | Preparation of N ‐(1‐Deoxy‐Α‐D‐Xylulos‐1‐Yl)‐Glutamic Acid via Aqueous Maillard Reaction Coupled with Vacuum Dehydration and Its Flavor Formation Through Thermal Treatment of Baking Process | Journal of Food Science | [10.1111/1750-3841.14733](https://doi.org/10.1111/1750-3841.14733) | Glu-xylose ARP preparation and flavour on baking (no S) | `xu2019.pdf` |
| 39 | 13 (4+18) | Zhang 2024 | Control Formation of Furans and Pyrazines Resulting from Dual Glycation Sites in Nα,Nε-Di(1-deoxy- d -xylulos-1-yl)lysine via Elevating Thermal Degradation Temperatures | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.4c07733](https://doi.org/10.1021/acs.jafc.4c07733) | di-xylulosyl-lysine ARP thermal degradation: furans/pyrazines (no S) | `zhang2024c.pdf` |
| 43 | 13 (1+24) | Xia 2023 | Simultaneously Enhanced Formation of Pyrazines and Furans during Thermal Degradation of the Glycyl- l -glutamine Amadori Compound by Selected Exogenous Amino Acids and Appropriate Elevated Temperatures | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.3c00085](https://doi.org/10.1021/acs.jafc.3c00085) | Gly-Gln ARP degradation + exogenous amino acids (no S) | `xia2023.pdf` |
| 48 | 13 (2+22) | Liu 2026 | Underlying substituted pyrazines formation in the glucose-alanine Maillard reaction: Carbon module labeling and carbonyl analysis strategy | Food Research International | [10.1016/j.foodres.2026.119409](https://doi.org/10.1016/j.foodres.2026.119409) | CAMOLA + carbonyl analysis: pyrazines in glucose-alanine (no S) | `liu2026.pdf` |
| 49 | 13 (1+24) | Wang 2008 | Comparison of 2-Acetylfuran Formation between Ribose and Glucose in the Maillard Reaction | Journal of Agricultural and Food Chemistry | [10.1021/jf802683a](https://doi.org/10.1021/jf802683a) | Wang & Ho: 2-acetylfuran from ribose vs glucose, reactivity and labelling (no S) | `wang2008c.pdf` |
| 50 | 13 (1+24) | Pripis‐Nicolau 2000 | Formation of Flavor Components by the Reaction of Amino Acid and Carbonyl Compounds in Mild Conditions | Journal of Agricultural and Food Chemistry | [10.1021/jf991024w](https://doi.org/10.1021/jf991024w) | alpha-dicarbonyls + amino acids incl. cysteine/methionine at mild (wine) conditions: S-heterocycles | `pripisnicolau2000.pdf` |
| 60 | 12.5 (0+25) | Totlani 2005 | Reactivity of Epicatechin in Aqueous Glycine and Glucose Maillard Reaction Models: Quenching of C2, C3, and C4 Sugar Fragments | Journal of Agricultural and Food Chemistry | [10.1021/jf050044x](https://doi.org/10.1021/jf050044x) | 13C/15N labelling: epicatechin quenching C2-C4 fragments (no S) | `totlani2005.pdf` |
| 80 | 11.5 (3+17) | Chen 2023 | Exogenous Threonine-Induced Conversion of Threonine-Xylose Amadori Compound to Heyns Compound for Efficiently Promoting the Formation of Pyrazines | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.3c02311](https://doi.org/10.1021/acs.jafc.3c02311) | Thr-xylose ARP -> Heyns conversion, pyrazines (no S) | `chen2023.pdf` |
| 82 | 11.5 (1+21) | Xu 2014 | Thermal reaction of anisaldehyde in the presence of l-cysteine, a model reaction of Chinese stew meat flavor generation | - | [10.7282/t3mp51kh](https://doi.org/10.7282/t3mp51kh) | thesis: anisaldehyde + cysteine model reaction (thiol/thiazolidine chemistry) | `xu2014.pdf` |
| 83 | 11.5 (1+21) | Zhang 2023 | Controlled Selective Formation of Amadori Compounds from α/ε Mono- or Di-glycation of Lysine with Xylose | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.3c00388](https://doi.org/10.1021/acs.jafc.3c00388) | selective mono/di-glycation of lysine by xylose: ARP formation (no S) | `zhang2023c.pdf` (renamed 2026-10-28: `zhang2023b.pdf` is the DMDS/DMTS paper, 10.1016/j.foodres.2023.113200) |
| 85 | 11.5 (1+21) | Weenen 1998 | Reactive intermediates and carbohydrate fragmentation in Maillard chemistry | Food Chemistry | [10.1016/s0308-8146(98)00074-0](https://doi.org/10.1016/s0308-8146(98)00074-0) | Weenen 1998: reactive intermediates and carbohydrate fragmentation (dicarbonyls) | `weenen1998.pdf` |
| 89 | 11.5 (4+15) | Elmore 2000 | Formation of 2-Alkyl-(2H)-thiapyrans and 2-Alkylthiophenes in Cooked Beef and Lamb | Journal of Agricultural and Food Chemistry | [10.1021/jf000124v](https://doi.org/10.1021/jf000124v) | Elmore & Mottram 2000: 2-alkylthiapyrans/alkylthiophenes in cooked beef/lamb, H2S + dienal route | `elmore2000.pdf` |
| 95 | 11.5 (2+19) | Buttery 1984 | Thiamin odor and bis(2-methyl-3-furyl) disulfide | Journal of Agricultural and Food Chemistry | [10.1021/jf00123a061](https://doi.org/10.1021/jf00123a061) | Buttery 1984: thiamine odour = bis(2-methyl-3-furyl) disulfide (MFT dimer identification) | `buttery1984.pdf` |
| 99 | 11 (1+20) | Deng 2024 | Promoted formation of pyrazines by targeted precursor addition to improve aroma of thermally processed methionine-glucose Amadori compound | Food Chemistry | [10.1016/j.foodchem.2024.142033](https://doi.org/10.1016/j.foodchem.2024.142033) | Met-glucose ARP + targeted precursors: pyrazines (sibling of deng2022 dossier) | `deng2024.pdf` |
| 100 | 11 (2+18) | Yu 2019 | Effective Mechanism of (−)-Epigallocatechin Gallate Indicating the Critical Formation Conditions of Amadori Compound during an Aqueous Maillard Reaction | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.9b00034](https://doi.org/10.1021/acs.jafc.9b00034) | EGCG and Ala-xylose ARP formation conditions (no S) | `yu2019.pdf` |
| 111 | 11 (2+18) | Van Lancker 2010 | Mechanistic Insights into Furan Formation in Maillard Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf102929u](https://doi.org/10.1021/jf102929u) | furan formation mechanisms, labelling (no S) | `vanlancker2010.pdf` |
| 112 | 11 (3+16) | Yaylayan 1999 | Origin of 2,3-Pentanedione and 2,3-Butanedione in d -Glucose/ l -Alanine Maillard Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf9902292](https://doi.org/10.1021/jf9902292) | Yaylayan 1999: 13C origin of 2,3-butanedione/pentanedione in glucose/alanine | `yaylayan1999.pdf` |
| 115 | 11 (1+20) | Boelens 1974 | Organic sulfur compounds from fatty aldehydes, hydrogen sulfide, thiols, and ammonia as flavor constituents | Journal of Agricultural and Food Chemistry | [10.1021/jf60196a036](https://doi.org/10.1021/jf60196a036) | Boelens 1974: S-heterocycles from fatty aldehydes + H2S/thiols/NH3 | `boelens1974.pdf` |
| 122 | 10.5 (2+17) | Chen 2026 | Exogenous xylose trapped threonine to regenerate Amadori compounds and accelerated pH decline impeding pyrazines formation upon thermal treatment of Amadori compounds revealed by sotope labeling | Current Research in Food Science | [10.1016/j.crfs.2026.101433](https://doi.org/10.1016/j.crfs.2026.101433) | isotope labelling: xylose trapping Thr, ARP regeneration (no S) | `chen2026.pdf` |
| 123 | 10.5 (2+17) | Tuong 2003 | Generation of Thiols by Biotransformation of Cysteine−Aldehyde Conjugates with Baker's Yeast | Journal of Agricultural and Food Chemistry | [10.1021/jf026198j](https://doi.org/10.1021/jf026198j) | cysteine-furfural conjugate -> 2-furfurylthiol by yeast (conjugate route to FFT) | `tuong2003.pdf` |
| 127 | 10.5 (1+19) | Yuan 2024 | Mechanism of pyrazines and thioethers formation promoted by high oxygen concentration in the methionine‐glucose Maillard reaction system | Journal of the Science of Food and Agriculture | [10.1002/jsfa.14094](https://doi.org/10.1002/jsfa.14094) | Met-glucose at high O2: pyrazines and thioethers mechanism | `yuan2024.pdf` |
| 137 | 10.5 (2+17) | Silwar 1989 | Gaschromatographisch-massenspektrometrische Untersuchung der im Modellsystem Cystein/Methionin/Furfural unter Röstbedingungen gebildeten Aromastoffe | European Food Research and Technology | [10.1007/bf01028064](https://doi.org/10.1007/bf01028064) | Silwar 1989: cysteine/methionine + furfural under roasting, aroma compounds (GC-MS) | `silwar1989.pdf` |
| 141 | 10 (4+12) | Yang 2017 | Preparation of 1-Amino-1-deoxyfructose Derivatives by Stepwise Increase of Temperature in Aqueous Medium and Their Flavor Formation Compared with Maillard Reaction Products | Food and Bioprocess Technology | [10.1007/s11947-017-2039-4](https://doi.org/10.1007/s11947-017-2039-4) | ARP derivatives by stepwise heating; flavour vs MRP | `yang2017.pdf` |
| 144 | 10 (1+18) | Deng 2022 | Mechanism of Pyrazine Formation Intervened by Oxidized Methionines during Thermal Degradation of the Methionine–Glucose Amadori Compound | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.2c06458](https://doi.org/10.1021/acs.jafc.2c06458) | oxidised methionine in Met-Glc ARP degradation: pyrazines (sibling of deng2022 dossier) | `deng2022b.pdf` |
| 146 | 10 (1+18) | Ma 2025 | Preparation, characterization, and thermal pyrolysis of three Maillard intermediates | Journal of Analytical and Applied Pyrolysis | [10.1016/j.jaap.2025.107274](https://doi.org/10.1016/j.jaap.2025.107274) | pyrolysis of three Maillard intermediates (Amadori/Heyns) | `ma2025.pdf` |
| 149 | 10 (3+14) | Wang 2022 | Controlled Formation of Pyrazines: Inhibition by Ellagic Acid Interaction with N-(1-Deoxy- d -xylulos-1-yl)-glycine and Promotion through Ellagic Acid Oxidation | Journal of Agricultural and Food Chemistry | [10.1021/acs.jafc.1c07391](https://doi.org/10.1021/acs.jafc.1c07391) | ellagic acid vs Gly-xylose ARP: pyrazines (no S) | `wang2022b.pdf` |
| 153 | 10 (2+16) | Fanun 2001 | Sucrose ester microemulsions as microreactors for model Maillard reaction | Colloids and Surfaces A Physicochemical and Engineering Aspects | [10.1016/s0927-7757(01)00786-5](https://doi.org/10.1016/s0927-7757(01)00786-5) | Garti group: sucrose-ester microemulsions as Maillard microreactors (check model; sibling of round-1 Yaghmur 2002) | `fanun2001.pdf` |
| 164 | 10 (3+14) | Jerić 2009 | Screening for glucose‐triggered modifications of glutathione | Journal of Peptide Science | [10.1002/psc.1159](https://doi.org/10.1002/psc.1159) | glucose-triggered modifications of glutathione (GSH glycation products) | `jeric2009.pdf` |
| 166 | 10 (2+16) | Yaylayan 2000 | Origin of Carbohydrate Degradation Products in l -Alanine/ d -[ 13 C]Glucose Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf000004n](https://doi.org/10.1021/jf000004n) | Yaylayan 2000: 13C glucose/alanine degradation products | `yaylayan2000.pdf` |
| 167 | 10 (2+16) | Yaylayan 2001 | Carbohydrate and Amino Acid Degradation Pathways in l -Methionine/ d -[13C] Glucose Model Systems | Journal of Agricultural and Food Chemistry | [10.1021/jf000986w](https://doi.org/10.1021/jf000986w) | Yaylayan 2001: methionine/[13C]glucose degradation pathways | `yaylayan2001.pdf` |
| 196 | 9 (1+16) | Xia 2022 | Structural diversity and concentration dependence of pyrazine formation: Exogenous amino substrates and reaction parameters during thermal processing of l-alanyl-l-glutamine Amadori compound | Food Chemistry | [10.1016/j.foodchem.2022.133144](https://doi.org/10.1016/j.foodchem.2022.133144) | Ala-Gln ARP pyrazine formation vs exogenous amino substrates (no S) | `xia2022b.pdf` |
| 201 | 9 (0+18) | Pfeifer 2010 | Investigation of Reactive α-Dicarbonyl Compounds Generated from the Maillard Reactions of l -Methionine with Reducing Sugars via Their Stable Quinoxaline Derivatives | Journal of Agricultural and Food Chemistry | [10.1021/jf1008988](https://doi.org/10.1021/jf1008988) | methionine + sugars: alpha-dicarbonyls as quinoxalines | `pfeifer2010.pdf` |
| 202 | 9 (0+18) | Wang 2024 | Rapid preparation of the Amadori rearrangement product of glutamic acid - xylose through intermittent microwave heating and its browning formation potential in microwave thermal processing | Food Research International | [10.1016/j.foodres.2024.114075](https://doi.org/10.1016/j.foodres.2024.114075) | Glu-xylose ARP microwave preparation | `wang2024.pdf` |
| 206 | 9 (3+12) | Lu 2026 | Formation of characteristic aroma and inorganic fugacious products during pyrolysis of Maillard intermediates with different reaction degrees | Journal of Analytical and Applied Pyrolysis | [10.1016/j.jaap.2026.107660](https://doi.org/10.1016/j.jaap.2026.107660) | pyrolysis of Maillard intermediates of different reaction degree | `lu2026.pdf` |

Reading the round-2 lists:

- **Class A is dominated by two old families that one hop missed.** (1) The **Ho-lab (Rutgers) series of 1988-1991** on
  thermal reactions of cysteine and glutathione: degradation in water (Zhang 1988), with 2,4-decadienal (1989), with glucose
  (1991) and with IMP (1991), plus Chen 2000/2002 (urea, carnosine in ribose-cysteine). These give volatile yields per
  precursor under stated T/pH/time, and several run cysteine and GSH side by side under one protocol. (2) **Mottram/Farmer and
  Hofmann/Schieberle siblings** of papers we already use: Farmer & Mottram 1992, Mottram 1995 nonaqueous, Hofmann &
  Schieberle 1998 (roasting), Schieberle & Hofmann 1996 (ribose-cysteine), Guth et al. 1995 (thiol/disulfide stability:
  a candidate sink step for FFT/MFT).
- The **thiamine** side gains Dreher 2003 (slow degradation at 35 C) and, in B, Buttery 1984 (MFT dimer as the thiamine odour).
- Class B adds labelling work from the Yaylayan, Wang & Ho and CAMOLA groups, mostly without sulfur; useful for the sugar-
  fragmentation trunk rather than the S lane.
- Reproduce: `r2.py`, `triage2.py`, `write2.py` and the cache are in the session scratchpad (`lit2/`, cache shared with `lit/`).
