# Census sweep for the seven EMPTY gaps (2026-11-01, read-only)

*Search record (what was searched, how, and which gaps stay empty). The papers to download live in the owner's private queue, not here.*

Method. Three channels.
- **OpenAlex census.** 102 structured queries, all `title_and_abstract.search`, `has_abstract:true`, years 1960-2026, run in three rounds (boolean, then narrow phrases). Round 3 was cut short when the free OpenAlex daily budget ran out (the last ~20 G1/G3/G6/G7 phrasings). Their content was covered with Crossref `query.bibliographic` instead (about 35 queries).
- **Patents.** Google Patents and USPTO pages, reached through web search and two fetched patent pages.
- **Chinese literature.** Web search on spkx.net.cn, spgykj.com, sciengine and CNKI, plus OpenAlex hits from Chinese journals.

Known-paper check. Each hit was compared with 1,507 DOIs gathered from docs/guides/LITERATURE_SEARCH_*.md, scratchpad/heldout_search.md, crossdomain/h2s searches, DOWNLOAD_QUEUE.md, data/keys/papers.yml, all of data/ and docs/reference/. Dossier stems in data/articles were checked too.

Raw files are in `scratchpad/census/` (the `oa/` and `cr/` JSON, plus `G*.tsv`, which flags each hit KNOWN or new).

DOI rule. "Crossref OK" means the DOI was resolved in api.crossref.org on 2026-11-01 and its title matches. Chinese-journal DOIs (10.7506, 10.13386) are registered with ChinaDOI, not Crossref. Those were checked by doi.org resolution only and are marked "ChinaDOI".

## Census counts (OpenAlex, unique works)

| gap | unique hits | already held (DOI) | on-topic after screening (approx.) | new plausible | fills the gap? |
|---|---|---|---|---|---|
| G1 pentose ARP split by in-situ OPD | 43 | 21 | ~20 | 3 (all weak) | **no** |
| G2 MFT / FFT heated-water loss | 96 | 36 | ~30 | 5 | **no** |
| G3 thiazolidine k vs pH and T, Ea | 126 | 11 | ~15 (most hits are pharma or biology) | 8 | **partly** (pH only; no Ea) |
| G4 O2 vs inert, 100-150 C | 280 | 19 | ~12 | 3 (weak) | **no** |
| G5 a_w 0.5-0.9 with deoxyosones or odorants | 272 | 9 | ~20 | 7 | **yes, one strong hit** |
| G6 pea / soy bound-amine kinetics | 139 | 10 | ~15 | 6 | **partly** |
| G7 free NF vs ribose MFT by labelling, >= 120 C | 17 | 8 | 9 | 0 | **no** |

Saturation signal. For G1, G2, G4 and G7 almost every on-topic hit was already held. Most were held through the same dozen labs: Jiangnan / Ho, DFA, Nestle, Reading, BTBU.

## Candidates (ranked: gap x confidence x new lab)

| # | citation | DOI | gap | what it likely reports (from abstract) | lab (new?) | conf. |
|---|---|---|---|---|---|---|
| 1 | Cui H, ..., Wang S 2026, LWT 258:120020, chestnut protein isolate + maltose / sucrose | 10.1016/j.lwt.2026.120020 (Crossref OK) | **G5** (+G6) | Full grid: a_w 0.53 / 0.75 / 0.97 x 100 / 150 / 200 C. LC-MS α-dicarbonyls (MGO ...), CEL / CML, GC-MS volatiles. Protein-bound amine. | Nankai University (NEW) | high |
| 2 | Wolf JC, Thompson DR, Reineccius GA 1977, J Food Process Preserv 1:271 | 10.1111/j.1745-4549.1977.tb00330.x (Crossref OK) | **G6** | Kinetics of available-lysine loss in heated soy protein isolate. Title says a no-loss phase reverts to a rapid-loss phase. | Minnesota (Reineccius; existing id `minnesota_reineccius`) | med-high |
| 3 | Bustillo Trueba P, ..., De Cooman L 2018, J Am Soc Brew Chem 76:272 | 10.1080/03610470.2018.1546094 (Crossref OK) | **G3** | Degradation rate of cysteinylated aldehydes (2-substituted TCAs) to free aldehyde at pH 4.4 / 5.2 / 6.0. Depends on the 2-substituent. Temperature series unknown. | KU Leuven brewing / Ghent campus (De Cooman, Aerts) (NEW; not Hendrickx) | med-high |
| 4 | Karpinski ZJ, Karny M, Kublik Z 1990, J Electroanal Chem 288:129 | 10.1016/0022-0728(90)80030-a (Crossref OK) | **G3** | Electrochemical (polarographic) study of TCA decomposition in water. This is the classical way to get ring-opening / hydrolysis rates vs pH. No abstract. | Univ. Warsaw, chemistry (NEW; Radomski 1989 lead is also Warsaw) | medium |
| 5 | Baert JJ, ..., De Cooman L 2015, JASBC 73:100 and 73:243 | 10.1094/asbcj-2015-0109-01 ; 10.1094/asbcj-2015-0531-01 (Crossref OK) | G3 | Extent of aldehyde + cysteine binding vs pH 4.4-6.0 (more binding at higher pH). Thiazolidine of furfural quantified in beer. | KU Leuven brewing (same as #3) | medium |
| 6 | Butvin P, Al-Ja'afreh J, Svetlik J, Havranek E 2000, Chem. Papers (Bratislava) | no Crossref DOI (only ChemInform abstract 10.1002/chin.200017131) | G3 | Solubility, stability and dissociation constants of (2RS,4R)-2-substituted TCAs in water. | Comenius Univ. Bratislava (NEW) | medium |
| 7 | Ames JM, Guy RCE, Kipping GJ 2001, J Agric Food Chem 49:4315 (not the corpus's ames2001 dossier, which is the cysteine / starch paper) | 10.1021/jf010198m (Crossref OK) | G5 | Effect of pH, temperature and moisture on volatiles from extrusion-cooked glycine / glucose. | Reading (DEV) | medium |
| 8 | Cremer DR, Eichner K 2000, Eur Food Res Technol 211:247 and 211:400 | 10.1007/s002170000192 ; 10.1007/s002170000196 (Crossref OK) | G5 | Companions of the cremer2000 dossier on disk. Strecker aldehydes in a low-moisture system vs pH. Strecker aldehydes from ARP + amino acids. | Muenster (Eichner; existing id) | medium |
| 9 | Ayllon-Parra N, del Castillo MD, ..., Ribas-Agusti A 2026/27, IFSET 115:104870 | 10.1016/j.ifset.2026.104870 (Crossref OK) | G6 | Maillard markers in fava-bean protein + oat blends, solid-state fermentation and low-moisture extrusion. | IRTA + CIAL Madrid (NEW) | medium |
| 10 | (authors not retrieved) 2020, Sci. Technol. Food Ind. 41(15), pork-hydrolysate glucose-cysteine model | 10.13386/j.issn1002-0306.2020.15.004 (ChinaDOI; not Crossref) | G5 | 3 x 3 factorial: 115 / 120 / 125 C x a_w 0.75 / 0.80 / 0.85. Browning, e-nose, volatile sulfur classes. Probably GC-MS class totals, not MFT. Chinese, with an English abstract (likely). | unknown Chinese lab (probably NEW) | low-med |
| 11 | Manzocchi E, ..., Noziere P 2023, J Dairy Sci 106:7675 | 10.3168/jds.2022-22902 (Crossref OK) | G6 | Lupin extruded with or without reducing sugars: MRP formation, CML. | INRAE Clermont-Theix ruminant nutrition (NEW; distinct from Berdague) | low-med |
| 12 | Lee Y, Kim Y 2026, npj Sci Food 10 | 10.1038/s41538-026-00804-9 (Crossref OK) | G6 | Hydrolysed vegetable proteins (mung bean / soy / pea) + reducing sugars: volatile formation. | Korea (may merge with `korea_kim_yun_park`) | low-med |
| 13 | Parliment TH, Stahl HD 1994, ACS Symp Ser 564:160; and 1995, Dev Food Sci 37:805 | 10.1021/bk-1994-0564.ch013 ; 10.1016/s0167-4501(06)80197-4 (Crossref OK) | G2 (formation side) | FFT generation in cysteine-pentose and coffee models. Possibly time courses with a late decline. No abstract. | Kraft General Foods (NEW) | low-med |
| 14 | Zweig I, Kamyshny A 2025, Geochim Cosmochim Acta 409:29 (dataset on Mendeley) | 10.1016/j.gca.2025.06.009 (Crossref OK) | G2 (O2 analogue) | Rate law for CH3SH + O2 in water. Protonated and thiolate pathways; orders ~2.2 / 1.6 in thiol and ~1 in O2. 25 C only. | Ben-Gurion Univ. (NEW, cross-domain) | low-med |
| 15 | Gajda T, ... 1990, J Chem Soc Dalton 3155 | 10.1039/dt9900003155 (Crossref OK) | G3 | Protonation constants of sugar-derived 2-(polyhydroxyalkyl)TCAs. Might give ring-chain information. | Szeged (NEW) | low |
| 16 | Zhang PJ, Zhao Z, Yu S, He S 2011, Int J Food Sci Technol 47:475 | 10.1111/j.1365-2621.2011.02866.x (Crossref OK) | G1 (weak) | Thermal degradation of xylose-Gly and xylose-Pro ARPs. Xylose regenerated (reversibility). No dicarbonyl trapping. | SCUT (existing id `scut_guan_yu`) | low |
| 17 | Hirsch J, Mossine VV, Feather MS 1995, Carbohydr Res 273:171 | 10.1016/0008-6215(95)00114-9 (Crossref OK) | G1 (weak) | Detection of dicarbonyl intermediates from Amadori degradation (trapping). Pentose coverage unknown. | Missouri (Feather; existing id) | low |
| 18 | Hellwig M, Henle T 2010, Eur Food Res Technol 230:903 | 10.1007/s00217-010-1237-3 (Crossref OK) | G1 (weak) | 3-deoxypentosone + lysine to formyline. Not an ARP split. | TU Dresden (NEW) | low |
| 19 | Hayase F et al. 1996, Biosci Biotechnol Biochem 60:1820 | 10.1271/bbb.60.1820 (Crossref OK) | G4 (weak) | O2 and transition metals in the advanced Maillard reaction of protein + glucose. Likely 37-50 C. | Meiji Univ. (NEW) | low |
| 20 | Haffenden LJW 2007, PhD thesis, "role of redox reactions in Maillard model systems" | no Crossref DOI | G4 (weak) | Oxidants / reductants and aroma / colour. | McGill (Yaylayan; existing id) | low |
| 21 | Tian H et al. 2004, Food Sci. (China) 25(10):33, oxidative stability of key thiols in reaction flavours | no DOI found | G2 | Ranks MFT as the least stable thiol; proposes degradation routes. Probably storage conditions, not 80-150 C. | BTBU (DEV) | low |
| 22 | Patent WO2005019165A1 (Givaudan; de Roos, Wolswinkel), cysteine-sugar Amadori meat flavour | patent page read (Google Patents) | G3 (weak) | Cys-Xyl Amadori conversion of 25-35% at 80-95 C, pH 5-5.4. Thiazolidine equilibrium not quantified. No rate constants. | Givaudan (DEV) | low |
| 23 | Patent EP1090557B1 / US6716468 (Quest, now Givaudan; Kerler et al.) | patent page read | G2 (weak) | Only one O2 test: MFT oxidised to disulfide at 20 C over 20 h. No heated kinetics. | Quest (DEV) | low |

Already held and re-found (not candidates): BTBU's Food Science papers (zhao2017b, Cao 2017, and the hou2017 pH paper); the Jiangnan pentose-ARP series; Hofmann 1996; Guth 1995; van Seeventer 2001; Blank 2002; the Xu supercritical CO2 vs N2 series; Yuan 2024 (high O2); Hartman 1984; Leahy 1989; Warmbier 1976; Cremer 2000 (Food Chem); Lee 2024 cakes.

## Gaps still EMPTY after this sweep

- **G1 EMPTY.** In-situ OPD trapping on a pentose ARP at pH 4.5-7.5 exists only in Jiangnan / Ho (already held). No second lab found.
- **G2 EMPTY.** No heated-water (80-150 C) lifetime for MFT, and none with O2, cysteine or carbonyl dependence. The nearest are analogues: Zweig 2025 (CH3SH + O2, 25 C) and the held DFA / Quest storage studies.
- **G4 EMPTY.** No O2-vs-inert contrast on Maillard odorants or dicarbonyls at 100-150 C beyond the held Glomb / Nestle / Jiangnan items.
- **G7 EMPTY.** Nothing besides the already queued Tressl 1994.
- **G3 partly.** There is pH-dependent stability from a new lab (KU Leuven brewing; Warsaw electrochemistry). Still no rate constants with Ea vs temperature for sugar thiazolidines.
- **G6 partly.** Wolf 1977 gives soy bound-lysine loss kinetics. Flavour-formation kinetics on pea / soy at cooking temperatures remain absent.
- **G5 filled at one point.** Nankai 2026 (a_w x T grid with α-DCs and volatiles).

## Patents and Chinese channels

- **Patents: EMPTY for G2, G3 and G7.** The flavour-house patents (Givaudan / Quest) give conversions and sensory data. None has a kinetic table or labelled results. No Nestle, Firmenich, Unilever or Ajinomoto patent with such data turned up.
- **Chinese.** The Food Science (China) cysteine-xylose papers come from BTBU (DEV) and are already held. The only new item is #10 (STFI 2020, a_w x T). The ChinaDOI prefixes are not in Crossref, so CNKI items cannot be Crossref-verified.

## Recommendation

Stop broad searching. For G1, G2, G4 and G7, 9 + 1 searches have converged: re-hits are almost entirely held papers. Further sweeps will most likely (~85 %) find nothing new; the remaining chance is a thesis or a non-indexed Chinese paper. Only two narrow follow-ups are worth doing:
1. Download #1, #2, #3 and #4 (about 1 h to triage).
2. One chained citation pass (citers and references) from Bustillo Trueba 2018 and Karpinski 1990, after the OpenAlex budget resets, looking for a temperature series and Ea for TCA ring opening (G3).

G2, G4 and G7 need the designed experiments already drafted (E3 brief), or data from the authors.
