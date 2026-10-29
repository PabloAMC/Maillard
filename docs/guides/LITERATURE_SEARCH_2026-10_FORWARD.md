# Literature search 2026-10: forward citation chaining from the sulfur-lane papers

Run 2026-10-28. This search answers Pablo's question: "What about articles that cite the ones we already have?" It
collects every paper that cites one of the corpus's quantitative sulfur-lane papers, ranks the ones the corpus does not
hold, and reads the top of that list for measured kinetics. The targets are the open path gaps in
[SULFUR_FLUX_MAP.md](SULFUR_FLUX_MAP.md) sections 4-5b. It records the search only. Nothing was added to `data/`, and no
paper was downloaded. It follows [LITERATURE_SEARCH_2026-10.md](LITERATURE_SEARCH_2026-10.md) (rounds 1-2),
[CITATION_CONTEXT_2026-10.md](CITATION_CONTEXT_2026-10.md),
[LITERATURE_SEARCH_2026-10_ROUND3.md](LITERATURE_SEARCH_2026-10_ROUND3.md),
[LITERATURE_SEARCH_2026-10_VARIABLES.md](LITERATURE_SEARCH_2026-10_VARIABLES.md) and
[LITERATURE_SEARCH_2026-10_PATHS.md](LITERATURE_SEARCH_2026-10_PATHS.md). Papers in those documents or in the corpus were
excluded before ranking.

## Bottom line

- **Forward chaining finds almost no modern kinetics for the sulfur lane that the older rounds missed.** 2 432 new citing
  works were ranked and about 120 read by title and abstract. **No citer reports a rate constant, an activation energy,
  or a multi-temperature or multi-pH time course for a cysteine (or glutathione) + pentose system in water at 90-145 °C.**
  The gaps it was meant to fill are (a) the Strecker sulfur split, (b) the TTCA / Cys-Amadori trunk timescale,
  (d) norfuraneol vs pH, (f) the cysteamine / 2-acetylthiazole route, (g) glutathione as its own thiol source and
  (h) MFT dimerisation. For these six, the forward pool holds identification, sensory and application papers only.
- **Why the pool is thin.** The modern Chinese-lab papers on cysteine + xylose intermediates (Cao 2017, Zhai 2019-2023,
  Hou 2017, Tang 2019, Zhao 2019, Zhou 2022-2023) are cited mostly by other papers from the same groups. About a third of
  those citers (217 of 706) are already in the corpus or in an earlier search document. Most of the rest apply the
  intermediates to protein hydrolysates, meat flavourings or extrudates without measuring rates. The 1990s Hofmann,
  Mottram and Cerny papers are cited for odorant identities, by aroma-extract and sensomics studies.
- **What the search does add** is quantitative, but it comes from outside the food literature or from a neighbouring
  system:
  - **Gap (e), furfural formation and decay.** Danon, Marcotullio & de Jong 2014 (Green Chem., open) review the published
    kinetics of pentose → furfural and of furfural loss in water. Danon, van der Aa & de Jong 2013 and Kim et al. 2011
    measure furfural decay and uncatalysed xylose dehydration, and Ershova et al. 2015 (open) adds the xylulose
    intermediate. All four work at pH < 3 or at 150-220 °C. They can bound the formation-vs-decay barrier structure
    (section 5's 64 vs 143 kJ/mol suspect), not set it at pH 5.
  - **Gaps (a)/(b), cysteine loss vs temperature.** Yu & Tang 2016 (Czech J. Food Sci., open): ascorbic acid + cysteine,
    125-155 °C, first-order reactant loss and an apparent Ea of 114 kJ/mol for browning. It is the only modern
    multi-temperature measurement involving cysteine loss that the forward pool contains. The partner is ascorbic acid,
    not a pentose.
  - **Gap (e), path level.** Cerny et al. 2021 (Food Chem.) adds a second labelling study to Poisson 2009: in coffee,
    furfuryl alcohol is an FFT precursor and furfural plays no major role.
- **Calibration.** About 58 % of the ranked citers have an abstract. I put **about 15 %** on the pool holding an aqueous
  90-145 °C rate law for gaps a, b, d, f, g or h that I missed because its abstract was absent or uninformative. Forward
  chaining from 83 seeds is fairly exhaustive for papers that cite us. The remaining risk is papers that cite none of our
  seeds: Chinese-language journals, theses, and the wine-chemistry literature on thiol oxidation.

## Method

- **Seeds: 83 papers, 76 with at least one citer.**
  - **Wave 1 (54 seeds).** The papers the brief named: both Hofmann & Schieberle 1998 JAFC model studies
    (10.1021/jf970990g, 10.1021/jf9705983), Hofmann 1998c, Hofmann 1995/1995b/1996/1997/2001, Cao 2017,
    Zhai 2019/2020/2021/2023 (three), Zhou 2023 and 2023c, Zhang 2025b, Kang 2026, Yaghmur 2002/2005,
    Cerny 2003/2004/2007/2008/2015, Kumazawa 2003, Gigl 2021, Du 2023, Whitfield 1999/2001,
    Mottram 1995/1995b/1996/2002, de Roos 1992/2005, Zheng & Ho 1994, Zheng 1997, Lee 2010, Shu 1985, Blank 2002,
    Guth 1995, Elmore 2002, Liu 2021 and Yiltirak 2026. It also includes the 2026-10-28 downloads: Larralde 1995,
    Garrett 1969, Edwards 1992, Wedzicha 1991, Kobayasi 1965, Mikami 2017, Poisson 2009 and Hicks 1975.
  - **Wave 2 (29 seeds).** Further quantitative sulfur-lane dossiers in `data/lit/extraction_dossiers/`, found by grepping
    for cysteine, thiol, glutathione, thiamine, TTCA and disulfide: Kreitman 2016 (both parts), Hou 2017, Tang 2019,
    Zhao 2017/2019, Yu 2010, Wang 2012, Chen 2000, Hartman 1984, Zhang 1988/1989, Meynier 1995, Madruga 1995,
    Farmer 1989, Mottram 1995c, Cerny 2013, Zhou 2022, Zhou 2025b, Zhang 2026, Wang 2026/2026b, Xu 2011/2013, Sun 2019b,
    Meng 2017, Zhang 2024b, Zheng 2022 and Mizutani 1960.
  - **Resolution.** DOIs came from `data/keys/papers.yml`, or from the dossier when the registry row lacks one. Seven
    seeds have no citer in OpenCitations: Kang 2026, Yiltirak 2026, Zhou 2025b, Zhang 2026, Wang 2026, Wang 2026b and
    Mizutani 1960. The 2026 papers are too new, and Mizutani has no DOI-linked citers.
- **Budget problems.** OpenAlex answered one request, then returned HTTP 429 with "Insufficient budget … $0 remaining;
  resets at midnight UTC" for this IP. The `cites:` filter was therefore not usable this round.
  - **Citers** came from OpenCitations (COCI) `index/v2/citations/doi:…`, anonymous, which returned 4 154 citation edges.
  - **Coverage check.** For 25 core seeds, Semantic Scholar's `/citations` (anonymous, with back-off) was compared against
    OpenCitations. Counts agree within about 10 % for most seeds. Semantic Scholar holds more citers for Zhang 2025b (29 vs
    12), Zhou 2023 (25 vs 18) and Cerny 2015 (16 vs 7). Its extra citers were all listed by title; only the Cannon 2019
    thesis was relevant.
  - **Metadata.** Titles, years and open-access status came from the Semantic Scholar `/paper/batch` endpoint, which needed
    up to 4 retries per call after 429s. Abstracts came from Semantic Scholar and, where it withholds them (mostly
    Elsevier), from the Europe PMC REST search by DOI. Author and venue details for the final list came from Crossref.
  - **Privacy and cache.** No key, e-mail or other identifier was sent. Responses are cached under `scratch/forward_cites/`
    (gitignored).
- **Pool.** The citation edges point to 2 684 distinct citing works.
  - **252 were excluded**: every DOI in `data/keys/papers.yml` and every DOI written in the five 2026-10 search documents.
    Of these, 112 are corpus papers and the rest appear in earlier documents (79 in rounds 1-2, 65 in CITATION_CONTEXT,
    39 in round 3).
  - **2 432 new citers remain**: 1 756 from 2010 or later, and 1 408 with an abstract.
- **Ranking.** Each citer scored 3 × (number of distinct seeds it cites) + 2 × (kinetics-regex hits, capped at 3) +
  1.5 × (system-regex hits, capped at 4) + 3 if the year is ≥ 2010.
  - The kinetics regex: kinetic, rate constant, activation energy, Arrhenius, time course, first/second order, half-life,
    T or pH dependence, modelling, quantitative.
  - The system regex: cysteine, glutathione, thiol, furanthiol, furfurylthiol, H2S, pentose, ribose, xylose, furfural,
    norfuraneol, Amadori, thiazolidine, deoxypentosone, disulfide, thiamine, cysteamine, thiazole.
  - 173 citers cite three or more seeds. All of them are reviews, book chapters or odorant-identification papers. The
    single most co-citing paper, a 2005 thiol review, cites 17 seeds.
- **Triage.** About 200 titles were scanned and about 120 abstracts read: the top of the ranking, every citer of the modern
  seeds with a relevant title, and every citer of Larralde, Garrett, Hicks and Mikami with a kinetics hit. The
  catalysis, solvent and biomass-pretreatment literature was excluded unless it worked in plain water.
  - **Classes.** **Q** = the abstract reports measured rates, time courses, or multi-T or multi-pH levels in an aqueous
    system. **M** = mechanism or labelling only. **B** = identification, sensory, application or method.
  - **Tags.** **[abstract]** marks what I read. **[inferred]** marks what I infer from the title or venue.
- **Counts over the triaged set.** Q 10, of which 9 are listed below. **None is Q for a cysteine/GSH + pentose aqueous pot
  at 90-145 °C.** Six are pentose or furfural chemistry under biomass or prebiotic conditions. One (Yu & Tang 2016) is
  cysteine + ascorbic acid. One (Deng 2021) is a methionine Amadori compound. One (Haleva-Toledo 1999) is from before
  2010. The tenth is Poisson et al. 2019, a coffee-bean time course (see gap e). M: about 8. B: the rest.

---

## Per-gap findings

### (a) Cysteine's Strecker sulfur → H2S vs inert

Nothing quantitative. The best citers of Kobayasi 1965, Zheng & Ho 1994 and Shu 1985 are:

| paper | DOI / access | what it has | class |
|---|---|---|---|
| Yu & Tang 2016, *Czech J. Food Sci.* 34:503 | [10.17221/86/2016-cjfs](https://doi.org/10.17221/86/2016-cjfs), **open** | [abstract] ascorbic acid + cysteine, 125-155 °C, 10-120 min. First-order loss of the reactants is measured; browning follows first-order kinetics with Ea 114.3 kJ/mol. Below 110 °C there is little or no reaction, and cysteine alone does not brown | Q, but for cysteine + an ascorbic-acid-derived dicarbonyl, not a pentose. A **cysteine-loss rate at four temperatures**, against which the model's shared 64 kJ/mol can be tested in a neighbouring system |
| Li, Schieberle & Steinhaus 2016, *Eur. Food Res. Technol.* 243:69 | [10.1007/s00217-016-2723-z](https://doi.org/10.1007/s00217-016-2723-z), paywalled | [abstract] SIDA for H2S (Na2³⁴S), methanethiol, ethanethiol and propanethiol. Applied to durian only | B (method). It is the assay the "measure the H2S split" experiment would use |
| Cannon 2019, Rutgers PhD thesis (Ho lab), "Quantitation of key odorants in arginine/cysteine-glucose Maillard reactions" | [10.7282/t3-nbkr-a666](https://doi.org/10.7282/t3-nbkr-a666), probably open (RUcore) [inferred] | [inferred] no abstract retrieved. Cites Kobayasi 1965 and de Roos 2005. Probably SIDA quantitation of thiols and H2S-derived odorants in cysteine + glucose; a time or pH series is possible | ? (confidence 0.3 that it holds a usable time course) |

### (b) The de Roos path (TTCA ⇌ Cys-Amadori → deoxyosones) and the trunk timescale

Nothing new and quantitative.

- The one hit, a Research Square preprint on TTCA/ARP storage stability vs temperature, pH and a_w
  (10.21203/rs.3.rs-7540317/v1), is the **preprint of Kang 2026**, which the corpus already holds.
- The Cao 2017 / Zhai / Hou / Tang citers outside the corpus are applications: chicken-carcass hydrolysate, mealworm,
  flaxseed, rice protein, MRI + oxidised fat (Food Res. Int. 2025, 116927; temperature screen only). Or they are
  Amadori-compound preparation papers (vacuum dehydration, Fru-Asp, theanine-glucose).
- **Closest analogue: Deng et al. 2021, JAFC 69:5167** ([10.1021/acs.jafc.1c02023](https://doi.org/10.1021/acs.jafc.1c02023),
  paywalled). [abstract] Thermal degradation of the methionine-glucose ARP, with and without added Met, at several
  temperatures and in buffer. Measured ARP loss rate, Met, and C6/GO/MGO dicarbonyls. Added Met does not change the ARP
  degradation rate. Q for a sulfur amino acid's ARP. It is a template rate for a Cys-ARP → deoxyosone step (Cys-ARP itself
  is not measured).
- Zhu et al. 2020, *Food Chem.* 326:126982
  ([10.1016/j.foodchem.2020.126982](https://doi.org/10.1016/j.foodchem.2020.126982)). [abstract] GSH competes with
  asparagine for glucose and forms the GSH-glucose ARP; the paper gives "quantitative evidence" from precursor and
  intermediate changes. M/Q-lite, probably one condition. Relevant to (g) more than (b).

### (c) Pentose 1,2- vs 2,3-enolisation, the pH law and phosphate catalysis

No partition data. Larralde 1995's 313 citers are almost all origin-of-life papers. Three of them measure pentose loss:

| paper | DOI / access | what it has | class |
|---|---|---|---|
| Nitta, Furukawa & Kakegawa 2015, *Orig. Life Evol. Biosph.* 46:189 | [10.1007/s11084-015-9472-z](https://doi.org/10.1007/s11084-015-9472-z), **open** | [abstract] Degradation rates of ribose and the other aldopentoses with and without silicate or phosphate, at high pH and moderate temperature. ³¹P NMR finds **no ribose-phosphate complex** | Q, at alkaline pH. An independent check on the phosphate-catalysis claim; transfer to pH 3-7 is weak |
| Furukawa, Horiuchi & Kakegawa 2013, *Orig. Life Evol. Biosph.* 43:353 | [10.1007/s11084-013-9350-5](https://doi.org/10.1007/s11084-013-9350-5), paywalled | [abstract] Relative stability of the four aldopentoses with and without borate. **Ribose is the least stable** | Q (relative rates). A ribose : xylose reactivity ratio to test against the Hofmann ribose vs xylose arms |
| Bains, Xiao & Yu 2015, *Life* 5:1054 | [10.3390/life5021054](https://doi.org/10.3390/life5021054), **open** | [abstract] First-order hydrolysis and decomposition constants and their T-dependence, compiled for 63 metabolites from the literature. Ribose is probably among them [inferred] | Q (secondary compilation). A cheap cross-check of the pentose-loss Ea |
| Danon, Marcotullio & de Jong 2014 (see gap e) | | Reviews the 1,2-enolisation, β-elimination and isomerisation routes from acyclic xylose, with kinetic data | M + Q (acid) |

### (d) Norfuraneol formation and decay vs pH

Nothing on norfuraneol itself. The citers of Mikami 2017 and Hicks 1975 are not kinetic.

- Haleva-Toledo, Naim, Zehavi & Rouseff 1999, JAFC 47:4140
  ([10.1021/jf9813788](https://doi.org/10.1021/jf9813788), paywalled, pre-2010). [abstract title] Effects of cysteine and
  N-acetylcysteine on the formation of furaneol, HMF and 5-methylfurfural and on browning, in buffered rhamnose or
  glucose + arginine. [inferred] It measures the cysteine sink on a 3(2H)-furanone at fixed T. That is the furaneol
  analogue of the norfuraneol + cysteine competition. M/Q, confidence 0.4 for usable time courses.
- A 2021 review of the stability of 4-hydroxy-2,5-dimethyl-3(2H)-furanone (ACS Food Sci. Technol.,
  10.1021/acsfoodscitech.1c00011) has no abstract available. It may tabulate stability vs pH (confidence 0.3).

### (e) Furfural formation vs temperature and the furfural → FFT branch

This is where forward chaining (from Garrett 1969) adds the most numbers. All of them come from biomass conditions.

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Danon, Marcotullio & de Jong 2014**, *Green Chem.* 16:39 (239 citations) | [10.1039/c3gc41351a](https://doi.org/10.1039/c3gc41351a), **open** | [abstract] Reviews the mechanism and **the published kinetic data for furfural formation from pentoses and for furfural disappearance in aqueous acid**, put on a common footing | Q (compendium). The quickest way to see whether furfural decay has a higher Ea than formation in the literature at all |
| **Danon, van der Aa & de Jong 2013**, *Carbohydr. Res.* 375:145 | [10.1016/j.carres.2013.04.030](https://doi.org/10.1016/j.carres.2013.04.030), paywalled | [abstract] Furfural degradation in dilute acid and saline solution: first/second-order fits and Arrhenius parameters. **Furfural decays faster with glucose present**, through an extra second-order furfural + sugar term (a possible Diels-Alder) | Q. A measured **furfural-decay rate law with a sugar-dependent term**. That bears directly on splitting `r_fur_decay` from the carbonyl-sink family |
| Kim, Lee, Park & Lee 2011, *React. Kinet. Mech. Catal.* 103:267 | [10.1007/s11144-011-0320-5](https://doi.org/10.1007/s11144-011-0320-5), paywalled | [inferred] no abstract. Xylose dehydration in high-temperature water without added acid, probably 180-220 °C, with rate constants and an Ea | Q, extrapolation needed |
| Ershova, Kanervo, Hellsten et al. 2015, *RSC Adv.* 5:66727 | [10.1039/c5ra10855a](https://doi.org/10.1039/c5ra10855a), **open** | [abstract] Uncatalysed and acid-catalysed xylose → furfural at 180-220 °C. Separate xylulose and furfural conversion runs, and a mechanistic kinetic model. Maximum furfural yield 48 % uncatalysed | Q, extrapolation needed |
| Cerny, Schlichtherle-Cerny, Gibe et al. 2021, *Food Chem.* 337:128008 | [10.1016/j.foodchem.2020.128008](https://doi.org/10.1016/j.foodchem.2020.128008), paywalled | [abstract] In-bean spiking: furfuryl alcohol → FFT dose-dependently. d2-FFA gives d2-FFT and **no labelled furfural**. Furfural "plays no major role" as an FFT precursor in roasting; the proposed route is FFA → furfuryl cation → S-furfuryl-cysteine → FFT | M, path level. With Poisson 2009 this makes **two labelling studies against free furfural + H2S as the main FFT route in a real matrix** |
| Poisson, Pittet, Schaerer et al. 2019, JAFC 68:4732 | [10.1021/acs.jafc.9b05436](https://doi.org/10.1021/acs.jafc.9b05436), paywalled | [abstract] Formation kinetics of key coffee odorants over 10 roasting time points, recombinate vs native beans, by SIDA | Q in a bean. Probably includes FFT [inferred]. Not an aqueous pot |

**Reading.** The biomass numbers sit at pH < 2 and/or 180-220 °C. They can rank the formation and decay barriers
against each other, which is what the section 5 suspect needs, but they cannot set either at pH 5 and 100-145 °C.
Cerny 2021 moves the prior a little further toward "the model's FFT route is wrong, not only its Ea". Its transfer from
bean to buffer is weak, as Poisson 2009's was. The experiment proposed in the PATHS document (Hofmann's furfural + Na2S
protocol at 100, 120 and 145 °C) remains the decisive move.

### (f) Cysteamine / 2-acetylthiazole route

Nothing. The citers of Hofmann 1995 are odorant identifications. One citer is a 2008 synthesis of diacetylthiazoles, and
another is a 1998 symposium paper on thiazolidines from cysteamine + carbonyls (10.1016/s0167-4501(98)80073-3,
mechanism only).

### (g) Glutathione as its own thiol source

Nothing at food temperatures.

- Lee 2010 and Zheng & Ho 1994 are cited by sensory papers on GSH Maillard products (beef stock, kokumi, yeast extract).
- Zhu 2020 (above) is the only one that follows a GSH-sugar Amadori compound.
- Dekker et al. 2022 (*Food Chem.* 371:131341, [10.1016/j.foodchem.2021.131341](https://doi.org/10.1016/j.foodchem.2021.131341))
  follows cysteinylated and glutathionylated polysulfides during fermentation. That is a time course, but at cellar
  temperature and yeast-driven. B for this gap.

### (h) MFT → bis(2-methyl-3-furyl) disulfide in water

Nothing.

- The Kreitman 2016 citers (Cu- and Fe-catalysed oxidation of H2S and thiols) are wine, catalysis and sensing papers.
- None measures furanthiol dimerisation at 100-140 °C.
- The stability papers in the pool (Sun 2020 on FFT loss to quinones in coffee brew; encapsulated flavourings) are
  qualitative.

---

## Ranked download list (new papers only; at most 12)

The confidence column gives my probability that the paper holds the number in a form the engine can use. File names
follow round 3.

| # | paper | DOI | year | gap | access | the number it would give | confidence |
|---|---|---|---|---|---|---|---|
| 1 | Danon, van der Aa & de Jong, *Carbohydr. Res.* 375:145 | [10.1016/j.carres.2013.04.030](https://doi.org/10.1016/j.carres.2013.04.030) | 2013 | e | paywalled | Arrhenius parameters for furfural decay, including a furfural + sugar second-order term: a measured decay barrier to set against formation (`danon2013.pdf`) | 0.7 (acid only) |
| 2 | Danon, Marcotullio & de Jong, *Green Chem.* 16:39 | [10.1039/c3gc41351a](https://doi.org/10.1039/c3gc41351a) | 2014 | e, c | **open** | Tabulated literature Ea and pre-factors for pentose → furfural and furfural loss, on one footing: tests whether decay Ea > formation Ea (`danon2014.pdf`) | 0.8 for a table; 0.2 for anything above pH 3 |
| 3 | Yu & Tang, *Czech J. Food Sci.* 34:503 | [10.17221/86/2016-cjfs](https://doi.org/10.17221/86/2016-cjfs) | 2016 | a, b | **open** | Cysteine and ascorbic acid loss rate constants at 125-155 °C and a browning Ea of 114 kJ/mol: a neighbouring-system Ea for cysteine destruction (`yu2016cjfs.pdf`) | 0.6 |
| 4 | Ershova, Kanervo, Hellsten et al., *RSC Adv.* 5:66727 | [10.1039/c5ra10855a](https://doi.org/10.1039/c5ra10855a) | 2015 | e | **open** | Uncatalysed xylose → furfural and furfural loss rate constants in water at 180-220 °C (`ershova2015.pdf`) | 0.7 (needs extrapolation) |
| 5 | Kim, Lee, Park & Lee, *React. Kinet. Mech. Catal.* 103:267 | [10.1007/s11144-011-0320-5](https://doi.org/10.1007/s11144-011-0320-5) | 2011 | e, c | paywalled | Ea for xylose dehydration in hot water without added acid (`kim2011.pdf`) | 0.6 [inferred] |
| 6 | Cerny, Schlichtherle-Cerny, Gibe et al., *Food Chem.* 337:128008 | [10.1016/j.foodchem.2020.128008](https://doi.org/10.1016/j.foodchem.2020.128008) | 2021 | e (path) | paywalled | Labelled FFA → FFT yield and the absence of furfural involvement: a path-level check on the 95-99 % furfural + H2S FFT branch (`cerny2021.pdf`) | 0.9 route; 0.05 kinetics |
| 7 | Deng, Cui, Hayat et al., JAFC 69:5167 | [10.1021/acs.jafc.1c02023](https://doi.org/10.1021/acs.jafc.1c02023) | 2021 | b (analogue) | paywalled | Degradation rate vs T of an S-amino-acid Amadori compound (Met-Glc ARP) and its dicarbonyl yields: a template for the Cys-ARP step (`deng2021.pdf`) | 0.5 |
| 8 | Nitta, Furukawa & Kakegawa, *Orig. Life Evol. Biosph.* 46:189 | [10.1007/s11084-015-9472-z](https://doi.org/10.1007/s11084-015-9472-z) | 2015 | c | **open** | Pentose degradation rate constants ± phosphate (alkaline): whether phosphate accelerates pentose loss at all (`nitta2015.pdf`) | 0.6 (wrong pH) |
| 9 | Haleva-Toledo, Naim, Zehavi & Rouseff, JAFC 47:4140 | [10.1021/jf9813788](https://doi.org/10.1021/jf9813788) | 1999 | d (analogue) | paywalled | Furaneol formation and loss with vs without cysteine in buffer: the size of a cysteine sink on a furanone (`halevatoledo1999.pdf`) | 0.4 |
| 10 | Cannon, Rutgers PhD thesis (Ho lab) | [10.7282/t3-nbkr-a666](https://doi.org/10.7282/t3-nbkr-a666) | 2019 | a, b | probably open [inferred] | SIDA levels of thiols and H2S-derived odorants in cysteine + glucose, possibly vs time or pH (`cannon2019.pdf`) | 0.3 |
| 11 | Furukawa, Horiuchi & Kakegawa, *Orig. Life Evol. Biosph.* 43:353 | [10.1007/s11084-013-9350-5](https://doi.org/10.1007/s11084-013-9350-5) | 2013 | c | paywalled | Relative degradation rates of ribose, arabinose, xylose and lyxose: a ribose : xylose reactivity ratio for the Hofmann arms (`furukawa2013.pdf`) | 0.5 |
| 12 | Bains, Xiao & Yu, *Life* 5:1054 | [10.3390/life5021054](https://doi.org/10.3390/life5021054) | 2015 | c | **open** | Compiled first-order decomposition constants vs T: a cross-check of the pentose-loss Ea (`bains2015.pdf`) | 0.3 |

**The most decisive item per gap:**

- (e): Danon 2013 + Danon 2014, then the three-temperature furfural + Na2S experiment from the PATHS document.
- (a)/(b): Yu & Tang 2016 is a neighbouring system only. Cao 2017, Zhang 2025b and Kang 2026, already in the corpus, remain
  the best timescale data.
- (c): no new partition data. Nitta 2015 tests the phosphate claim at the wrong pH.
- (d), (f), (g), (h): nothing found by forward chaining. These need the experiment or a different search route:
  Chinese-language journals (CNKI), theses, and the wine-chemistry thiol-oxidation literature searched by keyword
  instead of by citation.

## Verdict

Forward chaining from 83 corpus seeds reached 2 684 citing works. It confirms that the earlier rounds had already harvested
most of the modern quantitative sulfur-lane work: a third of the modern seeds' citers are already held, and the rest are
applications. It finds **no modern rate data for cysteine/GSH + pentose chemistry in water**. What it adds is acid and
hot-water furfural kinetics from the biomass literature, one cysteine + ascorbic acid Ea, and a second labelling study
against the furfural → FFT route.
