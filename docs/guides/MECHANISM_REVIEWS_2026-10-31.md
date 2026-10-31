# Mechanism reviews, 2026-10-31: what the literature says about the six core problems

*Six independent read-only reviews (one per core problem), each reading the dossiers and PDFs on disk and searching
for what is missing. The raw reports are in `docs/history/mechanism_reviews_2026-10-31/`. Numbers marked "derived"
are the reviewers' own arithmetic on keyed data; nothing here is a fit. This synthesis ranks what changes a decision
in `STATUS_AND_PLAN.md` sec. 5. Probabilities are the reviewers', kept where they were argued.*

## 1. The findings that change the plan

| # | problem | finding | P | consequence |
|---|---|---|---|---|
| 1 | deoxyosone hub | both deoxyosones are short-lived and quasi-steady next to their parent (Hofmann 1999, Zhai 2023 JAFC, Zhai 2021 pH 5.5; Martins M4: 1-DG ~1 min, 3-DG ~30 min); the engine's 3-DX lives ~35 h, so the ratio climbs | 0.65 | the hub step starts from Martins & van Boekel 2005 M4 as its null |
| 2 | hub / cysteine | free cysteine traps the aldehydic 3-DX, faster at high pH -- the pre-declared signature (vii) (+xylose / -xylose double ratio of TDP / DPO, measured 0.85 / 0.92 / 1.22 / 2.71 at pH 5.5 - 8) is reproduced at no pH 7 / 8 point by any arm (all ~0.73) | 0.45 | a cysteine-3-DX adduct with a free-amine factor enters the hub step |
| 3 | cysteine fate | the dilute fed-TTCA pots are kinetically far from equilibrium: adding xylose RAISES free cysteine (Zhai 2021 Fig. 6), which no equilibrium model allows; the same lab's 83 mM pot forms TTCA 8 - 30x faster per mol -- ring hydrolysis / recombination catalysed by solutes and buffer | 0.55 | B7b becomes catalysed rate laws k = k0 + k_aa [AA] + k_P [phosphate], dilute and concentrated pots fitted together, the xylose contrast must pass |
| 4 | thiols | the Cerny / Schieberle conflict most likely comes from the MFT LOSS law (the engine destroys 93 % of the MFT it makes at 100 C, so only an accelerating precursor -- norfuraneol -- can give Schieberle's rise); with an hours-long MFT lifetime (CF-013) a steady intact-C5 supply already gives the lag | 0.40 | a cheap dry check first: k_mft_decay x0.1 - 0.2 at the Cerny-passing corner |
| 5 | thiols | Liu 2023's "19 % fragment recombination" at 168 C is an EI-MS artefact: its MFT channels match the NIST spectrum's [M-H]+ / M+ (0.30 vs 0.28), as do furfural's and 2-methylfuran's | 0.8 | Liu's signature (thiol prereg sec. 11) is withdrawn; the "route mix changes with temperature" finding loses its high-T anchor |
| 6 | pyrazines | the lane uses a pyrazine-YIELD constant as the Strecker FLUX: fed-dicarbonyl pots turn ~1 % of the aminoketone into pyrazine (the rest to pyrroles, Adams 2008); Hofmann 2000 puts the true Strecker flux ~1.9 dex above B18's constant | 0.55 | pyrazines follow straight after the forced-input Strecker layer (roadmap step 2), with an aminoketone -> pyrrole sink |
| 7 | Strecker | a fed-donor yield is a branching share, not a rate (MGO is gone in ~15 min); the rate law is mass action on TOTAL amino acid with a shallow pH factor (<= 0.07 dex / pH); O2 opens a SEPARATE acid route (PA -> PAA < 1 %); lipid enals need O2 for the aldehyde itself | 0.5 | adopted in the forced-input prereg; a slope-1 free-amine law is refuted by > 3 dex |
| 8 | matrix | protein-bound lysine is the dominant amine in a real alt-protein product (~60 mmol/kg at 16 % pea isolate vs ~1 mmol/kg free amino acids; ~1.9x glycine's reactivity per amine): ~90 % of the sugar-amine flux | 0.40 (largest gap) | the composable network's amine classes (roadmap step 6); product numbers need designed data (P 0.85) |
| 9 | matrix | **code bug:** the B20 glycation arm returns lysine on fructoselysine decay; Nguyen 2016's model does not | -- | CF-072; fixed in its own PR |
| 10 | process | sec. 13 of the B7 prereg wrongly said Zhai 2021's pots were phosphate-buffered (they are water): phosphate catalysis is a live rival for the Davidek / Zhai tension | -- | corrected (dpo prereg, CORRECTION note) |

## 2. The cheapest discriminating checks (no new data)

1. The MFT-loss dry check (finding 4): ~30 min, one core.
2. Signature (vii) and the +-xylose free-cysteine contrast as must-pass tests of the hub step (findings 2, 3): declared.
3. The methylpyrazine ratio in Leahy Table I: methylpyrazine^2 / (pyrazine x 2,5-DMP) = 11.7 / 14.4 / 17.3 at 95 / 85 /
   75 C against 4 for one shared condensation constant -- methylpyrazine has an extra route (finding 6).
4. 2-mercapto-3-pentanone / MFT from ribose as a presence / absence signature across temperature (absent at 95 C,
   present at 145 C).
5. The Weenen 2001 vs Hofmann 2000 Strecker yield from MGO at pH ~3 (21 % vs ~1 %): a cross-laboratory test an
   aldehyde-loss term must explain.

## 3. Papers to download (deduplicated across the six reviews; DOIs verified against Crossref or PubMed)

"Abstract" = the reviewer read the abstract; "title only" = only the bibliographic record was verified, so the
content is NOT described. Already on disk and therefore left out: Kallen 1971 (dossier exists), Zamora 2015's
supplementary data, Koutsidis 2008a / b, Zhou 2024 (the pyrazine paper, 10.1021/acs.jafc.4c03706 -- the alanine / ARP
paper below is a different one).

**Tier 1 -- each decides a pending design choice**

| paper | DOI | verified | why |
|---|---|---|---|
| Smuda, Voigt & Glomb 2010, J. Agric. Food Chem. 58:6458 | 10.1021/jf100334r | abstract | amine-INDUCED beta-dicarbonyl cleavage of 1-DG (acetyl / formyl lysine amides): the amine-acceleration question the Davidek papers leave open |
| Wei, Cui, Hayat, Zhang & Ho 2022, J. Agric. Food Chem. 70:12164 | 10.1021/acs.jafc.2c04961 | abstract | cysteine traps the deoxyosones of a Gly-Xyl ARP as TTCA: the 3-DX trap behind signature (vii) |
| Zyzak, Richardson, Thorpe & Baynes 1995, Arch. Biochem. Biophys. 316:547 | 10.1006/abbi.1995.1073 | abstract | the Amadori compound decomposes faster in phosphate than HEPES, faster in air and at high pH (37 C): phosphate catalysis |
| Fife, Natarajan, Shen & Bembi 1991, J. Am. Chem. Soc. 113:3071 | 10.1021/ja00008a041 | title only | ring opening AND hydrolysis of thiazolidines as separate steps -- the two steps the engine lumps |
| Granvogl, Bugan & Schieberle 2006, J. Agric. Food Chem. 54:1730 | 10.1021/jf0525939 | abstract | Ile / Leu / Met / Phe + MGO (Val only in cocoa; corrected from the dossier): Strecker aldehydes vs amines, for the donor-generic layer |
| Low, Parker & Mottram 2007, J. Agric. Food Chem. 55:4087 | 10.1021/jf070044s | abstract | [2-13C]glycine into alkylpyrazines, two pathways with their shares |
| Hayashi & Namiki 1980, Agric. Biol. Chem. 44:2575 | 10.1271/bbb1961.44.2575 | abstract | the early C2 (Namiki) fragment from sugar + amine: the pre-Amadori nitrogen / dicarbonyl route |
| Warmbier, Schnickels & Labuza 1976, J. Food Sci. 41:528 | 10.1111/j.1365-2621.1976.tb00663.x | abstract | casein-glucose-glycerol: first-order available-lysine loss vs water activity (browning maximum a_w 0.4 - 0.5) |

**Tier 2 -- useful constants or rivals**

| paper | DOI | verified | why |
|---|---|---|---|
| Pesek & Frost 1975, Tetrahedron 31:907 | 10.1016/0040-4020(75)80099-8 | title only | thiazolidine decomposition in acid and base (a rate / pH profile, if any) |
| Schubert 1936, J. Biol. Chem. 114:341 | 10.1016/S0021-9258(20)65226-7 | title only | thiol-acid + aldehyde compounds; Schubert 1939 cites it as the companion with constants |
| Zhou, Huang, ..., Ho 2024, J. Agric. Food Chem. 72:5878 | 10.1021/acs.jafc.4c00021 | abstract | exogenous alanine speeds the ARP's reaction with 3-deoxyxylosone and lowers furfural: a 3-DX sink |
| Cui, Ma, ..., Ho 2022, J. Agric. Food Chem. 70:8409 | 10.1021/acs.jafc.2c03427 | abstract | glycylglycine catalyses its Amadori compound's degradation to deoxyosone (an amine factor on the split) |
| Voigt, Smuda, Pfahler & Glomb 2010, J. Agric. Food Chem. 58:5685 | 10.1021/jf100140h | abstract | O2-dependent lactic / glyceric acid from 1-DG + alanine |
| Hou, Zhou, Li, ..., Ho 2026, Food Chem. 499:147363 | 10.1016/j.foodchem.2025.147363 | abstract | cystine speeds ARP degradation and the retro-aldol of deoxyosones to GO / MGO |
| Thornalley, Langborg & Minhas 1999, Biochem. J. 344:109 | 10.1042/bj3440109 | abstract | GO / MGO / 3-DG from glucose and fructosyl-lysine depend on phosphate and trace metals (37 C) |
| Guan, Wang, Yu, Yu & Zhao 2012, J. Dairy Sci. 95:590 | 10.3168/jds.2011-4860 | abstract | phosphate speeds the reversible breakdown of early glucose + proline intermediates (bifunctional proton transfer) |
| Pripis-Nicolau et al. 2000, J. Agric. Food Chem. 48:3761 | 10.1021/jf991024w | abstract | six dicarbonyls x wine amino acids at pH 3.5: donor-generic Strecker |
| Hidalgo, Gallardo & Zamora 2005, J. Agric. Food Chem. 53:10254 | 10.1021/jf052240+ | abstract | 4-hydroxynonenal + Phe -> phenylacetaldehyde (17 %): the lipid-carbonyl Strecker donor |
| Zamora, Navarro, Gallardo & Hidalgo 2006, J. Agric. Food Chem. 54:6101 | 10.1021/jf061239n | abstract | epoxyalkenals -> alpha-keto acids (the acid route) |
| Scalone, Cucu, De Kimpe & De Meulenaer 2015, J. Agric. Food Chem. 63:5364 | 10.1021/acs.jafc.5b01129 | abstract | oligopeptides vs free amino acids in pyrazine formation, a_w 0.33, > 120 C |
| Van Lancker, Adams & De Kimpe 2010, J. Agric. Food Chem. 58:2470 | 10.1021/jf903898t | abstract | lysine dipeptides vs free amino acids in pyrazines (the 2012 follow-up is on disk) |
| Lu, Hao, Payne & Ho 2005, J. Agric. Food Chem. 53:6443 | 10.1021/jf050534p | abstract | glucose + glycine / di- / triglycine at 160 C, 0 - 100 % water in glycerol |
| Adams, Mottram, Parker & Brown 2001, J. Agric. Food Chem. 49:4333 | 10.1021/jf0100797 | abstract | volatile disulfides lost to ovalbumin, more after denaturation (registry cites it; no PDF) |
| Lo, Westwood, McLellan, Selwood & Thornalley 1994, J. Biol. Chem. 269:32299 | 10.1016/S0021-9258(18)31635-1 | search snippet only | MGO with lysine, cysteine and arginine (37 C): the dicarbonyl traps of protein amines |
| Shibamoto & Bernhard 1977, Agric. Biol. Chem. 41:143 | 10.1271/bbb1961.41.143 | abstract | pyrazine pathways in glucose-ammonia (O2, metals, ring closure by water loss) |

**Tier 3 -- title only, older, lower expected value:** van den Ouweland & Peer 1975 (10.1021/jf60199a045,
norfuraneol + H2S volatiles; in the registry, no PDF); Huber & Ledl 1990 (10.1016/0008-6215(90)84037-U, amine-bound
1,4-dideoxyhexodiuloses); Koehler & Odell 1970 (10.1021/jf60171a041); Shibamoto & Bernhard 1976 (10.1021/jf60206a045);
Eichner & Karel 1972 (10.1021/jf60180a025); Piloty & Baltes 1979 (10.1007/BF01186494).

## 4. Gaps no search filled (designed data or authors)

No time-resolved MFT isotope labelling at 95 - 100 C; no fed-DDP experiment; no H2S time course in a sugar + cysteine
pot; no DPO decay in water vs phosphate; no glycation kinetics on pea or soy protein; no thiol capture by a plant
protein at process temperature; no measured aminoketone-dimerisation or dihydropyrazine-oxidation constant; no
Strecker step measured in water (G1). Zhai 2021's 13C5-xylose label fraction in TTCA would give the formation constant
at 10 mM directly (only MS fragments are printed -- ask the authors).
