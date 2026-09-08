# Elmore & Mottram 1997 — EXTRACTION (alkanal + alpha-hydroxyketone or alpha-dione + ammonium sulfide in water, 140 C / 30 min; the 2-alkyl-3-thiazoline, 2-alkylthiazole and trialkylpyridine slate as percent of chromatogram area)
### The mechanism paper for the lipid–Maillard thiazolines: a fed-intermediate pot (no cysteine, no ribose) that makes every alkylthiazoline and alkylthiazole tentatively reported in pressure-cooked beef, and, with a 2,4-dienal, 2-pentylthiophene at 21 % of the chromatogram.

**Source on disk:** `data/articles/elmore1997.pdf` (owner's download, 2026-09-08). Born-digital ACS
PDF; the `pypdf` text layer is clean but drops inter-word spaces in the running text. Tables 1 and 2
were read from the text layer and checked against a rendering of p. 3597. Figure 1 (the formation
scheme) is a drawn mechanism, rendered and described in words below; Figure 2 is a mass-spectral
fragmentation scheme (not used). Tables 3-6 are LRI + m/z lists for the products (identification
data); not re-typed, only the LRIs of the compounds the repository names are quoted.

## 0. Identity

| field | value |
|---|---|
| Title | "Investigation of the Reaction between Ammonium Sulfide, Aldehydes, and alpha-Hydroxyketones or alpha-Dicarbonyls To Form Some Lipid-Maillard Interaction Products Found in Cooked Beef" |
| Authors | J. Stephen Elmore, Donald S. Mottram (Food Science and Technology, University of Reading) |
| Venue | J. Agric. Food Chem. 1997, 45 (9), 3595-3602 |
| DOI | 10.1021/jf970065u |
| Companion | Elmore, Mottram, Enser, Wood, J. Agric. Food Chem. 1997, 45, 3603-3607: the cooked-beef paper whose tentative thiazoline / thiazole identifications this paper confirms ("c" superscripts in Tables 3-4 mark the compounds found in beef) |
| Naming | "3-thiazoline" = 2,5-dihydrothiazole, ring S1-C2-N3=C4-C5; the 2-substituent is the aldehyde's chain minus its carbonyl carbon; 4- and 5-substituents come from the hydroxyketone or dione. "Alkane-alpha-dione" = 2,3-butanedione or 2,3-pentanedione. |

## 1. Why it matters

The roadmap's wishlist row "2-pentyl- and 2-hexyl-4-methylthiazole, alkylthiophenes,
2-pentylpyridine" (`tasks/roadmap_for_scientists.md`, the Programme 7 table) has no rule and no
number. This paper supplies the drawn mechanism for the thiazolines and thiazoles and shows the two
registry compounds directly: 4-methyl-2-pentylthiazole (`2_pentyl_4_methylthiazole`) at 0.7 % of
the chromatogram from hexanal + 1-hydroxy-2-propanone + (NH4)2S, and 2-hexyl-4-methylthiazole
(`2_hexyl_4_methylthiazole`) at 0.6 % from heptanal + 1-hydroxy-2-propanone + (NH4)2S (Table 1),
with the 3-thiazolines 18.3 % and 14.3 % beside them. It also gives the one clean dienal + H2S
observation in the group: (E,E)-2,4-nonadienal + acetoin + (NH4)2S makes 2-pentylthiophene at 21 %
of the chromatogram, which is the thiophene half of the Farmer 1990 / Mottram 2002 dienal scheme
(`farmer1990_extraction.md` §4, `mottram2002b_extraction.md` §4) run on a pure aldehyde. The pot is a
fed-intermediate system (H2S and NH3 supplied as ammonium sulfide, the carbonyl supplied), so it
fixes the ring-forming chemistry, not the rate at which cysteine + ribose supplies the partners.

## 2. Methods as they matter to a model

- **Charge:** 50 mL of 0.1 M ammonium sulfide in deionised water (made from 20 % aqueous (NH4)2S;
  5 mmol (NH4)2S = 10 mmol NH4+ + 5 mmol sulfide) + 0.005 mol aldehyde + 0.005 mol hydroxyketone or
  dione. Nominal concentrations 100 mmol/L aldehyde, 100 mmol/L carbonyl, 100 mmol/L (NH4)2S (200
  mmol/L N, 100 mmol/L S). No cysteine, no sugar, no buffer, no added salt.
- **pH:** not stated. An aqueous ammonium sulfide solution is alkaline; the pot is therefore
  not at the pH 5-6 of the sulfur lane's reference systems (flag 2).
- **Vessel and heating:** 100 mL autoclavable screw-top jar with PTFE-lined cap (about 50 mL
  headspace, air not excluded); autoclave, 140 C, 30 min. One pot per combination; no replicates
  stated.
- **Aldehydes used:** 2-methylpropanal, 3-methylbutanal, pentanal, hexanal, heptanal, octanal,
  nonanal, decanal, (E)-2-nonenal, (E,E)-2,4-nonadienal. Carbonyls: 3-hydroxy-2-butanone (acetoin),
  1-hydroxy-2-propanone (acetol), 1-hydroxy-2-butanone, 2,3-butanedione, 2,3-pentanedione. Not every
  aldehyde was run with every carbonyl; the combinations are the rows of Tables 1-2.
- **Work-up:** whole mixture extracted with 30 mL dichloromethane; extract stored at 4 C; 1 µL
  injected split 40:1 at 250 C; BPX5 50 m x 0.32 mm x 0.5 µm; 60 C (2 min) -> 280 C at 10 C/min
  -> 280 C (10 min); EI 70 eV, m/z 28-450; LRI against C6-C22 alkanes.
- **Quantification: percent of total chromatogram peak area (TIC).** No internal standard, no
  response factors, no calibration. "Total peak areas for all the reaction mixtures were similar"
  is the only cross-pot anchor. Unreacted aldehyde and carbonyl were each below 10 % of the total
  area except nonanal + acetol (35 %) and nonanal + 2,3-butanedione (33 %).
- **Identification:** LRI and full spectra (Tables 3-6), NIST where available; each 5-substituted
  3-thiazoline gives an E/Z pair "in approximately equal quantities", reported as the combined area.

## 3. Tables re-typed

### Table 1. "3-Thiazolines, Thiazoles, and Pyridines Isolated from the Reactions of alpha-Hydroxyketones with Ammonium Sulfide and Aliphatic Aldehydes." Percent of total chromatogram peak area.

Footnote a: for each 5-substituted 3-thiazoline two isomers (E and Z) were found in approximately
equal quantities; combined percentages are given. Blank = not listed.

**Reactions using 3-hydroxy-2-butanone (acetoin)**

| aldehyde | 3-thiazoline | % | thiazole | % | pyridine | % |
|---|---|---:|---|---:|---|---:|
| 2-methylpropanal | 4,5-dimethyl-2-isopropyl | 26.0 | 4,5-dimethyl-2-isopropyl | 0.6 | | |
| 3-methylbutanal | 4,5-dimethyl-2-isobutyl | 15.2 | 4,5-dimethyl-2-isobutyl | 0.2 | 3,5-diisopropyl-2-isobutyl | 0.8 |
| pentanal | 2-butyl-4,5-dimethyl | 18.2 | 2-butyl-4,5-dimethyl | 0.2 | 2-butyl-3,5-dipropyl | 2.4 |
| hexanal | 4,5-dimethyl-2-pentyl | 36.7 | 4,5-dimethyl-2-pentyl | 1.7 | 3,5-dibutyl-2-pentyl | 2.2 |
| heptanal | 4,5-dimethyl-2-hexyl | 41.8 | 4,5-dimethyl-2-hexyl | 0.3 | 3,5-dipentyl-2-hexyl | 9.7 |
| octanal | 4,5-dimethyl-2-heptyl | 33.5 | 4,5-dimethyl-2-heptyl | 0.5 | 3,5-dihexyl-2-heptyl | 22.6 |
| nonanal | 4,5-dimethyl-2-octyl | 31.4 | 4,5-dimethyl-2-octyl | 0.4 | | |
| decanal | 4,5-dimethyl-2-nonyl | 26.2 | 4,5-dimethyl-2-nonyl | 0.7 | | |
| (E)-2-nonenal | 4,5-dimethyl-2-hexyl | 0.9 | | | | |
| | 4,5-dimethyl-2-(1-octenyl) | 0.2 | | | | |
| | 4,5-dimethyl-2-octyl | 1.4 | | | | |
| (E,E)-2,4-nonadienal | 2-butyl-4,5-dimethyl | 1.7 | | | | |
| | 4,5-dimethyl-2-(1,3-octadienyl) | 0.7 | | | | |

**Reactions using 1-hydroxy-2-propanone (acetol)**

| aldehyde | 3-thiazoline | % | thiazole | % | pyridine | % |
|---|---|---:|---|---:|---|---:|
| 3-methylbutanal | 2-isobutyl-5-methyl | 1.8 | | | 3,5-diisopropyl-2-isobutyl | 2.8 |
| | 2-isobutyl-4-methyl | 6.6 | | | | |
| pentanal | 2-butyl-5-methyl | 2.3 | | | 2-butyl-3,5-dipropyl | 2.6 |
| | 2-butyl-4-methyl | 5.7 | 2-butyl-4-methyl | <0.1 | | |
| hexanal | 5-methyl-2-pentyl | 4.4 | | | 3,5-dibutyl-2-pentyl | 4.7 |
| | 4-methyl-2-pentyl | 18.3 | **4-methyl-2-pentyl** | **0.7** | | |
| heptanal | 2-hexyl-5-methyl | 3.8 | 2-hexyl-5-methyl | <0.1 | 3,5-dipentyl-2-hexyl | 39.6 |
| | 2-hexyl-4-methyl | 14.3 | **2-hexyl-4-methyl** | **0.6** | | |
| octanal | 2-heptyl-5-methyl | 1.5 | | | 3,5-dihexyl-2-heptyl | 43.5 |
| | 2-heptyl-4-methyl | 5.5 | 2-heptyl-4-methyl | 0.5 | | |
| nonanal | 5-methyl-2-octyl | 1.6 | | | | |
| | 4-methyl-2-octyl | 6.9 | 4-methyl-2-octyl | <0.1 | | |

**Reactions using 1-hydroxy-2-butanone**

| aldehyde | 3-thiazoline | % | thiazole | % | pyridine | % |
|---|---|---:|---|---:|---|---:|
| hexanal | 5-ethyl-2-pentyl | 15.2 | | | 3,5-dibutyl-2-pentyl | 11.3 |
| | 4-ethyl-2-pentyl | 16.0 | 4-ethyl-2-pentyl | 1.2 | | |
| heptanal | 5-ethyl-2-hexyl | 8.0 | 5-ethyl-2-hexyl | <0.1 | 3,5-dipentyl-2-hexyl | 19.5 |
| | 4-ethyl-2-hexyl | 8.3 | 4-ethyl-2-hexyl | 1.0 | | |
| octanal | 5-ethyl-2-heptyl | 4.8 | | | 3,5-dihexyl-2-heptyl | 48.1 |
| | 4-ethyl-2-heptyl | 7.6 | 4-ethyl-2-heptyl | 0.6 | | |
| nonanal | 5-ethyl-2-octyl | 10.5 | | | | |
| | 4-ethyl-2-octyl | 21.2 | 4-ethyl-2-octyl | 0.8 | | |

### Table 2. "Thiazolines, Thiazoles, and Pyridines Isolated from the Reactions of Alkane-alpha-diones with Ammonium Sulfide and Aliphatic Aldehydes." Percent of total chromatogram peak area (E/Z pairs combined).

**Reactions using 2,3-butanedione**

| aldehyde | 3-thiazoline | % | thiazole | % | pyridine | % |
|---|---|---:|---|---:|---|---:|
| pentanal | 2-butyl-4,5-dimethyl | 13.4 | 2-butyl-4,5-dimethyl | 12.8 | 2-butyl-3,5-dipropyl | 3.5 |
| hexanal | 4,5-dimethyl-2-pentyl | 15.7 | 4,5-dimethyl-2-pentyl | 16.8 | 3,5-dibutyl-2-pentyl | 8.1 |
| heptanal | 4,5-dimethyl-2-hexyl | 8.0 | 4,5-dimethyl-2-hexyl | 6.9 | 3,5-dipentyl-2-hexyl | 29.4 |
| octanal | 4,5-dimethyl-2-heptyl | 4.9 | 4,5-dimethyl-2-heptyl | 3.3 | 3,5-dihexyl-2-heptyl | 27.6 |
| nonanal | 4,5-dimethyl-2-octyl | 9.4 | 4,5-dimethyl-2-octyl | 2.6 | | |

**Reactions using 2,3-pentanedione**

| aldehyde | 3-thiazoline | % | thiazole | % | pyridine | % |
|---|---|---:|---|---:|---|---:|
| 3-methylbutanal | 4-ethyl-2-isobutyl-5-methyl | 6.0 | 4-ethyl-2-isobutyl-5-methyl | 3.1 | 3,5-diisopropyl-2-isobutyl | 0.7 |
| | 5-ethyl-2-isobutyl-4-methyl | 7.8 | 5-ethyl-2-isobutyl-4-methyl | 3.0 | | |
| pentanal | 2-butyl-4-ethyl-5-methyl | 10.4 | 2-butyl-4-ethyl-5-methyl | 6.0 | 2-butyl-3,5-dipropyl | 2.9 |
| | 2-butyl-5-ethyl-4-methyl | 10.1 | 2-butyl-5-ethyl-4-methyl | 4.4 | | |
| hexanal | 4-ethyl-5-methyl-2-pentyl | 14.8 | 4-ethyl-5-methyl-2-pentyl | 7.1 | 3,5-dibutyl-2-pentyl | 6.4 |
| | 5-ethyl-4-methyl-2-pentyl | 10.4 | 5-ethyl-4-methyl-2-pentyl | 5.2 | | |
| heptanal | 4-ethyl-2-hexyl-5-methyl | 3.8 | 4-ethyl-2-hexyl-5-methyl | 4.0 | 3,5-dipentyl-2-hexyl | 12.2 |
| | 5-ethyl-2-hexyl-4-methyl | 5.8 | 5-ethyl-2-hexyl-4-methyl | 4.2 | | |
| octanal | 4-ethyl-2-heptyl-5-methyl | 4.6 | 4-ethyl-2-heptyl-5-methyl | 1.3 | 3,5-dihexyl-2-heptyl | 29.1 |
| | 5-ethyl-2-heptyl-4-methyl | 4.9 | 5-ethyl-2-heptyl-4-methyl | 1.0 | | |
| nonanal | 4-ethyl-5-methyl-2-octyl | 10.4 | 4-ethyl-5-methyl-2-octyl | 4.2 | | |
| | 5-ethyl-4-methyl-2-octyl | 12.5 | 5-ethyl-4-methyl-2-octyl | 4.0 | | |
| decanal | 4-ethyl-5-methyl-2-nonyl | 18.7 | 4-ethyl-5-methyl-2-nonyl | 7.3 | | |
| | 5-ethyl-4-methyl-2-nonyl | 18.9 | 5-ethyl-4-methyl-2-nonyl | 5.7 | | |

### Numbers given only in the text ("Other Compounds" and the unsaturated-aldehyde pots), percent of total peak area

| pot | compound | % | note |
|---|---|---:|---|
| every acetoin pot | tetramethylpyrazine | 12 to 30 | from 2 x 3-amino-2-butanone (NH3 on acetoin) |
| every 2,3-butanedione pot | tetramethylpyrazine | 0.3 to 2.0 | |
| 2,3-butanedione + pentanal | 2-butyl-4,5-dimethylimidazole | 28 | imidazoles only in dione pots; falls with chain length |
| 2,3-butanedione + nonanal | 2-octyl-4,5-dimethylimidazole | 1.4 | |
| 2,3-pentanedione + 3-methylbutanal | the 2-isobutyl-ethylmethylimidazole | 27 | |
| 2,3-pentanedione + decanal | the 2-nonyl-ethylmethylimidazole | 10 | |
| acetol + pentanal | four 3-thiazolines of MW 227 (two pentanal units; tentatively 2-butyl-4-methyl-5-pentyl- and 2-butyl-5-methyl-4-pentyl-3-thiazoline, E/Z each) | 12 (sum) | present in every 1-hydroxy-2-alkanone pot, falling with aldehyde chain length |
| nonanal pots; decanal pots | 2-heptyl-2-undecenal; 2-octyl-2-dodecenal (aldol products) | about 20 | aldol share rises with chain length |
| (E)-2-nonenal + acetoin | tetramethylpyrazine | 16 | |
| (E)-2-nonenal + acetoin | 3,4-diheptyl-2(5H)-thiophenone, two isomers (structure by interpretation) | 18 and 11 | 2 x nonenal + 1 x H2S |
| (E,E)-2,4-nonadienal + acetoin | tetramethylpyrazine | 30 | |
| (E,E)-2,4-nonadienal + acetoin | **2-pentylthiophene** | **21** | "a major product"; no thiazoles or trialkylpyridines in either unsaturated pot |
| all pots | pyrazines, pyrroles, imidazoles, oxazoles, alkanethiols | trace to 30 | not itemised |

### LRIs (BPX5) of the compounds the repository names, from Tables 3-4

4-methyl-2-pentyl-3-thiazoline 1377; 4-methyl-2-pentylthiazole 1322; 2-hexyl-4-methyl-3-thiazoline
1483; 2-hexyl-4-methylthiazole 1426; 4,5-dimethyl-2-pentyl-3-thiazoline 1401 / 1412 (E/Z);
4,5-dimethyl-2-pentylthiazole 1384; 3,5-dibutyl-2-pentylpyridine 1918.

## 4. Routes and numbers the repository can use

| route or quantity | reactant -> product | mechanism as drawn (Figure 1) | measured numbers, units, conditions | evidence class |
|---|---|---|---|---|
| E97-A hydroxyketone thiolation | R2-C(=O)-CH(OH)-R3 + H2S -> R2-C(=O)-CH(SH)-R3 + H2O | the -OH is replaced by -SH ("substitution of the -OH group with an -SH"); no intermediate drawn | none (the mercaptoketone is not quantified; Farmer 1990 Table 2 rows 44-47 measure the same mercaptoketones in cysteine + ribose pots) | mechanism_drawn |
| E97-B aldimine | R1-CHO + NH3 -> R1-CH=NH + H2O | drawn as the imine of the alkanal | none | mechanism_drawn |
| E97-C 3-thiazoline closure | mercaptoketone + aldimine -> 2-R1-4-R2-5-R3-3-thiazoline + H2O | S of the mercaptoketone attacks the imine carbon giving R2-C(=O)-CH(R3)-S-CH(NH2)-R1; the NH2 condenses onto the ketone C=O with loss of water; ring S1(from H2S)-C2(aldehyde carbon, bears R1)-N3(from NH3)=C4(ketone carbon, bears R2)-C5(carbinol carbon, bears R3) | acetoin pots: 3-thiazoline 15-42 % of area for C4-C10 alkanals (Table 1); acetol pots: 4-methyl isomer 5.5-18.3 %, 5-methyl isomer 1.5-4.4 % (4-isomer favoured, more so for long chains); 1-hydroxy-2-butanone pots: 4-ethyl 7.6-21.2 %, 5-ethyl 4.8-15.2 % | mechanism_drawn + peak_area_only |
| E97-D thiazoline oxidation | 3-thiazoline -> thiazole | "[O]" arrow; the authors note the ammonium sulfide pot is reducing, which "discourages thiazole formation" | hydroxyketone pots: thiazole / 3-thiazoline area ratio 0.002-0.05 (e.g. hexanal + acetol 0.7 / 18.3 = 0.04); dione pots: about 1 (pentanal + butanedione 12.8 / 13.4; hexanal 16.8 / 15.7; heptanal 6.9 / 8.0) | within_study_ratio (area) |
| E97-E dione route to thiazole | R2-C(=O)-C(=O)-R3 + H2S -> R2-C(=O)-C(OH)(SH)-R3; either [H] -> mercaptoketone (joins E97-C) or + aldimine -> R2-C(=O)-C(R3)(OH)-S-CH(NH2)-R1 -> -H2O -> 5-hydroxy-3-thiazoline -> -H2O -> thiazole | as drawn: H2S adds to one carbonyl; the hemithioketal S attacks the imine; two dehydrations give the thiazole without an oxidation step | see E97-D: dione pots give thiazole and thiazoline in about equal amounts | mechanism_drawn |
| E97-F trialkylpyridine | 3 R-CH2-CHO + NH3 -> 2-(R-CH2)-3,5-R2-pyridine (+ 3 H2O + H2) | not drawn; cited to Shu 1985 and Hwang 1986 (three aldehydes condense with one ammonia); e.g. hexanal -> 3,5-dibutyl-2-pentylpyridine | 0.7-48 % of area, rising with chain length in every series (Tables 1-2); "absent" for nonanal and decanal because the C27-C30 pyridines did not elute (an analytical, not a chemical, zero) | peak_area_only |
| E97-G dienal + H2S -> 2-alkylthiophene | (E,E)-2,4-nonadienal + H2S -> 2-pentylthiophene | not drawn here (Farmer 1990 Fig 1 and Mottram 2002 Fig 3 draw it: 1,4-addition of H2S at C4, S onto the aldehyde carbon, dehydration; see those dossiers) | **21 % of total area** in the nonadienal + acetoin + (NH4)2S pot; the thiazolines from the same pot only 1.7 + 0.7 % | peak_area_only |
| E97-H enal + H2S -> thiophenone | 2 (E)-2-nonenal + H2S -> 3,4-diheptyl-2(5H)-thiophenone (two isomers) | not drawn; "formed from two molecules of (E)-2-nonenal and one molecule of hydrogen sulfide" | 18 % and 11 % of area | peak_area_only (structure by interpretation) |
| E97-I aminoketone -> pyrazine (competition for NH3) | 2 acetoin + 2 NH3 -> tetramethylpyrazine | via 3-amino-2-butanone (NH3 on the hydroxyketone) | 12-30 % of area in acetoin pots; 0.3-2.0 % in butanedione pots | peak_area_only |
| E97-J dione + 2 NH3 + aldehyde -> imidazole | 2,3-butanedione + 2 NH3 + R-CHO -> 2-R-4,5-dimethylimidazole | "both keto groups react with ammonia"; hydroxyketones do not give imidazoles because -OH is not replaced by -NH2 | 28 % (pentanal) falling to 1.4 % (nonanal) | peak_area_only |
| E97-K oxidative fragmentation of unsaturated aldehydes | (E)-2-nonenal -> heptanal fragment (2-hexyl-thiazoline 0.9 %); (E,E)-2,4-nonadienal -> pentanal fragment (2-butyl-thiazoline 1.7 %) | stated: "oxidative fragmentation of the double bonds occurred during the reaction" (air in the jar) | 0.9 % and 1.7 % of area | peak_area_only |

Chain-length bookkeeping (needed by any rule writer): an alkanal CnH2n+1-CHO gives a 2-(CnH2n+1)
substituent, i.e. the 2-alkyl has one carbon fewer than the aldehyde (hexanal -> 2-pentyl; heptanal
-> 2-hexyl). The trialkylpyridine from that alkanal is 2-(CnH2n+1)-3,5-bis(Cn-1H2n-1)-pyridine
(hexanal -> 2-pentyl-3,5-dibutyl).

## 5. Rule sketches (reactant -> product in words; controls are mine, SMILES mine)

**S1. alpha-hydroxyketone + H2S -> alpha-mercaptoketone (E97-A; the hydroxyketone analogue of R27,
which is dicarbonyl + H2S).** Reactant: R-C(=O)-CH(OH)-R'. Change: the carbinol O leaves as water,
S enters on that carbon.
- positive: acetoin `CC(O)C(C)=O` + `S` -> 3-mercapto-2-butanone `CC(S)C(C)=O` (the compound Farmer 1990 lists as "2-mercapto-3-butanone", row 44); acetol `CC(=O)CO` + `S` -> `CC(=O)CS` (species key `MP`).
- negative: 1-hexanol `CCCCCCO` + `S` -> no fire (no alpha-carbonyl; Farmer 1990's alkanol -> alkanethiol is a different, unproven route); hexanal `CCCCCC=O` + `S` -> no fire.

**S2. alkanal + NH3 -> aldimine (E97-B).** Positive: hexanal `CCCCCC=O` + `N` -> `CCCCCC=N`.
Negative: 2,3-butanedione (`DA`) + `N` -> not this rule (a ketone; the dione + NH3 route is E97-J).

**S3. alpha-mercaptoketone + aldimine -> 2-alkyl-3-thiazoline + H2O (E97-C; net of S-attack and
condensation).** Reactant pattern: HS-CH(R3)-C(=O)-R2 and R1-CH=NH.
- positive: `CC(S)C(C)=O` + `CCCCCC=N` -> 4,5-dimethyl-2-pentyl-3-thiazoline `CCCCCC1SC(C)C(C)=N1` (Table 1: 36.7 % with acetoin + hexanal); `CC(=O)CS` + `CCCCCC=N` -> 4-methyl-2-pentyl-3-thiazoline `CCCCCC1SCC(C)=N1` (18.3 %).
- negative: MFT (a thiol without a beta-ketone) + `CCCCCC=N` -> no fire; `CC(S)C(C)=O` + hexanal (no imine, no NH3) -> no fire.
- If the repository prefers one net three-component rule: R-CHO + NH3 + R2-C(=O)-CH(OH)-R3 + H2S -> 3-thiazoline + 3 H2O; controls as above with `CC(O)C(C)=O`, `CCCCCC=O`, `N`, `S`.

**S4. 3-thiazoline -> thiazole (E97-D, oxidation; terminal).** Positive: `CCCCCC1SC(C)C(C)=N1` ->
`CCCCCc1nc(C)c(C)s1`; `CCCCCC1SCC(C)=N1` -> 4-methyl-2-pentylthiazole `CCCCCc1nc(C)cs1`
(registry `2_pentyl_4_methylthiazole`); the heptanal analogue -> `CCCCCCc1nc(C)cs1`
(`2_hexyl_4_methylthiazole`). Negative: the thiazolidine `TTCA` -> no fire (saturated ring, no
C=N); 2-acetylthiazole `ACTZ` -> no fire (already aromatic).

**S5. alpha-dione + H2S + aldimine -> 2-alkylthiazole directly (E97-E; net, no oxidant).**
Positive: `CC(=O)C(C)=O` (`DA`) + `S` + `CCCCCC=N` -> `CCCCCc1nc(C)c(C)s1`. Negative: acetoin
+ `S` + `CCCCCC=N` -> this rule must not fire (hydroxyketones give the thiazoline, S3, and only
traces of thiazole, Table 1).

**S6. three alkanals + NH3 -> 2,3,5-trialkylpyridine (E97-F; terminal).** Positive: 3 x
`CCCCCC=O` + `N` -> 3,5-dibutyl-2-pentylpyridine `CCCCCc1ncc(CCCC)cc1CCCC`. Negative: 3 x
acetone -> no fire (ketone); 3 x benzaldehyde `O=Cc1ccccc1` -> no fire (no alpha-CH2). This rule
competes with S2/S3 for NH3 and with the aldol; Table 1 shows its share rising with chain length.

**S7. (E,E)-2,4-nonadienal + H2S -> 2-pentylthiophene (E97-G; see farmer1990 S2 for the drawn
mechanism).** Positive: `CCCC/C=C/C=C/C=O` + `S` -> `CCCCCc1cccs1`. Negative: hexanal + `S` -> no
fire; (E)-2-nonenal `CCCCCC/C=C/C=O` + `S` -> no thiophene (this paper: a 2(5H)-thiophenone from
two enals instead).

## 6. Flags

1. **Percent of total chromatogram area, no internal standard, one pot per combination.** Nothing
   here is a yield or a concentration. The thiazole / thiazoline ratio within one pot is the only
   number that survives as a ratio; cross-pot comparisons rest on "total peak areas were similar".
2. **pH not stated, no buffer.** The medium is 0.1 M ammonium sulfide in water, i.e. alkaline
   (my inference from the reagent; the paper gives no value). The sulfur lane's reference pots are
   pH 5-6. Thiazoline closure and the aminoketone -> pyrazine competition are both pH-sensitive;
   do not port the area shares to pH 5.
3. **Fed intermediates.** H2S and NH3 are supplied at 100 and 200 mmol/L; cysteine + ribose pots
   generate them at unknown, much lower, transient levels. The paper fixes what forms, not how fast.
4. **Air in a 100 mL jar with 50 mL liquid.** The authors invoke oxidative fragmentation of the
   unsaturated aldehydes (E97-K) and oxidation of thiazolines (E97-D); both depend on the headspace
   oxygen, which is not controlled.
5. **The 4- versus 5-methyl regiochemistry** (acetol pots) is decided by where the SH sits after
   E97-A; the 4-isomer dominates, more so for long aldehydes. A rule written on acetol
   (`CC(=O)CO`) gives the 4-methyl series; the 5-methyl series needs the acetol <-> 2-hydroxypropanal
   isomer, which the paper does not draw.
6. **MW 227 thiazolines** with two pentanal units (12 % of area in acetol + pentanal) are
   uncharacterised and larger than the "expected" thiazolines in that pot; the 1-hydroxy-2-alkanone
   systems are messier than Table 1 shows.
7. **Trialkylpyridines "absent" from nonanal and decanal pots is an elution artefact** (column
   held at 280 C; the C27-C30 pyridines did not elute); not evidence against the route.
8. **Registry match:** 4-methyl-2-pentylthiazole (0.7 %) and 2-hexyl-4-methylthiazole (0.6 %) are
   the two roadmap compounds; in this pot their 3-thiazolines are 26 and 24 times larger by area.
   Any rule for the thiazoles should carry the thiazoline as the primary product and the thiazole as
   its oxidation (S4), or the dione route (S5) if 2,3-butanedione is the partner.
9. **No cysteine-derived acetaldehyde run:** the Strecker aldehyde of cysteine (acetaldehyde) would
   give 2-methyl-thiazolines; the smallest aldehyde here is 2-methylpropanal. The paper's
   introduction cites Takken 1976 for acetaldehyde + pentanedione + H2S + NH3.
10. **No mass-unit number in the paper.** For a yield in mass units the group's dienal work
    (`mottram2002b_extraction.md` Table I) and the Whitfield 1988 ng table are the nearest.
