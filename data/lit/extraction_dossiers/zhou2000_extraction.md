# Zhou & Boatright 2000 — EXTRACTION (5-13C-2,4-decadienal and 15N-ammonia spiked into defatted soy flour slurry during isolate processing at room temperature; 2-pentylpyridine by isotope-dilution GC-MS; pH 4.5 / 7 / 9; eight amino acids)
### The isotope paper for 2-pentylpyridine: fixes the atom map (dienal C5 -> pyridine C2, ammonia N -> ring N), shows the reaction runs at room temperature in water at pH 9, and gives ppm levels in a real soy protein isolate.

**Source on disk:** `data/articles/zhou2000.pdf` (owner's download, 2026-09-08). Read from the `pypdf`
text layer, which is clean; Tables 1-3 came through row by row with their letters and were checked for
column count. Figure 1 (SPI process scheme), Figure 2 (synthesis of the labelled dienal), Figures 3-4
(mass spectra) are described from the text. **Figure 5 (pH dependence) is FIGURE-ONLY**: the text gives
the ordering pH 9 > 7 > 4.5 but no numbers; the layout re-extraction of page 5 found only rotated axis
text and no readable values.

## 0. Identity

| field | value |
|---|---|
| Title | "Precursors for Formation of 2-Pentyl Pyridine in Processing of Soybean Protein Isolates" |
| Authors | A. Zhou, W. L. Boatright (Animal Science Dept., University of Kentucky) |
| Venue | J. Food Sci. 65(7), 2000, 1155-1159 (MS 19991220) |
| Naming | 2-pp = 2-pentyl pyridine; SPI = soybean protein isolate; Stage I / II / III = sampling points of Figure 1 (after the pH 9 extraction; after the pH 4.5 precipitation; after lyophilisation, in that order as far as the text allows; Figure 1 not readable) |
| Companions | Boatright & Crum 1997 (2-pp has the largest flavour value of any SPI volatile); Zhou & Boatright 1999, J. Food Sci. 64, 852 (pH effect during SPI processing; 2-pp does not rise without the alkaline extraction); Boatright et al. 1998 (pro-oxidants FeCl3, CuCl2, UV raise 2-pp); Buttery 1977 (the Schiff-base / dihydropyridine mechanism they test) |

## 1. Why it matters

Programme 7 wants 2-pentylpyridine from the isolate's own lipid and nitrogen. This is the paper that
proves, in the soy matrix itself, that the carbon skeleton is 2,4-decadienal (5-13C label -> 2-13C in the
pyridine, M+ 149 vs 148) and the nitrogen is ammonia (15N -> 1-15N, same shift), that the reaction
needs no heat (it runs during a 1 h room-temperature alkaline extraction), that free ammonia rather than
ammonium is the reactive form (pH 9 > 7 > 4.5 in both buffer and slurry), and that in pH 7 buffer at
room temperature no amino acid (Arg, Lys, Asn, Gln included) gives 2-pp while in the flour slurry Arg,
Lys, Asn and Gln do, which the authors read as the flour supplying a deamidation route (possibly
enzymatic). The baseline 0.14-0.21 ppm of 2-pp in untreated isolates and the 2.0-2.6 ppm after a 5.64 mM
dienal spike are the only in-matrix numbers in the corpus. The dienal-limited / ammonia-rich picture of
the flour ("the amount of 2,4-decadienal was very limited, but ammonia existed") is the qualitative
stoichiometry the rule should carry.

## 2. Methods as they matter to a model

- **Flours:** Stressland, Edison, KS4694 (Purdue) and a lipoxygenase-1,2,3-null line (Kentucky);
  dehulled, ground, 20 mesh; defatted 3 x with 10 parts hexane; air-dried overnight.
- **SPI process (Figure 1):** flour : water 1 : 10 w/w (Nanopure); NaOH to pH 9.0, hold 1 h;
  centrifuge 1500 g 10 min; supernatant HCl to pH 4.5, hold 1 h; centrifuge; freeze overnight;
  lyophilise 2 d at room temperature. **No heating step anywhere.** Temperature is room temperature.
- **Spikes:** 2,4-decadienal in 2 mL chloroform into the slurry to **5.64 mM** final (chloroform alone
  gave no 2-pp increase). Ammonia (15N-labelled for the MS run; unlabelled for the tables) as ammonium
  hydroxide, slurry titrated with HCl to **8.05 mM** ammonia. Amino acids mixed with the dry defatted
  flour at **8 mM** (basis presumably the slurry). Buffer controls: dienal + ammonia or dienal + amino
  acid in aqueous buffer (pH 7 for Table 3; pH 4.5, 7, 9 for Figure 5), extracted directly with
  CHCl3/MeOH/H2O (5:10:4), solvent removed, taken up in methylene chloride.
- **Lipid extraction for analysis:** ~1 g lyophilised sample, 2 x 20 mL CHCl3/MeOH/H2O (5:10:4);
  rotary evaporation 50 C, N2 to dryness, taken up in "300 mL" methylene chloride (as printed; almost
  certainly 300 µL), stored at -15 C.
- **Internal standard and quantification:** d5/d6-2-pentylpyridine (made by the Tchitchibabine
  route from 2-picoline-d7) added as internal standard; HP G1800A GCD (EI); ratio of m/z 93 (analyte)
  to m/z 98 + 99 (IS), per Guth & Grosch 1990. **Isotope-dilution assay, one IS, one analyte.**
  Triplicate per treatment. Cool on-column injection (35 C injector, EC-5 30 m x 0.53 mm, 1.2 µm; 35 C
  -> 165 C at 10 C/min, 3 min, -> 210 C at 10 C/min, 10 min; m/z 35-250) to exclude injector-port
  artefacts; ammonia was flushed out with N2 before injection; the 2-pp odour was confirmed by sniffing
  before GC-MS.
- **Units:** "ppm" throughout Tables 1-3; the basis (per g lyophilised isolate, most plausibly, since
  ~1 g lyophilised sample was extracted; or per g slurry) is **not stated**.
- **Labelled dienal:** 5-13C-trans,trans-2,4-decadienal from 1-13C-hexanal (13C-DMF + pentyl-MgBr) and
  1-methoxy-1-Z-buten-3-yne (Barbier / Pippen-Nonaka route); mass spectra in Figure 3. The unlabelled
  spikes are presumably the same E,E isomer (not stated).
- **Statistics:** ANOVA (SAS), LSD at p < 0.05, Tukey-Kramer HSD; n = 3.

## 3. Tables re-typed

### Table 1. "Formation of 2-pentyl pyridine (ppm) when adding 2,4-decadienal." Mean ± SD, n = 3; letters compare stages within a row (a-b for dienal only, a'-b' for the parenthesised dienal + ammonia run).

Spike: 2,4-decadienal 5.64 mM (parenthesised row: + ammonia 8.05 mM). Room temperature.

| cultivar | Stage I | Stage II | Stage III |
|---|---:|---:|---:|
| Stressland | 2.330 ± 0.042 a | 2.546 ± 0.038 b | 2.600 ± 0.032 b |
| Stressland, dienal + ammonia | (4.426 ± 0.061 a') | (4.483 ± 0.066 b') | (4.550 ± 0.049 b') |
| LOX null (lipoxygenase 1, 2, 3 null) | 2.122 ± 0.038 a | 2.378 ± 0.029 a | 2.077 ± 0.026 a |
| Edison | 2.427 ± 0.051 a | 2.432 ± 0.036 a | 2.798 ± 0.032 b |
| KS4694 | 2.309 ± 0.035 a | 2.466 ± 0.029 a | 2.368 ± 0.040 a |

Authors: stage differences are "negligible compared to the increase in values by adding
2,4-decadienal"; 2-pp "was formed shortly after adding 2,4-decadienal, and remained during SPI
processing". (LOX null row: 2.378 carries letter a although it is outside the SD of 2.122; as printed.)

### Table 2. "Formation of 2-pentyl pyridine (ppm) when adding ammonia." Sampling at Stage I; mean ± SD, n = 3; letters compare columns within a row.

| cultivar | no addition | 2,4-decadienal only (5.64 mM) | ammonia only (8.05 mM) | 2,4-decadienal + ammonia |
|---|---:|---:|---:|---:|
| Stressland | 0.214 ± 0.004 a | 2.375 ± 0.038 b | 0.295 ± 0.004 c | 4.269 ± 0.056 d |
| LOX null | 0.138 ± 0.002 a | 2.015 ± 0.042 b | 0.189 ± 0.003 c | 4.189 ± 0.048 d |
| Edison | 0.172 ± 0.003 a | 2.591 ± 0.041 b | 0.194 ± 0.003 c | 4.593 ± 0.065 d |
| KS4694 | 0.198 ± 0.003 a | 2.229 ± 0.035 b | 0.236 ± 0.003 c | 4.568 ± 0.055 d |

Authors' reading: ammonia alone barely moves 2-pp (the flour's own dienal is limiting); dienal alone
gives ~2 ppm (the flour's own ammonia is enough for that); both together give more than the sum
(4.27 vs 2.375 + 0.295 - 0.214 = 2.46 for Stressland).

### Table 3. "Formation of 2-pentyl pyridine (ppm) when adding amino acids." Stressland, Stage I; mean ± SD, n = 3; letters compare the three slurry columns within a row; parenthesised = dienal + amino acid in pH 7 buffer.

Amino acid 8 mM; 2,4-decadienal 5.64 mM. Room temperature.

| amino acid | 2,4-decadienal only | amino acid only | 2,4-decadienal + amino acid (slurry) | (in pH 7 buffer) |
|---|---:|---:|---:|---:|
| arginine | 2.469 ± 0.029 a | 0.179 ± 0.003 b | 4.668 ± 0.056 c | (2.398 ± 0.033 a) |
| lysine | 2.544 ± 0.031 a | 0.184 ± 0.003 b | 3.819 ± 0.046 c | (2.466 ± 0.028 a) |
| aspartic acid | 2.444 ± 0.022 a | 0.191 ± 0.005 b | 2.491 ± 0.029 a | (2.461 ± 0.022 a) |
| asparagine | 2.512 ± 0.030 a | 0.216 ± 0.004 b | 4.025 ± 0.044 c | (2.529 ± 0.029 a) |
| glutamic acid | 2.169 ± 0.019 a | 0.188 ± 0.003 b | 2.292 ± 0.022 a | (2.206 ± 0.021 a) |
| glutamine | 2.306 ± 0.022 a | 0.200 ± 0.002 b | 3.835 ± 0.039 c | (2.332 ± 0.031 a) |
| glycine | 2.288 ± 0.028 a | 0.178 ± 0.003 b | 2.334 ± 0.025 a | (2.296 ± 0.028 a) |
| histidine | 2.299 ± 0.031 a | 0.169 ± 0.003 b | 2.318 ± 0.030 a | (2.278 ± 0.026 a) |

Each row is its own run (the "2,4-decadienal only" column varies 2.17-2.54 between rows).

### Figure 5 (FIGURE-ONLY): pH effect, buffer and Stressland slurry, dienal + ammonia

Text: "At pH 9, ammonia was the main form with some ammonium hydroxide, whereas at pH 4.5 the
ammonium ion was predominant"; "the highest level of 2-pp was found at pH 9, followed by pH 7, then
pH 4.5" in both buffer and slurry. No values printed. pKa of NH4+ ≈ 9.25 (ours): at pH 9 about 36 %
of total ammonia is NH3, at pH 7 about 0.6 %, at pH 4.5 about 0.002 %.

### Isotope evidence (Figure 4; text)

- Unlabelled 2-pp: M+ 148; ion 134 "after the first loss of methylene group"; base peak m/z 93 (the
  ion used for quantification).
- With 5-13C-2,4-decadienal: 2-13C-2-pentylpyridine, M+ 149, first-loss ion 135, base peak m/z 94.
- With 15N-ammonia: 1-15N-2-pentylpyridine, M+ 149, 135, base peak 94.
- Therefore: dienal C5 becomes pyridine C2 (the carbon bearing the pentyl), consistent with the
  Buttery 1977 mechanism the paper quotes: ammonia condenses with the aldehyde (Schiff base), ring
  closure (N onto C5), dihydropyridine, oxidation to 2-pp. Complete atom map (ours, from the label and
  the carbon count): dienal C1 -> pyridine C6 (H), C2 -> C5, C3 -> C4, C4 -> C3, C5 -> C2, C6-C10 ->
  pentyl; N from NH3.

### The isolate's own precursor levels

**Not measured.** The paper prints no 2,4-decadienal or ammonia concentration for the flour; it says
only that "the amount of 2,4-decadienal was very limited, but ammonia existed in soybean defatted
flours" and cites Arai et al. 1966 for free ammonia in the basic fraction of raw soybean. The
observable proxy is the no-addition 2-pp baseline: 0.138 (LOX null) to 0.214 ppm (Stressland), Table 2.

## 4. Routes and numbers the repository can use

Conditions for every row: defatted soy flour slurry 1 : 10 w/w in water (or aqueous buffer), room
temperature, pH 9 (1 h) then 4.5 (1 h) for the slurry; isotope-dilution GC-MS; n = 3; ppm (basis not
stated).

| route | reactant -> product | mechanism as drawn / stated | measured numbers (units, conditions) | evidence class |
|---|---|---|---|---|
| ZB-2PP | **2,4-decadienal + NH3 -> 2-pentylpyridine** | stated (Buttery 1977): Schiff base, ring closure, dihydropyridine, oxidation; not drawn; atom map fixed by 5-13C and 15N labels (C5 -> ring C2, N -> ring N) | dienal 5.64 mM alone: +2.0 to +2.4 ppm over baseline (Table 2: 2.015-2.591 vs 0.138-0.214); + ammonia 8.05 mM: 4.19-4.59 ppm; ammonia alone: +0.02 to +0.08 ppm | level_only; within_study_ratio; label-confirmed atom map |
| ZB-PH | same, pH 4.5 / 7 / 9 | NH3, not NH4+, is the nucleophile | ordering 9 > 7 > 4.5 in buffer and slurry; no numbers | figure_only |
| ZB-RT | same, at room temperature | "the synthesis of 2-pp is a spontaneous reaction from 2,4-decadienal and ammonium hydroxide" | Table 1: formed by Stage I (within the 1 h pH 9 hold), flat thereafter | level_only |
| ZB-AA-SLURRY | dienal + {Arg, Lys, Asn, Gln} in flour slurry -> more 2-pp; {Asp, Glu, Gly, His} -> none | flour "may provide the matrix for deamidation of amino acids, possibly by enzymatic means"; the Arg/Lys effect is unexplained | Table 3: Arg 4.668, Asn 4.025, Gln 3.835, Lys 3.819 vs dienal-only 2.3-2.5 ppm; Asp/Glu/Gly/His within 0.05-0.12 of dienal-only | level_only; within_study_ratio |
| ZB-AA-BUFFER | dienal + any of the eight amino acids in pH 7 buffer, room temperature -> no increase | "amino acids do not contribute to the synthesis in buffer solution at room temperature" | Table 3 parentheses: 2.206-2.529 vs 2.169-2.544 dienal-only | null result (level_only) |
| ZB-BASELINE | untreated isolate 2-pp | flour's own dienal (limiting) + own ammonia | 0.138-0.214 ppm across four cultivars; LOX null lowest | level_only |
| ZB-LOX | lipoxygenase-null flour | less enzymatic dienal | baseline 0.138 vs 0.172-0.214; with dienal spike the LOX-null line reaches the same ~2 ppm as the others | within_study_ratio |

Within-study ratios worth registering (Stressland, Stage I):
- (dienal + NH3) / dienal-only = 4.269 / 2.375 = 1.80; (dienal-only - baseline) / (NH3-only - baseline)
  = 2.161 / 0.081 = 27: the flour is ammonia-rich and dienal-poor for this reaction.
- Best amino acid (Arg) adds 2.20 ppm over dienal-only; ammonia at a similar molarity (8.05 vs 8 mM)
  adds 1.89 ppm: arginine in the slurry is at least as effective as ammonium hydroxide, unexplained.
- Conversion (ours, order of magnitude only, because the ppm basis is unstated): if ppm is µg per g of
  lyophilised isolate and the isolate is ~20 % of the flour mass, 2.16 ppm of 2-pp ≈ 14.5 nmol/g isolate
  ≈ 3 nmol per g of flour, against a dienal spike of 5.64 mM x 10 mL water per g flour = 56 µmol/g:
  about 5 x 10^-5 of the dienal. Even if the basis were per g slurry the conversion is < 10^-4. Treat as
  "well below 0.1 %", not as a number.

## 5. Rule sketches (repository suggestions, not the paper's)

Registry: **2-pentylpyridine has no key** in `data/keys/compounds.yml`; species `DECADIENAL`
(`CCCCC/C=C/C=C/C=O`), `Gln`, `Asn`, `Lys` exist in `data/species/structures.yml`; no arginine, no
ammonia species.

**S1. 2,4-decadienal + NH3 -> 2-pentylpyridine (net; atom map from this paper).** Ring = N + dienal
C1-C5; pentyl = C6-C10; N bonds C1 and C5; loses H2O and H2 (the oxidation of the dihydropyridine).
- positive: `CCCCC/C=C/C=C/C=O` + `N` -> `CCCCCc1ccccn1` + `O`
- label check the rule should reproduce (label at dienal C5): `CCCCC/[13CH]=C/C=C/C=O` + `N` ->
  `CCCCC[13c]1ccccn1` (2-13C-2-pentylpyridine). A SMIRKS that maps C5 onto the alkyl-bearing ring
  carbon passes; one that maps C1 onto it fails this paper.
- negative: hexanal `CCCCCC=O` + `N` -> no fire; DECADIENAL + `NCC(=O)O` (glycine) -> no 2-pp (Table 3, buffer and slurry); DECADIENAL + `NC(CCC(=O)O)C(=O)O` (Glu) -> no 2-pp.
- conditions to record: aqueous, 20-25 C, pH 9 >> 7 > 4.5; also 150-180 C (Du 2023, Zamora 2020, Kim 1998).

**S2. NH3 supply in a soy/pea slurry.** Two sources the paper distinguishes: (a) the flour's own free
ammonia (level not measured); (b) deamidation of Asn/Gln (and something from Arg/Lys) that happens in
the slurry but not in buffer. For the rule layer, S5 of `zamora2020_extraction.md` (Gln -> NH3) covers
(b) for Gln/Asn; arginine -> NH3 (arginase-like or alkaline hydrolysis to ornithine + urea -> NH3) would
be `proposed` with no source drawing it: `NC(CCCNC(N)=N)C(=O)O` -> `NC(CCCN)C(=O)O` + `NC(N)=O`; urea
-> 2 NH3 + CO2 (Zamora 2020 Table 1 shows urea is a working NH3 source at 180 C).

**S3. Speciation gate (not a SMIRKS).** The rule should carry "reactive N = NH3, not NH4+"; at pH 6-7
(the roadmap's cooks) only ~0.2-0.6 % of total ammonia is NH3.

## 6. Flags

1. **"ppm" with no basis** (per g lyophilised isolate? per g slurry?). Ratios within a table are safe;
   conversions to yield on the dienal are order-of-magnitude at best (§4).
2. **Figure 5 (pH) has no numbers**; the pH ordering is all that can be cited.
3. **The isolate's own 2,4-decadienal and ammonia were not measured**; the baseline 2-pp is the only
   proxy. The roadmap's "isolate carries 1-3 % lipid" cannot be joined to a dienal level from here.
4. **Amino-acid spikes were mixed into the dry flour** while dienal and ammonia went into the slurry;
   the "8 mM" for amino acids is nominal.
5. **"300 mL methylene chloride"** for the extract is printed twice; almost certainly 300 µL (a 1 g
   sample concentrated for GC). Does not affect the reported ppm.
6. **Room-temperature claim rests on the process itself having no heating step** (pH 9 hold 1 h, pH
   4.5 hold 1 h, lyophilisation); the authors excluded the 210 C injector as the site of formation by
   N2-flushing ammonia and by cool on-column injection. A kinetic rate at 25 C cannot be extracted: no
   time series inside the 1 h hold.
7. **Arg and Lys raise 2-pp in the slurry but not in buffer, and glycine does not**: the N-source rule
   "amide side chain only" (Kim & Ho 1998) is too narrow for the matrix; the paper leaves the Arg/Lys
   route unexplained ("possibly by enzymatic means").
8. **Isotope dilution with a d5/d6 IS and m/z 93 vs 98/99**: a good quantification for one analyte;
   no other compound (3-pentylpyridine, hexanal, 2-pentylfuran) is reported.
9. **Table 1 LOX-null letters** (2.122 a, 2.378 a, 2.077 a) are as printed.
10. **The dienal isomer of the unlabelled spike is not stated**; the labelled one is trans,trans.
