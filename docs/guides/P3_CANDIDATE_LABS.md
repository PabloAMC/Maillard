# Candidate laboratories for the factorial lever study (plan P3)

*2026-11-14. Compiled from public sources (Europe PMC, publisher and university pages). The owner decides and makes
every contact; nothing has been sent to anyone.*

**Purpose.** A shortlist from which the owner picks the one or two laboratories that would run the 2^3 factorial in
[FACTORIAL_LEVER_STUDY.md](FACTORIAL_LEVER_STUDY.md) (glucose vs ribose x cysteine vs thiamine x pH 5.5 vs 6.5, 15 %
pea protein isolate, sealed vials, 140 C 20 min, about 22 pots per laboratory, SIDA GC-MS for MFT, FFT, disulfides,
Strecker aldehydes, methional, pyrazines, hexanal/nonanal, DMHF, furfural).

**Criteria.** (1) An active group (papers 2021 - 2026) doing thermal-flavour / Maillard volatile chemistry with
quantitative GC-MS, ideally SIDA, and thiols above all; (2) plant-protein or meat-analogue flavour work is a plus;
(3) absent from the keyed record, meaning not in `docs/reference/contrast_benchmark/lab_ids.yml` / `split.yml` and not in
`docs/reference/lever_table/new_labs.yml`. A lab id found only in the screening logs (`field_checks/*_screening_log.json`,
decision "no") is still new. Close calls are flagged, not dropped; (4) an academic or institute group that takes
collaborations. EU groups are marked.

**Excluded on criterion 3** (good fits already in the record): TUM / Leibniz-LSB (Schieberle, Hofmann, Steinhaus, and
Dawid, who works in Hofmann's chair: `dfa_garching`); Reading (Parker, Elmore, Methven: `reading_mottram_parker`);
Hacettepe (Gokmen: `hacettepe_gokmen`, TEST); Ghent (De Meulenaer was merged into `ghent_de_kimpe`); Wageningen FQD
(Fogliano); BTBU (this includes its Key Laboratory of Green and Low-carbon Processing Technology for Plant-based Food,
Li Jian / Zhang Wei); Jiangnan (Chen Jie / Zeng Maomao, and Zhang / Cui); A*STAR SIFBI; FRC Nutrition Moscow; Graz
(`graz_esterbauer`).

## Candidates

| # | Group (PI, institution, country) | Relevant capability | Representative recent papers (DOI) | Why it adds a new laboratory | Caveats |
|---|---|---|---|---|---|
| 1 | **Chetschik group, Food Chemistry, Zurich University of Applied Sciences (ZHAW), Wädenswil, Switzerland** (Horizon-associated) | Sensomics with SIDA quantitation; quantitates key odour-active thiols; plant-based alternatives (egg) | 10.1021/acs.jafc.5c04304 (2025, key odorants of plant-based vs chicken egg, SIDA); 10.1021/acs.jafc.5c14203 (2026, first quantitation of key thiols in hemp, HPLC-MS/MS) | Not in any lab file; a SIDA laboratory outside the Munich group | Co-authors with Steinhaus (`dfa_garching`) on chocolate and hemp papers, and Blank (ex-Nestlé) is on the egg paper. Confirm that the method lead and the instrument are ZHAW's own. Cooked-matrix Maillard work is so far on chocolate and cocoa, not model pots |
| 2 | **Flores group, IATA-CSIC, Valencia, Spain** (EU) | Meat and meat-analogue aroma; GC-MS and GC-O; texturised pea protein; nitrogen and sulfur precursor effects on aroma | 10.1021/acs.jafc.3c08432 (2024, texturised pea protein fermentation for meat-analogue aroma); 10.1016/j.meatsci.2022.108896 (2022, nitrogen and sulfur precursors and sausage aroma) | `iata_csic_flores` appears only in the batch-3 screening log (decision "no"), not keyed | Does not use SIDA for thiols (it uses internal standards and GC-O). The FFT/MFT method would need to be set up or checked. Its main line is fermentation, not dry heating |
| 3 | **Peterson group, Flavor Research and Education Center, The Ohio State University, USA** | Flavour chemistry of pea protein isolates (volatile and non-volatile); long Maillard record (carbonyl trapping at Penn State) | 10.1016/j.foodres.2026.121033 (2026, aversive aroma compounds in pea protein isolates); 10.1111/1750-3841.70925 (2026, drivers of acceptance in pea protein isolates) | Not keyed. `pennsylvania_peterson` (Penn State era) appears only in a screening log | The same university as `ohio_state_barringer` (a different group); the owner should decide whether that counts as one laboratory. Recent output is about off-flavour and acceptance, not meaty thiols |
| 4 | **Cadwallader group, Food Science and Human Nutrition, University of Illinois Urbana-Champaign, USA** | Classic aroma-extract dilution plus SIDA quantitation of key odorants, including roasted and thermally generated ones | 10.1021/acs.jafc.5c01690 (2025, key odorants of peanut and sunflower seed butters); 10.1093/jas/skac305 (2022, Maillard products in pet foods) | Not in any lab file | No plant-protein matrix work found 2021 - 2026. Do not confuse with D. C. Cadwallader (NC State, Drake group) |
| 5 | **Fraunhofer IVV, Freising, Germany** (EU; sensory analytics and food process development; Buettner directs the institute) | Pea and lupin protein aroma profiles; GC-O tradition with SIDA in the Buettner school | 10.1002/fsn3.71889 (2026, pea protein isolate aroma after high-pressure processing); 10.3390/foods10061257 (2021, fermentation and the aroma of lupin protein) | Not in any lab file | It shares a town with TUM / Leibniz-LSB, and Buettner trained in the Schieberle school. Check that it is a separate method and instrument. No thermal Maillard model series found |
| 6 | **Zorn group, Institute of Food Chemistry and Food Biotechnology, Justus Liebig University Giessen, Germany** (EU) | Meat-like flavour from fungi and mycelium for meat alternatives; GC-MS-O identification | 10.1021/acs.jafc.3c00542 (2023, meat-like flavour from submerged-cultivated Laetiporus montanus) | Not in any lab file | Its focus is biotechnology, not sugar-amino-acid heating. Thiol quantitation by SIDA is not shown |
| 7 | **Fisk group, International Flavour Research Centre, University of Nottingham, UK** | Plant-based burger flavour: Maillard precursors (cysteine/cystine, dextrose/ribose) and metal ions; pulse volatiles | [eprints.nottingham.ac.uk/80016](https://eprints.nottingham.ac.uk/80016) (2025 MRes thesis, flavour pathways in plant-based meat analogues); 10.1016/j.foodchem.2022.133569 (2022, vetch volatiles by processing) | `nottingham_fisk` appears only in the screening logs | `new_labs.yml` records Fisk as a collaborator of the Jiangnan Zhang / Cui laboratory (`ho_network`), so under the collaboration rule a joint study would count as the Jiangnan laboratory. Accept only a study run by Nottingham alone. Its strength is headspace and APCI-MS; no SIDA for thiols found |
| 8 | **Kilcawley group, Teagasc Food Research Centre, Moorepark, Ireland** (EU) | Quantitative volatile GC-MS (lipid oxidation, dairy, grilled beef); pea-protein bread | 10.3390/foods14244239 (2025, volatile profile of beef by grilling temperature); 10.1016/j.foodchem.2026.151162 (2026, protein-enhanced pea bread) | Not in any lab file | Its strength is the off-note panel (hexanal, nonanal, 2-pentylfuran), not thiols. FFT/MFT would need a SIDA set-up |
| 9 | **Drake group, Dept. of Food, Bioprocessing and Nutrition Sciences, North Carolina State University, USA** | Pea protein aroma (GC-O, GC-MS, sensory); sulfur off-flavours in heated milk | 10.1016/j.foodchem.2022.134998 (2023, predominant aroma components of pea protein concentrates and isolates) | Not in any lab file | Mainly a sensory and dairy group; Maillard model systems are not its core work |
| 10 | **Salles / Guichard, CSGA (INRAE, CNRS, Institut Agro, Université Bourgogne), Dijon, France** (EU) | Faba-bean odour-active compounds (GC-O, GC-MS) | 10.1016/j.foodres.2022.112260 (2023, origins of volatile and odour-active compounds in faba bean fractions) | Not in any lab file (the INRAE units in the record are Clermont and Nantes) | The same parent institute as two keyed INRAE laboratories, though a different unit. Its work is on raw pulses, not Maillard heating |

## Ranked top 3

1. **ZHAW Chetschik (Switzerland).** It is the only new group found that runs SIDA on thiols and also works on plant-based flavour. Only check that its method is independent of Steinhaus's laboratory.
2. **IATA-CSIC Flores (Spain, EU).** It works on texturised pea protein and on sulfur-precursor aroma, so it is the closest fit to the matrix. The caveat is that its FFT/MFT quantitation would need a SIDA check.
3. **Ohio State Peterson (USA).** It studies pea-protein isolate flavour and has a long Maillard record. The owner must decide whether sharing a university with `ohio_state_barringer` matters.

The two-laboratory default could be ZHAW plus IATA-CSIC: different countries, methods and lineages, with both in or associated with Horizon Europe.
