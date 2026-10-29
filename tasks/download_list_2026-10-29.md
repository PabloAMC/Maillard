# Download list — missing supplements and patents (sulfur lane + trunk)

Audit date 2026-10-29. Nothing downloaded. Supplements on disk already (excluded): Kang2026, Thong2024,
Wang2022 (x2), Wang2026, Yiltirak2026, Zhang2025, brehm2020, zhou2023b (Zhou 2023 JAFC SI).
Content descriptions come from the dossiers (data/lit/extraction_dossiers/<key>_extraction.md);
"unknown" = the dossier does not say what the SI holds.
Confidence = my probability the file contains numbers the model can use (tables, not just plots).

## A. Supplements of papers we have (sulfur lane)

A1. Zhang 2024 FRI 182:114149 — https://doi.org/10.1016/j.foodres.2024.114149 (Elsevier "Appendix A" mmc)
  Holds: Fig. S1 dynamic MFT/MMFT/DMTS profiles, Fig. S2 single-factor temperature sweep, Fig. S3 RSM,
  Tables S1-S7 (S5 standard curves, S6 zero-order kinetics, S7 regressions; S6/S7 carry the k's whose units are missing).
  Why: only sulfur time series at 115 C with a quantified MFT sink (Cys vs GSH); temperature axis; resolves units.
  Already flagged in k3 as "single highest-value retrieval". Confidence usable numbers: 0.6 (tables S5-S7 yes; S1/S2 may be figures).

A2. Zhang 2026 Food Chem 148681 — https://doi.org/10.1016/j.foodchem.2026.148681
  Holds: Tables S6-S8 per-compound volatiles + GC-O FD at 130/140/150 C (GSH + glucose, pH 6.5, 8 points to 120 min per main text).
  Why: a fourth MFT temperature ladder, peptide system, paired with 19 printed k / Ea for non-volatiles. Confidence 0.75 (tables).

A3. Wang 2026b Food Chem 149913 — https://doi.org/10.1016/j.foodchem.2026.149913
  Holds: Table S1, 22-compound quantification (Xyl-Cys, 120 C, 0-150 min); main text prints only OAVs, Fig 4A broken axis.
  Why: "first external thiol-and-dimer time course" (MFT, FFT, dimers). Confidence 0.8.

A4. Zhai 2023 Food Chem 134420 — https://doi.org/10.1016/j.foodchem.2022.134420
  Holds: Table S2 MGO-GO-Cys validation model (MGO 10 + GO 10 + Cys 10 mM, pH 7, 100/120/140 C, 20-120 min) volatiles; Fig S2 Cys-alone degradation; Table S1 ROAV; Fig S1 map.
  Why: dicarbonyl + Cys -> sulfur volatiles at 3 T (tests the cysteine-Strecker / H2S step bypassing sugar); also bears on the Kang S4 provenance question. Confidence 0.6.

A5. Zhai 2023 JAFC (3c04166) — https://doi.org/10.1021/acs.jafc.3c04166
  Holds: Table S2 calibration curves of volatiles + main flavour compounds (ug/L) from thermal MGO-GO-Cys model; Fig S1/Table S1 TTCA/ARP identity; Fig S2 pathway.
  Why: decides whether Kang S4 numbers were externally calibrated (integrity item 1 in STRATEGY); MGO-GO-Cys data. ACS SI is free without subscription. Confidence 0.7.

A6. Zhai 2023c JAFC 71:16260 (3c06506) — https://doi.org/10.1021/acs.jafc.3c06506
  Holds: Fig S2 TTCA and ARP yields vs vacuum-dehydration time (Cys+Xyl 82.7 mM, 90 C); Fig S3 TTCA/ARP after spray drying; Table S2 calibration.
  Why: TTCA/ARP vs time under dehydration (water-activity axis of the Cys-pentose route). Figures, so digitisation. Confidence 0.4.

A7. Feng 2022 JAFC (2c02949) — https://doi.org/10.1021/acs.jafc.2c02949
  Holds: Table S1 calibration curves + LOD/LOQ (the absolute basis of the "second sulfur temperature ladder", RG-ARP 100-120 C pH 7, 120 min), Tables S2-S3 OAVs, Fig S2 Rib+GSH totals.
  Why: makes the second ladder auditable (absolute vs semi-quant); currently "unauditable". Confidence 0.7.

A8. Zhao 2019 Food Chem (2018.08.096) — https://doi.org/10.1016/j.foodchem.2018.08.096
  Holds: Table S1 pH, A420, sensory after heating for six 13C-TTCA pots (pH 5.5 phosphate, 140 C 30 min).
  Why: measured final pH per pot (pH lane). Endpoint only. Confidence 0.8.

A9. Zhang 2023b FRI 113200 — https://doi.org/10.1016/j.foodres.2023.113200
  Holds: Fig S2 final pH values, Fig S1 Met residuals with additives, Table S3 12 VSCs over 56 d storage at 50 C (Xyl-Cys-Met, pH 4.9).
  Why: measured final pH (lane 1) and a slow FFT-decay storage series. Confidence 0.5 (pH as figure; Table S3 numbers likely).

A10. Gigl 2021 JAFC (1c06163) — https://doi.org/10.1021/acs.jafc.1c06163
  Holds: Table S2 qHNMR time-dependency table (mmol/L and % free, 20 odorants incl. FFT, 10 time points, RT), Fig S3 FFT/NAC/MeFFT substitution time courses, Table S1 aroma model.
  Why: printed numbers for thiol loss to melanoidins (replaces +-4 pp digitisation). RT only. Confidence 0.85.

A11. Zhang 2022 FRI 111081 — https://doi.org/10.1016/j.foodres.2022.111081
  Holds: Tables S1-S6, Fig S1 (MFT/FFT per oil-MRP, Xyl-Cys pH 7.2, 120 C 90 min).
  Why: endpoint MFT/FFT with oil phases; HS-SPME response factor 1, so semi-quant. Confidence 0.6 numbers, low value.

A12. Yao 2025 Food Chem 141143 — https://doi.org/10.1016/j.foodchem.2024.141143
  Holds: Table S4 VSCs identified in four models (methional alone; Met+Glc; +-EGCG; 100 C, 60-120 min), Table S2 model calibration curves.
  Why: methional -> MeSH/DMDS/DMTS chain; the time courses themselves are Fig. 4 (main text). Confidence 0.4.

A13. Zhou 2022 JAFC, Ala-Xyl ARP + Cys — https://doi.org/10.1021/acs.jafc.2c07026
  Holds: Table S1 calibration equations, Fig S1 ARP characterisation, Fig S2 A420.
  Why: calibration basis only. Confidence 0.3; low value.

A14. Zhou 2025 LWT 117469 (coffee sulfury/roasty fading) — https://doi.org/10.1016/j.lwt.2025.117469
  Holds: Tables S1-S3, 86-compound x condition concentration matrix (storage of roasted coffee).
  Why: thiol (FFT etc.) loss in a real matrix over storage. Confidence 0.7 numbers, medium relevance.

A15. Li 2016 Free Radic Biol Med — https://doi.org/10.1016/j.freeradbiomed.2016.05.019
  Holds: Tables 2-7 (thiol + 1,2-benzoquinone rate constants vs pH), Table 8 (Gly), Table 9 (Lys).
  Why: thiol sink to quinones (polyphenol matrices); second-order k vs pH. Adjacent. Confidence 0.85.

A16. Stack 2018 Chem Res Toxicol (CC-BY) — https://doi.org/10.1021/acs.chemrestox.7b00239
  Holds: Figs S3/S4 Arrhenius plots (quinone + NAC/GSH). Adjacent, low value. Confidence 0.4.

## B. Supplements of papers we have (trunk: sugars, ARP, deoxyosones, dicarbonyls, furans, acids, pH)

B1. Sen & Gokmen 2022 Food Chem 133583 — https://doi.org/10.1016/j.foodchem.2022.133583
  Holds: Tables S1-S5 + Figs S1-S5, "the per-sample concentration time courses" (sucrose-rich nuts/seeds roasting; sugars, 3-DG, HMF ...).
  Why: tabulated multiresponse time courses (matrix HMF branch bound). Confidence 0.75.

B2. Piornos 2025 Food Chem 141532 — https://doi.org/10.1016/j.foodchem.2024.141532
  Holds: Fig S1 full concentration data panels a-s (12 curing runs, 65/78/90 C, green malt), Fig S2 predicted Int1, Table S1 moisture.
  Why: full multiresponse dataset of a Balagiannis/Parker fit + the moisture confound. Figures => digitise. Confidence 0.5.

B3. Xia 2022 JAFC (2c06639) — https://doi.org/10.1021/acs.jafc.2c06639
  Holds: Fig S5, 120 C 1-DG/3-DG scatter points + Arrhenius plot (Glc + Gly/diGly/triGly 0.2 M, pH 7.5).
  Why: trunk's own pot; second temperature for deoxyglucosone formation. Figure points. Confidence 0.5.

B4. Gursul Aktag 2020 Food Chem 126620 — https://doi.org/10.1016/j.foodchem.2020.126620
  Holds: Table S1 the 4 C data, S2 free amino acids, S3 pH and Brix, S4 comprehensive-model rate constants, Fig S1 mass balance.
  Why: pH-3.4 sugar/dicarbonyl/HMF storage network; S3 measured pH; S4 constants. Confidence 0.85 (tables).

B5. Nguyen 2016 Food Chem (2015.06.110) — https://doi.org/10.1016/j.foodchem.2015.06.110
  Holds: Fig S1 pH and mass balance; Fig S2 the ODEs; Table S1 pairwise comparisons (caseinate + Glc/Lac, pH 6.8, 120/130 C).
  Why: pH drift under phosphate + rate laws for printed k's. Confidence 0.5.

B6. Kocadagli & Gokmen 2016 JAFC (6b01862) — https://doi.org/10.1021/acs.jafc.6b01862
  Holds: Figs S1-S5 incl. S5 model-discrimination fit (3-DG-only HMF route rejected).
  Why: checks the HMF-from-fructose constraint; mostly fits/plots. Confidence 0.3.

B7. Goncuoglu Tas & Gokmen 2017 Food Chem (2016.11.159) — https://doi.org/10.1016/j.foodchem.2016.11.159
  Holds: Table S1 pH + moisture during hazelnut roasting, S2 amino acids, S3 colour, S4 Ea / k_b (authors disclaim the Ea).
  Why: measured pH/moisture trajectories in a roast; S4 of low trust. Confidence 0.8 numbers, medium value.

B8. Hamzalioglu 2026 JAFC (5c14296; PMC13220311 open) — https://doi.org/10.1021/acs.jafc.5c14296
  Holds: Fig S1 mass balances, Tables S1-S2 parameter correlation matrices (milk, 110-140 C).
  Why: correlation matrices let the printed k's be used as a joint posterior prior. Confidence 0.8, medium value.

B9. Zhou 2024 JAFC (4c03706) — https://doi.org/10.1021/acs.jafc.4c03706
  Holds: Table S1 Gly isotopomers, Fig S1 Gly-model pyrazine levels, Fig S2 Arrhenius plot (Ala + GO/MGO -> pyrazine, 100-120 C).
  Why: resolves the Ea / 110 C k inconsistency. Confidence 0.4.

B10. Cai 2024 FRI 114591 — https://doi.org/10.1016/j.foodres.2024.114591
  Holds: Table S1 per-timepoint means (5 responses, Glc + Pro, 100/120/140 C; 2-AP route), Table S2 fit performance, Fig S1 mass balance.
  Why: tabulated time courses for a 10-k glucose/proline model (Strecker/pyrroline branch, adjacent to trunk). Confidence 0.8.

B11. Quan 2020 Food Chem 127387 — https://doi.org/10.1016/j.foodchem.2020.127387
  Holds: Suppl. Fig 1 rate equations, Fig 2, Tables 1-2 validation (Lys/Asn/Trp + Glc, GO, MGO, 100/130 C).
  Why: would rescue the 88 refused constants (units/rate laws). Confidence 0.5.

B12. Hidalgo 2016 Food Chem (04.032) — https://doi.org/10.1016/j.foodchem.2016.04.032
  Holds: Suppl. Figs S-1..S-5, the k values (PEA/PAC from Phe + 2,4-decadienal, 100-180 C). Lipid-carbonyl Strecker; adjacent. Confidence 0.4.

B13. Mohsin 2018 Food Chem (2017.11.115) — https://doi.org/10.1016/j.foodchem.2017.11.115
  Holds: Table S-1 melanoidin C/N at 10 temperatures 130-200 C; Figs S-5/S-6 CHNO.
  Why: melanoidin stoichiometry check (MEL_C/N). Confidence 0.85, niche value.

B14. Hwang 1995 JAFC (jf00049a033) — https://doi.org/10.1021/jf00049a033
  Holds: 13-page supplementary material: abundances behind the Gly vs tested-amino-acid pyrazine split. Old ACS microfiche SI; may need ACS request. Confidence 0.4 (availability).

## C. Patents

Already in corpus (dossiers): US 10,076,126 (Nestec, Xyl+Cys 110 C 0.5-4 h MFT/FFT figs), US 8,920,862 (Nestec, 95 C 6 h SIDA endpoint),
US 9,943,096 + EP 2 943 072 B2 (Impossible, heme; GC-MS presence/0-S-M-L), US 2025/0221436 A1 (Area-%), US 3,904,655 (Lever, recipes only),
US 5,760,005 (NAC H2S floor). ON DISK BUT NO DOSSIER: data/articles/us10039306.pdf (Impossible, Vrljic et al. 2018, "Methods and compositions for consumables") — read before downloading anything else from the Impossible family.

Web search (2026-10-29) for new patents with quantitative Maillard data:

C1. US 9,005,689 B2, Nestec, priority 2009-08-13, "Process for the preparation of a flavor active composition" — https://patents.google.com/patent/US9005689B2/en
  Holds (Google Patents text): Glc + Glu + Cys 0.4 M each in 0.5 M phosphate pH 7, 100 C: FFT 0.4 ug/mol Cys at 30 min, 2.3 at 2 h;
  multistep (Cys 0.8 M pre-heated alone at pH 9 or pH 3, 100 C 6 h) FFT 14-23 ug/mol Cys; 2-acetyl-2-thiazoline 17 vs 37, ADHT 29 vs 114 ug/mol Cys;
  H2S from Cys ~4x higher at pH 9 than pH 3 after 5 h at 100 C; time courses in Figs 4, 5, 8 (figure only).
  Why: cysteine-alone degradation vs pH (H2S source) + dicarbonyl-pre-formed route, in phosphate, stated conditions. Confidence 0.6 (endpoints printed; time courses figures).

C2. US 2022/0256899 A1, Unilever (Conopco), priority 2018-10-25, "Process for the preparation of a flavouring" — https://patents.google.com/patent/US20220256899A1/en
  Holds: GSH yeast extract + nucleotide (IMP/GMP/CMP/UMP) yeast extract + lactic acid, 100 C 7 h: NMR before/after (GSH 0.27 -> <0.01 and 0.96 -> 0.26 mg/g; IMP 16.37 -> 15.55, 14.76 -> 13.29 mg/g), measured pH per run (4.3-6.7), time-temperature equivalence law t = 4096 e^(-bT), 90-160 C.
  Why: only quantitative GSH/IMP consumption data found in the flavour-patent space; matrix is yeast-extract paste (low water), so a qualitative check. No thiol numbers. Confidence 0.7 numbers, low-medium value.

C3. US 11,802,137 B2, Jiangnan University (Cui, Zhan et al.), filed 2020 — https://patents.google.com/patent/US11802137B2/en
  Holds: Cys 8 kg + Xyl 19.8 kg in 800 kg water (= 82.6 mM / 165 mM, pH 7.5, 100 C) intermediate (ARP/TTCA) HPLC content vs time (Fig 10) and A420 vs stage-1 time (Figs 1, 4, 6, 8); text gives only optimum times.
  Why: same charge as Zhai 2019 (likely the same lab data); figure only. Confidence 0.3.

C4. US 5,149,840 A, Unilever/Quest (now Givaudan), priority 1989 — https://patents.google.com/patent/US5149840A/en
  Holds: HDMF yields from rhamnose + Lys/Hyp in phosphate, 100 C 3 h (38 % direct, 61-63 % in mixture); xylose + Pro in glycerol/phosphate 65 C 7 h; no pH, no time course.
  Why: phosphate-catalysed furanone yields (trunk, furanone branch); rhamnose not a model sugar. Confidence 0.5 numbers, low value.

Searched and rejected (no quantitative Maillard data in text): US 6,716,468 (Quest; 20 ppm MFT in a spray-dried Xyl-Cys flavour, no conditions), US 2023/0180803 (Unilever, shear cell, sensory only),
US 11,903,401 (Jiangnan; ARP yields with EGCG, no cysteine), EP 2049509 / US 7,989,013 (Firmenich furyl thioalkanals, sensory), US 9,700,067 / 9,826,772 / 10,314,325 / 10,327,464 / 11,219,232 / 11,224,241 (Impossible family, same GC-MS presence style as US 9,943,096), US 12,408,684 (claims only).
No Givaudan, Firmenich, IFF, Symrise or Motif patent with a printed MFT/FFT concentration table was found.

## D. Ranking (expected value for the model)
1 A3 Wang 2026b Table S1 | 2 A1 Zhang 2024 SI | 3 A2 Zhang 2026 S6-S8 | 4 A5 Zhai 2023 JAFC S2 | 5 A4 Zhai 2023 FC S2
6 A10 Gigl 2021 Table S2 | 7 B1 Sen 2022 S1-S5 | 8 B4 Gursul 2020 S1-S4 | 9 A7 Feng 2022 S1 | 10 B10 Cai 2024 S1
11 A8 Zhao 2019 S1 (pH) | 12 C1 US 9,005,689 | 13 A9 Zhang 2023b (pH) | 14 B2 Piornos 2025 | 15 B3 Xia 2022 S5
16 B7 Goncuoglu Tas 2017 S1 | 17 B5 Nguyen 2016 | 18 B8 Hamzalioglu 2026 | 19 A14 Zhou 2025 LWT | 20 A6 Zhai 2023c
21 B11 Quan 2020 | 22 A15 Li 2016 | 23 C2 US 2022/0256899 | 24 A12 Yao 2025 | 25 B13 Mohsin 2018
Below cut: A11, A13, A16, B6, B9, B12, B14, C3, C4.

## E. Added 2026-10-29: the lipid, acrylamide, aroma and matrix ledgers

Items already above are not repeated (Cai 2024 = B10, Piornos 2025 = B2, Zhou 2024 = B9, Zhou 2022 = A13,
Goncuoglu Tas 2017 = B7, Quan 2020 = B11, Hidalgo 2016 = B12, Zhang 2022 = A11).

E1. Bi 2022 Food Chem 133044 — https://doi.org/10.1016/j.foodchem.2022.133044
  Holds: Tables S1-S4, Fig. S1. The SHIPPED pea-protein binding constants rest on Table S3, which the dossier quotes only
  from the text. Why: an integrity check on constants the matrix layer ships. Confidence numbers 0.8.
E2. Utz 2022 Foods 11:412 (open access) — https://doi.org/10.3390/foods11030412
  Holds: Tables S4 (SIDA concentrations of ~27 odorants in ten pea-protein isolates) and S5 (OAVs).
  Why: the off-flavour baseline of commercial pea isolates, by isotope dilution. Confidence 0.85.
E3. Zamora 2015 JAFC (5b02619) — https://doi.org/10.1021/acs.jafc.5b02619
  Holds: Supporting Figs. S1-S6 (O2, pH, dose and water series; k(T) time courses) of lipid-carbonyl Strecker degradation.
  Why: the only simulable time courses for that cross-route, with temperature. Figures, so digitisation. Confidence 0.5.
E4. Trikusuma 2020 Food Chem 126082 — https://doi.org/10.1016/j.foodchem.2019.126082
  Holds: Table S1 standard-addition slopes (matrix / water response factors for meat-analogue volatiles).
  Why: a measured matrix factor for meat analogues. Confidence 0.7.
E5. Conti 2025b FRI 116938 — https://doi.org/10.1016/j.foodres.2025.116938
  Holds: Supplementary Table 3, the numbers behind Fig. 2 (odorant release into saliva over time). Confidence 0.8, medium value.
E6. Goncuoglu 2026 Food Chem X 104019 (open access) — https://doi.org/10.1016/j.fochx.2026.104019
  Holds: Tables S3 (composition series) and S4 (screw speed and temperature series) of plant-protein extrusion volatiles.
  Why: the extrusion process axis (semi-quantitative). Confidence 0.7 numbers, medium value.
E7. Xin 2026b FRI 119010 — https://doi.org/10.1016/j.foodres.2026.119010
  Holds: Tables S1, S4-S7 (protein proximate, volatiles, OAVs, taste compounds) of plant-protein reaction flavours. Confidence 0.7.
E8. Jansson 2020 Int Dairy J 104653 — https://doi.org/10.1016/j.idairyj.2020.104653
  Holds: Table S2 standard curves, Fig. S1 the heating/cooling profile. Why: makes its calibrated volatiles usable. Confidence 0.6.
E9. Xiao 2025 Food Hydrocolloids 111119 — https://doi.org/10.1016/j.foodhyd.2025.111119
  Holds: Tables S1-S4, Fig. S1 (calibrated; content not described in the dossier). Confidence 0.4.
E10. Primary papers cited by the matrix code with no dossier: Barallat-Pérez 2023 (sealed lupin binding values) and
  Hall & Anderson 1983 (hexanal air-water partition, 7-9x from Meynier). Why: two shipped matrix inputs without a source
  on disk. Confidence 0.8 that each gives numbers.
E11. Grosch 1993, the MFT contrast Cerny 2015 cites (primary not on disk). Confidence 0.5.
Below cut (semi-quantitative relative areas, or no lane use): Singh 2021 S1, Wang 2025 S2, Wang 2026c S3, Du 2023 S1,
Zhao 2019b S1, Wang 2020 S1, Liu 2025 S1, Yang 2024 S1-S4, Miyazaki 2023 S1-2, Xu 2024 S1, Crowther 1979 (a thesis).
On disk, no download needed: us10039306.pdf (Impossible Foods) needs a dossier (queued, M3).

## F. One ranking across all lanes (2026-10-29)
Sulfur first (the next sub-models), then integrity of shipped constants, then the other lanes:
1 A3 | 2 A1 | 3 A2 | 4 A5 | 5 A4 | 6 A10 | 7 E1 (shipped pea constants) | 8 E10 (two unsourced matrix inputs)
9 B1 | 10 B4 | 11 A7 | 12 B10 | 13 A8 | 14 E2 | 15 C1 | 16 A9 | 17 E3 | 18 E4 | 19 B2 | 20 B3
21 B7 | 22 B5 | 23 B8 | 24 E6 | 25 E7 | 26 E5 | 27 A14 | 28 A6 | 29 B11 | 30 E8. Then the rest of D and E.
