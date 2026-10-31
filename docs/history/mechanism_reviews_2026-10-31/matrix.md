# Matrix gap: dilute pots vs alt-protein products (2026-10-31, read-only)

**The amine balance in a 16 % pea-isolate product:**
- **Bound lysine:** ≈ 0.47 mmol/g protein × ~130 g/kg ≈ **60 mmol/kg** (50-90).
- **The isolate's free amino acids:** ~1 mmol/kg (Jaeger 2023: "essentially none"; sugars 0.06-0.19 g/100 g).
- **Precursors are ADDED:** e.g. Cys 10 mM + sugar 20 mM (US 9,943,096).

**Per-amine reactivity, same Wageningen lab, pH 6.8:**
- casein-Lys: k3 = 1.5e-4 L/(mmol·min) at 120 C (Nguyen 2016);
- glycine: 7.8e-5 (Martins k_schiff taken to 120 C);
- **ratio ≈ 1.9** (mine). Bound Lys takes ~90 % of the sugar-amine flux. On globular pea or soy protein it should be lower.

## 1. Ranking
P = probability this gap is the largest shift versus a pot with the same free precursors.

| gap | P | size / sign | evidence on disk |
|---|---|---|---|
| (a) bound amines (Lys, Arg, N-termini) | **0.40** | 0.5-1 dex; sign varies by compound | ratio above; Lys/Arg trap the MGO/GO donors (Berk 2021); peptides beat free amino acids on pyrazines (Van Lancker 2012) |
| (c) binding / capture | **0.25** | 0.5-1.5 dex on FREE thiols, ~0 on formation | Xu 2022 (>98 % MFT headspace loss); Mottram 1996 (protein reduces MFT disulfide back to 45 % MFT) |
| (b) low a_w | **0.20** | 3-10x in the crust; splits unmeasured | Pereyra Gonzales 2010 (3.5x, dry basis) vs Bell 1995 (flat at fixed molality); MFT x77 from a_w 0.60 to 0.90 (Modi 2008) vs a maximum at 0.72 (Hartman 1984b) |
| (d) lipid crossover | **0.15** | 0.3-1 dex on S-heterocycles | Farmer, Whitfield, Zamora/Hidalgo |

**Caveats:**
- **Xu 2022 is not a binding constant.** 0.35 mM MFT meets only ~0.21 mM thiol + disulfide sites (mine). So the >98 % includes oxidative dimerisation and possibly an unheated blank.
- **Carried-in odorants:** the isolate powder already holds hexanal (OAV 6202) and 3-methylbutanal (OAV 10186) (Utz 2022). The engine has no input for that initial load.

## 2. What is on disk, and what needs new data

- **(a)** On disk: the reactivity ratio (Nguyen, Brands 2002c, Martins); Lys and Arg densities (Jaeger, Gorissen, Tang 2024); CML and furosine.
  - Missing: globular-protein glycation, pentose on protein, Arg trapping of MGO.
  - **Finding:** the B20 `r_flp_decay` returns lysine (LYSP), but Nguyen's SI returns none (`nguyen2016_extraction.md` ~l.422). The frozen k_glyc is 0.64x the printed k3, which looks like compensation.
- **(b)** Shapes on disk, sizes not: the sources fix different variables and the MFT data conflict.
- **(c)** Aldehydes partly covered (Snel 2023, Wang & Arntfield 2015). No rate for thiol capture by plant protein at process temperature exists anywhere.
- **(d)** Largely on disk.
- **Designed data are unavoidable for product-level numbers (P ≈ 0.85).** Minimal design:
  - isolate 0 / 8 / 16 % × a_w 0.98 / 0.75 / 0.5 (glycerol);
  - fixed Cys, ribose and glucose;
  - SIDA total and headspace-free, plus furosine.

## 3. Minimal model elements

- **(a)** LYSP competes with every sugar in the shared trunk at 1.9 × an accessibility factor of 0.4-1. Make the FLP-decay topology match its fit. Add an ARG dicarbonyl sink with Lo 1994 constants and a declared temperature band.
- **(b)** Mass action on water-phase molality plus a glass factor below Tg. Make water an explicit reactant in hydrolytic β-cleavage (B7), with Hartman and Modi as sign tests. Retire the dry-basis multiplier, which would double-count.
- **(c)** Charge the PSH and PSS pools together (the protein both releases and captures thiols). Give aldehydes Snel's K with a Langmuir cap. Report free and total separately.
- **(d)** Lipid carbonyls act as Strecker donors and thiol traps (Zamora/Hidalgo constants); carried-in volatiles become the initial state.

## 4. Papers not on disk

1. **Lo TW, Westwood ME, McLellan AC, Selwood T, Thornalley PJ.** *J Biol Chem* 1994, 269:32299, doi:10.1016/S0021-9258(18)31635-1.
   - Verified: Crossref record; summary seen in search results only.
   - Reported: MGO reacts with Lys reversibly and with Cys fast; on BSA the irreversible modification is mainly at Arg; 37 C.
2. **Warmbier HC, Schnickels RA, Labuza TP.** *J Food Sci* 1976, 41:528, doi:10.1111/j.1365-2621.1976.tb00663.x.
   - Verified: Crossref record plus UMN abstract.
   - Reported: casein-glucose-glycerol; first-order loss of available lysine; ~1 mol glucose per lysine lost; browning maximum at a_w 0.4-0.5.
3. **Warmbier, Schnickels, Labuza** (glucose:lysine ratio). *J Food Sci* 1976, 41:981, doi:10.1111/j.1365-2621.1976.tb14372.x.
   - Verified: Crossref record only; abstract not seen.
4. **Eichner K, Karel M.** *J Agric Food Chem* 1972, 20:218, doi:10.1021/jf60180a025.
   - Verified: Crossref record only; abstract not seen.
5. **Lu CY, Hao Z, Payne R, Ho CT.** *J Agric Food Chem* 2005, 53:6443, doi:10.1021/jf050534p.
   - Verified: PubMed 16076132.
   - Reported: glucose with Gly, Gly2 or Gly3 at 160 C, 0-100 % water in glycerol; volatiles rise as water falls; Gly ≈ Gly3 > Gly2.
6. **Scalone GLL, Cucu T, De Kimpe N, De Meulenaer B.** *J Agric Food Chem* 2015, 63:5364, doi:10.1021/acs.jafc.5b01129.
   - Verified: PubMed 25971942.
   - Reported: oligopeptides from whey hydrolysate raise pyrazines more than free amino acids; pyrazines highest at a_w 0.33 and above 120 C.
7. **Adams RL, Mottram DS, Parker JK, Brown HM.** *J Agric Food Chem* 2001, 49:4333, doi:10.1021/jf0100797.
   - Verified: PubMed 11559133.
   - Reported: volatile disulfides lost to ovalbumin at pH 8 and much more after denaturation; none at pH 2.
8. **Van Lancker F, Adams A, De Kimpe N.** *J Agric Food Chem* 2010, 58:2470, doi:10.1021/jf903898t.
   - Verified: PubMed 20121201.
   - Reported: Lys-dipeptides form more pyrazines than free amino acids.

**Gaps no search filled:** glycation kinetics on pea or soy; MFT or FFT capture by plant protein at process temperature.
