# The deoxyosone hub: why 3-DX / 1-DX is flat in time and falls with pH

**A correction first.** Prereg sec. 13 says the Zhai 2021 TTCA pots are 0.2 M phosphate. They are unbuffered water (`route_submodel.yml`, `buffer: "none (water)"`); so are Zhai 2019 and 2023. Davidek is in phosphate, so phosphate is still live.

## 1. Hypotheses (P = probability it is a major cause; they are not exclusive)

**H1. Both osones are in quasi-steady state; the engine's 3-DX is not (P 0.65).** If both are short-lived next to their parent (the Cys-ARP), the ratio stays flat at (f3/f1)(k1/k3). The engine's 3-DX lives about 35 h, so it accumulates and the ratio climbs.
- For:
  - Hofmann 1999: 3-DX and 1-DX peak at 10 min, then plateau at a ratio of about 3, measured with authentic quinoxaline standards.
  - Zhai 2023 JAFC Fig. 5: both peak with the ARP at 120 C.
  - Zhai 2021 Fig. 5 at pH 5.5: both plateau from 80 to 120 min.
  - Martins' M4: 1-DG lives about 1 min and 3-DG about 30 min.
- Against: B7T, the 3-DX exit, still failed the time signature (1 / 8). So a second source of 3-DX that grows in time must also be there (H5).

**H2. Cysteine traps the aldehydic 3-DX, and the trap is faster at high pH (P 0.45, mainly the pH 7 to 8 segment).**
- For:
  - Zhai 2019 Fig. 3b: Cys suppresses 3-DX to 0.20 of control but 1-DX only to 0.45.
  - Zhai 2021 Fig. 6: free Cys rises x4.8 from pH 5.5 to 8.
  - **New, from the Zhai 2021 Fig. 5 tables at 120 min.** Extra xylose re-traps the released Cys. It raises 3-DX/1-DX by x0.85 / 0.92 / 1.22 / 2.7 at pH 5.5 / 6 / 7 / 8. So removing Cys rescues 3-DX only at high pH, which matches the curvature of the measured segments.
- Against: the Zhai 2023 JAFC double ratio (0.53 to 2.6) conflicts.

**H3. The enolisation split itself shifts toward 1-DX at high pH (P 0.5 as a contributor).**
- Martins 2005b Table 3: pD 0.35 for k4 (to 3-DG) and 0.60 +- 0.52 for k7 (to 1-DG). That is -0.25 dex per pH unit in the split, against -0.49 measured.
- Nedvidek 1992 (trap in situ, so formation flux): 1-DG/3-DG rises x5 from pH 5 to 7.

**H4. The 1-DX sinks were calibrated in phosphate and applied to water (P 0.35).**
- `r_dpo_aa` is fitted on Davidek 2008 (phosphate).
- Davidek 2006b: beta-cleavage "favored in the presence of phosphate".
- Potman 1989: phosphate speeds Amadori formation and breakdown.
- Nothing on disk runs water against phosphate for a 1-deoxyosone.

**H5. The engine's pH-flat 3-DX source from xylose caramelisation is too large (P 0.4).** It supplies 39 to 45 % of the engine's 3-DX. Against it: at pH 5.5, adding 10 mM free xylose raises 3-DX only x1.24 but 1-DX x1.47 (Zhai 2021).

**H6. Oxygen diverts the ARP to the osone and away from 1-DX (P 0.15).** Gobert 2009: 1-DG 23.5 without air against 6.9 with air; glucosone 5.5 against 44.5. The Zhai headspaces are not stated.

**H7. OPD artefact (12 h post-hoc derivatisation) or pH drift in the unbuffered pots (P 0.1 each).**

## 2. Discriminating observations (on disk?)

- **H1:** ARP and both osones on one time grid. Partly on disk (Zhai 2023 JAFC, pH 7, coarse). Fed pentosone decay in water: no.
- **H2:** the +Xyl/-Xyl double ratio against pH. **Yes** (Zhai 2021); declare it as a signature. Fed 3-DX + Cys against pH: no.
- **H3:** in situ trapping on a pentose ARP across pH. Hexose only (Nedvidek, Martins 2005b).
- **H4:** 1-DX decay in water against phosphate. No.
- **H5:** free xylose at low pH. **Yes** (Zhai 2021).

## 3. Multiresponse nulls

- **Martins & van Boekel 2005, M4** (glucose + glycine, 0.1 M phosphate; k at 100 C in /min, Ea in kJ/mol):

  | step | k | Ea |
  |---|---|---|
  | k4: DFG -> 3-DG | 1.1e-2 | 97.1 |
  | k7: DFG -> 1-DG | 1.6e-2 | 107.3 |
  | k6: DFG -> MGO | 7.1e-3 | 124.5 |
  | k8: 1-DG -> acetic acid | 1.4 | 75.7 |
  | k5: 3-DG -> formic acid | 3.5e-2 | 29.6 |
  | k9: 3-DG + Gly | 8.1e-4 | 95.2 |

  - Use its pD values (2005b) only as bands: k5 and k8 are not monotone in pH.
  - Its topology already puts both osones in quasi-steady state.
- Other models to carry as nulls:
  - Martins 2003 Part II (DFG, 16 steps).
  - Brands & van Boekel 2002c (sugar + casein, 90-130 C, osones lumped).
  - Kocadagli & Gokmen 2016.
  - Knol 2005; Knol 2010 (Ea 84 for 3-DG -> formic acid).
  - Yu 2018 / 2020.
- Davidek, Hofmann, Weenen and Yaylayan (CAMOLA) give pathways, not models: use as branching priors.
- I found no published pentose or Cys-ARP multiresponse model.

## 4. Papers not on disk (DOIs checked against PubMed or Crossref records; content as the abstracts state it)

1. Wei, Cui, Hayat, Zhang, Ho. *JAFC* 2022, 70, 12164. doi:10.1021/acs.jafc.2c04961 (PubMed 36124743). Reported: Cys reacts with deoxyosones from the Gly-Xyl ARP to form TTCA, diverting them from short dicarbonyls (H2).
2. Zhou, Huang, ..., Ho. *JAFC* 2024, 72, 5878. doi:10.1021/acs.jafc.4c00021 (PubMed). Reported: alanine speeds an ARP + 3-deoxyxylosone reaction, which lowers furfural. That is a 3-DX sink proportional to the ARP (H1, H2).
3. Cui, Ma, ..., Ho. *JAFC* 2022, 70, 8409. doi:10.1021/acs.jafc.2c03427 (PubMed). Reported: glycylglycine catalyses ARP -> deoxyosone at 120 C (amine factor, H3).
4. Smuda, Voigt, Glomb. *JAFC* 2010, 58, 6458. doi:10.1021/jf100334r (PubMed 20429584). Reported: 1-DG + Boc-lysine gives acetyl, formyl, lactoyl and glycerinyl amides by beta-cleavage (amine-induced cleavage, H4).
5. Voigt, Smuda, Pfahler, Glomb. *JAFC* 2010, 58, 5685. doi:10.1021/jf100140h (PubMed 20441226). Reported: oxygen-dependent lactic and glyceric acid from 1-DG + alanine (H6).
6. Hou, Zhou, Li, ..., Zhang, Ho. *Food Chem* 2026, 499, 147363. doi:10.1016/j.foodchem.2025.147363 (PubMed 41344244). Reported: Cys oxidised to cystine speeds ARP degradation and the retro-aldol of deoxyosones to GO / MGO (H6).
7. Zyzak, Richardson, Thorpe, Baynes. *Arch Biochem Biophys* 1995, 316, 547. doi:10.1006/abbi.1995.1073 (PubMed 7840665). Reported: Amadori decomposition is faster in phosphate than HEPES, faster under air and at higher pH, and gives 3-DG (37 C; H4).
8. Thornalley, Langborg, Minhas. *Biochem J* 1999, 344, 109. doi:10.1042/bj3440109 (Crossref). Reported: GO / MGO / 3-DG formation from glucose and fructosyl-lysine depends on phosphate concentration and trace metals (37 C; H4, H6).

**Unverified content:** Huyghues-Despointes & Yaylayan 1996 (10.1021/jf9502921) and Beck, Ledl & Severin 1988 (10.1016/0008-6215(88)85058-4). The DOIs are on Crossref, but I could not read an abstract.

**Proposed next arm (no fit run):**
- Martins-like quasi-steady-state lifetimes in water: a 3-DX exit, and the 1-DX sink without the phosphate factor.
- Cut the caramelisation share.
- Add a Cys-3-DX trap with a free-amine factor.
- Declare the +Xyl double ratio beforehand as a held-out signature.
