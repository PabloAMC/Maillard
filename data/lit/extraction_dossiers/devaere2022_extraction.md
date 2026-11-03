# Devaere, De Winne, Dewulf, Fraeye, Šoljić, Lauwers, de Jong & Sanctorum 2022 — EXTRACTION (file stem `devaere2022`)
DOI: 10.3390/foods11131985 (printed p. 1). Foods **11**:1985 (Communication; received 5 Jun, accepted 1 Jul,
published 5 Jul 2022). Crossref 2026-11-03: "Improving the Aromatic Profile of Plant-Based Meat Alternatives: Effect
of Myoglobin Addition on Volatiles", Foods 11, article 1985, the same eight authors in the same order, issued
2022-07-05. **Matches the PDF and the brief.** No naming caveat.

### Soy-textured-protein burgers (25 % TSP, 15 % sunflower oil, 1 % methylcellulose, 1.5 % salt, 57.5 % water) with 0 / 0.5 / 1.0 % (m/m, on top) commercial bovine myoglobin, raw or oven-cooked 12 min at 250 °C; HS-SPME (DVB/CAR/PDMS, 45 °C) GC-MS, 40 volatiles as **raw peak areas (AU × 10³), no internal standard**, n = 3. **No sulfur compound detected, no sugar / cysteine / yeast extract in the recipe, no sensory, no core temperature, no kinetics.**

**Read 2026-11-03 from `data/articles/devaere2022.pdf`** (12 pp., MDPI, born-digital LaTeX). Prose from
`pdftotext -layout`; **Table 1 (pp. 5-7, landscape) read by eye from 200 dpi renders** (scratchpad
`ingest/devaere2022_work/p-05.png`, `p-06.png`, `p-07.png`); the text layer agreed with the renders cell for cell
apart from displaced letters on split cells (hexanal, 2-heptenal, 2-pentylfuran), which were taken from the renders.
Fig. 1 (PCA, p. 9) viewed at 110 dpi for content only; its printed variance shares are in the text. Class sums,
total volatiles and the text's area percentages re-computed with `bc` / `awk` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** KU Leuven Technology Campus Ghent (Devaere, De Winne, Dewulf, **Ilse Fraeye**: senior academic
author, supervision) with **Paleo b.v.** (Gosselies / Diest, Belgium; a fermentation-derived-myoglobin company):
corresponding author **Hermes Sanctorum** (Paleo consultant), de Jong (consultant), Šoljić and Lauwers (employees).
Funding VLAIO HBC.2021.0263. Data are part of patent application **PCT/EP2021/087884** [M]. **Not in
`data/keys/laboratories.yml`: NEW**, proposed id **`ku_leuven_ghent_fraeye`**, lab string
"KU Leuven Ghent (Fraeye) / Paleo (Sanctorum)". Same campus as `bustillotrueba2018_extraction.md` but a different
group (EFBT brewing, De Cooman): **not the same laboratory**. No other dossier by this group. Related heme dossiers:
`us9943096_extraction.md`, `ep2943072b2_extraction.md`, `us10039306_extraction.md` (Impossible Foods, LegH +
cysteine + sugar in buffer), `thong2024_extraction.md` (commercial burgers incl. leghemoglobin, sensory).
[!] **Industry conflict of interest**: four of eight authors are paid by the company that sells the ingredient.

## 1. What this paper answers for the model

| model need | what Devaere 2022 gives | tag |
|---|---|---|
| heme effect on lipid-oxidation volatiles in a plant matrix | Table 1: 3 Mb doses × raw / cooked, 11 aldehydes, 1-octen-3-ol, 2,3-octanedione, 2-pentylfuran | [M], peak area only |
| heme effect on sulfur / "meaty" volatiles | **none detected** (no S compound among 40); no S precursors added | [M] absence |
| heme effect on Maillard volatiles | 6 alkylpyrazines, 2-/3-methylbutanal, furfural, furfuryl alcohol, pyrrole, maltol | [M] |
| absolute concentrations (µg/kg) | **none**: no internal standard, no calibration, no LOD | [!] |
| time / temperature series, kinetics | **none**: one cooking condition (oven air 250 °C, 12 min); no core T, no weight loss | [!] |
| oxidation state of the lipid (PV, TBARS, O2) | **none** | [!] |
| pH, a_w, moisture after cooking | **none** | [!] |
| sensory | **none** (no panel; the "desirable" conclusion rests on GC-MS + PCA only) | [Q] |
| heme protein identity / purity / iron content | "commercial bovine Mb (Tebu-bio)" only; purity, form, met/oxy state not stated | [!] |

**Answers** (paper-specific questions).
1. **Recipe** [M]: RO water 57.5 % (45 % to hydrate TSP 30 min at 15 °C, 12.5 % to dissolve Mb), textured soy
   protein 25 % (Fuji Oil, Ghent; protein content, extrusion history and residual lipid not stated), **sunflower oil
   15 %** (Vandemoortele; regular vs high-oleic **not stated** [!], which sets the linoleate share the lipid lane
   needs), methylcellulose 50D 1 %, table salt 1.5 %, plus **bovine myoglobin 0 / 0.5 / 1.0 % m/m "on top"**
   (actual 0.498 / 0.990 % of the batter [Z]). **No sugar, no cysteine, no thiamine, no yeast extract, no
   flavouring** (confidence ~97 % from the methods list). Heme dose if the Mb were pure (17.6 kDa, 1 Fe): 1.0 % →
   0.57 mmol/kg batter, ~3.2 mg Fe/100 g, ~1.0 mM on the water phase; 0.5 % half of that [Z, upper bounds].
   Authors: 0.5-1.0 % mimics endogenous meat myoglobin [Q].
2. **Cooking** [M]: 30 g patties, Ø 6 cm (≈ 1.1 cm thick at density 1 [Z]), "baked for 12 min at 250 °C in an
   oven" (Rational Climaplus Combi CPC 61; convection / steam mode not stated). [!] **Abstract, table headings and
   discussion call this "grilled"**; the methods describe a combi oven. No core T, no surface T, no cooking loss, no
   turning, no resting. Raw samples were never heated except in the SPME vial (below).
3. **Volatiles and units** [M]: HS-SPME, 3 g in a 20 mL crimped vial, stored at 4 °C in the tray (duration not
   stated), 20 min equilibration + **40 min extraction at 45 °C** on 50/30 µm DVB/CAR/PDMS; GC-MS HP-1ms, full scan
   m/z 40-250; tentative identification by NIST08/Wiley275 + LRI. **Units: absolute peak area, AU × 10³, mean ± SD,
   n = 3**; no internal standard, no response factors, so only within-compound ratios between arms are meaningful.
4. **Sensory**: none. The conclusion that Mb improves the profile "in a desirable way" [V] is inferred from the shift
   towards aldehydes and pyrazines, not from a panel; the authors themselves say whether hexanal exceeds rancid
   thresholds needs "further research" [Q].
5. **Which volatiles move with myoglobin, how much, dose-dependent?** (fold changes [Z] from Table 1, sec. 4)
   - **Linoleate-type lipid oxidation products rise strongly, already in the raw batter** (confidence ~90 % that the
     sign is real): hexanal raw ×3.65 / ×4.53 (0.5 / 1.0 % Mb), cooked ×9.0 / ×11.9; 1-octen-3-ol raw ×3.7 / ×5.1,
     cooked ×7.1 / ×10.3; (E,E)-2,4-decadienal raw ×3.6 / ×7.2, cooked ×5.5 / ×5.9; 2,3-octanedione raw ×6.1 / ×8.6,
     cooked ×11.7 / ×18.1; 2-heptenal from n.d. to 911-1985; unsaturated-aldehyde class raw ×26 / ×47, cooked ×22 /
     ×38.
   - **Oleate-type products rise weakly**: nonanal raw ×1.04 / ×1.60, cooked ×1.23 / ×1.89; 2-decenal and
     2-undecenal from trace to tens; decanal **falls** to n.d.
   - **Dose response**: hexanal saturates (0 → 0.5 % gives 75 % raw / 73 % cooked of the 0 → 1.0 % increment);
     unsaturated aldehydes are closer to linear (56 % / 55 %) [Z]. Whether the hexanal saturation is chemistry or
     SPME fibre saturation (hexanal is the largest peak, 10⁴ AU × 10³) cannot be told; P(chemistry) ~40 %.
   - **Things that fall with Mb**: 2-pentylfuran raw ×0.67 / ×0.64 (cooked ×0.96 / ×0.93), 1-hexanol (raw ×0.66,
     cooked ×0.35), hexane and the other alkanes (×0.29-0.48), Strecker aldehydes in cooked samples (2- + 3-methyl-
     butanal ×0.38 / ×0.39). Hexane is a solvent residue of the TSP [Q], so a 50 % drop is hard to explain by
     chemistry; **competitive displacement on the Carboxen fibre** as the aldehyde load rises is the likelier cause
     (~60 %) and would bias every "decrease" in this table [!].
   - **Pyrazines**: unchanged by Mb in raw samples (carried in from the TSP [Q]); in cooked samples ×1.33 / ×2.24
     (total), with Tukey letters b / b / c (only 1.0 % differs from 0 %). Furfural cooked 16 → 106 → 148. Authors
     suggest Mb enhances Maillard browning [Q]; no browning was measured, and the falling Strecker aldehydes argue
     against a simple "more Maillard" reading [Z].
   - **Sulfur compounds: none detected in any arm** (confidence ~99 % that the paper has no sulfur data). The recipe
     has no added S precursor beyond soy protein methionine / cysteine; full-scan SPME without a sulfur-selective
     detector would miss ppt-level MFT / methanethiol anyway. **The claim that heme promotes meaty sulfur volatiles
     is not tested here**, neither confirmed nor refuted.
6. **Cooking × heme interaction** [Z]: without Mb, oven cooking lowers hexanal ×0.42 (2834 → 1201) and saturated
   aldehydes ×0.49; with Mb, cooked / raw is ×1.04 (0.5 %) and ×1.11 (1.0 %) (interaction p = 0.002 [M]). If the
   no-Mb fractional loss (58 %) also applied to the Mb arms, heme-dependent formation *during* cooking would be
   ~6.5 × 10³ AU × 10³ hexanal at 0.5 % and ~8.9 × 10³ at 1.0 % [Z, assumes equal loss]: comparable to the raw-batter
   effect. So heme acts both at ≤ 45 °C (batter, tray, vial) and during heating.
7. **Usable as a directional check of a heme-catalysis term the model lacks? Yes, for sign and rough fold only**
   (P(useful as sign check) ~85 %; P(useful for magnitude) ~10 %). The lipid lane (`src/kinetic_core/lipid.py`)
   decomposes a **declared** hydroperoxide pool (carrier + peroxide value) first-order with a Q10 assumption; it has
   **no initiation / propagation chemistry and no heme, iron or metal-catalysis term** (grep of `src/kinetic_core/`
   for heme / myoglobin / iron / catal: no hit; `data/lit/measurement_ledger/lipid.yml` has heme only in the
   Impossible Foods patent entries). Today the engine would predict **identical** lipid volatiles in all three Mb arms
   and **no** raw-batter change at all. Checks a future heme term should pass: (a) d[hexanal, 1-octen-3-ol,
   2-alkenals, 2,4-decadienals]/d[Mb] > 0 at 0-1 % Mb; (b) the effect exists without cooking (≤ 45 °C, ≤ ~1 h plus an
   unstated 4 °C hold), i.e. heme must create new hydroperoxide (initiation / propagation) or decompose it at low T,
   not only speed the hot step; (c) linoleate products respond far more than oleate products; (d) the heme-induced
   rise in hexanal is not accompanied by a rise in 2-pentylfuran (2-pentylfuran / hexanal raw 2.82 → 0.52 → 0.40,
   cooked 6.42 → 0.68 → 0.50 [Z]), which bears on any fixed 2-pentylfuran-per-hexanal constant. No absolute
   numbers can be scored (no IS, unknown vial time, unknown PV and oil type).

## 2. Conditions

- **Batter** [M]: as sec. 1 ans. 1. Grinding through a 4 mm plate (Minerva mincer) before Mb addition; Mb dissolved
  in the 12.5 % water and mixed in after grinding. Whether the 0 % arm received the plain 12.5 % water is not
  stated (presumed [Z]).
- **Heme protein** [M]: commercial bovine myoglobin, Tebu-bio (a distributor); purity, lot, oxidation state
  (met- / oxy-Mb), co-purified lipid and iron content **not stated** [!]. Not the fermentation-derived Mb the paper
  motivates.
- **pH / a_w / moisture**: none measured. **Buffer**: none (food matrix).
- **Heating** [M]: combi oven, 250 °C air, 12 min, 30 g patties; core T, mode, weight loss **not reported** [!].
- **Headspace / O2**: air during mixing and in the SPME vial (20 mL, 3 g sample, 45 °C, 60 min) [M/Z].
- **Replicates** [M]: "HS-SPME-GC-MS analysis performed in triplicate" [V]; whether three patties or three vials
  from one batch is not stated (P(analytical replicates of one batch) ~60 %) [!]. One batch per arm presumably.
- **Statistics** [M]: two-way ANOVA (%Mb, G, %Mb × G, SPSS 27); when the interaction is significant, one-way ANOVA +
  Tukey HSD over all six arms (letters); p < 0.05; PCA (Unscrambler X 10.5.1, mean-centred, unit variance).

## 3. Analytical methods

- **Extraction** [M]: Gerstel MPS2, 50/30 µm DVB/CAR/PDMS (Supelco), fibre conditioned 30 min at 270 °C;
  3 g ± 0.1 % in 20 mL vial, PTFE/silicone septum; 20 min at 45 °C (250 rpm, 5 s on / 2 s off), 40 min exposure at
  45 °C; 5 min desorption at 250 °C splitless (0.75 mm liner).
- **GC-MS** [M]: Agilent 7890A / 5975C, HP-1ms 30 m × 0.25 mm × 0.25 µm, He 1.2 mL/min; 35 °C (5 min) → 4 °C/min →
  215 °C → 7 °C/min → 250 °C (5 min); EI 70 eV, source 230 °C, quad 150 °C, scan m/z 40-250.
- **Identification** [M]: NIST08 / Wiley275 match + LRI (van den Dool & Kratz, C6-C16 alkanes) vs NIST WebBook
  values; **tentative** (no authentic standards stated). LRI gaps of 15-31 units for 2-heptenal (928 vs 951),
  2-octenal (1030 vs 1061), 2-decenal (1233 vs 1256), 2-undecenal (1334 vs 1350), decanal (1183 vs 1203) [M] [!].
- **Quantification** [M]: total-ion peak area, AU × 10³; **± = SD of n = 3**. No IS, no calibration, no LOD ("n.d.").
- Raw-ingredient (TSP) SPME runs were done but are "data not shown" [Q].

## 4. Data

### 4a. Table 1 (pp. 5-7): peak area, AU × 10³, mean ± SD, n = 3 [M]. Letters = Tukey groups (only where %Mb × G was significant). ANOVA p as printed (%Mb / G / %Mb × G; "-" = not analysed; NS).

| No. | compound | RMA | RMA+0.5%Mb | RMA+1.0%Mb | GMA | GMA+0.5%Mb | GMA+1.0%Mb | p %Mb / G / int |
|---|---|---|---|---|---|---|---|---|
| | **Acids** | 204.67 ± 37.24 | 163.53 ± 34.15 | 304.27 ± 50.84 | 235.46 ± 160.10 | 1093.37 ± 1199.55 | 568.79 ± 377.98 | NS/NS/NS |
| 1 | acetic acid | 179.22 ± 33.03 | 118.01 ± 49.49 | 208.00 ± 46.91 | 212.77 ± 172.16 | 309.64 ± 101.23 | 263.82 ± 173.09 | NS/NS/NS |
| 2 | octanoic acid | 5.52 ± 4.95 | 17.15 ± 8.81 | 34.27 ± 21.24 | 8.47 ± 7.44 | 101.15 ± 149.02 | 42.01 ± 17.96 | NS/NS/NS |
| 3 | nonanoic acid | 19.93 ± 8.77 | 28.36 ± 10.44 | 62.00 ± 50.16 | 14.23 ± 12.45 | 682.58 ± 1015.61 | 262.95 ± 214.72 | NS/NS/NS |
| | **Alcohols** | 3685.55 ± 524.65 | 4684.63 ± 591.67 | 5586.83 ± 163.30 | 1545.67 ± 228.85 | 3364.99 ± 174.73 | 4244.33 ± 318.35 | 0.000/0.000/NS |
| 4 | furfuryl alcohol | n.d. a | n.d. a | n.d. a | 68.35 ± 36.08 ab | 89.64 ± 35.55 bc | 155.34 ± 36.25 c | -/-/0.030 |
| 5 | 1-hexanol | 2976.13 ± 416.67 | 2089.81 ± 442.57 | 1950.16 ± 274.80 | 1116.23 ± 209.27 | 727.07 ± 96.73 | 386.26 ± 46.11 | 0.001/0.000/NS |
| 6 | 1-octen-3-ol | 709.42 ± 109.36 | 2594.82 ± 310.85 | 3636.68 ± 399.18 | 361.09 ± 24.91 | 2548.28 ± 164.25 | 3702.73 ± 301.48 | 0.000/NS/NS |
| | **Branched aldehydes** | 5.15 ± 8.91 a | 14.28 ± 24.73 a | 9.76 ± 16.90 a | 506.23 ± 21.30 c | 192.92 ± 48.22 b | 199.87 ± 41.69 b | -/-/0.000 |
| 7 | 3-methylbutanal | 5.15 ± 8.91 a | 9.66 ± 16.73 a | 5.73 ± 9.93 a | 143.27 ± 15.05 c | 66.47 ± 10.18 b | 63.01 ± 8.83 b | -/-/0.000 |
| 8 | 2-methylbutanal | n.d. a | 4.62 ± 8.00 a | 4.03 ± 6.98 a | 362.96 ± 9.83 c | 126.45 ± 38.42 b | 136.86 ± 34.70 b | -/-/0.000 |
| | **Cyclic aldehydes** | 1146.19 ± 107.85 | 1742.39 ± 299.87 | 2177.20 ± 300.31 | 699.86 ± 105.60 | 1167.18 ± 107.89 | 1368.04 ± 25.49 | 0.000/0.000/NS |
| 9 | furfural | n.d. a | n.d. a | n.d. a | 16.15 ± 27.97 a | 106.27 ± 52.19 b | 148.46 ± 35.45 b | -/-/0.005 |
| 10 | benzaldehyde | 1146.19 ± 107.85 | 1742.39 ± 299.87 | 2177.20 ± 300.31 | 683.71 ± 78.49 | 1060.91 ± 68.71 | 1219.58 ± 60.78 | 0.000/0.000/NS |
| | **Saturated aldehydes** | 3193.57 ± 343.90 a | 10,662.48 ± 1011.73 b | 13,330.93 ± 481.47 c | 1556.37 ± 173.13 a | 11,275.58 ± 921.93 b | 14,951.23 ± 41.23 c | -/-/0.002 |
| 11 | hexanal | 2834.01 ± 321.13 b | 10,340.95 ± 916.87 c | 12,848.83 ± 420.92 d | 1201.11 ± 176.28 a | 10,800.63 ± 911.84 c | 14,297.20 ± 11.47 d | -/-/0.002 |
| 12 | heptanal | 151.18 ± 16.64 ab | 114.66 ± 63.14 a | 165.32 ± 14.29 ab | 124.27 ± 7.53 a | 213.31 ± 7.45 bc | 253.48 ± 9.39 c | -/-/0.004 |
| 13 | nonanal | 198.36 ± 24.19 | 206.87 ± 44.04 | 316.78 ± 63.56 | 212.23 ± 8.77 | 261.64 ± 23.14 | 400.55 ± 26.92 | 0.000/0.012/NS |
| 14 | decanal | 10.02 ± 9.15 | n.d. | n.d. | 18.77 ± 2.25 | n.d. | n.d. | 0.000/NS/NS |
| | **Unsaturated aldehydes** | 64.91 ± 11.69 a | 1713.19 ± 333.87 b | 3028.08 ± 716.24 c | 111.56 ± 3.34 a | 2412.52 ± 257.87 bc | 4291.35 ± 375.27 d | -/-/0.046 |
| 15 | 2-heptenal | n.d. a | 910.59 ± 130.13 b | 1482.92 ± 264.78 c | n.d. a | 1135.85 ± 136.92 bc | 1984.58 ± 147.95 d | -/-/0.036 |
| 16 | 2-octenal | 34.27 ± 1.89 a | 310.70 ± 64.42 b | 563.49 ± 167.38 c | 29.13 ± 8.54 a | 396.22 ± 62.10 bc | 908.21 ± 79.95 d | -/-/0.010 |
| 17 | (E,E)-2,4-nonadienal | 8.68 ± 7.53 | 373.74 ± 113.42 | 719.18 ± 210.16 | n.d. | 338.64 ± 40.71 | 735.01 ± 87.60 | 0.000/NS/NS |
| 18 | 2-decenal | 0.85 ± 1.47 | 23.09 ± 8.46 | 55.54 ± 19.43 | 1.31 ± 2.27 | 30.08 ± 3.60 | 74.49 ± 8.48 | 0.000/NS/NS |
| 19 | (E,Z)-2,4-decadienal | n.d. a | 19.97 ± 7.63 ab | 48.12 ± 15.24 b | n.d. ab | 56.37 ± 4.00 c | 81.44 ± 9.77 c | -/-/0.004 |
| 20 | (E,E)-2,4-decadienal | 18.65 ± 4.85 a | 66.21 ± 26.82 ab | 133.86 ± 34.23 b | 79.56 ± 2.31 ab | 441.46 ± 44.42 c | 472.97 ± 57.88 c | -/-/0.000 |
| 21 | 2-undecenal | 2.47 ± 2.61 | 8.89 ± 4.20 | 24.96 ± 10.36 | 1.56 ± 1.37 | 13.90 ± 2.04 | 34.65 ± 3.66 | 0.000/NS/NS |
| | **Alkanes** | 9537.27 ± 3932.97 | 4757.54 ± 1085.68 | 4946.03 ± 605.80 | 9807.61 ± 521.13 | 4839.94 ± 633.10 | 3430.74 ± 447.39 | 0.000/NS/NS |
| 22 | 2-methylpentane | 577.36 ± 431.74 | 262.44 ± 111.86 | 290.99 ± 30.33 | 602.01 ± 74.63 | 321.95 ± 91.95 | 232.19 ± 28.72 | 0.020/NS/NS |
| 23 | 3-methylpentane | 1091.26 ± 558.45 | 521.94 ± 142.57 | 562.94 ± 91.92 | 1107.55 ± 86.87 | 531.33 ± 112.32 | 375.19 ± 63.76 | 0.001/NS/NS |
| 24 | hexane | 5502.04 ± 1770.46 | 2556.70 ± 468.23 | 2619.76 ± 338.67 | 6235.53 ± 306.44 | 2686.46 ± 304.45 | 1835.29 ± 209.84 | 0.000/NS/NS |
| 25 | methylcyclopentane | 1125.23 ± 673.07 | 687.57 ± 146.88 | 658.26 ± 54.99 | 807.21 ± 93.13 | 443.10 ± 75.67 | 338.51 ± 78.14 | 0.032/NS/NS |
| 26 | cyclohexane | 1159.83 ± 485.84 | 693.62 ± 196.39 | 760.50 ± 92.52 | 712.66 ± 63.45 | 395.53 ± 44.83 | 322.60 ± 47.96 | 0.013/0.003/NS |
| 27 | octane | 81.56 ± 20.14 a | 35.28 ± 30.63 a | 53.59 ± 6.93 a | 342.65 ± 6.78 b | 461.57 ± 23.92 c | 326.96 ± 30.81 b | -/-/0.000 |
| | **Ketones** | 625.06 ± 72.65 | 715.46 ± 134.18 | 796.79 ± 18.06 | 2708.99 ± 53.96 | 2488.64 ± 514.06 | 3340.68 ± 361.13 | 0.021/0.000/NS |
| 28 | 2-heptanone | 582.86 ± 66.40 | 458.72 ± 83.76 | 434.12 ± 59.12 | 2674.78 ± 53.78 | 2088.24 ± 518.36 | 2723.19 ± 399.77 | NS/0.000/NS |
| 29 | 2,3-octanedione | 42.21 ± 6.30 a | 256.74 ± 53.32 b | 362.67 ± 41.37 c | 34.21 ± 0.97 a | 400.41 ± 22.52 c | 617.49 ± 39.43 d | -/-/0.000 |
| | **Phenols** | 12.62 ± 3.70 | n.d. | 0.19 ± 0.33 | 185.83 ± 60.22 | 169.72 ± 31.93 | 310.78 ± 105.83 | NS/0.000/NS |
| 30 | guaiacol | 0.57 ± 0.99 a | n.d. a | 0.19 ± 0.33 a | 52.99 ± 32.95 ab | 94.46 ± 29.05 b | 205.96 ± 46.20 c | -/-/0.001 |
| 31 | 2-methoxy-4-vinylphenol | 12.05 ± 2.83 | n.d. | n.d. | 132.85 ± 50.37 | 75.26 ± 11.17 | 104.82 ± 61.59 | NS/0.000/NS |
| | **Pyrazines** | 300.22 ± 28.26 a | 270.90 ± 60.00 a | 285.85 ± 34.30 a | 3557.19 ± 948.44 b | 4717.95 ± 1057.69 b | 7969.41 ± 1309.12 c | -/-/0.001 |
| 32 | methylpyrazine | 21.31 ± 20.28 a | 51.44 ± 20.80 ab | 36.19 ± 31.53 a | 495.34 ± 123.11 bc | 802.87 ± 226.49 c | 1266.36 ± 305.34 d | -/-/0.006 |
| 33 | 2,5-dimethylpyrazine | 101.22 ± 4.32 a | 79.74 ± 17.63 a | 91.92 ± 2.71 a | 1264.96 ± 366.65 b | 1649.87 ± 388.37 b | 2938.08 ± 510.88 c | -/-/0.001 |
| 34 | 2-ethyl-6-methylpyrazine | 24.39 ± 4.90 a | 20.31 ± 2.23 a | 23.61 ± 1.16 a | 321.32 ± 100.02 b | 410.41 ± 112.05 b | 710.70 ± 145.65 c | -/-/0.005 |
| 35 | 2-ethyl-5-methylpyrazine | 57.26 ± 10.19 a | 41.63 ± 7.28 a | 51.90 ± 7.34 a | 647.55 ± 134.84 b | 850.09 ± 167.10 b | 1390.69 ± 158.40 c | -/-/0.000 |
| 36 | 2,5-dimethyl-3-ethylpyrazine | 89.18 ± 11.76 a | 77.79 ± 13.32 a | 79.78 ± 2.08 a | 701.87 ± 196.78 b | 853.54 ± 179.08 b | 1392.55 ± 198.91 c | -/-/0.002 |
| 37 | 2-methyl-3,5-diethylpyrazine | 6.87 ± 4.10 a | n.d. a | 2.46 ± 4.26 a | 126.15 ± 28.06 b | 151.17 ± 27.13 b | 271.03 ± 57.79 c | -/-/0.002 |
| | **Others** | 8461.00 ± 1206.46 b | 5849.08 ± 417.59 a | 5719.21 ± 224.49 a | 8937.64 ± 28.79 b | 8661.31 ± 152.26 b | 8847.58 ± 969.56 b | -/-/0.009 |
| 38 | pyrrole | n.d. | n.d. | n.d. | 280.33 ± 67.51 | 232.82 ± 117.28 | 270.24 ± 71.31 | NS/0.000/NS |
| 39 | 2-pentylfuran | 8001.63 ± 1252.79 b | 5354.32 ± 177.22 a | 5136.56 ± 258.42 a | 7706.66 ± 314.55 b | 7365.33 ± 194.02 b | 7196.63 ± 798.75 b | -/-/0.011 |
| 40 | maltol | 459.38 ± 88.24 | 494.76 ± 242.44 | 582.65 ± 118.56 | 950.65 ± 222.92 | 1063.15 ± 22.85 | 1380.71 ± 107.83 | 0.024/0.000/NS |
| | **Total volatiles** | 27,236.22 ± 6183.06 a | 30,573.47 ± 3231.45 ab | 36,185.15 ± 2222.29 bc | 29,852.42 ± 480.67 ab | 40,384.13 ± 1512.33 c | 49,522.78 ± 1936.45 d | -/-/0.035 |

RMA = raw, GMA = "grilled" (oven-cooked). LRI (unknown / reference) printed for each compound; not transcribed
except where flagged (sec. 3).

### 4b. Text percentages of total peak area [M] and my re-computation from Table 1 means [Z]

RMA: hydrocarbons 34.0 ± 6.7 % (35.0 [Z]), 2-pentylfuran 29.7 ± 2.2 % (29.4), alcohols 13.7 % (13.5),
unsaturated aldehydes 12.0 % [!] (0.24 [Z]: 64.91 / 27,236; the 12.0 % matches no class — saturated aldehydes are
11.7 %), hexanal 10.6 ± 1.3 % (10.4). Saturated aldehydes RMA+1.0 % 36.9 ± 1.0 (36.8), GMA+1.0 % 30.2 ± 1.3
(30.2). Hydrocarbons GMA 32.9 ± 2.1 (32.9), RMA+1.0 % 13.7 ± 1.1 (13.7), GMA+1.0 % 6.9 ± 0.7 (6.9). Pyrazines GMA
11.9 ± 3.0 (11.9), GMA+1.0 % 16.1 ± 2.3 (16.1). PCA: PC1 52 %, PC2 28 %, together 79.8 % [M].

### 4c. Derived Mb effects [Z] (ratio of means, 0.5 % / 1.0 % vs 0 % Mb, same cooking state)

| compound / class | raw | cooked |
|---|---|---|
| hexanal | 3.65 / 4.53 | 8.99 / 11.90 |
| saturated aldehydes | 3.34 / 4.17 | 7.24 / 9.61 |
| unsaturated aldehydes | 26.4 / 46.7 | 21.6 / 38.5 |
| (E,E)-2,4-decadienal | 3.55 / 7.18 | 5.55 / 5.94 |
| 1-octen-3-ol | 3.66 / 5.13 | 7.06 / 10.25 |
| 2,3-octanedione | 6.08 / 8.59 | 11.7 / 18.1 |
| nonanal | 1.04 / 1.60 | 1.23 / 1.89 |
| 2-pentylfuran | 0.67 / 0.64 | 0.96 / 0.93 |
| hexane | 0.46 / 0.48 | 0.43 / 0.29 |
| pyrazines (total) | 0.90 / 0.95 | 1.33 / 2.24 |
| 2-+3-methylbutanal | (trace) | 0.38 / 0.39 |
| total volatiles | 1.12 / 1.33 | 1.35 / 1.66 |

Cooked / raw at 0 / 0.5 / 1.0 % Mb: hexanal 0.42 / 1.04 / 1.11; saturated aldehydes 0.49 / 1.06 / 1.12;
pyrazines 11.8 / 17.4 / 27.9 [Z].

### 4d. Figures

Fig. 1 (PCA scores and correlation loadings): no model-relevant numbers beyond the printed variance shares; not
digitised.

## 5. Flags and internal consistency checks

1. **Arithmetic checks pass** [Z]: every class total equals the sum of its compounds (within 0.01), and the 11 class
   totals sum to the printed total volatiles in all six arms (within 0.02).
2. [!] **"Grilled" vs "baked in an oven"**: the cooking was a combi oven at 250 °C; every heading and the abstract
   say grilled. Heat-transfer mode matters for any thermal simulation.
3. [!] **No internal standard; peak areas on a competitive DVB/CAR/PDMS fibre.** Total volatile load rises ×1.66;
   hexane (an inert solvent residue) falls by half with Mb in raw and cooked samples. Every Mb-associated decrease
   (alkanes, 2-pentylfuran, 1-hexanol, decanal, possibly Strecker aldehydes) is suspect of fibre displacement.
4. [!] **Implausibly small SDs**: hexanal GMA+1.0 % 14,297.20 ± 11.47 (RSD 0.08 %), saturated aldehydes GMA+1.0 %
   ± 41.23 (0.3 %), "Others" GMA ± 28.79 (0.3 %), total volatiles GMA ± 480.67 (1.6 %), against 5-40 % RSD elsewhere.
   A typo or analytical (one-batch) replicates.
5. [!] **Raw samples incubate 60 min at 45 °C in air with heme in the SPME vial**, after an unstated 4 °C hold; part
   of the "raw" Mb effect may form in the vial.
6. [!] Text vs table: "2-methoxy-4-methylphenol" (text p. 8) vs 2-methoxy-4-vinylphenol (Table 1, No. 31); "All
   remaining volatile compounds are generally found to increase" with Mb [Q], yet 1-hexanol, decanal, 2-pentylfuran,
   alkanes and Strecker aldehydes fall; "short-chain carbohydrates (<C14)" (p. 4) means hydrocarbons; RMA
   "unsaturated aldehydes (12.0 %)" does not match Table 1 (0.24 %).
7. [!] "heptanal and nonanal (saturated aldehydes derived from oleic acid)" [Q]: heptanal is not a canonical oleate
   product (cf. the Frankel oleate slates in `lipid.yml`); heptanal also does not rise with Mb in raw samples.
8. [!] Saturated aldehydes "display the most pronounced increase" [Q]: largest absolute increase, but the
   unsaturated class rises ×47 vs ×4.2 (fold).
9. [!] Tentative identifications with LRI offsets of 15-31 for four 2-alkenals and decanal (sec. 3).
10. [!] Mb purity / form / iron content and sunflower-oil type unknown; n = 3 replication unit undefined; one cooking
    condition; industry authorship and patent tie.
11. "Data is contained within the article" [V], but raw-ingredient (TSP) SPME data are "not shown".

## 6. Existing records in the repo that cite this paper

grep 2026-11-03 (`foods11131985`, `devaere`, `myoglobin addition`, `Van Loo`, `Paleo`, `MDPI Burger`) across
`data/ docs/ src/ scripts/ tasks/ results/`: no ledger, registry, prior or code record. Hits: the benchmark
narrative, three research-corpus files and the alias line (not checked).

### 6a. `data/benchmarks/maillard_validation_benchmarks.md` sec. 3.5 (l. 486-499) and reference table l. 767

| field | repo value | paper | verdict |
|---|---|---|---|
| citation authors | "**Van Loo et al.**" (l. 488 and l. 767) | **Devaere et al.** (Devaere, De Winne, Dewulf, Fraeye, Šoljić, Lauwers, de Jong, Sanctorum) | **mismatch** — no author named Van Loo |
| title / journal / vol / DOI | as printed, Foods 2022;11(13):1985, 10.3390/foods11131985 | same | **matches** |
| matrix | "Textured soy protein (TSP) burgers; HS-SPME-GC-MS; grilled 12 min at 250 °C" | TSP burgers, HS-SPME-GC-MS; 12 min at 250 °C **in a combi oven** | **matches the paper's own label**; "grilled" is the paper's misnomer |
| hexanal 10.6 ± 1.3 % of total area (raw) | 10.6 ± 1.3 % | 10.6 ± 1.3 % [M] | **matches** |
| hexane "artefact of hexane-extracted SPI" | SPI | **TSP**; solvent origin is the authors' conjecture [Q] | **partial** (protein form wrong; origin unproven) |
| "Grilling (250 °C, 12 min) increased total saturated aldehydes by ~3-4× vs raw" | ×3-4 from grilling | grilling **lowered** them ×0.49 without Mb, ×1.06-1.12 with Mb; the ×3.3-4.2 is the **myoglobin** effect in raw samples | **mismatch** (effect misattributed) |
| "Upper bound for aldehyde generation during high-heat searing" | benchmark use | relative peak areas, no IS, one condition | **not supported** |
| myoglobin | not mentioned in key findings | the paper's only factor besides cooking | **omission** |

### 6b. `data/research_corpus/elicit_plant_based_cooking.md` l. 186, 235, 261, 277, 297, 315

Authors, year, title, DOI (l. 315): **match**. "Grilling at 250 °C for 12 minutes" (l. 235): matches the paper's
label (oven in methods). "Myoglobin ... increasing aldehydes and decreasing hydrocarbons while enhancing pyrazine
formation" (l. 261): **matches** the authors (pyrazine increase only after cooking). "Myoglobin addition ... to
catalyze lipid oxidation reactions and enhance Maillard browning" (l. 277): lipid part matches [Q]; **"enhance
Maillard browning" is the authors' speculation ("possibly"); no browning was measured**. "Primary study" (l. 186):
matches. `data/research_corpus/maillard_plant_based.md` l. 215 and `data/research_corpus/raw/Gemini_deep_research_D.md`
l. 148 list the paper (PMC link) as references only; title matches.

## 7. Proposed ledger entries (NOT written into the ledgers)

Two entries: lipid volatiles → `lipid.yml`; pyrazines / Strecker / furans → `aroma_other.yml`. **None for
`sulfur.yml`**: no sulfur compound was detected or reported, so there is no data set (absence under full-scan
SPME is not a measurement). `quant_class: other` because there is no internal standard (pure TIC peak area).

```yaml
# lipid.yml
  - id: devaere2022_tab1_myoglobin_dose_lipid_volatiles
    dossier: devaere2022_extraction.md
    paper: "Devaere, De Winne, Dewulf, Fraeye, Soljic, Lauwers, de Jong & Sanctorum 2022 Foods 11:1985"
    lab: "KU Leuven Ghent (Fraeye) / Paleo (Sanctorum)"
    system: "TSP 25 % + sunflower oil 15 % + methylcellulose 1 % + salt 1.5 % + water 57.5 %, bovine myoglobin 0 / 0.5 / 1.0 % on top; raw vs oven-cooked"
    medium: "30 g patties, combi oven 250 C air 12 min (core T not measured); HS-SPME DVB/CAR/PDMS 3 g in 20 mL vial, 45 C 20 + 40 min"
    temperature_c: [250]
    ph: {initial: [], measured_during: false, final: null}
    time_min: [12]
    species: [hexanal, heptanal, nonanal, decanal, 2-heptenal, 2-octenal, "(E,E)-2,4-nonadienal", 2-decenal, "(E,Z)-2,4-decadienal", "(E,E)-2,4-decadienal", 2-undecenal, 1-octen-3-ol, 1-hexanol, 2,3-octanedione, 2-heptanone, 2-pentylfuran, octanoic acid, nonanoic acid]
    kind: endpoint
    location: "Table 1 pp. 5-7 (n = 3, AU x 10^3, Tukey letters; dossier sec. 4a, fold changes sec. 4c)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "baked for 12 min at 250 °C in an oven"
    flags: [matrix_food, peak_area_only, ph_not_stated, air_headspace, condition_uncertain, internal_inconsistency, no_heme_species_in_engine, spme_fibre_competition_suspected, industry_authored]
    note: "Heme dose-response (0/0.5/1 % bovine Mb) of linoleate-type volatiles in a plant burger, raw AND cooked: hexanal x3.7/x4.5 raw, x9/x12 cooked; 2-alkenals, 2,4-dienals, 1-octen-3-ol, 2,3-octanedione rise; 2-pentylfuran and hexane fall (fibre displacement suspected). Directional check for a future heme / iron term (must act at <=45 C too). No IS, no PV, oil type and Mb purity not stated; 'grilled' in the paper = combi oven."
```

```yaml
# aroma_other.yml
  - id: devaere2022_tab1_myoglobin_dose_pyrazines_strecker
    dossier: devaere2022_extraction.md
    paper: "Devaere, De Winne, Dewulf, Fraeye, Soljic, Lauwers, de Jong & Sanctorum 2022 Foods 11:1985"
    lab: "KU Leuven Ghent (Fraeye) / Paleo (Sanctorum)"
    system: "TSP 25 % + sunflower oil 15 % + methylcellulose 1 % + salt 1.5 % + water 57.5 %, bovine myoglobin 0 / 0.5 / 1.0 % on top; raw vs oven-cooked; no added sugar or amino acid"
    medium: "30 g patties, combi oven 250 C air 12 min (core T not measured); HS-SPME DVB/CAR/PDMS 3 g in 20 mL vial, 45 C 20 + 40 min"
    temperature_c: [250]
    ph: {initial: [], measured_during: false, final: null}
    time_min: [12]
    species: [methylpyrazine, 2,5-dimethylpyrazine, 2-ethyl-6-methylpyrazine, 2-ethyl-5-methylpyrazine, 2,5-dimethyl-3-ethylpyrazine, 2-methyl-3,5-diethylpyrazine, 2-methylbutanal, 3-methylbutanal, benzaldehyde, furfural, furfuryl alcohol, pyrrole, maltol, guaiacol, 2-methoxy-4-vinylphenol]
    kind: endpoint
    location: "Table 1 pp. 5-7 (n = 3, AU x 10^3, Tukey letters; dossier sec. 4a, 4c)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "baked for 12 min at 250 °C in an oven"
    flags: [matrix_food, peak_area_only, ph_not_stated, air_headspace, condition_uncertain, no_heme_species_in_engine, spme_fibre_competition_suspected, industry_authored, "complements: devaere2022_tab1_myoglobin_dose_lipid_volatiles"]
    note: "Pyrazines carried in by the extruded TSP (raw ~300 AU, Mb-independent); cooking x12, and Mb raises cooked pyrazines x1.33/x2.24 while 2-/3-methylbutanal fall x0.38. No sulfur compound detected in any arm. Peak areas only; one cooking condition."
```

## 8. Supplementary material

No Supporting Information, Appendix or Table S is referenced (confidence ~97 %; the article ends at the references).
Raw-ingredient SPME data are "not shown" and not deposited. The same data and possibly more (fermentation-derived
Mb, other doses) may be in **patent application PCT/EP2021/087884** (Paleo b.v.), not on disk; its WO publication
number is not printed. Nothing model-relevant is known to sit only in SI.
