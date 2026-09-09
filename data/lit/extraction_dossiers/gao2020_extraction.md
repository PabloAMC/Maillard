# Gao et al. 2020 — EXTRACTION (yellow pea protein isolate at alkaline extraction pH 8.5 / 9.0 / 9.5: composition, free SH and S-S, solubility, beany-marker volatiles by HS-SPME, LOX activity through the AE-IEP process)
### The only pea-isolate paper on disk that prints free thiol and disulfide densities; its hexanal and LOX series are figures.

**Source on disk:** `data/articles/gao2020.pdf` (owner's download, 2026-09-08). Read from the
scratchpad text layer (`gao2020.txt`, clean); pypdf confirms two tables (Table 1 p. 4, Table 2
p. 6) and six figures. Repo status before this dossier: `data/species/protein_matrices.yml`
carries only beta-lactoglobulin and `results/validation/matrix_sites_prereg.md` §1 states that
"no dossier on disk gives the free thiol and disulfide content of pea or soy isolates"; the roadmap
names LOX-made hexanal as needing its own data programme. This paper supplies the first (Table 2)
and touches the second (Fig. 5, Fig. 6 — figures, with two LOX values printed in the text).

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Effect of alkaline extraction pH on structure properties, solubility, and beany flavor of yellow pea protein isolate" |
| Authors | Zili Gao, Peiyi Shen, Yang Lan, Leqi Cui, Jae-Bom Ohm, Bingcan Chen, Jiajia Rao (corresponding) — North Dakota State University / USDA-ARS Fargo |
| Venue | Food Research International 131 (2020) 109045; received 14 Jun 2019, revised 24 Jan 2020, accepted 27 Jan 2020, online 29 Jan 2020 |
| DOI | 10.1016/j.foodres.2020.109045 (PII S0963-9969(20)30070-3) |
| Naming | PPI = alkaline-extraction / isoelectric-precipitation isolate, neutralised to pH 7.0 and freeze-dried; "AE", "IEP", "PPI" in Fig. 6 = the three process stages (alkaline supernatant, pH-4.5 precipitate, final neutralised isolate), each lyophilised. "2-methoxy-3-isopropyl pyrazine" / "3-isopropyl-2-methoxypyrazine" / "2-methoxy-3-isoproply pyrazine" are the same marker. "SS bounds" (Table 2) = S-S bonds. |
| Material | Yellow field pea cv. Nette 2010, 2018 crop (Montana State University breeding programme); dehulled, milled to 250 um; flour 22.25 % protein. |
| Compound registry | hexanal -> `hexanal`; 1-octen-3-ol -> `1_octen_3_ol`; 1-pentanol, 3-methyl-1-butanol, 1-octen-3-one -> not in registry; 3-isopropyl-2-methoxypyrazine -> not in registry (registry has `3_isobutyl_2_methoxypyrazine` and group `methoxypyrazines`). |

## 1. Why it matters

Need (a). The matrix layer charges the sulfur lane's protein-disulfide pool and binds HMF to thiol
and amine pools from site densities in mmol/g, and refuses to run for any isolate that has no
sites on file. Table 2 of this paper is a direct measurement, by Ellman's reagent in 8 M urea, of
exposed free SH, total SH and S-S in a lab-made yellow-pea AE-IEP isolate at three extraction pHs:
free SH 12.25-14.72 umol/g, S-S 18.66-24.80 umol/g (per g isolate powder, 83-85 % protein). That is
0.0123-0.0147 mmol/g free thiol and 0.0187-0.0248 mmol/g disulfide — a third and a fifth,
respectively, of the BLG values the engine carries (0.0545 and 0.109 mmol/g). It says nothing about
amine sites (not measured) and nothing about aldehyde binding constants. Need (b). The paper asks
the right question for the "beany note before any heat" — how much hexanal and LOX activity survive
alkaline extraction at 8.5 / 9.0 / 9.5 into the isolate — but answers it in figures: the six
beany-marker peak areas (Fig. 5) and the LOX activities at the three process stages (Fig. 6) are
FIGURE-ONLY, with only two LOX values (43.7 and 18.6 U/g for the pH-9.0 route) and the qualitative
hexanal pattern printed in the text. No time course, no hexanal concentration in any unit.

## 2. Methods as they matter to a model

- **Flour.** Dehulled Nette 2010, Retsch ZM 200 at 8,000 rpm, 250 um screen, vacuum-cooled;
  22.25 % protein.
- **Isolate (AE-IEP).** 70.0 g flour in water 1:15 w/v (= 1,050 mL; **66.7 g flour/L, 14.8 g
  protein/L in the slurry**), pH "quickly adjusted" to 8.5, 9.0 or 9.5 with 1.0 M NaOH, stirred
  600 rpm for **60 min** (temperature not stated; flag 4), centrifuged 6000 rpm 20 min, filtered
  (Whatman 1), supernatant to **pH 4.5**, precipitate re-suspended in water and brought to **pH 7.0**
  with 1.0 M NaOH, freeze-dried 48 h, stored 4 C in glass. Three flour batches mixed before use.
- **Composition.** Kjeldahl N x 6.25 (not 5.7 as in the Wang papers); moisture and ash AOAC; lipid
  by accelerated solvent extraction (hexane, 10 MPa, 3 x 20 min); carbohydrate by difference.
  Extraction yield = PPI mass / flour mass; protein recovery = protein in PPI / protein in flour.
- **Free SH / total SH / S-S (verbatim core):** "bulk protein solution was prepared by dissolving
  75 mg of protein samples in 10 mL of Tris-Gly buffer (86 mM Tris, 90 mM glycine and 4 mM EDTA,
  pH 8.0) containing 8 M urea and then gently stirred overnight." Free SH: 1 mL + 4 mL buffer +
  0.05 mL Ellman's (2 mM DTNB), 20 min RT, A412. Total SH: 1 mL + 0.05 mL beta-mercaptoethanol +
  4 mL buffer 1 h, 12 % TCA precipitation, wash, redissolve in 10 mL buffer, Ellman's, A412.
  "SH (umol/g) = 73.53 x A x D / C", C = 7.5 mg/mL, D = 5 (free) or 10 (total), 73.53 = 10^6 /
  1.36 x 10^4 (epsilon 13,600 M^-1 cm^-1). "SS (umol/g) = (Total SH - free SH) / 2". The mass basis
  is **g of isolate powder** ("75 mg of protein samples"); the 8 M urea means "exposed free SH" is
  the unfolded protein's free cysteine, i.e. all non-disulfide cysteine, not a native-surface
  count (flag 2).
- **SEC-MALS-RI.** 3 mg/mL in 10 mM phosphate pH 7.0, mobile phase + 100 mM NaCl, dn/dc 0.185;
  peak 1 > 2000 kDa aggregates 10-18 %; peak 2 ~400 kDa (legumin); peaks 3-4 180-200 kDa
  (vicilin). Proportions Fig. 2B (FIGURE-ONLY).
- **SDS-PAGE (non-reducing).** Convicilin ~20 %, vicilin ~27 %, legumin ~45 % of PPI by ImageJ
  (text); "A light band at around 94 kDa that has previously been recognized as pea seeds
  lipoxygenase (LOX) also existed in our extracted PPI samples".
- **CD.** 0.1 mg/mL, 10 mM PBS pH 7.4, 25 C; beta-sheet dominant, no beta-turn; fractions in the
  Fig. 3 inset (FIGURE-ONLY).
- **Solubility.** 1.00 wt% PPI at pH 7.0 or 3.5, 2 h stirring RT, 9,100 g 10 min, Bradford on the
  supernatant vs total protein in 0.2 M NaOH. Text: pH 7.0 solubility 93.6 % -> 80.2 % as extraction
  pH 8.5 -> 9.5; the rest is Fig. 4 (FIGURE-ONLY).
- **Volatiles (verbatim core):** "the PPI sample (1.00 g) was mixed with 2.0 mL of 20% NaCl solution
  in 20 mL GC vials ... incubated at 60 C for 3 min ... while agitating at 250 rpm. Thereafter, the
  SPME fiber needle (DVB/CAR/PDMS ...) was inserted into the vial for 60 min to absorb volatiles at
  60 C and then transferred to the injector port (250 C) for 3 min of desorption." ZB-WAX 60 m x
  0.25 mm x 0.25 um; He 2 mL/min; 40 -> 85 C at 45 C/min -> 200 C at 9 C/min -> 250 C at 45 C/min,
  hold 3 min; EI 70 eV, m/z 40-350. So a **33 % w/v isolate slurry in ~13 % NaCl, 63 min at 60 C**,
  no internal standard, no calibration.
- **Quantification sentence (verbatim):** "The selected beany flavor markers were confirmed using the
  NIST 14 Library mass spectral database and their absolute peak areas in area counts were
  recorded." Six markers: 1-pentanol, 1-octen-3-ol, 3-methyl-1-butanol, hexanal, 1-octen-3-one,
  2-methoxy-3-isopropyl pyrazine; samples = pea flour and the three PPIs. Results = Fig. 5, peak
  areas, FIGURE-ONLY.
- **LOX activity (verbatim core):** "10.0 mL of phosphate buffer (pH 6.5) was added to 0.1 g PPI,
  and then stirred magnetically at room temperature (23 C) for 5 hrs. The mixture was then
  centrifuged at 9100g for 10 min and the supernatant was used as the crude enzyme extract ...
  linoleic acid (140 uL) and Tween 20 (140 uL) were mixed and emulsified into 8 mL of phosphate
  buffer (pH 6.5). Then, 1.1 mL of 0.5 M NaOH was added to clarify the solution, and the volume
  was brought to 50 mL with phosphate buffer (pH 6.5). The stock substrate solution was flushed with
  nitrogen ..., and was diluted (1:40, v/v) with 0.2 M sodium borate buffer (pH 9.0) before use. ...
  50 uL of the crude enzyme extract was added to a quartz cuvette contained 1.25 mL of the
  substrate with rapid mixing for 5 s, and the change in absorbance was recorded for 3 min. The
  unit of LOX activity was U/g, where U is defined as the numeric increase in absorbance per
  minute." Arithmetic: 140 uL x 0.902 g/mL = 126 mg = 0.450 mmol in 50 mL = 9.0 mM stock; 1:40 ->
  0.225 mM; in the 1.30-mL cuvette **0.216 mM linoleic acid in 0.2 M borate, pH ~9.0, 23 C**, with
  0.5 mg of isolate-equivalent extract (10 mg/mL x 50 uL). "U/g" is therefore A234 per min per gram
  of the lyophilised sample extracted, i.e. 1 U/g = 0.0005 A/min in the cuvette. Note the assay pH
  (9.0) against the paper's own citation that pea LOX peaks at pH 5.5 (Szymanowska 2009) and Zhang
  2020b's finding that pea carries only the pH-6.8 / 7.1 isoforms (flag 6).
- **Samples for LOX.** Lyophilised alkaline supernatant (AE stage, pH 8.5 / 9.0 / 9.5), lyophilised
  pH-4.5 precipitate (IEP stage), final PPI — each at "0.1 g" into 10 mL, so "U/g" is per gram of
  three different solids (flag 7).
- **Statistics.** CRD, two independent experiments, duplicate or triplicate measurements, ANOVA +
  Tukey p < 0.05.

## 3. Tables re-typed

### Table 1. "The yield and proximate composition analysis of pea protein isolate." (proximate on % wet basis)

| Extraction pH | Extraction yield (%) | Protein recovery yield (%) | Crude protein (%) | Lipid (%) | Carbohydrate (%) | Moisture (%) | Ash (%) |
|---|---|---|---|---|---|---|---|
| 8.5 | 12.93 +/- 0.91 b | 49.20 +/- 1.15 b | 84.67 +/- 0.13 a | 1.47 +/- 0.04 a | 4.08 +/- 0.18 a | 5.20 +/- 0.01 a | 4.58 +/- 0.01 a |
| 9.0 | 14.00 +/- 0.81 ab | 52.43 +/- 1.51 ab | 83.33 +/- 0.45 a | 1.46 +/- 0.02 a | 5.47 +/- 0.41 a | 5.17 +/- 0.06 a | 4.59 +/- 0.05 a |
| 9.5 | 15.36 +/- 0.51 a | 57.56 +/- 1.89 a | 83.40 +/- 0.73 a | 1.47 +/- 0.03 a | 5.27 +/- 0.73 a | 5.19 +/- 0.01 a | 4.68 +/- 0.04 a |

Letters: column-wise significance (p < 0.05). Protein by Kjeldahl N x 6.25. Row sums of the five
proximate columns: 100.00 / 100.02 / 100.01 (carbohydrate is by difference).

### Table 2. "The sulfhydryl group (SH) and disulfide bond (SS) contents of pea protein isolate." (umol per g isolate powder; 8 M urea)

| Extraction pH | Exposed free SH (umol/g) | Total SH (umol/g) | SS bonds (umol/g) |
|---|---|---|---|
| 8.5 | 14.72 +/- 0.23 a | 52.04 +/- 0.40 c | 18.66 +/- 0.20 c |
| 9 | 12.94 +/- 0.11 b | 55.27 +/- 0.70 b | 21.15 +/- 0.35 b |
| 9.5 | 12.25 +/- 0.50 b | 61.85 +/- 0.30 a | 24.80 +/- 0.15 a |

Letters: column-wise significance. Check of Eq. 4: (52.04 - 14.72)/2 = 18.66; (55.27 - 12.94)/2 =
21.165; (61.85 - 12.25)/2 = 24.80 — consistent. Converted to a protein basis with Table 1's protein
contents (0.8467 / 0.8333 / 0.8340): free SH **17.4 / 15.5 / 14.7 umol/g protein**; S-S **22.0 /
25.4 / 29.7 umol/g protein**; total half-cystine (free + 2 x S-S) 61.5 / 66.3 / 74.2 umol/g protein
= 0.75 / 0.80 / 0.90 g Cys / 100 g protein.

### Numbers printed only in the running text (Figs 1-6 otherwise FIGURE-ONLY)

| item | as printed | where |
|---|---|---|
| Subunit shares (SDS-PAGE, ImageJ) | convicilin ~20 %, vicilin ~27 %, legumin ~45 % of PPI; LOX band ~94 kDa present | §3.2.1 |
| SEC aggregates | peak 1 (> 2000 kDa) 10-18 % of protein | §3.2.2 |
| Solubility at pH 7.0 | 93.6 % (pH 8.5 extraction) -> 80.2 % (pH 9.5) | §3.3 |
| Hexanal (peak area) | "In pea flour, hexanal presented as the dominating aromatic compound ... PPI extracted at pH 9.5 showed similar level of hexanal as pea flour, whereas extraction at pH 8.5 and 9.0 greatly reduced the amounts of hexanal in PPI." | §3.4.1 |
| 1-Octen-3-one | "appeared in a small quantity in all PPI samples without any statistically difference" | §3.4.1 |
| Alcohols | 1-pentanol higher in flour than in PPI; 3-methyl-1-butanol higher in PPI than flour; all three alcohols highest in the pH-9.5 PPI | §3.4.1 |
| Methoxypyrazine | "Extraction process of PPI significantly reduced the content of 3-isopropyl-2-methoxypyrazine" | §3.4.1 |
| Overall | "PPI obtained at pH 9.0 contained the lowest amounts of beany odor associated compounds. However, increasing the alkaline extraction pH to 9.5 drastically escalated the content of beany odor in PPI." | §3.4.1 |
| LOX across stages | "the LOX activities in all samples maintained steady when processed from AE step at an alkaline pH to IEP extraction step where pH became 4.5 (p > 0.05)"; "with the PPI prepared via alkaline extraction at pH 9.0, the LOX activity rose up to 43.7 U/g at IEP step and dropped to 18.6 U/g in the final PPI product"; "the activity of LOX in the samples obtained at pH 9.0 always exhibited the lowest activity compared to the others over the course of extraction" | §3.4.2 |
| Cited pea LOX optimum | pH 5.5 (Szymanowska et al. 2009) | §3.4.2 |
| Lipid in Nette 2010 | "relatively low (~1.5%)" (this is the PPI lipid of Table 1, not the seed) | §3.4.2 |

## 4. Numbers the repository can use

| quantity | value | unit | conditions | source location | evidence class | fit |
|---|---|---|---|---|---|---|
| **Free thiol, pea AE-IEP isolate** | 14.72 +/- 0.23 / 12.94 +/- 0.11 / 12.25 +/- 0.50 (pH 8.5 / 9.0 / 9.5) | umol per g isolate powder | Ellman's in 8 M urea, pH 8.0, RT; = 0.0147 / 0.0129 / 0.0123 mmol/g powder; 0.0174 / 0.0155 / 0.0147 mmol/g protein | Table 2 | measured | not a `binding_constants.yml` type; **the `protein_matrices.yml` `free_thiol` site density for a pea isolate**, with `protein_basis: g_isolate_powder` and the protein fractions of Table 1. Choice of extraction pH is a declared band (0.0123-0.0147). |
| **Disulfide, pea AE-IEP isolate** | 18.66 +/- 0.20 / 21.15 +/- 0.35 / 24.80 +/- 0.15 | umol S-S per g powder | as above; 0.0187 / 0.0212 / 0.0248 mmol/g powder; 0.0220 / 0.0254 / 0.0297 mmol/g protein | Table 2 | measured | **the `disulfide` site density for a pea isolate** (charges `PROT_SS`); band 0.0187-0.0248 mmol/g powder |
| Total SH (after reduction) | 52.04 / 55.27 / 61.85 | umol per g powder | as above | Table 2 | measured | total cysteine; 0.75-0.90 g Cys / 100 g protein vs the 0.35 g / 100 g that Wang 2014 cites for a salt-extracted pea isolate (flag 3) |
| Amine site density | NOT MEASURED | — | — | — | — | the amine pool for pea remains unfilled by this paper |
| Protein content of PPI | 84.67 / 83.33 / 83.40 | % (N x 6.25, wet basis) | freeze-dried | Table 1 | measured | `protein_purity_fraction` for the rows above (N x 6.25; the Wang isolates used 5.7 — not interchangeable) |
| Lipid in PPI | 1.47 / 1.46 / 1.47 | % | hexane ASE | Table 1 | measured | substrate carried into the isolate (matters for any storage / LOX module) |
| Extraction yield; protein recovery | 12.93 / 14.00 / 15.36 %; 49.20 / 52.43 / 57.56 % | % | 1:15 w/v, 60 min, IEP 4.5 | Table 1 | measured | process datum |
| Solubility at pH 7.0 | 93.6 (pH 8.5) -> 80.2 (pH 9.5) | % | 1 wt%, 2 h, 9,100 g | §3.3 text | measured (two of three printed) | matrix datum |
| LOX activity, pH-9.0 route, IEP stage | 43.7 | U/g (A234 min^-1 per g lyophilised precipitate; assay pH 9.0, 23 C, 0.22 mM linoleic acid) | after AE 60 min + IEP 4.5 | §3.4.2 text (Fig. 6) | measured, level_only (printed in text) | none; unit conversion needs epsilon(234) and is not promoted (flag 8) |
| LOX activity, pH-9.0 route, final PPI | 18.6 | U/g PPI, as above | after neutralisation to 7.0 and freeze-drying | §3.4.2 text (Fig. 6) | measured, level_only | as above; **residual LOX in a finished pea isolate is non-zero** |
| LOX ratio IEP -> final PPI (pH-9.0 route) | 18.6 / 43.7 = 0.43 | x | same solid series | derived from the two printed values | within-study ratio | a loss factor for the neutralise-and-dry step, one route only |
| LOX at pH 8.5 and 9.5 routes, all stages | FIGURE-ONLY (pH 9.0 stated lowest) | U/g | | Fig. 6 | figure_only | none |
| Hexanal in flour and PPIs | FIGURE-ONLY | peak area counts | 1.00 g sample + 2 mL 20 % NaCl, 60 C, 63 min SPME | Fig. 5 | figure_only | none; qualitative order: flour ~ PPI(9.5) >> PPI(8.5), PPI(9.0) |
| 1-Pentanol, 1-octen-3-ol, 3-methyl-1-butanol, 1-octen-3-one, methoxypyrazine | FIGURE-ONLY | peak area counts | as above | Fig. 5 | figure_only | none |
| Subunit shares | convicilin ~20, vicilin ~27, legumin ~45 | % of PPI (ImageJ) | non-reducing SDS-PAGE | §3.2.1 | measured, approximate | matrix composition datum |
| Aggregates > 2000 kDa | 10-18 | % of protein | SEC-MALS, pH 7.0 | §3.2.2 | measured (range) | none |

Comparison the matrix layer can print (arithmetic, not a fit): BLG in `protein_matrices.yml` has
1 free thiol and 2 S-S per 18,362 Da = 0.0545 and 0.1089 mmol/g; this pea isolate has 0.23-0.27x
the free thiol and 0.17-0.23x the disulfide per gram of powder. Any pea run through the sulfur
lane's `PROT_SS` pool would therefore charge roughly a fifth of what the same g/L of BLG charges.

## 5. Flags

1. **Hexanal and every beany marker are FIGURE-ONLY** (Fig. 5), and would be peak areas in area
   counts without an internal standard or calibration even if read — not concentrations. The
   comparison is per gram of sample (1.00 g flour at 22 % protein vs 1.00 g PPI at 84 %), so "PPI
   at pH 9.5 showed similar level of hexanal as pea flour" is per gram of solid, not per gram of
   protein or per gram of lipid.
2. **"Exposed free SH" was measured in 8 M urea**, so it is the free cysteine of the unfolded
   protein (all non-disulfide Cys), which is the right quantity for a site density but not a
   native-surface count. Total SH was measured after beta-mercaptoethanol reduction and TCA
   precipitation; residual beta-ME is a known positive bias in that protocol. The basis is g of
   powder ("75 mg of protein samples"), converted to g protein above with Table 1.
3. **Cysteine disagreement with Wang 2014's cited value.** Total half-cystine here is 0.75-0.90 g /
   100 g protein; Wang 2014 cites 0.35 g / 100 g protein (Khattab 2009) for a salt-extracted pea
   isolate. Different cultivar, extraction and method (Ellman vs amino-acid analysis); a factor of
   2-2.6. The site density entered for pea should carry this as its uncertainty band, not pick one.
4. **Extraction temperature not stated** for the 60-min alkaline step (room temperature by
   implication), nor the temperature of the IEP and neutralisation steps; the LOX-relevant history
   of the slurry (time at each pH, temperature) is therefore only partly known: 60 min at pH
   8.5-9.5, then pH 4.5, then pH 7.0, then 48 h freeze-drying.
5. **LOX unit differs from Zhang 2020b's by 1000 and by basis.** Here 1 U = 1.0 A234/min per g of
   solid extracted (0.5 mg in the cuvette); in Zhang 2020b 1 U = 0.001 A234/min per mg of milk
   protein. In absorbance-rate terms the Zhang pea milk LOX-2 (2.16 A/min per mg protein) exceeds
   this final PPI (18.6 A/min per g = 0.0186 A/min per mg powder) by ~100x before any correction for
   assay pH (6.8 vs 9.0), substrate (2.4 mM vs 0.22 mM), reaction volume (30 mL vs 1.3 mL) or
   material (raw milk vs IEP isolate); with epsilon(234) = 25,000 and those volumes the molar rates
   are ~2.6 umol/min/mg vs ~0.001 umol/min/mg, a factor ~2,700. The two numbers are not comparable
   and must not be entered as one quantity; both are recorded in their printed units.
6. **Assay pH 9.0 for a pea LOX.** The substrate was diluted 1:40 into 0.2 M borate pH 9.0, the
   soybean LOX-1 optimum. The paper itself cites pea LOX optimum pH 5.5; Zhang 2020b finds pea has no
   LOX-1-type (pH 9) activity and assays pea at pH 6.8 / 7.1. The activities here are therefore
   likely far below the enzyme's capacity at the pH of the extraction slurry, and the stage-to-stage
   comparison assumes the pH-9 activity tracks the pH-6.8 one.
7. **"U/g" is per gram of three different solids** (lyophilised alkaline supernatant, lyophilised
   pH-4.5 curd, final PPI), so the AE -> IEP -> PPI series mixes bases; the printed 43.7 -> 18.6 is
   IEP curd -> PPI, the two most similar solids, and even there the neutralisation adds NaOH salt to
   the mass.
8. **Text vs figure wording.** "maintained steady when processed from AE step ... to IEP" and "rose
   up to 43.7 U/g at IEP step" sit in adjacent sentences; the Tukey letters resolving this are in
   Fig. 6 (FIGURE-ONLY).
9. **LOX extraction step is itself long:** 0.1 g PPI stirred 5 h at 23 C in pH 6.5 buffer before
   assay; with 1.47 % lipid in the PPI, some substrate turnover (and enzyme inactivation by its own
   hydroperoxides) can occur before the cuvette.
10. **The SPME step can make hexanal.** A 33 % w/v isolate slurry with active LOX (Fig. 6 shows it
    is non-zero in the final PPI) held 63 min at 60 C in 13 % NaCl: part of the hexanal signal may
    form in the vial. The pH-9.0 isolate having both the lowest LOX and the lowest hexanal is
    consistent with either in-vial formation or carry-over from processing; the paper cannot
    separate them.
11. **No time information on volatiles**, no concentration units, no isolate mass balance for
    hexanal; the "PPI extracted at pH 9.0 possessed the lowest beany flavor" claim is a ranking of
    peak areas.
12. **N x 6.25 here vs N x 5.7 in the Wang isolates**; a `protein_purity_fraction` taken from this
    paper is ~10 % higher than the same nitrogen would give under the Wang convention.
13. **Gao cites "Wang & Arntfield 2015" (FRI 77:1-9)** for the PPIa > PPIs aldehyde-binding result,
    which is in Wang & Arntfield 2014 (see `wang2014_extraction.md` §5).
