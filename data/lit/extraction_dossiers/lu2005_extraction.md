# Lu, Hao, Payne & Ho 2005 — EXTRACTION (glucose 40 mM + glycine, diglycine or triglycine 40 mM in glycerol / water 0-100 % (v/v), 160 °C, 1 h; volatiles by SBSE-GC-MS relative to an internal standard; remaining Gly, Gly-Gly, triGly and cyclo(Gly-Gly) by LC-MS/MS; peptides alone (no glucose) as a control; A440)

### A water-content x amine-form screen at one temperature and one time. Pyrazines rise ~x15 from 100 % to 10-25 % water; the peptide effect runs through peptide breakdown (triglycine → DKP + free glycine). Volatile numbers are relative areas, not concentrations; no pH, no time course.

**Read 2026-10-31 from `data/articles/lu2005.pdf`** (5 pp., J. Agric. Food Chem. 53:6443-6447; born-digital, clean text
layer). Tables 1-2 read from the text layer and checked against 130 dpi renders (scratchpad `dos7/lu-3.png`,
`lu-4.png`). Figs. 2-4 checked against the renders; Fig. 3 (peptides without glucose) exists only as a figure, read at
±3 mol %.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| what the model needs | what Lu 2005 gives | tag |
|---|---|---|
| peptide vs free amino acid as Maillard amine | volatile output glycine ≈ triglycine > diglycine (by pyrazines, relative areas) [Q]; at 0 % water diglycine gives the most 2,5-DMP [M] | [M] / [Q] |
| peptide breakdown at 160 °C | **Table 2 (with glucose, mol/mol %):** triglycine → DKP 47-81 % and free Gly 2.8-60 %, triGly left ≤ 0.07 %; diglycine → DKP 17-20 %, Gly 2.9-12.8 %, diGly left 0.7-28.5 % | [M] |
| does glucose consume the released glycine? | triGly without glucose (Fig. 3b) releases ~89-91 % Gly at 50-100 % water vs 41-60 % with glucose [F/M] → ~30-50 mol % of the released glycine reacts with glucose [Z] | [F] / [Z] |
| water-activity / water-content dependence | pyrazines (all three amines) highest at 10-25 % water, ~x5-15 above 100 % water; slight inhibition at 0 % (viscosity [Q]); browning A440 falls monotonically from 0 % to 100 % water (Fig. 4) | [M] / [F] |
| pyrazine pattern | trimethylpyrazine and 2,5-DMP dominate everywhere; methylpyrazine and 1-methylpyrroles appear only at ≤ 25 % water | [M] |
| rate constants, Ea, pH | **none.** One temperature (160 °C), one time (1 h), no pH printed or controlled | [!] |
| absolute amounts | **none for volatiles** (peak area / IS area × 10, compound-specific SBSE efficiency shown in Fig. 1); peptides quantified by external standard | [M] |

**Answer.** For the model's amine side: (i) at 160 °C peptides are not inert carriers — triglycine is almost fully
consumed in 1 h, mostly by cyclisation to cyclo(Gly-Gly) with release of glycine, and that released glycine is what
reacts (the authors' reading [Q], consistent with the no-glucose control [Z]); (ii) diglycine is more stable to
hydrolysis and gives fewer pyrazines except in near-dry glycerol. A peptide-bound amine therefore enters the Maillard
lane through a release step whose rate here is "near-complete in 1 h at 160 °C" for triGly and partial for diGly —
a bound, not a constant. Water content changes pyrazine output by an order of magnitude at fixed T, but in glycerol,
which is also a solvent change.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of Water Content on Volatile Generation and Peptide Degradation in the Maillard Reaction of Glycine, Diglycine, and Triglycine" |
| Authors | Chih-Ying Lu, Chi-Tang Ho* (Rutgers University); Zhigang Hao, Richard Payne (Colgate-Palmolive) |
| Venue | *J. Agric. Food Chem.* **2005**, 53 (16), 6443-6447; received 9 Mar 2005, revised 2 Jun 2005, accepted 5 Jun 2005, web 7 Jul 2005 |
| DOI | **10.1021/jf050534p** (printed; Crossref 2026-10-31: same title, Lu / Hao / Payne / Ho, vol. 53, 6443-6447). Matches the brief |
| Lab | Ho network (Rutgers) |

## 2. Conditions

| item | as printed |
|---|---|
| reactants | glucose 0.0002 mol + glycine, diglycine or triglycine 0.0002 mol (each heated separately), total 5 mL → **40 mM + 40 mM** [M/Z] |
| medium | glycerol (≥ 99.5 %) with **0, 10, 25, 50, 75, 100 % (v/v) water**; "0 %" = 1.05 % moisture by Karl Fischer [M] (Fig. 3 caption says 0.05 % [!]) |
| heating | oil bath **160 °C, 1 h**; vial 21 × 70 mm; sealing / headspace not stated [!] |
| pH / buffer | **not stated, not controlled** [!] (glycerol media; peptide pK2 9.8 / 8.25 / 7.91 quoted [C]) |
| internal standard | 4-methylpyrimidine 10 ppm, "surrogate IS for semiquantitation"; whether it was added before or after heating is not stated [!] |
| SBSE | 1 mL reaction product + 9 mL phosphate buffer pH 8 with 27 % NaCl; PDMS Twister 10 mm × 0.5 mm, RT, 1100 rpm, 1 h; thermal desorption, Varian 3800 / Saturn 2000 ion trap, ZB-5ms [M] |
| LC-MS/MS | 1 mL product in 200 mL MeOH/water + 0.5 % formic acid; silica column; SRM; external standards [M] |
| replication | volatiles in triplicate; non-volatiles "in duplicate", error bars "standard error of three replicates of analysis" [M] (two statements, ambiguous [!]) |
| peptides-only control | peptides heated without glucose, 160 °C, 1 h, 0 / 50 / 100 % water, duplicate (Fig. 3) [M] |

## 3. Quantitative results

**Table 1 (selected; relative amount = area / IS area × 10; G / Di / Tri = glycine / diglycine / triglycine with
glucose)** [M]:

| compound | 0 % G / Di / Tri | 10 % | 25 % | 50 % | 75 % | 100 % |
|---|---|---|---|---|---|---|
| pyrazine | 4.6 / 0.4 / 1.7 | 11.1 / 1.3 / 6.4 | 8.3 / 0.9 / 7.0 | 4.4 / 0.7 / 1.8 | 3.7 / 0.3 / 1.1 | 1.2 / ND / 0.4 |
| methylpyrazine | 33.7 / 3.0 / 4.7 | 40.3 / 6.8 / 11.0 | 12.7 / 1.4 / 12.7 | 1.6 / 0.1 / 0.7 | 0.9 / 0.0 / 0.2 | 0.2 / ND / trace |
| 2,5-dimethylpyrazine | 38.7 / 77.2 / 58.7 | 76.2 / 39.1 / 71.8 | 77.7 / 27.7 / 62.2 | 36.5 / 8.6 / 21.9 | 20.7 / 2.9 / 9.2 | 5.3 / ND / 1.1 |
| 2,3-dimethylpyrazine | 28.3 / 14.6 / 33.9 | 32.5 / 15.1 / 51.7 | 31.1 / 5.7 / 30.1 | 12.7 / 0.3 / 8.6 | 5.7 / 0.1 / 2.5 | 1.9 / trace / 0.6 |
| trimethylpyrazine | 46.3 / 63.4 / 80.6 | 68.0 / 41.8 / 110.4 | 121.9 / 23.3 / 91.6 | 107.5 / 5.5 / 58.5 | 55.2 / 1.4 / 25.7 | 12.0 / ND / 1.8 |
| tetramethylpyrazine | 10.0 / 7.4 / 21.6 | 18.3 / 11.3 / 38.0 | 40.3 / 6.3 / 33.7 | 37.1 / 0.5 / 22.9 | 12.7 / 0.1 / 4.0 | 1.9 / ND / 0.3 |
| furfural | 3.8 / 3.1 / 1.7 | 3.1 / 13.1 / 3.6 | 4.1 / 10.7 / 3.8 | 2.1 / 9.4 / 2.4 | 3.7 / 7.6 / 2.9 | 6.2 / 6.8 / 4.7 |
| 2-acetylpyrrole | ND / ND / ND | 2.3 / 0.4 / 1.5 | 2.1 / 0.4 / 2.5 | 1.5 / 0.2 / 1.3 | 1.2 / 0.1 / 0.6 | 0.4 / trace / 0.2 |

Full table (20 compounds incl. ethyl-substituted pyrazines, 1-methylpyrroles, 4,5-dimethyloxazole, 5-methylfurfural,
HMF) on p. 6444. No SD printed for Table 1 [!]. Furfural from glucose is an odd product (a pentose marker) —
uncommented by the authors [!].

**Table 2 (with glucose, mol/mol % of the peptide charged; LC-MS/MS)** [M]:

| | water % | 0 | 10 | 25 | 50 | 75 | 100 |
|---|---|---|---|---|---|---|---|
| diGly + Glc | DKP | 16.920 | 19.345 | 19.786 | 20.420 | 17.895 | 18.107 |
| | glycine | 2.940 | 6.036 | 9.698 | 12.756 | 10.667 | 9.387 |
| | diglycine | 0.693 | 1.488 | 3.787 | 12.139 | 23.841 | 28.505 |
| | triglycine | 0.009 | 0.010 | 0.011 | 0.014 | 0.016 | 0.023 |
| triGly + Glc | DKP | 47.019 | 59.563 | 67.766 | 74.496 | 79.469 | 80.656 |
| | glycine | 2.755 | 4.149 | 15.432 | 41.243 | 55.122 | 59.867 |
| | diglycine | 0.505 | 0.607 | 0.929 | 1.475 | 1.729 | 1.742 |
| | triglycine | 0.070 | 0.035 | 0.035 | 0.056 | 0.065 | 0.061 |

Mass balance at 100 % water [Z]: diGly system, diGly-equivalents accounted 51 % (rest Maillard products); triGly
system, glycine units accounted 75 % of 3 per molecule. Text values: diGly → DKP "about a 20 %" yield; triGly → DKP
50-80 % [M].

**Fig. 3 (peptides without glucose, 160 °C, 1 h; bars unlabelled, assigned diGly = dotted, triGly = hatched by
matching Table 2's DKP levels [!])** [F, ±3 mol %]: DKP 43 / 29 / 25 (diGly) and 74 / 80 / 91 (triGly) at 0 / 50 /
100 % water; glycine 16 / 29 / 16 (diGly) and 35 / 91 / 89 (triGly).

**Fig. 4 (A440)** [F]: all three systems ~0.51-0.56 at 0 % water, ~0.40 at 10 %, ~0.29-0.32 at 25 %, 0.12-0.18 at
50 %, 0.05-0.11 at 75 %, 0.05-0.09 at 100 %; diGly highest at ≥ 25 % water, triGly lowest.

## 4. Mechanism evidence and its kind

- **Kind 6 (within-study factors):** water content (6 levels) and amine form (3); a no-glucose control for peptide
  breakdown (Fig. 3). Ho network (Rutgers) with Colgate-Palmolive.
- **Mass balance:** partial (Table 2 non-volatiles), no Amadori / glycated-peptide species measured.
- DKP from triGly by N-terminal cyclisation and C-terminal cleavage (Rizzi 1989) [C/P]; "a high hydrolysis rate of
  the peptide bond may be one of the main factors that enhance volatile formation" [P].
- 2,5-DMP from pyruvaldehyde [C]; ethylpyrazines from amino-acid fragments [C].

## 5. Bearing on open questions

- **Pyrazine / Strecker precursors (evidence map 2.7):** LEAVES UNCHANGED for carbon / nitrogen origin. Adds that a
  peptide's pyrazine output tracks its free-glycine release [Q/Z]; useful for any peptide or protein-hydrolysate amine
  species, which the engine does not carry.
- **Evidence map 1.1 (amine identity, peptides) and 1.6 (water activity; "pyrazines: Reineccius only"):** adds a
  second pyrazine x water series (Ho network, glycerol, 160 °C). Directional only: relative areas, one time point,
  solvent changes with water content, no a_w measured. Does not lift "pyrazines on an a_w ladder" above declared band.
- **Matrix / peptide effects (STATUS_AND_PLAN triage 9, sensitivity page):** supports carrying peptide-bound
  glycine as a RELEASE-LIMITED amine (fast for triGly, slow for diGly at 160 °C) rather than as free amine [Z].
- **CF ids:** none touched directly (no register entry is about peptides or water content).

## 6. Disagreements found

1. **xia2022 (evidence map 1.1, Ho network / Jiangnan: "Glycine, Diglycine, and Triglycine Exhibit Different
   Reactivities")** — not compared here (that dossier not read for this note); same network, so not a cross-lab check
   [!].
2. **Internal:** "0 %" water = 1.05 % moisture (Methods, Fig. 2, Fig. 4) vs 0.05 % (Fig. 3 caption) [!]; replication
   for non-volatiles stated as duplicate and as three replicates [!]; Fig. 3 bars carry no legend [!].
3. **Reactivity ranking by pK2 fails** (authors: diGly pK2 8.25 and triGly 7.91 leave more amine unprotonated than
   Gly's 9.8 at any pH, yet Gly ≈ triGly > diGly) [Q]: a peptide's N-terminal amine reactivity cannot be set from pKa alone at 160 °C — relevant
   if a peptide amine is ever given a pKa-scaled rate.

## 7. Ledger-ready measurements

- **Table 2** (two readers): peptide breakdown products with glucose at six water contents — the only quantitative
  rows (external standards).
- **Fig. 3** (figure only, legend inferred): peptide-alone controls; pair with Table 2 for the "glycine consumed by
  glucose" difference.
- **Table 1** pyrazine columns at 100 % water (aqueous, the model's domain): relative areas G vs Di vs Tri — a
  ranking signature only.
