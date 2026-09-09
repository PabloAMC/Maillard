# Hofmann & Schieberle 2000 (b) — EXTRACTION (Strecker aldehyde by DIRECT OXIDATIVE degradation of the Amadori compound; argon vs air vs air + Cu²⁺)

**Source on disk:** `data/articles/hofmann2000b.pdf` (owner's download, 2026-08-28). Read-only extraction,
2026-09-07, from `pdftotext -layout`; all three tables' text layers were clean and are re-typed verbatim.
Companion to `hofmann2000_extraction.md` (JAFC 48:434), same group, same SIDA.

**Provenance codes:** **[M]** measured and printed · **[D]** derived by this extraction · **[FIG]** figure only,
no digitised numbers · **[NEG]** verified negative.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Aroma-Active Strecker-Aldehydes by a Direct Oxidative Degradation of Amadori Compounds" |
| Authors | Thomas Hofmann, Peter Schieberle (Deutsche Forschungsanstalt für Lebensmittelchemie, Garching) |
| Venue | J. Agric. Food Chem. 2000, 48 (9), 4301–4305 |
| DOI | 10.1021/jf000076e (received 18 Jan 2000, accepted 16 Jun 2000, web 11 Aug 2000) |
| Funding | DFG Schi 399/6-1 |

## 1. Why it matters

**This is the measurement the task brief asked for: an oxygen-dependent Amadori degradation step with Strecker
aldehyde yields under argon, air and air + Cu²⁺, in one table (Table 2), against the parent glucose/Phe pot under
the same three atmospheres.** The Amadori compound N-(1-deoxy-D-fructos-1-yl)-L-phenylalanine (ARP-Phe) gives
9× more phenylacetaldehyde under air than under argon, and 23× under air + Cu²⁺; the glucose/Phe pot gives
3.5× and 6.5×. The paper's mechanism (Fig. 4): O₂/metal oxidise the open-chain aminoketone of the ARP to an
iminoketone, which either hydrolyses to glucosone + Phe or (as the cyclic hemiketal) decarboxylates to the
Strecker aldehyde **without any free dicarbonyl**. Under argon the ARP instead eliminates Phe and gives
1-deoxyosone (Fig. 2). This is the first paper in the corpus that measures anything the model's "oxidant" pool
could be pinned to on the sugar path, and it also gives the ARP-Phe formation/decay series from glucose/Phe
(Table 3) — a fed-precursor time course for the trunk lane's Amadori node.

## 2. Methods as they matter to a model

| item | value | where |
|---|---|---|
| Precursors | ARP-Phe synthesised (glucose 300 mmol + Phe 400 mmol, MeOH/DMF reflux, malonic acid), **purity ≈ 90 %**, LC/MS 328 [M+1]⁺, ¹H NMR given; or D-glucose + L-Phe | p. 4301 |
| Charge and buffer | Table 2: "the precursors (1 mmol each) were heated (100 °C) in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0) for 120 min in a closed vial" — 0.1 mol/L ARP-Phe, or 0.1 mol/L glucose + 0.1 mol/L Phe. Experimental (p. 4301) states the same concentrations at 1/10 scale (0.1 mmol in 1 mL) and again writes "0.5 mmol/L" phosphate — ⚠ the same mmol/mol inconsistency as the 434 paper; the footnotes' 0.5 mol/L is the likely truth | p. 4301, Table 2 fn |
| Temperature / time / vessel | 100 °C, 120 min, closed vial (Table 2); reflux 30 min (Table 1, AEDA); reflux 10–300 min (Table 3, Fig. 1) | footnotes |
| **Atmosphere (Table 2)** | expt I: "under an atmosphere of argon"; expt II: "Argon was replaced by air oxygen" (text: "Flushing of the models with oxygen" — ⚠ air vs pure O₂ ambiguous; the footnote and expt III's "air atmosphere" point to **air**); expt III: "air atmosphere and in the presence of copper(II) ions (0.05 mmol CuSO₄)" = 5 mmol/L Cu²⁺ | Table 2 fn, p. 4302 |
| ARP-Phe degradation experiment (Figs 2, 3) | ARP-Phe 1.0 mmol + 1,2-diaminobenzene 1.2 mmol (in-situ dicarbonyl trap) in **2 mL** phosphate pH 7.0 ("0.5 mmol/L", sic) = **0.5 mol/L ARP-Phe**, **85 °C**, argon OR air + Cu²⁺ 0.05 mmol (25 mmol/L); aliquots at 20, 40, 80, 160, 320 min; quinoxalines by RP-HPLC | p. 4302 |
| Analytes and method | PA and PAA by SIDA ([¹³C₂] standards, GC-MS/CI, ether extraction at pH 3.0); Phe and ARP-Phe by amino-acid analyser (ninhydrin, 570/440 nm); glucosone, 1-deoxy-2,3-hexodiulose, 3-DG as quinoxalines by HPLC | p. 4302 |
| Quantification basis | **µmol product per mmol precursor** (= 0.1 mol %); Table 3 in absolute µmol per 10 mL / 1 mmol pot | Table 2 header |
| Replicates | **[NEG]** not stated; no uncertainties printed |

## 3. Every table, verbatim

### Table 1 — "Key Odorants Formed by Refluxing Solutions of N-(1-Deoxy-D-fructosyl)-L-phenylalanine (I) or Glucose/L-Phenylalanine (II)"
Footnotes: "a The flavor dilution (FD) factors were determined by the aroma extract dilution analysis. b A solution of N-(1-deoxy-D-fructosyl)-L-phenylalanine (1.0 mmol) in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0) was refluxed for 30 min. c A solution of glucose (1.0 mmol) and L-phenylalanine (1.0 mmol) in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0) was refluxed for 30 min."

| odorant | odor quality | FD factor I (ARP-Phe) | FD factor II (Glc/Phe) |
|---|---|---:|---:|
| phenylacetaldehyde | flowery | 8192 | 1024 |
| phenylacetic acid | honey-like | 1024 | 512 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone | caramel-like | 256 | 64 |

(FD factors are dilution steps, not concentrations: ARP-Phe over Glc/Phe = 8× for PA, 2× for PAA, 4× for HDMF, at
30 min reflux under uncontrolled atmosphere. Column II reproduces Table 1 of the 434 paper exactly.)

### Table 2 — "Concentrations (µmol/mmol) of Phenylacetaldehyde (PA) and Phenylacetic Acid (PAA) Generated upon Thermal Treatment of Either N-(1-Deoxy-D-fructosyl)-L-phenylalanine (ARP-Phe) or D-Glucose/Phenylalanine (Glc/Phe) under Various Conditions" **[M]** — THE OXYGEN TABLE
Footnotes: "a The precursors (1 mmol each) were heated (100 °C) in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0) for 120 min in a closed vial under an atmosphere of argon. b Argon was replaced by air oxygen. c Heating was performed in an air atmosphere and in the presence of copper (II) ions (0.05 mmol CuSO₄)."

| expt | atmosphere | ARP-Phe: PA | ARP-Phe: PAA | Glc/Phe: PA | Glc/Phe: PAA |
|---|---|---:|---:|---:|---:|
| I (a) | argon | 0.6 | 0.3 | 0.4 | 0.2 |
| II (b) | air | 5.5 | 3.0 | 1.4 | 1.8 |
| III (c) | air + Cu²⁺ 5 mmol/L | 13.8 | 7.6 | 2.6 | 3.9 |

Units µmol/mmol precursor (÷10 for mol %). **[D]** derived:

| quantity | ARP-Phe | Glc/Phe |
|---|---:|---:|
| PA, mol % (argon / air / air+Cu) | 0.06 / 0.55 / 1.38 | 0.04 / 0.14 / 0.26 |
| PAA, mol % (argon / air / air+Cu) | 0.03 / 0.30 / 0.76 | 0.02 / 0.18 / 0.39 |
| PA air / argon | **9.2** | **3.5** |
| PAA air / argon | **10** | **9** |
| PA (air+Cu) / air | 2.5 | 1.9 |
| PAA (air+Cu) / air | 2.5 | 2.2 |
| PA (air+Cu) / argon | 23 | 6.5 |
| PA / PAA (argon / air / air+Cu) | 2.0 / 1.8 / 1.8 | 2.0 / 0.78 / 0.67 |
| ARP-Phe over Glc/Phe, PA (argon / air / air+Cu) | 1.5 / 3.9 / 5.3 | |
| ARP-Phe over Glc/Phe, PAA (argon / air / air+Cu) | 1.5 / 1.7 / 1.9 | |

### Table 3 — "Formation of N-(1-Deoxy-D-fructosyl)-L-phenylalanine (ARP-Phe) from Glucose and L-Phenylalanine" **[M]**
Footnotes: "a A solution of glucose (1.0 mmol) and L-phenylalanine (1.0 mmol) was refluxed in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0). b The yield of ARP-Phe was calculated based on the amounts of L-phenylalanine reacted." Atmosphere not controlled (air, reflux).

| reaction time (min) | ARP-Phe (µmol) | phenylalanine unreacted (µmol) | yield of ARP-Phe (%) |
|---:|---:|---:|---:|
| 10 | 42 | 920 | 52.5 |
| 30 | 48 | 825 | 27.4 |
| 60 | 32 | 640 | 9.0 |
| 120 | 21 | 510 | 4.3 |
| 300 | 9 | 398 | 1.5 |

**[D]** The yield column reproduces exactly as ARP/(1000 − Phe unreacted): 42/80, 48/175, 32/360, 21/490, 9/602.
Phe consumed: 8.0 / 17.5 / 36.0 / 49.0 / 60.2 %. ARP-Phe peaks at 30 min at **4.8 mol % of glucose** and falls to
0.9 mol % by 300 min. Text: "After 300 min, about 60 % of the L-phenylalanine was degraded, but only 1.5 % of the
ARP-Phe was left" (sic — 1.5 % is the yield-on-Phe-reacted column). A first-order decay fitted to the 30 → 300 min
fall of ARP-Phe alone (48 → 9 µmol, 270 min) gives **k ≈ 0.0062 /min at reflux** — but the ARP is still being
formed over that window, so this is a LOWER bound on the ARP-Phe decay constant. For comparison, the trunk lane's
Martins 2005 step 4 + 6 + 7 (DFG-glycine → 3-DG, MG, 1-DG) sum to 0.0339 /min at 100 °C.

### Figures (no numbers in the text layer) **[FIG]**
- Fig. 1: PA and PAA vs time (0–500 min, reflux, air) from Glc/Phe and from ARP-Phe. Text: ARP-Phe gives PA > PAA
  at every time; Glc/Phe gives PA > PAA only for the first 60 min and PAA/PA = **1.7 at 500 min** (cf. "nearly twice"
  after ~8 h in the 434 paper, Fig. 2).
- Fig. 2 (ARP-Phe 0.5 M, 85 °C, argon, OPD trap): 1-deoxyosone the main dicarbonyl, "in particular after longer
  reaction times"; glucosone "very low"; Phe liberated **≈ 55 % after 320 min**.
- Fig. 3 (same, air + Cu²⁺ 25 mmol/L): **glucosone the main dicarbonyl**; 1-DG and 3-DG "only low yields"; Phe
  liberation "in the same order of magnitude" as under argon.

## 4. What the repo could take

Under the owner's rule (within-study ratios FIT; levels VALIDATE). Names follow `hofmann2000_extraction.md`.

**Fed-intermediate yields (levels — VALIDATE rows; 0.1 M precursor, 0.5 M phosphate pH 7.0, 100 °C, 120 min, closed vial):**

| row | precursor | atmosphere | product | mol % of precursor |
|---|---|---|---|---:|
| H00b-L1 | ARP-Phe | argon | PA | 0.06 |
| H00b-L2 | ARP-Phe | air | PA | **0.55** |
| H00b-L3 | ARP-Phe | air + Cu²⁺ 5 mM | PA | 1.38 |
| H00b-L4 | ARP-Phe | argon / air / air+Cu | PAA | 0.03 / 0.30 / 0.76 |
| H00b-L5 | glucose + Phe | argon / air / air+Cu | PA | 0.04 / 0.14 / 0.26 |
| H00b-L6 | glucose + Phe | argon / air / air+Cu | PAA | 0.02 / 0.18 / 0.39 |
| H00b-L7 | glucose + Phe, reflux, air | — | ARP-Phe at 10/30/60/120/300 min | 4.2 / 4.8 / 3.2 / 2.1 / 0.9 mol % of glucose |
| H00b-L8 | glucose + Phe, reflux, air | — | Phe consumed at 10/30/60/120/300 min | 8.0 / 17.5 / 36.0 / 49.0 / 60.2 % |

**Within-study ratios (FIT candidates):**

| row | ratio | value | note |
|---|---|---:|---|
| H00b-R1 | PA(air)/PA(argon) from ARP-Phe | **9.2** | the oxygen dependence of Amadori → Strecker aldehyde; the first such number in the corpus |
| H00b-R2 | PAA(air)/PAA(argon) from ARP-Phe | 10 | |
| H00b-R3 | PA(air)/PA(argon) from glucose + Phe | 3.5 | the parent pot is oxygen-sensitive too, less so |
| H00b-R4 | PAA(air)/PAA(argon) from glucose + Phe | 9 | matches the 434 paper's 4–5.5× on dicarbonyl donors in direction |
| H00b-R5 | PA(air + Cu)/PA(air), ARP-Phe / Glc-Phe | 2.5 / 1.9 | metal catalysis on top of oxygen |
| H00b-R6 | PA/PAA from ARP-Phe, any atmosphere | 1.8–2.0 | the ARP route's branch ratio is atmosphere-INdependent |
| H00b-R7 | PA/PAA from Glc/Phe, argon → air → air+Cu | 2.0 → 0.78 → 0.67 | the glucose pot's branch ratio flips with oxygen |
| H00b-R8 | PA from ARP-Phe over PA from Glc/Phe, air | 3.9 | fed Amadori vs parent pot at equal charge |
| H00b-R9 | ARP-Phe(30 min)/ARP-Phe(300 min), reflux | 5.3 | the shape of the Amadori transient |
| H00b-R10 | ARP-Phe yield on Phe reacted, 10 → 300 min | 52.5 → 1.5 % | the Amadori is the first sink of the amino acid, then not |
| H00b-R11 | FD(PA) ARP-Phe / Glc-Phe, 30 min | 8 | dilution steps, directional only |

**Directional claims a model must reproduce:**
1. **The Amadori compound has an oxygen-dependent route to the Strecker aldehyde that needs no free dicarbonyl**
   (Fig. 4; R1). In the model: an ARP → Strecker-aldehyde step whose rate scales with the oxidant pool. Under argon
   the ARP goes to 1-deoxyosone + amino acid instead (Fig. 2) — the Martins 2005 step 7 topology.
2. **Under air + metal the ARP's main dicarbonyl is glucosone, not 1-DG/3-DG** (Fig. 3 vs Fig. 2). A model whose
   ARP → dicarbonyl split does not respond to oxygen is refuted at the direction level.
3. Amino-acid liberation from the ARP (~55 % at 320 min, 85 °C, 0.5 M) is roughly the same with or without oxygen;
   oxygen changes WHERE the carbon goes, not how fast the amine returns.
4. In glucose/Phe, ARP-Phe holds > 50 % of the amino acid reacted at 10 min but only 1.5 % at 300 min (R10): the
   Amadori is a transient intermediate whose peak (4.8 mol % of glucose at 30 min reflux) is a level the trunk lane
   can be scored on.
5. The ARP route gives aldehyde over acid ≈ 2:1 regardless of atmosphere (R6); the dicarbonyl route (434 paper)
   gives acid over aldehyde under air. The two Strecker routes are distinguishable by their PA/PAA fingerprint.

## 5. Caveats

- **No kinetics.** Single temperature per experiment (100 °C for Table 2; reflux for Table 3; 85 °C for Figs 2–3);
  no rate constants or Ea printed. The ARP decay bound in §3 (Table 3) is this extraction's, not the authors'.
- **Oxygen charge is undefined.** "Argon replaced by air oxygen" in a closed 10 mL vial with an unstated headspace;
  neither the dissolved-O₂ nor the headspace volume is given, so R1–R4 are ratios between "no O₂" and "air-saturated
  closed vial", not a dose response. Whether expt II was flushed with pure O₂ (text) or air (footnote) is ambiguous.
- **ARP-Phe purity ≈ 90 %** (impurity unidentified; possibly Phe/glucose). The argon level from ARP-Phe (0.06 mol %)
  is close to the glucose/Phe argon level (0.04), so a 10 % contamination could explain part of the argon floor; the
  air and Cu²⁺ numbers are far above it.
- **Cu²⁺ 5 mmol/L (Table 2) and 25 mmol/L (Figs 2–3)** are single high doses; treat R5 as a direction.
- **No replicates, no SDs.**
- **Buffer "0.5 mmol/L" vs "0.5 mol/L"** typo as in the companion; the Fig. 2–3 experiment is at 0.5 M ARP-Phe in
  2 mL, five times the Table 2 concentration and at 85 °C, so its dicarbonyl pattern is not the Table 2 pot's.
- **Table 3's Phe-consumption series coincides with Table 3 of the 434 paper** (MGO/Phe, closed vial 98 °C): 8.0 /
  17.5 / 36.0 / 49.0 / 60.2 % here vs 9 / 18 / 36 / 49 / 61 % there. One series may have been reused; do not count
  the two as independent.
- Table 2 expt I–III are 120 min endpoints; Fig. 1 (the time course) is under uncontrolled atmosphere and not
  digitised, so the argon-vs-air comparison exists at one time only.
- The glucose/Phe rows of Table 2 (0.14 / 0.18 mol % at 120 min, air) are consistent with the 434 paper's 30-min
  values (0.042 / 0.029) and its Fig. 2 shape (PAA overtakes PA after 60 min); the two papers can be read as one pot.
