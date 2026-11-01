# Ayllón-Parra, del Castillo, Matas, Gou & Ribas-Agustí 2027 (online 2026) — EXTRACTION
DOI: 10.1016/j.ifset.2026.104870 (printed on p. 1; Crossref 2026-11-01: "Maillard reaction in fava bean protein and oat flour blends under solid-state fermentation and low-moisture extrusion", Ayllón-Parra, del Castillo, Matas, Gou, Ribas-Agustí, Innov. Food Sci. Emerg. Technol. 115, article 104870, issue dated 2027-01; online 2026-10-27. **Matches the brief**; the file stem says 2026, the volume says 2027.)
### Fava bean protein concentrate (62 % protein DM) + oat flour (17 %) at 50-50 and 80-20, two routes: (A) low-moisture twin-screw extrusion only (barrel up to 110 °C, < 1 min, then 70 °C oven 1.5 h); (B) hydrate to 60 % moisture, autoclave 121 °C 15 min, solid-state fermentation with *Pleurotus ostreatus* (28 °C, 14 d), oven 70 °C 2 h, then the same extrusion. Markers per stage: OPA free amino groups, NBT fructosamine, fluorescent AGEs, HMW (> 10 kDa) melanoidin index, colour; lysine and CML in the final extrudates only. **No a_w, no time course, no rate; the fungus confounds lysine; CML by acid hydrolysis without NaBH4.**

**Read 2026-11-01 from `data/articles/aylonparra2026.pdf`** (11 pp., born-digital). Prose and Tables 1-7 from
`pdftotext -layout`; **Fig. 1 (process table) read from a 130 dpi render of p. 3** (scratchpad `ing5/ap-03.png`).
Table S1 (sugar composition of the flours) is not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** IRTA Monells, Food Quality and Technology / Food Safety and Functionality (Gou; Ribas-Agustí), with
CIAL (CSIC-UAM) Madrid (del Castillo). **Not in `lab_ids.yml`: NEW** (proposed lead id `irta_gou`, partner CIAL
under `collaborations:`).

## 1. What this paper answers for the model

| question | answer | tag |
|---|---|---|
| which Maillard markers | free amino groups (OPA, Nα-acetyl-Lys eq. per mg soluble protein), fructosamine (NBT, 1-deoxy-1-morpholino-fructose eq.; LOD 50 uM), fluorescent AGEs (ex 360 / em 460, quinine eq.), HMW melanoidin index (> 10 kDa, A405, caramel E-150d eq.), L*a*b*, lysine and CML (6 N HCl 110 °C 24 h, ion-exchange AAA, ninhydrin 570 nm) | [M] |
| conditions | sec. 2; extrusion 110 °C (zones 40 / 60 / 80 / 100 / 110 x3 / die 110), < 1 min, 600 rpm, feed 3.32 g/min + water 1.67 g/min; autoclave 121 °C 15 min at ~62 % moisture; SSF 28 °C 14 d; ovens 70 °C | [M] |
| a_w | **not measured** (moisture only: 8.03 % raw, 11.41 % after extrusion A, 16.31 % after extrusion B) | [M] |
| pH | 1 : 10 water slurry at 30 °C, after each stage: 6.81 -> 6.85 (A); 6.81 -> 6.71 -> **7.95** (after SSF, ammonia [P]) -> 7.28 -> 7.16 (B) | [M] |
| a rate or a kinetic constant | **none** (stage end points, n = 3 batches) | — |
| final lysine / CML | lysine **71.15** (A) vs **59.29** (B) mg/g total amino acids (-17 %); CML **22.75** vs **49.27** ug/g DM (+117 %) | [M] |

**Answer.** Maillard markers are measured at stage end points of two food processes, not under controlled reaction
conditions. The model can use at most two plant-protein extrudate end points (CML 22.75 ug/g DM after a 110 °C,
< 1 min extrusion of an unfermented fava / oat blend; fructosamine below the NBT LOD in the same samples). Confidence
that anything here can enter a fit: ~0.05.

## 2. Conditions

- **Ingredients** [M]: fava bean protein concentrate (AGT Foods; 62 % protein DM), oat flour (Fazer; 17 % protein DM);
  blends 50-50 (40.43 % protein) and 80-20 (53.94 %). Sugars in Table S1 (not on disk).
- **Process A** [M]: extrusion of the blends (moisture 8.11 / 7.48 %) on a Process 11 twin-screw (11 mm, L/D 40),
  then 70 °C air oven 1.5 h.
- **Process B** [M]: 200 g hydrated to 60 % moisture (w.b.), gas-permeable lid, autoclave 121 °C 15 min; cubes
  inoculated with *P. ostreatus* M2191 (4-5 mg biomass), 28 °C 14 d; ground, 70 °C 2 h; then extrusion as A and
  70 °C 1.5 h.
- **Extracts** [M]: 500 mg in 20 mL 100 mM carbonate pH 10.8, 37 °C 30 min; soluble protein by Bradford; free amino
  groups, fructosamine and fluorescent AGEs normalised per mg SOLUBLE protein (extractability changes between stages
  acknowledged [Q]).
- **O2**: air throughout; not varied.

## 3. Quantitative results

**Process A** (Table 1, after mixing -> after extrusion; pooled over blends) [M]: moisture 8.03 -> 11.41 %; pH 6.81 ->
6.85; free amino groups 1.20 -> 3.76 umol/mg sol. protein; fructosamine ND -> ND; fluorescent AGEs 5.23 -> 48.53
nmol/mg; HMW melanoidins 7.21 -> 12.02 mg/g DM; L* 92.51 -> 82.09.

**Process B** (Tables 2-6) [M]: sterilisation — moisture 61.85 %, pH 6.71, free amino groups 1.15 -> 4.15 (50-50) and
1.23 -> 2.63 (80-20), fluorescent AGEs 4.47 -> 22.58 and 5.98 -> 37.12, HMW melanoidins ND; SSF — pH 7.95, free
amino groups 3.39 -> 29.05, fluorescent AGEs 29.85 -> 297.83, HMW melanoidins 21.33; drying — moisture 7.80 %, pH
7.28, fructosamine **6.91** umol DMF/mg, fluorescent AGEs 633.99, HMW melanoidins 28.79; extrusion — moisture 16.31 %,
pH 7.16, free amino groups 38.58 -> 22.79, fructosamine 6.92 -> 1.91, fluorescent AGEs 633.99 -> 285.53, HMW
melanoidins 29.30 -> 22.20 (50-50) and 28.29 -> 41.30 (80-20).

**Table 7** (final extrudates) [M]: lysine 71.15 +- 1.31 (A) vs 59.29 +- 0.49 (B) mg/g total AA; 50-50 63.55, 80-20
66.89; CML 22.75 +- 6.88 (A) vs 49.27 +- 6.22 (B) ug/g DM; 50-50 33.54, 80-20 38.47 (n.s.).

## 4. Mechanism evidence and its kind

Stage-wise end points with two-way ANOVA (stage x blend). No labelling, no mass balance, no controls separating
fungal metabolism from Maillard chemistry: the lysine loss in B "cannot be exclusively attributed" to the Maillard
reaction [V]. The fructosamine peak after the 70 °C drying of fermented material and its fall on extrusion is the
one sequential early -> late signal [M]. Laboratory: IRTA Monells (NEW).

## 5. Bearing on open questions

- **Atlas sec. 1.4 "plant-protein-bound amines"** — LEAVES essentially UNCHANGED: one CML end point for an
  unfermented fava / oat extrudate (22.75 ug/g DM), no rate, no a_w, no sugar data on disk.
- **Atlas sec. 1.4 "alkaline in-situ pH at cooking T"** — LEAVES UNCHANGED: pH 7.95 is after a 28 °C fermentation;
  the hot steps start from pH 6.7-7.3 (slurry pH, not in-situ).
- **Atlas sec. 1.4 "intermediate a_w with a measured a_w"** — LEAVES UNCHANGED (no a_w).
- **F05 (CML route; CF-072)** — LEAVES UNCHANGED; CML by acid hydrolysis without NaBH4 reduction [!] can include
  fructosyl-lysine converted to CML in hot HCl (the F05 marker note), and the B samples carry measurable fructosamine
  before extrusion, so part of B's +117 % CML may be analytical.

## 6. Disagreements found

None with the engine or another dossier (no comparable quantity). Internal [!]: the text gives post-drying moisture
7.20 %, Table 5 7.80 %; Table 2's "7 0.21" HMW melanoidin value after mixing is a typo for 7.21 (Table 1); CML by a
ninhydrin amino-acid analyser is unusual (LC-MS/MS is the norm in the F05 corpus), with no recovery or LOD printed.

## 7. Ledger-ready measurements

- Table 7 CML and lysine of extrudate A (unfermented fava / oat, 110 °C, < 1 min, moisture 8 -> 11 %): a plant-protein
  extrudate end point; low priority.
- Nothing else is keyable as kinetics.
