# Yang et al. 2024 — EXTRACTION (nine free fatty acids heated 80 C / 30 min in sealed vials; SPME-GC-Orbitrap volatiles; drawn hydroperoxide-isomer -> volatile maps for 6Z-/9Z-octadecenoic, linoleic, CLA, arachidonic and n-3 acids)
### An experimental "simulation" (model system), not a computation: the product maps are drawn hypotheses; the measured content is a presence/intensity heat-map (FIGURE-ONLY).

**Source on disk:** `data/articles/Yang2024.pdf` (owner's download, 2026-09-08; open access, LWT 214
(2024) 117083). No `yang2024_extraction.md` existed before this one. Read from the `pypdf` text layer;
the paper has NO numeric table. Fig. 1 (bar charts of log10 peak area) is FIGURE-ONLY. Figs 2-7 are
ChemDraw reaction maps, described in words from renderings of pages 3-5 and 8. Supplementary Figs
S1-S4 (fatty acid structures, CLA and n-3 maps) not on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Simulating fatty acid autoxidation and exploring the related volatiles formation mechanism" |
| Authors | Youyou Yang, Dapeng Liu, Weihai Xing, Chaohua Tang, Xiaohui Feng, Junmin Zhang (Institute of Animal Sciences, CAAS, Beijing) |
| Venue | LWT - Food Science and Technology 214 (2024) 117083 (PII S0023643824013665) |
| DOI | 10.1016/j.lwt.2024.117083 |
| What "simulating" means here | a bench model system: each pure free fatty acid heated in a sealed vial "which was the same as the heating procedure for stewed meat". There is no kinetic model, no quantum chemistry, no fitted parameters. |
| Naming | "HPOD" = hydroperoxide (their abbreviation, any fatty acid); "hyperoxide" in the text = hydroperoxide; CLA = (9Z,11E)-conjugated linoleic acid; LA = (9Z,12Z)-linoleic acid; ALA/EPA/DHA n-3; AA arachidonic; "n-9 FAs" = the two octadecenoic acids |

## 1. Why it matters

The paper heats nine pure free fatty acids under one mild condition and draws, for each, the
hydroperoxide isomers and the volatiles they should give. For the repo's three refused routes it
contributes: (i) an oleate (9Z-octadecenoic) map listing 8-, 9-, 10-, 11-HPOD -> alkanals heptanal to
decanal, with a positional-isomer control (6Z-octadecenoic -> 5-/6-/7-/8-HPOD -> nonanal to dodecanal)
that shows the alkanal chain length tracks the double-bond position; (ii) an LA map that assigns
2-pentylfuran and 2-butylfuran to 9- and 8-HPOD, 2-octenal / 2-heptenal to 11-/12-HPOD, and 1-octen-3-ol
to the 10-hydroperoxide via 2-octen-1-ol rearrangement; (iii) the statement (conclusion) that linoleic
acid and glycerol trilinoleate gave the same volatile pattern and intensities. Everything mechanistic
is a drawn proposal; the measured part is which volatiles appeared from which acid, shown only as bar
heights in Fig. 1.

## 2. Methods as they matter to a model

### 2a. What was MEASURED

- **Substrates (Sigma-Aldrich, purity not stated):** stearic acid (control), (9Z)-octadecenoic acid,
  (6Z)-octadecenoic acid, (9Z,11E)-CLA, (9Z,12Z)-LA, (9Z,12Z,15Z)-ALA, (5Z,8Z,11Z,14Z)-AA,
  (5Z,8Z,11Z,14Z,17Z)-EPA, (4Z,7Z,10Z,13Z,16Z,19Z)-DHA. GLA is mentioned in the text ("LA, GLA, and AA
  represented the n-6 FAs") but is not in the chemicals list. Glycerol trilinoleate is mentioned only in
  the conclusion. Free acids.
- **Charge and heating:** 10 mg fatty acid in a 20 mL glass vial, magnetic cap with PTFE-silicone
  septum (i.e. sealed, ambient air headspace); **80 C for 30 min**. No water, no catalyst, no light
  control stated, no replicate count stated.
- **SPME:** vial incubated 55 C / 20 min, extracted 55 C / 40 min, DVB/CAR/PDMS 50/30 µm; desorbed
  250 C / 3 min. GC Trace 1310 + Q-Exactive Orbitrap, VF-WAXms 60 m x 0.25 mm x 0.25 µm, He 1 mL/min,
  40 C (2 min) -> 230 C at 4 C/min (5 min); EI 70 eV, full scan 30-400 at 60 000 FWHM.
- **Identification:** NIST17, Wiley9 and a home library; HRF score > 95, match factor > 750, RI
  difference < 20 (home library) or < 50 (NIST) vs C7-C40 alkanes. Authentic standards were run for:
  heptanal, octanal, 2-heptanone, 2-octanone, 3-octanone, pentanal, 1-octen-3-ol, 1-octen-3-one,
  (2E,4E)-nonadienal, (2E,4E)-decadienal, 2-octenal, (2E,4E)-heptadienal, (2E)-nonenal, (2E)-heptenal,
  hexanal, nonanal, decanal, (2E)-undecenal, (2E)-decenal, dodecanal, 2-furanaldehyde, 2-pentylfuran,
  2-hexylfuran, 2-butylfuran. Furaldehyde identities were "verified by authentication standards".
- **Quantification: none.** Fig. 1 plots log10(A) (A = peak area, stacked per acid) for selected
  volatiles; no table, no internal standard, no units beyond area. FIGURE-ONLY.
- **Hydroperoxide isomers:** NOT measured. Every HPOD in Figs 2-7 is drawn from the radical-abstraction
  / resonance rule, not detected.

### 2b. What is DRAWN / PROPOSED (not measured)

Figs 2-7 and S2-S4: radical abstraction at the allylic / bis-allylic CH2, resonance, O2 addition to
give the HPOD set, O-O homolysis to the alkoxyl, "beta-scission reaction" / "homolytic cleavage" to
aldehydes and radicals, then radical + •OH / + H / + O2 follow-ups, retro-aldol of dienals, and
cyclisations to furans. No rate, no energy, no yield is attached to any arrow.

## 3. Tables re-typed

There are no numeric tables. The schemes:

### Fig. 1 — "The key volatiles derived from MUFA and PUFA" (FIGURE-ONLY)

Eight bar panels of log10(A): (A) alkanals/alkenals from (6Z)- vs (9Z)-octadecenoic acid — heptanal,
octanal, nonanal, decanal, undecanal, dodecanal, (2E)-heptenal, (2E)-octenal, (2E)-nonenal,
(2E)-decenal, (2E)-undecenal, 2,4-decadienal; (B) 2-pentylfuran, 2-hexylfuran, 2-octanone from the two
MUFAs; (C) C4-C10 alkanals/alkenals/dienals from CLA vs LA; (D) 2-butylfuran, 2-pentylfuran,
2-hexylfuran, furan-2-carbaldehyde from CLA vs LA; (E) ketones and 1-octen-3-ol from CLA vs LA; (F)-(H)
the same classes from AA, ALA, EPA, DHA. Bar heights are not read here. Presence/absence statements
from the text are recorded in §4.

### Fig. 2 — "The mechanism of autoxidation of cis-6 octadecenoic acid" (drawn; described)

(6Z)-octadecenoic acid -> H• abstraction at C5 or C8 + resonance -> HPOD isomers **5-, 6-, 7-, 8-HPOD**
-> alkoxyl radicals -> "homolytic cleavage" -> (left branch) undecyl-type radicals -> + •OH -> alcohol
-> dodecanal; a C13 radical chain -> (E)-3-tridecenal / (E)-3-dodecenal -> (E)-2-tridecenal /
(E)-2-dodecenal; decanal, undecanal and 2-decenal at the bottom; (right branch) alkoxyl -> (2E)-undecenal,
(2E)-dodecenal, nonanal, decanal. Text: "When the double bond was more adjacent to the carboxylic acid
group, the carbon length of the produced aldehydes was elongated from nonanal to dodecanal."

### Fig. 3 — "The mechanism of autoxidation of cis-9 octadecenoic acid" (drawn; described)

(9Z)-octadecenoic acid -> H• abstraction (C8 / C11) + resonance -> **8-HPOD (C9=C10), 9-HPOD (C10=C11),
10-HPOD (C8=C9), 11-HPOD (C9=C10)** (all four drawn with the OOH carbon and the shifted C=C) -> "- •OH,
beta-scission reaction" -> four alkoxyl radicals -> "homolytic cleavage" -> three product columns:
- left: a C9/C10 alkyl radical chain -> + •OH -> secondary/primary alcohols -> **decanal, nonanal,
  octan-2-one** (drawn via the alcohol, i.e. alcohol -> aldehyde/ketone oxidation);
- middle: a C8/C9 radical -> + •OOH -> hydroperoxide -> a conjugated dienyl-alkoxyl -> **2-hexylfuran**
  and **2-pentylfuran** (drawn from oleic acid; see flag 1);
- right: **2-undecenal, nonanal, 2-decenal, octanal**; then **2,4-decadienal** and **2,4-undecadienal**
  drawn below 2-decenal / 2-undecenal, and further down **2-nonenal, 2-octenal, 2,4-nonadienal,
  2-heptenal, hexanal, heptanal** (drawn as the retro-aldol / further-oxidation products of the dienals).
The figure does NOT connect a specific HPOD isomer to a specific aldehyde by arrows; the text says only:
"hyperoxide isomers including 8-, 9-, 10- and 11-hydroperoxide (HPOD) were produced, and further
decomposed into alkanals from heptanal to decanal."

### Fig. 4 — "The mechanism of autoxidation of LA" (drawn; described)

(9Z,12Z)-octadecadienoic acid -> H• abstraction + resonance -> + O2, H -> HPOD isomers drawn as **8-,
9-, 10-, 11-, 12-HPOD** plus two unlabelled structures at right (one with OOH at the 13 position, one
with OOH nearer the methyl end); text: "the types of LA oxidation extended to 9-, 10-, 11-, and
12-HPOD"; furans "through decomposition of 8- and 9-HPOD". -> alkoxyl radicals (three drawn, with the
scission bonds marked) -> "homolytic cleavage" -> three boxes:
- left box: a C10 alkyl radical -> + •OH -> 3-decenol / 3-nonenol -> **3-decenal / 3-nonenal** ->
  **2-decenal / 2-nonenal** -> **2,4-decadienal / 2,4-nonadienal** -> **hexanal, pentanal, 2-octenal,
  2-heptenal** (drawn as the fragmentation of the dienals); in the same box a C9 radical -> + O2, H ->
  hydroperoxide -> - •OH -> alkoxyl -> cyclisation -> **2-pentylfuran**; the "retro-aldol" arrow from
  2,4-decadienal -> **2-octenal + hexanal** is drawn between the left and middle boxes;
- middle box: a C8 radical -> + O2, H -> hydroperoxide -> - •OH -> alkoxyl -> **octanal** / **hexanal**
  and cyclisation -> **2-butylfuran**;
- right box: a short dialkoxyl fragment -> **furan-2-carbaldehyde**.
Text assignments: "2-octenal and 2-heptenal were produced through their 11- and 12-hydroperoxides,
respectively"; "hexanal and heptanal were derived from the decomposition of HPOD containing one C=C
bond"; "Autoxidation of LAs produced butyl furan and pentyl furan through decomposition of 8- and
9-HPOD"; "1-octen-3-ol was the main alcohol, which originated from the 10-hydroperoxides of the n-6
fatty acids and was generated from the rearrangement of 2-octen-1-ol" (this last step is not drawn in
Fig. 4).

### Fig. 5 — arachidonic acid; Fig. 6 — n-3 PUFA -> furan derivatives (3-furaldehyde, furfural,
5-methylfurfural, 5-ethyl-2-furaldehyde, 5-ethyl-2(5H)-furanone, 1-(2-furyl)propan-1-one);
Fig. 7 — (A) ALA -> 3,5-octadien-2-one; (B) CLA 13-hydroperoxide -> "LH" -> hexanal -> hexanoic acid;
(C) CLA 9-hydroperoxide (via rearrangement) -> "LH" -> octanoic acid. Not transcribed further; none
touches the three refused routes except 7B (13-OOH of CLA -> hexanal, drawn as a single arrow).

## 4. Routes and numbers the repository can use

No numbers are printed. "Present" below means a bar exists in Fig. 1 or the text says the compound was
detected from that acid at 80 C / 30 min.

| route | reactant -> product | mechanism as drawn (figure) | measured (presence only; intensities FIGURE-ONLY) | evidence class |
|---|---|---|---|---|
| OL-ALK | (9Z)-octadecenoic acid -> heptanal, octanal, nonanal, decanal (via 8-/9-/10-/11-HPOD) | Fig. 3, no isomer-to-product arrows | present (Fig. 1A); text: "alkanals containing the carbon number of 6-12 were the major products for MUFA oxidation" | mechanism_drawn (loose); figure_only |
| OL-POS | (6Z)-octadecenoic acid -> nonanal, decanal, undecanal, dodecanal (via 5-/6-/7-/8-HPOD) | Fig. 2 | present (Fig. 1A: undecanal, dodecanal appear only from the 6Z acid) | mechanism_drawn; figure_only; the positional control that supports "alkanal length tracks C=C position" |
| OL-ENAL | (9Z)-octadecenoic acid -> (2E)-decenal, (2E)-undecenal, (also (2E)-heptenal ... (2E)-nonenal, 2,4-decadienal drawn) | Fig. 3 right column | bars present in Fig. 1A for the 9Z acid | figure_only; see flag 1 |
| OL-PF? | (9Z)-octadecenoic acid -> 2-pentylfuran, 2-hexylfuran, 2-octanone | Fig. 3 middle column | bars present (Fig. 1B) | figure_only; chemically requires a diene: likely impurity (flag 1) |
| LA-PF | LA -> 2-pentylfuran (from 9-HPOD) and 2-butylfuran (from 8-HPOD) | Fig. 4: C9 / C8 alkyl radical -> O2 -> OOH -> alkoxyl -> cyclisation | present (Fig. 1D: 2-pentylfuran from both CLA and LA, 2-butylfuran from LA only, 2-hexylfuran from CLA only) | mechanism_drawn; figure_only |
| LA-HEX | LA -> hexanal, pentanal, 2-octenal, 2-heptenal | drawn as fragmentation of 2,4-decadienal / 2,4-nonadienal (left box) and retro-aldol 2,4-decadienal -> 2-octenal + hexanal; text also: 2-octenal from 11-HPOD, 2-heptenal from 12-HPOD, hexanal/heptanal "from HPOD containing one C=C bond" | present (Fig. 1C); "Butanal, pentanal, and heptanal, as well as 2,4-nonadienal can be detected during oxidation of LA oxidation, but not CLA" | mechanism_drawn (two inconsistent routes for hexanal: retro-aldol vs direct); figure_only |
| LA-DIEN | LA -> 2,4-decadienal (major), 2,4-nonadienal (minor) | Fig. 4 | text: 2,4-nonadienal "intensity was much less than that of 2,4-decadienal" | measured_ratio (qualitative, FIGURE-ONLY) |
| LA-OCTENOL | **LA 10-hydroperoxide -> 2-octen-1-ol -> 1-octen-3-ol (rearrangement)** | text only (§3.5); not drawn | 1-octen-3-ol present from LA (Fig. 1E) and from n-3/n-6 (Fig. 1H) | mechanism stated; figure_only |
| LA-TAG | glycerol trilinoleate vs linoleic acid | none | conclusion: "there is no difference in the pattern of volatiles and their intensities between linoleic acid and glycerol trilinoleate" (data not in the main text) | claim only |
| SFA-NULL | stearic acid -> nothing | none | "no volatiles were generated through oxidation of SFAs such as stearic acid" | measured_level (null) — the natural negative control for every rule |
| CLA-13 | CLA 13-OOH -> hexanal -> hexanoic acid; CLA 9-OOH -> octanoic acid | Fig. 7B, 7C, single arrows "LH" | present | mechanism_drawn |
| FURALD | n-3 / n-6 PUFA with > 2 C=C -> furfural, 3-furaldehyde, 5-methylfurfural, 5-ethyl-2-furaldehyde | Fig. 6 | present, confirmed with standards; not from LA/CLA/MUFA | mechanism_drawn; figure_only (outside the refused routes) |

## 5. Rule sketches (repository suggestions, not the paper's)

Keys: `hexanal`, `nonanal`, `heptanal`, `2_pentylfuran`, `1_octen_3_ol`, `e_2_octenal`,
`e_e_2_4_heptadienal`, `furan`, `2_acetylfuran`, `furfural` exist; **octanal, decanal, undecanal,
dodecanal, 2-butylfuran, 2-hexylfuran, (E)-2-decenal, (E)-2-undecenal, (E,E)-2,4-decadienal,
2-octen-1-ol, 5-methylfurfural-as-lipid-product** do not (`furfural` exists as a Maillard key; this
paper would add a lipid provenance to it).
Structures: no free-acid linoleate/oleate, no 10-HpODE; see miyazaki2023 §5 and cao2020 §5 for the
isomer SMILES.

This paper should be cited as a SUPPORTING witness (presence of products from a pure acid), not as the
mechanism source, because its arrows are loose and partly inconsistent. Specific uses:

**S1. Positional control for the oleate alkanal rules (cao2020 S1).** Add (6Z)-octadecenoic acid as an
extra positive control for the mono-ene B-scission rule (alkanal on the methyl side): the paper draws
its hydroperoxides as 5-/6-/7-/8-HPOD and reports undecanal and dodecanal only from this acid. With
Cao's geometry, the 8-OOH (C6=C7) B-scission gives undecanal (C8-C18) and the 7-OOH (C5=C6) gives
dodecanal (C7-C18); alkanal chain length = 19 - (OOH carbon number). A SMIRKS that hard-codes nonanal
would fail this control.
- positive: 7-hydroperoxy-5E-octadecenoic acid `CCCCCCCCCCCC(OO)/C=C/CCCC(=O)O` -> dodecanal `CCCCCCCCCCCC=O` + 7-oxoheptanoic acid `O=CCCCCCC(=O)O`
- positive: 8-hydroperoxy-6E-octadecenoic acid `CCCCCCCCCCC(OO)/C=C/CCCCC(=O)O` -> undecanal `CCCCCCCCCCC=O` + 8-oxooctanoic acid `O=CCCCCCCC(=O)O`
- negative control: stearic acid `CCCCCCCCCCCCCCCCCC(=O)O` (measured null in this paper).

**S2. 2-Octen-1-ol -> 1-octen-3-ol (allylic 1,3-transposition; new, net).**
`[CH3:8][CH2:7][CH2:6][CH2:5][CH2:4][CH1:3]=[CH1:2][CH2:1][OH:9]` -> `[CH3:8][CH2:7][CH2:6][CH2:5][CH2:4][CH1:3]([OH:9])[CH1:2]=[CH2:1]`
- positive: `CCCCC/C=C/CO` (2-octen-1-ol) -> `C=CC(O)CCCCC` (1-octen-3-ol, key `1_octen_3_ol`)
- negative: 1-octanol `CCCCCCCCO` (no allylic system); hexanal.
- caveat: Miyazaki 2023 (Fig. 3) instead derives 1-octen-3-ol from the 1-octen-3-hydroperoxide made by O2
  addition to the rearranged allyl radical, which is the more usual account; Yang's "rearrangement of
  2-octen-1-ol" is stated once without a drawing. Prefer Miyazaki's S4 for the repo; keep this as an
  alternative.

**S3. 2,4-Decadienal -> 2-octenal + hexanal (labelled "retro-aldol" in Fig. 4). Do NOT write as drawn.**
The arrow is not mass-balanced (C10 -> C8 + C6). A real retro-aldol of (E,E)-2,4-decadienal after
hydration at C5 gives hexanal + 2-butenal-type C4 fragment (C6 + C4), and after hydration at C3 gives
2-octenal + acetaldehyde (C8 + C2). If the repo wants a dienal-fragmentation rule it should come from a
paper that measures it (Chen 2017 cites Yu et al. 1994 for 2,4-decadienal -> pentanal, heptanal,
2-nonenal; also not on disk). Record Fig. 4's arrow only as a flagged claim.

## 6. Flags

1. **Products requiring a second double bond are drawn and reported from (9Z)-octadecenoic acid**
   (2-pentylfuran, 2-hexylfuran, 2,4-decadienal, 2-heptenal ... 2-nonenal, hexanal, heptanal from the
   dienals). Sigma oleic acid purity is not stated; the natural reading is a linoleate impurity. Do not
   cite this paper for oleate -> 2-pentylfuran or oleate -> 2,4-decadienal.
2. **Fig. 3 draws no isomer-specific arrows** for the oleate alkanals; the assignment "8-/9-/10-/11-HPOD
   -> heptanal to decanal" is a sentence. Cao 2020 (Fig. 6) is the drawn source; this paper corroborates
   presence.
3. **Fig. 4's hexanal / 2-octenal from 2,4-decadienal "retro-aldol"** is not mass-balanced as drawn
   (C10 -> C8 + C6). The LA map also gives two routes to hexanal (direct from 13-type HPOD in Fig. 7B;
   via dienal fragmentation in Fig. 4) without choosing.
4. **13-HPOD of LA is not labelled** in Fig. 4 (two unlabelled OOH structures); the text lists only 8-12.
   The repo's LOOH_13 rules are therefore not directly witnessed by this figure.
5. **No numbers at all** in the main text: Fig. 1 is log10(A) bars — FIGURE-ONLY. No replicate count,
   no internal standard, no SD. "Intensities" claims (e.g. AA ketones higher; 2,4-nonadienal << 2,4-
   decadienal) are qualitative.
6. **80 C / 30 min, 10 mg neat acid in a 20 mL sealed vial:** very mild, oxygen ample (~0.17 mmol O2 vs
   ~0.035 mmol acid), no matrix. Representative of stewing headspace, not frying.
7. **"Simulating"** in the title is a bench model system; nothing here is computed. No DFT; nothing to
   mark inadmissible.
8. **Trilinoleate = linoleic acid claim** (conclusion) has no data shown in the main text; supplementary
   not on disk.
9. Text says "2,4-denadienal obtained from the decomposition of 8-HPOD" (CLA) — typo for 2,4-decadienal;
   and "pentatonic acid" for pentanoic acid.
