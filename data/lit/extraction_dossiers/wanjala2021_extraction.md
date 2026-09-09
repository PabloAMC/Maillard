# Wanjala et al. 2021 — EXTRACTION (photosensitised linoleic acid hydroperoxide mixture + L-lysine in diethyl ether with BHT, 37 C / 30 min; GC-MS finds hexanal and 2-pentylfuran; drawn dioxetane / HNE-cyclisation route)
### A qualitative, hypothesis-driven note: lysine converts HpODE to hexanal and 2-pentylfuran with radicals scavenged; no amounts, no isomer resolution.

**Source on disk:** `data/articles/Wanjala2021.pdf` (owner's download, 2026-09-08; Scientific African 12
(2021) e00797, open access). Read from the `pypdf` text layer (7 pages, no tables); Schemes 1-5 are drawn
mechanisms, described in words from renderings of pages 2-6. Supplementary Figs S1 (TIC) and S2 (mass
spectra) are not on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Does lysine drive the conversion of fatty acid hydroperoxides to aldehydes and alkyl-furans?" |
| Authors | George W. Wanjala, Arnold N. Onyango, David Abuga, Calvin Onyango, Moses Makayoto (JKUAT; KIRDI, Nairobi) |
| Venue | Scientific African 2021, 12, e00797 (PII S2468227621001010) |
| DOI | 10.1016/j.sciaf.2021.e00797 |
| Naming | HPODE = hydroperoxy-octadecadienoic acid (free acid); 13-HPODE **1**, 9-HPODE **9**, 10-HPODE **10**, 12-HPODE **11**; HNE = 4-hydroxy-2-nonenal **18**; HPNE = 4-hydroperoxy-2-nonenal **14**; HEL = N-epsilon-(hexanoyl)lysine **7**; MDA = malondialdehyde **16**. Bold numbers are scheme numbers. |
| Prior hypothesis papers by the same group | Onyango 2016 (Oxid Med Cell Longev, singlet oxygen / dioxetane hypothesis); Onyango 2012 (Chem Phys Lipids 165:777); Onyango et al. 2010 (Food Res Int 43:925, dihydroperoxidation -> aldehydes); Wanjala et al. 2020 (J Chem, cholesterol-5-OOH + lysine) |

## 1. Why it matters

It is the one paper in the batch where an amine, not heat or metal, drives HpODE breakdown, and where
a radical scavenger is present — so the hexanal and 2-pentylfuran seen are claimed to arise by a
non-radical (dioxetane) route. For the repo's route (ii) it supplies a third, distinct 2-pentylfuran
mechanism: 9-/10-HPODE -> dioxetane -> (Z)-3-nonenal -> 4-hydroperoxy-2-nonenal -> 4-hydroxy-2-nonenal ->
cyclodehydration -> 2-pentylfuran (Scheme 3), with amino-acid catalysis of the last step cited to Adams
et al. 2011. It also says explicitly which HPODE isomers should give hexanal (12- and 13-) and which
should give 3-nonenal (9- and 10-) under the dioxetane mechanism (Scheme 2). It records no numbers, so
it can only anchor a mechanism_drawn rule, and one whose radical-free character is not proven.

## 2. Methods as they matter to a model

- **Hydroperoxide preparation:** linoleic acid 5 g (Sigma, purity not stated, "not purified before use")
  in 10 mL ethanol with 0.27 mM methylene blue, irradiated at 10 C with 366 nm UV (Funa SL-800G) from 25
  cm for 1 h. Dried over Na2SO4, silica-gel column (hexane/ethyl acetate 95:5). Hydroperoxydiene
  structure "confirmed by visualization on TLC by UV, and by their coloration with potassium iodide".
  **Isomers were NOT separated and their distribution was NOT measured**; the authors assume the
  photosensitised set 9-, 10-, 12-, 13-HPODE (Scheme 2, citing Minami et al. 2008). A trace of BHT was
  added at this stage. Conversion / yield of hydroperoxide from the 5 g acid is not stated.
- **Reaction:** hydroperoxide mixture **0.6 M in diethyl ether** with BHT + **L-lysine 0.3, 0.6 or 1.2 M**
  in sealed 10 mL vials; control vial without lysine. Either shaken 2 min or incubated **37 C / 30 min**
  with shaking. (Lysine's solubility in diethyl ether is not addressed; the mixture is presumably a
  suspension.) Volumes and headspace not stated.
- **Analysis:** 2 µL of the reaction mixture injected directly (liquid injection, not headspace) into
  GC-MS (Shimadzu QP2010 SE), BPX5 30 m x 0.25 mm x 0.25 µm, He 1 mL/min, 50 C (2 min) -> 150 C (1 min) at
  5 C/min; injector 240 C, EI 70 eV, scan m/z 35-1000. Identification by NIST match and retention
  time; hexanal and 2-pentylfuran standards were purchased (Sigma) — the text says identification was
  "based on comparison of retention time and mass spectra with those from the NIST Mass Spectral
  Library".
- **Quantification: none.** Results are a TIC (Fig S1) and two spectra (Fig S2). No areas, no internal
  standard, no replicate count, no lysine-concentration dependence reported, no time dependence.
- **Control:** "At a similar time, minimal product formation was detected for a reaction without lysine
  (not shown)."

## 3. Tables re-typed

No tables. Results in words: "The peaks at rt 4.45 min and 9.19 minutes were identified as hexanal 3
and 2-pentylfuran, respectively ... The mass spectrum for the peak at rt 8.46 (Fig S1 in SI) had the
NIST spectrum for 2-nonenal as a close match, but these spectra were not identical. Thus, the identity
of this peak requires further confirmation." (Z)-3-nonenal **13** was looked for and not detected.

### Scheme 1 — "Previously proposed pathways for the formation of hexanoyl-lysine as a product of the reaction of 13-HPODE 1 with lysine (RNH2)" (drawn; described; from Onyango 2016)

13-HPODE **1** + RNH2 <-> peroxide anion / ammonium ion pair -> loss of RNH2 with ring closure of the
peroxide oxygen onto C12 -> **1,2-dioxetane 2** spanning C12-C13 (the former 11E double bond is drawn
shifted; the dioxetane sits on the carbon bearing the OOH and its neighbour) -> dioxetane cleavage ->
**hexanal 3** (C13-C18) + the C1-C12 aldehyde (12-oxo-dodecadienoic acid, not numbered) -> hexanal +
RNH2 <-> carbinolamine **4** <-> Schiff base **5**; **5** + R'OOH -> peroxy-aminal **6** -> R''OH + HEL **7**;
alternatively **4** + 1O2 (from H2O2 / triplet carbonyl) -> **7**.

### Scheme 2 — "The expected lysine-catalysed conversion of different linoleic acid regioisomers to (Z)-3-nonenal 13 or hexanal 3" (drawn; described)

Linoleic acid **8** -> (photosensitised) -> 13-HPODE **1** (top; -> hexanal **3**) and, by routes A, B, C:
9-HPODE **9** (10E,12Z), 10-HPODE **10** (8E,12Z), 12-HPODE **11** (9Z,13E). **9** + RNH2 and **10** + RNH2 ->
the same **dioxetane 12** on C9-C10 (drawn with the 12Z double bond retained) -> **(Z)-3-nonenal 13**
(C10-C18) + 9-oxononanoic acid (not numbered). **11** -> dioxetane **2** (C12-C13) -> **hexanal 3**. Text:
"hydroperoxides 9 and 10 are expected to yield, via dioxetane 12, (Z)-3-nonenal (13), while 12-HPODE
(11), like 13-HPODE (1), affords hexanal (3)".

### Scheme 3 — "Suggested pathway for the formation of 2-pentylfuran (19) from dioxetane 12 via (Z)-3-nonenal (13) and 4-hydroxy-2-nonenal (18)" (drawn; described)

Dioxetane **12** -> triplet-excited (Z)-3-nonenal **13** -> energy transfer to 3O2 -> 1O2 + ground-state
**13** -> 1O2 ene reaction at the C3=C4 -> **4-hydroperoxy-2E-nonenal (HPNE) 14** (OOH on C4, new C2=C3).
Two branches: (a) **14** + lysine -> dioxetane **15** (C2-C3) -> **hexanal 3** + **malondialdehyde 16**
(+ NH3 from the lysine drawn); (b) **14** + lysine -> intermediate **17** (lysine alpha-carbon/amine adduct
on the peroxide oxygen, drawn as a cyclic transition structure) -> CO2 + 5-aminopentanal
H2N(CH2)3CH2CHO + **4-hydroxy-2-nonenal 18** -> cyclisation (loss of H2O) -> **2-pentylfuran 19**. Fe2+ is
drawn as an alternative reductant of **14** to **18**. Text: "HNE 18 then cyclizes to form 2-pentylfuran 19
[23, 9], and this cyclization can also be catalysed by lysine [24]" ([24] = Adams, Bouckaert, Van
Lancker, De Meulenaer, De Kimpe 2011, JAFC 59:11058, amino acid catalysis of 2-alkylfuran formation from
alpha,beta-unsaturated aldehydes).

### Scheme 4 — "The expected facile conversion of 10-HPODE 10 to octene radical 21 and 10-oxo-9-decenoic acid 22 during autoxidation" (drawn; described)

10-HPODE **10** -> alkoxyl at C10 (**20**) -> C10-C11 cleavage -> **2-octenyl (allylic) radical 21** +
**10-oxo-8-decenoic acid 22** (printed "10-oxo-9-decenoic"; the drawing has the C=C conjugated to the
CHO). Text: "formation of allylic radicals such as 21 is energetically favourable"; "12-HPODE 11 is
expected to be converted to 2-heptenal and an allylic radical". (Same C8 branch as Miyazaki 2023 Fig 3.)

### Scheme 5 — "Mechanism for the conversion of 13-HPODE 1 to HPNE 14 under autoxidative conditions" (drawn; described; after Schneider et al. 2001 and Onyango 2017)

13-HPODE **1** -> H-abstraction at C8 -> allylic radical **23** -> O2 at C8 -> 8-peroxyl **24** -> (AH) ->
**8,13-dihydroperoxide 25**; or **24** cyclises onto C9 -> dioxetanyl radical **26** (C8-C9) -> ring cleavage
-> carbon radical **27** = 4-hydroperoxy-2-nonenal-type C9 fragment with the radical at C2 (+ the C1-C9
oxo-acid, drawn as (CH2)5CH2OOH-bearing fragment) -> O2 -> peroxyl **28** -> (AH) -> 2,4-dihydroperoxy-
type **29**, or **28** cyclises -> peroxylactonyl radical **30** -> loses formyloxyl radical -> **HPNE 14**.
Text: antioxidants (BHT) trap **24**/**28** to the dihydroperoxides **25**/**29** and so limit HPNE by this
radical route; "the 13-hydroperoxy-group in 13-HPODE (1) remains intact during the auto-oxidative
conversion of the latter to HPNE (14) [26, 27]".

## 4. Routes and numbers the repository can use

Conditions for every row: HPODE mixture (isomers unresolved) 0.6 M in Et2O + trace BHT, lysine 0.3-1.2 M,
37 C / 30 min (or 2 min shaking), direct liquid injection GC-MS. No amounts.

| route | reactant -> product | mechanism as drawn (scheme) | measured numbers | evidence class |
|---|---|---|---|---|
| LYS-13-HEX | 13-HPODE (and 12-HPODE) + lysine -> hexanal + 12-oxo-dodecadienoic acid | amine-catalysed dioxetane on C12-C13, cleavage (Schemes 1, 2) | hexanal detected (rt 4.45 min); "a major product"; no amount | mechanism_drawn (hypothesis); measured_level (presence only) |
| LYS-9/10-3NON | 9-/10-HPODE + lysine -> (Z)-3-nonenal + 9-oxononanoic acid | dioxetane 12 on C9-C10 (Scheme 2) | 3-nonenal NOT detected | mechanism_drawn; null observation |
| LYS-PF | **(Z)-3-nonenal -> HPNE -> HNE -> 2-pentylfuran** (lysine, 1O2 from triplet carbonyl) | Scheme 3 | 2-pentylfuran detected (rt 9.19 min); HNE and MDA not detected ("poor volatility"); no amount | mechanism_drawn (hypothesis); measured_level (presence only) |
| LYS-HPNE-HEX | HPNE + lysine -> hexanal + malondialdehyde | dioxetane 15 (Scheme 3a) | not separately observable | mechanism_drawn |
| LYS-CTRL | HPODE mixture without lysine, same time | - | "minimal product formation ... (not shown)" | measured_level (qualitative null) |
| AUTOX-10-C8 | 10-HPODE -> 2-octenyl radical + 10-oxo-8-decenoic acid (radical route, suppressed here) | Scheme 4 | none | mechanism_drawn |
| AUTOX-13-HPNE | 13-HPODE -> 8-peroxyl -> dioxetanyl -> HPNE (radical route) or -> 8,13-dihydroperoxide (with antioxidant) | Scheme 5 | none | mechanism_drawn (after Schneider 2001) |
| OOH-RED | hydroperoxides -> alcohols promoted by lysine; alkoxyls -> alcohols favoured by antioxidants | text, citing Martin-Rubio et al. 2019 and Onyango et al. 2010 | none | claim (cited) |

Within-study ratio: none printable (TIC only, in SI).

## 5. Rule sketches (repository suggestions, not the paper's)

Keys: `hexanal`, `2_pentylfuran` exist in `data/keys/compounds.yml`; **3-nonenal, 4-hydroxy-2-nonenal,
4-hydroperoxy-2-nonenal, malondialdehyde, 9-oxononanoic acid (free acid; `ME_9_OXONONANOATE` is the
methyl ester), 2-nonenal** do not. `LOOH_13_ct` / `LOOH_9_ct` are the methyl esters of 13-/9-HPODE;
10-/12-HPODE have no structure entry. Free linoleic acid has no entry. Lysine exists as `Lys` in
`data/species/structures.yml`.

**S1. Amine-catalysed dioxetane cleavage of 13-HPODE -> hexanal (net; same products as R18b).**
This is mechanistically distinct from R18b (no alkoxyl, no radical) but the net transformation is
identical: C12-C13 bond breaks, C13 becomes CHO of hexanal, C12 becomes CHO of the C12 oxo-acid; the OOH
oxygens end up one on each aldehyde. If the repo wants to keep mechanisms apart, add a `conditions:`
"amine present, BHT present, 37 C" variant pointing to this dossier rather than a new SMIRKS.
- positive control: `CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)OC` (LOOH_13_ct) + `NCCCCC(N)C(=O)O` (Lys, catalyst, unchanged) -> `CCCCCC=O` (HEXANAL) + `O=C/C=C/C=C\CCCCCCCC(=O)OC`-type C12 oxo-dienoate (the paper does not draw the co-product's geometry; write it as methyl 12-oxo-dodeca-9,11-dienoate)
- second positive (12-HPODE free acid, predicted by Scheme 2 to give hexanal too): `CCCC/C=C/C(OO)C/C=C\CCCCCCCC(=O)O` -> `CCCCCC=O` + C12 fragment — NOTE this is the same 12-HPODE that Miyazaki 2023 found gives 2-heptenal (not hexanal) on heating; the two papers' mechanisms make opposite predictions for 12-HPODE, and Wanjala's mixture cannot tell which isomer gave the hexanal.
- negative control: methyl stearate `CCCCCCCCCCCCCCCCCC(=O)OC` (no OOH); tert-butyl hydroperoxide `CC(C)(C)OO` (the paper's own cited null: t-BuOOH + hexanal + lysine gave no HEL; also cannot form a dioxetane on an alkene).

**S2. 4-Hydroxy-2-nonenal -> 2-pentylfuran + H2O (cyclodehydration; net; amine-catalysed per Adams 2011).**
`O=[CH1:1][CH1:2]=[CH1:3][CH1:4]([OH:5])[CH2:6][#6:7]` -> `[cH:1]1[cH:2][cH:3][c:4]([CH2:6][#6:7])[o:5]1`
(the hydroxyl oxygen becomes the ring oxygen bonded to C1; the aldehyde oxygen leaves as water; the ring
is C1-C2-C3-C4-O).
- positive: `CCCCCC(O)/C=C/C=O` (4-hydroxy-2E-nonenal) -> `CCCCCc1ccco1` (2-pentylfuran, key `2_pentylfuran`)
- second positive: `CCCCC(O)/C=C/C=O` (4-hydroxy-2-octenal) -> `CCCCc1ccco1` (2-butylfuran)
- negative: 2-nonenal `CCCCCC/C=C/C=O` (no 4-OH: must not fire by THIS rule); hexanal; nonanal.
- note: Adams et al. 2011 (not on disk) is the primary for the amino-acid catalysis of 2-alkylfuran
  formation from 2-alkenals / 4-hydroxy-2-alkenals; if the repo wants this step as "established" it
  needs that dossier.

**S3. (Z)-3-nonenal -> 4-hydroperoxy-2E-nonenal (1O2 ene / autoxidation) -> 4-hydroxy-2E-nonenal (reduction).**
Two net steps: `CCCCC/C=C\CC=O` -> `CCCCCC(OO)/C=C/C=O` (O2 adds at C4, C=C shifts into conjugation);
`CCCCCC(OO)/C=C/C=O` -> `CCCCCC(O)/C=C/C=O` (OOH -> OH). Cited primaries for the first step: Gardner &
Hamberg 1993 (J Biol Chem 268:6971) and Noordermeer et al. 2000 (BBRC 277:112) — both plant / non-
enzymatic 3Z-alkenal -> 4-hydroxy-2E-alkenal. Negative: nonanal, 2-nonenal (already conjugated: 1O2 ene
at C2=C3 would not give the 4-OOH).

**S4. Chain of S1-type dioxetane on 9-/10-HPODE -> (Z)-3-nonenal + 9-oxononanoic acid (net).**
- positive: `CCCCC/C=C\C=C\C(OO)CCCCCCCC(=O)OC` (LOOH_9_ct) -> `CCCCC/C=C\CC=O` ((Z)-3-nonenal) + `O=CCCCCCCCC(=O)OC` (ME_9_OXONONANOATE)
- this is the R18b "side B" of LOOH_9 with the alkene geometry kept Z; the repo's R18b already yields
  ME_9_OXONONANOATE from LOOH_9_ct and a C9 3-enal; the only new content is that Wanjala predicts it
  and did NOT observe 3-nonenal (consumed onward to 2-pentylfuran, they argue), matching Miyazaki 2023's
  non-detection of 3-nonenal from 9-HpODE by heat.
- negative: methyl stearate (no OOH, no C=C). Oleate 9-OOH is NOT a safe negative: it has C10=C11
  adjacent to the OOH carbon and the same drawing would close a C9-C10 dioxetane; the paper is silent
  on mono-enes.

## 6. Flags

1. **Qualitative only.** No peak areas, no calibration, no replicates, no dependence on lysine
   concentration (three concentrations were run but no difference is reported), no time dependence.
   The TIC and spectra are in supplementary material not on disk.
2. **Isomer mixture, unresolved and unmeasured.** Whether hexanal came from 13-, 12- or (after
   isomerisation) 9-HPODE is not knowable from this experiment; the isomer assignments are Scheme 2's
   expectations.
3. **"Non-radical" is asserted, not demonstrated.** A "trace" of BHT is the only radical control; the
   paper itself concedes "it is perhaps not possible to completely prevent free radical reactions during
   the decomposition of hydroperoxides" and gives Fe2+ as an alternative reductant of HPNE. The lysine-
   free control is "not shown". Evidence class for the mechanism is mechanism_drawn (hypothesis), not
   established.
4. **Solvent/phase:** 0.3-1.2 M L-lysine in diethyl ether is a suspension; the reaction is heterogeneous
   and the "0.6 M hydroperoxide" concentration is nominal.
5. **Substrate purity** (Sigma linoleic acid, "not purified") and hydroperoxide yield not stated.
6. **Identification:** NIST match plus retention time; standards were bought but the text does not say
   they were co-injected. The rt 8.46 peak (2-nonenal-like) is unassigned.
7. **HNE, MDA, 5-aminopentanal, HEL** — none of the intermediates of Scheme 3 was detected; their
   absence is explained by volatility / reactivity.
8. **Scheme 4 product name** "10-oxo-9-decenoic acid" vs the drawing (conjugated 2-enal = 10-oxo-8-
   decenoic acid): naming slip.
9. **12-HPODE prediction conflict** with Miyazaki 2023 (hexanal here vs 2-heptenal there); see S1.
10. No DFT, no rates; nothing to mark inadmissible.
