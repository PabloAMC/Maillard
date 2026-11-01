# Mottram, Nobrega, Dodson & Elmore 1996 — EXTRACTION (cooked beef: SDE vs Tenax headspace give opposite thiol / disulfide pictures; furyl disulfides spiked into cooked beef before 2 h SDE come out ~60 % as free thiol; Table 3 reprints the egg albumin / casein / maltodextrin data of Mottram, Szauman-Szumski & Dodson 1996)

### An OBSERVATION-OPERATOR paper. SDE of a protein matrix is a reducing operator: it turns furyl disulfides into thiols. No rate, no temperature series, no pH; one endpoint per system.

**Read 2026-10-31 from `data/articles/mottram1996b.pdf`** (6 pp., Flavour Science: Recent Developments, RSC 1996,
pp. 413-418; an ABBYY OCR of a scan, text layer usable). Tables 1 and 2 checked against 110 dpi page renders
(scratchpad `dos7/mot-3.png`, `mot-4.png`); the OCR layer drops the MFT headspace cell of Table 1 ("-" on the render)
and garbles symbols (µg as "pg" / "jag"). No figures.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

**How this differs from `mottram1996_extraction.md`.** That dossier is Mottram, Szauman-Szumski & Dodson, *J. Agric.
Food Chem.* 1996, 44:2349 (doi 10.1021/jf960170o): disulfides heated with egg albumin, casein or maltodextrin in
water (SDE, 2 h), plus a commercial flavouring. **This** paper (Mottram, Nobrega, Dodson & Elmore, a symposium
chapter) adds two things the JAFC paper does not have: (i) a cooked-beef SDE vs Tenax-headspace comparison (Table 1),
and (ii) the two furyl disulfides spiked into cooked beef before SDE (Table 2). Its **Table 3 is a reprint** of the JAFC
Table 2 ("Adapted from Mottram et al." [V], ref 15), cell for cell identical for albumin / casein / maltodextrin
(checked: e.g. Exp 3 MFT disulfide 6 (2) / 272 (21) / 378 (24)); the JAFC "water only" column is omitted. **Do not
count Table 3 as a second measurement.** A third related text exists, Mottram, Nobrega & Dodson 1998, ACS Symp. Ser.
705, 78-84, "Extraction of thiol and disulfide aroma compounds from food systems" (cited as ref 9 by
adams2001) — **not on disk**; it may reprint Tables 1-2 again [C/!].

## 1. What this paper answers for the model

| what the model needs | what Mottram 1996b gives | tag |
|---|---|---|
| MFT / FFT loss law at process temperature | **nothing kinetic.** One endpoint after 2 h SDE (~100 °C) of cooked beef spiked AFTER cooking | [M] |
| reversibility: disulfide vs covalent sink | in cooked beef, ~79-80 % of each spiked disulfide is consumed and **58-60 % of the thiol units return as free thiol**; mixed and methyl disulfides add 4-5 %; **only 10-17 % of the units are unaccounted** (protein-bound, inferred) | [M] → [Z] |
| matrix contrast | beef (this paper) vs albumin / casein / maltodextrin (reprint): albumin ≈ beef in disulfide consumption, but albumin leaves 54-74 % unaccounted (mottram1996 dossier) vs 10-17 % in beef [Z] | [Z] |
| operator artefact | **SDE vs Tenax headspace on the same cooked beef:** SDE shows thiols (MFT, FFT, mercaptoketones) and almost no disulfides; headspace shows disulfides (dimethyl di/tri/tetrasulfide, furfuryl methyl disulfide) and no MFT, trace FFT (Table 1, relative areas) | [M] |
| which operator is "true" | **not tested.** The authors attribute the difference to redox reactions with protein "during extraction and/or analysis" [Q]; no spike was run through the headspace method | [Q] / [!] |
| pH, protein SH / SS titre | not printed. Text: "Albumin and meat protein contain much greater proportions of -SH groups than casein" [Q], no number | [Q] |

**Answers.**

1. **SDE of a protein-containing matrix is a reducing observation operator** [M/Q]. Spiked bis(2-furylmethyl)
   disulfide and bis(2-methyl-3-furyl) disulfide (500 µg each, both in 100 g cooked beef) come out of 2 h SDE as
   104 / 126 µg disulfide and 294 / 302 µg thiol; from water, 437 / 532 µg disulfide and < 1 µg thiol (Table 2).
   Molar thiol yields 58.3 % (FFT) and 59.9 % (MFT) of the thiol units charged [Z]. A model that compares engine
   thiol levels with SDE data from meat or protein-rich pots compares a thiol + disulfide sum, partly re-split.
2. **In cooked beef the disulfide is mostly a reservoir, not a sink** [Z]: thiol units accounted (thiol + 2 ×
   disulfide + mixed + methyl disulfides) **83 % (FFT) and 90 % (MFT)**. The authors' "binding to the protein" [P] is
   the 10-17 % remainder; no bound form was measured.
3. **No loss rate.** The beef was autoclaved (140 °C, 30 min) BEFORE the spike, so the spike saw only the SDE (2 h,
   ~100 °C). Volatiles are removed continuously during SDE, so the residence time of a molecule in the hot pot is
   unknown. Nothing here constrains `k_mft_decay` [!].

## 0. Identity

| field | value |
|---|---|
| Title | "Changes in Thiol and Disulfide Flavour Compounds Resulting from the Interaction with Proteins" |
| Authors | Donald S. Mottram, Ian C. Nobrega, Andrew T. Dodson, J. Stephen Elmore (University of Reading) — **four authors**; the brief listed three |
| Venue | *Flavour Science: Recent Developments* (eds. Taylor & Mottram), Royal Society of Chemistry, 1996, pp. 413-418 (section "Flavour Binding and Flavour Release") |
| DOI | **10.1533/9781845698232.7.413** (not printed; Crossref 2026-10-31: same title, Mottram / Nobrega / Dodson / Elmore, "Flavour Science", 413-418, 1996). Matches the PDF |
| Paper type | Symposium chapter; GC-MS, triplicates; no statistics beyond SD |

## 2. Conditions

| item | as printed |
|---|---|
| matrix | beef fillet (M. psoas major), fat trimmed, 100 g portions vacuum-packed, -20 °C; minced, in 500 mL bottles (PTFE-lined lids), **autoclaved 140 °C, 30 min** [M]. The Results text calls it "the boiled meat" [!] |
| spike | bis(2-methyl-3-furyl) disulfide (gift, flavour company) and bis(2-furylmethyl) disulfide (Aldrich), 1 mg/mL in ethanol, 0.5 mL **of each** into 100 g cooked minced beef, AFTER cooking, before SDE [M] |
| SDE (= the heating of the spike) | sample + 750 mL distilled water, modified Likens-Nickerson, **2 h**, n-pentane 27 mL + diethyl ether 3 mL; IS 1,2-dichlorobenzene 130 µg added AFTER extraction; Kuderna-Danish to ~0.5 mL [M] |
| headspace | cooked beef in a 250 mL flask at **60 °C**, oxygen-free N2 40 mL/min, **1.5 h**, Tenax TA 85 mg; IS 100 ng dichlorobenzene on the trap; thermal desorption 250 °C [M] |
| pH, buffer, O2 | not printed; SDE under air (not stated otherwise); headspace under N2 |
| GC-MS | HP5890 / 5972, BPX5 50 m; quantities of compounds 1-4 (MFT, its disulfide, FFT, its disulfide) against authentic standards; all others TIC area / IS with response factor 1 [M] |
| replication | triplicate analyses; SD in parentheses (Tables 1-2) [M] |

Concentrations [Z]: each disulfide 500 µg / 226.31 g/mol = 2.21 µmol in ~0.85 L (750 mL water + 100 g beef) ≈ 2.6 µM.

## 3. Quantitative results

**Table 1 (p. 415): cooked beef, relative TIC area × 10 000 (vs IS, RF 1), mean (SD), n = 3** [M]:

| compound | SDE | headspace (Tenax) |
|---|---|---|
| 3-mercapto-2-butanone | 20 (1) | - |
| 3-mercapto-2-pentanone | 13 (3) | - |
| 2-mercapto-3-pentanone | 8 (1) | - |
| 2-methyl-3-furanthiol | 4 (1) | - (render; OCR layer blank) |
| 2-furanmethanethiol | 62 (17) | tr |
| dimethyl disulfide | - | 356 (201) |
| dimethyl trisulfide | 9 (2) | 253 (82) |
| dimethyl tetrasulfide | - | 29 (14) |
| dipropyl disulfide | - | 7 (1) |
| 2-furanylmethyl methyl disulfide | 2 (1) | 23 (7) |
| 2-methyl-3-furanyl methyl disulfide | - | 7 (2) |
| 2-furanylmethyl methyl trisulfide | - | 6 (3) |
| 2-methyl-3-furanyl methyl trisulfide | - | tr |
| bis(2-furanylmethyl) disulfide | - | tr |
| 2-furanylmethyl 2-methyl-3-furanyl disulfide | - | tr |
| 4-methylthiazole / 4,5-dimethylthiazole | - / - | 11 (2) / 5 (1) |
| 2-acetylthiazole | 10 (3) | 26 (7) |
| 3-methylthiophene | - | 66 (32) |
| thiophene-2-carboxaldehyde / 3-methyl- / 5-methyl- | 7 (1) / 9 (3) / 2 (1) | 21 (7) / 5 (3) / 7 (2) |
| 2,3-dihydro-6-methylthieno[2,3c]furan | 5 (2) | 12 (0) |

The two columns are different operators on different subsamples (100 g + 750 mL water SDE vs chopped meat at
60 °C); relative areas are not comparable in absolute terms, only the thiol : disulfide pattern within each [Z].

**Table 2 (p. 416): µg recovered by SDE, mean (SD), n = 3** [M]:

| compound | added (µg) | water + disulfides | cooked beef + disulfides |
|---|---|---|---|
| bis(2-furanylmethyl) disulfide | 500 | 437 | 104 (27) |
| bis(2-methyl-3-furanyl) disulfide | 500 | 532 | 126 (27) |
| 2-furanmethanethiol | 1 (impurity) | < 1 | 294 (41) |
| 2-methyl-3-furanthiol | 1 (impurity) | < 1 | 302 (35) |
| 2-furanylmethyl 2-methyl-3-furanyl disulfide | — | 16 | 39 (6) |
| 2-methyl-3-furanyl methyl disulfide | — | — | 7 (1) |
| 2-furanylmethyl methyl disulfide | — | — | 3 (0) |

No SD for the water column [!]. Molar bookkeeping [Z] (MW disulfides 226.31, thiols 114.16, methyl disulfides
160.26; endogenous beef thiols ignored — Table 1 shows they are small relative areas, not µg):

| | disulfide left | thiol yield (of 2 × disulfide) | units accounted incl. mixed + methyl | without them |
|---|---|---|---|---|
| FFT, beef | 20.8 % | 58.3 % | 83.4 % | 79.1 % |
| MFT, beef | 25.2 % | 59.9 % | 90.0 % | 85.1 % |
| FFT disulfide, water | 87.4 % | ~0 | — | — |
| MFT disulfide, water | 106.4 % | ~0 | — | — |

Mixed disulfide 16 µg (0.071 µmol) forms from the two symmetrical disulfides in water alone (no thiol detected):
disulfide-disulfide scrambling during SDE without protein [M/Z].

**Table 3 (p. 417)** = reprint of mottram1996 Table 2 (see header). Not re-keyed here.

## 4. Mechanism evidence and its kind

- **Kind: within-study factor (matrix present / absent) + partial mass balance**, one laboratory (Reading, Mottram
  group). No labelling, no trapping of the bound form, no time course.
- Disulfide → thiol in the presence of meat: shown (K6, Table 2). Thiol-protein mixed disulfide as the sink:
  **[P]**, inferred from the 10-17 % gap; never measured.
- Mechanism (sulfhydryl-disulfide interchange with protein cysteine / cystine) [P], citing Jocelyn 1972 and
  Whitesides 1983 [C].
- The SDE vs headspace difference is attributed to protein redox during extraction [Q]; the alternative (oxidation of
  thiols on the Tenax trap or during 250 °C desorption, or at 60 °C under N2 purge) is not excluded [!].

## 5. Bearing on open questions

- **Evidence map 2.11 (thiol work-up artefacts, CAN):** SUPPORTS, and extends it from albumin to cooked meat. Same
  laboratory as mottram1996 (Reading), so it is not a second lab.
- **Evidence map 2.4, "disulfide as a reversible reservoir vs a terminal sink: CAN":** SUPPORTS (reservoir, ≥ 83 %
  of units recovered in beef). Note the contrast with albumin (54-74 % unaccounted in mottram1996): the sink share is
  matrix-dependent [Z].
- **MFT loss law (EX-050, `docs/history/thiol_rivals_2026-10-31/`):** LEAVES UNCHANGED (no rate). Relevant only to
  the live rival "a reversible MFT-cysteine mixed disulfide" (red team): protein-SH can both reduce MFT-S-S-MFT to MFT
  and hold MFT as a mixed disulfide [Q/P]; nothing quantitative.
- **CF-013 (rising vs falling thiol time courses):** LEAVES UNCHANGED.
- **CF-031 (lipid on thiols):** not touched (fat trimmed).
- **Engine (`src/kinetic_core/sulfur.py` `ch_protein_ss_*`, terminal PRB; `k_dimer_release` held at zero):** the
  engine has no protein-thiol + flavour-disulfide → flavour thiol channel; this paper shows that direction dominates
  in cooked beef during SDE [M/Z]. As in mottram1996: a qualitative gate, not a constant.

## 6. Disagreements found

1. **Matrix-dependent unaccounted share** [Z]: albumin (mottram1996, 67 g/L, SDE 2 h) leaves 54-74 % of thiol units
   unaccounted; cooked beef (this paper, ~12 % w/v) 10-17 %. Candidate causes: protein SH / SS titre and
   accessibility (albumin denatures during SDE; the beef was pre-cooked), protein-to-spike ratio, and the beef's own
   reductants (glutathione, cysteine). Not a contradiction of chemistry; a reason not to transfer one matrix's sink
   size to another.
2. **"Boiled meat" (Results) vs "autoclave at 140 °C for 30 min" (Experimental)** [!], same paper.
3. **Disulfides in cooked-meat literature measured by SDE are biased low and thiols high** when protein is present
   [Z from Table 2]; SDE-based thiol / disulfide ratios from meat (e.g. Table VII of grosch1992, boiled meat by SDE)
   inherit this. Direction known, size only for 2.6 µM spikes in pre-cooked beef.

## 7. Ledger-ready measurements

- **Table 2, beef column vs water column** (two readers): the cleanest number for an SDE operator in a protein matrix
  — disulfide 21-25 % left, thiol 58-60 % of units, 2 h SDE. Worth keying as an OPERATOR row (observation model),
  not a kinetic row.
- **Table 1** (pattern only): a qualitative operator contrast (SDE → thiols, headspace → disulfides) on the same
  cooked beef; relative areas with RF 1, not quantitative. Key only as a flag.
- Table 3: do NOT key (reprint of mottram1996 Table 2; key that one once).
