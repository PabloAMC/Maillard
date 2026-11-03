# Madruga & Mottram 1998 — EXTRACTION (file stem `madruga1998`)
DOI: **10.1590/S0103-50531998000300010** (not printed in the PDF; SciELO PDF `jbchs/v9n3/a10v09n3.pdf`).
Crossref 2026-11-03 (`api.crossref.org/works/10.1590/S0103-50531998000300010`): title "The effect of pH on the
formation of volatile compounds produced by heating a model system containing 5'-Imp and cysteine", authors Marta
Suely Madruga, Donald S. Mottram, *Journal of the Brazilian Chemical Society* **9**(3), published 1998-05; Crossref
carries no page range. **Matches the PDF** (running head "J. Braz. Chem. Soc., Vol. 9, No. 3, 261-271, 1998";
received 2 Nov 1997). No naming caveat: one 1998 Madruga paper on disk, stem is unambiguous.
**Not to be confused with** the companion **Madruga & Mottram 1998, Ciênc. Tecnol. Aliment. 18(4), "Sensory analysis
of a model system using 5'-IMP and cysteine at different pH"** (Crossref 10.1590/S0101-20611998000400007; not on disk,
not in the repo), nor with `madruga1995` (beef, HCl pH ladder) or `madruga1997` (beef + IMP / thiamine / cysteine).

### 5'-IMP 34.8 mg/mL + L-cysteine 12.1 mg/mL (≈ 100 + 100 mM [Z]) in pyrophosphate buffer at pH 3.0 / 4.5 / 6.0, 1 mL in N2-flushed flame-sealed Pyrex ampoules, autoclave 140 °C / 1 h; dynamic headspace (60 °C, 1 h, Tenax) → GC-FID/odour port and GC-MS on DB-5 and CP-Wax. **One table of 98 compounds coded +, nd or * per pH — presence/absence only. No internal standard, no peak area, no concentration, no precursor loss, no ribose, no hypoxanthine, no H2S, one T, one t.**

**Read 2026-11-03 from `data/articles/madruga1998.pdf`** (11 pp., born-digital Acrobat PDFWriter text with Fig. 1 as
an embedded 1-bit scan). Prose from `pdftotext -layout`. **Table 1 (pp. 264-267, rotated landscape) read cell by cell
from 200 dpi renders of pp. 4-7**; rows 85-98 re-checked at 300 dpi (`p7crop-07.png`). The text layer agrees with the
renders in every +/nd/* cell and every LRI checked. **Fig. 1 (p. 263) read from the 200 dpi render of p. 3**; relative
peak heights only, read error ±0.05 of full scale. Scratch: `.../scratchpad/ingest/madruga1998_work/` (`p-03..07.png`,
`table1.tsv` = my transcription, counts by `awk` [Z]).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Two affiliations printed: Madruga at DTQA, Universidade Federal da Paraíba (UFPB), João Pessoa; Mottram
at the Department of Food Science and Technology, University of Reading. No corresponding author marked. **The bench work
was done at Reading**: the acknowledgment thanks Mottram "for his supervision" and "the University of Reading where the
work was developed" [V]; methods cite Madruga's 1994 Reading PhD thesis (ref. 13) and the 8 synthesised standards come
from it; the apparatus (CERTOclav, Tenax, HP5988A, SGE odour port, pyrophosphate-buffered ampoules) is the Reading
group's (cf. `meynier1995_extraction.md`). Written up after Madruga moved to UFPB (received Nov 1997). **Lab id:
existing `reading_mottram`, ledger string "Reading (Mottram)"** (confidence ~90 %). Not "UFPB Joao Pessoa (Madruga)"
(that string, used for `madruga1997` and `desousafontes2024`, has no entry in `data/keys/laboratories.yml`; the
desousafontes2024 dossier proposes `ufpb_madruga`). Same-lab dossiers: `madruga1995`, `mottram1995`,
`meynier1995` (same buffer, same 140 °C / 1 h ampoule protocol, ribose + cysteine — the natural comparator),
`mottram1995b`, `mottram1996b`; same first author: `madruga1997`.

## 1. What this paper answers for the model

| model need | what Madruga & Mottram 1998 gives | tag |
|---|---|---|
| MFT / FFT / thiophenethiols vs pH from a nucleotide pentose source | detected (+) or not (nd) at pH 6.0 / 4.5 / 3.0; **no amounts** | [M] |
| furyl / thienyl disulfides vs pH | 12 disulfides: 0 / 8 / 12 detected at pH 6.0 / 4.5 / 3.0 | [M] |
| non-sulfur volatiles (furfural, norfuraneol, furans, ketones, pyrazines) | presence only; **no pyrazine at any pH**; norfuraneol only at pH 4.5 | [M] |
| quantities (absolute, IS-relative, peak area) | **none**; "total quantity ... increased significantly as the pH decreased" [V] with no number | [Q] |
| pH before / after heating | "decreased by less than 0.2 pH unit after heating" [V]; per-pot values not printed | [M] |
| IMP loss, IMP → inosine → ribose / R5P / hypoxanthine | **not measured**; IMP as ribose source is [P] | [P] |
| H2S, cysteine loss, browning | **not measured** | — |
| T or t series | **none** (140 °C, 1 h only) | — |
| odour descriptions | Fig. 1 sniffing summary (pH 3.0); 18 "meaty" compounds flagged in Table 1 (footnote e) | [M] |

**Answers to the paper-specific questions.**

1. **Exact pots** (confidence ~95 % on the printed numbers). Cysteine 0.121 g + 5'-IMP 0.348 g dissolved in 10 mL
   "appropriate pyrophosphate buffer (pH 3.0; 4.5, 6.0)" [M]; final 12.1 mg/mL Cys, 34.8 mg/mL IMP [M]. 1.0 mL aliquots
   into Pyrex ampoules, N2 blown over for 3 min, flame-sealed [M]; horizontal in a CERTOclav autoclave, **1 h at
   140 °C, 0.28 MPa (2.7 bar)** [M]; reactions in triplicate [M]. Molar [Z]: Cys 0.121/121.16 = **99.9 mM**; IMP
   0.348/348.21 (free acid) = **99.9 mM** — the masses are exactly 1.000 mmol of each on free-acid weights, i.e. an
   intended **equimolar 100 mM : 100 mM** pot (~85 %). If the reagent was the usual disodium salt (not stated) the IMP
   would be 88.7 mM (anhydrous, 392.2) or ~66 mM (7.5-hydrate) [Z, !]. Buffer molarity **not printed** (the same lab's
   Meynier 1995 pots used 0.2 M; ~70 % the same here). Ampoule volume and so headspace volume **not printed**. pH:
   initial 3.0 / 4.5 / 6.0, measured before and after, drop < 0.2 unit [M]; per-pot final values not printed.
   Come-up / cool-down not stated. The claim that the concentrations "were related to the amount of these compounds
   in meat muscle" [V] is loose: on `madruga1997`'s own "ten times normal" doses, 34.8 mg/mL IMP is ~13× and
   12.1 mg/mL Cys ~300× the native beef level [Z].
2. **Volatiles and units** (~98 %). 98 compounds (Table 1, transcribed in full in sec. 4a): MFT [19], FFT [29],
   2-methyl-4,5-dihydro-3-furanthiol [36], 2- or 3-thiophenethiol [37], **2-methyl-3-thiophenethiol [50]**,
   2-thiophenemethanethiol [54], **bis(2-methyl-3-furyl) disulfide [88]** and 11 other furyl/thienyl/oxoalkyl
   disulfides [87-98], three mercaptoketones [14, 26, 28], four thiophenones [33, 38, 40, 49], 14 thiophenes,
   dithiolanones/trithiolanes/dithianones/trithianes/tetrathiane, 8 bicyclics (thienothiophenes, kahweofuran),
   5 thiazoles, 14 furans incl. 2-furfural [15] and norfuraneol [21], 10 ketones incl. 2,3-pentanedione [8]; **no
   pyrazine, no 2-acetylthiazole, no thiazoline**. **Units: none** — the pH columns carry "+ compound gives clear MS;
   nd not detected; * compound gives no clear MS, detected mixed with other volatile" [V, footnote b]. The + / nd
   boundary is an MS-quality call, not a concentration threshold; no LOD.
3. **pH dependence** (direction ~85 %, magnitude 0 %). Detected counts 45 / 71 / 86 at pH 6.0 / 4.5 / 3.0 [Z].
   **MFT, FFT and 2-methyl-3-thiophenethiol: nd at pH 6.0, + at 4.5 and 3.0** [M]; the 2-methyl-3-furyl disulfides
   0 / 5 / 7 (pH 6.0 / 4.5 / 3.0); mercaptoketones 26 and 28 nd at pH 6.0; acyl-thiophenes [55, 63, 66, 68] only at
   pH 3.0; alkyl-/formyl-thiophenes, alkylthiazoles and trithiolanes at all three pH; trithianes/tetrathiane and
   short ketones more at pH 6.0 [M/Q]. So yes: **the same direction as ribose + cysteine** — Meynier & Mottram 1995
   (same lab, same buffer, 140 °C / 1 h ampoules; ribose 35 + Cys 50 mM final) measured MFT 15.2 → 1.0 and FFT
   8.5 → 1.5 µg/10 mg ribose from pH 4.5 to 6.0 (`meynier1995_extraction.md` Table 3), and Madruga & Mottram 1995
   beef MFT rises as pH falls. The paper has **no ribose + cysteine arm of its own** (it cites Farmer 1989 at pH 5.7
   and Mulders 1973 at pH 5.6 only for compound lists). The sulfur lane predicts the same sign for both MFT and FFT
   (acid-catalysed deoxyosone/furfural formation and thiolate-loss terms, `src/kinetic_core/sulfur.py` header,
   "pH, parts 2-3"), so this is a sign-agreement, not a test: (a) nd at pH 6.0 is uninterpretable as a fold-change
   without a detection limit; meynier1995's ~15× and ~6× falls would leave both thiols detectable, so the IMP pot is
   either poorer in thiols or the headspace/MS call is less sensitive (cannot tell); (b) **pH 3.0 is outside the
   lane's stated pH 4.5-7 window** (`sulfur.py`, "alkaline-pH wall" note); (c) the pentose supply is itself
   pH-dependent (answer 4). Where the IMP pot **differs** from meynier1995's ribose + Cys at matched pH: no
   pyrazines at pH 6.0 (meynier1995 methylpyrazine 0.9, dimethylpyrazine 0.3); thiazoles flat (meynier1995 thiazole
   and 2-acetylthiazole rise with pH; 2-acetylthiazole absent here); a thienothiophene [78] present at pH 6.0
   (meynier1995 thieno[2,3-b]thiophene nd at 6.0, possibly a different isomer); norfuraneol detected only at
   pH 4.5 (meynier1995: present 10.6 → 2.6 at all pH).
4. **Usability as a pH-series check for the sulfur lane** (calibrated). The engine has **no IMP, inosine,
   hypoxanthine, ribose-5-phosphate or nucleotide species** (`grep -i "inosin|nucleotid|\bIMP\b|r5p|hypoxanth"
   src/kinetic_core/` → no hits). Nothing in the paper measures the hydrolysis step: no IMP loss, no free ribose,
   R5P, inosine or hypoxanthine. Free ribose needs the **N-glycosidic** bond broken (inosine is not a reducing
   sugar); Shaoul & Sporns 1987 (`data/articles/shaoul1987.pdf`, pending; read here from its text layer only, not by
   eye — [C, !]) give IMP **phosphate** hydrolysis to inosine at 121 °C with half-lives ~32 / ~64 / ~71 min at
   pH 3.0 / 4.5 / 6.0, and say glycosidic cleavage becomes significant only at pH 3. So the ribose supply rises
   steeply toward pH 3 on its own, and an IMP pot's MFT pH slope = (release slope) × (sulfur-chemistry slope); this
   paper cannot separate them. P(usable as a quantitative check) ≈ 1 %; P(usable as a sign-only directional row
   that adds information beyond meynier1995) ≈ 10 %; P(useful as a qualitative constraint on a future IMP → ribose
   release step, e.g. "MFT appears at pH 4.5 from IMP within 1 h at 140 °C") ≈ 15 %. **Recommendation: ledger it as
   presence-only, `simulable: false`; do not build a check on it.**
5. **Proposed ledger entries** (sec. 7): **one** sulfur.yml entry, `madruga1998_tab1_imp_cys_ph_presence`
   (presence-only, `also_relevant: trunk` / `aroma_other`); no separate trunk.yml or aroma_other.yml entry, because
   the non-sulfur rows are the same data set and carry even less (furfural + at every pH; no pyrazine anywhere).
   **Duplication**: not a re-report of any ledgered pot (~95 %). `madruga1995` is HCl-acidified beef and
   `madruga1997` beef + 2.7 g IMP / 100 g at pH 5.6 / 30 min — different systems. The same IMP + Cys pots are very
   likely also behind (i) the companion CTA 1998 sensory paper (same title system, not on disk), (ii) Madruga's 1994
   Reading thesis (ref. 13) and possibly (iii) Mottram & Madruga 1994a/b ACS / Elsevier chapters (refs. 16-17,
   cited for meaty disulfides; ~40 % they show IMP + Cys data). None of these is in the repo; if any is ingested
   later it must be keyed as the same pots, not a second laboratory.

## 2. Conditions

- **Reactants** [M]: 5'-IMP 0.348 g + L-cysteine 0.121 g per 10 mL buffer (12.1 / 34.8 mg/mL). IMP form, supplier and
  purity **not printed** [!]; cysteine form (free base vs HCl) not printed (12.1 mg/mL matches free base 121.16 to
  0.1 % [Z]). No sugar added; IMP is the only carbon source besides cysteine.
- **Buffer / pH** [M]: pyrophosphate, pH 3.0 / 4.5 / 6.0; molarity and counter-ion not printed [!]. Pyrophosphate
  pKa ≈ 0.9 / 2.0 / 6.6 / 9.4 (25 °C, textbook) buffers poorly at pH 4.5 [Z]; at 100 mM IMP + 100 mM Cys the
  reactants' own groups (IMP phosphate secondary pKa ~6, Cys carboxyl ~1.9) are a large part of the titratable
  inventory [Z]. pH measured before and after (room temperature implied); drop < 0.2 unit [M].
- **Vessel / headspace / O2** [M]: Pyrex ampoules, 1.0 mL charge, N2 blown over 3 min, flame-sealed; ampoule volume not
  printed.
- **Heating** [M]: CERTOclav autoclave (Kelomat), ampoules horizontal, 140 °C, 1 h, 0.28 MPa (2.7 bar). Saturated steam
  at 140 °C is 0.36 MPa absolute / 0.26 MPa gauge, so 0.28 MPa is consistent with ~141 °C gauge [Z]. Come-up time not
  stated [!].
- **Water activity**: dilute aqueous (≈ 95 % water [Z]).
- **Replicates** [M]: reactions in triplicate; "no less than three headspace collections" per experiment. How the
  +/nd call combined replicates is not stated [!].

## 3. Analytical methods

- **Headspace** [M]: cooled ampoules broken, contents into a 250 mL conical flask with Dreschel head, **20 mL of the
  same pyrophosphate buffer added** "to dilute the mixtures"; flask at 60 °C, agitated; O2-free N2 40 mL/min for 1 h
  onto Tenax-GC (155 mm × 0.75 mm i.d. glass-lined steel trap); trap dried 5 min with N2. Whether one 1 mL ampoule or
  several were pooled per collection is not stated [!].
- **GC-FID / odour port** [M]: thermal desorption at 250 °C onto DB-5 (30 m × 0.32 mm × 1.0 µm) or CP-Wax 52CB
  (50 m × 0.32 mm × 0.21 µm), 0 °C 5 min → 60 °C 5 min → 4 °C/min to 250 °C (DB-5) / 220 °C (Wax); 1:1 split FID /
  sniffing port; four assessors.
- **GC-MS** [M]: HP 5988A / HP 5890, same columns; EI 70 eV, source 200 °C, m/z 29-290, 1.44 scans/s. LRI from C6-C20
  n-alkanes (external, run daily).
- **Identification** [M]: 56 compounds "MS+LRI" against authentic standards, 39 by literature MS only (refs. 1-8 of
  footnote d), 3 by interpretation ("ms": 67, 94, 97). My count of Table 1 reproduces 56 / 39 / 3 [Z]. The text says
  eight compounds "in italics" were synthesised (ref. 13); **no italics are visible in Table 1 at 300 dpi** [!], so
  which eight cannot be told.
- **Quantification: none.** No internal standard is mentioned (the beef papers `madruga1995` / `madruga1997` added
  65 ng dichlorobenzene; this one does not). No peak areas, no ±, no LOD.

## 4. Data

### 4a. Table 1 (pp. 264-267): 98 compounds, presence per pH [M] (read from renders; MS fragment column omitted)

Legend of cells as printed: **+** clear MS · **nd** not detected · **\*** no clear MS, co-eluting. "e" = "compounds
appeared to have meaty and similar aromas" (footnote e). LRI blank = not printed. ID: L = MS+LRI vs authentic; Mn =
literature MS (footnote d source n); ms = interpretation only.

| No | compound (as printed) | pH 6.0 | 4.5 | 3.0 | LRI DB-5 | LRI CP-Wax | ID |
|---|---|---|---|---|---|---|---|
| 1 | sulphur dioxide | + | + | + | – | – | M1 |
| 2 | carbon disulphide | + | + | + | – | – | M1 |
| 3 | 2-butanone | + | + | nd | 603 | | M2 |
| 4 | ethyl acetate | + | + | nd | 612 | | M1 |
| 5 | benzene | + | nd | nd | 664 | 930 | L |
| 6 | 2-pentanone | + | + | + | 689 | 974 | L |
| 7 | 3-pentanone | + | + | nd | 694 | | M1 |
| 8 | 2,3-pentanedione | + | nd | + | 725 | | L |
| 9 | 3-penten-2-one | nd | nd | + | 728 | | M1 |
| 10 | 2-methylthiophene | + | + | + | 777 | 1092 | L |
| 11 | 2,4-pentanedione | nd | + | + | 782 | | L |
| 12 | 2-hexanone | + | + | nd | 798 | 1095 | L |
| 13 | hexanal | + | + | + | 799 | 1097 | L |
| 14 | 3-mercapto-2-butanone | + | + | + | 817 | | L |
| 15 | 2-furfural | + | + | + | 829 | 1457 | L |
| 16 | 2-methylthiazole | + | + | + | 832 | 1250 | L |
| 17 | 2-methylcyclopentanone | + | nd | nd | 836 | 1366 | L |
| 18 | 2,4-dimethylfuran | + | nd | nd | | 1535 | M2 |
| 19 | **2-methyl-3-furanthiol** e | **nd** | **+** | **+** | 870 | 1294 | L |
| 20 | 2-ethylthiophene | + | + | + | 871 | 1170 | L |
| 21 | 4-hydroxy-5-methyl-3(2H)-furanone | nd | + | nd | 880 | | L |
| 22 | 2,5-dimethylthiophene | + | + | + | 882 | | L |
| 23 | 2-heptanone | + | + | + | 890 | 1180 | L |
| 24 | 2-butylfuran | + | + | + | 893 | | L |
| 25 | 3-heptanone | + | nd | nd | 895 | | M1 |
| 26 | 3-mercapto-2-pentanone | nd | + | + | 898 | 1347 | L |
| 27 | 2-acetylfuran | nd | + | + | 904 | 1498 | L |
| 28 | 2-mercapto-3-pentanone | nd | + | + | 910 | | L |
| 29 | **2-furylmethanethiol** e | **nd** | **+** | **+** | 913 | 1420 | L |
| 30 | 3 or 4-methyl-2-furfural | nd | nd | + | 920 | | M5 |
| 31 | dimethylsulphone | + | + | + | 925 | | M1 |
| 32 | 5-methyl-2-furfural | nd | + | + | 949 | 1572 | L |
| 33 | 4,5-dihydro-3(2H)-thiophenone | nd | + | + | 952 | | L |
| 34 | 1-(2-furyl)-2-propanone | nd | nd | + | 954 | 1512 | M2 |
| 35 | benzaldehyde | + | + | + | 958 | 1518 | L |
| 36 | 2-methyl-4,5-dihydro-3-furanthiol e | nd | nd | + | 975 | | M3 |
| 37 | 2- or 3-thiophenethiol e | + | + | + | 977 | 1554 | L |
| 38 | 4,5-dihydro-5-methyl-3(2H)-thiophenone | nd | + | + | 978 | | L |
| 39 | 5-methyl-2-acetylfuran | nd | nd | + | 985 | 1622 | L |
| 40 | 4,5-dihydro-2-methyl-3(2H)-thiophenone | + | + | + | 990 | 1520 | L |
| 41 | 2-pentylfuran | nd | + | nd | 994 | | L |
| 42 | 2,4,5-trimethylthiazole | + | + | + | 996 | 1390 | L |
| 43 | 2-formylthiophene | + | + | + | 1004 | 1687 | L |
| 44 | 4-ethyl-5-methylthiazole | + | + | + | 1005 | 1412 | L |
| 45 | 2-propionylfuran | nd | nd | + | 1019 | | L |
| 46 | 2-ethyl-1-hexanol | + | + | + | 1025 | 1490 | L |
| 47 | 1,3-dithiane | nd | nd | * | 1027 | | M5 |
| 48 | a dimethylfurfural | nd | + | + | 1040 | 1610 | M7 |
| 49 | 5-methyl-3(2H)-thiophenone | nd | * | * | 1042 | | M1 |
| 50 | **2-methyl-3-thiophenethiol** e | **nd** | **+** | **+** | 1065 | 1580 | L |
| 51 | 1-(2-furyl)-3-butanone | nd | nd | + | 1070 | 1640 | M2 |
| 52 | 3-methyl-1,2-dithiolan-4-one | nd | + | + | 1072 | 1664 | L |
| 53 | 2,5-dimethyl-4-ethylthiazole | nd | nd | + | 1074 | | L |
| 54 | 2-thiophenemethanethiol e | nd | nd | + | 1089 | 1689 | L |
| 55 | 2-acetylthiophene | nd | nd | + | 1090 | 1773 | L |
| 56 | a formylmethylthiophene | + | + | + | 1093 | 1712 | M5 |
| 57 | (E or Z) 3,5-dimethyl-1,2-dithiolan-4-one | nd | + | + | 1098 | 1568 | M5 |
| 58 | (E or Z) 3,5-dimethyl-1,2-dithiolan-4-one | nd | + | + | 1111 | 1597 | M5 |
| 59 | 2-formyl-5-methylthiophene | + | + | + | 1118 | 1783 | L |
| 60 | 2-formyl-3-methylthiophene | + | + | + | 1121 | 1815 | L |
| 61 | 1,2,4-trithiolane | * | * | * | 1127 | | M2 |
| 62 | (E) 3,5-dimethyl-1,2,4-trithiolane | + | + | + | 1140 | 1580 | L |
| 63 | 2-acetyl-5-methylthiophene | nd | nd | + | 1145 | 1874 | L |
| 64 | an ethylmethylfurfural | nd | + | + | 1147 | 1675 | M5 |
| 65 | (Z) 3,5-dimethyl-1,2,4-trithiolane | + | + | + | 1151 | 1603 | L |
| 66 | 2-acetyl-3-methylthiophene | nd | nd | + | 1161 | 1762 | L |
| 67 | 3-ethyl-1,2-dithiolan-4-one | nd | + | + | 1163 | | ms |
| 68 | 2-propionylthiophene | nd | nd | + | 1181 | 1840 | L |
| 69 | naphthalene | + | + | + | 1185 | 1730 | L |
| 70 | 3-ethyl-2-formylthiophene | + | + | + | 1203 | | L |
| 71 | 3,(5 or 6)-dimethyl-1,2-dithian-4-one | nd | nd | + | | 1853 | M5 |
| 72 | 2,3-dihydro-6-methylthieno[2,3c]furan | + | + | + | 1209 | | L |
| 73 | 3-methyl-1,2-dithian-4-one | nd | + | + | 1213 | 1896 | M4 |
| 74 | 2-(2-thienyl)-furan | nd | nd | + | 1227 | | M8 |
| 75 | benzothiazole | * | nd | * | 1231 | | M2 |
| 76 | a dimethylformylthiophene | + | + | + | 1240 | 1915 | M5 |
| 77 | 3-methyl-1,2,4-trithiane | + | nd | nd | 1253 | 1862 | L |
| 78 | a thienothiophene | + | + | + | 1257 | 1860 | M2 |
| 79 | a dihydrothienothiophene | + | + | + | 1284 | 1958 | M5 |
| 80 | a methylthienothiophene | nd | + | + | 1322 | 1943 | M5 |
| 81 | 1,2,4,5-tetrathiane | + | + | nd | 1337 | | M2 |
| 82 | a methyldihydrothienothiophene | nd | + | + | 1341 | 2003 | M5 |
| 83 | a methylthienothiophene | nd | + | + | 1364 | 1957 | M5 |
| 84 | a methyldihydrothienothiophene | nd | + | + | 1417 | 2017 | M5 |
| 85 | 2-ethylthienothiophene | nd | + | + | 1422 | 2012 | M2 |
| 86 | 3,4,6-trimethyl-1,2,5-trithiane | + | + | + | 1499 | | M2 |
| 87 | 2-methyl-3-furyl 1-methyl-2-oxopropyl disulphide e | nd | nd | + | 1505 | 2130 | L |
| 88 | **bis(2-methyl-3-furyl) disulphide** e | **nd** | **+** | **+** | 1540 | 2110 | L |
| 89 | 2-methyl-3-furyl 1-ethyl-2-oxopropyl disulphide e | nd | nd | + | 1574 | 2233 | L |
| 90 | 2-methyl-3-furyl 1-methyl-2-oxobutyl disulphide e | nd | + | + | 1596 | 2265 | L |
| 91 | 2-methyl-3-furyl 2-furylmethyl disulphide e | nd | + | + | 1640 | 2347 | L |
| 92 | bis(2-furylmethyl) disulphide e | nd | nd | + | 1683 | | L |
| 93 | 2-methyl-3-furyl 2- or 3-thienyl disulphide e | nd | + | + | 1699 | 2390 | M6 |
| 94 | 2- or 3-thienyl 1-methyl-2-oxobutyl disulphide e | nd | nd | + | 1732 | | ms |
| 95 | 2-methyl-3-furyl 2-methyl-3-thienyl disulphide e | nd | + | + | 1751 | 2395 | M6 |
| 96 | 2-furylmethyl 2- or 3-thienyl disulphide e | nd | + | + | 1789 | | M6 |
| 97 | bis(2- or 3-thienyl) disulphide e | nd | + | + | 1871 | | ms |
| 98 | 2-methyl-3-thienyl 2- or 3-thienyl disulphide e | nd | + | + | 1923 | 2715 | M6 |

(The MS fragment column of Table 1 is not transcribed; it is model-irrelevant. My spellings normalised "thiole" →
"thiol" in rows 29, 36, 37, 50, 54.)

### 4b. Counts by class, detected (+ or *) per pH [Z, my class assignment from Table 1 names]

| class | n | pH 6.0 | pH 4.5 | pH 3.0 |
|---|---|---|---|---|
| thiols (19, 29, 36, 37, 50, 54) | 6 | 1 | 4 | 6 |
| disulfides (87-98) | 12 | 0 | 8 | 12 |
| mercaptoketones (14, 26, 28) | 3 | 1 | 3 | 3 |
| thiophenones (33, 38, 40, 49) | 4 | 1 | 4 | 4 |
| thiophenes (incl. 74) | 14 | 9 | 9 | 14 |
| cyclic polysulfur (dithiolanones, trithiolanes, dithiane/ones, trithianes, tetrathiane) | 13 | 6 | 10 | 11 |
| bicyclic (72, 78-80, 82-85) | 8 | 3 | 8 | 8 |
| thiazoles (16, 42, 44, 53, 75) | 5 | 4 | 3 | 5 |
| furans incl. furfurals, norfuraneol | 14 | 3 | 8 | 11 |
| ketones (3, 6-9, 11, 12, 17, 23, 25) | 10 | 8 | 6 | 5 |
| SO2, CS2, dimethylsulphone | 3 | 3 | 3 | 3 |
| other (ethyl acetate, benzene, hexanal, benzaldehyde, 2-ethylhexanol, naphthalene) | 6 | 6 | 5 | 4 |
| **total** | **98** | **45** | **71** | **86** |

Patterns [Z]: 20 compounds only at pH 3.0; 5 only at pH 6.0 (benzene, 2-methylcyclopentanone, 2,4-dimethylfuran,
3-heptanone, 3-methyl-1,2,4-trithiane); 5 at pH 6.0 and 4.5 but not 3.0 (2-butanone, ethyl acetate, 3-pentanone,
2-hexanone, tetrathiane); 2 at 6.0 and 3.0 but not 4.5 (2,3-pentanedione, benzothiazole *); 2 only at 4.5
(norfuraneol, 2-pentylfuran). Heterocyclics = 98 − 22 = **76**, matching the Conclusion's "seventy-six" [Z].

### 4c. Fig. 1 (p. 263): GC-FID trace, pH 3.0 only, no y-axis scale [F]

Retention order matches the DB-5 LRI order (furfural 15 before MFT 19; 70-72 adjacent), so the trace is DB-5 (as the
text says), while the caption says the aromas summarised are from "the polar column effluent" [!]. Heights as a
fraction of the tallest peaks (±0.05): at or near full scale (tops level, so possibly detector-clipped, P ~60 %):
**19 (MFT), 29 (FFT)**, 36, 37, 62/65, 70/71, 79, 88, 93, and one unlabelled peak just before 19; 50 ≈ 0.70; 82 ≈ 0.55; 87 ≈ 0.54;
59, 91 ≈ 0.51; 97 ≈ 0.49; 95 ≈ 0.44; 89, 90 ≈ 0.40; 96 ≈ 0.33; 14 ≈ 0.15; 94, 98 ≈ 0.13; 35, 92 ≈ 0.11; 15
(furfural) ≈ 0.09. FID response, one run, one pH, no standard: **not ledgerable**; it only shows MFT and FFT were among
the largest headspace peaks at pH 3.0. Odour labels (e.g. MFT region "burnt, rubbery"/"sl. meaty") are descriptive.

## 5. Flags and internal consistency

1. [!] **Presence/absence only.** "+" means "clear MS", so nd vs + is an MS-quality call with no LOD; detection limits
   may differ between pH arms if co-elution differs. No IS, no area.
2. [!] **Abstract contradicts body on pyrazines**: abstract "pyrazines were completely inhibited at high pH" [V];
   body "Absence of pyrazines at low pH" and Conclusion "inhibited by acidic conditions"; **Table 1 has no pyrazine
   at any pH**. The abstract's "high" is an error (~90 %); the body's own claim (pH effect on pyrazines) is untested
   here since none was found even at pH 6.0.
3. [!] p. 263: "at pH 6.0 only one [heterocyclic disulphide] was detected" vs p. 268 and Table 1: **no disulfide at
   pH 6.0**; the one pH-6.0 compound in the meaty-thiol/disulfide group is the thiol 37. The p. 268 count (18 at pH
   3.0, 12 at pH 4.5, 1 at pH 6.0, across the 18 footnote-e compounds) reproduces from Table 1 [Z]; "four [disulfides]
   not detected at pH 4.5" = 87, 89, 92, 94 ✓.
4. [!] Text says 3-mercapto-2-pentanone [26] was "a trace peak" at pH 6.0; Table 1 has nd.
5. [!] Text "3-Methyl-1,2-dithian-4-one [52]": Table 1 row 52 is 3-methyl-1,2-**dithiolan**-4-one; the dithianone is 73.
   Text "2,4,6-trimethyl-1,3,5-trithiane [86]" vs Table 1 "3,4,6-trimethyl-1,2,5-trithiane". Text lists 22
   (2,5-dimethylthiophene) among "formyl thiophenes"; 45 (2-propionylfuran) among "alkylfurans"; calls 16/42/44/53
   "acylthiazoles" (they are alkylthiazoles).
6. [!] **Norfuraneol [21] LRI 880 on DB-5 is suspect**: the same lab's meynier1995 Table 3 gives 4-hydroxy-5-methyl-
   3(2H)-furanone LRI 1042 on a column whose other LRIs match this paper's DB-5 within ±8 (MFT 868/870, FFT
   909/913, 2-formylthiophene 1000/1004, 4,5-dihydro-2-methyl-3(2H)-thiophenone 987/990). Here 1042 is assigned
   to 5-methyl-3(2H)-thiophenone [49] (*, co-eluting). P(row 21 misassigned or its LRI misprinted) ≈ 70 %.
   Also 3-methyl-1,2-dithian-4-one: 1213 here vs 1312 in meynier1995 (isomer or misprint; unresolved).
7. [!] "Eight compounds (in italics) have been synthesised" — no italics survive in the PDF's Table 1.
8. [!] IMP salt form, buffer molarity, ampoule volume, per-pot final pH, come-up time not printed.
9. [!] "total quantity of volatile compounds increased significantly as the pH decreased" [V] — no data or test shown.
10. [!] 2,4-dimethylfuran [18] only at pH 6.0, explained by the authors as possibly lost in headspace capture [Q].
11. Counts check out: 98 rows; 56 L + 39 M + 3 ms = 98 ✓; 8 bicyclics ✓; 5 thiazoles ✓; 76 heterocyclics ✓.
12. Abstract says "over 90 volatile compounds"; Table 1 has 98 — consistent.

## 6. Existing records in the repo that cite this paper

`grep -rIn` over `data/ docs/ src/ scripts/ tasks/` for `S0103-5053`, `10.1590/s0103-50531998`, `madruga1998`,
`Braz. Chem`, `J Braz Chem Soc`, `Madruga & Mottram 1998`, `Madruga and Mottram 1998` (case-insensitive): **the only
hit is the pending queue line** in `data/keys/pdf_dossier_aliases.yml` l. 95:
"A5. Madruga & Mottram 1998 J Braz Chem Soc 9:261: volatiles from heated 5'-IMP + cysteine vs pH (MFT from nucleotide
ribose; E6, pH)".

| field | repo value | paper | verdict |
|---|---|---|---|
| authors / year | Madruga & Mottram 1998 | same | **match** |
| venue | J Braz Chem Soc 9:261 | 9(3):261-271 | **match** |
| system | heated 5'-IMP + cysteine vs pH | 5'-IMP + Cys, pH 3.0 / 4.5 / 6.0, 140 °C 1 h | **match** |
| "MFT from nucleotide ribose" | — | MFT detected at pH 4.5 / 3.0; the ribose route is the authors' [P], ribose never measured | **overstated** — presence only, mechanism not shown |

No ledger entry, no `papers.yml` entry, no conflicts-register, hold-out or directional-claims record cites it. The
Madruga records that do exist (`madruga1995_tab1_tab2_beef_ph_series`, `madruga1997_tab1_beef_imp_thiamine_cysteine`
in sulfur.yml; `madruga1995_tab2_lipid_volatiles_vs_ph_beef` in lipid.yml) are other papers with other pots (sec. 1,
answer 5) and were not re-checked here.

## 7. Proposed ledger entries (NOT written into the ledger)

One entry, `data/lit/measurement_ledger/sulfur.yml` (precedent for a presence-only Reading identification table:
`mottram1995_table1_dione_h2s_identification`). No trunk.yml / aroma_other.yml entry: the non-sulfur rows belong to the
same data set (cross-relevance carried as flags) and are uninformative alone (furfural at every pH, no pyrazine).

```yaml
  - id: madruga1998_tab1_imp_cys_ph_presence
    dossier: madruga1998_extraction.md
    paper: "Madruga & Mottram 1998 J Braz Chem Soc 9:261"
    lab: "Reading (Mottram)"
    system: "5'-IMP 34.8 mg/mL + L-cysteine 12.1 mg/mL (~100 + 100 mM if IMP free acid), 1 mL in N2-flushed flame-sealed Pyrex ampoules"
    medium: "pyrophosphate buffer (molarity not stated) at pH 3.0 / 4.5 / 6.0; CERTOclav autoclave 0.28 MPa"
    temperature_c: [140]
    ph: {initial: [3.0, 4.5, 6.0], measured_during: false, final: "decreased by less than 0.2 unit"}
    time_min: [60]
    species: [MFT, FFT, 2-methyl-3-thiophenethiol, 2-methyl-4,5-dihydro-3-furanthiol, 2-thiophenemethanethiol, bis(2-methyl-3-furyl) disulfide, MFT/FFT/thienyl mixed disulfides, mercaptoketones, thiophenones, thiophenes, trithiolanes, thiazoles, furfural, norfuraneol]
    kind: other
    location: "Table 1 (98 compounds, +/nd/* per pH; LRI DB-5 and CP-Wax; no amounts)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "Cysteine (0.121 g) and 5'-IMP (0.348 g) were dissolved in 10 mL of an appropriate pyrophosphate buffer"
    flags: [presence_only, buffer_molarity_not_stated, internal_inconsistency, no_imp_species_in_engine, "also_relevant: trunk", "also_relevant: aroma_other"]
    note: "Dynamic headspace (60 C, 1 h, Tenax), no internal standard; '+' = clear MS, no LOD. MFT, FFT, 2-methyl-3-thiophenethiol nd at pH 6.0, + at 4.5 and 3.0; disulfides 0/8/12 at pH 6.0/4.5/3.0; no pyrazine at any pH (abstract misstates 'inhibited at high pH'). Bench work at Reading (acknowledgment); same pots probably in the companion sensory paper (Cienc Tecnol Aliment 18(4) 1998) and Madruga's 1994 thesis - key once. IMP hydrolysis (ribose release) not measured, so the pH slope confounds release with sulfur chemistry; pH 3.0 outside the lane's 4.5-7 window. arms: pH 6.0 / 4.5 / 3.0 (same lab, same assay)."
```

## 8. Supplementary material

None referenced: no Supporting Information, appendix or "Table S" anywhere in the 11 pages. Everything model-relevant
is in Table 1. The methods and the eight synthesised reference compounds point to Madruga's PhD thesis (ref. 13:
"Studies on some factors affecting meat flavour formation", University of Reading, 1994, 206 pp.), which would be the
only place quantitative IMP + Cys data from these pots might exist (~30 %); not on disk.
