# El-Ghorab, Ashraf, Mohamed, Shaaban, El-massry & Farouk 2010 (J. Arab Soc. Med. Res. 5(2):131-139) — EXTRACTION (file stem `elghorab2010`; on disk `el-ghorab2010.pdf`)
### *"The Effect of pH on Flavor Formation and Antioxidant Activity of Amino Acid and Sugars Interaction Products"*
### Cysteine + ribose, cysteine + glucose and cysteine + beef fat, "0.1 M (0.5 g)" Cys + "0.07 M (0.5 g)" sugar in 50 mL water, **unbuffered, adjusted to label pH 4.5 / 7 / 10** with 0.1 N HCl / NaOH, **reflux at boiling, 1 h**. Measured: A420 browning, Folin total phenolics, DPPH and β-carotene antioxidant assays, a sensory profile. **No volatile analysis is described: no GC, no table of volatiles, no MFT number at any pH.** MFT appears once, in the Discussion, as a qualitative claim (present at pH 4.5, absent or traces at high pH).
DOI: none printed; none in the OpenAlex record (M4 case sec. 2.9). Journal web address printed: www.asmr.eg.net

**Read 2026-11-05 (first reader) from `data/articles/el-ghorab2010.pdf`** (9 pp., born-digital, text layer clean on
pp. 131-138; p. 139 is a non-Latin abstract whose text layer is unreadable glyphs and whose render was not read).
Prose from `pdftotext -layout`; **byline and Figs. 1-3 read from 110 dpi renders of pp. 131 and 134** (scratchpad
`m4pdf/el-ghorab2010/p-1.png`, `p-4.png`); Figs. 1 and 2 bar tops read against the gridlines: **read error ±0.02 OD
(Fig. 1), ±5 µg/mL (Fig. 2)**. Figs. 3-5 (% inhibition vs extract concentration) and Fig. 6 (sensory radar) were not
digitised: no engine observable.

**Reservation (M4, case `docs/reference/cases/M4_glucose_mft_protein_amine.md`, branch `m4-case-refresh`, commit
8660c005, sec. 2.9).** Role fixed before reading: D2 (second laboratory) only if MFT is reported at pH 4.5 / 7 / 10.
The paper's MFT statement and its Cys + glucose pots carry **RESERVED: unseen by the engine (M4); no engine run on
these pots until M4's design is frozen** (**[RESERVED M4]**). The role is not met (sec. 8). No engine run.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[Z]** derived by me · **[!]** flag. Provenance:
**first reader only, second reader pending.** Nothing here is in a ledger, fit table or benchmark.

## 0. Identity

| field | value as printed (p. 131) |
|---|---|
| Authors | **Ahmed H El-Ghorab**¹*, **Iram F Ashraf**², **Anjum F Mohamed**², **Hamdy A Shaaban**¹, **Khaled F El-massry**¹, **Amr Farouk**¹ |
| Affiliations | ¹Flavor and Aromatic Dept., National Research Center, Egypt; ²National Institute of Food Science and Technology, Agric. Uni. Faisal Abad, Pakistan; corresponding e-mail aghorab2@yahoo.com |
| Venue | *JASMR* (Journal of the Arab Society for Medical Research) **5** (2): **131-139** (2010) |
| DOI | **none printed** |

**Byline vs indexes [!].** Six authors on the PDF. The M4 case's candidate row (from OpenAlex) lists five: "El-Ghorab,
Ashraf, Mohamed, El-Massry & Farouk", without **Shaaban**. Other indexes transliterate the names differently (al-Ghurab,
Muhammad, Shaban, al-Misri, Mansur, per the brief). The byline above is the record; the key is `elghorab2010` after the
first author as printed.

**Laboratory.** National Research Center, Cairo, Flavor and Aromatic Dept. (El-Ghorab): the laboratory of
`elmassry2003_extraction.md` (ledger lab string "NRC Cairo (El-Ghorab)"), with co-authors from Faisalabad. Not DFA.

## 1. What this paper answers for the model

| model need | what El-Ghorab 2010 gives | tag |
|---|---|---|
| Glucose + Cys MFT at >= 3 pH (D2) | **none quantified**; one sentence: MFT "excellent presence" at pH 4.5, "absent at high pH or present in traces" | [Q] |
| Any volatile, any pH | **none**: no GC method, no table, no chromatogram | — |
| Browning, total phenolics, antioxidant activity vs pH | Figs. 1-5 and text values (sec. 4) | [M] / [F] |
| Rate constants, Ea, time course | **none** | — |

## 2. Conditions [M] (p. 132)

| item | as printed | derived [Z] |
|---|---|---|
| Charge | "0.1M (0.5g) cystiene amino acid, 0.07M (0.5g) ribose/glucose and ground beef fat sugar" weighed in a 50 mL beaker, mixed in distilled water, volume made to **50 mL** | **0.5 g Cys / 50 mL = 0.083 M (not 0.1); 0.5 g glucose = 0.056 M, 0.5 g ribose = 0.067 M (not 0.07)** [!]. Beef-fat amount not stated |
| Materials | the Methods list "glutamine and glutamic acid" as the proteins used; no such pot appears anywhere [!] | — |
| pH | adjusted to **4.5, 7 and 10** with 0.1 N HCl / NaOH before making to volume: **label pH, room temperature; no buffer; pH after heating not reported** | — |
| Heating | transferred to "flask of Reflux apparatus" and heated "at boiling temperature for one hour" | ~100 °C, open to the condenser (air headspace) [Z] |
| Work-up | "solvent-solvent" extraction into dichloromethane, rotary evaporation, stored at 5 °C | no internal standard, no analysis of the extract is described |
| Replicates | phenolics "in triplicate"; antioxidant values given as mean ± SD (n not stated) | — |

## 3. The assays [M] (pp. 132-133)

| assay | method | LOD / noise |
|---|---|---|
| Browning | A420 after "appropriate dilution" (dilution factor not stated), UV-1601 [C, Lertittikul et al.] | none stated; Fig. 1 bars carry small error bars |
| Total phenolics | Folin-Ciocalteu, 760 nm, gallic acid standard, µg/mL, triplicate | e.g. 409.13 ± 12.96 µg/mL |
| Antioxidant | DPPH (517 nm) and β-carotene / linoleic acid bleaching (470 nm), % inhibition at 10-500 µg/mL, BHT reference 91.135 ± 0.05 % | SD 0.2-0.5 % |
| Sensory | "according to El-massry, Farouk and El-Ghorab" [C]; odour intensity / acceptability, Fig. 6 | panel size not stated |
| **Volatiles** | **not described** (no GC, no MS, no identification, no quantification) | — |

## 4. Data

### 4.1 MFT and other volatiles: qualitative claims only [Q] [RESERVED M4]

- Discussion (p. 136): furan derivatives "especially 2-methyl-3-furanthiol" are "mostly dectected [sic] at acidic pH
  values between 4.5-5.5", followed by a citation to Nursten; the sentence is literature, not this paper's data [C].
- Discussion (p. 137): the glucose-cysteine system has a "pronounced flavor of meat at acidic pH (4.5) (fig. 2)" due to
  "excellent presence" of MFT, "which is absent at high pH or present in traces". **The figure cited (Fig. 2) is total
  phenolics, not a volatile measurement** [!]. No MFT value, peak area or detection method exists in the paper.
- Fat-cysteine pots: thiophene, trithiolane, thiazole, "thiophene faranmethyl thiol" [sic] and "2-pentyl 2(H) thiafuran"
  said to dominate (p. 136) [Q], again with no data.

### 4.2 Fig. 1 (p. 134): browning, optical density at 420 nm [F ±0.02]

| pot | pH 4.5 | pH 7 | pH 10 |
|---|---:|---:|---:|
| Cysteine / fat | 1.19 | 1.13 | 0.89 |
| Ribose / cysteine | 0.87 | 0.93 | 1.02 |
| **Glucose / cysteine [RESERVED M4]** | 0.80 | 0.90 | 0.76 |

**[!] The text contradicts the figure:** "browning degree was highest at pH 7 is 2.060" and "lowest ... at pH 4.5 in
fat model system is 0.692" (p. 133); the figure shows no bar above 1.2 and the fat pH 4.5 bar as the **highest**. The
text also says browning in fat pots rises with pH (the figure shows it falling). The dilution factor is not given, so
neither set is an absolute A420.

### 4.3 Fig. 2 (p. 134): total phenolics, µg/mL gallic acid equivalents [F ±5; text values M]

| pot | pH 4.5 | pH 7 | pH 10 |
|---|---:|---:|---:|
| **Glucose / cysteine [RESERVED M4]** | 408 (text **409.13 ± 12.96**) | 352 | 272 |
| Ribose / cysteine | 188 | 229 | 158 |
| Cysteine / fat | 12 | 66 | 52 (text "**67.47 ± 3.33**" at pH 10) |

**[!]** The text puts the fat-pot maximum (67.47) at pH 10 and orders fat phenolics 10 > 7 > 4.5; the figure has the
maximum at pH 7 (~66) and pH 10 at ~52. The text value matches the figure's pH 7 bar.

### 4.4 Antioxidant activity at 500 µg/mL (text, p. 133) [M]

| pot | highest | lowest |
|---|---|---|
| Cys + ribose | pH 7: **52.98 ± 0.38 %** | pH 4.5 (value not printed) |
| Cys + glucose [RESERVED M4] | pH 4.5: **92.62 ± 0.21 %** | pH 10: **36.76 ± 0.45 %** |
| Cys + fat | pH 10 ("significant"; value not printed) | — |
| BHT reference | 91.135 ± 0.05 % | — |

## 5. What it can and cannot test in the model

**Cannot:** anything on MFT, FFT, NF or any volatile (no measurement); any rate, Ea or pH law; A420 and phenolics are
not engine observables (melanoidin colour is not keyed to A420 in the engine; no species carries "total phenolics").

**Can:** nothing quantitative for M4. At most, a laboratory's unsupported qualitative statement that glucose + Cys MFT
falls from pH 4.5 to pH 10 at ~100 °C (reflux, 1 h), which would agree in sign with Hofmann 1998 Table 2's 5 → 7 limb
(`hofmann1998_reconciliation.md` sec. 2.2) if it were a measurement. It is not.

## 6. Conflicts and agreements with existing records

- **[!] Internal contradictions:** stated vs actual molarities (sec. 2); browning text vs Fig. 1 (sec. 4.2); phenolics
  text vs Fig. 2 (sec. 4.3); an MFT claim supported by a citation to a phenolics figure (sec. 4.1); Methods list
  glutamine and glutamic acid that are never used.
- **[!] Authorship record:** OpenAlex drops Shaaban (sec. 0).
- **Same laboratory as `elmassry2003_extraction.md`** (NRC Cairo): that paper is GSH + ribose volatiles with a GC
  method; this one has none.
- No conflict with any ledger value (nothing comparable is measured).

## 7. Provenance tiers (first reader only)

MFT statement: `qualitative_claim_without_data`, confidence none for any quantitative use. Figs. 1-2: `figure_digitised`,
low (text contradicts them). Text antioxidant values: `direct_measurement`, not an engine observable. **Proposed ledger
entries: none.**

## 8. M4 role (case 8660c005)

**Role fixed before reading** (case sec. 2.9, candidate table): "as Yeo 1991 (D2, second laboratory)", i.e. decisive
for D2 only if it is a glucose + Cys pot with MFT quantified at >= 3 pH; "peak areas at fixed design suffice for a
ratio". The brief restates it: D2 only if MFT is reported at pH 4.5 / 7 / 10.

| element of the criterion | what the paper has | verdict |
|---|---|---|
| a glucose + Cys pot at >= 3 pH, one design | yes: 0.083 M Cys + 0.056 M glucose (as computed from the printed masses), label pH 4.5 / 7 / 10, reflux 1 h | met |
| MFT reported at each pH (a number or a peak area) | **no**: no volatile analysis is described; MFT appears only as a qualitative Discussion claim ("excellent presence" at 4.5, "absent ... or ... traces" at high pH) with no value and no method | **not met** |

**Verdict: NOT MET.** No MFT number, peak area or detection method exists in the paper. D2 keeps one laboratory.
