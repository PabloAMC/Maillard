# Blank, Davidek, Devaud, Sagalowicz, Leser & Michel 2006 (Dev. Food Sci. 43:347) — EXTRACTION (file stem `blank2006`)
### *"Formation of 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol) in structured fluids"* — Weurman symposium proceedings, 4 pp.
### Xylose + glycine or L-leucine (**0.22 M each**) at **70 °C for 7 h**, either in **0.2 M phosphate pH 6.0 (the water control)** or in monoglyceride mesophases (lamellar, cubic, hexagonal; 2 mL buffer in 8 g Dimodan); residual xylose by HPAEC-ECD, norfuraneol (NF) by GC-MS against cyclotene. **No time course is shown**: aliquots were "taken at regular time intervals" [V] but only the **7 h endpoint** is plotted. Water control: **Xyl / Gly 63 % xylose left, NF 750 µg/mmol xylose; Xyl / Leu 80 %, 1520 µg/mmol**.

**Read 2026-11-05 (first reader) from `data/articles/blank2006.pdf`** (4 pp., book-chapter scan with an OCR text layer;
degree signs and µ are garbled in the text layer, read from 150 dpi renders `img/blank06-3.png`, `-4.png`). Bar values are
printed on the bars [M]; nothing needed digitising.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Imre Blank, Tomas Davidek, Stéphanie Devaud, Laurent Sagalowicz, Martin E. Leser, Martin Michel — Nestlé Research Center, Lausanne |
| Venue | W.L.P. Bredie & M.A. Petersen (Eds.), *Flavour Science: Recent Advances and Trends* (Developments in Food Science **43**), Elsevier 2006, pp. **347-350** |
| DOI | not printed; brief gives **10.1016/s0167-4501(06)80082-8** (Crossref not checked this session) |

**Laboratory.** **`nestle_blank`** (in `data/keys/laboratories.yml`; Blank and Davidek both authors). Same lab as
`blank1996` / `blank1997` / `blank2002` / `blank2003` / `davidek2002-2008`. The xylose HPAEC method is cited to Davidek,
Clety, Devaud & Blank 2003 JAFC 51:7259 [C] — that file is on the blinded list and was not opened.

## 1. What this paper answers for the model

| model need | what Blank 2006 gives | tag |
|---|---|---|
| a water-control NF time course | **No.** Only the 7 h endpoint is shown, for buffer and mesophases | [!] |
| xylose loss, 0.22 M Xyl + 0.22 M Gly or Leu, phosphate pH 6, 70 °C, 7 h | **63 % left (Gly), 80 % left (Leu)** | [M] |
| NF yield, same pots | **750 µg/mmol initial xylose (Gly), 1520 (Leu)** | [M] |
| mesophase effect | hexagonal phase: xylose 47 % (Gly) / 53 % (Leu) left, NF 2188 / 2535 µg/mmol; lamellar 45 % / 851; cubic 40 % / 1095 (Gly only) | [M] |
| pH / O2 | phosphate 0.2 M pH 6.0 (initial; not re-measured); Pyrex tube 30 mL with 10 mL → ~20 mL air headspace [Z] | [M] / [!] |

## 2. Conditions [M]

| item | buffer (control) | mesophases |
|---|---|---|
| reactants | D-xylose 0.330 g (2.2 mmol) + Gly 0.165 g or Leu 0.288 g (2.2 mmol) | Xyl 0.066 g (0.44 mmol) + Gly 0.033 g or Leu 0.057 g (0.44 mmol) |
| medium | 10 mL 0.2 M phosphate, pH 6.0 | 8 g monoglyceride (Dimodan U → hexagonal; HR → lamellar; 1:1 → cubic) + 2 × 1 mL 0.2 M phosphate pH 6.0 |
| aqueous concentration | **0.22 M Xyl, 0.22 M amino acid** [Z] | 0.22 M in the 2 mL aqueous domain [Z] |
| preparation | dissolved, heated | **monoglyceride at 116 °C 23 min; + amino acid in buffer (pre-heated 3 min at 116 °C), stirred at 116 °C ~22 min; + xylose in buffer (pre-heated 3 min, 116 °C), 5 min at 116 °C**; cooled |
| reaction | Pyrex tube 30 mL, **70 °C, up to 7 h**, aliquots at intervals | Pyrex tube, silicone bath **70 °C, up to 7 h** |
| replicates | duplicate | duplicate |

**Assays [M].** Xylose: HPAEC-ECD after 250× (buffer) or extraction + 25× (mesophase) dilution. NF: diethyl-ether
extraction (pH 4, NaCl saturated) with **cyclotene 372 µg internal standard**, micro-distillation, GC-MS (DB-Wax 60 m, EI);
**response factor to cyclotene not stated**, no LOD, no recovery. Mesophase samples need an extra back-extraction into
pH 8 phosphate — a different recovery path from the buffer samples [!].

## 3. Data [M] (bars, values printed; 7 h, 70 °C)

| system | xylose left (%) | NF (µg / mmol xylose) |
|---|---|---|
| **Xyl / Gly, buffer** | **63** | **750** |
| Xyl / Gly, lamellar | 45 | 851 |
| Xyl / Gly, cubic | 40 | 1095 |
| Xyl / Gly, hexagonal | 47 | 2188 |
| **Xyl / Leu, buffer** | **80** | **1520** |
| Xyl / Leu, hexagonal | 53 | 2535 |

**Arithmetic [Z].** NF molar yield per initial xylose = (µg/mmol) / 114.1: buffer Gly **0.66 mol %**, Leu **1.33 mol %**;
hexagonal Gly 1.92 %, Leu 2.22 %. Per xylose consumed: buffer Gly 1.8 %, Leu 6.7 %; hexagonal Gly 3.6 %, Leu 4.7 %.
Apparent first-order xylose loss over 7 h (buffer): Gly ln(1/0.63)/420 = **1.1 × 10⁻³ min⁻¹**, Leu **5.3 × 10⁻⁴ min⁻¹**
(second order in practice at 0.22 M each; one point, no order test).

## 4. What it can and cannot test

**Can:** one buffered, 70 °C pentose-amino-acid endpoint for xylose loss and NF yield with two amines (Gly faster
xylose loss, Leu higher NF yield per xylose); the Nestlé lab's NF yield scale at mild temperature (sub-1.5 mol %),
to set against Mikami's 13-16 mol % at 100 °C / pH 6.5 (dilute xylose, excess Lys) — **not a conflict**, very different
T, concentration and amine.

**Cannot:**
- **No time course** — what the brief hoped for is not in the paper; the intermediate aliquots are not reported.
- **NF absolute scale uncertain** (IS response factor not stated).
- **The mesophase contrast is confounded [!]**: mesophase pots spent ~5 min at 116 °C with xylose and amino acid present
  (and ~30 min at 116 °C with the amino acid) before the 70 °C step; the buffer control did not. Xylose loss at 116 °C
  over 5 min could be a sizeable part of the 15-23 % extra loss in the mesophases [P]. The "hexagonal phase accelerates"
  reading needs a pre-heated buffer control, which the paper does not have.
- No pH after heating; no O2 control; no temperature series.

## 5. Conflicts and agreements with the repo's records

1. No existing record cites this paper (grep 2026-11-05 of `data/ docs/ src/ tasks/ scripts/` for "blank2006", the DOI,
   "structured fluid", "Dimodan", "mesophase": "Blank 2006" hits only `davidek2006a/b`, different papers). **Related:**
   `us8920862_extraction.md` (a Nestlé patent, same lab family) heats cysteine + xylose in a Dimodan U / soybean-oil
   mesophase and reports MFT ×10 and FFT ×7.4 over the buffer arm (`lipid.yml`); its mesophase is also prepared hot
   (95 °C). The same preparation confound (sec. 4) applies there in principle — not checked against the patent text.
2. Consistent with the lab's `blank1997` statement that NF needs an amine and near-neutral pH (F13 atlas line 11); adds
   Leu ≥ Gly for NF per xylose at pH 6.

## 6. Proposed ledger entry (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: blank2006_fig1_fig2_xylose_norfuraneol_buffer_70c_7h
    dossier: blank2006_extraction.md
    paper: "Blank, Davidek, Devaud, Sagalowicz, Leser & Michel 2006 Dev Food Sci 43:347"
    lab: "Nestle Research (Blank)"
    system: "D-xylose 0.22 M + glycine or L-leucine 0.22 M (water control arm only)"
    medium: "0.2 M phosphate buffer pH 6.0, 10 mL in a 30 mL Pyrex tube"
    temperature_c: [70]
    ph: {initial: [6.0], measured_during: false, final: null}
    time_min: [420]
    species: [xylose, norfuraneol]
    kind: endpoint
    location: "Figs. 1-2 pp. 349-350 (values printed on bars; dossier sec. 3)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "was heated in a Pyrex tube (30 ml, 70 C, 7 h)"
    flags: [single_time_point, ph_start_only, nf_is_response_factor_not_stated, duplicates_only, no_time_course_despite_sampling]
    note: "Xyl left 63 % (Gly) / 80 % (Leu); NF 750 / 1520 ug per mmol xylose = 0.66 / 1.33 mol %. Mesophase arms not entered (confounded by 116 C preparation)."
```

## 7. Supplementary material

None. A full paper with the time courses may exist in the Nestlé series (not searched).
