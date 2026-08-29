# Blank, Fay, Lakner & Schlosser 1997 (`10.1021/jf960997i`) — single-paper extraction, Wave K5b, 2026-08-29

> Blank, I.†; Fay, L. B.†; Lakner, F. J.‡,§; Schlosser, M.*,‡
> "**Determination of 4-Hydroxy-2,5-dimethyl-3(2H)-furanone and 2(or 5)-Ethyl-4-hydroxy-5(or
> 2)-methyl-3(2H)-furanone in Pentose Sugar-Based Maillard Model Systems by Isotope Dilution
> Assays**"
> *J. Agric. Food Chem.* **1997**, 45 (7), 2642−2648.
> † Nestlé Research Center, Nestec Ltd., Vers-chez-les-Blanc, P.O. Box 44, CH-1000 Lausanne 26,
> Switzerland. ‡ Institut de Chimie Organique, BCh, Université de Lausanne, CH-1015 Lausanne.
> § Present address: Dept. of Biochemistry, University of Illinois, Urbana, IL 61801.
> Received for review December 27, 1996. Accepted March 24, 1997. Abstract published in Advance
> ACS Abstracts, **June 1, 1997**. Article ID `JF960997I`. Grant: Schweizerische Nationalfonds
> 20-36385-92.

**Source PDF:** `data/articles/blank1997.pdf` (7 pp., 259,214 bytes, on disk 2026-08-29 09:39).
**Read method:** full text layer, `pdftotext -layout` **and** `pdftotext -raw` (both runs cross-checked
against each other for every table cell, because the `-layout` pass mis-wrapped one cell of Table 3).
All five tables came out of the text layer intact. Figures 1−4 (a synthesis scheme and three EI mass
spectra) have no numeric text layer and nothing was digitised from them.

---

## 0. IDENTITY VERDICT — ★ THE FILE IS THE RIGHT PAPER, AND IT IS THE ONE THE REPO HAS BEEN MISSING

| Brief expected | Found |
|---|---|
| Blank, Fay, Lakner & Schlosser 1997, JAFC 45:2642−2648 | **EXACT MATCH** — title, all four authors, journal, volume, issue, year and page range all confirmed from the paper's own header and footer. |
| The isotope-dilution QUANTITATIVE companion to `blank1996.pdf` | **CONFIRMED.** It cites Blank & Fay 1996 as its own precursor and delivers precisely the "amounts … being quantified and parameters affecting this reaction … currently being studied" that the 1996 paper promised in its last Discussion sentence. |
| DOI `10.1021/jf960997i` | **NOW CERTAIN.** The article ID printed in the footer is `JF960997I`, and the first-page CCC line is `S0021-8561(96)00997-1`. This **upgrades `blank1996_extraction.md` §0a / §12 item 1 from "~85 % confidence, web lookup only, not CrossRef-verified" to VERIFIED FROM THE PRIMARY DOCUMENT.** |

**No wrong-file problem. This paper is exactly what it says it is.**

### 0a. WHAT IT DELIVERS — and it is the single richest DMHF file in the corpus

**Five tables, 39 quantified cells, all in µg per mmol of sugar, all by stable-isotope dilution
assay against synthesised [¹³C₂]HDMF and [²H₃]EHMF, all in duplicate with a stated maximum SD.**
That is: a **sugar series**, an **amino-acid series**, a **buffer/medium series**, a **pH ladder**,
and a **precursor-ratio dose−response**. There is no rate constant and no time course (§7), but
every other axis the DMHF channel needs is here.

### 0b. ★★ IT CONTRADICTS TWO OF THE FOUR CONSTRAINTS `blank1996_extraction.md` PROPOSED AS HOLD-OUTS

Read this before anything is declared. Both contradictions run the same way: **a GC-O `−` in the
1996 paper is not a zero, it is a non-detection at a sniffing port, and the 1997 IDA sees through it.**

| `blank1996_extraction.md` item | Its proposed role | What THIS paper measures | Verdict |
|---|---|---|---|
| **#7** — "HEMF **requires alanine**; absent with glycine, absent from pentose alone" | §10.2: **HOLD-OUT, structural**, "a model that emits HEMF from pentose+glycine fails" | Table 1: **EHMF = 0.3 (xylose/Gly), 0.7 (ribose/Gly), 1.3 (arabinose/Gly) µg/mmol.** Non-zero in every glycine system. | **★ CONTRADICTED AS AN ON/OFF SWITCH.** A model that emits HEMF from pentose+glycine is **right**, not wrong. The claim survives only in demoted form: a **preference of ≈10−25×**, not a truth table. **The proposed hold-out must be rewritten or it will fail a correct model.** |
| **#8** — "DMHF forms from **pentose alone**, no amino acid" | §10.2: **HOLD-OUT, structural**, "a zero-amino-acid positive control" | Table 1 footnote c: *"The control sample (without amino acid) resulted in **less than 0.01 µg of HDMF and EHMF per mmol of sugar**."* Table 3 rows 1−2: xylose in water **<0.01**, xylose in phosphate **<0.01**. | **★ CONTRADICTED AS A POSITIVE CONTROL.** The pentose-alone HDMF is **≥260× below** the xylose/glycine value and is reported only as an upper bound. It must be re-specified as a **CEILING (<0.01 µg/mmol)**, which is a *negative* control — the opposite role. |

The other two 1996 constraints **survive and are corroborated** — see §8.

---

## 1. SYSTEM DEFINITION — verbatim (Experimental Procedures, "Quantification Experiments", p. 2644)

> "Sample preparation was performed as recently described (Blank and Fay, 1996) with some
> modifications. **The precursors were allowed to react in water and phosphate- and
> malonate-buffered solutions (0.2 M, Na₂HPO₄, Na₂H₂C₃O₄) at different pHs (5, 6, 7) at 90 °C for
> 1 h.** After the reaction mixture was cooled rapidly, water (100 mL) and the labeled internal
> standards (**9.6−48.2 µg of [¹³C₂]HDMF and 10.0−50.2 µg of [²H₃]EHMF**) were added. The solution
> was saturated with NaCl, and **the pH was adjusted to 4.0 (aqueous HCl, 2 M) which is the pH
> optimum for the stability of HDMF in aqueous solutions (Hirvi et al., 1980)**. Neutral compounds
> were continuously extracted with Et₂O overnight using a rotation perforator. The organic phase
> was dried over Na₂SO₄ at +4 °C and concentrated to 0.5−1 mL. **All experiments were performed in
> duplicate.**"

Inherited from Blank & Fay 1996 (the "as recently described" base protocol, `blank1996_extraction.md`
§1.2): **5 mmol sugar + 5 mmol amino acid in 5 mL buffer**, i.e. 1.0 M : 1.0 M, sealed screw-cap
Pyrex tube, oil bath, magnetic stirring.

Reagents: D-xylose, D-ribose, D-arabinose, glycine, L-alanine **and 4-aminobutyric acid** (all ≥99 %,
Fluka). `[1,2-¹³C₂]`Ethanal and `[2,2,2-²H₃]`ethyl bromide from Dr. Glaser AG, Basel.

### 1.1 ⚠️ THE pH LABELS: WHAT IS AND IS NOT KNOWN

- **Table 4's pH 5 / 6 / 7 are buffer set-points in 0.2 mol/L phosphate.** The paper does **not**
  state whether pH was held during the hour, and gives **no final pH** for any run.
- The base protocol it inherits (Blank & Fay 1996) explicitly reports **drift from 6.0 to 5.0
  (xylose/Gly) and 5.3 (xylose/Ala)** under the identical buffer.
- Table 3 footnote c marks only the **water** rows as *"The pH was not controlled during the
  reaction."* The natural reading is that this footnote is there because the buffered rows **were**
  controlled — but **that is an inference, not a statement**, and the buffer capacity of 0.2 M
  phosphate against 1 M sugar + 1 M amino acid at 90 °C for 1 h is not demonstrated anywhere in
  either paper.
- **Carry this exactly as written: "pH 5/6/7, 0.2 M phosphate, set-point; hold not stated; the
  same buffer drifted 1.0 pH unit in the companion paper."** Same defect class as `k3` §C.11
  (Zhou 2023) and §1.2 of `blank1996_extraction.md`.

### 1.2 Analytics

- **GC-MS**: HP-5971 MSD + HP-5890 + HP-7673 autosampler. Splitless injector, **Carbowax capillary
  column**, chromatographic conditions per Blank et al. 1996b. Interface 220 °C, EI 70 eV,
  source ≈180 °C.
- **SIM ions** (Figure 2 legend, verbatim): *"the molecular ions **m/z 128 of HDMF, m/z 142 of EHMF
  (two tautomers), m/z 130 of [¹³C₂]HDMF, and m/z 145 of [²H₃]EHMF (two tautomers)**."*
- **Calibration** (p. 2645): *"The calibration curves were established with standard mixtures
  containing defined amounts of labeled and unlabeled compound in different ratios following the
  procedure described by Guth and Grosch (1990). **A good linearity was found in the concentration
  range of 3−50 µg/mL (r² = 0.999).** Samples for establishing the calibration curves and for
  quantifying HDMF and EHMF in the Maillard model reactions were **injected twice**."*
- **Qualitative MS** for the synthetic standards: Finnigan MAT 8430, EI 70 eV / CI 150 eV with
  ammonia, DB-FFAP 30 m × 0.25 mm × 0.25 µm, 50 °C (2 min) → 4 °C/min → 180 °C → 10 °C/min →
  240 °C (10 min).
- **NMR**: Bruker AC-250 (¹H, 250 MHz) and AMX-400 (¹³C, 100 MHz), CDCl₃. IR: Perkin-Elmer 1420.

**★ This is a real SIDA.** Internal standard added **before** workup, so the overnight ether
perforation's unknown recovery — the exact defect that made Blank & Fay 1996 §3.3 unusable — is
**corrected out**. The isotopomers of one compound share a response factor exactly. **These numbers
are ingestible in a way the 1996 "low mg/kg" was not.**

---

## 2. ★ TABLE 1 — verbatim (p. 2646). The sugar × amino-acid grid.

> **Table 1. Formation of HDMF and EHMF from Pentose Sugars in Maillard Model Systems Containing
> Glycine or L-Alanine^a**

| model system | HDMF^b | EHMF^b |
|---|---|---|
| arabinose/glycine | **5.1** | **1.3** |
| xylose/glycine^c | **2.6** | **0.3** |
| ribose/glycine | **3.6** | **0.7** |
| arabinose/L-alanine | **1.2** | **6.8** |
| xylose/L-alanine^c | **0.9** | **7.5** |
| ribose/L-alanine | **1.6** | **10.0** |

> ^a Equimolar amounts of precursors were used. Reaction conditions: **phosphate buffer (0.2 M,
> pH 6), 90 °C, 1 h**.
> ^b Data are means of **at least two assays, each of them injected twice; maximum SD ≤ 10 %**
> (in **µg/mmol of sugar**).
> ^c **The control sample (without amino acid) resulted in less than 0.01 µg of HDMF and EHMF per
> mmol of sugar.**

**Abstract restatement, verbatim:** *"The yields obtained were **2.6−5.1 µg of HDMF and 6.8−10 µg of
EHMF per mmol pentose**."* (Note the abstract quotes the HDMF range **from the glycine systems only**
and the EHMF range **from the alanine systems only** — it is a range of the *preferred* pairings,
not of the whole table.)

### 2.1 What Table 1 settles that Blank & Fay 1996 could not

1. **★ THE THREE PENTOSES ARE RESOLVED.** `blank1996_extraction.md` §11 records as a declared gap:
   *"The three pentoses are never resolved … No sugar-reactivity ordering can be extracted from
   this paper."* **That gap is now closed, and the lumping hid a genuine 2× spread — with a
   DIFFERENT ordering for each product:**

   | product | ordering | spread |
   |---|---|---|
   | **HDMF** (with glycine) | **arabinose (5.1) > ribose (3.6) > xylose (2.6)** | 1.96× |
   | **EHMF** (with alanine) | **ribose (10.0) > xylose (7.5) > arabinose (6.8)** | 1.47× |

   **The orderings are not the same and are close to reversed** (arabinose is best for HDMF and
   worst for EHMF; ribose is best for EHMF and middling for HDMF). Any model that treats "pentose
   reactivity" as one scalar cannot reproduce both columns.

2. **The amino-acid selectivity switch, quantified.** HDMF : EHMF ratio —

   | system | HDMF/EHMF |
   |---|---|
   | arabinose/Gly | 3.9 |
   | xylose/Gly | **8.7** |
   | ribose/Gly | 5.1 |
   | arabinose/Ala | 0.18 |
   | xylose/Ala | **0.12** |
   | ribose/Ala | 0.16 |

   **★ The xylose pair swings the selectivity ratio by 72× (8.7 → 0.12) on the amino acid alone.**
   That is the sharpest single number in the paper and it is a pure ratio, so it is immune to any
   absolute-scale error.

3. **The zero-amino-acid control is a CEILING, not a positive.** Footnote c, `<0.01 µg/mmol` for
   both products. See §0b.

### 2.2 [DERIVED — arithmetic only, flagged] mol % on sugar

MW(HDMF, C₆H₈O₃) = 128.13; MW(EHMF, C₇H₁₀O₃) = 142.15. Basis = 1 mmol sugar.
mol % = µg × 0.1 / MW.

| system | HDMF µg/mmol | **HDMF mol %** | EHMF µg/mmol | **EHMF mol %** |
|---|---|---|---|---|
| arabinose/Gly | 5.1 | **0.0040** | 1.3 | 0.00091 |
| xylose/Gly | 2.6 | **0.0020** | 0.3 | 0.00021 |
| ribose/Gly | 3.6 | **0.0028** | 0.7 | 0.00049 |
| arabinose/Ala | 1.2 | 0.00094 | 6.8 | **0.0048** |
| xylose/Ala | 0.9 | 0.00070 | 7.5 | **0.0053** |
| ribose/Ala | 1.6 | 0.0012 | 10.0 | **0.0070** |
| control, no AA | <0.01 | **<7.8 × 10⁻⁶** | <0.01 | <7.0 × 10⁻⁶ |
| **max in paper** (Table 5, xylose/Ala 1:4) | 1.6 | 0.0012 | **20.0** | **0.0141** |

**★ CROSS-CHECK ON A PRIOR DERIVATION.** `blank1996_extraction.md` §3.3 derived, from the 1996
paper's "low mg/kg" aside, that the expected scale was *"of order 10⁻³ − 10⁻² mol % on sugar,
spanning ~50×, on an ambiguous basis"* — and wrote it down explicitly so it could be sanity-checked
when the companion arrived. **The companion has arrived and the measured range is
7 × 10⁻⁴ − 1.4 × 10⁻² mol %.** The estimate was right, and the "reaction tube, 6.13 g" basis row of
that table (0.001−0.005 mol %) is the one that matched. Recording this because it validates the
method of that dossier, not because the derived numbers should now be ingested — **they still
should not; use the measured ones on this page instead.**

---

## 3. ★ TABLE 2 — verbatim (p. 2646). The 4-aminobutyric-acid experiment: a Strecker-null control.

> **Table 2. Formation of HDMF and EHMF in Maillard Model Systems Containing Xylose and Different
> Amino Acids^a**

| model system | HDMF^b | EHMF^b |
|---|---|---|
| xylose/4-aminobutyric acid | **0.4** | **0.1** |
| xylose/4-aminobutyric acid/glycine | **1.5** | **0.2** |
| xylose/4-aminobutyric acid/L-alanine | **0.7** | **3.2** |

> ^a Equimolar amounts of precursors were used. Reaction conditions: **phosphate buffer (0.2 M,
> pH 6.0), 90 °C, 1 h**.
> ^b Data are means of at least two assays, each of them injected twice; **maximum SD ≤ 10 %**
> (in µg/mmol of sugar).

**The design, verbatim (p. 2646):**
> "**4-Aminobutyric acid was employed to study the formation of 3(2H)-furanones by sugar
> fragmentation. As a γ-amino acid, it does not decompose by Strecker deamination.** The amounts of
> HDMF and EHMF formed were relatively low, i.e. 0.4 and 0.1 µg, respectively (Table 2). The data
> indicate that **the generation of HDMF by sugar fragmentation is favored compared to that of EHMF,
> most likely due to higher amounts of reactive C3 fragments formed, whereas EHMF requires both
> reactive C3 and C4 fragments.** As shown in Table 2, the amount of HDMF was significantly
> increased in the presence glycine (1.5 µg). Similarly, addition of L-alanine resulted in higher
> yields of EHMF (3.2 µg)."

**★ THIS IS A NON-ISOTOPIC, INDEPENDENT MEASUREMENT OF THE STRECKER/FRAGMENTATION SPLIT.** GABA is
a γ-amino acid: it accelerates the Maillard reaction as an amine but **cannot** donate a Strecker
aldehyde. So `xylose/GABA` is the fragmentation-only baseline, and the increment on adding a
Strecker-active amino acid at the same molarity is the Strecker channel.

**[DERIVED — arithmetic on Table 2 only]:**

| product | fragmentation-only (GABA) | + Strecker-active AA | Strecker increment | **Strecker share** | **fragmentation share** |
|---|---|---|---|---|---|
| **HDMF** (glycine added) | 0.4 | 1.5 | 1.1 | **73 %** | **27 %** |
| **EHMF** (alanine added) | 0.1 | 3.2 | 3.1 | **97 %** | **3 %** |
| HDMF (alanine added — cross-term) | 0.4 | 0.7 | 0.3 | 43 % | 57 % |

**★★ THE HDMF FIGURE IS 73/27, AND BLANK & FAY 1996 MEASURED 70/30 BY ¹³C LABELLING.** Two
completely independent methods — an isotopomer distribution and an amino-acid substitution — on the
same system, agreeing to within 3 percentage points. **This is the strongest corroboration in the
whole DMHF cluster** and it means the 70/30 split is safe to use as a branch prior.

Blank's own alternative framing of the same split (p. 2646, verbatim):
> "As indicated in Table 1, both furanones were also generated by sugar fragmentation, i.e.
> **0.9−1.6 µg of HDMF and 0.3−1.3 µg of EHMF** detected in the model reactions pentose/L-alanine
> and pentose/glycine, respectively. This, however, is **a minor reaction pathway where the amino
> acid indirectly contributes to the formation of these furanones by accelerating the Maillard
> reaction**, i.e. by favoring sugar fragmentation and condensation of the corresponding sugar
> fragments."

**[DERIVED] third estimate, from that framing:** if the HDMF seen in the *alanine* systems is all
fragmentation, the fragmentation share of HDMF in the glycine systems is 0.9/2.6 = **35 %**
(xylose), 1.2/5.1 = **24 %** (arabinose), 1.6/3.6 = **44 %** (ribose).

**Consolidated band for the fragmentation share of HDMF in a pentose/glycine system at 90 °C:**
**15 % (GABA-vs-xylose/Gly, 0.4/2.6) − 44 % (ribose alanine-baseline)**, with the two best-designed
estimates — the GABA increment (27 %) and the 1996 isotopomer measurement (30 %) — clustered
tightly near **30 %**.

### 3.1 The C3-fragment lead — flagged because it is the bridge to Wang & Ho 2008

> "Several reactive sugar fragments have been reported in the literature (review by Ledl and
> Schleicher, 1990). **First trials have shown that some of the well-known C3 fragments do generate
> HDMF, e.g. methylglyoxal, acetol, and dihydroxyacetone (data not shown). These results will be
> reported elsewhere.**"

**Blank's lab states, in 1997, that methylglyoxal generates HDMF — as an unpublished "data not
shown".** Eleven years later Wang & Ho 2008 (`10.1021/jf8012025`, dossiered in this wave) publish
exactly that measurement. The MGO → DMHF edge therefore has **two independent sources that agree on
its existence**, one of which is the same lab that built this pentose channel. **It is not an
isolated claim.**

---

## 4. ★ TABLE 3 — verbatim (p. 2647). The reaction-medium series. **The largest lever in the paper.**

> **Table 3. Effect of Reaction Medium on Formation of HDMF and EHMF from Xylose in Maillard Model
> Systems^a**

| model system | reaction medium | HDMF^b | EHMF^b |
|---|---|---|---|
| xylose | water^c | **<0.01** | **<0.01** |
| xylose | phosphate | **<0.01** | **<0.01** |
| xylose/glycine | water^c | **0.06** | **<0.01** |
| xylose/glycine | phosphate | **2.6** | **0.3** |
| xylose/glycine | malonate | **0.6** | **0.1** |
| xylose/L-alanine | water^c | **0.02** | **0.05** |
| xylose/L-alanine | phosphate | **0.9** | **7.5** |
| xylose/L-alanine | malonate | **0.3** | **0.8** |

> ^a Equimolar amounts of precursors were used. Reaction conditions: **pH 6, 90 °C, 1 h**.
> ^b Data are means of at least two assays, each of them injected twice; maximum SD ≤ 10 %
> (in µg/mmol of sugar).
> ^c **The pH was not controlled during the reaction.**

*(⚠️ The `pdftotext -layout` pass split the `xylose/glycine · water` row across two lines and put the
`<0.01` EHMF cell on its own line. The `-raw` pass renders it as one row, `xylose/glycine water^c
0.06 <0.01`, and the row arithmetic confirms that reading. The table above is the `-raw` reading.)*

**Authors' reading, verbatim (p. 2647):**
> "As shown in Table 3 the furanones were **preferentially formed in the phosphate-buffered
> systems**, i.e. 2.6 µg of HDMF in xylose/glycine and 7.5 µg of EHMF in xylose/alanine, **thus
> indicating a catalytic effect of phosphate, particularly in the presence of the amino acids**.
> The effect of malonate was less pronounced, yielding 0.6 µg of HDMF in xylose/glycine and 0.8 µg
> of EHMF in xylose/alanine. **The model systems without buffer gave rise to less than 0.1 µg of
> HDMF and EHMF per mmol of pentose.**"

**[DERIVED] catalysis factors — this is a HUGE effect and the repo has no term for it:**

| comparison | HDMF | EHMF |
|---|---|---|
| xylose/Gly: **phosphate / water** | 2.6 / 0.06 = **43×** | ≥30× (0.3 / <0.01) |
| xylose/Ala: **phosphate / water** | 0.9 / 0.02 = **45×** | 7.5 / 0.05 = **150×** |
| xylose/Gly: **phosphate / malonate** | 2.6 / 0.6 = **4.3×** | 3× |
| xylose/Ala: **phosphate / malonate** | 3× | 7.5 / 0.8 = **9.4×** |

**★ 0.2 M phosphate buys 43−150× over unbuffered water at the same nominal pH 6, and 3−9× over
0.2 M malonate at the same ionic strength and the same pH.** The malonate control is what makes
this a *specific* phosphate catalysis claim rather than a buffering or ionic-strength artefact.
Note the honest caveat the authors themselves flag: the water rows had **uncontrolled pH**, so part
of the 43−150× is a pH-drift effect and not all of it is catalysis. The **phosphate-vs-malonate**
comparison (3−9×) is the clean one — both buffered, both at pH 6.

**This is a phosphate-CONCENTRATION-free measurement: 0.2 M is the only level tested. There is no
phosphate dose−response in this paper** (same single-point limitation as Blank & Fay 1996 §1.2).

---

## 5. ★★ TABLE 4 — verbatim (p. 2647). The pH ladder. **The most model-relevant table here.**

> **Table 4. Effect of pH on Formation of HDMF and EHMF from Xylose in Maillard Model Systems
> Containing Glycine or L-Alanine^a**

| model system | pH | HDMF^b | EHMF^b |
|---|---|---|---|
| xylose/glycine | **5** | **2.6** | **0.3** |
| xylose/glycine | **6** | **2.6** | **0.3** |
| xylose/glycine | **7** | **3.1** | **0.7** |
| xylose/L-alanine | **5** | **0.3** | **2.0** |
| xylose/L-alanine | **6** | **0.9** | **7.5** |
| xylose/L-alanine | **7** | **2.5** | **13.5** |

> ^a Equimolar amounts of precursors were used. Reaction conditions: **phosphate buffer
> (0.2 mol/L), 90 °C, 1 h**.
> ^b Data are means of at least two assays, each of them injected twice; **maximum SD ≤ 10 %**
> (in µg/mmol of sugar).

*(The pH 5 and pH 6 xylose/glycine rows are **identical to two significant figures** — 2.6/0.3 in
both. This is not an OCR duplication: I verified it in both the `-layout` and the `-raw` passes and
the two rows are physically distinct lines in the PDF. It is a real, flat result.)*

**Authors' reading + the mechanistic reason they give, verbatim (p. 2647):**
> "**The formation of both furanones was favored at pH 7 resulting in 3.1 µg of HDMF and 13.5 µg of
> EHMF (Table 4).** This is in good agreement with the more intense caramel-like overall flavor of
> these samples compared to those prepared at pH 5. In general, 3(2H)-furanones are formed from
> 1-deoxyosones which are reactive intermediates in the Maillard reaction generated by degradation
> of Amadori compounds via 2,3-enolization (Hodge et al., 1972; Ledl and Schleicher, 1990). **This
> reaction is favored under neutral and slightly alkaline conditions as shown by Beck et al. (1988),
> i.e. the ratio 1-deoxyglucosone to 3-deoxyglucosone (formed via 1,2-enolization) was 20:1 at pH 7
> and 8:5 at pH 4.5.**"

### 5.1 ★ THE AMINO-ACID-SPECIFIC pH RESPONSE — a falsifiable two-slope structure

**[DERIVED] pH 5 → pH 7 fold-change:**

| system | product | pH5 | pH6 | pH7 | **fold pH5→pH7** |
|---|---|---|---|---|---|
| xylose/**glycine** | HDMF | 2.6 | 2.6 | 3.1 | **1.19×** |
| xylose/**glycine** | EHMF | 0.3 | 0.3 | 0.7 | 2.33× |
| xylose/**L-alanine** | HDMF | 0.3 | 0.9 | 2.5 | **8.3×** |
| xylose/**L-alanine** | EHMF | 2.0 | 7.5 | 13.5 | **6.75×** |

**★ The glycine system is essentially pH-INSENSITIVE over pH 5−7 (1.19×, flat from 5 to 6) while
the alanine system is 7−8× pH-SENSITIVE over the same range.** Both products, both systems, one
lab, one protocol, one table. A model with a single pH factor on a `Furanone_Formation` family
**cannot** produce a 1.2× and an 8.3× response from the same term. **This is the single best
falsification target in the paper.**

### 5.2 [CITED, NOT MEASURED] the Beck 1988 enolisation ratio

| quantity | value | source as this paper cites it |
|---|---|---|
| **1-deoxyglucosone : 3-deoxyglucosone at pH 7** | **20 : 1** | Beck, Ledl & Severin 1988, *Carbohydr. Res.* **177**, 240−243 |
| **1-deoxyglucosone : 3-deoxyglucosone at pH 4.5** | **8 : 5** (= 1.6 : 1) | same |

**That is a 12.5× swing in the 1-DG/3-DG partition over 2.5 pH units, and it is a HEXOSE
measurement being used to rationalise a PENTOSE result.** It is exactly the 2,3- vs 1,2-enolisation
partition the repo's deoxyosone fork needs. **Secondary here — it must be attributed to Beck 1988
and NOT to Blank 1997, and Beck 1988 has not been read.** Flagged for `k3` §C and as a fetch
candidate (§11).

---

## 6. ★ TABLE 5 — verbatim (p. 2647). The precursor-ratio dose−response. **A reaction order.**

> **Table 5. Effect of Amino Acid Concentration on the Formation of HDMF and EHMF from Xylose in
> Maillard Model Systems^a**

| model system | molar ratio xylose/amino acid | HDMF^b | EHMF^b |
|---|---|---|---|
| xylose/glycine | **1:1** | **2.6** | **0.3** |
| xylose/glycine | **1:2** | **3.2** | **0.4** |
| xylose/glycine | **1:4** | **4.2** | **0.5** |
| xylose/L-alanine | **1:1** | **0.9** | **7.5** |
| xylose/L-alanine | **1:2** | **1.2** | **12.5** |
| xylose/L-alanine | **1:4** | **1.6** | **20.0** |

> ^a Reaction conditions: **phosphate buffer (0.2 M, pH 6.0), 90 °C, 1 h**.
> ^b Data are means of at least two assays, each of them injected twice; **maximum SD ≤ 10 %**
> (in µg/mmol sugar).

**Authors' reading, verbatim (p. 2647):**
> "The amounts of HDMF and EHMF can also be influenced by the ratio of the precursors (Table 5). In
> general, **increasing concentrations of the amino acids favor the formation the furanones, e.g. a
> 4-fold increase of alanine resulted in 2.7-fold more EHMF.** One possible explanation is that
> **elevated concentrations of amino acids give higher levels of aldehydes, the Strecker deamination
> products**, required for the formation of HDMF and EHMF from pentose sugars (Blank and Fay, 1996).
> Alternatively, **higher amounts of amino acids may accelerate the Maillard reaction, producing
> reactive intermediates such as 1-deoxyosones which may decompose to smaller fragments and so
> generate furanones by condensation reactions** (Blank et al., 1996a)."

### 6.1 ★★ [DERIVED] REACTION ORDER IN AMINO ACID — and the model currently violates it

Fitting a power law `yield ∝ [AA]^n` on the 1:1 → 1:4 endpoints, `n = ln(fold) / ln(4)`:

| system | product | 1:1 | 1:4 | fold | **apparent order n in amino acid** |
|---|---|---|---|---|---|
| xylose/glycine | **HDMF** | 2.6 | 4.2 | 1.615 | **0.35** |
| xylose/glycine | EHMF | 0.3 | 0.5 | 1.667 | 0.37 |
| xylose/L-alanine | HDMF | 0.9 | 1.6 | 1.778 | 0.42 |
| xylose/L-alanine | **EHMF** | 7.5 | 20.0 | 2.667 | **0.71** |

*(Log-log linearity check on the three points, xylose/Ala EHMF: 7.5 → 12.5 → 20.0 gives step ratios
1.667 and 1.600 for two successive doublings — **very close to a clean power law**, n = 0.737 and
0.678 on the two segments. The glycine HDMF series 2.6 → 3.2 → 4.2 gives 1.231 and 1.313, n = 0.30
and 0.39. **Both series are sub-linear and reasonably power-law.**)*

**★ EVERY ORDER IS SUB-LINEAR: n ≈ 0.35 − 0.71, never 1.0.**

**Why this matters to the repo, concretely.** `src/reaction_templates.py::_furanone_generation`
(lines 2218−2274) emits the lumped step

```
pentose + glycine  -> DMHF + NH3 + CO2 + 2 H2O     (reaction_family="Furanone_Formation")
pentose + alanine  -> HEMF + NH3 + CO2 + 2 H2O     (reaction_family="Furanone_Formation")
```

as a **bimolecular elementary step**, i.e. **first order (n = 1) in amino acid by construction.**
Table 5 says the true order is **0.35 for the glycine/HDMF pairing** — a factor-of-3 error in the
*shape* of the amino-acid dependence, which will show up as soon as any benchmark varies the
sugar:amine ratio away from 1:1. **This is a structural mis-specification the model can be scored
against today, without touching a single barrier.** (Sub-linear order is exactly what a
saturating/rate-limited upstream step gives — the 1-deoxyosone pool, not the amine, is limiting.)

---

## 7. THE SYNTHESIS OF THE INTERNAL STANDARDS — recorded because it is what makes the numbers real

Three steps from one common unlabelled starting material (Figure 1, p. 2643), verbatim conditions:
**a**, LIDA, ¹³CH₃¹³CHO, (t-BuOCO)₂O, THF; **b**, LIDA, CD₃CH₂CHO, (t-BuOCO)₂O, THF; **c**, KMnO₄,
acetone/H₂O/CH₃COOH; **d**, oxalic acid, H₂O.

| # | compound | yield | note |
|---|---|---|---|
| **1** | 3-butyn-2-yl *tert*-butyl carbonate | **92 %**, distilled 36−38 °C/0.3 mmHg, n_D²⁰ 1.4195 | Eren & Keinan 1988 method |
| **2** | [1,2-¹³C₂]-3-hexyn-2,5-diyl bis(*t*-Bu carbonate) | **quantitative** (7.18 g), 1:1 diastereomers | from [1,2-¹³C₂]ethanal, 21.7 mmol |
| **3** | [1,2-¹³C₂]-3,4-hexanedione-2,5-diyl bis(*t*-Bu carbonate) | **84 % from ¹³CH₃¹³CHO** (6.38 g) | unlabelled analogue **88 %** |
| **4ab** | **[¹³C₂]HDMF** | **43 %** (1.02 g) after sublimation 50 °C/0.3 mmHg; **mp 77−79.5 °C** (lit. Re 1973: 77−79 °C) | 1:1 tautomers 4a/4b; unlabelled HDMF crude **85 %** |
| **5** | [7,7,7-²H₃]-3-heptyn-2,5-diyl bis(*t*-Bu carbonate) | 8.26 g, pale yellow oil | from CD₃CH₂CHO (Olah & Arvanaghi 1981 formylation) |
| **6** | [7,7,7-²H₃]-3,4-heptanedione-2,5-diyl bis(*t*-Bu carbonate) | **87 % from propanal-d₃** (7.43 g) | unlabelled **85 %** |
| **7ab** | **[²H₃]EHMF** | **74 %** (2.02 g) after flash chromatography (SiO₂, 1:1 hexane/EtOAc) | **tautomer ratio 7a/7b = 1:3**; unlabelled EHMF crude **88 %** |

**Design rationale, verbatim (p. 2645):** *"In the case of HDMF, **¹³C-labeling was indispensable to
avoid exchange of hydrogen and oxygen atoms during sample preparation** (Sen et al., 1991). On the
contrary, EHMF could be deuterated in the terminal position of the ethyl group."*
→ **A ²H-labelled HDMF standard would have back-exchanged and given wrong numbers.** Worth knowing
if anyone ever proposes a d-labelled furaneol assay.

**Reference MS for the standards (from the synthesis section, p. 2644):**
- **HDMF** (unlabelled), EI *m/z* (%): **128 (100), 85 (23), 57 (67), 43 (82), 29 (22)**; CI/NH₃
  148 (100, [M+NH₄]⁺). IR 3300 br, 1680, 1600, 1295, 1185 cm⁻¹.
- **[¹³C₂]HDMF**, EI *m/z* (%): **130 (100), 87 (12), 85 (13), 59 (32), 57 (41), 45 (47), 43 (45),
  29 (18)**; CI/NH₃ 148 (100). ¹³C-NMR δ 13.6 (d, J = 48 Hz), 16.5 (d, J = 39), 80.5 (d, J = 39),
  174.5 (d, J = 48).
- **EHMF** (unlabelled), minor tautomer EI 142 (54), 85 (12), 57 (100); major tautomer 142 (87),
  127 (34), 102 (32), 71 (64), 55 (25), 43 (100); CI/NH₃ 160 (100).
- **[²H₃]EHMF**, minor tautomer 145 (71), 89 (10), 85 (26), 60 (100), 57 (48), 29 (15); major
  tautomer 145 (96), 127 (41), 102 (37), 73 (77), 55 (28), 43 (100); CI/NH₃ 163 (100).

**Tautomer ratios differ between the labelled standards: HDMF 4a:4b = 1:1, EHMF 7a:7b = 1:3.**
(Compare `blank1996_extraction.md` §7: Re, Maurer & Ohloff 1973 report the natural HEMF **2a:2b =
1:2** — this synthesis gives **1:3**, a *different* ratio, and the SIM traces both tautomers
separately, so nothing hangs on it.)

---

## 8. CROSS-CHECK AGAINST `blank1996_extraction.md` — corroborate or contradict, item by item

| 1996 constraint | 1996 proposed role | 1997 evidence | **VERDICT** |
|---|---|---|---|
| **#9** norfuraneol ≫ DMHF, HEMF | §10.2 "**the single most valuable hold-out in the paper**", ordinal ceiling | p. 2646, verbatim: *"However, **in all samples analyzed, 4-hydroxy-5-methyl-3(2H)-furanone was the main reaction product (data not shown)** directly formed from pentose sugars (Feather, 1981)."* | **★ CORROBORATED, and now across two papers and six sugar/amino-acid systems.** Still **not quantified** — "data not shown", no number, in either paper. **The ceiling holds; the ratio does not exist.** |
| **#4/#6** 70/30 Strecker/fragmentation split (glycine) | §10.1 "FIT — as a PRIOR on the C1-donor branch ratio" | **73/27** from the GABA control, §3 above, by a wholly non-isotopic method | **★★ CORROBORATED INDEPENDENTLY.** Best corroboration in the cluster. |
| **#7** HEMF requires alanine (on/off) | §10.2 HOLD-OUT, structural truth table | **EHMF = 0.3 / 0.7 / 1.3 µg/mmol in the three glycine systems** | **★ CONTRADICTED as on/off.** Demote to a **≈10−25× preference**. See §0b. |
| **#8** DMHF forms from pentose alone | §10.2 HOLD-OUT, "zero-amino-acid positive control" | **<0.01 µg/mmol**, both products, both media (Table 1 fn c; Table 3 rows 1−2) | **★ CONTRADICTED as a positive.** Re-specify as a **ceiling / negative control**. See §0b. |
| **#3** C5 + C1 → C6 / C5 + C2 → C7, pentose backbone intact | §10.2 HOLD-OUT, atom-mapping | p. 2646: *"the quantitative results support the hypothesis (Blank and Fay, 1996) that HDMF and EHMF are mainly formed by **Strecker-assisted chain elongation of pentoses** … **formaldehyde (active C1) liberated from glycine is required for the formation of HDMF as is acetaldehyde (active C2 from L-alanine) for EHMF**."* | **CORROBORATED** — but this paper adds **no new isotope evidence**; it is quantitative support for the 1996 topology, not an independent atom-mapping. |
| **#14** pH drift 6.0 → 5.0/5.3 | §10.2 HOLD-OUT, weak | **Not repeated.** No final pH is reported anywhere in this paper. | **NEITHER.** Still a single 1996 observation. |
| **#12/#13** "low mg/kg" rough estimate | §10.2 "NEITHER. Excluded." | **SUPERSEDED.** Use Tables 1−5 instead. | **Retire item #12/#13 in favour of this paper's measured µg/mmol.** |

### 8.1 ⚠️ A NAMING CHANGE BETWEEN THE COMPANION PAPERS — do not treat as two compounds

| compound | `blank1996` calls it | `blank1997` calls it |
|---|---|---|
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone, Furaneol | **HDMF** | **HDMF** |
| 2(or 5)-ethyl-4-hydroxy-5(or 2)-methyl-3(2H)-furanone, Homofuraneol | **HEMF** | **EHMF** ← *changed* |
| 4-hydroxy-5-methyl-3(2H)-furanone, norfuraneol | **"HMF (3)"** ← *the trap* | *no abbreviation; always spelled out in full* |

**Same lab, same first author, one year apart, three different abbreviation conventions.
`HEMF` (1996) and `EHMF` (1997) are the same molecule.**

---

## 9. THE NAMING TRAP — what "HMF"/"HDMF"/"DMHF" denotes IN THIS PAPER

**There is no abbreviation table in this paper.** All abbreviations are defined inline, on first use,
in the title/abstract and the first paragraph of the Introduction. I read every one of the 60
occurrences of the string `HMF` in the text layer and classified them:

| token | occurrences | denotes, in this paper |
|---|---|---|
| **HDMF** | majority | **4-Hydroxy-2,5-dimethyl-3(2H)-furanone** = **Furaneol®** (registered trademark of Firmenich S.A., Geneva — stated on p. 2642) = the repo's `DMHF`. |
| **EHMF** | majority | **2(or 5)-Ethyl-4-hydroxy-5(or 2)-methyl-3(2H)-furanone** = **Homofuraneol** = the repo's `HEMF`. |
| **`HMF` standing alone** | **ZERO** | **The string "HMF" never appears as a free-standing abbreviation in this paper.** Every one of the 60 hits is inside `HDMF` or `EHMF`. |
| 5-hydroxymethylfurfural | **ZERO** | The other "HMF" is **not mentioned anywhere** in this paper. |
| **norfuraneol** | 1 mention | Written out in full as **"4-hydroxy-5-methyl-3(2H)-furanone"** (p. 2646), never abbreviated. |

**→ THIS PAPER IS TRAP-FREE.** Unlike Blank & Fay 1996, where `HMF (3)` = norfuraneol
(`blank1996_extraction.md` §7), and unlike Apriyantono & Ames 1993, where `HMFone` = norfuraneol.
**`HDMF` here means furaneol and nothing else.**

---

## 10. KINETICS — **NONE, and the reason is worth recording**

No rate constant. No activation energy. No Arrhenius plot. No half-life. **One temperature (90 °C)
and one time (1 h) in every single experiment in the paper, all 39 quantified cells.** There is no
second time point and no second temperature anywhere, so not even a two-point rate can be built.

**What there IS, and it is not nothing:** five orthogonal *parameter* axes at that one T,t —
sugar identity, amino-acid identity, buffer identity, pH, and precursor ratio. **That is enough to
constrain a channel's SELECTIVITY and its RESPONSE SURFACE completely, and its RATE not at all.**
Any Arrhenius parameter the repo assigns to a furanone family remains, after this paper, the repo's
own assumption — exactly as `research_round3_channels.md` §B.3 concluded. **This paper does not
change that; it changes what the channel's *shape* has to look like.**

---

## 11. CITED, NOT MEASURED — carried separately

| quantity | value | true source, as this paper cites it | page |
|---|---|---|---|
| **1-deoxyglucosone : 3-deoxyglucosone, pH 7** | **20 : 1** | **Beck, Ledl & Severin 1988**, *Carbohydr. Res.* **177**, 240−243 | 2647 |
| **1-deoxyglucosone : 3-deoxyglucosone, pH 4.5** | **8 : 5** | same | 2647 |
| pH **4.0** is the stability optimum of HDMF in aqueous solution | qualitative | **Hirvi, Honkanen & Pyysalo 1980**, *LWT* **13**, 324−325 | 2644 |
| HDMF/EHMF tautomer cyclisation procedure; HDMF mp 77−79 °C | mp | **Re, Maurer & Ohloff 1973**, *Helv. Chim. Acta* **56**, 1882−1894 | 2645 |
| ¹³C₂-HDMF previously synthesised; ²H-HDMF back-exchanges | qualitative | **Sen, Schieberle & Grosch 1991**, *LWT* **24**, 364−369 | 2642, 2645 |
| ²H-EHMF previously prepared, "considerably reduced yield" | qualitative | **Preininger & Grosch 1994**, *LWT* **27**, 237−244 | 2642, 2645 |
| IDA calibration-curve procedure | method | **Guth & Grosch 1990**, *LWT* **23**, 513−522 | 2645 |
| 3(2H)-furanones form from 1-deoxyosones via 2,3-enolisation | qualitative | **Hodge, Mills & Fisher 1972**, *Cereal Sci. Today* **17**, 34−40; **Ledl & Schleicher 1990**, *Angew. Chem. Int. Ed.* **29**, 565−594 | 2642, 2647 |
| norfuraneol forms **directly** from pentose sugars | qualitative | **Feather 1981**, *Prog. Food Nutr. Sci.* **5**, 37−45 | 2646 |
| Boc as an alcohol protecting group | method | **Greene & Wuts 1991**, p. 387 | 2645 |

**⚠️ Beck 1988's 20:1 / 8:5 is a HEXOSE deoxyosone partition quoted inside a PENTOSE paper.** It is
the most transferable-looking number in this reference list and therefore the most dangerous.
**Do not attribute it to Blank 1997. Beck 1988 has not been read.**

---

## 12. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured in this paper · **[C]** cited from elsewhere · **[D]** derived by me here.
Conditions for every **[M]** row: **90 °C, 1 h, 1 M sugar : 1 M amino acid unless a ratio is
stated, 0.2 M buffer, duplicate assays each injected twice, maximum SD ≤ 10 %.**

| # | quantity | value | conditions | class | verdict |
|---|---|---|---|---|---|
| 1 | **HDMF, arabinose/glycine** | **5.1 µg/mmol sugar** (0.0040 mol %) | pH 6, phosphate | **[M]** | ★ USE. Highest HDMF in the paper. |
| 2 | **HDMF, xylose/glycine** | **2.6 µg/mmol** (0.0020 mol %) | pH 6, phosphate | **[M]** | ★ USE. **The paper's reference cell** — it recurs in Tables 1, 3, 4 and 5. |
| 3 | **HDMF, ribose/glycine** | **3.6 µg/mmol** (0.0028 mol %) | pH 6, phosphate | **[M]** | ★ USE. |
| 4 | **HDMF, arabinose / xylose / ribose + L-alanine** | **1.2 / 0.9 / 1.6 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE. The fragmentation-channel baseline in Blank's own reading. |
| 5 | **EHMF, arabinose / xylose / ribose + glycine** | **1.3 / 0.3 / 0.7 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE. **These three cells are what contradict the 1996 on/off switch.** |
| 6 | **EHMF, arabinose / xylose / ribose + L-alanine** | **6.8 / 7.5 / 10.0 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE. |
| 7 | **Both furanones, pentose alone (no amino acid)** | **< 0.01 µg/mmol** | pH 6, phosphate **and** water | **[M]**, upper bound | ★ USE **as a CEILING**. ≥260× below the glycine system. |
| 8 | **Sugar ordering, HDMF** | **arabinose > ribose > xylose**, 1.96× spread | pH 6, Gly | **[M]** | ★ USE. Closes the 1996 lumping gap. |
| 9 | **Sugar ordering, EHMF** | **ribose > xylose > arabinose**, 1.47× spread | pH 6, Ala | **[M]** | ★ USE. **Different ordering from #8.** |
| 10 | **Amino-acid selectivity swing (xylose)** | HDMF/EHMF **8.7 (Gly) → 0.12 (Ala) = 72×** | pH 6, phosphate | **[D]** from #2/#5/#6 | ★ USE. Pure ratio, scale-free. |
| 11 | **HDMF, xylose/4-aminobutyric acid** (Strecker-null) | **0.4 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE. Fragmentation-only baseline. |
| 12 | **EHMF, xylose/4-aminobutyric acid** | **0.1 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE. |
| 13 | **HDMF, xylose/GABA/glycine** | **1.5 µg/mmol** | pH 6, phosphate | **[M]** | ★ USE with #11. |
| 14 | **EHMF, xylose/GABA/L-alanine** | **3.2 µg/mmol**; HDMF in that run **0.7** | pH 6, phosphate | **[M]** | ★ USE with #12. |
| 15 | **★★ Strecker / fragmentation split, HDMF** | **73 % / 27 %** | pH 6, 90 °C, from #11+#13 | **[D]** | ★★ USE. **Independently reproduces Blank & Fay 1996's 70/30 by a non-isotopic method.** |
| 16 | **Strecker / fragmentation split, EHMF** | **97 % / 3 %** | from #12+#14 | **[D]** | ★ USE. EHMF is almost purely Strecker-derived. |
| 17 | Fragmentation share of HDMF, alternative estimates | **15 % (GABA/xylose-Gly), 24−44 % (alanine-baseline)** | — | **[D]** | Use as the **band 15−44 %**, centre ≈30 %. |
| 18 | **Phosphate catalysis, vs water** | **43× (HDMF, Gly); 45× (HDMF, Ala); 150× (EHMF, Ala)** | pH 6 nominal; ⚠️ water rows pH-uncontrolled | **[D]** from Table 3 | ★ USE with the caveat. Confounded with pH drift. |
| 19 | **Phosphate catalysis, vs 0.2 M malonate** | **4.3× / 3× (HDMF); 3× / 9.4× (EHMF)** | pH 6, both buffered | **[D]** from Table 3 | ★★ USE. **The clean, non-confounded phosphate-specificity number.** |
| 20 | **Unbuffered water yields** | **HDMF 0.06 (Gly), 0.02 (Ala); EHMF <0.01 (Gly), 0.05 (Ala)** | pH uncontrolled | **[M]** | ★ USE as a near-floor. |
| 21 | **★ pH ladder, xylose/glycine** | HDMF **2.6 / 2.6 / 3.1**; EHMF **0.3 / 0.3 / 0.7** at pH **5 / 6 / 7** | 0.2 M phosphate | **[M]** | ★★ USE. |
| 22 | **★ pH ladder, xylose/L-alanine** | HDMF **0.3 / 0.9 / 2.5**; EHMF **2.0 / 7.5 / 13.5** at pH **5 / 6 / 7** | 0.2 M phosphate | **[M]** | ★★ USE. |
| 23 | **★★ Amino-acid-specific pH sensitivity** | pH5→7: **glycine 1.19× vs alanine 8.3× (HDMF)**; 2.33× vs 6.75× (EHMF) | — | **[D]** from #21/#22 | ★★ USE. **The best falsification target in the paper.** |
| 24 | **Dose−response, xylose/glycine 1:1/1:2/1:4** | HDMF **2.6 / 3.2 / 4.2**; EHMF **0.3 / 0.4 / 0.5** | pH 6, phosphate | **[M]** | ★ USE. |
| 25 | **Dose−response, xylose/L-alanine 1:1/1:2/1:4** | HDMF **0.9 / 1.2 / 1.6**; EHMF **7.5 / 12.5 / 20.0** | pH 6, phosphate | **[M]** | ★ USE. |
| 26 | **★★ Apparent reaction order in amino acid** | **n = 0.35 (HDMF/Gly), 0.37 (EHMF/Gly), 0.42 (HDMF/Ala), 0.71 (EHMF/Ala)** — **all sub-linear** | 90 °C, pH 6 | **[D]** from #24/#25 | ★★ USE. **The model's bimolecular template is n = 1 by construction and therefore already violates this.** |
| 27 | **Uncertainty** | **maximum SD ≤ 10 %**, n ≥ 2 assays × 2 injections | all tables | **[M]** | ★ USE as the σ for every cell. |
| 28 | Calibration linearity | **3−50 µg/mL, r² = 0.999** | SIDA | **[M]** | Method quality; supports #27. |
| 29 | Internal standard load | **9.6−48.2 µg [¹³C₂]HDMF; 10.0−50.2 µg [²H₃]EHMF** | added pre-workup | **[M]** | Method record. |
| 30 | Standard synthesis yields | HDMF chain **92 / quant. / 84 / 43 %**; EHMF chain **92 / quant. / 87 / 74 %** | — | **[M]** | Method record; not model-relevant. |
| 31 | Tautomer ratios of the standards | **[¹³C₂]HDMF 4a:4b = 1:1**; **[²H₃]EHMF 7a:7b = 1:3** | synthesis | **[M]** | Note: Re 1973 gives natural HEMF 2a:2b = **1:2**. |
| 32 | norfuraneol is the **main product** in every system | ordinal | all runs | **[M]**, "data not shown" | ★ USE as the CEILING. **No number exists.** |
| 33 | 1-DG : 3-DG = **20:1 (pH 7)**, **8:5 (pH 4.5)** | ratio | **hexose** | **[C]** Beck 1988 | **Secondary. Attribute to Beck, not Blank. Not read.** |
| 34 | HDMF stability optimum pH 4.0 | qualitative | aqueous | **[C]** Hirvi 1980 | Secondary. |
| 35 | methylglyoxal, acetol and dihydroxyacetone **do generate HDMF** | qualitative, **"data not shown"** | — | **[M]**, unpublished | ★ Note only. **Corroborates the Wang & Ho 2008 MGO edge from the same lab.** |
| 36 | **rate constant, Ea, half-life, time course, temperature series** | **NONE** | — | — | **ABSENT. §10.** |
| 37 | **final pH of any buffered run** | **NOT REPORTED** | — | — | **ABSENT. §1.1.** |
| 38 | **phosphate concentration series** | **NONE** — 0.2 M only | — | — | **ABSENT.** |

---

## 13. VERDICT — does this close the DMHF channel?

**It closes the SELECTIVITY and the RESPONSE SURFACE. It does not close the RATE.**

| requirement | met? | detail |
|---|---|---|
| **A quantitative anchor with a defensible basis** | **★★ YES.** | 39 SIDA cells in µg/mmol sugar, internal standard added pre-workup, SD ≤ 10 %, r² = 0.999. **The first ingestible absolute DMHF number in the corpus.** |
| **pH dependence** | **★★ YES**, and amino-acid-resolved | Table 4, three points × two amino acids × two products. |
| **Buffer / catalysis dependence** | **★ YES** | Table 3, phosphate vs malonate vs water. |
| **Precursor-ratio dependence** | **★★ YES**, and it yields a reaction order | Table 5, three ratios × two systems × two products. |
| **Sugar-identity dependence** | **★ YES** | Table 1, three pentoses, resolved for the first time. |
| **Branch fraction (Strecker vs fragmentation)** | **★★ YES**, independently of the isotope method | Table 2, GABA control. |
| **Kinetics (k, Ea, time course)** | **❌ NO** | One T, one t. §10. |
| **Absolute norfuraneol : DMHF ratio** | **❌ NO** | "data not shown" in both companion papers. |
| **Phosphate dose−response** | **❌ NO** | One concentration. |
| **Final pH / drift** | **❌ NO** | Not reported. §1.1. |

**One-line honest answer:** *this is the numbers paper `blank1996_extraction.md` §0a said had to be
fetched; it is now in the repo and read; it makes a pentose furanone channel CALIBRATABLE in
magnitude and in five parameter dimensions, it leaves the channel's RATE entirely to the repo's own
assumption, and it invalidates two of the four hold-outs the 1996 dossier proposed.*

---

## 14. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION.**

Nothing below has been written to any registry, benchmark file, config or declaration.

### 14.1 ⚠️ THE INDEPENDENCE PROBLEM — read before splitting anything

**The five tables are NOT 39 independent measurements.** Two cells are re-reported four times each:

- `xylose/glycine, pH 6, phosphate, 1:1` → **HDMF 2.6 / EHMF 0.3** appears in **Table 1, Table 3,
  Table 4 (pH 6 row) and Table 5 (1:1 row)**.
- `xylose/L-alanine, pH 6, phosphate, 1:1` → **HDMF 0.9 / EHMF 7.5** appears in the **same four
  tables**.

These are the same two experiments quoted four times, not sixteen results. **Any fit/hold-out split
that puts one instance of the reference cell in the fit set and another in the hold-out set has
leaked.** The safe unit of splitting is the **table**, or better, the **experimental axis**.

### 14.2 Proposed **FIT** candidates

| candidate | proposed role | why |
|---|---|---|
| **Table 1, the six sugar × amino-acid cells** (items #1−#6) | **FIT — magnitude.** Six absolute levels at one condition. | The only absolute DMHF/HEMF anchors the corpus has. If a magnitude is to be calibrated at all, it is calibrated here. Give every cell σ = 10 % (item #27), not a guessed σ. |
| **Item #15/#16, the 73/27 and 97/3 Strecker/fragmentation splits** | **FIT — as a branch PRIOR**, jointly with `blank1996_extraction.md` items #4/#5/#6. | Two independent methods agree to 3 pp. Ratio-form, so immune to absolute-scale error. **This is the safest fit target in the whole cluster.** |
| **Item #19, phosphate/malonate = 3−9.4×** | **FIT — only if a buffer/catalysis term is ever added.** | Clean, both arms buffered at the same pH. **Do NOT fit against item #18 (the water comparison): its pH was uncontrolled and the effect is confounded.** |
| **NO barrier or Ea may be fitted against this paper.** | — | **One temperature. §10.** A barrier fitted to a single-temperature dataset is a rate-constant fit wearing an Arrhenius costume; it will absorb every unmodelled factor in the system. This is the `thiol_addition_pentodiulose` failure mode (`barrier_constants.py:307`) in a new dress. **Record as a prohibited derivation.** |

### 14.3 Proposed **HOLD-OUT** candidates

| candidate | proposed role | why |
|---|---|---|
| **★★ Item #23 — the amino-acid-specific pH sensitivity (glycine 1.19× vs alanine 8.3× over pH 5→7)** | **HOLD-OUT, first choice.** | It is a **contrast of two slopes**, so it survives any error in the absolute scale, and no single pH factor on a shared reaction family can produce both. **Falsifiable in one run against a model that is already built.** |
| **★★ Item #26 — apparent reaction order in amino acid, n = 0.35−0.71** | **HOLD-OUT, second choice — and it is expected to FAIL today.** | `_furanone_generation` emits a bimolecular step, order 1 in amino acid by construction (§6.1). **Declaring this hold-out records a KNOWN, DIAGNOSED structural defect rather than discovering one later**, which is the honest way to carry it. |
| **Item #7 — pentose alone < 0.01 µg/mmol, both products** | **HOLD-OUT, structural CEILING.** ⚠️ **This REPLACES `blank1996_extraction.md` §10.2 item #8, which proposed the opposite sign.** | Cleanest negative control in the cluster: two products, two media, one bound. A model that emits appreciable DMHF from a pentose with no amine fails. |
| **Item #10 — the 72× HDMF/EHMF selectivity swing on xylose** | **HOLD-OUT, ratio.** | Scale-free; tests the amino-acid coupling directly. |
| **Items #8/#9 — the two DIFFERENT sugar orderings** | **HOLD-OUT, ordinal, hard.** | arabinose > ribose > xylose for HDMF **and** ribose > xylose > arabinose for EHMF, simultaneously. A single "pentose reactivity" scalar cannot satisfy both. |
| **Item #32 — norfuraneol is the main product in every system** | **HOLD-OUT, ordinal ceiling.** Carried jointly with `blank1996_extraction.md` #9. | Now supported by two papers. **Still unquantified — the hold-out can only ever be `norfuraneol > DMHF`, never a ratio.** |
| **Item #33 (Beck's 20:1 / 8:5)** | **NEITHER.** | Cited, hexose, unread. §11. |

### 14.4 ★ THE PAIRING THIS PAPER MAKES POSSIBLE — flagged for the synthesis

`blank1997` Table 4 and `wang2008` Figure 1 are **two pH ladders on two DIFFERENT DMHF formation
edges, with OPPOSITE SIGNS in the glycine system** (this paper: xylose/Gly rises 1.19× from pH 5 to
7; Wang & Ho: MG/Gly **falls 9×** from pH 3 to 8). **They must not be averaged into one pH term.**
Full treatment in `k5b_dmhf_synthesis.md` §4.

### 14.5 Prerequisites before any role is declared

1. **Correct `blank1996_extraction.md` §10.2 items #7 and #8** — or, better, leave that dossier as
   the record of what the 1996 paper said and let *this* dossier's §0b/§8 be the operative version.
   **What must not happen is item #7 being declared as a structural hold-out; it would fail a
   correct model.**
2. **Resolve defect D2 of `blank1996_extraction.md` §0c, which is now partly OBSOLETE.** That
   dossier states *"There is no pentose → DMHF step and no pentose → HEMF step anywhere in the
   network."* **That is no longer true**: `src/reaction_templates.py::_furanone_generation`
   (lines 2218−2274) emits both `pentose + glycine -> DMHF + NH₃ + CO₂ + 2 H₂O` and
   `pentose + alanine -> HEMF + NH₃ + CO₂ + 2 H₂O` as `Furanone_Formation`, gated at
   T ≥ 90 °C, `source_quality="literature"`, `barrier_uncertainty_kcal=6.0`. Both balance exactly.
   **The channel exists and is scorable against Table 1 today.** The 1996 dossier appears to have
   inspected `_furanone_and_mft_route` only.
3. **Do not use "0.2 M phosphate, pH 6" as a held condition** without carrying §1.1's caveat.

---

## 15. DECLARED GAPS — verbatim, for `k3` §C

> **"One temperature and one time.** 90 °C / 1 h in all 39 quantified cells. **No rate constant, no
> activation energy, no time course and no temperature series exist in this paper.** Five parameter
> axes are varied (sugar, amino acid, buffer, pH, precursor ratio) and none of them is time or
> temperature. Any Arrhenius parameter the repo assigns to a furanone family remains the repo's own
> assumption after reading this paper."

> **"The pH labels are buffer SET-POINTS and the paper never states whether they were held.** No
> final pH is reported for any run. The identical 0.2 mol/L phosphate buffer drifted 6.0 → 5.0/5.3
> over the same hour in the companion paper (Blank & Fay 1996). Only the *water* rows of Table 3
> carry an explicit 'pH was not controlled' footnote; the buffered rows carry no statement either
> way."

> **"The norfuraneol : HDMF ratio is 'data not shown' in BOTH companion papers.** Blank 1997 p. 2646
> states that 4-hydroxy-5-methyl-3(2H)-furanone 'was the main reaction product' in all samples
> analysed and gives no number. The ordinal ceiling norfuraneol ≫ DMHF is therefore supported twice
> over and quantified zero times."

> **"Blank & Fay 1996's 'HEMF requires alanine' is NOT an on/off switch.** The 1997 isotope-dilution
> assay measures EHMF at 0.3 (xylose), 0.7 (ribose) and 1.3 (arabinose) µg/mmol in the glycine
> systems where the 1996 GC-O table records '−'. The 1996 '−' is a non-detection at a sniffing port,
> not a zero. **Any hold-out asserting zero HEMF from pentose + glycine would fail a correct model.**"

> **"Blank & Fay 1996's 'HDMF forms from pentose alone' is a CEILING, not a positive control.** The
> 1997 control without amino acid gives < 0.01 µg/mmol of both furanones, in phosphate and in water
> — at least 260× below the xylose/glycine value, and reported only as an upper bound."

> **"The 20:1 / 8:5 1-deoxyglucosone : 3-deoxyglucosone ratio is Beck, Ledl & Severin 1988,
> Carbohydr. Res. 177:240−243, reached through this paper, and it is a HEXOSE measurement quoted
> inside a PENTOSE paper.** It has not been read. It must not be attributed to Blank 1997."

> **"No phosphate concentration series exists.** 0.2 mol/L is the only level, in every experiment in
> both companion papers. The 43−150× phosphate-vs-water effect is measured against an arm whose pH
> was uncontrolled; only the 3−9.4× phosphate-vs-malonate comparison is unconfounded."

> **"The 39 cells are not 39 independent measurements.** The `xylose/glycine, pH 6, phosphate, 1:1`
> and `xylose/L-alanine, pH 6, phosphate, 1:1` conditions each appear in FOUR of the five tables.
> Fit/hold-out splits must be made by axis or by table, never cell-by-cell."

---

## 16. WHAT TO FETCH NEXT — ranked, from this paper's reference list

| # | paper | why | confidence |
|---|---|---|---|
| 1 | **Beck, J.; Ledl, F.; Severin, T. 1988**, *"Formation of 1-deoxy-D-erythro-2,3-hexodiulose from Amadori compounds"*, *Carbohydr. Res.* **177**, 240−243 | The **pH-dependent 1-DG : 3-DG partition (20:1 at pH 7, 8:5 at pH 4.5)** — the enolisation fork the repo's deoxyosone node needs, and the mechanistic explanation Blank uses for his own pH ladder. **Currently a secondary, unread citation doing load-bearing work.** | Citation transcribed verbatim from the reference list. High. |
| 2 | **Blank, I.; Devaud, S.; Fay, L. B. 1996a**, *"New aspects on the formation of 3(2H)-furanones through the Maillard reaction"*, in *Flavour Science: Recent Developments* (Taylor & Mottram, eds.), RSC, Cambridge, pp. 188−193 | Cited **three times** here as the source of the **sugar fragmentation−condensation** route, and it is the companion the 1997 paper leans on for everything the 1996 JAFC paper does not cover. Book chapter, no DOI. | High on identity, unknown availability. |
| 3 | **Sen, A.; Schieberle, P.; Grosch, W. 1991**, *LWT* **24**, 364−369 | The original HDMF stable-isotope-dilution assay. Needed if these numbers ever need method scrutiny. Already ranked #5 in `blank1996_extraction.md` §12. | High. |
| 4 | **Hirvi, T.; Honkanen, E.; Pyysalo, T. 1980**, *LWT* **13**, 324−325 | *"Stability of 2,5-dimethyl-4-hydroxy-3(2H)-furanone … in aqueous buffer solutions."* **This is a DMHF DEGRADATION paper with a pH axis** — i.e. a candidate source for the DMHF **sink** the model does not have. Pairs directly with Shu & Ho 1988. | High on identity; **★ raised in priority by this wave**, see the synthesis. |
| 5 | **Preininger, M.; Grosch, W. 1994**, *LWT* **27**, 237−244 | The ²H-EHMF standard and the Emmentaler OAV work. Method only. | High. |

---

## 17. SOURCES USED IN THIS EXTRACTION BEYOND THE PDF

**None.** Every statement above is from `data/articles/blank1997.pdf` itself, from
`data/lit/extraction_dossiers/blank1996_extraction.md`, from
`data/lit/extraction_dossiers/research_round3_channels.md` §B, or from reading
`src/reaction_templates.py` and `src/barrier_constants.py` in the working tree (read-only; nothing
in `src/` was modified). No web lookup was performed and none was needed — the DOI is printed in the
document.

*End of dossier. Nothing outside this file was created or modified.*
