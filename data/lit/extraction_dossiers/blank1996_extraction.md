# Blank & Fay 1996 (10.1021/jf950439o) — single-paper extraction, 2026-08-29

> Blank, I.\*; Fay, L. B.
> "**Formation of 4-Hydroxy-2,5-dimethyl-3(2H)-furanone and 4-Hydroxy-2(or 5)-ethyl-5(or 2)-methyl-3(2H)-furanone through Maillard Reaction Based on Pentose Sugars**"
> *J. Agric. Food Chem.* **1996**, 44 (2), 531−536. Nestlé Research Center, Nestec Ltd.,
> Vers-chez-les-Blanc, CH-1000 Lausanne 26, Switzerland.
> Received for review July 14, 1995. Accepted November 22, 1995. Abstract published in
> Advance ACS Abstracts, January 15, 1996. Article ID `JF950439O`.

**Source PDF:** `data/articles/blank1996.pdf` (6 pp.).
**Read method:** full text layer, `pdftotext -layout`. Tables 1 and 2 came out of the text layer
intact and are transcribed verbatim below. **Figures 1, 2, 4 (GC-MS/MS spectra) and Figures 3
and 5 (the mechanism schemes) have NO text layer** — they are line-art/spectra only. Every
mechanistic step below is therefore taken from the **running text and the figure legends**,
which spell the scheme out completely; nothing was digitised off a figure and nothing was
inferred from a drawing I could not read.

---

## 0. VERDICT UP FRONT — read this before anything else

### 0a. ⚠️ THE BRIEF'S PREMISE IS HALF WRONG. THIS IS THE RIGHT PAPER, BUT NOT THE RIGHT *KIND* OF PAPER.

The DOI, title, authors, journal, year and volume/pages **all match the brief exactly**. The file
is not mis-filed. But the brief describes it as *"isotope-**dilution quantification** of
DMHF/HEMF from pentoses at ~90 °C with **pH and phosphate varied**"*, and **none of that is in
this paper**:

| The brief expected | What the paper actually contains |
|---|---|
| Isotope-**dilution** quantification (SIDA) | **NO.** Isotope **labelling** for mechanism only. No internal standard, no calibration, no response factor, no absolute concentration. |
| Yield tables in mol % | **NO YIELD TABLE OF ANY KIND EXISTS.** The paper has exactly two tables: GC-O odour intensities (`+`/`++`/`+++`) and MS fragmentation *m/z* lists. |
| pH varied | **NO.** One buffer, pH 6.0 initial, uncontrolled drift to 5.0 / 5.3. |
| Phosphate varied | **NO.** One phosphate concentration, 0.2 mol/L, in every experiment. |
| Time courses / kinetics | **NO.** One time point, 1 h. One temperature, 90 °C. |
| ~90 °C | **YES** — 90 °C, correct. |

The authors say so themselves, in the last sentence of the Discussion (p. 535):

> "**As a continuation of this study, the amounts of HDMF and HEMF are being quantified and
> parameters affecting this reaction are currently being studied.**"

**The paper the brief actually described is the COMPANION, which is NOT in `data/articles/`:**
Blank, I.; Fay, L. B.; Lakner, F. J.; Schlosser, M., *"Determination of 4-Hydroxy-2,5-dimethyl-3(2H)-furanone
and 2(or 5)-Ethyl-4-hydroxy-5(or 2)-methyl-3(2H)-furanone in Pentose Sugar-Based Maillard Model
Systems by **Isotope Dilution Assays**"*, **J. Agric. Food Chem. 1997, 45, 2642−2648, DOI
`10.1021/jf960997i`** (~85 % confidence on the DOI/pagination — taken from a web bibliographic
lookup during this extraction, **NOT CrossRef-verified and NOT read**; the title, author list,
journal and year are corroborated by three independent listings). That is the paper carrying the
IDA quantification and the parameter study. **It should be fetched.** See §12.

### 0b. WHAT THIS PAPER DOES DELIVER, AND IT IS NOT NOTHING

1. **A complete, isotope-evidenced formation route for DMHF and HEMF from PENTOSES**, stated by
   the authors as a numbered scheme with every step named (Figure 3 legend, p. 534). The
   precursor question the brief asked — *acetylformoin, or 1-deoxyosone + Strecker aldehyde?* —
   is answered **BOTH, in series**: 1-deoxyosone **+** Strecker aldehyde **→ acetylformoin →**
   DMHF. See §5. This is a **route the model does not have**.
2. **Three isotopomer-distribution measurements** (70/30; 4/12/84; 15/47/38 %) which are genuine
   *numbers* and constitute a **branch-fraction constraint** on the C1-unit source
   (Strecker vs. sugar fragmentation), amino-acid-resolved. See §4.3. These are the only
   quantitative results in the paper.
3. **An ordinal magnitude constraint**: norfuraneol ("HMF", **3**) is the **major** volatile of the
   pentose/amino-acid system, while DMHF and HEMF are **minor peaks** formed "by a side
   reaction". See §3.3 / §7.
4. **A confirmed amino-acid on/off switch**: HEMF requires alanine (absent with glycine, absent
   from pentose alone); DMHF forms with glycine, with alanine, **and from pentose alone**. §3.2.

### 0c. THE THREE REPO CITATION DEFECTS THIS EXTRACTION FOUND — flagged, not fixed

I was told to touch nothing but this file. These are **for the orchestrator**. All three are
places where the repo cites `10.1021/jf950439o` for something the paper does not say.

| # | Where | The repo's claim | What the paper says | Severity |
|---|---|---|---|---|
| **D1** | `src/barrier_constants.py:325` (`furanone_amino_acid_reduction`); `src/reaction_templates.py` (Wave P item 6 comment); `AUDIT.md:571`; `tasks/audit_remediation.md:2842`; `tests/unit/test_chemistry_soundness.py:287` | *"the accepted mechanism names the reductant and it is the amino acid (Blank & Fay 1996 …)"* | **Blank & Fay 1996 do NOT name the amino acid as the reductant of the acetylformoin step.** They name the alternatives explicitly (p. 534): *"The reduction may occur **either by a dismutation or by a reaction with further enoloxo compounds**, as recently reported by Schieberle (1992)."* The abstract phrase the repo quotes — "reduction of the resulting acetylformoin-type intermediates" — names **no reductant at all**. In this paper the amino acid's demonstrated role is the **Strecker-aldehyde donor**, i.e. a *carbon* donor, not a hydride donor. | **HIGH** — the citation is doing load-bearing work: it is the stated justification for deleting the `[HH]` pool gate. *(Partial defence, stated for fairness: the Strecker step itself does reduce the dicarbonyl, and the co-citation **Kerler et al. 2010** genuinely does carry a Strecker-active-amino-acid claim — I have not read it. So the CLAIM may survive on Kerler alone; the **attribution to Blank & Fay 1996 does not.**)* |
| **D2** | `tasks/audit_remediation.md:1713` | *"DMHF/HEMF from pentose + amino acid is CONFIRMED with the exact **C5+C1 / C5+C2 stoichiometry the model uses**"* | The C5+C1 / C5+C2 stoichiometry **is** what the paper confirms (§5). But **the model does not use it.** `_furanone_and_mft_route` gives the pentose 1-deoxyosone exactly two fates — `Furanone_Cyclisation` → norfuraneol, and `Deoxyosone_Reduction` → 1,4-dideoxypentodiulose. **There is no pentose → DMHF step and no pentose → HEMF step anywhere in the network**, and `test_chemistry_soundness.py` pins that absence (`assert [...] == ["Furanone_Cyclisation"]`). The repo's only DMHF step is the **hexose** one. | **HIGH** — a "CONFIRMED" verdict recorded against a route that is not implemented. |
| **D3** | `src/barrier_constants.py:325` step stoichiometry: `hexose 1-deoxyosone + amino acid -> DMHF + Strecker aldehyde + CO2 + NH3 + H2O` | The Strecker aldehyde is a **co-product** that leaves. | In Blank & Fay's mechanism the Strecker aldehyde is a **reactant that is consumed into the furanone ring** (aldol condensation, step *g*). The repo's step and the cited paper's mechanism are **opposite in the aldehyde's fate**. | **MEDIUM** — for a *hexose* the repo's intact-C6 topology is independently right (Wang & Ho CAMOLA), so the *hexose* step is not wrong; but Blank & Fay 1996 is the wrong citation to hang on it, and the *pentose* case genuinely does incorporate the aldehyde. |

---

## 1. SYSTEM DEFINITION — verbatim (Experimental Procedures, p. 532)

### 1.1 Materials

> "D-Xylose, D-ribose, D-arabinose, glycine, and L-alanine of highest purity (>99%) were obtained
> from Fluka (Buchs, Switzerland). [2-¹³C]Glycine and [1-¹³C]-D-xylose were from Cambridge Isotope
> Laboratories (Andover, MA), and [3-¹³C]-L-alanine was from Tracer Technologies (Somerville, MA).
> **The isotopic content of the labeled compounds was 99%.**"

Reference compounds: HDMF (Furaneol) from Aldrich; HEMF (homofuraneol) from Givaudan.
Buffer salt: disodium hydrogen phosphate dihydrate (Merck).

### 1.2 The ONE reaction condition — there is no other in the paper

> "In a 15 mL Pyrex tube, **5 mmol of pentose** (xylose, ribose, or arabinose) and **5 mmol of amino
> acid** (glycine or alanine) were dissolved in **5 mL of phosphate buffer (0.2 mol/L Na₂HPO₄,
> pH 6.0)**. The tube was sealed with a screw cap and heated at **90 °C for 1 h** in an oil bath
> while stirring with a magnetic stirrer. The reaction was stopped by rapid cooling with tap water.
> **During the reaction, the pH dropped to 5.0 (xylose/glycine) and 5.3 (xylose/alanine).**"

Derived, arithmetic only: **1.0 mol/L pentose and 1.0 mol/L amino acid, 1:1**, sealed, unstirred
headspace, no pH control. **⚠️ The pH label "6.0" is INITIAL, not held** — the same defect class
already recorded for Zhou 2023 in `k3_final_parameter_inventory.md` §C.11. The *operating* pH of
these runs is a **6.0 → 5.0 (Gly) / 6.0 → 5.3 (Ala) drift**, and any use of "pH 6" for this system
must carry that.

### 1.3 Workup (matters, because it is a hot-acid step downstream of a labile furanone)

> "Then, 40 mL of water was added to the dark brown reaction mixture, which was then saturated
> with 16 g of NaCl. **The pH was adjusted to 4.0 (aqueous HCl, 2 mol/L)**, and the neutral
> compounds were continuously extracted with 50 mL of diethyl ether overnight using a rotation
> perforator … The organic phase was separated, **dried over sodium sulfate at 4 °C**, and
> concentrated to 0.5 mL using a Vigreux column … and a microdistillation device according to the
> procedure of Bemelmans (1979)."

> "Experiments with the labeled compounds were performed **in the same way** as described above."

**Note for any future absolute-yield reader:** an overnight ether perforation of a 45 mL aqueous
phase at pH 4 with no internal standard has **unknown and uncorrected recovery**. This is exactly
why §3.3's "low mg/kg" statement cannot be an anchor, and exactly why the 1997 companion had to be
done as an **isotope-dilution** assay.

### 1.4 Analytics

- **GC-O**: Carlo Erba Mega 2, on-column, effluent split 1:1 FID / sniffing port. Columns OV-1701
  (medium polarity) and FFAP (high polarity), both 30 m × 0.32 mm, 0.25 µm film. Program 50 °C
  (2 min) → 6 °C/min → 180 °C → 10 °C/min → 240 °C (10 min). LRI by van den Dool & Kratz (1963).
- **GC-MS**: HP-5971 MSD + HP-5890, Carbowax 30 m × 0.32 mm × 0.25 µm. Program 20 °C (0.5 min) →
  30 °C/min → 100 °C → 4 °C/min → 145 °C (10 min) → 70 °C/min → 220 °C (2.5 min). Splitless 250 °C,
  interface 220 °C, EI 70 eV, source ~180 °C.
- **GC-MS/MS**: Finnigan TSQ-700 + HP-5890, FFAP column. Splitless 280 °C; 60 °C (1 min) →
  10 °C/min → 200 → 30 °C/min → 240 °C (2 min). EI 70 eV, source 150 °C. CID of the molecular ion,
  daughter spectra 20–200 Da, **collision energy 10 eV (laboratory frame)**, argon at **1.1 mTorr**.

Rationale given for MS/MS (p. 533): *"in such complex matrices, some interferences can contaminate
the mass spectra, particularly at the low mass ranges (20−100 Da) where characteristic fragments of
HDMF and HEMF are observed."*

---

## 2. TABLE 1 — verbatim (p. 532). **This is the paper's only "results" table, and it is ORDINAL.**

> **Table 1. Sensory Contribution of HDMF and HEMF Estimated by GC-Olfactometry^a in Maillard
> Systems Based on Pentoses and either Glycine (Gly) or L-Alanine (Ala) and in Model Reactions
> Containing Only Pentose^b**

| compd | LRI OV-1701 | LRI FFAP | pentose/Gly^c | pentose/Ala^c | pentose^c |
|---|---|---|---|---|---|
| HDMF **1** | 1235 | 2040 | **++** | **+** | **+** |
| HEMF **2a** | 1310 | 2090 | **−** | **+++** | **−** |
| HEMF **2b** | 1325 | 2190 | **−** | **+**^d | **−** |

> ^a GC-O data are presented in terms of **odor intensities perceived at the sniffing port**
> (+: weak; +++: intense).
> ^b The model mixtures were heated in a phosphate buffer (0.2 mol/L, pH 6) at 90 °C for 1 h.
> ^c The pentoses used were **ribose, xylose, and arabinose**.
> ^d The tautomer 2b was sensorily detectable only in highly concentrated samples.

**Three things to notice, because they limit this table hard:**
1. It is **odour intensity at a sniffing port**, not concentration. `+++` vs `+` is **not** a
   concentration ratio — it convolves the concentration with the (very different) odour thresholds
   of 2a and 2b, which the authors say differ ("the odor activity of 2a was significantly
   higher than that of 2b (Blank, unpublished results, 1995)").
2. **The three pentoses are LUMPED into one column.** Ribose, xylose and arabinose are not
   resolved anywhere in this paper. There is **no sugar-reactivity ordering** to be had here, and
   any attempt to read one into it would be inventing data.
3. The `+`/`−` grid **is** a clean set of on/off structural facts. That part is usable.

---

## 3. THE QUANTITATIVE CONTENT, SUCH AS IT IS

### 3.1 Identification claim (p. 533)

> "To the best of our knowledge, **this is the first time that the furanones 1 and 2 have been
> reported in Maillard reaction systems based on pentoses.**"

Basis: LRI match on two columns of different polarity + GC-MS + interference-free GC-MS/MS
daughter spectra against authentic reference compounds. **Identification is solid.**

### 3.2 The amino-acid on/off structure (p. 533, Table 1)

> "**HDMF was mainly responsible for the caramel-like note** of the reaction mixture containing
> pentose/glycine or only a pentose (Table 1). In the system pentose/alanine, **both HDMF (1) and
> HEMF (2) contributed** to the caramel-like, sweet character. **Both furanones were represented
> by minor peaks in the gas chromatograms.**"

| system | HDMF | HEMF |
|---|---|---|
| pentose + glycine | formed (strongest GC-O) | **NOT DETECTED** |
| pentose + alanine | formed | formed (2a intense; 2b only in concentrates) |
| pentose alone (no amino acid) | **formed** | **NOT DETECTED** |

Abstract, verbatim: *"HEMF was detected in the system pentose/alanine. HDMF was formed in both
pentose/glycine and pentose/alanine systems **as well as directly from pentoses**."*

### 3.3 ⚠️ THE ONLY CONCENTRATION STATEMENT IN THE PAPER (p. 533) — and it is not an anchor

> "**Rough estimates of the amounts of HDMF and HEMF formed from pentoses were in the low
> milligrams per kilogram (parts per million) range.**"

That is the entire quantitative yield content. **It has no method, no basis of reference (per kg of
*what*?), no per-system breakdown, no uncertainty, and the authors label it a "rough estimate"
themselves.** It is a one-clause aside inside a paragraph whose next sentence defers quantification
to a future paper.

**[DERIVED — DO NOT INGEST]** For scale only, so the orchestrator knows what order of magnitude the
1997 companion is likely to report, here is my arithmetic with the basis ambiguity carried
explicitly. Reaction tube contents ≈ 0.751 g xylose + 0.375 g glycine + ~5.0 g buffer ≈ **6.13 g**;
after workup dilution ≈ **61 g** aqueous. "Low mg/kg" read as 1–5 mg/kg, MW(HDMF) = 128.13,
5 mmol sugar basis:

| basis assumed | HDMF mass | HDMF | **mol % on sugar** |
|---|---|---|---|
| reaction tube, 6.13 g | 6–31 µg | 0.05–0.24 µmol | **≈ 0.001 – 0.005 mol %** |
| diluted aqueous phase, 61 g | 61–307 µg | 0.5–2.4 µmol | **≈ 0.01 – 0.05 mol %** |

So the honest statement is **"of order 10⁻³ – 10⁻² mol % on sugar, spanning ~50×, on an ambiguous
basis, from a rough estimate with an uncorrected extraction recovery."** **This must NEVER enter a
benchmark file, a `mol %` field, or a fit target.** It is written here only so that (a) a future
reader does not re-derive it and think it is a finding, and (b) when the 1997 companion arrives
its numbers can be sanity-checked against the right order of magnitude. Recording it as a fit
target would repeat, exactly, the `cys_ribose_140C_Hofmann1998` fabrication that Wave S2b
retracted (`barrier_constants.py:307`).

### 3.4 The isotopomer distributions — **the paper's only real numbers**

**(i) Xylose + [2-¹³C]glycine** (p. 533, and the abstract):

> "The intensity of the molecular ions at *m/z* 128 and 129 revealed that **about 70% of HDMF was
> labeled**, thus indicating that HDMF was preferentially, but not exclusively, produced by
> incorporation of the labeled carbon of glycine. **The remaining 30% might be formed by sugar
> fragmentation.**"

Abstract restates it: *"The presence of ¹²C-HDMF, which was **approximately 30% of the total HDMF
amount** found in xylose/glycine, indicates that HDMF is **partly formed by sugar fragmentation**."*

**(ii) [1-¹³C]xylose + [2-¹³C]glycine and [1-¹³C]xylose + [3-¹³C]alanine** (p. 535) — the
double-labelling experiment, and the most informative numbers in the paper:

> "The total HDMF amount detected in the model reaction based on [1-¹³C]-D-xylose and [2-¹³C]glycine
> was composed of **unlabeled (4%), singly labeled (12%), and doubly labeled (84%)** HDMF. The
> corresponding figures for the system [1-¹³C]-D-xylose/[3-¹³C]-L-alanine were **15%, 47%, and 38%**,
> respectively. The data indicate that HDMF was generated by **different types of retro-aldol
> fragmentation of the pentose sugar**. **In the presence of glycine, HDMF was preferably formed
> involving the Strecker reaction.**"

Tabulated, with the interpretation the authors give:

| system | unlabelled HDMF | singly ¹³C | doubly ¹³C | reading |
|---|---|---|---|---|
| [1-¹³C]xylose + [2-¹³C]Gly | **4 %** | **12 %** | **84 %** | intact-C5 sugar backbone **+** Strecker-derived C1 dominates, **84 %** |
| [1-¹³C]xylose + [3-¹³C]Ala | **15 %** | **47 %** | **38 %** | Strecker channel **much weaker** with alanine; the C1 unit comes largely from sugar fragmentation instead |

**★ This is a measured, amino-acid-resolved BRANCH FRACTION at the C1-donor fork**, and it is the
single most model-relevant number in the paper. It is a *labelling* percentage, not a yield, so it
is immune to the response-factor / recovery caveat that kills §3.3 — the isotopomers of one
compound share a response factor exactly. Two independent measurements of the glycine case (70 %
labelled in (i); 84 + 12 = 96 % carrying at least one label in (ii)) are **consistent in sign and
roughly in size**, but they are **not the same quantity** (different labelled reactant sets) and
must not be averaged.

### 3.5 TABLE 2 — verbatim (p. 534). Mass-spectral, not kinetic.

> **Table 2. Mass Spectral Data of Nonlabeled and ¹³C-Labeled HDMF**

| compd | fragmentation pattern (*m/z*) | | | | | technique | ref |
|---|---|---|---|---|---|---|---|
| HDMF | 128 | 85 | 72 | 57 | 43 | GC-MS | Sen et al. (1991) |
| HDMF^a | 128 | — | 72 | — | — | GC-MS/MS | **this work** |
| ¹³CH₃-HDMF^b | 129 | 85/86 | 73 | 57/58 | 43/44 | GC-MS | Tressl et al. (1993) |
| [2,5-¹³CH₃]-HDMF | 130 | 86 | 74 | 58 | 44 | GC-MS | Sen et al. (1991) |

> ^a Parent scan of *m/z* 72 obtained by GC-MS/MS using commercially available HDMF.
> ^b MS data reported for a 1:1 mixture of two singly labeled isotopomers [2-¹³CH₃]HDMF and
> [5-¹³CH₃]HDMF.

**Three of the four rows are CITED, not measured here** (Sen 1991 ×2, Tressl 1993 ×1). Only the
GC-MS/MS parent-scan row is this work's own.

---

## 4. THE ISOTOPE EVIDENCE, EXPERIMENT BY EXPERIMENT

All labelling experiments used **xylose** as the pentose. Ribose and arabinose were never labelled.

### 4.1 Xylose + [3-¹³C]-L-alanine (p. 533)

Result: *"unlabeled HDMF and **monolabeled HEMF** tautomers."*
Evidence for the label position in HEMF: daughter ions **m/z 127 [M − ¹³CH₃]⁺** and **m/z 114
[M − CH₂ − ¹³CH₃]⁺** place the ¹³C **in the ethyl group**; the shifted series **58, 72, 87, 100**
(vs 57, 71, 86, 99 unlabelled) confirms incorporation; **[M − CO − CH₃]⁺ at m/z 100 (labelled) vs 99
(unlabelled)** and **[CO−CH₃]⁺ at m/z 43 for BOTH** show *"the methyl group of the labeled HEMF does
**not** bear the ¹³C atom."*
Conclusion, verbatim: *"the data strongly support the presence of a ¹³CH₃CH₂ group attached to C-2
of HEMF, i.e. **2-([2-¹³C]ethyl)-4-hydroxy-5-methyl-3(2H)-furanone**."*
→ **The C-3 carbon of alanine becomes the terminal carbon of HEMF's ethyl group.** Alanine enters
as **acetaldehyde**, a C2 unit.

### 4.2 Xylose + [2-¹³C]glycine (p. 533−534)

Result: mixture of monolabelled and unlabelled HDMF, **~70 % labelled** (§3.4(i)).
Position: daughter spectrum of *m/z* 129 gives **86, 73, 58, 44**; the ion pairs **43/44, 57/58,
85/86 of nearly equal intensity** correspond to *"the symmetric structure of HDMF which may occur as
a 1:1 mixture of two isotopomers."* The singly-labelled **[M − 2 × CO]⁺ at m/z 73** excludes C-3 and
C-4 as the label site (the unlabelled counterpart *m/z* 72 was shown by parent scan to come
exclusively from the molecular ion).
Conclusion: *"the ¹³C atom might be located in the methyl group, by analogy to HEMF preferentially
at C-2, i.e. **4-hydroxy-2-[¹³C]methyl-5-methyl-3(2H)-furanone**."*
→ **The C-2 carbon of glycine becomes a methyl carbon of HDMF.** Glycine enters as **formaldehyde**,
a C1 unit.

### 4.3 The double-labelling confirmations — the hypothesis test (p. 534−535)

The authors state the prediction *first*, then test it. This is the strongest evidence in the paper.

**(a) [1-¹³C]-D-xylose + [3-¹³C]-L-alanine → doubly labelled HEMF.**
> "the double-labeled HEMF tautomers were detected by GC-MS/MS … the fragment at *m/z* 100 in
> Figure 4A indicated the incorporation of one ¹³C into the C-2 ethyl group. On the other hand, the
> fragment at *m/z* 44, represented by the ion **[CO−¹³CH₃]⁺**, suggested the second labeling
> position in the **C-5 methyl group**. In concordance with the proposed formation pathway (Figure
> 3), the data strongly support the presence of **2-[2-¹³C]ethyl-4-hydroxy-5-[¹³C]methyl-3(2H)-furanone**."

**(b) [1-¹³C]-D-xylose + [2-¹³C]glycine → ¹³C₂-HDMF labelled in BOTH methyls.**
> "The mass spectrum shown in Figure 4B indicated the incorporation of two ¹³C atoms, i.e. *m/z* 130
> and 74. The fragment at *m/z* 44 ([CO−¹³CH₃]⁺) suggested that **both methyl groups were labeled,
> because the fragment at m/z 43 was not at all present** in the mass spectrum … thus suggesting
> that both methyl groups were labeled in the ¹³C₂-HDMF produced, i.e.
> **4-hydroxy-2,5-[¹³C]dimethyl-3(2H)-furanone**."

→ **Sugar C-1 becomes one methyl; the Strecker aldehyde carbon becomes the other methyl.** The
pentose C5 skeleton is carried through **intact** and is *extended*, not fragmented and reassembled,
in the dominant channel.

---

## 5. ★ THE MECHANISM, AS THE AUTHORS STATE IT (Figure 3 + text, p. 534). THIS IS THE DELIVERABLE.

### 5.1 The one-line stoichiometry the paper puts in display type (p. 534)

> **C5 + C1(2) → C6(7)**
>
> "C5 represents the pentose, **indicating that the pentose carbon chain remains intact**. The C1
> and C2 units may be **formaldehyde and acetaldehyde, the Strecker degradation products of glycine
> and alanine**, respectively."

### 5.2 Figure 3 legend, verbatim — the step labels

> **Figure 3.** Hypothetical formation pathway of HDMF (1) and HEMF (2) from pentoses (explanation
> in the text): **a**, early stage of Maillard reaction including Amadori rearrangement; **b**,
> degradation via 2,3-enolization; **c**, keto/enol or vinylogous keto/enol tautomerization; **d**,
> cyclization; **e**, dehydration; **f**, Strecker reaction; **g**, aldol condensation; **h**,
> reduction.

Note the word **"Hypothetical"** in the authors' own legend, and **"proposed"** in the abstract.
This is a **PROPOSED** mechanism supported by isotope evidence — not a measured one.

### 5.3 The route in full, in the authors' words (p. 534)

> "The formation mechanism proposed in Figure 3 is based on the formation of the **Amadori compound
> of the pentose** and its subsequent decomposition via **2,3-enolization** to form a **C5
> 1-deoxydiketose (4)**. This compound and other reductones react with the amino acids forming the
> corresponding aldehydes by **Strecker degradation, i.e. formaldehyde and acetaldehyde**. Further
> enolization of the 1-deoxyosone (4) may lead via compound (5) to the intermediate (6) which
> **reacts with the Strecker aldehydes via an aldol-type condensation**. The **chain elongation by
> the Strecker aldehydes** results in the intermediates **7a and 7b containing six and seven carbon
> atoms**, respectively. **Enolization and dehydration give rise to the cyclic intermediates
> acetylformoin (8a) and the corresponding ethyl homologue (8b).** The target molecules HDMF (1)
> and HEMF (2) are formed by **reduction, enolization, and water elimination** of the intermediates
> (8a) and (8b). **The reduction may occur either by a dismutation or by a reaction with further
> enoloxo compounds**, as recently reported by Schieberle (1992) for the formation of HDMF from (8a)
> in a model reaction based on hexoses."

**As a node-and-edge list (this is what a channel implementation needs):**

| # | node → node | transformation | carbons | paper's step letter |
|---|---|---|---|---|
| 1 | pentose + amino acid → **Amadori** | condensation + Amadori rearrangement | C5 + AA | a |
| 2 | Amadori → **C5 1-deoxydiketose (4)** = **1-deoxypentosone** | 2,3-enolisation | C5 | b |
| 3 | (4) + amino acid → **Strecker aldehyde** (HCHO from Gly, CH₃CHO from Ala) + aminoketone | Strecker degradation, dicarbonyl-driven | C1 / C2 released | f |
| 4 | (4) → (5) → **(6)** | further enolisation / vinylogous keto-enol tautomerisation | C5 | c |
| 5 | **(6) + Strecker aldehyde → 7a (C6) / 7b (C7)** | **aldol-type condensation — CHAIN ELONGATION** | C5+C1 / C5+C2 | g |
| 6 | 7a/7b → **acetylformoin (8a)** / **ethyl-acetylformoin (8b)** | enolisation + dehydration + cyclisation | C6 / C7 | d, e |
| 7 | **8a → HDMF (1)**, **8b → HEMF (2)** | **reduction** + enolisation + water elimination | C6 / C7 | h |

**★ ANSWER TO THE BRIEF'S KEY QUESTION.** The two candidate routes the brief posed —
*"acetylformoin? 1-deoxyosone + Strecker aldehyde?"* — **are not alternatives. They are steps 5−6
and step 7 of one route.** Acetylformoin **is** the immediate DMHF precursor (as Schieberle 1992
established for hexoses), and for **pentoses** acetylformoin is itself made by **1-deoxypentosone
+ Strecker aldehyde**. The model already holds the *input* node (`DPO` = 1-deoxypentosone,
`species_sulfur.py:99`) and already generates Strecker aldehydes; it holds **neither**
acetylformoin **nor** any pentose furanone product.

Definition given in the Introduction (p. 531), so the species is unambiguous:
> "acetylformoin **[2,4-dihydroxy-2,5-dimethyl-3(2H)-furanone]** is a key intermediate in the
> formation of HDMF from hexoses" *(attributed to Schieberle 1992).*

### 5.4 The SECOND, parallel route: sugar fragmentation (p. 535). Also stated by the authors.

> "The **C1 unit formaldehyde may originate from both Strecker degradation AND sugar fragmentation**,
> i.e. by **retro-aldol cleavage of 1-deoxyosones** (Ledl and Schleicher, 1990). This explains the
> presence of **unlabeled HDMF** in the model reactions containing only pentoses,
> xylose/[2-¹³C]glycine, and xylose/[3-¹³C]-L-alanine. **HDMF can be generated by recombination of
> pentose fragmentation products which may also include the fragments C2, C3, and C4 formed by
> retro-aldol reaction.**"

So DMHF has **two parallel channels from a pentose**, and §3.4 **quantifies their split**:

- **Channel A (Strecker chain-elongation, C5+C1):** dominant with glycine — 84 % doubly labelled.
- **Channel B (fragment recombination, C1 from retro-aldol of the 1-deoxyosone; possibly also
  C2/C3/C4 recombination):** ~30 % with glycine (§3.4(i)); **dominant-ish with alanine**, where the
  doubly-labelled share falls to 38 % and singly-labelled rises to 47 %. It is the **only** channel
  in the pentose-alone system.

Conclusions section, verbatim (p. 535):
> "The labeling experiments suggest that the furanones 1 and 2 are **mainly formed via Maillard
> reaction of pentoses with the amino acids glycine and alanine**, respectively (Figure 5). The C5
> moiety of the sugar is prolonged by a C1 or C2 unit, i.e. formaldehyde or acetaldehyde formed by
> Strecker degradation of the corresponding amino acids. **Alternatively, HDMF can also be produced
> without direct interaction of glycine.**"

### 5.5 The Strecker-aldehyde framing the authors close on (p. 535)

> "The results presented in this paper show the **manifold role of Strecker aldehydes**. They may
> contribute to flavor on their own (methional, phenylacetaldehyde, etc.) **or participate in the
> formation of other key odorants**. As shown in our laboratory, alanine and glycine are actively
> involved in the formation of sensory relevant 3(2H)-furanones (this paper) **and alkylpyrazines
> (Amrani-Hemaimi et al., 1995)**. **Although the furanones are formed by a side reaction**, these
> compounds may significantly contribute to the overall flavor due to their low threshold values."

*(Amrani-Hemaimi, Cerny & Fay 1995, JAFC 43, 2818−2822 — the same lab, the same tracer approach.
`k3_final_parameter_inventory.md` §A.10 already lists an "Amrani-Hemaimi on/off switch". This paper
is its furanone sibling and the two share the Strecker-aldehyde-as-building-block premise.)*

---

## 6. KINETICS — **NONE. NOT PARTIAL. NONE.**

No rate constant. No activation energy. No Arrhenius plot. No half-life. No plateau. No
concentration–time series. **One temperature (90 °C) and one time (1 h) in every experiment in the
paper**, including the labelled ones. There is no second time point anywhere from which even a
crude two-point rate could be built, and **no absolute concentration** for such a rate to be
expressed in. Nothing in this paper can be made into a kinetic parameter by any transformation.
Any wave that claims otherwise is over-reading it.

---

## 7. CITED, NOT MEASURED — carried separately so it is never mistaken for this paper's data

| quantity | value | true source, as this paper cites it | page |
|---|---|---|---|
| HDMF retronasal odour threshold, water | **160 µg/kg** | Huber (1992), *Perfum. Flavor.* 17(4), 15−19 | 531 |
| HEMF retronasal odour threshold, water | **20 µg/kg** | Huber (1992) | 531 |
| **Norfuraneol (HMF, 3)** retronasal odour threshold, water | **8300 µg/kg** | Huber (1992) | 533 |
| HEMF tautomer equilibrium **2a : 2b** | **1 : 2** | Re, Maurer & Ohloff (1973), *Helv. Chim. Acta* 56, 1882 | 531 |
| acetylformoin is the key HDMF intermediate from **hexoses** | qualitative | Schieberle (1992), ACS Symp. Ser. 490, 164−175 | 531, 534 |
| HDMF from 6-deoxysugars (rhamnose) | qualitative | Hodge 1963; Shaw & Berry 1976 | 531 |
| HDMF from fructose, from hexoses, from hexose-phosphates | qualitative | Mills & Hodge 1976; Shaw 1968; Fagerson 1969; Schieberle 1992 | 531 |
| pentoses degrade via 2,3-enolisation to **norfuraneol** | qualitative | Feather (1981) | 531 |
| retro-aldol cleavage of 1-deoxyosones gives formaldehyde | qualitative | Ledl & Schleicher (1990) | 535 |
| 6,7-dideoxyheptose as a speculative HEMF intermediate | **speculation, explicitly** | Huber (1992) — *"Huber (1992) **speculated** on the existence of…"* | 531 |

**⚠️ The three odour thresholds are Huber 1992, a trade-magazine article (*Perfumer & Flavorist*),
reached through this paper.** They are **secondary** here. If the repo wants a furaneol threshold
for §A.7, it should be taken from a primary source, and it must **not** be attributed to Blank &
Fay 1996. Note the *ratio* form survives provenance better (`k3` §A.7.3): **HEMF : HDMF : norfuraneol
= 1 : 8 : 415** on Huber's own scale, single-source and therefore internally consistent.

**The one directional fact this section supports and the repo can use** (p. 533):
> "**HMF (3), which was the major volatile in the pentose/amino acid systems**, did not contribute
> much to the caramel-like aroma."

→ **In a pentose + Gly/Ala system at 90 °C, norfuraneol ≫ DMHF and norfuraneol ≫ HEMF by mass.**
Combined with "minor peaks" (§3.2) and "formed by a side reaction" (§5.5), this is a **hard ordinal
ceiling on any pentose DMHF channel the model might add**, and it is exactly the constraint that
`k3` §B3.7 (Nedvidek's DPO partition) already warns about from the other direction.

---

## 8. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured in this paper · **[C]** cited from another paper · **[P]** proposed /
hypothetical by the authors · **[D]** derived by me during this extraction.

| # | quantity | value | conditions | class | verdict |
|---|---|---|---|---|---|
| 1 | DMHF **route topology** from pentose: Amadori → 2,3-enolisation → 1-deoxypentosone → (+Strecker aldehyde, aldol) → **acetylformoin** → (reduction) → **DMHF** | route, 7 steps | 90 °C, pH 6→5, 0.2 M PO₄ | **[P]**, isotope-supported **[M]** | **★ USE as topology.** The paper's central result. |
| 2 | HEMF route: identical, with **acetaldehyde** in place of formaldehyde → ethyl-acetylformoin → HEMF | route | same | **[P]**/**[M]** | **★ USE as topology.** |
| 3 | **C5 + C1 → C6** / **C5 + C2 → C7** stoichiometry; pentose backbone **intact** | exact | same | **[M]** (double labelling) | **★ USE.** Falsifiable atom-mapping constraint. |
| 4 | DMHF label distribution, [1-¹³C]xylose + [2-¹³C]Gly | **4 % / 12 % / 84 %** (un-/single-/double-) | 90 °C, 1 h | **[M]** | **★ USE — the paper's best number.** Branch fraction, response-factor-immune. |
| 5 | DMHF label distribution, [1-¹³C]xylose + [3-¹³C]Ala | **15 % / 47 % / 38 %** | 90 °C, 1 h | **[M]** | **★ USE.** Pairs with #4 as an amino-acid contrast. |
| 6 | ¹³C-labelled fraction of DMHF, xylose + [2-¹³C]Gly | **~70 % labelled / ~30 % unlabelled** | 90 °C, 1 h | **[M]** | USE, but **do not merge with #4** — different labelled-reactant set. |
| 7 | HEMF requires **alanine**; absent with glycine; absent from pentose alone | on/off | 90 °C, 1 h, GC-O | **[M]** ordinal | **★ USE — clean switch.** |
| 8 | DMHF forms from **pentose alone**, no amino acid | on/off | 90 °C, 1 h, GC-O | **[M]** ordinal | **★ USE — a zero-amino-acid positive control.** |
| 9 | **norfuraneol ≫ DMHF, HEMF** ("major volatile" vs "minor peaks", "side reaction") | ordinal | 90 °C, 1 h | **[M]** ordinal | **★ USE as a CEILING on any pentose DMHF channel.** |
| 10 | GC-O intensity grid (Table 1) | `++`/`+`/`+`; `−`/`+++`/`−`; `−`/`+`/`−` | 90 °C, 1 h | **[M]** ordinal | USE as on/off only. **Not a concentration ratio.** |
| 11 | LRI: HDMF 1235 (OV-1701) / 2040 (FFAP); HEMF 2a 1310/2090; 2b 1325/2190 | — | — | **[M]** | Identification support only. |
| 12 | HDMF / HEMF amounts | **"low mg/kg (ppm)", rough estimate** | 90 °C, 1 h | **[M]**, self-labelled rough | **⚠️ DO NOT INGEST.** No method, no basis, no recovery, no breakdown. §3.3. |
| 13 | implied mol % on sugar | **~10⁻³ – 10⁻² mol %, ambiguous basis, ~50× span** | as above | **[D]** | **⚠️ DO NOT INGEST — my arithmetic, not a measurement.** Scale check only. |
| 14 | initial pH / final pH | **6.0 → 5.0 (Gly)**, **6.0 → 5.3 (Ala)** | 90 °C, 1 h, 0.2 M Na₂HPO₄ | **[M]** | **USE as a drift constraint**, and as the label caveat (§1.2). |
| 15 | phosphate | **0.2 mol/L Na₂HPO₄**, one level | — | **[M]** | **No phosphate series exists.** Single point. |
| 16 | loading | **1.0 M pentose, 1.0 M amino acid, 1:1** | 5 mmol / 5 mmol / 5 mL | **[M]**/**[D]** | Context for #4−#9. **No dose–response.** |
| 17 | acetylformoin = 2,4-dihydroxy-2,5-dimethyl-3(2H)-furanone | structure | — | **[C]** Schieberle 1992 | USE for the species definition. |
| 18 | odour thresholds HDMF 160 / HEMF 20 / norfuraneol 8300 µg/kg water, retronasal | — | — | **[C]** Huber 1992 | **Secondary.** Ratio form 1 : 8 : 415 is safer. |
| 19 | HEMF tautomer ratio 2a:2b = 1:2 | — | — | **[C]** Re 1973 | Secondary. |
| 20 | activation energy, rate constant, half-life, time course, temperature series | **NONE** | — | — | **ABSENT. §6.** |

---

## 9. VERDICT — channel, or yield snapshot?

**Neither, exactly. It is a ROUTE paper with a BRANCH-FRACTION anchor and NO yield anchor.**

Scored against the brief's own test — *"formation route + at least one quantitative anchor"*:

| requirement | met? | detail |
|---|---|---|
| **Formation route, fully specified** | **✅ YES, and better than expected.** | Seven steps, every one named by the authors, every intermediate identified, the C5+C1/C5+C2 topology confirmed by **double** labelling with the prediction stated before the test. Acetylformoin is pinned as the immediate precursor. A second parallel channel (fragment recombination) is named and its share measured. |
| **At least one quantitative anchor** | **⚠️ YES, but not the kind that closes a channel.** | The anchors are **isotopomer branch fractions** (#4, #5, #6) and **ordinal on/off facts** (#7−#10). They constrain **where the carbon comes from and in what proportion** — they do **not** constrain **how much furanone is made**. |
| **A yield the model could be scored against** | **❌ NO.** | Item #12 is a self-declared rough estimate on an unstated basis. There is **no mol %, no ppb, no µg/kg with a basis, no internal standard, no recovery correction, and no per-system number.** |
| **Kinetics** | **❌ NO.** | §6. Single T, single t. |
| **pH / phosphate dependence** | **❌ NO.** | Single buffer, single concentration, uncontrolled drift. |

**So: this paper ENABLES a pentose DMHF/HEMF channel TOPOLOGICALLY and CONSTRAINS ITS INTERNAL
BRANCHING, but it CANNOT CALIBRATE IT.** A channel added on this paper alone would have a
correct, isotope-evidenced structure and a **completely unconstrained magnitude** — with only the
ordinal ceiling `norfuraneol ≫ DMHF` (#9) preventing it from being arbitrarily large. **That is a
real and useful state to be in** — it is strictly better than the repo's current position, which is
that the pentose furanone channel does not exist at all (defect **D2**) — but it is **not** the
"quantified furaneol channel" the brief was hoping for. **The magnitude has to come from the 1997
IDA companion (§0a), and that paper is not in the repo.**

**One-line honest answer:** *this is the mechanism paper, not the numbers paper; the numbers paper
is Blank, Fay, Lakner & Schlosser 1997 and it must be fetched before any pentose furanone channel
can be calibrated or scored.*

---

## 10. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION. DO NOT EDIT THE DECLARATION.**

Nothing below has been written to any registry, benchmark file, or config. These are proposals.

### 10.1 Proposed **FIT** roles — *none for magnitude*

| candidate | proposed role | why |
|---|---|---|
| Items #4/#5 (isotopomer branch fractions) | **FIT — as a PRIOR on the C1-donor branch ratio**, if and only if a pentose DMHF channel is implemented with both Channel A (Strecker) and Channel B (fragment recombination). Target: with glycine, Channel A ≥ ~70−84 % of DMHF flux; with alanine, Channel A ≈ 38−47 %. | Response-factor-immune, amino-acid-resolved, and it is a *ratio*, so it survives the absence of an absolute scale. It is the only thing here that can legitimately move a parameter. |
| **NO barrier should be fitted against this paper.** | — | There is no magnitude to fit to. Fitting a barrier against item #12 would reproduce the `thiol_addition_pentodiulose` / `cys_ribose_140C_Hofmann1998` failure mode verbatim (`barrier_constants.py:307`, Wave S2b) — a constant fitted to a number that is not a measurement. **Record as a prohibited derivation.** |

### 10.2 Proposed **HOLD-OUT** roles

| candidate | proposed role | why |
|---|---|---|
| #7 — HEMF **requires alanine** (zero with glycine, zero with pentose alone) | **HOLD-OUT, structural.** Three-cell truth table; a model that emits HEMF from pentose+glycine fails. | Free of any calibration. Pure topology. Falsifiable in one run. |
| #8 — DMHF forms from **pentose alone** | **HOLD-OUT, structural.** A zero-amino-acid positive control for the fragment-recombination channel. | The repo has very few zero-precursor controls (`k3` §C.4 complains about exactly this absence in the extrusion lane). |
| #9 — **norfuraneol ≫ DMHF and ≫ HEMF** in pentose/Gly and pentose/Ala at 90 °C | **HOLD-OUT, ordinal ceiling.** | Directly guards the failure mode `k3` §B3.7 names: over-funnelling the DPO pool. If a new DMHF channel inverts this ordering, the channel is mis-sized and the test says so **without needing an absolute number**. **This is the single most valuable hold-out in the paper.** |
| #3 — pentose backbone **intact**, C5+C1/C5+C2 | **HOLD-OUT, atom-mapping.** Pins that all five sugar carbons reach DMHF and the sixth comes from the Strecker aldehyde. | Same class as the existing Wang & Ho intact-C6 CAMOLA pin (`test_chemistry_soundness.py`), and it is its **pentose counterpart** — the repo currently pins the hexose case only. |
| #14 — pH drift 6.0 → 5.0 / 5.3 | **HOLD-OUT, weak.** Only if the model has a pH-drift prediction for 1 M sugar + 1 M amino acid at 90 °C. | Two points, one system. Low information, but free. |
| #12 / #13 (the "low mg/kg") | **NEITHER. Excluded from both roles.** | §3.3. |

### 10.3 Proposed **PREREQUISITE** — do this before either role is declared

1. **Fetch Blank, Fay, Lakner & Schlosser 1997, JAFC 45, 2642−2648 (`10.1021/jf960997i`, DOI
   unverified).** Without it there is no magnitude anchor for this channel, and a channel with a
   free magnitude and one ordinal ceiling is a **liability in a validation panel**, not an asset.
2. **Resolve defects D1, D2, D3 (§0c) first.** D2 in particular: `audit_remediation.md:1713`
   currently records this route as "CONFIRMED … the model uses", and it is not implemented. That
   record should be corrected *before* anyone reads it as evidence that the channel already exists.
3. If the channel is implemented, it needs **acetylformoin as a new explicit species**
   (2,4-dihydroxy-2,5-dimethyl-3(2H)-furanone, C₆H₁₀O₄) and its **ethyl homologue** (C₇H₁₂O₄).
   `species_sulfur.py` has `DPO` (the input) and `NF` (the competitor) already; the intermediate
   and the two products are all missing.

---

## 11. DECLARED GAPS — verbatim, for `k3` §C

> **"There is no yield table in this paper.** The only statement of amount is *'Rough estimates of
> the amounts of HDMF and HEMF formed from pentoses were in the low milligrams per kilogram (parts
> per million) range'* (p. 533) — a self-declared rough estimate with **no method, no basis of
> reference, no recovery correction, no internal standard, no uncertainty and no per-system
> breakdown**, published in the same paragraph in which the authors defer quantification to a
> future study. **It is not a measurement and it must never become a fit target or a benchmark
> row.**"

> **"One temperature, one time, one pH, one phosphate concentration, one loading.** 90 °C / 1 h /
> pH 6.0 initial / 0.2 mol/L Na₂HPO₄ / 1 M : 1 M. **No Arrhenius data, no time course, no
> dose–response, no pH series and no phosphate series exist in this paper.** The brief's premise
> that pH and phosphate were varied here is **wrong**; that is the 1997 companion."

> **"The three pentoses are never resolved.** Table 1's columns lump ribose, xylose and arabinose
> under one heading, and **every labelling experiment used xylose only**. **No sugar-reactivity
> ordering can be extracted from this paper**, and any ribose-vs-xylose claim sourced to it is
> invented."

> **"The odour thresholds are Huber 1992, a trade-magazine secondary source reached through this
> paper.** 160 / 20 / 8300 µg/kg water, retronasal. **They are not this paper's measurements** and
> must not be attributed to Blank & Fay 1996."

> **"The mechanism is labelled 'Hypothetical' by the authors in their own figure legend and
> 'proposed' in their own abstract.** The isotope evidence constrains the carbon topology tightly
> and the double-labelling prediction was stated before it was tested — but **no intermediate
> between the Amadori compound and the product was isolated or quantified**, acetylformoin
> included. Steps 4−7 of §5.3 are inferred, not observed."

> **"Blank & Fay 1996 do NOT name the amino acid as the reductant.** They name *'either a
> dismutation or a reaction with further enoloxo compounds'* (p. 534). Any repo text attributing
> the amino-acid-as-reductant claim to `10.1021/jf950439o` is an over-read of the abstract's
> reductant-free phrase 'reduction of the resulting acetylformoin-type intermediates'."

---

## 12. WHAT TO FETCH NEXT — ranked

| # | paper | why | confidence in the citation |
|---|---|---|---|
| 1 | **Blank, Fay, Lakner & Schlosser 1997**, *"Determination of 4-Hydroxy-2,5-dimethyl-3(2H)-furanone and 2(or 5)-Ethyl-4-hydroxy-5(or 2)-methyl-3(2H)-furanone in Pentose Sugar-Based Maillard Model Systems by **Isotope Dilution Assays**"*, JAFC 1997, 45, 2642−2648, `10.1021/jf960997i` | **The paper the brief actually described.** Same lab, same 90 °C / 1 h system, with SIDA quantification and the "parameters affecting this reaction" study the 1996 paper promises. **This is the only route to a magnitude anchor for the pentose furanone channel.** | Title/authors/journal/year: high (~95 %). DOI + pagination: **~85 %, web lookup only, not CrossRef-verified, not read.** |
| 2 | **Schieberle 1992**, *"Studies on the formation of furaneol in heat-processed foods"*, ACS Symp. Ser. 490, 164−175 | The primary source for **acetylformoin as the key HDMF intermediate** and for the reduction step's mechanism (dismutation / enol-oxo). Cited twice here and load-bearing for the last step of the route. Book chapter, no DOI given. | High on identity; availability unknown. |
| 3 | **Amrani-Hemaimi, Cerny & Fay 1995**, JAFC 43, 2818−2822, *"Mechanisms of formation of alkylpyrazines in the Maillard reaction"* | The sibling tracer study from the same lab, already named in `k3` §A.10 as an "on/off switch". Shares the Strecker-aldehyde-as-C-donor premise. | High. |
| 4 | **Huber 1992**, *Perfum. Flavor.* 17(4), 15−19 | The actual source of the HDMF/HEMF/norfuraneol thresholds. **Trade magazine** — likely to be a dead end for a primary threshold, in which case the repo should source furaneol's threshold elsewhere and say so. | High on identity, low on usefulness. |
| 5 | **Sen, Schieberle & Grosch 1991**, *LWT* 24, 364−369 | The **stable-isotope-dilution assay for HDMF** itself — the method the 1997 companion applies. Useful if the companion's numbers need method scrutiny. | High. |

---

## 13. SOURCES USED IN THIS EXTRACTION BEYOND THE PDF

The companion paper's identity (§0a, §12 item 1) was established by web search during this
extraction and is **not** verified against CrossRef. Listings consulted:

- [Determination of … in Pentose Sugar-Based Maillard Model Systems by Isotope Dilution Assays — ResearchGate](https://www.researchgate.net/publication/248744102_Determination_of_4Hydroxy25-dimethyl-32_H_-furanone_and_2or_5Ethyl4-hydroxy-5or_2-methyl-32_H_-furanone_in_Pentose_Sugar-Based_Maillard_Model_Systems_by_Isotope_Dilution_Assays)
- [Same title — Semantic Scholar](https://www.semanticscholar.org/paper/Determination-of-and-2(or-5)-Ethyl-4-hydroxy-5(or-Blank-Fay/f90ff430befe3c3a86ce84c09f15cfaa80df9ae7/figure/2)
- [Blank & Fay 1996 (this paper) — ACS Publications](https://pubs.acs.org/jafcau/article-abstract/44/2/531/69763/Formation-of-4-Hydroxy-2-5-dimethyl-3-2H-furanone)

*End of dossier. Nothing outside this file was created or modified.*
