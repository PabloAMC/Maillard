# van den Ouweland & Peer 1975 — EXTRACTION
DOI: 10.1021/jf60199a045 (not printed on the pages; Crossref 2026-11-01: "Components contributing to beef flavor. Volatile compounds produced by the reaction of 4-hydroxy-5-methyl-3(2H)-furanone and its thio analog with hydrogen sulfide", Van den Ouweland & Peer, J. Agric. Food Chem. 23(3), 501-505, 1975-05)
### Fed intermediate: norfuraneol (4-hydroxy-5-methyl-3(2H)-furanone, 0.5 g in 30 mL water) or its thiophenone analogue + H2S (15 g, ~100x molar), glass-lined autoclave, 95-100 °C, 4 h; no amine, no buffer, no pH. Qualitative identification only (GC-IR-MS of trapped peaks against synthesised references); NO yields. **The aromatic furanthiol identified is printed as 4-mercapto-2-methylfuran (SH NOT next to the methyl), not 2-methyl-3-furanthiol (MFT).** The MFT-pattern products are the dihydrofuran (3-mercapto-2-methyl-4,5-dihydrofuran) and the thiophenes.

**Read 2026-11-01 from `data/articles/vandenouweland1975.pdf`** (5 pp., JAFC pp. 501-505; Acrobat Capture OCR;
p. 501 starts with the tail of an unrelated sesquiphellandrol paper). The rotated Tables I-II, Figs. 1-4 and every
structure were read from 130 / 220 dpi renders (scratchpad `ing1101/vo-*.png`, `vohi-*.png`); the text layer of
both tables is unusable. The GC traces (Figs. 2-3) are recorder plots without a quantitative scale: no peak was
digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a
figure · **[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived
by me, never printed · **[!]** flag.

**Laboratory.** Unilever Research Duiven, Zevenaar, The Netherlands (van den Ouweland, Peer). `lab_ids.yml` has
`unilever_vlaardingen` (Potman & van Wijk), a different Unilever site and group; the measurement ledger already
uses the string "Unilever (van den Ouweland)" for their patent US 3,904,655 (`us3904655_extraction.md`). Proposed:
a new id `unilever_duiven_van_den_ouweland`, merged with the patent (same authors, same chemistry), kept apart
from Vlaardingen.

## 1. What this paper answers for the model

The question: a fed-intermediate test of the norfuraneol (NF) route to MFT (engine `k_nf_mft`; atlas F23 free-NF
sub-route; E6).

| model need | what the paper gives | tag |
|---|---|---|
| does free NF + H2S make MFT? | **not as printed.** Table I peak 3 is "4-Mercapto-2-methylfuran" (m/e 114, base peak), its structure drawn with SH on the ring carbon remote from the methyl, confirmed against a synthesised reference (footnote a). No 2-methyl-3-furanthiol is named | [M] / [!] |
| MFT-pattern products (SH next to the methyl) | 3-mercapto-2-methyl-4,5-dihydrofuran (peak 4, M 116, "roasted meat"), 3-mercapto-2-methylthiophene (peak 10, M 130, "roasted meat"), 3-mercapto-2-methyl-2,3- and -4,5-dihydrothiophene (peaks 11, 17), 3-mercapto-4-hydroxy-2-methyl-2,3-dihydrothiophene (14) | [M] |
| yields / rate / mol % | **none**: no amounts, no internal standard; Figs. 2-3 are FID traces in "recorder response" | [!] |
| NF conversion | NF (peak 12) remains the largest peak of Fig. 2 after 4 h with ~100x H2S [F, qualitative] | [F] |
| conditions relevant to food | none: no amine, no buffer, 100x H2S "a much higher concentration ... than ... a natural beef broth" [Q] | [Q] |
| temperature / pH dependence | none (one condition) | — |

**Answers.**

1. **As a fed-intermediate test of NF -> MFT, the paper shows NF + H2S -> sulfur-substituted 2-methylfurans /
   thiophenes at 95-100 °C, but NOT MFT itself as printed** [M]. Confidence that the print names no MFT: ~95 %.
   Whether peak 3 is in fact MFT misnamed: P ~0.3 [Z]. Against misnaming: the name, the drawn structure and the
   footnote (structure confirmed with a synthesised reference) agree, and the paper names the 3-mercapto-2-methyl
   pattern correctly for peaks 4, 10, 11 and 17 in the same table. For it: both isomers are chemically accessible
   (thionation at the C-3 carbonyl gives the 4-mercapto isomer, at the C-4 enol carbon via the open dione of Fig. 1
   gives MFT [Z]), and the paper has often been cited as an MFT source. The MS (m/e 114 100, 113 50, 85 60, 71 45)
   gives [M-H]+/M+ 0.50, against 0.28 for MFT's NIST EI spectrum (the value used in CF-086): mild evidence against
   MFT [Z], instrument-dependent.
   **Added 2026-11-01 (resolve2_sulfur and its red team; CF-092). P ~0.3 is kept.** Two further readings: (i) the
   paper's own **Fig. 4** (p. 505) POSTULATES the NF -> MFT pathway (H2S replaces the enol OH next to the methyl ->
   4-mercapto-5-methyl-3(2H)-furanone -> reduction / dehydration -> 2-methyl-3-furanthiol and its dihydro /
   tetrahydro forms), the pathway Hofmann 1998 Fig. 1 prints; peak 3's 4-SH assignment contradicts the paper's own
   mechanism. So the paper is a legitimate citation for the PROPOSED NF -> MFT pathway, not for MFT IDENTIFIED. (ii)
   Hofmann & Schieberle 1998 (intro) state that van den Ouweland and Peer were unable to identify MFT; Whitfield &
   Mottram 1999 attribute MFT to them twice (intro and discussion). Against "peak 3 is MFT": four later labs make MFT
   from NF and none reports 2-methyl-4-furanthiol. For it: Fig. 4's regiochemistry; peak 3's IR lists no S-H band
   (~2550), unlike peaks 7, 8 and 10, a small flag on the identification as a whole. Not a misnaming in the print
   (P <= 0.05). Unknown: whether peak 3's reference compound was synthesised by a route that fixes the
   regiochemistry (not printed).
2. **No rate, no yield, no share**: nothing here bounds `k_nf_mft` or the free-NF share (E6). The paper is a
   qualitative licence that NF + H2S, without any amine, reaches the 2-methyl-3-thio skeleton (as dihydrofuran and
   thiophene), at a 100x H2S excess.
3. **Authors' own reading [Q]**: the meaty compounds are "probably not formed via a Maillard-type reaction" (NF from
   ribose 5-phosphate by dephosphorylation-dehydration) and none had yet been identified in meat (1975).

## 2. Conditions

| item | as printed / derived |
|---|---|
| substrate | 4-hydroxy-5-methyl-3(2H)-furanone (NF, compound 1, R = H) 0.5 g in water 30 mL = ~146 mmol/L [Z: M 114.1]; or 4-hydroxy-5-methyl-3(2H)-thiophenone (2, R = H), "exactly the same way" |
| sulfur | H2S 15 g = ~440 mmol, ~100x molar over NF [Z] |
| amine / buffer / pH | none / none / not stated |
| vessel | glass-lined autoclave (closed; headspace and O2 not stated) |
| T / time | "4 hr at 95-100°" |
| work-up | poured into 100 mL ice water, five CH2Cl2 extractions (125 mL), concentrated at atmospheric pressure to 10 mL |
| analysis | packed GLC (600 x 0.4 cm, 10 % Apiezon L-1 % Carbowax 20M), FID with 10 % split, 60 -> 220 °C at 5 °C/min (4 °C/min for the thiophenone run); trapped peaks for MS and IR (Copier & van de Maas method); odour at the exhaust |

## 3. Quantitative results

None in the kinetic sense. Printed per peak: name, structure, MS (m/e and relative intensity), IR bands, odour
(Tables I-II) [M]. Table I (NF + H2S), peaks 1-17: 2-methylthiophene; 2-methyl-2,3(or 2,5)-dihydrothiophene;
**4-mercapto-2-methylfuran** (114: 100, 113: 50, 85: 60); 3-mercapto-2-methyl-4,5-dihydrofuran (116: 70);
1-methylthio-2-butanone; 3-hydroxy-2-methyltetrahydrofuran; 4-mercapto-3-oxotetrahydrofuran;
4-mercapto-2-methyl-2,3-dihydrothiophene; 4-mercapto-5-methyl-3-oxotetrahydrofuran; **3-mercapto-2-methylthiophene**
(132: 9, 130: 100, 129: 50); 3-mercapto-2-methyl-2,3-dihydrothiophene; NF (starting material);
4-mercapto-2-methylthiophene; 3-mercapto-4-hydroxy-2-methyl-2,3-dihydrothiophene;
4-mercapto-2-methyl-4,5-dihydrothiophene; 4-hydroxy-5-methyl-3(2H)-thiophenone; 3-mercapto-2-methyl-4,5-
dihydrothiophene. Footnote a (structures checked against synthesised references) on peaks 1, 3, 4, 5, 8, 9, 10,
12, 13, 15, 16, 17. Table II (thiophenone + H2S), 16 peaks, including 2,3-pentadione, 2,3-dimercaptobutane,
4-mercapto-5-methyl-3(2H)-furanone, cis / trans-3-mercapto-2-methyltetrahydrothiophene, 4-mercapto-5-methyl-3(2H)-
thiophenone [M].

## 4. Mechanism evidence and its kind

- **Fed intermediate (kind 2), qualitative**, one laboratory, one condition, no controls (no NF-free blank, no
  amine arm). Fig. 1: NF ⇌ open 2,4-diketone, which with H2S closes to the thiophenone [P]. Fig. 4: H2S replaces
  the enol OH, then reduction / dehydration to the mercapto-dihydro and mercapto-aromatic furans / thiophenes [P].
- Identification is by MS + IR against synthesised references for the footnoted peaks [M]; the rest by spectra only.

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| E6 (MFT 99.9 % via free NF in the engine; measured 7 %) | **LEAVES UNCHANGED** | no yield; and no MFT by name, so not even a qualitative positive for NF -> MFT |
| F23 free-NF sub-route (atlas) | weakly SUPPORTS that free NF + H2S reaches the 2-methyl-3-thio skeleton (dihydrofuran, thiophene) without an amine | qualitative, 100x H2S |
| engine `k_nf_mft` (NF + H2S -> MFT) | LEAVES UNCHANGED | no rate; and the paper's direct product is the dihydro / thiophene, not MFT |
| CF-086 / F24 | not touched | no labelling |

## 6. Disagreements found

1. **With the repo's legacy descriptions** (`docs/history/AUDIT_legacy_lane_2026-08.md` "describes a synthesis of MFT
   from norfuraneol"; `docs/validation/isotope_topology_evidence.md` row 12): the print names
   4-mercapto-2-methylfuran, not MFT. Quantity: identity of the aromatic furanthiol from NF + H2S. Candidate causes:
   (a) the citation trail treats the paper as an MFT source without reading Table I; (b) a misnaming in the paper
   (P ~0.3, sec. 1).
2. **With Whitfield & Mottram 1999** (`whitfield1999_extraction.md`; Reading; NF 50 mM + H2S ~97 mM, 0.5 M phosphate
   pH 4.5, 140 °C, 60 min): Whitfield measures MFT from NF + H2S and finds **none** of the di- or tetrahydro
   derivatives of MFT or 2-methyl-3-thiophenethiol; van den Ouweland finds the dihydro derivatives (peaks 4, 11, 17)
   and no MFT by name. Candidate causes: H2S excess (~100x here, ~2x there: a strongly reducing medium favours the
   dihydro products), no buffer vs phosphate pH 4.5, 95-100 °C / 4 h vs 140 °C / 1 h, packed vs capillary GC
   (co-elution of MFT with peak 3 or 4), identification method.

## 7. Ledger-ready measurements

1. **Table I identities** (NF + H2S, 95-100 °C, 4 h, no amine): a qualitative "operating" row for NF + H2S ->
   3-mercapto-2-methyl-4,5-dihydrofuran and 3-mercapto-2-methylthiophene; an "absent as printed" note for MFT,
   with the misnaming caveat. Two readers: re-check the structure of peak 3 against the name.
2. Nothing numeric is worth keying (no amounts).
