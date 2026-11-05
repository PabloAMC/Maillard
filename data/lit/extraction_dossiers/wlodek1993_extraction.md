# Wlodek, Rommelspacher, Susilo, Radomski & Höfle 1993 — EXTRACTION
DOI: 10.1016/0006-2952(93)90632-7 (PII printed as the PDF title; Crossref 2026-11-05: "Thiazolidine derivatives as source of free L-cysteine in rat tissue", Wlodek, Rommelspacher, Susilo, Radomski, Höfle, Biochemical Pharmacology 46(11), 1917-1928, 1993-12. **Matches the brief.**)
### Fourteen thiazolidines, among them FIVE cysteine-PENTOSE rings (L-xylose, D- and L-arabinose, D-ribose, D-lyxose), 1.5 and 2.5 mmol/L in 0.1 M phosphate **pH 7.4, 37 °C, ONE time point (30 min)**: "non-enzymatic hydrolysis" read as DTNB-reactive SH (Fig. 2, bar chart, mean + SD, n = 5-10). The rest of the paper is rat-mitochondria metabolism. **One pH, one temperature, no time course, no K.**

**Read 2026-11-05 from `data/articles/wlodek1993.pdf`** (12 pp., Acrobat Capture scan with a usable OCR layer for prose;
figures and structure drawings garbled in the OCR). Prose from `pdftotext -layout`; **Fig. 2 (p. 1923) digitised from a
300 dpi render** (scratchpad `m1c/f2.png`); Figs. 1, 3-5 read at 120 dpi for content only. Read error on Fig. 2:
**+-1.5 nmol/mL** on bar tops, +-2 nmol/mL on SD caps. Reader: one reader (Claude, first reader). A second reader is
required before any value here enters a ledger, a band or a fit row; none is proposed.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by
me · **[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Lidia Wlodek (Medical Academy, Krakow), Hans Rommelspacher* and Rudy Susilo (Free University Berlin), Jan Radomski (Warsaw University, Chemistry), Gerhard Höfle (GBF Braunschweig) |
| Title | "Thiazolidine derivatives as source of free L-cysteine in rat tissue" |
| Venue | Biochemical Pharmacology 46(11), 1917-1928 (1993); received 13 Nov 1992, accepted 13 Aug 1993 |
| DOI | 10.1016/0006-2952(93)90632-7 (PII on the PDF; Crossref) |
| Sugar thiazolidines | made by Radomski & Temeriusz 1989 (ref. 40; `radomski1989_extraction.md`, method A) or a modified Weitzel 1959 (ref. 41); stereochemistry by polarimetry + 1H NMR (ref. 40). Listed (Structure 1): 6 CGlu(D), 7 CAra(D) 2(S)-D-arabino, 8 CAra(L) 2(R)-L-arabino, 9 CFuc(D), 10 CFuc(L), 11 CGal(D), 12 CRib(D) 2(S)-D-ribo, 13 CXyl(L) 2(S)-L-xylo-tetrahydroxybutyl |

## 1. What this paper answers for the model

| model need | what Wlodek 1993 gives | tag |
|---|---|---|
| pentose-thiazolidine release rate vs pH (6-8) | **one pH only (7.4)**, 0.1 M phosphate, 37 °C, 30 min: DTNB-SH released from CXyl(L) **40 / 55 nmol/mL** (1.5 / 2.5 mmol/L charge), CAr(D) 38 / 55, CRib(D) 16 / 40, CAr(L) 12 / 28, CLyx(D) 16 / 21 | [F] |
| the same as fractions and apparent k | 0.8-2.7 % of the ring in 30 min; apparent first-order 2.7e-4 - 9.0e-4 /min (sec. 3) | [Z] |
| "cysteine measured"? | **total DTNB-reactive thiol**, read 1 min after adding 2 mmol/L DTNB; not free cysteine specifically. A ring opened to its open-chain thiol form, and cystine-free cysteine, both read | [M] / [!] |
| temperature law | none (37 °C only) | -- |
| K | none (one time point; equilibrium not tested) | -- |
| bis-Cys (2 Cys : 1 aldehyde) pool | **methyl-djenkolic acid (MDA)** from the ACETALDEHYDE ring CA only, identified by 1H NMR of its bis-DNP derivative (Table 2), in buffer with or without mitochondria, pH 7.4, 37 °C. "The lower the pH the more MDA is formed from CA" [Q], "preliminary experiments", no numbers. **No bis-Cys species is reported for any sugar ring** | [M] / [Q] |
| pentose rings and enzymes | sugar thiazolidines are not oxidised by liver or brain mitochondria (Fig. 1, A250 flat) | [Q] / [F] |

**Answer.** The paper's chemistry for the model is one bar chart at pH 7.4 and 37 °C. It places five pentose rings'
SH release at 0.8-2.7 % in 30 min, between the stable formaldehyde ring CF (0.3-0.4 %) and the labile glucose and
acetaldehyde rings (3.1-3.8 %). It says nothing about pH dependence, temperature dependence or equilibrium. The bis-Cys
pool named in the queue description exists here only for acetaldehyde. Confidence that nothing else in the paper bears
on thiol supply in a heated cysteine + pentose pot: ~0.9.

## 2. Conditions [M]

| item | as printed |
|---|---|
| Non-enzymatic hydrolysis (Fig. 2) | TD 1.5 or 2.5 mmol/L, **0.1 M phosphate buffer pH 7.4, 37 °C, 30 min, 1 mL final volume** (Fig. 2 legend). Methods: after 30 min, 1 mL of 2 mM DTNB added, stirred, **A412 read 1 min later** (Zeiss PMQ photometer). "Values are means and SD of 5-10 experiments" |
| pH | buffer label 7.4; no measured pH; the drift of 1.5-2.5 mmol/L solute in 0.1 M phosphate is negligible [Z] |
| Atmosphere, vessel | not printed (air presumably) |
| DTNB calibration | extinction coefficient and standard **not printed**; result "nmoles SH / 1 mL incub. med." |
| LOD / LOQ | not printed. [Z] The smallest bar (Th, 1.6 nmol/mL, SD ~1) suggests a working floor of ~1-2 nmol/mL |
| Mitochondrial work | 3 mmol/L TD, 0.25 mL mitochondria (6-7 mg protein/mL) in 0.1 M phosphate pH 7.4, 37 °C, 15 / 30 / 60 min, A250 (Fig. 1); SH in mitochondria with 2.5 mmol/L TD, corrected for non-enzymatic hydrolysis (Figs. 3-4); GSH 0.1 mmol/L or BSA 2 mg/mL with 0.1-0.5 mmol/L TD (Fig. 5); NAC / NFC deacylase in cytosol, pH 7.5, 30 °C (Table 3) |
| MDA | from cys 10 g + acetaldehyde 2 g in 6 M HCl, then pH 5.0 (synthesis); identified as product of CA in buffer by TLC (R_f 0.175, solvent IV) and 1H NMR of the bis-DNP derivative (Table 2: most signals doubled, equal intensity) |

## 3. Quantitative results

### Fig. 2 (p. 1923), DTNB-SH released, nmol/mL, 30 min, 37 °C, pH 7.4 [F]

Axis: 0 -> 100 nmol/mL, 8.97 px per unit at 300 dpi; bars read at their tops, SD caps as printed above them.

| ring | sugar / aldehyde | class | 1.5 mmol/L: SH (SD) | 2.5 mmol/L: SH (SD) | fraction released 1.5 / 2.5 [Z] | apparent k, /min, 1.5 / 2.5 [Z] |
|---|---|---|---|---|---|---|
| CGlu(D) | D-glucose | hexose | 56.7 (10) | 79.8 (8.4) | 3.8 / 3.2 % | 1.3e-3 / 1.1e-3 |
| CA | acetaldehyde | aliphatic | 50.3 (3.8) | 78.0 (5.4) | 3.4 / 3.1 % | 1.1e-3 / 1.1e-3 |
| **CXyl(L)** | **L-xylose** | **pentose** | **40.1 (5.2)** | **55.4 (5.9)** | **2.7 / 2.2 %** | **9.0e-4 / 7.5e-4** |
| **CAr(D)** | **D-arabinose** | **pentose** | **38.0 (4.4)** | **54.5 (5.9)** | **2.5 / 2.2 %** | **8.5e-4 / 7.3e-4** |
| CFuc(D) | D-fucose | 6-deoxyhexose | 25.9 (2.9) | 57.4 (4.2) | 1.7 / 2.3 % | 5.8e-4 / 7.7e-4 |
| CFuc(L) | L-fucose | 6-deoxyhexose | 22.5 (5.0) | 50.1 (4.9) | 1.5 / 2.0 % | 5.0e-4 / 6.7e-4 |
| **CRib(D)** | **D-ribose** | **pentose** | **16.3 (3.7)** | **39.9 (3.6)** | **1.1 / 1.6 %** | **3.6e-4 / 5.4e-4** |
| **CAr(L)** | **L-arabinose** | **pentose** | **12.3 (4.2)** | **27.5 (5.8)** | **0.8 / 1.1 %** | **2.7e-4 / 3.7e-4** |
| CGal(D) | D-galactose | hexose | 12.3 (2.2) | 26.6 (2.8) | 0.8 / 1.1 % | 2.7e-4 / 3.6e-4 |
| CIA | indole-3-acetaldehyde | aliphatic | 14.6 (5.4) | 26.8 (2.7) | 1.0 / 1.1 % | 3.3e-4 / 3.6e-4 |
| **CLyx(D)** | **D-lyxose** | **pentose** | **15.5 (4.5)** | **20.8 (5.2)** | **1.0 / 0.8 %** | **3.5e-4 / 2.8e-4** |
| CP | pyruvate | keto acid | 8.0 (2.1) | 14.7 (3.1) | 0.5 / 0.6 % | 1.8e-4 / 2.0e-4 |
| CF | formaldehyde | unsubstituted TC | 5.0 (1.8) | 9.9 (2.3) | 0.3 / 0.4 % | 1.1e-4 / 1.3e-4 |
| Th | (cysteamine + formaldehyde) | no carboxyl | 1.6 (1.0) | 3.1 (1.4) | 0.1 / 0.1 % | -- |

Apparent k = -ln(1 - f) / 30 min, treating the 30-min point as first-order release with no re-formation and the
DTNB reading as free cysteine; both assumptions are untested (sec. 1). Order as printed in the text: least stable
CGlu(D) and CA, most stable Th [Q].

Derived [Z]:
- **Concentration scaling.** Going from 1.5 to 2.5 mmol/L (x1.67), SH rises x1.38 (CXyl(L)), x1.43 (CAr(D)), x2.45
  (CRib(D)), x2.24 (CAr(L)), x1.34 (CLyx(D)). First-order release predicts x1.67, release to a dissociation
  equilibrium x~1.29 (sqrt). The pentose rings scatter on both sides within their SDs; the bar chart cannot choose.
- **The model's sugar is D-xylose; TTCA is the D-xylose + L-cysteine ring.** CXyl(L) here is the L-xylose + L-cysteine
  ring, a diastereomer of TTCA, not TTCA. The five pentose rings span x3.2 in apparent k (2.8e-4 to 9.0e-4 /min), so
  the sugar's configuration alone moves the rate by about half a decade at 37 °C, pH 7.4.

### Other results (not bearing on thiol supply) [M] / [Q]
- Fig. 1: mitochondrial A250 rise (thiazoline formation) for CF > CP, CIA > CA; none for Th or the sugar rings.
- Table 3: cytosolic deacylase, apparent K_m (mol/L) and V_max (umol/mL x min x mg protein): liver NAC 6.5 x 10^? (exponent
  lost in the scan) / 11.016; liver NFC 1.7 x 10^? / 0.062; brain NFC 2.4 x 10^? / 0.160 [!].
- Figs. 3-5: SH levels in mitochondria and in GSH / BSA solutions with TDs; CA and CGlu(D) lower SH at > ~0.6-0.8
  mmol/L, raise it below [Q].

## 4. Mechanism evidence and its kind

| claim | kind | strength |
|---|---|---|
| sugar rings release a thiol non-enzymatically at pH 7.4, 37 °C | DTNB thiol after 30 min, one time point | moderate for "some thiol appears"; weak for "free cysteine" |
| the open ring is the route (Scheme 2: ring opening -> Schiff base -> hydrolysis -> cys + aldehyde) | [P] scheme | none for the sugar rings |
| CA -> MDA (2 Cys : 1 acetaldehyde), equilibrating with CA | 1H NMR identity; TLC | strong for identity (acetaldehyde only); none for size |
| more MDA at lower pH | "preliminary experiments" [Q] | none |

## 5. What the paper can and cannot test in the model

| model part | can test? | why |
|---|---|---|
| R-rel's release pH law (n_r) | **no** | one pH (7.4) |
| `k_ttca_cys` level at process temperature | **no** | 37 °C only |
| K (C1, R-K) | **no** | single time point, equilibrium untested |
| P1a open-chain member | no | DTNB would read it, but cannot tell it from free cysteine |
| a bis-Cys pool for pentoses | no | MDA shown for acetaldehyde only; and 2 Cys : 1 aldehyde cannot make Zhang's 1 : 1 deficit (the case's own note) |
| a same-class level of pentose-ring release near neutral pH, 37 °C | yes, as a band-type datum | five pentose rings, 2.7e-4 - 9.0e-4 /min apparent, DTNB caveat |

## 6. Conflicts with existing records

1. **Karpinski 1990 (`karpinski1990_extraction.md`), the same L-arabinose ring (AThz = CAr(L)).** Karpinski: 4.8e-5 /min
   at pH 7.4, 25 °C (argon, voltammetric free cysteine, RSD up to 30 %). Wlodek: 2.7e-4 / 3.7e-4 /min at pH 7.4, 37 °C
   (air, DTNB) [Z]. Ratio 5.6-7.7 over 12 K implies an apparent Ea of ~110-130 kJ/mol [Z], against ~25-45 kJ/mol from
   Karpinski vs Belikov 1967 (`belikov1967_extraction.md`, 50 °C, sugar release). The formaldehyde ring tells the same
   story: CF here 1.1-1.3e-4 /min vs Karpinski's Thz 3.6e-5 /min at pH 7.4, 25 °C (x3-3.6). Likely sources: DTNB reads
   thiol that is not free cysteine (open chain, during the 1-min read), different atmosphere, different assay. Note
   that Radomski (a co-author here) and Karpinski were both at Warsaw University chemistry; whether the AThz samples
   share a batch is not stated. Not a conflict of kind; recorded for whoever sets a release band.
2. **Relative stability, pentose vs formaldehyde ring.** Karpinski: AThz / Thz = 1.0 (pH 6.0) and 1.3 (pH 7.4).
   Wlodek: CAr(L) / CF = 2.5 (1.5 mmol/L) and 2.8 (2.5 mmol/L) at pH 7.4 [Z]. Within a factor ~2, outside Karpinski's
   30 % RSDs.
3. **CF-076's release band** (k_r(100 °C, pH 4.5) in [2e-4, 8e-4] /min, declared; 0.9-3.2e-4 feasible). Wlodek's
   pentose rings at 37 °C, pH 7.4 already sit at 2.7e-4 - 9.0e-4 /min [Z]. Through Karpinski's pH shape (release
   peaks near pH 6, falls toward 7.4 and 4.5) this does not directly contradict a pH 4.5 value, and the DTNB caveat
   applies. Recorded, not resolved; not entered in the register (this branch does not edit `docs/reference/`).
4. **Queue description vs paper.** The download queue described "non-enzymatic hydrolysis of sugar TCAs; a bis-Cys
   thioaminal pool". The bis-Cys species (MDA) is for acetaldehyde only.

## 7. Flags

1. Fig. 2's x-axis shows **CLyx(D)**, which is not in Structure 1's list (6-13); its synthesis is not described.
2. The DTNB reading is total thiol 1 min after a 2 mmol/L DTNB addition. A ring that opens on that time scale (Radomski
   1989: relaxation ~0.2 /min at 37 °C for one sugar ring, `radomski1989_extraction.md`) could add to the reading
   during the minute; the low values (<= 3.8 %) say any such addition is small.
3. "0.1 MM phosphate" in the Fig. 3 legend and "0.1 mM phosphate" in Fig. 5 are typesetting variants of 0.1 M (the
   Methods print 0.1 M) [!].
4. Table 3 exponents are lost in the scan.
5. No sugar is measured; the "hydrolysis" is inferred from thiol alone.

## 8. Ledger-ready measurements (none proposed)

None is proposed. If an amendment to the M1c case later admits a same-class release level near neutral pH, the
candidates are Fig. 2's five pentose bars, with a second reader first and the DTNB caveat attached.

## M1c case criterion (frozen bd641ae8)

The frozen case (`docs/reference/cases/M1c_thiol_supply_structure.md` on `origin/m1c-case`, sec. 4.3) names wlodek1993
**decisive only if at >= 90 °C**, on: *a pentose TCA's release rate at >= 2 pH between 6 and 8, cysteine measured*;
otherwise *a band source for R-rel's release pH law*.

**Verdict: NOT MET.**
- **pH: one value, 7.4** (0.1 M phosphate label). The criterion needs >= 2 pH values in 6-8.
- **Temperature: 37 °C only.** The decisive role needs >= 90 °C.
- **Release rate: a single 30-min point**, not a rate. Five pentose rings (L-xylose, D-/L-arabinose, D-ribose,
  D-lyxose) release 0.8-2.7 % of their ring thiol in 30 min (12-55 nmol/mL from 1.5-2.5 mmol/L), apparent first-order
  2.7e-4 - 9.0e-4 /min [Z]. None is the D-xylose ring (TTCA).
- **Cysteine measured: partly.** Total DTNB-reactive thiol, not free cysteine; no sugar measured.
- **The fallback band role (R-rel's release pH law) is also not filled**: a pH law needs >= 2 pH values, and the paper
  has one.
- The bis-Cys pool the queue described is methyl-djenkolic acid from acetaldehyde, not from any sugar.
