# Davidek, Devaud, Robert & Blank 2006a — EXTRACTION
### *"Sugar Fragmentation in the Maillard Reaction Cascade: Isotope Labeling Studies on the Formation of Acetic Acid by a Hydrolytic β-Dicarbonyl Cleavage Mechanism"*
### Glucose + glycine (0.1 M each), 0.2 M phosphate, pH 4 / 6 / 8, 120 °C (and 90 °C, labelling only); 3-DG + glycine; **1-deoxyglucosone (1-DG) WITHOUT amine, pH 5-8**; 2,4- / 2,3-pentanedione + glycine pH 6 / 8 / 10; 13C-glucose isotopomers. Acetic acid by isotope dilution.

**Read 2026-10-31 from `data/articles/davidek2006.pdf`** (10 pp., JAFC 54:6667-6676; born-digital, clean text
layer; figures checked against 150 / 300 dpi page renders). The PDF file is named `davidek2006.pdf`; the measurement
ledger and the fit table call this paper **davidek2006a** (its companion, JAFC 54:6677, is `davidek2006b`), so this
dossier keeps that key and `data/keys/pdf_dossier_aliases.yml` maps `davidek2006 -> davidek2006a`. The keyed numbers
live in `data/lit/fit_rows/ph_acid_davidek_acids.yml` (ledger `davidek2006a_*` in `data/lit/measurement_ledger/trunk.yml`;
two-reader verification `data/lit/source_verification/ph_acid_davidek_acids.yml`); this dossier records what the paper
is and what it answers.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[Z]** derived by me, never printed ·
**[!]** flag.

## 1. What this paper answers for the model

Question (B7a, `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 11): the engine's DPO beta-cleavage to
acetic acid (`r_dpo_aa`) carries a base factor fitted on Davidek 2008's xylose + **glycine** pots; it drains DPO
hardest at high pH and cancels the enolisation split's pH law in the amine-free TTCA pots. The proposed next arm
(B7a-2) makes that cleavage AMINE-ASSISTED. Does the cleavage need, or run faster with, an amine?

| need | what Davidek 2006a gives | tag |
|---|---|---|
| (a) amine-free vs amine | **The one amine-free deoxyosone experiment in the Davidek series is here (Fig. 3):** synthetic 1-DG heated in phosphate buffer with **no glycine** (Methods: the deprotected 1-DG solution "was diluted five times with phosphate buffer"; no amine added) gives **0.80-0.84 mol acetic acid / mol 1-DG within 60 min at 120 °C, pH 6.5-8**. No paired 1-DG + glycine run, so amine *acceleration* is untested; amine *requirement* is excluded for the hexose analogue | [M] (end labels) / [F] |
| (a) other amine-bearing comparators | glucose + glycine 1.00 / 16.9 / 38.4 mol % at pH 4 / 6 / 8 after 240 min (Table 1); Amadori compound DFG "about 4, 28, and 56 mol %" at pH 4 / 6 / 8 after 4 h, attributed to ref 21 (Davidek 2005, Ann NY Acad Sci 1043:73, not on disk) | [M] / [C] |
| (b) pH, amine-free | 1-DG (Fig. 3), 120 °C: 60-min yield **0.20 (pH 5), 0.56 (pH 6), 0.80-0.84 (pH 6.5 / 7 / 8, indistinguishable)**; 5-min yield ~0.03 / ~0.28 / ~0.54-0.59. Steep between pH 5 and 6.5, **flat above 6.5** ("no significant differences ... at pH 6.5, pH 7, and pH 8") | [M] / [F] / [Q] |
| (b) pH, with glycine | 2,4-pentanedione + glycine (Fig. 5A, 120 °C): 240-min acetic acid ~52 / 89 / ~91 mol % at pH 6 / 8 / 10 (89 printed for pH 8); 20-min ~13 / ~40 / ~60 mol %, i.e. only ~3x per 2 pH units from 6 to 8 and < 2x from 8 to 10 | [M] / [F] |
| (c) labelling | 13C-glucose + glycine (Table 2): acetic acid from C-1 67.0-70.3 %, C-2 68.7-71.4 %, C-3 4.3-8.6 %, C-4 4.5-11.0 %, C-5 19.4-25.0 %, C-6 18.6-24.2 % (pairs C1/C2 ~70 %, C3/C4 ~8-10 %, C5/C6 ~20 %); "pH had only a negligible effect on the origin" | [M] |
| (d) conditions / method | 0.1 M precursors, 0.2 M Na phosphate (NaH2PO4 / Na2HPO4), pH set before heating, not held, not re-measured; 0.4-0.5 mL in 1.5 mL screw-cap vials, silicone bath; acetic acid by SIDA ([2H3]acetic acid, SPME-GC/MS or ether extract GC/MS); formic acid against [2H3]acetic acid (not SIDA) | [M] |

**Answer.** For the hexose 1-deoxyosone, the beta-cleavage to acetic acid **does not require an amine**: without any
amine, 1-DG gives ~0.8 mol/mol acetic acid in 60 min at pH >= 6.5 and 120 °C. Confidence that an amine-free 1-DG
pathway to acetic acid exists at these yields: **~90 %** (the Methods are explicit that 1-DG was diluted into
phosphate buffer only; residual risk: unstated carry-over of anything from synthesis / deprotection). Transfer to the
pentose DPO (1-deoxypentosone, the engine's species) is by analogy, **~75 %** (same 1-deoxy-2,3-diulose -> 2,4-diulose
-> beta-cleavage chemistry; Davidek 2008 finds 77-89 % of pentose acetic acid from C-1/C-2, the same signature).
Whether an amine *accelerates* the cleavage is **not answerable from this paper**: no matched precursor was run with
and without glycine. Implications for B7a-2:

1. **A purely amine-multiplied law (zero rate without amine) contradicts Fig. 3.** If B7a-2 is run, the
   amine-free term must stay non-zero; an *additive* amine term (k0 * f(pH) + k_N * [free amine] * f(pH)) is the
   form this paper allows. Confidence the pure multiplicative form is wrong for hexoses: ~85 %.
2. **The amine-free pH law is a saturating step, not a log-linear base factor.** 60-min yields 0.20 -> 0.56 ->
   0.80-0.84 from pH 5 -> 6 -> 6.5 and flat to pH 8 [M]. A base factor of the engine's kind, if fitted to
   glycine pots up to pH 8, would over-predict the pH 6.5 -> 8 rise in an amine-free pot (if its slope is steep
   there). Look at `r_dpo_aa`'s pH factor between 6.5 and 8 against this curve before adding an amine term.
3. **The cleavage is downstream of 1-DG, not 1-DG itself.** 1-DG is gone within 5 min at pH >= 6 (8 % left at
   pH 5 after 5 min, by HPLC [M]) while acetic acid keeps rising for 60 min (pH 6: 0.28 -> 0.56) [M/F]: the authors
   assign the slow tail to 1-deoxy-2,4-hexodiulose, acetylformoin and C4 fragments (Fig. 7). A one-step
   DPO -> AA reaction compresses this two-time-scale release.
4. **Glycine does not show up as a free-amine law in the beta-diketone test.** With glycine present, 2,4-pentanedione
   cleavage rises only ~3x from pH 6 to 8 and < 2x from 8 to 10 in the 20-min yields [F]; a term proportional to
   free glycine (pKa ~9.6) would rise ~100x from 6 to 8. Caveat: pH not held, and 0.2 M phosphate barely buffers
   at pH 10, so the pH 10 pot drifts down as ~0.09 M acid forms [Z]. Weak evidence (~60 %) against a
   free-amine-order term.
5. **Phosphate is not the discriminating variable here**: every Davidek pot and the Zhai 2021 TTCA pots alike sit
   in 0.2 M phosphate (`data/lit/fit_rows/route_submodel.yml`), so phosphate general-base catalysis cannot explain the
   Davidek / Zhai tension (Davidek 2006b's introduction says the beta-cleavage is "favored in the presence of
   phosphate buffer" [C], citing Hayami 1961 / Weenen 1998).

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Tomas Davidek, Stephanie Devaud, Fabien Robert (Nestle Research Center, Lausanne), Imre Blank* (Nestle Product Technology Center Orbe) |
| Title | as in the header |
| Venue | *J. Agric. Food Chem.* **2006**, 54 (18), 6667-6676 |
| DOI | **10.1021/jf060667q** (printed p. 6667; Crossref 2026-10-31: same title, authors, volume 54, issue 18, pp. 6667-6676) |
| Dates | received 9 Mar 2006; revised 1 Jul 2006; accepted 9 Jul 2006; published on web 9 Aug 2006 |
| Content | Tables 1 (acetic acid vs time / pH), 2 (13C shares), 3-5 (MS fragments of tetroses, tetrulose, quinoxaline); Figs. 1 (alpha-cleavage hypothesis), 2 (1-DG deprotection), 3 (1-DG -> acetic acid), 4 (beta-cleavage scheme), 5 (pentanediones), 6, 7 (scheme), 8. Supporting information (tetrose reference spectra) not on disk |

## 3. Methods [M]

| experiment | as printed |
|---|---|
| Glucose / glycine | glycine 0.35 mmol + D-glucose (or 13C isotopomer) 0.35 mmol in phosphate buffer (3.5 mL, 0.2 mol/L) -> 0.1 M each; 0.4 mL into 1.5 mL screw-cap vials; silicone bath **90 °C 15 h or 120 °C 4 h**; pH 4, 6, 8. Aliquot 0.3 mL + [2H3]acetic acid 250-600 µg, acidified to pH 3 (H3PO4), SPME-GC/MS. "Each sample was prepared at least in duplicate and injected twice" |
| Tetrose / C4H6O3 identification | glycine + glucose 0.5 mmol each in 5 mL buffer (0.2 M, pH 6 or 8) in a Pyrex tube; methoxime-TMS (tetroses) and o-phenylenediamine (quinoxaline) derivatives, GC-MS, not quantified |
| Pentanediones | glycine 0.5 mmol + 2,3- or 2,4-pentanedione 0.5 mmol in 5 mL buffer (0.2 M) -> 0.1 M each; 0.4 mL vials; 120 °C, 10-240 min; pH 6, 8, 10. SIDA for acetic acid; [2H5]propanoic acid (93 µg) and [13C3]acetone (170 µg) as standards |
| 3-DG (2) / glycine | glycine 0.05 mmol + 3-DG 0.05 mmol in 0.5 mL buffer (0.2 M) -> 0.1 M each; 120 °C 4 h; pH 4, 6, 8; SIDA |
| **1-DG (1), no amine** | protected 1-DG (4) 0.1 mmol in 2 mL water, Dowex 50W x8 H+ 2.5 h (complete deprotection by HPLC, Fig. 2B), resin filtered off, **"diluted five times with phosphate buffer (0.2 mol/L)"** -> nominal ~10 mM 1-DG in ~0.16 M phosphate [Z]; 0.5 mL in 1.5 mL vials, 120 °C up to 60 min; pH 5, 6, 6.5, 7, 8. Acetic acid: 0.4 mL + [2H3]acetic acid 61 µg, pH 3, ether extract, GC-MS (m/z 60 / 63). 1-DG: OPD quinoxaline, RP-HPLC 320 nm against 2-methylquinoxaline |
| Quantification | SIDA for acetic acid, acetone, propanoic acid; "[2H3]Acetic acid was used as internal standard for the quantification of formic acid" (not SIDA) [!]; natural 13C (1.10 %) corrected. CV < 10 % (Tables 1, 2, Fig. 3), < 5 % (Fig. 5), < 3 % (Fig. 2) |

[!] pH: set before heating, temperature of setting not stated, never re-measured; 0.2 M phosphate cannot hold pH 10
or pH 4 (outside pKa2 7.2 +- 1); acid production of up to ~0.09 M (pentanedione pH 8 / 10) is half the buffer.
[!] Methods give "90 °C, 15 h or 120 °C, 4 h", yet Table 1 prints a pH 4 cell at 480 min.

## 4. Results [M] / [F]

**Table 1 (glucose + glycine, 120 °C, acetic acid mol % of glucose; ND = not determined):**

| min | pH 4 | pH 6 | pH 8 |
|---|---|---|---|
| 15 | ND | ND | 21.4 +- 0.1 |
| 30 | ND | 9.7 +- 0.1 | 28.3 +- 0.2 |
| 60 | 0.28 +- 0.02 | 12.6 +- 0.1 | 33.0 +- 0.1 |
| 120 | 0.48 +- 0.04 | 15.1 +- 0.1 | 37.0 +- 0.2 |
| 240 | 1.00 +- 0.01 | 16.9 +- 0.2 | 38.4 +- 0.4 |
| 480 | 1.69 +- 0.17 | ND | ND |

**3-DG + glycine, 4 h, 120 °C (text p. 6670):** acetic acid 0.8 / 1.5 / 2.4 mol %, formic acid 0.9 / 1.5 / 1.7 mol %
at pH 4 / 6 / 8. Authors: hydrolytic alpha-cleavage of 3-DG "cannot be considered as a major pathway" to formic acid.

**Fig. 3, 1-DG without amine, 120 °C** (mol acetic acid / mol 1-DG; 60-min end labels printed; earlier points read
by me from a 300 dpi render, read error ~ +-0.02):

| min | pH 5 | pH 6 | pH 6.5 / 7 / 8 cluster |
|---|---|---|---|
| 0 | ~0.0 (one t = 0 marker cluster ~0.02) | | |
| 5 | ~0.03 | ~0.28 | ~0.54-0.59 |
| 10 | ~0.07 | ~0.41 | ~0.62-0.65 |
| 15 | ~0.085 | ~0.43 | ~0.70-0.76 |
| 30 | ~0.13 | ~0.48 | ~0.76-0.79 |
| 60 | **0.20** | **0.56** | **0.84 / 0.82 / 0.80** (labels; which marker is which pH is ambiguous [!]) |

1-DG decay: "Already after 5 min of heating, only 8% of 1 was found by RP-HPLC in the pH 5 sample, and no 1 could be
detected anymore at pH >= 6" (data not shown) [M/Q].

**Fig. 5 (with glycine, 120 °C, mol % of precursor):** 2,4-pentanedione -> acetic acid ~ acetone (near-equimolar,
"constantly increased"); 240 min: ~52 (pH 6), **89 (pH 8, printed)**, ~91 (pH 10); 2,3-pentanedione -> acetic acid
**11 % (pH 8, 4 h, printed)**, ~9-13 % at all pH, with propanoic acid at similar levels; no propanal or acetaldehyde
detected. Authors assign 2,3-pentanedione's acids to an oxidative alpha-cleavage (detailed in davidek2006b).

**Table 2 (13C shares of singly labelled acetic acid, %, glycine present):**

| label | pH 6, 90 °C 15 h | pH 6, 120 °C 4 h | pH 8, 90 °C 15 h | pH 8, 120 °C 4 h |
|---|---|---|---|---|
| [1-13C] | 68.9 | 70.3 | 69.8 | 67.0 |
| [2-13C] | 70.3 | 70.5 | 71.4 | 68.7 |
| [3-13C] | 4.3 | 7.7 | 8.1 | 8.6 |
| [4-13C] | 4.5 | 7.4 | 7.3 | 11.0 |
| [5-13C] | 25.0 | 19.4 | 21.9 | 22.1 |
| [6-13C] | 24.2 | 18.6 | 22.4 | 19.6 |

Mechanism (Fig. 7) [Q]: 1-DG isomerises to 1-deoxy-2,4-hexodiulose (5); OH- attack at C-2 releases acetic acid
(C-1/C-2) + a C4 enediol (-> tetrulose, erythrose, threose: identified for the first time, tetrulose 5-10x the
aldotetroses by peak height); acetylformoin and 2-hydroxy-3-oxobutanal (= 1-deoxy-2,3-tetrodiulose, OPD quinoxaline)
give acetic acid from C-5/C-6 and C-3/C-4. Glyceric acid (attack at C-4) is "significantly" lower than acetic acid
(data in davidek2006b). "Does not require the intermediacy of any oxidizing species" [V].

## 5. Flags

- [!] **No amine-free vs amine pair on the same precursor.** 1-DG was run only without amine, pentanediones only
  with glycine, glucose / 3-DG only with glycine.
- [!] 1-DG concentration is nominal (deprotection yield, resin losses not printed); the readers' adjudication
  (2026-10-30) left Fig. 3 un-keyed for this reason and for the overlapping pH 6.5 / 7 / 8 markers.
- [!] The DFG yields (4 / 28 / 56 mol %) sit in a sentence citing ref 21 (Davidek et al. 2005, Ann NY Acad Sci
  1043:73-79, 10.1196/annals.1333.009): their method is in that paper, which is not on disk (download candidate:
  title "The Effect of Reaction Conditions on the Origin and Yields of Acetic Acid Generated by the Maillard Reaction").
- [!] Abstract gives "90-120 °C, pH 6-8"; the experiments span pH 4-10.
- Fit table spot-check (2026-10-31): all 13 Table 1 rows and all 6 3-DG rows in
  `data/lit/fit_rows/ph_acid_davidek_acids.yml` match the PDF (values, sd, conditions, volumes); the row text
  "1.0 +- 0.01" drops the printed trailing zero of "1.00 +- 0.01" (same value).
