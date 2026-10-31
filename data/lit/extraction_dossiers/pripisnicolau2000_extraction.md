# Pripis-Nicolau, de Revel, Bertrand & Maujean 2000 — EXTRACTION
### *"Formation of Flavor Components by the Reaction of Amino Acid and Carbonyl Compounds in Mild Conditions"*
### Wine conditions: 14 amino acids x 6 carbonyls (glyoxal, methylglyoxal, diacetyl, pentane-2,3-dione, acetoin, acetol), 20 mmol/L each, 12 % hydro-alcoholic, pH 3.5 (H3PO4 / NaOH), 25 °C in the dark, 4 weeks. **Strecker aldehydes are IDENTIFIED, never quantified.** The only numbers are cysteine loss (pH 3.5 vs 8, three carbonyls, 25 °C) and 2-acetyl-2-thiazoline / 2-acetylthiazole from cysteine + MGO (2 mM, 10 / 25 °C). **Not a second laboratory for donor-generic Strecker magnitudes.**

**Read 2026-10-31 from `data/articles/pripis-nicolau2000.pdf`** (6 pp., JAFC 48:3761-3766; born-digital text layer,
clean). All six pages read; Figs. 4 and 6 checked against 150 / 400 dpi renders (scratchpad `pn-*.png`,
`fig4crop.png`, `fig6crop.png`) and read by me: read error ~ +-2 percentage points on Fig. 4, ~ +-10 mg/L
(cysteine) and ~ +-0.02 mg/L (products) on Fig. 6, time ~ +-0.05 log10 units. Figs. 1-3 are schemes, Fig. 5 is
two GC-FID chromatograms (no scale usable as amounts).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

The question it was fetched for (Strecker forced-input prereg sec. 3 "Cross-laboratory hold-outs"; evidence map
2.6; `docs/guides/MECHANISM_REVIEWS_2026-10-31.md` row: "six dicarbonyls x wine amino acids at pH 3.5: donor-generic
Strecker"): is this a usable SECOND LABORATORY for donor identity factors, yields, time courses and pH of the
Strecker step, and does it agree with Hofmann 2000 / 2000b?

| model need | what Pripis-Nicolau 2000 gives | tag |
|---|---|---|
| Strecker aldehyde YIELD per donor (PA, 3-MB, 2-MB, 2-MP, methional) | **none.** Table 1 lists products by GC-MS identity and odour intensity only (*, **, ***); no amount, no area, no internal-standard ratio for any aldehyde | [M] (absence) |
| donor ranking for the aldehyde | **none.** Table 1 gives the SAME products and the same odour score ("**" moderate) for Val / Leu / Ile / Phe with all four dicarbonyls, listed as one block | [M] |
| which donors | glyoxal, MGO, diacetyl, pentane-2,3-dione (α-dicarbonyls) + acetoin, acetol (α-hydroxyketones). **No 3-DG, no glucosone, no sugar.** Overlap with Hofmann 2000 Table 4: MGO and GO only | [M] |
| amino-acid identity factor | **none quantitative.** PA (and a "low level" of benzaldehyde) from Phe; 2-methylpropanal / 3-methylbutanal / 2-methylbutanal from Val / Leu / Ile; methional, MeSH, DMDS from Met | [Q] |
| Strecker time course | **none** for any aldehyde; extracts analysed at 7, 14, 30 days | [M] |
| pH law | only for CYSTEINE loss (pH 3.5 vs 8; Fig. 4) and for 2-acetylthiazole / thiazoline (more at pH 3.5; Fig. 5, chromatograms) | [F] / [Q] |
| temperature | 25 °C vs 10 °C, cysteine + MGO 2 mM only (Fig. 6) | [F] |
| amino-acid consumption | HPLC (OPA / IDA derivatisation) was run, but only cysteine is shown | [M] |
| acid (PAA) | not analysed | [M] |

**Answer: NO, it is not a usable second laboratory** for the forced-input Strecker layer (confidence ~0.97). No
quantity in it can be set beside Hofmann 2000 / 2000b, so no dex agreement or disagreement can be computed. What
it adds is an existence statement: dicarbonyls (and, after 4 weeks, even the hydroxyketones, presumably via their
redox partners) give Strecker aldehydes at **25 °C, pH 3.5, 20 mM, in 12 % ethanol** [Q]. The repo's own
`docs/guides/CITATION_CONTEXT_2026-10.md` already records a citer's note that the aldehyde concentrations were not
reported; this dossier confirms it from the PDF.

What differs from Hofmann 2000, for the record:

| condition | Pripis-Nicolau 2000 | Hofmann 2000 (Table 4 / Table 3) |
|---|---|---|
| temperature | 25 °C (10 °C in Fig. 6b) | 98-100 °C |
| time | 7-30 days (aldehydes); 0-7 d (Cys, Fig. 4); 0-~14 d (Fig. 6) | 10-500 min |
| pH / buffer | 3.5, H3PO4 (1/3) + 1 M NaOH, strength not given; Fig. 4 also pH 8 | 7.0, phosphate 0.1-0.5 mol/L |
| solvent | 12 % hydro-alcoholic (ethanol by context of wine; the solvent is named only "hydro alcoholic") | water |
| concentrations | 20 mmol/L each (stoichiometric); 2 mM in Fig. 6 | 100 mmol/L each |
| vessel / O2 | "stored at 25 °C in the dark", vessel and headspace not stated | closed vial, air vs argon |
| analysis | solvent extraction (ethyl acetate, ether, ether/hexane), GC-FID / FPD / NPD / MS; SIM quantification only for 2-acetylthiazole and 2-acetyl-2-thiazoline (IS octan-3-ol) | SIDA GC-MS |

## 2. Conditions [M]

- **Model mixtures:** carbonyl + amino acid "in stoichiometric conditions (20 mmol/L)" in a 12 % hydro-alcoholic
  solution, adjusted to pH 3.5 with H3PO4 (1/3) and 1 M NaOH; stored at 25 °C in the dark; flavour followed over a
  4-week storage period. 14 amino acids (Arg, Phe, Lys, Cys, Met, Pro, Val, Leu, Ile, Gly, Ser, Ala, GABA, Glu).
- **Extraction:** 50 mL of mixture extracted with ethyl acetate (2 x 5 mL), ether (2 x 5 mL) or ether / hexane 50/50
  (5 + 2 x 2 mL) after 7, 14 and 30 days.
- **GC:** FPD on Carbowax 20M; NPD / FID on HP-5; MS (HP 5972, EI 70 eV) on BP-21. 2-Acetylthiazole (m/z 99, 127)
  and 2-acetyl-2-thiazoline (m/z 101, 129) quantified in SIM against octan-3-ol (m/z 58, 83). No calibration or
  response factor is printed.
- **Amino acids:** RP-HPLC with OPA / 2-sulfanylethanol derivatisation and a second derivatisation with iodoacetic
  acid for sulfur amino acids; fluorescence 356 / 445 nm.
- **Kinetic study (Fig. 4):** cysteine + MGO, diacetyl or acetoin at pH 3.5 and pH 8, 25 °C. Concentration not
  restated; the general recipe is 20 mmol/L [Z: inferred, the text contrasts Fig. 6 as "a low concentration (2 mM)"].
- **Fig. 6:** cysteine + MGO "at a low concentration (2 mM)", 10 and 25 °C, pH 3.5 (by context).
- **Replicates:** not stated; no error bars.

[!] Fig. 6's cysteine axis starts at ~355 mg/L [F], which is 2.93 mM as free cysteine (121.16 g/mol) [Z], not 2 mM.
It equals 2.02 mM if the mass is of cysteine hydrochloride monohydrate (175.6 g/mol) [Z]. The paper does not say
which salt was used; treat the mg/L axis as relative (fraction of initial) when keying.

## 3. Quantitative results

**Fig. 4: cysteine remaining (% of initial) vs time (log axis), 25 °C** [F] (times read on the log axis; the
first plotted time is ~56 min, i.e. ~1 h; later points ~3, ~5, ~20-30, ~48, ~72-100 h and ~7 d):

| carbonyl, pH | ~1 h | ~3 h | ~5 h | ~20-30 h | ~48 h | ~7 d |
|---|---:|---:|---:|---:|---:|---:|
| acetoin, 3.5 | 99 | 98 | — | 98-97 | 96 | 87 |
| acetoin, 8 | 86 | 54 | — | 38-29 | 27 | 27 |
| MGO, 3.5 | 14 | — | 14 | 11 | 10.5 | 10.5 |
| MGO, 8 | ~3 | ~3 | ~3 | ~3 | ~3 | ~3 |
| diacetyl, 3.5 | 53 | — | 44 | 43-39 | 34 | 33 |
| diacetyl, 8 | ~6 | ~6 | ~6 | ~5 | ~5 | ~5 |

Text [Q]: decomposition of cysteine faster at pH 8 than at 3.5, with MGO and diacetyl; slow with acetoin,
"particularly at pH 3.5". Plateaus, not completion: 10-14 % of cysteine remains with MGO at pH 3.5 for a week [F].

**Fig. 6: cysteine + MGO, 2 mM, pH 3.5** [F] (left axis cysteine mg/L; right axis products mg/L):

| t | Cys 25 °C | thiazoline 25 °C | 2-AT 25 °C | Cys 10 °C | thiazoline 10 °C | 2-AT 10 °C |
|---|---:|---:|---:|---:|---:|---:|
| 0 | ~355 | 0 | 0 | ~355 | 0 | 0 |
| 10 min | ~43 | 0.42 | 0.10 | ~120 | — | — |
| 20-30 min | ~37-40 | 0.49 | 0.10 | ~90 | 0.36 / 0.35 | 0.08 / 0.07 |
| 35-50 min | ~38 | 0.58 | 0.13-0.14 | ~82 | 0.36 | 0.08 |
| ~60 min | ~40 | **0.65 (max)** | 0.17 | ~78 | — | — |
| ~120-180 min | ~40 | — | — | ~72 / 57 | 0.65 / **0.67 (max, ~180 min)** | 0.13 |
| ~2500-5000 min | ~40 | 0.10-0.11 | 0.33-0.40 | ~52 | 0.32-0.29 | 0.12-0.15 |
| ~20000 min (~14 d) | ~30 | 0.11 | **0.58** | ~47 | 0.22 | 0.39 |

Text [Q]: the thiazoline shows an "intermediate" curve (max 60 min at 25 °C; 180 min at 10 °C, "3 h later") and is
"oxidized in 2-acetylthiazole which continued to increase regularly" [P: the oxidation is inferred from the shapes].

**Derived** [Z]:
- 2-Acetylthiazole at ~14 d, 25 °C: 0.58 mg/L = 4.6 µM = **0.23 % of a 2 mM cysteine charge**; thiazoline peak 0.65
  mg/L = 5.0 µM = 0.25 %. So ~87 % of the cysteine is gone in 10 min while < 0.5 % appears as the two measured
  heterocycles: the fast cysteine sink is NOT these products (the text names 2-acetylthiazolidine, the
  cysteine-MGO thiazolidine, as present but unquantified).
- Second-order lower bound on Cys + MGO consumption at 25 °C, pH 3.5 (equimolar, first point only):
  **>= 0.34 L mmol^-1 min^-1** from Fig. 6a (2 mM, 10 min); the 20 mM Fig. 4 point (~1 h) gives only >= 0.0055
  (the reaction is already at its plateau by then). Lower bounds, not rates.
- The plateau is not a single 1:1 equilibrium: an apparent K = [adduct]/[Cys][MGO] of ~35 L/mmol at 2 mM vs ~4 L/mmol
  at 20 mM (if every lost cysteine were a 1:1 adduct). Candidate causes: other sinks (H2S, acetaldehyde route, Fig.
  2b), MGO hydrate / oligomer speciation, or OPA-HPLC partially reading the adduct as cysteine.

## 4. Mechanism evidence and its kind

Laboratory: **Faculté d'Œnologie, Bordeaux (INRA / Université Victor Segalen; de Revel group) with Reims.**
- **Product identification only** (GC-MS against purchased standards, retention / spectra): Strecker aldehydes,
  H2S, MeSH, CS2, DMDS, methional, pyrazines (pyrazine, 2,5- / 2,6-DMP, TMP, diethyl-dimethylpyrazines),
  2-methylthiazole, 2-acetylthiazole, 2-acetyl-2-thiazoline, 2-acetylthiazolidine, trimethyloxazole and two
  dimethylethyloxazoles [Q].
- **Within-study factor (kind 6):** pH 3.5 vs 8 on cysteine loss (Fig. 4) and on acetylthiazole / thiazoline
  (Fig. 5: more at pH 3.5, chromatogram heights only); 10 vs 25 °C (Fig. 6).
- **Proposed, not shown [P]:** cysteine decarboxylation to H2S + acetaldehyde + NH3 "in the presence of glyoxal"
  (Fig. 2b); oxazoles from α-aminoketone + acetaldehyde (Fig. 3); methional -> MeSH / DMDS by retro-Michael (Fig. 1);
  thiazoline -> thiazole oxidation.
- No labelling, no trapping, no mass balance, no fed intermediate beyond the carbonyl itself.

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| Forced-input Strecker layer is SINGLE-LAB (prereg secs. 3, 10.9; EX-049 / EX-051) | **LEAVES UNCHANGED** | no Strecker aldehyde amount; no 3-DG / glucosone; 25 °C, pH 3.5. Cannot be a cross-lab check of k_D, s_D, phi_D or any signature (T, O, G, I groups) |
| Donor-identity ranking (evidence map 2.6 "CAN (ranking)") | **LEAVES UNCHANGED** | no ranking: one odour score for all four dicarbonyls |
| CF-073 (Weenen vs Hofmann, x14-20) | **LEAVES UNCHANGED** | no yield |
| CF-074 / CF-075 | **LEAVES UNCHANGED** | Hofmann-internal |
| Identity factors f_Leu, f_Ala (prereg I1 / I2) | **LEAVES UNCHANGED** | identification only |
| Sulfur lane: Cys + MGO chemistry (`r_cys_actz`, `r_mgo_cys_strecker`, `MGO_CYSTEINE_STRECKER` fitted, no measured constant) | **CONTEXT, not a constant** | Cys + MGO is fast at 25 °C (>= 0.34 L/(mmol min) lower bound, pH 3.5) but the fast sink is an adduct pool (thiazolidine) with < 0.5 % to the thiazoline / thiazole; any cysteine-loss rate taken from this paper is NOT a Strecker rate. pH 8 faster than 3.5 for Cys loss; the acetylthiazole / thiazoline the reverse [F / Q]. 25 °C only: no transfer to 100 °C without an Ea |
| Cysteine fate / reversible carbonyl adducts (CF-076 context) | **LEAVES UNCHANGED** (analogue) | consistent with fast, partly reversible Cys-carbonyl adducts at room temperature (cf. Zhai 2019's TTCA before heating), but MGO is not a sugar; nothing transfers to K or k_f of TTCA |

## 6. Disagreements found

1. **`docs/guides/MECHANISM_REVIEWS_2026-10-31.md` (literature row): "six dicarbonyls x wine amino acids at pH 3.5:
   donor-generic Strecker"** vs the PDF: four α-dicarbonyls + two α-hydroxyketones, and NO quantified Strecker
   output. Cause: the row was written from the abstract (it says so). Not a usable donor-generic dataset.
2. **`data/lit/extraction_dossiers/jousse2002_extraction.md` (R8): Jousse 2002 calibrated a prefactor to
   Pripis-Nicolau's "cysteine + methylglyoxal fully consumed < 90 min at 25 °C"** vs Figs. 4 / 6: ~86-87 % consumed
   by the first point (10 min at 2 mM; ~1 h at 20 mM) and a 10-14 % cysteine plateau for a week [F]; and the
   consumed cysteine is not a Strecker product (< 0.5 % to the measured heterocycles). Any constant calibrated as
   "Strecker of cysteine" on this curve measures adduct formation. Cause: chemistry (adduct vs Strecker), plus a
   "fully" that the figure does not show.
3. **Internal [!]:** Fig. 6 initial cysteine ~355 mg/L vs "2 mM" (2.93 mM as free base); consistent with 2 mM only as
   the HCl monohydrate (sec. 2).

## 7. Ledger-ready measurements

- **Fig. 6a / 6b** (Cys, 2-acetyl-2-thiazoline, 2-acetylthiazole; 2 mM Cys + MGO, pH 3.5, 25 / 10 °C, 10 min-14 d):
  the only time-resolved Cys + MGO -> acetylthiazoline -> acetylthiazole series at low temperature; useful as a
  shape check of the sulfur lane's `r_cys_actz` topology (thiazoline intermediate) if the lane ever runs at 25 °C.
  Key the mg/L axes, flag the 355 mg/L start.
- **Fig. 4** (Cys % remaining, MGO / diacetyl / acetoin x pH 3.5 / 8, 25 °C): a within-study pH factor for cysteine
  + carbonyl adduct formation; diagnostic only (semi-quantitative, 20 mM inferred).
- Nothing for the Strecker lane: no aldehyde was quantified.
