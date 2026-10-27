# Guth, Hofmann, Schieberle & Grosch 1995 — EXTRACTION
### *"Model Reactions on the Stability of Disulfides in Heated Foods"*
### A disulfide paper, not a thiol paper. At 100 °C for 2 h in water, MFT-MFT and FFT-FFT split back into their thiols (86 % and ≥ 98 % of the disulfide gone when a thiol trap is present, 62 % and 92 % of the theoretical two thiols recovered). Without the trap the MFT-MFT pool is conserved (95.5 %) and only scrambles. Nothing in it measures the loss of a thiol. One temperature, no pH, no time series. Read 2026-10-28.

**Full read of `data/articles/guth1995.pdf`** (5 pages, clean text layer). Tables 1–4 were checked
on a 200-dpi raster of p. 2200; the text layer renders "h/" as "w" in Table 4's time column
("12 w100 °C" is **12 h/100 °C**) and drops superscripts from the deuterium labels. Figures 1–7
are structures, two mass spectra, one mass chromatogram and one EPR trace; no kinetic figure.

**Provenance:** **[M]** measured, printed as a number · **[Q]** qualitative · **[Z]** derived by me
from the paper's own numbers, never printed · **[C]** outside the paper (general chemistry or
another dossier) · **[!]** flag.

---

## 0. IDENTITY — matches (with an author-list correction)

| field | value as printed |
|---|---|
| Authors | **Helmut Guth, Thomas Hofmann, Peter Schieberle, Werner Grosch*** (four authors; the brief's "Guth & Grosch" is short by two) |
| Title | "Model Reactions on the Stability of Disulfides in Heated Foods" |
| Venue | *J. Agric. Food Chem.* **1995**, 43 (8), **2199–2203** |
| DOI | **10.1021/jf00056a044**. Not printed in the PDF (only the code JF950138E). Crossref REST API (`api.crossref.org/works/10.1021/jf00056a044`, 2026-10-28) returns this title, JAFC 43(8):2199–2203, Aug 1995, authors Guth, Hofmann, Schieberle, Grosch. Matches |
| Affiliation | Deutsche Forschungsanstalt für Lebensmittelchemie and Institut für Lebensmittelchemie der TU München, Garching |
| Dates | received 7 Mar 1995, accepted 1 Jun 1995; ASAP 1 Jul 1995 |
| PDF character | 5 pages: Table 1 (MS calibration), Tables 2–4 (model experiments, one end point each), Table 5 (MS of adducts), Figs. 1–7 (no kinetics). No replicates or SD stated anywhere [!] |

**[!] The thiol-oxidation companion is a different paper.** The thiol → disulfide stability study
cited here as "Hofmann et al., 1995, submitted" appeared as Hofmann, Schieberle & Grosch 1996,
JAFC 44:251–255 (`hofmann1996_extraction.md`, already in the corpus; diethyl ether, 6 °C). The
pairing "stability of thiols / their disulfides" in the brief spans both papers; **this one holds
only the disulfide half.**

---

## 1. THE ONE-PARAGRAPH ANSWER

Sealed vials, **100 °C, 2 h** (Table 4: 2 or 12 h), disulfides at **µg scale** in **50 mL water**
(unbuffered; **pH never stated**) or **1–5 mL benzene**. (i) **Scrambling without a trap**
(Table 2): a 1:1 mix of MFT-MFT and [²H₆]MFT-MFT gives **21.5 % mixed disulfide in water and
12.5 % in benzene** after 2 h, with the total disulfide mass **conserved at 95.5 %** in both. The
FFT-FFT pair scrambles too, but in water **only 45 % of the disulfide mass survives** (no fate
given). (ii) **With 4-vinylpyridine as thiol trap** (Table 3): in water **FFT-FFT ≥ 98 %, MFT-MFT
86 % and cystine ≥ 97.6 % disappear**, giving **FFT 52 µg (92 % of two thiols per disulfide), MFT
112 µg (62 %) and cysteine 18 %**. In benzene nothing happens (≤ 12 % loss, thiol < 1 µg).
(iii) D₂O labelling puts one D in the trapped MFT (yield 63 %): the authors propose **hydrolysis
RS–SR + H₂O → RSH + RSOH, with sulfenic-acid disproportionation** supplying the >50 % excess.
(iv) In benzene with a **hydrogen donor** (1,4-cyclohexadiene, 1,4-hexadiene, BHT) thiyl radicals
abstract H and **regenerate the thiols** (FFT 50–58 µg from 522 µg FFT-FFT in 12 h; MFT 30 µg from
210 µg MFT-MFT in 2 h), and with cyclohexadiene **22 % (FFT) and 38 % (MFT)** of the sulfur leaves
as non-thiol, non-disulfide products, identified for FFT as four **1-(2-furfurylthio)cyclohexadiene**
isomers. Conclusion printed: "Disulfides are not stable during heating of foods. They are mainly
hydrolyzed with formation of thiols." **No thiol was ever the starting material, and no step has
a time course, a second temperature, a pH or an oxygen variable.**

---

## 2. METHODS [M] (pp. 2199–2200)

| item | as printed |
|---|---|
| Thiols / disulfides | MFT, FFT commercial (Aldrich); MFT-MFT and FFT-FFT synthesised (Hofmann et al. 1995 procedure); [²H₆]MFT-MFT, [²H₄]FFT-FFT, [²H₃]MFT, [²H₂]FFT labelled standards |
| Loading | disulfide in diethyl ether pipetted into a vial, ether blown off under N₂, residue taken up in water or benzene (+ reactants), **vial sealed and heated**. Headspace gas and its oxygen **not stated** [!] |
| Water volume | **50 mL** (Tables 2–3) [Z: 100 µg MFT-MFT = 0.44 µmol → **8.8 µM** of each disulfide; 180 µg → 15.9 µM] |
| Benzene volume | 1 mL (Tables 2–3), 0.5–5 mL (Table 4) |
| Temperature / time | **100 °C, 2 h** (Tables 2, 3; Table 4 nos. 3–4); **12 h** (Table 4 nos. 1–2). Come-up not stated |
| pH | **not stated; unbuffered water** [!]. 4-Vinylpyridine (pyridine-type base, pKa ≈ 5.6 [C, general]) at 500 µg/50 mL = 95 µM will push unbuffered water to pH ≈ 7–8 [Z, order of magnitude], so the trapped and untrapped water runs need not share a pH |
| Thiol trap | 4-vinylpyridine 100 µg (benzene) / **500 µg** (water); excess removed afterwards with L-cysteine ("cf. Experimental Procedures", where the step is **not described** [!]) |
| Work-up (water) | internal standards added after reaction; ether extraction 3 × 15 mL, Na₂SO₄, concentrated to 1 mL, HRGC-MS (1 µL) |
| Quantitation | isotope dilution, HRGC-MS (FFAP, CP-Wax 51); calibration factors Table 1; cystine/cysteine by amino-acid analysis as the pyridylethyl derivative |
| EPR | 3.1 µmol disulfide + 10 µmol DMPO in 3 mL benzene or water, sealed 10-mL vial, 100 °C, **1 h**, ice-quenched; Bruker ESP 300, 9.75 GHz |

---

## 3. RESULTS — tables as printed [M] (p. 2200, checked on the raster)

### Table 1. Thin-film capillaries, selected ions, calibration factors

| analyte | capillary | ion (m/z) | IS | IS ion (m/z) | calib factor |
|---|---|---|---|---|---|
| 2-furfurylthiol (FFT) | FFAP | 81 | [²H₂]FFT | 83 | 0.83 |
| bis(2-furfuryl) disulfide (FFT-FFT) | FFAP | 227 | [²H₄]FFT-FFT | 231 | 0.83 |
| [β-(4-pyridyl)ethyl] 2-furfuryl thioether (PEFTE) | FFAP, CP-Wax 51 | 220 | [²H₂]PEFTE | 222 | 0.83 |
| 2-methyl-3-furanthiol (MFT) | FFAP | 115 | [²H₃]MFT | 118 | 0.98 |
| bis(2-methyl-3-furyl) disulfide (MFT-MFT) | FFAP | 227 | [²H₆]MFT-MFT | 233 | 0.98 |
| [β-(4-pyridyl)ethyl] 2-methyl-3-furyl thioether (PEMFTE) | FFAP, CP-Wax 51 | 220 | [²H₃]PEMFTE | 223 | 0.98 |

### Table 2. Thermolysis of mixtures of disulfides and their labelled analogues (after 2 h)

| no. | reaction system disulfides (µg) | starting disulfides (µg) | product disulfide (µg) |
|---|---|---|---|
| | **I. In benzene (1 mL) at 100 °C** | | |
| 1 | MFT-MFT (100) plus [²H₆]MFT-MFT (100) | MFT-MFT (83) plus [²H₆]MFT-MFT (83) | [²H₃]MFT-MFT (25) |
| 2 | FFT-FFT (183) plus [²H₄]FFT-FFT (181) | FFT-FFT (121) plus [²H₄]FFT-FFT (122) | [²H₂]FFT-FFT (73) |
| | **II. In water (50 mL) at 100 °C** | | |
| 3 | MFT-MFT (100) plus [²H₆]MFT-MFT (100) | MFT-MFT (73) plus [²H₆]MFT-MFT (75) | [²H₃]MFT-MFT (43) |
| 4 | FFT-FFT (100) plus [²H₄]FFT-FFT (100) | FFT-FFT (56) plus [²H₄]FFT-FFT (23) | [²H₂]FFT-FFT (11) |

Footnote: "At the end of the reaction the disulfides were quantified by mass chromatography."
Text: mixed MFT disulfide "formed in yields of 12.5 % (benzene) and 21.5 % (water)" [M]
(= product ÷ total initial disulfide, [Z] 25/200 and 43/200).

### Table 3. Formation of thiols by thermolysis of disulfides with 4-vinylpyridine (after 2 h)

| no. | disulfide (µg) | 4-vinylpyridine (µg) | starting disulfide (µg) | thiol (µg) |
|---|---|---|---|---|
| | **I. In benzene (1 mL) at 100 °C** | | | |
| 1 | FFT-FFT (150) | 100 | FFT-FFT (132) | FFT (< 1) |
| 2 | MFT-MFT (150) | 100 | MFT-MFT (145) | MFT (< 1) |
| | **II. In water (50 mL) at 100 °C** | | | |
| 3 | FFT-FFT (56) | 500 | FFT-FFT (< 1) | FFT (52) |
| 4 | MFT-MFT (180) | 500 | MFT-MFT (25) | MFT (112) |
| 5 | cystine (4.1) | 500 | cystine (< 0.1) | cysteine (0.74) |

Footnotes: thiols measured as their pyridylethyl thioethers (PEFTE, PEMFTE) against
[²H₄]PEFTE / [²H₃]PEMFTE; cystine/cysteine by amino-acid analysis. Text: disulfides "degraded by
100 % (FFT-FFT) or 86 % (MFT-MFT)"; cysteine "in a yield of 18 %". D₂O repeat of no. 4:
[²H₁]PEMFTE "the major product (yield 63 %)".

### Table 4. Thermolysis of disulfides in benzene containing a hydrogen donor

| no. | disulfide (µg) | hydrogen donor | vol (mL) benzene | time/temp | disulfide at end (µg) | thiol at end (µg) |
|---|---|---|---|---|---|---|
| 1 | FFT-FFT (522) | 1,4-cyclohexadiene (0.5 mL) | 0.5 | 12 h/100 °C | FFT-FFT (359) | FFT (50) |
| 2 | FFT-FFT (522) | BHT (1 mg) | 1 | 12 h/100 °C | FFT-FFT (312) | FFT (58) |
| 3 | MFT-MFT (210) | 1,4-cyclohexadiene (0.1 mL) | 5 | 2 h/100 °C | MFT-MFT (100) | MFT (30) |
| 4 | MFT-MFT (210) | 1,4-hexadiene (0.1 mL) | 5 | 2 h/100 °C | MFT-MFT (180) | MFT (30) |

Text: for no. 1 the recovered FFT-FFT + FFT "amounted only to 78 %" of the initial; for no. 3 "the
loss (38 %) was even greater". The no. 1 deficit is attributed to thiyl + cyclohexadienyl radical
coupling, four 1-(2-furfurylthio)cyclohexadiene isomers (M⁺ 192, Table 5, Fig. 5). Not quantified.

### EPR [Q]
DMPO alone: no signal. MFT-MFT + DMPO in benzene, 100 °C, 1 h: six weak unassigned lines. In
water: the **DMPO/OH adduct** (a_N = a_H = 14.9 G), read as thiyl-DMPO adducts hydrolysed during
heating (hypothesis, Fig. 7). So both an ionic (hydrolysis) and a radical (homolysis) route are
proposed in water.

---

## 4. DERIVED QUANTITIES [Z] (mine, not printed)

MW: MFT = FFT = 114.16; MFT-MFT = FFT-FFT = 226.30; cystine 240.30; cysteine 121.16.

**Mass balance and thiol stoichiometry, water + trap, 100 °C, 2 h (Table 3):**

| disulfide | initial (µmol) | left (µmol) | disappeared | thiol found (µmol) | thiol / (2 × initial) | thiol / disulfide lost | apparent k_loss (min⁻¹) |
|---|---|---|---|---|---|---|---|
| FFT-FFT | 0.247 | < 0.0044 | ≥ 98.2 % | 0.456 | 92 % | ≥ 1.84 | **≥ 0.034** (≥ 5.6 × 10⁻⁴ s⁻¹) |
| MFT-MFT | 0.795 | 0.110 | 86.1 % | 0.981 | 62 % | 1.43 | **0.0165** (2.7 × 10⁻⁴ s⁻¹) |
| cystine | 0.0171 | < 0.0004 | ≥ 97.6 % | 0.0061 | 18 % | ≤ 0.37 | ≥ 0.031 |

k_loss = ln(initial/left)/120 min, assuming first order with the trap holding the thiol at zero.
**One end point each, so first order is an assumption, not a finding.** Cystine: 82 % of its
sulfur is unaccounted for.

**Without a trap (Table 2), 100 °C, 2 h:**

| system | disulfide mass recovered | mixed disulfide / total initial | fraction of the 50 % statistical equilibrium | k_scramble (min⁻¹) = −ln(1 − f)/120 |
|---|---|---|---|---|
| MFT, water | **95.5 %** (191/200) | 21.5 % | 0.43 | **4.7 × 10⁻³** |
| MFT, benzene | 95.5 % (191/200) | 12.5 % | 0.25 | 2.4 × 10⁻³ |
| FFT, benzene | 86.8 % (316/364) | 20.1 % | 0.40 | 4.3 × 10⁻³ |
| FFT, water | **45 %** (90/200) | 5.5 % (12 % of survivors) | 0.11 (0.24) | ~1–2 × 10⁻³ |

- k_scramble is a **lower bound on the S–S cleavage rate per disulfide** (re-pairing with the same
  label is invisible). In water it sits within a factor ~3.5 of the trapped MFT-MFT loss rate,
  consistent with one cleavage step that is reversible unless the thiol is removed.
- **Net decay of the MFT-MFT pool in water without a trap: ≤ ~5 % in 2 h at 100 °C**, i.e.
  k_net ≲ 4 × 10⁻⁴ min⁻¹, at the precision of a single unreplicated isotope-dilution run.
- **FFT-FFT in water without a trap loses 55 % of its mass in 2 h** (k_net ≈ 6.7 × 10⁻³ min⁻¹). The
  labelled and unlabelled halves disagree (56 vs 23 µg from 100 + 100) [!], and there is no
  independent internal standard for this table (Table 1's IS for FFT-FFT is the reactant itself),
  so this number is the least secure in the paper.
- Benzene + H donor: (FFT-FFT + FFT) mass balance 78 % (no. 1, 12 h), 71 % (no. 2, BHT, 12 h); (MFT-MFT
  + MFT) 62 % (no. 3, cyclohexadiene, 2 h), 100 % (no. 4, hexadiene, 2 h).
- **Arrhenius transfer for context only (not a prior):** with the model's shared disulfide barrier
  (Zhang 2026, 122.2 kJ/mol, `ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL`) the factor 100 → 145 °C is
  ×69, putting the trapped MFT-MFT cleavage at ~1.1 min⁻¹ (log₁₀ ≈ +0.06) and the untrapped
  scrambling at ~0.3 min⁻¹ (log₁₀ ≈ −0.5) at 145 °C. B17's fit returned log₁₀ k_dimer_release =
  −8.26 at 145 °C, **unidentified and not shipped** (`results/validation/kinetic_core_b17_ship_rule.md`).
  The barrier is borrowed, the pH is unknown, and the trap shifts pH, so this only says "fast
  compared with anything B17 explored", not a value.

---

## 5. WHAT THE MODEL CAN TAKE FROM IT

**Context.** The sulfur lane destroys MFT/FFT by (a) oxidative dimerisation `ch_dimer_mft/fft`
(2 RSH + OX → disulfide, thiolate-weighted), (b) the thiolate-mediated non-disulfide sink
`ch_thiolate_loss_*` (B2.1, Kumazawa 2003 shape), (c) the unassigned `ch_mft_decay`/`ch_fft_decay`,
and (d) disulfide decay `ch_mftd_decay`/`ch_fftd_decay` (`k_dimer_decay`). B17 added the reverse
step `ch_dimer_release_*` (`k_dimer_release`, inert at zero, B17 not shipped). The unidentified
thiol sink is declared to act through the thiolate and has no measured constant.

**Verdict on G1 (quantified thiol loss kinetics): NO.** No thiol is heated as a starting material,
there is one time point per experiment, one temperature, no pH, no oxygen control. **It cannot
identify or bound the thiol sink itself.** What it does bound sits on the disulfide side:

1. **Disulfide formation is not a terminal sink in hot water.** MFT-MFT and FFT-FFT cleave in
   water at 100 °C on a timescale of tens of minutes to hours (k ≥ 5 × 10⁻³ min⁻¹ untrapped,
   1.6–3.4 × 10⁻² min⁻¹ trapped) and return thiol. A model in which `ch_dimer_*` is irreversible
   (release at zero) will **over-retain sulfur as disulfide** at ≥ 100 °C in water. This supports
   making `ch_dimer_release_*` non-zero, as an order-of-magnitude lower bound at 100 °C; it does not
   supply a value with a pH or temperature dependence.
2. **Bound on `k_dimer_decay` (MFTD):** net loss of the MFT-MFT pool ≤ ~5 % in 2 h at 100 °C in
   water without a trap → **k_dimer_decay(MFTD, 100 °C, pH unknown) ≲ 4 × 10⁻⁴ min⁻¹**, one run.
   For FFTD the same table gives ~55 % in 2 h (≈ 6.7 × 10⁻³ min⁻¹), fate unknown, flagged above.
   The model shares one `k_dimer_decay` across MFTD and FFTD; this paper (weakly) says the two
   differ by > 10×.
3. **The sink cannot be "thiol → disulfide" at 100 °C in water.** If the unidentified thiol sink
   were oxidation to the disulfide, the disulfide would hand the thiol back (point 1). So an
   irreversible sink at ≥ 100 °C must be a different chemistry. Guth offers one **candidate
   mechanism, qualitative only**: thiyl radicals coupling with **allylic radicals of 1,4-dienes**
   (thioethers; 22–38 % of sulfur in benzene with cyclohexadiene), and the authors suggest
   linoleic acid can play the donor. That is relevant only to lipid-containing matrices, it is
   measured in benzene, and it competes with H-abstraction that regenerates the thiol. **It says
   nothing about a thiolate dependence**; the paper varies no pH.
4. **Mechanism note for the thiolate claim:** the hydrolysis route proposed (RS–SR + H₂O → RSH +
   RSOH) is ionic, and pyridine-trapped runs were probably at higher pH than untrapped ones. This is
   compatible with, but does not test, a thiolate/hydroxide dependence of disulfide cleavage.

**Not usable:** any temperature or pH dependence (none); any number as a fit target for a thiol
(none); the benzene data as aqueous constants; the EPR as a rate.

**Caveats.** No replicates or SD; pH and headspace oxygen unknown; micromolar loadings (8.8–16 µM
disulfide) far below flavour-pot sulfur; the trap's own chemistry (Michael addition, base) is part
of the measured rate; Table 2 FFT-FFT has no independent internal standard.
