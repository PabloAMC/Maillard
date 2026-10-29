# Li, Schieberle & Steinhaus 2016/2017 — EXTRACTION (stable isotope dilution assays for H2S, methanethiol, ethanethiol, propane-1-thiol via 3-buten-2-one derivatisation, GC-GC-MS(CI); applied to durian pulp)
### A METHOD paper. It is the H2S method the designed experiment cites; it contains NO measurement of H2S formed by heating (cysteine, sugar + cysteine, meat or any Maillard model). Its only food data are four concentrations in raw durian pulp.

**Source on disk:** `data/articles/li2016b.pdf` (owner's download; 11 pages, Eur Food Res Technol (2017)
243:69-79). Read from the `pdftotext -layout` text layer in the scratchpad (2026-10-29); the text layer is clean.
Figures 1-7 are images (reaction scheme, MS-CI spectra, one chromatogram); no data are lost in them that matter
to a model (Fig. 3 is a single chromatogram, no numbers). Repo status before this dossier: cited only in
`docs/guides/DESIGNED_TIMECOURSE_PROTOCOL_DRAFT.md` (analyte table, H2S row, as "Li, Schieberle & Steinhaus
2016"); no dossier, no ledger entry (checked by grep).

## 0. Identity

| field | value |
|---|---|
| Title | "Development of stable isotope dilution assays for the quantitation of the food odorants hydrogen sulphide, methanethiol, ethanethiol, and propane-1-thiol and application to durian (Durio zibethinus L.) pulp" |
| Authors | Jia-Xiao Li, Peter Schieberle, Martin Steinhaus (corresponding) — Deutsche Forschungsanstalt für Lebensmittelchemie, Freising |
| Venue | Eur Food Res Technol (2017) 243:69-79; received 9 Mar 2016, revised 17 May 2016, accepted 1 Jun 2016, online 10 Jun 2016 |
| DOI | 10.1007/s00217-016-2723-z |
| Repo key | `li2016b` (year of online publication; the issue is 2017). The protocol draft cites it as 2016. |
| Lineage | GC-GC-MS system and durian screening from Li, Schieberle & Steinhaus 2012 (J Agric Food Chem 60:11253, ref. [6]); 3-buten-2-one thiol addition in water from Khatik 2006 (ref. [38]); Na/phenanthrene sulphide synthesis from Takata 2003 (ref. [42]); SAFE from Engel 1999 (ref. [43]) |

Abbreviations used below: ISC = isotopically substituted compound (the labelled standard); IUC = isotopically
unmodified compound (the analyte); AC = area count; MVK = 3-buten-2-one (methyl vinyl ketone).

## 1. Why it matters

The designed time-course protocol (`docs/guides/DESIGNED_TIMECOURSE_PROTOCOL_DRAFT.md`, analyte table) names this
paper as the H2S method, because "H2S supply over time ... caps both thiols" (item 5 of what the experiment must
identify). The sulfur lane has H2S formation from cysteine alone (Zheng & Ho 1994, `zheng1994_extraction.md`,
in the B2.3 fit) but no H2S time course in any heated cysteine + sugar system, and nothing on H2S loss. This
dossier records (a) the method in enough detail to judge whether it can deliver that observable, and (b) the
answer to "does this paper itself contain heated-system H2S data": **no** (section 6).

**Correction for the protocol draft (not edited here):** the draft's H2S row says "headspace, ³⁴S-sulfide SIDA".
The method is **not a headspace method**. It is an in-solution derivatisation of total sulphide (H2S + HS-) with
MVK, followed by solvent extraction, SAFE distillation and GC-GC-MS. The authors list "no headspace sampling
equipment is needed" as an advantage of the derivatisation route over headspace SIDA.

## 2. Labelled standards and how they were made

| standard | source | how made / purity | stock solution |
|---|---|---|---|
| sodium (³⁴S)sulphide (ISC for H2S) | synthesised from ³⁴S8 (Aldrich, 99.5 % isotopic purity) | ³⁴S powder 144 mg (4.2 mmol S) + phenanthrene 300 mg (1.7 mmol, electron-transfer catalyst) added to a dispersion of sodium 654 mg (9.4 mmol) in liquid paraffin with anhydrous 1,2-dimethoxyethane 15 mL; stirred under argon, 70 °C, 4 h; hexane 30 mL added, centrifuged, solid washed with hexane (2 × 30 mL), dried: **282 mg raw product, yellow powder; purity 52 % by SIDA/GC-GC-MS of its MVK derivative; yield 44 %**; stored under argon | 46 µg raw product per mL in phosphate buffer pH 7.4; **exact concentration determined in triplicate against the unlabelled Na2S stock: 0.30 µmol/mL** |
| sodium sulphide nonahydrate (IUC reference) | Aldrich | commercial | 0.28 µmol/mL in phosphate buffer pH 7.4 |
| sodium (²H3)methanethiolate (ISC for MeSH) | synthesised from (²H6)dimethyl sulphide (Aldrich) | Na 0.80 g (34.8 mmol) in condensed NH3 15 mL (dry ice/2-propanol), argon, 15 min; (²H6)DMS 1.36 g (20 mmol); 30 min; NH3 evaporated: **2.30 g crude; purity 54 %; yield 97 %** (on sodium: 2 Na per cleavage) | 41.8 µg raw product per mL in phosphate buffer pH 9.0 |
| sodium methanethiolate (IUC) | Aldrich | commercial | 0.29 µmol/mL in phosphate buffer pH 9.0 |
| (²H5)ethanethiol, (²H7)propane-1-thiol | C/D/N Isotopes | commercial | in ethanol (concentration not printed) |
| ethanethiol, propane-1-thiol | Aldrich | commercial | in ethanol (concentration not printed) |

Checks [Z]: 282 mg × 0.52 / 79.95 g/mol (Na2³⁴S) = 1.83 mmol of a possible 4.2 mmol → 44 % (matches). 46 µg/mL ×
0.52 / 79.95 = 0.299 µmol/mL, matching the 0.30 µmol/mL determined by reverse isotope dilution. Why ³⁴S: "Due to
the acidity of hydrogen sulphide, stable isotopic substitution is only possible in the sulphur atom"; ³⁴S8 was the
only commercial isotopologue of sufficient purity. The paraffin/sulphur/diatomaceous-earth route (H2S trapped in
alkali) was tried and judged ineffective at the price of ³⁴S8 (data not shown). The Na/phenanthrene route avoids
oligosulphide by-products (per ref. [42]).

Reference derivatives (for peak assignment), all synthesised and NMR-confirmed (¹H/¹³C, COSY/HSQC/HMBC):

| derivative | of | synthesis | purity / yield | MS-CI ions |
|---|---|---|---|---|
| 4,4′-sulphanediyldibutan-2-one | H2S | Na2S·9H2O 2.40 g (10 mmol) + DCM 12 mL + MVK 3.50 g (50 mmol) + acetate buffer pH 5.0 10 mL, argon, 3 h at −60 °C then 15 h ambient; flash chromatography | > 99 % (GC-FID), 75 % | m/z 105 (100, loss of MVK), 175 (7, M+H+); ³⁴S analogue m/z 107 / 177 |
| 4-(methylsulphanyl)butan-2-one | MeSH | NaSMe 0.35 g + MVK 0.42 g in phosphate pH 9.0, overnight, RT | 92 %, 84 % | 119 (M+H+); d3: 122 |
| 4-(ethylsulphanyl)butan-2-one | EtSH | EtSH 0.75 g + MVK 0.70 g in water, 30 min RT | 97 %, 81 % | 133 (M+H+); d5: 138 |
| 4-(propylsulphanyl)butan-2-one | PrSH | PrSH 0.91 g, same procedure | 96 %, 85 % | 147 (M+H+); d7: 154 |

## 3. Sample preparation — how H2S is trapped and stabilised

**Trapping principle:** H2S/HS- is converted in the aqueous sample, in a closed vial, to the non-volatile,
non-oxidisable bis-adduct 4,4′-sulphanediyldibutan-2-one by Michael addition to an excess of MVK. The labelled
³⁴S-sulphide is added before derivatisation and equilibrated with the sample, so every later loss (extraction,
SAFE, concentration, cyclisation of the derivative) is compensated by the ratio.

**H2S workup (verbatim sequence, "Quantitation of hydrogen sulphide"):**

1. 20-mL septum-sealed crimp-top vial: phosphate buffer pH 7.4 (8 mL) + ³⁴S-sulphide stock (2 mL = 0.6 µmol;
   0.2-2 mL for calibration) + either unlabelled Na2S stock (2 mL; 0.2-2 mL for calibration) or durian pulp
   (1-10 g).
2. Homogenise by magnetic stirring, ambient temperature, **30 min** (standard-analyte equilibration).
3. Add MVK **200 µL** (≈ 173 mg, ≈ 2.5 mmol [Z], i.e. > 400-fold over the 0.6 µmol standard), stir
   **overnight** (ambient).
4. Extract with dichloromethane (2 × 25 mL); dry over anhydrous Na2SO4.
5. Remove non-volatiles by SAFE at 40 °C.
6. Concentrate distillate to **1 mL** (Vigreux column).
7. Inject 0.5-2 µL into GC-GC-MS.

Whether the MVK was added through the septum (vial kept closed) is not stated; the vial is described as
septum-sealed. Total time per sample [Z]: ~0.5 h + overnight (~15-18 h) + extraction/SAFE/concentration (a few
hours) + one GC-GC-MS run (target elutes at ~24 min in the second dimension) → roughly one working day plus one
night per batch; samples can be run in parallel batches.

**Why pH 7.4:** "selected as optimum between the high volatility of hydrogen sulphide (pKA 7.0) at low pH and
the susceptibility of the open-chain derivative to intramolecular cyclisation at high pH."

**Thiols:** same workup, but phosphate pH 9.0 for methanethiol (pKa 10.3) and plain water for ethanethiol and
propane-1-thiol (handled as free thiols; b.p. 35 and 68 °C). **H2S and methanethiol need separate aliquots**
(different buffers); ethanethiol and propanethiol share a water workup.

## 4. Instrument

- **GC-GC-MS (quantitation):** two-dimensional heart-cut system as in ref. [6]. First dimension FFAP,
  30 m × 0.32 mm i.d., 0.25 µm (Agilent J&W); second dimension DB-5, 30 m × 0.25 mm i.d., 0.25 µm. Heart-cut
  window 1 min around the target, cold-trap refocusing, re-chromatography on DB-5; ion-trap MS in **chemical
  ionisation, methanol reagent gas**. Both ovens 40 °C start, 6 °C/min; 4,4′-sulphanediyldibutan-2-one elutes
  at ~27 min (FFAP) and ~24 min (DB-5).
- **GC-MS (EI, for reference spectra):** HP 5890 GC + Finnigan MAT 95, 70 eV.
- **Quantifier ions:** m/z 105 / 107 (H2S derivative, ³²S / ³⁴S fragment after MVK loss — the weak M+H+ 175/177
  is not used); 119 / 122 (MeSH); 133 / 138 (EtSH); 147 / 154 (PrSH).
- NMR: Bruker Avance 400 (structure confirmation only).

## 5. Calibration, recovery, LOD/LOQ, precision, matrix effects

### 5a. Calibration (all printed)

AC(ISC)/AC(IUC) = m × n(ISC)/n(IUC) + b, molar ratios **1:10 to 10:1**, linear regression in Excel 2010;
analyte n(IUC) = m × n(ISC added) / [AC(ISC)/AC(IUC) − b].

| analyte | m (slope) | b (intercept) | R² | calibration range, absolute (per workup) | equivalent in durian (1-10 g) |
|---|---:|---:|---:|---|---|
| H2S | 0.9113 | −0.0311 | 1.000 | ~2 to 200 µg H2S (with 0.6 µmol ³⁴S standard) | ~200 µg/kg to ~200 mg/kg |
| methanethiol | 0.7878 | −0.0015 | 0.9997 | not printed in absolute units | ~300 µg/kg to ~300 mg/kg |
| ethanethiol | 1.2457 | −0.0943 | 0.9981 | not printed | ~1 mg/kg to ~1 g/kg |
| propane-1-thiol | 1.4562 | −0.1624 | 0.9999 | not printed | ~250 µg/kg to ~250 mg/kg |

Natural-abundance correction used when determining the ³⁴S stock concentration (not in the calibration
equation itself, which absorbs it into m and b): AC(IUC) = AC(m/z 105)/0.9502;
AC(ISC) = AC(m/z 107) − 0.0443 × AC(m/z 105) (³⁴S NA 0.0421, ³²S NA 0.9502). No m/z 105 signal from the ISC.

Check [Z]: 2 µg H2S = 0.059 µmol and 200 µg = 5.9 µmol, i.e. 1:10 and 10:1 against the 0.6 µmol standard —
consistent with the stated ratio range.

### 5b. Recovery

**Not reported.** No spike-recovery experiment, absolute or relative, for any of the four analytes. The design
relies on SIDA compensation of losses. (The ³⁴S stock concentration determination in triplicate is the only
replicate method-validation number besides the durian CVs.)

### 5c. LOD / LOQ

**No LOD or LOQ is stated for any analyte.** The lower end of the calibrated range is the only quantitative
floor: **~2 µg H2S per workup (≈ 0.06 µmol), i.e. ~200 µg/kg at 10 g sample**. The authors argue (no
experiment) that the limit of quantitation could be pushed below the 10 µg/kg odour threshold of H2S in water
"by e.g. increasing the sample size to 50 g, reducing the amount of added standard to 20 %, and further
concentration of the SAFE distillate to 200 µL." [Z]: with 0.12 µmol standard and the 1:10 ratio floor, that is
~0.4 µg per workup, ~8 µg/kg at 50 g — a projected, not demonstrated, LOQ.

### 5d. Precision

- Calibration linearity R² 0.9981-1.000 (table above); number of calibration points not printed.
- Durian: mean of **triplicate workups**; CV 6 % (H2S), 13 % (MeSH), 20 % (EtSH), 22 % (PrSH) (Table 1).
- ³⁴S-sulphide stock: triplicate, 0.30 µmol/mL (no SD printed).
- No inter-day precision, no inter-laboratory comparison, no reference material.

### 5e. Matrix effects

- Durian pulp: in the 200 µg/kg-200 mg/kg range "peak shapes were excellent and no matrix interferences were
  visible in the specific ion chromatograms" (Fig. 3). ~4 g pulp + 1 mL final distillate gave the optimum
  (³²S/³⁴S ratio ~1).
- **Blanks without MVK** were run in parallel to show the four derivatives are not natively present in durian;
  none was "in relevant amounts". The authors note that 4-(methylsulphanyl)butan-2-one (the methanethiol
  derivative) **occurs naturally in some thermally treated foods** (mushrooms, roasted pistachios, refs. [46,
  47]); a no-MVK blank is therefore obligatory in heated samples.
- No test in any other matrix (no heated model system, no meat, no buffer series, no high-thiol or high-amine
  matrix).

## 5f. Pitfalls the authors name

1. **Volatility and oxidation losses:** the high volatility of H2S and short-chain thiols "and their tendency
   towards the oxidative formation of disulphanes ... may lead to substantial losses during the analytical
   workup" (Introduction; refs. [2, 33, 34]). This is the stated reason for SIDA + derivatisation.
2. **Multiple MVK products from H2S** (Fig. 1): 4-sulphanylbutan-2-one (mono-adduct; negligible only if MVK is
   in high excess), the target bis-adduct, and cyclic aldol products (1-(4-hydroxy-4-methyltetrahydro-2H-
   thiopyran-3-yl)ethanone diastereomers; 1-(4-methyl-3,6-dihydro-2H-thiopyran-3-yl)ethanone) formed by
   intramolecular aldol reaction **under basic conditions**. Preparative reaction in DCM/phosphate pH 7.4 gave a
   mixture of the bis-adduct and both cyclic diastereomers; in acetate pH 5.0 it gave the bis-adduct only.
   Countermeasure: "3-buten-2-one was applied in high excess and strong basic conditions were avoided".
3. **pH trade-off** (pKa 7.0 volatility vs cyclisation) — section 3.
4. **Unassigned by-products** in one-dimensional GC-MS of the derivatised Na2S standard ("none of these could be
   assigned to any of the structures depicted in Fig. 1") — the reason for heart-cut GC-GC.
5. **Natural ³⁴S in the analyte** (4.4 % of the m/z 105 area appears at m/z 107) — corrected as above.
6. **Low purity of the synthesised standards** (52 %, 54 %): the standard's concentration must be determined by
   reverse isotope dilution, not by weighing.
7. **Native derivative in heated foods** (methanethiol adduct) — section 5e.

**Not discussed:** artefactual formation of H2S during workup (e.g. from cysteine, thiazolidines, other
sulphur precursors during the overnight step or SAFE at 40 °C); stability of underivatised H2S in a stored
sample; consumption of MVK by other nucleophiles (thiols, amines) in a sample; behaviour at pH other than 7.4.
None of these is tested.

## 6. Every quantitative result

### 6a. Table 1 — durian pulp ('Monthong', tree-ripened, organic, Thailand; optimum eating ripeness)

"Concentrations and OAVs of hydrogen sulphide and aroma-active alkane-1-thiols in durian pulp." Mean of
triplicate workups; CV = SD/mean; odour threshold = orthonasal, in water (DFA in-house database, ref. [8]).

| odorant | odour quality | conc. (µg/kg) | CV (%) | SD [Z] (µg/kg) | µmol/kg [Z] | odour threshold (µg/kg) | OAV |
|---|---|---:|---:|---:|---:|---:|---:|
| hydrogen sulphide | rotten egg | 3340 | 6 | ~200 | 98 | 10 | 330 |
| methanethiol | rotten, cabbage | 26,700 | 13 | ~3500 | 555 | 0.59 | 45,000 |
| ethanethiol | rotten, onion | 45,400 | 20 | ~9100 | 731 | 0.095 | 480,000 |
| propane-1-thiol | rotten, durian | 3070 | 22 | ~680 | 40 | 0.49 | 6300 |

OAVs re-derived [Z]: 334, 45,250, 477,900, 6270 — all match the printed rounding. Molar masses used: H2S 34.08,
MeSH 48.11, EtSH 62.13, PrSH 76.16 g/mol. Basis: per kg fresh pulp (fresh weight; moisture not stated). Raw,
unheated fruit: **not simulable** by the kinetic core.

### 6b. Method numbers

Recorded in sections 2 and 5: synthesis yields/purities, stock concentrations (0.28, 0.29, 0.30 µmol/mL),
calibration m/b/R², calibration ranges, NA correction factors.

### 6c. Model-system, recovery and heated-system experiments

**None.** The only "experiments" besides durian are the three method-development runs (unlabelled Na2S in pH 7.4
buffer; ³⁴S-Na2S in pH 7.4 buffer; ³⁴S-Na2S + durian pulp at several masses), reported qualitatively.

### 6d. Answer: does the paper contain ANY measurement of H2S formed by heating?

**No.** There is no heated cysteine system, no sugar + cysteine system, no meat or Maillard model, no time
course and no fixed-time heated sample. The only heating in the whole paper is in syntheses (70 °C, 4 h for
Na2³⁴S) and in SAFE (40 °C, a workup step). Every mention of thermal H2S is a literature citation in the
Introduction:

- Introduction: "In food, hydrogen sulphide may be formed either enzymatically or thermally from cysteine or
  cysteine-containing peptides."
- Introduction: "The thermal formation of hydrogen sulphide from cysteine follows the Strecker degradation path
  [14–17]." Refs. [14] = Zheng & Ho 1994 (ACS Symp Ser 564:138; already dossiered, `zheng1994_extraction.md`,
  ledger `zheng1994_table1_cysteine_h2s_release_k`, in the B2.3 fit), [15] Hofmann & Schieberle 1995
  (ribose/cysteine), [16] Hofmann & Schieberle 1997 (glucose/rhamnose + cysteine), [17] Kobayasi & Fujimaki 1965
  (already dossiered, `kobayasi1965_extraction.md`).
- Introduction: "Under Maillard reaction conditions, methanethiol may also be formed from cysteine [15, 16]."
- Results: 4-(methylsulphanyl)butan-2-one "has been reported in some foods, such as mushrooms [46] and roasted
  pistachios [47], after thermal treatment" (native derivative; not an H2S datum).

The paper adds nothing to H2S formation or loss kinetics. It supplies a method only.

## 7. What the method implies for the designed experiment's H2S analyte

Model predictions to compare against (final H2S, `results/validation/sulfur_flux_attribution.md`, B60
calibration), converted [Z] at 34.08 g/mol, and their position against the demonstrated calibration window
(0.06-6 µmol per workup with 0.6 µmol ³⁴S standard; window scales with the amount of standard added):

| reference pot | model final H2S (mM) | µg/mL | µmol in a 1 mL aliquot | µmol in a 10 mL aliquot | verdict |
|---|---:|---:|---:|---:|---|
| hofmann1998 ribose+Cys pH 3, 145 °C 20 min | 0.901 | 30.7 | 0.90 | 9.0 | 1 mL is mid-window |
| cao2017 xylose+Cys 90 °C 300 min | 0.242 | 8.2 | 0.24 | 2.4 | 1 mL fine |
| food-like ribose+Cys 140 °C 30 min | 0.23 | 7.8 | 0.23 | 2.3 | 1 mL fine |
| hofmann1998 pH 5, 145 °C 20 min | 0.0631 | 2.15 | 0.063 | 0.63 | 1 mL at the floor; use 5-10 mL or ~0.1 µmol standard |
| zhai2019 xylose+Cys 90 °C 160 min | 0.0382 | 1.30 | 0.038 | 0.38 | needs 5-10 mL or reduced standard |
| yiltirak ribose+Cys 120 °C 60 min | 0.0287 | 0.98 | 0.029 | 0.29 | needs 5-10 mL or reduced standard |
| zhang2025b xylose+Cys pH 6.5, 120 °C 50 min | 0.00931 | 0.32 | 0.009 | 0.093 | 10 mL + ~0.1 µmol standard |
| hofmann1998 pH 7, 145 °C 20 min | 3.1e-5 | 0.001 | 3e-5 | 3e-4 | **below even the projected LOQ** (~0.012 µmol, 50 g / 20 % standard) by ~40×: a "<LOQ" result is the expected outcome; a quantifiable H2S at pH 7 would falsify the model |

Implications [Z — my inference from the method; none of it is tested in the paper]:

1. **Sensitivity is adequate** for the model's 0.01-0.9 mM range with 1-10 mL aliquots: those levels are
   0.3-31 mg/L, i.e. 0.1-9× the durian H2S (3.3 mg/kg, 0.098 mmol/kg) that the method measured at CV 6 %; at
   10 mL even 0.3 mg/L gives ~3 µg per workup, above the ~2 µg calibrated floor. The ³⁴S standard
   should be dosed to the expected level (ratio within 1:10-10:1, ideally ~1), e.g. 0.1 µmol for pH 5-6 pots,
   0.6 µmol for pH 3 pots. The pH 7 prediction (~nM-µM) is not measurable by this method as published.
2. **Sample volume:** 1-10 mL of reaction mixture per time point; the published vial (20 mL crimp) holds 8 mL
   buffer + 2 mL standard + ≤ 10 g sample, so a 5-10 mL reaction fill in the same 20-mL vial is compatible.
3. **A sealed hot pot cannot be sampled with this method; it must be quenched closed, spiked through the septum,
   and derivatised in the reaction vial.** H2S is volatile (pKa 7.0; dimensionless gas/water partition ~0.4 at
   25 °C for undissociated H2S, standard Henry's-law data, not from this paper): in a 20-mL vial with 10 mL fill
   at pH 5, ~30 % of the sulphide sits in the headspace after cooling, and it is lost the moment the vial is
   opened or an aliquot is withdrawn. The method's own logic (equilibrate standard for 30 min in a closed
   septum vial, then derivatise) only compensates losses that occur after the standard has equilibrated.
   Therefore: run the H2S time points in septum-crimped vials that double as derivatisation vials; quench on
   ice; inject ³⁴S-sulphide stock (in pH 7.4 phosphate) through the septum, stir 30 min, inject MVK through
   the septum, stir overnight; only then open. This also fixes the "separate vials per time point" design.
4. **Buffer/pH at derivatisation:** the method was validated at pH 7.4 only. The designed pots are in 0.2 M
   phosphate at pH 5-7 (cooled). The preparative bis-adduct synthesis worked cleanly at pH 5.0 (acetate), so
   derivatising at pH 5-7 is plausible, but untested in this assay. Adding pH 7.4 buffer via the septum
   (as the method does with 8 mL) would dilute and pressurise a filled vial; a partner lab must validate one or
   the other.
5. **Workup time and throughput:** ~1 day + overnight per batch, plus SAFE and heart-cut GC-GC-MS(CI). For ~130
   vials (or ~75 in the tight design) this is the most labour-intensive analyte in the table. A cheaper
   alternative for a first pass is the protocol's own fallback, methylene blue (ref. [31] of the paper), which
   is not isotope-corrected.
6. **Matrix-specific risks not covered by the paper** (must be validated in a heated ribose + cysteine matrix by
   spiking ³⁴S-sulphide into a heated blank and into an unheated ribose + cysteine blank):
   - MVK is consumed by cysteine and other thiols/amines in the sample (30 mM Cys × 10 mL = 0.3 mmol vs ~2.5 mmol
     MVK): still ≥ 8-fold excess, but the excess condition that suppresses the mono-adduct is weaker than in
     durian. SIDA compensates if the ³⁴S standard reacts identically, which it should.
   - Artefactual H2S from cysteine, thiazolidines or Cys-Amadori compounds during 30 min + overnight at room
     temperature and SAFE at 40 °C: not tested. Likely small at room temperature, but the unheated
     ribose + cysteine (t = 0) vial is the control that measures it.
   - Native MVK adducts in heated Maillard samples (MVK itself, and 4-(methylsulphanyl)butan-2-one, are thermal
     products): run a no-MVK blank per condition (the paper's own blank design).
   - Oxidation of H2S to polysulphides/elemental sulphur in an air headspace during heating is model chemistry,
     not a workup loss: the method measures free sulphide at quench, not total H2S produced.
7. **What it measures:** free sulphide (H2S + HS-) present at quench, in liquid + headspace of the closed vial.
   It does not measure H2S that left the system (open vessels), nor sulphur bound in thiols, polysulphides or
   heterocycles. In the model this is the H2S pool after r_h2s_loss and the thiol-forming sinks — the right
   observable for "H2S supply over time", provided the vial stays closed from fill to derivatisation.

## 8. Ledger entries (`data/lit/measurement_ledger/sulfur.yml`)

- `li2016b_text_h2s_thiol_sida_method_validation` — kind other, digitised_in_dossier (sections 2-5).
- `li2016b_tab1_durian_h2s_thiols` — endpoint, digitised_in_dossier, simulable false (section 6a).
- No heated-system entry: the paper has none (section 6d).
