# Charles-Bernard, Kraehenbuehl, Rytz & Roberts 2005 (Part 1) — EXTRACTION
### *"Interactions between Volatile and Nonvolatile Coffee Components. 1. Screening of Nonvolatile Components"*
### Medium-roast spray-dried coffee extract fractionated by SPE, dialysis (3.5 kDa), Sephadex G-25 SEC and DEAE anion exchange; each fraction reconstituted in 10 mM acetate pH 5.2 and dosed with a five-volatile mix (FFT, EtSH, DMS, pyrrole, 2,3-pentanedione); loss read as SPME-GC-MS headspace area relative to a matrix-free blank, at 1 h, 24 h, or every 2 h for 18.5 h, at 25 C. Plus an O-phenylhydroxylamine consumption assay ("electrophilic reactivity") and amperometric oxidizability of the anion-exchange fractions. **One temperature (25 C), one pH (5.2), no MFT, no product measured, no reversibility test.**

**Read 2026-10-31 from `data/articles/charles-bernard2005a.pdf`** (9 pp., JAFC 53:4417-4425; born-digital text
layer, clean; Tables 1-3 re-typed from it). **Every figure (1-6) was checked against 300 dpi renders, and Fig. 4 and
Fig. 2's ion-exchange rows against 600 dpi crops.** Figure values are my readings. Read error: about ±1.5 percentage
points on the bar charts (Figs. 1, 2, 3, 6), ±2 points on Fig. 4 above 10 %, worse below 5 %. Fig. 5: ±0.005 mmol/g
and ±30 nA·s. The paper says the intergrid distance of every figure is its 5 % LSD. That is ~8.3 points on Figs. 2, 3
and 6 and ~12.5 points on Figs. 1 and 4 [Z from F], so the paper's own noise is larger than my read error.

**DOI:** printed `10.1021/jf048021q` (p. 4417 footer and `JF048021Q`). Crossref 2026-10-31 returns the same title,
vol. 53, pp. 4417-4425, and the authors Charles-Bernard, Kraehenbuehl, Rytz and Roberts.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown · **[Z]**
derived by me, never printed · **[!]** flag.

---

## 0. How Part 1 relates to Part 2 (`charlesbernard2005_extraction.md`, JAFC 53:4426)

| aspect | Part 1 (this paper) | Part 2 | same? |
|---|---|---|---|
| laboratory | Nestlé Research Center, Lausanne; Charles-Bernard, Kraehenbuehl, **Rytz** (statistics), Roberts | same group, without Rytz | **same lab** |
| coffee | blend of Santos 50 % + Colombia 30 % (both Arabica) + Ivory Coast Robusta 20 %; CTN 85 medium roast; countercurrent extraction 3 x 12 min at 120 C; spray-dried; **extraction yield 22.7 %** [M] | Arabica 80 / Robusta 20, CTN 85, yield 22.7 %, preparation "as described earlier (7)" | **almost certainly the same extract** (blend, roast and yield all match). Part 2 adds no fractionation |
| matrix in the pot | 1 % t.s. whole brew, OR fractions at their relative-yield concentration, OR at equal mass (1 % or 0.13 %) | 1 % t.s. whole brew only | Part 1 is the superset |
| buffer, pH, T | 10 mM acetate pH 5.2; 25 C in the autosampler; mixing at RT | 10 mM acetate pH 5.2; 25 C | **same** |
| thiols | FFT and EtSH (+ DMS, pyrrole, 2,3-pentanedione), at ppm: FFT "4 mg/mL" [!] (read 4 mg/L = **35 uM** [Z]), EtSH 2 mg/L = **32 uM** [Z] | ten thiols at 3.2-9.7 uM (FFT 3.9, EtSH 5.4 uM) | **thiol charge ~6-9x higher in Part 1** |
| measurand | SPME (PDMS/DVB 65 um, 1 min) headspace area of the sample / area of a matrix-free blank vial measured at the same time; duplicates; DB-Wax; same oven programme | identical construction | **same method** |
| kinetics | every 2 h, 0.5-18.5 h (Fig. 4); no rate constants printed | every 2-4 h to ~15 h; k(obs) printed (Table 2, unit misprinted) | same design; Part 2 prints constants |
| electrophile probe | **O-phenylhydroxylamine consumption by HPLC**, mmol reacted / g dry (Fig. 5) | hydroxylamine dose that stops protecting thiols ("8-10 mmol/g") | **different probe, different read-out**; the two are never compared in either paper (sec. 6, D1) |
| O2 | glovebox only for the amperometry samples; the one O2 statement is cited (ref 5) [C] | air vs N2 measured (Table 3) | Part 2 owns the O2 evidence |

**Bottom line on the pair:** same laboratory, same extract, same pot (1 % t.s., pH 5.2, 25 C) and same headspace-ratio
method. Part 1 asks **which fraction** removes thiols; Part 2 asks **how**. Part 1's whole-brew EtSH kinetics reproduce
Part 2's EtSH constant to within ~25 % (sec. 3.4), so the two can be treated as one data set measured on one material.

---

## 1. What this paper answers for the model

The model's questions: which matrix fractions bind or destroy thiols, and how (EX-050's MFT loss law and the matrix
layer: `MELE` / `k_thioether` in `src/kinetic_core/species_sulfur.py` and `parameters_sulfur.py`; protein and
melanoidin binding of thiols, STATUS_AND_PLAN sec. 4 item 7).

| need | what Part 1 gives | tag |
|---|---|---|
| (a) **which fraction destroys thiols** | Every size class does, given enough time: both dialysis halves (<3.5 and >3.5 kDa, at relative-yield concentration) leave **FFT = EtSH = 0 at 24 h** (Fig. 1 caption). At equal mass (1 % t.s., 1 h), EtSH reactivity is U-shaped in MW: SEC-C (medium MW) is the most inert (EtSH ~70 % left); the phenol-rich low-MW SEC-I and SEC-J leave 0 %; the high-MW SEC-A/B leave ~37-38 %; whole brew 16 % (Figs. 2, 3). FFT is 0 at 1 h in every SEC fraction | [M] caption, [F] |
| (b) **charge** | Strong cation exchange removes most of the fast sink: the SCX and Dowex eluates (1 % t.s., 1 h) leave **FFT 73 / 71 %** and **EtSH 57 / 61 %**, against FFT 0 / EtSH 16 % for the whole brew. SAX (acids removed) leaves EtSH 0 %. DEAE: AEC-A (neutral + cationic, **78 % of the solids**) carries most of the activity; AEC-C (13 %) is potent; the anionic, darkest AEC-B and AEC-D are weak | [F], [Q] |
| (c) **pure components** (at brew concentrations) | Chlorogenic acid (1.0 g/L) and caffeic acid (concentration not printed [!]) take **FFT to 0 % at 24 h** and leave EtSH, pyrrole and pentanedione unaffected. Caffeine 0.5, trigonelline 0.02, N-methylpyridinium iodide 0.03, larch and acacia arabinogalactan 3 g/L: no significant effect | [M] conc., [Q] |
| (d) **the O2 dependence of the CGA + FFT loss** | "completely suppressed in the absence of air", cited to the authors' 2003 proceedings (ref 5); not shown here | [C] |
| (e) **class attribution** | Pearson r over 13 samples (24 h losses), Table 3: melanoidins (by difference) 0.41 FFT / 0.83 EtSH; protein 0.41 / 0.83; carbohydrates 0.41 / 0.73; chlorogenic acids 0.33 / 0.63; **overall concentration 0.39 / 0.80**. Every class is collinear with concentration (r 0.97-0.98), and FFT is censored at 100 % loss in 10 of 13 samples. **The correlations cannot apportion the sink between classes** | [M], [Z] |
| (f) **protein** | Authors: "might play a minor role". The r values match the melanoidins' r exactly. No native-protein test | [Q], [Z] |
| (g) **reversibility** (covalent vs disulfide vs physical) | **Not tested.** No release, dilution, reducing-agent or product measurement. The authors infer a chemical reaction from the pseudo-first-order decay to zero [Q]. Thiol/disulfide interchange with matrix protein is excluded by an unshown experiment: "No decrease of dimethyl disulfide was observed (data not shown)" [V]. A ~10-12 % instantaneous headspace drop (fit intercepts, sec. 3.4) may be physical partition [Z] | [Q], [Z] |
| (h) **matrix electrophile density** | Whole brew consumes **~0.29 mmol O-phenylhydroxylamine per g dry** at 13 h, still rising slowly (Fig. 5); the fractions sum to ~0.21 (72 %). **About 30x below Part 2's "8-10 mmol/g"**; the reagent charge itself caps this assay near 0.7-0.9 mmol/g [Z] | [F], [Z] |
| (i) **order of the sink in matrix loading** | AEC-A at 7.8 vs 1.3 g/L (6.0x): the EtSH pseudo-first-order k rises only **~3.1x** (apparent order ~0.6 in matrix) [Z from F]. One pair of loadings, with a lag phase | [Z] |
| (j) **a loss constant at process temperature** | **None.** One temperature (25 C), one pH (5.2), 0.01-1 % solids. No constant is printed. The only kinetic data are Fig. 4 curves, which I converted to 25 C constants myself (sec. 3.4) | — |
| (k) MFT | **Not studied.** FFT and EtSH are the only thiols | [M] |

**Answer to the brief's last question: Part 1 gives NO loss constant usable at process temperature.** Its kinetic
content is figure-only, at 25 C, pH 5.2, in reconstituted coffee fractions. It prints no rate constant, no
temperature series and no Ea. The constants in sec. 3.4 are my fits to figure readings, valid **only at the
conditions tested**. They cannot set λ(110 C) for MFT, and they do not change EX-050's INCONCLUSIVE verdict.

What it does add, weighed:
1. **A same-lab, same-extract direct measurement that undercuts the 8-10 mmol/g site density** (sec. 6, D1). It
   supports the Part 2 dossier's Flag 1 and the hofmann2001 capacity bound (>= 0.028 mmol/g). P(the thiol-reactive
   electrophile density of 1 % brew solids is < 1 mmol/g) ~0.7. The O-phenylhydroxylamine assay sees carbonyls, not
   necessarily the thiol-reactive Michael acceptors, so this is not settled.
2. **The fast thiol sink is in the cation-exchange-retained, neutral-to-cationic, coloured material**, not in the
   carbohydrates, caffeine, trigonelline or (for EtSH) the free chlorogenic acids. Weight is moderate. Cation exchange
   also removes metals and K+ and may acidify the eluate, and neither the pH nor the metal content of the
   reconstituted fractions is reported (sec. 2, A3). Part 2's DTPA null (1/k_rel 0.5-1.0) argues against the metal
   reading.
3. **FFT is 8-20x more labile than EtSH in every fraction** [Z, sec. 3.4]. In the pure CGA and caffeic acid pots it
   is the only thiol attacked. This agrees with Part 2 (FFT > 6.5x EtSH, a floor).
4. **Equal-mass comparisons break the paper's "run in parallel" claim** (sec. 6, D3). Per gram, the most
   thiol-active fraction (AEC-A) has the lowest O-phenylhydroxylamine reactivity and the lowest oxidizability.
   Neither assay is a proxy for `MELE`.

---

## 2. Conditions

| item | as printed / derived |
|---|---|
| matrix | spray-dried medium-roast extract (CTN 85), reconstituted at **1 % t.s. (10 g/L)** "to mimic a brew", or its fractions; purified components at brew concentrations [M] |
| fraction concentrations (interaction tests) [M] | dialysis permeate **2.0** g/L, retentate **7.8** g/L [! likely swapped: yields are 76.5 / 19.5 %, so relative-yield values would be ~7.65 / 1.95]; SEC-A 0.90, SEC-D 1.0, SEC-E 1.2, SEC-F 1.7 g/L; AEC-A 7.8, AEC-B 0.41, AEC-C 1.3, AEC-D 0.15 g/L. **Figs. 2 and 3 are at equal mass, 1 % t.s. per fraction** (text p. 4421). Fig. 4 adds AEC-A at 0.13 % (equal mass with AEC-C) |
| volatiles (final) [M] | DMS 1 mg/L, FFT "4 mg/mL" [!], EtSH 2 mg/L, pyrrole 10 mg/L, 2,3-pentanedione 10 mg/L. My molarities [Z]: DMS 16 uM, FFT 35 uM (reading mg/L, as the "ppm level" sentence requires), EtSH 32 uM, pyrrole 149 uM, pentanedione 100 uM. Stable 24 h in buffer [Q] |
| buffer, pH | 10 mM acetate, pH 5.2 (label). **The pH of the reconstituted fractions is never measured** [!]. Amperometry: 100 mM acetate pH 5, degassed. Electrophile assay: no buffer stated |
| mixing, vessel | 2 mL aroma solution + 2 mL matrix solution (each at 2x), 4 mL amber vials, stirred **15 min at RT**; 800 uL moved to 2 mL amber silanized vials (1.2 mL air headspace); **equilibrated 24 h at 25 C** in the autosampler, or 1 h for the screens (Figs. 2, 3); kinetics: blank and sample vials alternated from t = 0, sampled every 2 h, 0.5-18.5 h |
| temperature | 25 C (autosampler) for every interaction read. RT for mixing and for fraction concentration. Dialysis 4 C, 5 d. Extraction (industrial) 120 C. **No temperature factor** |
| O2 | air headspace, sealed vial; my estimate ~10 umol O2 per vial against ~0.03 umol thiol, so O2 is not limiting [Z]. No anaerobic arm for the fractions |
| water activity | dilute solutions (~0.01-1 % solids), a_w ~1 |
| SPME-GC-MS | PDMS/DVB 65 um, 1 min headspace equilibration, 1 min absorption (chosen "so as to not saturate the fiber"); linearity of each compound checked alone and in the mix; desorption 240 C 5 min (purge off 3 min); HP 5973, DB-Wax 30 m, 35 C (3 min) to 170 C at 4 C/min, to 220 C at 20 C/min; scan 29-300 amu. **No internal standard, no isotope dilution** |
| measurand | "relative headspace concentration" = sample area / blank area at the same sampling time, %; "loss" = 100 % minus that. Duplicates at least; one-way ANOVA pooled variance; LSD (5 %) shown as the grid spacing |
| electrophile assay | 1.1 % t.s. sample diluted 9:1 with 10 mg/mL O-phenylhydroxylamine (Fluka); reagent consumption by HPLC (C18, water / 0.1 % TFA / CH3CN), every 30 min for 3 h, then every 2 h; reported as mmol reacted electrophile / g dry, **assuming 1:1**; curves fitted in TC2D, **no constants printed**. My dose [Z]: coffee 9.9 g/L, reagent 1.0 g/L = 9.2 mM (free base) or 6.9 mM (HCl salt; form not stated), i.e. **0.93 or 0.69 mmol reagent per g coffee**. Temperature and pH not stated |
| oxidizability assay | flow-injection amperometry, glassy carbon vs Ag/AgCl, 200-800 mV in 100 mV steps; 0.01 % w/v sample in degassed 100 mM acetate pH 5, prepared in a glovebox |
| Fig. 5 normalisation | "normalized to the yield of individual fraction" [V]. I read it as equal-mass measurements multiplied by the fraction's yield, i.e. per g of starting coffee. That reading is ~0.75 likely; the alternative is relative-yield dosing (same numbers, same use) |

**Fractionation artefacts that bear on every fraction result [!][Z]:**
- **A1. Fractions are freeze-dried and dried over P2O5.** Freeze-drying is known to change melanoidin retention
  (the authors cite Andriot 2004, ref 19). Recombining the four AEC fractions restores the thiol kinetics (Fig. 4,
  star) but only ~72 % of the electrophile and ~70 % of the oxidizability (Fig. 5).
- **A2. Eluent residues.** AEC fractions come from pyridinium acetate (0.05-1 M) evaporated to "oily residues". SAX-2
  was eluted with 1 M HCl, SCX-2 and Dowex-2 with 0.5 M NH4OH. Residual salt, pyridine or ammonium is not quantified,
  and 10 mM acetate may not hold pH 5.2 against it.
- **A3. H+-form cation exchange.** It removes K+, Mg2+ and transition metals along with the cationic organics, and it
  acidifies the eluate. The loss of thiol reactivity in SCX-1 and Dowex-1 is therefore confounded among three
  causes: "positively charged melanoidins", metal removal and pH. Part 2's DTPA (0.5 mmol/g) did not slow the loss,
  which weakens the metal reading. pH is untested.
- **A4. Two concentration bases.** Figs. 2 and 3 use equal mass (a 1.2 %-yield fraction is dosed at ~80x its natural
  level). Figs. 1, 4 and 6 use relative-yield doses. Rankings across the two bases are not comparable.

---

## 3. Quantitative results

### 3.1 Table 1 (p. 4418) — composition, % w/w [M], re-typed

Footnotes: a after combustion; b by difference; c after hydrolysis (carbohydrates as anhydrosugars, amino acids as
peptides); d used in the correlation study. Lipids are said to be subtracted but have no column.

| sample | yield % | ash | melanoidins | carbohydrates | proteins | caffeine | chlorogenic acids | organic acids |
|---|---|---|---|---|---|---|---|---|
| light brew^d | — | 14.5 | 33.6 | 11.8 | 6.0 | 7.7 | 12.4 | 13.9 |
| medium brew^d | — | 12.8 | 36.9 | 13.9 | 7.0 | 6.6 | 7.3 | 15.4 |
| dark brew^d | — | 12.6 | 43.7 | 15.3 | 5.0 | 5.8 | 1.8 | 15.7 |
| dialysis permeate | 76.5 | 15.3 | 40.2 | 6.1 | 5.3 | 8.0 | 7.8 | 17.3 |
| dialysis retentate | 19.5 | 4.5 | 53.5 | 31.8 | 8.6 | 0.2 | 0.5 | 0.9 |
| SEC-A^d | 9.2 | 3.9 | 36.6 | 51.3 | 6.3 | 0.0 | 0.0 | 1.9 |
| SEC-B | 5.6 | 9.0 | 39.7 | 31.3 | 18.0 | 0.0 | 0.0 | 2.0 |
| SEC-C | 2.9 | 8.7 | 48.2 | 22.0 | 19.0 | 0.0 | 0.1 | 1.9 |
| SEC-D^d | 9.9 | 15.7 | 35.2 | 12.0 | 9.1 | 0.0 | 0.0 | 28.0 |
| SEC-E^d | 11.9 | 17.3 | 38.1 | 4.7 | 3.5 | 0.0 | 0.0 | 26.5 |
| SEC-F^d | 17.3 | 28.1 | 35.7 | 2.8 | 6.3 | 0.0 | 0.0 | 22.8 |
| SEC-G | 20.8 | 7.0 | 49.7 | 2.6 | 3.1 | 25.6 | 7.6 | 3.4 |
| SEC-H | 4.1 | 10.4 | 22.7 | 1.3 | 1.1 | 7.2 | 51.7 | 5.2 |
| SEC-I | 7.6 | 11.5 | 50.2 | 1.5 | 1.0 | 0.4 | 30.1 | 5.3 |
| SEC-J | 1.3 | 6.1 | 75.5 | 3.6 | 1.2 | 0.5 | 11.2 | 1.8 |
| SEC-K | 1.2 | 8.5 | 79.7 | 3.7 | 2.0 | 0.3 | 3.4 | 2.5 |
| AEC-A^d | 77.9 | 12.7 | 34.4 | 11.2 | 6.6 | 7.7 | 7.4 | 20.0 |
| AEC-B^d | 4.1 | 1.4 | 55.2 | 3.4 | 8.9 | 0.4 | 18.1 | 12.5 |
| AEC-C^d | 13.0 | 0.5 | 74.4 | 5.3 | 11.0 | 0 | 1.6 | 7.3 |
| AEC-D^d | 1.5 | 0.5 | 84.8 | 2.2 | 9.2 | 0 | 0.4 | 2.8 |

Checks [Z]: the medium-brew row sums to 99.9 %; the SEC yields sum to 91.8 % (printed "92%"); the AEC yields sum to
96.5 % with Table 1's 4.1 % for AEC-B. **The Methods print AEC-B's yield as "0.1%"** [!], but AEC-B's dose (0.41 g/L)
matches 4.1 %. Eleven samples carry footnote d. The "13 selected samples" are presumably those eleven plus the two
dialysis fractions; the paper does not list them.

### 3.2 Table 2 (p. 4418) — SPE yields, % [M]

Fraction 1 = eluted, fraction 2 = reversibly retained. Ph 57 / 45; CH 60 / 36; C2 59 / 44; C8 65 / 33; C18 52 / 51;
CN 93 / 7; cellulose 98 / 3; polyamide 90 / 15; **Dowex 50Wx8 43 / 17** (60 % recovered, "partially recovered" [Q]);
SCX 79 / 23; SAX 79 / 26.

### 3.3 Table 3 (p. 4423) — Pearson r, composition vs volatile loss at 24 h, 13 samples [M]

| | FFT | EtSH | pyrrole | DMS | 2,3-pentanedione |
|---|---|---|---|---|---|
| overall concentration | 0.39 | 0.80 | 0.95 | 0.97 | 0.88 |
| melanoidins | 0.41 | 0.83 | 0.97 | 0.98 | 0.92 |
| ash | 0.34 | 0.76 | 0.90 | 0.94 | 0.83 |
| protein | 0.41 | 0.83 | 0.92 | 0.95 | 0.81 |
| caffeine | 0.37 | 0.75 | 0.90 | 0.95 | 0.81 |
| organic acids | 0.30 | 0.75 | 0.88 | 0.95 | 0.83 |
| carbohydrates | 0.41 | 0.73 | 0.95 | 0.83 | 0.87 |
| chlorogenic acids | 0.33 | 0.63 | 0.72 | 0.79 | 0.56 |

Text [M/Q]: melanoidins are 34-85 % w/w of the samples and 0.01-0.45 % t.s. in the pot, and correlate with overall
concentration at r = 0.98. Ash, protein, caffeine and organic acids correlate with it at r = 0.97. For EtSH, above
0.20 % melanoidins the 24 h loss is 100 %; below it, loss and melanoidin content correlate at r = 0.96.
**[Z] The thiol rows of this table are censored, not weak.** FFT is lost completely in 10 of 13 samples. A partial
correlation on concentration, which would apportion the classes, was not computed.

### 3.4 Figures (all [F]; my readings)

**Fig. 1 (24 h, relative-yield doses):** caption "FFT = EtSH = 0 in all samples" [M] (light, medium and dark brews;
dialysis HMW and LMW). DMS: light 15, medium 9, dark 7, HMW 88, LMW 20 %. Pyrrole: 66 / 69 / 58 / 88 / 82 %.
2,3-pentanedione: 89 / 89 / 79 / 102 / 93 %. LSD ~12.5 points.

**Fig. 2 (1 h, equal mass 1 % t.s.), % of blank, eluted (fr. 1) / retained (fr. 2):**

| solid phase | EtSH fr. 1 | EtSH fr. 2 | FFT fr. 1 | FFT fr. 2 |
|---|---|---|---|---|
| whole medium brew | 16 | — | 0 | — |
| Ph | 12 | 27 | 0 | 0 |
| CH | 10 | 23 | 0 | 0 |
| C2 | 18 | 23 | 0 | 0 |
| C8 | 10 | 15 | 0 | 0 |
| C18 | 9 | 15 | 0 | 0 |
| CN | 10 | n.d. | 0 | n.d. |
| cellulose | 0 (no bar) | n.d. | 0 | n.d. |
| polyamide | 7 | 34 | 0 | n.d. |
| **SCX** | **57** | no bar [!] | **73** | no bar |
| **Dowex 50Wx8** | **61** | **42** | **71** | no bar |
| **SAX** | 0 (no bar) | 21 | 0 | 0 |

"n.d." = footnote a, "Fraction 2 not determined".
- SCX fraction 2 (both thiols) and Dowex fraction 2 (FFT) have no visible bar and no n.d. footnote. Each is either
  0 % (complete loss) or unplotted [!].
- Polyamide carries the n.d. footnote yet shows a retained EtSH bar (34 %) [!].
- The reversed-phase "0" FFT entries agree with the text ("FFT completely disappeared in all of these samples").
- LSD ~8.3 points.

**Fig. 3 (1 h, equal mass 1 % t.s.), EtSH % of blank:** SEC-A 37, B 38, **C 70**, D 45, E 50, F 30, G 27, H 12,
**I 0, J 0**. SEC-K not shown. Caption: FFT = 0 in all samples [M].

**Fig. 4 (relative-yield doses unless noted; every 2 h from 0.5 h; % of blank), EtSH:**

| matrix (g/L) | 0.5 | 2.5 | 4.5 | 6.5 | 8.5 | 10.5 | 12.5 | 14.5 | 16.5 | 18.5 h |
|---|---|---|---|---|---|---|---|---|---|---|
| medium brew (10) | 64 | 24 | 10 | 2.5 | ~1.5 | ~1 | 0 | 0 | 0 | 0 |
| AEC-ABCD recombined | 72 | 28 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| AEC-A (7.8) | 92 | 33 | 7 | 2 | 2 | ~0.5 | 0 | 0 | 0 | 0 |
| **AEC-A (1.3, = 0.13 %)** | 95 | 79 | 51 | 34 | 23 | 14 | 8 | 6 | 0 | 0 |
| AEC-C (1.3) | 76 | 69 | 62 | 65 | 46 | 40 | 40 | 29 | 25 | 20 |
| AEC-B (0.41) | **106** | 96 | 95 | 91 | 89 | 84 | 79 | 74 | 70 | 69 |
| AEC-D (0.15) | 87 | 83 | 83 | 89 | 83 | 78 | 80 | 78 | 76 | 75 |

**Fig. 4, FFT:** medium brew, AEC-A (7.8) and ABCD: **0 from the first point** (ABCD 5 % at 0.5 h). AEC-A (1.3): 57 %
at 0.5 h, 0 at 2.5 h. AEC-C (1.3): 56, 10, then 0. AEC-B (0.41): **off-scale above 110 % at 0.5 h** [!], then 78, 58,
44, 33, 20, 10, then 0 from 14.5 h. AEC-D (0.15): 75, 60, 43, 38, 24, 17, 16, 8, 5, 3. LSD ~12.5 points.

**My pseudo-first-order fits to these readings [Z]** (ln(%/100) on t; free intercept; through-origin k in brackets;
25 C, pH 5.2):

| series | k (h^-1) | k (s^-1) | intercept | r² | t1/2 |
|---|---|---|---|---|---|
| EtSH, medium brew 10 g/L | **0.53** [0.55] | **1.5e-4** | 90 % | 0.99 | 1.3 h |
| EtSH, ABCD | 0.78 [0.71] | 2.2e-4 | 131 % | 0.95 | 0.9 h |
| EtSH, AEC-A 7.8 g/L | 0.65 [0.58] | 1.8e-4 | 141 % | 0.995 | 1.1 h |
| EtSH, AEC-A 1.3 g/L | 0.21 [0.19] | 5.7e-5 | 123 % | 0.99 | 3.4 h |
| EtSH, AEC-C 1.3 g/L | 0.075 [0.085] | 2.1e-5 | 88 % | 0.96 | 9.3 h |
| EtSH, AEC-B 0.41 g/L | 0.024 [0.019] | 6.6e-6 | 106 % | 0.98 | 29 h |
| EtSH, AEC-D 0.15 g/L | 0.008 [0.018] | ~2e-6 | 88 % | 0.75 | ~80 h (poor) |
| FFT, AEC-B 0.41 g/L (from 2.5 h) | 0.19 [0.16] | 5.4e-5 | 143 % | 0.96 | 3.6 h |
| FFT, AEC-D 0.15 g/L | 0.17 [0.18] | 4.8e-5 | 99 % | 0.97 | 4.0 h |
| FFT, AEC-C 1.3 g/L (2 points) | ~0.86 | ~2.4e-4 | 86 % | — | 0.8 h |
| FFT, medium brew | > ~4-5 (if "0" means < ~6 %; the detection floor is not stated) | > ~1e-3 | — | — | < ~10 min |

What the fits say [Z]:
- **Same-lab replicate.** Whole-brew EtSH 1.5e-4 s^-1 against Part 2 Table 2's 1.19e-4 s^-1: within 25 %. Part 1's
  thiol charge was ~6x higher, and the constant is not lower, so there is no sign of site depletion.
- **Matrix order.** AEC-A at 6.0x the loading gives 3.1x the k (mid-curve two-point slopes give 3.4x), an apparent
  order of ~0.6-0.7 in matrix. A k_obs = k2[E] law with [E] proportional to solids predicts 6x.
- **Lag phase.** AEC-A (both loadings) and ABCD have intercepts of 123-141 %: slower loss in the first ~2 h, then
  faster. The two-point slopes for AEC-A 1.3 g/L are 0.09 h^-1 at 0.5-2.5 h and 0.20-0.25 h^-1 afterwards. This fits
  a reactant that has to form first, for example an oxidation-generated electrophile [P]. The paper does not discuss
  it.
- **Intercepts below 100 %** (brew 90 %, AEC-C and AEC-D 88 %) mean ~10-12 % disappears before the first read. That
  may be fast physical partition (reversible) or fast chemistry; the data cannot separate them. Intercepts **above
  100 %** at t = 0.5 h (AEC-B: EtSH 106 %, FFT > 110 %) show that the headspace ratio also carries a matrix effect on
  partition (salting-out or residual eluent salt), of order ±10 %.
- **FFT / EtSH k ratio:** ~8 (AEC-B), ~10-20 (AEC-D), ~11 (AEC-C).

**Fig. 5, top: O-phenylhydroxylamine consumed (mmol / g dry, normalised to yield):**

| sample | 0.5 h | 1 h | 2 h | 3.5 h | 5.2 h | 8 h | 10-11 h | 13 h |
|---|---|---|---|---|---|---|---|---|
| medium brew | 0.105 | 0.13 | 0.17 | 0.21 | 0.23 | 0.26 | 0.275-0.295 | **0.29** |
| "total calculated" (sum of fractions) | 0.11 | 0.12 | 0.14 | 0.16 | 0.17 | 0.19 | 0.20 | **0.21** |
| AEC-A | 0.08 | 0.07 | 0.10 | 0.11 | 0.11 | 0.13 | 0.13-0.15 | 0.145 |
| AEC-C | 0.026 | 0.03 | 0.035 | 0.04 | 0.040 | — | — | — |
| AEC-B | 0.008 | — | — | 0.01 | — | 0.015 (6.6 h) | — | — |
| AEC-D | ~0.002 | — | — | — | — | — | — | — |

My per-gram-of-fraction values (value / yield) [Z]: AEC-A 0.19 mmol/g at 13 h (0.145 at 5 h); AEC-C 0.31 at 5 h;
AEC-B ~0.37 at 6.6 h; whole brew 0.23 at 5 h and 0.29 at 13 h.

**Fig. 5, bottom: amperometric signal at 800 mV (nA·s, normalised to yield):** brew ~2180, "total calculated" ~1515,
AEC-A ~1035, AEC-C ~330, AEC-B ~110, AEC-D ~40. Per gram of fraction [Z]: AEC-A ~1330; AEC-C ~2540; AEC-B ~2680;
brew 2180.

**Fig. 6 (24 h, relative-yield doses), % of blank:** EtSH: AEC-A 0, AEC-B 68, AEC-C 20, AEC-D 74, SEC-A 74, SEC-D
68, SEC-E 66, SEC-F 33. FFT: AEC-A, B, C 0; AEC-D 3; SEC-A 0; SEC-D 16; SEC-E 13; SEC-F 0. My melanoidin doses
(% t.s.) [Z]: SEC-A 0.033, SEC-D 0.035, SEC-E 0.046, SEC-F 0.061, AEC-B 0.023, AEC-C 0.097, AEC-D 0.013.

---

## 4. Mechanism evidence and its kind

Laboratory: **Nestlé Research Center, Lausanne** (Kraehenbuehl / Roberts group; Rytz statistics).

| claim | kind of evidence | strength |
|---|---|---|
| the coffee matrix (not the buffer) removes thiols | matched matrix-free blank at every time point; blank mix stable 24 h | good (within-study control) |
| the fast thiol sink is in cation-exchange-retained, coloured material | **within-study factor**: fractionation by charge and size, at equal mass and at relative yield | moderate; confounded by pH, metals and eluent residues (A2, A3) |
| "melanoidins are the main" actors | **correlation** over 13 fractions; melanoidins defined **by difference** (everything not otherwise assigned, including colourless low-MW unknowns, which the authors acknowledge) | weak: collinear with concentration, thiol losses censored |
| proteins minor; no thiol / disulfide interchange | authors' judgement on equal r values; DMDS not lost (data not shown) | weak / [Q] |
| carbohydrates inert | no effect of arabinogalactans at 3 g/L; carbohydrate-rich fractions (SEC-A, SEC-D, retentate) low on DMS | moderate for DMS; for thiols SEC-A still kills FFT at 24 h |
| chlorogenic / caffeic acids attack FFT only, needing O2 | pure-component pots (24 h); O2 part cited [C] (ref 5) | moderate (existence); O2 part not shown here |
| thiols "chemically react" (not physical) | **fitted time courses** (pseudo-first-order decay to zero) [Q] | weak for mechanism. Decay to zero excludes a fast physical equilibrium, but not slow reversible addition or disulfide formation |
| the electrophile and oxidizability assays track thiol activity | yield-weighted comparison of four AEC fractions | **not supported on an equal-mass basis** (sec. 6, D3) |
| CROSSPY (pyrazinium) cations as thiol binders | cited from Hofmann & Schieberle (ref 3); "no method to quantify this family" [V] | [C] / [P] |

**No isotope labelling, no fed intermediate, no trapping, no product identified, no mass balance.** By the plan's
evidence ranking ("How a mechanism is identified"), this is rank 6 (one factor varied in one study: the matrix
fraction) plus rank 7 (time courses).

**Reversibility, stated plainly:** the paper cannot tell covalent from disulfide from physical binding.
- Its evidence against physical binding is kinetic only.
- Its evidence against thiol/disulfide interchange is an unshown DMDS null on roasted-coffee solids. The authors
  argue that roasting consumed the protein SH ([P], not measured).
- It measures no disulfide, no released thiol and no bound thiol.
- Cross-lab context: Gigl 2021 (TUM, qNMR, `gigl2021_extraction.md`) finds three things for FFT on the HMW coffee
  melanoidin at room temperature: an instantaneous non-covalent (π-π) drop; difurfuryl disulfide appearing from
  60 min; and only ~5 % (48 h) to 11 % (96 h) covalent incorporation. Part 2 finds the 25 C loss 15-45x slower under
  N2. **Read together, the 25 C coffee thiol "sink" is mostly oxidative (disulfide / radical) plus a physical
  fraction, with covalent addition a minority.** Part 1 neither supports nor contradicts this; it simply cannot see
  it.

---

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| **EX-050 / MFT loss law** (`docs/history/thiol_rivals_2026-10-31/`: A0 rejected, no rival admitted, k_mft_decay to become a declared band λ(110 C) <= ~0.02-0.03 /min) | **LEAVES UNCHANGED** | no MFT, 25 C only, no temperature or O2 factor, no cysteine. Nothing transfers to 95-145 C |
| Triage item 3 (thiol losses anchored on Cerny) | LEAVES UNCHANGED | different system (coffee solids, 25 C) |
| EVIDENCE_MAP 2.12, "MFT loss: NONE at process temperature" | LEAVES UNCHANGED | adds a 25 C FFT / EtSH matrix-sink row, not a process-temperature lifetime |
| EVIDENCE_MAP 1.3 O2 row (charlesbernard2005) | LEAVES UNCHANGED | Part 1's only O2 statement is cited (ref 5); Part 2 holds the air / N2 contrast |
| **`MELE` site density 8-10 mmol/g** (`parameters_sulfur.py` MELANOIDIN_SITE_DENSITY_MMOL_PER_G; Part 2 dossier Flag 1; hofmann2001 capacity >= 0.028 mmol/g) | **SUPPORTS the flag that 8-10 mmol/g over-states the density** | same lab, same extract: a direct consumption assay gives ~0.29 mmol/g at 13 h, flattening at ~31-42 % of a reagent charge of ~0.7-0.9 mmol/g. Caveat: it is a carbonyl probe, not a thiol-reactive-site titre (sec. 6, D1) |
| **order of the sink in matrix** (`k_thioether` order 2: k_obs = k2[MELE], [MELE] proportional to solids) | **weakly CONTRADICTS strict proportionality** | AEC-A 6x loading gives a 3.1x k (apparent order ~0.6) [Z from F]; one pair, lag phase, figure-read. A designed loading series would test it |
| reversible conjugation (Stack 2018 K ~1e2-1e3 M^-1 in the engine) | **tension, LEAVES UNCHANGED** | complete EtSH / FFT loss within 10 h needs K[E] >= ~100. With [E] ≈ 2.9 mM (0.29 mmol/g x 10 g/L) and K <= 1e3, an equilibrium sink would stop at 26-78 % free [Z]. Complete loss therefore needs [E] ~0.1 M (Part 2's number), a larger K, or a mostly irreversible / oxidative channel (Part 2's N2 effect, Gigl's disulfide). Part 1 cannot choose |
| identifying `MELE` with "melanoidins" / electrophile assays | **weakly CONTRADICTS** | on an equal-mass basis, the electrophile and oxidizability assays rank AEC-A last while thiol activity ranks it first (D3) |
| STATUS_AND_PLAN sec. 4 item 7 / triage 9 (protein binding of thiols, e.g. Xu 2022 pea protein > 98 % MFT bound) | **LEAVES UNCHANGED** | roasted-coffee protein is 5-7 % of the solids and denatured; the paper tests no native protein. Nothing here transfers to plant-protein matrices (cf. mottram1996: protein reduces flavour disulfides) |
| "LMW fractions are more thiol-reactive" (Hofmann 2002, cited) and the engine's use of MW > 3000 material | **partly CONTRADICTS a monotonic MW law** | at equal mass, EtSH activity is U-shaped (HMW SEC-A/B active, medium SEC-C inert, phenol-rich SEC-I/J most active); at 24 h both dialysis halves remove all thiol. Per gram, the HMW SEC-A/B fractions have ~0.55x the whole-brew EtSH activity at 1 h [Z: ln(1/0.37) / ln(1/0.16)] |
| CF-006(b) (hofmann2001 / hofmann2002 one paper; k_thioether from it) | LEAVES UNCHANGED | — |
| CF-001 … CF-076 otherwise | **none touched** | no conflict-register entry concerns the coffee-matrix thiol sink. D1 below is a candidate new entry (not added: dossier-only brief) |

---

## 6. Disagreements found

**D1 (NEW) — matrix electrophile density, Part 1 vs Part 2 (same lab, same extract).**
- Quantity: electrophilic sites per g dry medium-roast coffee solids.
- Part 1: **~0.29 mmol/g** O-phenylhydroxylamine consumed at 13 h, slowly rising (Fig. 5) [F].
- Part 2: **8-10 mmol/g**, the hydroxylamine dose above which thiol protection stops increasing [M].
- Ratio ~30. Neither paper mentions the other's number.
- Candidate causes:
  - **Method.** Part 1 consumes a non-saturating reagent charge (~0.7-0.9 mmol/g [Z]), so it cannot read above that
    cap. Part 2 reads a protection plateau, which measures the dose needed, not the sites consumed, and oximation needs
    excess reagent.
  - **Chemistry.** Both probes are carbonyl nucleophiles. Neither titrates Michael acceptors or quinones specifically.
  - **Kinetics.** Part 1 is not at an endpoint.
- Direction: supports the Part 2 dossier's reading that 8-10 mmol/g is an upper bound by one to two orders.
  P(true thiol-reactive density < 1 mmol/g) ~0.7.

**D2 (NEW, within paper) — the 1 h screen vs the kinetic series.**
- Fig. 2 (whole brew, 1 h screen): EtSH **16 %**.
- Fig. 4 (whole brew, kinetic series): 64 % at 0.5 h and 24 % at 2.5 h, i.e. **~45-50 % at 1-1.25 h** by log-linear
  interpolation [Z]. Part 2's constant gives 65 % at 1 h.
- The residual differs ~3x, or the clocks are offset by ~2 h.
- Candidate causes: a different time origin (the screens count 1 h after a 15 min stir, while the series' "0.5 h" may
  count from autosampler start); different vials or handling; a different brew preparation.
- Consequence: **the 1 h screen values (Figs. 2, 3) are rankings within their figure, not points on the kinetic
  clock.**

**D3 (NEW) — "electrophilic reactivity and oxidizability run in parallel with the thiol binding activity" (p. 4422).**
- This holds only on yield-weighted values, which mostly track mass (AEC-A is 78 % of the solids).
- Per gram of fraction [Z]: electrophile AEC-C 0.31 > AEC-A 0.15 (at 5 h); oxidizability AEC-C ~2540 > AEC-A ~1330.
- Yet at equal mass (0.13 %) AEC-A destroys EtSH ~2.8x faster than AEC-C (k 0.21 vs 0.075 h^-1).
- Candidate causes: the assays probe carbonyls and redox-active phenols, not the thiol-reactive species; or my
  reading of "normalized to the yield" (~0.75) is wrong. If the plotted values were measured at relative-yield doses
  instead, the per-gram conversion is the same, so the conclusion stands either way.

**D4 (NEW, minor, within paper).** The text says FFT losses were 100 % "for all samples except three, which present
less than 0.05 % t.s. melanoidins". But SEC-A (0.033 %) and AEC-B (0.023 %) also sit below 0.05 % and still lose 100 %
of FFT (Fig. 6) [Z]. **Melanoidin dose is not a threshold for FFT.**

**D5 (minor, within paper).** The abstract quotes "R = 0.83-0.98, depending on the volatile compound" for melanoidin
and protein. Table 3 gives FFT 0.41 (both) and protein vs pentanedione 0.81.

**D6 (typos and inconsistencies) [!].**
- FFT "4 mg/mL" (read mg/L).
- AEC-B yield "0.1%" in the Methods vs 4.1 % in Table 1.
- Dialysis doses "2.0 g/L permeate / 7.8 g/L retentate" vs yields 76.5 / 19.5 %: probably swapped. Fig. 1's
  "relative yield" caption requires ~7.65 / 1.95.
- The caffeic acid concentration is not printed.
- The kinetic sampling is printed as "every 2 h", while the Methods' "after 1, 24, or every 2 h between 0.5 and
  18.5 h" is consistent with Fig. 4. No conflict.

**Agreement worth recording (not a disagreement).**
- Whole-brew EtSH k: Part 1 1.5e-4 s^-1 [Z from F] vs Part 2 1.19e-4 s^-1 [M]. Same lab, same extract, 25 C, pH 5.2,
  6x different thiol charge: a repeatability check, within 25 %.
- FFT is faster than EtSH in both papers. Part 1's first datapoint (0.5 h) would tighten Part 2's floor (> 7.7e-4
  s^-1) by ~1.3-2x only if the detection floor were known. It is not, so do not tighten it.

**Against the engine.** `k_thioether`'s note says a real-matrix pot should charge `MELE` lower and carry a larger
k2. Part 1 agrees on the density. But its matrix-loading pair (order ~0.6) and its equal-mass ranking say that **no
single site-pool concentration proportional to solids reproduces the fraction data**. The coffee sink is better
treated as an empirical pseudo-first-order term at the tested loading than as a mechanistic bimolecular one.

---

## 7. Ledger-ready measurements

Registry: `2_furfurylthiol` and `chlorogenic_acid` are keyed. **Ethanethiol, dimethyl sulfide, pyrrole and
2,3-pentanedione are not** (`data/keys/compounds.yml`), and there is no melanoidin or coffee-fraction key. Every row:
reconstituted medium-roast coffee (CTN 85) or its fractions, 10 mM acetate pH 5.2 (label), 25 C, air headspace,
SPME-GC-MS area ratio to a matched blank, n >= 2, no IS. **Two readers for every [F] row.**

1. **Fig. 4, EtSH and FFT time courses, 7 matrices + AEC-A at 0.13 %.** The paper's only kinetic data, and the only
   **matrix-loading pair** (AEC-A 7.8 vs 1.3 g/L). Key as `matrix_sink_time_series` at 25 C (not process
   temperature). Flags: lag phase, intercepts ≠ 100 %, LSD 12.5 points, zeros below an unstated floor.
2. **Fig. 5 top, O-phenylhydroxylamine consumption vs time (brew + 4 fractions + sum).** The direct electrophile
   number that bears on `MELE` (D1). Record the reagent charge (0.69-0.93 mmol/g [Z]) next to it. It is a lower bound,
   not a titre.
3. **Fig. 2, SPE eluted / retained, EtSH and FFT at 1 h, equal mass.** The cation-exchange effect (FFT 0 → 71-73 %).
   Flag the pH / metal confound (A3) and D2.
4. **Fig. 3, SEC EtSH at 1 h, equal mass**, with Table 1 compositions. The U-shaped MW dependence.
5. **Table 1 (printed)**: composition of every fraction. Covariates for any reuse of rows 1-4 and 6.
6. **Fig. 6, 24 h EtSH and FFT for 8 fractions at relative-yield doses.** Endpoint data; mostly censored for FFT.
7. *(lower priority)* Fig. 5 bottom, amperometric oxidizability; Fig. 1 DMS / pyrrole / pentanedione (not thiols;
   relevant only if the matrix layer ever carries sulfides or diketones).

**What to request from the authors:**
- (i) the numeric Fig. 4 and Fig. 5 data and the TC2D fit constants;
- (ii) the pH of each reconstituted fraction;
- (iii) the detection floor behind every "0";
- (iv) the DMDS experiment ("data not shown");
- (v) whether the O-phenylhydroxylamine assay was ever run to saturation, which would settle D1.
