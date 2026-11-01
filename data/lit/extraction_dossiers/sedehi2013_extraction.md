# Sedehi, Takano, Blasic, Sullivan & De Haan 2013 — EXTRACTION (glyoxal and methylglyoxal + ammonium sulfate, glycine, serine, arginine, methylamine in D2O, 0.2-1.5 M, pH 2-7, 277-298 K, 1H NMR of reactant losses and imidazole / formic acid formation; pH-dependent second-order rate laws and Ea)

### The only pH- AND temperature-resolved rate laws for alpha-dicarbonyl + amine / ammonium consumption in water on disk. Atmospheric conditions (≤ 25 °C, molar, D2O). The rate constants are referred to the computed dehydrated-aldehyde fraction; only the product k·f_Ald is measured.

**Read 2026-10-31 from `data/articles/sedehi2013.pdf`** (8 pp., Atmos. Environ. 77:656-663; born-digital, clean text
layer). Table 1 and eqs 1-4 checked against 130 dpi renders (scratchpad `dos7/sed-4..6.png`); Fig. 2's lower panel
re-read at 300 dpi (`sedzoom-4.png`) to check the Table 1 fits against the plotted lines (glyoxal + glycine line
through (2.2, -5.2) and (5.7, -2.6): slope 0.74, intercept -6.83, matching Table 1 within read error ±0.1 dex).
**Supplement (Tables S1-S4, Figs. S1-S6, incl. the Arrhenius plots Fig. S5 and per-run Table S2) NOT on disk.**

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| what the model needs | what Sedehi 2013 gives | tag |
|---|---|---|
| rate law for GO / MGO + amine in water | **Rate = k f_Ald [Ald]tot [Am]tot** (eq 1), first order in each reactant below ~1 M; k reported per DEHYDRATED aldehyde, per TOTAL amine | [M] |
| k(pH), glyoxal | log k = slope·pH + intercept (M⁻¹ s⁻¹, room T, Table 1): **Gly 0.741 pH − 6.84; Ser 0.998 pH − 7.59; Arg 0.994 pH − 7.69; MeNH2 0.972 pH − 8.33; ammonium sulfate 1.05 pH − 7.45** (AS fit from Yu 2011 data); "accurate to within a factor ≈2" | [M] |
| k(pH), methylglyoxal | **Gly 0.262 pH − 3.40; Ser 0.421 pH − 3.91; Arg 0.422 pH − 4.09; MeNH2 0.334 pH − 3.22; AS 0.834 pH − 5.91** | [M] |
| unified laws | glyoxal: **Rate = (7 ± 6) × 10¹ M⁻¹ s⁻¹ f_Ald [Glx]tot f_Am [Am]tot** (eq 3, f_Am = deprotonated fraction); methylglyoxal: **Rate = 10^[(0.36 ± 0.06) pH − (3.6 ± 0.3)] M⁻¹ s⁻¹ f_Ald [MG]tot [Am]tot** (eq 4) | [M] |
| Ea | **amino acids 50 ± 11 kJ/mol** (all aldehyde + amino-acid pairs); MG + Arg 55 ± 24, MG + Gly 54 ± 33, MG + Ser 40 ± 9, glyoxal + Ser 51 ± 10; **ammonium sulfate 18 ± 8** (MG + AS 15 ± 5, glyoxal + AS 21 ± 15); 277-298 K only. **No Ea for glyoxal + glycine** | [M] |
| what the loss produces | glyoxal + AS: imidazole + formic acid, ≥ 75 % (pathway 1) at pH 4; glyoxal + amines at pH > 4: mainly N-containing oligomers (pathway 2; at pH 6 imidazole + formate only 18 % of loss for MeNH2); glyoxal + Gly / Ser: ~half the glyoxal loss is self-reaction (pathway 3-like); **MG + amines: mostly amine-catalysed aldol condensation** (MG loss > 20x amino-acid loss, methylimidazole ~1 % yield), except MG + Gly at pH 3.5 (mostly N-methylimidazole) | [M] / [Q] |
| Strecker products (aldehyde, CO2, aminoketone) | **not measured, not discussed** | [!] |
| temperature range | **277-298 K**; nothing above 25 °C | [M] |

**Answers.**

1. **Form.** For glyoxal the consumption is first order in carbonyl and in UNPROTONATED amine (slope ~1 per pH unit
   for Ser, Arg, MeNH2, AS; 0.74 for Gly, where self-reaction dilutes it). For methylglyoxal it is first order in
   TOTAL amine with a weak pH slope (0.26-0.42): the amine (as aminium) catalyses MG aldol condensation, which does not
   need a free amine. This is a structural statement the engine can use: an MG sink that scales with [amine], not a
   first-order constant (`k_mgo_mel`) [Z].
2. **Size at food temperature** (transfer, all assumptions in sec. 5). The measured combination k·f_Ald (f_Ald 7.0 ×
   10⁻⁴ glyoxal, 2.7 × 10⁻³ MG, both COMPUTED [!]) at 25 °C, carried to 100 °C with Ea 50 (54 for MG + Gly), in the
   engine's units L/(mmol·min) [Z]:

   | pair | pH 5 | pH 6 | pH 7 | pH 8 (extrapolated) |
   |---|---|---|---|---|
   | glyoxal + glycine | 1.8e-6 | 9.8e-6 | 5.4e-5 | 3.0e-4 |
   | glyoxal + ammonium (Ea 21) | 1.5e-6 | 1.6e-5 | 1.8e-4 | 2.0e-3 |
   | MG + glycine | 1.1e-4 | 1.9e-4 | 3.5e-4 | 6.4e-4 |
   | MG + ammonium (Ea 15) | 1.0e-5 | 6.8e-5 | 4.6e-4 | 3.2e-3 |

   Ea ±11 kJ/mol moves the 100 °C values by x0.4-2.4 (glyoxal + Gly); to 140 °C multiply the Ea-50 rows by a further
   ~4.8 [Z]. These are TOTAL consumption constants (all pathways), not Strecker constants.
3. **Against the engine's Strecker constants** [Z]: `k_go_ak` (glyoxal + glycine → aminoacetaldehyde, B18) is 2.9e-7
   L/(mmol·min) at 100 °C, pH 8; Sedehi's total glyoxal + glycine consumption transported to the same point is
   ~3e-4, i.e. **~1000x larger**. `k_mgo_ak` 3.0e-8 vs ~6e-4 (~2 × 10⁴ x). Not a contradiction (different observables),
   but it says the Strecker branch is ≤ 0.1 % of glyoxal + amine consumption at 100 °C if both are right, and that
   the amine-driven dicarbonyl sinks the B18 note calls "conditional_on_b13_dicarbonyl_sinks" are of real size: at 0.1
   M glycine, pH 6, 100 °C the glyoxal + Gly consumption is ~1e-3 /min (3 % of the shipped `k_go_sink` 0.0326 /min)
   and MG + Gly ~0.019 /min (t½ ~36 min) [Z].

## 0. Identity

| field | value |
|---|---|
| Title | "Temperature- and pH-dependent aqueous-phase kinetics of the reactions of glyoxal and methylglyoxal with atmospheric amines and ammonium sulfate" |
| Authors | Nahzaneen Sedehi, Hiromi Takano, Vanessa A. Blasic, Kristin A. Sullivan, David O. De Haan* (University of San Diego) |
| Venue | *Atmospheric Environment* **77** (2013) 656-663; received 22 Mar 2013, revised 23 May 2013, accepted 28 May 2013 |
| DOI | **10.1016/j.atmosenv.2013.05.070** (printed; Crossref 2026-10-31: same title, five authors, vol. 77, 656-663). Matches the brief |
| Lab | De Haan group (San Diego); same lab as haan2009 |

## 2. Conditions

| item | as printed |
|---|---|
| solvent | **D2O** (99.9 % D), 1 % v/v acetonitrile as NMR internal standard [M] |
| reactants | glyoxal (from trimer dihydrate, overnight hydrolysis), methylglyoxal (40 % aq.), methylamine (40 % aq.), glycine, serine, arginine, ammonium sulfate [M] |
| concentrations | **0.2-1.5 M** initial (e.g. glyoxal 0.5 M + arginine 0.4 M, Fig. 1) [M] |
| pH | adjusted DOWN with acetic acid-d6, UP with **sodium phosphate**; ISFET sensor; pH 2.5-6 runs tracked through the run by formic / acetic acid NMR shifts and averaged; three runs at pH > 6 use the initial pH only [M]. No pD correction stated in the main text [!] |
| temperature | **277-298 ± 1 K** (Table S2, not on disk) [M] |
| analysis | 1H Wet1D NMR (Varian 400 / 500 MHz), continuous over ≥ 16 h; analysis within ~10 min of mixing; glyoxal "initial" set at t = 6 h after dimer-monomer re-equilibration [M] |
| detection | product (imidazole, formic acid) detection limit ~200 µM, rate constants down to ~10⁻⁶ M⁻¹ s⁻¹ measurable [M] |
| f_Ald | from ab initio ΔG of dehydration (4.3 and 3.5 ± 1 kcal/mol) → f_Ald (7.0 ± 5.7) × 10⁻⁴ (glyoxal), (2.7 ± 2.2) × 10⁻³ (MG) at 298 K [M]. **Computed, not measured** [!] |

## 3. Quantitative results

**Table 1 (log k vs pH, M⁻¹ s⁻¹, room temperature; Ea kJ/mol)** [M]:

| dicarbonyl | reactant | log k(pH) | Ea |
|---|---|---|---|
| MG | arginine | 0.422 pH − 4.09 | 55 ± 24 |
| MG | glycine | 0.262 pH − 3.40 | 54 ± 33 |
| MG | serine | 0.421 pH − 3.91 | 40 ± 9 |
| MG | other amino acids | eq 4 | — |
| MG | methylamine | 0.334 pH − 3.22 | — |
| MG | ammonium sulfate | 0.834 pH − 5.91 | 15 ± 5 |
| MG | OH radical | 9.4 (cited) | 13 ± 6 (cited) |
| glyoxal | arginine | 0.994 pH − 7.69 | — |
| glyoxal | glycine | 0.741 pH − 6.84 | — |
| glyoxal | serine | 0.998 pH − 7.59 | 51 ± 10 |
| glyoxal | methylamine | 0.972 pH − 8.33 | — |
| glyoxal | ammonium sulfate | 1.05 pH − 7.45 (fit to Yu 2011, as 2nd order) | 21 ± 15 |
| glyoxal | OH radical | 9.5 (cited) | 13 ± 1 (cited) |

Footnote c: Ea "calc'd from methylglyoxal loss rates or formic acid production rates" [M]. The pH ranges behind each
fit are only in Figs. 2-3: glyoxal + Gly ~2.2-5.7, + Ser ~2.4-5.4, + Arg ~4-6, + MeNH2 ~4-6, + AS ~1.8-8.3 (incl.
Yu / Nozière data) [F].

**Other printed numbers** [M]: glyoxal + AS formic-acid fit slope 1.05 ± 0.12; imidazole-derived fit (Yu + Nozière)
slope 1.34 ± 0.08; new formic-acid data within x4 of Yu's fit; upper limit 25 % for pathways 2-3 in glyoxal + AS at pH
4; glyoxal + MeNH2 glyoxal-loss slope 0.97 ± 0.06; at pH 6 imidazole / formate production only 18 % of reactant loss
(MeNH2); glyoxal : arginine loss ratio 0.9 ± 0.5, arginine imidazole yield ~0.005; glyoxal + Ser rate constants ~400x
the De Haan 2009a values, glyoxal + Gly ~250x (pH drift); glyoxal + Gly pH slope 0.74 ± 0.08; MG + AS methylimidazole
pathway ≥ 70 %, slope 0.94 ± 0.15; MG loss / amine loss 1.6 ± 0.3 for Arg and Gly at low pH (self-reaction 20-30 %),
> 20x for the other MG + amine runs; MG loss slope 0.36 ± 0.08; MG loss varies < x6 across the four amines.
Kua 2011 computation [C]: C-N bond formation Ea 14.1 kJ/mol, dehydration steps 55 and 61 kJ/mol.

**Atmospheric part (not used)** [M]: aerosol f_Ald 0.046 (glyoxal) and 0.09 (MG); relative contributions at pH 5.5.

## 4. Mechanism evidence and its kind

- **Kind 4-5 (directly measured reactant losses and products, with a stoichiometry comparison):** reactant loss vs
  imidazole + formic acid formation partitions the loss among pathways 1 (imidazole), 2 (imine oligomers), 3 (aldol)
  [M]. **Kind 6:** pH and temperature varied within study, four amines plus ammonium. De Haan group, San Diego.
- Rate-determining step: low Ea with AS read as C-N bond formation; ~50 kJ/mol with amino acids read as later
  dehydration (via Kua 2011 computations [C]) [P].
- MG aldol catalysed by protonated primary amines [P from rate pattern: weak pH slope, no pKa ordering].

## 5. Bearing on open questions — and the transfer to a 100-140 °C food pot

**Transfer assumptions (every one is needed to use these numbers in the engine):**

1. **Use k·f_Ald, never k alone.** f_Ald is computed (ab initio) at 298 K; the measured quantity is the product.
   Multiply Table 1's k by the paper's own f_Ald (7.0e-4 / 2.7e-3) to recover the total-basis constant; this undoes
   the computed number exactly at 298 K and does not import it [Z]. (Repo rule: no computed constants enter.)
2. **Temperature.** Arrhenius from 277-298 K extrapolated to 373-413 K (75-115 K beyond the data). The apparent Ea on
   totals already includes the T-dependence of the hydrate equilibria over 277-298 K; assumed unchanged above
   (dehydration is endothermic, so f_Ald rises with T and the true high-T rate is more likely ABOVE the extrapolation
   than below) [Z]. Also assumed: no change of rate-determining step. Glyoxal + glycine has no own Ea: borrow the
   amino-acid class 50 ± 11.
3. **pH.** Glyoxal + glycine data span pH ~2.2-5.7; pH 6-8 food pots are extrapolations of the slope. pH in D2O read
   on an ISFET without a stated pD correction (pD ≈ reading + 0.4 if uncorrected) → ±0.4 pH units of systematic offset
   [!]. A solvent kinetic isotope effect (D2O vs H2O) is not addressed [!].
4. **Concentration regime.** 0.2-1.5 M here vs µM-mM dicarbonyl and 10-500 mM amine in food. First-order in each
   reactant is claimed for these conditions (two-step mechanism, first step rate-determining unless [Glx][Am] > 1.2
   M²) [M/C]; the low-concentration limit is the better-behaved one for a food pot [Z].
5. **Buffer.** Upward pH adjustment used sodium phosphate; phosphate catalysis was not tested [!] — the phosphate
   question the engine carries (evidence map 1.2) is not answered and may be confounded in the higher-pH runs.
6. **Observable.** These are consumption constants (imidazoles, imine oligomers, aldol oligomers). The model's
   Strecker steps (`k_go_ak`, `k_mgo_ak`) are a sub-branch; Sedehi bounds them from ABOVE only.

**Per question:**
- **Pyrazine / Strecker precursors (evidence map 2.6-2.7; `parameters_pyrazine.py` flag
  "conditional_on_b13_dicarbonyl_sinks"):** SUPPORTS that dicarbonyl sinks with amine are real and pH-shaped; gives a
  form (k[Ald][Am], pH slope 0.74 glyoxal + Gly, 0.26 MG + Gly) and an order of magnitude after transfer. Does not
  measure Strecker yields; LEAVES the Strecker constants UNCHANGED.
- **CF-073 (Strecker aldehyde from MGO + Phe: Weenen ~21 % vs Hofmann ~1 %):** LEAVES UNCHANGED. Sedehi's MG law
  (eq 4) is first order in MG and in amine, as a Strecker step is, so a competing amine-catalysed aldol sink of this
  form predicts a concentration-INDEPENDENT Strecker share — it cannot explain the x14-20 gap between 26 mM and 100 mM
  pots [Z]. It does say that most MG consumed by an amine at pH 4-6 (25 °C) goes to aldol products, so a Strecker
  share of a few % (Hofmann) is the less surprising of the two values — weak, transfer-laden.
- **Engine `k_mgo_mel` (MG → melanoidin carbon, first order, amine-independent, fitted):** CONTRADICTS the form at
  25 °C — MG consumption is first order in total amine, weakly pH-dependent [M]. Whether the form holds at 100-140 °C is
  untested.
- **Engine `k_go_sink` (0.0326 /min, barrier fixed to zero, DISPUTED_SINK_KEYS):** an amine-dependent glyoxal sink of
  ~1e-3 /min at 0.1 M Gly, pH 6, 100 °C [Z] covers ~3 % of it; at pH 8 (extrapolated) ~3e-2 /min, comparable. A
  glyoxal sink with an amine and pH term is a candidate structure for the disputed sink.
- **Thiol losses / observation operators:** not touched.

## 6. Disagreements found

1. **Same lab, x250 (glycine) and x400 (serine): Sedehi 2013 vs haan2009 Table 1** — the 2009 constants assumed a
   constant pH 6; formic acid production dropped the pH by 1-2 units. Explained by an unmonitored pH (operator /
   conditions). Use Sedehi; do not use haan2009's k [M].
2. **Glyoxal + AS: imidazole-based vs formic-acid-based pH slopes** (1.34 vs 1.05) and the Nozière 2009 UV-209 nm
   data (later attributed to imidazole and imidazole-2-carboxaldehyde, Yu 2011 [C]) lie on the imidazole fit;
   Nozière's own rate law (second order in glyoxal, exp(2.5 pH)) is a different functional form (see
   noziere2009_extraction.md sec. 6).
3. **Engine Strecker constants (B18, Zhou 2024 fit) vs total amine-driven consumption** — see sec. 1 answer 3 (~10³
   for glyoxal); not a contradiction, a bound.

## 7. Ledger-ready measurements

- **Table 1 rows glyoxal + Gly, MG + Gly, glyoxal + AS, MG + AS** (two readers) with their Ea where printed, keyed as
  k·f_Ald (total basis) at 298 K with the pH law, flagged "atmospheric, ≤ 25 °C, D2O, 0.2-1.5 M".
- **Eqs 3 and 4** (the aggregated laws) — candidate priors for an amine-dependent dicarbonyl sink, as a declared band
  after the transfer (never as a fit row: no Strecker observable).
- **Table S2** (per-run k, T, pH) would be the real keying target: SI not on disk; fetch if this lane is used.
