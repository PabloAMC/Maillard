# Charles-Bernard 2005 — EXTRACTION (reconstituted coffee brew model, 1 % total solids in 0.01 M acetate pH 5.2, 25 C, air vs nitrogen; ten volatile thiols followed by SPME-GC-MS over ~15 h; pseudo-first-order k_obs per thiol plus a seven-additive inhibition screen reported as k_rel)

### The paper DOES print "8-10 mmol/g dry coffee" and calls it the number of electrophilic sites of the matrix — but it is a hydroxylamine dose at which a thiol-protection curve levels off, read on whole 1 % coffee solids at pH 5.2, not a titre on the MW > 3000 melanoidin fraction the repository applies it to, and it sits ABOVE the stoichiometric ceiling for any carbonyl-per-sugar-residue count, so the recast that produces `k_thioether` = 5.01e-4 L mmol^-1 min^-1 reproduces arithmetically but rests on a site density the paper cannot support as a stoichiometric concentration.

**Source on disk:** `data/articles/charles-bernard2005.pdf` (8 pp., J. Agric. Food Chem. 2005, 53 (11), 4426-4433).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/charles-bernard2005.txt`, 453 lines). **Tables 1, 2 and 3 came through clean**
and are re-typed in full below; the only damage is cosmetic (Table 2's footnote "In acetone/H2O
(3:1 v/v)" is broken across two lines by the subscript 2, and Table 3's columns interleave with
Figure 7's caption). Figures 1-11 are images. **Figure 3 is the dose-response curve from which the
8-10 mmol/g number is read, and it is figure-only**: the number reaches the page only as a sentence
in the Results and again in the Discussion. There is no supplementary material. Repo status before
this dossier: this paper is cited by `src/kinetic_core/parameters_sulfur.py` (the `k_thioether`
`source_anchor`, the `site_density_midpoint_9_mmol_per_g` flag and the
`charles_bernard_table2_unit_erratum_values_are_per_second` flag), by
`src/kinetic_core/species_sulfur.py` (the `MELE` species note), and by
`k3_final_parameter_inventory.md` secs. A.4 and F — but it has **no extraction dossier**, so the
one number the whole bimolecular recast turns on has never been checked against the page.

## 0. Identity

| field | value |
|---|---|
| Title | "Interactions between Volatile and Nonvolatile Coffee Components. 2. Mechanistic Study Focused on Volatile Thiols" |
| Authors | Marielle Charles-Bernard, Deborah D. Roberts, Karin Kraehenbuehl (corresponding) — Nestlé Research Center, P.O. Box 44, Vers-Chez-les-Blanc, CH-1000 Lausanne 26, Switzerland |
| Venue | J. Agric. Food Chem. 2005, 53 (11), 4426-4433. Received 26 November 2004, revised 17 March 2005, accepted 20 March 2005, web 26 April 2005 |
| DOI / article ID | 10.1021/jf048020y (printed as `JF048020Y`) |
| Part 1 of the pair | Bernard, M.; Kraehenbuehl, K.; Roberts, D. D. "Interactions between volatile and nonvolatile coffee components. 1. Screening of the nonvolatile components." J. Agric. Food Chem. 2005, 53, 4417-4425 — **not on disk**. Part 1 is where "melanoidins are the main responsible components for thiols degradation" is established; **this paper does no fractionation at all** |
| Naming | `k(obs)` = pseudo-first-order observed rate constant; `k(ref)` = k(obs) in coffee without additive; `k(rel)` = k(obs)/k(ref); Table 3 reports **1/k(rel)**, called "rate attenuation". "t.s." = total solid |
| Lineage | thiol-melanoidin covalent binding from Hofmann & Schieberle 2002 (ref 4, the CROSSPY pyrazinium work); quinone [1,4]-addition from Cilliers & Singleton 1989/1990 (refs 12, 13); relative thiolate nucleophilicities from Dmuchovsky 1966 (ref 23) |
| Companions on disk | `hofmann2002_extraction.md` (the 30 C model systems and the 80 C brew hold-out this paper corroborates) |

## 1. Why it matters

`src/kinetic_core/species_sulfur.py` declares `MELE`, the matrix electrophile site pool, as a
zero-atom `role='site'` species and justifies its existence in one sentence: it is

> "a titrated site density (Charles-Bernard 2005: 8-10 mmol per g dry coffee solids, by
> hydroxylamine titration)"

and the same file gives the reason the pool exists at all — "**both are DEPLETABLE, which is what
makes the thiol sink a bimolecular, saturating channel rather than an eternal first-order drain**."

`src/kinetic_core/parameters_sulfur.py` then uses that density to recast the corpus's
pseudo-first-order thiol-loss constant into `k_thioether`, order 2, `k_ref` = 5.01e-4, with the
`source_anchor` reading "recast bimolecular on Charles-Bernard's titrated site density of 8-10
mmol per g dry solids (p. 4428)" and the flags `bimolecular_recast_from_pseudo_first_order` and
`site_density_midpoint_9_mmol_per_g`. `THIOL_CHANNELS` describes the channel as "second (thiol x
electrophile SITE)".

**So one number from one sentence of this paper carries the entire order-2 structure of the sulfur
lane's dominant 25-30 C thiol sink.** This dossier's job is to say whether the page supports it.
The short answers, expanded in section 4 and Flags 1-5:

1. **Is the density printed?** Yes, twice, in the running text — p. 4428 ("This experiment is an
   indirect method to measure the number of electrophilic sites of the matrix (in this case 8-10
   mmol/g dry coffee)") and p. 4430 ("dramatically stabilized when 8-10 mmol of hydroxylamine/g dry
   coffee was added"). The repository's page cite is correct.
2. **By what method?** **Not a titration.** It is the hydroxylamine dose above which the protective
   effect on aliphatic thiols stops increasing, read off Figure 3, which is an image. The authors
   themselves call it "an indirect method". Only four hydroxylamine doses were run — 3.6, 7.2, 14.4
   and 21.6 mmol/g (Table 1) — so **there is no experimental point anywhere in 8-10 mmol/g**; the
   value is an interpolation between the second and third dose.
3. **On what material?** A reconstituted brew of a whole medium-roast coffee extract (Arabica 80 %
   / Robusta 20 %, CTN 85, extraction yield 22.7 %) at **1 % total solids in 0.01 M acetate,
   pH 5.2, 25 C**. Whole soluble solids, no fractionation. The paper states in the same paragraph
   that "these electrophilic sites include also reducing sugars, which present moderate to low
   reactivity toward thiols."
4. **Is the recast what the paper supports?** The arithmetic is faithful and I reproduce it exactly
   (section 4). The chemistry is not: the density is applied to Hofmann's **isolated MW > 3000
   melanoidin fraction at 12.5 g/L**, a material from which every reducing sugar has been removed
   by the cutoff, and the resulting site concentration (112.5 mmol/L) exceeds the thiol
   concentration by a factor of ~250 in Hofmann's pot and ~2000 in this one — so **the "depletable"
   pool depletes by at most 0.4 % in any system in the corpus and the order-2 form is numerically
   indistinguishable from the pseudo-first-order form it replaced.** The recast buys no behaviour;
   it buys a stated mechanism.

The paper's second contribution to the lane is corroborative and is real: **Table 2 is an
independent laboratory's pseudo-first-order thiol-loss constant in a real coffee matrix at 25 C**,
and Table 3 is a mechanistic partition (nucleophilic addition vs oxygen/radical) that no other
paper in the corpus supplies. `MELE`'s note claims the pool is "a LUMP OVER AT LEAST TWO MEASURED,
CHEMICALLY DISTINCT CHANNELS"; this paper is the source of the claim that the addition happens on
**oxidised** species and needs air, which is what makes `OX` and `MELE` co-required rather than
independent.

What this paper does NOT give: any temperature series (so **no activation energy** — the `ea` on
`k_thioether` comes from Stack 2018, not from here), any pH series, any absolute concentration-time
number (all figure-only), any identification of the electrophile, any measurement on
**2-methyl-3-furanthiol** or **methanethiol** (neither is in the thiol mixture, despite methanethiol
being named in the abstract), and any fractionation of the matrix.

## 2. Methods as they matter to a model

- **Matrix.** Coffee brew models from a blend of **Arabica 80 % + Robusta 20 %**; the medium-roasted
  extract **CTN 85** with **extraction yield 22.7 %**. Coffee stock solutions at **2.5 % total
  solid in 0.01 M acetate buffer, pH 5.2**. Preparation "as described earlier (7)" — i.e. in Part 1,
  which is not on disk.
- **The pot the kinetics were run in.** Coffee stock (2.5 % t.s.) diluted 4:1 with additive stock,
  stirred 1 h at 25 C, then diluted **1:1** with the aroma mixture. **Final coffee = 1 % t.s.**,
  i.e. **10 g dry coffee solids per L**. 800 uL into a 2 mL amber silane-treated glass vial,
  equilibrated 1 h at 25 C before headspace analysis. Headspace-to-liquid ratio 1.2 mL : 0.8 mL.
- **Thiols, final concentrations at t0 [umol/L]** (printed in Methods): 2-methyl-2-propanethiol
  (2M2P) **3.7**; 3-mercapto-3-methylbutyl formate (MMBF) **4.9**; 2-butanethiol (2BT) **3.7**;
  ethanethiol (EtSH) **5.4**; propanethiol (PropSH) **4.4**; butanethiol (BuSH) **3.7**;
  pentanethiol (PentSH) **3.2**; 2-furfurylthiol (FFT) **3.9**; benzylthiol (BnSH) **3.4**;
  thiophenol (PhSH) **9.7**. Sum = **46.0 umol/L (mine)**. Thiol stock prepared in a nitrogen
  glovebox at twice final concentration, then diluted 1:1. Chosen "to be in the linear range of the
  SPME fiber and at a ppm level."
- **The measured quantity is a RATIO, not a concentration.** For every time point a blank vial (A:
  thiols + buffer + additive, no coffee) and a coffee vial (B) were prepared at time zero, run
  alternately on the autosampler, and "the ratio of the two integration surfaces was plotted as a
  function of additive concentration and/or time". 100 % in every figure is the matched blank. So
  **k(obs) is a matrix-attributable loss constant, already net of fibre drift, partitioning and any
  loss in buffer.** That is a good property and should be said when the number is compared with
  Hofmann's, which is a direct concentration decay.
- **Kinetic treatment, printed verbatim.** "The data were treated assuming pseudo first-order
  kinetics. For each volatile compound, the ln of concentration was expressed as a function of time
  **[s]**. The slope of the curve gave -k(obs), the observed rate constant." **The time base is
  seconds.** This is the printed sentence that settles the unit erratum (Flags 2).
- **Sampling cadence.** "vials containing blank samples A and coffee samples B were prepared at time
  zero and put alternatively on the autosampler so that the headspace was sampled every 2-4 h in
  intact 'aged' vials." Kinetics therefore run over roughly 15 h (Figures 2 and 5 both report a
  15 h endpoint). **Aromatic and benzylic thiols "were already undetectable at the first
  datapoint"** — so their constants are bounds, not fits (Flags 3).
- **Analysis.** SPME-GC-MS. Varian CP-820 autosampler; Hewlett-Packard 5973; DB-Wax 30 m x
  0.25 mm i.d. x 0.25 um film, 0.9 mL/min constant flow. **PDMS/DVB fibre 65 um, 1 min
  equilibration**; desorption 5 min at 240 C, unsplit for the last 2 min. Oven 35 C (3 min),
  35-170 C at 4 C/min, 170-220 C at 20 C/min, 220 C for 10 min. Scan mode, 29-300 amu. **No
  internal standard and no isotope dilution** (contrast Hofmann 2002's SIDA).
- **Anaerobic trials.** "entirely prepared in a glovebox (Easy Box EB 80-1 spez., MecaPlex) and with
  previously degassed buffer (**Ar bubbling for 1 h**)." Note the glovebox atmosphere is nitrogen
  and the degassing gas is argon; the paper labels the condition "under N2".
- **Blank stability control.** "At room temperature (RT), the mixture containing the nine thiols in
  the working buffer with air as headspace was shown to be stable over 24 h." (Ten thiols are
  listed and Table 2 has ten rows — Flags 8.)
- **Additive referencing is not uniform.** Na2SO3 and Na2S2O3 "both decreased the headspace
  concentration of the blank thiol mixture even upon short time equilibration (1 h). The results
  with these additives are therefore expressed **relative to a blank thiol mixture solution without
  additive**." All other additives are referenced to blank + additive. **The two reducing-agent rows
  of Table 3 therefore have a different denominator from every other row** (Flags 6).
- **Temperature.** 25 C throughout, with some steps at unspecified "RT". **One temperature only.**
- **pH.** 5.2, set by 0.01 M acetate, additive stocks re-adjusted to 5.2 with HCl or NaOH. Never
  re-measured; the Discussion refers to "coffee pH (~5)" and "the coffee beverage (pH ~ 5)".
- **Replication.** "Each data point was measured in duplicate, and the error bars are deviations
  from average." **n = 2.** No confidence intervals on any rate constant anywhere in the paper.

## 3. Tables re-typed

### Table 1 (p. 4427). "Preparation of Coffee Brews with Additives"

Footnote a: "Concentrations of stock solutions and final concentrations in sample for interaction
measurement."

| additive | expected action | additive stock solutions [mg/mL] | additive final concentrations [mmol/g dry coffee] |
|---|---|---|---|
| hydroxylamine x HCl | nucleophilic competitor | 25, 50, 100, 150 | **3.6, 7.2, 14.4, 21.6** |
| ascorbic acid (Na salt) | radical scavenger | 2, 5, 10, 20, 50 | 0.1, 0.5, 1, 2, 5 |
| caffeic acid | oxygen scavenger | 1, 2, 5 | 0.05, 0.11, 0.28 |
| DTPA | metal chelator | 3.9-19.6 | 0.1, 0.5 |
| Na2S2O3·5H2O | reducing agent | 2.5, 25, 125 | 0.1, 1, 5 |
| Na2SO3 | nucleophile and reducing agent | 1.3, 6.3, 12.6 | 0.1, 0.5, 1 |

**This is the table that constrains the site-density claim.** The four hydroxylamine doses are
3.6, 7.2, 14.4 and 21.6 mmol/g. The claimed level-off at "8-10 mmol/g" falls **between the second
and third dose, with no measurement in it.**

### Table 2 (p. 4430). "Rate Constants of Thiol Degradation in the Presence of Reconstituted Coffee Brew, t.s. 1 %, pH 5.2; See Also Figure 2"

Header exactly as printed: **`rate constant kref [mol-1 s-1]`** (on this unit see Flags 2 — the
values are per second). "relative nucleophilicity" column is credited to reference (24) in the
table but the text attributes the quantity to Dmuchovsky, reference (23) (Flags 9). Footnote a:
"In acetone/H2O (3:1 v/v)". Footnote b: "In water (23)."

| flavor compounds | type | rate constant kref [mol-1 s-1] *(as printed)* | pKa | relative nucleophilicity (24) |
|---|---|---|---|---|
| PhSH/PhS- | aromatic | > 7.70 x 10^-04 | 8.6,^a 6.5^b | — |
| FFT | primary/benzylic | > 7.70 x 10^-04 | 11.3^a | 1400 |
| BenzSH | primary/benzylic | > 7.70 x 10^-04 | 11.8^a | 1200 |
| PentSH | primary | 1.83 x 10^-04 | — | — |
| BuSH | primary | 1.42 x 10^-04 | 12.6^a | 1000 |
| PropSH | primary | 1.32 x 10^-04 | — | — |
| EtSH | primary | 1.19 x 10^-04 | — | — |
| 2BT | secondary | 8.11 x 10^-05 | 12.9^a | 380 |
| MMBF | tertiary | 2.12 x 10^-05 | — | — |
| 2M2P | tertiary | 1.02 x 10^-06 | 13.1^a | 340 |

### Table 3 (p. 4430). "Rate Attenuation (1/kRel) of Thiol Losses Observed in the Presence of Various Additives"

Footnote a: "Relative to aerobic reference without hydroxylamine." Footnote b: "Relative to
anaerobic without hydroxylamine."

| additive | amount (mmol/g) | EtSH | PropSH | BuSH | PentSH | FFT |
|---|---|---|---|---|---|---|
| aerobic | — | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| anaerobic | — | 17.9 | 19.8 | 15.9 | 16.4 | 45.1 |
| hydroxylamine | 10 | 10.5 | 11.4 | 10.8 | 13.4 | 4.0 |
| anaerobic + hydroxylamine^a | 21 | 64.8 | 82.5 | 101.3 | 236.1 | 92.6 |
| anaerobic + hydroxylamine^b | 21 | 3.6 | 4.2 | 6.4 | 14.4 | 4.1 |
| Na2SO3 | 1 | 51.5 | 61.7 | 18.4 | 62.7 | 67.9 |
| ascorbic acid | 1 | 9.0 | 7.8 | 4.4 | 3.7 | 23.6 |
| caffeic acid | 0.28 | 0.3 | 0.3 | 0.3 | 0.3 | 1.0 |
| DTPA | 0.5 | 0.5 | 0.7 | 0.7 | 1.0 | 1.0 |

**The "10 mmol/g" in the hydroxylamine row matches none of Table 1's four doses (3.6, 7.2, 14.4,
21.6)**; the "21" of the combined rows is evidently 21.6 rounded. Flags 7.

### Numbers printed in the running text (everything else in this paper is figure-only)

| quantity | value | where |
|---|---|---|
| **electrophilic site density of the matrix** | **8-10 mmol/g dry coffee** ("This experiment is an indirect method to measure the number of electrophilic sites of the matrix (in this case 8-10 mmol/g dry coffee)") | Results, p. 4428 |
| the same, restated | "the aliphatic thiols (e.g., ethanethiol) were dramatically stabilized when **8-10 mmol of hydroxylamine/g dry coffee** was added 1 h prior to the contact between thiols and matrix" | Discussion, p. 4430 |
| hydroxylamine level-off, aliphatic thiols | "for all primary, secondary, and tertiary thiols it leveled off **above 8-10 mmol hydroxylamine/g coffee solids**" | Results, p. 4428 |
| hydroxylamine, FFT exception | "2-Furfurylthiol, a benzylic thiol, is further stabilized even above this concentration"; "still degraded even in the presence of **7-21 mmol** of hydroxylamine/g dry coffee" | Results p. 4428, Discussion p. 4430 |
| hydroxylamine preincubation | "almost completed after only 1 h"; extending 1 h -> 24 h did not significantly increase stabilisation except for FFT | Results, p. 4428 (Fig. 4) |
| anaerobic residual loss | "In the absence of oxygen, linear aliphatic thiols are only **20-30 % degraded after 15 h**" | Discussion, p. 4429-4430 |
| anaerobic attenuation, aliphatic | "a **15-20 times** attenuation factor of the reaction rate as compared to the reference system under air" | Discussion, p. 4430 |
| anaerobic attenuation, FFT | "the degradation of FFT is slowed **more than 40 times**" | Discussion, p. 4430 |
| ascorbic acid threshold | "The presence of **> 1 mmol ascorbic acid/g dry coffee**... strongly decreased the magnitude of the interactions"; "leveled off above **0.5 mmol/g dry coffee**" | Results, p. 4429 |
| the synergy numbers quoted in text | "for FFT 1/krel = **45, 4, and 93**; EtSH 1/krel = **18, 11, and 65**" (N2; hydroxylamine; N2 + hydroxylamine) | Discussion, p. 4431 |
| benzylic thiols at t = 1 h | FFT and benzylthiol "were undetectable after only 1 h exposure to the coffee matrix"; "Aromatic and benzylic thiols were already undetectable at the first datapoint" | Results p. 4428, Discussion p. 4430 |
| blank stability | thiol mixture in buffer under air stable over **24 h** at RT | Methods, p. 4428 |
| Na2S2O3 | "the thiol mixture was unaffected or slightly destabilized depending on concentration" | Discussion, p. 4430 |
| extraction yield of the extract used | **22.7 %** (CTN 85, medium roast) | Methods, p. 4427 |

**Figure-only in this paper (not typed as numbers).** Fig. 1 (pathway scheme); **Fig. 2**
(degradation kinetics of the thiols in 1 % brew — the raw data behind Table 2); **Fig. 3 (the
hydroxylamine dose-response from which 8-10 mmol/g is read)**; Fig. 4 (1 h vs 24 h preincubation);
Fig. 5 (air vs N2 at 15 h); Fig. 6 (Na2SO3 and Na2S2O3 at 1 h); Fig. 7 (ascorbic acid); Fig. 8
(caffeic acid and DTPA); Fig. 9 (NH2OH 21 mmol/g and N2, separate and combined, at 1 h); Fig. 10
(proposed mechanism); **Fig. 11 (rate constant vs log P for the linear aliphatic thiols — the
"direct correlation" claimed in the abstract has no printed slope, intercept or r^2 anywhere)**.

### Arithmetic on the printed numbers (all mine)

**1. The site concentration this paper actually measures.** 1 % t.s. = **10 g dry solids/L**, so
8-10 mmol/g gives **[E] = 80-100 mmol/L in the pot, midpoint 90 mmol/L = 0.090 mol/L**. Against a
total thiol charge of 46.0 umol/L, the sites outnumber the thiols by **~2000-fold**. Nothing in this
paper's own pot could ever be limited by site depletion.

**2. Reproducing the repository's recast, exactly.** Hofmann's model system is 12.5 g/L melanoidin
(MW > 3000). At the 9 mmol/g midpoint, [E] = 12.5 x 9 = **112.5 mmol/L = 0.1125 mol/L**. Hofmann's
Fig. 6 pseudo-first-order constant 9.4e-4 s^-1 divided by 0.1125 mol/L = **8.36e-3 L mol^-1 s^-1**,
and 8.36e-3 x 60 / 1000 = **5.01e-4 L mmol^-1 min^-1** — which is `k_thioether`'s `k_ref` to three
figures. **The registry's arithmetic is faithful; I could not find an error in it.** (Using Table
2's 9.8e-4 s^-1 instead gives 8.71e-3 L mol^-1 s^-1 = 5.23e-4 L mmol^-1 min^-1, which is where the
"8.6e-3 L/(mol*s)" in the `QUINONE_THIOL_K2_CONTEXT` comment sits — between the two.)

**3. The pool cannot deplete.** In Hofmann's system FFT is 438 umol/L against [E] = 112.5 mmol/L, so
1:1 thioether addition consumes **0.39 %** of the pool at complete thiol conversion. In this paper's
pot, 46.0 umol/L against 90 mmol/L is **0.051 %**. `MELE`'s stated purpose — "both are DEPLETABLE,
which is what makes the thiol sink a bimolecular, saturating channel rather than an eternal
first-order drain" — is **not delivered at this site density**: to within 0.4 %, `k_thioether * [E]`
is a constant and the channel is arithmetically identical to the pseudo-first-order drain it
replaced. The recast changes the lane's *statement of mechanism*, not its numbers. That is worth
having, but it must not be reported as a modelling improvement in the outputs.

**4. The stoichiometric ceiling, and 8-10 mmol/g is above it.** A site density of 9 mmol/g requires
a mean mass per electrophilic site of 1/0.009 = **111 g/mol**. For comparison: a hexose residue in a
polysaccharide is 162 g/mol (ceiling **6.2 mmol/g** if every single residue carried a free
carbonyl); free glucose is 180 g/mol (**5.6 mmol/g**); 5-hydroxymethylfurfural is 126 g/mol
(**7.9 mmol/g**); chlorogenic acid is 354 g/mol (**2.8 mmol/g**). Coffee brew solids are dominated
by polysaccharides, melanoidins, chlorogenic acids and caffeine, none of which can supply one
thiol-reactive electrophile per 111 g. **8-10 mmol/g therefore cannot be a stoichiometric count of
electrophilic sites; it is a saturating dose of a reagent**, and hydroxylamine oximation is an
equilibrium that needs large excess to run to completion, quite apart from the hydroxylamine
consumed by everything in the matrix that is not a thiol-reactive site. My best reading is that the
true thiol-reactive site density is **one to two orders of magnitude below** the quoted figure, and
therefore that the true bimolecular constant is **one to two orders of magnitude above**
5.01e-4 L mmol^-1 min^-1. I state that as a bound and a direction, not as a number (Flags 1). It is
consistent with the direction of `parameters_sulfur.py`'s own carried-but-not-operative note, which
records that at literature quinone-thiol constants of 5e5-7e5 M^-1 s^-1 "a NANOMOLAR electrophile
pool already reproduces the coffee-matrix loss rate."

**5. The unit erratum, proven three ways from this paper alone (mine).** (i) The Methods sentence
prints the time base: ln(concentration) against time **in seconds**, slope = -k(obs) — a
first-order construction, so the unit is s^-1. (ii) Take EtSH under nitrogen: k = 1.19e-4 / 17.9
(Table 3) = 6.65e-6 s^-1; over 15 h = 54000 s that is kt = 0.359, leaving 69.8 %, i.e. **30.2 %
degraded**, against the printed "linear aliphatic thiols are only 20-30 % degraded after 15 h."
PropSH gives 30.2 % by the same route, BuSH 38 %, PentSH 45 %. (iii) Read as a genuine second-order
constant, 1.19e-4 mol^-1 s^-1 acting on a 4-46 umol/L thiol pool would produce a loss of order
1e-9 per second — nothing would happen in 15 h and Figure 2 would be flat. **The values are s^-1.**
The registry's flag `charles_bernard_table2_unit_erratum_values_are_per_second` is confirmed
independently of the inventory's half-life argument.

**6. Where "> 7.70e-4" comes from (mine).** 4 x ln2 / 3600 s = **7.702e-4 s^-1** exactly. The three
undetectable thiols (PhSH, FFT, BenzSH) were gone by the first 1 h datapoint, and the printed bound
is precisely the constant that puts four half-lives in that hour, i.e. that leaves 1/16 = 6.25 % of
the initial. **The bound is a detection-floor construction, not a fit**, and it is the same number
for all three thiols because it depends only on the sampling time. This matters: FFT's constant here
is a *floor set by the experiment's time resolution*, so the repository's use of it as "the
corroborating lower bound" is exactly the right reading and it must never be read as an estimate.

**7. Half-lives from Table 2 (mine, first order, at 25 C, 1 % t.s., pH 5.2).**
PhSH / FFT / BenzSH **< 15.0 min**; PentSH **63.1 min**; BuSH **81.4 min**; PropSH **87.5 min**;
EtSH **97.1 min**; 2BT **142 min**; MMBF **545 min (9.1 h)**; 2M2P **189 h**. The full spread from
the fastest bound to the slowest is **> 755-fold**, and within the primary aliphatics alone (EtSH to
PentSH) it is **1.54-fold**.

**8. Table 3's two anaerobic+hydroxylamine rows are consistent for four thiols and inconsistent for
FFT (mine).** Footnote b's value times the anaerobic value should reproduce footnote a's value:
EtSH 17.9 x 3.6 = 64.4 vs 64.8; PropSH 19.8 x 4.2 = 83.2 vs 82.5; BuSH 15.9 x 6.4 = 101.8 vs 101.3;
PentSH 16.4 x 14.4 = 236.2 vs 236.1. All four close to within 1 %. **FFT: 45.1 x 4.1 = 184.9 against
a printed 92.6 — off by a factor of 2.05.** Either FFT's footnote-b entry should be 2.05, or its
footnote-a entry should be ~185. It is suspicious that the printed 4.1 is within rounding of the
aerobic hydroxylamine value 4.0 in the row above. **FFT is the repository's thiol, so this is the
one row where the inconsistency bites** (Flags 5).

**9. The two stabilisations are sub-multiplicative (mine).** For EtSH the separate attenuations are
17.9 (no O2) and 10.5 (hydroxylamine), whose product is 188, against a measured combined 64.8 — the
combination recovers **34 %** of the product. If the FFT footnote-a value is taken at face value:
45.1 x 4.0 = 180 against 92.6, i.e. **51 %**. The paper describes this as "an addition of the two
stabilization effects"; it is neither additive (28.4 and 49.1 predicted) nor multiplicative.
**Read: the nucleophilic-addition and oxidation channels overlap — removing oxygen removes part of
the electrophile pool, because the electrophiles are themselves oxidation products.** That is the
paper's own conclusion and it is the mechanistic warrant for `MELE` and `OX` being co-required in
`species_sulfur.py` rather than independent multipliers.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Of the ten thiols in this study, exactly
**one** is keyed: `2_furfurylthiol` (FFT). `methanethiol` and `2_methyl_3_furanthiol` are keyed but
are **not studied in this paper**. **Ethanethiol, propanethiol, butanethiol, pentanethiol,
2-butanethiol, 2-methyl-2-propanethiol, 3-mercapto-3-methylbutyl formate, benzylthiol and
thiophenol are all absent from the registry**, as are hydroxylamine, ascorbic acid, caffeic acid,
DTPA and "melanoidins". Every row below shares: reconstituted coffee brew, **1 % total solids
(10 g/L) in 0.01 M acetate, pH 5.2, 25 C, air headspace in a sealed 2 mL vial (0.8 mL liquid),
unstirred during measurement, duplicate, ~15 h**, and is expressed against a matched
coffee-free blank.

| quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **matrix electrophile site density** | **8-10** (midpoint 9) | **mmol per g dry coffee solids** | whole 1 % t.s. brew, pH 5.2, 25 C, hydroxylamine preincubated 1 h, aerobic | — (a saturating dose, not a rate) | Results p. 4428, restated Discussion p. 4430 | **level_only** — the sentence is printed, the curve behind it (Fig. 3) is figure_only, and there is no measured dose in 8-10 (Table 1: 3.6, 7.2, 14.4, 21.6). **NOT a titre.** |
| the same, as a concentration in this paper's pot | 80-100 (midpoint 90) | mmol/L | as above | — | derived from 1 % t.s. (mine) | derived_assumption |
| the same, applied to Hofmann's 12.5 g/L MW>3000 fraction | 112.5 | mmol/L | **a different material** (Flags 1) | — | derived (mine) | derived_assumption |
| k_ref, PhSH | **> 7.70e-4** | **s^-1** (header prints `mol-1 s-1`; Flags 2) | as above | pseudo-first order in thiol | Table 2 p. 4430 | measured_rate — **lower bound only**, set by the 1 h sampling floor |
| k_ref, **FFT** | **> 7.70e-4** | **s^-1** | as above | pseudo-first order in thiol | Table 2 | measured_rate — **lower bound only** |
| k_ref, BenzSH | **> 7.70e-4** | **s^-1** | as above | pseudo-first order | Table 2 | measured_rate — lower bound only |
| k_ref, PentSH | 1.83e-4 | s^-1 | as above | pseudo-first order | Table 2 | measured_rate |
| k_ref, BuSH | 1.42e-4 | s^-1 | as above | pseudo-first order | Table 2 | measured_rate |
| k_ref, PropSH | 1.32e-4 | s^-1 | as above | pseudo-first order | Table 2 | measured_rate |
| k_ref, **EtSH** | **1.19e-4** | **s^-1** | as above | pseudo-first order | Table 2 | measured_rate — **the best-determined constant in the paper** (a full aliphatic decay curve, not a bound) |
| k_ref, 2BT | 8.11e-5 | s^-1 | as above | pseudo-first order | Table 2 | measured_rate |
| k_ref, MMBF | 2.12e-5 | s^-1 | as above; note the ester hydrolyses (ref 22) | pseudo-first order | Table 2 | measured_rate (a lumped hydrolysis + addition constant, per the paper's own reading) |
| k_ref, 2M2P | 1.02e-6 | s^-1 | as above | pseudo-first order | Table 2 | measured_rate |
| half-lives from Table 2 | <15.0 / 63.1 / 81.4 / 87.5 / 97.1 / 142 / 545 / 11330 | min, for (PhSH,FFT,BenzSH) / PentSH / BuSH / PropSH / EtSH / 2BT / MMBF / 2M2P | as above | — | derived from Table 2 (mine) | derived_assumption (arithmetic only) |
| **oxygen dependence of the sink, EtSH** | **17.9x attenuation under N2** | dimensionless (1/k_rel) | 15 h, otherwise as above | ratio of two pseudo-first-order constants | Table 3 | **within_study_ratio** |
| oxygen dependence, PropSH / BuSH / PentSH | 19.8 / 15.9 / 16.4 | dimensionless | " | " | Table 3 | within_study_ratio |
| **oxygen dependence, FFT** | **45.1x** | dimensionless | " | " | Table 3 | **within_study_ratio** — the benzylic thiol is 2.5x more oxygen-dependent than the aliphatics |
| **nucleophile-blockable fraction, EtSH** | **10.5x attenuation at 10 mmol/g hydroxylamine** | dimensionless | 1 h preincubation, aerobic | " | Table 3 | within_study_ratio |
| the same, PropSH / BuSH / PentSH | 11.4 / 10.8 / 13.4 | dimensionless | " | " | Table 3 | within_study_ratio |
| **the same, FFT** | **4.0x only** | dimensionless | " | " | Table 3 | within_study_ratio — **FFT is the thiol least protected by blocking electrophiles**, i.e. the one the repository models is the one for which the `MELE` channel is least dominant |
| combined N2 + hydroxylamine (vs aerobic) | 64.8 / 82.5 / 101.3 / 236.1 / 92.6 for EtSH/PropSH/BuSH/PentSH/FFT | dimensionless | 21 mmol/g, anaerobic | " | Table 3 footnote a | within_study_ratio (**FFT's row is internally inconsistent**, Flags 5) |
| Na2SO3 attenuation | 51.5 / 61.7 / 18.4 / 62.7 / 67.9 | dimensionless at 1 mmol/g | different denominator (Flags 6) | " | Table 3 | within_study_ratio (weakened) |
| ascorbic acid attenuation | 9.0 / 7.8 / 4.4 / 3.7 / 23.6 | dimensionless at 1 mmol/g | aerobic | " | Table 3 | within_study_ratio |
| caffeic acid | 0.3 / 0.3 / 0.3 / 0.3 / 1.0 | dimensionless at 0.28 mmol/g | aerobic | " | Table 3 | within_study_ratio — **values below 1 mean caffeic acid ACCELERATES the loss ~3x** for the aliphatics |
| DTPA | 0.5 / 0.7 / 0.7 / 1.0 / 1.0 | dimensionless at 0.5 mmol/g | aerobic | " | Table 3 | within_study_ratio — also accelerating or neutral |
| activation energy of any of it | **absent** | — | one temperature only (25 C) | — | — | **not present in this paper** |
| thiol pKa values | 8.6 / 6.5 (PhSH), 11.3 (FFT), 11.8 (BenzSH), 12.6 (BuSH), 12.9 (2BT), 13.1 (2M2P) | — | acetone/H2O 3:1, except 6.5 in water | — | Table 2, from refs 23/24 | level_only (**not measured here**; borrowed) |
| relative thiolate nucleophilicity | 1400 (FFT), 1200 (BenzSH), 1000 (BuSH), 380 (2BT), 340 (2M2P) | dimensionless | maleic-anhydride addition, Dmuchovsky 1966 | — | Table 2 | level_only (borrowed) |
| rate vs log P for linear aliphatics | — | — | — | — | Fig. 11 | **figure_only** (no slope, intercept or r^2 printed) |
| every concentration-time course | — | % of blank vs h | — | — | Figs. 2, 4, 5, 6, 7, 8, 9 | **figure_only** |

### Does the paper support the repository's recasting? The four-part answer.

**(a) The number is printed and the page cite is right.** "8-10 mmol/g dry coffee" appears verbatim
on p. 4428 and again on p. 4430. `species_sulfur.py`'s `MELE` note and `parameters_sulfur.py`'s
`source_anchor` both quote it correctly.

**(b) "By hydroxylamine titration" over-states the method, and the repository should soften that
phrase.** There is no titration, no endpoint, no equivalence point, no analytical determination of
consumed hydroxylamine. There is a protection-versus-dose curve on four doses (Figure 3) whose
plateau the authors locate by eye between the second and third dose, and which they themselves label
"an indirect method". A defensible re-wording for the `MELE` note: *"a saturating hydroxylamine dose
above which thiol protection stops increasing (Charles-Bernard 2005, 8-10 mmol per g dry coffee
solids on whole 1 % brew solids, read from a four-point figure) — an upper bound on, not a
measurement of, the thiol-reactive site density."*

**(c) The material is the wrong one for where the number is used.** The density is measured on
**whole soluble coffee solids**, and the paper says in the same paragraph that "these electrophilic
sites include also reducing sugars, which present moderate to low reactivity toward thiols." The
repository multiplies it by **12.5 g/L of Hofmann's MW > 3000 melanoidin fraction**, from which
every reducing sugar and every other small electrophile has been removed by the cutoff. The transfer
therefore requires assuming that a high-MW polymer fraction has the same electrophile density per
gram as the whole extract, **after the paper has told you that a chunk of the whole-extract number
is small molecules.** This is the single largest unstated assumption in the `k_thioether` entry, and
it biases [E] **upward**, hence k2 **downward**, by an unknown factor.

**(d) The number is above its own stoichiometric ceiling.** Section 3 arithmetic 4: 9 mmol/g demands
one electrophile per 111 g of dry solids, below the mass of a single hexose residue. Combined with
(c), the working conclusion is that **`k_thioether` = 5.01e-4 L mmol^-1 min^-1 is an under-estimate
of the true bimolecular constant by one to two orders of magnitude, and `[MELE]0` is an
over-estimate by the same factor.** Because the lane only ever uses the product, the shipped
predictions are unaffected — which is also why the error has been invisible. It becomes visible the
moment anything in the lane changes the site concentration independently (dilution series, a
solids-loading axis, the 80 C brew hold-out where `MELE` is supposed to be partly pre-consumed).
**Anything that leans on `MELE` depleting is leaning on nothing at the current density.**

### What this paper adds that Hofmann 2002 does not

Two things, both real. First, **an independent laboratory, an independent matrix (a real brew, not a
model melanoidin solution), an independent method (SPME headspace ratio against a matched blank, not
SIDA) and an independent pH (5.2, not 6.0), reaching the same order of magnitude for the thiol-loss
constant**: FFT > 7.70e-4 s^-1 here against Hofmann's 9.4e-4 s^-1 at 30 C. That is a genuine
cross-laboratory corroboration of the *rate*, and it is what the `source_anchor` claims. Second,
**the oxygen dependence** — 15.9-19.8x for the aliphatics and 45.1x for FFT — which is the only
measurement in the corpus that ties the electrophile pool to the oxidant state, and therefore the
only warrant for `OX` and `MELE` being co-required.

## 5. Flags

1. **"8-10 mmol/g" is a saturating reagent dose, not a site titre, and it exceeds the
   stoichiometric ceiling for a carbonyl count on coffee solids.** One site per 111 g/mol is below
   a hexose residue (162 g/mol). Hydroxylamine oximation is an equilibrium requiring excess, and the
   dose is consumed by everything in the matrix, not only by thiol-reactive sites. **Treat the
   number as an upper bound on the site density**, hence `k_thioether` as a lower bound on the
   bimolecular constant, and say so in the parameter's note. The registry's product
   `k_thioether * [MELE]0` is what the data constrain and it is unaffected.
2. **The unit on Table 2's header is wrong and this paper proves it internally.** The header prints
   `[mol-1 s-1]`; the values are **s^-1**. Three independent proofs are in section 3 arithmetic 5
   (the Methods sentence stating a first-order ln-vs-seconds construction; the anaerobic 15 h
   degradation reproducing the printed "20-30 %"; the absurdity of a second-order reading at
   micromolar thiol). The existing registry flag
   `charles_bernard_table2_unit_erratum_values_are_per_second` is correct and now rests on the
   primary paper.
3. **Three of the ten constants are detection-floor bounds, not measurements, and FFT is one of
   them.** PhSH, FFT and BenzSH all print "> 7.70 x 10^-4", which is exactly 4 ln2 / 3600 s — the
   constant that empties the vial to 1/16 by the first sampling time. **The repository's thiol
   therefore has no measured rate in this paper, only a floor.** Never quote 7.70e-4 as FFT's
   constant; quote it as "faster than the 1 h sampling could resolve".
4. **No temperature series, so no activation energy, so nothing here licenses any rate transfer.**
   `k_thioether`'s `ea` comes from Stack 2018 and its `rate_transfer="not_licensed"` is right. This
   paper is a single point at 25 C, and Hofmann's is a single point at 30 C; two points 5 K apart
   from two laboratories in two matrices cannot make a slope.
5. **Table 3's FFT column is internally inconsistent by a factor of 2.05.** Footnote b's value times
   the anaerobic value reproduces footnote a for EtSH, PropSH, BuSH and PentSH to within 1 %, but
   for FFT gives 184.9 against a printed 92.6. One of the two FFT entries is wrong and the paper
   gives no way to tell which; the printed 4.1 may be a transcription of the 4.0 in the row above.
   **This is the repository's thiol.** Do not use FFT's combined-condition attenuation without
   saying so.
6. **The two reducing-agent rows have a different denominator from every other row of Table 3.**
   Na2SO3 and Na2S2O3 depress the *blank* headspace, so their results are referenced to a blank
   without additive while all other additives are referenced to blank + additive. The Na2SO3
   attenuations (18.4-67.9x, the largest in the table) therefore contain an unremoved
   blank-suppression term and are not comparable with the hydroxylamine or anaerobic rows.
7. **Table 3's hydroxylamine dose, "10 mmol/g", matches none of the four doses in Table 1** (3.6,
   7.2, 14.4, 21.6). The "21" of the combined rows is plainly 21.6 rounded, but 10 is not a rounding
   of any of them. Either a fifth dose was run and not tabulated, or the entry is nominal. **Do not
   pair Table 3's hydroxylamine column with a Table 1 dose without asking the authors.** It is also
   the number closest to the claimed 8-10 mmol/g plateau, so the coincidence is worth resolving.
8. **Minor internal inconsistencies.** The Methods list ten thiols and Table 2 has ten rows, but the
   stability-control sentence says "the nine thiols". The abstract and conclusion name
   **methanethiol** among the thiols the paper explains, but methanethiol is in neither the mixture
   nor Table 2 — its behaviour is imported from the introduction's references, not measured here.
9. **Table 2's "relative nucleophilicity" column is credited to reference (24)** — Pascual et al.'s
   EPR study of coffee radicals — **while the text credits the same quantity to reference (23)**,
   Dmuchovsky et al.'s maleic-anhydride study, which is evidently the correct source. Treat the
   column as Dmuchovsky 1966.
10. **What this paper does NOT contain**: any activation energy or temperature dependence; any pH
    series; any water-activity or solids-loading series (a single 1 % t.s.); any fractionation of
    the matrix (that is Part 1, not on disk); any identification of the electrophile — Figure 1 and
    Figure 10 are hypotheses, and the words "quinone" and "CROSSPY" appear only in the introduction
    as literature; any measurement on 2-methyl-3-furanthiol or methanethiol; any disulfide or other
    product measurement (**no mass balance anywhere: the thiol that disappears is never accounted
    for**); any confidence interval on any rate constant; any absolute concentration; any
    supplementary material.
11. **What to request from the authors**: (i) the numeric dose-response behind Figure 3, and in
    particular whether any dose between 7.2 and 14.4 mmol/g was ever run — the "8-10" claim depends
    entirely on this; (ii) the hydroxylamine dose actually used in Table 3's third row; (iii) the
    resolution of the FFT inconsistency in Table 3; (iv) the numeric decay curves behind Figure 2,
    with the number of points per fit and the r^2 of each pseudo-first-order regression, none of
    which is printed; (v) whether the hydroxylamine consumed by the matrix was ever measured
    directly — that would convert the level_only site density into a real titre and is the single
    experiment that would repair `MELE`; (vi) Part 1 (53:4417-4425) for the molecular-weight
    fractionation, which is what would license or refuse the transfer of the whole-solids density to
    a MW > 3000 fraction.
12. **Registry gaps against `data/keys/compounds.yml`**: `2_furfurylthiol` is present and is the
    only studied compound that is. **Absent: ethanethiol, propanethiol, butanethiol, pentanethiol,
    2-butanethiol, 2-methyl-2-propanethiol, 3-mercapto-3-methylbutyl formate, benzylthiol,
    thiophenol, hydroxylamine, ascorbic acid, caffeic acid, DTPA, and any melanoidin entry.** The
    aliphatic homologous series EtSH-PropSH-BuSH-PentSH is the paper's cleanest quantitative result
    (four constants spanning 1.54x with a matched matrix) and **none of its four members is keyed**,
    so it cannot currently be used as a benchmark. `MELE` and `OX` are network-local species names
    in `species_sulfur.py`, not registry ids, and correctly so — neither is a molecule.
