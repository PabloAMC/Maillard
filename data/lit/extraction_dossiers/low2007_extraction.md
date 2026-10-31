# Low, Parker & Mottram 2007 — EXTRACTION
### *"Mechanisms of Alkylpyrazine Formation in a Potato Model System Containing Added Glycine"*
### Potato cakes (drum-dried flake + 1 % w/w [2-13C]glycine solution, 1 : 1.3), 1.3 % w/w glycine (dry basis), 180 °C impingement oven, 15 and 30 min. 13C isotopomer shares (unlabelled / singly / doubly labelled) of 22 alkylpyrazines and two alpha-diketones by dynamic headspace GC-MS. **No absolute concentrations in this paper.**

**Read 2026-10-31 from `data/articles/low2007.pdf`** (8 pp., JAFC 55:4087-4094; born-digital, clean text layer; Table 1
checked against a 110 dpi render of p. 4089, every cell matches). Figures 1-5 are reaction schemes; Fig. 2 prints the
expected isotopomer proportions in square brackets, but those are scheme labels, not measurements. No row of this paper is
keyed anywhere in the repository yet.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[Q]** qualitative statement or
interpretation by the authors · **[C]** cited from another paper · **[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Context: the pyrazine lane (`src/kinetic_core/parameters_pyrazine.py`, B18) carries two Strecker steps. Glyoxal + Gly gives
the aminoacetaldehyde pool (AKG); methylglyoxal + Gly gives aminoacetone + CO2 + formaldehyde (AKM). Three aminoketone
condensations share one fast constant: AKG + AKG -> pyrazine, AKM + AKM -> 2,5-dimethylpyrazine, AKG + AKM -> methylpyrazine.
The lane has no dihydropyrazine species and no aldehyde-addition step. MECHANISM_REVIEWS sec. 1 finding 6 and the pyrazine
review (`docs/history/mechanism_reviews_2026-10-31/pyrazines.md`) ask whether the lane's precursor routes are right and
whether methylpyrazine has an extra route (the Leahy "mixing test").

| need | what Low 2007 gives | tag |
|---|---|---|
| (a) which pyrazines take glycine carbon | **essentially none** (<= 7 % labelled, any time): methyl-, ethyl-, 2,5-dimethyl-, isobutyl-, 2-methyl-6-isobutylpyrazine, and pyrazine itself (text; pyrazine is not in Table 1) | [M] / [Q] |
| | **About half labelled, singly, at a methyl:** 2,6-dimethyl- (SL 42.9 / 45.3 % at 15 / 30 min), 2,3-dimethyl- (78.8 / 70.4 %), trimethyl- (80.9 / 71.6 %), 2-ethyl-3,5-dimethyl- (73.9 / 68.3 %) | [M] |
| | **Doubly labelled:** tetramethylpyrazine (DL 54.6 / 46.6 %, SL 39.0 / 43.3 %); 2,5-diethyl-3-methyl- (DL 21.8 / 26.7 %) | [M] |
| | **Partly labelled:** the ethylmethyl isomers (2-ethyl-6-methyl- SL 21.5 / 26.6 %, 2-ethyl-5-methyl- 8.0 / 10.5 %) and propyl- / methylpropylpyrazines (label in the propyl chain) | [M] |
| (b) label in the alpha-diketones | 2,3-butanedione SL 56.3 / 42.4 %, DL 2.9 / 6.8 %, "singly labeled at C1"; 2,3-pentanedione SL 20.9 / 34.1 %, label at C1 (or C2) | [M] |
| (c) the two pathways | **"x + x"**: two alpha-aminocarbonyls (from alpha-dicarbonyl / hydroxycarbonyl + amino acid, Strecker) condense to a dihydropyrazine, which is oxidised (or dehydrated) to the pyrazine. **"x + x + y"**: a Strecker aldehyde y (formaldehyde from Gly, acetaldehyde from Ala, methylpropanal from Val, ...) adds to the dihydropyrazine before aromatisation and becomes an extra alkyl substituent; only one y per ring | [Q] |
| (d) their relative contributions (authors' arithmetic) | x + x + y up to **~45 % of 2,6-DMP** (x + x <= 55 %); **~50 % of 2,3-DMP**; **~50 % of trimethylpyrazine**; tetramethylpyrazine more doubly labelled than the 1 : 2 : 1 that x + x alone would give. x + x + y with formaldehyde gives **25 % of 2-ethyl-6-methyl-** and **9 % of 2-ethyl-5-methylpyrazine**, and ~70 % of 2-ethyl-3,5-dimethyl- (labelled at the 3-methyl). 2,5-DMP comes only from x + x | [Q] / [M] inputs |
| (e) conditions | 180 °C, impingement oven, 15 or 30 min; cakes 3 mm x 73 mm, ~18 g; 15 cakes per batch milled together; triplicate analyses of the pooled powder | [M] |
| (f) water | **not measured or stated.** Dough = flake : 1 % glycine solution 1 : 1.3 w/w, so ~56 % water at the start plus the flake's own moisture [Z]; end-of-bake moisture and product temperature not reported | [!] |
| (g) quantification | dynamic headspace (7 g sample + 93 g water, 37 °C, 1 h, N2 purge onto Tenax), GC-MS EI. Isotopomer shares by the Hwang 1993 ion-cluster method, corrected for natural 13C and [M-H]+. 1,2-dichlorobenzene IS (used for concentrations in the companion paper, not here). No SDs printed. The "change" column (+ / - / nc) is quoted from the companion, Low et al. 2006 | [M] |

**Answer.**

1. **The lane's three species are pure "x + x" products whose carbons come from sugar fragments.** Methylpyrazine, pyrazine
   and 2,5-DMP carry no glycine carbon, even at 400x glycine excess (2,5-DMP 3.7 - 6.3 % SL). Amrani-Hemaimi 1995 found the
   same for [3-13C]alanine (`amranihemaimi1995_extraction.md`). The aminoketone's carbons come from the dicarbonyl, not the
   amino acid; the lane's structure agrees. Confidence: ~0.85 for these three species in a sugar-amino-acid matrix at 160 - 180 °C;
   aqueous 95 °C (Leahy) is untested.
2. **Formaldehyde, the lane's Strecker co-product, is not a dead end.** It alkylates dihydropyrazines, and in this
   glycine-rich potato it makes about half of 2,6-DMP, 2,3-DMP and trimethylpyrazine. The lane does not carry these
   species, so nothing breaks today. **If the lane adds 2,3- / 2,6-DMP, TMP or TeMP, it needs a dihydropyrazine
   intermediate with a branch: aromatisation vs Strecker-aldehyde addition.** One shared condensation constant, with the
   ring formed and oxidised in one step, cannot produce these isotopomer patterns. Confidence that aldehyde addition carries
   >= 30 % of 2,3-DMP / TMP in a Gly-rich matrix: ~0.7. The shares below are model-dependent inferences, and Fig. 1 marks the
   x + x + y route to 2,5-DMP as possible but "not observed in practice".
3. **Methylpyrazine's extra route (the Leahy mixing ratio 11.7 - 17.3 vs 4) is not formaldehyde addition to the bare
   dihydropyrazine, at least not glycine-derived formaldehyde.** Methylpyrazine is 95 - 100 % unlabelled although 2,6-DMP
   takes ~45 % from labelled HCHO. Caveats:
   - Sugar-derived HCHO is invisible to the label, so a sugar-HCHO route is not excluded.
   - Leahy's pot is lysine + glucose in water at 75 - 95 °C, not this matrix.

   Confidence this paper disfavours a glycine-HCHO route for methylpyrazine: ~0.7. Its bearing on Leahy: weak (~0.55).
   2,6-DMP is not detected in Leahy's lysine-glucose pots (`leahy1989a_extraction.md`), so leaving 2,6-DMP out of the mixing
   ratio is not an artefact there.
4. **The dicarbonyl supply depends on the amino acid.** Glycine carbon ends up in 2,3-butanedione at C1 (~50 % singly
   labelled). Glycine raises butanedione and lowers 2,3-pentanedione (by 50 - 60 %, quoted from Low 2006) [Q]. Glycine
   addition also cut pyrazine, methyl-, ethyl-, propyl- and isobutylpyrazine "in almost equal proportions (approximately
   70%)" [Q, from Low 2006]: competition for a shared precursor (glyoxal, or the unsubstituted dihydropyrazine). Neither
   effect is in the lane, whose dicarbonyl supply comes from the trunk with no amine identity. Both point to the
   composable-network amine classes (roadmap step 6), not to step 2.
5. **No rate, no yield, no temperature or water series.** The 15 vs 30 min contrast changes yields (absolute numbers are in
   Low 2006) but "labeling patterns were similar irrespective of cooking time" [Q]. Nothing here calibrates a constant;
   everything is a scale-free carbon-origin signature.

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Mei Yin Low, Jane K. Parker, Donald S. Mottram* (Department of Food Biosciences, University of Reading) |
| Title | as in the header |
| Venue | *J. Agric. Food Chem.* **2007**, 55 (10), 4087-4094 |
| DOI | **10.1021/jf070044s** (printed p. 4087; article id "JF070044S" p. 4094) |
| Dates | received 6 Jan 2007; revised 16 Mar 2007; accepted 25 Mar 2007; published on web 21 Apr 2007 |
| Content | Table 1 (13C isotopomer shares, 2 diketones + 22 pyrazines, 15 / 30 min, label positions); Figs. 1-5 (x + x and x + x + y schemes for di-, tri- / tetra-, ethylmethyl- / methylpropyl-, ethyldimethyl-, diethylmethylpyrazines) |
| Companion (holds the concentrations) | ref 5: Low, Koutsidis, Parker, Elmore, Dodson & Mottram 2006, JAFC 54:5976-5983, "Effect of citric acid and glycine addition on acrylamide and flavor in a potato model system" (26 pyrazines, butanedione, pentanedione; 180 °C, 15 - 60 min). **Not on disk; DOI not verified here.** The PhD thesis (Low 2006, Reading) carries the multi-response kinetics that `balagiannis2015_extraction.md` row 7 describes |

## 3. Methods [M]

| item | as printed |
|---|---|
| Matrix | drum-dried potato flake (Maris Piper : Pentland Dell 1 : 1, McCain), prepared after Elmore et al. 2005 |
| Glycine | unlabelled or [2-13C]glycine (99 %) as a 1.0 % w/w aqueous solution; flake : solution 1 : 1.3 w/w; "glycine concentration in the product of 1.3% w/w (dry weight) which was a 400-fold increase in glycine over that in the untreated potato flake". [Z] ~173 mmol/kg dry added; native ~0.43 mmol/kg dry. Native alanine "about 5.5 times more than glycine" (p. 4092), i.e. ~2.4 mmol/kg dry [Z] |
| Forming | mechanical dough maker and roller, 3 mm sheet, 75 mm cutter, docked; raw cakes 3 mm x 73 mm, ~18 g |
| Baking | 180 °C, Impinger II moving-band impingement oven (conveyor off), 15 or 30 min; 3 rows of 5 cakes under the upper-deck jets. Cooled, 15 cakes pooled and milled ("average heating conditions in the oven"), 3 portions analysed |
| Volatiles | 7 g powder + 93 g deionised water, 37 °C, 1 h, oxygen-free N2 at 40 mL/min onto 85 mg Tenax TA; IS 100 ng 1,2-dichlorobenzene added to the trap |
| GC-MS | HP 5890 / 5972, ZB-Wax 60 m x 0.25 mm x 0.25 µm, EI 70 eV, m/z 29-400; LRI from C6-C25 alkanes |
| Label calculation | per pyrazine, % unlabelled (UL), singly (SL), doubly (DL) labelled; "No triply labeled pyrazine species were found"; label position from fragment ions (Porter 1985 fragmentation); corrected for natural 13C and loss of one H |
| Not stated | end moisture, product core temperature, oven air speed, replicate SD, pH of the dough |

## 4. Results [M]

### 4.1 Table 1 (13C isotopomer shares, %; "change" = concentration change with glycine, from Low 2006)

| compound | change 15 | UL 15 | SL 15 | DL 15 | change 30 | UL 30 | SL 30 | DL 30 | label position (as printed) |
|---|---|---|---|---|---|---|---|---|---|
| 2,3-butanedione | + | 40.9 | 56.3 | 2.9 | (blank) | 50.8 | 42.4 | 6.8 | singly labeled at C1 |
| 2,3-pentanedione | − | 76.8 | 20.9 | 2.3 | (blank) | 57.8 | 34.1 | 8.1 | singly labeled at C1 (or C2) |
| methyl- | − | 95.2 | 3.7 | 1.1 | − | 100 | 0.0 | 0.0 | |
| 2,5-dimethyl- | nc | 96.0 | 3.7 | 0.3 | nc | 93.3 | 6.3 | 0.4 | |
| 2,6-dimethyl- | nc | 56.4 | 42.9 | 0.7 | nc | 52.4 | 45.3 | 2.3 | singly labeled at a methyl substituent |
| ethyl- | − | 93.8 | 6.0 | 0.2 | − | 95.9 | 1.4 | 2.8 | |
| 2,3-dimethyl- | + | 19.8 | 78.8 | 1.4 | + | 27.0 | 70.4 | 2.7 | singly labeled at a methyl substituent |
| 2-ethyl-6-methyl- | − | 76.6 | 21.5 | 1.9 | − | 71.2 | 26.6 | 2.2 | singly labeled at methyl substituent |
| 2-ethyl-5-methyl- | nc | 91.4 | 8.0 | 0.6 | − | 88.6 | 10.5 | 1.0 | singly labeled at methyl substituent |
| trimethyl- | + | 15.5 | 80.9 | 3.6 | + | 22.2 | 71.6 | 6.3 | singly labeled at 2- or 3-methyl substituent |
| propyl- | − | 25.8 | 18.6 | 55.7 | − | 45.3 | 38.7 | 16.0 | at least one 13C atom in propyl substituent |
| 2,6-diethyl- (e) | − | 86.6 | 13.4 | 0.0 | − | 90.5 | 8.2 | 1.3 | |
| 2-ethyl-3,6-dimethyl- | + | 94.9 | 4.1 | 1.0 | − | 90.3 | 7.6 | 2.1 | singly labeled at 3-methyl substituent |
| isobutyl- | − | nd | nd | nd | − | 94.1 | 2.4 | 3.5 | |
| 2,3-diethyl- | − | nd | nd | nd | − | 87.9 | 9.8 | 2.3 | |
| 2,5-diethyl- (e) | nc | 83.3 | 14.0 | 2.7 | nc | 81.2 | 15.3 | 3.6 | |
| 2-ethyl-3,5-dimethyl- | + | 18.5 | 73.9 | 7.7 | + | 23.1 | 68.3 | 8.6 | singly labeled at 3-methyl substituent |
| 2-methyl-6-propyl- (e) | − | nd | nd | nd | − | 49.8 | 43.5 | 6.7 | |
| tetramethyl- | + | 6.4 | 39.0 | 54.6 | + | 10.1 | 43.3 | 46.6 | singly at a methyl; doubly at 2,5- or 2,6-methyl substituents |
| 2-methyl-3-propyl- | nc | nd | nd | nd | − | 42.5 | 36.3 | 21.1 | |
| 2-methyl-6-isobutyl- (f) | nc | 100 | 0.0 | 0.0 | − | 97.1 | 2.4 | 0.5 | |
| 2,3-diethyl-5-methyl- (e) | nc | 90.1 | 8.7 | 1.2 | nc | 83.8 | 13.7 | 2.4 | singly labeled at methyl substituent |
| 2,6-diethyl-3-methyl- (e) | nc | 90.3 | 8.5 | 1.2 | − | 86.5 | 11.6 | 1.9 | singly labeled at methyl substituent |
| 2,5-diethyl-3-methyl- (e) | + | 33.1 | 45.1 | 21.8 | + | 31.6 | 41.8 | 26.7 | singly at methyl; doubly at methyl and an ethyl substituent |

(e) tentative identification by literature MS and LRI; (f) tentative, MS only; nd not detected; footnote h "1.3% w/w dry
wt." 2-Ethyl-3-methylpyrazine: not determinable (co-elutes with trimethylpyrazine). Unsubstituted pyrazine: not in the
table; the text calls it unlabelled [!].

### 4.2 The authors' pathway arithmetic [Q], with my check [Z]

The authors average 15 and 30 min. They assume every formaldehyde is 13C (glycine in large excess), that 2,3-butanedione
is ~50 % singly labelled and not limiting, and that the C2 / C3 dicarbonyls are unlabelled. Let L be the labelled fraction
(SL + DL) and assume the x + x + y product is always labelled. Then the x + x + y share is f = (L_obs − L_xx) / (1 − L_xx).

| pyrazine | x + x precursors | L_xx assumed | L_obs (avg) | authors' x + x + y share | my f [Z] |
|---|---|---|---|---|---|
| 2,6-DMP | MGO + MGO (minor amination site) | ~0 (2,5-DMP 4 - 7 %) | 45.6 % | "as much as 45%" | 0.43 - 0.46 |
| 2,3-DMP | glyoxal + 2,3-butanedione | ~0.5 | 76.7 % (SL 74.6) | "approximately 50%" | 0.49 - 0.54 |
| trimethyl- | MGO + butanedione | ~0.5 | 81.2 % (SL 76.3) | "about 50%" | 0.53 - 0.62 |
| tetramethyl- | butanedione x 2 | 1 : 2 : 1 (UL : SL : DL) | ~8 : 41 : 51 | "higher percentage of the doubly labeled species" (both routes) | not unique |
| 2-ethyl-6-methyl- / 2-ethyl-5-methyl- | x + x + y with HCHO | -- | SL 24 / 9 % | 25 % / 9 % via glycine-HCHO; x + x + y with acetaldehyde (Ala) judged larger | -- |
| 2-ethyl-3,5-dimethyl- | x + x + y with HCHO on 2-ethyl-5-methyl-dihydropyrazine | -- | SL 71 % | "Approximately 70%"; 2-ethyl-3,6-DMP mostly from Ala-acetaldehyde (SL 4 - 8 %) | -- |

[!] f is a LOWER bound on aldehyde addition if some formaldehyde is sugar-derived (unlabelled). It is an UPPER bound if
glycine carbon also enters the x + x dicarbonyls; butanedione shows glycine carbon does enter some of them. The two biases
pull opposite ways, so the shares are uncertain at roughly +-15 points [Z].

Regiochemistry [Q]: MGO amination prefers the aldehydic site, which favours 2,5- over 2,6-DMP in x + x. HCHO adds to
2-methyl-3,6-dihydropyrazine "at carbons 3 or 6 with equal preference" (the 2,6- and 2,3-DMP x + x + y yields are similar).
Ethyl addition (acetaldehyde) prefers carbon 6. Methyl from glycine adds preferably at carbon 3 of 2,5- / 2,6-disubstituted
dihydropyrazines (2-ethyl-5,6-DMP not formed).

### 4.3 Qualitative results quoted from the companion (Low 2006) [C]

Glycine at 1.3 % raised 2,3-DMP > 2-fold, tri- and tetramethylpyrazine "more than doubled", raised 2-ethyl-3,5-DMP and
2,5-diethyl-3-methyl-, and cut pyrazine / methyl- / ethyl- / propyl- / isobutyl- by ~70 %. 2,5-DMP was the largest pyrazine
with or without glycine, followed by 2,6- and then 2,3-DMP. 2,3-Pentanedione fell 50 - 60 %; 2,3-butanedione rose.
Without glycine, 2-ethyl-3,6-DMP was 3 - 5x 2-ethyl-3,5-DMP.

## 5. Flags

- [!] **No concentrations, no SDs, no yield.** Every quantity is an isotopomer share. Concentrations, and the + / − / nc
  calls, live in Low 2006 (not on disk).
- [!] **Matrix conditions are open:** no water content, product temperature or pH; impingement oven with surface browning,
  so the cake is not isothermal. Pooled 15 cakes, so the shares are batch averages.
- [!] **Glycine is 400x its native level**, so the shares describe a glycine-flooded potato, not a native one. The authors
  say so: the x + x + y share is sensitive to the amino-acid mix.
- [!] Pathway shares are the authors' inferences from the label pattern under stated assumptions (sec. 4.2), not
  measurements. The x + x / x + x + y mechanisms are literature hypotheses (Amrani-Hemaimi 1995, Chiu 1990, Weenen 1994,
  Rizzi 1972 [C]).
- [!] Pyrazine itself is called unlabelled in the text but has no Table 1 row.
- [!] The 30-min "change" cells for the two diketones are blank in the table.
- [!] Amrani-Hemaimi's elution order of 2-ethyl-5- / 2-ethyl-6-methylpyrazine is disputed here (the authors use the
  "generally accepted" DB-Wax order: 2-ethyl-6- first). Any cross-paper comparison of those two isomers inherits this.

## 6. Supplementary material

**None.** There is no "Supporting Information Available" paragraph and no SI reference. The data the model would need
next, the absolute pyrazine and diketone concentrations with and without glycine at 15 - 60 min, are in the companion Low
et al. 2006 (JAFC 54:5976, ref 5) and the Low 2006 PhD thesis (Reading). The companion is the paper to download if the
pyrazine lane ever wants potato-matrix yields or the glycine effect as a number.

## 7. Keying recommendation

- **Not worth two-reader keying now.** The lane has no isotopomers and none of the labelled species (2,3- / 2,6-DMP, TMP,
  TeMP). Table 1 can serve only as scale-free signatures for a future dihydropyrazine + aldehyde-addition arm. If that arm is
  pre-registered, key ~8 rows (2,5-, 2,6-, 2,3-DMP, TMP, TeMP, methylpyrazine, 2,3-butanedione, 2-ethyl-3,5-DMP; 15 and 30
  min) as HOLD-OUT signatures, never as fit rows. That is ~30 min for two readers; the text layer is clean (raster-checked).
- Usable today without keying, as a one-line signature: methylpyrazine and 2,5-DMP carry no amino-acid carbon (this paper
  plus Amrani-Hemaimi 1995).
