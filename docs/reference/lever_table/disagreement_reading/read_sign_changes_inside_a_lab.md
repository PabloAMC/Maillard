# Sign changes inside one laboratory: reading of eight cells

*Written 2026-11-12. Literature and contrast files only; no engine output. Effects are in dex per declared step
(per doubling of time or dose, per pH unit, per 10 C, or level against level). "Endpoint" slopes below are my own
recomputations from the keyed arms; "keyed" values are the ones in the cell list.*

## What was read

| source | read in full | read in part | not read (contrast file only) |
|---|---|---|---|
| zhou2023c | text layer, Fig. 1 heat map as image (p. 3) | | |
| feng2022 | methods; Tables 1-2 thiol rows; results pp. 9098-9101 | discussion | supplement |
| hofmann1998 | methods; Tables 1-8; results text | | references |
| cerny2007 | abstract; Table 2; MFT row of Table 3 | discussion | |
| grosch1992 | Tables III-IV and the text around them (pp. 185-187) | rest of chapter | |
| Zhang2024 | | secs. 2.3, 2.6, 3.2 | Fig. 2d not viewed (values digitised by readers) |
| li2021 + supplement | | main text (abstract, methods 2.6-2.7, secs. 3.1-3.2); Table S2 benzaldehyde and benzeneacetaldehyde rows, all six blocks | Li et al. 2021b (the extraction method) |
| Theng2024 | Table 1 benzaldehyde row (pp. 6330-6331) | methods | |
| leahy1989 | methods; Table III text; results and discussion pp. 81-89 | Table I (via the contrast file) | |
| us10076126, thomas2014, ding2023, martin2001, zamora2006, monforte2021, hidalgo2013, zamora2015, boonbumrung2023 | | | contrast files only |

## 1. 2-Furfurylthiol (FFT) x time (Ho network: feng2022 +0.12, +0.20 against zhou2023c -0.50)

**Designs.** Both are Jiangnan / Ho papers, 20 mmol/L Amadori product + 20 mmol/L GSH, 120 C, 30 - 120 min,
unbuffered (NaOH-set initial pH), HS-SPME DVB/CAR/PDMS with 1,2-dichlorobenzene internal standard and
calibration curves. They differ in the Amadori product and pH: feng2022 uses the ribose-**glutathione** ARP
(RG-ARP, so the sulfur source is also built into the ARP) at pH 7.0; zhou2023c uses the **alanine-xylose** ARP at
pH 7.5 (zhou2023c Methods, p. 14313).

**What the data show.** feng2022 Table 2 (ARP + GSH): FFT 14.31, 30.45, 39.35, 44.63 ug/L at 30, 60, 90, 120 min,
monotone rise. zhou2023c Fig. 1 (ARP + GSH): FFT 124, 383, 270, 116 - a **peak at 60 min**. The table keeps the
longest monotone run (60 - 120 min, falling) and drops the rising 30 -> 60 min step (+0.49 dex per doubling, the
same sign as feng2022). Over the full 30 - 120 min range zhou2023c's FFT is flat (124 -> 116: -0.01 per doubling).
zhou2023c's levels are 4 - 13 times feng2022's at the same times, which fits a system that runs further and passes
its maximum earlier; the paper itself shows deoxypentosones largely consumed in the ARP-GSH model (Fig. 4) and says
organic acids drop the pH rapidly (p. 14317), but prints no thiol loss mechanism.

**Explanation.** (a) a time course read before (feng2022) against after (zhou2023c) a maximum; (e) the
longest-run rule picks the post-peak segment and discards the rising piece; behind the different peak time, (b)
different ARP (ribose-GSH against xylose-alanine), pH 7.0 against 7.5, and much higher thiol levels in zhou2023c.
Measurement is the same method family in both, so it cannot flip the sign.

**Other laboratory.** us10076126 (Nestle/Firmenich, 110 C, 30 - 240 min, deep-eutectic, glycerol or phosphate
media): every FFT ladder rises monotonically to 240 min (contrast file). It sits on the pre-peak side like feng2022;
no peak inside its range at 110 C.

**Would confirm / refute.** A zhou2023c-type run sampled at 15 - 60 min (rise expected), or feng2022 extended past
120 min (fall expected if the peak is a matter of reaction extent). zhou2023c prints no SDs in Fig. 1, so the 60-min
peak cannot be tested against replicate noise from the paper.

**Confidence.** Segment mechanics: high. Why zhou2023c peaks earlier: low - medium.

**Keying.** zhou2023c Fig. 1 values checked against the page image: correct. feng2022 Tables 1-2: correct. Note:
the keyed +0.20 for feng2022 ARP + GSH FFT is not reproduced by an unweighted log-log fit (+0.25); it is
reproduced in direction and roughly in size by an SD-weighted fit dominated by the 90-min arm (SD 0.13). Not a
keying error; flag it as a generator property to check.

## 2. 2-Methyl-3-furanthiol (MFT) x pH, 4.5 - 6 (DFA Garching: hofmann1998 -0.34 and +0.33)

**Designs.** One paper, one method (SIDA with [2H3]-MFT, 0.5 mol/L phosphate, 145 C, 20 min, autoclave), but the
two contrasts come from **different precursor systems**:
- -0.34: cysteine (3.3 mmol) + ribose (10 mmol), Table 2: MFT 55.3, 19.8, 2.5 ug at pH 3, 5, 7.
- +0.33: hydroxyacetaldehyde + mercapto-2-propanone (1 mmol each), Table 8: MFT 15.5, 268.1, 311.5 ug at pH 3, 5, 7.

**Explanation.** (b) The pH effect depends on which step limits MFT formation. From ribose/cysteine the paper
attributes MFT to norfuraneol (NF, formed "in very high yields from ribose/cysteine", Table 5) plus H2S, a route
favoured by acid. From the preformed C2 + C3 fragments the step is an aldol condensation (Fig. 3), favoured by
less acid conditions. The same paper shows a third behaviour: with glucose or rhamnose MFT peaks at pH 5 (Table 2:
0.3, 1.9, 0.4 and 0.1, 0.8, 0.1), not keyed into this cell. So inside one laboratory the sign of "MFT per pH unit"
is a property of the precursor system, and the hydroxyacetaldehyde / mercaptopropanone pots are a mechanistic
intermediate experiment, not a sugar + amino acid reaction. Both are strongly buffered, so (c) does not apply.

**Other laboratory.** cerny2007 (+0.21; xylose + cysteine + thiamin, 0.5 mol/L phosphate, 145 C, 20 min, SPME peak
areas): MFT 415, 336, 371, 539, 391 at pH 4, 5, 5.5, 6, 7 (Table 2), i.e. **flat** within 0.2 dex; endpoints pH 4
-> 7 give -0.01 per pH unit. The keyed +0.21 is the 5 -> 6 monotone run of a flat, noisy ladder (SD stated only as
<= 23 %, about 0.1 dex per arm; the authors warn pH may change headspace partition, text below Table 2). The system also
differs (b): 10 - 25 % of MFT carbon comes from thiamin (Table 3). So cerny2007 differs for a related reason
(another precursor system) plus a segment pick (e) and peak-area measurement (d), not for the reason inside
hofmann1998.

**Would confirm / refute.** A pH ladder of ribose/cysteine run in the cerny2007 set-up with SIDA; if it falls like
hofmann1998, the cerny2007 flatness belongs to thiamin and SPME, not to the laboratory.

**Confidence.** Within-lab explanation: high (the paper states both trends). Cross-lab: medium.

**Keying.** hofmann1998 Tables 2 and 8 and cerny2007 Table 2 checked: correct. Flag for the table, not the keying:
the +0.33 contrast is a precursor-intermediate system pooled into a sugar/amino-acid pH cell.

## 3. MFT x precursor dose, sulfur donor varied (DFA Garching: grosch1992 -0.29, -0.37, -0.35 against +0.30 ... +0.42)

**Designs.** One paper (Grosch and Zeiler-Hilgart 1992), 0.2 mol/L Na-pyrophosphate pH 5.7, ribose 670 umol in 500
mL, reflux 2 h, isotope dilution. The cell pools **two different levers**:
- the four positive contrasts (+0.30 to +0.42, midpoint 18) are **thiamine 6 vs 30 umol** at fixed cysteine 170,
  340, 680, 1360 umol (Table III rows);
- the three negative ones are **cysteine dose** ladders: 30 umol thiamine column 325 -> 2720 umol (Table III: 19,
  14, 10, 4.8, 2.0 ug), 6 umol thiamine column 340 -> 1360 (2.8, 1.6, 0.5), and 0.1 mol/L pyrophosphate 325 -> 680
  (Table IV: 8.8, 7.2, 3.5).

**What the data show.** The full 30-umol-thiamine cysteine ladder is 5.2, 14, 18, 18, 19, 14, 10, 4.8, 2.0 ug at 0,
85, 170, 255, 325, 340, 680, 1360, 2720 umol: **rise then fall**, maximum at 170 - 325 umol. The authors say so:
cysteine "in a concentration range of 170 to 325 umol ... produced the maximum amount of MF-SH, while a higher
concentration was inhibitory" (p. 186). The table keeps the falling run and drops the rise (including the 0-dose
arm). Thiamine is the main MFT carbon source here: ribose + cysteine without thiamine gave only traces (Table V,
text p. 186), so more thiamine always means more MFT.

**Explanation inside the laboratory.** (e) two levers (thiamine amount, cysteine amount) keyed into one "sulfur
donor" cell; (a) the cysteine ladder is non-monotone and only its falling segment survives. Note the 325 and 340
umol arms (19 and 14 ug, adjacent doses, the 325 arm reprinted from Table II) differ by 0.13 dex, a fair gauge of
pot-to-pot scatter.

**Other laboratories.** All positive, and all on rising segments: Zhang2024 cysteine 0.5 -> 62.5 mg/mL rises
(0.45 -> 10.96) and is flat from 62.5 to 75 mg/mL (10.96, 10.92); thomas2014 is a **thiamine** ladder (same lever as
grosch's positive rows); ding2023 rises from 0.25 to 1.0 % cysteine with 0.25 % xylose. The fall in grosch1992 begins
at about 0.5 mol cysteine per mol ribose and about 11 per thiamine; ding2023 still rises at about 5 mol cysteine per
mol xylose (my conversion of 1.0 % and 0.25 % w/v), so a cysteine : sugar ratio alone does not predict the turnover.
A plausible (b): in grosch1992 the MFT precursor is scarce thiamine (60 umol/L) and excess cysteine competes with
or diverts its intermediates, while in ding2023 the carbon comes from xylose in a protein hydrolysate. Not tested.

**Would confirm / refute.** Split the cell by lever (thiamine vs cysteine); keep the rising piece of the grosch
cysteine ladder (0 -> 325 umol is about +0.1 per doubling from 85 umol on, more from 0). A cysteine ladder extended
past 5 mol/mol sugar in a thiamine-poor xylose system would show whether ding2023 also turns over.

**Confidence.** Within-lab (lever pooling + segment): high. Why the turnover appears only in grosch1992: low.

**Keying.** grosch1992 Tables III - IV checked against the text layer: correct. Zhang2024, thomas2014, ding2023
not re-checked here (digitised or previously checked).

## 4. MFT x time (Ho network: feng2022 +0.10, +0.23 against zhou2023c -0.85)

**Designs.** As cell 1.

**What the data show.** feng2022 MFT: ARP alone 4.45, 5.95, 6.62, 7.22; ARP + GSH 11.27, 21.88, 28.63, 32.15 ug/L
(30 - 120 min), monotone rise. zhou2023c ARP + GSH: 59.1, 64.1, 57.1, 8.1 (Fig. 1): essentially flat to 90 min,
then a seven-fold drop at 120 min. Unlike FFT, the MFT fall is not only a segment choice: the 30 -> 120 min endpoint
slope is -0.43 per doubling, and the 30 -> 60 piece is only +0.04.

**Explanation.** (a) zhou2023c's MFT has plateaued by 30 min and collapses late, feng2022's is still rising at 120
min: the two runs sit at different points of reaction extent; (b) different ARP (ribose-GSH against xylose-alanine),
pH and level (zhou2023c 3 - 9 times higher at 30 - 90 min). MFT is the least stable of the thiols (hofmann1998
introduction), so late loss is expected once formation stops, but zhou2023c gives no SD and a single 120-min point
carries the whole fall: (d) cannot be excluded for that one arm. (e) applies less here than for FFT.

**Other laboratory.** us10076126: MFT rises monotonically to 240 min at 110 C in all three media (contrast file);
same side as feng2022.

**Would confirm / refute.** Replicated zhou2023c 90 - 150 min points; bis(2-methyl-3-furyl) disulfide in the same
pots (not reported in Fig. 1) would show whether the drop is oxidation.

**Confidence.** Reaction-extent difference: medium. Mechanism of the 120-min collapse: low.

**Keying.** Correct (Fig. 1 image; feng2022 Tables 1-2).

## 5. Benzaldehyde x precursor dose, sugar varied (A*STAR theng2024 -0.02 ... +0.11; NUS li2021ph -0.34 ... +0.26)

**Designs.** theng2024: nutritional-yeast hydrolysate, protein : sugar 1 : 0.5, 1, 2, initial pH 6.5 unbuffered,
100 C 2 h, HS-SPME, internal standard, "relative concentration (ppm)". li2021ph: enzymatic pork-trimming
hydrolysate, xylose 0, 0.5, 1.5, 2.5, 3.5 g/100 mL, initial pH 5.5 or 4.5 unbuffered (final pH down to 4.86 and
4.35), 90/95/100 C 60 min, **GC-FID raw peak areas, no internal standard** (sec. 2.6).

**theng2024 inside.** The four sugars' benzaldehyde values all carry the same Duncan letter "a" (Table 1:
glucose 0.347, 0.334, 0.319; fructose 0.303, 0.316, 0.308; xylose 0.383, 0.503, 0.506; ribose 0.338, 0.437 +- 0.252,
0.395); the no-sugar NYH control is 0.311 and the unheated-hydrolysate value 0.369. Benzaldehyde is present before
heating and does not respond to sugar. The "sign change" is between effects of -0.02 and +0.11 dex with arm CVs of
10 - 60 %: (f) in the sense of no detectable effect, i.e. noise around zero.

**li2021ph inside.** Six ladders (Table S2). Every one rises from 0 to 0.5 g/100 mL xylose (+0.21 to +0.38 dex),
then wanders: peaks at 0.5 - 2.5 and falls at 3.5 in five of six. The table drops the 0-dose arm and keeps whichever
monotone run is longest, so the keyed slopes are pieces of the upper ladder: +0.26 (pH 5.5, 90 C, 1.5 -> 3.5
rising), -0.10, +0.07, -0.15, -0.34, -0.28. Endpoint 0.5 -> 3.5 slopes are -0.12 to +0.09 per doubling (mean about
-0.04): flat. A measurement caveat (d): the same unheated hydrolysate reads 11.23 (pH 5.5 group) and 5.44 (pH 4.5
group), a 0.31 dex difference with no heating, which says block-to-block scatter in raw FID areas is of the order of
the keyed effects. Furfural rises 27 - 40 fold with xylose (sec. 3.5), so fibre competition at high xylose could
depress benzaldehyde, but the extraction method (Li et al. 2021b) was not read: hypothesis only.

**Explanation.** theng2024: (f) no effect. li2021ph: (a) non-monotone ladder (consistent first-dose rise, flat or
falling after), (e) zero-dose arm dropped and segment picked, (d) uncorrected FID areas with ~0.3 dex block scatter.
Across laboratories there is no conflict: above the first addition, benzaldehyde is insensitive to sugar dose in
theng2024, li2021ph and boonbumrung2023 (0.02 in every arm, at the floor).

**Would confirm / refute.** Keep the 0 -> 0.5 piece as its own contrast (it would be uniformly positive); a
replicate of li2021ph blocks on one day with an internal standard.

**Confidence.** medium - high.

**Keying.** All 30 li2021ph benzaldehyde arms checked against Table S2 and the 12 theng2024 arms against Table 1:
correct.

## 6. Phenylacetaldehyde x pH, 4.5 - 6 (NUS li2021ph: 15 contrasts, -0.59 ... +0.38; other lab zamora2006 -0.05)

**Design.** li2021ph, as cell 5; the pH lever is initial pH 5.5 against 4.5, unbuffered, two arms per contrast,
one contrast per temperature x xylose block.

**What the data show.** The 15 contrasts have mean -0.04 dex and SD 0.29 dex. No pattern with xylose; a weak one with
temperature (90 C: 5/5 negative, mean -0.18; 95 C: mean -0.05; 100 C: 3/5 positive, mean +0.10). Phenylacetaldehyde
is not detected in the unheated hydrolysate, so all of it forms during heating. The scatter is larger than the
printed SDs imply (arm CVs 3 - 25 %, about 0.05 - 0.1 dex per two-arm contrast), and one arm is an outlier: B90-35
reads 25.61 +- 0.89 against 6.66 - 8.18 at neighbouring xylose and temperature, and alone makes the -0.59. The pH
lever is also weaker than labelled: the pH 5.5 pots drop to as low as 4.86 and the pH 4.5 pots to 4.35 (sec. 3.2),
which shrinks the step but cannot flip a sign.

**Explanation.** (d) between-block scatter of raw FID areas without internal standard, larger than the printed SDs,
plus one outlier arm; (c) unbuffered drift (shrinks the effect); possibly (b) a temperature interaction, too weak to
call. The paper itself says pH "might not significantly affect the volatile profiles" (sec. 3.5).

**Other laboratory.** zamora2006 (citrate-buffered, 37 C, lipid-oxidation-driven Strecker of phenylalanine): 5.1 ->
3.6 over pH 3 -> 6, about -0.05 per pH unit. That agrees with li2021ph's mean (-0.04); the sign changes are inside
li2021ph's noise, not a between-laboratory difference.

**Would confirm / refute.** Re-measuring B90-35; a replicate pH pair within one block with an internal standard.

**Confidence.** medium (that it is noise around a small negative or zero effect); low for any temperature
interaction.

**Keying.** All 15 pairs checked against Table S2: correct. Flag B90-35 (25.61) as a suspect value in the paper,
keyed as printed.

## 7. Phenylacetaldehyde x temperature (NUS li2021ph: -0.69 ... +0.75; other labs +0.12 ... +0.58)

**Design.** li2021ph, 90 / 95 / 100 C (5 C steps), 60 min, per pH x xylose block.

**What the data show.** Most ladders are V- or A-shaped at 95 C (e.g. pH 5.5, 0.5 g xylose: 5.31, 4.22, 13.96; pH
4.5, 0.5 g: 6.72, 15.85, 8.54). The longest-run rule, with ties broken toward the first segment, keeps the 90 -> 95
step in seven of ten contrasts, and a 5 C step scaled to 10 C doubles its noise. The **90 -> 100 C endpoints are
positive in 9 of 10 ladders** (+0.04 to +0.43 per 10 C, mean +0.23 without the outlier); the one negative endpoint
(-0.54) is the B90-35 outlier arm of cell 6.

**Explanation.** (e) segment pick on 3-arm zig-zag ladders (the 95 C arm decides the sign); (d) block scatter and
the outlier arm. With endpoints, li2021ph agrees in sign with every other laboratory (monforte2021 +0.58, hidalgo2013
+0.12, zamora2015 +0.23, +0.12).

**Would confirm / refute.** Recompute li2021ph temperature contrasts as 90 vs 100 C or as a fit over all three arms;
the sign conflict should vanish.

**Confidence.** high.

**Keying.** Correct. Two keyed slopes (+0.31 for pH 5.5 / 3.5 g, -0.69 for pH 4.5 / 3.5 g) are not reproduced by an
unweighted fit (+0.43, -0.54); like cell 1 they are consistent with SD weighting. Not a keying error.

## 8. Pyrazines x sugar identity, glucose vs fructose (Minnesota leahy1989 +0.36 and -0.47; Reading martin2001 -0.28)

**Design.** leahy1989: 0.1 M sugar + 0.1 M amino acid in 0.1 M borate pH 9.0 (pH fell with time; samples
readjusted to 9.0 before analysis), 95 C 2 h, purge-and-trap, NPD, 2-methoxypyrazine internal standard. The two
contrasts differ only in the **amino acid**: lysine (glucose 13.1, fructose 5.7 ppm total pyrazines: +0.36) and
asparagine (glucose 0.74, fructose 2.2: -0.47) (Table III, text p. 85-87).

**Explanation.** (b) amino acid x sugar interaction, and partly an aggregation effect of "total pyrazines". The
pyrazine mix differs by amino acid: unsubstituted pyrazine is 56 % / 42 % of the total in lysine-glucose /
lysine-fructose but 19 % / 1 % with asparagine; dimethylpyrazines are 3 % / 9 % with lysine but 22 % / 70 % with
asparagine (text p. 87). Fructose favours the dimethylpyrazines in both amino acids (also in the Table I rate
constants, e.g. asparagine, one dimethylpyrazine: k 0.087 vs 0.487 at 95 C, contrast leahy1989-38), while glucose favours unsubstituted
pyrazine. So the total flips with whichever pyrazine dominates, and that is set by the amino acid (lysine also yields
~18 times more pyrazines). The authors' own reading is reaction extent: fructose reacts faster at first and glucose
overtakes later (Kato et al.), so the faster lysine systems are past the crossover at 2 h and the slow asparagine
systems are not (p. 89) - (a) with an amino-acid-dependent clock. Borate complexes sugars and the pH drifted, but
both apply to both amino acids, so neither explains the reversal.

**Other laboratory.** martin2001 (fried potato slices, asparagine infused, 180 C 2 min): GA 910 vs FA 1730 (Fig. 2),
-0.28, the same sign as leahy1989's asparagine pots. The two laboratories agree for asparagine; the lysine contrast
is the outlier, for the reason above.

**Would confirm / refute.** Per-compound contrasts (not totals) should keep one sign per pyrazine across amino
acids; a lysine time course at 95 C sampled at 15 - 60 min should show fructose ahead early.

**Confidence.** medium (the composition data are in the paper; the reaction-extent account is the authors'
hypothesis).

**Keying.** leahy1989 Table III totals match the text; martin2001 not re-checked (digitised, read off a bar chart).

## Summary

| cell | explanation | confidence |
|---|---|---|
| FFT x time | a, e (b behind the peak timing) | high (segment); low - medium (why) |
| MFT x pH (4.5 - 6) | b (precursor system); cross-lab also d, e | high within; medium across |
| MFT x precursor dose (sulfur donor) | e (two levers pooled), a (rise-then-fall cysteine ladder); b across labs | high within; low (mechanism) |
| MFT x time | a, b; d possible for the single 120-min arm | medium |
| Benzaldehyde x sugar dose | theng2024 f (no effect); li2021ph a, e, d | medium - high |
| Phenylacetaldehyde x pH | d (block scatter, one outlier arm), c (shrinks only); b weak | medium |
| Phenylacetaldehyde x temperature | e, d | high |
| Pyrazines x sugar identity | b (amino acid sets pyrazine mix and reaction extent), a | medium |

Overall pattern:
1. None of the eight within-laboratory sign changes needs a laboratory effect: each comes from a second factor the
   study itself varied (precursor system, amino acid, a second dosed precursor) or from where on a non-monotone curve
   the table's segment falls.
2. The longest-monotone-run rule (and the dropped zero-dose arm) is again the commonest proximate cause (cells 1, 3,
   5, 7); with endpoints or kept pieces, cells 5 - 7 become agreement or noise around zero.
3. Where a laboratory quantifies by raw FID or SPME peak area without an internal standard (li2021ph), block-to-block
   scatter (~0.3 dex) reaches the size of the effects, so its two-arm contrasts carry more than the 0.01 - 0.03 dex
   repeatability; no keying errors were found in any arm checked.
