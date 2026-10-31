# Zyzak, Richardson, Thorpe & Baynes 1995 — EXTRACTION
### *"Formation of Reactive Intermediates from Amadori Compounds under Physiological Conditions"*
### Model Amadori compound Nα-formyl-Nε-fructoselysine (fFL, 0.4 or 5.0 mM) and glycated rat-tail collagen, **pH 7.4, 37 °C, days to weeks**; phosphate 10 / 50 / 100 / 200 mM vs 0.2 M HEPES; air vs nitrogen + chelators. Half-lives (Table I), product yields (Table II: tetroses, pentoses, hexoses, 3-DG). **One pH and one temperature only** -- the abstract's "and/or pH" is not supported by any pH series in the paper.

**Read 2026-10-31 from `data/articles/zyzak1995.pdf`** (8 pp., Arch. Biochem. Biophys. 316:547-554). **Scanned page
images, NO text layer** (pdftotext returns nothing): everything below was read from 130 / 180 dpi page renders, Tables
I and II from 300 dpi crops (scratchpad `dos3/zzT1-4.png`, `zzT2-6.png`).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Question (hub / cysteine-fate rivals, `docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 1 findings 3 and 10;
`results/validation/kinetic_core_dpo_partition_prereg.md` sec. 13 CORRECTION and sec. 14 H4 / H6): is the Amadori
compound's decomposition buffer-catalysed (phosphate vs HEPES; phosphate-concentration law), how much does air
accelerate it, what pH / temperature dependence, and what products with what yields?

| need | what Zyzak 1995 gives | tag |
|---|---|---|
| (a) phosphate-concentration law (N2 + chelators) | fFL t1/2 **57.1 +- 2.0 / 40.9 +- 2.3 / 24.4 +- 2.1 / 16.8 +- 0.6 d** at 10 / 50 / 100 / 200 mM phosphate (triplicates) -> k = 0.0121 / 0.0169 / 0.0284 / 0.0413 d^-1; linear fit **k = 0.0105 + 0.157 [Pi]/M d^-1**; k(200) / k(10) = 3.4; log-log slope 0.41 | [M] / [Z] |
| (b) the same under air (no chelators) | t1/2 **10.7 +- 0.2 / 9.8 +- 1.2 / 5.8 +- 0.1 / 4.2 +- 0.2 d** -> k = 0.065 / 0.071 / 0.120 / 0.165 d^-1; k = 0.055 + 0.56 [Pi]/M d^-1; k(200) / k(10) = 2.5 | [M] / [Z] |
| (c) phosphate vs HEPES (0.2 M each) | HEPES t1/2 **~200 d (N2), ~56 d (air)** (Fig. 2B labels, n not stated) vs 16.8 / 4.2 d in phosphate: phosphate **~12x (N2) and ~13x (air)** faster. HEPES (k 0.0035 d^-1, N2) is **3x below** the phosphate law's zero-phosphate intercept (0.0105) | [M] / [F] / [Z] |
| (d) air vs N2 | "~4 times faster under oxidative compared to anti-oxidative conditions at identical phosphate concentrations" [V]; Table I ratios 5.3 / 4.2 / 4.2 / 4.0 at 10-200 mM. **Confound:** N2 pots also contain 1 mM phytic acid + 1 mM DTPA, air pots none | [M] / [!] |
| (e) pH | **None.** pH 7.4 throughout; the abstract's "accelerated by increasing the phosphate concentration and/or pH" has no pH data behind it in the text, tables or figures | [!] |
| (f) temperature | **None.** 37 °C throughout | — |
| (g) products (% of initial fFL; Table II) | **tetroses 19-31 %, hexoses (glucose + mannose, reversal) 7.5-22 %, pentoses 1.8-5.4 %, 3-DG 1.0-1.6 %** at ~half-life; 3-DG maximal "up to 5% yield within 1 week", then declines [Q]; glucosone not detected (its own t1/2 "<=12 h" under both conditions [M]); 1-DG tentatively < 1 % (L-rhamnitol as model standard) | [M] / [Q] |
| (h) product distribution vs phosphate | "phosphate catalyzes the decomposition reaction, but does not significantly affect the product distribution" [V]; Table II, N2: 10 mM / day 49 vs 200 mM / day 13 at similar loss (40.5 vs 48.2 %): tetroses 23.6 vs 19.0, hexoses 21.8 vs 19.9, 3-DG 1.1 vs 1.6 % | [M] / [Q] |

**Answer.** At 37 °C and pH 7.4, **phosphate is a strong, concentration-dependent catalyst of Amadori-compound
decomposition, and HEPES is not** (confidence ~90 % in the direction, ~80 % that the size, x3.4 from 10 to 200 mM and
x12 vs 0.2 M HEPES, is right within a factor 1.5; one research group, one compound, triplicates). The catalysis acts on
all channels together (reversal, rearrangement to 3-DG, retro-aldol), so it sits on a common early step (enolisation /
ring opening) [Z, consistent with the authors' "increases the proportion of the open chain" and "promotes the
enolization" [C, refs 22-23]]. Air adds a ~4x faster oxidative channel (CML, LL; and the product spectrum shifts to
fewer hexoses). Weighed for the model:

1. **Phosphate vs water is a live, large rival for the Davidek (0.2 M phosphate) / Zhai (water) tension** (dpo prereg
   sec. 13 CORRECTION). If the 37 °C law transfers, an Amadori compound in 0.2 M phosphate decomposes ~4-12x faster
   than in water (water is closer to HEPES or to the zero-phosphate intercept; which one is unknown). Transfer to
   100-120 °C and to a glycine-xylose ARP is **unmeasured**: P(phosphate still >= 2x at 100-120 °C) ~0.6 [Z, my
   judgement: general-base catalysis usually persists but the uncatalysed path gains at high T].
2. **The law is not first-order in phosphate.** Linear in [Pi] with a large intercept, or ~0.4 order in a log-log fit;
   and HEPES lies below the intercept. Either HEPES actively slows the reaction or the low-phosphate end is curved.
   A `k = k0 + kP [Pi]` form (mechanism review finding 3) fits the four points (k0 / kP = 0.067 M at 37 °C under N2,
   i.e. 0.2 M phosphate multiplies the uncatalysed rate by ~4, by ~12 against HEPES) [Z].
3. **3-DG is a minor, transient product at 37 °C** (~1-1.6 % at the half-life; peak <= 5 %), glucosone does not
   accumulate, and Amadori **reversal** to the aldoses + amine is a major route (hexoses ~20 % under N2). Reversal is
   not in a typical hot-pot Maillard model and is relevant to the TTCA question: a reversed ARP releases free aldose
   that cysteine traps as a thiazolidine (see wei2022_extraction.md, sec. 1 point 2).
4. **O2 (hub H6):** air accelerates decomposition ~4x, but the products under air are CML / LL (oxidative cleavage)
   plus more tetroses and pentoses, and **no glucosone** (detected only with added Cu2+ in Kawakishi 1991 [C]). This
   is weak evidence AGAINST H6's "O2 diverts the Amadori compound to the osone" at physiological trace-metal levels
   (P(H6) unchanged-to-lower; the air / chelator confound applies).

## 2. Identity

| field | value as printed |
|---|---|
| Authors | David V. Zyzak, John M. Richardson, Suzanne R. Thorpe, John W. Baynes* (University of South Carolina) |
| Title | as in the header |
| Venue | *Arch. Biochem. Biophys.* **1995**, 316 (1), 547-554 (January 10) |
| DOI | **10.1006/abbi.1995.1073** -- NOT printed on the PDF (1995 scan); Crossref 2026-10-31: same title, authors, vol. 316, issue 1, pp. 547-554 |
| Dates | received 8 Aug 1994 |
| Funding | NIH NIDDK DK-19971 |
| Content | Table I (half-lives vs phosphate, N2 / air), Table II (product yields); Figs. 1 (amino-acid chromatograms), 2 (fFL decay, phosphate vs HEPES, air / N2), 3 (collagen FL decay), 4 (GC-MS SIM of sugars), 5 (product kinetics, 0.2 M phosphate, N2), 6 (pathway scheme) |

## 3. Methods [M]

| item | as printed |
|---|---|
| fFL | synthesised from Nα-formyllysine + glucose (Finot & Mauron 1969); 0.4 or 5.0 mM (which concentration per Table I row is not stated [!]) |
| Aerobic pots | under air, 0.2 M phosphate pH 7.4 (and the phosphate series), dark, 37 °C, Teflon-lined screw-cap tubes, 2 drops toluene (antibacterial) |
| Anaerobic pots | + phytic acid and DTPA 1 mM each, gassed 5 min with dithionite-scrubbed N2, sealed under N2 |
| HEPES | 0.2 M, pH 7.4, air and N2 (Fig. 2B) |
| Collagen | rat-tail collagen glycated 4 weeks in 0.2 M phosphate + 250 mM glucose (anti-oxidative), washed; 50 mg wet / 4 mL buffer (10 mM or 0.2 M phosphate), N2; FL as furosine by GC-MS |
| Analytics | fFL as Nε-mannitol- + glucitol-lysine after NaBH4 reduction and deformylation (HPLC amino-acid analyser, OPA fluorescence); sugars as alditol acetates (NaBH4 / NaBD4), GC-MS SIM (m/z 145 tetroses / pentoses / hexoses, 231 3-DG, 210 myo-inositol IS) |
| Discrimination | alditol acetates cannot tell aldose from ketose or glucose from glucosone; NaBD4 labelling: > 95 % of tetrose / pentose / hexose peaks are monocarbonyl sugars; tetroses and pentoses are aldose + ketose mixtures, hexoses exclusively aldoses (O-methyloxime check) |

## 4. Results

**Table I (t1/2, days; fFL means +- SD of triplicates; collagen mean +- range of duplicates or single) [M]:**

| | phosphate (mM) | nitrogen | air |
|---|---|---|---|
| fFL | 10 | 57.1 +- 2.0 | 10.7 +- 0.2 |
| fFL | 50 | 40.9 +- 2.3 | 9.8 +- 1.2 |
| fFL | 100 | 24.4 +- 2.1 | 5.8 +- 0.1 |
| fFL | 200 | 16.8 +- 0.6 | 4.2 +- 0.2 |
| glycated collagen | 10 | 75.1 +- 3.2 | nd |
| glycated collagen | 200 | 20.0 (single) | nd |

**Fig. 2 (0.2 M buffers) [F / M labels]:** phosphate t1/2 ~5 d (air, composite of three experiments, 95 % CI shown),
~17 d (N2); HEPES ~56 d (air), ~200 d (N2). First-order (log-linear) decay over ~30 d.

**Table II (% of initial FL) [M]:**

| condition | tetritols | pentitols | hexitols | 3-DG | % loss of FL |
|---|---|---|---|---|---|
| fFL 10 mM Pi / N2 / day 49 | 23.6 +- 1.4 | 2.4 +- 0.2 | 21.8 +- 1.1 | 1.1 +- 0.2 | 40.5 +- 1.4 |
| fFL 200 mM Pi / N2 / day 13 | 19.0 +- 0.2 | 1.8 +- 0.1 | 19.9 +- 1.2 | 1.6 +- 0.1 | 48.2 +- 1.7 |
| fFL 10 mM Pi / air / day 13 | 30.6 +- 6.4 | 5.4 +- 0.6 | 12.4 +- 1.2 | 1.0 +- 0.4 | 61.3 +- 1.1 |
| fFL 200 mM Pi / air / day 5 | 19.7 +- 4.7 | 4.8 +- 1.9 | 7.5 +- 4.3 | 1.1 +- 0.5 | 58.4 +- 2.8 |
| collagen 10 mM Pi / day 35 | 7.3 | 1.6 | 16.7 | 1.5 | 28 |
| collagen 200 mM Pi / day 28 | 22.5 | 2.1 | 17.1 | 0.6 | 58 |

[!] In one row the products exceed the loss (10 mM / N2: 48.9 % products vs 40.5 % loss; the other rows sum to
57-97 % of the loss): either the % loss is low there or a fFL can yield more than one counted sugar; not discussed.

**Fig. 5 (fFL, 0.2 M phosphate, N2; % of initial fFL) [F]** (points at ~0, 5, 13, 29, 49 d): fFL 100 / ~70 / ~54 /
~25 / ~12; total sugars 0 / ~25 / ~44 / ~56 / ~60; lysine 0 / ~8 / ~16 / ~37 / ~48; tetritols 0 / ~8 / ~20 / ~28 /
~21; hexitols 0 / ~15 / ~22 / ~24 / ~31; pentitols 0 / ~1 / ~2 / ~4 / ~6; 3-DG 0 / ~2.5 / ~2 / ~0.5 / ~0. Text:
sugars or lysine "typically accounted for >75% of fFL decomposition" [V]; gradual browning signals other products.

**Mechanism (Fig. 6) [P]:** under +-O2, reversal (mannose, glucose, lysine), rearrangement (3-DG), reverse aldol
(tetroses, tetruloses, pentoses, pentuloses); under +O2 also CML, LL and (proposed) glucosone. Shown: products and
their yields; the reverse-aldol and double-bond-migration steps are inferred from the products, not observed.
Phosphate's role (acid-base catalyst, pKa near 7) is argued from refs 21-23 [C] (Watkins et al. 1987 J. Biol. Chem.
262:7202, phosphate catalysis of glycation; Burton & McWeeny 1963; Thornalley et al. 1984).

## 5. Supplementary material

**None.** A 1995 ABB article: no Supporting Information is mentioned anywhere in the paper. The aldose / ketose
discrimination data "will be published separately in detail (D. V. Zyzak, Ph.D. thesis, in preparation)" [V] -- not
supplementary material, and not needed for the rate questions above.

## 6. Flags

- [!] **No pH or temperature series**; the abstract's pH claim is unsupported in the body.
- [!] Air vs N2 is also chelators-absent vs chelators-present (1 mM phytate + 1 mM DTPA): the ~4x "oxidative" factor
  includes any trace-metal catalysis.
- [!] HEPES half-lives are figure labels (~56, ~200 d) without SD or n; HEPES only at 0.2 M (no HEPES concentration
  series), so "HEPES is inert" is not tested, only "HEPES is far slower than phosphate".
- [!] fFL concentration per experiment (0.4 vs 5.0 mM) not stated per row.
- [!] Fig. 3's caption assigns open circles to 10 mM and filled to 0.2 M phosphate; the legend and text say the
  opposite (open = 0.20 M, t1/2 ~20 d). This dossier follows the legend and text.
- [!] Hexose secondary-amine Amadori compound (fructoselysine), 37 °C, days: transfer to a glycine-xylose ARP at
  100-120 °C (minutes) is an assumption, not evidence.
