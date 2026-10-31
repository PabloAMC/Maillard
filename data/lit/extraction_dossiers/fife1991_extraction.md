# Fife, Natarajan, Shen & Bembi 1991 — EXTRACTION (thiazolidine RING OPENING to an iminium ion and its HYDROLYSIS to aldehyde + amino thiol measured as SEPARATE steps; full pH-rate profiles pH 1-12; general-ACID catalysis of ring opening (acetic acid) and general-BASE catalysis of the hydrolysis (acetate, cacodylate, PHOSPHATE, imidazole, N-ethylmorpholine, carbonate; Brønsted β 0.29); 30 / 50 / 80 / 90 °C, μ 0.5 M KCl. Substrates are 2-aryl/styryl 1,3-thiazolidines from p-(dimethylamino)cinnamaldehyde, cinnamaldehyde and p-(dimethylamino)benzaldehyde with 2-aminoethanethiols — NO cysteine, NO 4-carboxylate, NO sugar)

### THE MECHANISTIC PAPER THE ENGINE'S LUMP HIDES. Ring opening (C-S cleavage, iminium + thiol) is FAST and reversible; the slow, rate-determining step is the water / OH⁻ attack on the iminium, and THAT step is general-base catalysed by buffer bases including phosphate. No constant transfers to TTCA (wrong aldehydes, wrong amine, wrong ring substituents), but the rate-law FORM k = k0(pH) + Σ k_B [B] and the size of the buffer term relative to k0 are measured.

**Source on disk:** `data/articles/fife1991.pdf` (9 pp.; Acrobat Capture OCR of a scan; p. 1 starts with the end of the
preceding osazone paper, p. 9 with the following cellulose paper). Prose from the `pdftotext -layout` layer; **Tables
I-II, eqs 1-4, 9, 11, 13 and Figs. 1-6 read from 110 / 220 dpi page renders** (scratchpad `pdf3/fife-*.png`,
`fifehi-*.png`). The OCR layer garbles every exponent in the prose (e.g. "5.3 X IO-' M"); all exponents below were
read from the renders.

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[F]** read from a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed / hypothesised, not
measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question (plan step B7b; `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 12;
`docs/history/mechanism_reviews_2026-10-31/cysteine.md` H1, P 0.55): the engine's `k_ttca_cys` (~0.028 min⁻¹, fitted,
no pH or catalyst term; forward rate by detailed balance with Cao's K ~60 L/mol) lumps ring opening and hydrolysis.
Why are dilute (10 mM) fed-TTCA pots in water far from equilibrium while 0.4 M pots in 0.2 M phosphate equilibrate
in 30 min?

| what the model needs | what Fife 1991 gives | tag |
|---|---|---|
| ring opening and hydrolysis as separate steps | **yes, both measured.** Ring opening (iminium appearance, 480 / 525 nm, stopped flow) and hydrolysis (aldehyde appearance at 400 nm = iminium decay at 480 nm) are separate first-order processes; "the rate-determining step in the overall reaction must be hydrolysis of the iminium ion intermediate at all pH values" (for I-III) | [M] / [Q] |
| ring-opening rate constants | II (N-butyl), 30 °C: k₁' = **1.4 s⁻¹** (pH-independent), k_H = **7 × 10⁵ M⁻¹ s⁻¹** (H₃O⁺), K_a = 3 × 10⁻⁵ M (eq 1, fit pH 3.5-6.5). III (N-phenyl), 30 °C: k₁' = **5.7 × 10⁻³ s⁻¹** pH 4-10; k(H₂O)/k(D₂O) = 2.25 (pD 6.15 / 6.65: 2.51 / 2.55 × 10⁻³ s⁻¹); II opens 250× faster than III | [M] |
| hydrolysis rate constants (zero buffer) | Table I below (50 °C; IV also 80 °C; V 90 °C) | [M] |
| pH profile | full log k vs pH 1-12 for I, II (Fig. 1), III (Fig. 2), IV (Fig. 4), V (Fig. 5); up to five slope changes; apparent pK 6.3 (II) / 7.1 (III) is K_eq, a composite of ring opening + protonation | [M] / [F] |
| buffer / general catalysis | **ring opening:** general ACID, acetic acid k_HA = **0.046 M⁻¹ s⁻¹** (III, 30 °C, pH 4.1-5.1, Fig. 3); none with Tris at pH 8. **hydrolysis:** general BASE, Table II (I: phosphate **0.050 M⁻¹ s⁻¹**); Brønsted β = **0.29** (II, Fig. 6); **no** buffer catalysis for the N-acetyl V | [M] |
| temperature | no Ea. Only IV at two temperatures (k_OH'' 0.057 → 0.36 M⁻¹ s⁻¹, 50 → 80 °C) | [M] |
| which thiazolidines | I 2-(p-(dimethylamino)styryl)-1,3-thiazolidine (NH); II its N-butyl; III its N-phenyl; IV 2-styryl-N-butyl (cinnamaldehyde); V 2-(p-(dimethylamino)phenyl)-N-acetyl. All from 2-aminoethanethiols — **no cysteine / 4-carboxylate, no aliphatic or sugar aldehyde** | [M] |
| equilibrium constant with free thiol + aldehyde | **none.** Only the ring-chain constants: K_eq = **5.3 × 10⁻⁷ M** (II, ring + H₃O⁺ ⇌ iminium-SH), K₂ = 1.2 × 10⁻¹⁰ M (thiol, pK₂ 9.9), K₀' = 2.3 × 10⁻⁴ (zwitterion; "rapid reclosure of the thiazolidine ring via attack of the thiol anion") | [M] |

**Answers.**

1. **Can any constant anchor `k_ttca_cys`? No** (confidence ~95 %). Every substrate carries a conjugated aryl / styryl
   group at C-2 that stabilises the iminium by orders of magnitude relative to TTCA's polyhydroxyalkyl C-2, none has
   the 4-carboxylate, and the temperatures are 30-90 °C, not 90-145 °C. The authors themselves show the regime
   flips with iminium stability: I-IV hydrolyse with fast pre-equilibrium opening; V (carbocation destabilised by
   N-acetyl) has rate-determining, H₃O⁺-catalysed opening with NO buffer catalysis. Which regime TTCA is in is not
   decided by this paper [Q/Z].
2. **The engine's lump has the wrong shape, in the direction H1 needs** [Q]. "Ring opening ≠ hydrolysis": the fast
   C-S ring-chain step (seconds at 30 °C here; minutes-hours for the sugar / penicillamine thiazolidines in
   Roberts 1987, Nagasawa 1981) is NOT the slow release of free thiol + aldehyde. The overall release is k_hyd ×
   (fraction open), and the fraction open falls at high pH (K₀' 2.3 × 10⁻⁴; "ring opening is only favorable at pH
   <10 where the thiol group is un-ionized"). This resolves the cysteine review's "against" column: fast
   tautomerism is compatible with slow hydrolysis.
3. **Rate law form for B7b: k = k0(pH) + Σ k_B,i [B_i]**, linear in each buffer base, with k_B following a Brønsted
   line of small slope (β 0.29: a 3-unit pKa change moves k_B only ×7.4) [M/Z]. So ANY base present at 0.1-0.4 M —
   phosphate, carboxylates, amines — contributes; the catalyst is the BASE form, so phosphate's activity at pH 5.5 is
   set by its small HPO₄²⁻ fraction.
4. **How big is the phosphate term? Modest at pH 5.5, large at pH 7** [Z, transfer-laden]. HPO₄²⁻ fraction at pH 5.5
   with the printed pKa 7.05: 0.027, so 0.2 M phosphate gives [B] ≈ 5.5 mM.
   - Compound I, 50 °C: k_B·[B] = 0.050 × 0.0055 = 2.7 × 10⁻⁴ s⁻¹ against k0(pH 5.5) ≈ 6 × 10⁻⁴ s⁻¹ → **× 1.5**.
     At pH 7 (0.094 M HPO₄²⁻): 4.7 × 10⁻³ s⁻¹ → **× 9**.
   - Compound II (phosphate not run; Brønsted line through acetate, k_B ≈ 1.2 × 10⁻³ M⁻¹ s⁻¹ at pKa 7.05):
     6.6 × 10⁻⁶ s⁻¹ against k0(pH 5.5) ≈ 6 × 10⁻⁵ s⁻¹ [F] → **× 1.1**.
   - Cao vs Zhai needs ≥ ~4× at face value (Cao ≥ 0.013 min⁻¹ at 90 °C vs Zhai ≤ 0.0033 min⁻¹ at 100 °C), ≥ ~8-10×
     once the 10 °C is allowed for. **On Fife's numbers phosphate GENERAL-BASE catalysis at pH 5.5 alone gives
     1.1-1.5×, not 10×** — it would have to be ~10× more effective relative to the water reaction for TTCA than for
     these aryl thiazolidines. Mild evidence against phosphate as the whole explanation (consistent with the review's
     split: phosphate 0.10 of H1's 0.55). It says nothing about amino-acid catalysis by the 0.4 M reactants
     themselves beyond the Brønsted slope (a carboxylate of pKa ~2.3 — assumed, not in the paper — would be ~0.2×
     acetate per mol).
5. **Temperature:** IV's apparent OH⁻ term gives Ea ≈ **58 kJ/mol** [Z: 0.057 → 0.36 M⁻¹ s⁻¹ over 323 → 353 K; it
   folds in the T-dependence of K_eq]. The only temperature information; not a TTCA Ea.

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanism of Thiazolidine Hydrolysis. Ring Opening and Hydrolysis of 1,3-Thiazolidine Derivatives of p-(Dimethylamino)cinnamaldehyde" |
| Authors | Thomas H. Fife, R. Natarajan, C. C. Shen, Ramesh Bembi — Dept. of Biochemistry, University of Southern California, Los Angeles |
| Venue | *J. Am. Chem. Soc.* 1991, **113**(8), 3071-3079; received June 5, 1990. Code 0002-7863/91/1513-3071$02.50/0 |
| DOI | 10.1021/ja00008a041 (not printed; Crossref returns this title, these four authors, JACS 113(8) 3071-3079, 1991; checked 2026-10-31). **Matches the brief.** |
| Funding | NIH, NSF |
| Paper type | Physical-organic kinetics: UV-vis and stopped-flow spectrophotometry, pseudo-first-order, extrapolation to zero buffer |
| Related | Fife & Natarajan 1986 JACS 108:2425 (O,S-thioacetal ring opening, ref 10); Pesek & Frost 1975 Tetrahedron 31:907 (ref 28); Luhowy & Meneghini 1979 JACS 101:420 (alkaline thiazolidine hydrolysis, ref 29a) — none on disk |

## 2. Methods as they matter to a model

- **Substrates** [M]: I-IV by refluxing the aldehyde (p-(dimethylamino)cinnamaldehyde or cinnamaldehyde) with
  2-aminoethanethiol, 2-(butylamino)ethanethiol or 2-(phenylamino)ethanethiol in dry benzene (Dean-Stark); V by
  N-acetylation of 2-(p-(dimethylamino)phenyl)-1,3-thiazolidine.
- **Kinetics** [M]: 15 µL acetonitrile stock into 3 mL buffer (substrate µM-scale; concentration not printed);
  pseudo-first-order ≥ 4 half-lives; ionic strength 0.5 M with KCl throughout; k0 by extrapolation of k_obsd to zero
  buffer; K_w 5.5 × 10⁻¹⁴ (50 °C), 2.5 × 10⁻¹³ (80 °C). Durrum D-110 stopped flow for ring opening.
- **Detection** [M]: aldehyde at 400 nm (I-III), 300 nm (IV), 360 nm (V); iminium at 480 nm (I, II) / 525 nm (III);
  protonated aldehyde at 280 nm below pH 2.

## 3. Numbers printed

**Table I** (hydrolysis = aldehyde formation, H₂O, 50 °C, μ 0.5 M KCl) [M]:

| compd | k₁, s⁻¹ | k₂, s⁻¹ | k_OH', M⁻¹ s⁻¹ | k_OH'', M⁻¹ s⁻¹ | pK_app,1 | pK_app,2 |
|---|---:|---:|---:|---:|---:|---:|
| I | 3.6 × 10⁻³ | 6.0 × 10⁻⁴ | 13000 | 0.12 | 3.0 | 6.2 |
| II | 1.7 × 10⁻⁴ | 3.5 × 10⁻⁵ | 9100 | 2.0 | 1.8 | 6.3 |
| III | | 7.6 × 10⁻⁴ | 170 | 0.014 | | 7.1 |
| IV | | | | 0.057 | 5.0 | |
| IV (80 °C) | | | | 0.36 | 4.5 | |
| V (90 °C) | 1.0 × 10⁻³ (= k_H K_a) | | | | 4.5 | |

k₁: water attack on the dicationic (protonated dimethylamino) iminium at low pH; k₂: water on the cationic iminium
(un-ionised SH); k_OH': OH⁻ on the cationic iminium (kinetically = water on the zwitterion); k_OH'': OH⁻ at high pH.
V: k_obsd = k_H a_H K_a / (K_a + a_H), k_H = 26 M⁻¹ s⁻¹, pK_a 4.5 (eq 2); k_obsd = 0.001 s⁻¹ pH 1-4, slope -1 above.

**Table II** (general-base catalysis of hydrolysis, 50 °C, μ 0.5 M) [M], 10³ k_B in M⁻¹ s⁻¹:

| compd | acetate (pKa 4.70) | cacodylate (6.20) | phosphate (7.05) | imidazole (6.80) | N-ethylmorpholine (7.80) | carbonate (10.20) |
|---|---:|---:|---:|---:|---:|---:|
| I | 5.6 | 23 | 50 | — | — | — |
| II | 0.25 | 1.9 | — | 1.25 | 2.30 | 9.8 |

"the observed catalytic effect of the buffer in the hydrolysis of I-IV was due to the base species" [V].

**Other printed constants** [M]: II ring-opening eq 1 at 30 °C (k₁' 1.4 s⁻¹, k_H 7 × 10⁵ M⁻¹ s⁻¹, K_a 3 × 10⁻⁵ M);
III k₁' 5.7 × 10⁻³ s⁻¹, k_HA(acetic acid) 0.046 M⁻¹ s⁻¹; K_eq (II) 5.3 × 10⁻⁷ M; K₂ 1.2 × 10⁻¹⁰ M; K₀' 2.3 × 10⁻⁴;
k₃ (water on the zwitterion of II, calculated) 4 s⁻¹ at 50 °C; IV pK_a 4.6 by spectrophotometric titration (30 °C).

## 4. Figures (read from renders, ± 0.1 log unit) [F]

- Fig. 1 (50 °C): I — log k0 ≈ -2.4 (pH 1-2), ≈ -3.2 (pH 4-9), OH⁻ rise above 9. II — ≈ -3.8 (pH 1-2.5), ≈ -4.45
  (pH 3.5-5), rise through pH 6.3 to ≈ -3.0 (pH 7-10), rise above 10.
- Fig. 2: III ring opening (30 °C) flat at ≈ -2.25 pH 4-11, steep below 4; III iminium hydrolysis (50 °C) ≈ -3.1 (pH
  3.5-6.5), ≈ -4.1 (pH 8-11), rise above 11.
- Fig. 3: k_obsd linear in total acetate 0.1-0.5 M at pH 4.10 / 4.65 / 5.10 (slopes fall with pH: acid catalysis).
- Fig. 4 (IV, 80 °C): ≈ -2.55 (k0 ≈ 2.8 × 10⁻³ s⁻¹ ≈ 0.17 min⁻¹) pH 5-10, ≈ -4.0 at pH 3, rise above 11.
- Fig. 5 (V, 90 °C): -3.0 pH 1-4, slope -1 above pH 4.5 (≈ -4.0 at pH 5.4, ≈ -4.8 at pH 6.1).
- Fig. 6: Brønsted plot for II, slope 0.29.

Read once, not two-reader keyed; nothing here is a FIT row.

## 5. What the repository takes

No constant enters the engine (wrong substrates). Context for the B7b design (`cysteine.md` H1), never a prior:
(i) thiazolidine → thiol + aldehyde is two steps, fast reversible C-S ring opening to an iminium and a slower,
rate-determining iminium hydrolysis [M/Q]; (ii) the hydrolysis is general-base catalysed, linear in base
concentration, Brønsted β 0.29, phosphate the strongest catalyst measured for I [M]; (iii) the ring-opening step is
general-acid catalysed only near the H₃O⁺ / water crossover (pH ~4.5) [M/Q]; (iv) on these substrates 0.2 M
phosphate at pH 5.5 adds only 10-50 % to k0 [Z] — weigh this against the phosphate share of H1.

## 6. Supplementary material

**None for this paper.** The "Supplementary Material Available" paragraph on p. 3071 belongs to the PRECEDING
osazone paper (crystallographic tables). Fife et al. print no supplementary statement; all kinetic data are in
Tables I-II and Figs. 1-6. Verdict: nothing held elsewhere.

## 7. Flags

1. **Analogue only.** Conjugated aryl / styryl C-2, 2-aminoethanethiol ring, no carboxylate — the iminium is far more
   stable than TTCA's; the regime (pre-equilibrium vs rate-determining opening) can differ (V shows it flips).
2. **μ = 0.5 M KCl in every run**: salt effects, relevant to the 0.4 M vs 10 mM question, were never varied.
3. **Substrate concentration not printed** (µM-scale by the injection volumes [Z]); no bimolecular re-formation step
   from free thiol + aldehyde was studied.
4. **OCR exponents are wrong throughout the text layer**; use the renders (Tables I-II are clean).
5. **Fig. 2 / Fig. 4 slopes** are drawn with theoretical unit slopes; the reading in sec. 4 is of plateaus.
6. **Registry:** none of I-V is a compound key (not needed).
