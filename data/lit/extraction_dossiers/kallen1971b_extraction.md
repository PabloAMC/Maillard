# Kallen 1971b — EXTRACTION (L-cysteine + formaldehyde -> thiazolidine-4-carboxylic acid (TC): full pH-rate profile of FORMATION, pH 1-13, 25 °C, ionic strength 1.0 M (KCl); rate-determining step changes at pH ~6, from amine attack on formaldehyde (acid side, NO buffer catalysis) to carbinolamine dehydration to the iminium (alkaline side, general-acid catalysed, Brønsted α 0.66); phosphate catalysis SATURATES; morpholine nucleophilic catalysis; hemithioacetal pre-equilibrium. Formation rates only; the reverse rate follows from the companion equilibria paper)

### THE RATES THAT kallen1971_extraction.md SAID WERE MISSING. Combined with Kallen 1971's K_app (same lab, same 25 °C / I = 1.0 M), they give the TC hydrolysis constant k_r ≈ 1 × 10⁻⁶ s⁻¹ (pH 5.6-6.6), about 120-270x SLOWER than the relaxation-based estimate in kallen1971_extraction.md arithmetic 3, and a k_r(pH) that RISES with pH from 2.8 to 5.6. Formaldehyde, 25 °C only: a mechanism and pH-shape analogue for TTCA, not a constant.

**Read 2026-11-01 from `data/articles/kallen1971b.pdf`** (13 pp.: JACS pp. 6236-6248; the last page continues into an
unrelated Brossi et al. alkaloid paper). Acrobat Capture OCR. Prose from the `pdftotext -layout` layer. **Tables I-V,
Schemes I-II, eqs 2-10 and Figs. 1-8 checked against 140 dpi page renders**, Figs. 6 and 7 against 300 dpi crops
(scratchpad `dos6/kal-0*.png`, `k6-05.png`, `k7-06.png`). The OCR layer mis-reads several Table I entries (pH 4.47 as
"4.41", 4.77 as "4.71", k at pH 5.00 as "0.21", K1(app) 17.8 as "11.8"); every number below is from the renders.
Figure values read by me: ±0.1 log unit on Figs. 3 / 6, ±0.1 × 10³ on Fig. 7.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the author · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Questions: CF-076 (dilute vs concentrated Cys + Xyl ⇌ TTCA exchange), the E1 / E2 / E4 rivals of
`docs/history/cys_exchange_2026-10-31/` (prereg bands: E1 k0, c_aa, c_P; E2 / E4 acid-limb slope `a`), the engine's
`k_ttca_cys` (0.0286 min⁻¹ at 100 °C, Ea 64.1 kJ/mol, pH- and solute-flat).

| model need | what Kallen 1971b gives | tag |
|---|---|---|
| forward k_f, pH law | second-order k_f (total cysteine, formaldehyde hydrate F) from pH 1.00 (0.0013 M⁻¹ s⁻¹) to 9.41 (1.50), peak ~4 M⁻¹ s⁻¹ near pH 7-8.3 (Table I, Fig. 3); per amine-free-base cysteine: k1a 1.11 × 10³, k1b 2.92 × 10³ M⁻¹ s⁻¹ (Table II) | [M] / [F] |
| reverse k_r | **not printed.** With Kallen 1971's K_app: k_r ≈ 9 × 10⁻⁷ s⁻¹ (pH 5.57), 9.3 × 10⁻⁷ (6.56), 3.2 × 10⁻⁷ (7.45); 2.3 × 10⁻⁷ (4.47); 3 × 10⁻⁸ (2.84) at 25 °C | [Z] |
| K | not here (companion paper); K1 for the N-hydroxymethyl carbinolamine 26.8 ± 4.0 M⁻¹ (Table II) | [M] |
| rate-determining step | amine attack on unhydrated formaldehyde below pH ~6; acid- and water-catalysed carbinolamine dehydration to the cationic imine above; ring closure fast at all pH ("does not appear to become rate determining") | [Q] / [M] |
| buffer catalysis | **none** by cyanoacetate, formate or acetate, 0.05-0.4 M, pH 2.3-5.0 [M]. General ACID catalysis of dehydration above pH 5.5 (Table III: H₃O⁺ 1.4 × 10⁸, acetic acid 9.3 × 10⁴, phosphate monoanion 1.2 × 10⁴, MES 2.5 × 10², HFP 1.4, phosphate dianion 0.17, water 0.12/55.5 M⁻¹ s⁻¹; α 0.66 ± 0.07). **Phosphate catalysis saturates** at pH 6.12-6.62 (Fig. 7): x3.5-5 at 0.2-0.3 M | [M] / [F] |
| nucleophilic (amine) catalysis | morpholine (secondary amine) "more than three orders of magnitude more effective" than general acids of similar pKa, via the morpholine-formaldehyde iminium (Table IV); imidazole is NOT a nucleophilic catalyst and inhibits by binding formaldehyde (Table V) | [M] / [Q] |
| temperature / Ea | **none.** 25 °C throughout; morpholine series at 10 °C; the 10 °C k0 values in Table IV footnotes cannot be combined with the 25 °C profile into an Ea (sec. 6) | [!] |
| hemithioacetal (thiol trap of the aldehyde) | forms first, fast and reversibly: k 7.6 × 10⁴ (species 2) and 1.8 × 10⁵ M⁻¹ s⁻¹ (species 12) for thiolate attack (eq 2); kinetically competent but not rate-limiting | [M] |

**Answers.**

1. **`k_ttca_cys`: no bound** (confidence ~90 %). TC's K is ~10⁴ above Cao's TTCA K and the temperature is 25 °C; how
   that gap splits between k_f and k_r for an aldose is not in the paper. Arrhenius transfer of TC's k_r
   (5.4 × 10⁻⁵ min⁻¹, pH 5.6) to 100 °C gives 1.4 × 10⁻³ / 9.8 × 10⁻³ / 0.18 min⁻¹ for Ea 40 / 64 / 100 kJ/mol [Z],
   a band that brackets both the engine's 0.0286 and the dilute-pot ceiling ~1 × 10⁻³ min⁻¹. Uninformative as a number.
2. **E2 / E4 acid limb: CONTRADICTED in sign on the analogue** [Z, transfer]. In E2's notation
   k_r = k_r45 · 10^(−a (pH − 4.5)), TC gives a ≈ **−0.54** between pH 2.8 and 5.6 (k_r rises with pH), a ≈ 0 between
   5.6 and 6.6, a ≈ +0.5 between 6.6 and 7.5. E4's feasible region (a 0.4-0.6, `verdict.yaml`) has the opposite sign
   over the in-situ pH of Cao (4.2-4.7) and Zhai 2021 (4.7-5.5). Mechanism: the release step is C-N cleavage from the
   carbinolamine, whose concentration falls with acidity (iminium + H₂O ⇌ carbinolamine + H⁺). Fife 1991's acid limb is
   ring OPENING (C-S) or an N-acyl ring, not release from an N-H ring. Transfer assumption: TTCA's release is limited by
   the same C-N step at pH 4-5.5.
3. **E1 catalysed exchange: disfavoured at the size it needs** [Z, transfer]. (i) No carboxylate catalysis below pH 5
   (Cao's pot): the attack step limits and is not catalysed. (ii) Phosphate catalysis is real above pH 5.5 but
   SATURATES at x3.5-5 (half-saturation 17-29 mM total phosphate, Fig. 7 C₃); E1's linear k0(1 + c_P[HPO₄²⁻]) has no
   ceiling. Before saturation the measured slope is c_P ≈ 270 per M HPO₄²⁻ (pH 6.1-6.6) [Z], the top of E1's band.
   (iii) Amino-acid ammonium as a general acid (Brønsted line through H₃O⁺ and water, slope −0.62; pKa 8.9-9.6):
   c_aa ≈ 0.1-0.3 M⁻¹ at pH 6.0, 0.3-0.8 at 6.5, 0.9-2.6 at 7.0, 2.3-6.4 at 7.4, ±1 log for the scatter of points
   about the line [Z]. E1 was feasible only with c_aa ≥ 20-30 M⁻¹; E4's ~3 M⁻¹ is reachable only near pH 7.
   The one mechanism here that can be large is **nucleophilic catalysis by a secondary amine** (morpholine); whether a
   primary amino acid, or TTCA's own secondary ring nitrogen, does this for an aldose is untested anywhere on disk [P].
   That would be a concentration-dependent rate of exactly the kind CF-076 needs: an unlisted rival, not evidence for it.
4. **Microscopic reversibility check of the catalyst species**: Kallen's rate law eq 4, written on the predominant
   protonated cysteine, is v = k_cat′[RNH₃⁺][F][A⁻], so the buffer BASE form appears; for phosphate monoanion as acid,
   A⁻ is HPO₄²⁻. E1's use of [HPO₄²⁻] is consistent with this and with Fife 1991's general-base catalysis of iminium
   hydration (the reverse of the dehydration catalysed here) [Z].

## 2. Identity

| field | value |
|---|---|
| Title | "The Mechanism of Reactions Involving Schiff Base Intermediates. Thiazolidine Formation from L-Cysteine and Formaldehyde" |
| Author | Roland G. Kallen, Department of Biochemistry, School of Medicine, University of Pennsylvania |
| Venue | *J. Am. Chem. Soc.* **93**(23), 6236-6248, November 17, 1971; received August 3, 1970 |
| DOI | not printed. Crossref (2026-11-01): **10.1021/ja00752a040**, same title, author, JACS 93, 6236-6248 (1971). Matches the brief |
| Erratum | **Crossref lists "Additions and Corrections" to this paper: JACS 95, 6510 (1973), DOI 10.1021/ja00800a602. NOT on disk; its content is unknown** [!] |
| Companion | Kallen 1971, JACS 93, 6227 (equilibria; on disk, `kallen1971_extraction.md`) |
| Funding | NIH GM 13,777 |

## 3. Conditions

- **Reactants:** L-cysteine 2.5-10 mM (polarimetry), lower for UV; formaldehyde in large excess (0.01-0.40 M, Table I
  last column), pseudo-first-order; constants expressed against formaldehyde HYDRATE F [M].
- **Medium:** water, KCl to I = 1.0 M, EDTA present ("containing buffer, EDTA, and KCl"; concentration not printed
  here), pH set with KOH; buffers 0.01-0.1 M (HCl, formate, acetate, DABCO, phosphate) per Table I; Figs. 7 / 8 at
  pH held by 0.03 M MES [M]. pH measured (glass electrode); no drift reported; D₂O runs with pD = meter + 0.408.
- **Temperature:** 25 ± 0.1 °C; morpholine series (Table IV) at 10 °C [M].
- **Methods:** optical rotation (546.1 / 589.3 nm), UV absorbance below 240 nm (210 nm), titrimetry by alkali uptake
  at pH 5-6; at least three half-times; correlation ≥ 0.98. Zero-buffer rates by extrapolation (four or more buffer
  concentrations, duplicate) [M].
- **O₂ / headspace:** not stated; stock cysteine made daily because neutral stocks precipitate "probably cystine" [Q].

## 4. Quantitative results (from renders)

**Table I** (25 °C, I 1.0 M; k_obsd/[F] per TOTAL cysteine; K1(app) and k2 from saturation) [M]:

| pH | buffer, M (free base %) | k_obsd/[F], M⁻¹ s⁻¹ | K1(app), M⁻¹ | K1(app)/α_RNH2T, M⁻¹ | k2obsd, s⁻¹ |
|---|---|---:|---:|---:|---:|
| 1.00 | HCl | 0.0013 | | | |
| 1.60 | HCl | 0.0040 | | | |
| 2.84 | formate 0.01 (13) | 0.014 | | | |
| 3.35 | formate 0.01 (37) | 0.023 | | | |
| 4.47 | acetate 0.1 (40) | 0.12 | | | |
| 4.59 | acetate 0.1 (47) | 0.16 | | | |
| 4.60 | formate 0.02 (91) | 0.12 | | | |
| 4.77 | acetate 0.1 (57) | 0.19 | | | |
| 5.00 | acetate 0.1 (74) | 0.27 | | | |
| 8.32 | DABCO 0.1 (13) | 3.33 | | | |
| 8.43 | DABCO 0.1 (16) | | 4.4 | 26.7 | 0.67 |
| 8.97 | DABCO 0.1 (39) | | 6.8 | 25.7 | 0.35 |
| 9.41 | DABCO 0.1 (64) | 1.50 | | | |
| 9.47 | DABCO 0.1 (68) | | 9.6 | 26.9 | 0.21 |
| 9.95 | DABCO 0.1 (86) | | 14.8 | 31.0 | 0.12 (k2 0.12) |
| 10.50 | phosphate 0.1 (10) | | 15.5 / 17.8 | 22.2 / 25.3 | 0.12 / 0.10 |
| 11.08 | phosphate 0.1 (30) | | 19.3 | 21.8 | 0.12 |
| 11.57 | phosphate 0.1 (56) | | 20.0 | 20.9 | 0.17 |

Averages: K1 26.8 ± 4.0 M⁻¹; k2 0.126 ± 0.026 s⁻¹.

**Table II** (25 °C, I 1.0 M, water activity 1.0) [M]: K1 = [>NCH₂OH]/([>NH][F]) 26.8 ± 4.0 M⁻¹; k1a (COOH,
NH₂ species) 1.11 ± 0.13 × 10³ and k1b (COO⁻, NH₂) 2.92 ± 0.11 × 10³ M⁻¹ s⁻¹; k−1a 41.8 s⁻¹, k−1b 1.05 × 10² s⁻¹
(= k1/K1); k2 (water) 0.12 ± 0.012 / 0.126 ± 0.026 s⁻¹; k2′ (H₃O⁺) 1.40 ± 0.05 × 10⁸ M⁻¹ s⁻¹.

**Table III** (general-acid catalysis of N-hydroxymethylcysteine dehydration, 25 °C) [M]: catalyst, pH, concentration
range, n, fraction HA, pKa′, k2″ (M⁻¹ s⁻¹): H₃O⁺ 7.5-9.5, n 34, −1.7, 1.4 × 10⁸; acetic acid 5.72, 0.001-0.20 M, 15,
0.07, 4.60, 9.3 × 10⁴; MES 6.44-7.02, 0.01-0.30, 18, 0.20-0.39, 6.13, 2.5 × 10²; phosphate monoanion 6.12-6.62,
0.001-0.3, 47, 0.43-0.68, 6.50, 1.2 × 10⁴; HFP 9.26-9.74, 0.05-0.40, 22, 0.35-0.50, 9.22, 1.4; phosphate dianion
10.47, 0.05-0.30, 12, 0.90, 11.42, 0.17; HOH 9.5-13, 8, 15.7, 0.12/55.5. Imidazolium ~6 × 10³ (text).

**Fig. 7 (phosphate, pH 6.12 / 6.31 / 6.62, per amine free base, ×10³ M⁻¹ s⁻¹)** [F]: zero-buffer ~0.7 / ~0.65 /
~0.37; plateau at 0.2-0.3 M ~3.5 / ~2.25 / ~1.6 → **x5 / x3.5 / x4.3**. Printed fit coefficients C₁ 6.24, 6.22, 8.99;
C₂ 6.54, 8.55, 4.24 × 10⁻³; C₃ 2.9, 2.62, 1.72 × 10⁻² M; derived k1 2695, k2′ 1.87 × 10⁸, k2″ 1.22 × 10⁴ [M]. Text:
"two-threefold increases in rate from 0 to 0.025 M buffer" [M].

**Fig. 6 (log k per amine free base, no buffer)** [F]: ~3.05 at pH 1-3; ~3.4 at 4.5-5.2; 3.05 (5.7), 2.9 (6.1-6.3),
2.75 (6.4), 2.6 (6.6), 2.35 (7.0), 2.13 (7.4), 1.9 (7.7), 1.5 (7.9), 1.35 (8.4), ~1.0 (9.0), ~0.5 above pH 10.
**Fig. 3 (log k per total cysteine)** [F]: −2.9 (pH 1) rising to ~0.6 at pH 7-7.5, dip ~0.4 at 8.5-9, ~0.5 at 10.

**Table IV** (morpholine-catalysed, 10 °C, F_T 0.025 M) [M]: k*/a_H × 10⁹ (M⁻² s⁻¹) average 4.3 ± 0.7 (pH 6.38),
5.4 ± 0.6 (7.22), 5.9 ± 1.6 (8.35), 16.1 ± 5 (9.25). **Table V** (imidazole, pH 7.79, 25 °C): k_obsd 2.0 × 10⁻² s⁻¹
with no imidazole falling to 0.52 × 10⁻² at 0.30 M; k_obsd/[F_free] constant 0.60-0.73 M⁻¹ s⁻¹ [M].

**Derived [Z]** (microscopic pKa from Kallen 1971 Table II; Fig. 6 readings; Kallen 1971 Table IV K_app):

| pH | α free-base amine | k_f, M⁻¹ s⁻¹ (total Cys) | K_app, M⁻¹ | k_r, s⁻¹ | k_r, min⁻¹ |
|---|---:|---:|---:|---:|---:|
| 2.84 | — | 0.014 [M] | ~4.9 × 10⁵ | 3 × 10⁻⁸ | 1.8 × 10⁻⁶ |
| 4.47 | — | 0.12 [M] | ~5.3 × 10⁵ | 2.3 × 10⁻⁷ | 1.4 × 10⁻⁵ |
| 5.57 | 4.7 × 10⁻⁴ | 0.59 | 6.53 × 10⁵ | 9.0 × 10⁻⁷ | 5.4 × 10⁻⁵ |
| 6.56 | 4.5 × 10⁻³ | 2.2 | 2.32 × 10⁶ | 9.3 × 10⁻⁷ | 5.6 × 10⁻⁵ |
| 7.45 | 3.2 × 10⁻² | 4.2 | 1.30 × 10⁷ | 3.2 × 10⁻⁷ | 1.9 × 10⁻⁵ |

Check against Kallen 1971 Fig. 4 (TC 7.5 mM hydrolysis, 0.01 M phosphate): the relaxation k_r + k_f·2[Cys]eq predicts
t½ ≈ 71 min at pH 5.39 (read 1.2-1.4 h) and ≈ 29 min at pH 6.19 with x1.7 phosphate catalysis (read ~0.4 h). Both
agree; the relaxation is ~130x (pH 5.4) to ~210x (pH 6.2) faster than k_r because re-formation dominates it [Z].

## 5. Mechanism evidence and its kind

- Kallen laboratory (Penn; Jencks school). Kind 6 (one factor varied inside one study: pH 1-13, buffer identity and
  concentration, formaldehyde concentration, D₂O) plus kind 4 (the hemithioacetal intermediate seen directly as a fast
  first UV phase with isosbestic points, Fig. 2; the carbinolamine by saturation kinetics, Fig. 5) [M].
- The change in rate-determining step is shown three ways: the bell in the pH profile, saturation in formaldehyde above
  pH 8, and saturation in phosphate at pH 6.1-6.6 [M/Q]. The sulfonium-ion path is excluded by a 10¹⁰ rate deficit
  (footnote 49) [Q].
- Ring closure fast: shown indirectly (morpholine-catalysed path goes faster through the same iminium) [Q].
- No temperature series, no TC hydrolysis kinetics measured directly here, no sugar.

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| CF-076 E1 (catalysed exchange) | **CONTRADICTS the linear form; weakly against the size** | Catalysis saturates (x3.5-5, phosphate, pH 6.1-6.6); none by carboxylates below pH 5; ammonium general-acid c_aa ≤ ~3 M⁻¹ (pH ≤ 7, central). Transfer from formaldehyde assumed |
| CF-076 E2 / E4 (in-situ pH, acid limb) | **CONTRADICTS the sign of `a` below pH 5.6** | TC k_r rises with pH (a ≈ −0.54) from 2.8 to 5.6; E4's feasible a 0.4-0.6 |
| CF-076 E3 (second Cys sink) | LEAVES UNCHANGED | hemithioacetals are fast and reversible (not a sink) |
| CF-076, unlisted rival | **NEW lead [P]** | amine nucleophilic catalysis (morpholine x10³ over Brønsted) is concentration-dependent; untested for primary amines or TTCA's ring N with an aldose |
| CF-011 (phosphate catalysis, Xie vs Jiangnan) | LEAVES UNCHANGED for the Amadori step; **against phosphate catalysis of the ring exchange at Cao's in-situ pH 4.2-4.7** | attack-limited, no buffer catalysis below pH 5 |
| evidence map 2.12 "cysteine fate / TTCA: catalysed hydrolysis vs K(pH)" | sharpens | K(pH) and k_r(pH) both pH-shaped on the analogue; catalysis bounded |
| triage item 3 (thiol losses, cysteine fate) | supports a structured exchange law (two steps, saturable catalysis) over a lumped k |
| engine `k_ttca_cys` | LEAVES UNCHANGED | no transferable constant |

## 7. Disagreements found

1. **With `kallen1971_extraction.md` arithmetic 3 and its sec. 5 row "relaxation k ≈ 1.5-5 × 10⁻⁴ s⁻¹; k_f ~10³
   M⁻¹ s⁻¹"** [Z]. That dossier treated the approach-to-equilibrium relaxation as the hydrolysis constant ("with
   [Cys]eq/[TC] ~ 1 % ... ≈ the hydrolysis constant"). For A + B ⇌ C near equilibrium the relaxation is
   k_r + k_f([A]eq + [B]eq) = k_r(1 + 2[TC]/[Cys]eq) ≈ 130-210 k_r here. Measured k_f + K give k_r ≈ 1 × 10⁻⁶ s⁻¹
   (~120-270x lower than 1.5-5 × 10⁻⁴) and k_f ≈ 0.6-2.7 M⁻¹ s⁻¹ per total cysteine (the "10³" is the free-base-amine constant).
   Cause: arithmetic, not data.
2. **With `docs/history/cys_exchange_2026-10-31/prereg.yaml` E1 k0 band basis** ("fast end Kallen 1971 TC relaxation
   0.009-0.03 /min at 25 C x Arrhenius"): inherits (1). TC's k_r at 25 °C is 5.4 × 10⁻⁵ min⁻¹; the fast end so derived
   would be ≤ 0.18 min⁻¹ at 100 °C (Ea 100), not 1.0. The band's slow end and the verdicts are not affected by this
   alone (E1 failed on c_aa, not k0).
3. **With the E2 / E4 `a` band [0, 1] (basis Fife 1991)**: sign, as in sec. 1 (2). Candidate cause: Fife measured ring
   opening and N-acyl / aryl substrates; Kallen the N-H ring's overall formation (and by detailed balance release).
4. **With Fife 1991 on linearity of buffer catalysis**: Fife's general-base catalysis of iminium hydrolysis is linear to
   0.5 M (aryl iminiums, slow hydration); Kallen's saturates by 0.1 M (attack takes over). Chemistry (which step limits),
   not method.
5. Internal: k2′ 1.40 × 10⁸ (Table II) vs 1.87 × 10⁸ (Fig. 7 fit); Table IV footnote k0 at 10 °C (pH 6.38 / 7.22 /
   8.35 / 9.25: 1.32 × 10⁻³ / 8.83 × 10⁻³ / 8.53 × 10⁻² / 0.228 s⁻¹) against 25 °C gives 25/10 °C ratios from ~12 to
   ~0.4, physically incoherent as an Ea (pK shifts with T, formaldehyde basis unclear): no Ea is taken [!].
6. The 1973 "Additions and Corrections" (not on disk) may revise any of the above [!].

## 8. Ledger-ready measurements

1. **Table I** k_obsd/[F] vs pH (pH 1.00-9.41, total cysteine, 25 °C, I 1.0 M): the forward pH law of a cysteine
   thiazolidine, printed, two readers trivial.
2. **Table II** k1a, k1b, K1, k2, k2′: the microscopic constants.
3. **Table III** general-acid k2″ with pKa′: the Brønsted line (catalyst class anchors for E1-type bands).
4. **Fig. 7** zero-buffer and plateau per pH (two readers): the saturation ceiling of phosphate catalysis.
5. **"No buffer catalysis ... 0.05-0.4 M ... pH 2.3-5.0"** (carboxylate buffers): a null result worth a row.
6. Derived k_r(pH) (sec. 4) only as a [Z] note beside Kallen 1971's K_app rows, never as a measurement.
