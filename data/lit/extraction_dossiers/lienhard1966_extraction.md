# Lienhard & Jencks 1966 — EXTRACTION (equilibria and kinetics of thiol addition to aliphatic aldehydes and ketones → hemithioacetals / hemithioketals in water, 25 °C, ionic strength 1 M (KCl), 10⁻⁴ M EDTA, pH 0-5; UV carbonyl absorbance at 278-286 nm; acetaldehyde, isobutyraldehyde, pivaldehyde, pyruvate, acetone × 2-methoxyethanethiol, ethanethiol, methyl mercaptoacetate, 2-mercaptoethanol, benzenethiol)

### HEMITHIOACETAL FORMATION IS A FAST, MODERATELY FAVOURABLE PRE-EQUILIBRIUM: **K_hemi (vs unhydrated carbonyl) = 32-60 M⁻¹ for acetaldehyde, 15-18 M⁻¹ for isobutyraldehyde, 3.8-4.8 M⁻¹ for pivaldehyde, 5.2 M⁻¹ for pyruvate, < 0.2 M⁻¹ for acetone** (25 °C, 1 M KCl), and **nearly independent of thiol pKa** (2× over pKa 7.8-10.3). Rate law v = k_H⁺·a_H⁺·[RSH][C=O] + k_RS⁻·[RS⁻][C=O]; **k_RS⁻ ≈ 1.4-2.8 × 10⁷ M⁻¹ min⁻¹ for acetaldehyde regardless of thiol pKa**; no uncatalysed thiol pathway. Unhydrated carbonyl is the reactive species. **No aromatic or furanic aldehyde, no cysteine, no H2S, no temperature dependence.**

**Source on disk:** `data/articles/lienhard1966.pdf` (14 pp., ACS scan with OCR text layer; the PDF also contains the end of
the preceding article and the start of the following one, about a Grignard addition — ignore those pages). Tables I and II
re-read from a 200-dpi render (`scratch/lit8/lien-04.png`, `lien-05.png`); Fig. 2 not digitised (its lines are eq 7 with
Table II constants).

## 0. Identity

| field | value |
|---|---|
| Title | "Thiol Addition to the Carbonyl Group. Equilibria and Kinetics" |
| Authors | Gustav E. Lienhard, **William P. Jencks** — Brandeis University / Harvard |
| Venue | J. Am. Chem. Soc. **88** (17), 3982-3995 (article starts p. 3982), 5 Sep 1966. Received 4 May 1966 |
| DOI | not printed (1966); the brief's **10.1021/ja00969a017** matches title, journal, volume, issue, first page. **Identity verified.** |

## 1. Why it matters

Brief item 8: equilibrium and rate constants for thiol + aldehyde hemithioacetals, pKa dependence, T, pH; any data for
aromatic/furanic aldehydes. Relevant to furfural + H2S/HS⁻ → FFT (an addition to the carbonyl), to thiol sequestration of
aldehydes (Strecker aldehydes, methylglyoxal, mercaptoacetaldehyde) and to whether such adducts need explicit species.

## 2. Methods (pp. 3983-3984)

- 25 °C; ionic strength 1 M with KCl; 10⁻⁴ M EDTA; thiols distilled under N2; aldehydes freshly distilled.
- **Equilibria:** decrease of carbonyl n→π* absorbance (278-286 nm) at fixed carbonyl (0.02-0.05 M) and 3+ thiol
  concentrations (0.014-2.0 M), mostly at pH 4.6 (≤ 0.1 M acetate) or in HCl; benzenethiol by its own 236 nm band.
  K_obsd = [adduct]/([thiol][hydrated + unhydrated carbonyl]); K_hemi = K_obsd (1 + K_hyd).
- **Kinetics:** pseudo-first-order in excess thiol; first reading ~5 s after mixing; ±5 %; pH 0-4.3 (HCl, formate,
  acetate); acetaldehyde K_hyd = 0.85 (1 N KCl, 25 °C; measured here), 0.99 in D2O.

## 3. Numbers

### Table I (p. 3985), equilibria, 25 °C, μ = 1 M

| carbonyl (M) | thiol (M) | pH | n | K_obsd (M⁻¹) | K_hyd | K_hemi (M⁻¹) |
|---|---|---|---|---|---|---|
| acetaldehyde (0.030) | 2-methoxyethanethiol (0.030-0.25) | 4.6 | 5 | 17.4 ± 1 | 0.85 | **32** |
| acetaldehyde (0.030) | same, in D2O | 4.6 | 5 | 36.4 ± 4 | 0.99 | 72 |
| acetaldehyde (0.032) | methyl mercaptoacetate (0.018-0.27) | 4.6 | 6 | 32.7 ± 3 | 0.85 | **60** |
| acetaldehyde (0.034) | ethanethiol (0.014-0.058) | 4.6 | 3 | 19.6 ± 3 | 0.85 | **36** |
| isobutyraldehyde (0.025) | 2-methoxyethanethiol (0.083-0.23) | 2.0 | 5 | 12.8 ± 1 | 0.24 | 16 |
| isobutyraldehyde (0.025) | same (0.02-0.33) | 4.3 | 6 | 11.8 ± 1 | 0.24 | 15 |
| isobutyraldehyde (0.02-0.20) | benzenethiol (10⁻⁴) | 1 N HCl | 7 | 14.8 | 0.24 | 18 |
| isobutyraldehyde (0.01-0.38) | benzenethiol (10⁻⁴) | 4.6 | 10 | 14.8 | 0.24 | 18 |
| pivaldehyde (0.025) | 2-methoxyethanethiol (0.06-0.25) | 4.6 | 3 | 4.8 ± 0.7 | < 0.1 | 4.8 |
| pivaldehyde (0.025) | 2-mercaptoethanol (0.05-0.97) | 0.4 N HCl | 4 | 4.2 ± 0.5 | < 0.1 | 4.2 |
| pivaldehyde (0.025) | 2-mercaptoethanol (0.026-2.0) | 4.6 | 8 | 3.8 ± 0.5 | < 0.1 | 3.8 |
| pyruvate (0.02) | methyl mercaptoacetate (0.05-0.28) | 5.2 | 5 | 5.2 ± 0.2 | not hydrated | 5.2 |
| acetone (0.05) | 2-methoxyethanethiol (0.52) | 1.5 | 1 | < 0.2 | not hydrated | < 0.2 |

Text (p. 3989): K for thiol + acetone is ≥ 100× smaller than for acetaldehyde, which is ~20× smaller than for formaldehyde
(hydrate basis, Kallen & Jencks). D2O raises K_hemi 2.27× (K_H2O/K_D2O = 0.44).

### Table II (p. 3986), rate constants, 25 °C, μ = 1 M

| aldehyde | thiol | pK'a (RSH, 25 °C, μ 1 M) | k_H⁺ (M⁻² min⁻¹) | k_RS⁻ (M⁻¹ min⁻¹) |
|---|---|---|---|---|
| acetaldehyde | methyl mercaptoacetate | 7.83 | 84 | **1.41 × 10⁷** |
| acetaldehyde | 2-methoxyethanethiol | 9.50 | 530 | **2.82 × 10⁷** |
| acetaldehyde | ethanethiol | 10.25 | 1 960 | **2.83 × 10⁷** |
| pivaldehyde | 2-mercaptoethanol | 9.51 | 200 / 228 | 5.8 / 6.4 × 10⁶ |
| pivaldehyde | 2-methoxyethanethiol | 9.50 | 180 | 5.6 × 10⁶ |

(pivaldehyde: from disappearance / appearance.) Rate law (eqs 6-7): k_obsd = (k_H⁺ a_H⁺ + k_RS⁻ K'_RSH/(K'_RSH + a_H⁺))
× (K_hemi[RSH] + 1 + K_hyd)/(K_hemi + K_hemi K_hyd). Minimum of the pH-rate profile at pH ~1.5-3 (Fig. 2).
General acid catalysis (acetaldehyde + 2-methoxyethanethiol): dichloroacetic 28 M⁻² min⁻¹ (pK'a 1.3), difluoroacetic 25
(1.15), phosphoric 23 (1.7), cyanoacetic 2.4 (2.4); Brønsted α ≈ 0.7 (0.4 through H3O⁺ and strong carboxylic acids);
no catalysis by formate/acetate ions. Isotope effect on k_H⁺: k_H2O/k_D2O = 0.59. Pivaldehyde: k_dehyd 0.75, k_hyd
0.57 min⁻¹ — dehydration becomes rate-limiting above the measured pH range; also rate-limiting for formaldehyde and
methylglyoxal hydrates (cited).

### Arithmetic (mine)

1. **Bound fraction.** At 25 °C, 33 mM of an alkanethiol binds acetaldehyde K_obsd[RSH] = 17-33 × 0.033 = 0.6-1.1 → 36-52 %
   of total acetaldehyde as hemithioacetal; a branched Strecker-type aldehyde (isobutyraldehyde) 28 %; pivaldehyde 13 %.
2. **Speed.** At pH 7 with a thiol of pKa 8.5 (thiolate fraction ~3 %), k_RS⁻ 2.8e7 → pseudo-first-order for 33 mM thiol
   ≈ 2.8e7 × 0.033 × 0.03 ≈ 3e4 min⁻¹ (t½ ~ 1 ms); dehydration of the aldehyde hydrate then limits (~1 min⁻¹ at 25 °C).
   **Hemithioacetal formation is a pre-equilibrium on every engine time scale.**
3. **pKa insensitivity of K:** the equilibrium shifts with R by almost exactly as much as the thiol's acidity (p. 3989), so
   a single K per carbonyl suffices across alkyl thiols (cysteine was not tested).

## 4. Numbers the repository can use

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| RSH + acetaldehyde ⇌ hemithioacetal | K_hemi | 32-60 | M⁻¹ | 25 °C, μ 1 M, pH 4.6 | Table I | measured |
| RSH + isobutyraldehyde | K_hemi | 15-18 | M⁻¹ | same, pH 1-4.6 | Table I | measured |
| RSH + pivaldehyde | K_hemi | 3.8-4.8 | M⁻¹ | same | Table I | measured |
| RSH + pyruvate | K_hemi | 5.2 | M⁻¹ | pH 5.2 | Table I | measured |
| RSH + acetone | K | < 0.2 | M⁻¹ | pH 1.5 | Table I | upper bound |
| RS⁻ + acetaldehyde | k_RS⁻ | 1.4-2.8e7 | M⁻¹ min⁻¹ | 25 °C | Table II | measured (extrapolated from pH ≤ 4.3) |
| H⁺-catalysed path | k_H⁺ | 84-1 960 | M⁻² min⁻¹ | 25 °C | Table II | measured |
| T dependence, aromatic/furanic aldehydes, cysteine, H2S | — | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** order-of-magnitude K for alkanethiol adducts of aliphatic aldehydes (tens of M⁻¹ unbranched, a few M⁻¹ for
α-branched) and the fact that the adducts equilibrate in milliseconds-to-minutes; a pH-independent K.
**Cannot:** furfural (conjugated aromatic aldehyde, no data here), H2S/HS⁻ as the nucleophile, cysteine (which goes on to
thiazolidines), any temperature dependence (one temperature; no ΔH), conditions above pH 5 for equilibria.

## 6. Flags

1. The PDF bundles fragments of the adjacent JACS articles (p. 3982 top; final pages on a benzophenone Grignard reaction).
2. k_RS⁻ is inferred from the basic limb at pH ≤ 4.3 assuming thiolate attack; no measurement near neutral pH.
3. μ = 1 M KCl throughout; activity effects relative to Maillard buffers unassessed.

## 7. Use in the model

- **No furanic data**: nothing here can be used directly for furfural + H2S → FFT or for furfural sequestration by
  cysteine. By analogy, an aliphatic-aldehyde hemithioacetal with K ~ 10-60 M⁻¹ at 25 °C binds a substantial fraction of an
  aldehyde at tens of mM thiol but < 1 % at sub-mM thiol; at 100-145 °C the adducts (exothermic, ΔH unknown) will be weaker.
  An engine treatment as a fast pre-equilibrium (not a kinetic step) is supported; a declared K with a declared ΔH would be
  an assumption, not an anchor.
- **Mercaptoacetaldehyde (MAA, Kobayasi 1965)** is itself a thiol + aldehyde; this paper's K's suggest its intermolecular
  hemithioacetal/dithiane equilibria are significant at tens of mM, which reduces free MAA — relevant if B59 adds MAA →
  H2S + acetaldehyde. Not quantified here.
- **Methylglyoxal:** the GSH + methylglyoxal hemithioacetal (K_d 3 mM, Szajewski 1980 eq 5) is the closer analogue for the
  `r_cys_actz` / MGO pool than these simple aldehydes.
- No bearing on the furfural barrier swap or B61.
