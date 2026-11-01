# Nikolantonaki & Waterhouse 2012 — EXTRACTION (4-methyl-1,2-benzoquinone (Q4MeC, 1 mM) + 4 mM of each of ten nucleophiles in model wine (12 % ethanol, 5 g/L tartaric acid, pH 3.5), 10 °C: pseudo-first-order quinone decay constants K by UV at 400 nm, ten replicates; products by HPLC-HRMS. Thiols: 2-furanmethanethiol (2FMT = FFT), 3-sulfanylhexanol, 4-methyl-4-sulfanylpentan-2-one, GSH, H₂S; also SO₂, ascorbate, phloroglucinol, Met, Phe. NO cysteine, NO MFT)

### THE ONLY MEASURED RATE ON DISK FOR A MAILLARD ODORANT THIOL (FFT) WITH AN ELECTROPHILE. K(2FMT) = 0.0837 s⁻¹ at 4 mM, i.e. k₂ ≈ 21-24 M⁻¹ s⁻¹ [Z], 4.5x slower than GSH and H₂S, 14x faster than the tertiary thiol 4MSP. One electrophile (a quinone, absent from sugar + cysteine pots), one pH (3.5), one temperature (10 °C), 12 % ethanol; thiol-quinone adducts re-aromatise and are taken as irreversible. No Ea, no pH law, no reverse rate.

**Read 2026-11-01 from `data/articles/nikolantonaki2012.pdf`** (8 pp., the ASAP version paginated A-H with volume
"XXXX"; the issue version is JAFC 60(34), 8484-8491). Born-digital, clean text layer. **Table 1 and Fig. 2 checked
against a 130 dpi render** (scratchpad `dos6/nik-4.png`). Fig. 2 read ±0.1 ln unit.

**Legend.** **[V]** verbatim · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative claim ·
**[C]** cited · **[P]** proposed · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| model need | what Nikolantonaki 2012 gives | tag |
|---|---|---|
| thiol + electrophile forward rate | K (s⁻¹, 4 mM Nu, 1 mM Q4MeC): 4MSP 0.0060, 3SH 0.0578, **2FMT 0.0837**, SO₂ 0.3343, AA 0.3471, GSH 0.3808, H₂S 0.4188; Phl 0.0064; Met 0.0002, Phe 0.0005 (Table 1, means of 10) | [M] |
| second-order k₂ | not printed; K/[Nu]: 2FMT 21 (K/4 mM) to 24 (K/3.5 mM, mean Nu over the run) M⁻¹ s⁻¹; 3SH 14-17; 4MSP 1.5-1.7; GSH 95-109; H₂S 105-120; SO₂ 84-96 | [Z] |
| reversibility | not measured; products are thiol-catechol adducts (mono- and di-adducts, catechol dimers), HRMS; SO₂ mostly REDUCES the quinone back (89 % 4MeC, 11 % sulfonates); ascorbate reduces it completely | [M] / [Q] |
| pH law, buffer, Ea | **none** (pH 3.5, 10 °C only) | — |
| amino acids | Met, Phe "essentially zero" at pH 3.5; no Strecker aldehyde (contrast with Rizzi's pH 7.17 phosphate result) [Q/C] | [M] |

**Answers.**

1. **MFT loss law: no bound on λ** (confidence ~95 %). Transfer needs (i) an electrophile of quinone reactivity in a
   sugar + cysteine pot (none shown on disk), (ii) MFT ≈ FFT in nucleophilicity (MFT is a heteroaryl thiol, FFT a
   furfuryl (benzylic-type) thiol: different pKa class, not in the paper), (iii) an Ea from 10 to 100-145 °C (none),
   (iv) a pH law from 3.5 to 5-6.5 (none). λ = k₂[E]: at 10 °C, the engine's 1.6 × 10⁻³ s⁻¹ would need ~70-80 µM of a
   Q4MeC-like quinone [Z].
2. **Structure of the loss law, in kind** [M → Z]: thiol loss to an electrophile is second order (A2), and a sacrificial
   nucleophile present in excess protects the odorant thiol: GSH, SO₂ and H₂S react 4.5-5x faster than 2FMT per mole and
   ascorbate reduces the quinone [Q, the authors' conclusion]. In a pot with 30-100 mM cysteine (and H₂S), a µM odorant
   thiol takes a negligible share of any electrophile flux: the A3-type "suppressed by free cysteine" via competition
   (cysteine not measured here; GSH stands in).
3. **Steric order among odorant thiols** [M]: primary (2FMT) > secondary (3SH, x0.69) > tertiary (4MSP, x0.07). MFT's
   thiol is on an aromatic ring carbon, a fourth class not covered.
4. **CF-076, `k_ttca_cys`**: nothing.

## 2. Identity

| field | value |
|---|---|
| Title | "A Method To Quantify Quinone Reaction Rates with Wine Relevant Nucleophiles: A Key to the Understanding of Oxidative Loss of Varietal Thiols" |
| Authors | Maria Nikolantonaki, Andrew L. Waterhouse* — Department of Viticulture and Enology, University of California, Davis |
| Venue | *J. Agric. Food Chem.* 2012, **60**(34), 8484-8491 (Crossref); PDF on disk is the ASAP version (pp. A-H); received May 8, revised July 18, accepted August 6, 2012 |
| DOI | **10.1021/jf302017j** (printed on every page; Crossref 2026-11-01: same title, authors, JAFC 60(34) 8484-8491). Matches the brief |
| Funding | American Vineyard Foundation |
| Related | Nikolantonaki et al. 2010 (ref 11, thiols + catechin / epicatechin under wine oxidation; not on disk) |

## 3. Conditions

- Model wine: 12 % v/v ethanol, 5 g/L L-tartaric acid, pH 3.50; < 1 % acetonitrile from the quinone stock; thiol
  stocks in ethanol (< 2 % extra) [M]. Q4MeC 1 mM (from 23 mM stock made with periodate resin under argon, ~95 %
  conversion, used within 30 min); nucleophile 4 mM, added before the quinone, thermostated cell, **10 °C**; dead time
  0.5 s; ≥ 4 half-lives; 10 replicates [M]. Thiol stocks standardised with DTNB.
- Without nucleophile Q4MeC decays zero order (not significant on the timescale) [Q].
- Products: Polymerx RP-1 HPLC, LTQ-Orbitrap negative ESI, < 2 ppm.
- O₂: quinone made under argon; reaction headspace not stated.

## 4. Quantitative results

**Table 1** [M] (regression ln(Q/Q₀) = intercept − K·t; relative K vs 3SH):

| Nu | regression | K, s⁻¹ | relative |
|---|---|---:|---:|
| Met | −0.0055 − 0.0002x | 0.0002 | 0.003 |
| Phe | −0.0082 − 0.0005x | 0.0005 | 0.009 |
| Phl | −0.1623 − 0.0064x | 0.0064 | 0.110 |
| 4MSP | −0.0258 − 0.0060x | 0.0060 | 0.100 |
| 3SH | −0.2295 − 0.0578x | 0.0578 | 1.000 |
| 2FMT | −0.3409 − 0.0837x | 0.0837 | 1.440 |
| SO₂ | −0.3549 − 0.3343x | 0.3343 | 5.800 |
| AA | −0.5987 − 0.3471x | 0.3471 | 6.000 |
| GSH | −0.6169 − 0.3808x | 0.3808 | 6.600 |
| H₂S | −0.7535 − 0.4188x | 0.4188 | 7.200 |

**Fig. 2** [F] (ln K, mean ± SD): 4MSP ≈ −5.2 ± 0.2, 3SH ≈ −2.75 ± 0.2, 2FMT ≈ −2.3 ± 0.15, SO₂ / AA / GSH / H₂S ≈
−1.1 to −0.9 ± 0.15; Met ≈ −9.7, Phe ≈ −8.9, Phl ≈ −7.6 (wide SD). Letters a-e mark significance groups.

**Table 2** [M]: 2FMT products C₁₂H₁₂O₃S (mono-adduct, three isomers 2:6:1), C₁₇H₁₆O₄S₂ (di-adduct), C₁₉H₁₈O₅S
(dimeric catechol + thiol, low yield); analogous sets for 3SH, 4MSP, H₂S (incl. S-bridged dimers), GSH (mono 16:2:1,
di); Met and Phe adducts not detected.

## 5. Mechanism evidence and its kind

- UC Davis (Waterhouse). Kind 6 within one study (nucleophile identity at fixed medium) and product identification by
  HRMS (kind 4-like, end-point only; no time-resolved intermediates, no labelling) [M].
- Proposed [P]: Michael-type 1,4-addition at three ring positions, then oxidation of the thio-catechol by more quinone
  and a second addition; catechol-dimer quinones via redox exchange (Scheme 2).

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| MFT loss law (A0 first-order, A2 electrophile pool, A3 cysteine protection) | SUPPORTS A2's form and A3's competition in kind; no bound on λ | analogue electrophile, 10 °C, pH 3.5 |
| evidence map 2.12 "MFT loss: O₂ / carbonyl dependence" | SUPPORTS that an oxidation-born electrophile (quinone) consumes thiols fast; O₂ dependence not tested | |
| CF-031 (lipid vs thiols) | LEAVES UNCHANGED | no lipid |
| CF-076, `k_ttca_cys` | LEAVES UNCHANGED | |

## 7. Disagreements found

1. **Abstract vs Table 1** [!]: the abstract (and the Results paragraph) assign K = 0.3343, 0.3471, 0.3808 to "SO₂,
   GSH, AA"; Table 1 gives AA 0.3471 and GSH 0.3808. This dossier follows Table 1.
2. **Fig. 2 vs Table 1 for phloroglucinol** [!]: Fig. 2 plots ln K ≈ −7.6, Table 1's K 0.0064 is ln −5.05. The thiol
   points agree with Table 1 within ~0.2; Fig. 2 is probably a mean of per-replicate ln K, which a noisy low rate drags
   down. Irrelevant to the thiol numbers.
3. **The regressions' intercepts imply a fast initial phase** [Z][!]: ln(Q/Q₀) at t = 0 is −0.34 (2FMT), −0.62 (GSH),
   −0.75 (H₂S), i.e. 29-53 % of the quinone gone before the fitted line starts, while K × 0.5 s (dead time) accounts for
   only −0.04 to −0.21. Either Q₀ is nominal or the true initial rate exceeds K. If the latter, the printed K
   UNDERSTATES the thiol-quinone rate (by an unknown factor) and the [Z] k₂ above are lower bounds.
4. No overlap with other dossiers on the same quantity.

## 8. Ledger-ready measurements

1. **Table 1, 2FMT row** (K 0.0837 s⁻¹, 4 mM, 1 mM Q4MeC, 10 °C, pH 3.5, model wine, n = 10): the one odorant-thiol
   electrophile constant on disk; key with the intercept flag (sec. 7.3).
2. Table 1 thiol rows (3SH, 4MSP, GSH, H₂S) as the steric / class series around it.
