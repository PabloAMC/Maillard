# De Haan, Corrigan, Smith, Stroik, Turley, Lee, Tolbert, Jimenez, Cordova & Ferrell 2009 — EXTRACTION (glyoxal + glycine, serine, aspartic acid, ornithine, arginine: products by NMR / ESI-MS / AMS; 1,3-disubstituted imidazoles whose "odd" ring carbon comes from a SECOND glyoxal, with formic acid released; bulk NMR kinetics at 1 M, 295 K, D2O; drying-droplet and aerosol-uptake experiments)

### Mechanism paper for the glyoxal + amino acid → imidazole + formic acid sink (product identity, stoichiometry, no Strecker-formaldehyde needed). Its rate constants are SUPERSEDED by the same lab (sedehi2013: x250 for glycine, x400 for serine) because the pH was assumed, not monitored.

**Read 2026-10-31 from `data/articles/haan2009.pdf`** (7 pp., Environ. Sci. Technol. 43:2818-2824; born-digital, clean
text layer). Table 1, Schemes 1-3 and eqs 1-3 read from the text layer; the ESI-MS, FTIR, kinetics and AMS figures
(Figs. 1-4) were not needed for numbers (all kinetic numbers are in Table 1 and the text) and were not re-rendered.
**SI (Tables S1, Figs. S1-S4: NMR assignments, 1/[AA] plots, product time courses, AMS spectra) NOT on disk.**

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| what the model needs | what De Haan 2009 gives | tag |
|---|---|---|
| what glyoxal + amino acid makes in water | **1,3-disubstituted imidazoles** (zwitterionic, "extremely stable") with Gly, Ser, Asp, Orn; with Gly and Ser "the only products detected" by ESI-MS / 13C / 2-D NMR [Q], plus browning (melanoidin-like) and a minor diimine; arginine instead forms 1:1 side-chain adducts (two non-aromatic five-membered rings, 3:1) | [M] / [Q] |
| where the odd ring carbon comes from | **a second glyoxal molecule**, with C-C rupture and loss of formic acid (Scheme 2) — not formaldehyde from a Strecker degradation: products carrying a third R group (m/z 275 Ser, 359 Asp, 356 Orn) were not observed [M → Q] | [M] / [P] |
| stoichiometry | glycine: imidazole formed at 60 ± 20 % of the glycine loss rate (signal-matched 2 H each) → ~1.2 ± 0.4 imidazole per 2 Gly; diimine 15 ± 5 % → 0.3 ± 0.1 per 2 Gly; total 1.5 ± 0.5 per 2 Gly, "within uncertainty" of 1:2 [M] | [M] |
| order | 1/[reactant] linear in time (equimolar) → second order; diimine and imidazole yields proportional throughout (diimine a side-product) [M/Q]; drying-droplet AMS gives overall order ~1.6 [M] | [M] |
| rate constant | Table 1 (295 K, 1 M glyoxal): Gly 0.12 ± 0.04, Ser 0.09 ± 0.01, Arg (4.5 ± 0.6) × 10⁻⁴ M⁻¹ s⁻¹ on REACTIVE forms (f_glyoxal 7 × 10⁻⁴ computed; f_AA at an ASSUMED pH 6) — **superseded** (sec. 6) | [M] / [!] |
| temperature | "The temperature dependence of glyoxal-amino acid reactions is unknown" [V] | [!] |
| amine site selectivity | ornithine: 80 % at the α-amine (pKa 8.69) vs 20 % side chain (10.76), by 1H NMR integration [M] | [M] |

**Answer.** For a food pot this paper settles the PRODUCT and STOICHIOMETRY of the glyoxal + primary-amino-acid sink
at 25 °C: 2 glyoxal + 2 amino acid → 1,3-di(carboxyalkyl)imidazole + formic acid (+ water), irreversible, with browning
as a side channel. That is a glyoxal sink that consumes amine without making a Strecker aldehyde or an aminoketone,
and a formic-acid source (relevant to the trunk's acid budget, CF-069) [Z]. It gives no usable rate: use sedehi2013's
pH-monitored constants (same lab), and only after the transfer listed there.

## 0. Identity

| field | value |
|---|---|
| Title | "Secondary Organic Aerosol-Forming Reactions of Glyoxal with Amino Acids" |
| Authors | David O. De Haan* (University of San Diego; CIRES Boulder), Ashley L. Corrigan, Kyle W. Smith, Daniel R. Stroik, Jacob J. Turley, Frances E. Lee (San Diego), Margaret A. Tolbert, Jose L. Jimenez, Kyle E. Cordova, Grant R. Ferrell (Boulder) |
| Venue | *Environ. Sci. Technol.* **2009**, 43 (8), 2818-2824; received 12 Dec 2008, revised 17 Feb 2009, accepted 19 Feb 2009, web 12 Mar 2009 |
| DOI | **10.1021/es803534f** (printed; Crossref 2026-10-31: same title, ten authors, vol. 43, 2818-2824; Crossref files the first author under family name "Haan", hence the PDF stem `haan2009`). Matches the brief |
| Lab | De Haan group (San Diego) with Tolbert / Jimenez (Boulder); same lab as sedehi2013 |

## 2. Conditions

| item | as printed |
|---|---|
| reagents | 15N-enriched glycine (Cambridge Isotopes), amino acids, glyoxal trimer dihydrate (Sigma-Aldrich), used as received [M] |
| product study | 0.1-1.0 M glyoxal + amino acid in water; 1 µL droplets dried on ATR diamond (FTIR); 150 µL dried in vials to brown solids, redissolved 1 mg/mL (ESI-MS) or in D2O (NMR) [M] |
| bulk kinetics | **1H NMR in D2O at 295 K**: glyoxal 1 M + glycine 1 M, arginine 0.86 M or serine 0.48 M (solubility limits); two reactions per amino acid; initial loss rates from linear fits (all data for Gly and Ser; first 7 days for Arg) [M] |
| pH | **not measured during kinetics**; f_AA computed for "experimental pH 6" [M] — the later paper shows pH fell 1-2 units (sedehi2013) [!] |
| drying droplets | 0.03-10 mM glyoxal + amino acid, nebulised (~1.5 or ~5 µm), N2 chamber > 70-80 % RH, HR-ToF-AMS / Q-AMS after ~10-20 min [M] |
| aerosol uptake | solid glycine aerosol (84 nm mode), 50 % RH, 45 ppb gas-phase glyoxal; SMPS + Q-AMS [M] |

## 3. Quantitative results

**Table 1 (295 K, D2O, no drying)** [M]:

| amino acid | [glyoxal] M | [AA] M | f_glyoxal | f_AA | AA initial loss (mM d⁻¹) | k (M⁻¹ s⁻¹) |
|---|---|---|---|---|---|---|
| glycine | 1.0 | 1.0 | 7 × 10⁻⁴ | 1.7 × 10⁻⁴ (pH 6) | 1.9 ± 0.6 | 0.12 ± 0.04 |
| serine | 1.0 | 0.48 | 7 × 10⁻⁴ | 7.1 × 10⁻⁴ (pH 6) | 2.0 ± 0.2 | 0.09 ± 0.01 |
| arginine | 1.0 | 0.86 | 7 × 10⁻⁴ | 1 | 24 ± 3 | (4.5 ± 0.6) × 10⁻⁴ |

f_glyoxal: monohydrate fraction from ab initio thermodynamics (ref 21) [M, computed]. Rate law (eq 1): -d[AA]/dt =
k [Glyoxal]tot f_glyoxal [AA]tot f_AA. Glyoxal loss "should be the same" as AA loss, not measured [Q].

**Recomputation from the printed loss rates** [Z]: serine 0.097 (printed 0.09) and arginine 4.6 × 10⁻⁴ (printed 4.5 ×
10⁻⁴) reproduce; **glycine gives 0.185, not 0.12** (1.9 mM/d = 2.2 × 10⁻⁸ M/s over 1 × 7e-4 × 1 × 1.7e-4) — a x1.5
inconsistency inside Table 1 [!]. f_AA values reproduce from pKa 9.78 (Gly) and 9.15 (Ser) at pH 6 [Z].

**Total-basis constants (what was measured)** [Z]: glycine 1.9 mM/d at 1 M + 1 M → 2.2 × 10⁻⁸ M⁻¹ s⁻¹ (1.3 × 10⁻⁹
L/(mmol·min)) at 295 K and an unknown, falling pH. Serine 4.8 × 10⁻⁸ M⁻¹ s⁻¹; arginine 3.2 × 10⁻⁷ M⁻¹ s⁻¹.

**Other numbers** [M]: ornithine α : side-chain 80 : 20; glycine imidazole 60 ± 20 % and diimine 15 ± 5 % of the
glycine loss rate; arginine products 1 : 2 = 3 : 1, α-amine products ≥ 30x less abundant; drying 4 mM droplets → imidazole
(m/z 185) within ~10 min; at 32 µM the m/z 86 fragment is 2000x lower (overall order ~1.6); with f_glyoxal = 1
predicted reaction times within x6 of the ~30 min observed; solid glycine aerosol + 45 ppb glyoxal: +20 % volume in
30 min, glycine AMS peaks halved, uptake coefficient 4 × 10⁻⁴; SMPS wall-loss constant 0.0034 min⁻¹.

## 4. Mechanism evidence and its kind

- **Kind 4-5:** products identified (1H / 13C / 2-D NMR, ESI-MS, HR-ToF-AMS; 15N-glycine available for N assignment),
  product-to-reactant stoichiometry from matched NMR signals (a partial mass balance on glycine). **Exclusion test**
  for the Strecker-formaldehyde route: the masses it would require for Ser, Asp, Orn are absent (kind 3-like negative
  result).
- Mechanism (Scheme 2): nucleophilic attack of the unprotonated amine on glyoxal monohydrate → imine; diimine side
  product; imine + imine (80 % of the time, by NMR) → zwitterion → ring → dehydration triggers C-C rupture → imidazole
  + formic acid [P beyond the observed products].
- Order (second, equimolar) from 1/[AA] linearity; diimine not an intermediate (constant product ratio) [M/Q].
- Laboratory: De Haan (San Diego) + Tolbert / Jimenez (Boulder).

## 5. Bearing on open questions — and the transfer to a 100-140 °C food pot

**Transfer assumptions:**
1. **Use the PRODUCT and STOICHIOMETRY, not the rate.** The k values are invalidated by the pH drift (sedehi2013:
   true k x250 for glycine, x400 for serine) [M in sedehi2013].
2. The product claim was made at 295 K and 0.1-1 M. At 100-140 °C and mM glyoxal, competing glyoxal channels
   (Strecker on the amine, Cannizzaro, hydrate-driven oligomers, the engine's unassigned `k_go_sink`) change the
   partition; the imidazole route's share is not transferable without a high-T product study [Z].
3. **Aldehyde route excluded only for these conditions**: with added formaldehyde at 80 °C (cited ref 34 [C]) the
   odd carbon can come from formaldehyde; in a food pot, Strecker aldehydes (RCHO) are present and could supply it,
   giving 2-substituted imidazoles — not excluded at Maillard temperatures [Z].
4. The 15N-glycine was bought for NMR, but no 15N quantitation is reported [!].

**Per question:**
- **Pyrazine / Strecker precursors (evidence map 2.6-2.7):** adds a documented NON-Strecker glyoxal + amine sink
  (imidazole + formic acid) that removes both glyoxal and amino acid without producing an aminoketone [M/Z]. Supports
  the B18 flag "conditional_on_b13_dicarbonyl_sinks" as a real dependency; LEAVES the Strecker constants unchanged.
- **Engine `k_go_ak` ("glyoxal + glycine → aminoacetaldehyde + CO2 + HCHO"):** the competing imidazole channel is
  absent from the engine [Z]. Whether it matters at 100 °C is unknown (no Ea here; sedehi2013's class Ea 50 ± 11).
- **CF-069 (acetic / formic acid budget of the trunk):** a glyoxal-derived formic acid source with amines; LEAVES
  UNCHANGED (no rate at process temperature) but names a candidate formate route [Z].
- **Thiol losses / observation operators:** not touched.

## 6. Disagreements found

1. **haan2009 Table 1 vs sedehi2013 (same laboratory):** glycine k x250 and serine k x400 lower here; cause: pH
   assumed 6, while formic acid production lowered it by 1-2 units (the later paper's diagnosis) — operator /
   conditions. Use sedehi2013.
2. **Internal:** glycine k recomputes to 0.185 from the printed loss rate and f values, printed 0.12 (x1.5) [!].
3. **Odd-carbon origin vs older literature** (ref 33 [C]: formaldehyde from glycine's Strecker degradation): this paper
   argues the second glyoxal for its conditions (no 3-R products) — a mechanism disagreement settled here at 25 °C,
   not at Maillard temperatures.

## 7. Ledger-ready measurements

- **Stoichiometry rows** (two readers): glycine imidazole 60 ± 20 % and diimine 15 ± 5 % of glycine loss (text, p.
  2822); ornithine α : side-chain 80 : 20.
- **Product identity** (Schemes 2-3): a qualitative row "glyoxal + primary amino acid → 1,3-disubstituted imidazole +
  formic acid, 25 °C" for any future glyoxal-amine sink.
- **Do NOT key Table 1's k** except as a superseded value with the sedehi2013 correction attached.
