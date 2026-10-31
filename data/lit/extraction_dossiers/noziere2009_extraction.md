# Nozière, Dziedzic & Córdova 2009 — EXTRACTION (glyoxal 0.05-1 M in 0.5-4 M ammonium sulfate, chloride, fluoride, phosphate, pH 4.6-8.2 (unbuffered or phosphate / borate buffers), 298 K; LC-HRMS products; UV kinetics of a 209 nm product: second order in glyoxal, k_II rising exponentially with NH4+ activity and pH)

### An "iminium pathway" rate for glyoxal at molar ammonium and 25 °C only. No temperature dependence. The absolute k_II scales with an ASSUMED absorption cross-section (crotonaldehyde's), and the 209 nm product was later attributed to imidazole / imidazole-2-carboxaldehyde. Usable for food pots only as a form and an order of magnitude, with heavy caveats.

**Read 2026-10-31 from `data/articles/noziere2009.pdf`** (7 pp., J. Phys. Chem. A 113:231-237; born-digital, clean text
layer). Eqs 3-13, Table 1 and the printed k_II values read from the text layer and checked against 130 dpi renders
(scratchpad `dos7/noz-4.png`, `noz-5.png`); Figs. 3-4 (log k_II vs a_NH4+ and vs pH) read at ±0.2 dex. Table 2
(HRMS masses, proposed structures) read from the render. No SI statement.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| what the model needs | what Nozière 2009 gives | tag |
|---|---|---|
| rate law, glyoxal + NH4+ | product formation (209 nm) **first order in glyoxal in 20 wt % (NH4)2SO4** (k_I = (6.2 ± 1.6) × 10⁻⁶ s⁻¹), **second order in glyoxal at ≥ 25 wt % and in all other salts**: rate = k_II [product-equivalent glyoxal]² (eqs 5-6) | [M] |
| dependence on ammonium | **exponential in activity**, not linear in concentration: k_II = (1.7 ± 0.8) × 10⁻⁴ exp(0.94 a) (sulfate, eq 7); (2.0 ± 1.0) × 10⁻⁴ exp(0.74 a) (chloride, eq 8); (0.5 ± 0.3) exp(1.28 a) ((NH4)2HPO4, eq 9); all salts, common pH: A exp((1.5 ± 0.8) a) (eq 10) | [M] |
| pH dependence | k_II = f(a) exp((2.5 ± 1.0) pH) (eq 11) = [H⁺]^-1.086 [M]; **abstract prints exp((2.5 ± 0.2) pH)** [!] | [M] |
| combined law | **k_II (M⁻¹ s⁻¹) = (2 ± 1) × 10⁻¹⁰ exp(1.5 a_NH4+) exp(2.5 pH)** (eq 12), 298 K | [M] |
| temperature / Ea | **none** (298 ± 2 K only) | [!] |
| counterion / phosphate | k_II(HPO4²⁻ salt) > F⁻ > SO4²⁻ ≈ Cl⁻, "mostly ... attributed to different pH" [Q]; phosphate-salt solutions were at pH 8.0-8.2 | [M] / [Q] |
| absolute scale | k_II ∝ 1 / σ(209 nm); σ ASSUMED = crotonaldehyde's 37 150 M⁻¹ cm⁻¹; k_II uncertainty 50 %, k_I 25 % [M] | [M] / [!] |
| what the product is | 209 nm band: "π → π* transition of a C=C or C=N bond", second-generation product of two glyoxal molecules [P]; HRMS: mostly acetal oligomers (Brønsted pathway) plus N-containing m/z 192.0503 [M]; the 209 / 207 nm band was later assigned to imidazole + imidazole-2-carboxaldehyde (Yu 2011, cited in sedehi2013) [C] | [P] / [C] |

**Answer.** For a food pot the paper gives a FORM (two glyoxal molecules per product, an ammonium-catalysed path
whose rate climbs about one decade per pH unit) and an order of magnitude at molar ammonium, 25 °C. Its own
combined law cannot be taken to dilute ammonium: eq 12 at a_NH4+ = 0, pH 6 still returns 6.5 × 10⁻⁴ M⁻¹ s⁻¹, while
the blank buffers without ammonium gave no 209 nm product [M/Z]. The exponential-in-activity form is a fit over a =
0.3-2.6 M, not a mechanism. Confidence ~0.9 that no constant from this paper should enter the engine; the
De Haan-group laws (sedehi2013, first order in glyoxal and in free NH3) are the transferable description of the same
chemistry.

## 0. Identity

| field | value |
|---|---|
| Title | "Products and Kinetics of the Liquid-Phase Reaction of Glyoxal Catalyzed by Ammonium Ions (NH4+)" |
| Authors | Barbara Nozière* (Applied Environmental Sciences), Pawel Dziedzic, Armando Córdova* (Organic Chemistry), Stockholm University |
| Venue | *J. Phys. Chem. A* **2009**, 113 (1), 231-237; received 3 Sep 2008, revised 6 Nov 2008, web 9 Dec 2008 (Crossref "issued" 2008) |
| DOI | **10.1021/jp8078293** (printed; Crossref 2026-10-31: same title, Nozière / Dziedzic / Córdova, vol. 113, 231-237). Matches the brief |
| Related | Nozière & Córdova 2008 JPCA 112:2827 (amino-acid catalysis, ref 9); Nozière, Dziedzic & Córdova 2007 GRL (ref 8) — not on disk |

## 2. Conditions

| item | as printed |
|---|---|
| glyoxal | 40 % aqueous (~8.8 M), Aldrich; **0.05-1 M, generally 0.1 M**, added by microsyringe [M] |
| salts | (NH4)2SO4 20-43 wt % (1.6-3.5 M, saturated), NH4Cl 2-4 M, NH4F 0.5-2 M, (NH4)2HPO4 0.5-1.33 M [M] |
| buffers | Fixanal: pH 7 **KH2PO4 / Na2HPO4**; pH 8 and 9 Na2B4O7 / HCl; pH 10 Na2B4O7 / NaOH [M] |
| pH (Table 1, ± 0.2) | unbuffered sulfate 4.9-5.8; sulfate in pH 7 buffer 6.1-6.2, in pH 8 buffer 7.0-7.2, pH 9 7.7, pH 10 8.1; NH4F 6.3-8.1; (NH4)2HPO4 8.0-8.2; NH4Cl 4.6-5.7 [M] |
| activity coefficients | γNH4+ 0.42-0.46 (sulfate, AIM model), 0.61-0.64 (chloride, AIM); fluoride and phosphate ASSUMED 0.60 (dilute NH4Cl) [M] |
| temperature | **298 ± 2 K** only; 4 mL bulk, stirred, dark [M] |
| analysis | UV-vis (Agilent 8453), 1 mm cells, 190-1100 nm, absorbance at 209 nm; LC-HRMS (Bruker MicrOTOF) for products [M] |
| blanks | buffers without ammonium (pH 7-10): no 209 nm product [M] |
| scale | "over a hundred experiments" [M] |

## 3. Quantitative results

**Printed rate constants** [M]: k_I(20 wt % (NH4)2SO4) = (6.2 ± 1.6) × 10⁻⁶ s⁻¹; k_II((NH4)2SO4 25 wt %) = (8.0 ±
4.0) × 10⁻⁴; k_II(NH4Cl 2 M) = (6.2 ± 2.0) × 10⁻⁴; k_II(NH4F 0.5 M) = (3.7 ± 2.0) × 10⁻³; k_II((NH4)2HPO4 0.5 M) =
(1.2 ± 0.6) M⁻¹ s⁻¹. Check [Z]: eq 7 at a = 1.68 (2 M sulfate → 4 M NH4+, γ 0.42) gives 8.2 × 10⁻⁴, matching the
25 wt % value.

**Eqs 7-12** as in sec. 1. Eq 12 checks against the paper's own atmospheric numbers [Z]: a = 4, pH 4 → 1.8 × 10⁻³ (printed
2 × 10⁻³); pH 6 → 0.26 (printed 0.26).

**Fig. 4 (log k_II vs pH, fixed salt; re-read at 300 dpi, `nozfig4-5.png`)** [F, ±0.15 dex]: sulfate 35 wt % (blue):
0.002 (pH 5.5, unbuffered), **0.14 (pH 6.1, phosphate pH 7 buffer)**, **0.0055 (pH 7.0, borate pH 8 buffer)**, 2.6
(pH 7.6, borate pH 9), 5.5 (pH 8.1, borate pH 10); dashed fit ~1.15 dex / pH. NH4F 0.5 M (pink): 0.0035 (pH 6.4),
0.009 (pH 6.7, phosphate buffer), 0.04-0.09 (pH 7.5-7.7), 0.12-0.37 (pH 8.1); on its line. NH4Cl 4 M (green):
0.0012-0.002 (pH 4.8), 0.075 (pH 5.7). The authors' anomaly ("ammonium sulfate solutions at pH = 7", cause "not
clear" [Q]) is the sulfate series in the pH 8 BORATE buffer, ~x20 below the line (Fig. 3 blue "pH ~7.0-7.2",
2-6 × 10⁻³). The sulfate series in PHOSPHATE buffer (pH 6.1) sits ~x10 ABOVE the line, while the fluoride series in
the same phosphate buffer (pH 6.7) is on its line [F/Z]. Borate complexes diols (glyoxal hydrates) — a candidate cause
the authors do not raise [Z].

**Atmospheric estimates (not used)** [M]: a = 4, pH 4 / 6 → k_II 2 × 10⁻³ / 0.26; first-order 7 × 10⁻⁷ / 9 × 10⁻⁵ s⁻¹
at 3.5 × 10⁻⁴ M glyoxal; uptake γ_r 0.6-7 × 10⁻³.

## 4. Mechanism evidence and its kind

- **Kind 4 (product formation monitored) + kind 6 (a_NH4+, pH, counterion varied within study)** + product
  identification by HRMS (mostly acetal oligomers identical to acid catalysis; some N-containing ions). Nozière /
  Córdova, Stockholm.
- Two pathways (Scheme 1): Brønsted-acid (NH4+ as proton donor → acetal oligomers) and iminium (NH3 + protonated
  glyoxal → iminium → heterocycles / oligomers) [P].
- Second order in glyoxal for the 209 nm product → a product built from two glyoxal molecules [Q/P]; consistent with
  haan2009's NMR finding that the 1,3-disubstituted imidazole's odd carbon comes from a second glyoxal [Z, cross-paper].
- The iminium path is "a lower limit" for ammonium-catalysed glyoxal consumption [Q].

## 5. Bearing on open questions — and the transfer to a 100-140 °C food pot

**Transfer assumptions needed (all unverified):**
1. **No temperature dependence measured.** Any use above 25 °C needs a borrowed Ea (sedehi2013's ammonium value 18-21
   ± 8-15 kJ/mol is the nearest, same chemistry, another lab).
2. **Ammonium regime.** a_NH4+ 0.3-2.6 M here; food NH4+ is 0-50 mM (deamidation, Strecker). The exp(1.5 a) form must
   be replaced by a law that vanishes at zero ammonium; the paper gives none [!]. Sedehi's k[Glx][NH4+]tot f_NH3 does.
3. **Glyoxal order.** Second order in glyoxal at 0.1 M; at food glyoxal (µM-mM) a first-order step can become
   rate-limiting (the paper saw the first-order regime at 20 wt % sulfate) [M/Z].
4. **Absolute scale** tied to an assumed σ(209 nm) and to a product later re-identified; rescale if Yu 2011's
   imidazole cross-sections are used [!].
5. **pH**: 4.6-8.2, compatible with food pots; the slope (~1.1 dex / pH) agrees in size with sedehi2013's AS slope
   (1.05-1.34).
6. **Ionic strength** 1.6-10 m; salting effects on glyoxal hydration not separated [!].

**Per question:**
- **Pyrazine / Strecker precursors, ammonium as an amine:** LEAVES UNCHANGED (no Strecker observable; the iminium path
  makes imidazoles, a dicarbonyl sink that does NOT feed pyrazines) [Z].
- **Phosphate catalysis (evidence map 1.2):** LEAVES UNCHANGED. The authors attribute the (NH4)2HPO4 excess mostly
  to pH (eq 10) [Q]; but in Fig. 4 the sulfate series in phosphate buffer sits ~x10 above its pH line while the
  fluoride series in the same buffer does not [F]. Ambiguous, one point each, 298 K.
- **Engine `k_go_sink` (disputed):** an ammonium path is a candidate contributor only at high ammonium; negligible at
  food NH4+ by this paper's own numbers (eq 12 is not extrapolable; sec. 1).
- **CF ids:** none touched.

## 6. Disagreements found

1. **Abstract vs eqs 11-12:** pH coefficient 2.5 ± 0.2 (abstract) vs 2.5 ± 1.0 (eqs 11-12) [!].
2. **Functional form vs sedehi2013 / Yu 2011 (same reaction, glyoxal + ammonium sulfate):** Nozière: second order in
   glyoxal, exponential in a_NH4+, ~[H⁺]^-1.09. Sedehi (formic acid and glyoxal loss, NMR): first order in glyoxal and
   in total ammonium x free-NH3 fraction, slope 1.05 per pH. In sedehi2013's Fig. 2 Nozière's points lie on the
   imidazole-based fit (slope 1.34 ± 0.08, with Yu's imidazole data) [C]. Candidate causes: different observables
   (one product's UV band vs reactant loss), the molar glyoxal here (0.1 M; Sedehi cites a regime change when
   [Glx][Am] > 1.2 M²), and the assumed cross-section. For a food pot prefer Sedehi's form.
3. **Buffer anomalies** [F/Z]: the low-k "pH 7" sulfate series is in borate buffer (measured 7.0-7.2), not in the
   phosphate pH 7 buffer; the phosphate-buffered sulfate point (pH 6.1) is high by ~x10. The authors attribute
   neither to the buffer; for a food model, buffer identity is an uncontrolled factor in this data set.

## 7. Ledger-ready measurements

- **Printed k_II values** (four salts) and **eqs 7-12**, keyed as "298 K, molar NH4+, assumed σ(209 nm)" — context
  rows only; no fit should read them.
- **Table 1** (pH and γNH4+ per solution) — needed to re-express any of Fig. 3-4 on a free-NH3 basis if this route is
  ever compared with sedehi2013.
