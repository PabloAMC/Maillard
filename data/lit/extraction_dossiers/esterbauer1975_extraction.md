# Esterbauer, Zöllner & Scholz 1975 — EXTRACTION (glutathione + 21 α,β-unsaturated carbonyls (2-alkenals, 4-hydroxy-2-alkenals, enones, acrylyl derivatives): forward second-order k₁, reverse first-order k₂ and dissociation K, 66 mM phosphate pH 7.4, 20 °C; pH-rate profile 3.5-10.5 for GSH + crotonaldehyde; buffer catalysis for GSH, cysteine and ethyl thioglycolate)

### THE ONE PAPER IN THIS SET WITH FORWARD, REVERSE AND K FOR A THIOL TRAP. Michael addition of a thiol to an enal is REVERSIBLE: simple 2-alkenal adducts dissociate with t½ 4-6 h at 20 °C, 4-hydroxyalkenal adducts in 6-19 days. Rate ∝ thiolate fraction, buffer (general-acid) catalysed for thiols without a neighbouring NH₃⁺, NOT for cysteine (which reacts ~60x faster than GSH unbuffered). One temperature (20 °C), no Ea. Thiol is GSH (cysteine and ethyl thioglycolate in one table), not MFT.

**Read 2026-11-01 from `data/articles/esterbauer1975.pdf`** (8 pp., Z. Naturforsch. 30c, 466-473; scanned and
re-typeset by the publisher, "Download Date 7/3/15"). Prose from the `pdftotext -layout` layer; **Tables II-IV and
Fig. 1 checked against 150 dpi renders, Fig. 1 against a 300 dpi crop** (scratchpad `dos6/est-3.png`, `est-4.png`,
`estf1-4.png`). Fig. 1 values read by me ±0.1 log unit.

**Legend.** **[V]** verbatim · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative claim ·
**[C]** cited · **[P]** proposed · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

The model questions: the MFT loss law (engine `k_mft_decay`, first order, 0.095 min⁻¹ at 110 °C; thiol-rivals A2
"second order in a growing electrophile pool", A3 "suppressed by free cysteine"; the patent bound λ(110 °C) ≲
0.02-0.05 min⁻¹), the CF-076 E3 "second cysteine sink", and CF-031 (lipid vs meaty thiols).

| model need | what Esterbauer 1975 gives | tag |
|---|---|---|
| forward k (thiol + electrophile) | GSH, pH 7.4, 20 °C, 66 mM phosphate: acrolein 1.21 × 10², MVK 31.9, 4-hydroxy-2-alkenals 1.09-2.44, 4-keto-2-pentenoic acid 3.42, crotonaldehyde 0.785, 2-pentenal 0.471, 2-hexenal 0.330, cyclohexenone 0.336, citral 3.23 × 10⁻², ethyl crotonate 3.10 × 10⁻³, mesityl oxide 2.30 × 10⁻³, acrylamide 0.013, acrylic acid 0.001 M⁻¹ s⁻¹; crotonate, 2-pentenoate no reaction (Tables II, III) | [M] |
| reverse k (adduct → thiol + electrophile) | first order: 2-alkenals 3.0-4.7 × 10⁻⁵ s⁻¹ (t½ 4-6 h); acrolein 1.76 × 10⁻⁶ (4.6 d); 4-hydroxyalkenals 4.1-12.6 × 10⁻⁷ (6-19 d); MVK 1.19 × 10⁻⁷; 3-methyl-3-buten-2-one 9.0 × 10⁻³ (calculated) | [M] |
| K (dissociation, (GSH)(Carb)/(Add)) | 1.40 × 10⁻⁸ M (acrolein) to 1.65 × 10⁻² M (ethyl crotonate); 2-alkenals 4.2 × 10⁻⁵-1.4 × 10⁻⁴ M; 4-hydroxyalkenals 2.5-8.9 × 10⁻⁷ M | [M] |
| pH law | GSH + crotonal, pH-stat, no buffer (Fig. 1): k₁ ≈ 4 × 10⁻⁴ (pH 3.5), 1.6 × 10⁻³ (4.0), 9 × 10⁻³ (5.0), 1.35 × 10⁻² (5.4), 2.4 × 10⁻² (6.0), 3.4 × 10⁻² (6.5), 0.10 (7.0), 0.29 (7.5), 0.60 (8.0), 1.3 (8.5), 1.7 (9.0), 1.9 (10.0), 2.3 (10.5) M⁻¹ s⁻¹; "increases with pH approx. parallel with α_SH" (pKa_SH 8.56) | [F] / [Q] |
| buffer catalysis | crotonal, pH 7.4, 20 °C (Table IV): GSH none 0.27 → 0.1 M phosphate 1.25 (x4.6); 0.066 M phosphate 0.78; Tris 0.83, triethanolamine 0.53, borate 0.66, NH₄Cl 0.48, acetate 0.25, formate 0.27, 0.5 M NaCl 0.28. Ethyl thioglycolate 0.10 → 2.00 (x20). **Cysteine 16.5 → 16.7 (none)** | [M] |
| Ea | **none** (20 ± 1 °C only) | — |
| mechanism of the pH / buffer law | eq (2)-(3): thiolate adds (k_a), the enolate intermediate is protonated by H₃O⁺, water or HX; below pH 4.5 the addition limits (β = 1), above pH 6.0 the proton transfer limits (β < 1, buffer-dependent); cysteine's NH₃⁺ protonates the intermediate intramolecularly | [Q] / [P] |

**Answers.**

1. **MFT loss law, structure: a Michael-type thiol capture is a REVERSIBLE partition, not a first-order irreversible
   decay** [M → Z, transfer]. At 20 °C the simple-enal adducts revert in hours; at 100-145 °C (Ea unknown) they would
   revert in minutes or less, so the capture approaches its equilibrium: free / total thiol = K/(K + [E]). With the 20 °C
   2-alkenal K (4 × 10⁻⁵-1.4 × 10⁻⁴ M) and [E] = 1 mM, 88-96 % would be bound; at 0.1 mM, 40-70 % [Z]. A pot-independent
   first-order λ (A0) cannot represent this; the "reversible sink" family (unlisted in the thiol-rivals prereg, cf. Guth
   1995's reverting disulfide) can. Transfer assumptions: MFT behaves as an alkyl thiolate nucleophile; the pot holds
   Michael acceptors at 0.1-1 mM (not shown for any Maillard pot on disk); K's temperature dependence unknown
   (addition is exothermic in general, so K_diss rises with T: binding weaker at process temperature).
2. **MFT loss law, pH law** [F → Z]: the rate follows the thiolate fraction: ×~2.5 per pH unit from 5 to 6.5 and
   ×~6-8 per unit from 6.5 to 8.5 for GSH. For any thiol with pKa above the pot pH, Michael capture slows ~10x (GSH)
   from pH 7 to pH 5. The engine's `k_mft_decay` has no pH term; this predicts one (direction: faster at higher pH),
   for MFT only if its thiol pKa is above ~6 (not in this paper).
3. **MFT loss law, buffer** [M → Z]: thiols without a neighbouring NH₃⁺ (GSH, ethyl thioglycolate; MFT is one) are
   accelerated x4.6-20 by 0.1 M phosphate at pH 7.4; the effect needs pH > 6 (proton-transfer-limited regime) [Q].
   At pH 5 (the Hofmann / Schieberle 0.5 M phosphate recipes) the addition limits and buffer should matter little.
4. **A3 (cysteine protects MFT) has a competition basis** [M → Z]: cysteine reacts with crotonal 16.5 M⁻¹ s⁻¹ unbuffered,
   61x GSH and 165x ethyl thioglycolate; in 0.1 M phosphate still x8 over ethyl thioglycolate. At 30-100 mM cysteine
   against µM MFT, cysteine takes > 99.9 % of any shared electrophile flux. Competition is a mechanism for A3's
   "suppressed by free cysteine" alongside the redox one; the paper does not test MFT.
5. **λ magnitude: no bound** without an Ea. λ = k₁[E]: the engine's 0.095 min⁻¹ (1.6 × 10⁻³ s⁻¹) at 20 °C would need
   ~13 µM acrolein or ~2-5 mM of a 2-alkenal [Z]; at 110 °C the requirement drops by the unknown Arrhenius factor.
6. **CF-076 E3 (second cysteine sink)** [Z]: cysteine + crotonal 0.99 L mmol⁻¹ min⁻¹ at pH 7.4 (20 °C) sits at the top of
   E3's declared k_s band (10⁻⁴-10⁻¹ L mmol⁻¹ min⁻¹); at pH 5.5 the thiolate law gives ~1/15 of it. But E3's electrophiles
   in the Zhai pots are α-dicarbonyls, which take cysteine by hemithioacetal / thiazolidine chemistry (Schubert 1936,
   Kallen 1971b), not Michael addition: the constants do not transfer. Leaves E3 unchanged.
7. **`k_ttca_cys`, E1 / E4**: nothing (no thiazolidine).

## 2. Identity

| field | value |
|---|---|
| Title | "Reaction of Glutathione with Conjugated Carbonyls" |
| Authors | Hermann Esterbauer, Helmward Zöllner, Norbert Scholz — Institut für Biochemie, Universität Graz |
| Venue | *Z. Naturforsch.* **30c**, 466-473 [1975]; received January 17 / April 1, 1975 |
| DOI | not printed. Crossref (2026-11-01): **10.1515/znc-1975-7-808**, same title, Esterbauer / Zöllner / Scholz, Z. Naturforsch. C 30(7-8) 466-473. Matches the brief |
| Funding | Fonds zur Förderung der wissenschaftlichen Forschung, Wien |
| Related | Esterbauer 1970, Mh. Chem. 101, 782 (forward method, ref 29; not on disk) |

## 3. Conditions

- **Forward** [M]: 66 mM phosphate pH 7.4, 20 ± 1 °C; GSH + carbonyl 0.05 + 0.05 mM (acrolein, MVK), 0.1 + 2.0 mM (most),
  0.1 + 1.0, 10 + 10, 5 + 5, 0.1 + 40, 0.2 + 8.0 mM (per reaction, Table II caption); 10 % ethanol for reactions 5, 12,
  13 ("and 24" — no reaction 24 exists [!]). GSH by DTNB at 412 nm, or carbonyl by UV 222 nm (ε 13,600). SD 5-10 % for
  k₁, 3-5 % for k₂, 10-15 % for K; ≥ 3 determinations.
- **Second-order check (Table I)**: 4-hydroxypentenal, GSH 0.1-0.3 mM × aldehyde 0.1-2.0 mM: k₁ 2.08-2.30, mean 2.19 ±
  0.08 M⁻¹ s⁻¹ [M].
- **Reverse and K** [M]: GSH + carbonyl equilibrated (0.1 + 0.2 to 10 + 20 mM) in 66 mM phosphate pH 7.4 (from 266 mM
  stock), N₂-degassed for experiments 6-15 and 20 (stock solutions only for volatile carbonyls); then DTNB in excess
  traps released GSH; k₂ = v/(a₀ − a) from the initial slope over 0.5-5 h. k₂ independent of DTNB and of adduct and
  crotonaldehyde concentration; p-chloromercuribenzoate gave reagent-dependent rates and was not used [Q].
- **pH profile (Fig. 1)**: pH-stat with HCl / NaOH, no buffer; GSH 0.1 mM (pH ≥ 7) or 1.0 mM (pH < 7), crotonal 20x;
  20 ± 1 °C.
- **Buffer table (Table IV)**: pH 7.4, 20 °C, thiol 0.1 mM, crotonal 2.0 mM.

## 4. Quantitative results

**Table II** (pH 7.4, 66 mM phosphate, 20 °C) [M]; * = calculated from the other two:

| # | carbonyl | k₁, M⁻¹ s⁻¹ | k₂, s⁻¹ | K, M | % free GSH (0.1 + 0.1 mM) |
|---|---|---|---|---|---|
| 1 | acrolein | 1.21 × 10² | 1.76 × 10⁻⁶ | 1.40 × 10⁻⁸ | 1.2 |
| 2 | crotonaldehyde | 7.85 × 10⁻¹ | 3.01 × 10⁻⁵ | 4.16 × 10⁻⁵ | 46.8 |
| 3 | 2-pentenal | 4.71 × 10⁻¹ | 3.44 × 10⁻⁵ | 7.38 × 10⁻⁵ | 56.0 |
| 4 | 2-hexenal | 3.30 × 10⁻¹ | 4.70 × 10⁻⁵ | 1.42 × 10⁻⁴ * | 67.0 |
| 5 | citral | 3.23 × 10⁻² | 5.00 × 10⁻⁶ * | 1.63 × 10⁻⁴ | 70.5 |
| 6 | 4-hydroxy-2-pentenal | 2.19 | 5.60 × 10⁻⁷ | 2.96 × 10⁻⁷ | 5.2 |
| 7 | 4-hydroxy-2-hexenal | 1.56 | 4.10 × 10⁻⁷ | 2.63 × 10⁻⁷ * | 5.0 |
| 8 | 4-hydroxy-2-heptenal | 1.83 | 4.57 × 10⁻⁷ | 2.50 × 10⁻⁷ * | 4.9 |
| 9 | 4-hydroxy-2-octenal | 1.74 | 7.49 × 10⁻⁷ | 4.89 × 10⁻⁷ | 6.7 |
| 10 | 4-hydroxy-2-nonenal | 1.09 | 9.60 × 10⁻⁷ | 8.85 × 10⁻⁷ * | 8.5 |
| 11 | 4-hydroxy-2-decenal | 1.96 | 9.50 × 10⁻⁷ | 4.85 × 10⁻⁷ * | 4.7 |
| 12 | 4-hydroxy-2-undecenal | 1.47 | 9.36 × 10⁻⁷ | 6.38 × 10⁻⁷ | 7.7 |
| 13 | 4-hydroxy-2-dodecenal | 2.44 | 1.26 × 10⁻⁶ | 5.42 × 10⁻⁷ * | 7.1 |
| 14 | 4-hydroxy-4-isopropyl-2-pentenal | 3.57 × 10⁻² | 7.76 × 10⁻⁷ | 2.31 × 10⁻⁵ | 38.0 |
| 15 | 4-ethoxy-2-pentenal | 1.83 | 1.75 × 10⁻⁵ | 9.60 × 10⁻⁶ | 26.4 |
| 16 | methyl vinyl ketone | 3.19 × 10 | 1.19 × 10⁻⁷ | 3.60 × 10⁻⁹ | 0.6 |
| 17 | 3-methyl-3-buten-2-one | 6.00 × 10⁻¹ | 9.00 × 10⁻³ * | 1.50 × 10⁻² | 99 |
| 18 | mesityl oxide | 2.30 × 10⁻³ | 1.12 × 10⁻⁵ | 4.80 × 10⁻³ | 98 |
| 19 | 2-cyclohexen-1-one | 3.36 × 10⁻¹ | 5.16 × 10⁻⁸ | 1.53 × 10⁻⁵ * | 32 |
| 20 | 4-keto-2-pentenoic acid | 3.42 | 2.32 × 10⁻⁶ | 6.07 × 10⁻⁷ | 7.5 |
| 21 | ethyl crotonate | 3.10 × 10⁻³ | 5.75 × 10⁻⁵ | 1.65 × 10⁻² | 99 |
| 22 | crotonate, 2-pentenoate | no measurable reaction | | | |

**Table III** (acrylyl derivatives, k₁, M⁻¹ s⁻¹) [M]: acrolein 121, MVK 31.9, ethyl acrylate 0.243, acrylamide 0.013,
acrylic acid 0.001. "aldehyde > ketone > ester > amide > carboxylate", ~five orders of magnitude; one alkyl on the double
bond ~x100 slower, two ~x10⁴ [Q].

**Table V** (derived by the authors) [M]: adduct-formation t½ (0.1 + 0.1 mM): acrolein 1.4 min, MVK 5.2 min,
4-hydroxyalkenals 1.1-2.5 h, 2-alkenals 3.6-8.5 d; reverse t½ acrolein 4.6 d, MVK 60.7 d, 4-hydroxyalkenals
6.3-19.2 d, crotonal / pentenal / hexenal 4.0-6.3 h, 3-methyl-3-buten-2-one 1.2 min.

**Fig. 2** [M]: example reverse run (2.5 mM GSH + 5.0 mM 4-HNE, 3.5 h, then DTNB), k₂ = 10 × 10⁻⁷ s⁻¹.

## 5. Mechanism evidence and its kind

- Graz laboratory (Esterbauer). Kind 6 (one factor varied within the study: carbonyl structure, pH with pH-stat, buffer
  identity and concentration, thiol identity); kind 3 (DTNB as a trap driving the reverse reaction to completion,
  checked for reagent independence) [M/Q].
- Second-order forward / first-order reverse orders shown (Table I; reverse independent of adduct concentration) [M].
- The intramolecular-proton-transfer explanation of cysteine's buffer independence is proposed [P].
- No product identification beyond stoichiometry; no temperature series; all forward rates at pH 7.4.

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| MFT loss law (evidence map 2.12 "first-order loss law and its O₂ / carbonyl dependence") | **CONTRADICTS a pot-independent irreversible first-order law IF Michael capture dominates; SUPPORTS A2's carbonyl dependence in kind** | reversible, K-set partition; rate ∝ [E], thiolate fraction, buffer (pH > 6). Analogue only (GSH, 20 °C) |
| thiol-rivals A3 (cysteine suppresses loss) | **SUPPORTS a competitive mechanism** | cysteine ×61-165 faster per mol than GSH / ethyl thioglycolate unbuffered |
| CF-031 (lipid lowers vs raises meaty thiols) | SUPPORTS feasibility of "lipid aldehydes as thiol sinks" for UNSATURATED aldehydes only | 2-alkenals and 4-hydroxyalkenals capture thiols; saturated aldehydes (hexanal) are not Michael acceptors (not tested here); 2-alkenal capture reverts in hours |
| CF-076 E3 | LEAVES UNCHANGED | wrong electrophile class for the Zhai pots |
| `k_ttca_cys`, E1 / E4 | LEAVES UNCHANGED | |

## 7. Disagreements found

1. **With Doorn & Petersen 2002 (`doorn2002_extraction.md`)**: 4-HNE + GSH k₁ 1.09 M⁻¹ s⁻¹ here (20 °C, 66 mM phosphate)
   vs 1.33 ± 0.05 there (23 °C, 50 mM phosphate). Agreement within 22 %, a cross-lab check passed (Doorn cites it).
2. Internal [Z][!]: MVK reverse t½ from k₂ = 1.19 × 10⁻⁷ s⁻¹ is 67 d, Table V prints 60.7 d (k₂ would be 1.32 × 10⁻⁷);
   other rows agree. Table II caption cites a reaction 24 that does not exist.
3. With the engine: `k_mft_decay` is pH-, buffer- and electrophile-flat and irreversible; this paper's class chemistry is
   none of these (analogue; sec. 1 transfer assumptions).

## 8. Ledger-ready measurements

1. **Table II** k₁, k₂, K for the 2-alkenals, 4-hydroxyalkenals, acrolein, MVK (GSH, pH 7.4, 20 °C): a class table for
   any thiol-capture band (A2-type k₂, reversibility), printed.
2. **Fig. 1** k₁(pH) for GSH + crotonal (13 points, two readers): the thiolate pH law.
3. **Table IV** buffer and thiol identity effects: phosphate catalysis for NH₃⁺-free thiols, none for cysteine.
