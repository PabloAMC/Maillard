# Karpinski, Karny & Kublik 1990 — EXTRACTION
DOI: 10.1016/0022-0728(90)80030-a (not printed in the PDF; Crossref 2026-11-01: "Electrochemical studies of the decomposition of thiazolidine-4-carboxylic acids in aqueous solutions", Karpinski, Karny, Kublik, J. Electroanal. Chem. Interfacial Electrochem. 288(1-2), 129-141, 1990-08. **Matches the brief.**)
### Thiazolidine-4-carboxylic acid (Thz = the formaldehyde ring, Kallen's "TC"), its 2-propyl derivative (PrThz) and **the L-arabinose ring, 2-(L-arabino-tetrahydroxybutyl)-Thz (AThz), a pentose-cysteine thiazolidine and a C-2 diastereomer family member of TTCA**, decomposing to free L-cysteine in deaerated buffers at pH 4.5-8.0, **25 °C only**, followed by normal-pulse and square-wave voltammetry on a hanging mercury drop (ring compound and cysteine measured together, mass balance <= 5 %). Pseudo-first-order release constants (Table 1), plus air-saturated arms. **No ring pKa, no temperature series, no Ea, no equilibrium constant.**

**Read 2026-11-01 from `data/articles/karpinski1990.pdf`** (13 pp., scanned with an OCR layer that garbles exponents
and names). Prose from `pdftotext -layout`; **Table 1 (p. 139) and Figs. 7-8 (p. 137) read from 200 dpi renders**
(scratchpad `ing5/kp-09.png`, `kp-11.png`): the Table 1 header is "10^5 k_a / s", so every entry is x 10⁻⁵ s⁻¹
(the OCR layer prints "10' k,/s"). Fig. 8 read error +-0.05 in ln c.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Chemistry, University of Warsaw (Kublik; electroanalysis). The thiazolidines were made
by Jan Radomski (acknowledged; ref 23 = Radomski & Temeriusz 1989, `radomski1989_extraction.md`, same department,
which proposed `warsaw_temeriusz`). **Not in `lab_ids.yml`: NEW.** Under merge rule 2 (one institute, interlinked
authors: Karpinski & Radomski co-authored ref 19) the conservative choice is ONE Warsaw lab id for both dossiers.

## 1. What this paper answers for the model

CF-076's provisional admission (in-situ pH with Kallen-shaped K(pH) and k_r(pH)) declares k_r(100 °C, pH 4.5) in
[2e-4, 8e-4] min⁻¹ and rests on two unmeasured quantities: phi, the factor that transfers TC's release rate to TTCA
("not identified by any data on disk", cluster-5 verdict sec. 3), and TTCA's ring-N pKa.

| model need | what Karpinski 1990 gives | tag |
|---|---|---|
| release rate of a CYSTEINE-PENTOSE thiazolidine | **AThz: 0.3 × 10⁻⁵ s⁻¹ at pH 6.0, 0.08 × 10⁻⁵ s⁻¹ at pH 7.4** (25 °C, argon) = 1.8 × 10⁻⁴ and 4.8 × 10⁻⁵ min⁻¹ [Z units]; RSD "reached 30%" | [M] |
| phi (pentose ring vs TC, same lab, same method) | Thz 0.3 (pH 6.0) and 0.06 (pH 7.4): **AThz / Thz = 1.0 and 1.3** -> phi ~0.6-1.9 within the 30 % RSDs [Z] | [M] / [Z] |
| k_r(pH) shape | PrThz (best data, RSD < 10 %), 10⁵ k_a/s: pH 4.5 **2.1**, 5.2 **4.7**, 6.0 **7.3**, 6.7 **4.3**, 7.4 **1.6**, 8.0 **0.03** — "a maximum at pH near 6" [V]. AThz and Thz only at 6.0 and 7.4 (x3.75 and x5 slower at 7.4) | [M] |
| what is measured | overall release: cysteine appears as the ring compound disappears, sum = c0 within 5 %; no open-ring thiol intermediate detected (LOD ~1 uM) | [M] |
| ring opening vs hydrolysis | the authors conclude ring opening is rate-determining [P]; the data cannot separate it from a fast ring-chain pre-equilibrium followed by slow hydrolysis (Fife 1991) [Z][!] | [P] / [Z] |
| ring pKa | **none** | — |
| temperature / Ea | **none**: 298 +- 0.2 K only | — |
| O2 | air speeds loss and makes cystine: Thz pH 6.0 0.3 -> **1.3**; PrThz pH 6.0 7.3 -> **14.8**; PrThz pH 7.4 1.6 -> **2.3** (10⁻⁵ s⁻¹) | [M] |
| concentration dependence | PrThz k_a independent of c0 from 10⁻⁵ to below 10⁻³ M; "almost a factor of two lower" at 10⁻³ M | [M] |

**Answers.**
1. **Rate constants vs pH / T:** Table 1 above (25 °C only). Mechanism: release of L-cysteine identified by CV / OSW
   matching authentic cysteine; open-chain species below detection; ring opening called rate-determining [P].
2. **Ring pKa: not given** (confidence ~98 %).
3. **The new datum: phi ≈ 1 for a pentose ring** (same lab, same method, same day-type experiment): AThz releases
   cysteine as fast as unsubstituted TC at pH 6.0 and 7.4 [M / Z]; the 2-propyl ring is 24x / 20x faster. A
   polyhydroxyalkyl C-2 substituent therefore behaves like H, not like an alkyl group. Confidence that phi(25 °C,
   pH 6-7.4) lies in 0.5-2: ~0.75 (30 % RSDs on two ratios, ~10 % conversion in 8 h, and AThz's waves overlap
   cysteine's).
4. **Consequence for the CF-076 band** [Z, transfer-laden]: AThz at pH 4.5 is not measured. With PrThz's pH shape
   (k(4.5) / k(6.0) = 0.29), AThz(25 °C, pH 4.5) ≈ 5.2 × 10⁻⁵ min⁻¹. Reaching the declared k_r(100 °C, pH 4.5) =
   2-8 × 10⁻⁴ min⁻¹ needs **Ea ≈ 17-34 kJ/mol**. With the Ea 64 kJ/mol used for TC in the cluster-5 transfer, the
   pentose ring would release **~9 × 10⁻³ min⁻¹ at 100 °C, 12-46x above the band**. If Karpinski's pH-6 phosphate
   inflated the rate by Kallen 1971b's saturating x3.5-5, the needed Ea rises to ~35-52 kJ/mol. The verdict's own
   escape ("phi ~0.1 or Ea near 40") loses its phi branch: what remains is a low Ea (~20-50 kJ/mol) for release.

## 2. Conditions

- **Compounds** [M]: Thz, PrThz, AThz from L-cysteine + aldehyde (ref 23); AThz printed as
  "2-L-arabino-tetrahydroxybutyl" on p. 130 (the abstract's OCR is garbled).
- **Solutions** [M]: stock 1-5 mM in deaerated 0.1 M Na2HPO4 or 0.05 M CH3COOH, kept under argon; kinetic runs
  mostly at 10⁻⁵-10⁻⁴ M (PrThz 10⁻⁵-10⁻³ M); "phosphate buffer" named for pH 6.0, 6.7, 7.4 runs (figure captions);
  the pH 4.5 / 5.2 buffer is not named (acetate likely, from the stock) [!]; ionic strength not stated; pH a set
  buffer value, drift not reported.
- **Temperature** [M]: water-jacketed cell, 298 +- 0.2 K.
- **O2** [M]: deaerated, argon throughout; air-saturated arms flagged in Table 1.
- **Analysis** [M]: NP voltammetry (50 ms pulses) and Osteryoung square-wave (10 Hz, 25 mV) at an SMDE (0.0167 cm²);
  calibration per pH for each Thz and cysteine; ln c vs t linear (r > 0.99 for PrThz, pH 4.5-7.4); NP and OSW k_a
  agree within 5 %. Thz / AThz runs: ~8 h, ~10 % conversion, RSD up to 30 %.

## 3. Quantitative results

**Table 1** (first-order decomposition constants, 25 °C, 10⁵ k_a in s⁻¹) [M]:

| compound | pH | 10⁵ k_a / s⁻¹ | k_a, min⁻¹ [Z] |
|---|---|---:|---:|
| Thz | 7.4 | 0.06 | 3.6e-5 |
| Thz | 6.0 | 0.3 | 1.8e-4 |
| Thz | 6.0 + O2 (air-saturated) | 1.3 | 7.8e-4 |
| PrThz | 8.0 | 0.03 | 1.8e-5 |
| PrThz | 7.4 | 1.6 | 9.6e-4 |
| PrThz | 7.4 + O2 | 2.3 | 1.4e-3 |
| PrThz | 6.7 | 4.3 | 2.6e-3 |
| PrThz | 6.0 | 7.3 | 4.4e-3 |
| PrThz | 6.0 + O2 | 14.8 | 8.9e-3 |
| PrThz | 5.2 | 4.7 | 2.8e-3 |
| PrThz | 4.5 | 2.1 | 1.3e-3 |
| AThz | 7.4 | 0.08 | 4.8e-5 |
| AThz | 6.0 | 0.3 | 1.8e-4 |

**Fig. 8** (PrThz, pH 6.0) [F]: deaerated ln c from ~ -10.0 to ~ -11.6 over ~3-21 × 10³ s (slope ~ -9 × 10⁻⁵ s⁻¹,
consistent with 7.3 × 10⁻⁵ within the read error); with air, steeper. Fig. 4 (CV, pH 6.0, 5 × 10⁻⁴ M PrThz): "all
PrThz had decomposed" in 6 h [M].

**Text** [M]: decomposition "much slower at pH ca. 3 and ca. 8 than at pH from 5 to 7"; with O2, Thz at pH 6.0 and
PrThz at pH 7.4 are linear only for the first 3 h, then accelerate.

## 4. Mechanism evidence and its kind

Measured product (cysteine identified electrochemically; mass balance ring + Cys = c0 within 5 %); absence of a
detectable open-chain thiol (an upper bound on its concentration, not a rate); within-study pH and O2 factors.
Kind: measured intermediate (negative) + mass balance + within-study factor. The "ring opening is rate-determining"
conclusion is [P]: a low steady-state open fraction is equally expected when ring-chain opening is fast and
ring-favoured and the hydrolysis is slow (Fife 1991, `fife1991_extraction.md`). Laboratory: Warsaw (NEW).

## 5. Bearing on open questions

- **CF-076, the provisional admission (E2c / E4c family; P 0.45)** — **CONTRADICTS one of its two escape routes and
  narrows the other.** Phi ~0.1 is not supported: a pentose-cysteine ring releases like TC at 25 °C (phi ~1). The
  admission then requires release Ea ~17-34 kJ/mol at face value (to ~52 with a phosphate-catalysis allowance), not
  the 64 kJ/mol used for TC. Not a rejection: Ea has never been measured for any thiazolidine on disk.
- **CF-076, the kallen1971 correction (k_r 120-270x slower than the relaxation)** — **SUPPORTS** it. Karpinski's TC
  (Thz) at pH 6.0, 25 °C is 1.8 × 10⁻⁴ min⁻¹, x3.3 the corrected Kallen value (5.4 × 10⁻⁵ min⁻¹ at pH 5.6; plateau to
  6.6), and 50-170x below the old relaxation reading (1.5-5 × 10⁻⁴ s⁻¹). Independent lab, independent method.
- **CF-076, k_r(pH) shape** — SUPPORTS a bell with its maximum near pH 6 for PrThz and a fall above 6.6 (Thz x5, AThz
  x3.75 from 6.0 to 7.4), as in Kallen's shape. It does not give the acid limb for the pentose ring.
- **CF-076 remedy (1), the TTCA ring-N pKa** — LEAVES UNCHANGED (not measured).
- **F04 unknown "K(pH) and K(T) for any cysteine-sugar thiazolidine"** — LEAVES UNCHANGED for K; adds the first
  measured RELEASE rate of a cysteine-pentose ring (two pHs, 25 °C).
- **F04 unknown "inert vs air for the thiazolidine itself"** — **FILLED at 25 °C** (not at process temperature): air
  speeds ring loss x1.4-4.3 and yields cystine. Bears on the cysteine-sink rival in Wei 2022 / Zhai 2021 (headspace O2
  oxidises thiol): O2 acts on the released cysteine and possibly on an open-chain intermediate.
- **Atlas sec. 1.4 "O2 vs inert at >= 100 C"** — LEAVES UNCHANGED (25 °C).

## 6. Disagreements found

1. **CF-076 declared band vs Karpinski-transferred release (NEW).** Quantity: k_r(100 °C, pH 4.5) for a pentose-Cys
   thiazolidine. Band 2-8 × 10⁻⁴ min⁻¹; Karpinski AThz (25 °C, pH-shape and Ea 64 transferred) ~9 × 10⁻³ min⁻¹ [Z],
   x12-46. Candidate causes: release Ea much lower than 64 kJ/mol (chemistry); AThz's acid limb steeper than PrThz's
   (the pH 4.5 value is a transfer); phosphate catalysis in Karpinski's pH 6 buffer (conditions); L-arabino vs
   D-xylo configuration (chemistry, probably small); the Zhai dilute pots not in the regime the band assumes.
2. **Position of the release maximum among 2-substituted rings (NEW).** PrThz (2-alkyl): maximum near pH 6.0;
   k(4.5) / k(6.0) = 0.29. Bustillo Trueba 2018 (`bustillotrueba2018_extraction.md`; 2-aryl, 2-furyl and
   2-(methylthio)ethyl, 20 °C): fastest at pH 4.4-5.2, x1.6-4.7 SLOWER at pH 6.0. Candidate causes: substituent
   (alkyl vs aryl / heteroaryl; ring pKa); buffer (phosphate / acetate here vs Carmody borate-citrate-phosphate);
   method (voltammetric cysteine vs UPLC-UV ring loss); concentration (10⁻⁵-10⁻⁴ M vs 0.5 mM); 25 vs 20 °C.
3. **PrThz pH 7.4 -> 8.0 drops x53** (1.6 -> 0.03) [!] — steeper than one ionisation would give; the authors note the
   literature conflicts above pH 7 (Pesek vs Riemschneider) but do not discuss their own step.

## 7. Ledger-ready measurements

- Table 1 AThz (pH 6.0, 7.4) and Thz (pH 6.0, 7.4) — the pentose-ring release rate and its TC reference from one lab;
  two readers; carry the 30 % RSD.
- Table 1 PrThz pH ladder (4.5-8.0, deaerated) — a 2-alkyl release pH profile.
- Table 1 air-saturated arms (Thz 6.0, PrThz 6.0 and 7.4) — the only O2 contrast on thiazolidine stability on disk.
