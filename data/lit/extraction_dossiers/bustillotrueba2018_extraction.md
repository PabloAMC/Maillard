# Bustillo Trueba, De Clippeleer, Van der Eycken, Guevara Romero, De Rouck, Aerts & De Cooman 2018 — EXTRACTION
DOI: 10.1080/03610470.2018.1546094 (printed on the cover page and p. 272; Crossref 2026-11-01: "Influence of pH on the Stability of 2-Substituted 1,3-Thiazolidine-4-Carboxylic Acids in Model Solutions", J. Am. Soc. Brew. Chem. 76(4), 272-280, 2018-10-02. **Matches the brief.**)
### Three cysteine-aldehyde thiazolidines — of furfural (FURF-CYS), benzaldehyde (BENZ-CYS) and methional (METH-CYS) — at 100 mg/L (~0.5 mM; ring loss by UPLC-PDA) and 2 mg/L (~10 uM; free aldehyde by HS-SPME on-fibre PFBHA GC-MS), Carmody borate / citrate / phosphate buffers at pH 2.0 / 4.4 / 5.2 / 6.0 / 9.0, prepared under N2, **20 °C only**, 0-24 h, n = 3. **No sugar ring, no temperature series, no Ea, no K, no measured ring pKa** (5-6 is cited from Butvin 1999).

**Read 2026-11-01 from `data/articles/bustillotrueba2018.pdf`** (10 pp.: T&F cover page + pp. 272-280; born-digital).
Prose from `pdftotext -layout`; **Fig. 4 (% remaining / % released) and Fig. 6 (mass-balance bars) read from 170 dpi
renders of pp. 277-278** (scratchpad `ing5/bt-07.png`, `bt-08.png`); the remaining-% values below come from the Fig. 6
black bars (cleaner than the Fig. 4 lines). Read error +-3 percentage points.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** KU Leuven, Technology Campus Ghent, Laboratory of Enzyme, Fermentation and Brewing Technology (EFBT;
Aerts, De Cooman), with Ghent University / University College Ghent brewing science (De Clippeleer) and KU Leuven
LOMAC (Van der Eycken). **Not in `lab_ids.yml`: NEW** (proposed id `kuleuven_ghent_efbt`). [!] The existing
text_aliases "Ghent" (-> `ghent_de_kimpe`) and "Leuven" (-> `leuven_hendrickx`) would both mis-assign this group
if it is ever named bare in the cited text.

## 1. What this paper answers for the model

| model need (CF-076) | what Bustillo Trueba 2018 gives | tag |
|---|---|---|
| ring loss / aldehyde release vs pH | five pH values, three 2-substituted rings, 0-24 h at 20 °C (Figs. 4, 6) | [F] |
| release = ring loss? | yes early: mass balance (ring % + aldehyde %) "about 100%" [V] at 0-1 h; drifts at 10-24 h | [F] / [Q] |
| a rate constant | **none printed**; initial apparent first-order constants derivable from the 0-5 h bars [Z] | [Z] |
| shape vs pH | **bell, maximum at pH 4.4-5.2**: slower at 6.0 (x1.6-4.7), slow at 2.0, nearly stable at 9.0 | [F] / [Z] |
| T series / Ea | **none** (20 °C) | — |
| K (formation equilibrium) | **none** | — |
| ring-N pKa | **cited only**: "pKa only 5 up to 6" depending on the C-2 substituent (Butvin et al. 1999) | [C] |
| a sugar (polyhydroxyalkyl) ring | **none** | — |

**Answers.**
1. **Stability vs pH (20 °C, initial 0-5 h apparent k, h⁻¹)** [Z from Fig. 6]: BENZ-CYS 0.20 (pH 4.4), 0.17 (5.2),
   0.063 (6.0), 0.037 (2.0), < 0.005 (9.0); FURF-CYS 0.11, 0.074, 0.023, 0.035, < 0.005; METH-CYS 0.077, 0.074, 0.047,
   0.012, < 0.006. In min⁻¹ at pH 4.4: 3.3 × 10⁻³ (BENZ), 1.8 × 10⁻³ (FURF), 1.3 × 10⁻³ (METH). Order BENZ > FURF >
   METH, as the authors predict from the C-2 substituent [Q].
2. **No T series, no Ea** (confidence ~99 %).
3. **Does it bound CF-076's k_r / K(pH)? Not directly** (confidence ~85 %). No sugar ring, one temperature, no K. What
   it does say:
   - the release maximum of 2-substituted Cys rings sits at pH ~4.4-5.2, about one unit more acid than Kallen's TC
     (rising 4.5 -> 5.6, flat to 6.6). That is the direction CF-076's admission needs (TTCA ring pKa <= 5.9), but for
     aryl / alkyl substituents, which `karpinski1990_extraction.md` shows release 20-24x faster than a pentose ring;
   - 2-substituted rings release 100-240x faster than TC at room temperature (1.3-3.3 × 10⁻³ min⁻¹ at 20 °C, pH 4.4, vs
     TC 1.4 × 10⁻⁵ min⁻¹ at 25 °C, pH 4.5, cluster-5 verdict sec. 3) [Z]. Not a bound on TTCA: Karpinski's pentose ring
     releases like TC, not like these.
   - the acid limb is real: pH 2.0 is x3-6 slower than pH 4.4. The authors' explanation (S-protonation blocking ring
     opening) is "personal communication" [P].

## 2. Conditions

- **Adducts** [M]: synthesised from aldehyde (11 mmol in MeOH) + L-cysteine (10 mmol in water), RT, 1-2 h, closed
  flask; precipitate washed with ether, dried; 1H-NMR in DMSO-d6 (4 isomers: C-2, C-4). Eight made; only FURF-CYS,
  METH-CYS, BENZ-CYS detected well enough to follow.
- **Buffers** [M]: Carmody universal series (boric acid, citric acid, trisodium phosphate; 0.2 / 0.05 / 0.1 M
  component stocks; buffer strength 0.025-0.03 mol/L per pH unit), pH 2.0, 4.4, 5.2, 6.0, 9.0; label pH, not
  re-measured. N2-sparged 20 min; samples made in an anaerobic workstation.
- **Concentrations** [M]: 100 mg/L for UPLC (FURF-CYS 0.50, BENZ-CYS 0.48, METH-CYS 0.49 mM [Z]); a 50-fold dilution
  (2 mg/L, ~10 uM) in 20 mL headspace vials for the aldehyde. [!] Ring loss and aldehyde release are measured at
  two concentrations 50x apart, then summed into one "mass balance".
- **Time / T** [M]: separate vials per time point, 0, 1, 5, 10, 24 h at 20 °C.
- **Analysis** [M]: UPLC-PDA at 227 nm (HSS C18, 45 °C, water / MeOH 55 : 45 + acetic acid, pH 5.5); HS-SPME with
  on-fibre PFBHA, GC-MS, deuterated internal standards, a calibration per pH (R 0.9955-1.000).

## 3. Quantitative results

**Remaining ring, % of t = 0** (Fig. 6 black bars; pH 2.0 / 4.4 / 5.2 / 6.0 / 9.0) [F]:

| adduct | 1 h | 5 h | 10 h | 24 h |
|---|---|---|---|---|
| FURF-CYS | 94 / 88 / 87 / 95 / 98 | 84 / 58 / 69 / 89 / 98 | 75 / 40 / 49 / 85 / 97 | 53 / 16 / 23 / 68 / 97 |
| BENZ-CYS | 94 / 75 / 75 / 91 / 99 | 83 / 37 / 43 / 73 / 98 | 72 / 19 / 21 / 54 / 97 | 42 / 6 / 3 / 21 / 93 |
| METH-CYS | 97 / 84 / 83 / 89 / 98 | 94 / 68 / 69 / 79 / 97 | 91 / 60 / 64 / 79 / 97 | 85 / 40 / 60 / 62 / 90 |

**Released aldehyde at 24 h, % of full release** (Fig. 4B) [F]: furfural 29 / ~98 (SD very large) / 67 / 20 / 3;
benzaldehyde 11 / 76 / 86 / 33 / 3; methional 5 / 42 / 38 / 43 / ~2.

**Mass balance** (Fig. 6) [F]: ~95-114 % at 0-1 h; at 24 h 53-115 % (BENZ pH 2.0 53 %, pH 6.0 54 %; FURF pH 4.4
114 %).

## 4. Mechanism evidence and its kind

Within-study pH factor with a two-sided measurement (ring loss and aldehyde gain), anaerobic. Kind: within-study
factor + partial mass balance. Mechanism (Fig. 5: N-protonation, C-S ring opening to the iminium, hydration,
release) is [P], drawn not shown. The METH-CYS plateau at ~60 % (pH 5.2-6.0, 10-24 h) may be an approach to
equilibrium; but the 10 uM aldehyde vials release about the same fraction (~40 %) as the 0.5 mM ring vials lose,
which an equilibrium at a single K would not give (dilution favours dissociation) [Z] — so no K is derivable.
Laboratory: KU Leuven Ghent EFBT (NEW).

## 5. Bearing on open questions

- **CF-076 (k_r(pH) shape)** — SUPPORTS a bell-shaped release with an acid limb and a high-pH fall; its maximum
  (pH 4.4-5.2) is more acid than Kallen's TC and Karpinski's PrThz (pH ~6): weak support for a lower ring-N pKa in
  2-substituted rings (the admission's pKa <= 5.9), for aryl / alkyl substituents only.
- **CF-076 (k_r magnitude, phi)** — LEAVES UNCHANGED: no sugar ring; Karpinski 1990 is the relevant phi datum.
- **CF-076 remedy (1), ring pKa** — LEAVES UNCHANGED; a lead: Butvin, Al-Ja'Afreh, Svetlik & Havranek 1999,
  "Solubility, Stability, and Dissociation Constants of (2RS,4R)-2-Substituted Thiazolidine-4-Carboxylic Acids in
  Aqueous Solutions", Chem. Pap. 53:315-322 [C] (no Crossref DOI for the journal article; only a ChemInform abstract,
  10.1002/chin.200017131) — the measured 2-substituted ring pKa values.
- **CF-011 (phosphate)** — LEAVES UNCHANGED: phosphate is in every Carmody buffer; not varied.
- **F04 unknown "competition of other carbonyls for cysteine"** — LEAVES UNCHANGED (20 °C stability only; furfural,
  benzaldehyde, methional rings are aroma-relevant reservoirs but no process-temperature data).
- **F26 / Strecker aldehyde release** — adds that Cys-bound Strecker aldehydes (methional, benzaldehyde) release at
  room temperature on hours at pH 4.4-6.0: a possible storage-operator effect on aldehyde measurements in Cys pots.

## 6. Disagreements found

1. **Release pH maximum vs Karpinski 1990 (NEW).** Bustillo: fastest at pH 4.4-5.2, k(4.4)/k(6.0) = 1.6-4.7 (20 °C).
   Karpinski (PrThz, 25 °C): maximum near pH 6, k(4.5)/k(6.0) = 0.29. Candidate causes: C-2 substituent (aryl,
   heteroaryl, thioether-alkyl vs n-propyl) shifting the ring pKa; buffer (Carmody borate-citrate-phosphate vs
   phosphate / acetate); method (UPLC ring loss vs voltammetric cysteine gain); concentration (0.5 mM vs 10-100 uM);
   20 vs 25 °C (cannot move a pH maximum).
2. **Mass balance across two concentrations** [!]: the paper's "mass balance" adds % ring lost at 0.5 mM to % aldehyde
   released at ~10 uM. For a reversible dissociation these need not agree; the near-100 % early balance says the
   systems are far from equilibrium at 0-1 h, not that the method is validated.

## 7. Ledger-ready measurements

- Fig. 6 / Fig. 4A ring-remaining ladders (3 adducts x 5 pH x 4 times) — a 20 °C pH profile of 2-substituted Cys-ring
  stability under N2; key as % remaining (two readers), with the concentration noted.
- Fig. 4B aldehyde release at 24 h (pH 4.4-6.0) — for the methional / benzaldehyde storage-release check.
- Not keyable: any K, Ea or ring pKa (none measured).
