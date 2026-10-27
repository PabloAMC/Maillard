# Goldberg, Kishore & Lennen 2002 — EXTRACTION
### *"Thermodynamic Quantities for the Ionization Reactions of Buffers"*
### The evaluated pK, ΔrH° and ΔrCp° at 298.15 K, I = 0 for 64 buffers. Transcribed: acetate, citrate, L-cysteine, diethanolamine, glycine, glycylglycine, phosphate, phthalate, pyrophosphate, Tricine. Not in GKL: formic acid, asparagine, lysine, sulfide/H2S. Includes an audit of every `[GKL]` citation in `ph_state.py` / `parameters_sulfur.py`. Read 2026-10-27.

**Full read of the front matter (Sections 1–4), Section 8 (summary table) and Tables 7.2, 7.22,
7.23, 7.24, 7.31, 7.33, 7.51, 7.52, 7.56, 7.65 of `data/articles/Goldberg2002.pdf`** (141-page
PDF with a text layer, extracted with `pdftotext -layout`). Each per-table "Selected values" line was
checked against the Section 8 summary (pp. 359–360); they agree for every entry below.

**Provenance:** **[M]** printed value · **[Z]** derived by me from printed values, never printed ·
**[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Robert N. Goldberg, Nand Kishore, Rebecca M. Lennen |
| Title | "Thermodynamic Quantities for the Ionization Reactions of Buffers" |
| Venue | *J. Phys. Chem. Ref. Data* **31** (2), **231–370** (2002) |
| DOI | 10.1063/1.1416902 |
| Affiliation | Biotechnology Division, NIST, Gaithersburg MD |
| Dates | received 21 June 2001, accepted 16 Sept 2001, published 24 April 2002 |
| Character | Evaluated compilation (no new measurements). Per buffer: selected values, a table of literature values, those values adjusted to 298.15 K and I = 0, and a comments paragraph. |

---

## 1. WHAT IS TABULATED AND HOW TO USE IT

- **Quantities [M]:** pK, ΔrG°, ΔrH°, ΔrCp° for each **macroscopic** ionisation at T = 298.15 K,
  p = 0.1 MPa. Standard state: hypothetical ideal solution of unit **molality**, I = 0. The pK
  values are pK_m. The authors state they "limited the scope to the macroscopic ionization
  constants" (p. 233). No microscopic constants are given (relevant to cysteine's SH/NH3+ pair).
- **Temperature dependence (Eq. 12, Clarke & Glew), R = 8.31451 J K⁻¹ mol⁻¹:**
  R ln K = −ΔrG°/Tr + ΔrH°(1/Tr − 1/T) + ΔrCp°[(Tr/T) − 1 + ln(T/Tr)] + (Tr/2)(∂ΔrCp°/∂T)p[...].
  **No (∂ΔrCp°/∂T) value is selected for any buffer** (p. 234), so in practice the constant-ΔrCp°
  truncation is the GKL prescription. In pK form:
  pK(T) = pK(Tr) + ΔrH°/(R ln10)·(1/T − 1/Tr) − ΔrCp°/(R ln10)·[ln(T/Tr) + Tr/T − 1].
- **Intended range [M]:** "the temperature range that is generally encountered in the vast
  majority of biochemical studies" (p. 232). Using it at 120–180 °C is an extrapolation GKL do not
  endorse. **[!]**
- **Ionic strength [M]:** GKL adjust to I = 0 with extended Debye–Hückel, ln γ = −A z² √I/(1 + B√I),
  **B = 1.6 kg^½ mol^−½** (p. 234). Not Davies. They warn the I = 0 extrapolation is the main
  uncertainty, and that adjustments from I > 0.1 M can carry errors "substantially larger" than
  ±0.05–0.10 in pK.
- **Evaluation codes [M] (p. 235)**, one letter per quantity in the order pK, ΔrH°, ΔrCp°:

| code | pK | ΔrH° (kJ/mol) | ΔrCp° (J K⁻¹ mol⁻¹) |
|---|---|---|---|
| A | ±0.003 | ±0.2 | ±10 |
| B | ±0.03 | ±0.5 | ±20 |
| C | ±0.1 | ±2.0 | ±40 |
| D | ±0.3 | ±5.0 | ±80 |
| U | larger than D, or not assessable | | |

A two-letter code (e.g. "CD") means ΔrCp° was not determined. "≈" before a value is GKL's.

---

## 2. SELECTED VALUES AT 298.15 K AND I = 0 [M]

Transcribed from each table's "Selected values" line; units pK (–), ΔrG° and ΔrH° kJ mol⁻¹,
ΔrCp° J K⁻¹ mol⁻¹.

| GKL table (page) | buffer | reaction | pK | ΔrG° | ΔrH° | ΔrCp° | evaluation |
|---|---|---|---|---|---|---|---|
| 7.2 (p. 238) | acetate | CH3COOH = H⁺ + CH3COO⁻ | 4.756 | 27.147 | −0.41 | −142 | AAA |
| 7.22 (p. 265) | citrate | (1) H3L = H⁺ + H2L⁻ | 3.128 | 17.855 | 4.07 | −131 | AAB |
| | | (2) H2L⁻ = H⁺ + HL²⁻ | 4.761 | 27.176 | 2.23 | −178 | AAB |
| | | (3) HL²⁻ = H⁺ + L³⁻ | 6.396 | 36.509 | −3.38 | −254 | AAB |
| 7.23 (p. 268) | L-cysteine | (1) H3L⁺ = H⁺ + H2L | 1.71 | 9.76 | ≈ −0.6 | — | CD |
| | | (2) H2L = H⁺ + HL⁻ | 8.36 | 47.72 | 36.1 | ≈ −66 | BBD |
| | | (3) HL⁻ = H⁺ + L²⁻ | 10.75 | 61.36 | 34.1 | ≈ −204 | BBD |
| 7.24 (p. 270) | diethanolamine | HL⁺ = H⁺ + L | 8.883 | 50.705 | 42.08 | 36 | AAA |
| 7.31 (p. 281) | glycine | (1) H2L⁺ = H⁺ + HL± | 2.351 | 13.420 | 4.00 | −139 | AAA |
| | | (2) HL± = H⁺ + L⁻ | 9.780 | 55.825 | 44.2 | −57 | AAB |
| 7.33 (p. 287) | glycylglycine | (1) H2L⁺ = H⁺ + HL± | 3.140 | 17.923 | 0.11 | −128 | AAB |
| | | (2) HL± = H⁺ + L⁻ | 8.265 | 47.177 | 43.4 | −16 | AAB |
| 7.51 (p. 321) | phosphate | (1) H3PO4 = H⁺ + H2PO4⁻ | 2.148 | 12.261 | −8.0 | −141 | AAA |
| | | (2) H2PO4⁻ = H⁺ + HPO4²⁻ | 7.198 | 41.087 | 3.6 | −230 | AAA |
| | | (3) HPO4²⁻ = H⁺ + PO4³⁻ | 12.35 | 70.49 | 16.0 | −242 | BCA |
| 7.52 (p. 326) | phthalate | (1) H2L = H⁺ + HL⁻ | 2.950 | 16.838 | −2.70 | −91 | AAA |
| | | (2) HL⁻ = H⁺ + L²⁻ | 5.408 | 30.869 | −2.17 | −295 | AAA |
| 7.56 (p. 333) | pyrophosphate | (1) H4P2O7 = H⁺ + H3P2O7⁻ | 0.83 | 4.74 | −9.2 | ≈ −90 | CDD |
| | | (2) H3P2O7⁻ = H⁺ + H2P2O7²⁻ | 2.26 | 12.90 | −5.0 | ≈ −130 | BCD |
| | | (3) H2P2O7²⁻ = H⁺ + HP2O7³⁻ | 6.72 | 38.36 | 0.5 | −136 | BBC |
| | | (4) HP2O7³⁻ = H⁺ + P2O7⁴⁻ | 9.46 | 54.00 | 1.4 | −141 | DCC |
| 7.65 (p. 351) | Tricine | (1) H2L⁺ = H⁺ + HL± | 2.023 | 11.547 | 5.85 | −196 | AAB |
| | | (2) HL± = H⁺ + L⁻ | 8.135 | 46.435 | 31.37 | −53 | AAA |

**Not in GKL [M, by absence from the 68-entry table of contents and a full-text search]:**
formic acid/formate, L-asparagine, L-lysine, sulfide/hydrogen sulfide, thiazolidine-4-carboxylic
acid. Other GKL entries that touch this model: ammonia (Table 7.7: pK 9.245, ΔrH° 51.95, ΔrCp° 8,
AAA), carbonate (Table 7.20), L-histidine (7.39), 2-mercaptoethanol (7.43: pK 9.75, ΔrH° 26.2, no
ΔrCp°).

### 2.1 L-cysteine, the literature behind Table 7.23 [M]

GKL adopt pK1 from Borsook et al. 1937 alone; pK2 and pK3 are averages of Borsook 1937, Wrathall
et al. 1964 and Coates et al. 1969; ΔrH° for (2) and (3) are **Wrathall 1964 calorimetry**; the
ΔrCp° values (≈ −66, ≈ −204) are **GKL's own fit to Coates 1969's pK(T)**. Coates et al. measured
pK at I = 0 from 278.15 to 368.15 K, which allows an in-paper check of the Clarke–Glew
extrapolation:

| T (K) | Coates pK2 [M] | Coates pK3 [M] | GKL Eq. 12, (2) [Z] | GKL Eq. 12, (3) [Z] | constant-ΔH 26 / 43 from 8.33 / 10.28 (model B9) [Z] |
|---|---|---|---|---|---|
| 278.15 | 8.84 | 11.20 | 8.82 | 11.21 | 8.66 / 10.82 |
| 298.15 | 8.37 | 10.70 | 8.36 | 10.75 | 8.33 / 10.28 |
| 338.15 | 7.69 | 10.08 | 7.64 | 10.12 | 7.79 / 9.39 |
| 368.15 | 7.28 | 9.76 | 7.23 | 9.84 | 7.46 / 8.85 |

GKL's selected set reproduces Coates to within 0.05 (pK2) and 0.08 (pK3) at 95 °C. The model's B9
pair is off by +0.18 (pK2) and **−0.91 (pK3)** at 95 °C. Comments [M]: cysteine purity can bias pK
by ≈0.04, and cysteine is light-sensitive. Ionic-strength 3.0 M data (GKL code 72GRA/WIL) give pK2 8.78 and
pK3 10.71, i.e. **I matters by ~0.4 for pK2**.

### 2.2 Phosphate comments [M]

ΔrH°(1) = −8.0 is the mean of all studies except Bjerrum & Unmack 1929; the calorimetric mean is
−8.2. ΔrH°(2) = 3.6 (all-study mean 3.69, calorimetric mean 3.57). **[!] Typo in GKL:** the
comments paragraph writes ΔrCp°(1) "= 141"; the selected-values line and Section 8 both print −141.

---

## 3. AUDIT OF THE MODEL'S CONSTANTS AGAINST GKL

Code read 2026-10-27 on branch `b51-ph-state-everywhere` (uncommitted working tree), **not edited**.
`ph_state.TitratableGroup` stores (pKa, ΔH kJ/mol, ΔCp J K⁻¹ mol⁻¹) per step. "Cited" is what the
code's source string attributes to `[GKL]`. Values cited to `[CRC]` are checked too, where GKL has
them.

| model constant | value in code | cited to | GKL value (table) | verdict |
|---|---|---|---|---|
| `FORMIC` ΔH | −0.13 | [GKL] | **formic acid not in GKL** | **MISATTRIBUTED**: the number cannot come from GKL |
| `FORMIC` ΔCp | −142 (acetic, flagged as a carried class value) | GKL 7.2 by analogy | acetate −142 (7.2) | MATCH as a carried value; the carry is disclosed |
| `ACETIC` pK / ΔH / ΔCp | 4.756 [CRC] / −0.41 / −142 | GKL Table 7.2 | 4.756 / −0.41 / −142 (7.2, AAA) | MATCH ×3 |
| `CYSTEINE_B9` pK | 1.92 / 8.33 / 10.28 | [CRC] | 1.71 / 8.36 / 10.75 (7.23) | not cited to GKL; pK3 differs by **0.47**. GKL's literature has no 10.28 (Borsook 10.78, Wrathall 10.76, Coates 10.70) |
| `CYSTEINE_B9` ΔH | 0 / 26 / 43 | [GKL] | ≈ −0.6 / **36.1 / 34.1** (7.23) | **MISMATCH** (steps 2 and 3). The code's own string already says so |
| `CYSTEINE_GKL` pK | 1.710 / 8.36 / 10.75 | GKL Table 7.23 | 1.71 / 8.36 / 10.75 | MATCH ×3 |
| `CYSTEINE_GKL` ΔH | −0.6 / 36.1 / 34.1 | GKL Table 7.23 | ≈ −0.6 / 36.1 / 34.1 | MATCH ×3 |
| `CYSTEINE_GKL` ΔCp | −139 / −66 / −204 | GKL 7.23; step 1 carried from glycine | — (not printed) / ≈ −66 / ≈ −204 | MATCH ×2. Step 1 is correctly disclosed as glycine's (7.31) −139 |
| `TTCA_GROUP` ΔH | 0 / 40 | [GKL, amine class] | **GKL prints no class value.** Nearest: glycine NH3⁺ 44.2, glycylglycine 43.4, diethanolamine (secondary) 42.08, Tricine (N-substituted glycine, secondary) 31.37 | **NOT PRINTED**: 40 is an inference, not a GKL number |
| `TTCA_GROUP` ΔCp | −139 / −57 | GKL Table 7.31 (as class values) | glycine −139 / −57 (7.31) | MATCH to Table 7.31. **[!]** The "class" is the code's assumption: diethanolamine's secondary ammonium is **+36** and Tricine's −53 |
| `PHOSPHATE` pK | 2.148 / 7.198 / 12.35 | [CRC] | 2.148 / 7.198 / 12.35 (7.51) | MATCH (also GKL) |
| `PHOSPHATE` ΔH | **−7.9** / 3.6 / 16.0 | [GKL] | **−8.0** / 3.6 / 16.0 (7.51) | **MISMATCH step 1** (0.1 kJ/mol; inside grade A ±0.2; moves pK1 at 145 °C by 0.005). Steps 2 and 3 match |
| `PHOSPHATE` ΔCp | −141 / −230 / −242 | GKL Table 7.51 | −141 / −230 / −242 | MATCH ×3 |
| `CITRATE` pK | 3.128 / 4.761 / 6.396 | [CRC] | 3.128 / 4.761 / 6.396 (7.22) | MATCH (also GKL) |
| `CITRATE` ΔH | 4.07 / 2.23 / −3.38 | [GKL, MS] | 4.07 / 2.23 / −3.38 (7.22) | MATCH ×3 |
| `CITRATE` ΔCp | −131 / −178 / −254 | GKL Table 7.22 | −131 / −178 / −254 | MATCH ×3 |
| `PYROPHOSPHATE` pK / ΔH / ΔCp | 0.83/2.26/6.72/9.46; −9.2/−5.0/0.5/1.4; −90/−130/−136/−141 | GKL Table 7.56 | identical (7.56) | MATCH ×12. **[!]** grades CDD/BCD/BBC/**DCC**: pK4 is ±0.3 |
| `PHTHALATE` pK / ΔH / ΔCp | 2.950/5.408; −2.70/−2.17; −91/−295 | GKL Table 7.52 | identical (7.52, AAA/AAA) | MATCH ×6 |
| `PhDrift.arp_group` carboxyl | pK 3.00, ΔH 0 | [CRC] class value | not in GKL | n/a |
| `PhDrift.arp_group` ammonium ΔH | **45.0** | [GKL, secondary-amine class] | **GKL prints no class value.** The closest structural analogue in GKL (an N-substituted glycine, secondary ammonium) is **Tricine, 31.37** (7.65); diethanolamine 42.08 (7.24) | **NOT PRINTED**: 45 is not a GKL number and exceeds every secondary-amine entry transcribed here |
| `PhDrift.arp_group` ΔCp | −139 / −57 | "glycine's class values" | glycine −139 / −57 (7.31) | MATCH to Table 7.31 as a carry; same caveat as TTCA |
| `parameters_sulfur.PKA_25C["thiol"]`, `["cysteine_amine"]` | 8.33 / 10.28 | "standard aqueous physical chemistry" | 8.36 / 10.75 (7.23) | not cited to GKL here, but `CYSTEINE_B9` cites the same pair to GKL. pK3 differs by 0.47 |
| `parameters_sulfur.DELTA_H_IONISATION_KJ_MOL["thiol"]`, `["cysteine_amine"]` | 26.0 / 43.0 | same as above (and `CYSTEINE_B9` → [GKL]) | **36.1 / 34.1** (7.23) | **MISMATCH**: the B2 pair. Off by −10.1 and +8.9 kJ/mol |
| `parameters_sulfur.PKA_GKL_CYSTEINE` | thiol (8.36, 36.1, −66); amine (10.75, 34.1, −204) | GKL Table 7.23 | identical | MATCH ×6 |
| `parameters_sulfur` `h2s_1` | 7.05 / 22.0 | not GKL | **no sulfide in GKL** | n/a. See `suleimenov1997_extraction.md`: pK1 6.99 (s.v.p.), ΔH 22.7 kJ/mol, ΔCp −265 |

Table numbers cited in the code: **7.2, 7.22, 7.23, 7.31, 7.51, 7.52, 7.56**. Each is the right
table for the buffer named. The pre-registration (`results/validation/kinetic_core_b51_prereg.md`)
cites the same seven.

**Formula check.** `ph_state.pk_at` implements GKL Eq. 12 with constant ΔrCp°, the sign
conventions correct (R in kJ, ΔCp converted from J). [Z] With GKL's Table 7.51 values it gives
phosphate pK2 = **7.633 at 145 °C** (7.017 with constant ΔH), the figure the code comments quote.
R differs from GKL's 8.31451 in the sixth figure, which is negligible.

**Summary of discrepancies:** (i) FORMIC ΔH −0.13 cited to GKL, but formic acid is not in GKL;
(ii) PHOSPHATE ΔH1 −7.9 vs GKL −8.0; (iii) the B9 / `parameters_sulfur` cysteine ΔH pair 26/43 vs
GKL 36.1/34.1 (already flagged in the code and `tasks/lessons.md`); (iv) TTCA 40 kJ/mol and ARP
45 kJ/mol cited to a GKL "amine class" that GKL does not print; (v) the ARP/TTCA ΔCp "class values"
are glycine's numbers, correctly transcribed, but GKL's secondary amines span −53 to +36.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor:**
- pK, ΔH and ΔCp at 25 °C, I = 0, for every buffer the panels declare (phosphate, citrate,
  acetate, phthalate, pyrophosphate) and for cysteine. These are evaluated equilibrium constants
  with stated uncertainties.
- The **form** of the temperature extrapolation (Eq. 12 with constant ΔCp). For cysteine it is
  validated in-paper to 95 °C (§2.1).
- A grade-based uncertainty band per constant. Pyrophosphate pK4 (D, ±0.3) and cysteine pK1 (C, ±0.1)
  should carry wider bands than phosphate pK2 (A, ±0.003).

**Cannot anchor:**
- **Anything above ~95 °C.** GKL target the biochemical range. The highest-T data behind any
  transcribed value are Coates 1969's cysteine to 368 K. At 120–180 °C the omitted (∂ΔCp/∂T) term
  is unknown, and Suleimenov 1997 shows ΔCp of H2S ionisation goes from −265 at 25 °C to −321 at
  150 °C. A constant-ΔCp error is plausibly a few tenths of a pK unit at 145 °C for high-|ΔCp|
  steps (phthalate 2, citrate 3, phosphate 2/3). [Z, order of magnitude only]
- **High ionic strength.** Values are at I = 0. GKL use extended Debye–Hückel (B = 1.6) only to
  extrapolate *to* I = 0; nothing here validates Davies at I ≈ 0.5–1.5 M. Cysteine at 3 M
  (72GRA/WIL) shows pK2 shifted +0.4.
- **Microscopic constants.** Cysteine steps 2/3 are macroscopic. Mapping step 2 to the thiol and
  step 3 to the ammonium, as `CYSTEINE_GKL` does, is an approximation outside GKL's scope.
- **Formic acid, asparagine, lysine, sulfide, Amadori compounds, thiazolidines:** absent. The B54+
  "measured asparagine ionisation enthalpy" need is **not** met by GKL.
- **Class values.** GKL never assigns a ΔH or ΔCp to a functional-group class. Any "[GKL, X class]"
  citation is an inference and should be labelled as such.

**Caveats.** All numbers were read from the PDF's text layer and cross-checked between the
per-table line and Section 8. The literature tables (§2.1 Coates data) were read from the same
layer. The [Z] Clarke–Glew values were computed in the `maillard` env with the formula above, not by
running the model.
