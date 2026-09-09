# Brands & van Boekel 2001 — EXTRACTION (monosaccharide–casein reaction network at 120 C, pH 6.8 phosphate)

**Source on disk:** `data/articles/brands2001.pdf` (9 pp., born-digital ACS PDF; text layer clean). Read
method: `pdftotext -layout` for all body text, plus 300-dpi rasters of pp. 3-5 for Figures 1-6, which are
the **only** place any number lives — the paper has **no tables, no rate constants and no activation
energies**. Every concentration below is digitised by this extraction from the figures (±3 % of the axis
span, i.e. ±4 mmol/L on the sugar panels, ±0.2 mmol/L on the acid panels, ±0.4 mmol/L on the lysine
panels, ±0.1 A420 on the browning panels). Read-only extraction, 2026-09-07.

## 0. Identity

| field | value |
|---|---|
| Title | "Reactions of Monosaccharides during Heating of Sugar-Casein Systems: Building of a Reaction Network Model" |
| Authors | Carline M. J. Brands and Martinus A. J. S. van Boekel (Wageningen University) |
| Venue | J. Agric. Food Chem. 2001, 49 (10), 4667-4675 |
| DOI as printed | `10.1021/jf001430b` (p. 4667 footer, "10.1021/jf001430b CCC: $20.00"); manuscript id JF001430B |
| Dates | received Nov 28 2000; revised Jun 12 2001; accepted Jun 15 2001; web Aug 31 2001 |
| Companion | the kinetic analysis promised here is Brands & van Boekel 2002, JAFC 50, 6725-6739 (not on disk as such; `brands2002b.pdf` is the melanoidin extinction-coefficient paper) |

## 1. Why it matters to the model

This is the same laboratory, the same buffer (0.1 M phosphate, pH 6.8), the same tubes and the same
120 C as the Martins & van Boekel 2005 glucose/glycine backbone, but with **casein lysine** as the amine
and with **fructose, galactose and tagatose** run alongside glucose. It is therefore a direct test of the
trunk's sugar-side steps (glucose <-> fructose isomerisation, glucose -> formic + acetic acid, Amadori
build-up and breakdown) at 120 C in water, with and without an amine, plus the first same-lab
fructose-side dataset. The paper establishes the reaction network (Figure 10) that the 2002 companion fits;
it fits nothing itself.

## 2. Methods as they matter to a model

- **Charges:** sodium caseinate 3 % w/w (spray-dried, 90 % protein) + 150 mmol/L monosaccharide (glucose,
  fructose, galactose or tagatose) in 0.1 M phosphate buffer pH 6.8; "molar ratio of sugar to lysine
  residues of 10:1", i.e. ~15 mmol/L lysine residues (Fig. 1/2 right panels start at 15.0 / 14.5 mmol/L).
- **Vessel:** screw-capped glass tubes (Schott, 16 x 160 mm); volume per tube not stated. Oil bath 120 C.
- **Times:** 0, 2, 5, 10, 15, 20, 30, 40 min; "reported heating times include the heating-up period of
  ~2-3 min". Ice-water quench. "Heat-treated and analyzed in at least duplicate."
- **pH:** initial 6.7 (abstract) / buffer 6.8 (methods); NOT controlled during heating; fell 0.3 (glucose)
  and 0.4 (fructose) units in 40 min (Fig. 3).
- **Sugar-only controls:** glucose and fructose "heated in the absence of protein (remaining conditions
  kept unchanged)" (Fig. 6).
- **Isolated Amadori:** 150 mM glucose + 3 % casein incubated 65 C, 15 h; glycated protein separated on
  Sephadex G25; then heated at 120 C (Fig. 7; figure not digitised here, text: acetic ~1.5x formic).
- **Analytics:** sugars and organic acids after Sephadex G25 desalting, HPLC ION-300 (0.0025 M H2SO4, 85 C),
  RI for sugars, UV 210 nm for acids; total acid by titration to pH 8.3 with 0.1 N NaOH; available lysine
  by OPA fluorescence; Amadori as furosine (8 M HCl, 110 C, 23 h) x 3.1; Heyns via CML after periodate
  (not detected); HMF, furfuryl alcohol, HHMF, DDMP by RP-HPLC (UV 280 / 220 nm); methylglyoxal as
  OPD quinoxaline (qualitative only); browning A420 after 4x dilution in 16 % SDS; protein-bound
  melanoidins from A420 with epsilon = 500 L/(mol cm) "(unpublished results)" (later 477 ± 50, Brands 2002b).
- **Quantification basis:** mmol/L of solution; melanoidins as "sugar units incorporated"; mass balance
  as % of initial sugar.

## 3. Data (all digitised from figures; no printed tables exist)

**Figure 1** — "Glucose-casein solutions heated at 120 C: glucose; fructose; formic acid; acetic acid;
lysine residues; Amadori compound." (mmol/L)

| t (min) | glucose | fructose | formic | acetic | lysine res. | Amadori |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 150 | 0 | 0 | 0 | 15.0 | 0 |
| 2 | 150 | 0 | 0 | 0 | 15.0 | 0.1 |
| 5 | 139 | 5 | 0.35 | 0.1 | 13.8 | 0.6 |
| 10 | 128 | 13 | 0.75 | 0.75 | 12.7 | 1.0 |
| 15 | 121 | 19 | 0.95 | 1.55 | 11.0 | 1.0 |
| 20 | 110 | 24 | 1.1 | 2.15 | 10.5 | 0.8 |
| 30 | 102 | 30 | 2.8 | 3.5 | 8.7 | 0.6 |
| 40 | 93 | 35 | 3.3 | 5.3 | 8.7 | 0.5 |

**Figure 2** — "Fructose-casein solutions heated at 120 C" (same legend). Amadori = 0 at every time;
Heyns not detected.

| t (min) | fructose | glucose | formic | acetic | lysine res. |
|---:|---:|---:|---:|---:|---:|
| 0 | 146 | 0 | 0 | 0 | 14.5 |
| 2 | 142 | 0 | 0.45 | 0 | — |
| 5 | 141 | 0 | 1.2 | 0 | 13.7 |
| 10 | 133 | 4 | 2.0 | 0.45 | 12.2 |
| 15 | 121 | 9 | 2.5 | 1.4 | 11.1 |
| 20 | 108 | 10 | 2.9 | 2.15 | 9.7 |
| 30 | 102 | 15 | 4.1 | 4.1 | 8.7 |
| 40 | 97 | 17 | 5.0 | 5.85 | 7.3 |

**Figure 3** — "pH and total amount of acids as found by titration and HPLC in heated glucose-casein (A)
and fructose-casein (B) systems." (acids in mmol/L)

| t (min) | A: pH | A: titration | A: HPLC | B: pH | B: titration | B: HPLC |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 6.72 | 0 | 0 | 6.70 | 0 | 0 |
| 5 | 6.70 | 1.0 | 0.5 | 6.65 | 2.8 | 1.1 |
| 10 | 6.67 | 3.0 | 1.5 | 6.61 | 5.5 | 2.5 |
| 15 | 6.61 | 7.2 | 2.5 | 6.55 | 9.8 | 3.9 |
| 20 | 6.58 | 9.2 | 3.2 | 6.49 | 13.2 | 5.1 |
| 30 | 6.53 | 14.2 | 6.3 | 6.38 | 21.0 | 8.3 |
| 40 | 6.46 | 17.2 | 8.5 | 6.30 | 25.2 | 10.8 |

Titrated acid is ~2x the formic + acetic found by HPLC in both systems: "other organic acids were formed
but were not identified" (lactic, glycolic, saccharinic acids proposed).

**Figure 4** — "Browning of total system, protein fraction, and sugar fraction expressed in absorbance
units measured at 420 nm and concentration of protein-bound melanoidins in heated glucose-casein (A) and
fructose-casein (B) systems." (A420 after 4x SDS dilution; melanoidin mmol/L = A420(protein)/500 x 1000)

| t (min) | A: total | A: protein | A: sugar fr. | A: melanoidin mM | B: total | B: protein | B: sugar fr. | B: melanoidin mM |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0.2 | 0.2 | 0 | 0.4 | 0.3 | 0.3 | 0 | 0.6 |
| 15 | 0.75 | 0.6 | 0.1 | 1.2 | 1.15 | 0.9 | 0.05 | 1.8 |
| 20 | 1.35 | 1.1 | 0.25 | 2.2 | 1.75 | 1.15 | 0.6 | 2.3 |
| 30 | 2.3 | 1.95 | 0.35 | 3.9 | 3.5 | 2.8 | 0.65 | 5.6 |
| 40 | 3.3 | 2.65 | 0.7 | 5.3 | 4.75 | 3.65 | 1.05 | 7.3 |

**Figure 5** (mass balance, % of initial sugar, 40 min): glucose-casein: glucose ~62, fructose ~23, Amadori
~1, total acids ~10, brown ~3, missing ~1; fructose-casein: fructose ~66, glucose ~12, acids ~17, brown ~4.
Text: "almost negligible amount of missing compounds after 40 min"; "between 10 and 30 min more compounds
were missing".

**Figure 6** — "Glucose solutions (top) and fructose solutions (bottom) heated without casein at 120 C:
glucose; fructose; formic acid; acetic acid; absorbance at 420 nm."

| t (min) | Glc-only: glucose | fructose | formic | acetic | A420 | Fru-only: fructose | glucose | formic | acetic | A420 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 152 | 0 | 0 | 0 | 0 | 147 | 0 | 0 | 0 | 0 |
| 2 | 153 | 0 | 0 | 0 | 0 | 149 | 0 | 0.1 | 0 | 0 |
| 5 | 142 | 6 | 0.15 | 0 | 0 | 140 | 0 | 0.9 | 0 | 0.05 |
| 10 | 135 | 15 | 0.4 | 0 | 0 | 135 | 3 | 3.2 | 0.5 | 0.15 |
| 15 | 129 | 22 | 0.55 | 0 | 0 | 126 | 10 | — | 1.75 | 0.55 |
| 20 | 112 | 30 | 1.05 | 0.2 | 0.05 | 115 | 11 | 5.45 | 3.0 | 0.85 |
| 30 | 105 | 38 | 1.7 | 1.1 | 0.15 | 112 | 21 | 6.3 | 5.25 | 1.65 |
| 40 | 93 | 45 | 2.4 | 2.05 | 0.5 | 99 | 17 | 6.55 | 6.55 | 2.45 |

**Other quantities in the text:** HMF and furfuryl alcohol "0-40 uM"; HHMF and DDMP identified by
spectrum only; methylglyoxal detected as OPD quinoxaline (no number); isolated protein-bound Amadori
heated at 120 C gave formic and acetic acid with acetic "~1.5 times higher", and "no sugars were formed";
mannose and psicose not detected; galactose/tagatose systems "more rapid" than glucose/fructose (no numbers).

**Reaction scheme proposed (Figure 10, "Reaction network model for sugar-casein reactions"):** aldose <->
ketose (via 1,2-enediol anion); aldose -> Cn (unidentified sugar fragments) + acids; ketose -> Cn + acids;
aldose + lysine-R -> Amadori -> acids (formic via 3-deoxyaldoketose C1-C2 cleavage; acetic via
1-deoxy-2,3-diketose C2-C3 cleavage or via triose/methylglyoxal) + AMP; ketose + lysine-R -> AMP directly
(no Heyns detected); Cn + lysine-R -> AMP; AMP -> melanoidins. Figures 8-9 give the chemistry
(beta-elimination, alpha-dicarbonyl cleavage, retro-aldolisation; Amadori via 1,2-enaminol, 2,3-enolisation
to 1-deoxy-2,3-diketose). No step is assigned a rate constant in this paper.

## 4. What the repo could take

Nothing here is a printed rate; the FIT-eligible items are **within-study ratios and rates derived from
digitised time courses**, all at one temperature (120 C), so no activation energy can come from this paper.

**Derived first-order loss rates at 120 C (this extraction, ln(c0/c40)/40 min; heating-up 2-3 min included
in t, so true isothermal rates are ~5-8 % higher):**

| system | k_loss (/min) | comparison |
|---|---:|---|
| glucose + casein | 0.012 | Martins 2005 trunk at 120 C: k_glc_fru 1.6e-3 x exp(Ea 122.6 kJ/mol, 100->120 C) = 0.012 /min for the isomerisation step alone; plus k_schiff x [amine]. Same order. |
| glucose alone | 0.012 | identical to with-casein: casein does not change net glucose loss at 120 C over 40 min |
| fructose + casein | 0.010 | |
| fructose alone | 0.010 | |

**Initial glucose -> fructose rate (0-10 min):** 13 mmol/L per 10 min from 150 mM = 0.0087 /min (with casein);
15 mM per 10 min = 0.010 /min (alone). Trunk `k_glc_fru` extrapolated to 120 C = 0.012 /min — agreement within
1.4x. Candidate FIT row: `k_glc_fru @ 120 C, pH 6.8 phosphate, no amine = 0.010 /min (±20 %, digitised)`.

**Fructose -> glucose:** in fructose-only, glucose reaches 21 mM at 30 min from 147 mM: initial rate
~3-10 mM per 10 min, i.e. 0.002-0.007 /min. Trunk `k_fru_glc` (9.2e-3 at 100 C, Ea 93.4) extrapolates to
0.043 /min at 120 C, which would convert ~80 % of a fructose-only charge in 40 min; Brands' fructose-only
run loses only 33 % of fructose in total. **The trunk's step-3 constant is a glucose-system fit that
does not transfer to a fructose-rich charge** — a validation claim, not a fit row.

**Within-study ratios (all 40 min, 120 C):**

| ratio | value | reading |
|---|---:|---|
| acetic / formic, glucose + casein | 5.3 / 3.3 = 1.6 | trunk (Martins, 100 C, 4 h, glycine) has ~5; casein at 120 C gives less acetic-dominance |
| acetic / formic, glucose alone | 2.05 / 2.4 = 0.85 | without amine, formic >= acetic |
| acetic, glucose+casein / glucose alone | 5.3 / 2.05 = 2.6 | the amine (via Amadori 2,3-enolisation) is the main acetic-acid route from glucose |
| formic, glucose+casein / glucose alone | 3.3 / 2.4 = 1.4 | formic acid is mostly sugar-only chemistry |
| acetic / formic, fructose + casein | 5.85 / 5.0 = 1.2 | |
| acetic / formic, fructose alone | 6.55 / 6.55 = 1.0 | |
| total acid (fructose alone) / (glucose alone) | 13.1 / 4.45 = 2.9 | ketose degrades to acids ~3x faster than aldose; text: "ketoses seemed to be more reactive in the sugar degradation reactions" |
| A420 fructose alone / glucose alone | 2.45 / 0.5 = 4.9 | caramelisation browning matters for fructose, not glucose (text cites Buera 1987) |
| melanoidin (fructose-casein) / (glucose-casein) | 7.3 / 5.3 = 1.4 | |
| lysine lost, glucose-casein / fructose-casein | 6.3 / 7.2 mM (42 % / 50 %) | "about equal or somewhat higher" for ketose |
| Amadori peak, glucose-casein | 1.0 mM at 10-15 min, falling to 0.5 at 40 min | 7 % of initial lysine; Amadori never accounts for lysine loss |
| titrated acid / HPLC (formic+acetic) | 2.0 (glucose), 2.3 (fructose) | half the acid is unidentified — the trunk's FA + AA underestimate total acid by ~2x |

**Directional claims:** (i) acetic-acid formation from glucose has a lag (~10 min at 120 C without amine,
~5 min with casein) while formic has none — the 1-deoxy-2,3-diketose intermediate is real; (ii) casein does
not accelerate net glucose loss but does redirect it toward acetic acid; (iii) no Heyns product from
fructose-casein, and fructose-casein acid formation equals fructose-alone acid formation ("acid formation via
the Maillard reaction is apparently not significant" for ketoses); (iv) glucose alone: fructose is 45 mM at
40 min vs 35 mM with casein — the amine competes for the glucose.

## 5. Caveats

- **No table, no fitted constant, no Ea.** All numbers above are my digitisation of 300-dpi rasters; use
  them as ±3 %-of-axis, and as VALIDATION unless re-digitised at higher resolution.
- Single temperature (120 C); heating-up 2-3 min inside the reported times; pH uncontrolled (drifts to
  6.3-6.5 by 40 min).
- Matrix is casein (protein-bound lysine, 15 mM), not a free amino acid; the melanoidin axis uses epsilon =
  500 L/(mol cm), a value the same group later measured as 477 ± 50.
- Amadori is furosine x 3.1 (a milk-derived conversion factor); Heyns not detectable by the method.
- Mass balance is only "almost" closed at 40 min and has a dip at 10-30 min; half of the titratable acid is
  unidentified.
- Ajandouz 2008 Table 3 credits this paper with "Ea 121 kJ/mol, 110-150 C" for glucose-casein browning;
  no such number exists here (it is Brands 2002 or Morales 1998).
- Figure 7 (isolated Amadori time course) was not digitised; only the text ratio (acetic ~1.5x formic) is
  recorded.
