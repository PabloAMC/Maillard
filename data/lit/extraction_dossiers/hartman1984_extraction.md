# Hartman, Carlin, Scheide & Ho 1984 — EXTRACTION
### *"Volatile Products Formed from the Thermal Degradation of Thiamin at High and Low Moisture Levels"*
### Qualitative. Very concentrated thiamine HCl (40 g in 95 g water or propylene glycol, water pH 2.3) at 135 °C for 30 min: 2-methyl-3-furanthiol is "medium" in water and absent in PG, its disulfide "extra large" in PG and "small" in water. Size classes only (XS–XL); no concentrations, no time course, one temperature, no H2S. Read 2026-10-28.

**Full read of `data/articles/hartman1984.pdf`** (4 pages, text layer). Table I was checked on a
220-dpi raster of p. 1016; the text layer misaligns the PG/water columns in several rows (e.g. it
puts 2,5-dimethyl-3-furanthiol's "XS" in the wrong column). **[!] The PDF also carries the end of
the preceding article (p. 1015, top half: CAS registry list that includes hydrogen sulfide, carbon
disulfide and dimethyl disulfide) and the start of the following one (p. 1018, Saura-Calixto et al.,
almond sugars).** Those registry numbers do not belong to this paper; H2S is not in its Table I.

**Provenance:** **[M]** measured, printed · **[Q]** qualitative · **[Z]** derived by me, never
printed · **[C]** outside the paper · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Guy J. Hartman, James T. Carlin, Jürgen D. Scheide, Chi-Tang Ho* |
| Title | "Volatile Products Formed from the Thermal Degradation of Thiamin at High and Low Moisture Levels" |
| Venue | *J. Agric. Food Chem.* **1984**, 32 (5), **1015–1018** |
| DOI | **10.1021/jf00125a016**. Not printed in the PDF. Crossref REST API (`api.crossref.org/works/10.1021/jf00125a016`, 2026-10-28): same title, JAFC 32(5):1015–1018, Sep 1984, authors Hartman, Carlin, Scheide, Ho. Matches |
| Affiliation | Dept. of Food Science, Cook College, NJ Agricultural Experiment Station, Rutgers, New Brunswick NJ |
| Dates | received 2 Mar 1984, accepted 25 Jun 1984; NJAES Publication No. D-10205-1-84 |
| PDF character | Table I (36 compounds × 2 solvents, size classes), Fig. 1 (mechanism scheme for 3-methyl-4-oxo-1,2-dithiane). No numbers beyond the recipe; no replicates stated |

The paper's "Hartman et al., 1984" citations are the companion *Lebensm.-Wiss. Technol.* paper
(thiamine/cystine/ascorbate/MSG, high vs low moisture, in press); the a_w series in the same system
is `hartman1984b_extraction.md`.

---

## 1. THE ONE-PARAGRAPH ANSWER

**40 g thiamine HCl + 95 g water or propylene glycol (PG)**, sealed in a 0.3-L stainless cylinder,
**135 °C for 30 min** with stirring; **water solution pH 2.3** (PG pH not given). Volatiles by
vacuum steam distillation, ether extraction, GC-MS. Table I lists 36 compounds with a size class
per solvent (XS, S, M, L, XL; footnote defines them, but **no scale or detector response is
given**). Water gives more compounds and more of them: thiophenes, pentanones, most thiazoles and
all furanthiols appear **only in water**. **2-Methyl-3-furanthiol: M in water, absent in PG.
bis(2-Methyl-3-furyl) disulfide: XL in PG, S in water.** The authors note "it is not known why"
the disulfide is favoured at low moisture. 2,5-Dimethyl-3-furanthiol (XS) and the PG acetals
(dioxolanes, L–XL) are PG-only. H2S is discussed from the literature (Dwivedi & Arnold 1971) as
the likely partner for thiophene formation but **was not measured**. The authors' summary: water
plays a "hydrolytic as well as amphoteric role".

---

## 2. METHODS [M] (pp. 1015–1016)

| item | as printed |
|---|---|
| Thiamine | reagent-grade thiamin hydrochloride (Sigma), **40 g** [Z: 0.119 mol] |
| Solvent | **95 g water** or **95 g propylene glycol** [Z: ~1.25 mol thiamine per kg water, ~1 M by volume; ~10⁴× the 1.2 × 10⁻⁴ M of Dwivedi & Arnold 1971] |
| pH | "The pH of the water solution of thiamin was **2.3**" (unbuffered, thiamine HCl's own acidity); PG not given |
| Vessel | **0.3-L Hoke SS-DOT sample cylinder**, loaded at atmospheric pressure and sealed (air headspace, not stated) |
| Heating | **135 °C for 30 min**, magnetic stirring; heating medium and come-up not stated [!] |
| Isolation | diluted to 95 g PG + 695 g water; vacuum steam distillation at 0.3 mmHg, flask at 4 °C, four-stage cold-finger traps (dry ice, dry ice/acetone); NaCl-saturated, extracted with 1 L ether, back-washed 3×, Na₂SO₄, concentrated to 1 mL (spinning band) |
| Analysis | HP-5985 GC/MS, DB-1 30 m (1 µm film), 40 °C/2 min → 250 °C at 4 °C/min, 1 µL splitless |
| Quantitation | **none**: "S = small; XS = extra small; M = medium; L = large; XL = extra large" [Q]. No internal standard, no replicate |

---

## 3. RESULTS — Table I as printed (p. 1016, checked on the raster)

"Volatile Products of Thiamin in Water or Propylene Glycol". Blank = not reported in that solvent.

| class / compound | PG | water | previously identified |
|---|---|---|---|
| **carbonyls** | | | |
| acetic acid | | XL | |
| ethyl acetate | M | M | |
| n-propyl acetate | | S | |
| acetol acetate | S | | |
| 2-pentanone | | S | |
| 3-pentanone | | S | |
| 2,3-pentanedione | | S | |
| 5-chloro-2,3-pentanedione | S | | |
| 5-chloro-2-pentanone | S | | |
| **furans** | | | |
| 2-methyl-4,5-dihydrofuran | | S | van der Linde et al. (1979) |
| 2-methyltetrahydrofuran-3-one | | XL | van der Linde et al. (1979) |
| 2-methyl-3(2H)-furanone | | M | van der Linde et al. (1979) |
| 2-methyltetrahydro-3-furanthiol | | S | |
| **2-methyl-3-furanthiol** | | **M** | van der Linde et al. (1979) |
| 2-methyl-4,5-dihydro-3-furanthiol | | L | van der Linde et al. (1979) |
| 2,5-dimethyl-3-furanthiol | XS | | |
| **bis(2-methyl-3-furyl) disulfide** | **XL** | **S** | |
| 2-methyl-3-furyl 2-methyl-4,5-dihydro-3-furyl disulfide | S | S | |
| bis(2-methyl-4,5-dihydro-3-furyl) disulfide | | S | van Dort et al. (1984) |
| 2-methyl-4,5-dihydro-3-furyl 2-methyltetrahydro-3-furyl disulfide | | S | |
| 2,5-dimethyl-4,5-dihydro-3-furyl 2-methyl-4,5-dihydro-3-furyl trisulfide | | S | |
| **thiophenes** | | | |
| 2-methylthiophene | | M | van der Linde et al. (1979) |
| 2-methyl-4,5-dihydrothiophene | | L | Dwivedi and Arnold (1973) |
| 2-methyltetrahydrothiophene | | M | |
| 2-methyltetrahydrothiophen-3-one | | M | Dwivedi and Arnold (1973); van der Linde et al. (1979) |
| 2-acetylthiophene | | L | |
| **thiazoles** | | | |
| 4-methylthiazole | M | L | |
| 4,5-dimethylthiazole | | S | Dwivedi and Arnold (1973); van der Linde et al. (1979) |
| 5-ethyl-4-methylthiazole | | S | van der Linde et al. (1979) |
| 5-(2-hydroxyethyl)-4-methylthiazole | | M | Dwivedi and Arnold (1973); van der Linde et al. (1979) |
| 5-(2-chloroethyl)-4-methylthiazole | XS | S | |
| **dioxolanes** | | | |
| 4-methyl-1,3-dioxolane | L | | |
| 2,4-dimethyl-1,3-dioxolane | L | | |
| 2-ethyl-4-methyl-1,3-dioxolane | XL | | |
| **miscellaneous sulfur compounds** | | | |
| dimethyl disulfide | | S | |
| methyl thiol acetate | | S | |
| 3-mercapto-2-pentanone | | M | van der Linde et al. (1979) |
| 3-methyl-4-oxo-1,2-dithiane | | S | Hartman et al. (1984) |
| 1-methylbicyclo[3.3.0]-2,8-dioxa-4-thiaoctane | | M | van Dort et al. (1984) |
| 1,3-dimethylbicyclo[3.3.0]-2,8-dioxa-4-thiaoctane | | L | van Dort et al. (1984) |

The MS-reference column and the printed partial mass spectra (for compounds without a library
reference) are omitted here; they identify, they do not quantify.

Text [Q]: carbonyls, furans and thiazoles "greater in number and quantity" in water; thiophenes
"only identified in the water system", possibly because they need H2S; pentanones only in water,
consistent with water being needed to form 3-mercapto-5-hydroxy-2-pentanone (the model's HMP); the
chlorinated pentanones (PG only) "most likely unreacted synthetic intermediates of thiamin";
3-methyl-4-oxo-1,2-dithiane minor here, major in the companion system with cystine and ascorbic
acid (proposed route via 3-mercapto-5-hydroxy-2-pentanone, Fig. 1).

---

## 4. DERIVED QUANTITIES [Z]

None that carry units. Two ordinal readings, both at 135 °C, 30 min, pH 2.3 (water):
- **Monomer/dimer ordering flips with solvent:** water MFT (M) > MFT-MFT (S); PG MFT (absent) <
  MFT-MFT (XL). With Guth et al. 1995 (`guth1995_extraction.md`: disulfides hydrolyse back to thiols in
  hot water) this is consistent with disulfide cleavage being water-dependent, but the paper offers
  no test.
- **pH 2.3 still yields MFT, its dihydro analogue (L) and the HMP-derived pentanones.** Dwivedi &
  Arnold 1972 (`dwivedi1972_extraction.md`) bound all non-thiazole sulfur to ≤ 4 % of degraded
  thiamine at pH 3.5–6, 121 °C. No contradiction is testable here (no amounts); at ~1 M thiamine a
  sub-percent channel still gives GC-visible peaks.

---

## 5. WHAT THE MODEL CAN TAKE FROM IT

**Context.** Thiamine (THI) feeds HMP (`r_thi_hmp`, → MFT) and MeSH (`r_thi_mesh`) and, from B50,
H2S through `thiamine_h2s_share(pH)` (zero below pH 6.25). The sulfur lane has no water-activity
or solvent term.

**Verdict: (c) qualitative only.** No rate, no yield, no concentration, one temperature.
- **G5 (thiamine degradation with H2S at ≥ 2 temperatures): NO.** One temperature, H2S not measured.
- **G1 (thiol loss): NO.**
- **Qualitative checks it supports:** (i) MFT forms from thiamine alone at **pH 2.3, 135 °C** in water
  (the HMP route is open at low pH; whether `k_thi_hmp` carries the right pH shape is not testable
  from size classes); (ii) the **low-moisture direction pushes MFT into its disulfide** (or prevents
  the disulfide's hydrolysis); a model calibrated in water would under-predict MFT-MFT and
  over-predict free MFT at low moisture; (iii) the thiophenes require water (proposed H2S
  involvement).
- **Contradictions with the model:** none testable.

**Not usable:** the size classes as numbers (no scale); the H2S, CS2, DMDS registry numbers on p. 1015
(another paper's).

**Caveats.** ~1 M thiamine in a sealed cylinder (self-buffered at pH 2.3) is far from any food;
no replicate; no blank isolate reported; extraction at 4 °C under vacuum loses H2S and other
highly volatile species.
