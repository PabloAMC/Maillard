# Buttery, Ling & Mon 1986 — EXTRACTION (a quantitative method for 2-acetyl-1-pyrroline in rice: 200 g rice, 2 h Likens-Nickerson with an acid solvent phase, collidine internal standard, measured relative recovery)

### The paper that puts a number on what Buttery 1983 could not measure: the relative recovery of 2-acetyl-1-pyrroline through steam distillation is 28.0 % (SD 4.9 %), so the 1983 levels are multiplied by 3.57 — and four re-measured rice levels, in ppb.

**Source on disk:** `data/articles/buttery1986.pdf` (J. Agric. Food Chem. 1986, 34 (1), 112-114; the
PDF's first page carries the tail of the preceding article). The `pdftotext -layout` text layer is
clean for prose and lost only the fraction inside the recovery sentence and the body of Table I,
both of which were recovered from 300-dpi rasters of printed pages 113 and 114
(`scratchpad/img/bu86eq3-2.png`, `bu86t1b-3.png`). Figures 1, 2 and 3 (a capillary chromatogram of
basmati steam volatile oil; found-versus-added recovery plot; two basic-fraction chromatograms) are
**FIGURE-ONLY** and were not read for values.

## 0. Identity

| field | value |
|---|---|
| Title | "Quantitative Analysis of 2-Acetyl-1-pyrroline in Rice" |
| Authors | Ron G. Buttery, Louisa C. Ling and Thomas R. Mon (Western Regional Research Center, ARS, USDA, Albany CA) |
| Venue | J. Agric. Food Chem. 1986, 34 (1), 112-114. Received for review June 10 1985; accepted October 30 1985. "This article not subject to U.S. Copyright" |
| Registry number printed | 2-acetyl-1-pyrroline **99583-29-6** — **different from the 85213-22-5 printed in Buttery 1983 and stored in `data/species/desirable_targets.yml`**; see flag 5 |
| Naming | "AP" = 2-acetyl-1-pyrroline; "TMP" = collidine (2,4,6-trimethylpyridine), the internal standard. "ppb = parts per billion (10^9)", footnote a of Table I |
| Predecessor | Buttery, Ling, Juliano & Turnbaugh 1983, JAFC 31:823 (`buttery1983_extraction.md`), on disk — the method here is built on its isolation step and explicitly corrects its numbers |
| Repo registry | `2_acetyl_1_pyrroline` in `data/keys/compounds.yml`; target row in `data/species/desirable_targets.yml` |

## 1. Why it matters

Two things, both about how a 2-acetyl-1-pyrroline number should be read.

**It fixes the level scale.** `buttery1983_extraction.md` Table I is the repository's only
2-acetyl-1-pyrroline level table and its authors called it "meant only to give a general idea".
This paper measures the one thing that was missing — the compound's recovery through the isolation,
**28.0 % relative to the internal standard, SD 4.9 %** — and re-runs three of the same varieties with
it applied. Any level benchmark the engine is ever scored against should come from this table, not
from the 1983 one.

**It is the repository's only evidence on the compound's stability, and it is the reason the engine's
missing loss step matters.** The engine has no decay channel for 2-acetyl-1-pyrroline. This paper
prints, without a rate: that the compound is unstable; that its recovery loss "seems to occur during
the steam distillation step" and "may be partly due to the instability"; that GC injector
temperatures above 150-170 C "seemed to cause some decomposition"; and — the fact that matters most
for a kinetic model — that "**what the method actually measures is not the concentration in the rice
at any particular time but the total amount of 2-acetyl-1-pyrroline produced from the given weight
of rice during the heating (2-h) period**". A measured level of this compound is a time-integrated
production, not a state, in every source the repository has.

What it is not: it is an analytical-method paper. It contains no rate constant, no barrier, no
temperature series and no Maillard model system.

## 2. Methods as they matter to a model

- **Standards.** Authentic 2-acetyl-1-pyrroline synthesised as in Buttery et al. 1983; a standard
  solution at **290 ppm in benzene**, stored at -20 C, "seemed quite stable over several months".
  Collidine (2,4,6-trimethylpyridine, Eastman No. 4815) as a **30.0 ppm water solution**, "quite
  stable at room temperature". Volatile-free water made by boiling off about 10 % of the volume;
  antifoam made by concentrating 20 mL of GE AF 60 silicone antifoam emulsion in 600 mL water down
  to 200 mL. All glassware oven-heated at 120 C for several hours.
- **Isolation.** 6 L water in a 12 L round-bottom flask with a large magnetic stirrer on a 1300 W
  two-circuit 115 V mantle; **200 g of rice** added gradually while stirring; 50 mL antifoam and
  **5.00 mL of the 30 ppm collidine solution (= 150 µg of internal standard)** added; a
  Likens-Nickerson head (Kontes K-523010-0000) attached; the solvent flask charged with **80 mL
  water + 2.0 mL concentrated sulfuric acid and 120 mL diethyl ether**, stirred. Mantle 100 V upper
  / 80 V lower to boiling, then 80 V / 70 V; solvent mantle 80 V. **2 h.** Explicit warning: "Care
  should be taken to prevent any 'burning' of the rice on the bottom of the 12-L flask as this
  produces interfering compounds (**alkylpyrazines**)."
- **The method's one idea.** The acid sits in the solvent flask, so the basic volatiles are stripped
  out of the steam distillate and locked into the aqueous acid layer continuously, while the
  non-basic bulk of the rice volatiles stays in the ether: "A convenient and efficient way to extract
  the basic fraction is to use both dilute acid and diethyl ether in the solvent flask in the steam
  distillation continuous-extraction apparatus." Figure 1's chromatogram of the whole volatile oil
  is what this avoids — 2-acetyl-1-pyrroline is peak 62 among a very large number of components.
- **Recovery of the basic fraction.** Aqueous acid layer separated, covered with 120 mL ether,
  15 g NaHCO3 added in portions over about 5 min, shaken, ether layer separated, dried over 40 g
  anhydrous Na2SO4, filtered through a washed cotton plug, concentrated through a 17 cm Vigreux over
  a ~70 C bath to about 1 mL, transferred to a graduated 2.5 mL micro test tube and taken to about
  **0.05 mL** over a ~50 C bath. **One ether extraction only**, deliberately: "To minimize
  manipulation time, only one extraction with ether was used. Extracting several consecutive times
  would undoubtably improve the transfer of material but would add considerably to the analysis
  manipulation time. A recovery factor relative to the internal standard seemed to be sufficient to
  compensate."
- **Choice of internal standard, stated as criteria.** Collidine was chosen because it must (1) have
  related properties to 2-acetyl-1-pyrroline — basic, similar water solubility, similar volatility;
  (2) be stable; (3) elute close to it; (4) be commercially available. Measured
  **ether/water partition coefficients for 2-acetyl-1-pyrroline and collidine were in the ratio
  0.74 / 1.0**, "the collidine being slightly less soluble in water".
- **Gas chromatography.** FID. Main column: laboratory-made **150 m x 0.5 mm i.d. Pyrex glass
  capillary, wall coated with Carbowax 20M**, 15 psi He, 50 C for 30 min after injection then
  1 C/min to 170 C; 2 µL of the ether concentrate injected into a **170 C injector** connected
  directly to the column, no split. Alternative commercial capillary: J&W 60 m x 0.326 mm i.d.,
  0.25 µm Durabond wax, 50 to 230 C at 4 C/min, 22 psi, 0.2 µL split 1/50 — retention times
  2-acetyl-1-pyrroline **26 min**, collidine **29 min**. Alternative packed column: 1.2 m x 3 mm o.d.
  Pyrex, 80-100 mesh Chromosorb G with 10 % Carbowax 20M, isothermal 90 C, 6 psi He, 2 µL, 170 C
  injector — retention times **22 min** and **26 min**.
- **The calculation, printed in full.**

      concn (ppb) = (area of AP peak / area of TMP peak) x 150 x 3.57 x 5

  with "The factor 150 is the number of micrograms of internal standard added, and the factor 5 is
  needed to convert to parts per billion (ppb)", and 3.57 the reciprocal of the recovery.
  **The arithmetic checks (mine):** 5.00 mL x 30.0 µg/mL = 150 µg of collidine; 200 g of rice is
  0.2 kg, so µg per 200 g x 5 = µg/kg = ppb; and 100 / 28.0 = 3.571.
- **The recovery factor itself.** "The relative recovery factor was determined by adding known
  quantities of both 2-acetyl-1-pyrroline and collidine to 6 L of water in the 12-L flask and
  carrying them through the whole process. This factor for the conditions described in the
  Experimental Section was determined to be **28.0 %** (with a standard deviation of **4.9 %** for
  **nine GLC determinations and three steam distillation continuous-extraction procedures**). It is
  necessary to multiply the results based on the internal standard then by the factor 100/28.0 which
  equals **3.57**." Note the recovery was measured **in water, without rice** — see flag 3.
- **Method test.** A common American rice containing **less than 10 ppb** 2-acetyl-1-pyrroline was
  spiked with known amounts and carried through; found is plotted against added in Figure 2 with a
  found = added line. The paper's verdict is prose only: "there is some scatter around the line, but
  the method seems accurate enough for the purpose, considering the small concentrations involved
  and the unstable nature of 2-acetyl-1-pyrroline." **No regression, slope, intercept, r2 or number
  of spike levels is printed**, and Figure 2's axes run to about 3000 ppb. The figure is FIGURE-ONLY
  and no point is read from it here.

## 3. Tables re-typed

### Table I. "Concentration of 2-Acetyl-1-pyrroline Found in Some Samples of Cooked Rice"

Footnote a as printed: "ppb = parts per billion (10^9)."

| rice sample | concn, ppb |
|---|---:|
| Malagkit Sungsong (brown) | 760 |
| basmati 370 (brown) | 610 |
| IR 841-76-1 (brown) | 560 |
| Texas long grain (Bluebell var.; polished) | 6 |

That is the paper's only table. The accompanying caution, printed immediately below it: "It should
be pointed out that we have not demonstrated that the method is accurate for the low value found for
Texas long-grain rice. However, it would be expected to be of the right order."

**Reconciliation with Buttery 1983 (mine).** The 1983 paper's brown-rice column is in ppm by weight;
multiplying it by 1000 and then by this paper's 3.57 recovery factor gives what the older data
should have said:

| variety (brown) | 1983 Table I, ppb (= ppm x 1000) | x 3.57 (mine) | 1986 Table I, ppb | ratio |
|---|---:|---:|---:|---:|
| Malagkit Sungsong | 200 | 714 | 760 | 1.06 |
| basmati 370 | 170 | 607 | 610 | 1.005 |
| IR 841-76-1 | 200 | 714 | 560 | 0.78 |
| Texas long grain (polished, 1983 milled) | <8 | <29 | 6 | — |

This is the arithmetic behind the paper's own sentence: "In previous analyses of these rice
varieties by the authors (Buttery et al., 1983), the low efficiency of recovery of
2-acetyl-1-pyrroline in the steam distillation continuous-extraction step had not been taken into
account. When the recovery factor is applied to the old figures, however, they do agree fairly well
with the present analyses." Two of the three agree to within 6 %; the third differs by 22 %, which
is inside the scatter one sample per variety with no replicates can carry.

### Figures (FIGURE-ONLY, no values read)

| figure | what it shows | printed caption details |
|---|---|---|
| 1 | capillary GLC of the **total** steam volatile oil from basmati rice; 2-acetyl-1-pyrroline is peak 62 | "using the 150 m X 0.5 mm i.d. Pyrex capillary wall coated with Carbowax 20-M"; time axis marked at 60 and 120 min |
| 2 | amount found by the method against amount added to a bland rice, with the found = added line drawn | axes to about 3000 ppb; no fitted line, no statistics printed |
| 3 | GLC of the **basic fraction** from basmati rice (A) against Texas long-grain rice (B) | "AP = 2-acetyl-1-pyrroline. TMP = collidine (2,4,6-trimethylpyridine) internal standard" |

## 4. Kinetic numbers the repository can use

**Registry mapping:** 2-acetyl-1-pyrroline -> `2_acetyl_1_pyrroline`. Collidine
(2,4,6-trimethylpyridine) is **not** in `data/keys/compounds.yml`; neither is hexanol. The
"interfering compounds (alkylpyrazines)" the burning warning names are the registry's `pyrazines`
class.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **relative recovery of 2-acetyl-1-pyrroline through the whole isolation, against collidine** | **28.0 (SD 4.9)** | % | spiked into 6 L water in the 12 L flask, no rice; 2 h Likens-Nickerson with acid solvent phase; nine GLC determinations over three isolations | p. 113, right column | measured_ratio (recovery), **within_study_ratio** |
| the multiplier that follows | 100 / 28.0 = 3.57 | — | as above | p. 113 | derived (the authors' own arithmetic) |
| ether/water partition coefficient, 2-acetyl-1-pyrroline / collidine | 0.74 / 1.0 | — | not otherwise specified | p. 113, right column | within_study_ratio |
| 2-acetyl-1-pyrroline in cooked rice | 760 / 610 / 560 / 6 | ppb (µg/kg of rice) | 200 g rice in 6 L water, 2 h atmospheric steam distillation, collidine internal standard, recovery-corrected | Table I, p. 114 | level_only — and see flag 1 on what the level means |
| the bland test rice used for spiking | < 10 | ppb | same method | p. 114, left column | level_only (a detection-limit statement) |
| Texas long grain, accuracy caveat | "we have not demonstrated that the method is accurate for the low value ... However, it would be expected to be of the right order" | — | — | p. 114 | level_only (order of magnitude) |
| retention times, 2-acetyl-1-pyrroline / collidine | 26 / 29 min (60 m Durabond wax, 50-230 C at 4 C/min); 22 / 26 min (1.2 m packed Carbowax, isothermal 90 C) | min | as printed | p. 113 | identification, not kinetic |
| any rate constant, half-life, barrier, temperature series or pH | **NOT PRESENT** | — | — | — | — |
| Figures 1, 2, 3 | chromatograms and the found-vs-added plot | — | — | Figs 1-3 | figure_only |

### What this paper says about the compound's stability — the whole of it

The engine carries no loss step for 2-acetyl-1-pyrroline. Everything this paper contributes toward
one is listed here, verbatim or near-verbatim, so that its **qualitative** character is unmistakable.

1. "The chemical analysis of 2-acetyl-1-pyrroline in rice by conventional methods is, however,
   difficult because of the presence of other interfering compounds and **the instability of this
   compound**." (Introduction.)
2. "The major difference in the recovery of 2-acetyl-1-pyrroline and the internal standard seems to
   occur **during the steam distillation step**. The reasons for this are not known although they
   **may be partly due to the instability of 2-acetyl-1-pyrroline**." (p. 113.) The recovery loss is
   72 % over 2 h at boiling in water — but it is a lumped loss over distillation, partition,
   neutralisation, single ether extraction and two concentration steps, and the authors themselves
   decline to attribute it. **It is not a decay rate and must not be read as one.**
3. "the method seems accurate enough for the purpose, considering the small concentrations involved
   and **the unstable nature of 2-acetyl-1-pyrroline**." (p. 114.)
4. "**Because of the instability of 2-acetyl-1-pyrroline it is advisable to use injector temperatures
   of 150-170 C. Higher injector temperatures seemed to cause some decomposition.**" (p. 114.) A
   bracket, not a rate: intact at 170 C on the timescale of an injection, decomposing above it.
5. In the standard solution the compound is stable: **290 ppm in benzene at -20 C, "quite stable
   over several months"** — consistent with Buttery 1983's observation that it is more stable
   dilute and in solution than neat.
6. The formation statement that governs how any measured level should be read: "The
   2-acetyl-1-pyrroline is **formed from the rice during the cooking process**, and what the method
   actually measures is **not the concentration in the rice at any particular time but the total
   amount of 2-acetyl-1-pyrroline produced from the given weight of rice during the heating (2-h)
   period. For this reason it is important to keep the degree of heating and the heating period the
   same for the different samples.**"

**The conclusion for the engine.** There is no measured loss rate for 2-acetyl-1-pyrroline anywhere
in this paper. What there is: a bracket at 150-170 C on a seconds timescale (item 4); a lumped 72 %
loss over a two-hour aqueous boil that the authors will not attribute (item 2); stability in cold
dilute organic solution (item 5). If a wave wants a decay step, these bound it loosely from both
sides and none of them can parameterise it. **A refusal is the correct answer today, and the wishlist
entry is a measured aqueous half-life at 100-140 C.**

## 5. Flags

1. **A measured level of this compound is a production integral, not a state.** The authors say so
   explicitly (section 4, item 6): the 2 h steam distillation *is* the cook, and Table I reports the
   total made over that 2 h from 200 g of rice, not a concentration at a time. A benchmark row that
   compares an engine concentration at a stated time against these 760 / 610 / 560 / 6 ppb is
   comparing two different quantities unless the engine integrates the same schedule. The 1983
   paper's own vacuum-versus-atmospheric comparison (a factor of about 1/10,
   `buttery1983_extraction.md` flag 3) is the size of the discrepancy this can produce.
2. **The recovery factor has an SD but the method has no accuracy statement.** 28.0 % with SD 4.9 %
   is 17.5 % relative, from nine determinations over three isolations — so the 3.57 multiplier
   carries about +/- 18 % before anything else. The method test (Figure 2) is reported in prose only,
   with no slope, no intercept, no r2, no number of spike levels and no replicate count. **The paper
   supports "this method is right to tens of percent", not better**, and Table I's three digits
   (760, 610, 560) over-state it.
3. **The recovery was measured in water, not in rice.** "adding known quantities of both
   2-acetyl-1-pyrroline and collidine to 6 L of water in the 12-L flask". The method test with a
   bland rice (Figure 2) is the only matrix check and it is not quantified. Whether 200 g of rice
   solids changes the recovery — by binding, by pH, or by the compound's own instability in a
   different medium — is not measured.
4. **Two printed inconsistencies in the column descriptions.** The Methods give the main column as
   "150-m length **0.5-mm i.d.** Pyrex glass capillary" and Figure 1's caption agrees, but the
   Results say "Figure 3 shows GLC analysis using the laboratory-constructed **0.66-mm i.d.** Pyrex
   capillary column"; and the packed column is "1.2-**m** length by 3-mm o.d." in the Methods and
   "1.2-**mm** length X 3-mm o.d." in the Results. Both are typographical and neither affects a
   number, but they are recorded because they are the kind of thing a reader might otherwise treat
   as two different columns.
5. **Two different CAS registry numbers for the same compound, in two papers by the same authors.**
   This paper prints **99583-29-6**; Buttery 1983 prints **85213-22-5**, which is what
   `data/species/desirable_targets.yml` stores. The repository should keep 85213-22-5 (it is the
   number in general use for 2-acetyl-1-pyrroline and the one the earlier paper assigned) but should
   know that a search on the 1986 number will not find the 1983 data, and vice versa. Neither paper
   comments on the discrepancy.
6. **A method-borne artefact that is directly a Maillard signal.** "Care should be taken to prevent
   any 'burning' of the rice on the bottom of the 12-L flask as this produces interfering compounds
   (**alkylpyrazines**)." Any rice volatile dataset produced this way can carry pyrazines that came
   from the flask bottom rather than from the cook, which matters if a future benchmark ever pairs
   2-acetyl-1-pyrroline with pyrazines from the same isolation.
7. **What the paper does NOT contain.** No rate constant of any kind; no half-life; no activation
   energy; no time course (one 2 h endpoint); no temperature series (one boiling-water isolation and
   one injector-temperature bracket); no pH; no water activity; no precursor experiment and no
   mention of proline; no Maillard model system; no odour threshold (that is in
   `buttery1983_extraction.md`); no statistics on Figure 2; no replicate levels in Table I; no
   milled/brown pairs (three of the four samples are brown, one is polished, so the 1983 paper's
   milling ratio cannot be re-measured here).
8. **Registry gaps against `data/keys/compounds.yml`:** collidine / 2,4,6-trimethylpyridine (the
   internal standard) is absent, as is hexanol. The target compound itself is present.
9. **What to request.** (i) An aqueous decay measurement for 2-acetyl-1-pyrroline at 100-140 C with
   time — the single number that would let the engine carry a loss step instead of refusing one.
   (ii) The recovery re-measured **in a rice matrix** rather than in water, which is what flag 3
   leaves open. (iii) Buttery, Ling & Juliano 1982 (Chem. Ind. 958), still not on disk, for the
   vacuum-isolation numbers that both this paper and the 1983 one refer to but neither reprints.
