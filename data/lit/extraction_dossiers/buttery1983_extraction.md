# Buttery, Ling, Juliano & Turnbaugh 1983 — EXTRACTION (2-acetyl-1-pyrroline in ten rice varieties by 2 h atmospheric steam distillation continuous extraction, and its odour threshold in water on a 16-judge panel)

### The source of the repository's 2-acetyl-1-pyrroline odour threshold: a Guadagni-plot determination printed as a VOLUME ratio, 0.1 part (mL) per 10^9 parts (mL) of water, plus a ten-variety level table that the authors themselves say is a semi-quantitative peak-area estimate.

**Source on disk:** `data/articles/buttery1983.pdf` (J. Agric. Food Chem. 1983, 31 (4), 823-826; the
PDF's first page carries the tail of the preceding article's reference list and its last page the
head of the following one). The text layer is an OCR layer with glyph damage in prose ("parts/lOg"
for "parts/10^9", "Malagkit"/"Malagkit" inconsistencies, a dropped fraction on page 825).
**Tables I and II were verified against 200-400 dpi rasters of printed pages 824 and 825**
(`scratchpad/img/bu83t1-2.png`, `bu83t2-3.png`, `bu83t2h2-3.png`, `bu83thr2-3.png`); every number
below matches the raster, and the dropped fraction on page 825 was recovered from the raster as
**1/10** (`scratchpad/img/bu83fracz2-3.png`). Table III (panel odour descriptions) and Table IV (a
rank order with no numbers) are re-typed in full. Figure 1 is a synthesis scheme with no numbers.

## 0. Identity

| field | value |
|---|---|
| Title | "Cooked Rice Aroma and 2-Acetyl-1-pyrroline" |
| Authors | Ron G. Buttery, Louisa C. Ling and Jean G. Turnbaugh (Western Regional Research Center, ARS, USDA, Berkeley CA) and Bienvenido O. Juliano (International Rice Research Institute, Manila) |
| Venue | J. Agric. Food Chem. 1983, 31 (4), 823-826. Received December 13 1982; accepted March 21 1983 |
| Registry numbers printed | 2-acetyl-1-pyrroline **85213-22-5**; 2-acetylpyrrole 1072-83-9; 2-(1-hydroxyethyl)pyrrolidine 63848-93-1 |
| Naming | "ppm" in Table I is explicitly defined by footnote a as parts by **weight** per 10^6 parts of rice, dry weight. "parts/10^9 parts of water" in Table II and the threshold sentence is explicitly a **volume** ratio, "part (mL) ... per 10^9 parts (mL)" |
| Companion | Buttery, Ling & Juliano 1982, Chem. Ind. (London) 958 (the first identification of 2-acetyl-1-pyrroline in rice; **not on disk**) and Buttery, Ling & Mon 1986, JAFC 34:112, the quantitative method (`buttery1986_extraction.md`, on disk) |
| Threshold method lineage | Guadagni & Buttery 1978, J. Food Sci. 43:1346 (the procedure) and Guadagni, Maier & Turnbaugh 1973, J. Food Sci. 38:1277 (the plotting), both cited and **neither on disk** |
| Repo registry | `2_acetyl_1_pyrroline` in `data/keys/compounds.yml`; the target row is in `data/species/desirable_targets.yml` at line 266 |

## 1. Why it matters

`data/species/desirable_targets.yml` carries 2-acetyl-1-pyrroline with
`odour_threshold_ug_per_kg: 0.1` and names this paper as the source. That single number is what makes
the compound the highest-potency member of the desirable list and therefore what drives every odour
activity value the engine computes for it; it is also the threshold that
`kinetic_core_b19_prereg_draft.md` section 2 assumes exists when it lists 2-acetyl-1-pyrroline as
the proline target to build toward. Until now the number has been documented only through the
target-file comment. `k3_final_parameter_inventory.md` line 592 recorded the opposite — "**NO
THRESHOLD EXISTS IN EITHER 2-AP PAPER** ... GAP — 2-AP's whole claim to being a top odorant is an
OAV argument the repo cannot score" — and line 1152 and line 1644 opened a retrieval request for
Buttery, Ling & Juliano 1982 as "the **only** route to a 2-AP odour threshold". **That entry is
wrong and this dossier closes the request:** the threshold is printed in this paper, on page 825,
with its panel size, its raw percent-correct series in Table II and its plotting method named. The
1982 Chem. Ind. note is not needed for it.

The paper matters a second way that the target file does not use. It prints the compound's
**instability**, in prose and with a mechanism — the pure liquid darkens on storage at -20 C under
vacuum, it "could not be gas chromatographed" on the authors' packed silicone or Carbowax columns,
and it "may be more stable in dilute solution". The engine has no loss step for 2-acetyl-1-pyrroline;
this and `buttery1986_extraction.md` are where the evidence for one would come from.

And it matters a third way, negatively, which section 5 flag 2 states: the ten-variety level table
is not a quantitative measurement, the authors say so twice, and the same authors withdrew it three
years later by a factor of 3.57 (`buttery1986_extraction.md`).

## 2. Methods as they matter to a model

**Two different isolation methods appear in this paper and they differ by a factor of ten.**

- **The method used for Table I (atmospheric).** 500 g of rice + 6 L of water in a 12 L flask,
  Likens-Nickerson steam distillation continuous extraction head, 125 mL freshly distilled diethyl
  ether in a 250 mL solvent flask, **2 h at atmospheric pressure**. "**Cooking occurs during the
  isolation process.**" Ether extract dried over Na2SO4 and concentrated to 0.15 mL on a warm water
  bath through low-hold-up columns. Volatile oil yield "of the order of 5 parts per million (ppm)".
- **The authors' earlier method (vacuum), Buttery et al. 1982.** Vacuum steam distillation continuous
  extraction of **already cooked** rice. Page 825, verbatim: "This vacuum isolation method probably
  gives a more accurate figure for the amount actually present in the cooked rice. **This figure was
  about 1/10 that found by using the atmospheric isolation method.**" So Table I's numbers are about
  ten times the amount present in cooked rice by the authors' own preferred method.
- **What Table I therefore measures**, in the authors' words: "The method of isolation of the
  volatiles used in the present work gives a figure for the **total 2-acetyl-1-pyrroline produced by
  the rice during the isolation period**. In the normal cooking of a rice for food, much of this
  would be lost and the amount left in the rice would be considerably less."
- **Basic fraction (used for some samples only).** Ether concentrate into 50 mL hexane, extracted
  with 3 N HCl (3 x 25 mL), acid washed with ether (1 x 50 mL), neutralised with excess NaHCO3 under
  100 mL ether with ice cooling, aqueous re-extracted 2 x 50 mL ether, dried, concentrated to
  0.01 mL. "GLC analysis of this basic fraction gave an even better separation."
- **Analysis.** 150 m x 0.64 mm i.d. Pyrex glass capillary coated with Carbowax 20M into a modified
  Consolidated 21-620 cycloidal MS through a Llewellyn-Littlejohn single-stage silicone rubber
  membrane separator; 50 C for 30 min after injection, then 50 to 170 C at 1 C/min, 2 h final hold.
  2-Acetyl-1-pyrroline elutes at **Kovats index 1320**, "just before hexanol".
- **Quantification: peak area, no internal standard, no response factor, one sample per variety.**
  Page 824, verbatim: "**The figures in Table I are meant only to give a general idea** of the
  variation of 2-acetyl-1-pyrroline with the different varieties. The GLC peak area method was used,
  making the usual assumption that all compounds, in the complex mixture, have the same response in
  the flame ionization detector. **Only one sample from each variety was studied and no attempt was
  made to study the variation within any variety.** The difference found between varieties such as
  Malagkit Sungsong (0.09 ppm) and Basmati 370 (0.06 ppm) **cannot be considered meaningful** for
  this type of study. The factor of about 10 times difference, between the more aromatic Asian rice
  varieties and the American Calrose and Texas Long Grain varieties, **is certainly meaningful**."
- **Identification.** Mass spectrum and GLC retention against a synthetic sample. The printed
  spectrum: "molecular ion at m/e 111 (5), other major ions at 43 (100), 41 (50), 42 (24), 83 (11),
  69 (11), 68 (8), 55 (2), 52 (0.9), 54 (0.2), 67 (0.2)". IR (CCl4): "major absorption maxima at
  1695, 1620, 1435, 1370, 1340, 1250, 1080, 1000, 975, and 940 cm-1 in the 2000-600 cm-1 region".
- **Synthesis (Figure 1), for the record.** 2-acetylpyrrole (1.7 g) in 50 mL methanol, hydrogenated
  over 5 % rhodium on alumina (2.0 g), room temperature, 10 psi H2, 15 h, stirring -> 1.8 g crude
  2-(1-hydroxyethyl)pyrrolidine; that (1.8 g) into stirred silver carbonate on Celite (16 g) in
  100 mL benzene under nitrogen, refluxed 15 h; filtered, concentrated to 5 mL; isolated by GLC on a
  2 m x 0.64 cm o.d. aluminium column packed with 15 % Amine 220 on 60-80 mesh Chromosorb P.
  **Overall yield from 2-acetylpyrrole: 10 %.** Stored in sealed 3 mm o.d. Pyrex tubes under vacuum
  at -20 C.
- **Threshold determination.** Page 825, verbatim: "An odor threshold of 2-acetyl-1-pyrroline was
  determined in water solution by using established procedures [e.g., Guadagni and Buttery (1978)]
  **with a trained panel consisting of 16 judges**. As in previous threshold determinations the odor
  judges were presented with **two Teflon squeeze bottles, one containing the solution and the other
  odor-free water**. The task for each judge was to determine which of the coded bottles contained
  the odorant. Table II lists the results obtained. **Plotting the data as outlined by Guadagni et
  al. (1973) gave a threshold of 0.1 part (mL) of compound per 10^9 parts (mL) of water.**" Teflon
  bottles and tubing throughout (Buttery et al. 1981). Odour **quality** evaluations used 100 mL
  opaque Pyrex flasks.
- **What the threshold procedure is.** A two-alternative forced choice against odour-free water,
  repeated at eight concentrations with 16 to 113 total judgments per concentration; the threshold
  is read off a plot of percent-correct against log concentration, by the Guadagni 1973 construction.
  The paper does not print the criterion percentage, the plot, or a confidence interval.

## 3. Tables re-typed

### Table I. "Concentration of 2-Acetyl-1-pyrroline Found in Cooked Rice Varieties in Terms of Dry Weight of Rice"

Header as printed: "2-acetyl-1-pyrroline concn, ppm^a", split into "milled rice" and "brown rice".
Footnote a: "ppm = parts (weight) of compound per million (10^6) parts of rice (dry weight)." Blank
cells are blank in the printed table (no brown-rice sample was run).

| variety | milled rice (ppm) | brown rice (ppm) |
|---|---:|---:|
| Malagkit Sungsong | 0.09 | 0.2 |
| IR841-76-1 | 0.07 | 0.2 |
| Khao Dawk Mali 105 | 0.07 | 0.2 |
| Milagrosa | 0.07 | — |
| Basmati 370 | 0.06 | 0.17 |
| Seratus Malam | 0.06 | — |
| Azucena | 0.04 | 0.16 |
| Hieri | 0.04 | 0.1 |
| Texas Long Grain | <0.008 | — |
| Calrose | <0.006 | — |

Sample provenance, from Experimental: 1981 crops of Azucena, IR841-76-1 (a line derived from Khao
Dawk Mali 105) and Milagrosa (Philippines), Basmati 370 (Pakistan), Hieri (Japan), Khao Dawk Mali
105 (Thailand), Seratus Malam (Indonesia), and the 1982 crop of Malagkit Sungsong (Philippines), all
received as brown rice and **milled in the laboratory removing about 10 % of the outer layers**;
Calrose (milled) from local markets in Berkeley CA, 1982; Texas Long Grain (milled) from Comet Rice
Mills, Houston TX, "probably Labelle variety". Most varieties were obtained through IRRI, Manila.

Note the abstract and page 824 disagree slightly with the table on Calrose: the abstract says "less
than 0.006 parts per million (ppm) for Calrose", the table says "<0.006", and page 824 says "the
Calrose the least amount with 0.006 ppm" (no "less than"). **The table's "<0.006" is the form to
carry**; it is a detection limit, not a measurement. Milled-to-brown ratios, where both exist
(mine): 0.09/0.2 = 0.45, 0.07/0.2 = 0.35, 0.06/0.17 = 0.35, 0.04/0.16 = 0.25, 0.04/0.1 = 0.40 —
i.e. milling to remove about 10 % of the outer layers removes 55-75 % of the compound.

### Table II. "Odor Threshold Determination of 2-Acetyl-1-pyrroline"

Header exactly as printed: "concn, parts/10^9 parts of water | % correct judgments | total no.
judgments".

| concn, parts/10^9 parts of water | % correct judgments | total no. judgments |
|---:|---:|---:|
| 7 | 100 | 16 |
| 3.5 | 94 | 16 |
| 0.9 | 94 | 16 |
| 0.35 | 95 | 19 |
| 0.18 | 86 | 52 |
| 0.09 | 75 | 81 |
| 0.045 | 65 | 113 |
| 0.023 | 53 | 94 |

Total judgments across the eight rows (mine): 407. The series behaves as a forced-choice curve
should: 53 % at the lowest level is chance (50 %), and the printed threshold of 0.1 falls between
the 0.09 row (75 % correct) and the 0.18 row (86 %). **The 0.1 in the threshold sentence is not a
row of this table**; it is the value read off the Guadagni plot of these eight points, and 0.09 —
the concentration at which the panel was 75 % correct, i.e. halfway between chance and certainty —
is its nearest measured neighbour. That coincidence is worth stating because it makes the printed
threshold legible as a criterion rather than as a fitted extrapolation.

### Table III. "Most Often Used Panel Odor Descriptions of a 0.05-ppm Solution of 2-Acetyl-1-pyrroline and of Cooked Malagkit Sungsong Rice"

Two independent panels; percentages are "% of judges using description".

| odour description | 0.05-ppm solution of 2-acetyl-1-pyrroline in water (22 judges) | cooked Malagkit Sungsong rice (23 judges) |
|---|---:|---:|
| popcorn | 82 | 60 |
| cooked oatmeal | 45 | 56 |
| cooked rice | 23 | 48 |
| sweet | 23 | 30 |
| nutty | 17 | 14 |

The 0.05 ppm solution is **500 times the threshold** (mine, 0.05 ppm = 50 parts per 10^9 against a
threshold of 0.1).

### Table IV. "Ranking of Cooked Rice Samples in Terms of Those Having the Greatest Popcorn-like Aroma (at the Top of the Table) to Those Having the Least (at the Bottom of the Table)"

No numbers are printed; the table is an ordered list with a "greatest popcorn aroma" arrow at the top
and "least popcorn aroma" at the bottom. Order as printed, top to bottom:

Malagkit Sungsong · Milagrosa · Khao Dawk Mali 105 · IR841-76-1 · Basmati 370 · Seratus Malam ·
Azucena · Hieri · Calrose · Texas Long Grain

The panel was 21-23 judges, presented with groups of two or three cooked samples with common samples
so the groups could be overlapped into one ranking. The rank order agrees with Table I's
concentration order except that Milagrosa (0.07) outranks Khao Dawk Mali 105 (0.07) and
IR841-76-1 (0.07) — a tie in the table — and that Calrose (<0.006) outranks Texas Long Grain
(<0.008), which is a reversal inside the two detection limits.

### The two difference tests (prose, page 826, not tabulated)

| test | judgments | correct | printed verdict |
|---|---:|---:|---|
| Calrose vs Malagkit Sungsong, matching an unknown to a labelled control | 41 | 83 % | "highly significant data that the panel can tell the difference" |
| same test after adding 25 mL of a 0.05-ppm 2-acetyl-1-pyrroline solution to each Calrose sample (and 25 mL odour-free water to each Malagkit Sungsong sample) | 40 | 62 % | "only slightly better than pure chance where the correct sample would be matched 50% of the time" |

This is the paper's actual causal argument: doping the bland variety with the compound abolishes a
difference the panel could otherwise make.

## 4. Kinetic numbers the repository can use

**Registry mapping:** 2-acetyl-1-pyrroline -> `2_acetyl_1_pyrroline` (CAS 85213-22-5, which matches
the target file's `cas`). 2-acetylpyrrole, 2-(1-hydroxyethyl)pyrrolidine, 1-pyrroline, hexanol and
2-acetyl-1,4,5,6-tetrahydropyridine are not in `data/keys/compounds.yml`.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **odour threshold of 2-acetyl-1-pyrroline in water** | **0.1 part (mL) of compound per 10^9 parts (mL) of water** | volume ratio (v/v) | water solution, 16-judge trained panel, two-bottle forced choice against odour-free water, Guadagni et al. 1973 plot of the Table II series | p. 825, right column | **threshold** |
| the percent-correct series the threshold was read from | 100 / 94 / 94 / 95 / 86 / 75 / 65 / 53 % at 7 / 3.5 / 0.9 / 0.35 / 0.18 / 0.09 / 0.045 / 0.023 parts per 10^9 | % correct, and parts per 10^9 by volume | as above; 16 to 113 judgments per level, 407 in total | Table II, p. 825 | threshold (the raw data behind it) |
| 2-acetyl-1-pyrroline in cooked rice, ten varieties, milled | 0.09 down to <0.006 | ppm by weight, dry rice | 500 g rice + 6 L water, 2 h atmospheric Likens-Nickerson, FID peak area, one sample per variety | Table I, p. 824 | level_only — and see flag 2 |
| the same, brown rice, five varieties | 0.2 / 0.2 / 0.2 / 0.17 / 0.16 / 0.1 | ppm by weight, dry rice | as above | Table I | level_only |
| the vacuum-vs-atmospheric isolation ratio | atmospheric x **1/10** = the vacuum figure, which "probably gives a more accurate figure for the amount actually present in the cooked rice" | — | the authors' two methods on the same material | p. 825, left column | within_study_ratio |
| milled / brown ratio | 0.25 to 0.45 (mine, five pairs) | — | ~10 % of the outer layers removed | Table I | within_study_ratio |
| aromatic Asian varieties / American varieties | "about 10 times", declared meaningful by the authors | — | Table I | p. 824 | within_study_ratio |
| synthesis yield, 2-acetylpyrrole to 2-acetyl-1-pyrroline | 10 % overall | % | two steps, Figure 1 | p. 823-824 | level_only (not a Maillard number) |
| descriptor profile at 0.05 ppm | popcorn 82 %, cooked oatmeal 45 %, cooked rice 23 %, sweet 23 %, nutty 17 % (22 judges) | % of judges | water, opaque flasks | Table III | level_only (sensory) |
| difference test, doped Calrose | 83 % correct (41 judgments) falls to 62 % (40 judgments) when 25 mL of 0.05 ppm 2-acetyl-1-pyrroline is added to the bland variety | % correct | 100 mL opaque flasks | p. 826 | within_study_ratio (sensory) |
| any rate constant, barrier, time course or temperature dependence | **NOT PRESENT** | — | — | — | — |
| Figure 1 | a synthesis scheme, no numbers | — | — | Figure 1 | figure_only |

### Is the repository's stored 0.1 ug/kg a sound reading of what is printed? — YES, with one declared assumption

The printed quantity is a **volume ratio**: 0.1 mL of compound in 10^9 mL of water. The stored
quantity in `data/species/desirable_targets.yml` is a **mass ratio**, `odour_threshold_ug_per_kg:
0.1`. The conversion (mine):

    0.1 mL compound in 1e9 mL water
      = 0.1 x rho grams of compound in 1e9 g = 1e6 kg of water        (water taken as 1.000 g/mL)
      = 1e5 x rho micrograms per 1e6 kg
      = 0.1 x rho  ug/kg

with rho the density of liquid 2-acetyl-1-pyrroline in g/mL. So **0.1 ug/kg is exactly right if
rho = 1.00 g/mL**, and the stored number is the printed number under that single assumption, which
the target file's comment already states ("printed as a volume ratio, 0.1 ug/kg on unit density").

The assumption is sound, for three reasons that should travel with it.

1. **Neither density is printed.** The paper gives no density for 2-acetyl-1-pyrroline (C6H9NO,
   MW 111.14) and none is anywhere on disk. The number 1.00 g/mL is the repository's, not the
   authors'.
2. **The error it can introduce is small and bounded by chemistry.** For the stored 0.1 ug/kg to be
   wrong by a factor of two, the compound would have to be a liquid of density 0.5 or 2.0 g/mL. A
   small oxygen- and nitrogen-bearing five-membered heterocyclic ketone cannot be either; liquids of
   this composition sit close to water. The realistic error is at most some tens of percent, and
   it is one-sided (a density above 1 makes the true threshold slightly higher than 0.1 ug/kg, i.e.
   the stored value is slightly conservative for odour activity).
3. **It is far inside the determination's own resolution.** Table II's points are spaced by a factor
   of about two, the criterion percentage is not printed, no confidence interval is given, and the
   two nearest measured levels are 0.09 (75 % correct) and 0.18 (86 %). A threshold read off eight
   points on a two-alternative forced choice is not a two-significant-digit quantity. **The honest
   statement is "0.1 ug/kg in water, order of magnitude, panel of 16"**, and a density correction
   would be lost inside that.

What is **not** sound, and is a separate matter from the density, is treating this water threshold as
a threshold in any other medium. `k2_matrix_and_thresholds.md` is the repository's place for that
distinction; nothing in this paper measures a threshold in oil, in a rice matrix, or at any
temperature other than the panel booth's.

## 5. Flags

1. **The threshold is a volume ratio and the repository stores a mass ratio.** The two coincide only
   because the compound's density is taken as water's. Section 4 states the arithmetic and declares
   the assumption; the target file's note already carries it. Any future re-expression of the
   threshold (in mol/L, say — 0.1 ug/kg / 111.14 g/mol = **0.90 nmol/L**, mine, and inheriting the
   same density assumption) must carry it too.
2. **Table I is not a quantitative measurement, and the same authors reduced it threefold three
   years later.** The authors call it "meant only to give a general idea", they assume equal FID
   response for every compound in a complex mixture, they ran one sample per variety with no
   replicates, and they say a 0.09-vs-0.06 difference "cannot be considered meaningful".
   `buttery1986_extraction.md` then measured the steam-distillation recovery of the compound as
   **28.0 %** and multiplied its own results by 3.57 — and states that when the recovery factor is
   applied to "the old figures", i.e. these, "they do agree fairly well with the present analyses".
   **Table I's numbers are uncorrected for a recovery the same laboratory later measured at 28 %.**
   Do not use Table I as a level benchmark; use the 1986 table.
3. **And in the other direction, Table I over-reports what is in cooked rice by about ten.** The
   authors' vacuum method, which they call more accurate for the cooked-rice content, gave "about
   1/10" of the atmospheric figure, because the 2 h atmospheric isolation is itself a cook: "the
   method ... gives a figure for the total 2-acetyl-1-pyrroline produced by the rice during the
   isolation period". Flags 2 and 3 push in opposite directions (a 3.57x recovery correction up, a
   10x isolation-artefact correction down) and neither has an uncertainty. **The safe reading of
   Table I is: an ordering of varieties spanning about tenfold, and no absolute level.**
4. **The compound is unstable, in three separate printed statements, and the engine has no loss step
   for it.** (i) Neat, sealed under vacuum at -20 C, it "slowly turned to a red color which became
   darker the longer the storage"; the authors speculate "a conjugated pyrroline polymer ... formed
   by condensation of the carbonyl groups with the 5-positions of other molecules". (ii) "The
   compound showed considerable instability to general gas chromatography conditions and could not
   be gas chromatographed by using the authors' silicone or Carbowax 20M **packed** columns. This
   may explain why it was not detected in the earlier studies of rice volatiles." (iii) "For this
   reason 2-acetyl-1-pyrroline may be more stable in dilute solution. In the authors' experience
   this seemed to be true with dilute water solutions." **None of this is a rate.** There is no
   half-life, no temperature, no concentration dependence and no time. It is a qualitative warning
   that a second-order self-condensation exists and that the compound survives better dilute — which
   is the regime the engine works in — and it is the reason `buttery1986_extraction.md` exists.
5. **The threshold has no confidence interval, no criterion percentage and no plot.** Table II's
   407 judgments are printed, the Guadagni 1973 construction is named, and the plot is not shown.
   The two papers that define the procedure (Guadagni & Buttery 1978; Guadagni, Maier & Turnbaugh
   1973) are **not on disk**; `k3_final_parameter_inventory.md` line 1655 already lists Guadagni,
   Buttery & Turnbaugh 1972 as a retrieval item because four of six aqueous thresholds in the
   repository divide by it. Carry the threshold as an order of magnitude.
6. **The panel size is stated twice with two different numbers, for two different tests.** The
   threshold panel is "a trained panel consisting of **16 judges**" — and Table II's "total no.
   judgments" column runs to 113, so judges were presented repeatedly at the lower levels. The
   odour-quality panels are separate and larger: 22 judges for the solution, 23 for the rice
   (Table III), 21-23 for the ranking (Table IV), 41 and 40 judgments for the two difference tests.
   Only the 16 belongs with the threshold.
7. **What the paper does NOT contain.** No formation rate, no barrier, no time course, no
   temperature series, no pH; no precursor experiment (the "Possible Origin" section is speculation
   about proline and hydroxyproline, and ends "the authors were unable to carry this out in the
   laboratory by normal chemical means"); no measurement in a Maillard model system of any kind; no
   threshold in any medium other than water; no density, vapour pressure or partition coefficient;
   no stability rate. The engine's route to the compound (proline Strecker, then acylation of
   1-pyrroline) comes from `hofmann1998b_extraction.md`, not from here, and this paper cannot
   corroborate it.
8. **Registry gaps against `data/keys/compounds.yml`:** the target compound is present as
   `2_acetyl_1_pyrroline`. Absent: 2-acetylpyrrole, 1-pyrroline, 2-acetyl-1,4,5,6-tetrahydropyridine
   and 2-acetyl-2-thiazoline — the last two are the "cracker-like" comparators the paper names
   ("2-Acetyl-1-pyrroline seems to be the most potent of the 'cracker-like' group of odor compounds
   which includes 2-acetyl-1,4,5,6-tetrahydropyridine, 2-acetylpyrazine, 2-acetyl-2-thiazoline"),
   and the repository does carry 2-acetylthiazole and 2-acetylfuran, so the family is half-named.
9. **A correction to record in the inventories.** `k3_final_parameter_inventory.md` line 592 states
   "NO THRESHOLD EXISTS IN EITHER 2-AP PAPER" and lines 1152 and 1644 open a retrieval for Buttery,
   Ling & Juliano 1982 as the only route to it. Both are superseded by this dossier: the threshold
   is printed here, with its panel and its raw series. The 1982 Chem. Ind. note remains
   un-retrieved and is still the first identification of the compound in rice, but it is no longer
   needed for the threshold.
10. **What to request.** (i) Guadagni, Maier & Turnbaugh 1973 and Guadagni & Buttery 1978, to
    recover the threshold criterion and to put an interval on the 0.1 — they are cited by four
    other repository thresholds as well. (ii) A published density for 2-acetyl-1-pyrroline, which
    would close flag 1 outright. (iii) Any measurement of the compound's decay with time and
    temperature in dilute aqueous solution: the engine's missing loss step needs a rate, and neither
    this paper nor the 1986 one supplies one.
