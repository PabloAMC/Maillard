# van Seeventer, Weenen, Winkel & Kerler 2001 — EXTRACTION (a ribose/cysteine process flavouring made at 130 C, then stored in the dark at 50 C for 24 h in 0.5 M phosphate pH 5.0 under air and under argon; five character-impact compounds followed at 0/4/8/24 h by GC-FID against maltol; losses reported as ZERO-ORDER % per day, plus one mass balance that does not close and one H-D exchange experiment)

### THE PAPER OF THIS CLUSTER THAT CARRIES THE SINK NUMBER, and it is a declared HOLD-OUT so the objective has legitimately never seen it: **2-methyl-3-furanthiol disappears at 59 % per day at 50 C and cannot be found afterwards either as the thiol or as any (mixed) disulfide, while 2-mercapto-3-butanone, 2-furfurylthiol and 2,5-dimethyl-3-furanthiol all close their balances in the same experiment.** That is a measured, non-oxidative, MFT-specific sink at a temperature 95 C below the sulfur module's reference, with air ≈ argon. And a second number, mine, that nobody has used: the "zero-order" loss is **not** zero-order across pots — a 10-fold larger thiol pool loses **at least 15 times more per day**, so the sink is not a fixed-capacity drain.

**Source on disk:** `data/articles/vanseeventer2001.pdf` (4 pp., J. Agric. Food Chem. 2001, 49 (9),
4292-4295). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/vanseeventer2001.txt`, 293 lines). **Table 1 — the paper's only table — came
through clean** and is re-typed in full below, including both footnotes. Figures 1, 3 and 5 are the
paper's quantitative plots and are **figure_only**: Figure 1 (loss rate against cysteine
concentration for the three thiols), Figure 3 (H-D exchange in MFT over 64 days) and Figure 5 (H-D
exchange in bis(3-furanyl) disulfide). Figures 2 and 4 are a proposed mechanism and a synthesis
scheme. There is no supplementary material. Repo status before this dossier: this paper is named in
`src/kinetic_core/parameters_sulfur.py` **`NO_MEASURED_RATE["k_oligomer"]`** (the one channel the
module declares structurally and holds at exactly 0.0), in `PROHIBITED_DERIVATIONS` ("van
Seeventer's 59 %/day bolted on as THE MFT sink"), in the `acid_catalysed_C5_oligomerisation` entry of
`THIOL_CHANNELS`, and in `src/kinetic_core/species_sulfur.py` as the note on the `OLG` species. In
`docs/reference/FIT_HOLDOUT_DECLARATION.md` D.4 its Table 1 is a **★ HOLD-OUT** ("a third
temperature, a third mechanism, and **zero order** — it tests the *functional form*, not just the
magnitude"), while D.3 makes its **precursor conversion (55 % / 75 %) a FIT row**. **It has never had
a dossier of its own.**

## 0. Identity

| field | value |
|---|---|
| Title | "Stability of Thiols in an Aqueous Process Flavoring" |
| Authors | Paul B. van Seeventer, Hugo Weenen, Chris Winkel (corresponding), Josef Kerler — Quest International, Huizerstraatweg 28, 1411 GP Naarden, The Netherlands. (Present addresses at publication: van Seeventer at Zuivelfabriek De Kievit bv, Meppel; Weenen at TNO Nutrition and Food, Zeist.) |
| Venue | J. Agric. Food Chem. 2001, 49 (9), 4292-4295. Received 14 March 2001, revised 1 June 2001, accepted 1 June 2001, web 10 August 2001 |
| DOI / article ID | 10.1021/jf010348t (printed as `JF010348T`) |
| Naming | MFT = 2-methyl-3-furanthiol; FFT = 2-furfurylthiol; **MB** = 2-mercapto-3-butanone in the abstract and Table 1 but "3-mercapto-2-butanone" in the Results text (Flags 2); HDF = 4-hydroxy-2,5-dimethyl-3(2H)-furanone (the registry's `hdmf` / furaneol); sotolone = 3-hydroxy-4,5-dimethyl-2(5H)-furanone; DMFT = 2,5-dimethyl-3-furanthiol |
| Lineage | the process flavouring is **Hofmann & Schieberle's own** ribose/cysteine model (refs 1 and 2 = the Hofmann 1995 dissertation and this cluster's `hofmann1995` paper), reproduced "according to the procedure of Hofmann and Schieberle (2)" |
| Companions on disk | this cluster's `hofmann1995_extraction.md` (the pot this paper stores) and `hofmann2001_extraction.md` (the 30 C and 80 C covalent sink); `zhang2024b_extraction.md` and `zhou2023_extraction.md` (the 115-120 C dimerisation); `charlesbernard2005_extraction.md` (the 25 C ladder) |

## 1. Why it matters

**What it contributes to the thiol-sink question. This is the paper that carries a number, and it
carries three.**

1. **A thiol loss rate at a temperature other than 145 C — at 50 C, in the sulfur lane's own
   reference chemistry.** Table 1 entry 1: in the real ribose/cysteine process flavouring, **MFT
   59 %/day, FFT 28 %/day, MB 14 %/day, HDF < 10 %/day, sotolone not determined**. In the
   reconstituted mixture with the matrix stripped out (entry 2), **all three thiols are > 90 %/day**.
   That is a 145 C-free measurement of exactly the quantity the refused waves were trying to find,
   in a pot made from ribose and cysteine at pH 5.
2. **A mass balance that does not close, stated in the paper's own words** (p. 4294): "The amount of
   MB and FFT at the beginning was equal to the total amount of both thiols and (mixed) disulfides
   at the end of storage. This was absolutely not the case for MFT. After storage a large amount of
   MFT was not detected, neither as thiol nor as (mixed) disulfide, using GC-analysis." The same
   experiment with **DMFT closes**. So: the sink is real, it is specific to MFT, it is **not** the
   disulfide, and it is **not** oxidative — "there was almost no difference between storage of the
   aqueous model process flavoring under air or argon atmosphere". The conclusions put it plainly:
   "The instability is not due to disulfide formation, but appears to result from electrophilic
   coupling reactions."
3. **A number nobody has taken from it (mine): the loss is NOT zero-order across pots.** The paper
   fits zero order *within* each run and the module's channel entry repeats that ("ZERO in thiol").
   But the process flavouring's thiols sit "in the order of ... around 5 µM" (Table 1 footnote a)
   and the reconstituted mixtures at **50 µM** — ten times higher — and the losses per day go the
   wrong way for a zero-order sink. **MFT: 59 %/day of 5 µM = 2.95 µM/day, against > 90 %/day of
   50 µM = > 45 µM/day. A factor of at least 15 across a 10-fold change in the pool.** A truly
   zero-order, capacity-limited drain would have given the same µM/day in both. Even entry 4, with
   150 mM cysteine restored, gives 69 %/day of 50 µM = 34.5 µM/day, still **12x** the process
   flavouring's absolute rate. Either the order is nearer first than zero, or the intact Maillard
   matrix protects the thiol by a factor of about ten and a half — and the paper argues the latter
   ("an anti-oxidative effect of the matrix"). Either reading is new information for a sink
   objective, and the two are separable by an experiment nobody has run.

**And what it does NOT contribute, which must be said just as plainly.** Table 1 is a **★ HOLD-OUT**
and this dossier does not license fitting to it. `NO_MEASURED_RATE["k_oligomer"]` records the
consequence the module pre-registered: the model predicts no oligomerisation loss at 50 C and will
therefore fail this row, and "that failure is informative and is reported as such". Nothing in this
dossier changes that. What the three refused sink structures could legitimately have taken from this
paper *before* fitting anything is the **structural** brief — a sink that (a) is not the disulfide,
(b) is not oxygen-driven, (c) discriminates sharply between MFT and DMFT (a 5-methyl group turns it
off), (d) is suppressed by residual cysteine with a **minimum at about 50 mM cysteine for MFT and a
monotone dependence for MB and FFT**, and (e) is acid-catalysed and electrophilic at **C-5**. The
B17/B25 record reaches candidate (c) of its own accord — "the unsaturated-carbonyl adducts ... an
IRREVERSIBLE sink" — and that is the same *kind* of object as this paper's electrophilic coupling,
but it is not the same object: this paper's electrophile is **the protonated thiol itself**
(Figure 2), i.e. the sink is a self-reaction and needs no partner pool at all. **That is a fourth
structure the pre-registrations never named**, and it is the only one of the four whose rate would
be second order in the thiol and would therefore *rise* with the fed-thiol loadings of the 145 C
panel rather than fighting them.

**What the code already claims and whether it holds.** Every claim verifies against the print:
"59 % of initial per DAY for MFT, 28 % for FFT" (Table 1 entry 1); "air ~ argon (so it is NOT
oxidative)" (p. 4293, but **"data not presented here"** — Flags 3); "the MFT mass balance FAILS to
close as thiol + disulfide while MB, FFT and DMFT all close" (p. 4294); "85 % H-D exchange at C-5
against 10 % at C-4" (p. 4295). The **precursor conversion FIT row (55 % / 75 %)** also verifies, by
my arithmetic: residual cysteine **15 mM** from a 33 mM charge is **54.5 % converted** and residual
ribose **25 mM** from a 100 mM charge is **75 % converted**.

What this paper does NOT give the repository: any rate constant with a unit of inverse time; any
activation energy; any temperature other than 50 C for the storage (the 130 C is the *making* of the
pot, with no thiol data); any pH other than 5.0; any identification of the missing MFT; any
quantification of the disulfides that did form; any yield of anything from the 130 C cook; and any
air-versus-argon number.

## 2. Methods as they matter to a model

- **The process flavouring (entry 1).** **180 mmol D-ribose + 60 mmol L-cysteine in 1800 mL of
  0.5 M phosphate buffer, pH 5.0** — i.e. **ribose 100 mM, cysteine 33.3 mM (mine)**, matching
  Table 1 footnote a. Heated in a **2 L autoclave from room temperature to 130 C in 10 min, held at
  130 C for 20 min**, then rapidly cooled with tap water. The character-impact components come out
  "in the order of concentrations of around **5 µM**".
- **Storage of the process flavouring.** **100 mL portions, in the dark, at 50 C, in closed bottles,
  under both air and argon.** Sampled at **0, 4, 8 and 24 h**.
- **The reconstituted mixtures (entries 2-5).** MFT, MB, FFT, HDF and sotolone from an ethanolic
  stock into 0.5 M phosphate pH 5.0 to **50 µM each** — "approximately 10-fold higher than that
  found in the aqueous model of Hofmann and Schieberle". Ribose and/or cysteine added "typically in
  the range of **25-250 mM**". **5 mL portions in closed 15 mL glass containers, in the dark, at
  50 C, under air only.** Sampled at 0, 4, 8 and 24 h, **in duplicate**.
- **Quantification.** Internal standard **maltol** — 50 µg into the 100 mL process-flavouring
  sample, 500 µg into the 5 mL reconstituted sample (a 100-fold difference in the standard-to-sample
  ratio; Flags 6). Extraction with **dichloromethane/diethyl ether 7:3 v/v** (20 mL from the process
  flavouring, 2.0 mL from the reconstituted), dried over Na2SO4, concentrated at room temperature
  and 200 mbar to 2 mL then under argon to 1.0 mL for the process flavouring. **GC-FID**, HP-5
  50 m x 0.32 mm x 1.05 µm, injector 225 C, detector 250 C, 65 -> 120 C at 3 C/min then to 250 C at
  40 C/min. **Recoveries >= 70 %, RSD <= 10 %**, validated in triplicate. HRGC-MS on a Finnigan MAT
  TSQ 70, EI at 70 eV. **Cysteine and ribose by capillary electrophoresis** (HP 3D CE).
- **The kinetic model, in the authors' words.** "All plots of relative area against storage time
  showed a **linear decrease, indicating zero-order or pseudo-zero-order kinetics** for the compounds
  investigated. In all storage experiments described in this report, zero-order models were used on
  the basis of **visual assessment** and the coefficient of correlation (R) obtained from regression
  analysis. The relative decrease in the starting concentration of the analytes (% per day) was
  calculated by **linear regression, using the least-squares method**." **No R value is printed
  anywhere**, and no standard error on any rate (Flags 1).
- **The H-D exchange experiment.** MFT and, separately, bis(3-furanyl) disulfide at **0.01 mM each
  in CH3OD containing DCl**, at **room temperature**, followed by 1H NMR (JEOL 400 EX, 399.65 MHz,
  TMS reference) over several days — Figure 3's caption says **64 days**. **This is a methanolic,
  acidified, sub-micromolar system and not the storage matrix**, and the authors say why: MFT's low
  solubility in water. The disulfide was used in place of 3-furanthiol because pure 3-furanthiol
  could not be isolated.
- **The mass-balance experiment.** "storage of a mixture of the three thiols MB, FFT, and MFT in
  phosphate-buffered solution (0.5 M; pH 5.0) **in the absence of cysteine**", then a comparison of
  the initial thiol against the final thiol **plus (mixed) disulfides**. The DMFT repeat is marked
  "data not presented". **No numbers are printed for this experiment at all** — it is stated as an
  equality and an inequality (Flags 4).

## 3. Tables re-typed

### Table 1. "Decrease in Concentration (% per day) of Five Character-Impact Components of an Aqueous Model Process Flavor (Entry 1) and Reconstituted Models (Entries 2-5) During Storage at 50 °C: Influence of Ribose and Cysteine"

Column head exactly as printed: **"decrease in concentration (% per day)"**.

| compound | 1^a model process flavoring | 2^b no addition | 3^b + 250 mM ribose | 4^b + 150 mM cysteine | 5^b + 150 mM cysteine + 250 mM ribose |
|---|---:|---:|---:|---:|---:|
| 2-mercapto-3-butanone (MB) | 14 | > 90 | > 90 | 26 | 25 |
| 2-methyl-3-furanthiol (MFT) | **59** | > 90 | > 90 | 69 | 42 |
| 2-furfurylthiol (FFT) | **28** | > 90 | > 90 | 39 | 43 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone (HDF) | < 10 | < 10 | < 10 | 77 | 66 |
| 3-hydroxy-4,5-dimethyl-2(5H)-furanone (sotolone) | n.d.^c | < 10 | < 10 | 22 | < 10 |

Footnote a, exactly as printed: "According to the procedure of Hofmann and Schieberle (2), a
phosphate-buffered (0.5 M, pH 5.0) solution of D-ribose (100 mM) and L-cysteine (33 mM) was heated
in an autoclave from room temperature to 130 °C in 10 min and was kept at 130 °C for 20 min.
Character-impact components were formed in the order of concentrations of around 5 µM."
Footnote b: "50 µM of each character-impact component was added to an aqueous phosphate-buffered
solution (0.5 M; pH 5.0)." Footnote c: "Not determined."

**Read the columns carefully.** Column 1 is a **real Maillard matrix at ~5 µM thiol** with its own
residual cysteine (15 mM) and ribose (25 mM) and all its Maillard products. Columns 2-5 are **clean
buffer at 50 µM thiol** with nothing but what was added. They are not five levels of one experiment;
they are one experiment (column 1) and a four-level experiment (columns 2-5), and only 2-5 are
mutually comparable.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| abstract summary | MFT "59 % decrease/24 h", FFT "28 % decrease/24 h", 2-mercapto-3-butanone "14 % decrease/24 h", 2,5-dimethyl-4-hydroxy-3(2H)-furanone "max. 10 % decrease/24 h" | p. 4292 abstract |
| charge of the process flavouring | 180 mmol ribose + 60 mmol cysteine in 1800 mL | Methods |
| **residual cysteine in the fresh process flavouring** | **15 mM** | p. 4293 |
| **residual ribose in the fresh process flavouring** | **25 mM** | p. 4293 |
| air versus argon | "almost no difference ... (data not presented here), probably due to an anti-oxidative effect of the matrix" | p. 4293 |
| effect of ribose alone | "ribose addition has no effect on stability (compare entry 3 with entry 2)" | p. 4293 |
| effect of cysteine | "a stabilizing effect through the addition of cysteine was clearly present (entry 4)" | p. 4293 |
| effect of both | "an even larger stabilizing effect on MFT and the furanones than with cysteine alone" | p. 4293 |
| **MFT's cysteine optimum** | "a clear **minimum** at a cysteine concentration of **around 50 mM**, which is different from what is observed for MB and FFT" | p. 4294, Figure 1 |
| explanation offered for the MFT optimum | ribose lowers the free cysteine by forming a **thiazolidinecarboxylic acid** | p. 4294 |
| **the mass balance** | MB and FFT: initial thiol = final thiol + (mixed) disulfides. MFT: "absolutely not the case ... a large amount of MFT was not detected, neither as thiol nor as (mixed) disulfide" | p. 4294 |
| the DMFT control | DMFT's balance closes ("data not presented") | p. 4294 |
| **H-D exchange in MFT** | **85 %** at the **5-position**, **10 %** at the **4-position**, after several days in CH3OD/DCl at room temperature (Figure 3 caption: 64 days) | p. 4295 |
| H-D exchange in bis(3-furanyl) disulfide | "high exchange rate at the 2- and 2'-positions ... remarkable when compared to that of the other positions" — **no percentages printed** | p. 4295, Figure 5 |
| dimers not observed in the NMR experiment | "except for the formation of the corresponding disulfide, no evidence for the formation of dimers was found, as was expected according to the mechanism in Figure 2" | p. 4295 |
| the coupling-to-other-thiols branch is ruled out | "Because no loss of the other thiols was observed in the reconstituted mixture experiments, we must conclude that this part of the proposed mechanism is not significant." | p. 4294 |
| method validation | recoveries >= 70 %, RSD <= 10 %, triplicate | Methods |

**Figure-only in this paper:** every point of Figure 1 (% per day against cysteine concentration,
three thiols) except the qualitative "minimum around 50 mM"; every point of Figures 3 and 5 (H-D
exchange against time and position) except the printed 85 % and 10 %. Per house rule they are not
typed as numbers.

### Arithmetic on the printed numbers (all mine)

**1. Absolute loss rates, and the order test.** Zero order means a constant amount per unit time.
Converting the printed percentages onto the stated pool sizes:

| pot | pool | MFT | FFT | MB |
|---|---:|---:|---:|---:|
| entry 1, process flavouring | ~5 µM | 2.95 µM/day | 1.4 µM/day | 0.70 µM/day |
| entry 2, clean buffer | 50 µM | **> 45 µM/day** | > 45 µM/day | > 45 µM/day |
| entry 4, + 150 mM cysteine | 50 µM | 34.5 µM/day | 19.5 µM/day | 13 µM/day |
| entry 5, + cysteine + ribose | 50 µM | 21 µM/day | 21.5 µM/day | 12.5 µM/day |

**Entry 1 against entry 2 is a 10x change in the pool and a >= 15x change in the absolute rate.**
Zero order across pots is therefore refused by the paper's own numbers. Two readings survive and the
data cannot separate them: (i) the reaction is nearer first order (or, for a self-coupling, second
order) in the thiol and the "zero order" is only the near-linearity of an exponential over its first
half-life; (ii) the intact Maillard matrix protects the thiol, over and above its 15 mM of residual
cysteine, by a further factor of about **11.7** (entry 1's 2.95 µM/day against entry 4's 34.5, both
per unit pool: 0.59/day against 0.69/day is only 1.17x, so on a *fractional* basis entries 1 and 4
nearly agree and the whole discrepancy is the pool size — which is reading (i)). **On the fractional
basis the process flavouring and the 150 mM-cysteine reconstitution agree to within 17 % for MFT,
39 % for FFT and 86 % for MB.** That is the sharpest thing in this dossier: **the fractional rate,
not the absolute rate, is what transfers**, which is first-order behaviour and not zero-order
behaviour, and the module's `THIOL_CHANNELS` entry describing this channel as "ZERO in thiol"
inherits an artefact of how the authors chose to report.
*Caveat carried:* entries 1 and 4 differ in matrix as well as in pool, so this is suggestive, not
decisive. It is exactly the kind of claim a hold-out is supposed to test rather than absorb.

**2. Half-lives implied at 50 C, if first order (mine).** From entry 1: MFT t½ = ln2/(-ln(1-0.59))
days = 0.78 day = **19 h**; FFT **1.9 days**; MB **4.6 days**. From entry 2, all three < 0.3 day.
These are *conversions of the printed percentages under an assumed order*, not measurements.

**3. Precursor conversion in the 130 C cook (mine) — the FIT row.** Ribose 100 mM charged, 25 mM
residual: **75 % converted**. Cysteine 33 mM charged, 15 mM residual: **54.5 % converted**. Both
after 10 min ramp + 20 min hold at 130 C in 0.5 M phosphate at pH 5.0. This is the reactant-side
constraint D.3 declares FIT, and it verifies to the two figures the declaration quotes.

**4. What the 5 µM tells us about the 130 C cook's yield (mine).** ~5 µM of each character-impact
compound from 33 mM of cysteine is a molar yield of about **0.015 %** per compound on the cysteine,
or about **0.005 %** on the ribose. The figure is loose ("in the order of ... around 5 µM", one
number for five compounds) but it is the only yield this paper implies and it is the right order for
the Hofmann 1995 pot it reproduces.

**5. HDF and sotolone move the opposite way from the thiols on cysteine addition.** HDF is
< 10 %/day in entries 1, 2 and 3 and **77 %/day in entry 4** — cysteine *destabilises* the furanone
by a factor of at least 7.7 while *stabilising* every thiol. Sotolone does the same (< 10 -> 22).
This is a clean within-study crossover on one additive and it is the only place in this paper where
a non-thiol is destroyed; the lane's furanone chemistry has no such term.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Keyed: `2_methyl_3_furanthiol`,
`2_furfurylthiol`, `hdmf` (= HDF, 4-hydroxy-2,5-dimethyl-3(2H)-furanone). **Not keyed:**
2-mercapto-3-butanone (MB), 3-hydroxy-4,5-dimethyl-2(5H)-furanone (sotolone),
2,5-dimethyl-3-furanthiol (DMFT), bis(3-furanyl) disulfide, maltol, thiazolidinecarboxylic acid, and
both reactants. See Flags 8. Note that **sotolon is explicitly out of scope** for the sulfur wave per
`sulfur.py`'s OUT_OF_SCOPE block.

**EVERY ROW BELOW EXCEPT THE LAST TWO IS PART OF A DECLARED ★ HOLD-OUT** (`FIT_HOLDOUT_DECLARATION.md`
D.4, "van Seeventer 2001 Table 1, 50 °C zero-order MFT/FFT"). They are recorded here so that the
hold-out can be *scored*, not fitted. Shared conditions: **0.5 M phosphate, pH 5.0, dark, 50 C,
0-24 h, GC-FID against maltol, recoveries >= 70 %, RSD <= 10 %.**

| step / quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| MFT loss, process flavouring | **59** | % per day | ~5 µM MFT in the 130 C ribose/cysteine pot, air or argon | **zero order** by the authors, by visual assessment | Table 1 entry 1 | **measured_rate** (★ HOLD-OUT) |
| FFT loss, process flavouring | **28** | % per day | as above | zero order | Table 1 entry 1 | measured_rate (★ HOLD-OUT) |
| MB loss, process flavouring | **14** | % per day | as above | zero order | Table 1 entry 1 | measured_rate (★ HOLD-OUT) |
| HDF loss, process flavouring | < 10 | % per day | as above | zero order | Table 1 entry 1 | **threshold** |
| sotolone, process flavouring | not determined | — | — | — | Table 1 entry 1 | — |
| MFT / FFT / MB loss, clean buffer | **> 90 / > 90 / > 90** | % per day | 50 µM each, no additive, air | zero order | Table 1 entry 2 | **threshold** (lower bounds) |
| the same + 250 mM ribose | > 90 / > 90 / > 90 | % per day | 50 µM each | zero order | Table 1 entry 3 | threshold |
| MFT / FFT / MB + 150 mM cysteine | **69 / 39 / 26** | % per day | 50 µM each | zero order | Table 1 entry 4 | measured_rate (★ HOLD-OUT) |
| MFT / FFT / MB + 150 mM cysteine + 250 mM ribose | **42 / 43 / 25** | % per day | 50 µM each | zero order | Table 1 entry 5 | measured_rate (★ HOLD-OUT) |
| HDF + 150 mM cysteine | **77** (against < 10 with no additive) | % per day | 50 µM | zero order | Table 1 entries 2 and 4 | measured_rate — the furanone crossover |
| sotolone + 150 mM cysteine | 22 (against < 10) | % per day | 50 µM | zero order | Table 1 entries 2 and 4 | measured_rate |
| **MFT mass balance** | **does NOT close**: after storage MFT is absent "neither as thiol nor as (mixed) disulfide" | — | 50 µM three-thiol mix, no cysteine, pH 5.0, 50 C | — | p. 4294 | **level_only** (a stated inequality with no number attached) |
| MB and FFT mass balance | **closes** as thiol + (mixed) disulfides | — | same experiment | — | p. 4294 | level_only |
| DMFT mass balance | **closes** ("data not presented") | — | same experiment with DMFT for MFT | — | p. 4294 | level_only |
| air versus argon | "almost no difference" — **data not presented** | — | process flavouring, 50 C, 24 h | — | p. 4293 | **level_only** — a stated null with no data behind it |
| **residual cysteine after the 130 C cook** | **15** | mM (from 33 mM charged) | 0.5 M phosphate pH 5.0, 130 C, 10 min ramp + 20 min hold | — | p. 4293 | **measured_rate** basis for the FIT conversion row |
| **residual ribose after the 130 C cook** | **25** | mM (from 100 mM charged) | as above | — | p. 4293 | as above |
| cysteine conversion | **54.5** | % | as above | — | derived (mine) | derived_assumption — **this is D.3's FIT row** |
| ribose conversion | **75** | % | as above | — | derived (mine) | derived_assumption — **D.3's FIT row** |
| character-impact concentration in the pot | ~5 | µM each | after the 130 C cook | — | Table 1 footnote a | level_only (one figure for five compounds) |
| implied molar yield on cysteine | ~0.015 | % | as above | — | derived (mine) | derived_assumption |
| H-D exchange in MFT | **85 at C-5, 10 at C-4** | % of H exchanged | 0.01 mM in CH3OD/DCl, room temperature, 64 days | — | p. 4295 | **measured_rate**-adjacent, but in a **different solvent and matrix** — carry as level_only for the storage system |
| H-D exchange in bis(3-furanyl) disulfide | high at 2 and 2', low elsewhere; **no numbers** | — | same solvent | — | p. 4295, Fig. 5 | **figure_only** |
| loss rate against cysteine concentration, three thiols | MFT has a **minimum near 50 mM**; MB and FFT do not | — | 50 µM each, 25-250 mM cysteine, 50 C | — | Figure 1 | **figure_only** |
| absolute loss, entry 1 vs entry 2 | >= 15x more per day on a 10x larger pool | — | — | — | derived (mine) | **within_study_ratio (mine)** — the order test, section 3 item 1 |

### Can any of this be put on the same basis as a shipped constant? Deliberately, no.

The module's position is that it should not be, and this dossier agrees for the reason the
declaration gives: **the functional form is the thing being tested.** Three further obstacles, each
of which would have to be cleared before a rate from this paper could be transported:

- **The unit is % per day of a pool whose size is stated only to one significant figure** ("around
  5 µM") in the one entry that is a real Maillard matrix. A rate constant of any order requires the
  pool; entry 1's pool is approximate and is quoted once for five different compounds.
- **There is no barrier and no second temperature.** 50 C is the only storage temperature in the
  paper. Pairing it with the 145 C panel to extract an activation energy is already a named
  `PROHIBITED_DERIVATION` and this paper gives no reason to reopen that — indeed it gives a reason to
  keep it shut, since its own sink is expressly not the one operating at 115-120 C (the disulfide) or
  at 25-30 C (the covalent thioether to melanoidin).
- **The mechanism experiment is in methanol.** The 85 %/10 % H-D exchange, which is the only direct
  evidence for the proposed C-5 electrophilic route, was run in CH3OD/DCl at 0.01 mM because MFT is
  too insoluble in water — a different solvent, a different acidity, a 5000-fold lower
  concentration and 64 days instead of one. It supports the *regiochemistry*; it does not measure
  the storage reaction.

## 5. Flags

1. **Not one error bar, correlation coefficient or replicate value is printed.** The Methods say
   zero-order models were chosen "on the basis of visual assessment and the coefficient of
   correlation (R) obtained from regression analysis" and **no R appears in the paper**. The
   validation figures (recovery >= 70 %, RSD <= 10 %) are for the analytical method, not for the
   fitted rates. Every percentage in Table 1 is a point estimate of unstated precision, and three of
   the fifteen thiol cells are "> 90", which is a censored observation and not a rate.
2. **The paper names the same compound two ways.** The abstract and Table 1 print
   **2-mercapto-3-butanone**; the Results text on p. 4293-4294 prints **3-mercapto-2-butanone** for
   the same abbreviation MB. These are different molecules. The introduction's list of the
   ribose/cysteine odorants gives "2-mercapto-3-butanone (MB)", and the Materials section bought
   "2-mercapto-3-butanone (10 % solution in triacetin)" — so **the purchased compound is
   2-mercapto-3-butanone** and the Results text is the typo. This matters because 3-mercapto-2-
   butanone is the isomer Cerny 2003 and Whitfield 1999 work with, and the two literatures could be
   crossed by this slip.
3. **The air-versus-argon comparison is asserted, not shown.** "there was almost no difference
   between storage of the aqueous model process flavoring under air or argon atmosphere (**data not
   presented here**)". The reconstituted experiments (entries 2-5), which carry all the mechanistic
   contrasts, were run **under air only**. The module's channel entry states "air ~ argon (so it is
   NOT oxidative)" as if it were a measurement; it is the authors' summary of unpublished data. The
   *mass-balance* argument for a non-oxidative sink is independent of it and is stronger.
4. **The mass-balance experiment has no numbers at all.** It is the most important result in the
   paper for this repository and it is reported as two sentences and one "data not presented". No
   initial concentration, no final thiol, no final disulfide, no recovery, no replicate count, no
   statement of how much "a large amount" is. **Nothing quantitative can be built on it**; what it
   licenses is a structural claim, and that is how this dossier records it.
5. **The proposed mechanism is not confirmed by the experiment run to confirm it.** Figure 2
   predicts MFT dimers/oligomers; the NMR experiment found "no evidence for the formation of dimers
   ... except for the formation of the corresponding disulfide", and the authors explain the failure
   by the change of solvent. The oligomer (`OLG`) species the module carries is therefore named after
   a **hypothesis whose direct test was negative in the only system where it was tried.** The
   evidence for the sink is the missing mass, not the found product.
6. **The internal-standard ratio differs 100-fold between the two experiment types** (50 µg maltol
   into 100 mL versus 500 µg into 5 mL) and the extraction ratio differs too (20 mL solvent for
   100 mL sample versus 2.0 mL for 5 mL). Both are self-consistent, but entry 1 and entries 2-5 are
   different analytical methods as well as different chemistries, which is a second reason not to
   read entry 1 against entry 2 as a concentration series. My order test in section 3 carries this
   caveat.
7. **What to request from the authors**: (i) the air-versus-argon data; (ii) the mass-balance
   numbers — initial MFT, final MFT, final disulfides, and the same for MB, FFT and DMFT; (iii) the
   numeric points behind Figure 1, especially the position and depth of the MFT minimum; (iv) the
   regression outputs (R, intercept, and the 0/4/8/24 h points) that would let the zero-order claim
   be re-tested as a first- or second-order fit — **this is the single most valuable request in the
   cluster**, because the functional form is the thing the module holds this row out to test;
   (v) whether the missing MFT can be found in the non-volatile residue.
8. **Registry gaps against `data/keys/compounds.yml`**: **2-mercapto-3-butanone** and
   **2,5-dimethyl-3-furanthiol** are the two that matter — MB is one of the three thiols carrying
   the hold-out's rates and DMFT is the negative control that makes the MFT-specific sink
   interpretable. Also absent: sotolone (out of scope by declaration), bis(3-furanyl) disulfide,
   maltol (the internal standard).
9. **What this paper does not contain**: any second storage temperature; any pH series; any
   activation energy; any rate constant with a unit of inverse time; any identification or
   quantification of the MFT degradation products; any measurement of the disulfides that did form;
   any headspace or partition data; any oxygen concentration; any light exposure (everything is
   dark); any water-activity variation; and any supplementary material.
