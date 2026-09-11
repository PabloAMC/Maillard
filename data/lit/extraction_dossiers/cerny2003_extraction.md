# Cerny & Davidek 2003 — EXTRACTION (ribose + cysteine, 0.5 mol/L potassium phosphate pH 5.00, 95 C, 4 h, sealed 2 mL vials; five CAMOLA isotope pots read by HS-SPME-GC-MS; fifteen compounds' isotopomer shares in per cent, and NOT ONE rate, yield or concentration)

### THE PAPER THE SULFUR LANE LEANS ON MOST IS A LABELLING PAPER AND NOTHING ELSE: it fixes the topology of `r_ddp_mft`, `r_nf_mp3p` and `r_fur_fft` with isotopomer shares that are ratios by construction, and it contains **no thiol loss rate, no binding plateau, no mass balance and no measured sink partner** — so on the question that refused three sink structures it is silent, with one exception that is a citation rather than a measurement: it quotes van Seeventer (ref 41) for the negative that disulfide formation is *not* the cause of 2-methyl-3-furanthiol instability at 50 C, which is the same negative the B17 variant (b) refusal reached from the other side.

**Source on disk:** `data/articles/cerny2003.pdf` (8 pp., J. Agric. Food Chem. 2003, 51 (9), 2714-2721).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/cerny2003.txt`, 455 lines). **Tables 1, 2 and 3 came through clean** and are
re-typed in full below; the column alignment of Table 2 is intact and was checked row by row against
the running text (the text quotes 2-methyl-3-furanthiol as "approximately 1:1", 2-methyl-3-
(methylthio)furan as "3:1:3:1", thiazole as "56 % unlabeled / 44 % singly labeled" and the
disulfides as "approximately 1:2:1", and the table reproduces every one). Figures 1-6 are structure
drawings and proposed mechanisms — **no figure in this paper carries a datum**, so nothing here is
figure-read. There is no supplementary material. Repo status before this dossier: **42 files in the
repository match "Cerny 2003"** — `src/kinetic_core/sulfur.py` (the `r_ddp_mft` note),
`src/kinetic_core/species_sulfur.py` (the DDP and MP3P notes),
`src/kinetic_core/parameters_sulfur.py` (line 1928), `docs/reference/FIT_HOLDOUT_DECLARATION.md`
(where T2/T3 are a declared **HOLD-OUT**), every B2.x fit and hold-out generator and report, the B8
panels, two isolate benchmarks, `k3_final_parameter_inventory.md` and `cerny2004_extraction.md` —
and **it has never had a dossier of its own.** Its numbers have reached the code only through the
inventory and through the one-line notes quoted above.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction" |
| Authors | Christoph Cerny (corresponding; by publication at Firmenich SA, Meyrin-Geneva) and Tomas Davidek — Nestle Research Center, Vers-chez-les-Blanc, Lausanne, and Centres R&D Nestle S.A.S., Amiens |
| Venue | J. Agric. Food Chem. 2003, 51 (9), 2714-2721. Received 13 November 2002, revised 21 February 2003, accepted 2 March 2003, web 29 March 2003 |
| DOI / article ID | 10.1021/jf026123f (printed as `JF026123F`) |
| Method acronym | CAMOLA — "carbohydrate module labeling" (ref 25) / "carbon module labeling" (ref 26); a 1+1 mixture of unlabelled and [13C5] precursor, read as isotopomer shares of the molecular ion |
| Naming | "1,4-dideoxyosone" of ribose = 5-hydroxy-2,3-pentanedione = laurencione; "4-hydroxy-5-methyl-3(2H)-furanone" = norfuraneol (the repository's `NF`); compound numbers 1-29 refer to Figures 1-6 |
| Companions on disk | `cerny2004_extraction.md` (in-situ branching), `cerny2007_extraction.md` (the pH ladder and the thiamine/xylose split), `whitfield2001_extraction.md` and this cluster's `whitfield1999_extraction.md` (the norfuraneol pots this paper's Table 3 argues against), `hofmann1998b_extraction.md` and this cluster's `hofmann1995_extraction.md` (the same ribose/cysteine system at 145 C) |

## 1. Why it matters

**What it contributes to the thiol-sink question: nothing quantitative, and this dossier says so
plainly.** The paper prints no rate constant, no half-life, no yield in mol %, no concentration, no
mass balance and no time series. Every number in it is a **share of a molecular-ion signal among
that same molecule's isotopomers** — a within-study ratio in the strictest sense, and one that is
blind to how much of anything was made or lost. Three sink structures were refused on an objective
that never saw this paper; **had it been read first, it would not have changed the objective by one
row**, because it has no row to give. The one sentence in it that touches the sink is a *citation*,
not a measurement (p. 2717): "According to van Seeventer (41) disulfide formation is not the cause
of 2-methyl-3-furanthiol instability in aqueous solutions of reacted cysteine/ribose kept at 50 C."
That is the third paper of this cluster speaking, and it is the same negative the B17 variant (b)
refusal reached from the model side. What Cerny adds to it, from his own data, is structural and is
in section 4: his disulfide isotopomer shares are **binomial on a well-mixed free-thiol pool**
(26 : 48 : 22 against the 25 : 50 : 25 of random pairing, Table 2 rows 13 and 15), so whatever
disulfide exists at 95 C comes from the *free* thiol after it is formed, not from a bound or
sequestered pool — which is a constraint on where a sink may sit, not a number for one.

**What it does contribute, and what the code already takes from it.** The repository's sulfur
network hard-codes three topology claims out of this paper:

| code site | the claim as the code states it | this paper's source |
|---|---|---|
| `sulfur.py` `r_ddp_mft` (DDP + H2S -> MFT) | "Cerny 2003 T2: 49 % unlabelled / 46 % 13C5 with no fragment pattern, 'pathways via ribose fragmentation were not relevant' => ~93 % of MFT carries the intact pentose skeleton at 95 C" | Table 2 row 7 (49 / 1 / 0 / 1 / 3 / 46) and p. 2715 |
| `sulfur.py` `r_nf_mp3p` (NF + H2S -> MP3P) | "2-mercapto-3-pentanone is 96 % unlabelled from Cerny's NF spike — the sharpest single-species NF-route marker in the corpus" | Table 3 row 16 (96 % unlabelled / 4 % 13C5) |
| `species_sulfur.py` DDP note | "1,4-dideoxyosone + RCHO + CO2 + NH3. Cerny 2003 T2 makes this the ..." | Figures 2 and 6, the proposed route; **the route itself is proposed, not demonstrated** (Flags 2) |
| `sulfur.py` `r_fur_fft` | furfural as the dominant FFT precursor | this paper's fifth pot: FFT from a 2-furaldehyde + [13C5]ribose + cysteine pot is **92 : 8 unlabelled : 13C5** at equimolar loading, p. 2718 |

The intact-skeleton claim is the load-bearing one and it holds: **MFT, FFT and 3-mercapto-2-
pentanone are each essentially only unlabelled or fully 13C5**, with the intermediate mass channels
(1 to 4 labels) summing to 5, 4 and 4 % respectively. The paper's own words (p. 2715): "pathways via
ribose fragmentation were not relevant."

Two further findings the lane has *not* taken up. (a) **Norfuraneol is demoted as an MFT
intermediate**: in the pot where norfuraneol competes against [13C5]ribose head to head, MFT comes
out **93 % labelled**, i.e. from the ribose, "only to a small extent from 4-hydroxy-5-methyl-3(2H)-
furanone" (Table 3, p. 2718). The `FIT_HOLDOUT_DECLARATION` records this as "the NF <= 7 % ceiling"
and holds it out. (b) **2-Furaldehyde is a far more efficient FFT precursor than ribose at equimolar
loading** (92 : 8), which is the qualitative parent of the "60x ribose" figure the `r_fur_fft` note
carries from elsewhere — this paper supports the *ordering* and does not print a factor.

What this paper does NOT give the repository: any rate, any barrier, any absolute or relative yield,
any time course, any temperature other than 95 C, any pH other than 5.00, any headspace-to-solution
partition, and any statement of how much of the fed precursor reacted.

## 2. Methods as they matter to a model

- **Pots.** Five, all in the same buffer and regime (Table 1): **potassium phosphate 0.5 mol/L,
  pH 5.00, 95 C, 4 h**, reactants dissolved "to a total amount of 500 mg", in **silanized 2 mL glass
  vials, septum-closed**, heated in a Reacti-Therm metal block. The charges in µmol are re-typed in
  section 3. The molar ratio cysteine : other precursor is **1 : 3** in every pot except the
  cysteine/norfuraneol pot, which is **2 : 3**.
- **The 1+1 labelling trick.** Where a pot carries both `ribose` and `[13C5]ribose`, they are
  **72.5 µmol each** — a nominal 1 : 1 mixture. If the five-carbon skeleton survives, the product
  can only be all-12C or all-13C5 and the two must appear in equal amounts; any intermediate mass
  channel is the signature of fragment recombination. This is the whole logic of the paper and it is
  why the shares are *shares*, never yields.
- **[13C5]ribose enrichment is 98 %** (Cambridge Isotope Laboratories), not 100. The residual 2 %
  propagates into every "13C5" share as a small deficit and into the intermediate channels as a
  small excess; the paper does not correct for it (Flags 4).
- **Analysis.** HS-SPME-GC-MS, **duplicate**. At least 1 h equilibration at 20 C, then a
  PDMS-DVB fibre (65 µm) exposed **60 min at 20 C** to the headspace, **without agitation**;
  5 min desorption at 250 C through a 0.75 mm i.d. liner. GC 6890A / MSD 5973, HP-PONA
  50 m x 0.20 mm x 0.50 µm, 35 -> 240 C at 6 C/min then 10 min isothermal. **EI at 70 eV, scan
  m/z 28-350.** No internal standard, no isotope-dilution quantification, no calibration —
  **this instrument configuration cannot produce a concentration and the paper never claims one.**
- **How the shares are corrected.** Table 2 footnote a: the M+ + 1 and M+ + 2 signals are corrected
  by subtracting the natural abundances of **13C (1.10 %), 33S (0.76 %) and 34S (4.20 %)**, and "the
  loss of hydrogen frequently observed with the molecular ion in EI-MS was also corrected in the
  labeled molecular ions by the ratio (M+ - 1)/M+". So the shares are corrected ion intensities,
  compound by compound, and are internally comparable within a row and **not across rows**.
- **Identification.** By mass spectra and retention indices against authentic compounds "analyzed in
  our laboratory or in literature data" (both table footnotes b).
- **Regime.** The authors place the conditions inside the Council of Europe process-flavour
  guidelines (ref 27: <= 180 C, <= 24 h for temperatures <= 110 C, pH <= 8.0).
- **Sampling of the headspace, not the pot.** Every share is measured on what the fibre absorbed
  from the headspace of a sealed vial at 20 C. Isotopomers of one compound share a partition
  coefficient to a very good approximation, so **a share is safe where a level would not be** — this
  is the methodological reason the paper's numbers are ratios and the reason they may be used as
  such.

## 3. Tables re-typed

### Table 1. "Model Reactions" — footnote a: "Reaction in phosphate buffer (0.5 mol/L; pH 5.00) at 95 C (4 h)."

Header exactly as printed: **amount (µmol)**, columns A-E.

| reactant | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|
| cysteine | 48.3 | 48.3 | 48.3 | 48.3 | 48.3 |
| ribose | 144.9 | | 72.5 | | |
| [13C5]ribose | | | 72.5 | 72.5 | 72.5 |
| 4-hydroxy-5-methyl-3(2H)-furanone | | 72.5 | | 72.5 | |
| 2-furaldehyde | | | | | 72.5 |

(Blank cells are blank in the printed table: the reactant is absent from that pot.)

### Table 2. "Proportion of Isotopomers from the Reaction between Ribose, [13C5]Ribose, and Cysteine" (pot C)

Header exactly as printed: **"proportion of labeled carbon atoms in the molecule^a (%)"** over
columns **0^c, 1, 2, 3, 4, 5, 6, 10**; footnote c: "Number of 13C atoms in the molecule".

| no. | compound^b | m/z (M+) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 10 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | furan | 68 | 49 | 0 | 0 | 1 | 50 | | | |
| 2 | 2-methylfuran | 82 | 51 | 0 | 0 | 1 | 9 | 39 | | |
| 3 | thiazole | 85 | 56 | 44 | 0 | 0 | | | | |
| 4 | 2-methylthiophene | 98 | 45 | 7 | 0 | 0 | 3 | 44 | | |
| 5 | 3-mercapto-2-butanone | 104 | 49 | 0 | 0 | 1 | 50 | | | |
| 6 | 2-furaldehyde | 96 | 47 | 0 | 0 | 2 | 5 | 45 | | |
| 7 | 2-methyl-3-furanthiol | 114 | 49 | 1 | 0 | 1 | 3 | 46 | | |
| 8 | 3-mercapto-2-pentanone | 118 | 49 | 1 | 0 | 1 | 2 | 47 | | |
| 9 | 2-furfurylthiol | 114 | 48 | 0 | 1 | 1 | 2 | 47 | | |
| 10 | 2-methyl-3-(methylthio)furan | 128 | 36 | 11 | 1 | 0 | 2 | 36 | 13 | |
| 11 | 3-thiophenethiol | 116 | 95 | 1 | 2 | 1 | 1 | | | |
| 12 | 2-methyl-3-thiophenethiol | 130 | 44 | 2 | 6 | 14 | 10 | 24 | | |
| 13 | bis(2-methyl-3-furyl) disulfide | 226 | 26 | 0 | 0 | 0 | 3 | 48 | 0 | 22 |
| 14 | (2-methyl-3-furyl) (2-oxo-3-pentyl) disulfide^d | 230 | 27 | 0 | 1 | 0 | 3 | 49 | 0 | 20 |
| 15 | bis(2-furfuryl) disulfide | 226 | 28 | 0 | 0 | 0 | 0 | 54 | 0 | 18 |

Footnote a as quoted in section 2. Footnote b: identification by MS and retention index against
authentic compounds. Footnote d: "Mass spectrum of the unlabeled compound m/z (%) 230 (M+, 43),
187 (13), 145 (27), 114 (34), 113 (92), 85 (15), 81 (14), 43 (100)."

**Row sums (mine).** Rows 1-9 and 11-12 sum to 100 within rounding except row 2 (100), row 4 (99),
row 6 (99) and row 11 (100). Row 10 sums to 99; rows 13, 14 and 15 sum to 99, 100 and 100. The
blank cells are blanks in the print, i.e. that isotopomer channel does not exist for that molecule
(furan has 4 carbons, thiazole 3, and so on), **not zeros that were omitted**.

### Table 3. "Proportion of Labeling in Compounds from the Reaction between 4-Hydroxy-5-methyl-3(2H)-furanone, [13C5]Ribose, and Cysteine" (pot D)

Header exactly as printed: **"no. of 13C atoms in the labeled molecule"**, **"unlabeled compound
(%)"**, **"13C-labeled compound (%)"**.

| no. | compound^a | m/z (M+) | no. of 13C atoms in the labeled molecule | unlabeled compound (%) | 13C-labeled compound (%) |
|---|---|---:|---|---:|---|
| 5 | 3-mercapto-2-butanone | 104 | 5 | 94 | 6 |
| 6 | 2-furaldehyde | 96 | 5 | 2 | 98 |
| 7 | 2-methyl-3-furanthiol | 114 | 5 | 7 | 93 |
| 8 | 3-mercapto-2-pentanone | 118 | 5 | 42 | 58 |
| 16 | 2-mercapto-3-pentanone | 118 | 5 | 96 | 4 |
| 9 | 2-furfurylthiol | 114 | 5 | 0 | 100 |
| 10 | 2-methyl-3-(methylthio)furan | 128 | 5, 6 | 6 | 84, 10 |
| 11 | 3-thiophenethiol | 116 | 4 | 98 | 2 |
| 13 | bis(2-methyl-3-furyl) disulfide | 226 | 5, 10 | 1 | 15, 84 |

In this pot the **unlabelled** carbon comes from norfuraneol (fed unlabelled) and the **13C5** carbon
from [13C5]ribose. Read the "unlabeled (%)" column as the **norfuraneol share** of that compound and
the "13C-labeled (%)" column as the **ribose share**.

### Numbers printed in the running text (everything else in this paper is a structure drawing)

| quantity | value | where |
|---|---|---|
| MFT isotopomer ratio, pot C | unlabelled : 5-times-labelled "approximately 1:1"; m/z 115-118 "practically not observed" | p. 2715 |
| FFT and 2-furaldehyde, pot C | "only unlabeled or completely labeled compound in a ratio of 1:1" | p. 2715 |
| thiazole, pot C | "56 % unlabeled (m/z 85) and 44 % singly labeled (m/z 86)"; m/z 58 shows loss of H-[12C]CN and H-[13C]CN, label at **C-2** | p. 2716 |
| 2-methyl-3-(methylthio)furan, pot C | unlabelled : singly : 5-times : 6-times = **3 : 1 : 3 : 1**; the 2-methylfuran moiety is "exclusively from ribose", **half the thiomethyl carbon from cysteine** | p. 2716 |
| the three disulfides, pot C | "isotopomer ratio of approximately 1:2:1 (unlabeled/[13C5]/[13C10])", read as oxidation of the corresponding thiols | p. 2716 |
| 2-methyl-3-thiophenethiol | ratio **not determinable** — the peak coeluted with another compound | p. 2716 |
| 2-mercapto-3-pentanone, pot C | **not detected at all** (the paper calls this "surprisingly") | p. 2715 |
| FFT from the 2-furaldehyde pot (pot E) | unlabelled : 13C5 = **92 : 8** (m/z 114 vs 119), at **equimolar** 2-furaldehyde and ribose | p. 2718 |
| pot E, why so few compounds | "only a few compounds were detected due to saturation of the PDMS-DVB fiber with unreacted 2-furaldehyde" | p. 2718 |
| the sink sentence | "According to van Seeventer (41) disulfide formation is not the cause of 2-methyl-3-furanthiol instability in aqueous solutions of reacted cysteine/ribose kept at 50 C" | p. 2717 |
| the disulfide caveat | "To verify whether the disulfides represent true reaction products or are formed as artifacts during SPME, additional experiments would be necessary." | p. 2717 |
| literature comparison, mercaptoketones | ribose + cysteine at **140 C / 30 min** (ref 30) and **up to 145 C / 20 min** (ref 8 = Hofmann 1995) both gave 3-mercapto-2-pentanone **and** 2-mercapto-3-pentanone, "with 3-mercapto-2-pentanone dominating" | p. 2715 |

**No concentration, no yield, no rate appears anywhere in this paper.** There is no figure with an
axis.

### Arithmetic on the printed shares (all mine)

**1. The intact-skeleton fraction.** For a 1 : 1 unlabelled/[13C5] ribose charge, an intact skeleton
gives share(0) + share(5) = 100 and every intermediate channel 0. Summing Table 2:

| compound | share(0) + share(max) | intermediate channels (1-4) | reading |
|---|---:|---:|---|
| 2-methyl-3-furanthiol (7) | 49 + 46 = **95** | 1 + 0 + 1 + 3 = **5** | intact; the code's "~93 %" is conservative against this |
| 3-mercapto-2-pentanone (8) | 49 + 47 = **96** | **4** | intact |
| 2-furfurylthiol (9) | 48 + 47 = **95** | **4** | intact |
| 2-furaldehyde (6) | 47 + 45 = **92** | **7** | intact (the 4-label channel at 5 % is the largest of the three) |
| 3-mercapto-2-butanone (5) | 49 + 50 = **99**, but the maximum is **4 labels** | 1 | intact minus **one carbon**: a C4 product from a C5 chain, no recombination |
| furan (1) | 49 + 50 = **99** at 4 labels | 1 | as above |
| 2-methylfuran (2) | 51 + 39 = 90 at 5 labels, **plus 9 % at 4 labels** | 10 | mostly intact, with a real one-carbon-loss channel |
| 2-methylthiophene (4) | 45 + 44 = 89 | **10**, of which 7 at one label | mostly intact |
| 2-methyl-3-thiophenethiol (12) | 44 + 24 = 68 | **32** | **fragment-derived in large part — and the paper says the peak is coeluted, so do not use this row** |
| 3-thiophenethiol (11) | 95 unlabelled, max 4 labels | 5 | **from cysteine, not ribose** |

**2. The disulfides are binomial on the free thiol pool (mine).** Random pairing of two equally
abundant thiol isotopomers gives 25 : 50 : 25 for (0 : 5 : 10) labels. Table 2 measures
**26 : 48 : 22** for bis(2-methyl-3-furyl) disulfide (row 13) and **28 : 54 : 18** for
bis(2-furfuryl) disulfide (row 15), and 27 : 49 : 20 for the mixed disulfide (row 14). The
agreement with the binomial is close for row 13 (chi-square on 100 arbitrary units is small) and
poorer for row 15, where the homodimer channels are 6 to 7 points light and the 13C5 channel is
4 points heavy. **What this licenses: the disulfide draws from a well-mixed pool of free thiol
after the thiol is formed.** What it does not licence: any statement about how much disulfide there
is. There is no quantity here — only a composition.

**3. The norfuraneol ceiling on MFT (mine, from Table 3).** In pot D, norfuraneol (72.5 µmol,
unlabelled) and [13C5]ribose (72.5 µmol) compete at **equimolar loading** for the same cysteine
(48.3 µmol). MFT comes out 7 % unlabelled. Reading the shares as the two routes' contributions at
equal precursor charge, **the norfuraneol route supplies 7 % of the MFT and the ribose route 93 %,
i.e. a route ratio of 13.3 : 1 in ribose's favour at 95 C, pH 5** (mine). The same arithmetic on the
other rows gives: 2-mercapto-3-pentanone **24 : 1 in norfuraneol's favour**; 3-mercapto-2-butanone
**15.7 : 1 in norfuraneol's favour**; 3-mercapto-2-pentanone **1.38 : 1 in ribose's favour** (58/42);
FFT **norfuraneol contributes 0**; 3-thiophenethiol 98 % unlabelled but its carbon is cysteine's
(Table 2 row 11), so that row measures nothing about the two sugars. These are the sharpest
route-partition numbers in the paper and they are all within-study ratios at one condition.

**4. Where the paper contradicts its own Table 2 across pots (mine).** 2-Mercapto-3-pentanone is
**not detected** in pot C (ribose + cysteine alone) but is a named product of pot D and pot B
(norfuraneol + cysteine). Both facts are printed. Together they say: **norfuraneol makes
2-mercapto-3-pentanone and ribose does not**, at 95 C and pH 5. That is a qualitative branch test
the lane's `r_nf_mp3p` note already leans on, and it is stronger than the 96 % figure alone because
it is a presence/absence contrast in the same laboratory and the same week.

**5. What the 98 % enrichment does to the shares (mine).** With 98 % enrichment, a five-carbon
product of "labelled" origin has (0.98)^5 = **0.904** probability of being fully 13C5, so about
**10 % of the labelled molecules land in the 1-4 channels by isotopic impurity alone**, split mostly
into the 4-label channel. Against a labelled share of ~47, that predicts roughly **4-5 points** of
intermediate signal with no fragmentation whatever. **The measured intermediate channels for MFT,
3-mercapto-2-pentanone and FFT are 5, 4 and 4 points.** They are therefore consistent with *zero*
fragmentation, and the code's "~93 % intact" is if anything an under-statement. This is the single
most useful piece of arithmetic in the dossier for the `r_ddp_mft` note; it is mine, and it assumes
the 98 % is per-carbon, which the vendor specification usually means but which the paper does not
state (Flags 4).

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Keyed: `2_methyl_3_furanthiol`,
`2_furfurylthiol`, `bis_2_methyl_3_furyl_disulfide`, `furan`, `2_methylthiophene`, `furfural`
(= 2-furaldehyde), `norfuraneol`, `hydrogen_sulfide`, `mercapto_2_propanone`. **Not keyed and
appearing in this paper's tables:** 3-mercapto-2-pentanone, 2-mercapto-3-pentanone,
3-mercapto-2-butanone, 2-methyl-3-(methylthio)furan, 3-thiophenethiol, 2-methyl-3-thiophenethiol,
bis(2-furfuryl) disulfide, the mixed (2-methyl-3-furyl)(2-oxo-3-pentyl) disulfide, thiazole,
2-methylfuran, and both reactants (ribose, cysteine) — see Flags 8.

Every row below shares: cysteine 48.3 µmol, the co-reactant(s) at 72.5 or 144.9 µmol, potassium
phosphate **0.5 mol/L at pH 5.00**, total charge 500 mg, sealed silanized 2 mL vial, **95 C for
4 h**, duplicate, HS-SPME at 20 C for 60 min, EI-MS.

| step / quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| MFT skeleton origin | 49 unlabelled / 46 13C5; intermediate channels 5 total | % of the M+ isotopomer signal | pot C, 95 C, pH 5, 4 h | none — this is a composition, not a rate | Table 2 row 7, p. 2716 | **within_study_ratio** |
| FFT skeleton origin | 48 / 47; intermediates 4 | % | pot C | none | Table 2 row 9 | within_study_ratio |
| 3-mercapto-2-pentanone skeleton origin | 49 / 47; intermediates 4 | % | pot C | none | Table 2 row 8 | within_study_ratio |
| 2-furaldehyde skeleton origin | 47 / 45; intermediates 7 | % | pot C | none | Table 2 row 6 | within_study_ratio |
| 3-mercapto-2-butanone: one carbon lost | 49 / 50 with the maximum at **4** labels | % | pot C | none | Table 2 row 5 | within_study_ratio |
| 3-thiophenethiol is cysteine-derived | 95 % unlabelled, maximum 4 labels | % | pot C | none | Table 2 row 11 | within_study_ratio |
| thiazole is one-carbon-from-ribose | 56 / 44 single label, at **C-2** | % | pot C | none | Table 2 row 3 + p. 2716 | within_study_ratio |
| 2-methyl-3-(methylthio)furan | 36 / 11 / 1 / 0 / 2 / 36 / 13; ring from ribose, **half the SMe carbon from cysteine** | % | pot C | none | Table 2 row 10 + Fig. 4 | within_study_ratio |
| disulfide composition, MFT dimer | 26 / 48 / 22 at 0 / 5 / 10 labels (binomial 25/50/25) | % | pot C | none | Table 2 row 13 | within_study_ratio |
| disulfide composition, FFT dimer | 28 / 54 / 18 | % | pot C | none | Table 2 row 15 | within_study_ratio |
| **MFT: norfuraneol route share** | **7 % norfuraneol / 93 % ribose** at equimolar charge | % | pot D, 95 C, pH 5, 4 h | none | Table 3 row 7 | **within_study_ratio** — this is the "NF <= 7 % ceiling" the hold-out declaration names |
| 2-mercapto-3-pentanone: norfuraneol route share | 96 % norfuraneol / 4 % ribose | % | pot D | none | Table 3 row 16 | within_study_ratio |
| 3-mercapto-2-butanone: norfuraneol route share | 94 / 6 | % | pot D | none | Table 3 row 5 | within_study_ratio |
| 3-mercapto-2-pentanone: split route | 42 norfuraneol / 58 ribose | % | pot D | none | Table 3 row 8 | within_study_ratio |
| FFT: norfuraneol route share | **0 / 100** | % | pot D | none | Table 3 row 9 | within_study_ratio |
| 2-furaldehyde: norfuraneol route share | 2 / 98 | % | pot D | none | Table 3 row 6 | within_study_ratio |
| MFT disulfide in pot D | 1 / 15 / 84 at 0 / 5 / 10 labels | % | pot D | none | Table 3 row 13 | within_study_ratio |
| **FFT: furfural route share** | **92 % from 2-furaldehyde / 8 % from ribose at equimolar loading** | % | pot E, 95 C, pH 5, 4 h | none | p. 2718 running text | within_study_ratio |
| 2-mercapto-3-pentanone from ribose alone | **not detected** | — | pot C | — | p. 2715 | **threshold** (a verified absence at this method's detection limit, which is not stated) |
| route ratio MFT ribose : norfuraneol | **13.3 : 1** | — | pot D, equimolar | — | derived from Table 3 (mine) | within_study_ratio (mine) |
| route ratio 2-mercapto-3-pentanone norfuraneol : ribose | **24 : 1** | — | pot D | — | derived (mine) | within_study_ratio (mine) |
| intermediate-channel budget from 98 % enrichment | ~4-5 points expected with zero fragmentation | % | any C5 product | — | derived (mine) | derived_assumption (mine) |

### Can any of this be put on the same basis as a shipped constant? No, and here is exactly why.

A share of an isotopomer signal has no dimension of amount and no dimension of time. It constrains
**which of two parallel routes carried the atoms**, and nothing else. Concretely:

- It licenses a **flux ratio at one condition**: at 95 C and pH 5 with equimolar precursors, the
  ribose-direct route to MFT carries 93 % of the flux and the norfuraneol route 7 %. In the shipped
  network that is a constraint on `k_ddp_mft`(+`k_ddp_mft_hs`) against `k_nf_mft` **evaluated in
  Cerny's own pot**, which is precisely how `FIT_HOLDOUT_DECLARATION.md` line 43 says to use it:
  "the NF <= 7 % ceiling must be **evaluated at Cerny's conditions**."
- It does **not** license any transport to 145 C, the sulfur module's `T_REF_S_K`. The paper itself
  supplies the reason: at 140-145 C (refs 8 and 30) both pentanone isomers appear, and at 95 C only
  one does. **The route mix is measured to change with temperature**, which is the declaration's
  stated ground for holding this paper out.
- It does not license any statement about the level of anything, because the analysis is
  headspace SPME with no internal standard.

## 5. Flags

1. **This paper contains no rate, no yield, no concentration and no mass balance.** It cannot be
   made into a scored row of any objective that measures amounts. Anyone reaching for it to explain
   a missing thiol sink should stop here: the strongest sink-relevant sentence in it is a citation of
   van Seeventer 2001 (this cluster's third paper), and the strongest sink-relevant datum of its own
   is a *composition* (the binomial disulfides) with no amount attached.
2. **The 1,4-dideoxyosone route is PROPOSED, not demonstrated.** The paper's own language throughout
   Figures 2 and 6 is conditional: "is proposed as an intermediate", "if proven", "could be a
   suitable precursor", "further experiments are necessary to verify whether the new pathway
   proposed in Figure 2 can be considered to be the main reaction mechanism". The 1,4-dideoxyosone
   of ribose was **never detected in this study** — its prior identification is in a different
   system (ref 32, xylose + hydrolysed wheat protein). The repository's `DDP` species and the
   `r_ddp_mft` reaction are therefore built on a *hypothesis this paper advances*, supported by the
   negative evidence that norfuraneol cannot be the intermediate. That is a defensible modelling
   choice and it should be labelled as one wherever it appears; the note in `sulfur.py` currently
   states the isotope shares (which are measured) and the route name (which is not) in one breath.
3. **The disulfides may be SPME artefacts and the authors say so.** p. 2717: solvent extraction was
   avoided, which removes the Hofmann 1996 artefact mechanism, "However ... To verify whether the
   disulfides represent true reaction products or are formed as artifacts during SPME, additional
   experiments would be necessary." Any use of rows 13-15 as evidence that a disulfide pool exists
   in the pot must carry that sentence. Note that an SPME-borne artefact would *still* be binomial
   on the free thiol pool, so my section-3 arithmetic does not distinguish the two.
4. **The 98 % enrichment is uncorrected.** No isotopic-purity correction is described anywhere. My
   estimate that isotopic impurity alone accounts for 4-5 points of intermediate-channel signal in a
   C5 product assumes the 98 % is per-carbon uniform labelling; if it is instead 98 % of molecules
   fully labelled, the effect is about 2 points. The paper does not say which, and the difference
   matters for exactly one claim — whether the residual intermediate channels are zero or merely
   small.
5. **One row is unusable and the paper says so**: 2-methyl-3-thiophenethiol (row 12, shares
   44/2/6/14/10/24) coeluted with another compound, so "an unambiguous determination of the
   isotopomer ratio ... was not possible". Its 32 points of intermediate channel are the largest in
   Table 2 and must not be read as evidence of fragmentation.
6. **Duplicate, with no dispersion reported.** "All samples were analyzed in duplicate" and not one
   error bar, standard deviation or range appears in either table. Every share in this dossier is a
   point estimate of unknown precision. For shares near 50 the two-pot logic is self-checking (the
   unlabelled and labelled channels must be equal), and rows 7, 8 and 9 pass that check to within
   3 points — which is the only precision estimate available and it is mine.
7. **What to request from the authors**: (i) the duplicate values behind every share, or a stated
   precision; (ii) the total ion counts or any quantification at all, which would turn the route
   shares into route fluxes; (iii) the promised follow-up — "CAMOLA experiments with 4-hydroxy-5-
   methyl-3(2H)-furanone/cysteine ... Further experiments are currently under way" (p. 2718) — which
   would settle whether the mercaptoketones come from an intact norfuraneol chain; (iv) whether the
   disulfides survive a non-SPME sampling; (v) the isotopic-purity convention behind "98 %
   enrichment".
8. **Registry gaps against `data/keys/compounds.yml`.** Nine of this paper's fifteen Table 2
   compounds have no id: 3-mercapto-2-pentanone, 2-mercapto-3-pentanone, 3-mercapto-2-butanone,
   2-methyl-3-(methylthio)furan, 3-thiophenethiol, 2-methyl-3-thiophenethiol, bis(2-furfuryl)
   disulfide, the mixed (2-methyl-3-furyl)(2-oxo-3-pentyl) disulfide, and thiazole. The first three
   are the α-mercaptoketones the network already carries internally as `MP3P` / `MP` and they are
   the compounds the paper's sharpest route contrasts are about; **`3-mercapto-2-pentanone` and
   `2-mercapto-3-pentanone` are the two that most deserve keys**, because the presence/absence
   contrast between pot C and pot D is a topology test the lane could score if the species were
   nameable. Neither ribose nor cysteine is in the registry (it carries no Maillard reactants at
   all).
9. **What this paper does not contain**: any temperature other than 95 C; any pH other than 5.00;
   any time other than 4 h (no time course); any buffer strength other than 0.5 mol/L; any
   oxygen/argon contrast; any water-activity variation; any quantity of any kind; any figure with a
   numeric axis; any supplementary material.
