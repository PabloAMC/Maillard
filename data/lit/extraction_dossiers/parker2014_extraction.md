# Parker, Bishara, Baines & Mottram 2014 — EXTRACTION (cysteine 25 mM + ribose or ascorbic acid 25 mM, 0.2 M pyrophosphate pH 5.5 or unbuffered water, 140 °C, 30 min, one time point; volatiles by dynamic headspace and SPE GC-MS, semi-quantitative; NO H2S, NO thiazolidine, NO time course measured)

### A FOUR-PAGE PROCEEDINGS CHAPTER WITH ONE MEASURED CONTRAST AND ONE PROPOSAL. Measured (semi-quantitatively, one time, 140 °C, pH 5.5): norfuraneol is ~5× higher with ribose than with ascorbic acid, yet MFT is detected only with ascorbic acid and the MFT disulfides are far higher there. Proposed, not measured: the ribose-cysteine thiazolidine ties up the cysteine that would supply H2S. Neither H2S nor the thiazolidine was measured.

**Source on disk:** `data/articles/parker2014.pdf` (4 pp., book pp. 211-214, born-digital). Prose, Table 40.1 and
references read from the `pdftotext -layout` text layer; Fig. 40.1 (m/z 113 extracted-ion chromatograms) and
Fig. 40.2 (scheme) read from 150 dpi page renders (scratchpad `img/parker-3.png`, `img/parker-4.png`).

**Legend.** **[V]** verbatim (short quote) · **[M]** a number or result printed in this paper · **[F]** read from a
figure · **[Q]** qualitative statement · **[C]** cited from another paper · **[P]** proposed, not measured ·
**[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The questions put to it: the 100 °C MFT lag in ribose + cysteine and Cerny 2003's isotope constraint (fed norfuraneol
is a minor MFT source at 95 °C) — `results/validation/kinetic_core_thiol_submodel_prereg.md` secs. 13-16, where
mechanism (2) "slow H2S supply as a co-limiter" (prior 0.20) and a "ribose-cysteine thiazolidine reservoir" are
candidates.

| what the model needs | what Parker 2014 gives | tag |
|---|---|---|
| time course of H2S | **none. H2S was not measured** (not even detected/reported as a volatile) | — |
| time course of MFT | **none.** One time point (30 min at 140 °C). MFT detected in CAB (ascorbic, buffered) duplicates; **not detected** in CRB (ribose, buffered) on m/z 113 | [M] / [F] |
| MFT disulfides | peaks 4-7 (2-methyl-3-furyl 1-methyl-2-oxopropyl disulfide, bis(2-methyl-3-furyl) disulfide, 2-methyl-3-furyl 2-thienyl disulfide, 2-methyl-3-furyl 2-methyl-3-thienyl disulfide) "relatively high abundance" in CAB vs "trace levels" in CRB | [M] / [F] |
| norfuraneol | most abundant SPE compound in both buffered pots; **five-fold higher in CRB than CAB** | [M] |
| thiazolidine (TTCA-type, ribo) | **not measured.** Proposed to sequester cysteine in CRB | [P] |
| furfural | ~40 % (CAU) and ~60 % (CRU) of total DHE peak area in the unbuffered pots | [M] |
| compound counts | >150 identified in total; >120 in each buffered pot; 58 in CAU; 30 in CRU | [M] |
| quantitation | semi-quantitative against an internal standard (1,2-dichlorobenzene for DHE; 3-chlorophenol for SPE); **no response factors, no absolute amounts printed** | [M] |

**Answer, measured vs proposed.**
- **Measured** [M]: at 140 °C / 30 min / pH ~5.5 (0.2 M pyrophosphate), equimolar 25 mM, NF(CRB)/NF(CAB) ≈ 5 while
  MFT(CRB) is below detection and MFT(CAB) is detected; disulfides likewise CAB ≫ CRB. So **MFT is anti-correlated
  with norfuraneol across the two precursors.**
- **Proposed** [P]: "the H2S required for the formation of MFT is tied up in the semi-stable" ribose-cysteine
  thiazolidine [V], which "is readily formed" (cited to de Roos 1992, ACS Symp. Ser. 409, 203-216 [C]); ascorbic acid
  reaches 1-deoxypentosone by Vernin 1998's route [C], bypassing the thiazolidine and leaving cysteine free to degrade
  to H2S. (Chemically it is cysteine, the H2S precursor, that would be tied up, not H2S.)

**Bearing on the thiol sub-model (my reading).**
1. **Consistent with Cerny 2003's constraint** [Z], confidence ~0.7: if norfuraneol were the main MFT precursor, the
   pot with 5× more NF should make more MFT, not less. Two readings remain and Parker cannot separate them:
   (a) **H2S / free-cysteine supply limits MFT** in ribose + cysteine (the authors' proposal; mechanism (2) and the
   thiazolidine-reservoir rival of sec. 13), and (b) **a non-NF intact-C5 route dominates MFT** and ascorbic acid feeds
   it better (mechanism (1) of sec. 14, via 1-deoxypentosone). Both predict MFT not tracking NF. (b) gets slightly
   more support from the fact that the ascorbic pot makes *less* NF yet *more* MFT, which a pure H2S limit on a
   shared NF route would not produce unless H2S differed by >5× [Z].
2. **On the 100 °C lag: no information** (confidence ~0.95). One temperature (140 °C), one time (30 min), no H2S.
   The thiazolidine-reservoir idea is a qualitative prior for R4-type arms, not a rate.
3. **Weak support for "the thiazolidine reservoir delays sulfur supply"**: a proposal, plausibly correct in sign
   (cf. `cao2017_extraction.md`: ~89 % of 0.4 M cysteine is in TTCA by 30 min at 90 °C; `zhai2019_extraction.md`),
   but Parker adds no measurement of it. Do not cite it as evidence that the reservoir limits H2S.

**Confounds the chapter does not address** [Z]: ascorbic acid is a reductant and its degradation products
(dehydroascorbic acid, 2,3-diketogulonic acid, C5 fragments after decarboxylation) are not ribose's; its pot also lacks
ribose's Amadori route entirely. The MFT contrast therefore varies the C5 precursor chemistry and the redox state
together with thiazolidine formation. Ascorbate also forms no thiazolidine with cysteine in the way an aldose does —
which is the authors' point, but it is one of several differences.

**Digitisable data?** Fig. 40.1 is a stacked m/z 113 extracted-ion trace (duplicate CRB and CAB, offset baselines,
full scale 60,000 counts, several peaks overrun the trace above). Peak-height ratios against the internal-standard
peak 2 could be read, but single-ion, no response factors, overlapping offsets: **not worth two-reader keying**. The
usable content is the sign (MFT CAB > CRB, MFT in CRB below detection) and the printed NF ratio (5×).

## 0. Identity

| field | value |
|---|---|
| Title | "Comparison of Ribose and Ascorbic Acid in Model Process Reactions Containing Cysteine" (Chapter 40) |
| Authors | Jane K. Parker, Sandra Bishara (Dept. of Food and Nutritional Sciences, University of Reading); David A. Baines (Baines Food Consultancy, Thornbury, Bristol); Donald S. Mottram (Reading) |
| Venue | V. Ferreira & R. Lopez (Eds.), *Flavour Science* (symposium proceedings, cited elsewhere in the repo as the Weurman symposium; the chapter itself does not name it), Academic Press / Elsevier, pp. 211-214. Footer prints "© 2014 2013 Elsevier Inc."; other repo documents list it as "Parker 2013" (`docs/guides/CITATION_CONTEXT_2026-10.md`). Not the same paper as `parker2013_extraction.md` (J. Sci. Food Agric. 2013) |
| DOI as printed (p. 211 footer) | **10.1016/B978-0-12-398549-1.00040-4** ✔ matches the brief; Crossref returns this title, these four authors, *Flavour Science* pp. 211-214, 2014 (checked 2026-10-31) |
| Paper type | Short conference chapter; single-time-point model reactions; volatile profiling |
| Same lab | Mottram / Parker, Reading: `parker2013_extraction.md`, `balagiannis2010_extraction.md`; Whitfield & Mottram 2001 (NF + cysteine / H2S at pH 6.5, cited as ref 1) |

## 2. Methods as they matter to a model

- **Pots** (Table 40.1) [M]:

| code | cysteine (mM) | ascorbic acid (mM) | ribose (mM) | buffer (M) | initial pH | final pH |
|---|---|---|---|---|---|---|
| CAB | 25 | 25 | – | 0.2 | 5.4 | 5.6 |
| CRB | 25 | – | 25 | 0.2 | 5.5 | 5.6 |
| CAU | 25 | 25 | – | – | 3.0 | 3.0 |
| CRU | 25 | – | 25 | – | 4.7 | 3.8 |

- **Buffer** [M]: "pyrophosphate buffer (pH 5.5)" in the text; Table 40.1 gives 0.2 M. Counter-ion not stated.
  Unbuffered = de-ionised water. (Pyrophosphate, not orthophosphate: relevant if buffer catalysis is modelled per
  buffer species.)
- **Heating** [M]: 20 mL in airtight screw-top reagent bottles, **140 °C, 30 min**, then allowed to cool (heating medium,
  come-up and cool-down not stated). Internal standard 1,2-dichlorobenzene added before heating.
- **DHE** [M]: 5 mL heated sample + 5 mL water, 250 mL flask at 60 °C, N₂ 40 mL/min for 1 h onto Tenax; thermal
  desorption; GC-MS (HP5890 II / HP5972), VF-5MS 30 m × 0.25 mm × 0.25 µm.
- **SPE** [M]: 20 mL + 1 µL 3-chlorophenol (1 mg/mL in 5 % MeOH); GC-MS on DB-Wax 30 m.
- **Replication** [M]: duplicates shown in Fig. 40.1; no SDs printed.

## 3. Results as printed

| quantity | value | conditions | where | tag |
|---|---|---|---|---|
| NF, CRB / CAB | **5** (five-fold, ribose higher) | SPE, 140 °C 30 min, pH 5.5 | p. 213 | [M] |
| NF rank | most abundant SPE compound in both CRB and CAB | same | p. 213 | [M] |
| MFT (peak 1, m/z 113) | present in CAB duplicates; absent in CRB duplicates | DHE, same | Fig. 40.1, p. 213 | [M] / [F] |
| MFT disulfides (peaks 4-7) | CAB relatively high; CRB trace | DHE, same | Fig. 40.1 | [M] / [F] |
| furfural share of total DHE area | ~40 % CAU; ~60 % CRU | unbuffered, pH 3.0 / 4.7→3.8 | p. 212 | [M] |
| volatiles identified | >150 total; >120 in each of CRB, CAB; 58 CAU; 30 CRU | DHE | p. 212 | [M] |
| buffered-pot classes | up to 60 thiophenes, 10 thiazoles, 20 cyclic sulfides, 25 MFT- or FFT-moiety sulfides | DHE | p. 213 | [M] |
| bis(2-methyl-3-furyl) disulfide odour threshold | 0.02 ng/kg | — | p. 213 | [C] (unreferenced) |

## 4. What the repository takes

At most two directional claims at 140 °C, 30 min, pH 5.5 (0.2 M pyrophosphate), cysteine : C-source 25 : 25 mM:
(i) norfuraneol is ~5× higher from ribose than from ascorbic acid; (ii) MFT and its disulfides are higher from
ascorbic acid than from ribose (MFT below detection with ribose). Together they are a sign test that MFT does not
track norfuraneol across C5 precursors — compatible with Cerny 2003's isotope result. No FIT row (no absolute amount,
one time point, no response factors); the engine has no ascorbic-acid lane, so claim (ii) is not scoreable today. The
thiazolidine-limits-H2S statement is a **proposal** and may be cited only as such.

## 5. Flags

1. **"H2S is tied up in the thiazolidine" is imprecise** — the thiazolidine binds cysteine; any H2S limit follows from
   less free cysteine. And H2S was never measured.
2. **Semi-quantitative only.** Internal-standard ratios without response factors; Fig. 40.1 is one ion (m/z 113),
   which also matches other furyl sulfides. "Five-fold" is the only ratio printed.
3. **Pyrophosphate, 0.2 M, and the buffered/unbuffered contrast is confounded with pH** (unbuffered ascorbic pot at
   pH 3.0; unbuffered ribose pot drifts 4.7 → 3.8). The buffered-vs-water comparison here says nothing about buffer
   catalysis at equal pH.
4. **The cited mechanism for ascorbic acid → 1-deoxypentosone (Vernin 1998) is cited, not tested here.**
5. **Year ambiguity:** chapter printed © 2014 (book) / 2013; repository documents use both "Parker 2013" and
   "Parker 2014" for this DOI. The dossier stem follows the file on disk (`parker2014.pdf`).
6. **Registry:** `2_methyl_3_furanthiol`, `bis_2_methyl_3_furyl_disulfide`, `norfuraneol`, `furfural`,
   `hydrogen_sulfide`, `cysteine`, `ribose` are in `data/keys/compounds.yml`; ascorbic acid and the ribose-cysteine
   thiazolidine are not.
