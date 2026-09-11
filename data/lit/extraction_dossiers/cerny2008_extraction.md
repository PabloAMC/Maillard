# Cerny & Guntz-Dubini 2008 — EXTRACTION (thiamine hydrochloride 25.00 g + cysteine 3.00 g + xylose 11.00 g in 750 g potassium phosphate 0.5 mol/L pH 5.0, stirred pressure reactor, 145 C / 45 min; isolation and MS + NMR identification of 5-hydroxy-3-mercapto-2-pentanone — a STRUCTURE paper with no kinetics and no quantification)

### THE EXISTENCE PROOF FOR THE ENGINE'S `HMP` NODE AND NOTHING MORE: it establishes that 5-hydroxy-3-mercapto-2-pentanone is a real, isolable compound in a thiamine + cysteine + xylose pot at 145 C — which licenses `HMP` as a species — but it prints no concentration, no yield, no rate and no time course, so all three of the engine's HMP constants (`k_thi_hmp`, `k_hmp_mft`, `k_hmp_mp2p`) stay exactly as fitted as they were.

**Source on disk:** `data/articles/cerny2008.pdf` (4 pp., J. Agric. Food Chem. 2008, 56 (22), 10679-10682).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/cerny2008.txt`, 191 lines). The text layer is clean and complete: the running
text, both column blocks, the Materials and Methods, and all twelve references came through legible.
**This paper contains NO printed table.** Every numbered display is a figure: Figure 1 (proposed
formation pathway of 2-methyl-3-furanthiol from thiamine via 5-hydroxy-3-mercapto-2-pentanone,
redrawn from reference 6), Figure 2 (EI mass spectrum of the compound), Figure 3 (13C NMR spectrum),
Figure 4 (1H NMR spectra, A isolated and B synthesized), Figure 5 (EI mass spectrum of the
trimethylsilylated oxime). The spectral assignments quoted below are printed **in the running text**,
not read off the figures; the spectra themselves are figure-only. There is no supplementary material.
Repo status before this dossier: `HMP` is carried as a species in
`src/kinetic_core/species_sulfur.py` (line 128, molar mass 134.20 at line 447) and reaches the
network through four rows in `src/kinetic_core/sulfur.py` (`r_thi_hmp`, `r_hmp_mft`, `r_hmp_mp2p`,
`r_hmp_decay`), but **no extraction dossier for its identification paper existed** and none of its
constants carries a literature note.

## 0. Identity

| field | value |
|---|---|
| Title | "Identification of 5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine, Cysteine, and Xylose" |
| Authors | Christoph Cerny (corresponding; phone +41 22 780 22 11, christoph.cerny@firmenich.com) and Renée Guntz-Dubini — Firmenich SA, Corporate R&D Division, P.O. Box 148, CH-1217 Meyrin, Geneva 2, Switzerland |
| Venue | J. Agric. Food Chem. 2008, 56 (22), 10679-10682. Received for review 9 June 2008, revised 2 October 2008, accepted 5 October 2008, web 5 November 2008 |
| DOI / article ID | 10.1021/jf801762c (printed as `JF801762C`) |
| Paper type | **Identification / structure elucidation.** Not a kinetic paper, not a quantification paper, not a comparative-yield paper. The abstract states the whole scope: "We have identified the compound in a thermally treated mixture of thiamine, cysteine, and xylose and characterized it by MS and NMR." |
| Compound numbering in the paper | **1** = thiamine; **2** = 5-hydroxy-3-mercapto-2-pentanone (open chain); **3** = the cyclic ketal formed by cyclisation of **2**; **4** and **5** = two thiamine-derived aroma compounds named only by number in the Results. **My reading**, from the intro's list of van der Linde's products and from Figure 1's caption: **5** = 2-methyl-3-furanthiol and **4** = 4,5-dihydro-2-methyl-3-furanthiol. This is an inference from context, not printed (Flags 5). |
| Lineage | the compound was **postulated** as a key thiamine intermediate by Van der Linde et al. 1979 (ref 6) and by Güntert et al. 1993 (ref 7), who "consider [it] unstable and very reactive" and did not identify it; Matsukawa et al. 1948/1951 (refs 8, 10) reported the acetylated compound and a 1951 identification in refluxed aqueous thiamine with no MS/NMR data; Onural 1991 (ref 9) reported a synthesis "but the author does not give any details on the procedure, yields, or analytical data" |
| Predecessor pot | ref 11 = Cerny, C., "Formation of aroma compounds in the Maillard reaction of xylose, cysteine and thiamine", Recent Highlights in Flavor Chemistry & Biology, 2008, pp. 261-264. The present paper says the model reaction here used "similar reaction conditions" to that study. |
| Companions on disk | `cerny2007_extraction.md` (Cerny & Briffod 2007, the [13C5]xylose + cysteine + thiamin pH ladder at **145 C / 20 min** — the same laboratory, same three precursors, same temperature, and the source of the engine's thiamine-vs-xylose isomer split), `cerny2004_extraction.md`, `cerny1994_extraction.md` |

## 1. Why it matters

The engine's thiamine route (`src/kinetic_core/sulfur.py`, the block headed "THE THIAMINE ROUTE",
lines 511-529) is four reactions wide and passes entirely through one node:

| row in `sulfur.py` | rate key | what this paper says about it |
|---|---|---|
| `r_thi_hmp`: THI -> HMP + FRAG_C 7 + FRAG_N 4 | `k_thi_hmp` (order 1) | **the product exists and is isolable** — nothing about the rate |
| `r_hmp_mft`: HMP -> MFT | `k_hmp_mft` (order 1) | Figure 1 draws this step, but as a **proposal taken from reference 6**, not measured here |
| `r_hmp_mp2p`: HMP -> MP2P | `k_hmp_mp2p` (order 1) | 3-mercapto-2-pentanone is named as a main SPME peak in the **predecessor** pot (ref 11), not measured here |
| `r_hmp_decay`: HMP -> FRAG_C 5 + FRAG_S | `k_thiol_decay` (shared) | the compound is **stable enough to isolate**, which argues against a fast unassigned sink (Flags 4) |

In `src/kinetic_core/parameters_sulfur.py` the step table at lines 1956-1964 gives `k_thi_hmp` and
`k_hmp_mft` an **empty note string** — they are the only two steps in the thiamine block with no
literature remark of any kind — and all three are routed to the `thiol_assembly` barrier family
(lines 2401-2402), i.e. they carry a family Ea, not their own. None of the three appears in a
`MEASURED_SULFUR` row. So the whole thiamine lane's kinetics is fitted, and this paper **does not
change that**: it supplies zero rates and zero concentrations.

What it does supply is the one thing a species registry needs and the engine did not have on file:
**a primary identification.** Before this paper the compound was, in the authors' own account of the
literature, a postulate — proposed as an intermediate by two groups, never isolated, "considered
unstable and very reactive", with "no scientific publication which reports on the spectroscopic
properties". This paper isolates 20 mg of it from an actual thiamine + cysteine + xylose reaction at
145 C, matches its EI-MS and its 1H and 13C NMR against an independently synthesized reference, and
confirms the elemental composition **C5H10O2S, M+ m/z 134** — which is exactly the molar mass
(134.20) and the C5/S1 formula `species_sulfur.py` already carries. That is an independent check on
the species record, and it is the only one available.

Two findings bear on the network as it is built, and both are cautionary rather than enabling:

**(a) The node is not one compound.** The isolated material is an equilibrium mixture of the open
chain **2** and its cyclic ketal **3**, at roughly 69:31 by 1H integration in the isolated sample and
79:21 in the synthesized reference (section 3). The engine's `HMP` is a single state variable, and
its two productive fates (`r_hmp_mft` to the furanthiol, `r_hmp_mp2p` to the mercaptopentanone) are
plausibly not fed by the same tautomer. Nothing here resolves that; it is a declared lump, and this
dossier is the first place it is written down.

**(b) The node is invisible to the method the sibling pot uses.** The Results state that GC-MS of the
diethyl ether extract gave two major peaks — sulfurol and this compound — and that **neither
corresponded to the volatiles seen before by HS-SPME** in the predecessor pot, the authors' reading
being that "the sulfurol and the unknown compound are not volatile enough or too polar to be detected
by SPME under the experimental conditions." The engine's thiamine-vs-xylose isomer diagnostic comes
from Cerny 2007 Table 4, an **HS-SPME** measurement in the same laboratory at the same 145 C. So the
intermediate that the engine routes the entire thiamine flux through is, by this paper's own
statement, **not measurable by the method that produced the data the branch is calibrated against.**
That is the single most useful sentence in this paper for the repository, and it is the reason no
amount of re-reading Cerny 2007 will ever produce an `HMP` concentration.

What this paper does NOT give: any concentration, any yield, any time course, any temperature series,
any pH series, any rate, any barrier, any mass balance, any comparison of routes, and any measurement
at all of MFT or of 3-mercapto-2-pentanone in this pot.

## 2. Methods as they matter to a model

- **The pot.** "Thiamine hydrochloride (25.00 g), cysteine (3.00 g), and xylose (11.00 g) were
  dissolved in potassium phosphate buffer (750 g, 0.5 mol/L, pH 5.0) and thermally reacted in a
  stirred pressure reactor (Glass, type 2, 1-L, Büchi, Uster, Switzerland) at **145 °C for 45 min**."
  Note the buffer is given as a **mass** (750 g), not a volume, so every concentration below is mine
  and approximate (Flags 2). **Stirred**, and a sealed pressure reactor — 145 °C is above the boiling
  point, so this is a closed system under autogenous pressure, which matters for H2S retention.
- **pH.** Set to 5.0 by 0.5 mol/L potassium phosphate at t = 0. **Never re-measured**, and no drift is
  reported. Treat 5.0 as an initial pH. (Cerny 2007, same laboratory, runs a five-point pH ladder
  4.0-7.0 at 145 C and shows the sulfur volatiles move strongly across it.)
- **Work-up.** After cooling to room temperature the product was saturated with sodium chloride and
  extracted twice with diethyl ether (300 mL + 200 mL); the organic phase was dried over anhydrous
  sodium sulfate and concentrated to approximately 1.5 mL on a rotary evaporator. **Solvent
  extraction, not headspace** — this is the methodological difference from the SPME pots, and it is
  why this compound was seen at all.
- **Isolation.** Flash chromatography on silica gel (20 × 2.5 cm), eluent cyclohexane / ethyl acetate
  / ethanol 8 + 1 + 1 (v/v); target fractions combined, a spatula of silica added, concentrated to
  dryness, and **submitted to flash chromatography a second time** with the same eluent. "The purest
  fraction after solvent evaporation yielded **20 mg**." That 20 mg is a mass of purified isolate
  after two columns with no recovery correction — it is a floor on what was formed, **not a yield**
  (Flags 1).
- **GC-MS.** Agilent GC 6890A coupled to an MSD 5973; HP-5MS capillary 30 m × 0.25 mm, film 0.25 µm;
  oven 50 °C for 1 min, then 10 °C/min to 250 °C; 1 µL injected, split 1:30; EI at 70 eV, scan range
  m/z 15-400. **No internal standard is mentioned, no calibration, no quantification of any kind.**
- **NMR.** Bruker DPX 400, in deuterochloroform, tetramethylsilane as internal standard.
- **Derivatisation (structure confirmation only).** The deuterated solvent was blown off under
  nitrogen; hydroxylamine hydrochloride in pyridine (500 µL, 2.5 % w/v) added; the closed vial heated
  30 min at 70 °C; hexamethyldisilazane (900 µL) and trifluoroacetic acid (100 µL) added; shaken and
  left 15 min; analysed by GC-MS. This makes the trimethylsilylated **oxime**.
- **Independent synthesis of the reference.** From 5-(triethylsilyloxy)-3-thioacetyl-2-pentanone
  (provided by Dr. Roger Snowden), itself obtained from 3-chloro-5-hydroxy-2-pentanone by
  triethylsilylation of the hydroxyl and substitution of the chloride by a thioacetoxy group. The
  silyl ether and thioacetate were then cleaved: 26.5 mg of the protected compound in 2.50 mL
  tetrahydrofuran added dropwise at 4 °C to aqueous lithium hydroxide (5 %); after 60 min at 0-4 °C
  added to aqueous ammonium chloride (8.9 %) at 0-4 °C; warmed to 20-25 °C; saturated with sodium
  chloride and extracted with diethyl ether (5 × 0.50 mL); the aqueous phase adjusted to pH 4-5 with
  0.5 mol/L sulfuric acid and extracted with dichloromethane (6 × 2 mL); dried over sodium sulfate;
  solvent evaporated. **No yield is reported for this synthesis either.**
- **Reagent sources.** Dichloromethane, hydroxylamine hydrochloride, lithium hydroxide, pyridine,
  sodium chloride: Merck. Ammonium chloride, ethyl acetate, trifluoroacetic acid: Acros. Cysteine,
  hexamethyldisilazane, anhydrous sodium sulfate, sulfuric acid, **thiamine hydrochloride, xylose**:
  Fluka. Diethyl ether, ethanol: Carlo Erba. Cyclohexane: Riedel de Haën. Silica gel (32-63, 60 Å):
  Brunschwig. All reagents analytical grade.
- **Replication.** **None reported.** One reaction, one isolation, one set of spectra.
- **Reference temperature.** 145 °C = 418.15 K. This is the same temperature as Cerny 2007 (which
  ran 20 min against this paper's 45 min) and matches the sulfur lane's own 145 C fit rows, so no
  temperature transport would be needed — if there were a number here to transport.

## 3. Tables re-typed

**There is no table in this paper.** The four printed tables the house format expects do not exist:
the paper's entire data content is five figures plus spectral assignments quoted in the running text.
Everything below is transcribed from the running text and is printed, not figure-read; the spectra in
Figures 2-5 are figure-only and are not typed as numbers.

### 3.1 13C NMR assignments of the open-chain compound 2 (Results, p. 10682, printed in text)

Recorded in CDCl3, TMS internal standard, Bruker DPX 400. Introduced as "The following 13C NMR
signals (δ/ppm) agree with the chemical structure of 2".

| assignment | δ / ppm |
|---|---|
| C, C-2 | 206.1 |
| CH2, C-5 | 59.9 |
| CH, C-3 | 44.7 |
| CH2, C-4 | 36.4 |
| CH3, C-1 | 27.5 |

### 3.2 13C NMR assignments of the cyclic ketal 3 (Results, p. 10682, printed in text)

Introduced as "In addition, ca. **30 %** of ketal 3, formed by cyclization of 2, is present, as
suggested by the following 13C NMR signals (δ/ppm)".

| assignment | δ / ppm |
|---|---|
| C, C-2 | 103.4 |
| CH2, C-5 | 64.9 |
| CH, C-3 | 45.8 |
| CH2, C-4 | 34.5 |
| CH3 | 24.3 |

The C-2 shift moves from 206.1 (ketone carbonyl) to 103.4 (ketal carbon) — that pair is the whole
structural argument for the cyclisation.

### 3.3 1H NMR (Figure 4; the two singlets and the integration are printed in text)

| quantity | value | as printed |
|---|---|---|
| methyl singlet of the open-chain molecule 2 (CH3, C-1) | 2.35 ppm | "The singlet signals at 2.35 ppm and 1.52 ppm correspond to the methyl protons of the open chain molecule 2 (CH3, C-1) and of the cyclic form 3 (CH3), respectively." |
| methyl singlet of the cyclic form 3 (CH3) | 1.52 ppm | same sentence |
| open-chain share, **isolated** compound | **69 %** | "Integration of the peak areas indicates a share of 69 % open chain form in the isolated compound and 79 % for the synthesized reference." |
| open-chain share, **synthesized** reference | **79 %** | same sentence |
| number of possible isomers | **six** | "Because of the presence of one chiral carbon in 2 and two chiral carbons in 3, in total six isomer compounds can be present. Information on the ratio of the optical isomers was not obtained from the 1H NMR spectra." |

### 3.4 EI mass spectrum of the free compound (Figure 2; the assignments are printed in text)

| m/z | assignment as printed |
|---|---|
| 134 | "corresponds to the molecular ion M+" |
| 116 | "The loss of water ... is characteristic for alcohols" |
| 90 | "can be explained by a McLafferty rearrangement with neutral loss of acetaldehyde" |
| 43 | "results from the separation of the acetyl group" |

**Relative intensities are not printed anywhere** — only these four assignments. The full spectrum is
Figure 2 and is figure-only.

### 3.5 EI mass spectrum of the trimethylsilylated oxime (Figure 5; assignments printed in text)

| m/z | assignment as printed |
|---|---|
| 365 | "corresponds to the molecular peak (M+)" |
| 350 | "indicates the loss of CH3" |
| 249 | "The base peak ... originates most likely from a McLafferty rearrangement with loss of the neutral enol fragment trimethylsilyloxyethene" |
| 73, 103, 147 | "common fragments for trimethylsilylated compounds" |

### 3.6 Everything else printed as a number in this paper

| quantity | value | where |
|---|---|---|
| thiamine hydrochloride charged | 25.00 g | Materials and Methods, Model Reaction |
| cysteine charged | 3.00 g | same |
| xylose charged | 11.00 g | same |
| potassium phosphate buffer | 750 g, 0.5 mol/L, pH 5.0 | same |
| reactor | 1-L stirred glass pressure reactor (Büchi, type 2) | same |
| reaction temperature and time | **145 °C, 45 min** | same |
| ether extraction | 300 mL + 200 mL, twice, after NaCl saturation | same |
| extract concentrated to | ~1.5 mL | same |
| flash column | silica gel, 20 × 2.5 cm; cyclohexane/EtOAc/EtOH 8+1+1 v/v; run **twice** | Isolation |
| **purest fraction after evaporation** | **20 mg** | Isolation — a mass of isolate, **not a yield** (Flags 1) |
| protected precursor taken into the deprotection | 26.5 mg in 2.50 mL THF | Synthesis |
| ketal share by 13C | ca. 30 % | Results |
| open-chain share by 1H | 69 % (isolated) / 79 % (synthesized) | Results |
| possible isomers | 6 | Results |

### Arithmetic on the printed charges (all mine)

**1. Molar charges.** Thiamine hydrochloride (thiamine chloride hydrochloride, M = 337.27 g/mol):
25.00 / 337.27 = **74.1 mmol**. L-cysteine (M = 121.16): 3.00 / 121.16 = **24.8 mmol**. D-xylose
(M = 150.13): 11.00 / 150.13 = **73.3 mmol**. **Molar ratio thiamine : cysteine : xylose =
3.0 : 1.0 : 3.0** — cysteine is the minor component by a factor of three, which is the opposite of
the usual sulfur-lane model pot and matters if anyone tries to compare this pot's sulfide supply
with the corpus.

**2. Approximate concentrations.** Treating 750 g of aqueous 0.5 mol/L phosphate as ~0.750 L (the
density of that buffer is above 1.0, so these are **upper** estimates and are approximate — the paper
gives a mass, not a volume): thiamine ~**98.8 mmol/L**, cysteine ~**33.0 mmol/L**, xylose ~**97.7
mmol/L**, phosphate 500 mmol/L. The phosphate is 5x the thiamine, so unlike the Knol-type pots this
one is genuinely buffered.

**3. The 20 mg as a floor, not a yield.** HMP is C5H10O2S, M = 134.19 g/mol (which reproduces
`species_sulfur.py`'s 134.20). 20 mg / 134.19 = **0.149 mmol**, against 74.1 mmol thiamine charged =
**0.20 %** of the thiamine on a molar basis. **This number must not be used as a yield.** It is what
survived (i) an aqueous-to-ether partition of a polar hydroxy-thiol, (ii) concentration to 1.5 mL on
a rotary evaporator, (iii) two silica flash columns, and (iv) selection of "the purest fraction" —
each of which discards material by an unmeasured amount. It is a lower bound on formation with an
unknown and probably large multiplier. It is recorded in section 4 as `level_only` with that caveat
attached, and it is the only number in this paper that could be mistaken for a yield.

**4. Formula check against the registry.** M+ = 134 with loss of water at 116 and loss of an acetyl
at 43, plus five 13C carbons at 206.1 / 59.9 / 44.7 / 36.4 / 27.5, is consistent with C5H10O2S and
with the engine's carbon count of 5 and sulfur count of 1 for `HMP`. The TMS-oxime at M+ 365 is the
same skeleton plus one oxime nitrogen and two trimethylsilyl groups, and 134 - 1 (OH proton) + 73
(TMS) = 206 on the hydroxyl, then + 15 (=N-) + 73 - 2 on the ketone gives 365 by my count; the
authors do not print the arithmetic and I do not lean on mine. **The species record is confirmed.**

## 4. Kinetic numbers the repository can use

**There are none.** This section exists to record that fact precisely, and to place the few printed
quantities in their correct evidence classes so that a later wave does not mistake any of them for a
rate or a yield.

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** `2_methyl_3_furanthiol`,
`2_furfurylthiol`, `furfural`, `2_3_butanedione`, `mercapto_2_propanone` and `thiamine_availability`
are keyed. **5-hydroxy-3-mercapto-2-pentanone is NOT in the registry**, and neither is
3-mercapto-2-pentanone, 2-mercapto-3-pentanone, sulfurol (5-(2-hydroxyethyl)-4-methylthiazole),
xylose, cysteine, or the cyclic ketal. The engine's internal keys `HMP`, `MP2P`, `MP3P`, `THI` are
network-local names and are not registry ids.

Every row below shares: thiamine hydrochloride ~98.8 mmol/L + cysteine ~33.0 mmol/L + xylose
~97.7 mmol/L in 0.5 mol/L potassium phosphate at **initial** pH 5.0, stirred 1-L glass pressure
reactor, **145 °C, 45 min**, single run, no replicate, no internal standard, no calibration.

| quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **5-hydroxy-3-mercapto-2-pentanone is formed and is isolable** from thiamine + cysteine + xylose | — (presence) | — | 145 °C, 45 min, pH 5.0, sealed stirred reactor | **none — nothing is fitted in this paper** | Results p. 10682; Abstract | **threshold** (a presence/absence result: the compound is above the detection and isolation threshold of a solvent-extraction GC-MS method) |
| mass of purified isolate, "purest fraction" after two flash columns | 20 | mg | as above, from 25.00 g thiamine hydrochloride charged | — | Materials and Methods, Isolation | **level_only** — **NOT a yield**; recovery is uncorrected and unreported (Flags 1) |
| the same mass as a molar share of charged thiamine | 0.20 | mol % (mine) | as above | — | derived from the printed 20 mg and 25.00 g (mine) | **derived_assumption** — a **floor** on formation, upper multiplier unknown |
| ketal 3 share of the isolated material, by 13C | ca. 30 | % | **in CDCl3 at NMR probe temperature**, not in the pot | — | Results, 13C paragraph | **within_study_ratio** (a solution equilibrium in an organic solvent after isolation) |
| open-chain 2 share, isolated material, by 1H integration | 69 | % | as above | — | Results, 1H paragraph | **within_study_ratio** |
| open-chain 2 share, synthesized reference, by 1H integration | 79 | % | as above | — | Results, 1H paragraph | **within_study_ratio** |
| number of possible stereoisomers (1 chiral C in 2, 2 in 3) | 6 | — | — | — | Results | **level_only** (a structural count, not a measurement) |
| 13C shifts of 2 | 206.1 / 59.9 / 44.7 / 36.4 / 27.5 | ppm | CDCl3, TMS | — | Results | **level_only** (spectroscopic identity) |
| 13C shifts of 3 | 103.4 / 64.9 / 45.8 / 34.5 / 24.3 | ppm | CDCl3, TMS | — | Results | **level_only** |
| 1H methyl singlets, 2 and 3 | 2.35 / 1.52 | ppm | CDCl3, TMS | — | Results | **level_only** |
| EI-MS ions of the free compound | 134 (M+), 116, 90, 43 | m/z | EI 70 eV | — | Results / Figure 2 assignments | **level_only** (no intensities printed) |
| EI-MS ions of the TMS-oxime | 365 (M+), 350, 249 (base peak), 73, 103, 147 | m/z | EI 70 eV, after NH2OH + HMDS/TFA | — | Results / Figure 5 assignments | **level_only** |
| Figures 2, 3, 4A, 4B, 5 (the spectra themselves, all intensities and integrals as drawn) | — | — | — | — | Figs. 2-5 | **figure_only** |
| Figure 1, the thiamine -> HMP -> MFT pathway | — | — | — | — | Fig. 1, "(6)" | **figure_only** and, additionally, a **proposal taken from reference 6**, not a result of this paper (Flags 3) |

### Can anything here be put on the same basis as the sulfur lane's constants? No, and here is why.

**(a) There is no rate and no barrier.** One temperature, one time, one run, one endpoint. Nothing in
this paper constrains `k_thi_hmp`, `k_hmp_mft`, `k_hmp_mp2p` or the shared `k_thiol_decay`, either
individually or as a ratio. The paper does not measure MFT, does not measure 3-mercapto-2-pentanone,
and does not measure thiamine loss, so not even a within-study branch ratio can be formed.

**(b) The 20 mg cannot be turned into a fed-intermediate yield.** A `fed_intermediate_yield` class
requires a known charge of the intermediate and a measured product; here the compound is the
*product* of an uncharacterised extraction chain and no downstream product was measured in the same
run. Classing it as a yield would be exactly the "a peak-area ratio is not a yield" error, one step
further removed.

**(c) The 69:31 / 79:21 ketal equilibrium does not transfer to the pot.** It is measured in
deuterochloroform at NMR temperature on isolated material, after the compound has been removed from
water. The ring-chain equilibrium of a 5-hydroxy ketone is solvent- and water-activity-dependent by
construction; in a 0.5 mol/L aqueous buffer at 145 °C it will be something else. The number is worth
recording because it says the node is a lump, and worth nothing as a partition coefficient.

**(d) What DOES transfer is a structural licence and a methodological warning.** The species record
is confirmed (C5H10O2S, 134). The `r_hmp_decay` sink should not be assumed fast (section 5, Flags 4).
And the branch this paper's compound sits on cannot be calibrated from any HS-SPME dataset, this
laboratory's included.

## 5. Flags

1. **The 20 mg is not a yield and must never be entered as one.** It is the mass of the purest
   fraction after a two-column flash purification of an ether extract that was itself concentrated to
   1.5 mL by rotary evaporation. Recovery through those four steps is unmeasured, uncorrected and
   unmentioned. The 0.20 mol % of charged thiamine that it corresponds to (mine) is a **floor**, and
   the true formation could be an order of magnitude above it or more. Class it `level_only`, carry
   the derived percentage as `derived_assumption`, and do not let it into any objective row.
2. **The pot's concentrations are mine and approximate.** The buffer is charged as **750 g**, not as
   a volume, and the reactor is a 1-L vessel filled with 750 g of liquid plus 39 g of solutes. My
   ~98.8 / ~33.0 / ~97.7 mmol/L are computed by treating 750 g as 0.750 L and are upper estimates
   (0.5 mol/L phosphate has a density above 1.0). The **molar ratio 3 : 1 : 3** is robust to that
   assumption; the absolute concentrations are not.
3. **Figure 1 is a citation, not a finding.** Its caption reads "Proposed formation pathway of
   2-methyl-3-furanthiol from thiamine via 5-hydroxy-3-mercapto-2-pentanone **(6)**" — reference 6 is
   Van der Linde et al. 1979. This paper demonstrates that the postulated intermediate **exists**; it
   does not demonstrate that MFT comes from it, and it measures no MFT at all. The engine's
   `r_hmp_mft` row therefore still rests on the 1979 proposal, and this dossier does not upgrade it.
4. **The stability claim cuts against a fast `r_hmp_decay`, but softly.** Güntert et al. (ref 7)
   "consider [the compound] unstable and very reactive" and failed to identify it; this paper's
   conclusion is "It was stable enough to be isolated and characterized." That is a statement about
   survival through a room-temperature work-up and two silica columns, **not** about its half-life at
   145 °C in an aqueous pot, and it is not evidence for a value of `k_thiol_decay`. What it does say
   is that the compound is not so transient that carrying it as a discrete pool species is
   unphysical — which is a defence of the network's topology, not of any number in it.
5. **Compounds 4 and 5 are never named in the text.** The Results say the predecessor pot's SPME run
   "revealed 3-mercapto-2-pentanone, 2-furfurylthiol, 4, and 5 as main peaks", and Figure 1's caption
   names 2-methyl-3-furanthiol as the pathway's endpoint. My reading — **5** = 2-methyl-3-furanthiol,
   **4** = 4,5-dihydro-2-methyl-3-furanthiol (the two compounds van der Linde is cited for in the
   introduction) — is an inference from context. Read Figure 1 from the page image before treating
   either identification as fixed.
6. **The `HMP` node is a lump of at least two species and possibly six.** Open chain **2** plus
   cyclic ketal **3**, at ~69:31 in CDCl3, with one chiral centre in **2** and two in **3** giving
   six possible isomers and no information on the optical ratio. The engine's single `HMP` state is a
   declared lump. If `r_hmp_mft` and `r_hmp_mp2p` turn out to be fed by different tautomers, the
   fixed branch ratio the two fitted constants imply is an approximation with no support here.
7. **No replication anywhere.** One reaction, one isolation, one NMR sample, one synthesis. No error
   bar, no standard deviation, no n, in the whole paper.
8. **This pot is not the corpus's usual sulfur pot.** Cysteine is the **minor** precursor (1 part
   against 3 parts each of thiamine and xylose), the phosphate is 0.5 mol/L, the pH is 5.0, and the
   vessel is a sealed stirred pressure reactor at 145 °C. Compare Cerny 2007 in the same laboratory:
   same three precursors and same 145 °C, but 20 min and a pH ladder. Do not merge the two pots'
   conditions.
9. **What this paper does not contain**: any concentration; any yield; any time course; any
   temperature series; any pH series; any rate constant; any activation energy; any mass balance; any
   measurement of MFT, of 3-mercapto-2-pentanone, of 2-mercapto-3-pentanone or of thiamine loss; any
   MS intensity; any quantitative NMR against a standard; any table; any replicate; any supplementary
   material.
10. **What to request from the authors**: (i) whether HMP was ever quantified in this or the
    predecessor pot, by any method, and against what standard — the 20 mg is the only quantity in the
    paper and it is a purification artefact; (ii) the isolation recovery, or a spiked-recovery
    experiment, which would turn the 20 mg into a real yield; (iii) the ring-chain equilibrium in
    **water** rather than CDCl3, and if possible at temperature — this is what decides whether the
    `HMP` lump is defensible; (iv) any stability measurement of the isolated compound in aqueous
    buffer at 100-145 °C, which is the only thing that would put a number on `r_hmp_decay`; (v) the
    identity of compounds **4** and **5** in the compound numbering (Flags 5); (vi) the raw GC-MS
    trace of the ether extract, to see whether MFT was present at all in this 45 min pot.
11. **Registry gaps against `data/keys/compounds.yml`**: `2_methyl_3_furanthiol`,
    `2_furfurylthiol`, `furfural` and `mercapto_2_propanone` are present.
    **5-hydroxy-3-mercapto-2-pentanone is absent** — the very compound this paper identifies, and the
    one the engine routes its whole thiamine flux through. Also absent: **3-mercapto-2-pentanone**
    (engine `MP2P`), **2-mercapto-3-pentanone** (engine `MP3P`), **sulfurol**
    (5-(2-hydroxyethyl)-4-methylthiazole, one of only two major peaks in this extract and a species
    the engine does not carry at all), **thiamine** as a reactant (the registry has only
    `thiamine_availability`), **cysteine** and **xylose**. If a benchmark row is ever built on a
    thiamine pot, at least thiamine, cysteine and xylose would need keying.
12. **A species the engine does not carry appears here as a major product.** Sulfurol
    (5-(2-hydroxyethyl)-4-methylthiazole) is one of the **two** major GC-MS peaks in the ether
    extract, described as "a well-known thiamine degradation compound (6, 10, 12)". The engine's
    `r_thi_hmp` and `r_thi_mesh` send thiamine to HMP and methanethiol with the rest going to
    `FRAG_C`/`FRAG_N`; sulfurol is a named, abundant, nitrogen- and sulfur-bearing thiamine product
    that is currently inside that fragment lump. No rate exists for it here either, but it should be
    on the record as a candidate third thiamine fate rather than an anonymous fragment.
