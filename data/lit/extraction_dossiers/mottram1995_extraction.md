# Mottram, Madruga & Whitfield 1995 — EXTRACTION (2,3-butanedione and 2,3-pentanedione, 2.5 % v/v in ethanol, saturated with gaseous H2S for 40 min at -15 C, then warmed to room temperature, nitrogen-purged and held 30 min in air; GC-MS + GC-odour-port + trapped-peak 1H NMR of thirty sulfur compounds — an IDENTIFICATION paper with no quantification of any kind)

### THE STRUCTURAL CHARTER FOR THE MERCAPTOKETONE SPECIES AND FOR THE DISULFIDE SINK, AND IT SAYS THE SINK IS BUILT TOO NARROW: every thiol in a mixed pool pairs with every other, so a pot containing 2-methyl-3-furanthiol, 2-furfurylthiol and two mercaptopentanones makes at least ten distinct disulfides — while the engine's `ch_dimer_mft` and `ch_dimer_fft` can only make the two homodimers, and the mixed disulfides this paper reports are precisely the ones with the meaty odour.

**Source on disk:** `data/articles/mottram1995.pdf` (5 pp., J. Agric. Food Chem. 1995, 43 (1), 189-193).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/mottram1995.txt`, 388 lines) for the running text, which came through legible
apart from routine OCR damage to units and italics (noted where it matters). **Table 1's two data
columns — LRI and the MS/NMR block — were LOST by the text extractor for compounds 1-12 and the MS
block was lost for all thirty rows**, so page 3 (journal p. 191) was rendered at 250 and 500 dpi with
`pdftoppm` and Table 1 was read from the page image. It is re-typed in full in section 3 from that
render. Figure 1 (the structures of the thirty compounds) is an image and is figure-only. There is no
supplementary material. Repo status before this dossier: this paper is **not cited anywhere in
`src/kinetic_core/`** — neither `sulfur.py`, `parameters_sulfur.py`, `species_sulfur.py` nor
`parameters_dicarbonyl.py` names it — and it has no extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Some Novel Meatlike Aroma Compounds from the Reactions of Alkanediones with Hydrogen Sulfide and Furanthiols" |
| Authors | Donald S. Mottram (corresponding; fax +44 734 310080) and Marta S. Madruga — Department of Food Science and Technology, The University of Reading, Whiteknights, Reading RG6 2AP, United Kingdom; **Frank B. Whitfield** — CSIRO Food Research Laboratory, P.O. Box 52, North Ryde, New South Wales 2113, Australia |
| Venue | J. Agric. Food Chem. 1995, 43 (1), 189-193. Received for review 28 June 1994, accepted 11 October 1994; abstract published in Advance ACS Abstracts, 15 November 1994 |
| Article ID | `JF940343E` (no DOI printed on the page) |
| Title check | the task brief guessed "Mottram & Nobrega?" — **the authors are Mottram, Madruga and Whitfield**. The title is otherwise exactly as given in the brief. |
| Paper type | **Identification / structure elucidation / sensory description.** No kinetics, no quantification, no time course, no temperature series. The word "yield" does not appear; abundance is stated as "major", "large amount", "small amounts" — GC peak language throughout. |
| Underlying thesis | Madruga, M. S., "Studies on some factors affecting meat flavour formation", Ph.D. Thesis, The University of Reading, 1994 — the source of the isolations that prompted this work |
| Companions on disk | `whitfield1988_extraction.md`, `whitfield2001_extraction.md` (same third author; the Discussion cites Whitfield, Mottram & Shaw 1993, the **4-hydroxy-5-methyl-3(2H)-furanone = norfuraneol** + cysteine/H2S system, as the place compounds 7, 9-12 and 20-22 were also found), `mottram2002_extraction.md`, `mottram2002b_extraction.md`, `farmer1990_extraction.md`, `schutte1972_extraction.md` |

## 1. Why it matters

The engine's sulfur lane carries four mercapto-carbonyl species and two disulfides, and this paper is
the structural source for what most of them are:

| engine species (`species_sulfur.py`) | this paper's compound | how the engine makes it |
|---|---|---|
| `MP` 1-mercapto-2-propanone | — (the C3 analogue; not in this paper) | `r_mgo_mp`: methylglyoxal + H2S, second order (`k_mgo_mp`) |
| — **no species** for 3-mercapto-2-butanone | **1**, LRI 815 | **no route exists** |
| `MP2P` 3-mercapto-2-pentanone | **2**, LRI 899 | `r_hmp_mp2p`: from thiamine's HMP only (`k_hmp_mp2p`, first order) |
| `MP3P` 2-mercapto-3-pentanone | **3**, LRI 904 | `r_nf_mp3p`: norfuraneol + H2S (`k_nf_mp3p`, second order) |
| `MFTD` bis(2-methyl-3-furyl) disulfide | **19**, LRI 1547 | `ch_dimer_mft`: 2 MFT + oxidant |
| `FFTD` bis(2-furfuryl) disulfide | **24**, LRI 1701 | `ch_dimer_fft`: 2 FFT + oxidant |

**Finding 1 — the dicarbonyl-plus-sulfide route to the mercaptoketones is one step, and the engine
only has it at C3.** This paper makes 3-mercapto-2-butanone from 2,3-butanedione + H2S and **both**
mercaptopentanone isomers from 2,3-pentanedione + H2S, in dilute ethanol, in one operation. The
engine has exactly the same chemistry at C3 (`r_mgo_mp`, methylglyoxal + H2S -> 1-mercapto-2-propanone,
a `MEASURED_SULFUR` step sized on Hofmann 1998 T7) but reaches `MP2P` only through thiamine and
`MP3P` only through norfuraneol. There is **no `DA + H2S`** row in `sulfur.py`, and there cannot be
one as the module stands, because **`DA` is in `TRUNK_ONLY_KEYS`** (`src/kinetic_core/species.py`
line 234) and is therefore absent from the sulfur state vector altogether. `species.py`'s own note on
`DA` already declares the gap in as many words: *"the 3-mercapto-2-butanone precursor once a sulfur
wave adopts it"*. This paper is the chemistry that adoption would be based on.

**Finding 2 — the disulfide sink is combinatorial, and the engine's is not.** `ch_dimer_mft` makes
`MFTD` from two MFT; `ch_dimer_fft` makes `FFTD` from two FFT. There is no cross term. This paper
puts four thiols in one pot and gets, in addition to the three homodimers, **every mixed pair**:
compound **20** (MFT-mercaptobutanone), **21** and **22** (MFT-mercaptopentanones), **26**, **27**,
**28** (FFT-mercaptoketones), **10**, **11**, **12** (mercaptoketone-mercaptoketone crosses) and
**30** (MFT-FFT). And the odour verdict runs the wrong way for the engine's convenience: the
Abstract's own conclusion is that "disulfides containing the 2-methyl-3-furyl group had meaty aromas,
whereas those without this group were sulfurous or onion-like" — so the **mixed** MFT disulfides
(20, 21, 22, 30) are aroma-active meaty compounds, not aroma loss. That is the same argument
`species_sulfur.py` already makes for `MFTD` ("NOT AROMA LOSS ... 15.6x MORE POTENT than its own
monomer"), extended to a set of species the engine cannot form. A pot in which MFT is consumed by
pairing with a mercaptoketone instead of with itself will, in the current engine, show that MFT as
simply gone.

**Finding 3 — the disulfides appear in 30 minutes at room temperature after the H2S is purged
out.** The procedure removes excess H2S with nitrogen and then holds the ampoules at room temperature
for 30 min before analysis; the disulfides are then "the major components formed in the reaction".
No oxidant is added and no metal is reported. `THIOL_CHANNELS`'s `oxidative_dimerisation` entry
defends its oxidant gate with "METAL-FREE thiol autoxidation in water is NEGLIGIBLE -- <= 2e-6 1/s
(Ngamchuea 2016)", which would give a half-life of about four days, not thirty minutes. The tension is
real but it is **not** a refutation, for four reasons set out in Flags 4: this is ethanol not water,
the thiols are at ~0.2-0.3 mol/L rather than at aroma levels, the ampoules were opened and purged in
air, and nothing here is quantified so "major peak" is a peak-area statement. It belongs on the
record as a bound on the argument, not as a number.

**Finding 4 — a citation the HMP node should carry.** The Discussion states that
5-hydroxy-3-mercapto-2-pentanone "is the intermediate for a number of thiols including
2-methyl-4,5-dihydro-3-furanthiol and 2-methyl-3-furanthiol **as well as the mercaptoketones 1-3**"
(citing Güntert 1993b, Hartman 1984, van der Linde 1979). The engine's `r_hmp_mp2p` sends HMP to
`MP2P` (compound **2**) only. If this reading of the older literature is right, HMP also feeds
compound **1** (3-mercapto-2-butanone, which the engine has no species for) and compound **3**
(`MP3P`, which the engine sources from norfuraneol) — and the second of those would blunt the
isomer-split diagnostic `species_sulfur.py` rests on (Cerny 2007: `MP2P` 77-90 % thiamine-derived,
`MP3P` 94->95 % xylose-derived). This is a **citation inside a discussion**, not a result of this
paper, and it does not override an isotope measurement. It is flagged so that nobody meets it later
and mistakes it for new evidence.

What this paper does NOT give: any concentration, any yield, any rate, any barrier, any time course,
any temperature series, any pH, any odour threshold measured here, any peak area, any response
factor, any mass balance, and any measurement at a cooking temperature. Its highest temperature is
room temperature.

## 2. Methods as they matter to a model

- **The pots — seven of them, lettered A-H in Table 1's footnote.** Aliquots (2 mL) of a solution of
  **2,3-butanedione (Aldrich) in ethanol, 2.5 % v/v**, were placed in **5 mL glass ampoules in an
  ice-salt mixture (temperature -15 °C)**. Hydrogen sulfide from a lecture bottle was passed through
  the solutions for **40 min at a flow of approximately 30 mL/min** (the text layer garbles the unit
  to "mumin"; the render reads mL/min). **2-Methyl-3-furanthiol (50 µL, Aldrich)** was then added to
  one reaction mixture and **a similar amount of 2-furylmethanethiol (Aldrich)** to another; H2S was
  passed through each for a **further 5 min**. A third mixture received no further treatment. The
  ampoules were allowed to reach room temperature, **then nitrogen was blown through each reaction
  mixture until excess hydrogen sulfide was removed**, and finally the ampoules were **held at room
  temperature for 30 min before analysis**. A similar set was prepared with **2,3-pentanedione** and
  with a **1:1 mixture of 2,3-butanedione and 2,3-pentanedione**. Separately, a mixture of
  2-methyl-3-furanthiol and 2-furylmethanethiol (**50 µL each in 2 mL ethanol**) was allowed to stand
  at **room temperature for 1 h** before analysis — that is mixture **H**, and it contains no
  dicarbonyl and no added H2S.
- **The seven systems, as printed in Table 1's footnote.** A: 2,3-butanedione + H2S. B:
  2,3-pentanedione + H2S. C: 2,3-butanedione + 2,3-pentanedione + H2S. D: 2,3-butanedione +
  2-methyl-3-furanthiol + H2S. E: 2,3-pentanedione + 2-methyl-3-furanthiol + H2S. F: 2,3-butanedione
  + 2-furylmethanethiol + H2S. G: 2,3-pentanedione + 2-furylmethanethiol + H2S. H:
  2-methyl-3-furanthiol + 2-furylmethanethiol.
- **Solvent, not water.** Every reaction is in **ethanol**. There is no buffer, no pH, and no pH is
  reported or reportable. This is the single largest obstacle to using anything here quantitatively:
  H2S dissociation, thiol pKa, and the whole thiolate chemistry the engine's `thiolate` factor
  expresses do not carry over from ethanol to an aqueous pot.
- **Temperature.** **-15 °C during the sulfiding, room temperature (unstated, ~20-25 °C) for the
  30 min hold.** Nothing here is at a cooking temperature. No temperature is varied. **No rate or
  barrier of any kind can be extracted, and none may be transported to 100-145 °C.**
- **Gas chromatography.** Split/splitless injection, **1 µL** of each reaction mixture; fused-silica
  capillary **DB-5, 30 m × 0.32 mm i.d., 1 µm film** (J&W Scientific), in a Hewlett-Packard HP5890;
  helium at **2 mL/min**; oven **60 °C for 5 min, then 4 °C/min to 250 °C**. The column effluent was
  **split equally between a flame ionisation detector and an odour port**. A **C8-C20 n-alkane**
  solution was run before each new reaction mixture to give **linear retention indices (LRI)**.
- **Odour-port evaluation.** **Four individuals**, experienced in aroma evaluation, who "provided
  descriptions for each aroma detected". No intensity scale, no dilution series, no AEDA, no FD
  factors, **no threshold measured here**. The descriptions in Table 1 are free-text and are
  consensus-free — they are listed, not scored.
- **GC-MS.** Hewlett-Packard **HP5988A** mass spectrometer with an HP5890 GC and an HP Chemstation;
  same column and conditions as above; source **250 °C**; ionising voltage **70 eV**; scan range
  **m/z 29-290**, **1 scan/s**. Note the low-mass cutoff at 29 and the **high-mass cutoff at 290**:
  the trisulfide 15/16 molecular ion at m/z 266 is inside it, but anything above 290 would not have
  been seen at all.
- **NMR of trapped GC peaks.** Effluent split between FID and a collection port; when a component of
  interest eluted, a **100 µL Microcap (1.42 mm o.d.)** was pushed into a tight-fitting PTFE sleeve
  in the collection port; **the tube was surrounded by solid carbon dioxide** during collection.
  **Three collections into one tube, between 5 and 50 µg**, from separate GC analyses. One end sealed
  in a microflame, **7 µL of 99.96 % deuteriobenzene** added, contents concentrated at the sealed end
  by spinning in a centrifuge, open end sealed. Spectra on a **Bruker CXP 100** with a 1H microprobe
  head, the Microcap inserted into a PTFE sleeve in an inverted 5 mm NMR tube. **1000 to 12 000 scans
  (1-12 h)** depending on sample size. **Spectral analyses assumed first-order principles.** NMR was
  acknowledged to B. H. Kennett.
- **Replication.** **None reported.** One preparation per system; the only repetition is the three
  GC collections pooled per NMR tube.
- **Quantification.** **None anywhere.** No internal standard, no calibration, no response factors,
  no peak areas printed, no relative amounts printed. Abundance appears only as the words "major",
  "large amount", "small amounts", "smaller peaks", "a small but significant" ion.

## 3. Tables re-typed

There is exactly **one** table in this paper. It is re-typed here in full from the 250/500 dpi page
render (the text layer lost both data columns). Figure 1 is a structure drawing and is figure-only.

### Table 1. "Spectral Data, Linear Retention Indices (LRI), and Odor Descriptions for Sulfur Compounds Obtained from the Reaction of Alkanediones, 2-Methyl-3-furanthiol, 2-Furylmethanethiol, and Hydrogen Sulfide"

Column headers as printed: `no.` / `compound` / `LRI` / `spectral data: MS (m/z, %) and 1H NMR
(100 MHz, C6D6, TMS)` / `odor port description` / `reaction mixture a`. Molecular ions are printed in
**bold** in the original; they are marked **bold** here. Two LRI values separated by a slash mean the
two diastereoisomers were resolved as separate GC peaks. NMR rows appear only for the compounds that
were trapped; a blank NMR line means no spectrum was recorded, not a null result.

| no. | compound | LRI | MS (m/z, %) | 1H NMR (100 MHz, C6D6, TMS) | odour port description | reaction mixture |
|---|---|---|---|---|---|---|
| 1 | 3-mercapto-2-butanone | 815 | 43 (100), 61 (75), **104** (47), 60 (25), 35 (9), 105 (3), 106 (2) | δ 1.09 (d, J = 8 Hz, 3H), 1.40 (d, J = 10 Hz, 1H), 1.75 (s, 3H), 2.68 (m, 1H) | burnt bread, burnt cereal, burnt vegetables, sulfury, burnt hamburger | A, C, D, F |
| 2 | 3-mercapto-2-pentanone | 899 | 43 (100), 41 (53), 75 (50), 74 (38), **118** (18), 57 (16), 45 (14) | δ 0.72 (t, J = 7 Hz, 3H), 1.40 (d, J = 9 Hz, 1H), 1.73 (m, 2H), 1.95 (s, 3H), 2.60 (bm, 1H) | sulfury, hydrogen sulfide | B, C, E, G |
| 3 | 2-mercapto-3-pentanone | 904 | 57 (100), 61 (59), **118** (23), 43 (16), 41 (9) | δ 0.97 (t, J = 3 Hz, 3H), 1.17 (d, J = 8 Hz, 3H), 1.40 (d, J = 10 Hz, 1H), 2.12 (m, 2H), 2.60 (bm, 1H) | sulfury, rotten meat, hydrogen sulfide | B, C, E, G |
| 4 | bis(1-methyl-2-oxopropyl) sulfide | 1359 | 43 (100), **174** (40), 87 (28), 59 (27), 103 (17), 115 (12), 128 (5), 71 (6), 175 (3), 176 (2) | — | chopped onions, sulfury, rotten egg, cooked cabbage | A, C |
| 5 | bis(1-methyl-2-oxobutyl) sulfide | 1460 | 57 (100), 43 (33), 113 (26), 153 (18), **202** (14), 164 (13), 135 (9) | — | sulfury, cooked vegetables | B, C |
| 6 | 2-[(1-methyl-2-oxopropyl)thio]-3-pentanone | 1498 | 43 (100), 57 (23), 117 (5), 71 (4), 85 (2), **188** (2), 129 (1), 142 (1), 159 (1) | — | sulfury, freshly cut spring onion | C |
| 7 | bis(1-methyl-2-oxopropyl) disulfide | 1476/1489 | 43 (100), 71 (18), 72 (16), **206** (14), 103 (17), 59 (15), 93 (10), 119 (9), 163 (8), 207 (2), 208 (1) | δ 1.08 (d, J = 8 Hz, 6H), 1.86 (s, 3H), 1.90 (s, 3H), 3.01 (m, 2H) | fried onion, chopped onion | A, C, D, F |
| 8 | bis(1-ethyl-2-oxopropyl) disulfide | 1606/1615 | 43 (100), 85 (36), 73 (21), **234** (15), 107 (15), 117 (14), 191 (7) | δ 0.70 (t, J = 7 Hz, 3H), 0.72 (t, J = 7 Hz, 3H), 1.71 (bm, 4H), 1.92 (s, 6H), 2.96 (m, 2H) | sulfury, burnt onion | B, C, E, G |
| 9 | bis(1-methyl-2-oxobutyl) disulfide | 1640/1649 | 57 (100), 117 (27), 43 (16), 85 (15), **234** (11), 175 (6), 133 (6) | δ 0.97 (t, J = 9 Hz, 3H), 1.00 (t, J = 9 Hz, 3H), 1.16 (d, J = 8 Hz, 3H), 1.18 (d, J = 7 Hz, 3H), 2.25 (bm, 4H), 3.09 (bm, 2H) | roasted, sulfury, hydrogen sulfide | B, C, E, G |
| 10 | 3-[(1-methyl-2-oxobutyl)dithio]-2-pentanone | 1629/1631 | 57 (100), 43 (96), 85 (42), 117 (27), **234** (21), 177 (8), 191 (4), 235 (2), 236 (2) | δ 0.71 (t, J = 7 Hz, 3H), 0.98-1.02 (t, J = 8 Hz, 3H), 1.18 (d, J = 9 Hz, 3H), 1.71 (bm, 2H), 1.92 (s, 3H), 2.25 (bm, 2H), 2.96 (bm, 1H), 3.09 (bm, 1H) | rotten egg, burnt, onion | B, C, E, G |
| 11 | 3-[(1-methyl-2-oxopropyl)dithio]-2-pentanone | 1544/1550 | 43 (100), 71 (10), 57 (9), 85 (7), **220** (5), 101 (3), 117 (3), 129 (20), 143 (2) — **the "129 (20)" is printed exactly so and breaks the row's descending order; see Flags 2** | δ 0.70-0.72 (t, J = 7 Hz, 3H), 1.10 (d, J = 8 Hz, 3H), 1.70 (bm, 2H), 1.88 (s, 3H), 1.93 (s, 3H), 2.96 (bm, 1H), 3.10 (bm, 1H) | oniony, freshly cut spring onion, sulfury | C |
| 12 | 2-[(1-methyl-2-oxopropyl)dithio]-3-pentanone | 1562/1567 | 57 (100), 43 (95), 59 (19), 71 (13), 103 (11), 85 (10), 117 (10), **220** (7), 163 (4), 129 (3), 145 (2), 177 (2) | δ 0.96-1.00 (t, J = 8 Hz, 3H), 1.08-1.11 (d, J = 7 Hz, 3H), 1.18 (d, J = 8 Hz, 3H), 1.88 (s, 3H), 2.26 (bm, 2H), 2.96 (bm, 1H), 3.10 (bm, 1H) | roast, chopped onion, sulfury | C |
| 13 | bis(1-methyl-2-oxopropyl) trisulfide | 1731 | 43 (100), 59 (18), 93 (14), **238** (10), 71 (8), 195 (4), 103 (3), 123 (2), 240 (1) | — | freshly cut onion | A, C |
| 14 | bis(1-ethyl-2-oxopropyl) trisulfide | 1877 | 43 (100), 85 (15), 105 (10), 57 (11), 73 (11), 149 (8), **266** (7), 117 (7) | — | sulfury, oniony, rotten egg | B, C |
| 15 | bis(1-methyl-2-oxobutyl) trisulfide | 1899 | 57 (100), 43 (12), **266** (7), 149 (6), 117 (6), 85 (6) | — | sulfury, burnt, hydrogen sulfide | B, C |
| 16 | 3-[(1-methyl-2-oxobutyl)trithio]-2-pentanone | 1888 | 57 (100), 43 (85), 149 (16), **266** (14), 117 (10), 85 (13), 73 (10), 105 (8) | — | sulfury, burnt, cooked cabbage | B, C |
| 17 | 3-[(1-methyl-2-oxopropyl)trithio]-2-pentanone | 1805 | 43 (100), 59 (8), 71 (6), 85 (6), 135 (6), 105 (4), **252** (4), 117 (4), 209 (3), 149 (2), 167 (2) | — | sulfury, freshly cut onion, fried | C |
| 18 | 2-[(1-methyl-2-oxopropyl)trithio]-3-pentanone | 1820 | 57 (100), 43 (67), 59 (18), 85 (7), 117 (7), 135 (6), **252** (4), 149 (3), 177 (1), 209 (1) | — | fried onions, sulfury, freshly cut spring onion | C |
| 19 | **bis(2-methyl-3-furyl) disulfide** | 1547 | 113 (100), **226** (44), 43 (30), 114 (14), 69 (13), 85 (12), 45 (13), 115 (7), 155 (5), 227 (5), 228 (4), 183 (2) | — | **beefy, meaty, meat soup** | D, E, H |
| 20 | 3-[(2-methyl-3-furyl)dithio]-2-butanone | 1510 | 113 (100), 43 (79), **216** (36), 114 (29), 59 (18), 45 (14), 81 (13), 173 (11), 217 (6), 218 (5) | — | **meaty, burnt meat, sulfury** | D |
| 21 | 3-[(2-methyl-3-furyl)dithio]-2-pentanone | 1580 | 43 (100), 113 (92), **230** (43), 114 (34), 145 (27), 85 (15), 81 (14), 187 (13), 231 (6), 232 (5) | — | **meaty, boiled meat, roast beef, sulfury** | E |
| 22 | 2-[(2-methyl-3-furyl)dithio]-3-pentanone | 1593 | 113 (100), 57 (78), **230** (37), 114 (32), 43 (29), 173 (14), 85 (12), 145 (7), 231 (5), 232 (4) | — | **roast meat, boiled meat** | E |
| 23 | bis(2-furylmethyl) sulfide | 1416 | 81 (100), 53 (15), **194** (12), 113 (8), 43 (7), 195 (2) | — | overcooked stew, burnt, sulfury | F, G, H |
| 24 | **bis(2-furylmethyl) disulfide** | 1701 | 81 (100), 53 (14), **226** (7), 69 (5), 112 (4), 193 (2) | — | spring onion, freshly cut onion, fried onion, roasted, rubbery | F, G, H |
| 25 | bis(2-furylmethyl) trisulfide | 1879 | 81 (100), 53 (20), 43 (10), 45 (8), 161 (6), 113 (3), **258** (1) | — | sulfury, hydrogen sulfide | F, G, H |
| 26 | 3-[(2-furylmethyl)dithio]-2-butanone | 1589 | 81 (100), 43 (14), 53 (13), 45 (5), **216** (2), 185 (2) | — | sulfury, onion | F |
| 27 | 3-[(2-furylmethyl)dithio]-2-pentanone | 1671 | 81 (100), 43 (14), 53 (11), **230** (2), 113 (1), 196 (1) | — | rubbery, oniony | G |
| 28 | 2-[(2-furylmethyl)dithio]-3-pentanone | 1686 | 81 (100), 53 (13), 45 (4), **230** (2), 113 (1) | — | roasted, oniony | G |
| 29 | 2-methyl-3-[(2-furylmethyl)thio]furan | 1501 | 81 (100), 53 (21), **194** (17), 113 (13), 45 (11), 126 (7), 195 (1) | — | bland meat | H |
| 30 | **2-methyl-3-[(2-furylmethyl)dithio]furan** | 1649 | 81 (100), 113 (18), 53 (18), 43 (14), **226** (9), 45 (12), 85 (11), 162 (1) | — | **meat, burnt meat, roast meat, roast coffee** | H |

Footnote a, printed in full: "Reaction mixtures in which compounds were found: A, 2,3-butanedione +
hydrogen sulfide; B, 2,3-pentanedione + hydrogen sulfide; C, 2,3-butanedione + 2,3-pentanedione +
hydrogen sulfide; D, 2,3-butanedione + 2-methyl-3-furanthiol + hydrogen sulfide; E, 2,3-pentanedione
+ 2-methyl-3-furanthiol + hydrogen sulfide; F, 2,3-butanedione + 2-furylmethanethiol + hydrogen
sulfide; G, 2,3-pentanedione + 2-furylmethanethiol + hydrogen sulfide; H, 2-methyl-3-furanthiol +
2-furylmethanethiol."

### Statements of relative abundance in the running text (words, never numbers)

| statement | where |
|---|---|
| butanedione + H2S: "three major GC peaks corresponding to 3-mercapto-2-butanone (1) and the corresponding disulfide, bis(1-methyl-2-oxopropyl) disulfide (7)" — the disulfide's two diastereoisomers resolved into separate peaks | Results |
| "Small amounts of the corresponding monosulfide 4 and trisulfides 13 were also found, although diastereoisomers were not evident" | Results |
| pentanedione + H2S: "Symmetrical disulfides 8 and 9 and the unsymmetrical disulfide 10 ... were **the major components formed in the reaction**", each resolved into two diastereoisomeric pairs; "the corresponding three trisulfides 14-16 were also formed, but **only one related sulfide**, bis(1-methyl-2-oxobutyl) sulfide (5), was produced" | Results |
| MFT + pentanedione: "A **large amount** of compound 19 was also formed, together with the two symmetrical disulfides 8 and 9, and the mixed disulfide 10 ... **Related sulfides and trisulfides were not found.**" | Results |
| FFT + butanedione: "The chromatogram ... showed other **large peaks** corresponding to bis(2-furylmethyl) disulfide (24) and compound 7, with **smaller peaks** containing bis(2-furylmethyl) sulfide (23) and its trisulfide homologue 25" | Results |
| mixture H: "2-methyl-3-[(2-furylmethyl)dithio]furan (30) and **a small amount** of the corresponding monosulfide 29 were produced together with the two symmetrical disulfides 19 and 24" | Results |
| the sensory conclusion | "The aromas detected for the volatiles containing the 2-methyl-3-furyl moieties were described as 'meaty, beefy, boiled meat, roast meat'. However, 'sulfury, burnt, oniony, rubbery' notes were mainly used to describe the compounds containing the 2-furylmethyl group" |
| novelty claim | "Compounds 4-6, 8, 13-18, and 26-29 are reported here for the first time." |
| where the mercaptoketone disulfides had also been seen | "We have recently found compounds 7, 9-12, and 20-22 in a heated model system containing **4-hydroxy-5-methyl-3(2H)-furanone and cysteine or hydrogen sulfide** (Whitfield et al., 1993)" — i.e. in a **norfuraneol** system |
| Güntert's thresholds (cited, not measured here) | "The odor threshold values of such compounds were found to be in the **low micrograms per kilogram range**, and most of those containing the 2-methyl-3-furyl group had meaty aromas" (Güntert et al. 1993b) |

### Diagnostic fragment rules the paper states in words (useful for reading anyone else's spectra)

- m/z **43** (CH3CO+) marks compounds formed from the **3-mercapto-2-alkanones** (1 and 2), i.e. an
  acetyl group.
- m/z **57** (CH3CH2CO+) marks compounds formed from **2-mercapto-3-pentanone** (3), i.e. a propionyl.
- Unsymmetrical compounds from a mixture of the two mercaptopentanones "showed abundant ions at both
  m/z 57 and 43".
- m/z **113** marks compounds containing the **2-methyl-3-furyl** group.
- m/z **81** is the intense fragment for compounds containing the **2-furylmethyl** moiety, whose
  molecular ions are "very small".
- 3-mercapto-2-pentanone: base peak 43, major fragment **75** (C2H5CHSH+) from fission at C2.
  2-mercapto-3-pentanone: base peak **57**, strong ion **61** (CH3CHSH+) from fission at C3. Both have
  small molecular ions at m/z **118**, and "agreed with published spectra (Hartman et al., 1984)".
- "All showed molecular ions and M + 2 ions arising from the 34S isotope."

### Arithmetic on the printed charges (all mine, and all approximate)

**1. What is in a 2 mL ampoule.** 2.5 % v/v gives **50 µL** of dione in 2 mL. Converting to moles
needs densities the paper does **not** print, so these are mine with the literature density stated:
2,3-butanedione (ρ ≈ 0.990 g/mL, M = 86.09) → 49.5 mg → **0.575 mmol**, i.e. **~288 mmol/L**;
2,3-pentanedione (ρ ≈ 0.957 g/mL, M = 100.12) → 47.9 mg → **0.478 mmol**, i.e. **~239 mmol/L**.
The added thiols, 50 µL each and both C5H6OS with M = 114.17: 2-methyl-3-furanthiol (ρ ≈ 1.09 g/mL)
→ **~0.48 mmol, ~239 mmol/L**; 2-furylmethanethiol (ρ ≈ 1.13 g/mL) → **~0.49 mmol, ~247 mmol/L**.
**These are order-of-magnitude figures only** — the densities are not from this paper (Flags 3).

**2. The sulfide charge is a vast excess and is not a concentration.** 40 min at ~30 mL/min is
~1200 mL of H2S **passed through** the solution. At -15 °C and 1 atm the molar volume is
22.414 × 258.15/273.15 = 21.18 L/mol, so ~**57 mmol** of H2S was swept through against ~0.5-0.6 mmol
of dione — a ~100-fold molar excess **delivered**, of which only the dissolved fraction reacts and
that fraction is unknown. The further 5 min after the thiol addition delivers ~7 mmol more. **No
sulfide concentration can be computed**, and therefore no order in H2S can be read off this paper —
which is the reason its chemistry cannot be turned into a rate law even in principle.

**3. Molecular ions confirm the assignments the engine cares about, by my arithmetic.** 3-mercapto-2-butanone
C4H8OS = 104 (printed as the M+ of compound 1, bold at 47 %). Both mercaptopentanones C5H10OS = 118
(compounds 2 and 3, both bold) — which reproduces `species_sulfur.py`'s **118.20** for `MP3P` and
`MP2P`. The mercaptobutanone homodisulfide C8H14O2S2 = 206 (compound 7); the mercaptopentanone
homodisulfides C10H18O2S2 = 234 (compounds 8, 9, 10); bis(2-methyl-3-furyl) disulfide C10H10O2S2 =
226 (compound 19) and bis(2-furylmethyl) disulfide, the same formula, also **226** (compound 24) —
so **`MFTD` and `FFTD` are isobaric at nominal mass**, and are told apart only by their base peaks
(113 against 81) and their retention (1547 against 1701). Worth knowing before anyone reads a mass
chromatogram from a mixed pot.

**4. The combinatorics the engine does not have (mine).** With four thiols in one pot — MFT, FFT,
3-mercapto-2-pentanone and 2-mercapto-3-pentanone — the number of distinct disulfides is
4 + C(4,2) = **10**. This paper reports the homodimers 8, 9, 19, 24 and the crosses 10, 21, 22, 27,
28, 30 — that is ten, i.e. **the complete set**, plus the butanone series on top of it. The engine
forms 2 of the 10.

## 4. Kinetic numbers the repository can use

**There is no rate, no barrier, no yield and no concentration in this paper.** Every printed number
is a retention index, a mass-spectral relative intensity, an NMR shift or coupling, or a charge in
the recipe. This section records them in their correct classes so that a later wave does not mistake
a relative MS intensity for an abundance or an LRI for a measurement of anything chemical.

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Keyed: `2_3_butanedione` (with alias
"diacetyl"), `2_methyl_3_furanthiol`, `2_furfurylthiol`, `bis_2_methyl_3_furyl_disulfide`,
`mercapto_2_propanone`, `hydrogen_sulfide`, `dimethyl_disulfide`, `dimethyl_trisulfide`, `furfural`.
**Not keyed:** 2,3-pentanedione, 3-mercapto-2-butanone, 3-mercapto-2-pentanone, 2-mercapto-3-pentanone,
bis(2-furylmethyl) disulfide (the engine's `FFTD`), and **all 24 other compounds in Table 1**.

Every row below shares: 2 mL of a 2.5 % v/v alkanedione solution in **ethanol** (no water, no
buffer, no pH), 5 mL glass ampoule, gaseous H2S bubbled 40 min at ~30 mL/min at **-15 °C**, N2 purge,
then **30 min at room temperature in air**, single preparation, no replicate, no internal standard,
no quantification.

| quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **2,3-butanedione + H2S gives 3-mercapto-2-butanone** (compound 1) | — (presence) | — | ethanol, -15 °C then RT | **nothing is fitted in this paper** | Results; Table 1 row 1 | **threshold** (presence, above the GC-MS detection limit) |
| **2,3-pentanedione + H2S gives BOTH mercaptopentanone isomers** (2 and 3) | — (presence, both) | — | as above | — | Results; Table 1 rows 2, 3 | **threshold** |
| the two isomers 2 and 3 are formed **together** from one dione, and no split between them is stated | — | — | as above | — | Results ("gave a mixture of 3-mercapto-2-pentanone (2) and 2-mercapto-3-pentanone (3)") | **level_only** — **the paper prints NO ratio**; do not invent one |
| every thiol pair present forms its mixed disulfide (10 of 10 possible for the four thiols) | — (presence, complete set) | — | as above | — | Table 1 rows 10, 21, 22, 27, 28, 30 + homodimers 8, 9, 19, 24 | **threshold** (the structural finding; see section 3 arithmetic 4) |
| mono-, di- and **tri**sulfides all form; disulfides are "the major components" | — | — | as above | — | Results | **peak_area_only** — a chromatogram-eye statement, no areas printed |
| bis(2-methyl-3-furyl) disulfide (`MFTD`) forms whenever MFT is present | — (presence) | — | mixtures D, E, H | — | Table 1 row 19 | **threshold** |
| MFT + FFT alone, no dicarbonyl, RT, 1 h: gives the **mixed** disulfide 30 and the two homodimers 19 and 24 | — (presence) | — | ethanol, RT, 1 h, **no added oxidant, no H2S** | — | Results, mixture H | **threshold** — the cleanest statement in the paper that disulfides form fast at RT (Flags 4) |
| linear retention indices, DB-5, C8-C20 alkanes | 815 to 1899 (thirty values, section 3) | LRI units | DB-5, 30 m × 0.32 mm, 1 µm, He 2 mL/min, 60 °C (5 min) then 4 °C/min to 250 °C | — | Table 1, LRI column | **level_only** (an identification aid; not a chemical measurement) |
| mass-spectral relative intensities | as re-typed in section 3 | % of base peak | EI 70 eV, source 250 °C, scan m/z 29-290, 1 scan/s | — | Table 1, MS column | **level_only** — **relative intensities within one spectrum, NOT abundances between compounds** |
| 1H NMR shifts and couplings for compounds 1, 2, 3, 7, 8, 9, 10, 11, 12 | as re-typed in section 3 | ppm / Hz | 100 MHz, C6D6, TMS, on 5-50 µg of GC-trapped material | — | Table 1, NMR column | **level_only** |
| odour-port descriptions | free text (section 3) | — | four assessors, no intensity scale, no dilution series | — | Table 1, odour column | **level_only** — descriptions, **not thresholds and not intensities** |
| the sensory rule: 2-methyl-3-furyl disulfides are meaty, the rest are sulfurous/onion | — | — | as above | — | Abstract and Discussion | **level_only** (a qualitative verdict from four assessors) |
| odour thresholds "in the low micrograms per kilogram range" | — | µg/kg | — | — | Discussion, **citing Güntert et al. 1993b** | **level_only** and **not measured here** — do not attribute to this paper |
| Figure 1 (the thirty structures) | — | — | — | — | Fig. 1 | **figure_only** |

### Can anything here be put on the same basis as the sulfur lane's constants? No — and the reason is structural, not a matter of effort.

**(a) There is no time resolution and no temperature axis.** One sulfiding at -15 °C, one 30 min hold
at an unstated room temperature, one endpoint. Nothing constrains `k_dimer_mft`, `k_dimer_fft`,
`k_nf_mp3p`, `k_hmp_mp2p`, `k_thiol_decay` or any other constant, individually or as a ratio.

**(b) The solvent is ethanol.** The lane's pH machinery — `neutral_h2s` against `hs_anion` tags, the
`thiolate` factor on `ch_dimer_*`, `k_thiolate_loss` — is aqueous acid-base chemistry. It has no
meaning in ethanol, and a rate measured there could not be transported even if one existed.

**(c) The sulfide is a gas stream, not a pool.** ~57 mmol of H2S swept through 2 mL of solution
gives no concentration and no order. Every H2S-dependent step in `sulfur.py` is second order in
(carbonyl × H2S); this paper cannot test that, confirm it, or size it.

**(d) "Major peak" is not an abundance.** No areas, no response factors, no standard. The house rule
that a peak-area ratio is not a yield applies here twice over: there are not even peak areas, only
adjectives.

**What DOES transfer** is structure and topology: three mercapto-carbonyl identities with their
mass spectra and retention indices, the one-step dicarbonyl + sulfide route at C4 and C5, the
combinatorial disulfide set, the isobaric `MFTD`/`FFTD` warning, and the sensory verdict that mixed
2-methyl-3-furyl disulfides are meaty. Those licence **network edges and species records**, not
constants.

## 5. Flags

1. **Table 1's data columns were lost by the text extractor and were read from a page render.**
   `pdftotext -layout` returned the compound names, the LRI values for rows 13-30 and most odour
   descriptions, but **dropped the entire MS/NMR column and the LRI values for rows 1-12**. Section 3
   is transcribed from `pdftoppm` renders of journal page 191 at 250 dpi and, for one ambiguous row,
   500 dpi. Two odour descriptions were additionally garbled in the text layer ("su12xh$drogen" for
   row 2, "au1fu79 hy rogen sulfide" for row 3, "sulkfuakmt" for row 16) and are given in section 3
   as read from the render.
2. **One printed intensity is internally inconsistent and is transcribed as printed.** Row 11 reads
   "43 (100), 71 (10), 57 (9), 85 (7), **220** (5), 101 (3), 117 (3), 129 (20), 143 (2)". Every other
   row lists ions in descending intensity, so a value of 20 between two 3s and a 2 breaks the
   pattern; the 500 dpi render shows "(20)" unambiguously. It is either a typesetting error for "(2)"
   or a departure from the ordering. **I have not corrected it**, and nobody should silently read it
   as 2.
3. **All molar concentrations in this dossier are mine and rest on densities the paper does not
   print.** The recipe gives volumes (2 mL, 2.5 % v/v, 50 µL) and a gas flow, never a molarity. My
   ~288 / ~239 / ~239 / ~247 mmol/L for butanedione, pentanedione, MFT and FFT use literature
   densities from outside this paper and should be treated as order-of-magnitude. The **ratios** are
   sturdier than the absolutes.
4. **The room-temperature disulfide formation is a real observation but a weak constraint on the
   engine's oxidant gate.** Mixture H — two thiols in ethanol, one hour at room temperature, no H2S,
   no added oxidant — yields the mixed disulfide 30 and both homodimers; mixtures A-G give disulfides
   as "the major components" within 30 min after the H2S was purged out with nitrogen.
   `THIOL_CHANNELS`'s `oxidative_dimerisation` entry justifies making dimerisation first order in an
   explicit oxidant pool by citing metal-free autoxidation in **water** at <= 2e-6 1/s. Four reasons
   this paper does not overturn that: (i) the solvent is **ethanol**, not water; (ii) the thiols are
   at ~0.24 mol/L, four to six orders above aroma concentrations, and a second-order channel scales
   as the square; (iii) the ampoules were purged with a gas stream and handled in air, so dissolved
   O2 is present and unmeasured; (iv) trace metals in the reagents are neither excluded nor measured,
   and nothing is quantified, so "major peak" cannot be converted into a rate. **Record it as a
   qualitative caution on the gate, not as evidence against it, and do not derive any constant from
   it.**
5. **`DA` cannot currently reach the sulfur lane at all.** `species.py` line 234 puts `DA` in
   `TRUNK_ONLY_KEYS`, so it is absent from the sulfur state vector; `network.py` makes it
   (`r_odg_da`, Kocadagli step 12) and sinks it at a rate the source gives as zero (`r_da_sink`).
   Adopting this paper's C4 chemistry would mean moving `DA` out of `TRUNK_ONLY_KEYS`, adding a
   3-mercapto-2-butanone species, and adding a `DA + H2S` step — three structural changes, none of
   which this paper supplies a rate for. **The C3 analogue `k_mgo_mp` is the only sized member of the
   family** (Hofmann 1998 T7), and borrowing it across two carbon numbers is not licensed.
6. **The disulfide sink is structurally incomplete and the missing members are aroma-active.** The
   engine can form `MFTD` and `FFTD` only. This paper reports six mixed disulfides involving MFT or
   FFT (20, 21, 22, 26, 27, 28) plus the MFT-FFT cross (30), and the four containing the
   2-methyl-3-furyl group are described as meaty, boiled meat, roast meat, roast coffee. In a pot
   containing mercaptoketones, MFT that leaves as a mixed disulfide is scored by the engine as MFT
   destroyed, and its aroma contribution is lost from the prediction. **No number here fixes this**;
   it is a topology flag.
7. **`MFTD` and `FFTD` are isobaric (both M+ 226)** and differ in the table only by base peak (113
   vs 81) and LRI (1547 vs 1701). Any benchmark row built from a single-ion trace at m/z 226 in a pot
   containing both thiols would conflate them.
8. **The HMP-to-mercaptoketone claim in the Discussion is a citation, not a result.** "An important
   thiamin degradation product is 5-hydroxy-3-mercapto-2-pentanone, and this very reactive compound
   is the intermediate for a number of thiols including 2-methyl-4,5-dihydro-3-furanthiol and
   2-methyl-3-furanthiol as well as the mercaptoketones 1-3" cites Güntert 1993b, Hartman 1984 and
   van der Linde 1979. If taken at face value it would give HMP a route to `MP3P` as well as to
   `MP2P`, which would compete with `r_nf_mp3p` and blunt the Cerny 2007 isomer diagnostic. **It is
   third-hand here and must not be adopted from this paper.** See `cerny2008_extraction.md` for the
   HMP identification itself.
9. **No replication, no error bar, no n, no standard, anywhere in the paper.** One preparation per
   system; the odour panel is four people producing free text.
10. **Nothing is at a cooking temperature.** -15 °C for the sulfiding and room temperature for the
    hold. Any use of this paper to justify a rate, a barrier, or a branch ratio at 100-145 °C would
    be an extrapolation across 120-160 °C from a measurement that is not a rate to begin with.
11. **What this paper does not contain**: any concentration; any yield; any peak area; any response
    factor; any rate or barrier; any time course; any temperature or pH variation; any aqueous
    system; any threshold measured here; any isomer ratio between compounds 2 and 3; any mass
    balance; any statement of how much of the dione was consumed; any supplementary material.
12. **What to request from the authors**: (i) the FID peak areas behind the "major"/"small" language,
    which would at least give a `peak_area_only` isomer ratio for compounds 2 and 3 from one dione —
    the single most useful missing number, since the engine sources those two isomers from two
    entirely different precursors; (ii) whether the 30 min room-temperature hold was in air or under
    nitrogen, and whether the reagents were metal-screened (Flags 4); (iii) the room temperature
    itself; (iv) any repeat of mixture H with a shorter hold, which would bound the room-temperature
    dimerisation timescale; (v) whether the mercaptoketone was ever observed to survive at all once
    the H2S was purged, i.e. the thiol-to-disulfide split at the endpoint.
13. **Registry gaps against `data/keys/compounds.yml`**: `2_3_butanedione`,
    `2_methyl_3_furanthiol`, `2_furfurylthiol`, `bis_2_methyl_3_furyl_disulfide` and
    `hydrogen_sulfide` are present. **Absent: 2,3-pentanedione** (the C5 dione, which has no engine
    species either), **3-mercapto-2-butanone** (engine has no species), **3-mercapto-2-pentanone**
    (engine `MP2P`), **2-mercapto-3-pentanone** (engine `MP3P`), **bis(2-furylmethyl) disulfide**
    (engine `FFTD` — the only one of the engine's two disulfides that is not keyed), and the
    twenty-four mixed and poly-sulfides. If the disulfide sink is ever scored against a real
    chromatogram, at least `FFTD` and the two mercaptopentanones need keying.
