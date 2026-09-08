# Sägesser et al. 2024 — EXTRACTION (Roquette Nutralys F85M pea isolate and Solae Alpha 8 soy CONCENTRATE as high-moisture-extrusion references against two Auxenochlorella microalgal biomasses: composition, free amino acids, amino-acid classes, minerals, pH, NSI, zeta — every number a figure)
### Names two commercial reference powders to the lot and measures exactly what Programme 7 asks for on them (free amino acids by RP-HPLC, true protein, ash, water, minerals, pH) — and prints none of it as a number: Figs 1, 3 and 5 are raster images, the amino-acid table is in the companion paper Sägesser et al. 2023, which is not on disk.

**Source on disk:** `data/articles/Sagesser2024.pdf` (owner's download, 2026-09-08; 9 pages, no
supplement on disk; "Data will be made available on request"). Read from the scratchpad text
layer (`Sagesser2024.txt`, clean). pypdf on pp. 4-7 confirms that Fig. 1 (composition), Fig. 2
(extrudates), Fig. 3 (protein mass distribution, NSI, amino-acid classes and free amino acids),
Fig. 4 (interaction potential) and Fig. 5 (conductivity, zeta, minerals, pI) are each a single
embedded JPEG with no text layer — FIGURE-ONLY throughout; the one table (Table 1) is the dry-feed
recipe list. Repo status before this dossier: Nutralys F85M is already named in
`data/lit/binding_constants.yml` (Snel 2023 ketone-binding records, `protein_material: "Pea
protein isolate (PPI, Nutralys F85M, Roquette)"`); nothing on disk gives its composition. No
dossier cites this paper or Sägesser 2023.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Structurability of microalgae, soy and pea protein for extruded high-moisture meat analogues" |
| Authors | Corina Sägesser, Johanna Meta Kallfelz, Samy Boulos (co-corresponding), Joseph Dumpler, Lukas Böcker, Thomas Mair, Laura Nyström, Alexander Mathys (corresponding) — Sustainable Food Processing Laboratory and Laboratory of Food Biochemistry, ETH Zurich |
| Venue | Food Hydrocolloids 156 (2024) 110290; received 8 Mar 2024, revised 10 May, accepted 10 Jun, online 10 Jun 2024; open access CC BY |
| DOI | 10.1016/j.foodhyd.2024.110290 |
| Products | **PPI** = "pea protein isolate Nutralys F85M from Roquette (Lestem [sic; Lestrem], France; lot code: W007T)"; **SPC** = "Soy protein concentrate Alpha 8 from Solae LLC (Missouri, USA; lot code: R880002119)" — a CONCENTRATE, not an isolate (flag 1); **MAB1** = AlgaRich Vegi P55, Daesang (lot SCA24001) and **MAB2** = Golden Chlorella, Alver AG (lot C2B07-C318234), both heterotrophic Auxenochlorella protothecoides; maltodextrin 6 (Nutricia) and canola oil as fillers |
| Nitrogen factor | **sample-specific**, from Sägesser et al. 2023 (Bioresource Technology 390, 129849): protein = (total N - non-protein N) x a sample-specific factor; the values themselves are not printed here |
| Naming | "true protein" = the above; "NSI" = nitrogen solubility index; "native pH" = pH of the powder in deionised water (10 % w/v); S1-S3, SM1-SM3, P = extrudate recipes of Table 1 |
| Compound registry | none (no volatiles; amino-acid classes and minerals only) |

## 1. Why it matters

Programme 7 wants, per commercial isolate, the protein content with its nitrogen factor, the
protein-bound lysine, the free amino acids and the free sugars. This paper is the only one on
disk that (a) names the exact commercial pea isolate (Nutralys F85M, the product the repo's
Snel 2023 binding constants were measured on) and a commercial soy concentrate to the lot, (b)
determines true protein with a sample-specific nitrogen factor rather than 6.25, (c) measures
free amino acids explicitly (RP-HPLC after PITC derivatisation, norleucine internal standard,
triplicate) "to determine the number of amino acids bound in peptides and proteins", and (d)
measures ash, water (Karl Fischer), minerals (ICP-MS) and the powder's own pH. It is therefore
the right paper for the programme — and it delivers nothing quotable: every result is a bar or
pie in a raster figure, and the amino-acid profile is in Sägesser 2023. What can be carried is
qualitative and process-side: the two conventional powders have low free amino acids relative
to the microalgae, comparable cysteine, a native pH more than one unit above the microalgae
and about three units above their isoelectric points, and a defined high-moisture extrusion
history (80-140 C barrels, 50 C die, 2 kg/h, 500 rpm) that produced a 33 % protein extrudate
from pure PPI. The action item is to obtain Sägesser 2023 (flag 2).

## 2. Methods as they matter to a model

- **Composition (Fig. 1, FIGURE-ONLY).** Explicit protein and carbohydrate rather than N x 6.25 and
  by-difference: "The true protein content of the raw materials has been determined previously
  via multiplication of proteinaceous nitrogen, derived by the difference of total nitrogen minus
  non-protein nitrogen, with a sample-specific nitrogen-to-protein conversion factor by Sägesser
  et al. (2023)." Total carbohydrate = hydrolysed neutral sugars (HPAEC-PAD, Dionex ICS-5000+,
  two-step H2SO4 hydrolysis per Manns 2014; sorbitol internal standard; calibrated with
  arabinose, xylose, glucose, galactose, rhamnose, mannose, fructose) + uronic acids
  (colorimetric, Eurofins). **This is total structural carbohydrate after acid hydrolysis, not
  free sugar; free sugars were not measured** (flag 3). Fat by Weibull-Stoldt (4 M HCl 20 min
  100 C, Soxhlet; external lab AVS). Ash: 1 g at 500 C overnight. Water: Karl Fischer. Fig. 1
  caption: "carbohydrates include dietary fibres"; text: "Other unidentified components in
  Fig. 1 include organic acids, inorganic substances, nucleic acids, lignin, pigments, vitamins
  and metabolites."
- **Amino-acid profile.** "had been analysed chromatographically after hydrolysis by Sägesser
  et al. (2023)" — not in this paper. Shown in Fig. 3b as classes (aliphatic, aromatic, basic,
  acidic, amidic, hydroxyl, sulfur ...) of the protein-bound amino acids, "the difference between
  the sum of the embedded AA to 100 % are the free AA" — FIGURE-ONLY.
- **Free amino acids (verbatim core):** "free amino acids were quantified in this research by
  Reversed Phase High Performance Liquid Chromatography (RP-HPLC) to determine the number of
  amino acids bound in peptides and proteins. Prior to the analysis, intact cells were disrupted
  via High-Pressure-Homogenisation (HPH) of 0.1 % (w/v) sample suspensions ... three passages at
  1000 bar with cooling to 20 C ... filtered through a 0.45 um hydrophilic filter and 15 uL of
  0.5 mM NorLeu in 0.1 M HCl was added as internal standard to 500 uL filtrate before it was
  dried under vacuum at 35 C ... For derivatization, 200 uL of a 7:1:1:1 mixture of
  EtOH/water/NEt3/PITC was added ... 20 min at room temperature ... dissolved in 150 uL 5 mM
  phosphate buffer (pH 7.4, 5% ACN) ... Agilent 1200 series LC-system ... Pico Tag column ... as
  described by Kwanyuen and Burton (2010). Calibration curves were prepared with the amino acid
  standard mix AAS 18 as well as NorLeu, Gln, Asn and Trp". So: the free amino acids of a 1 g/L
  aqueous suspension, filtrate only (what dissolves at 20 C), corrected for NorLeu recovery.
  Results: Fig. 3b, FIGURE-ONLY; text: "MAB contain higher amounts of free amino acids as shown
  in Fig. 3b" (flag 4).
- **Protein size.** SDS-PAGE (NuPAGE 4-12 % Bis-Tris, MES) after the same HPH; Coomassie
  intensity -> cumulative mass distribution (Fig. 3a). Text: MAB proteins "20 % smaller in
  molecular weight" than PPI / SPC.
- **NSI.** From Sägesser 2023: N in the supernatant of a 10,000 g centrifugation over total N.
  Fig. 3a; text: MAB "18-68 % higher solubility" than PPI / SPC; HMEC has been done with NSI up
  to 64 %.
- **Zeta potential and pI.** 0.05 % w/v suspensions, Zetasizer Nano ZS; pI by HCl titration to
  ~0 mV. Fig. 5; text: "the native pH was approximately 3 pH units above the isoelectric point"
  for all; zeta "similar for PPI, SPC and MAB1".
- **pH and conductivity.** 10 % (w/v) suspensions in deionised water; pH found independent of
  dilution from 0.05 to 10 % ("presence of buffering agents"), so taken as the dough / extrudate
  pH. Fig. 5a; text: PPI and SPC "differed by more than 1 pH unit" above the MAB.
- **Minerals.** ICP-MS (iCap RQ, KED mode) after HNO3 microwave digestion; P, Na, K, Mg, Ca the
  relevant ones; chloride and non-phosphate anions not covered. Fig. 5b-f; text: "In PPI and
  MAB, P was present in excess to neutralize divalent cations. In SPC, similar amounts of
  divalent cations and P were present."
- **Interaction-potential arithmetic (Fig. 4).** Sum over sidechains of (count x bond energy):
  S-S 277 kJ/mol, ion pairs 38-73, cation-pi 12-25, amide bridges 12-14, H-bonds ~4 kJ/mol;
  hydrophobicity as summed sidechain mass. Declared as indicative only.
- **Extrusion.** Thermo Process 16 twin-screw, L/D 40, 8 barrels; powder barrel 1, water barrel 3,
  oil barrel 4; 2 kg/h, 500 rpm; barrels 2-8 at 80 / 80 / 120 / 130 / 130 / 140 / 140 C, end
  plate 140 C, cooling die 50 C (28 x 6 x 360 mm); vacuum-packed, 4 C. Feed water set to give
  **23.0 % protein in every extrudate** except P (pure PPI, **33 % protein**: "It was not possible
  to reach a protein content as low as 23 % using PPI alone since it cannot absorb the additional
  water"). SME and STE per recipe: Fig. 2, FIGURE-ONLY.
- **Statistics.** Triplicates, mean +/- SD, no tests (n = 3).

## 3. Tables re-typed

### Table 1. "Dry feed composition of tested recipes" (the only table)

| product name | dry feed composition |
|---|---|
| S1 | 100 % SPC |
| S2 | 88 % SPC + 12 % MD |
| S3 | 90 % SPC + 6 % MD + 4 % oil |
| SM1 | 50 % SPC + 50 % MAB1 |
| SM2 | 50 % SPC + 50 % MAB1-HPH |
| SM3 | 50 % SPC + 50 % MAB2 |
| P | 100 % PPI |

### Everything else: FIGURE-ONLY, with the sentences the text prints about PPI and SPC

| figure | content | what the text says about PPI / SPC (verbatim or close) |
|---|---|---|
| Fig. 1 | composition bars: protein, fat, carbohydrate (incl. fibre), ash, water, other — four powders, triplicate SD | "Compared to SPC, PPI contains a comparable amount of fat, substantially more protein and less potential passive fillers than MAB." No number. |
| Fig. 2 | extrudate dry-matter composition, SME, STE, cutting force, tensile strength, yield strength, Young's modulus, photos | "PPI showed an even lower SME than this blend [SPC + MAB]"; extrudate P "had slightly lower strength, stiffness and number of fibres, but higher hardness" than S1. |
| Fig. 3a | cumulative protein mass distribution (SDS-PAGE); NSI | MAB proteins "20 % smaller"; NSI of MAB "18-68 % higher" than PPI / SPC. |
| Fig. 3b | amino-acid classes as % of total protein; free amino acids as the remainder to 100 % | "MAB contain higher amounts of free amino acids"; PPI and SPC "contain more aliphatic and aromatic but fewer basic and amidic amino acids" than MAB; "The amount of cysteine was comparable in all raw materials"; "MAB1 contained the highest number of cysteine residues in its proteins". |
| Fig. 4 | theoretical interaction potential (covalent, electrostatic, H-bond energies; hydrophobic sidechain mass) | "the covalent and electrostatic interaction potential of all raw materials is comparable"; MAB1 has 15 %, MAB2 50 % fewer hydrophobic sidechains than PPI or SPC. |
| Fig. 5a | conductivity of 10 % suspensions; zeta at native pH | PPI and SPC pH "higher" than MAB by "more than 1 pH unit"; zeta "similar for PPI, SPC and MAB1", MAB2 less negative; MAB conductivity higher, MAB1 "substantially higher". |
| Fig. 5b-f | ICP-MS minerals; pI and mineral pies per powder | "the total amount of minerals was comparable in all raw materials"; native pH ~3 units above pI for all; P in excess over divalent cations in PPI, balanced in SPC ("Salt bridges were only expected to exist in SPC extrudates"). |

Microalgae (one line): both A. protothecoides biomasses (~50 % protein by selection criterion,
10.2-10.5 x 10^9 intact cells per g dry matter, reduced 83 % to 1.7 x 10^9 by HPH) gave weaker,
less fibrous extrudates at 50:50 with SPC; the paper attributes this to smaller proteins, higher
NSI, more free amino acids, fewer hydrophobic sidechains, lower pH and higher ionic shielding,
not to fat, protein content or intact cells.

## 4. Numbers the repository can use

| product | quantity | value +/- sd | unit as printed | mmol per g protein | method | source | evidence class |
|---|---|---|---|---|---|---|---|
| Nutralys F85M (Roquette, lot W007T) | identity of the pea isolate on which Snel 2023's ketone binding constants (`binding_constants.yml`) were measured; the "F85" grade name implies a nominal ~85 % protein (manufacturer's naming, not a measurement here) | — | — | — | — | §2.1 | identity only |
| Alpha 8 (Solae, lot R880002119) | soy protein CONCENTRATE used as the soy reference | — | — | — | — | §2.1 | identity only |
| PPI, SPC | true protein content (sample-specific N factor); fat; total carbohydrate (after hydrolysis); ash; water | FIGURE-ONLY | g/100 g (Fig. 1 bars) | — | Dumas / Kjeldahl-type N with NPN subtraction (2023); Weibull-Stoldt; HPAEC-PAD + uronic acids; 500 C; Karl Fischer | Fig. 1 | figure_only |
| PPI, SPC | lysine, arginine, cysteine, methionine per g protein | NOT IN THIS PAPER (Sägesser 2023); classes only, FIGURE-ONLY | % of total protein by class | — | chromatography after hydrolysis (2023) | Fig. 3b | figure_only / external |
| PPI, SPC | free amino acids (each and total) | FIGURE-ONLY; qualitatively LOWER than in both microalgal biomasses | % of total protein (remainder to 100 % in Fig. 3b) | — | RP-HPLC, PITC, Pico-Tag, NorLeu IS, 0.1 % w/v HPH suspension filtrate, triplicate | Fig. 3b; §3.3 | figure_only (direction stated in text) |
| PPI, SPC | free sugars | NOT MEASURED (the sugar profile is total carbohydrate after H2SO4 hydrolysis) | — | — | — | §2.2 | — |
| PPI, SPC | cysteine | "comparable in all raw materials" (MAB1 highest) | — | — | 2023 profile | §3.4, §3.5 | level_only, qualitative |
| PPI, SPC | pH of a 10 % (w/v) suspension in deionised water | FIGURE-ONLY; > 1 pH unit above the microalgae; ~3 units above the powder's pI; independent of dilution 0.05-10 % | — | — | pH meter | Fig. 5a; §2.4, §3.4 | figure_only (relations stated in text) |
| PPI, SPC | isoelectric point | FIGURE-ONLY | pH | — | zeta titration with HCl | Fig. 5c-f | figure_only |
| PPI, SPC | NSI | FIGURE-ONLY; MAB 18-68 % higher | % | — | 10,000 g supernatant N / total N (2023) | Fig. 3a | figure_only |
| PPI, SPC | minerals P, Na, K, Mg, Ca | FIGURE-ONLY; P > divalent cations in PPI, P ~ divalent cations in SPC | mg/g (Fig. 5b) | — | ICP-MS | Fig. 5b-f | figure_only |
| PPI (recipe P) | protein content of the pure-PPI high-moisture extrudate | 33 | % (wet) | — | recipe arithmetic | §2.5 | measured / declared; PPI "cannot absorb the additional water" needed for 23 % |
| SPC recipes (S1-S3, SM1-SM3) | protein content of the extrudates | 23.0 | % (wet) | — | recipe | §2.5 | declared |
| all recipes | extrusion heat history | barrels 80 / 80 / 120 / 130 / 130 / 140 / 140 C, end plate 140 C, die 50 C; 2 kg/h; 500 rpm; L/D 40 | — | — | — | §2.5 | declared; a usable time-temperature envelope for an HMEC-cook spec (residence time not printed) |
| MAB1, MAB2 | intact cells per g dry matter | 10.2 +/- 1.1 and 10.5 +/- 0.7 x 10^9; after HPH 1.7 +/- 0.5 x 10^9 (-83 %) | cells/g DM | — | Neubauer count | §3.1 | measured (microalgae only) |

Nothing in this table can enter `protein_matrices.yml`. The product identity is the one carry:
any future entry for a "pea isolate, Nutralys F85M" can bind Snel 2023's binding constants and
Sägesser 2023's composition to the same named material.

## 5. Flags

1. **The soy reference is a concentrate (Alpha 8, Solae), not an isolate.** Its protein-bound
   lysine per g PROTEIN is comparable to an isolate's, but per g POWDER it carries more
   carbohydrate (soluble sugars, oligosaccharides or fibre depending on the concentrate process)
   — which is why the paper needed maltodextrin and oil to macro-match SPC to the microalgae.
   Any use of "SPC" numbers as "soy isolate" numbers must be per g protein.
2. **The amino-acid table is in Sägesser et al. 2023, Bioresource Technology 390, 129849
   ("A novel approach for the protein determination in food-relevant microalgae").** That paper
   holds, for the same lots, the total amino-acid profile, the sample-specific nitrogen factor,
   the NPN and the NSI. It is not on disk; it is the item to fetch for Programme 7 — it would
   give Nutralys F85M's lysine per g true protein, which no other dossier has.
3. **"Sugar profile" here is total carbohydrate after two-step sulfuric-acid hydrolysis** (neutral
   sugars + uronic acids), i.e. starch, fibre and oligosaccharides counted as their monomers. It
   says nothing about free glucose / sucrose / maltose; the free-sugar charge of these two powders
   remains unmeasured (Jaeger 2023 measured it on other products: 0.19 / 0.06 g/100 g DM).
4. **Free amino acids were measured but are unreadable**: they appear only as the gap between the
   embedded-amino-acid classes and 100 % in Fig. 3b, a raster image, with the direction given in
   one sentence. The method also defines them as what passes a 0.45 um filter from a 1 g/L
   suspension at 20 C after 3 x 1000 bar homogenisation — a water-extractable free pool, fine for
   an isolate, and one that would include any small peptides' N-termini only if they were
   derivatised and co-eluted (PITC derivatises peptides too; the Pico-Tag calibration is for free
   amino acids).
5. **Every composition, pH, mineral, NSI, zeta and texture value is FIGURE-ONLY.** Figs 1-5 are
   embedded JPEGs; no supplementary file is on disk; "Data will be made available on request."
   Nothing was read off a figure for this dossier.
6. **The interaction-potential numbers (Fig. 4) are arithmetic on the amino-acid profile**, not
   measurements — the paper says so ("only provide an indication of reactivity. They may deviate
   significantly from actual values"). Do not cite them as disulfide or ion-pair counts.
7. **Protein basis differs from every other isolate dossier.** True protein = (N - NPN) x
   sample-specific factor, versus N x 6.25 in Jaeger 2023, Gorissen 2018, Gao 2020. For pea and
   soy the sample-specific factor is ~5.4-5.7, so the same powder reads 8-14 % lower in protein
   here than under 6.25; when Sägesser 2023's per-g-protein lysine is eventually entered it will
   be correspondingly HIGHER than a 6.25-based figure and must be labelled with its basis.
8. **"Lestem, France" is a misprint for Lestrem** (Roquette's site); the lot code W007T is as
   printed.
9. **Residence time and product temperature in the extruder are not printed**, so the heat history
   in §4 is a barrel-setpoint envelope, not a time-temperature profile; SME / STE per recipe are
   in Fig. 2 (FIGURE-ONLY).
