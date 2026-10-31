# Davidek, Robert, Devaud, Arce Vera & Blank 2006b — EXTRACTION
### *"Sugar Fragmentation in the Maillard Reaction Cascade: Formation of Short-Chain Carboxylic Acids by a New Oxidative α-Dicarbonyl Cleavage Pathway"*
### Glucose + glycine (0.1 M each), 0.2 M phosphate pH 8, 90 °C 15 h under argon / air / 18O2 (acetic acid by isotope dilution; C4-aldonic acids by HPAEC); 17O-NMR of 1-DG, 3-DG and glucose + glycine at pH 10 in 17O-water; 13C / 18O labelling of erythronic, threonic and glyceric acid

**Read 2026-10-31 from `data/articles/davidek2006b.pdf`** (8 pp., JAFC 54:6677-6684; born-digital, clean text
layer). Ledger `davidek2006b_tab1_glc_gly_acetic_atmosphere` (`data/lit/measurement_ledger/trunk.yml`); keyed rows in
`data/lit/fit_rows/ph_acid_davidek_acids.yml`. Companion paper: `davidek2006a_extraction.md`.

**Legend.** **[V]** verbatim · **[M]** printed number · **[Q]** qualitative statement by the authors · **[C]** cited ·
**[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

Question (B7a-2, `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 11): is the deoxyosone
beta-cleavage to acetic acid amine-assisted, and how does acid formation depend on pH, amine and oxygen?

| need | what Davidek 2006b gives | tag |
|---|---|---|
| (a) amine-free vs amine | **Nothing quantitative.** Every quantified pot has glycine. The 17O-NMR runs (pH 10, 90 °C, 3 h) of 1-DG and 3-DG are described in Methods as "glycine (0.35 mmol) and D-glucose, 2, or 1, each 0.35 mmol", i.e. **with glycine** [!]; Results call only sample C "glycine-induced", so the text is ambiguous. NMR peaks only, no amounts | [M] / [!] |
| (b) pH | one quantified pH (8). The beta-cleavage is said to occur "in a wide pH range (pH 3-11)" and to be "favored in the presence of phosphate buffer" [C: Hayami 1961, Weenen 1998]; no data here | [C] |
| (c) labelling | 18O2 at pH 8: acetic acid 22.9 +- 0.3 unlabelled + 10.6 +- 0.2 mono-18O = 33.5 +- 0.5 mmol/L, so **~32 % carries O2-derived oxygen** [Z]; only one 18O per molecule. C4-aldonic acids (erythronic, threonic) carry glucose C-3..C-6 (C-1 of the acid = glucose C-3), 1/3 18O-labelled at the carboxyl; glyceric acid from glucose C-4..C-6, 1/3 labelled; lactic acid ~20 % labelled. No 18O exchange of acetic acid with water under the conditions (checked) | [M] / [Q] |
| (d) yields | acetic acid (mmol/L, 100 mM glucose): **argon 30.5 +- 0.4, air 32.3 +- 3.3, 18O2 33.5 +- 0.5** (Table 1), "about 32 mol %" independent of atmosphere; C4-aldonic acids ~0.7 / 2.2 / 5.1 mmol/L under argon / air / 18O2 | [M] |

**Answer.** This paper bears on B7a-2 only indirectly. (i) **Anoxic acetic acid = oxic acetic acid** (30.5 vs 32.3
vs 33.5 mmol/L): the hydrolytic (non-oxidative) route carries essentially all the yield when O2 is excluded, so the
engine is right to make DPO / 1-DG beta-cleavage oxygen-free [M]. (ii) Under O2 about one third of the acetic acid
comes through the oxidative alpha-cleavage of 1-DG (Baeyer-Villiger type, mixed anhydride; acetic + erythronic acid),
yet the total does not rise: the two routes compete for the same pool (the authors' reading) [Q/Z]. (iii) No amine
dependence can be read from it. Confidence that the paper contains no amine-free acid yield: ~95 % (all pages read;
the only amine-free candidates, the 17O-NMR runs, are unquantified and probably contain glycine).

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Tomas Davidek, Fabien Robert, Stephanie Devaud, Francia Arce Vera (Nestle Research Center, Lausanne), Imre Blank* (Nestle PTC Orbe) |
| Venue | *J. Agric. Food Chem.* **2006**, 54 (18), 6677-6684 |
| DOI | **10.1021/jf060668i** (printed p. 6677; Crossref 2026-10-31: same title, authors, 54(18):6677-6684) |
| Dates | received 9 Mar 2006; revised 1 Jul 2006; accepted 9 Jul 2006; published on web 9 Aug 2006 |
| Content | Tables 1 (acetic acid by atmosphere), 2-3 (MS fragments of aldonic / glyceric acid TMS derivatives); Figs. 1, 4, 5, 6 (schemes), 2 (HPAEC chromatogram), 3 (17O NMR spectra) |

## 3. Methods [M]

| item | as printed |
|---|---|
| Glucose / glycine | glycine 0.3 mmol + glucose 0.3 mmol in phosphate buffer (3 mL, 0.2 mol/L, pH 8) -> 0.1 M each; 10 mL flask; argon (3 vacuum / argon cycles, 10 min Ar bubbling) or air or 18O2 (93-94 %), or 10 % H2 18O; 90 °C, 15 h, atmospheric pressure; "at least in duplicate" |
| Acetic acid | 0.3 mL + [2H3]acetic acid 600 µg, pH 3 (H3PO4), ether extract (1 mL), GC-MS SIDA (method 1 of davidek2006a); each injected twice; CV <= 10 % |
| C4-aldonic acids | HPAEC (CarboPac PA1, NaOH gradient to 300 mM), pulsed electrochemical detection (printed as "electron capture detector" [!]), external calibration; 10x dilution |
| Identification | TMS derivatives, GC-MS (HP-PONA); 13C isotopomers [1,2-13C2], [3-13C], [6-13C]glucose under air |
| 17O NMR | glycine 0.35 mmol + glucose, 3-DG or 1-DG 0.35 mmol in 10 % 17O-water (3.5 mL), **pH 10 (NaOH)**, 90 °C 3 h; 48.82 MHz, 70-120 k scans; no acetic acid / water exchange at pH 10 |

## 4. Results [M]

Table 1 (90 °C, 15 h, pH 8): see section 1 (d). 18O2 MS of acetic acid: m/z 43 (100), 45 (95), 47 (45), 60 (50), 62 (35).
A second 18O2 run at **pH 6** is reported in Results only (m/z 60 (70), 62 (35)) [!: Methods mention pH 8 only].
17O NMR: acetic acid (δ 282 ppm) from 1-DG; formic acid (δ 289 ppm) from 3-DG; both in glucose + glycine; no 17O in
carbonyl or alcohol groups. Authors' conclusion [Q]: "at least two newly described reaction mechanisms": hydrolytic
beta-dicarbonyl cleavage as the major pathway, oxidative alpha-dicarbonyl cleavage "induced by oxidizing species" as
a minor one; "about 30-33 mol % acetic acid ... of which only 1/3 was formed by the oxidative reaction pathway";
hydrolytic alpha-dicarbonyl cleavage of 1-DG "can be ruled out" under their conditions (dilute aqueous, 90-120 °C,
pH 6-10).

## 5. Flags

- [!] The "1/3 oxidative" share comes from the 18O2 pot; argon gives 30.5 mmol/L with no oxidant, so in the absence
  of O2 the hydrolytic route compensates fully. A model that adds an O2-dependent acetic acid source on top of the
  hydrolytic one would over-predict oxic pots.
- [!] HPAEC detector named "electron capture detector" with a gold electrode and a pulse waveform: an amperometric
  (PAD) detector, misnamed in print.
- Fit table spot-check (2026-10-31): the three rows (argon 30.5 +- 0.4, air 32.3 +- 3.3, 18O2 33.5 +- 0.5 =
  22.9 + 10.6) and the pot conditions in `data/lit/fit_rows/ph_acid_davidek_acids.yml` match the PDF.
