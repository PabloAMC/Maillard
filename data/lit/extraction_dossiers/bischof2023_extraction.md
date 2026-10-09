# Bischof et al. 2023 — EXTRACTION (beef loin, 1H NMR, glucose, nucleotides and free amino acids over 28 days)

**Source on disk:** `data/articles/Bischof2023.pdf` (the publisher's open-access PDF). Table 1 read from
the PDF and checked by eye on 2026-09-14. **Every number in the earlier version of this dossier (written
the same day from an automated fetch of the HTML) was wrong**, not only the dry-/wet-aged columns it
flagged as suspect: the fetch returned the α-glucose row under "glucose" for one breed and the β-glucose
row for the other, and for IMP, inosine, hypoxanthine, leucine and methionine it returned numbers that
appear nowhere in Table 1. The table below is the PDF's. Written for the cultivated-tissue composition
box (`src/cultivated_tissue_invariance.py`), not for the core fit.

| field | value |
|---|---|
| Title | "NMR-based comparison of the metabolome of beef from Simmental and black-and-white young bulls during wet- and dry-aging" |
| Authors | G. Bischof, F. Witte, N. Terjung, E. Januschewski, V. Heinz, A. Juadjur, M. Gibis |
| Venue | European Food Research and Technology 2023, 249, 2113-2124 |
| DOI | 10.1007/s00217-023-04283-0 |

## 1. Methods

M. longissimus thoracis et lumborum from seven Simmental (S) and seven Black-and-White (BW) young bulls,
twenty months old; each loin split and dry-aged (D) and wet-aged (W), sampled at 0, 7, 14, 21 and 28 days
(the day-0 value is shared by D and W). 200 mg of meat, methanol / water extraction, 1H NMR, 30
metabolites quantified against standards. Table 1: "Concentration of metabolites (mean ± standard
deviation in µmol/g wet meat)", one row per breed × aging type. Glucose is printed as **two rows, α
glucose and β glucose** (the NMR anomers); free glucose is their sum.

## 2. Findings that matter (Table 1; µmol per g wet meat, mean ± SD; n = 7 bulls per breed)

Glucose, with the two anomer rows and their sum (SDs of the sum added linearly, an upper bound):

| row | S d0 | S D d28 | S W d28 | BW d0 | BW D d28 | BW W d28 |
|---|---|---|---|---|---|---|
| α glucose | 1.83 ± 0.91 | 2.71 ± 1.74 | 2.46 ± 1.36 | 2.44 ± 0.67 | 4.23 ± 0.63 | 3.53 ± 0.84 |
| β glucose | 2.60 ± 1.70 | 3.55 ± 2.38 | 3.21 ± 1.72 | 4.16 ± 1.48 | 5.78 ± 0.79 | 4.85 ± 1.37 |
| **glucose (α + β)** | **4.43 ± 2.61** | 6.26 ± 4.12 | 5.67 ± 3.08 | 6.60 ± 2.15 | **10.01 ± 1.42** | 8.38 ± 2.21 |
| glucose 6-phosphate | 0.98 ± 0.71 | 0.87 ± 0.61 | 0.80 ± 0.50 | 1.81 ± 0.77 | 1.47 ± 0.23 | 1.21 ± 0.34 |

Nucleotides and the amino acids the box or the sulfur lane cares about:

| compound | S d0 | S D d28 | S W d28 | BW d0 | BW D d28 | BW W d28 |
|---|---|---|---|---|---|---|
| IMP | 2.94 ± 0.32 | 1.27 ± 0.38 | 1.02 ± 0.23 | 3.19 ± 0.42 | 1.20 ± 0.22 | 1.00 ± 0.20 |
| inosine | 0.91 ± 0.11 | 1.35 ± 0.23 | 1.31 ± 0.35 | 1.12 ± 0.14 | 1.54 ± 0.16 | 1.30 ± 0.14 |
| hypoxanthine | 1.68 ± 0.41 | 4.37 ± 1.01 | 4.32 ± 1.03 | 1.64 ± 0.68 | 4.46 ± 1.17 | 4.81 ± 1.44 |
| leucine | 0.33 ± 0.05 | 1.30 ± 0.39 | 1.48 ± 0.37 | 0.29 ± 0.07 | 1.44 ± 0.50 | 1.53 ± 0.52 |
| isoleucine | 0.21 ± 0.03 | 0.75 ± 0.21 | 0.86 ± 0.22 | 0.18 ± 0.04 | 0.82 ± 0.27 | 0.91 ± 0.30 |
| valine | 0.39 ± 0.04 | 1.30 ± 0.37 | 1.52 ± 0.38 | 0.33 ± 0.07 | 1.35 ± 0.44 | 1.53 ± 0.51 |
| methionine | 1.18 ± 0.11 | 1.80 ± 0.27 | 1.85 ± 0.23 | 1.08 ± 0.10 | 1.72 ± 0.19 | 1.70 ± 0.23 |

Ribose, cysteine and thiamine are not among the 30 metabolites quantified. The methionine row is five
to twenty times the GC-MS and CE-TOFMS values for free methionine in beef (Koutsidis 2008a/b 0.05 to
0.40 mmol/kg; Muroya 2019 0.04 to 0.32 µmol/g); an NMR assignment at 1D resolution is the likely reason,
and the box does not use it.

## 3. What the repo takes

Beef-side ranges for the composition box, in mM in tissue water (µmol/g ÷ 0.75), spanning mean − SD to
mean + SD across both breeds and both aging types over 0 to 28 days:

- **glucose (free, α + β)**: 1.82 to 11.43 µmol/g → **2.4 to 15 mM** (means alone: 4.43 to 10.01 µmol/g,
  5.9 to 13.3 mM). The earlier dossier's "1.2 to 7.5 mM" was the anomer rows read as totals, which is
  why the box needed a cited review value as its upper corner; the primary sum covers Koutsidis 2008b's
  GC-MS range (9.8 to 13.7 mM) and Aliani 2013's cited point (11 mM) without it.
- **leucine (free)**: 0.22 to 2.05 µmol/g → 0.29 to 2.7 mM (earlier: 0.64 to 2.4, from the wrong rows).
- **IMP**: 0.80 to 3.61 µmol/g → 1.1 to 4.8 mM, inside Muroya 2019's 0.10 to 10 mM (earlier: 1.7 to 8.5,
  from the wrong rows).

Fourteen animals, two breeds, one laboratory. The day-0 glucose difference between breeds (4.4 vs
6.6 µmol/g) is of the same size as the ageing change within a breed; the Simmental glucose SDs are
about half the mean, so the low corner of the range is set by animal-to-animal spread, not by ageing.
