# Kocadağlı & Gökmen 2016 (JAFC, NaCl) — DUPLICATE. This is `kocadagli2016jafc_extraction.md`, already on disk.

### THIS IS NOT A THIRD PAPER. `data/articles/Kocada2016.pdf` is the JAFC glucose ± NaCl caramelization paper, and it already has a full 478-line dossier at `data/lit/extraction_dossiers/kocadagli2016jafc_extraction.md` — the source of the trunk's furanic and dicarbonyl constants. This file records the identity check and stops.

**Source on disk:** `data/articles/Kocada2016.pdf` (810,015 bytes; ACS "Just Accepted" manuscript,
38 pages).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Kocada2016.txt`, 2011 lines). The identity check below was made against the
ACS cover page and the first page of the manuscript body, both of which are clean in the text layer.

## 0. Identity

| field | value | verdict against the two existing dossiers |
|---|---|---|
| file on disk | `data/articles/Kocada2016.pdf` | **the same file** `kocadagli2016jafc_extraction.md` §0 names as its source |
| Title | "Effect of sodium chloride on α-dicarbonyl compounds and 5-hydroxymethyl-2-furfural formations from glucose under caramelization conditions – A multiresponse kinetic modelling approach" | **identical** to the title in `kocadagli2016jafc_extraction.md` §0 |
| Authors | Tolgahan Kocadağlı, Vural Gökmen (corresponding, vgokmen@hacettepe.edu.tr) — Food Quality and Safety (FoQuS) Research Group, Department of Food Engineering, Hacettepe University, 06800 Beytepe Campus, Ankara, Turkey | identical |
| Venue | J. Agric. Food Chem., **Just Accepted Manuscript**, web publication date 30 July 2016 | identical |
| DOI as printed | `10.1021/acs.jafc.6b01862` | **identical** |
| **Is this a third Kocadağlı 2016 paper?** | **NO.** | **It is `kocadagli2016jafc_extraction.md` under the file name the brief gave it.** |
| The other 2016 Kocadağlı dossier on disk | `kocadagli2016foodchem_extraction.md` = "Multiresponse Kinetic Modelling of Maillard Reaction and Caramelisation in a Heated Glucose/Wheat Flour System", Food Chemistry 211:892-902, DOI `10.1016/j.foodchem.2016.05.150`, file `data/articles/Kocadagli2016.pdf` | **a genuinely different paper**, different file, different DOI, different matrix |

**The file-name trap, which the existing dossier already documents and this check confirms.** The
two file names are the opposite way round from the naïve guess:

| file on disk | is actually | DOI | dossier |
|---|---|---|---|
| `Kocada2016.pdf` (shorter stem) | the **JAFC** glucose ± NaCl caramelization paper | `10.1021/acs.jafc.6b01862` | `kocadagli2016jafc_extraction.md` |
| `Kocadagli2016.pdf` (longer stem) | the **Food Chemistry** glucose/wheat-flour paper | `10.1016/j.foodchem.2016.05.150` | `kocadagli2016foodchem_extraction.md` |

`kocadagli2016jafc_extraction.md` §0 carries this as a "⚠️ WRONG-FILE WARNING, REPORT UPSTREAM",
and the warning was correct: a reader who assumed the longer stem was the JAFC paper would have
swapped them. **The brief that produced this dossier made the same assumption in the other
direction — it treated `Kocada2016.pdf` as a possible third paper — and the answer is that it is
the JAFC one.**

## 1. Why this matters, and why nothing more is written here

The brief asked for the identity to be settled first, and it is settled: **there is no third
Kocadağlı 2016 paper.** Writing a second extraction of the same PDF would put two independent
transcriptions of the same 102-cell Table 1 into `data/lit/extraction_dossiers/`, which is exactly
the way a transcription error becomes a disagreement between sources. The existing dossier already
holds:

- **§3, Table 1 re-typed in full** — 18 elementary steps × 3 temperatures × 2 systems (glucose and
  glucose-NaCl), each with its 95 % HPD, in `min⁻¹ × 10³`, with §3.1 listing the 13 cells whose
  HPD equals or exceeds the estimate and which must be refused.
- **§4, Table 2 re-typed in full** — the reparameterised Arrhenius fit, **T_b = 180 °C**, from which
  `src/kinetic_core/parameters_furanic.py` takes `k_glc_tdg`, `k_tdg_ddg`, `k_ddg_hmf`,
  `k_fru_int`, `k_int_hmf`, `k_fru_odg` and `k_tdg_mgo` (steps 3, 4, 5, 6, 7, 8 and 11), and
  `src/kinetic_core/parameters_dicarbonyl.py` takes `k_glc_g`, `k_g_go`, `k_odg_da`, `k_go_sink`
  and `k_da_sink` (steps 9, 10, 12, 15 and 17) — the whole B7 and B13 constant set.
- **§5**, an independent three-point Arrhenius refit of the paper's own Table 1 against its
  published Table 2, with the discrepancies named step by step.
- **§6**, the NaCl effect re-derived as within-study ratios, and §6.1 the mole conversions to HMF.
- **§8**, the verified negatives, and **§9**, a per-item USE / USE-Q / REFUSE verdict.

Nothing in the identity check contradicts any of it.

## 2. The one substantive point the brief raised, answered

`results/validation/kinetic_core_b21_prereg.md` records that **the trunk's α-dicarbonyl levels in
WATER are an open question**, and the brief asked whether any α-dicarbonyl concentration this paper
prints would be high-value against that. **It would not, and the reason is the matrix.**

- The system is **not aqueous**. Glucose (0.1 M) and glucose-NaCl (0.1 M each) solutions were made
  up in water only as a pipetting convenience — 0.5 mL per glass tube — and then **frozen at −80 °C
  and freeze-dried** before heating, expressly "to observe caramelization conditions during
  heating". The authors state that freeze-drying "led glucose and glucose-NaCl systems to two
  different types of amorphous states, **which was not characterized in this study**". The pot that
  is heated is a dry amorphous glass at 160, 180 and 200 °C for up to 30 min in PTFE-sealed tubes,
  in duplicate. There is **no water activity, no moisture content and no pH** anywhere in the paper.
- The model's units are **µmol per tube**, an absolute amount and not a concentration — the paper
  says so ("The amount of reactants and products were expressed as μmol"). Without a volume, and
  in a solid, they cannot be turned into mmol/L.
- **Every α-dicarbonyl amount-versus-time datum is figure-only** (Figures 3 and 4). The text layer
  preserves the axis labels — "3-deoxyglucosone, µmol", "1-deoxyglucosone, µmol",
  "3,4-dideoxyglucosone, µmol", "glucosone, µmol", "glyoxal, µmol", "methylglyoxal, µmol",
  "diacetyl, µmol", "HMF, µmol", "glucose/fructose, µmol" — but not the curves. Per house rule they
  are not typed as numbers, and the existing dossier already records this in §8.
- Three of them (3,4-dideoxyglucosone, 1-deoxyglucosone, glucosone) are **semi-quantitated against
  3-deoxyglucosone's calibration curve**, which the existing dossier flags as
  `absolute_concentration: false`.
- The only absolute numbers the paper prints are the initial glucose charges — **47.1 ± 0.66 µmol**
  (glucose) and **56.4 ± 0.77 µmol** (glucose-NaCl), a 19.7 % mismatch the existing dossier §1
  already flags — and the HMF mole conversions of its §6.1.

So this paper is the **source** of the glass constants whose transfer to water B21 was written to
fix, not evidence about water. It cannot close the B21 question, and `kinetic_core_b13_prereg.md`
§5's wishlist item stands unchanged: what is needed is *an isothermal aqueous glucose-amine pot with
glucosone, glyoxal and diacetyl quantified against time*. On disk, the papers that speak to that are
`hamzalioglu2026_extraction.md` (aqueous, milk, 110-140 °C), `quan2020_extraction.md` and
`xia2022_extraction.md` (levels in water) — and, newly, `goncuoglu2016_extraction.md`, whose
glyoxal sink is measured at three temperatures in a roasted hazelnut and whose authors found **no
glucosone at all**.

## 3. Flags

1. **Do not add a second transcription of this paper.** `kocadagli2016jafc_extraction.md` is the
   record. If a number from it is ever doubted, re-read the PDF and amend that file in place rather
   than creating a rival dossier.
2. **The file-name mapping is a live hazard.** Anything that ingests by file-name stem alone will
   swap `Kocada2016.pdf` and `Kocadagli2016.pdf`. Both existing dossiers carry the warning; this
   file adds a third copy of it so that a search for "kocadagli2016nacl" also lands on it.
3. **This dossier deliberately contains no re-typed table, no section 4 kinetic table and no
   arithmetic**, because everything it would contain is already in
   `kocadagli2016jafc_extraction.md` and duplicating it would create a second authority for the same
   102 cells.
