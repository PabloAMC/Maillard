# Schieberle 2005 (Ann. N.Y. Acad. Sci. 1043:236-248) — EXTRACTION (file stem `schieberle2005b`)
### *"The Carbon Module Labeling (CAMOLA) Technique: A Useful Tool for Identifying Transient Intermediates in the Formation of Maillard-Type Target Molecules"*
### **Not the paper the queue expected.** The queue wanted Schieberle, Fischer & Hofmann 2003 (10th Weurman Symposium proceedings, pp. 447-452, the first CAMOLA chapter); this file is Schieberle's 2005 short review of CAMOLA, which cites that chapter as its ref. 13. **No cysteine, no MFT, no thiol anywhere.** Its CAMOLA data are on 4-hydroxy-2,5-dimethyl-3(2H)-furanone (HDMF, "4-HDF") from 1 : 1 [12C6] / [13C6]-glucose + proline (intact C6 under roasting; **~60 % from C3 + C3 recombination in aqueous buffer, 145 °C, 20 min**) and on two bitter pyrrolizines; plus a fed hydroxyacetone + 2-oxopropanal → HDMF pH series (pH 3 / 5 / 7: **0.05 / 0.28 / 1.10 mol %**). **The 2003 Weurman chapter remains not obtained.**
DOI: 10.1196/annals.1333.029 (printed on p. 236)

**Read 2026-11-05 (first reader) from `data/articles/schieberle2005b.pdf`** (13 pp., publisher PDF with a clean text
layer; 7.6 MB because the figures are images). Prose and Tables 1-3 from `pdftotext -layout` (the tables are short and
came through intact). Figs. 1-10 are a flow chart, structures, schemes and mass spectra (Figs. 5 and 8); **no
isotopomer intensities were digitised from the spectra**: the shares below are the ones printed in the text.

**Reservation.** This paper is not one of the case's candidates (the candidate was the 2003 chapter, role "evidence for
D3 only (DFA); never decisive for D1 / D2"). It has no Cys + glucose pot and no thiol, so no value here falls in a pot
the case names as decisive; nothing is reserved. **If** Table 3 is ever used as a band source for rival R3's fragment
recombination, that needs an amendment registered before any M4 engine result (case sec. 2.9, last paragraph). No
engine run.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[Q]** qualitative statement
by the author · **[C]** cited (another paper's data) · **[P]** proposed scheme · **[Z]** derived by me · **[!]** flag.
Provenance: **first reader only, second reader pending.** Nothing here is in a ledger, fit table or benchmark.

## 0. Identity

| field | value as printed |
|---|---|
| Author | **Peter Schieberle** (sole author), Technical University of Munich, D-85748 Garching; Chair for Food Chemistry, Lichtenbergstrasse 4 |
| Venue | *Ann. N.Y. Acad. Sci.* **1043**: **236-248** (2005), © New York Academy of Sciences |
| DOI | **10.1196/annals.1333.029** (printed p. 236) |
| Kind | "short state-of-the-art review" (abstract) with some unpublished data ("Details ... will be published elsewhere") |

**File-name note.** `data/articles/schieberle2005.pdf` is a different paper (Schieberle & Hofmann, Woodhead 2005,
`schieberle2005_extraction.md`). This file's stem `schieberle2005b` is kept as its key (author + year + suffix); no alias
needed. **The expected paper**, Schieberle, Fischer & Hofmann 2003, "The carbohydrate modul labelling technique ..."
(Weurman 2002 proceedings, Beaune, pp. 447-452; no DOI known), **is not on disk**: it stays on the queue.

**Laboratory.** DFA / TU München, Schieberle: `dfa_hofmann_schieberle` in `data/keys/laboratories.yml`.

## 1. What this paper answers for the model

| model need | what Schieberle 2005b gives | tag |
|---|---|---|
| Glucose + Cys CAMOLA on MFT (D3 evidence) | **none**: no cysteine, no MFT | — |
| Glucose fragment recombination in water at 145 °C (R3's premise, in a proline pot) | HDMF from 1 : 1 [12C6] / [13C6]-glucose + Pro: low water → 1 : 1 [12C6] : [13C6] (intact); aqueous buffer, 145 °C, 20 min → [13C3] the main isotopomer, "about 60%" of HDMF from two C3 fragments | [M] / [Q] |
| Fed C3 + C3 → C6 furanone vs pH | hydroxyacetone + 2-oxopropanal → HDMF: 59.6 / 364.5 / 1450.5 µg = 0.05 / 0.28 / 1.10 mol % at pH 3 / 5 / 7 (Table 3, "Data from Frauendorfer & Schieberle", then submitted) | [M] / [C] |
| Amine (Pro) pyrrolizines from fragments | 5-AMDP 34× 7-AMDP; isotopomers 1 : 2 : 1 (C3 + C3) and 1 : 1 : 1 : 1 (C2 + C4) | [M] |
| Rate constants, NF, MFT, FFT | **none** | — |

## 2. Conditions [M]

| experiment | conditions as printed | missing |
|---|---|---|
| HDMF CAMOLA, roasting (Fig. 5A) | glucose + [13C6]-glucose (1 + 1) with L-proline, **15 min at 160 °C, 10 % water** | amounts, vessel |
| HDMF CAMOLA, aqueous (Fig. 5B) | the same "in aqueous buffer in a autoclave at 145°C" [V], **20 min** | buffer, pH, amounts [!] |
| Table 3 | 2-oxopropanal + hydroxyacetone in aqueous buffer, pH 3.0 / 5.0 / 7.0 | amounts, T, time, vessel not printed; the mol % is consistent with a 1 mmol basis (sec. 4.3) [Z] |
| 5- / 7-AMDP CAMOLA | 1 : 1 [12C6] / [13C6]-glucose + proline in phosphate-buffered saline, **145 °C, 1 h, autoclave** | amounts, pH |
| C3 / C4 + C2 feeding (AMDP) | C3 compounds + Pro, **150 °C, 1 h, autoclave**; erythrose + acetaldehyde + Pro | amounts; "data not shown" |
| Fig. 7 bitter compounds | glucose + proline, **150 °C, 1 h**, aqueous | — |

## 3. The assay

GC/MS isotopomer analysis (m/z 128 / 131 / 134 for HDMF; m/z 163-169 for AMDP); quantification of AMDP and HDMF "by
stable isotope dilution assays" is the group's stated practice (Fig. 1 concept) but **no assay detail, LOD, recovery or
replicate noise is printed for any number in this paper** [!]. Table 1 / 2 food values are from earlier work [C].

## 4. Data

### 4.1 HDMF CAMOLA, glucose + proline [M] / [Q] (pp. 241-242)

| condition | isotopomer pattern | author's reading |
|---|---|---|
| roasting, 160 °C, 15 min, 10 % water | **1 : 1 [12C6] : [13C6]**; [13C3] "only as a very minor compound" | intact glucose skeleton |
| aqueous buffer, 145 °C, 20 min | **[13C3] the main isotopomer** | statistical weighting: m/z 131 should be 2× m/z 128 or 134 for C3 + C3, so **"about 60%"** of HDMF by recombination of two C3 fragments |

### 4.2 Table 3 (p. 243): HDMF from 2-oxopropanal and hydroxyacetone vs pH [M] ("Data from Frauendorfer & Schieberle")

| pH | HDMF µg | mol % |
|---|---:|---:|
| 3.0 | 59.6 | 0.05 |
| 5.0 | 364.5 | 0.28 |
| 7.0 | 1450.5 | 1.10 |

Title as printed: "... formed from 2-oxopropanol [sic] and hydroxyacetone"; the text says 2-oxopropanal. Author's
reading [Q]: formation "drastically lowered" at pH 3 vs 7, so a retro-aldol step (Fig. 6) is "most probable" [P].

### 4.3 Arithmetic [Z]

1. Table 3 on HDMF (128.13 g/mol) and a 1 mmol basis: 0.047 / 0.284 / 1.13 mol % vs printed 0.05 / 0.28 / 1.10 (the
   printed pH 7 % implies 1409 µg, 3 % below 1450.5). pH 7 / pH 3 = **24×** (+1.39 dex), pH 5 / pH 3 = 6.1×.
2. **[!] Two of the three µg values coincide exactly with cells of Hofmann & Schieberle 1998** (`hofmann1998_reconciliation.md`):
   **364.5 µg** = Table 6, furan-2-aldehyde at pH 5 (from hydroxyacetaldehyde + 2-oxopropanal), and **59.6 µg** = Table 9,
   1-mercapto-2-propanone at pH 5 (glucose + Cys); and 1450.5 is within 0.5 % of Table 6's FA at pH 7 (1443.1). The
   mol % column here is consistent with HDMF on a 1 mmol basis (item 1), which argues the numbers are HDMF's; but two
   exact coincidences in one three-row table are improbable by chance. A transcription from the 1998 tables cannot be
   excluded. Do not use Table 3 without Frauendorfer & Schieberle's primary paper.
3. The "about 60 %" is the author's statistical estimate; the spectra (Fig. 5) were not digitised here.

### 4.4 Pyrrolizines, glucose + proline, PBS, 145 °C, 1 h [M] (pp. 243-244)

5-AMDP : 7-AMDP = **34 : 1** by amount. 5-AMDP isotopomers m/z 163 : 166 : 169 = **1 : 2 : 1** (C3 + C3); 7-AMDP m/z
163 : 165 : 167 : 169 = **1 : 1 : 1 : 1** (C2 + C4). Glyceraldehyde + hydroxyacetone (1 : 1) + Pro make ~30× more 5-AMDP
than glucose; erythrose + acetaldehyde ~20× more 7-AMDP (data not shown) [M].

### 4.5 Tables 1-2 (food HDMF, context) [C]

Cocoa, unroasted / roasted: HDMF 17 / 990 µg/kg; 3-hydroxy-4,5-dimethyl-2(5H)-furanone 11 / 10; 3-methylbutanal 1100 /
26,400; 3-methylbutanoic acid 14,200 / 17,300. HDMF in foods 820 (barley malt) to 131,000 µg/kg (coffee), OAV 14-2183
(threshold in cellulose). No engine observable.

## 5. What it can and cannot test in the model

**Cannot:** anything on MFT, FFT, NF, Cys or the which-carbon question; any rate constant (one time each).

**Can (evidence, weak):**
- That glucose fragments recombine into C6 furanones in water at 145 °C in an amine pot (60 % of HDMF with Pro), but
  stay intact under low-water roasting. This bears on R3's premise (fragment recombination at Hofmann's 145 °C) only by
  analogy: different amine (Pro), different product (HDMF, a C6), no thiol.
- A fed C3 + C3 → HDMF pH law rising ~24× from pH 3 to 7, the same direction as Hofmann 1998 Table 6's C2 + C3 → NF
  (0.13 / 0.78 / 3.18 mol %), **subject to flag 4.3 item 2**.
- The engine has no HDMF-from-fragments step; HDMF is not an M4 species.

## 6. Conflicts and agreements with existing records

1. **[!] Possible transcription of Hofmann & Schieberle 1998 cells into Table 3** (sec. 4.3 item 2).
2. **Agreement:** Hofmann 1998 Table 6 (fragment → NF, FA rising with pH) and this Table 3 (fragment → HDMF rising with
   pH) give the same direction; same laboratory, so not independent.
3. **Nuance vs `k5b_dmhf_synthesis.md` row B11**, which credits this paper with a route starting "2 MG → (Cannizzaro)
   → MG + acetol". This paper's Fig. 6 starts from hydroxyacetone (acetol) + 2-oxopropanal and names no Cannizzaro step;
   Cannizzaro is Cannon & Ho's reading (`cannon2023_extraction.md`), not Schieberle's [!].
4. **`docs/validation/isotope_topology_evidence.md`** lists this paper from its abstract only ("no compound-specific
   data in the abstract"). It does hold compound-specific shares (sec. 4.1, 4.4); none concerns a sulfur compound.
5. Compound keys not checked (no data rows written).

## 7. Provenance tiers (first reader only)

Sec. 4.1 / 4.4 shares: `direct_measurement` reported in a review without methods, confidence low-medium. Table 3:
`secondary` (cited as another paper's data, then unpublished), confidence low (flag 4.3 item 2). **Proposed ledger
entries: none.**

## 8. M4 role (case 8660c005)

**Role fixed before reading** (case sec. 2.9, candidate table) for the expected paper, Schieberle, Fischer & Hofmann
2003: "evidence for D3 only (DFA); never decisive for D1 / D2"; the measurement that would count: "a glucose + Cys
CAMOLA on MFT".

| element of the criterion | what this file has | verdict |
|---|---|---|
| it is the 2003 Weurman chapter | **no**: Schieberle 2005, Ann. N.Y. Acad. Sci. 1043:236 (cites the 2003 chapter as ref. 13) | not the candidate |
| a glucose + Cys CAMOLA on MFT | **no**: CAMOLA on HDMF and two pyrrolizines, glucose + proline only; no cysteine, no MFT | **not met** |

**Verdict: NOT MET** (and not the candidate paper). It bears on M4 only weakly (sec. 5). **The 2003 chapter remains
not obtained**; its role in the case is unchanged.
