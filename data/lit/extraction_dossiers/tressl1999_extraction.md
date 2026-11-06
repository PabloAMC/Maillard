# Tressl & Rewicki 1999 (Flavor Chemistry: 30 Years of Progress, pp. 305-325) — EXTRACTION (file stem `tressl1999`)
### *"Heat Generated Flavors and Precursors"* (chapter 26; a review of the authors' [13C]-sugar labelling work)
### A review chapter, no Methods section and no new pots. One sentence on MFT: its formation pathway from D-glucose "was established" in experiments "analogous" to the [1-13C]-glucose / Cys ones, **with no data, no figure, no conditions and no citation.** **No MFT isotopologue is printed anywhere.** The one positional result on a glucose + Cys thiol is **Fig. 14: 2-furfurylthiol (FFT) and furfural from [1-13C]-D-glucose / L-cysteine are each ~1 : 1 unlabelled : singly labelled** (C-1 / C-6 label equilibration, then C1-C5 cleavage via 3-deoxyglucosone), conditions not stated.
DOI: 10.1007/978-1-4615-4693-1_26 (not printed in the PDF; Crossref 2026-11-05: "Heat Generated Flavors and Precursors", Roland Tressl and Dieter Rewicki, in Flavor Chemistry, Springer US, 1999, pp. 305-325)

**Read 2026-11-05 (first reader) from `data/articles/tressl1999.pdf`** (21 pp., a scan with an Acrobat Paper Capture
OCR layer that garbles every scheme). Prose from `pdftotext -layout`; the text was grepped for every MFT spelling
("furanthiol", "methyl-3-fur", "mercapto-2-methylfuran", "thiol"): one hit for MFT (p. 318). **Figs. 12, 13 and 14 were
read from 100 dpi renders of pp. 317-318** (scratchpad `m4pdf/tressl1999/p-13.png`, `p-14.png`). The other figures
(Figs. 1-11, 15-17: schemes) were checked through their captions and the prose; none names MFT. No figure carries a
digitisable quantity; all numbers below are printed percentages in the schemes.

**Reservation (M4, case `docs/reference/cases/M4_glucose_mft_protein_amine.md`, branch `m4-case-refresh`, commit
8660c005, sec. 2.9).** The case fixed this paper's role before it was read: decisive for the which-carbon gap only if it
reports MFT from [1-13C]- or [6-13C]-glucose. The glucose + Cys label values below (Fig. 14) carry
**RESERVED: unseen by the engine (M4); no engine run on these pots until M4's design is frozen** (**[RESERVED M4]**).
The role is not met (sec. 7); release is the coordinator's call. No engine run.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this chapter · **[Q]** qualitative statement
by the authors · **[C]** cited by the authors · **[P]** proposed scheme, not shown by data here · **[Z]** derived by me
· **[!]** flag. Provenance: **first reader only, second reader pending.** Nothing here is in a ledger, fit table or
benchmark.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Roland Tressl**¹ and **Dieter Rewicki**² (byline p. 305) |
| Affiliations | ¹Technische Universität Berlin, Institut für Biotechnologie, Seestraße 13, 13353 Berlin; ²Freie Universität Berlin, Institut für Organische Chemie, Takustraße 3, 14195 Berlin |
| Venue | *Flavor Chemistry: 30 Years of Progress*, edited by Teranishi et al., Kluwer Academic / Plenum Publishers, New York, **1999**, pp. **305-325**, chapter **26** (the number printed at the head of p. 305) |
| DOI | **10.1007/978-1-4615-4693-1_26** (not printed; Crossref matches title, both authors, pages 305-325, 1999) |
| Kind | review; the newest work cited is "1998c, submitted" |

**Laboratory.** TU Berlin, Tressl group (with Rewicki, FU Berlin): the laboratory of `tressl1994_extraction.md`. Not in
`data/keys/laboratories.yml`. **Not** the TU Berlin Kanzler / Rohn or Kroh groups (different institute and decade).

## 1. What this chapter answers for the model

| model need | what Tressl 1999 gives | tag |
|---|---|---|
| Which glucose carbon MFT loses (positional label on MFT) | **none**: one sentence asserting the pathway "was established", no isotopologues, no conditions, no reference | [Q] |
| Positional label on another glucose + Cys thiol | **FFT and furfural, [1-13C]-glucose / Cys: each 50 % singly labelled** (Fig. 14) | [M] |
| Route proposed for that FFT | 3-deoxyglucosone (3-deoxyaldoketose) → C-1 / C-6 equilibration → C1-C5 (vinylogous retro-aldol) cleavage → C5 + H2S (Strecker) | [P] |
| Cys fate in [1-13C]-glucose / Cys | Strecker (mercaptoacetaldehyde → acetaldehyde + H2S), transamination (→ 3-mercaptopropionic acid), β-elimination (pyruvate + H2S → 2-mercaptopropionic acid) (Fig. 13) | [P] / [C] |
| MFT sources | thiamine (Werkhoff 1990 [C]); only a minor part of MFT from Cys in [34S]-Cys + thiamine + sugar pots (Bolton 1994 [C]); NF + H2S gives MFT (van den Ouweland & Peer 1975 [C]) | [C] |
| Rate constants, yields, pH law | **none** | — |

## 2. Conditions

**None printed for any experiment.** The chapter states no concentrations, buffer, pH, temperature, time, vessel or
headspace for the experiments behind Figs. 12-14. The glucose + Cys labelling work it cites is Tressl, Kersten, Nittka
& Rewicki 1994 (ACS Symp. Ser. 564, pp. 224-235): equimolar, aqueous, initial pH 5, 160 °C, 1.5 h according to
`tressl1994_extraction.md`. **Fig. 14 has no citation**, so its conditions cannot be assumed to be those [!].

## 3. The assay

Not described. Isotopomer shares were obtained by GC-MS (and 1H-NMR for one pyrazine, Weenen & Tjan 1992 [C]) in the
cited works. No LOD, no replicate noise anywhere in the chapter.

## 4. Data (all printed in schemes; no error stated)

### 4.1 Fig. 14 (p. 318): [1-13C]-D-glucose / L-cysteine [M] [RESERVED M4]

| product | label pattern as printed |
|---|---|
| 2-furylmethanethiol (FFT) | **50 % singly labelled** (caption: "1 : 1, unlabeled to singly labeled") |
| 2-formylfuran (furfural) | **50 % singly labelled** |

Authors' reading [Q]: the same precursors make furfural and FFT; the 1 : 1 ratio "indicate[s] an equilibration of the
label between the 1- and 6-position" with subsequent vinylogous retro-aldol cleavage. The scheme draws 3-deoxyglucosone →
an equilibrating hexosulose pair → loss of C-1 or C-6 as formaldehyde [P].

### 4.2 MFT (p. 318) [Q]

The only MFT sentence on glucose: in "analogous experiments" the formation pathway of MFT from D-glucose "was
established" [V, two fragments]. **No isotopologue share, no figure, no conditions, no reference.** The next sentences
are about thiamine and Bolton et al. 1994 ([34S]-Cys; "only a minor part" of MFT from Cys) and NF + H2S (van den
Ouweland & Peer 1975).

### 4.3 Other label results printed in the schemes (not M4 species; context) [M]

- Fig. 12 ([1-13C]-glucose / amino acid, via 3-deoxyglucosone and retro-aldol): aminocarbonyl pools A* 76 % / A 24 %,
  B* 41 % / B 59 %; pyrazine and pyrrolizine isotopomers unlabelled / singly / doubly 6/36/58 (×2), 12/51/37,
  16/54/30, 15/55/30, 36/47/17, 38/47/15, 34/49/17.
- Weenen & Tjan 1992 [C]: [1-13C]-glucose / Asn 2,5-dimethylpyrazine, 26 % 13C in the methyls, 20 % in ring CH (1H-NMR).
- Simon 1962 / Simon & Heubach 1965 [C]: [1 (or 6)-14C]-glucose → C6-reductone, label 72 % ring / 28 % exocyclic methyl.
- Fig. 9 ([1-13C]-glucose / Pro, GABA, Ile): cyclopentenones, pyranones, furanones "100% singly labeled", i.e. intact
  C6 via 1-deoxyglucosone (or diacetylformoin).
- [1-13C]-lactose / GABA: maltol 97 % of the aroma, 100 % singly labelled [M, p. 314].
- p. 308: from ketoses, "up to 50%" of Amadori (not only Heyns) rearrangement products, depending on the amino acid
  [C, Heyns et al. 1961]: glucose-fructose interconversion.

## 5. What it can and cannot test in the model

**Cannot:** the which-carbon question for MFT (no data); any rate, yield, Ea or pH law; anything needing conditions.

**Can (evidence only, after release):** a glucose + Cys **FFT** positional share (50 % from [1-13C]-glucose), which bears
on part E / F25-c (hexose FFT skeleton), not on part G. It says the FFT C5 can carry glucose C-1 half the time, which a
simple C-1-loss route (Tressl 1994's furfuryl alcohol) does not give. Conditions unknown, so it cannot be put against
an engine pot.

## 6. Conflicts and agreements with existing records

1. **[!] CONFLICT with `tressl1994_extraction.md` (same laboratory).** Tressl 1994 Table IV: furfuryl alcohol (and
   2-hydroxymethylthiophene) from [1-13C]-glucose + Cys is **unlabelled**, and 100 % singly labelled from
   [6-13C]-glucose: C-1 lost, no equilibration. Tressl 1999 Fig. 14: furfural and FFT from [1-13C]-glucose + Cys are
   **50 % labelled**: C-1 and C-6 equilibrated. Same carbon skeleton (furfuryl), opposite label result. Either the
   conditions differ (Fig. 14 cites no source) or one of the two is a scheme, not a measurement. Tressl 1994 reports no
   FFT from glucose at all. Unresolved; do not use either share for the furfuryl skeleton without the primary data.
2. **The M4 case sec. 2.1 summary of Tressl 1994** ("every glucose-derived C5 sulfur heterocycle ... is glucose C2-C6")
   is contradicted for FFT by this chapter's Fig. 14 (flag 1).
3. **Agreement with Bolton 1993 / 1994** (`bolton1993_extraction.md`): Cys a minor MFT sulfur source when thiamine is
   present [C].
4. **The case's which-carbon gap stands** (case sec. 2.1 last row, sec. 4.3): no positional label on MFT exists on disk;
   this chapter's "was established" points to no datum.

## 7. Provenance tiers (first reader only)

Fig. 14 shares: `review_scheme_value` (printed in a review figure with no source cited and no conditions), confidence
low. Fig. 12 / 9 values: `review_secondary` (summaries of Tressl 1993 / 1995 and Weenen & Tjan 1992), confidence
medium for the label shares, none for conditions. **Proposed ledger entries: none.**

## 8. M4 role (case 8660c005)

**Role fixed before reading** (case sec. 2.9, candidate table): "decisive for the which-carbon gap only if it reports
MFT from [1-13C]- or [6-13C]-glucose". The measurement that would count: "MFT isotopologues from a positional glucose
label".

| element of the criterion | what the chapter has | verdict |
|---|---|---|
| a positional glucose label ([1-13C] or [6-13C]) with Cys | yes for FFT and furfural ([1-13C], Fig. 14, 50 % singly labelled each) | met for FFT only |
| MFT isotopologues reported from it | **no**: one sentence ("the formation pathway ... was established"), no share, no figure, no citation | **not met** |

**Verdict: NOT MET.** No MFT isotopologue is printed. The which-carbon gap remains open; the designed experiment (case
sec. 4.4, arm 1) is still the only way to close it.
