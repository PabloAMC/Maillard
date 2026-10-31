# Schubert 1936 — EXTRACTION (thioglycolic acid anilide + formaldehyde / butyraldehyde / chloral / pyruvic acid / quinone / isatin give 1-hydroxyalkyl thioethers (hemithioacetals); L-cysteine + formaldehyde / butyraldehyde / chloral / benzaldehyde / furfural give crystalline condensation products (thiazolidines), cysteine + pyruvic acid an addition compound; preparative chemistry with qualitative dissociation tests; NO rate, K, pH series or temperature series)

### THE COMPANION schubert1939_extraction.md HOPED FOR, AND IT HAS NO CONSTANTS EITHER. What it adds is a trap experiment (iodoacetate) showing a cysteine-benzaldehyde thiazolidine releases free cysteine in bicarbonate at room temperature within ~0.5 h to overnight, and the statement that these compounds, formaldehyde's included, absorb iodine "as readily as cysteine". Formaldehyde's thiazolidine releases cysteine with k_r ~10⁻⁶ s⁻¹ (Kallen), so the iodine titration measures ring-chain OPENING, not release.

**Read 2026-11-01 from `data/articles/schubert1936.pdf`** (10 pp., J. Biol. Chem. 114, 341-350; publisher re-issue
with CC BY footer; iText text layer of a scan). Prose from `pdftotext -layout`, which is clean; p. 344 checked against
an 80 dpi render (scratchpad `dos6/sch-04.png`). No figures or tables.

**Legend.** **[V]** verbatim · **[M]** number printed · **[F]** figure · **[Q]** qualitative claim · **[C]** cited ·
**[P]** proposed · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| model need | what Schubert 1936 gives | tag |
|---|---|---|
| forward / reverse rate constants | **none.** Preparations: "a crystalline product begins to separate in 10 to 30 minutes" (cysteine-HCl 10 g + potassium acetate 6 g in water or 50 % ethanol + aldehyde, room temperature) | [M] (protocol time) |
| K | **none** | — |
| pH law | qualitative, nitroprusside test for free -SH: negative in bicarbonate for all cysteine compounds except chloral's; positive in weak ammonia for benzaldehyde, furfural, chloral; "very faint" for butyraldehyde; none for formaldehyde and phenylglyoxal. "these compounds dissociate to different extents, depending on the nature of the aldehyde and on the alkalinity" | [Q] |
| buffer catalysis, Ea | none | — |
| trap of released cysteine | benzaldehyde-cysteine compound 1 g in 12 mL water + 1 g NaHCO₃ + 1 g iodoacetic acid: "After standing a half hour benzaldehyde begins to separate"; overnight, no thiazolidine left; S-carboxymethylcysteine isolated (yield not printed) | [M] / [Q] |
| trap of the thiol (iodine) | cysteine-aldehyde compounds "absorb iodine just as readily as cysteine", one I atom per molecule, cystine formed; acetylated compounds absorb iodine "only very slowly" | [Q] |
| second trap (phenylhydrazine) | benzaldehyde phenylhydrazone forms "readily" from the benzaldehyde compound in 2 M potassium acetate | [Q] |
| thiol + pyruvic acid / quinone | hemithioacetal-type addition compounds (cysteine-pyruvic acid compound gives no nitroprusside test even in weak ammonia); quinone + thioglycolic anilide gives a red adduct "immediately" in alcohol | [Q] |

**Answers.**

1. **`k_ttca_cys`, E1 / E4: no bound** (confidence ~99 %; all ten pages read). No sugar, no temperature, no rate.
2. **Release to an irreversible trap (benzaldehyde-TCA, bicarbonate pH ~8, room temperature)**: free aldehyde appears
   within ~0.5 h and conversion is complete overnight [M/Q]. As a release rate this says only k_r(25 °C, pH ~8, aryl
   C-2) ≳ 10⁻⁵ s⁻¹ (completion in ≤ ~16 h) [Z]. Aryl thiazolidines are the LEAST stable class (Pesek 1975), so this
   does not transfer to TTCA in either direction.
3. **Iodine titration measures ring-chain opening, not release** [Z]. Formaldehyde's TC is in the list that absorbs
   iodine "as readily as cysteine", yet its release constant at 25 °C is ~1 × 10⁻⁶ s⁻¹ (t½ ~9 days,
   `kallen1971b_extraction.md`). The fast iodine uptake therefore reports C-S opening to the iminium-thiol, which iodine
   drains; it gives no rate for C-N cleavage to free cysteine + aldehyde. `schubert1939_extraction.md` says the same for
   the sugar compounds ("ring-chain opening ... not slow").
4. **"Buffering of thiol concentration"** [Q, the author's phrase]: the reservoir concept the engine's TTCA ⇌ Cys + Xyl
   implements, stated in 1936; no numbers.
5. **MFT loss law**: nothing quantitative. Thiol + quinone adds "immediately" (red product) [Q]; consistent in kind with
   Nikolantonaki 2012's fast thiol-quinone constants.

## 2. Identity

| field | value |
|---|---|
| Title | "Compounds of Thiol Acids with Aldehydes" |
| Author | Maxwell P. Schubert, Laboratories of The Rockefeller Institute for Medical Research, New York |
| Venue | *J. Biol. Chem.* **114**, 341-350 (1936); received February 21, 1936 |
| DOI | not printed. Crossref (2026-11-01): **10.1016/S0021-9258(20)65226-7**, same title, Schubert, JBC 114, 341-350, 1936. Matches the brief |
| Related | Schubert 1935 JBC 111, 671 (thiols + methylglyoxal / phenylglyoxal; Crossref 10.1016/S0021-9258(18)75013-8; not on disk); Schubert 1939 (on disk) |

## 3. Conditions

- Room temperature throughout; aldehyde in excess; solvents water, 50 % ethanol or ethanol with potassium acetate;
  products on ice for hours, recrystallised [M]. Concentrations: cysteine-HCl 10 g in 25-200 mL (≈0.3-2.5 M as
  charged, [Z]); no pH measured, no buffer beyond acetate / bicarbonate / ammonia.
- Analysis: elemental S, N, Cl; melting points; iodine titration in alcohol + 2 mL 4 M potassium acetate ("end-points
  are quite sharp", "titrated to a sharp end-point in a few minutes" for the anilide compounds) [M/Q].

## 4. Quantitative results

Only composition and iodine values [M]: thioglycolic anilide compounds, iodine value (equivalents I per mole) 0.95
(formaldehyde), 1.07 (butyraldehyde), 1.00 (chloral), 1.3 (pyruvic acid); thioglycolic anilide itself 0.97. Melting
points and analyses for all products (pp. 345-350). No kinetic number anywhere.

## 5. Mechanism evidence and its kind

- Rockefeller Institute (Schubert; Michaelis lineage). Kind 3 (trapping: iodoacetate, iodine, phenylhydrazine) and
  derivatisation (acetylation of the OH / NH stabilises the compounds against iodine and nitroprusside) [Q].
- Proposed structures [P]: thiazolidine ring (V) for cysteine + aldehydes (independently found by H. T. Clarke, footnote);
  1-hydroxyalkyl thioether (I) for the anilide adducts, preferred over a "loose molecular compound" (II) because
  acetylation gives (III) [Q].

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| CF-076 (dilute vs concentrated exchange) | LEAVES UNCHANGED | no rate, no concentration series |
| cys_exchange prereg E1 k0 fast-edge basis ("Schubert 1939 ring opening as fast as free Cys to I₂ at RT") | **CONTRADICTS its use as a RELEASE rate** | formaldehyde's TC titrates as fast, yet releases with k_r ~10⁻⁶ s⁻¹; iodine reads C-S opening |
| evidence map 2.12 cysteine fate | LEAVES UNCHANGED | |
| MFT loss law | LEAVES UNCHANGED | qualitative thiol-quinone addition only |

## 7. Disagreements found

1. **With `docs/history/cys_exchange_2026-10-31/prereg.yaml` (E1 `k0_100` basis)**: the iodine evidence was cited for
   the fast edge of a release constant; on Kallen's numbers for the formaldehyde compound it bounds only ring-chain
   opening (sec. 1, answer 3). `schubert1939_extraction.md` drew the correct distinction; the prereg did not carry it.
   Cause: interpretation, not data.
2. None with other dossiers on numbers (there are none).

## 8. Ledger-ready measurements

None. If a qualitative row is ever wanted: the iodoacetate trap (benzaldehyde-TCA → benzaldehyde + S-carboxymethylcysteine,
bicarbonate, room temperature, onset ~0.5 h, complete overnight), as the one timed release observation.
