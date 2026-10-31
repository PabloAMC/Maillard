# Schubert 1939 — EXTRACTION (L-cysteine + six reducing sugars incl. D-xylose and D-arabinose condense 1:1 with loss of water to thiazolidines at room temperature; preparative note; qualitative thiol tests; NO rate constant, NO equilibrium / dissociation constant, NO pH, temperature or buffer series)

### A THREE-PAGE PREPARATIVE NOTE, NOT A PHYSICAL-CHEMISTRY PAPER. The brief expected dissociation constants; there are none. The constants, if any, are in the companion Schubert 1936 (J. Biol. Chem. 114:341, thiol acids + aldehydes), cited as footnote 1 and not on disk. What this note does add is the first (and only, in the corpus) qualitative evidence that a **xylose-** and **arabinose-cysteine** thiazolidine releases its full thiol content **as fast as free cysteine** to an irreversible trap (iodine) at room temperature.

**Source on disk:** `data/articles/schubert1939.pdf` (3 pp., J. Biol. Chem. pp. 601-603; publisher re-issue with a
CC BY open-access footer). **Typeset text layer (iText re-issue of a scan).** Prose read from the `pdftotext -layout`
layer and **checked against 110 dpi page renders** (scratchpad `pdf2/sch-1.png` ... `sch-3.png`); Table I was read
from the render because the text layer garbles its headers and molecular formulae.

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[Q]** qualitative statement by the
author · **[C]** cited from another paper · **[P]** proposed, not measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question put to it (plan step B7b, `results/validation/kinetic_core_dpo_partition_prereg.md` secs. 7 and 12,
blocking item B2): the engine's TTCA ring opening `k_ttca_cys` (~0.028 min⁻¹, fitted, no pH factor; reverse step
`r_ttca_form` by detailed balance with Cao's K ~60 L/mol) cannot reconcile Cao 2017 (`cao2017_extraction.md`: 0.4 M,
0.2 M phosphate pH 5.5, apparent K 58-65 L/mol, equilibrium within 30 min at 90 °C) with Zhai 2021
(`zhai2021_extraction.md`: 10 mM TTCA in unbuffered water, natural pH 5.5, TTCA loss ≤ 0.0033 min⁻¹ at 100 °C).

| what the model needs | what Schubert 1939 gives | tag |
|---|---|---|
| rate of ring opening / formation, any aldose-Cys thiazolidine | **none timed.** Formation: preparations stand 48 h at room temperature, then 50-70 h more after pyridine; these are isolation protocols, not kinetics | [M] (protocol) |
| equilibrium constant K = [TCA]/([Cys][aldose]) | **none.** No free-Cys fraction, no titration curve | — |
| pH dependence | **qualitative, one sign:** "in a solution containing excess sodium bicarbonate none of these compounds gives a positive thiol test with sodium nitroprusside", while free cysteine "gives a strong test"; in dilute NH₄OH / (NH₄)₂SO₄ the glucose and lactose compounds give "fair" tests and the rest (incl. xylose, arabinose) "very faint" | [Q] |
| temperature dependence | **none** (room temperature only; recrystallisation from hot water survives for Man / Gal) | — |
| buffer / general-acid catalysis | **not tested** | — |
| release of the thiol to an irreversible trap | "Aqueous solutions of all of these compounds absorb iodine as rapidly as does free cysteine and in amount equivalent to the cysteine they contain"; in bicarbonate they "react rapidly with sodium iodoacetate to give S-carboxymethylcysteine" (galactose compound, ~5 h at room temperature, then worked up; yield not printed) | [Q] / [M] |

**Answer to B2: nothing quantitative** (confidence ~99 % that the note holds no rate, K, pH series or buffer
comparison: all three pages read). What it does say, weighed:

1. **The ring-chain opening of a pentose-cysteine thiazolidine is not slow at room temperature** [Z, from Q]. The
   open-chain (iminium / Schiff base) form carries the free thiol; iodine oxidises only that form, so complete,
   cysteine-equivalent iodine uptake "as rapidly as" free cysteine means the closed ring re-supplies the open form on
   the titration timescale (seconds to minutes, not timed) in the compounds' own solutions ("acid to litmus", pH not
   stated; ring nitrogen pKa ~6 by analogy with Kallen 1971's 6.24, so below that). At 90-100 °C the step can only be
   faster. This is the same direction as Roberts 1987's C-2 mutarotation (`roberts1987_extraction.md`) and Nagasawa
   1981's epimerisation (`nagasawa1981_extraction.md`), now for **xylose and arabinose themselves**.
2. **At mildly alkaline pH the closed form dominates strongly** [Q]: no detectable thiol in excess bicarbonate
   (pH ~8 [Z]) at whatever concentration the spot test used (not stated). Same sign as Kallen 1971 (K_app rises ~20×
   from pH 5.6 to 7.5) and Roberts' pyridine remark. Hexose (glucose) and lactose adducts dissociate more than the
   pentose ones in ammonia [Q]: an ordering, not a number.
3. **Caveat that limits (1).** The trap captures the open thiol, so (1) bounds the **C-S ring opening**, not the
   hydrolysis of the open iminium to free cysteine + xylose; the engine's `r_ttca_cys` lumps both. A fast opening
   with a slow (possibly buffer-catalysed) hydrolysis remains possible.

**Bearing on the Cao-vs-Zhai contradiction (my reading, confidence ~0.6):** the note weakens "the uncatalysed ring
opening is intrinsically slow in water" as the explanation, and so shifts weight toward (a) the **equilibrium
position** differing between the two pots (K at 10 mM, pH 5.5, 100 °C in Zhai would have to be ≳ 10⁴ L/mol for
TTCA to stay ≥ 90 % closed, ≳ 2 × 10⁵ for the 98 % Zhai reports at 20 min, if equilibrium were reached [Z: K =
(1-α)/(cα²), c = 0.01 M, α = 0.1 / 0.02; at K = 60, α ≈ 0.77], against 60 L/mol in Cao), (b) a slow
**hydrolysis** step after the opening, or (c) the analytics re-forming TTCA on cooling (Zhai). It cannot rank them,
and gives no number for any.

**Can anything anchor the engine's constant?** No. No rate, no K, no conditions beyond "room temperature". **Keying:
nothing worth keying** (Table I is elemental analysis, melting points and isolated yields: 20 % for xylose, 65 % for
arabinose; isolated yields after alcohol precipitation are not equilibrium yields).

**Leads (none on disk; DOIs checked against Crossref 2026-10-31):**
- **Schubert 1936**, *J. Biol. Chem.* 114:341-350, "Compounds of thiol acids with aldehydes" — DOI
  10.1016/s0021-9258(20)65226-7. Footnote 1 of this note; the likeliest place for the dissociation constants the
  brief expected (aldehyde-cysteine thiazolidines; Roberts 1987 cites it for the pH shift of the equilibrium).
- **Pesek & Frost 1975**, *Tetrahedron* 31:907-913, "Decomposition of thiazolidines in acidic and basic solution" —
  DOI 10.1016/0040-4020(75)80099-8 (cited by Nagasawa 1981, ref 10a). **The most direct lead for a ring-opening /
  hydrolysis rate and its pH profile.**
- Kallen 1971b, *J. Am. Chem. Soc.* 93:6236-6248 — DOI 10.1021/ja00752a040 (already listed in
  `docs/guides/LITERATURE_SEARCH_2026-10-29_phosphate_thiazolidine.md`; rates and buffer catalysis for the
  formaldehyde analogue).
- Ratner & Clarke 1937, *J. Am. Chem. Soc.* 59:200-206, "The action of formaldehyde upon cysteine" — DOI
  10.1021/ja01280a050.

## 0. Identity

| field | value |
|---|---|
| Title | "The Combination of Cysteine with Sugars" |
| Author | Maxwell P. Schubert — Laboratories of The Rockefeller Institute for Medical Research, New York |
| Venue | *J. Biol. Chem.* 1939, **130**(2), 601-603. Received July 21, 1939 (volume and issue not printed in the PDF; page numbers are) |
| DOI | 10.1016/s0021-9258(18)73532-1 (not printed; Crossref returns this title, this author, J. Biol. Chem. 130(2) 601-603, 1939-10; checked 2026-10-31). **Matches the brief.** |
| Paper type | Preparative note: synthesis, elemental analysis, qualitative reactivity |
| Companion | Schubert 1936, J. Biol. Chem. 114:341 (aldehydes; footnote 1). Not on disk |

## 2. Methods as they matter to a model

- **Sugars** [M]: D-arabinose, **D-xylose**, D-glucose, D-mannose, D-galactose, lactose. "No compound of fructose and
  cysteine has been isolated" [Q] (ketose: no thiazolidine, consistent with Cys reacting at the aldehyde).
- **Stoichiometry** [Q]: 1:1, "with elimination of a mole of water"; structure assigned as thiazolidine by analogy with
  the aldehyde adducts; the galactose adduct does not react with 2,4-DNPH in acetic acid (no free carbonyl) [Q]; the
  peracetylated mannose adduct (hexaacetate, N-acetylated) no longer absorbs iodine [Q].
- **Preparation** [M]: 20 mM (= **millimoles**, 1939 usage) each of cysteine hydrochloride and sugar shaken with 7 cc
  water (10 cc for lactose) — nominally ~2.9 mol/L each [Z: 20 mmol / 7 mL, before the solute volume], strongly
  acidic (Cys·HCl, no base); clear in ~15 min; 48 h at room temperature; then 20 mmol pyridine; 50-70 h more; 200 cc
  absolute ethanol; recrystallised from water + 150 cc ethanol (xylose and lactose: two or three more cycles).
- **Iodoacetate trap** [M]: 10 mmol galactose adduct + 10 mmol iodoacetic acid + 30 mmol NaHCO₃, ~5 h at room
  temperature (volume not printed), acidified, product S-carboxymethylcysteine. Yield not printed.

## 3. Numbers printed (none is a kinetic or thermodynamic constant)

Table I (read from the render) [M]:

| sugar | N found (%) | S found (%) | m.p. with dec. (°C) | recrystallised from water (cc) | yield (%) |
|---|---:|---:|---:|---|---:|
| calc. C₈H₁₅O₆NS (pentose) | 5.5 | 12.6 | | | |
| D-arabinose | 5.6 | 12.3 | 153 | 25 hot | 65 |
| D-xylose | 5.6 | 12.4 | 133 | 13 cold | 20 |
| calc. C₉H₁₇O₇NS (hexose) | 4.9 | 11.3 | | | |
| D-glucose | 5.1 | 11.5 | 167 | 35 hot | 40 |
| D-mannose | 4.9 | 11.2 | 171 | 50 hot | 63 |
| D-galactose | 4.9 | 10.8 | 138 | 45 hot | 58 |
| calc. C₁₅H₂₇O₁₂NS (lactose) | 3.1 | 7.2 | | | |
| lactose | 3.1 | 6.4 | 130 | 20 cold | 18 |

Author's caveat [Q]: melting points shift by up to 10 °C with drying / humidity. Mannose-adduct hexaacetate m.p.
150-152 °C, N 3.0 / S 6.0 % (calc. 2.6 / 6.0).

## 4. Figures

None.

## 5. What the repository takes

Nothing quantitative: no FIT row, no constant, no pot, no keyed table. Two qualitative statements may be cited as
context for the TTCA ring-chain step, never promoted to a prior: (i) xylose- and arabinose-cysteine thiazolidines
deliver their full thiol to iodine as fast as free cysteine at room temperature [Q] — the C-S ring opening is not
slow; (ii) in excess bicarbonate no free thiol is detectable [Q] — the closed form dominates at mildly alkaline pH.

## 6. Flags

1. **"mM" means millimoles**, not mmol/L (1930s J. Biol. Chem. usage); the synthesis is at ~2-3 mol/L, not 20 mM.
2. **No pH, concentration or time is given for any thiol test.** Do not back-compute a K from "very faint" or a rate
   from "as rapidly as".
3. **Trap ≠ dissociation.** Iodine and iodoacetate capture the open-chain thiol; they bound the ring opening, not the
   hydrolysis to free cysteine + sugar.
4. **The brief's expectation (dissociation constants) belongs to Schubert 1936, not this note.**
5. **Registry:** `cysteine`, `xylose` and `ttca` (the xylose-cysteine thiazolidine, engine species TTCA) are in
   `data/keys/compounds.yml`; the arabinose, hexose and lactose adducts are not registered.
