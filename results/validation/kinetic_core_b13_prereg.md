# Wave B13 pre-registration — the dicarbonyl trio on the trunk lane (programme step W4, part 2)

*Written 2026-09-07 before any B13 number was scored. Owner's instruction: "implement all the backlog."
Modules: `src/kinetic_core/parameters_dicarbonyl.py`, `network.DICARBONYL_REACTIONS`; licence:
FIT_HOLDOUT_DECLARATION.md Amendment 22.*

## 1. What is added, and from where

Three species and five steps on the TRUNK lane only, with the constants Kocadağlı & Gökmen 2016
(JAFC 64:6446) fitted to their amine-free glucose glass at 160 / 180 / 200 C, re-referenced from the
paper's T_b = 180 C exactly as wave B7 did for the furanic channel (`kocadagli2016jafc_extraction.md`
sec. 4, Table 2, glucose system):

| step | k_b at 180 C (x1e-3 /min, HPD) | Ea (kJ/mol, HPD) | carried as |
|---|---:|---:|---|
| glucose -> glucosone | 0.069 +/- 0.005 | 125.9 +/- 4.9 | measured |
| glucosone -> glyoxal + C4 | 737 +/- 58.9 | 93.8 +/- 6.4 | measured |
| 1-deoxyglucosone -> diacetyl + C2 | 12.2 +/- 1.12 | 150.8 +/- 8.8 | measured |
| glyoxal -> unassigned | 32.6 +/- 8.83 | 0 (FIXED by the authors) | declared, flagged |
| diacetyl -> unassigned | 0 +/- 0 | (blank) | zero, as a prediction |

Why these three: glyoxal is the CML precursor and methylglyoxal (already on the trunk) the CEL
precursor, and the panel refuses both AGE rows today for want of the species; diacetyl is the
buttery odorant and the 3-mercapto-2-butanone precursor Yiltirak 2026 quantifies; glucosone is the
only route to glyoxal in the source.

## 2. What does not change

- The sulfur network keeps exactly the topology wave B9 was fitted on (`network.REACTIONS`); the five
  steps run only when the trunk integrates on its own (`network.TRUNK_REACTIONS`). A request for
  glyoxal, glucosone or diacetyl that resolves to another lane is REFUSED by name, not answered with
  a silent zero. The acrylamide network builds its own set and is untouched.
- The three species are appended at the end of the trunk table: every existing index is unchanged.
- No fit row is added or read. The B1 pairs, the furanic block and B12's terms are untouched.
- Every panel row: the only trunk row (Steinhagen 2021, glucose alone at 121 C) gains the tiny
  glucose -> glucosone drain (1e-6 /min at 121 C) -- expected to move HMF by far less than 0.1 %.

## 3. Validation, honestly

No independent measurement of glyoxal, glucosone or diacetyl exists in the corpus with a charge the
trunk can run: Gürsul Aktağ 2020's juices carry sucrose and the juice's own amines at 27-37 C and its
own barriers are mostly negative; Kocadağlı 2016 Food Chem (wheat flour) prints constants, not
concentrations; Lee 2022 / 2024 (cake, glucose +/- leucine) hold the right data in figures under a
baking temperature profile and are the next validation set once digitised. B13 therefore ships the
trio as DECLARED, ONE-SOURCE, EXTRAPOLATED constants, flagged as such on every run, and the wishlist
names the measurement that would replace each (a glyoxal loss rate at two temperatures; a diacetyl
loss rate; glucosone in a glucose/glycine solution at 100-145 C).

## 4. Pre-registered checks

- The trunk's carbon balance closes with the five steps (construction-time invariant).
- Every existing prediction at the reference conditions reproduces within 0.1 % (pinned by test).
- A glucose/glycine pot at 120 C for 60 min predicts glyoxal and diacetyl as positive numbers on the
  trunk and is refused for them on the sulfur lane.
- The CML / CEL rows stay refused (the lysine adduct steps are not in this wave) -- the refusal text
  now says the precursor exists and the adduct step does not.

## 5. Outcome (2026-09-07, after the checks)

All four checks held. The trunk's carbon balance closes with the five steps (import-time invariant).
Every existing trunk prediction at the reference conditions reproduces within 1e-3 (test pin; the
Steinhagen 2021 row's HMF is unchanged at 1.46e3 ug/L, 11.93x). A glucose/glycine pot at 120 C for
60 min answers glyoxal 21 ug/L, glucosone 184 ug/L and diacetyl as positive numbers on the trunk, and a
ribose/cysteine pot is refused for glyoxal by name. The CML / CEL rows stay refused. **One finding the
wave did not expect:** at 120 C the source's own constants make glucosone ACCUMULATE above glyoxal
(its onward step runs at 0.016 /min there against a 0.033 /min glyoxal sink at every temperature),
because the glyoxal sink carries the authors' fixed-zero barrier. That is a testable prediction and
the first item the wishlist names. Scorecard 4/39 and 3/38, directional 18/30 and envelope coverage
5/33 unchanged; the envelope's median width moved 1.369 -> 1.338 dex because the two B12 trunk bands
now join the draw table (they are inert on every panel row today).
