# L1 scoping: do the Maillard levers move the lipid aldehydes? (2026-11-14, literature only)

**Status: scoping note, not a sec. 7 case.** Nothing here changes the engine; the engine is paused (STATUS_AND_PLAN
row 6, STEP5_RECORD). It records what the corpus on disk says, so the case can be written when the pause ends.

## Why

In the P3 engine dry run (`results/validation/p3_engine_dry_run.md`) the lipid lane predicts hexanal, nonanal and
2-pentylfuran at exactly 0 dex for every lever: sugar identity, sulfur source and pH. The zero is structural. The lane
(`src/kinetic_core/lipid.py`) decomposes a pre-existing hydroperoxide pool and has no initiation term. Its only
coupling to the Maillard lanes, the aldehyde-lysine `COVALENT_SINK`, is inert by ruling
(FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2). `lane_coupling_verdict` considers lysine only, and never looks at
cysteine, H2S or pH.

## What the corpus on disk says (ratios: with the lever / without; dex = log10)

| source | system | lever | hexanal | others |
|---|---|---|---|---|
| farmer1992 | 4 lipids +/- cysteine + ribose, pH 5.6, 140 C 1 h | Cys + ribose | -0.08 to -0.49 | nonanal -0.31 to -1.11; 2,4-decadienal -0.23 to -1.30; 2-pentylfuran +0.34 to +0.51 (3 of 4) |
| elmore2002 | methyl linoleate +/- cysteine + ribose, pH 5.5, 140 C 30 min | Cys + ribose | -0.23 (not significant) | 2-enals and dienals 5x to >250x lower; 2-pentylfuran +0.74 |
| mandin1999 | linoleic acid +/- glucose or methionine, pH 6, ~100 C 2 h | glucose; methionine | +0.24; ~0 | decadienal +0.18; -0.6 |
| xin2026 | extruded plant-protein analogue, 8 carbohydrates + control | carbohydrate identity | spread <= 0.38 over all arms | nonanal spread 0.22 |
| madruga1995 | beef, pH 4.0 -> 5.6, 140 C 30 min | pH | ~ -0.14 per pH unit | 2-pentylfuran -0.52 overall |
| desousafontes2024 | goat hydrolysate + xylose + Cys + thiamine, pH 4 vs 6, 121 C 1 h | pH | +0.15 (opposite sign to Madruga) | nonanal +0.25; 2-pentylfuran +0.98 |
| zamora2015b | 2,4-decadienal -> hexanal | pH | yield peaks near pH 8 | |
| zhang1989 | 2,4-decadienal + cysteine, 180 C | cysteine | cysteine also FORMS hexanal (278.6 mg/mol vs 45.4 with glutathione) | |
| zhang2021, du2023, du2024 | aldehyde + cysteine or H2S | thiol sink | thiazolidines, 2-alkylthiophenes (figures only) | |
| fischer2021 | pea isolate, unheated | pH (release) | the isolate carries ~5 ug/g hexanal; release +59 % at pH 4.5 vs 6.5 | |

Nothing on disk measures Maillard reaction products acting as antioxidants on lipid oxidation in a heated pot.

## Reading (calibrated; the probability that the lever moves hexanal by more than 0.3 dex in the P3 pot)

- **Sulfur source (cysteine vs thiamine): P ~ 0.30, cysteine lower.** Thiazolidines and H2S trapping are specific to
  the thiol. Nonanal P ~ 0.45, cysteine lower; 2-pentylfuran P ~ 0.3, sign leaning cysteine-higher. A zero is not
  defensible here.
- **Sugar (glucose vs ribose): P ~ 0.15, sign unknown.** Xin 2026 bounds the spread across carbohydrates at 0.38 dex.
- **pH (6 vs 7): P ~ 0.12, sign unknown.** The two on-disk pH series disagree in sign.

Two things in this pot lower all three probabilities. The isolate's carried-over hexanal dominates, and the protein's
amines already make up most of the carbonyl sink, so swapping one additive changes less.

## What the case would be, when the pause ends

- A **thiol / H2S carbonyl sink**, carbonyl-class specific: saturated aldehydes ~2x, enals and dienals 5x - 80x,
  ketones and 1-octen-3-ol untouched (Farmer 1992's pattern). Its structure comes from elmore1997, zhang1989 and
  du2023/2024, with a rival in which cysteine also forms hexanal (zhang1989).
- Until it is built, the shipped zeros should carry declared bands, not point zeros:
  - sulfur source on hexanal [-0.5, 0] dex;
  - sugar and pH +/- 0.3 dex.

  This is a presentation change, made only when the pause ends.
- `lane_coupling_verdict` should name cysteine as a shared species alongside lysine.
- **The hold-out conflict comes first.** farmer1992 and elmore2002 are designated hold-outs in their dossiers, so a
  band set from them would be fitted to data read for scoring. Either the band comes from the mechanism papers alone,
  with those two kept as the test, or they are formally re-designated before use.
- P3 itself measures hexanal, nonanal and 2-pentylfuran on all three levers, so it is the designed test of this.
