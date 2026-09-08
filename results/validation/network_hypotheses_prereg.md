# Pre-registration: the hypothesis layer (2026-09-08)

## 1. What it is

A read-only layer that applies the literature's reaction rules (`data/lit/reaction_rules.yml`, one
SMIRKS per transformation a dossier draws, each with its anchor and with positive and negative
controls) to the engine's species (`data/species/structures.yml`) and places every proposed step
against the engine's reactions: **modelled** (the engine has the step, with a rate), **mechanism
known** (the source draws it; the engine has no step), or **proposed** (analogous to a cited rule, on
reactants no source shows). It writes `results/validation/network_hypotheses.{json,md}`.

It produces no rate and no concentration. The engine never reads it (`tests/unit/test_network_hypotheses.py`
asserts that nothing under `src/kinetic_core` imports it). Its job is to show the frontier: which
steps the chemistry allows that the model does not have, so that a refusal can say "no rate" rather
than "impossible", and so that the sink-structure work has a candidate list it did not write by hand.

## 2. Why now

The thiol-sink re-tuning (`kinetic_core_b16_prereg.md` section 6) found that the model removes the
thiols too fast at 100 °C and at 140 °C alike and that no single Arrhenius sink can serve both. The
next step needs the candidate sinks laid out: what the thiols can react with in the pot, which of
those steps the engine has, and which it lacks. A hand-curated table from the dossiers already read
(`docs/validation/thiol_sink_candidates.md`) is written alongside; this layer is the reproducible
version of it, and it generalises to the other lanes.

## 3. What runs

Six charges (`src/network_hypotheses/report.CHARGES`): each lane's reference pot at depth 2, and a
depth-1 probe of MFT and FFT with every carbonyl and thiol partner the pots hold. Depth is the number
of rule applications; products above 40 heavy atoms are dropped; a run that exceeds 400 products stops.

## 4. What counts as success, declared before the run

- **T1, positive control.** From pentose + cysteine the rules must rediscover the engine's own steps:
  cysteine to hydrogen sulfide is placed *modelled*, and the 1-deoxypentosone is reached within two
  steps. From the hydroperoxides, hexanal and pentane are reached and placed *modelled*.
- **T2, no nonsense.** Every proposed step conserves carbon except where the rule's source loses CO2
  (the Strecker rule, cysteine thermolysis, the acrylamide step); no product fails RDKit sanitisation.
- **T3, the frontier is visible.** The thiol probe places at least one thiol-removing step as
  *mechanism known* or *proposed* that the engine lacks (the engine's sinks are disulfide formation,
  the methanethiol coupling, an inert thioether with a lumped electrophile pool, protein disulfide
  exchange, oligomer and first-order decay; it has no thiol plus carbonyl step).
- **T4, controls.** Every rule's positive and negative controls pass (the unit test).

Failure of T1 or T2 means the rules are wrong and the artifact does not ship. Failure of T3 means the
rule set is too small to say anything about the sinks, and the artifact ships with that stated.

## 5. What it will not do

It will not rank candidates, estimate rates, or feed the fit. A candidate step becomes a model step
only through a pre-registered wave with a measured or declared rate.

## 6. Outcome (2026-09-08, first run; `results/validation/network_hypotheses.md`)

26 rules, 128 proposed steps over the six charges: 22 placed *modelled*, 64 *mechanism known*,
42 *proposed*; 30 products are engine species, 4 are registry compounds, 103 are neither.

- **T1 passed.** Pentose + cysteine: cysteine to hydrogen sulfide is *modelled* (`r_cys_h2s`) and the
  1-deoxypentosone is reached at depth 2 through the cysteine Amadori compound. Glucose + glycine:
  all four steps *modelled* (the Amadori rule matched through the Schiff base, `r_schiff > r_amadori`).
  The hydroperoxides: pentane, hexanal and the four oxo-esters and dienals all *modelled* against the
  lipid lane's position table.
- **T2 passed.** Every step conserves carbon except the three declared CO2 losses; no product failed
  sanitisation (the unit test).
- **T3 passed.** The thiol probe places as *mechanism known* or *proposed*, and the engine lacks:
  the mixed MFT-FFT disulfide and the thiol-cysteine mixed disulfides (R12); MFT and FFT addition to
  2,4-decadienal and to acrylamide (R13, the acrylamide + cysteine analogue is modelled); the
  hemithioacetals of both thiols with furfural, HMF, methylglyoxal, glyoxal, hexanal and the pentose
  (R14, proposed); the forward thiazolidine step from cysteine and the pentose (R15: the engine
  charges TTCA and models only its ring opening). These are the candidate sinks the hand-curated
  table (`docs/validation/thiol_sink_candidates.md`) lists with the literature's numbers.
- **T4 passed.** All 26 rules' positive and negative controls pass.

Two things the run also showed. The Amadori compound of cysteine itself (the first product of
pentose + cysteine by R01) is not an engine species, so the enolisation steps that follow are placed
*mechanism known* although the engine lumps them as pentose + cysteine to the deoxypentosones; the
placement is right about the species, not wrong about the chemistry. And beyond the first
generation the walk was restricted to known compounds and the sink rules made terminal, because the
unrestricted run produced adducts of adducts no source names (140 steps, most of them noise).

The artifact ships. It is read by `maillard explain` ("possible, not modelled") and by the wishlist
("reachable by a cited rule" against each refused compound).

**Addendum, 2026-09-08 evening.** Rule R28 (two alpha-aminoketones to a pyrazine, anchored to Zhou 2023's
mechanism figure) and a Strecker charge (the small dicarbonyls with alanine and cysteine). The Strecker
rule now puts the amine on the carbon that was the aldehyde, which is what the pyrazine condensation
needs (the first version put it on the substituted carbon). With them the layer reaches
2,5-dimethylpyrazine, trimethylpyrazine and tetramethylpyrazine within two steps, and `maillard
explain 2,5-dimethylpyrazine` now says the refusal means no rate, not no route. Nonanal from the
oleate hydroperoxides, 2-pentylfuran and 1-hexanol still have no cited rule: no dossier on disk draws
their route, and a rule without a source is not written. 27 rules, 201 steps, T1 to T4 unchanged.

**Addendum, 2026-09-08 late evening: the lipid routes.** Six lipid papers were read the same night
(`miyazaki2023`, `cao2020`, `chen2017`, `yang2024`, `choe2006`, `wanjala2021` extraction dossiers) and
four rules entered the table: R29 and R30, the two scissions of a mono-ene allylic hydroperoxide
(Cao 2020's drawn map for the oleate 8-, 9-, 10- and 11-hydroperoxides: nonanal, octanal, decanal and
their oxo-esters; 2-decenal and 2-undecenal with the alkane-ended esters; Chen 2017 corroborates
9-OOH to nonanal and 2-decenal on the free acid); R31, the furyl route of a conjugated-diene
hydroperoxide (Miyazaki 2023's proposed route from the 13-hydroperoxide to 2-pentylfuran plus
9-oxononanoate, and its mirror from the 9-hydroperoxide to hexanal; the 12-hydroperoxide gives no
2-pentylfuran in their table, and the rule's shape refuses it); R32, the 10-hydroperoxide to
1-octen-3-ol (Miyazaki 2023, Figure 3). The species these rules name that the engine does not hold
(the four oleate isomers the lipid lane lumps, the linoleate 10-hydroperoxide, the C8 to C11
aldehydes, the ester-side fragments, methyl stearate as every rule's negative control) live in
`data/species/literature_structures.yml`, loaded beside the engine's table and never read by the
engine. With them the layer reaches nonanal from an oleate charge and 2-pentylfuran and 1-octen-3-ol
from the linoleate charge, all placed *mechanism known, not modelled*; `maillard explain` now says
"no rate, not no route" for the three, and `maillard explain 1-hexanol` still says no route: none of
the six papers draws an aldehyde-to-alcohol step (the alcohols they see come from hydroperoxide
reduction or the allyl radical), so the rule is not written. 31 rules, 219 steps; T1 to T4 unchanged;
all 31 rules' controls pass. Four papers give four different parents for 2-pentylfuran (the
13-hydroperoxide here; the singlet-oxygen 10-hydroperoxide in Choe & Min's review; the
9-hydroperoxide in Yang 2024's map; 4-hydroxynonenal in Wanjala 2021): R31 is anchored to the one
paper that measured isomer by isomer, and the others are recorded on the rule.

One correction fell out of the controls: the registry's 1-octen-3-ol SMILES had seven carbons
against its own eight-carbon InChI (`data/species/off_flavour_targets.yml`); fixed at the source and
the registry regenerated (lesson in `tasks/lessons.md`).

*Later the same night: Yao 2024 (methyl linoleate at 180 °C) was read; it detects and quantifies the
dihydrofuran intermediate R31 assumes and draws R32's scission for the 10-hydroperoxide, so both rules
now cite it as corroboration beside Miyazaki 2023. Its twenty barriers are DFT and are recorded as
inadmissible. Rules and steps unchanged: 31 rules, 219 steps.*

