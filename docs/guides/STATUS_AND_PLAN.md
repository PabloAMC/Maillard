# Status and plan

*The one plan for this repository. Rewritten 2026-11-02 for readability; the content is the 2026-11-01 plan as
refined on 2026-11-02. The previous text is kept verbatim in
[docs/history/STATUS_AND_PLAN_2026-11-01.md](../history/STATUS_AND_PLAN_2026-11-01.md). Older records that cite
"STATUS_AND_PLAN sec. 5", "rule (a)", "B0", "roadmap row 2" or "triage item 3" refer to that file.*

**Contents.**
1. Who this is for and what they need
2. Where the model stands
3. Why it fails: five root problems
4. The plan at a glance
5. The workstreams
6. Timeline: what runs in parallel
7. How a change ships
8. Checkpoints and odds
9. Not in the plan
10. Glossary of codes

---

## 1. Who this is for and what they need

The users are alternative-protein scientists who formulate meaty flavour on pea, soy or similar proteins. They ask
questions like these:

- "If I add 0.5 % yeast extract, or swap ribose for glucose, does the meaty note go up or down, and by how much?"
- "What does raising the pH, adding fat, or cooking longer do to the meaty, roasted and off-notes?"
- "Which change is worth testing in the lab first?"

So the tool must:

1. **Get the direction and rough size of each lever right** for the key odorants. The odorants are the meaty thiols
   (MFT, FFT, H2S), Strecker aldehydes, pyrazines and lipid aldehydes.
2. **Accept real recipes and real cooking.** Ingredients such as protein isolate, yeast extract and oil. Cooking such
   as a pan-fried patty, not only a sealed vial of buffer.
3. **Say honestly how sure it is,** for a laboratory it has never seen.
4. **Improve on a lab's own data** (the `calibrate` verb, which exists).

Absolute levels that hold across laboratories within 3x are NOT a near-term goal. Laboratories disagree by 10x to
1000x on the same pot; no model removes that without designed data.

**Owner decisions, 2026-11-02.** Work continues without new laboratory data; the designed-experiment brief is paused.
The plan's priority is mechanisms that make the model better, not more process.

## 2. Where the model stands

Shipped calibration B60. Numbers from `results/validation/core_panel_scores.md`, `contrast_scores.md` and EX-056.

| question | today |
|---|---|
| absolute level within 3x, out of sample | 10 / 23 rows (median miss 4.4x; the thiol levels refused by H0 and the levels the falsified pH balance moves by S3 are out of the count); the "typical yield" baseline does better on the rows both answer (5 / 8 against 4 / 8) |
| ratio between two pots of one study | LOSES to "no change": 19 out-of-sample pairs in 4 studies, median error 0.30 dex against 0.15 for "no change", direction right 2 / 10 (after the S3 refusal of 2026-11-07; 37 pairs, 0.31 against 0.28, 12 / 25 before) |
| meaty sulfur odorants | 10 / 40 rows within 3x; thiol leave-one-lab-out NOT ADEQUATE (90 % intervals cover 0.34) |
| acrylamide, lipid aldehydes | the null model wins |
| pyrazines | about 3 orders of magnitude too low |
| Strecker aldehydes, methional, 2-acetyl-1-pyrroline | absent or refused |
| a real recipe with both sulfur and acrylamide chemistry | refused (the lanes do not combine) |

What works today: the experiment ranking (`rank`), per-lab calibration (`calibrate`), honest refusals, and an
evidence base of about 1,800 two-reader-verified fit rows, a 32-family pathway atlas and 20 ranked engine errors
([ENGINE_VS_LITERATURE.md](../reference/pathway_atlas/ENGINE_VS_LITERATURE.md)).

## 3. Why it fails: five root problems

Every task in this plan exists to fix one of these.

| # | root problem | example | fixed by |
|---|---|---|---|
| R1 | **Wrong mechanisms in the meaty sulfur chemistry** | the model destroys thiols far faster than any pot; MFT runs 99.9 % through a route labelling puts at ~7 % | workstream W2, steps M1, M2 |
| R2 | **Missing chemistry** | no Strecker aldehydes; glucose makes almost no MFT; pyrazines far too low (a dicarbonyl-supply and amine-identity deficit, not a sink: M6 case) | W2, steps M3, M4 (supply and amine classes), M6 (the sinks) |
| R3 | **Lab pots do not transfer to food** | constants measured in phosphate buffer at fixed pH; foods are unbuffered, fatty, drier, with protein as the main amine | W2 steps M1b (the pH balance and the pentose-side pH laws), M4, M5; W5 |
| R4 | **The software cannot represent a whole recipe, and fits are slow** | four separate lanes that do not combine; one fit evaluation takes 2.5 - 21 s | W3 |
| R5 | **Uncertainty is not honest for a new laboratory** | fitted constants absorb laboratory offsets; out-of-sample coverage 0.34 | W1 (H0), W4 |

## 4. The plan at a glance

Six workstreams. W1 measures progress. W2 is the core: it fixes the chemistry. W3 and W4 give the engine and its
uncertainty. W5 turns the engine into a tool for a formulator. W6 is the owner's outside-data requests.

```mermaid
graph TD
    W1["W1 Yardstick<br/>contrast benchmark, H0"]
    M0["M0 Misfit inventory<br/>+ seen-data checks"]
    M1["M1 run (not separated); M1.1 measurement model;<br/>M1c supply structure; M1b pH balance + pentose pH"]
    M2["M2 MFT routes, in one refit<br/>with M1c (never at instant supply)"]
    M3["M3 Strecker aldehydes"]
    W3["W3 One network (A2)<br/>+ fast engine (H3)"]
    M4["M4 Glucose to MFT,<br/>protein amines"]
    M5["M5 Buffer, fat,<br/>water activity"]
    M6["M6 Pyrazines"]
    JR["Joint refit"]
    H1["H1 Hierarchical model spec"]
    H2["H2 Tool bake-off"]
    H4["H4 Full posterior"]
    W5a["W5 Yeast extract,<br/>crust, heme"]
    LT["Lever table<br/>(the deliverable)"]

    W1 --> M1 --> M2 --> JR
    M0 --> M4
    M0 --> M5
    M0 --> M6
    M0 --> JR
    M3 --> W3
    M3 --> H2 --> W3
    W3 --> M4 --> M5 --> M6 --> JR
    W1 --> H1
    M2 --> H1 --> H4
    W3 --> H4
    JR --> H4 --> LT
    W3 --> W5a --> LT
```

Read it as three lines of work that meet at the end:
- **the chemistry line** (M1 - M6, then one joint refit);
- **the engine and uncertainty line** (one network, fast engine, hierarchical posterior);
- **the product line** (ingredients, cooking, the lever table).

## 5. The workstreams

### W1. Yardstick: know whether a change helps

| task | what | state | effort |
|---|---|---|---|
| W1.1 | **Contrast benchmark scorer** on the frozen split v2: does the engine get within-study contrasts (pH, time, cysteine +/-, phosphate) right, against two null models | **done** (PR #48, 2026-11-03): prereg v6 frozen, scored once (`results/validation/kinetic_core_contrast_benchmark_score_v1.md`). Direction not distinguishable from a coin (weighted S 0.43, 90 % 0.19 - 0.77); log-ratio error not distinguishable from predicting no change (median 0.66 vs 0.41 dex). PROVISIONAL (3 labs) | -- |
| W1.2 | **H0: between-laboratory error out of sample** from the thiol leave-one-lab-out (EX-056). Sets how wide intervals must be for a new lab | **scored** (EX-059, PR #53): calibrated, not shown better than a null; tau_ABS 0.964 dex, D2 branch (i). **Applied** 2026-11-05 (absolute thiol levels refused for a new laboratory, ratios with omega 0.419 in their intervals; in review). **Lifting the refusal** needs absolute MFT / FFT from at least two new laboratories, re-scored under H0's rule (H0' is triggered by the first batch that keys >= 2 new laboratories, H0 prereg sec. 9; wording aligned 2026-11-06): besides W6 replies and the WUR blind test, a literature route is open. A forward-citation sweep of the SIDA thiol papers (723 citing works) found Vauthey 2000 (Nestec, 10.1021/jf991254a; cysteine + ribose with an aqueous control) and Lopez 2019 (10.1021/acs.jafc.9b04153; a mushroom-protein matrix, likely not simulable); Vauthey's roles were declared before its data were read (M1c amendment A1). A qualifying laboratory counts toward that trigger; once two have qualified, a pre-registered H0 amendment re-scores tau_ABS; owner sign-off, since shipped outputs move | -- |
| W1.3 | **Profile one fit evaluation** (Jacobian, pH solve, integrator) before the next long refit | **done** (PR #55): ~80 % of pot time is the Python rate rebuild + pH solves; the FD Jacobian ~90 % of fit time | -- |
| W1.4 | **Literature tooling: `scripts/tools/openalex.py`** (added 2026-11-05): `verify` checks every `data/keys/papers.yml` entry against OpenAlex (authors, year, venue, DOI resolution, `is_retracted`; report only, never writes `data/`); `sweep` runs forward / backward citations from seeds, deduplicated against the corpus (~1 req/s, no account); `labs` gives institution ids as a first pass for laboratory assignment (a person confirms group moves). `verify` would have caught an invented first author and the crossed Bork files; laboratory identity is what H0's tau and H1's offsets rest on | run `verify` now and before each benchmark re-score; `sweep` per case (sec. 7 item 2) and before each W6 request; `labs` feeds H0 / H1 | 1 - 2 h build; verify ~15 min unattended |

**What the benchmark can and cannot judge (2026-11-04).** W1.1 answers 4 of the 65 held-out cells; the rest are
refused (an operation the engine lacks 23, inputs missing 17, a species it does not model 12, a precursor it cannot map
7, other 2). On what it does answer, it would detect a 70 % per-step agreement in direction only 43 % of the time (its detection
curve). So it cannot show on its own that a W2 step helped:
- **A step's verdict rests on its own decisive data** (sec. 7, item 4), never on the benchmark score.
- **The benchmark is the side-effect record** (sec. 7, item 7): every held-out row, before and after each shipped step.
- **The January re-score is descriptive.** It reports coverage first, then the direction and log-ratio scores on the items
  answerable both times. It is read as evidence about the plan as a whole only where coverage has grown, for instance
  the Strecker and methional items after M3 (today the engine predicts no change on all 20 methionine items) or the yeast-extract
  pots after W5.1.

### W2. Chemistry: fix the mechanisms (the core)

Ordered by expected gain for plant-based odorants. Each step names its evidence and its decisive data; how a step ships is sec. 7. "E" numbers
refer to the ranked errors in ENGINE_VS_LITERATURE sec. 3 (and its amendments in secs. 5 and 7).

| step | mechanism | what to build | evidence | decisive data (sec. 7, item 4) | effort |
|---|---|---|---|---|---|
| **M0** | **Misfit inventory and seen-data checks** (added 2026-11-05: the misfit protocol, [OPERATING_PROCEDURE](OPERATING_PROCEDURE.md) sec. 4, applied across the whole model instead of case by case) | One hand-kept inventory (`docs/reference/misfit_inventory.md`, draft PR) of every open structural misfit: size and sign, the protocol step reached (integrity, conditions, measurement, mechanism, unresolved), candidate mechanisms already on record, the shared structure it touches, the cheapest check on already-seen data, the data that check spends, and its leverage on the lever table. Misfits are grouped by shared structure (one piece explaining >= 2 misfits is the protocol's best test of a mechanism: it predicts beyond the misfit it came from). Every misfit's leading candidate mechanisms are checked (owner, 2026-11-05), in the compute waits, cheapest discriminating check first and grouped so one check tests a shared piece once; the first is the M3 argon diagnostic (every Strecker rival ~5x high on Hofmann 2000b Table 2's argon PA + PAA). **A check never ships anything:** it nominates rivals for a sec. 7 case, tested on data the check did not spend (sec. 7, item 2). Each check writes its pass rule before computing, and every outcome is recorded, failures included (no selective reporting across many checks). Expected yield: most checks narrow the candidates rather than settle them (P ~0.25 that a given check leaves one clear candidate); P ~0.5 that at least one shared piece explains >= 2 misfits. Its output orders and shapes the M4 - M6 cases and the joint refit | The registers' engine misfits (CF-012 residual, CF-084, CF-095, CF-122, the Bork norfuraneol decay, ...), the diagnostics report, the scorecards, the slice profile, the M1 / M1c / M3 records | -- (a triage, not a verdict) | inventory ~3 h; ~1 - 2 h per check |
| **M1** | **Thiol supply and loss, judged as combinations** (E10, E12; rescoped 2026-11-04 from its [case](../reference/cases/M1_thiol_loss_oxidant_ph_drift.md) and CF-107) | **M1 run (EX-061, 2026-11-05): NOT SEPARATED** under the frozen rule, but artefact-driven (a linearised floor sensitivity; evaluated exactly the reservoir favours C1). The data: a reservoir is real (the incumbent loses ~99 % of cysteine and pentose within 5 min at 120 C), C1's is 3 - 5x too slow and loses the FFT shape; before M2 fix the measurement model (no between-lab tau on C/C0 series, floors from the assay, exact densities) and free the reservoir's release rate ([review](../../results/validation/kinetic_core_m1_results_review.md)). **Then (2026-11-05) the supply screen** (`kinetic_core_supply_structure_signature_screen.md`; plug-in runs, no fits): no supply structure (C1; a pH-dependent exchange, whose slope Cao's drift bounds at m <= 0.18; a flat K plus a 1:1 drain) reproduces Zhai 2021 Fig. 4b (TTCA lost ~10x faster at pH 7 - 8 than at 5.5) or Zhang 2025b's mass balance (~19 mM of each reactant in no measured pool); a base-favoured TTCA exit passes the first, nothing inside non-Zhang bands the second. Supply moves to **M1c**; no Cys / Xyl label swap in Zhang 2025b (checked on the PDF). **M1.0 done (EX-060, 2026-11-04):** G1 holds at B60, more strongly than at B9 (at most 0.009 % of the oxidant used; phi = 1 adds at most 0.08 %), so no oxidant source enters M1. In Cao 2017 free cysteine and xylose fall within 30 min to a plateau (free Cys 0.06 - 0.08 of 0.4 mol/L, 15 - 20 % of the charge; corrected 2026-11-05) but the TTCA / cysteine-Amadori reservoir still holds about all of both after 5 h (1.06 / 1.09 of the charge), where R0 has no reservoir and burns the pentose in 1.8 s: the measured intermediates require C1. On Schieberle 2000, R0's MFT 360 / 30 and 720 / 30 ratios miss by -1.09 and -1.42 dex (peak at 2.6 h). **M1.0 as planned:** does B27's gate G1 still hold at B60 (the dimer step uses < 1 % of the oxidant, so no oxidant SOURCE can move the disulfide share; P ~0.85 it holds); the engine at B60 on the decisive data below. Then three combinations, each with `k_ttca_cys` in the structure (Cao 2017's K, not frozen) and the frozen pot-independent sink `k_mft_decay` freed or banded: **C1** slow cysteine-Amadori / TTCA supply (B59b); **C2** C1 + a loss that free cysteine slows (EX-050's form); **C3** C2 + the disulfide side: Guth's release declared, not fitted, and the deoxypentosones as an oxidant POOL with a declared pH law (not B27's source). Out, with reasons: B17a (Stack's measured K(T) makes it inert at process temperature); B27 as built (G1 arithmetic); B11's O2 (no within-study O2 arm at process temperature; Xu 2010 at 140 C finds air vs N2 changes neither thiol nor disulfide) | CF-107: four refused sink waves (B16, B17b, B17a, B25) are one fit optimum in a degenerate state, 77 % of its cost from TTCA rows on the frozen `k_ttca_cys`, best starts unconverged; against instant supply and a frozen fast sink, every added sink structure was switched off. Whole-model refusals are therefore not evidence against C1 - C3. B16 with its barrier ceiling lifted (converged) still turns MFT over at 6 h: the 100 C rise is a supply problem. Guth 1995 measures disulfide release far faster than B17b's fitted floor | Schieberle 2000's own ratios (MFT 360 / 30 and 720 / 360 min; measured still rising) separate the incumbent from C1 - C3; free cysteine and sugar over time (Cao 2017, Zhang 2025b, Zhai 2019) separate C1 from C2; the disulfide share's pH slope inside Zhou 2023, Zhang 2024's cystine vs cysteine arms and Whitfield 1999's fed pots judge C3 | M1.0 ~3 - 4 h; C1 - C3 2 - 3 days |
| **M1.1** | **Measurement model and fit driver, before any thiol refit** (new 2026-11-05, from M1's review and the premortem of the re-sequence) | **Superseded in part by sec. 7 item 5 (2026-11-06):** omega leaves the units' covariance; units no longer favour or veto one by one (one pooled T); 'means differ by more than omega' no longer substitutes for adequacy; the rest stands. The parked M2 prereg must conform before it freezes. No between-lab tau on C/C0 or ratio series (tau on absolute levels only); floors from each assay's reported LOD / LOQ, never a data point (only Zhai 2021 prints one: elsewhere a declared transfer, run as a sensitivity); H0's within-lab omega (0.419 dex: the worse of H0's two sets governs, prereg sec. 3.6 / D2 and A1 item 3; the reset set gives 0.390) in every unit's primary covariance; exact or Monte-Carlo densities on every form that can flip a sign (PSIS k-hat < 0.7, ESS >= 100); a unit favours or vetoes only where the favoured rival is adequate (posterior-predictive p > 0.05) or the means differ by more than omega, else "both inadequate"; Holm across guard groups; unidentified coordinates integrated over their band; central-difference polish (the forward-difference Jacobian is up to 7 % off, bake-off); asserted convergence. The Cao 2017 pot-input audit is done (2026-11-05): the keying is correct (free-base cysteine, Na phosphate set before mixing, 25 C label); a reported 0.8-unit gap compared the engine's in-situ (90 C) pH with Cao's cooled readings. On the cooled scale the engine starts at 5.50 like Cao and, under C1, reads ~0.2 low at 30 min, all of it 7.9 mM of cysteine carboxyl with its amine declared destroyed (CF-084 explanation 1); Cao stays a descriptive cross-check, not a control. Every pH comparison states its scale (cooled or in situ) | M1's results review (three MAJOR findings); the M2 draft's adversarial review (BLOCKER 1); the JAX bake-off's gradient check | -- (infrastructure: tested against M1's records, which it must reproduce under the old settings) | 8 - 12 h |
| **M1c** | **Thiol supply structure** (new 2026-11-05, from the supply screen) | A base-favoured TTCA exit, plus an account of Zhang 2025b's ~19 mM pool (open-chain or Schiff-base TTCA, cystine, the cysteine Amadori compound onward). Pool members enter only with a measured or same-class constant or a declared band, never a free pool sized to the gap (that would be CF-107 again). Items 1 - 5 frozen before any newly found paper is read: the screen's pots (Cao 2017, Zhai 2021, Zhang 2025b) are declared seen, so the decisive data must be others, or the verdict is development-grade and the decider is designed data (W6). Screened first on four signatures: (a) Cao's flat K under drift, (b) Zhai's pH series, (c) Zhang's mass balance, (d) the FFT shape (Schieberle; B60's FFT rows; C1 lost it). Runs at the measured pH trajectory where one exists, whatever M1b(i) finds. **Stop rule:** if no structure passes (a), (b) and (d) within ~25 h, supply is "not identifiable on disk": the R0 / C1 envelope stays under H0's relative-only ruling, and the designed-data request goes out | The screen record; Cao 2017 (K, pool closure, drift); Zhai 2021 Fig. 4b; Zhang 2025b (development data now: its blank pot was seen under the engine on 2026-11-05); M1's review; the literature search of 2026-11-05 (ring-chain pH law; norfuraneol decay at pH 5 - 6.5; the missing pool; assay LODs) | To be named in the case before the search's papers are read; P ~0.6 that nothing unseen and decisive exists on disk | case 3 - 4 h; build + screen 10 - 18 h; fit 4 - 6 h + 6 - 12 h compute |
| **M1b** | **pH: the pot's balance (CF-084) and the pentose-side steps** (E16, E2; rescoped 2026-11-04 from its [case](../reference/cases/M1b_sugar_step_ph_sulfur_lane.md)) | **Re-sequenced 2026-11-05:** part (i) runs first and on its own, timeboxed at ~14 h, judged on CF-084's own controls (Meynier, Hou, Zhou; Cao's trajectory is NOT a control, since TTCA and the cysteine Amadori compound are titratable and their budget is M1c's question); it does not gate M1c; its stop rule: if neither rival passes both controls, "balance not identified" and the measured-trajectory operator stays; shipping it moves pH everywhere, so the owner signs off first. Part (ii)'s chosen rival enters the joint refit (its norfuraneol decay and 3-deoxypentosone removal move M2's decisive pot-D share), or is declared deferred with that effect stated. **E16 (2026-11-05): undecidable on disk**, four blockers (the O2 status; unbuffered pH drift; no temperature series for norfuraneol alone; norfuraneol self-reaction with its own fragments at 20 mM). Bork 2022 / 2023 (Kanzler group, a new laboratory) set the joint refit's acid-pH decay band (4.4 - 6.1e-3 /min at 130 C, ~20 mM, unbuffered, initial pH 5, O2 not stated) and show the shipped pH-free decay 3.1 - 4.3x too slow there (a CF candidate). Bork 2024 / 2024b and Kanzler 2026 are being extracted (P ~0.15 they make E16 decidable: they would need norfuraneol alone at a second temperature or a known O2 status); Wietstock 2016 bears on M3's O2 budget. (i) The charge-balance pH scores pots only once its two controls pass, with CF-084's two rivals built and compared: the released ammonia survives (less acid), or the acid forms after the Amadori product is gone. (ii) pH laws on the pentose-side steps for the rival the data need: base-catalysed norfuraneol decay (E16; B61's prereg, its fit stopped and owed on B64), the enolisation split carried to pentoses (E2), a pH-dependent pentose + amine entry, or a base-favoured 3-deoxypentosone removal (EX-045, EX-046). Not built: `ph.trunk_terms` on its own, which slopes only the three hexose Amadori steps and is inert in pentose pots; E2's hexose part waits with E2 | The shipped engine's norfuraneol and furfural are pH-free in glycine-pentose pots (identical to 14 digits on Meynier's ladder; the pentose + amine entry is off, the amine-free pentose steps and the norfuraneol decay carry no pH). CF-084: the balance fails both controls. EX-001 - EX-003: Martins' slopes on the acids (acetic transferred, formic did not) | For E16, fed norfuraneol loss vs pH does not exist on disk: each laboratory fed it at one pH (Mikami 2017 and Nakamura 2020, one laboratory, pH 8, 100 C, and Mikami is B61's anchor; Bork 2022 / 2023, initial pH 5, 130 C) (corrected 2026-11-05); the 3- / 1-deoxypentosone ratio vs pH inside Zhai 2021 and Hou 2026 (R3En, the structure scored in EX-033, has the wrong sign; an R3En result, not B60's; corrected 2026-11-05); acetic acid vs pH (Davidek 2006a, 2008); the pH controls (Meynier's drift <= 0.2, Hou's trajectory). Meynier's norfuraneol ladder stays held-out (item 7 only). The split and the removal rivals cannot be separated without an in-situ pentose trap across pH (W6) | 2 - 3 days |
| **M2** | **MFT route shares** (E6 + E11), in ONE refit with M1c's supply structure and M1b(ii)'s chosen rival | Per-route barriers bounded by measured classes, fitted together with the supply structure the data favour (the steps share rows, and their errors stack). **Never at instant supply** (half of CF-107's compensator; the drafted design put P(anything ships) at ~0.1 there). The route design is drafted and parked (`kinetic_core_m1m2_supply_routes_refit_prereg.md`, v2: it supersedes the 2026-11-01 route prereg, moves Cerny 2007's 145 C labelled shares in-sample, holds the furfural barriers at the shipped values in every rival, and needs a synthetic power check before it freezes: under omega the shipped engine misses pot D by 6.9 sd but Hofmann by only 0.6 sd) | labelling in one laboratory (Cerny, Nestlé); fed-intermediate comparisons in two (de Roos 2005 / Laarhoven at Givaudan, reported only: `simulable: false`; Hofmann) (corrected 2026-11-05) | Cerny pot D norfuraneol share (measured 0.07, engine 0.98) and MP3P share; Hofmann 145 C fed pots; Zhai 2021's fed-TTCA 140 C MFT time-course shapes (E11; shared 120 C denominators with B60 declared) | 2 days + ~1 day on >= 8 cores |
| **M3** | **Strecker aldehydes** (E3; rescoped 2026-11-04 from its [case](../reference/cases/M3_strecker_aldehydes.md)) | **Phase 1 (2026-11-05): undecidable on disk, not frozen; moved to W6.** The forced-input layer was built and checked (EX-049 reproduced bit-exact; mass balances to ~1e-14; screen: no rival killed). Its power check, under the prereg's own rule that an unresolved density reads as a tie, gave every contrast zero power; amendment v2.1, made before any decisive data were read, added the Cu-law-free sub-units U3a / U3b (the prereg's fallback) and an adaptive-tempering SMC density estimator (validated against an exact grid and a 32-dimension known answer; reusable for other ~30 - 35-coordinate held-out densities). Still no decisive contrast reaches usable power (best: the route pair crediting R1 at 0.065; transfer and decay law ~0): P > 0.9 a run ends "not separated". The Hofmann 2000b Figs. 2 - 3 keyed rows were never read and stay unused. Open diagnostic (M0): every rival over-predicts Table 2's argon PA + PAA ~5x (median 0.39 - 0.48 vs 0.09 mM), so the shared structure (route A's rate law carried from fed-dicarbonyl pots, or the anoxic donor's efficiency) is suspect; a case for it follows only if the seen-data check nominates a rival. As designed before: first the prereg EX-051 asked for and that does not yet exist: one ARP-Phe decay shared by the fed-ARP and glucose pots, and two routes (dicarbonyl + amino acid; the Amadori product oxidised directly). Rivals: the current engine (it refuses); the claim; the dicarbonyl route only, with unmeasured or O2-made osones; both routes with a pot-specific decay (EX-049's optima); the claim with the decay moving with pH, buffer and order. Then feed the layer from the network; that phase shares the pyrazine steps (`r_go_ak`, `r_mgo_ak`), so it moves shipped pyrazines: owner first | Fed intermediates: Hofmann 2000b, Weenen 2001, Deng 2022 (the only second laboratory); measured ARP decay (Martins 2003, Deng 2021). The whole-model refusals (B22, B22b, B24, B29) are weak evidence: compensating errors named in the case. Cremer (v1 TEST) and Balagiannis (v2 TEST) carry no route marker and no donor: item 7 only | The fed ARP's early aldehyde shape (no lag if the Amadori route runs, a lag if only dicarbonyls do); fed vs in-pot ARP decay (no amino-acid ARP decay measured at 100 C and <= 0.1 M phosphate reaches EX-049's 0.058 - 0.069 /min; a linear phosphate law carrying Deng 2021's 0.1 M arm to Hofmann's 0.5 M is weak support at most, not inconsistent only at Deng's pH 7.5, CF-120; corrected 2026-11-05); the fed ARP's O2 mass balance (the dicarbonyl-only rival needs >= ~0.5 - 0.7 mol free osone per mol ARP; needs Hofmann 2000b Figs. 2 - 3 keyed); PA / PAA ratios. The TUM Becker contrast (F16#14, TEST in v1 and v2) measures no donor or ARP, sits at pH 5.2 outside the layer's scope, and is read only in item 7. The dicarbonyl-only and two-route rivals rest on one laboratory (DFA) | prereg 3 - 4 h; keying Hofmann 2000b Figs. 2 - 3 2 - 3 h; layer 2 - 3 days |
| **M4** | **Glucose to MFT; amine classes** (E1, E8; rescoped 2026-11-04 from its [case](../reference/cases/M4_glucose_mft_protein_amine.md)) | **Re-planned 2026-11-06 (independent premortem; case on PR #90):** on disk part G is development grade. D3's label shares falsify C2 + C3 recombination as the hexose MFT route (the incumbent's route and R3) and R4-as-written is falsified, but no rival is credited (the falsify-never-credit rule, sec. 7 item 2); the claim is restated as "a hexose C5 entry into 1-DPO, carbon origin (C-1 retained vs lost) undecided"; part A (amine classes, absorbing case M0a, PR #91) is an engine-structure requirement with no decisive data on disk. So M4 is NOT built on the network now: its structure goes into W3 at parity and off by default (one cysteine-sugar Amadori class module shared with M1c; shared 1-DG / glucosone / MGO pools; amine classes as species with their Amadori chains; one pH module on a declared scale; an O2 / headspace hook; R1 / R2 / R5 as declared-zero steps), a reduced forced-input model sizes the designed experiment's arms, the case freezes as "development grade, decider W6", and items 6 - 9 wait for designed data (decisive_measurements entry glucose_mft_carbon_origin). As designed before: On the one network (W3): glucose + cysteine -> the cysteine-glucose Amadori compound -> 1-deoxyglucosone with cysteine released -> the shared pentose MFT steps (part G); amine classes as species (free amino acid, cysteine, peptide, protein lysine) with their competition declared (part A). Rivals: the incumbent; the claim; a route that bypasses norfuraneol; C2 + C3 recombination; the hexose -> C5 entry switched on alone. The amine-catalysed enolisation (E8) waits: switching on M1's cysteine-Amadori route already retires the steps it would act on, and its development evidence at process temperature is too thin to fit | Hofmann Table 5 (one lab, 145 C): norfuraneol / MFT 7.1 in the glucose pot vs 2754 for ribose; glucose MFT vs pH peaks at pH 5, the opposite of ribose; labelling rejects C2 + C3 recombination by a measured share; de Roos: the Cys-glucose Amadori compound forms yet gives < 0.05 mg MFT per mol Cys at ~100 C. The B8 - B51.1 fits are whole-model, weak evidence; the red-team probe of the hexose entry inherits E6 / E11 | Within-pot shapes that cancel scale: glucose's MFT / norfuraneol ratio against ribose's (390x its norfuraneol ratio; separates the claim from the bypass), and the pH profile of glucose MFT. The plan's earlier ribose / glucose end-point ratios (10.4, 4.3) are moved by E6, E11 and the H2S budget too and the ribose arm is a B60 fit row: side evidence. No fed hexose intermediate to MFT or positional glucose label exists | 1 week (after W3) |
| **M5** | **Buffer, fat, water activity** (E5, E7; rescoped 2026-11-04 from its [case](../reference/cases/M5_buffer_fat_water_activity.md)) | Buffer: a per-step phosphate factor, linear with an intercept, against no factor / saturating / pH-only / step-specific rivals (B53 never ran; the form in EX-058 rested on TEST sources, so it is not settled). Fat: the thiols' loss to lipid, as partition (the plan's earlier claim) against lipid-electrophile capture of MFT (CF-031) and both. Water activity: B54's design (whole-phase concentration basis, Tg gate) against the shipped table, concentration only and step-specific | Davidek's pH-stat (x2.8 - 3.6 at held pH) rejects the pH-only rival; Bell's 0.5 vs 0.2 M ratio (1.74) rejects saturation with K <= 0.1 M; Mottram 2002: MFT x2.2 - 2.6 with FFT unchanged. The lipid fall is MFT-selective (x13 in Elmore; x9 to x2.2 in Yiltirak with fresh oil), which partition alone barely produces: the data favour capture | Buffer: dose ladders at held pH (Bell, Zyzak, Guan); fat: the MFT / FFT selectivity of the lipid fall (scale-free) and Zhang 2022's phase-resolved pot; water: Hofmann dry vs wet (ribose FFT x8, glucose x0.5). Potman and Kallen (v2 E5 TEST), Kim 1996, Zamora 2015, Hidalgo 1993 (v2 E7 TEST) and Xu 2011 (E12 TEST): item 7 only. No measured oil-water partition coefficient for MFT or FFT, and no MFT water-activity series beyond Modi 2008's peak areas | 1 week |
| **M6** | **Pyrazines: the sinks, and where the deficit really is** (E13, E14, CF-095; rescoped 2026-11-04 from its [case](../reference/cases/M6_pyrazines.md)) | A pH-flat dicarbonyl + amino-acid adduct sink; a finite aminoketone condensation with a measured-class constant; aminoacetone + MGO -> pyrrole as a competing sink. Rivals: the incumbent; the claim; the adduct sink only; a ring that never passes through free aminoketone; the deficit is supply and amine identity. Both sink structures LOWER pyrazines, while B18 T3 / T4 put the 2.9-dex deficit on the trunk's dicarbonyl supply in water and on glycine standing for lysine: the deficit is M4 / W3 supply and amine classes, not M6 | Piloty: 82 - 86 % amino-acid consumption vs 1 - 12 % Strecker products; Zhou 2024's zero-order rise from t = 0, where the 5-aminolevulinic-acid analogue implies a ~115 h lag (a lag < 30 min needs a constant ~5e4x the analogue): the measured analogue rejects the finite condensation as built, unless the ring bypasses free aminoketone. No aminoketone has ever been measured in a pot | Zhou 2024's early time course; Adams' pyrrole and DMP against MGO excess and pH; Piloty across designs; EX-021's 2.3-dex gap shared with M3. Leahy's pH ratios are B18 calibration inputs and the four E13 TEST studies measure the mixing ratio, not the sinks: item 7 only (their scores are still owed) | 3 - 4 days |
| **JR** | **Joint refit** | One refit after M1 - M6, not one per step | -- | all of the above, together | ~1 day + compute |

**Already decided, not in the steps:** E4 (pyrazine pH law) is a declared band, no change (owner, 2026-11-01). CF-077
(methylpyrazine cross-pair factor) is merged (PR #37).

**Waiting for evidence:** E15, E17 - E20, and E2's hexose part (E2 and E16 enter W2 as rivals inside M1b, 2026-11-04). Each enters W2 when a labelling, fed-intermediate,
trapping or within-study result discriminates its rivals, for instance from a W6 author reply.

### W3. Engine: one network, built to be fast

| task | what | state | effort |
|---|---|---|---|
| A3 | One explicit configuration object instead of global switches | **done** (PRs #34, #36) | -- |
| **A2 + H3** | **One composable network** (each chemistry module adds its species, reactions, rate laws and pH terms; one integrator) **written so it compiles with exact gradients** in **JAX + diffrax** (chosen 2026-11-05 by the sulfur-lane bake-off, [jax_bakeoff_sulfur_lane.md](../../results/validation/jax_bakeoff_sulfur_lane.md): rates to 1e-12 and trajectories to ~5e-6 dex of the engine, four stiff corners pass; the exact Jacobian 2.5x faster than the fits' finite-difference one at equal tolerance, 6.6x at rtol 1e-6; one forward solve 2.5x slower, ~75 % of it the pH root-find; the fits' forward-difference Jacobian is up to 7 % off its own central difference). The lanes become modules of one recipe. Agreement with today's lanes pinned by tests: rates to 1e-12 relative, trajectories to 1e-6 dex (bit-identical cannot hold across integrators). Forward mode for fits of ~20 parameters; fits at rtol 1e-6 once a check shows panel residuals move < 1e-3 dex; one focused step on the pH root-find. CasADi only if H4 moves off HMC. Then the chemistry of M4 moves onto it. **Network infrastructure (2026-11-05, PR #74, `a2-network-infra`, open; `results/validation/jax_network_infra.md` on that branch):** forced inputs, the fit driver and the density bridge reproduce the engine (fit cost to 3e-8 relative), but at 5 coordinates the JAX fit is ~1.2x slower than the engine's polished fit (plus a one-time compile), and batch draws give no CPU gain at equal tolerance (~2.4x only at rtol 1e-6). Porting sulfur early would save ~0 - 5 serial h of M1b / M1c density compute against 9 - 18 h of porting, so sulfur was to stay last; the lever is FD-dominated fits with many coordinates (P ~0.5 that it saves more than 15 serial h). **Sulfur port started 2026-11-05 (owner),** in parallel with M1b(i) / M1c: structure derived from the NumPy lane where it is data, custom terms hand-ported and listed, at the shipped default (off-by-default rivals stay NumPy-only). A parity test in CI (panel pots plus random in-band parameter sets, and identical reaction and parameter sets) replaces the one-writer rule for `sulfur.py`: any later edit to the NumPy lane updates the JAX lane in the same PR | contract in [A2_MODULE_CONTRACT.md](A2_MODULE_CONTRACT.md); its requirements include each case's "needs from the network" (`docs/reference/cases/` M1 - M6: cysteine-hexose Amadori and 1-DG species, amine classes as species with declared competition, an oil phase with phase-tagged observables, buffer and water-activity hooks, aminoketone / pyrrole / adduct pools, isotopologue pools, the forced-input driver); H3 was planned after A2, now merged into it so the engine is rewritten once | 75 - 115 h (bake-off estimate: module contract 20 - 30, lanes and switches 30 - 45, thermal and pH-stat modes 6 - 10, fit driver 10 - 15, pH step 8 - 15) |

Why together: A2 rewrites the lanes anyway. Writing them once in a form with exact gradients makes every later refit
faster (the bake-off measured 2.5 - 6.6x on the Jacobian, which is ~90 % of fit time per the W1.3 profile) and makes the posterior (W4) affordable. Starts when M3's Strecker layer needs the shared trunk: M1 / M2 still reshape the sulfur and trunk lanes, and a port started earlier would chase them.

### W4. Uncertainty: the hierarchical Bayesian fit

This was "inference track H". It gives intervals that hold for a NEW laboratory. It does not fix a wrong structure
(P ~0.15 that it alone makes the thiol lane beat the nulls), and it can hide one inside laboratory offsets, so it
runs after the chemistry steps it depends on, with tests that catch that.

| task | what | depends on | effort |
|---|---|---|---|
| H0 | between-lab error out of sample (listed under W1.2) | -- | -- |
| **H1** | **Specify the hierarchical model,** pre-registered: laboratory offsets on observation levels and on reported conditions (heat-up lag, pH drift, buffer); on a rate constant only where >= 3 labs inform it. Masking tests: (a) lab offsets regressed on recorded pH, T, buffer, method -- a trend names a missing mechanism; (b) offsets must not change within-lab signature verdicts; (c) an offset > 0.3 dex on a physical constant is a structural error | H0, M2 | 4 - 6 h |
| **H2** | **Tool bake-off.** The engine half was pulled forward and done on the sulfur lane (2026-11-05, see A2 + H3: JAX + diffrax chosen). Left for the Strecker layer: NumPyro; simulation-based calibration; a deliberately wrong structure must be caught; GPU benchmarked here before any Colab run | M3 | 5 - 8 h |
| H3 | the fast engine: merged into W3 (A2 + H3) | -- | -- |
| **H4** | **Full posterior** on each lane that passes its structure tests; compared with the cheap alternative (Laplace + H0's offset) by held-out log density | W3, H1, JR | 1 - 2 days + 6 - 24 h per run |

Rules kept from track H: a separate `maillard-jax` environment (the core lock unchanged); Colab returns posterior draws
only, and every artifact is made locally from them; nothing private goes to Colab.

### W5. The tool for formulators: recipes, cooking, and the deliverable

Standalone modules to the A2 contract. They move no shipped constant. Each says in its output that it is not
validated on a plant-based product.

| task | what | sources | effort |
|---|---|---|---|
| W5.1 | **Ingredients to precursors.** Yeast extract to cysteine, glutathione, IMP -> ribose (nucleotide breakdown kinetics) and thiamine; protein isolates to their amine and cysteine content. Removes the scorer's "precursor unmapped" refusals | 14 yeast-extract papers selected 2026-11-02 (composition: Münch 1997, Kuroda 1997, Jacob 2019, Tomé 2021; nucleotide kinetics: Matoba 1988, Shaoul 1987, Davídek 1972, Kuchiba 1990, Kavitha 2007; routes: Peer 1968, Madruga 1998, Alim 2019 labelling; checks: Münch 1998, Alim 2018) | 2 - 3 days |
| W5.2 | **Crust.** A heat and moisture history for a pan-fried or grilled patty, feeding the existing thermal-history and water-activity terms | 8 crust papers selected 2026-11-02 (Parker 2012, Balagiannis 2019, Chan 2020, Purlis 2009, Bassama 2011, Sprague 2011, Oroszvári 2005, Kondjoyan 2016). No plant-based patty profile exists; it is checked on meat and potato data | 3 - 5 days |
| W5.3 | **Heme and iron catalysis** of thiol and lipid oxidation (leghemoglobin products). Builds on M1's oxidation chemistry | Impossible Foods patents already on file | after M1 |
| **W5.4** | **The lever table: the deliverable.** Per odorant and per lever (cysteine or yeast extract, sugar type, thiamine, pH, fat, temperature and time, cooking method): direction, size with interval, and whether the held-out benchmark supports it | JR, H4 (Laplace intervals as the fallback) | 2 - 3 days |

### W6. Outside data, without a laboratory (owner)

- Send the 5 drafted author data requests (`data/articles/OUTREACH_DRAFTS.md`, laptop only).
- Add requests for existing data that would test a W2 step (ENGINE_VS_LITERATURE sec. 3.1, the "Or" alternatives):
  Cerny's raw 145 C isotopologue data (M2); DFA's raw Table 4 and Quest's recovery (M3, CF-073); Zhou 2024's
  supplementary ion data (M6); the BTBU / CAAS MFT yields per glucose consumed (M4).
- For M1b: whether any group holds in-situ trapped deoxypentosones across pH (none is on disk; ENGINE_VS_LITERATURE
  sec. 3.1, E2). It would separate the split rival from the removal rival.
- Literature keying only where a W2 step's test lacks held-out rows (E3, E10, E18 are one laboratory short).
- For M1c (added 2026-11-05): a designed cysteine + xylose time course (free cysteine, xylose, TTCA and pH, cooled,
  at two or three pH values) from any group able to run it; M1c's case names the design.

- For M3 (added 2026-11-05, phase 1 undecidable on disk): a fed 1-deoxyglucosone + phenylalanine pot, and a fed
  ARP-Phe pot with the ARP measured under argon and under air and the headspace stated.

**What would decide each open question is recorded in one place:** `docs/reference/decisive_measurements.yml` (PR #79),
rendered as the data wishlist's last section. Every step that ends "undecidable on disk" adds or updates its entry
there (question, rivals, why the data on disk cannot decide, the smallest design, what each outcome means, what it
unlocks), so nothing is lost while W6 waits.

Before each W6 request, run a targeted forward-citation sweep (W1.4), so outreach asks only for data the literature
does not already hold.

**W6 is optional and comes last (owner, 2026-11-05).** No step waits on it; each step's stop rule covers its absence. Without it:
M1c most likely ends "supply not identifiable on disk" (P ~0.6 that no unseen decisive data exist) and the thiol lane
stays relative-only under H0 with the R0 / C1 supply envelope carried as uncertainty; M2 most likely reads "conflict:
errors stack"; absolute thiol levels for a new laboratory stay refused until second-laboratory absolute data arrive
(WUR blind test or a reply); M1b(ii)'s split and removal rivals stay unseparated. The odds that M1c + M2 noticeably
improve the sulfur lane fall from ~0.25 to ~0.15. A reply that arrives later re-opens the step it names.

## 6. Timeline: what runs in parallel

Each row is a span of calendar days, sized from the effort column at about 5 working days a week; no holiday time is
assumed. Columns are workstreams. Cells in one row run in parallel. Re-sequenced 2026-11-04: M4 needs the one network,
so it now follows A2 + H3 instead of running beside it, and the M1 / M1b / M2 block takes its full ~8 working days.
Re-sequenced again 2026-11-05, after M1's verdict, the supply screen and an independent premortem: M1.1 (measurement
model) first; M1b(i) and M1c next, neither gating the other; M2 only in the joint refit with M1c's structure; M3's
first phase (prereg EX-051, keying Hofmann 2000b; Strecker files only) runs in M1c's compute waits; A2 + H3 stays
after M3 (the lanes are still being reshaped, more so now), its module contract is drafted early, and it ports the
sulfur module last. The block grows by ~40 - 70 h; with M3 phase 1 in the waits, A2 + H3 slips by ~0 - 15 h only.
Re-dated 2026-11-05 (evening): M3 phase 1 closed early (undecidable on disk) and every JAX lane but sulfur is already
ported, so the M4 - M6 cases (no engine runs) and the network join move into November, M4 follows in late November,
and the end moves from 02-05 to ~01-23. No holidays are assumed; treat 02-05 as the buffer.

| dates | W1 yardstick | W2 chemistry | W3 engine | W4 uncertainty | W5 product | W6 owner |
|---|---|---|---|---|---|---|
| 11-02 - 11-08 | benchmark scorer (done); H0 (done, applied); profile (done) | M1 run (EX-061); supply screen; M1.1; Cao input audit; M1c case items 1 - 5 frozen; M1b and M2 preregs drafted; M3 phase 1 built, power-checked and closed (undecidable on disk, 11-05); M0 inventory started | JAX bake-off; contract v2.2; lipid, trunk + acrylamide, trunk lane and network infrastructure ported (PRs #58, #65, #71, #74); sulfur port started (11-05) | -- | -- | -- (W6 last) |
| 11-09 - 11-20 | -- | M1b(i) (timeboxed); M1c screen, build and fit; M0's checks (every leading candidate) in the compute waits; the cases for M4, M5, M6 written, ordered and shaped by M0's groups (sec. 7 items 1 - 5; no engine runs) | sulfur port lands (parity test in CI); the lanes joined into one network with exact gradients | H2's rest (NumPyro, simulation-based calibration, a wrong structure caught) on the network | -- | -- |
| 11-21 - 12-05 | -- | joint M1c + M2 (+ M1b(ii)) refit, with M1's C2 (a sink that free cysteine slows) in its rival set by amendment before its freeze; M5 (moved up 2026-11-06: M4's build waits for designed data); a case for any shared structure M0 nominates | A2b (declared-pH modules beside the inline lane, one titratable table and one pH module across lanes) and M4's structure, at parity and off by default | H1 spec (after M2) | W5.1 yeast extract (standalone) | -- |
| 12-06 - 12-19 | -- | M6; M4 frozen as development grade (its comparison waits for designed data) | -- | -- | W5.2 crust (standalone); W5.1 / 5.2 joined to the network | -- |
| 12-20 - 01-09 | -- | joint refit | -- | H4 posterior (sulfur lane first) | W5.3 heme | -- |
| 01-10 - 01-23 | benchmark re-scored (coverage first; descriptive, see W1) | -- | -- | H4 (the other lanes) | W5.4 lever table | -- |

Hard dependencies: M1.0 before M1's combinations; M1.1 before any thiol refit; M1c's structure, M2 and M1b(ii)'s
chosen rival in one refit (same rows; errors that stack are judged together, CF-107), never M2 at instant supply;
M1b(i) does not gate M1c (M1c runs at the measured pH trajectory where one exists); M0 before the M4 - M6 cases are
written; W3 before M4 and H4; M2 and H0 before H1; the joint refit before H4; H4 (or its fallback) before the lever
table. (M3 no longer gates H2 or W3: phase 1 closed undecidable on 11-05.) While A2 + H3 is built, a chemistry step
that edits a lane being ported updates its JAX module in the same PR, enforced by the lane's parity test (owner,
2026-11-05, for sulfur); otherwise the chemistry line writes the next three cases, which need no engine runs.

## 7. How a change ships

**Principle.** A change ships when the evidence favours it over the current model and over the other explanations,
judged on the data that measure that mechanism most directly. The rule is written so that the decision does not
depend on who makes it or on what they hope to find. It replaces the 2026-11-01 admission gate (>= 2 held-out
laboratories before any work) and the owner-approval step.

**The case.** One page in the PR. Items 1 - 5 are written BEFORE any result of the change is computed.

1. **The claim.** Which mechanism, and which steps of the network change.
2. **All the evidence, by kind.** The kinds rank by how directly they show a mechanism, and the rank is fixed:
   isotope labelling; fed intermediates; trapping; measured intermediates and lifetimes; mass balance; one factor varied
   inside one study; fitted time courses. The case lists every source on disk that bears on the claim, including those
   against it. A source left out needs a stated reason. Before items 1 - 5 freeze, the case also records a citation
   sweep of its decisive sources (forward and backward, W1.4), keyed to the case's own decisive criterion rather than
   its topic, with every unseen citing work listed or excluded with a reason (titles / abstracts only, so blinding
   holds). A source found after the freeze enters only by amendment, declared before its data are read.
   Data that an M0 check or any earlier diagnostic used to NOMINATE this claim are listed as spent: they count as
   development evidence and cannot be its decisive data (item 4). Spent or seen data may still FALSIFY a rival whose
   prediction is constant-free (structural: nothing was tuned on them), but never CREDIT one (2026-11-06, from the
   M4 premortem).
3. **The rivals.** The current model, the proposed change, and at least one other explanation of the same evidence.
4. **The decisive data.** The measurements on which the rivals predict different things, chosen from those closest to
   the mechanism: route shares from labelling, branching ratios, ratios that cancel scale, the sign of a within-study
   contrast. End-point levels that other errors also move are not decisive. If the rivals predict the same thing on
   every available measurement, the data cannot decide, and the case says so.
5. **Tolerances from the data** (revised 2026-11-06 after the judge calibration, PR #99, an independent review and a
   confirmation run). Each comparison uses the noise of the measurement at the level at which the compared value
   varies: the SD of independent runs where printed, integrated over its sampling uncertainty (n - 1 df; an SD whose replicate level is unstated counts as
   repeatability, a floor; an unstated n is taken as 2); else the
   same laboratory class's measured run-to-run reproducibility of the same observable on a repeated design (a run
   offset and a per-point term where the repeat shows both); an instrument's repeatability or a reading error is a
   floor, never the whole noise of a reaction-system observable. A model-discrepancy term is a variance component (per
   point, or shared within a run where its estimate shows that structure), estimated on the same observable and
   laboratory class from rows that are neither decisive nor fitted with a shared constant, as the residual variance
   pooled over the rivals' fits, integrated over its estimation uncertainty, the
   same for every rival, and redrawn in the calibration. Where no run-to-run noise (a printed run SD or a same-observable reproducibility) or no such discrepancy
   estimate exists, zero is not assumed: only values in which a
   run-level term cancels (within-run ratios, shapes, signs of within-study contrasts) are decisive for credit, and
   absolute levels are falsification-only or reported; a falsification on them stands only under the largest noise any
   listed source supports, at least the lane's within-laboratory term (S1 - S4 preamble). A within-run value's per-point
   term is the printed within-run replicate SD or a same-observable within-run repeat; with only an instrument floor, it
   is a variance component integrated over a prior declared with items 1 - 5 and bounded below by the floor, and the
   false-credit check below decides whether it can credit. Across laboratories: the measured between-laboratory spread of
   the same observable and value basis (H0's tau for thiol-lane absolute levels; else a measured spread of that same observable listed in the S1 record); with none, a
   cross-laboratory comparison is not decisive. No term estimated on another observable is transferred, and no bound is
   used as a standard deviation. No round-number thresholds. Item 5 governs the comparisons of items 6, 9 and M1.1; the lane
   terms of S2 (b), S3 and S4 are those rules' own declared terms.

   **The decision statistic (M1.1).** Rivals a and b (either may be the incumbent) are compared by ONE pre-declared
   statistic T: the pooled log predictive density difference over the decisive units, each density conditioned on the
   spent rows (p(y_dec | y_spent); the spent rows' own likelihood ratio never enters T). a is credited only if T exceeds
   the 95th percentile of T under b as truth AND a is adequate on every decisive unit (posterior-predictive p >= 0.05
   after Holm across the decisive units, each p under p(y_dec,u | y_spent)); b by the mirror rule; if both meet their
   conditions, neither is credited; otherwise the pair is undecided and both are carried (S3). A rival falsified under
   S3 is neither credited nor used as a null. A rival left with no unfalsified rival to serve as its null is carried, not
   credited (a falsification, above all on spent data, never credits the survivor), and S3's refusal clause applies to
   outputs that depended on a falsified incumbent. Each null is a parametric bootstrap from the null rival's spent-conditioned
   posterior, its truths drawn outside the scoring ensemble, redrawing the decisive data, the conditioning rows, the
   forced inputs, the replicate SDs, the discrepancy component and the form index, with >= 1000 replicates and the
   threshold's Monte Carlo error reported; the power check also reports false credit under a generator that adds, beyond the
   noise the scoring model uses, the same-observable reproducibility term (else the lane's within-laboratory term), and if the lower end of its 95 % Wilson
   interval exceeds 0.05 the case cannot credit. Forms: an uncertain experimental fact is a discrete nuisance coordinate; an
   unsourced modelling choice keeps required agreement; forms that only hedged omega are dropped; each form's kind, and a nuisance coordinate's weights (equal
   unless sourced), are declared with items 1 - 5. Unit contributions and
   leave-one-unit-out verdicts are printed; a credit rests on a unit when leaving that unit out drops T below its recalibrated threshold; such a credit
   stands only if that unit's own T_u exceeds the 95th percentile of T_u under the null rival.
   **A verdict needs resolved densities and a complete null:** a density is resolved when the shared estimator accepts it
   (`src/kinetic_core/predictive_densities.py`: k-hat < 0.7, ESS >= max(100, its precision requirement)); a null that does
   not redraw every listed input is not the rule's null; a bootstrap replicate whose density the estimator does not
   accept is re-estimated once where fresh draws are supplied (the verdict records whether re-estimation was attempted
   and how many replicates it retried); any still unresolved is never read as estimated but bounded: in testing a's
   credit every such replicate's T, T_u and T_-u is set to +infinity in both nulls, in testing b's to -infinity, and a
   credit stands only if it survives that bound; every crossing above is strict: the null's weight at or beyond T must
   be below 0.05 (a null whose weighted unresolved fraction reaches 0.05 therefore cannot credit); each null's
   unresolved fraction is reported; the judge calibration's 0.70 - 0.90 power used the prototype's adequacy,
   conditioned on the spent AND decisive rows, so the spent-only adequacy above is uncalibrated except in synthetic
   worlds (added 2026-11-07, PR #105: adversarially reviewed and independently confirmed with edits, 2026-11-07;
   applies to a case only through that case's power check under this text); a pair lacking either is 'not computable', reported so and
   carried undecided (S3); an unresolved density is re-estimated (larger draw sets, or the reusable SMC
   estimator), never read as a tie; every power check reports which cells are computable and the cost of making the
   rest computable, before the run. This does not reopen M1b part (i) (its amendment A1 is frozen); A1's would-be kills enter S3 as falsifications only once their sec. 3.8 confirmation fits
   pass S3's test.

Then:

6. **The comparison.** Every rival is scored the same way, on the same data, under the same measurement model.
   Constants come from measurements or same-class analogues where they exist. Only converged fits count. The result
   is a comparison between rivals (likelihood on the decisive data), not a pass mark (item 5's adequacy gate is a
   condition on credit).
7. **The side effects.** The full held-out benchmark, before and after, every row. A row that gets worse is expected
   when the old model relied on a compensating error. The case names the suspected compensating error and registers it
   (a CF entry). A loss does not block shipping. A large loss with no candidate explanation is a reason to look
   further before shipping.
8. **The adversarial review.** An independent reviewer whose task is to build the strongest case AGAINST the change.
   The reviewer sees the evidence and the data, not the author's conclusion. Every objection is answered in the record,
   or the change waits.
9. **The decision.** The change ships if it is credited under M1.1 against every rival not falsified under S3 (at
   least one) and no objection is left
   unanswered. The decision and its basis go into the experiment register (an EX entry). A shipped change is reversed
   when later evidence favours another rival.

**Symmetry: the incumbent is not privileged (owner, 2026-11-06; revised after an independent review the same day).**
Items 1 - 9 make a change beat the current model, but the current model never had to pass the same test, so when the
data could not decide, its errors persisted by default. Rules S1 - S4 close that gap. Order: S1's record must exist
before S2 or S4 is applied to any constant. Every S-rule change is a case with its own adversarial review and
confirmation (item 8), and never upgrades a prior without an artifact in `results/validation/` (AGENTS.md). The
owner-first clauses that still apply are: amendments to H0's tau, E4 (the pyrazine pH band), W6 (outreach), numbers
already sent outside, and deleting scripts; otherwise S-applications ship under item 9's authority ("the working session
applies this rule and ships"). The lane's within-lab / between-lab terms below are, for the thiol lane, H0's omega
(0.419 dex) / tau_ABS (0.964 dex); a lane without H0 uses the pooled between-laboratory spread of the reaction class in
the S1 record, and its absolute outputs are refused as in the thiol lane until that lane has its own H0.

- **S1. Evidence tiers, per coordinate** (log k_ref and Ea separately), in one hand-kept record with a citation per
  row (paper key, ledger row or fit report): measured (a direct measurement of the step: every channel the engine's
  constant lumps present in the measured system, at the engine's rate law and concentration regime, under comparable
  conditions); analogue (a same-class measurement transferred through a measured law, with its band); bounded
  (one- or two-sided, stated; a measurement of a SUBSET of a lumped step's channels is a one-sided bound, not a
  value: a lower bound on the lumped constant when the channels add); fitted-identified (a converged fit whose marginal profile is pinned with every other unidentified coordinate
  profiled over its band; any band-edge, flat-profile or band-spanning-CI criterion overrides a reported stderr;
  re-derived for the coordinates whose rows a structural change moves); fitted-unidentified; structural (a declared structural assumption);
  physical constant; conflicting measurements; withdrawn. The record also lists the incumbent's structural choices
  that never had a case. The existing vocabularies (`provenance_tier`, the sulfur evidence `grade`) are mapped onto it
  or retired. Upgrades need an artifact and a second reader; downgrades need a citation, and a downgrade that moves a
  shipped output ships as an S4 change, review included. A CI gate fails on a shipped coordinate without a tier; the
  record is complete (the gate green) before any S2 / S4 application, and the record file and the gate are named in the
  PR that creates them. Each never-cased structural choice the record lists either gets a case scheduled in sec. 6 or a
  named rival carried under S3.
- **S2. A direct measurement replaces a contradicted fitted-identified constant, under strict conditions** (a
  fitted-unidentified coordinate has no point to contradict: it is carried under S4). It applies only when
  (a) the measurement is "of the step" in S1's sense (all lumped channels; the engine's rate-law order and
  concentration range at the outputs it serves; otherwise the measurement defines a rival structure, items 1 - 9, or a
  band, S4 — a rival when it implies a different rate law or species, a band when it only constrains the value);
  (b) the 95 % interval of (measurement minus fitted value) excludes zero, with the reported noise, the within-lab
  between-paper spread, the between-lab term and the fitted coordinate's own profile interval combined in quadrature;
  the between-lab term is the measured spread of that step or its reaction class where >= 2 laboratories measured it,
  else the lane's between-lab term (above); no other source; a transfer
  across T, pH, buffer, O2, concentration or matrix adds a term only through a measured law; if a condition the step
  depends on is unstated, the band is undefined and the measurement enters under S4, not S2; (c) the case lists EVERY
  measurement of the step, agreeing ones included, pooled with the between-lab term; one laboratory's contradiction
  among agreeing laboratories does not trigger S2; an audit pass covers every coordinate tiered fitted-* in the S1
  record for the lane, each contradiction's z is converted to a two-sided p, all are judged together (Holm, as in M1.1)
  with the number expected by chance reported, and no coordinate enters S2 outside a full-lane pass; (d) the measurement was read from the source by two
  readers, passes OpenAlex `verify`, and is neither a quote of a quote nor an inferred bound (patents count as bounds
  only). S2 is an explicit exception to item 2's spent-data rule (a direct measurement of the parameter was not tuned
  to it). Side effects: (e) if any held-out row is lost beyond the lane's within-lab term, the named compensating error is a
  specific step, and a forward run shows that moving it
  within its own tier band restores the lost rows' sign and rough size; if no such step exists, S2 waits; (f) rows
  lost beyond the lane's within-lab term that are decisive-type data (within-study ratios, shapes, label shares, multi-laboratory) DO
  arbitrate and block; (g) every constant co-fitted with the replaced one drops to fitted-unidentified until
  re-derived; (h) replacements in one lane ship in one batch case with one before / after. If S2 waits under (e) or is
  blocked under (f), the coordinate's tier becomes "conflicting", a CF entry is opened, and it is carried under S4 as
  a band spanning the fitted value and the measurement's interval: it never ships at the fitted point alone.
  Measurements of one step are "conflicting" when they disagree beyond (b)'s band among themselves (e.g. Kanzler 2021
  vs 2026), else they are pooled under (c); a conflict band spans the measurements' own intervals, with a CF entry.
- **S3. Undecided rivals are carried, not resolved by default.** A rival is carried if it was named in a case before
  results and has a mechanism with at least one evidence kind from item 2's list above "fitted time courses", or uses
  only measured / analogue-tier constants; it is falsified when its converged best fit is inadequate
  (posterior-predictive p < 0.05 after Holm across the datasets tested, M1.1) on any data, spent included. The envelope is the union
  of the carried rivals' own envelopes, unweighted ("never credit"). Tolerance: the lane's within-laboratory term for
  relative outputs, its between-laboratory term for absolute ones. Rivals that disagree in sign make the lever's direction refused; rivals that agree
  in sign show its size with the spread. If the incumbent itself is falsified with no rival credited (e.g. C2 + C3
  recombination as the hexose MFT route, M4), the outputs that depend on it (those whose sensitivity to the
  falsified structure's switch exceeds the lane's within-lab term) are refused, never shown at the incumbent's point.
  The refusal is lifted only by an S3 change of its own (items 6 and 9 replaced as below) that ships the unweighted
  union of the remaining carried rivals' envelopes. A rival enters the union only if its converged best fit has passed
  this rule's test on the data that falsified the incumbent; the union credits no rival. Until then the refusal stands.
  *(Lift clause added 2026-11-07. Its wording follows M1b part (i) prereg sec. 21.8.6 item 2 (the B1 u B13 envelope
  "only through its own reviewed case") and that item's confirmation edit, sec. 21.9 (only after converged S3 fits of
  B1 and B13); pending independent confirmation.)* The reviewer (item 8) checks both missing and padded rivals. Carried rivals must be runnable for
  prediction: until W3 runs them in JAX, they run in the NumPy engine behind their switches, and the cost is stated.
- **S4. Unidentified coordinates are bands, and fits may move them only on decisive-type evidence.** A band comes from
  measured bounds, a same-class analogue, or the reaction-class distribution of measured constants (as B55's
  `ea_class_prior`); a search bound is never a band; with no class, the output is refused rather than computed over an
  invented band. Bands are sampled jointly within the region the decisive-type data allow (common random numbers
  across the arms of a contrast); held-out coverage is reported before and after (over-coverage flags an inflated
  band; an over-coverage flag sends the band back to its source for review). A fit (the joint refit, M2, H4) may move a fitted-unidentified coordinate only if every coordinate whose sensitivity
  on those rows is not negligible next to the lane's within-lab term is free in the same fit (CF-107's logic), and
  its contraction counts only if the profile interval from a refit on the decisive-type rows alone (item 4) lies
  inside the band; M2 v2's held furfural barriers are checked against this before M2 freezes; every fit reports per-coordinate
  contraction, and a coordinate contracted only by level rows stays a band.

For S3 / S4 changes, items 6 and 9 are replaced by: band sources verified; before / after listed; held-out coverage
reported; review. These rules move shipped outputs: intervals widen and some answers become refusals (on current
shares, S4 will likely widen or refuse most absolute thiol outputs and may reach relative ones). That is the intended,
honest consequence; sec. 8's odds and W5.4 are re-estimated in the PR that first applies them.

**Why this resists persuading ourselves.** The ranking of evidence is fixed in advance. The rivals and the decisive
data are named before results exist. Tolerances come from the measurements. The decision is a comparison, not a
threshold someone picked. Evidence against must be listed. The reviewer's job is to disagree.

The working session applies this rule and ships. The owner sees a summary of what shipped and why.

Working rules that stay: run through `./scripts/docker_maillard.sh run`; nothing writes into `data/` at runtime; one
writer per shared file (registers, ledger, atlas); stage files by explicit path; check new CF / EX ids against
`origin/main` and open branches; a "held-out" claim is checked against the split files. The full list of lessons is
`tasks/lessons.md`.

## 8. Checkpoints and odds

| date | done means | P |
|---|---|---|
| 11-07 | benchmark scored; H0 scored and applied; profile done (re-dated 2026-11-05: done but for H0's application). M3 phase 1 reached its verdict early, on 11-05: undecidable on disk, moved to W6 | ~0.8 |
| 11-20 | M0's checks run on every leading candidate, the misfits they leave grouped by shared structure; the M4 - M6 cases frozen (items 1 - 5); the one network runs with exact gradients, sulfur included (rates to 1e-12, trajectories to 1e-6 dex of today's lanes) | ~0.5 |
| 12-05 | M1c, M1b(i) and M2 have verdicts (or M1c's stop rule fired); H1 spec | ~0.45 |
| 12-12 | M4 frozen as development grade (decider W6), its designed experiment sized; W3 carries M4's structure at parity; M5 has a verdict | ~0.4 |
| 12-19 | M6 done | ~0.35 |
| 01-09 | the joint refit done | ~0.3 |
| 01-23 | a posterior on the sulfur lane; the first lever table (Laplace intervals if H4 is not done); the held-out benchmark re-scored, coverage reported first (was 02-05; moved earlier 2026-11-05 evening) | ~0.3 |

Other odds (lowered 2026-11-04: six earlier sink structures were refused, B11 - B27; lowered again 2026-11-05: no
supply structure on disk passes the screen, and P ~0.6 that no unseen decisive supply data exist): M1c + M2 noticeably
improve the sulfur lane out of sample, ~0.25. The sulfur lane beats both null models after M1 + M2, ~0.3. Ranked relative predictions with honest intervals for most levers by the end of January, ~0.35. Absolute levels
within 3x across laboratories for most rows, ~0.15 (not a goal).

Re-estimated 2026-11-07 with the first S3 application, as sec. 7 asks (the outputs that depend on M1b part (i)'s
falsified pH balance are refused, per pot: `results/validation/ph_balance_s3_sensitivity.md`; on balance-solved pots,
33 % of the decided concentration levels move beyond tau, 48 % of the pH levels beyond the pH-level term and 46 % of
the pH-lever ratios beyond omega or in sign, between B0 and the carried rivals' screen-passing extremes; prereg sec.
21.10). The within-study contrasts that remain lose to "no change" (sec. 2). Until an S3 change ships the union of the carried rivals' envelopes (S3's lift
clause: B1 and B13 need converged S3 fits on CF-084's controls and the CO2 ledger audit, CF-134), W5.4's pH, buffer and
dose rows on the sulfur lane are refused wherever the balance moves them; recipes with a held (clamped) pH are
answered. Ranked relative predictions with honest intervals for most levers by the end of January: ~0.35 -> ~0.3. The
01-23 checkpoint stays ~0.3 (its lever table already waited on the envelope).

## 9. Not in the plan

- New laboratory experiments and the designed-experiment brief (paused by the owner).
- Broad literature searches and the consolidation wave V. Reading is limited to the papers named above and W6.
- New engine switches without a case (sec. 7).
- The E6 + E11 refit queue on its own branch as a separate exercise: it becomes M2.

## 10. Glossary of codes

| code | meaning | where |
|---|---|---|
| W1 - W6 | the workstreams of this plan | sec. 5 |
| M1 - M6 (and M1b), JR | the chemistry steps and the joint refit | sec. 5, W2 |
| E1 - E20 | ranked engine mechanism errors | [ENGINE_VS_LITERATURE.md](../reference/pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in secs. 5 and 7 |
| CF-nnn | a registered conflict between sources or between data and model | `docs/reference/conflicts_register.yml` |
| EX-nnn | a registered experiment (fit, test, review) and its verdict | `docs/reference/experiment_register.yml` |
| A2 | the one composable network (W3) | [A2_MODULE_CONTRACT.md](A2_MODULE_CONTRACT.md) |
| A3 | the explicit engine configuration (done) | `src/kinetic_core/config.py` |
| H0 - H4 | the uncertainty stages (W1.2, W4); H3 is merged into A2 | sec. 5 |
| B-numbers (B1 - B64) | earlier fit waves of the shipped calibration (B60 is shipped) | `results/validation/kinetic_core_b*` |
| TEST / DEVELOPMENT | the frozen split of held-out vs development evidence | `docs/reference/contrast_benchmark/` |
| lane | one of the four separate networks today (trunk, sulfur, acrylamide, lipid) | README |
| MFT, FFT | 2-methyl-3-furanthiol, 2-furfurylthiol: the key meaty odorants | [GLOSSARY.md](GLOSSARY.md) |
| norfuraneol (NF), TTCA, DPO, 1-DG | intermediates: a furanone, the cysteine-sugar ring compound, 1-deoxypentosone, 1-deoxyglucosone | [GLOSSARY.md](GLOSSARY.md) |
| dex | a factor of 10 on a log scale (0.3 dex ~ 2x) | -- |
