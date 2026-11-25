# B65: cysteine as a thiol redox partner in the sulfur lane — pre-registration

*Written 2026-11-26, before any feasibility run. Frozen at the commit that adds this file; changes after that are
dated amendments below, never edits. Engine-free: nothing here changes a shipped constant or structure.*

## 1. Why

**Owner decision (2026-11-26).** The engine pause (STEP5_RECORD) is lifted for two named changes only:

- (A) cysteine reducing thiol disulfides (bis(2-methyl-3-furyl) disulfide, MFTD; bis(2-furfuryl) disulfide, FFTD)
  back to the free thiol by thiol–disulfide exchange;
- (B) cysteine raising the yield of 2-methyl-3-furanthiol (MFT) on thiamine's route.

Both are to be built from the literature, not fitted to a partner's development data.

**What is already on record.**

- **CF-089.** M0 checks CK-42 / CK-42b, `results/validation/m0_checks_batch2.md` sec. 2:
  - Grosch & Zeiler-Hilgart 1992 shows a ×7 cysteine effect on thiamine-MFT that is not hydrogen sulfide release.
  - In that dilute system the effect follows the cysteine concentration (an optimum near 0.34–0.68 mM at two thiamine
    charges), not the Cys : thiamine ratio.
  - A first-order cysteine capture cannot make its high-dose fall.
  - The register's remedy names the form g(total free thiol) on the thiamine → MFT steps. It records that Grosch is
    spent: it nominated the form, so it cannot test it.
- **B17** (`kinetic_core_b17_prereg.md` sec. 6). A generic disulfide release was not shipped. The engine's dimer share
  is oxidant-limited (0.04–0.9 % against a measured 6.5–9.6 %), not release-limited.
- **EX-050** (`docs/history/thiol_rivals_2026-10-31/`). A slow cysteine exchange was rejected by the red team.
- **`CYSTINE_SULFIDE_EXCHANGE_DEFERRED`** (`parameters_sulfur.py`). The engine has no cystine species. Cystine enters
  only as the oxidant charge of Zhang 2024's Fig. 1 arm.
- **The new tension this wave must face.**
  - Zhang 2024 Fig. 2 (thiamine + xylose + methionine, 115 °C, 60 min, pH 4.9): MFT rises monotonically with cysteine
    from 0.5 to 62.5 mg/mL (about 4 to 516 mM) and then plateaus. There is no fall.
  - Grosch's dilute system falls below its cysteine-free value by about 5 mM.
  - A concentration law built on Grosch alone predicts strong suppression in concentrated systems.

Literature verified at source for this wave: Grosch 1992 (Tables II, III, VI, VIII), Bolton 1993/1994, Nogami 1967,
Nogami 1968, Szajewski 1980, Keire 1992, Lees 1993, Rothwarf 1992, and US 6,358,549 B2. The verification notes are
private working files. Every number used below is cited to its primary source.

## 2. Questions

- **Q-A.** Is there a structure in which cysteine keeps MFT free (as opposed to dimerised) that reproduces the
  measured dimer-share signatures with constants inside cited bands?
- **Q-B.** Is there a structure for cysteine's effect on thiamine-MFT that reproduces BOTH the dilute optimum-and-fall
  (Grosch) AND the concentrated monotone rise (Zhang Fig. 2) with constants inside cited or declared bands? Or do the
  two regimes need a medium or temperature term that no measurement supports?

## 3. The reduced module (engine-free, isothermal)

The module carries only what the questions need. It is a standalone ODE in the style of
`docs/history/thiol_rivals_2026-10-31/`.

- **Species.**
  - THI (thiamine), HMP (5-hydroxy-3-mercapto-2-pentanone), MFT, MFTD.
  - Cys; CSSC (cystine), new; MSC (MFT–S–S–Cys mixed disulfide), new, rival A1 only.
  - OX (oxidant equivalents, with a supply term), H2S.
- **Shared steps (all rivals).** Each uses the engine's shipped constant at the pot temperature, with its shipped band.
  - THI → HMP (k_thi_hmp).
  - HMP → MFT (k_hmp_mft).
  - HMP → loss (k_thiol_decay).
  - MFT → loss (k_mft_decay; plus the thiolate loss, k_thiolate_loss, with the engine's thiolate factor).
  - 2 MFT + OX → MFTD (k_dimer_mft, measured barrier 122.2 kJ/mol).
  - Cys → H2S (k_cys_h2s, Zheng & Ho 1994, measured).
- **Oxidant.**
  - OX(0) is declared in [0.01, 10] mM.
  - Supply (oxygen ingress) is declared in [0, 1] mM/min.
  - Both are wide on purpose: B17 found the engine oxidant-limited, and no measurement of oxidant supply in these pots
    exists on disk.
  - Cystine charged as an additive enters as CSSC (rival A1) or as OX (rival A0, the engine's convention).
- **The ribose route** (norfuraneol + H2S → MFT) is included at shipped constants for pots that contain a pentose. In
  Grosch Table II it is a minor source: ribose + Cys gives < 0.1 µg.

## 4. Rivals

### Q-A: the free/dimer split

| rival | structure | new constants and bands |
|---|---|---|
| **A0** (engine) | 2 MFT + OX → MFTD; Cys + OX → fragments (k_cys_ox); no exchange; cystine = OX | none |
| **A1** (fast exchange) | Cys + OX → ½ CSSC (cysteine oxidation goes to cystine, not fragments). MFTD + Cys ⇌ MFT + MSC; MSC + Cys ⇌ MFT + CSSC. Equilibrium-controlled: forward k at least 10³ × the slowest competing step, reverse from K. | K1, K2 (dimensionless) each in [0.1, 10]. Source: thiol–disulfide equilibrium constants for low-molecular thiol pairs of similar pKa are about 0.3–3 (Szajewski & Whitesides 1980, K^SH ≈ 1; Keire 1992 Cys/GSH; Lees 1993, thiophenol vs mercaptoethanol 0.31 at pD 7). The band is widened ×3 each side for an aryl thiol at 115–150 °C. Speed: Nogami 1967/1968, Cys + thiamine alkyl disulfides, apparent Ea 68–80 kJ/mol at fixed pH, about 10⁵–10⁶ L mol⁻¹ min⁻¹ extrapolated to 150 °C. Equilibrium control holds across that band. |
| **A2** (negative control) | as A1 but with K1, K2 in [10³, 10⁴] (exchange drives MFT into disulfides) | must be infeasible on S2; if it is not, the signatures have no power for Q-A |

### Q-B: cysteine on the thiamine route

| rival | structure | new constants and bands |
|---|---|---|
| **B0** (engine) | no cysteine term | none |
| **B1** (CF-089 form, concentration law) | k_hmp_mft × g(c), with c = total free thiol (mM): g = (1 + a·c/(c + K_p)) / (1 + (c/K_d)²) | a in [0, 10]; K_p in [0.03, 3] mM; K_d in [0.3, 30] mM. Shape from CK-42/CK-42b; the bands span Grosch's nominated optimum and fall, widened ×3. T-independent. |
| **B1T** (B1 with temperature) | as B1, with K_d(T) = K_d(100 °C)·exp(E_d/R·(1/373 − 1/T)) | E_d in [−150, +150] kJ/mol (declared; no source). If only B1T is feasible, the result is "needs an unmeasured temperature term", recorded as an open hypothesis, not shipped. |
| **B2** (oxidant scavenging, mechanistic) | HMP + OX → loss (k_hmp_ox) and the shared MFT + OX dimerisation compete with Cys + OX (k_cys_ox, shipped). Cysteine protects by consuming oxidant, so the boost saturates when Cys ≫ oxidant supply. There is no fall term. | k_hmp_ox in [0, 10]·k_dimer_mft (declared); OX band as sec. 3 |
| **B3** (B2 + inhibition of thiamine breakdown) | as B2, with k_thi_hmp / (1 + c/K_i) | K_i in [1, 1000] mM. Arnold 1971: "no thiamine lost" with 8.3 mM Cys at pH 5–7 in boiling water, qualitative; CF-089 explanation 2 gives the detection caveat. |
| **B4** (negative control) | as B2 but cysteine adds oxidant (sign reversed) | must be infeasible on S1/S3 |

The rivals are not exclusive. Every feasible one is reported. Under the symmetric sec. 7, undecided rivals are
carried as uncertainty; none is chosen on score.

## 5. Signatures (mechanism evidence; scale-free where possible)

Tolerances are written now. "Dex" is log10. Read error for Zhang figures is ±3–5 % of full scale (dossier). Grosch SDs
come from Table II.

| id | source (status) | conditions | signature | tolerance | discriminates |
|---|---|---|---|---|---|
| S1 | Zhang 2024 Fig. 1a (in-sample fit row in B2.3) | thiamine 44.5 mM, xylose 99.9 mM, 115 °C, 60 min, pH 4.9; Cys 123.8 / GSH 48.8 / cystine 62.4 mM | MFT ratio to control: Cys 19, GSH 15, cystine 29 | each within ±0.3 dex; ordering cystine ≥ Cys ≥ GSH not required (read error) | B (all) |
| S2 | Zhang 2024 Fig. 1a | as S1 | dimer share by mass, MFTD/(MFT + MFTD): control 0.36, Cys 0.079, GSH 0.063, cystine 0.35 | Cys and GSH < 0.15; cystine > 0.2; control > 0.2 | A0/A1/A2 |
| S3 | Zhang 2024 Fig. 2d (not a fit row) | Met:VB1:Xyl 1:3:3 w/w (absolute thiamine not printed; declared band 10–100 mM), Cys 0.5–75 mg/mL, 115 °C, 60 min, pH 4.9 | MFT rises at least ×9 from 0.5 to 62.5 mg/mL, with no step falling by more than 0.1 dex up to 62.5 mg/mL | as stated | B1 vs B2/B3 |
| S4 | Grosch 1992 Table III (spent; direction only) | thiamine 60 µM, ribose 1.34 mM, 0.2 M pyrophosphate pH 5.7, 100 °C, 2 h | optimum at 0.34–0.68 mM Cys; MFT at 5.44 mM below the Cys-free value | optimum within one printed grid step; fall below Cys-free | B1/B2/B3 |
| S5 | Grosch 1992 Table III, 6 µmol thiamine row (spent) | as S4, thiamine 12 µM | optimum at about 0.68 mM (concentration, not ratio) | within one grid step | ratio vs concentration |
| S6 | Grosch 1992 Table VIII (spent) | as S4 | in-pot MFT loss with 0.65 mM Cys ≤ loss without (0.034 vs 0.044 /min) | direction | B2 (protection) |
| S7 | Bolton 1993/1994 (panel) | bouillon, 34S-Cys, 120 °C, 1 h | at least 80 % of MFT sulfur from thiamine with Cys present | ≥ 0.8 | excludes any rival routing Cys sulfur into thiamine-MFT |
| S8 | US 6,358,549 B2 (new) | FFTD + Cys / GSH / cystine (1:10 w/w), 100 °C, steam-stripped | thiol recovered: Cys 87 % (pH 7.6), GSH 80 %, cystine 2.8 % | direction only (open, stripped system; not an equilibrium yield) | A1 vs A0 |

**Power check.** Before scoring, each signature is run on synthetic data generated by each rival at the tolerance's
noise. A signature that cannot tell any two rivals apart is dropped and listed.

**Feasible.** A rival is feasible if at least one constant set inside its bands meets every applicable signature.
Sampling: Latin hypercube, at least 20,000 draws per rival, then a local polish of the best 20. The result is reported
as the fraction of draws passing each signature, and as the joint pass.

## 6. Decision rule (written now)

- A rival infeasible on any signature is rejected for this module and recorded.
- **Q-A.**
  - If A1 is feasible and A0 is not: implement A1 in the engine (cystine species; exchange as a fast equilibrium with K
    bands; oxidant supply as a declared band). Re-run the Zhang Fig. 1 fit row this forces (CYSTINE_SULFIDE_EXCHANGE
    deferral note). Judge it by the acceptance rule: evidence audit, then leave-one-paper-out adequacy.
  - If both are feasible: no choice on score. A1 is carried as a rival.
  - If neither is feasible: the oxidant supply is the binding unknown, and it goes to the experiment ranking.
- **Q-B.**
  - If a T-independent rival (B1, B2 or B3) is feasible on S1, S3, S4 and S5: implement it with bands, judged as above.
  - If only B1T is feasible: record "a temperature or medium term is required and unmeasured". The cysteine effect is
    carried as a declared band on the thiamine route, a widening of the interval, not a point. The decisive
    measurement is a cysteine series at two temperatures in one medium (added to `docs/reference/decisive_measurements.yml`).
  - If none is feasible: the structure list is incomplete. Record it and stop.
- Expected outcomes, stated now:

  | outcome | probability |
  |---|---|
  | A1 feasible, A0 infeasible | 0.55 |
  | both A feasible | 0.25 |
  | neither A feasible | 0.20 |
  | some T-independent B feasible | 0.30 |
  | only B1T feasible | 0.35 |
  | no B feasible | 0.35 |

## 7. Development data from a partner (private overlay)

A partner laboratory's unpublished blind-test data (development split) bear on both questions. They are not used to
choose a rival or size a constant here. After the public verdict, each feasible rival is run on those pots, and the
directions are recorded in a private file outside the repository. The public record says only whether the directions
agreed. This keeps the partner's data out of the public repository (lesson 2026-11-24) and out of the design.

## 8. What this wave will not do

- No fit to TEST contrasts.
- No change to any shipped constant.
- No new reaction in the engine before the decision rule says so and the owner approves the implementation.
- No use of Grosch 1992 as evidence for adequacy: it is spent.
- No sizing of the cysteine effect from the partner's data.

## Amendments

**2026-11-26: not frozen; the wave closes at the screen.** An independent premortem
([kinetic_core_b65_premortem.md](kinetic_core_b65_premortem.md)) found the draft unfit to freeze:
- A0 and B2 used `k_cys_ox` "as shipped", which is zero (B11 did not ship), so A0 could not win.
- B1's ceiling on S3 is x1.75 by arithmetic.
- S1 and S3 are confounded by the pentose route.
- S2 may be set during sampling.
- Cystine is charged as Cys + OX in the engine, not OX alone.
- EX-050's rejected exchange was cysteine + xylose (TTCA), not thiol-disulfide.

As the 2026-11-05 lesson requires, a plug-in screen ran before any redesign
([kinetic_core_b65_plugin_screen.md](kinetic_core_b65_plugin_screen.md): 88 variants, 57 pots, shipped values and band
corners only, no fit).

**Outcome, under sec. 6's rule.**
- **Q-A.**
  - A1 (exchange at equilibrium) is refuted: it drives the dimer to ~1e-9 of MFT, against measured shares.
  - A0 with `k_cys_ox` on its band meets S2 only with an oxidant charge near 10 mM.
  - A measurement model (dimer formed in the sampling vial) is not excluded.
  - So the binding unknown is where the dimer forms. It goes to `decisive_measurements.yml`
    (`thiol_dimer_in_pot_vs_sampling`).
- **Q-B.**
  - No single structure inside cited bands reproduces S1 and S3.
  - Only combinations pass, and they need the pentose route about +1.9 dex outside its fitted band.
  - None also reproduces Grosch's dilute optimum and fall.
  - Zhang's cystine arm is unreachable under the engine's cystine convention.
  - The structure list is incomplete. This goes to `decisive_measurements.yml` (`cysteine_boost_thiol_yield_concentrated`).
- **No engine change is proposed.** The lift of 2026-11-26 produced no admissible step. The pause is back in force for
  these two changes until a measurement named in the two entries exists.
- **Pots seen by the screen:** listed in the screen record, sec. 6. They cannot be decisive for these questions later.
