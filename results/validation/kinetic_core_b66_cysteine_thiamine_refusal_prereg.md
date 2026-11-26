# B66: MFT on thiamine's route is refused where cysteine is also charged (pre-registration)

*2026-11-26. Written before any engine file is edited and before the panel is re-scored. Owner approval:
STEP5_RECORD amendment of 2026-11-26 (second entry). Nothing is fitted; no shipped constant changes.*

## 1. Why

B65 looked for a mechanistic step inside cited bands that makes cysteine raise 2-methyl-3-furanthiol (MFT) on
thiamine's route, and found none ([plug-in screen](kinetic_core_b65_plugin_screen.md)). Its outcome line then kept the
incumbent ("no engine change"). Under the symmetric sec. 7 the incumbent is not privileged, and it is contradicted:

| source | system | MFT with cysteine / without | engine as shipped |
|---|---|---|---|
| Grosch & Zeiler-Hilgart 1992, Table III (doi 10.1021/bk-1992-0490.ch014) | thiamine 60 µM + ribose, 100 °C, pH 5.7, 2 h | x2.7 to x3.7 at 0.17 - 0.65 mM Cys; x0.38 at 5.4 mM | x1.27 at 0.65 mM (screen row 1: 5.44 mM / 0 = 1.91) |
| same, Table II | no ribose | x7.0 | — |
| Zhang et al. 2024, Fig. 1 (doi 10.1016/j.foodres.2024.114149) | thiamine 44.5 mM + xylose, 115 °C, pH 4.9, 1 h | x19 | x1.20 |

Every measured ratio is beyond omega / 2 (0.21 dex) of the engine's, in two published laboratories. (A partner's
private development data agree in direction; they are not in this repository.) The mechanism is undecided (`decisive_measurements.yml`: `cysteine_boost_thiol_yield_concentrated`). So the engine
should not show MFT at the shipped value in such a pot. This is the thiamine-H2S refusal's form
(`src/kinetic_core/thiamine_h2s_refusal.py`, PR #158): a carried rival, and the outputs it moves are refused.

## 2. The rule

- **Applies to:** a sulfur-lane pot whose charge holds thiamine (`THI` > 0) and cysteine (`Cys` > 0). Cystine counts,
  because the engine charges it as cysteine plus oxidant.
- **Rival:** the rate of `r_hmp_mft` (5-hydroxy-3-mercapto-2-pentanone -> MFT) multiplied by **3.7**: Grosch's
  largest ratio with ribose present, the smallest of the published maxima. The rival is never shipped. It is on only
  inside the rule's own second integration.
- **Refused:** each requested output whose log10 spread between the shipped point and the rival exceeds omega / 2,
  the thiamine rule's per-arm threshold, for the same reason (an arm judged alone must be safe in any ratio).
- **Reason text:** names the direction where it is known (higher at cysteine up to about 10 times thiamine), the
  published range (x2.7 to x19), the fall at a large cysteine excess in dilute solution (Grosch), and what lifts it.
- **Switch:** config key `sulfur.cysteine_thiamine_mft_refusal`, shipped True; False where the model's own numbers
  are computed (the rule's own integrations, S3's rival points, fits).
- **Not judged:** a pot without thiamine or without cysteine; every other lane. There the engine is bit-identical.

Why multiply `r_hmp_mft` and not thiamine's cleavage: the probe (sec. 4) shows both move MFT equally, and the cleavage
also moves 3-mercapto-2-pentanone by the same factor, for which no measurement exists.

## 3. Declared choices and their cost

- **Any cysteine triggers it.** The lowest measured dose is 0.17 mM (cysteine : thiamine 2.8 : 1), which already gave
  x2.7. Below that nothing is measured. A trace of cysteine therefore over-refuses. Stated, accepted.
- **No ratio or concentration domain.** At a large cysteine excess Grosch's MFT falls to x0.38 (−0.42 dex), which is
  also beyond omega / 2. The sign of the error changes; the refusal does not.
- **Glutathione is not covered.** One laboratory (Zhang 2024: x15). Recorded here as known and uncovered.
- **The value 3.7 does not matter much.** Any rival above x1.62 crosses the threshold where thiamine's route carries
  the MFT. A larger value would refuse the same outputs.

## 4. The probe that preceded this text (read-only, 2026-11-26)

Shipped sulfur lane, isothermal, pH held; `r_hmp_mft` x3.7:

| pot | thiamine route's share of the MFT flux | MFT | MFT dimer | FFT, H2S, 3-mercapto-2-pentanone |
|---|---|---|---|---|
| Zhang Fig. 1 + Cys (115 °C, pH 4.9) | 0.80 | +0.48 dex | +0.96 | 0.00 / 0.00 / −0.01 |
| Grosch + 0.65 mM Cys (100 °C, pH 5.7) | 0.77 | +0.48 | +0.95 | 0.00 / 0.00 / −0.01 |

The branch of the intermediate that reaches MFT is 1.5 - 1.8 %, so the multiplier is not capped.

## 5. The free-thiol / dimer split

In a pot the rule applies to, MFT and its dimer are both refused, so the split is not shown there. Elsewhere the
evidence on the split is mixed: Grosch 1992 (standards added right after heating) finds the dimer below about 4 % of
MFT, as the engine does; Zhang 2024's control has 36 %; and where the dimer is read by headspace sampling, oxidation
during sampling is not excluded (`decisive_measurements.yml`: `thiol_dimer_in_pot_vs_sampling`). A second
carried rival (an oxidant charge of 10 mM, the screen's row 4) is put to the premortem, and is decided after a count
of the shipped outputs it would remove. It is not part of this change unless that record says so.

## 6. What is expected, written before the run

| outcome | probability |
|---|---|
| Switch off is bit-identical to origin/main on every sulfur-lane panel pot | 0.97 |
| MFT and the dimer are refused in every panel pot holding thiamine and cysteine; no other compound is | 0.80 |
| Between 4 and 20 scored panel rows become refused | 0.70 |
| The frozen-fit drift gate and the JAX parity test are unchanged | 0.90 |
| A directional claim that was scored becomes not evaluable | 0.50 |

## 7. Checks before merge

1. Switch off: bit-identical to origin/main on the sulfur-lane panel (fixture test, as PR #158).
2. S3's and the thiamine rule's own integrations do not pick the rule up (both pass the switch off).
3. Frozen-fit drift gate, JAX parity, the seven gates, unit and scientific tiers.
4. The count of refused panel rows, by study, in the PR.

## 8. What would lift it

A measured cysteine series at fixed thiamine in one medium, with thiamine and cysteine followed over time, with and
without pentose (`cysteine_boost_thiol_yield_concentrated`), or a step that reproduces Grosch's and Zhang's dose
curves inside cited bands.

## Amendments

**2026-11-26, after the independent premortem** ([kinetic_core_b66_premortem.md](kinetic_core_b66_premortem.md)),
before any commit. Sections 1 - 8 stand as written, except that a partner's private numbers were removed from secs. 1,
4 and 5 (they do not belong in a public repository). What changed in the rule, and why:

- **Donors (sec. 2, sec. 3).** The rule applies to thiamine with cysteine, glutathione or TTCA. Sec. 3 was wrong
  that glutathione has one laboratory: Grosch Table II no. 7 against no. 3 is 23 against 5.2 µg, x4.4, beside Zhang's
  x15. TTCA opens to cysteine in the pot. "Cystine counts" was empty: the front door has no cystine alias.
- **Factor (sec. 2, sec. 3).** 7.0, not 3.7: Grosch Table II no. 4 against no. 1, the largest ratio measured where
  thiamine is the pot's only MFT source. "The value does not matter much" was wrong for mixed pots: at 3.7 an output
  is refused only where thiamine's route carries more than about 23 % of the MFT, at 7.0 above about 10 %. Zhang's
  x19 is not used as the factor, because the sugar route also makes MFT in that pot. On the tracked evidence the two
  values refuse the same outputs.
- **The test against the incumbent (sec. 1).** A with / without ratio is judged by omega (0.42 dex), not omega / 2.
  Against omega: x2.7 is inside (0.33 dex from the engine); x3.7 (0.46) and x4.4 are just outside; x7.0 (0.85) and
  Zhang (1.2) are far outside; the fall to x0.38 has the other sign. The clause that applies is S3's "incumbent
  falsified with no rival credited: the outputs sensitive to it are refused" (STATUS_AND_PLAN, S3). The multiplier is
  a sensitivity probe, not a carried rival in S3's sense.
- **No direction is stated in the reason text.** The effect follows the cysteine concentration, not its ratio to
  thiamine (CF-089; Grosch's 12 µM thiamine series peaks at 57 : 1), and in Bolton 1994's pot the engine is already
  1.3 dex above the measured level.
- **"No TEST contrast is read"** means no new one: Grosch's contrasts were already spent.
- **Sec. 5 is closed: no second rule for the dimer split now.** A lane-wide oxidant charge of 10 mM would push an arm
  beyond omega / 2 in 5 of the 17 evaluable sulfur directional claims, for an output that no scored row and no
  evaluable claim reads. The best method on disk (Grosch, standards added after heating) agrees with the engine,
  10 mM is far above oxygen saturation with no source, and the sampling question is already registered.

**2026-11-26, outcome against sec. 6** (sec. 6 is not edited).

| forecast | probability given | result |
|---|---|---|
| Switch off is bit-identical to origin/main on every sulfur-lane panel pot | 0.97 | **Held.** All 25 pots, both rival switches off (fixture test). |
| MFT and the dimer are refused in every panel pot with thiamine and cysteine; no other compound is | 0.80 | **Held, on one pot.** Bolton 1994 is the only such pot; its one target, MFT, is refused (0.79 dex at x7.0). |
| Between 4 and 20 scored panel rows become refused | 0.70 | **Missed.** No scored row changes: Bolton's MFT was already not evaluable (absolute level). The rows refused for an absolute level go from 27 to 26. |
| The frozen-fit drift gate and the JAX parity test are unchanged | 0.90 | **Half.** No frozen fit moved. The two new config keys changed the stored config digests, so the JAX trunk references were rebuilt (numbers unchanged), as the 2026-11-18 lesson warns. |
| A directional claim that was scored becomes not evaluable | 0.50 | **Happened, twice.** WANG12-01 (an agree) and MOD08-PH-01 (a miss): headline 10 of 19 to 9 of 17; 8 of 14 without pH and water activity; 1 of 3 on them; 58 independent claims not evaluable. |

The miss on the third row came from not looking at the panel before forecasting: one pot holds both precursors, and
its row was already refused. The cultivated-tissue surface charges thiamine and cysteine in every pot, so MFT is
refused there outright; no number was shipped there before either (absolute levels were already refused).

**Unrelated finding.** Re-running the chance baseline of the mechanism signatures
(`kinetic_core_mechanism_signatures_chance_r3en`) changed one pass rate from 0.195 to 0.19 (one draw in 200). The same
re-run on unmodified origin/main (1d16bcdb) gives 0.19 as well, so the stored value no longer reproduces there either.
This change did not cause it. The record is re-stamped at the value both runs give.
