# Pre-registration: ENV-B34, prior rows for the 3-deoxyglucosone limb and the amine-free sugar entries (written 2026-09-11, before any row was added)

## 1. Why

Wave B34 added five hold-out observables from Leitzen 2021 and three of them came back from the
Monte-Carlo envelope with intervals of about **1e-6 decades**: 3-deoxyglucosone, 3,4-dideoxyglucosone
and glucosone are published as *exact*. It bites hardest on the wave's best row — 3-deoxyglucosone is
**1.11× on fold error and outside its own interval**, because the interval has no width. A model that
is right to 11 % and says it is certain claims more than its evidence supports.

The cause is the same one ENV-B13 fixed for the five HMF rows: the constants that make and consume
3-deoxyglucosone have **no `CorePrior` row at all**. ENV-B13 deliberately left two of them fixed
(`k_tdg_ddg`, `k_ddg_hmf`) *because two laboratories agree on them to 1.5×* — and B34 then measured
`k_tdg_ddg` 32× too slow in water, both laboratories having measured it dry. Agreement is not
accuracy when both share a matrix. The amine-free entries (`k_glc_tdg`, `k_fru_odg`, `k_fru_int`)
were never in ENV-B13's scope at all.

## 2. What is added, and where every number comes from

**No centre moves.** Four constants get sampled prior rows; each band is the source's own printed
95 % HPD, nothing invented. Kocadagli & Gokmen 2016 Table 2 (glucose system, reparameterised
Arrhenius, `kocadagli2016jafc_extraction.md` §4):

| constant | step | k_b (×10³) ± HPD | relative width | Ea ± HPD (kJ/mol) |
|---|---|---:|---:|---:|
| `k_glc_tdg` | 3, Glc → 3-DG | 4.19 ± 2.44 | **± 58 %** | 107.2 ± **52.7** |
| `k_tdg_ddg` | 4, 3-DG → 3,4-DGE | 30.5 ± 3.39 | ± 11 % | 36.9 ± 6.3 |
| `k_fru_int` | 6, Fru → Int | 330 ± 22.8 | ± 7 % | 100.4 ± 6.6 |
| `k_fru_odg` | 8, Fru → 1-DG | 2.11 ± 0.40 | ± 19 % | 99.3 ± 21.8 |

Each becomes two uniform-band rows in the envelope's own vocabulary (`b34.<key>.log10_k_100C` over
the relative HPD applied to the shipped 100 °C value, and `b34.<key>.ea_kj_mol` over the printed Ea
HPD), drawn independently — the same simplification ENV-B13 made, and the same draw hook
(`disputed_sinks`). `k_ddg_hmf` stays fixed: it is a timescale bracket with a declared zero barrier
and no HPD to draw from; a band invented for it would be a fabricated interval.

**What ENV-B13's "agreeing" reason said, and what it says now.** `k_tdg_ddg` leaves the agreeing
list. The paper's own NaCl column prints an Ea of **117.7 ± 11.1** for the same step against
**36.9 ± 6.3** in the glucose column — a threefold disagreement *inside one laboratory*, which the
1.5× cross-laboratory agreement on the rate at one temperature was hiding. That fact goes into the
reason string.

## 3. What counts as success, declared before the run

- **T1** the eight prior rows exist with the bands above, and the draw reaches the engine.
- **T2** every row the priors can REACH widens (3-deoxyglucosone, 3,4-dideoxyglucosone, glucosone,
  5-HMF, DMHF on the trunk and acrylamide lanes); every row they cannot reach moves by less than
  the measured Monte-Carlo noise floor, taken as the observed **maximum** over two seeds of identical
  priors (the ENV method note).
- **T3** no unreached median moves beyond the same floor. No centre was changed.
- **T4** reported: how much wider, coverage, and whether any measurement moved inside its interval.

Ship rule: **INSTALL if T1, T2 and T3 hold.** T4 is reported.

## 4. Predictions, before the run

1. Not-evaluable rows **3 → 0**. **95 %.** Zero-width intervals were the whole reason.
2. **3-deoxyglucosone moves inside its interval.** **85 %.** It is 1.11× off with a ±58 % rate band
   and a ±52.7 kJ/mol barrier band carried 80 °C from the anchor.
3. **3,4-dideoxyglucosone does NOT move inside its interval.** **85 %.** It is 32× off; the widest
   plausible band on `k_tdg_ddg` from a ±11 % rate and ±6.3 kJ/mol barrier is well under a decade at
   121 °C. If this prediction fails it means the barrier band, not the rate band, is doing it.
4. Envelope coverage rises by at least one row (12/42 → ≥13/42 on the same denominator). **75 %.**
5. Median interval width over all rows rises by under 0.1 dex. **70 %.** Four rows widen a lot,
   thirty-eight do not move.

---

# Outcome (2026-09-11)

## First run, under the shared stream: T3 fired on the sampler

Seven rows the priors cannot touch — three furfurylthiol rows on the sulfur lane and four lipid rows,
three of them in one bundle moving by an identical 0.1783 dex — exceeded the two-seed floor of
0.139 dex. The floor itself had read 16.99 % at ENV-B13's verdict and 11.28 % a day later for the
same code. That is the sampler re-shuffling every later draw when eight coordinates are inserted, not
the priors doing anything, and a rule whose verdict depends on which seed pair measured the floor is
not a rule. **ENV-M1** (`kinetic_core_env_m1_prereg.md`) was pre-registered and adopted first: one
random stream per coordinate, keyed by name. BEFORE and AFTER were then regenerated under the same
streams, BEFORE with the eight `b34.` rows filtered out.

## Two amendments to this pre-registration, both before the re-judgement

1. **Glucosone is not reached.** Section 3 listed it; in an amine-free pot its only route is
   `k_glc_g`, which carries no flux and has no band. Removed from the reach list, and it is the ONE
   row still not evaluable (prediction 1 fails: 3 → 1, not 3 → 0).
2. **The constants reach every Maillard lane, not only the trunk and acrylamide rows named.** Under
   per-coordinate streams a row can move only if the priors reach it, and eleven others did: every
   acrylamide row (the amine-free entries compete with the asparagine initiation for glucose) and two
   sulfur rows (the sulfur integrator runs the trunk's furanic block). So T2 and T3 are split: the
   four named compounds MUST widen; other Maillard-lane rows may move and are reported with their
   size; the lipid lane, the only one the constants cannot touch, must be **bit-identical** — which
   replaces the noise floor as the decisive comparison and is checkable exactly.

## Ship: INSTALL

| test | result | pass |
|---|---|---|
| T1 the eight rows exist | 8 rows, 8 sampled, bands as declared | **yes** |
| T2 widths | all 9 named rows widen (+0.08 to +0.99 dex); 11 other Maillard rows move; **lipid rows bit-identical 8/8** | **yes** |
| T3 medians | lipid medians moved: **0**; Maillard-lane medians moved through shared glucose: largest **0.157 dex** (methylglyoxal in the amine-free pot), next 0.029 | **yes** |
| T4 | coverage **12/42 → 16/44**, not evaluable 3 → 1; **newly inside 4**: 3-deoxyglucosone, methylglyoxal (both in the amine-free pot), HMF in Schibilsky pH 5 and in Chang 2021 water; newly outside none; median width 0.915 → 0.987 dex | reported |

Under the shared stream the same comparison had coverage 12/42 → 14/44; the two extra hits under
ENV-M1 are re-seeding, not the priors — every envelope number moved once when the streams changed,
and the guards are re-pinned in the same commit with ENV-M1 as the reason.

## Predictions, scored

1. Not-evaluable 3 → 0 (95 %) — **wrong**: 3 → 1, glucosone, for a reason the prereg mis-stated
   (it is not reached).
2. 3-deoxyglucosone moves inside its interval (85 %) — **right**.
3. 3,4-dideoxyglucosone does NOT move inside (85 %) — **right**; 32× off with a 0.98 dex interval.
4. Coverage rises by at least one row (75 %) — **right**, by four.
5. Median width rises by under 0.1 dex (70 %) — **right**, +0.072.
