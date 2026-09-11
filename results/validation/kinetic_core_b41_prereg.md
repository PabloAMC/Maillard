# Pre-registration: wave B41, the pH term on the formic-acid exit only (written 2026-09-11, BEFORE any run)

## 1. Why

B40 put a pH term on both 3-deoxyglucosone exits and the hold-out answered precisely. With the
term on, the fed pot's peak landed at 25.5 min (P1 held), Leitzen's 3-DG stayed at 1.36× (the
2–6× I feared was wrong), 3,4-DGE went 32× → 6× and HMF improved on five pots — and
**methylglyoxal on Leitzen went 1.28× → 33×.** Every gain came from the term on the formic-acid
exit, which is Martins' own k6 measured at two pH values on the very step the trunk carries. The
loss came from the term on `k_tdg_mgo`, which was *declared* from Martins' k5 — a lumped
"3-DG → fragments" step — and transferred to Kocadagli's amine-free methylglyoxal route. The
hold-out rejected that transfer and kept the other. B41 keeps what the data kept.

## 2. The declaration

`THREE_DEOXY_EXIT_PH` carries **`k_tdg_fa` only**: 10^(0.77·(pH − 6.8)), band (0.65, 0.89), from
Martins 2003 Table 3 k6 at 100 and 120 °C. `k_tdg_mgo` keeps no pH term, and the reason is on
record: the one hold-out that measures methylglyoxal in an amine-free pot below pH 6.8 says the
transfer is wrong by 26×.

## 3. The fit

B39's generator, rows, sigmas, bounds and starts, with the one-exit term on. Leitzen never read.

## 4. Predictions

- **P1.** Fed rows fit and the fed-3-DG peak stays inside 20–60 min (B40 put it at 25.5 with two
  exits slowed; with one, it moves earlier and I predict it stays above 20).
- **P2.** Leitzen: 3,4-DGE under 10×, 3-DG at or under 3×, HMF within 3 % of 11.9× or better,
  **methylglyoxal within 3 % of 1.28×** — the row that failed B40, now decisive.
- **P3.** No row now within 3× leaves the band; the HMF rows B40's candidate brought in (Schibilsky
  pH 8, Chang water) stay in, since they came from the formic-acid term.
- **P4.** At least three of the five coordinates PINNED; the reverse pair weakens relative to B39
  because the second exit no longer trades against them.

## 5. Ship rule

SHIP if P1, P2 and P3 hold. P4 reported. Frozen pair under `_b41_baseline/`. If it ships: the term's
flag defaults on with the one-exit table, B39's literals become the B41 optimum (`SHIPPED_B39`
True), the envelope gains the five `b39.` Laplace rows and drops ENV-B34's printed band on
`k_tdg_ddg`, the panel headline is re-pinned with the reason, and B39's and B40's records stand as
the two steps that got here.

## 6. Outcome (written 2026-09-11, after the run) — SHIP

Artifacts: `kinetic_core_b41_fit_report.*`, `kinetic_core_b41_ship_rule.*`, frozen pair under `_b41_baseline/`.

| prediction | result |
|---|---|
| P1 fed rows fit, peak inside 20–60 min | **HELD**: peak at **23 min**, cost 0.81 on twelve rows, χ²_red 0.12; every maximum within 0.1 dex, every Zhang ratio within 0.06 dex |
| P2 Leitzen, never read by the fit | **HELD on all four rows**: 3,4-DGE 32.4× → **6.6×**; 3-DG 1.11× → **1.27×**; HMF 11.9× → **9.3×**; methylglyoxal 1.28× → **1.01×** — the row that failed B40 now improves |
| P3 no row leaves the band | **HELD**, and two enter: within-3× **10/45 → 12/45**, out-of-sample **9/44 → 11/44** (Schibilsky pH-8 HMF 3.6× → 1.6×, Chang-water HMF 3.2× → 2.6×) |
| P4 ≥ 3 pinned | HELD: σ k_tdg_ddg 0.06, k_ddg_tdg 0.18, k_ddg_dgal 0.24, k_dgal_ddg 0.26, k_ddg_hmf 0.06 dex; verdicts {'log10_k_tdg_ddg_100C': 'PINNED', 'log10_k_ddg_tdg_100C': 'PINNED', 'log10_k_ddg_dgal_100C': 'PINNED', 'log10_k_dgal_ddg_100C': 'WEAK', 'log10_k_ddg_hmf_100C': 'PINNED'} |

**Installed.** `trunk_conditions.THREE_DEOXY_EXIT_PH_TERM` defaults on with the formic-acid exit
alone; the rejected `k_tdg_mgo` declaration is kept beside it as `THREE_DEOXY_EXIT_PH_REJECTED_B40`.
`parameters_dicarbonyl.SHIPPED_B39` is True and `FROZEN_B39` carries this report's optimum. The
envelope gains five `b39.` rows at the fit's own Laplace σ and drops ENV-B34's printed band on
`k_tdg_ddg`. The panel headline is re-pinned in `tests/scientific/test_core_headline_guards.py`
with this wave as the reason.

**What the three waves together established.** The 3-deoxy limb's error was never one slow step. It
was a one-way step where the chemistry runs both ways, a missing epimer, an exit measured at pH 6.8
and applied at pH 5, and a 3,4-DGE → HMF rate carried from a 160–200 °C glass with a barrier fixed
to zero. Fed pots on a small network pinned all five constants where the level-scored fits could
not (B38's prediction), and a hold-out the fit never read rejected the one declaration that did not
transfer and accepted the rest, row by row.
