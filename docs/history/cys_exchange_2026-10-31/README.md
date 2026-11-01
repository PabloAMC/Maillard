# Catalysed cysteine exchange (CF-076): rivals, power check, feasibility (2026-10-31; EX-054)

This is an archived analysis record, not a generator. `prereg.yaml` was written at 13:44; the first run was at 13:52. No red team was run, because nothing here moves the engine; one is due before any exchange law changes.

**Erratum (2026-11-01; this record is otherwise unchanged).** The k_r bands in `prereg.yaml` rest on two misreadings:
- E1 `k0_100` basis: "Kallen 1971 TC relaxation 0.009-0.03 /min at 25 C" is the approach-to-equilibrium relaxation of
  Kallen's Fig. 4, not the hydrolysis constant; with the companion kallen1971b, TC k_r is 5.4e-5 /min at pH 5.6 (t1/2
  ~9 d), 120-270x slower (kallen1971 dossier, arithmetic 3, corrected 2026-11-01). "Schubert 1939 (ring opening as fast
  as free Cys to I2 at RT)" measures C-S ring opening, not cysteine release. The band [1e-5, 1] /min survives
  numerically on a rewritten basis, so E1's fate here (it failed on c_aa) does not move on this account.
- `c_aa` [0.3, 30] and `c_P` [4, 300] (Fife 1991, general base, linear) are superseded by kallen1971b: ammonium
  Brønsted c_aa 0.1-0.3 /M at pH 6 and ~1-3 at pH 7 (+-1 log), no catalysis below pH ~5.5, saturation at x3.5-5. E2 /
  E4's acid limb (a in [0, 1], k_r falling with pH) has the wrong sign below pH 5.6 on the N-H analogue, where k_r
  RISES with pH.

The corrected bands were re-run on 2026-11-01 (cluster 5 of the conflict-resolution programme, EX-055; record in
`docs/history/conflict_resolution_2026-11-01/cluster5/`). That re-run SUPERSEDES the verdicts below: with Kallen's
k_r(pH) sign, in-situ pH + K(pH) is feasible and E4c is admitted provisionally (P 0.45), E1c is rejected, and E0 is
now rejected in two laboratories (Zhang 2025b's blank added as T8).

**Outcome: no rival is admitted.**
- **E0, the engine's exchange: REJECTED, provisionally (single lab).** It releases cysteine 2 - 30x faster than TTCA is lost in Zhai 2021's dilute pots, which allow k_r(100 C, 10 mM) <~ 1e-3 /min. The T4 gate has power 0.93 / 0.99 against E1 / E2.
- **E3, E0 plus a second cysteine sink:** rejected on the same gate.
- **E2, in-situ pH + K(pH):** rejected (max |z| 5.7).
- **E1, catalysed exchange (k0 + k_aa[AA] + k_P[P]):** rejected by the rule (4.15). Feasible on T4 - T7 only at its band corner: amino-acid catalysis >= 20 - 30 /M, Ea at the 40 kJ/mol floor.
- **E4, E1 + E2:** fails only T3 (3.54). T3 also fails the no-exchange control (3.71), so it looks like a shared-model failure; without T3, E4 is feasible (catalysis ~3 /M). Reading it as inconclusive is POST HOC.
- **Controls:** NC1 (instant exchange) and NC2 (none) fail as required.

**Notes:**
- The conflict is inside ONE laboratory: Zhai 2021 dilute vs Zhai 2019 concentrated. Cao 2017 (T7) passes any fast-enough exchange and adds little.
- Kang 2026 Fig. S3 at 100 C is a digit-for-digit copy of Zhai 2021 Fig. 6 at pH 7, so it is not independent data. Kang's measured pH drift (5.5 / 7 / 8 -> 5.1 / 4.9 / 4.7) set E2's band.
- Ea sits at its band floor in every feasible point, so temperature dependence is part of the conflict.

**Next:**
- Designed: TTCA at 10 vs 100 mM in water at 100 C, x +/- 0.2 M phosphate, x +/- 0.2 M non-reacting carboxylate, with initial free cysteine over 2 - 10 min and the final pH.
- Ask Zhai for the 13C5-xylose label fraction in TTCA.
- An unlisted rival to test: TTCA re-forming after the quench in the concentrated pots (a storage / re-assay check).
