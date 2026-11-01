# Pyrazine formation: rival structures, power check, scoring, red team (2026-10-31; EX-052)

This is an archived analysis record, not a generator. `prereg.yaml` was written before any run; `analysis_README.md` is the investigator's own summary; `red_team.md` is the independent review, which saw the prereg, code and data, not the conclusions.

**Outcome after the red team:**
- **ENGINE BUG (code mismatch, CF-077).** `network.py` gives the three condensations (`r_akg_pz`, `r_akm_dmp`, `r_ak_mpz`) one shared `k_cond` in plain mass action, so the engine's mixing ratio R = MPZ^2 / (PZ * DMP) is 1. Its docstrings (`network.py`, `parameters_pyrazine.py`) promise the statistical cross-pair factor, MPZ = 2 sqrt(PZ * DMP), i.e. R = 4. Methylpyrazine is therefore under-predicted 2x against the engine's own design, and 3.5x / 7x against Leahy's lysine-glucose / lysine-ribose pots.
- **P0, the engine: REJECTED, two labs.** Reineccius (S1, the mixing bound), Ho network (S2, calibrated rows with R >= 1.9) and Ghent (S7, the sign of the pH 4 vs 7 contrast on fed MGO).
- **P1, one shared constant with the factor 2: INCONCLUSIVE by the rule** (one lab rejects it with a pre-registered test). P(wrong) ~0.9: Lys-Rib z >= 4.9 at the top of the response-factor band.
- **P1x, pair-specific constants: inconclusive, leaning rejected as the sole explanation.** It has an unlisted sibling P1y (a pair-specific dihydropyrazine fate). The red team's post hoc label test, label(MPZ) = mean of label(PZ) and label(DMP) for every pairing rival, fails on Zhou 2024 at 30 min (0.371 vs 0.489) even under heavy [M-H]+ leakage.
- **P2b, an extra MPZ route: provisional,** admitted mostly by its flexibility.
- **P2a_aa, Strecker-aldehyde addition: rejected** by 13C labels at 180 C (two labs).
- **The negative control fails, as required.**
- **Process findings:** S3's leak band and the `near_loq` band were set after looking at the data they decide (circular); the zhou2023c ARP + Gly row went unscored; S5a's pyrrole mass balance assumes free aminoacetone.

P(the engine's pyrazine structure is wrong in a way that changes predictions) ~0.95.

**Minimum engine change (its own PR, owner approval since shipped numbers move):**
1. The statistical factor on the cross pair (2 k_cond on `r_ak_mpz`, or k_cond / 2 on the homo pairs).
2. A declared one-sided MPZ band of +0.3 to +1.1 dex until pairing and the extra route are separated.
3. Leahy's pH slope below 7 not applied to the fed-dicarbonyl Strecker steps (a band instead; S7).
