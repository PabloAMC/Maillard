# Thiol loss and cysteine fate: rival structures, power check, scoring, red team (2026-10-31; EX-050)

This is an archived analysis record, not a generator: the scripts are standalone reduced modules (A2 contract sec. 3), run once by an investigation agent. `prereg.yaml` was written at 12:50, before the first run at 12:58. `red_team.md` is the independent review, which saw the prereg, code and data but not the conclusions.

**Outcome after the red team:**
- **(a) MFT loss: INCONCLUSIVE.** The engine's loss law (A0) is rejected; no rival is admitted.
  - The patent's spiked-standard time course (Nestle, cysteine 0.25 M) needs lambda(110 C) <~ 0.02 - 0.05 /min, against the engine's 0.095. One source, no replicates.
  - Van Seeventer's 50 C storage data do NOT test a processing loss law (storage and processing sinks are different channels), so "two labs" does not hold.
  - Cerny pot D is a ROUTE problem: the loss law moves it by at most 0.3 dex against a 1.1 dex miss.
  - P(engine k_mft_decay too fast by >= 3x at <= 110 C in a cysteine-rich pot) ~0.75; as a pot-independent constant up to 145 C ~0.45.
  - Not listed, and live: the van Seeventer loss minimum near 50 mM cysteine; a reversible MFT-cysteine mixed disulfide; a shared redox buffer (cysteine + reductones).
- **(b) Cysteine fate: the engine's exchange (B0) is rejected provisionally; "slow exchange" (B1) is REJECTED.**
  - B1's ceiling was derived from the same xylose contrast it was scored on.
  - The same lab's concentrated pots (Zhai 2019) and Cao 2017 need exchange >= 2.6x faster.
  - With the xylose-consumption test included, the joint score ties (B0 -32.1, B1 -32.4, B4 -31.9).
  - What survives (P ~0.75): one pH-flat k_r and K cannot fit both the dilute and the concentrated pots. This is registered as a conflict. Catalysed exchange and in-situ pH in unbuffered pots are the next test.

**Minimum engine change:** none structural. `k_mft_decay` is to be carried as a declared single-lab band with lambda(110 C) <= ~0.02 - 0.03 /min (Ea not taken from storage data), in the next evidence-audit PR after A3. `k_ttca_cys` is untouched until catalysed exchange is tested.
