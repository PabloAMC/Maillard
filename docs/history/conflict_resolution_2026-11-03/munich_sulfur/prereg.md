# Pre-registration: Munich sulfur cluster (CF-016, CF-065), 2026-11-03

Written before the page renders and the engine computation below were run. Exploration done before writing it:
text-layer greps of the five candidate PDFs and the dossiers of munch1997, munch1998, hofmann1995, schieberle1996,
schieberle1998, schieberle2000 (to find which tables exist). Builds on cluster 4 sec. 9 (CF-016) and cluster 3 sec. 2
(CF-065) of 2026-11-01; does not repeat their operator / speciation work.

## CF-016 (furfural in the ribose 10 mmol + Cys 3.3 mmol pot: Hofmann 67.5 vs Schieberle 5.2 ug per pot)

**H1a, a third print exists.** Test: page reads of every on-disk table that reports this pot (munch1997 Tables 4-5,
munch1998 compound list, hofmann1995 (AEDA), schieberle1996 Table 3, schieberle2000 Table III). Decision: a furfural
LEVEL (not an FD factor) for this pot within 0.2 dex of one print sides with it. Falsified if no table prints one.

**H1b, Hofmann's ribose cell is a 10x slip (true 6.75).** Test outside the cell: Hofmann 1998's own running text on the
ribose / glucose furfural ratio (table: 67.5 / 6.9 = 9.8). A slip in the ribose cell predicts the text gives ~1x or
contradicts the table. Decision: text ratio within 1.3x of 9.8 -> the Hofmann single-cell slip is rejected
(P <= 0.05). A column-wide scale error is NOT tested by this (it is cluster 4's "Hofmann high").

**H1c, Schieberle's cell is a unit slip (52 is ug per pot, not per 100 mmol).** No discriminating test is known in
advance; recorded as an observation only (whether the 52 is a cell new to the 1998 chapter, i.e. absent from the
1996 Table 3 that the rest of the aqueous column reprints).

**H1d, FFT closure.** Using Munch 1997 T5 expt 1 (furfural 1 mmol + Cys 1 mmol -> FFT 25.3 ug), the second-order
scaling to the ribose pot's Cys 33 mM gives the time-mean furfural needed to make the pot's FFT (12.1 ug) by the
furfural route. Decision: if needed / printed > 10 for BOTH prints, closure does not discriminate (and the FUR route
supplies < 10 % of in-situ FFT under either print: a model finding, not a verdict).

## CF-065 (FFT from furfural + sulfide: Hofmann 0.48 mol % vs Zhai ~24 %)

**H2a, order in H2S.** Munch 1997 T5 expts 2 / 3 (H2S 0.25 vs 1 mmol, furfural 1 mmol). n = log(FFT3/FFT2) / log 4, with
the +/-10 % Hofmann replicate bound on each. Decision: n in [0.8, 1.25] -> linear (the second-order law cluster 3 used
to transfer Hofmann's constant is supported by a row outside the conflict); n < 0.7 -> saturating (cluster 3's
transfer is invalid at 10 mM); n > 1.3 -> threshold / fixed sink. Either non-linear outcome reopens cluster 3's
1.1 dex condition correction.

**H2b, fraction of furfural ending as FFT.** Computed directly per furfural and per H2S charged. If it is < 1 % at
both doses, neither reagent is depleted by FFT formation and the yield is rate-limited, not supply-limited.

**H2c, cysteine acts only as an H2S donor.** The engine's measured Cys -> H2S constant (k_cys_h2s, Zheng 1997,
extrapolated to 145 C, pH 5) predicts R = (integral of Cys-derived H2S) / (integral of fed H2S) over 20 min, the same
for any H2S-only route. Observed: FUR route R_FUR = 25.3 / 550.8 (Munch T5 1 / 3); NF route R_NF = 0.05 / 0.19 mol %
(Hofmann T4). Decision: both within 0.3 dex of R_pred -> H2S-only donor supported; R_NF matches but R_FUR > 0.3 dex
below -> a cysteine-specific furfural sink (thiazolidine) indicated; both above R_pred -> an H2S-independent Cys route
or k_cys_h2s too low at 145 C. The fed-H2S headspace partition (cluster 3) applies to both arms and cancels to first
order. Falsifier of the "donor-only" engine structure: R_FUR / R_NF outside [0.5, 2].

Expected (before running): n ~1 (P 0.75); R_pred between 0.02 and 0.3 (P 0.7); CF-016 stays open (P 0.85).
