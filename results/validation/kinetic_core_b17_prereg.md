# Pre-registration: wave B17, the sink structure (written 2026-09-08, not yet run)

## 1. Why

B16 (`kinetic_core_b16_prereg.md` section 6) showed that no re-tuning of the thiol sinks as first-order
Arrhenius steps serves both 100 °C and 145 °C: the reference pot at 100 °C keeps making MFT and FFT for
twelve hours while every fitted version peaks near six, the 145 °C fed pots break when the sinks are
weakened, and Wang 2022 shows the sinks are too strong at 140 °C as well. The candidate table
(`docs/validation/thiol_sink_candidates.md`) says why: every measured sink has a partner and saturates
or reverses, and the model's dominant sink has neither. B17 changes the FORM of one sink, not its
numbers.

## 2. The variants, in the order they are run

**(b) The disulfide gives the thiol back.** One new step and one new fitted coordinate:

    ch_dimer_release_mft   MFTD -> 2 MFT        k_dimer_release   (log10 k at 145 °C, free)
    ch_dimer_release_fft   FFTD -> 2 FFT        k_dimer_release   (shared)

with the release sharing the dimerisation's own barrier (122.2 kJ/mol, Zhang 2026's peptide dimer,
`MEASURED_EA_OVERRIDES`) so the equilibrium constant, not the two rates, carries the temperature
dependence. The dimer decay (`k_dimer_decay`, 2e-10 per minute, effectively zero) is left as it is:
the dimer becomes a reservoir instead of a grave. What the data bear on it: Kumazawa 2003's apparent
loss rate at pH 6 halves when the cook doubles (0.073 to 0.042 per minute), which a reversible or
oxidant-limited sink produces and a first-order one cannot; Zhou 2023's dimer holds 6.5 to 9.6 % of
MFT as thiol equivalents across pH; Zhang 2024's dimer share responds to the oxidant, not the thiol;
Gigl 2021's dimer reaches about 21 % of the initial thiol at 96 h.

**(a) A saturable covalent sink on a pool browning makes**, run only if (b) fails T1: the matrix
electrophile pool `MELE`, today a charged input that is zero in every pot, produced from the sulfur
lane's caramelisation and Amadori decay at a declared yield per unit of browning carbon, so the
thioether channel (Hofmann 2002's measured rate and Stack 2018's measured equilibrium) runs with a
pool the pot itself makes and exhausts. Its one new coordinate is the site yield. Its assumption,
the sites per browning carbon, has no direct measurement; Hofmann 2002's plateau (about 80 % of
400 µg FFT bound by 12.5 g/L melanoidin) is the anchor and the band is wide. This is why (a) runs
second.

## 3. What runs

The objective is B16's: B9's 54 primary-evidence rows, the seven Schieberle 2000 Table IV
within-study ratios at 100 °C, and Zhai 2021's three TTCA rows (64 rows). Free coordinates: B9's 23
plus `k_dimer_release` (variant b) or the site yield (variant a). Two starts, B8's protocol, the
600-evaluation budget, quick mode for the search and a careful evaluation of the optimum. Then the
Laplace covariance at the optimum and the slice profile of the new coordinate.

Hold-outs, never in the objective: the four Hofmann 1998 pH-5 level bundles, the four Yiltirak 2026
pots, Wang 2022's 140 °C series (shape only), Bolton 1994, and the four returned Hofmann pH-3 and
pH-7 rows.

## 4. What counts as success, declared before the run

- **T1, the reference pot at 100 °C.** MFT and FFT at 360 and 720 min within 0.5 dex of Table IV
  as within-study ratios to the 30 min point, and both still rising between 6 and 12 h (no peak
  before 12 h). B16 failed this with a peak at 6 h.
- **T2, the 145 °C fed pots do not break.** Every B9 fit row within 0.3 dex of its B9 residual.
  B16 broke the fed-ribose row by 1.9 dex.
- **T3, the dimer shares.** Zhou 2023's and Zhang 2024's dimer-to-thiol ratios within 0.3 dex.
- **T4, another laboratory.** Yiltirak 2026's median fold error below 20 (B9: 115; B16's
  weakened sinks reached 14 while breaking T2).
- **T5, the shape at 140 °C.** Wang 2022's MFT and FFT decline from 30 to 180 min by less than one
  decade (the pot declines gently; B9 loses three decades).
- **T6, identification.** The Laplace covariance identifies the new coordinate (sigma below one
  decade), and the slice profile is not bound-limited.

Ships as B17 if T1, T2 and T6 pass; T3 to T5 are reported. If T1 or T2 fails, the variant is kept
as a record like B16 and the other variant runs. If both fail, the outcome section says so and the
experiment in the introduction's section 8 is the only remaining route.

## 5. What it will not do

It will not move the formation steps, the pH structure or the barrier bands; it will not read any
hold-out; it will not ship a variant that passes T1 by breaking T2.

## 6. Outcome

*(not yet run; the generator `generate_kinetic_core_b17_fit.py` is to be derived from B16's, with the
new reaction and coordinate added to the sulfur network and the B23 vector as B10 added the route
barriers)*
