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

## 6. Outcome (2026-09-08, variant (b) run the same night) — DO NOT SHIP, and a sharper diagnosis

**What was built.** `ch_dimer_release_mft` / `ch_dimer_release_fft` (MFTD → 2 MFT, FFTD → 2 FFT) in the
sulfur network on one shared constant `k_dimer_release`, inert (zero) unless a B17 report supplies it,
so every earlier wave reproduces bit for bit (the B11 discipline); its barrier is the dimerisation's
own measured 122.2 kJ/mol. Generator `generate_kinetic_core_b17_fit.py`: B16's 64-row objective, the
constant appended to the vector as B11 appended its oxygen consumers (24 free), two starts, the
600-evaluation budget; Laplace `--wave b17`; ship rule `generate_kinetic_core_b17_ship_rule.py`.

**The fit.** Both starts land on B16's optimum: cost 931.1 and 931.3 (B16: 931.0), with the release
constant driven to the floor of its band (log10 −8.26 and −8.09 per minute at 145 °C: off). The cost
slice along the coordinate is flat to four decimals over two decades either side; the Laplace σ is
4.8 × 10⁴ dex. The data neither want nor can see a disulfide release.

- **T1 failed**, as in B16: MFT at 100 °C peaks at six hours (110.6 → 97.4 µg/L from 6 to 12 h) where
  the pot keeps rising; the 360 and 720 min ratios are 0.85 and 0.97 dex low. FFT within 0.35 dex.
- **T2 failed**, as in B16: the fed-ribose MFT row moves +1.90 dex; 16 B9 rows move more than 0.3 dex.
- **T3 failed, and this is the finding.** The model's dimer shares are 0.04 / 0.35 / 0.90 % of the
  free thiol at pH 6 / 7 / 8 in Zhou 2023's pot (measured 8.6 / 6.5 / 9.6 %) and 0.39 % in Zhang
  2024's cysteine arm (measured 8.7 %; a later re-read of the paper's text, `zhang2024b_extraction.md`,
  finds this share is a figure read-off and the printed text supports only the ordering cystine >
  cysteine = glutathione, so Zhou 2023's shares carry the comparison): ten to two hundred times too little disulfide, with both
  dimerisation constants already on the upper edge of their bands. The dimerisation is not
  rate-limited; it is oxidant-limited. The ambient oxidant pool the lane charges (the B11 reservoir
  shipped inert) runs out, so the disulfide channel cannot hold the 7 to 10 % of the thiol the two
  laboratories find there, and making that channel reversible returns nothing because there is
  nothing to return.
- **T4 passed**: Yiltirak's four pots, median fold 13.6 (B9: 115), the same gain B16 showed.
- **T5 failed**: Wang 2022's 140 °C pot, MFT falls 1.32 dex from its 30-minute peak (the pot declines
  gently); FFT 0.11 dex.
- **T6 failed**: unidentified, slice flat.

**Variant (a)** (the saturable thioether sink on a pool browning makes) is the next run by this
pre-registration's own order; it was not run tonight. What T3 adds to its brief: the oxidant supply
is now a named suspect for the dimer share, which is a different quantity from the missing thiol,
and the two should not be conflated. A pot that holds 7 to 10 % of its thiol as disulfide with the
ambient oxidant exhausted needs either a larger reservoir (B11's vessel plumbing exists and ships
inert) or a second oxidant; that is a question for the laboratory experiment in the introduction's
section 8, where the disulfides are quantified in the same run.

Kept as a record: `kinetic_core_b17_fit_report.json`, `kinetic_core_b17_laplace_covariance.json`,
`kinetic_core_b17_ship_rule.{json,md}`, the two members. The engine keeps reading B9; the release
steps stay in the network at zero.

## 6b. Outcome, variant (a) (2026-09-09) — DO NOT SHIP, and the reason is structural

**What was built.** Three parallel branches of the deoxyosone decay, `ch_mele_from_dpo` / `_tdp` /
`_ddp` (DPO, TDP, DDP → 5 FRAG_C + 1 MELE), at rate `k_mele_site` = yield × `k_osone_decay` with the
carbonyl-sink family's barrier, inert (zero) unless a report supplies `mele_site_log10_yield`, so
every earlier wave reproduces bit for bit. The sites feed the measured thioether channel
(`k_thioether`, Hofmann 2002; K(T), Stack 2018), whose constants the fit cannot move. Generator
`generate_kinetic_core_b17a_fit.py`: B16's 64-row objective, the yield appended to the vector as B11
appended its oxygen consumers (24 free), band log10 −4 to 0.2, two starts, the 600-evaluation
budget; Laplace `--wave b17a`; ship rule `generate_kinetic_core_b17a_ship_rule.py` (T6 on the yield;
Zhang 2024's dimer share reported as figure-derived after `zhang2024b_extraction.md`, Zhou 2023
alone deciding T3).

**What happened.** Start 0 (B9's optimum, yield at the band centre) converged in 467 evaluations to
B16's optimum, cost 930.98, with the yield at log10 0.048 (1.1 sites per osone); start 1 (B8's
perturbation) reached 933.24 with the yield at log10 −3.96, its floor. The same cost at the two ends
of the band: the cost slice along the yield is flat to the fourth significant figure
(930.978 to 930.979 over ±1 decade) and the Laplace sigma is 8.6e4 decades. T1, T2, T3 and T5 fail
exactly as they did for variant (b) (the reference pot's MFT peaks at 6 h and falls; the fed-ribose
row moves +1.91 dex; the dimer shares are 0.03 to 0.93 % against Zhou's 6.5 to 9.6 %); T4 passes
(Yiltirak median fold 13.4). Verdict by the pre-registered rule: DO NOT SHIP.

**Why it could not have worked as pre-registered, and what that teaches.** Two facts, both visible
in the frozen numbers, not in the fit's noise.

1. *The pool has no source at the optimum.* B16's optimum (which this fit re-finds) carries
   `k_osone_decay` at log10 −8.59 per minute at 145 °C with a 123.5 kJ/mol barrier: the osone decay
   the sites were tied to is switched off. The fitted site constant is 2.9e-9 per minute; over a
   twelve-hour cook at 100 °C the pot makes of the order of 1e-9 of its sugar into sites. A yield per
   osone decayed can be anything when nothing decays. The pre-registration tied the sites to "the
   caramelisation and Amadori decay at a declared yield per unit of browning carbon"; the lane's
   browning carbon flows through `r_osone_decay_*`, and the earlier waves had already driven that
   flux to zero to hold the 145 °C fed pots. The right source, if the idea is tried again, is the
   flux the lane actually carries at the optimum (the ARP and pentose steps themselves, or the
   melanoidin pool the trunk makes), not a sink that is dead.
2. *The measured equilibrium releases the thiol at cook temperature.* Stack 2018's conjugation has
   K = 5.64 M⁻¹ at 19.4 °C and ΔH = −28.5 kJ/mol, so K = 0.45 M⁻¹ at 100 °C and 0.167 M⁻¹ at 145 °C.
   Even a pool equal to the whole sugar charge of the reference pot (0.1 M) would hold under 5 % of
   the thiol at 100 °C and under 2 % at 145 °C. A reversible thioether on this equilibrium cannot be
   the missing sink at cooking temperature whatever its site density; the channel earns its place at
   25 to 80 °C, where it was measured. Hofmann 2002's 80 % plateau at 80 °C is a covalent,
   effectively irreversible binding to melanoidin (CROSSPY-type), a different object from the
   quinone adduct Stack measured, and the lane lumps the two under one MELE.

**What the record now says about W7.** Both pre-registered structures are refused, for different
reasons: (b) because the model makes almost no disulfide (oxidant-limited); (a) because its source
is dead at the optimum and its equilibrium is too weak when hot. What the two share is the
diagnosis that the thiol loss at 100 °C is neither the disulfide nor a reversible adduct. The
reading of 2026-09-09 (`farmer1990_extraction.md`, `whitfield1988_extraction.md`,
`mottram2002b_extraction.md`; backlog "lipid-Maillard II") supplies the third candidate the
pre-registration did not name: the unsaturated-carbonyl adducts (thiophenes, thiapyrans, Michael
adducts) that halve the thiols when a lipid supplies electrophiles, an IRREVERSIBLE sink on a pool
the pot's own sugar fragments and any carried lipid make. And Xu 2010 (`xu2010_extraction.md`) says
the oxidant behind the disulfide share is internal to the pot, not the headspace. A variant (c),
irreversible addition to an electrophile pool sourced from the flux the lane carries, with the
dimer share left to an internal oxidant, is the next thing to pre-register; it is not run tonight.

Kept as a record: `kinetic_core_b17a_fit_report.json`, `kinetic_core_b17a_laplace_covariance.json`,
`kinetic_core_b17a_ship_rule.{json,md}`, the two members. The engine keeps reading B9; the site steps
stay at zero.

