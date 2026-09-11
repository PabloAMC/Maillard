# Pre-registration: wave B42, two papers and two fed half-lives at 50 °C (written 2026-09-11, BEFORE the probes ran)

## 1. Why

Two papers arrived against the reading list: Zamora et al. 2015 (thermal breakage of unsaturated
aldehydes; `zamora2015b_extraction.md`) and Gobert & Glomb 2009 (glucose degradation with lysine at
50 °C; `gobert2009_extraction.md`). The first turned out to be a hexanal *source*, not the hexanal
loss law the reading list hoped for, and is recorded as a candidate step. The second prints two
**fed half-lives at 50 °C** — 3-deoxyglucosone 40 h, glucosone 8 h, each 42 mM with 42 mM lysine at
pH 7.4 — which are the lowest-temperature aqueous measurements of those pools on disk, and a direct
test of two barriers this model carries with a caveat.

## 2. The two probes

Both pots: 42 mM of the fed compound + 42 mM lysine, water, pH 7.4, 50 °C, integrated to 48 h (3-DG)
and 8 h (glucosone); the model's half-life is read off the fed compound's own series.

## 3. Predictions

- **P1, fed 3-DG.** The model's half-life at 50 °C is **under 4 h** against the printed 40 h — at
  least tenfold too short. The reason is named in advance: the formic-acid exit carries Martins'
  30 kJ/mol barrier (flagged since B1 as conflicting with Knol 2010's 84), which barely slows the
  step between 100 and 50 °C, and B41's pH term *raises* it at pH 7.4.
- **P2, fed glucosone.** The model's half-life at 50 °C is **under 30 min** against the printed 8 h —
  more than tenfold too short — and it sends most of the glucosone to glyoxal, where the paper finds
  0.07 %. The reason: B21's aqueous glucosone → glyoxal barrier is 4.2 kJ/mol, consistent with zero,
  fitted at 110–140 °C. A barrier near zero cannot be right if the step is fast at 120 °C and slow at
  50 °C.
- **P3.** Nothing moves: no constant, no artifact but this record and the two dossiers.

## 4. What follows if the predictions hold

Neither probe is a fit. What each licenses is a pre-registered wave with its own hold-out: for P1,
adopting Knol 2010's measured 84 kJ/mol on the formic-acid exit in place of Martins' 30, judged on
the Martins fit rows (which are at 100 °C and would be unaffected at the reference) and on Leitzen;
for P2, giving `k_g_go` a barrier from the two temperatures now on disk (Hamzalioglu at 110–140 °C,
Gobert at 50 °C) and re-judging B21's glyoxal rows. Both are named in `EXPERIMENTS.md` and neither
is done here.

## 5. Outcome (written 2026-09-11, after the probes)

**Both predictions held, by far more than they claimed.**

| probe | printed (Gobert & Glomb 2009, 50 °C, pH 7.4) | model | factor |
|---|---|---|---:|
| fed 3-DG 42 mM, half-life | **40 h** | **≈ 1.2 h** with lysine charged (acrylamide lane), **≈ 0.5 h** without (trunk lane) | 30–80× too short |
| fed glucosone 42 mM, half-life | **8 h** | **≈ 3 min** | ≈ 160× too short |
| glyoxal from fed glucosone at 8 h | **0.07 %** of the charge | **69 %** at the peak (0.1 h), gone by 8 h | ≈ 1000× |

One limitation surfaced on the way: charging lysine routes a pot to the acrylamide lane, which does
not compose the dicarbonyl block, so the glucosone probe had to run amine-free on the trunk lane;
the 3-DG probe was run both ways and the lane changes the half-life by 2.5×, not the conclusion.

**What this says, in order of certainty.**

1. **B21's aqueous glucosone → glyoxal barrier cannot be near zero.** Hamzalioglu's 0.33 /min at
   120 °C and Gobert's 8-hour half-life at 50 °C are two temperatures on the same step; a barrier
   that connects them is of the order of 100 kJ/mol, not 4. The model's glyoxal at 50 °C is three
   decades too high because the constant fitted at 110–140 °C was carried down 70 °C flat. This is
   the second time a barrier declared "consistent with zero" has failed a low-temperature hold-out
   (the glass `k_ddg_hmf` was the first, in B41). Next wave: `k_g_go`'s barrier from the two
   temperatures now on disk, with B21's own glyoxal rows and the pyrazine total as the hold-outs.
2. **The 3-DG pool's lifetime is still too short, now at 50 °C by 30–80×,** after B41 fixed its
   120 °C behaviour on fed pots. The exit that dominates at 50 °C is Martins' formic-acid step at
   30 kJ/mol, which the parameter table has flagged since B1 as conflicting with Knol 2010's
   84 ± 14. A barrier of 84 would slow that exit about 15× more between 100 and 50 °C. Next wave:
   Knol's barrier on the formic-acid exit, judged on the Martins fit rows (at the 100 °C reference,
   unaffected) and on Leitzen and the fed pots.
3. **Zamora 2015 is a source, not a sink.** 2,4-decadienal breaks to hexanal (11.5 %) with a 21 kJ/mol
   barrier and a half-life of tens of minutes at 120–200 °C. A candidate lipid-lane step; the hexanal
   loss the sixteen papers showed remains without a law on disk.

**Nothing moved.** No constant, no artifact but this record and the two dossiers.
