# Pre-registration: ENV-B13, the eight trunk constants the envelope treats as certain (written 2026-09-10, before anything ran)

## 1. Why

ENV-B18 found that the Monte-Carlo envelope had no prior row for the pyrazine step. Checking the
rest of the table found something larger: **the whole dicarbonyl and furanic-sink block is absent
from it.** Eight constants — the two sugar entries, the two osone steps, and the four sinks — appear
in no prior row, so they contribute exactly zero to every published interval. The envelope asserts
them with certainty.

That would be defensible if they were well determined. Two of them are declarations their own
authors flagged, and a second laboratory has now refuted both:

| constant | what ships | what a second laboratory measures | apart by |
|---|---|---|---:|
| `k_da_sink` | **0, with a blank barrier** — carried as "a prediction the data may reject" | 130e-3 /min in roasted hazelnut at 160 °C | the prediction is **rejected** |
| `k_go_sink` | 32.6e-3 /min with the **barrier FIXED TO ZERO by its authors** | 18 / 61 / 290e-3 /min over 150–170 °C | rate agrees 1.87x; **the zero barrier is refuted** |
| `k_odg_da` | 12.2e-3 /min at 180 °C, barrier 150.8 | 895e-3 /min at 160 °C | **466x** |
| `k_hmf_self` | 8.97e-7 /min, **barrier zero by declaration**, from 0.9 % lost in 7 days at 5 °C | 21e-3 /min at 160 °C (hazelnut); 0.111 /min at 180 °C (Gökmen) | **2.3e4** and **1.2e5** |

Five panel rows score hydroxymethylfurfural. Their intervals today contain no contribution from a
sink constant that two other laboratories put five orders of magnitude away.

## 2. What this wave does, and what it refuses to do

**It does not move a centre.** Which laboratory is right for a plant-protein cook is not settled by
either of them: one is an aqueous amine-free glass, the other a dry, lipid-rich, whole-tissue nut,
and the two disagree on the dicarbonyl ORDER as well as the rates. Installing either number as the
new centre would be choosing a matrix by preference.

**It gives each disputed constant a declared band spanning the disagreement, and integrates across
it.** That is what an envelope is for. A band that spans two laboratories is an honest statement of
what is known; a point with no band is not.

| constant | band, and its basis |
|---|---|
| `k_da_sink` | log10 k flat from an effective zero to the measured 130e-3 /min. "Somewhere between nothing and what the second laboratory measured" is exactly the state of knowledge. |
| `k_go_sink` | the RATE is left alone, since the two agree to 1.87x. The BARRIER is banded 0 to 150.8 kJ/mol — from the refuted zero to the steepest barrier the trunk itself carries. The 20 °C window is too narrow to fit a credible barrier and this wave does not pretend otherwise. |
| `k_odg_da` | log10 k flat across the 466x, at the common temperature. |
| `k_hmf_self` | log10 k flat across the span from the shipped value to Gökmen's, which is the widest of the three. Its barrier stays zero by declaration: one temperature each, no Arrhenius. |

The four undisputed constants (`k_glc_g`, `k_g_go`, `k_tdg_ddg`, `k_ddg_hmf`) get prior rows too,
but **fixed and unsampled**, with the reason recorded — three of them are the trunk's first
cross-laboratory agreement, inside a factor of two, and the fourth has one determination. Listing
them makes the table complete, so the next person checking for a missing block finds a row rather
than a silence.

## 3. What counts as success, declared before the run

- **T1** all eight constants appear in the prior table; the four disputed ones sampled, the four
  agreeing ones fixed with a stated reason; every band's endpoints traceable to a printed number.
- **T2 decisive.** Every hydroxymethylfurfural interval gets WIDER. Not one row anywhere may narrow.
- **T3** the centre draw still reproduces the deterministic prediction exactly, so no median moves.
- **T4 reported.** How much wider, and whether any measurement moves from outside its interval to
  inside.

Ship rule: **install if T1, T2 and T3 hold.** T4 is reported.

## 4. Predictions, before the run

1. T2 holds and the HMF intervals widen by **more than a decade**. **75 %.** A flat band across five
   orders of magnitude on the only sink a compound has should dominate everything else in its
   interval.
2. At least one HMF measurement moves from outside its interval to inside. **60 %.** The furanic
   channel prints its own warning that HMF is expected to be over-predicted for want of this sink;
   if that is right, a band reaching the second laboratory's much faster sink should reach down
   towards the measurements.
3. This makes the model look worse and is still correct. **95 %.** A wider interval on a row that
   already misses is not an improvement in accuracy, and publishing a narrow interval that omits a
   known five-decade disagreement is not honesty.
