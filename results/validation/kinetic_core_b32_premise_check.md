# Wave B32, killed by its own premise check before a pre-registration was written (2026-09-10)

## What was proposed

B29's oxygen axis left one finding standing: the trunk sends **54.6 %** of the fed-Amadori pot's
Strecker flux through the oxidative entries against **64.1 %** of the sugar pot's, where Hofmann &
Schieberle 2000b demand the opposite order — their Amadori pot is the more oxygen-sensitive (an
air/argon ratio of 9.2× against 3.5×). The proposal was to refit the branching at the Amadori
compound, `k_ama_g` against `k_ama_odg` and `k_ama_mgo`, on two within-study ratios already in the
corpus.

**Nothing was fitted. No constant moved.** Three probes were run first, and the third refutes the
proposal outright.

## Probe 1: `k_ama_g` cannot do it, at any value

| shift on `k_ama_g` | ×    | fed-Amadori share | sugar share | in the right order? |
|---:|---:|---:|---:|:--|
| −1.50 | 0.032 | 13.3 % | 17.5 % | no |
| −0.50 | 0.316 | 41.1 % | 49.0 % | no |
| **0.00** | **1** | **54.6 %** | **64.1 %** | **no** |
| +0.50 | 3.16 | 62.6 % | 74.5 % | no |
| +1.50 | 31.6 | 66.9 % | 82.4 % | no |
| +3.00 | 1000 | 67.4 % | 83.6 % | no |

Five decades, and the gap never closes — **it widens**. B21's declared transfer band on this constant
is ±0.5 dex, so even the licensed range was never the question. The reason is structural: the Amadori
step is SHARED. The sugar pot reaches the Strecker aldehyde through Glc → AMA → G, the very same
constant, so raising it lifts both pots together.

## Probe 2: the one asymmetric lever is structurally dead

`k_glc_g`, glucose → glucosone, is the only oxidative entry the sugar pot has and the fed pot does
not. Moving it over **six decades changes the shares by nothing at all** — 54.6 % and 64.1 % to four
figures at every value. Its rate constant is 5.34e-8, and at 100 °C it carries no flux. The
hypothesis that this entry was inflating the sugar pot's share was wrong.

## Probe 3, which is the answer: THE TWO POTS ARE NOT BRANCHING DIFFERENTLY

After 120 minutes at 100 °C the fed pot has 0.124 mM of Amadori compound left and no glucose. The
sugar pot has 2.13 mM of Amadori compound and **73.9 mM of its glucose still unreacted**. They are
not two branchings; they are the same reaction at two very different extents. Run the sugar pot on:

| pot | minutes | Amadori left | glucose left | oxidative share |
|---|---:|---:|---:|---:|
| fed Amadori | 120 | 0.124 | 0 | **54.6 %** |
| glucose + Gly | 120 | 2.128 | 73.94 | 64.1 % |
| glucose + Gly | 480 | 1.131 | 41.91 | 62.7 % |
| glucose + Gly | 1920 | 0.149 | 5.92 | 62.5 % |
| glucose + Gly | 7680 | 7.2e-5 | 0.0029 | **55.4 %** |

**At matched extent the two shares are 54.6 % and 55.4 %.** There is no branching asymmetry in this
model to refit. A fit would have moved a measured constant to absorb a difference in how far the
reaction has run, and it would have reported the resulting agreement as chemistry.

## What this says about B29's finding, which is now sharper

The finding survives and improves. The trunk's oxygen sensitivity enters at **one node that both pots
share**, so the model can only separate the two pots by REACTION EXTENT — and the extent effect runs
the wrong way for Hofmann's data, since the pot that is further along is the *less* sensitive one.
Hofmann's Amadori pot is the more sensitive. No re-parameterisation of the existing topology can
produce that, because no constant in it is on a route the Amadori pot has and the sugar pot does not.

**What would produce it is a route the Amadori compound has and the sugar pot does not reach**, which
is a structural claim needing a source, not a fit. That is the shape of the remaining gap, and it is
a more useful statement than "the branching is wrong" — which was my proposal, and which this check
shows is not a thing that can be true of this model.

## A defect in the measurement B29 reported, found on the way

B29's "oxidative share" is `(air − argon) / air` on the summed Strecker aldehyde. Under argon the
**non-oxidative product goes UP**, not down — AKM rises by 81 % in the fed pot and 93 % in the sugar
pot — because the Amadori compound that is no longer drained to glucosone flows down the deoxyosone
route instead. That is sound competitive branching, but it means the quantity is a NET difference
between two competing channels and not a share of anything. The share of the aldehyde that is
oxidative by construction is `AKG / (AKG + AKM)`: **74.9 %** and **81.3 %**. The order is wrong on
that measure too, so nothing above depends on which is used — but B29's number should be read as what
it is, and the ship rule's wording is corrected accordingly.

## Verdict

**NOT PRE-REGISTERED AND NOT RUN.** The premise fails. Three probes cost minutes; the wave they
prevented would have fitted a constant to an extent difference and shipped it as a branching fix.
