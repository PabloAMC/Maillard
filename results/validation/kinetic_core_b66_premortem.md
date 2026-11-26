# B66 premortem (independent, read-only), 2026-11-26

*An independent agent read the draft pre-registration
([kinetic_core_b66_cysteine_thiamine_refusal_prereg.md](kinetic_core_b66_cysteine_thiamine_refusal_prereg.md)), the
B65 records, the thiamine-H2S rule (PR #158) and the code, and ran two short probes. It assumed the change had gone
wrong and ranked the causes. This is its report, condensed; a partner's private numbers quoted in the draft are left
out. V = verified in code, data or a run; I = inferred.*

## Ranked causes

| # | What goes wrong | P | Fix taken |
|---|---|---|---|
| 1 | The draft carried a partner's private numbers and pot into a public record (V) | 0.9 if merged as written | Removed from the prereg and the test; "two published laboratories" |
| 2 | The blast radius is not where the draft expected: no scored level row changes, two independent directional claims drop out, and pinned guards and quoted counts break (V for the spreads, I for the counts) | 0.85 | Guards re-pinned with a dated comment, quotes updated, the miss named in the PR, the forecast miss recorded |
| 3 | The rule was narrower than its evidence on the donor: glutathione has two laboratories (Grosch Table II, x4.4; Zhang, x15), and TTCA releases cysteine in the pot; neither was refused (V) | 0.6 | `applies()` on thiamine with cysteine, glutathione or TTCA |
| 4 | "The value 3.7 does not matter" is wrong for mixed pots: 3.7 refuses only above a 23 % thiamine share of MFT; Zhang's x19 is confounded by the sugar route (V on tracked evidence: no difference) | 0.35 | Factor 7.0 (Grosch Table II, thiamine the only MFT source); cut-off about 10 %, declared |
| 5 | The reason text claimed "probably higher" below a cysteine : thiamine ratio; the record says the effect follows concentration, and in Bolton's pot the engine is already above the measurement (V) | 0.5 | No direction and no ratio rule in the text |
| 6 | The falsification was argued with the per-arm threshold and the wrong S3 clause (a with / without ratio is judged by omega; x2.7 is inside it) | 0.3 | Clause and comparison corrected in the amendment |
| 7 | Two new config keys re-pin every stored config digest (lesson 2026-11-18) | 0.4 red CI | Digest-bearing artifacts regenerated before pushing; the step-5 look files left alone |
| 8 | Three per-arm rules judged separately could add up; the spreads seen are bimodal (V) | 0.05 | One sentence in the module |

## Blast radius (rule on)

| Surface | Pots with thiamine and a donor | Newly refused |
|---|---|---|
| Scored panel | Bolton 1994 | MFT (0.54 dex at x3.7); the row was already not evaluable, so no scored row changes |
| Directional, independent | WANG12-01 pot C; MOD08-PH-01 both arms | MFT: one agree and one miss become not evaluable |
| Directional, other arms with both | PH-03, CER07-PH-01 / 02, WANG12-02 | none (they read FFT or furfural, 0.00 dex) |
| Within-study contrasts (36 pairs) | none | none |
| P3 dry run | none | none |
| Cultivated-tissue surface | every box pot | MFT; one unit test updated |

## Answers to the draft's open questions

- **Was the constant fitted on cysteine + thiamine pots?** No for the level (V): the anchor is thiamine alone (Hofmann
  1998). The pair pots enter only as a route share and within-pot ratios. So the pair pots are the right ones to
  refuse and thiamine-alone pots are not implicated.
- **Is `r_hmp_mft` the right object?** Yes. Thiamine's cleavage would also move 3-mercapto-2-pentanone and the derived
  H2S release, and Arnold 1971 reports no extra thiamine loss with cysteine. Grosch's FFT does not rise with cysteine,
  so leaving FFT shown is consistent. The interaction with the other two rules is clean (all own-integrations pass
  every switch off; caches are keyed by the config digest).
- **A second rule for the dimer split?** Not now (P about 0.8 that deferring is right): see the prereg's amendment.
