# Wave B30 (W8) gating test: DO NOT SHIP -- and nothing was fitted

*T0 GATES THE WAVE: the model must at least FALL from pH 4.5 to 6.5 on the one pot measured at both. Nothing is fitted unless it does.*

| arm | pH 4.5 | pH 6.5 | ratio 4.5 / 6.5 |
|---|---:|---:|---:|
| as shipped | 0.0383 | 0.07767 | 0.493 |
| thiolate loss off | 0.03843 | 0.09901 | 0.388 |
| hydrosulfide branch off | 0.0383 | 0.07767 | 0.493 |
| **measured** | 0.150 mol % | < 0.001 | **>= 150** |

The model gives MORE thiol at the higher pH, where the measurement collapses. The sign is wrong and the ratio is out by about 304x. No slope was fitted.

## Attribution

- The thiolate loss carries almost none of it: switching it off moves the ratio by -0.10.
- **The hydrosulfide branch is not active on this pot at all** (True): switching it off changes the answer by nothing to four figures. the two-branch sulfide mechanism was built for the deoxypentosone route (r_ddp_mft_hs) and the furfural route (r_fur_fft_hs). The NORFURANEOL route has no hydrosulfide partner -- r_nf_mft and r_nf_mp3p are single steps in neutral H2S -- so on the one pot in the corpus with a measured pH PAIR, the lane's pH mechanism is structurally absent.
- a hydrosulfide branch on the norfuraneol steps pushes the SAME way as the one that already exists: more hydrosulfide at higher pH means faster addition means MORE thiol at pH 6.5, where the measurement wants at least 150x less. The collapse is not in the nucleophile. It is in the substrate or in the sulfide budget, and neither is modelled.
