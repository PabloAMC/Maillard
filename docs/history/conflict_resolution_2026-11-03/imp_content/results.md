# CF-041 (IMP in cultured / cultivated meat): results

Pre-registration: `prereg.md` (this directory). Script for T4:
`/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/7ae75f75-b2b3-4540-b5f1-5f645399920f/scratchpad/conflicts/imp_content/t4_remap.py`.
Pages read by eye: Joo2022.pdf p. 5 (methods) and p. 7 (text + Fig. 2); Kim2024.pdf p. 3-4 (methods) and p. 10 (Table 4).

## Units and basis as printed

| source | sample | basis printed | IMP cultured | IMP conventional (same paper) |
|---|---|---|---|---|
| Joo 2022 | 2D cell sheet, scraped, 3x PBS wash + 500 x g spin; 1 g pellet in perchloric acid; HPLC | "mmol/kg of the samples" (as-is; no drying step) | chicken 0.81, cattle 1.98 mmol/kg (text) | Fig. 2B by eye: chicken ~24, cattle ~14 mmol/kg |
| Kim 2024b | 3D construct, 1e7 cells/mL in gelatin hydrogel, ~90 % moisture; 5 g as-is, extract freeze-dried, 1H NMR | "Contents (mg/kg)" (as-is) | pig 0.11 mg/kg = 0.32 nmol/g | pork leg 792.19 mg/kg = 2.28 umol/g |

No per-protein, per-dry or per-cell basis in either paper. Gap as printed: 1.98 umol/g / 0.32 nmol/g = 10^3.8.

Two dossier errors found: (1) joo2022 sec. 2 says cattle TM IMP "roughly 2.5 to 3" read from the figure; Fig. 2B
shows ~14 (cattle) and ~24 (chicken) mmol/kg. (2) joo2022 omits the chicken CMT value (0.81 mmol/kg, printed in
the text). The kim2024b dossier matches Table 4 exactly.

## Tests

- **T1 (Kim's scale) -- passes.** Kim's pork 2.28 umol/g sits at the bottom edge of the beef references
  (koutsidis2008a/b 2.5-6.3 mmol/kg; muroya2019 D14 3.4 umol/g) and just below chicken breast (aliani2013 Table 2,
  2.7-5.3 umol/g); pork leg is expected lower than breast. A unit slip in Kim's table (H1a) would move the pork
  control too, and it does not. P(H1a) ~0.05.
- **T2 (Joo's scale) -- supports H1b, below the pre-registered 2x margin.** Joo's TM IMP exceeds every reference
  IMP by 2.2x (cattle, vs 6.3) to 4.5x (chicken, vs 5.3), and exceeds the WHOLE adenine-derived pool that IMP is
  drawn from: chicken breast total base-containing pool 8.3-15 umol/g (aliani2013 Tables 1-2), beef ATP at D0 6.5
  umol/g and IMP + inosine + Hx at D1 ~8.5 umol/g (muroya2019). Joo's TM is 1.6-1.7x above the pool, beyond its
  own error bar (chicken 24, bar to ~20). On a wet basis that is not possible; a dry basis at ~75 % moisture
  (24/4 = 6, 14/4 = 3.5 umol/g) would put both TMs inside the references. The margin is below the 2x I
  pre-registered, so I call this support, not proof. P(Joo's absolute scale inflated 1.6-5x, all samples) ~0.65;
  of which P(undisclosed dry basis) ~0.4, calculation / dilution factor ~0.25.
- **T3 (living-state IMP) -- passes at the low cell fraction only.** muroya2019 D0 (pre-rigor, ATP intact): IMP 78
  nmol/g. Kim's construct at 1e7 cells/mL x ~3 pL per myoblast is ~3 % cells at seeding (my estimate, not printed;
  up to ~10-30 % after 16 d growth). 78 x 0.03-0.3 = 2.3-23 nmol/g vs Kim 0.32: 7x at 3 % (within the decade),
  1.5 decades at 30 %. So scaffold dilution plus a no-post-mortem state explains Kim's corner without a slip if
  cultured cells hold IMP at or below pre-rigor muscle; the upper half of the cell-fraction range needs cells
  ~10x below pre-rigor muscle, plausible (no AMP-deaminase activation in fed, oxygenated culture) but not tested
  here. Also: 0.11 mg/kg from 5 g is ~0.55 ug IMP on the NMR, near a cryoprobe's quantification floor, so Kim's
  value carries roughly a factor 2-3 uncertainty (my judgement; the paper prints only SEM 41.7 pooled across
  columns).
- **Joo's own ratio (not a pre-registered test; it follows from T2).** CMT/TM = 0.034 (chicken), 0.14 (cattle).
  Scaled onto reference meat (2.5-6.3 umol/g beef, 2.7-5.3 chicken): cattle CMT ~0.35-0.9 umol/g, chicken
  ~0.09-0.18 umol/g. Still 5-10x above muroya D0, consistent with partial ATP -> IMP breakdown during Joo's
  scrape-wash-spin handling (~30 min without glucose). That handling step is the "harvest state" both papers'
  discussions invoke; neither measured it.
- **T4 (box effect on the frozen artifact) -- as predicted.** Remapping each frozen draw's cultivated IMP from
  [0.0003, 3.0] to [0.0003, 1.2] mM by the same quantile (exact for a same-seed regeneration; the recomputed
  all-declared ranking matched the stored one in 200/200 draws per programme): candidate naive ranking changed 0/200,
  engine arms untouched (IMP is in UNRANKABLE), every aggregate and both verdicts (T1, T1) unchanged.
  `naive_ranking_all_declared` changes in 43/200 draws per programme; IMP tops it in 130 instead of 115 draws.

## Verdict

The four decades are not one slip: decomposed (log10) as Joo's inflated absolute scale ~0.2-0.65, Kim's scaffold
dilution ~0.5-1.5, harvest state (Joo's pellet partly post-mortem, Kim's construct near living-state) ~1-2, Kim
near the NMR floor +/-0.4; the sum (2-4.5) covers the 3.8 observed. The two numbers measure different things (a
washed cell pellet vs a hydrogel product), both honestly labelled "cultured meat".

- Status **explained**, P ~0.6 that this decomposition is right in kind (P ~0.25 that a slip in Joo's absolute
  scale is larger than 5x, P ~0.15 that cultured myotubes simply hold far less IMP than pre-rigor muscle and
  dilution plays a minor role -- both still "explained" in kind, different weights).
- What would change my mind: Joo's raw data or a corrigendum giving TM inside 2-8 umol/g wet (H1b falls); a
  measured cell fraction in Kim's construct near 100 % (dilution fails); any third cultured-muscle IMP on a wet
  cell-pellet basis (the decisive measurement; it is already in the box's gap_map experiment).

## Recommendation

1. Joo values enter as ratios to its own TM, not as absolute mM. Box upper corner of cultivated IMP 3.0 -> 1.2 mM
   (Joo cattle ratio 0.14 x beef 2.5-6.3 umol/g at 75 % moisture; rounded up). Lower corner 0.0003 kept, labelled
   "hydrogel construct, product basis, near NMR floor". Effect on the frozen `measured_box_prediction` artifact: zero
   on every aggregate and verdict; 43/200 per-draw all-declared rankings move. Do not regenerate it now; the change
   lands at its next declared regeneration. `STUB_BOX` and `cultivated_tissue_invariance.json` are untouched.
   This is not tuning: no fit or score uses IMP.
2. Dossier fixes (joo2022 sec. 2-3): TM values ~24 / ~14 mmol/kg read from Fig. 2B; chicken CMT 0.81; the TM
   exceeds the adenine pool (cite this record). Ledger notes of both entries: replace "~4 decades" with the basis
   and pointer to CF-041.
3. Red team should attack: my bar reading of Fig. 2B (render at higher dpi); the ~3 pL cell volume and the 3-30 %
   cell fraction; whether a dry basis is ruled out anywhere in Joo or in its method source (Ismail et al. 2020a);
   whether 1.2 mM should instead stay at 3.0 given only P ~0.65 on H1b (the box is a sensitivity range, so
   keeping 3.0 is the conservative alternative and changes nothing that is scored).
