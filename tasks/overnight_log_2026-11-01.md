# Overnight 2026-11-01: which reactions happen under which conditions

**Owner's goal for the night:** "figure out what reactions happen in which conditions"; address the key problems; delegate; stay on track.

## Plan (in order of value)

1. **The shared vocabulary.** docs/reference/pathway_atlas/families.yml: 32 reaction families, condition bins, status and evidence vocabularies. (Done, main thread.)
2. **Literature pathway-condition atlas.** Six agents by family group (early; hub A; hub B; Strecker / N-heterocycles; sulfur; other lanes), one YAML per family. Every cell is sourced (dossier + location + lab + evidence kind); unknown regimes are explicit.
3. **Engine regime map.** A generator runs the SHIPPED engine over the same condition grid (pH x T x sugar x amine x O2 x buffer) and records flux per family (an engine-reaction -> family map). Overnight compute, <= 4 workers.
4. **Engine vs atlas comparison (morning).** Every cell where the engine's dominant route differs from the literature's, or where the engine runs a route the literature shows absent: a ranked list of engine mechanism errors by condition.
5. **A3, phase 1, in an isolated worktree off main.** A frozen Config object with a from_globals bridge, a job-carried config applied in workers (spawn-safe), and a fork-vs-spawn equivalence test. Bit-identical; all gates + full suite; a DRAFT PR for the owner (it does not merge overnight). Nearly every engine fix found on 10-31 / 11-01 waits on A3.
6. **Housekeeping:** cluster 5 (cysteine exchange) and its register edits; PR #33 merge when CI is green; the ingest branch's full suite and PR.

**Not tonight:** the contrast benchmark (tomorrow; it uses the atlas); the thiol leave-one-lab-out (its draft is parked in scratchpad/lolo_draft).

**Compute:** <= 9 cores overnight (the regime map at 4 workers; A3's suites; the ingest suite).

## Checkpoints
- 01:20 vocabulary written; premortem launched; atlas agents launching (read-only on all but their own files).
- 01:40 premortem adopted, three changes:
  1. atlas cells keyed on 3 axes (pH, T, water) with tags; markers per family; status capped by evidence kind; a verbatim `quote:` per cell; a validator and a 10 % red team in the morning;
  2. the engine grid replaced by the engine-reaction -> family map tonight, recipe-paired runs tomorrow, and the freed cores given to the thiol leave-one-lab-out;
  3. A3 phase 1 narrowed, with a switch-coverage test, fork-vs-spawn with all switches flipped, and a 05:00 stop rule.
  
  Compute capped at 8 cores; staging by explicit path; separate worktrees: A3 isolated, LOLO in scratchpad/lolo_wt.
- 01:45 launched: 6 atlas agents (F01-05, F06-10, F11-15, F16-21, F22-28, F29-32); the engine family map; A3 phase 1 (isolated worktree, draft PR); the LOLO finish + launch (lolo_wt, 5 workers). Pending: cluster 5; PR #33 CI.
- 01:35 engine family map done (docs/reference/pathway_atlas/engine_family_map.yml; checker passes; read-only).
  - **Coverage:** 180 reaction keys; 27 / 32 families have an engine reaction, 25 a shipped one.
  - **Families missing:**
    - no reaction: F09 Namiki, F12 alpha-oxidative cleavage, F18 lipid Strecker, F20 extra pyrazine routes, F32 DHA;
    - none shipped: F17 the O2 Strecker acid route, F21 pyrroles / pyridines;
    - partial: F02 (no Heyns arm), F04 (TTCA formation off), F23 (no hexose intact-C5 entry), F30 (no crossover into Maillard routes).
  - **Condition-dependence red flags in the CODE** (each cites a line):
    - the sulfur lane computes pH, but its 29 trunk steps get no pH term (TRUNK_PH_TERMS=False);
    - opposite pH signs for the same chemistry: hexose enolisation rises with pH in the trunk, pentose r_pent_tdp falls with [H+], and every "base" step uses cysteine's pKa even for the alanine ARP;
    - furans: no trunk pH term; r_ddg_hmf Ea = 0; the DMHF steps are T-independent;
    - the water-activity term sits on k_schiff only, and glycation has no a_w term;
    - O2: OX = 1.0 in every sulfur-lane pot, inert-gas pots included;
    - ALL sulfur formation steps share one Ea (64.1), so the F23 / F24 route mix CANNOT shift with T (a structural reason for the route-mix-vs-T failure);
    - fitted barriers on their bounds (20 / 260 kJ/mol);
    - no buffer / phosphate catalysis term anywhere.
- 02:xx atlas F11 - F15 done: 67 cells, every quote grepped in its dossier (scratchpad/validate_atlas_f11_f15.py).
  - **F11 beta-cleavage:** the main acetic-acid route at pH 6 - 8, 37 - 120 C (fed 1-DG in two labs, plus 18O2).
  - **F12 alpha-oxidative cleavage:** operating, never dominant; no data at 100 - 140 C.
  - **F13 furfural:** from the intact C5, falling 39x to > 1500x from pH 4.5 to 6.5; the norfuraneol pH optimum flips between 100 and 140 C (candidate conflict).
  - **F14 amine-free degradation:** substantial at 120 C in phosphate, needs phosphate at pH 7.
  - **F15 redox:** single-lab and qualitative.
  - **Engine mismatches (agent's arithmetic, to verify):**
    - r_tdg_fa sends ~8x too much 3-DG to formic acid;
    - acid yield is flat in pH (fed 1-DG gives 4x less at pH 5);
    - r_aa_frag destroys acetic acid against the Davidek plateau (CF-069);
    - no O2 C3-acid route; no phosphate catalysis;
    - norfuraneol decay has no O2 / metal dependence;
    - Strecker-coupled reduction runs only with cysteine;
    - no furfural -> furfuryl alcohol step.
- atlas F06 - F10 done: 75 cells, quotes grepped (scratchpad/atlas_hub_check.py); cells tagged formation_flux vs standing_level.
  - **F07 1-deoxyosone:** DOMINATES formation at pH >= 6 in phosphate: Beck 20 : 1 in-situ; Davidek labelling + fed (>= 67 - 70 % of DFG at pH 8); Desclaux fit 51 : 1. Three labs, hexose / phosphate only.
  - **F06 3-deoxyosone:** 20 - 40 % of formation below pH ~5, ~5 % at pH 7, yet still the larger STANDING level at pH 7.
  - **F08 oxidative osone:** shown at 37 - 50 C only.
  - **F09 Namiki C2:** shown only in ethanol with alkylamines; every aqueous cell unattributed.
  - **F10 retro-aldol:** absent or minor from fed hexose deoxyosones at 37 - 50 C.
  - **KEY engine mismatches:**
    - dpo_share = 0.30 comes from a STANDING ratio but is used as a FORMATION split; flux data put the 1-deoxyosone share at 0.5 - 0.98;
    - k_ama_tdg / k_ama_odg / k_ama_mgo share one pH exponent (0.69), so the hexose split cannot move with pH (data: 0.25 - 0.45 dex / pH);
    - r_arp_tdp is labelled acid-catalysed while Martins' k4 rises with pH;
    - no O2 lever beyond k_ama_g / k_glc_g; r_g_go barrier ~0;
    - no pentose glyoxal source and no Namiki step;
    - r_dpo_c2c3 rests on a hexose analogue the Glomb lab refuted.
- atlas F29 - F32 done: 69 cells, quotes grepped. This agent used `not_stated` for unreported pH / sugar / buffer / water, which is outside the vocabulary: add `not_stated` to every bin in families.yml (morning validator).
  - **F29 acrylamide:** operates 120 - 200 C, pH 4 - 8, down to a_w 0.34 (2 labs kinetics, 3 labs foods); no labelling on disk; fructose > glucose (3 labs).
  - **F30 lipid crossover:** lipid in cysteine + pentose pots (120 - 145 C, pH 5 - 6.5) cuts MFT 2 - 143x (4 labs); the dienal + NH3 -> 2-pentylpyridine route is dominant (labelling, 2 labs).
  - **F31 melanoidins:** labelling in 3 labs (sugar + amine enter together, 70 - 125 C); covalent thiol binding shown only at 25 - 30 C and O2-dependent.
  - **F32 DHA cross-links:** one study only.
  - **Engine mismatches:**
    - acrylamide starts only from glucose (no fructose / sucrose / HMF / 3-APA); a phosphate pH law applied to every buffer (citrate at pH 4 ~7x lower); irreversible amine scavenging vs ~90 % reversible;
    - **the LIPID LANE is a direct sum with its covalent sink off, so the engine predicts NO MFT change with lipid, against a 2 - 143x drop in 4 labs** (high impact for plant-based meat);
    - the melanoidin C/N is fixed at 8 vs a measured 6.6 - 21;
    - the thiol sink has no O2 dependence (and CF-081);
    - DHA is absent.
  - **Register points:** Yiltirak 2026 contradicts half of CF-031; the friedman1990 dossier is a browning paper, not Friedman's lysinoalanine work; Xu 2011 vs Farmer 1992 on Maillard antioxidants (no CF).
- atlas F01 - F05 done: 96 cells, quotes checked. Normalised the vocabulary: ph not_measured -> not_stated; protein_bound_nterm added to the amine bin.
  - **F01 Schiff base:** a fast reversible pre-equilibrium, K 1.1 - 2.1 /M at 37 - 54 C (2 labs).
  - **F02 Amadori:** operating at 37 - 130 C, Ea 83 - 97 kJ/mol (4 labs), ~0.5 dex / pH; Heyns absent from fructose + casein.
  - **F03 amine release:** DOMINANT (fed DFG returns 71 - 96 % of its glycine, 3 labs); reversion to sugar is minor at 90 - 120 C; phosphate speeds both.
  - **F04 thiazolidine:** fast, reversible, K ~60 /M (3 labs; the rest analogues).
  - **F05 glycation:** in every regime; its laws rest on 1 - 2 labs.
  - **Engine mismatches:**
    - r_schiff is irreversible, with no pH / phosphate term;
    - no phosphate catalysis of Amadori decay (4 labs show 1.3 - 12x);
    - no reversion step, no Heyns step;
    - the fed pentose ARP decays release no amine while the same lab measures 25 - 40 % glycine release, and the disabled r_pama route assumes full return (the two pentose routes contradict each other);
    - TTCA formation is off by default and has no pH law;
    - CF-072 is two-sided (three other papers DO return lysine);
    - no lactose step, no pH / a_w term for glycation.
- atlas F22 - F28 (sulfur) done: 156 cells, quotes grepped; every `dominant` cell rests on labelling / fed. Two empty python3 calls outside the wrapper (no code ran).
  - **F22 cysteine -> H2S:** thermolysis rises 9 - 17x from pH 5 to 9 at 80 - 110 C, Ea 123 - 135 (2 labs within 1.5x).
  - **F23 MFT from an intact C5:** the MAIN MFT carbon route at 95 - 168 C (Cerny, CAAS, BTBU labelling); hexose MFT is also a C5; at 95 C, pH 5, free norfuraneol carries only 7 %.
  - **F24 MFT from C2 + C3:** <= 1 - 6 % in labelled pots at every T (Liu's 19 % is an EI-MS artefact, ~6 %).
  - **F25 FFT:** furfural + H2S dominant at pH 4 - 6, 95 - 145 C; FFT collapses by pH 7.
  - **F26:** mercaptoketones collapse while thiazoles / thiophenones rise from pH 4.5 to 7.
  - **F27:** thiamine supplies 54 - 98 % of MFT at dilute loadings (3 labs); glutathione a weak donor.
  - **F28 thiol losses:** unattributed in hot water; O2 / cysteine effects shown only at <= 50 C.
  - **Engine mismatches:**
    - hexoses reach MFT only via C2 + C3 (the C5 entry r_odg_dpo is off), while labelling says intact C5;
    - Cerny pot D: too much MFT runs through free norfuraneol;
    - **FFT formation RISES with pH in the engine (r_fur_fft_hs) while the measured yield FALLS (a sign error)**;
    - k_mft_decay has no cysteine / O2 / matrix term;
    - disulfides are terminal;
    - 3-mercapto-2-pentanone forms only from thiamine;
    - no buffer catalysis (Mottram 2002: 2 - 8x);
    - H2S release above 110 C extrapolates ~20x high against the one 160 C point.
  - **Unregistered conflicts:** Arnold vs Grosch (cysteine effect on thiamine); Hincelin vs Bolton (xylose); Liu 2023's 19 % (an EI-MS artefact).
- atlas F16 - F21 done: 98 cells, quotes grepped. ALL 32 FAMILIES MAPPED (~560 cells).
  - **F16:** the dicarbonyl Strecker aldehyde is dominant in fed pots at 98 - 170 C, pH 3 - 7, and needs no O2 (2 labs on ranking; yields 20x apart, CF-073).
  - **F17:** in air the O2 route carries 77 - 82 % of the Strecker ACID; the fed Amadori compound gives x9 more aldehyde in air (one lab).
  - **F18:** lipid Strecker 37 - 180 C, all quantitative cells from Seville.
  - **Engine mismatches:**
    - B18's Strecker steps carry the Leahy pyrazine pH term (x0.091 at pH 5) while the fed-dicarbonyl aldehyde peaks at pH 5 / is largest at pH 3.2;
    - no 3-DG / glucosone donors, no Strecker acid, no amine exit;
    - B29's O2 factors are unfitted (argon refused);
    - lipid is a direct sum (x1).
  - **Possible conflict:** Zamora 2013 vs Hidalgo 2005, x10 - 30 (with CF-078).
- next: the atlas validator (sole writer of atlas files now), then a 10 % red team and the engine-vs-atlas comparison in parallel.
- atlas validator: 561 cells; 53 -> 7 failures (mechanical fixes only, one cap-breaking dominant downgraded); new gate scripts/ci/pathway_atlas_gate.py. Launched: the fix-up of the last 7 + CF entries; the 10 % red team; the engine-vs-literature comparison.
- cluster 5 done (scratchpad/cluster5/verdict.md).
  - **Bands corrected:** the TC k_r at 25 C is 5.4e-5 /min (pH 5.6), so the kallen1971 dossier's value was the relaxation, 120 - 270x fast; k_r rises with pH up to 5.6.
  - **Rejected:**
    - the engine exchange E0 (+E3), now in TWO labs (z 33 on the dilute gate; Zhang 2025b blank);
    - catalysed exchange (corrected bands), secondary-amine catalysis, and a quench / storage operator.
  - **Feasible:** in-situ pH with Kallen-shaped K(pH) and k_r(pH). E4c admitted PROVISIONALLY (z 1.45); the TTCA ring-N pKa sits at its band floor, Zhai 2019's pH is latent, and Zhai 2021 Fig. 7a at face value rejects it.
  - **CF-076: EXPLAINED, provisional (P 0.45);** k_r(100 C, pH 4.5) ~2 - 8e-4 /min, 40 - 100x below the engine.
  - **CF-011:** explained in part (phosphate P 0.5; ELSD quantitation P 0.3).
  - **Pending:** the register / dossier edits (kallen1971 correction; evidence-map K4; CF-076 / 011 YAML) go to the next single-writer pass after the atlas fix-up. No engine change before a red team and the pKa.
- atlas red team (98 cells, seed 20261101): WRONG 10 / 98 (6 / 68 random = 8.8 %, CI 4 - 18 %); any issue 60 % (mostly tags / bins).
  - **`absent` 5 / 16 wrong (31 %):** non-detects without a detection limit, and product capture.
  - **`dominant` 2 / 25.**
  - **Fit to drive changes:** the labelling / fed cells of F23 / F24, F16 / F17, F19, F01 - F06.
  - **NOT fit until re-audited:** `absent` cells as switch-offs; F25 dominant outside Cerny's fed-furfural pot; F31 dominant; F30 quantities; cross-file T / water bins.
  - The comparison agent was told to apply these rules.
  - **Pending fix-up (after the comparison):**
    - the gate recomputes T bins from T_C;
    - `absent` requires a detection_basis, else it becomes `not_detected` (a new status);
    - fix the 10 wrong cells; add Zhao 2019b's labelled FFT cell;
    - a README limitation: many quotes are dossier paraphrase, so the gate checks the dossier, not the paper.
- engine-vs-literature comparison committed. Two agent leads corrected by engine probes (the FFT pH 'sign error' is a 10 - 70x slope error; dpo_share 0.30 is latent, not shipped). Top errors: no hexose intact-C5 MFT route; the enolisation split is pH-frozen; no Strecker aldehydes + the pyrazine pH law on the Strecker step; no buffer catalysis; MFT 99.9 % via free NF; lipid inert. Next: a red team of the ranking; the atlas fix-up.
- red team of the engine-error ranking: E1 / E10 band-now fixes REFUTED (B58 already said do-not-ship for TTCA; E1 needs a new Glc + Cys step); E2 weakened; E3 a refusal. Order: E4 (owner approval) -> E6+E11 refit -> E10 via B59b -> E1 -> E2; E5 / E7 / E3 after A3. Amendment in ENGINE_VS_LITERATURE sec. 5.
- recount on the fixed atlas: settled 31 -> 25, unsettled 217 -> 233, unknown 1800 -> 1790; all 22 pooled claims survive (P12 / P13 firmer); E18 narrowed to >= 145 C, E14 narrowed; no mismatch changes kind.
- 02:36 A3 phase 1 green and bit-identical: draft PR #34 (42 ConfigKeys; 9 tests incl. the AST coverage scan and spawn-vs-fork with all keys flipped; unit 1195 / scientific 108 / 14 gates pass; core panel byte-identical). NOT merged (owner review). An independent code review launched. Not covered: script-level statics (TDP_EXIT_STATIC etc.) still need fork, so no persistent pool yet.
- A3 PR #34 independent review: APPROVE WITH CHANGES. applied() correct (all bindings restored, nesting OK, not thread-safe); bit-identity holds. Issues:
  1. MEDIUM: build-phase keys (shipped_b39, shipped_b22) can be installed but give an inconsistent engine (tables built at import); make install refuse non-default build keys.
  2. MEDIUM: the AST coverage scan misses whole classes of new switch (mode strings, argument-built holders, plain numbers / dicts, setattr by variable, in-place dict edits); fix with "every uppercase global read by a function is a key or on an explicit allowlist" (~60).
  3. LOW-MEDIUM: hash(Config) is not stable across processes; add digest() before caching.
  4. LOW: bands declared but not checked; types not checked.
  5. LOW: the joint fit still builds its Config in the worker.
  6. LOW: key granularity (whole route dataclasses) vs the A2 contract's field-level keys; decide before phase 2.
  7. NIT: a test leaves sys.path modified.
  
  The spawn test is meaningful for 25 of 42 keys only (17 have no reader in the probe).
  
  **Existing bug found:** build_choice_robustness.py:99 sets only parameters_lipid.Q10_ASSUMPTION; the rates use it, but lipid.py:580 REPORTS q10_used from the other binding (a reporting mismatch).
  
  Owner decisions before phase 2: items 1 and 6.
- 03:xx next-phase plan written, premortemed (circularity, power, re-tuning compensations, A3 blocking, the owner bottleneck) and amended into STATUS_AND_PLAN. Started: the benchmark's frozen split.
- 03:45 E4 corrected: the red team misread Leahy (14 = pyrazine; MPZ = 215, shipped 170 matches); E4 expected REFUSE -> declared band + designed experiment; the order becomes CF-077 then A3 phase 2. Both preregs committed (E6 + E11 runs after EX-056).
