# Wolf 1977a — findings (2026-11-01)

Paper: Wolf, Thompson & Reineccius 1977, "Initial losses of available lysine in model systems", J. Food Sci. 42(6):1540-1544,
DOI 10.1111/j.1365-2621.1977.tb08420.x (Crossref-verified: title, authors, volume, issue, pages match the PDF; the DOI is
not printed on the PDF). Dossier: `data/lit/extraction_dossiers/wolf1977a_extraction.md`. No alias line needed (the file
stem `wolf1977a` matches the dossier stem, the same convention that covers `wolf1977.pdf`; dossier_coverage_gate PASS
with the main-checkout articles dir).

Computations (stdlib, run through ./scripts/docker_maillard.sh; nothing in the repo read at runtime or written):
- `wolf1977a_rivals.py` -> `wolf1977a_rivals.json`: shipped B20 constants (k_glyc 10^-4.707 L/(mmol min), Ea 96.8;
  k_flp_decay 10^-1.856 /min, Ea 97.1; k_flp_cml 10^-3.317, Ea 113) in Wolf's Table 2 / Table 3 pots, with and without the
  lysine return; plus each structure refitted to Table 2 and transferred to Table 3. Glucose: 1:1 into FLP plus an
  amine-free sink calibrated to Wolf's printed 4 % -> 2.5 % in 30 min at 105 C. Moisture (not printed) scanned
  w = 0.08 / 0.10 / 0.12 g water / g dry.
- `wolf1977a_surface.py`: Wolf's Case 1-3 surfaces evaluated; the engine's k_glyc·[G]_water at the same conditions.
- `wolf1977a_design_endpoints.py`: rival k_app = -ln(C/C0)/t at Wolf's design end points (80 C / 360 min, 125 C / 1.5 min).

## 1. CF-072 verdict

**Wolf 1977a does NOT decide CF-072. P(it discriminates the two rival structures beyond between-lab error) = 0.10.**

Rivals: R (shipped) FLP -> 3-DG + LYSP; N (Nguyen 2016) FLP -> products, no lysine back. Observable = FDNB available lysine.

| Wolf observable | R predicts | N predicts | observed | R - N (dex) | decides? |
|---|---|---|---|---|---|
| Table 2, 105 C, 4 % glc, 54 min, fraction left (shipped constants, w 0.08-0.12) | 0.32-0.42 | 0.05-0.14 | 0.43 (reslurry 0.52) | +0.5 to +0.8 | only with k_glyc pinned; see below |
| same, each structure refitted to Table 2 (k_glyc free) | rms 0.06 dex; k_flp_decay -> 0 (collapses to N) | rms 0.06 dex; k_glyc 0.3-0.5 dex below shipped | — | ~0 | no |
| same, R with Wolf's own phase-one k (0.025 /min) + shipped k_flp_decay | steady free fraction 0.46 | 0.26 at 54 min | 0.43 | +0.25 | weak (one lab, SD 0.2-1.0 mg/g) |
| Table 3, 95 C, 1 % glc, 75 min (shipped) | 0.72-0.79 | 0.62-0.71 | 0.68 | +0.05 to +0.07 | no |
| order expt, 105 C, 18 min (initial slope) | k_glyc·G | k_glyc·G | 0.024 / 0.026 /min | 0 | no |
| design end points 125 C / 1.5 min, k_app, 3.68 % glc | 0.20-0.30 | 0.22-0.32 | 0.31 (Case 1) | -0.03 to -0.04 | no |
| design end points 80 C / 360 min, k_app | 0.00021-0.0029 | 0.00033-0.0060 | 0.00072-0.0017 | -0.2 to -0.3 (signs of misfit differ by glucose) | no; both miss the glucose slope by ~0.7 dex |
| plateau depth vs glucose (4 % 0.43, 1 % 0.69) | k_d/(k_d + k_glyc G): 0.46 / 0.66 with Wolf's k | no plateau | 0.43 / 0.69 | — | R vs bare N: yes in sign; R vs N + inaccessible pool: confounded (10 C cooler, glucose at 1 % unprinted, C0 differs 0.12 dex) |

Reading:
- The plateau with sugar still in excess (2.5 % glucose at 30 min; ~4 mol glucose lost per mol lysine) **rejects a BARE
  no-return structure** (P 0.8). It does not reject N plus an inaccessible lysine pool or N plus an FDNB response to
  Amadori lysine — and the engine already applies `available_fraction(amine_band)` to matrix lysine, so the engine's N
  variant is "N + pool" by construction.
- Between R and N+pool the paper is ~neutral: likelihood ratio R : N+pool ~1 (every observable within 0.07 dex once
  k_glyc is free within its between-lab band). R : bare N ~2.
- Tension carried from the companion: at 125 C, Wolf 1977a's k (0.30 Case 1, 0.57 Case 2 at 4 %) with the shipped
  k_flp_decay (0.10 /min at 125 C) gives R a steady free fraction 0.15-0.25 and a second-pass loss ~30 %; Wolf 1977
  sees 0. N + pool (p 0.5) with k 0.57 predicts pass 1 to 0.66 (obs 0.52) and a ~17 % second-pass loss. Neither fits;
  the 2-min lag-corrected times are not printed.

**What would decide CF-072** (in order of cost):
1. Nguyen 2016 Fig. S1B lysine mass balance (already the remedy's decider; on disk as SI).
2. A sugar-chase on plateau samples, one lab: heat plateau material after removing free glucose (re-slurry, dialyse or
   ultrafilter, re-dry to the same a_w). R predicts available AND total free lysine (LC-MS/MS after enzymatic hydrolysis,
   not acid: acid regenerates lysine from fructoselysine) rising at ~k_flp_decay (>= +0.2 dex in 60 min at 105 C);
   N (+pool) predicts no rise. Thompson 1976's "phase two" rise (98 % confidence; cited, not on disk) is the only in-lab
   hint; download Thompson, Wolf & Reineccius 1976, Trans. ASAE (DOI unknown; check Crossref before listing).
3. A glucose step at the plateau (add glucose): R predicts a deeper new plateau (k_d / (k_d + k_glyc G')); a fixed pool
   predicts none.
4. Furosine (FLP) alongside FDNB over the plateau: R predicts a steady FLP pool (~(1 - f)·L0) with continuing 3-DG
   production; N + pool predicts FLP peaking then falling.

### Proposed register text (do NOT apply without owner review)

Append to CF-072 `tests:`
```yaml
      - "wolf1977a (J Food Sci 42:1540; Crossref-verified 2026-11-01), FDNB lysine in soy / glucose / cellulose at a_w 0.65, 95-125 C: does NOT discriminate. Shipped constants: return structure leaves 0.32-0.42 at 105 C / 54 min / 4 % glucose, no-return 0.05-0.14, observed 0.43 (+0.5 to +0.8 dex between rivals) -- but with k_glyc free (a no-return refit needs it 0.3-0.5 dex below shipped, inside the between-lab band) both structures fit Table 2 at rms 0.06 dex, and Table 3 (95 C, 1 %), the order test and the 125 C design end points separate them by <= 0.07 dex. The plateau with 2.5 % glucose left rejects a BARE no-return network (P 0.8) but not no-return + an inaccessible lysine pool (the engine's available_fraction) or + FDNB counting Amadori lysine. LR return : no-return+pool ~1; return : bare no-return ~2. P(Wolf 1977a decides CF-072) 0.1. Scratch: wolf1977a_findings.md."
      - "Decider (proposed): a sugar-chase on plateau material -- free glucose removed, re-dried to the same a_w, reheated; the return structure predicts available and LC-MS total-free lysine rising at ~k_flp_decay (>= +0.2 dex in 60 min at 105 C), the no-return structure none. Cheaper first: Nguyen 2016 Fig. S1B lysine balance."
```
Status unchanged (`error_confirmed`: the code-vs-Nguyen transcription mismatch stands regardless of the chemistry).

New conflict proposed (rival-independent; for the owner to register or merge into F05):
```yaml
  - id: CF-XXX
    kind: magnitude
    scope: corpus
    sources: [wolf1977a]
    quantity: "protein-glycation barrier and glucose order: k_glyc (Ea 96.8 kJ/mol, first order in glucose) vs FDNB lysine loss in soy at a_w 0.65"
    observation: >-
      Wolf 1977a's response surface (Case 1, pH 7, a_w 0.65) gives E 132-138 kJ/mol at 2-4 % glucose and a k ratio of only
      2.3 (80 C) to 6.6 (125 C) for an 11.5x glucose step. The engine's k_glyc·[G]_water (moisture 8-12 g/100 g, not
      printed) matches Wolf at 125 C (+-0.1 dex) but over-predicts by +0.24 to +0.42 dex at 100 C and +0.5 to +0.7 dex at
      80 C. Rival-independent (both CF-072 structures share the initial slope).
    explanations:
      - "Plateau-capped 80 C / 360 min end points (k = -ln(C/C0)/t taken after the loss stopped) inflate Wolf's E and flatten the glucose slope (P ~0.5)."
      - "A real barrier difference between free glycine (Martins 2005, the engine's source) and protein-bound lysine in a low-moisture solid."
      - "Water-phase glucose concentration not proportional to dry-basis glucose (dissolution / glass / moisture shifts at 125 C)."
    tests:
      - "Any soy or pea glycation TIME COURSE at >= 3 temperatures, 80-130 C, with moisture printed and the plateau excluded from the rate fit."
    status: open
    affects:
      - {file: src/kinetic_core/parameters_glycation.py, symbol: EA_GLYC_KJ_MOL}
```

## 2. Disagreements with existing dossiers / ledger

1. `wolf1977_extraction.md` sec. 3 Table 1 row 5 and ledger `wolf1977_tab1_tab2_soy_available_lysine_plateau_reprocessing`
   (note: "Table 1 row 5 is from Wolf 1977a (secondary_source)"): **Wolf 1977a prints no 125 C / 2 min / 4 % glucose value
   (5.0 +- 0.7)**. Its tables are 105 C (Table 2) and 95 C (Table 3); the 125 C design raw data are not printed. The
   citation in the companion cannot be verified against 1977a; flag `secondary_source_unverified` (likely Wolf 1978 or
   unpublished).
2. `wolf1977_extraction.md` sec. 1 answer 1 / sec. 7: "The kinetic constants ... are in the companion papers": true only
   as regression surfaces (no per-condition k table); the companion's >= 0.33 /min bound at 125 C agrees with 1977a's
   Case 1 0.30 /min (0.04 dex) — consistent.
3. Product and composition: Promine D "soybean protein concentrate" (1977a) vs "soy protein isolate" (companion dossier and
   ledger `system`); 15 % vs 20 % protein. Ledger `system` field for the companion should not be used to describe 1977a.
4. Lysine per protein: 6.3 (1977a Table 2), 4.8 (Table 3), ~7.3 implied by the order test, 4.8-5.1 (companion) g/100 g
   protein. Batch spread ~0.15 dex in one lab and assay: the floor for any within-lab comparison across batches.
5. `wolf1977_extraction.md` sec. 4: "sugar exhaustion and local sugar depletion (by Wolf 1977a / b, [C])" — 1977a excludes
   bulk exhaustion at 4 % only (2.5 % left at 30 min). It does not test local depletion, and at 1 % glucose (Table 3)
   exhaustion is not excluded (glucose unprinted; ~4 mol glucose lost per mol lysine at 4 %).
6. Engine vs Wolf (new, rival-independent): Ea 96.8 vs 132-138 kJ/mol; glucose order 1 vs 0.35-0.8 (sec. 1 above).
7. Internal to 1977a: Case 3 does not reproduce Cases 1-2 (centre E 19.9 vs 32.9 kcal/mol); the text's 5.3x / 25x glucose
   effects are not reproduced by any printed case; E is printed as "Kcal/mol·K" and "Cal/mole".

## 3. Lane notes

- E5 (buffer catalysis): Wolf's pH arms switch buffer species (acetate 5.32; phosphate 6.0, 7.0; borate 8.0, 8.68) at
  0.1 M x 1.5 mL / g = ~0.15 mmol buffer / g dry, i.e. ~1-2 M in the water phase at a_w 0.65 [Z, moisture assumed]. Do not
  read Figs. 4-5 as a pH law; borate complexes glucose.
- a_w: a_w acts only as a_w x glucose (Case 2 k_100 term 0.00104·G·a_w; Fig. 6 crossing at a_w ~0.43); no maximum near
  0.65. Supports a weak sugar-dependent a_w term for F05 (d), not a bell.

## 4. Proposed ledger entries (aroma_other.yml, beside the companion; NOT applied)

```yaml
  - id: wolf1977a_tab2_tab3_soy_available_lysine_plateau_time_course
    dossier: wolf1977a_extraction.md
    paper: "Wolf, Thompson & Reineccius 1977 J Food Sci 42:1540"
    lab: "Minnesota (Reineccius)"
    system: "15 % soy protein (Promine D, called a concentrate) + microcrystalline cellulose + glucose 4 % (Table 2) or 1 % (Table 3); Table 2 arm re-slurried at 30 min, freeze-dried, re-equilibrated, reheated"
    medium: "slurried in 0.1 M sodium phosphate pH 7.0 (1 : 1.5 w/v), freeze-dried, a_w 0.65 over saturated salt (hygrometer-checked), sealed retort pouches (<= 0.07 g / cm2)"
    temperature_c: [105, 95]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 15, 30, 42, 45, 54, 60, 75]
    species: ["available lysine (FDNB, mg/g sample)", "glucose (% , one value: 4 -> 2.5 at 30 min, 105 C)"]
    kind: time_series
    location: "Table 2 (105 C; continuous 0 / 30 / 42 / 54 min, reslurry 42 / 54) and Table 3 (95 C, 0-75 min), printed mean +- SD; glucose in text p. 1542; dossier sec. 3.5-3.6"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "glucose levels were set at 1%, water activity at 0.65 and pH at 7.0"
    flags: [matrix_food, low_moisture_no_aw_term, ph_start_only, "condition_not_stated: moisture content", "condition_not_stated: glucose in the 1 % series", "condition_not_stated: n per time point", internal_inconsistency, "also_relevant: matrix"]
    note: "Phase three with sugar in excess: 105 C, 4 % glucose 9.4 -> 6.5 -> 4.0 -> 4.0 mg/g (0 / 30 / 42 / 54 min), re-slurry 3.6 / 4.9 (n.s.); glucose 4 -> 2.5 % at 30 min (not exhausted). 95 C, 1 % glucose 7.2 -> 6.7 -> 4.9 -> 5.4 -> 5.0 -> 4.9 (plateau ~0.69). C0 differs 0.12 dex between the series at the same protein level. Simulable only with a declared moisture band. Does not decide CF-072 (findings 2026-11-01). arms: 4 vs 1 % glucose (also 105 vs 95 C), continuous vs re-slurried (same lab, same FDNB assay)."
  - id: wolf1977a_order_test_105c_two_protein_levels
    dossier: wolf1977a_extraction.md
    paper: "Wolf, Thompson & Reineccius 1977 J Food Sci 42:1540"
    lab: "Minnesota (Reineccius)"
    system: "7.5 or 15 % soy protein (Promine D) + 4 % glucose + microcrystalline cellulose"
    medium: "0.1 M phosphate pH 7.0 slurry, freeze-dried, a_w 0.65, sealed pouches"
    temperature_c: [105]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [18]
    species: ["available lysine (FDNB)"]
    kind: rate_constant
    location: "p. 1541-1542 text (losses 2.3 / 4.1 mg/g; first-order k 0.024 / 0.026 min-1); dossier sec. 3.2"
    status: text_only
    quant_class: calibrated
    simulable: false
    condition_quote: "The specific reaction rates for the first order curves were 0.024 and 0.026"
    flags: [ph_start_only, "condition_not_stated: C0", "condition_not_stated: moisture content", low_moisture_no_aw_term, matrix_food]
    note: "First order in lysine (equal slopes at two protein levels; zero / second order rejected at 99 %). Two-point (18 min + controls) slopes; C0 not printed (0.026 over 18 min implies ~11 mg/g at 15 % protein, above Tables 2-3)."
  - id: wolf1977a_response_surface_e_k100_soy_glycation
    dossier: wolf1977a_extraction.md
    paper: "Wolf, Thompson & Reineccius 1977 J Food Sci 42:1540"
    lab: "Minnesota (Reineccius)"
    system: "15 % soy protein (Promine D) + glucose 0.32-3.68 % + microcrystalline cellulose; central composite pH x a_w x glucose, two replicates, two temperatures"
    medium: "0.1 M buffer slurries (acetate pH 5.32; phosphate 6.0, 7.0; boric acid 8.0, 8.68), freeze-dried, a_w 0.33-0.94 over saturated salts (hygrometer-checked), sealed pouches; 80 C / 360 min and 125 C / 1.5 min (lag-corrected)"
    temperature_c: [80, 125]
    ph: {initial: [5.32, 6.0, 7.0, 8.0, 8.68], measured_during: false, final: null}
    time_min: [1.5, 360]
    species: ["available lysine (FDNB): first-order E and k_100"]
    kind: rate_constant
    location: "p. 1542 response-surface equations Cases 1-2 (Case 3 unusable); Table 4 regression summary (p. 1543); Figs. 3-6 drawn from Case 2; dossier sec. 3.3-3.4"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "The two processing temperatures were 125°C and 80°C."
    flags: [derived_by_dossier, ph_start_only, "condition_not_stated: moisture content", "condition_not_stated: composition-specific C0 (center-point controls only)", buffer_species_confounded_with_ph, endpoint_derived_rate_constant, internal_inconsistency, matrix_food]
    note: "Case 1: E = 32.9 + 1.94 G - 1.64 aw G - 1.32 G^2 + 1.50 pH G kcal/mol; k_100 = 0.0124 + 0.0043 G min-1 (coded: G = glucose% - 2, pH - 7, (aw - 0.65)/0.15). Centre E 137.7 kJ/mol; residual SD k_100 0.003 min-1 (0.1 dex). Two-temperature Arrhenius per composition, endpoint k = -ln(C/C0)/t; 5 / 96 points discarded as outliers. Held-out TEST and Ea band for protein glycation in a low-moisture plant matrix, NOT a fit target. Engine k_glyc·G agrees at 125 C, over-predicts +0.3 to +0.7 dex at 80-100 C. Text's 5.3x / 25x glucose effects not reproduced by the printed equations."
```

Companion ledger entry `wolf1977_tab1_tab2_soy_available_lysine_plateau_reprocessing`: proposed flag change — replace
`secondary_source` with `"secondary_source_unverified: Table 1 row 5 cited to Wolf 1977a, not printed there"`.
