"""
src/barrier_constants.py — Centralised FAST-mode heuristic barrier constants.

These values are FITTED calibration constants (kcal/mol) for each Maillard
reaction family: literature ranges were used as starting points, but many
entries were subsequently adjusted so the FAST ranking reproduces the
benchmark panel (several inline comments below are calibration rationales,
e.g. "must stay competitive enough to recover the ... sulfur balance").
Treat them as model parameters with literature-informed priors — NOT as
independently anchored literature barriers.

They are used by both `pipeline.py` and `run_pipeline.py` for the
instant FAST-mode rankings.  Update this single file when new data
is available — both call-sites import from here.

Starting-point sources (see audit caveats)
------------------------------------------
* Yaylayan & Huyghues-Despointes 1994 — AUDIT 2026-08-26: that review covers
  Amadori rearrangement products; cited here for Schiff ΔG‡, off-topic.
* Martins & van Boekel 2003 (Amadori kinetics)
* Hofmann & Schieberle 2000 (Strecker degradation)
* Wedzicha 1984 — AUDIT 2026-08-26: paper is a kinetic model of the
  SULPHITE-inhibited Maillard reaction, not cysteine thermolysis.
* Hodge 1953; Nursten 2005 (overall Maillard kinetics)
* data/Gemini_Deep_Research/maillard_meat.md, data/Gemini_Deep_Research/maillard_plant_based.md (project literature reviews)

AUDIT 2026-08-26 (forensics) — the sulfur-branch values (thiol_addition,
thiol_addition_hexose, thiol_oxidation, aminoketone_condensation,
strecker_degradation, thiohemiacetal_formation, thiol_dehydration) were all set
in commit 2ea7d12 by a joint fit whose targets included the now-quarantined
Mottram1994/Farmer1999 benchmarks. Post-quarantine, only thiol_addition retains
a literature constraint (Hofmann1998, admissible window [28.10, 28.85] kcal/mol);
it was RE-CENTRED to the Hofmann-only optimum 28.60 on 2026-08-26 (owner-approved).
Treat the rest of the sulfur branch as free parameters pending re-anchoring — see
tasks/audit_remediation.md.

CORRECTION 2026-08-27 (Wave S2c) TO THE PARAGRAPH ABOVE — the sentence "only
thiol_addition retains a literature constraint (Hofmann1998 …)" IS FALSE, and the
paragraph is kept verbatim because it records what the repo believed, not because it
is true.  Wave S2b settled the provenance of that "constraint": cys_ribose_140C_Hofmann1998's
MFT 342 ppb and FFT 200 ppb are NOT from 10.1021/jf9705983 or from any other paper.
They are interior points of two invented bands (MFT `~0.02–0.05` mol %, FFT
`~0.01–0.03` mol %) in data/benchmarks/maillard_validation_benchmarks.md §1.3, an
abstract-reconstructed table committed in c7efbbc — the same commit that created the
benchmark JSON.  On the file's declared (and unattested) 10 mM basis with MW 114.17 the
arithmetic closes exactly: 0.0300 mol % → 342.5 → 342 ppb; the geometric mean of the FFT
band, 0.017321 mol % → 197.7 → 200 ppb.  Confidence ~90%.
**THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.**  Every value in it —
thiol_addition included — is an estimate or a fit to a repo-internal number.  The
[28.10, 28.85] "Hofmann window" is a window around the repo's own guess.

UPDATED 2026-08-28 (Wave W) — READ BOTH HALVES, THEY SAY DIFFERENT THINGS NOW.
The BENCHMARK half of that paragraph is superseded: the sulfur branch has three absolute,
stable-isotope-dilution literature anchors, `hofmann1998_{ribose,glucose,fructose}_cysteine
_145C_20min_pH5`, taken from Table 1 of the real 10.1021/jf9705983 after the full text was
obtained.  The paper also confirms from the primary source that 342 and 200 ppb appear
nowhere in it, so everything above about the provenance of the retired values stands.
The CONSTANTS half is UNCHANGED and is the part that matters here: **no barrier in this
file is anchored to a literature measurement, and Wave W fitted nothing.**  The new anchors
were deliberately not used to refit anything — the model misses them by 12.27x, 29.58x and
14.46x, and closing that by solving four free barriers against six rows would recreate
exactly the circularity Wave H produced and Wave S2c had to retract.  What changed is that
these estimates can now be SCORED against real measurements instead of being unfalsifiable.
They score badly.  The [28.10, 28.85] "Hofmann window" is still a window around the repo's
own guess and is still not a literature window.

AUDIT 2026-08-27 (Wave H, sulfur refit) — Wave G1 replaced the fabricated one-step
MFT shortcut with the accepted 1-deoxyosone → norfuraneol → MFT route, which moved
the Hofmann constraint OFF `thiol_addition` (that key now labels only the demoted
shortcut and the lumped `Thiol_Addition_H2` step) and ONTO
`thiol_addition_norfuraneol`.  That constant was refit against Hofmann1998 alone
(28.60 → 26.85; see its comment).  Measured at the same time, and recorded so the
next person does not re-run the search: `thiohemiacetal_formation` has EXACTLY ZERO
derivative on this benchmark over its whole defensible range; `thiol_dehydration` buys
under 0.01 dex; and `furanone_cyclisation`, which only becomes identifiable once
`thiol_addition_norfuraneol` drops below it, turns out to sit exactly at its own optimum
(achievable gain 0.0000 dex). All three keep their incumbents.  The refit's own headline
is that it does not work: no barrier value in any defensible range gets MFT closer
than ~5.6x under, because the deficit is in the volatile-budget ALLOCATION (furfural,
unmeasured in this benchmark, takes ~78% of a total budget that is itself the right
order of magnitude), not in the barriers.  Full profile:
results/validation/sulfur_barrier_refit_hofmann.md.

AUDIT 2026-08-27 (Wave I fix 8) — THE PARAGRAPH ABOVE IS NOW STALE, and is kept
verbatim because it records what was true when the refit ran, not because it still
describes the model.  Wave I fixed red-team finding H4: the `[HH]` reducing-equivalent
token that gates `Thiol_Addition_Norfuraneol` had exactly one producer reachable from
a cysteine/sugar system, the pyrazine aromatisation, so the whole sulfur branch was
running behind a pyrazine bottleneck it has no business being behind.  Giving the token
its own source (`2 cysteine -> cystine + 2[H]`, src/reaction_templates.py
`_thiol_reductant_pool`) moved Hofmann1998 MFT from 61.25 ppb to 345.04 ppb against a
342 ppb reference (5.58x UNDER -> 1.01x) and FFT from 61.44 to 187.49 ppb against 200.
NO BARRIER WAS TOUCHED: `thiol_addition_norfuraneol` is still the 26.85 that Wave H
fitted.  Consequences the next owner must not miss:
  * the "5.6x under, deficit is in the allocation" diagnosis was an artefact of the
    H4 coupling, so the Wave H profile in results/validation/sulfur_barrier_refit_hofmann.md
    was computed on a network that no longer exists and must be RE-RUN before anyone
    cites it again;
  * 26.85 was selected as the conservative edge of an indifference band on a saturated
    profile.  On the fixed network the profile is no longer saturated, so this value is
    now an UNREFIT constant sitting near a benchmark it happens to reproduce.  That
    agreement is a coincidence of two independent choices, NOT a calibration — do not
    report it as one.  OPEN OWNER ITEM: re-run scripts/generators/refit_sulfur_barriers_hofmann.py.

AUDIT 2026-08-26 — cross-file inconsistency: data/lit/arrhenius_params.yml
disagrees with this table on several families while both claim anchoring
(mutarotation 22.6 vs 5.0 kcal/mol; thiol_addition 29 kJ/mol there vs
28.85 kcal/mol here — same numerals, different units, suspected kJ/kcal
transcription collision). Values left UNCHANGED pending human re-anchoring;
see tasks/audit_remediation.md.
"""

import json
import yaml
import math
from pathlib import Path
from typing import Any, Dict, Tuple, Optional, Sequence

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
ARRHENIUS_FILE = ROOT / "data" / "lit" / "arrhenius_params.yml"
REFINEMENT_PATCH_FILE = ROOT / "data" / "lit" / "refinement_surrogate_patches.json"
COMPUTATIONAL_PRIORS_FILE = ROOT / "data" / "lit" / "computational_priors.json"

# Exact Mapping: normalized reaction family name → barrier in kcal/mol.
# Replaces fragile substring matching.
FAST_BARRIERS: Dict[str, Tuple[float, str]] = {
    # ── Sugar prerequisite ──────────────────────────────────────────
    "mutarotation":         ( 5.0,  "Ring opening is near-barrierless (hemiacetal ⇌ open-chain)"),
    "ring_opening":         ( 5.0,  "Ring opening is near-barrierless (hemiacetal ⇌ open-chain)"),

    # ── Core Maillard cascade ───────────────────────────────────────
    # AUDIT 2026-08-27 (Wave I fix 10, red-team H5) — SCHIFF/AMADORI ORDERING
    # INVERSION between this table and data/lit/arrhenius_params.yml. BOTH
    # VALUES UNCHANGED; this is an annotation, not a fix. At 150 C:
    #     HERE (screening lane):  schiff 15.0 kcal -> t1/2 2.59e-4 s
    #                             amadori 23.0 kcal -> t1/2 5.26 s
    #                             => Schiff ~2.0e4 x FASTER
    #     YAML (Cantera lane):    schiff A=1.5e11, 97.0 kJ -> t1/2 4.36 s
    #                             amadori A=1.0e11, 59.0 kJ -> t1/2 1.33e-4 s
    #                             => Amadori ~3.3e4 x FASTER
    # The two tables disagree about which of the first two steps of the entire
    # cascade is rate-determining, and disagree by ~6.6e8 in the ratio.
    # AUTHORITY, stated plainly so nobody has to guess:
    #   * THIS TABLE (FAST_BARRIERS) is authoritative for the SCREENING /
    #     recommend lane — src/recommend.py, src/pipeline.py, and every
    #     benchmark prediction through src/benchmark_validation.py. These are
    #     FITTED calibration constants, as the module docstring says.
    #   * data/lit/arrhenius_params.yml is authoritative for the CANTERA EXPORT
    #     lane — get_arrhenius_params -> src/cantera_export.py ->
    #     scripts/run_cantera_kinetics.py. Those claim literature anchoring.
    # No code path reads both, which is how the inversion survived.
    #
    # ---------------------------------------------------------------------
    # RESOLVED 2026-08-27 (Wave S3) — THE DATA DECIDED. THIS TABLE HAS THE
    # ORDERING BACKWARDS. VALUES STILL UNCHANGED; see below for why.
    # ---------------------------------------------------------------------
    # Wave S3 fitted the trunk to Martins' measured glucose/glycine
    # trajectories at 80/100/120 C — which measure the Amadori compound (DFG)
    # DIRECTLY — plus an experiment on ISOLATED DFG. Both are in
    # data/lit/timeseries/; the fit is
    # scripts/generators/generate_trunk_rate_calibration.py and the result is
    # results/validation/trunk_rate_calibration_refit.{json,md}.
    #
    # The two steps can only be compared as PSEUDO-FIRST-ORDER rates, because
    # the condensation is bimolecular and the rearrangement is not. At the
    # experiment's own 200 mmol/L glycine and 100 C:
    #     k_schiff * [Gly]  = 3.39e-3 /min      (condensation)
    #     k_amadori         = 1.51e-1 /min      (rearrangement)
    #     => AMADORI IS ~45x FASTER. The condensation is rate-determining.
    # A profile likelihood over that ratio rejects EVERY value from 0.1 to
    # 1e6 outside a narrow band around 45, at 95% confidence.
    #
    # SCORECARD:
    #   * THIS TABLE claims the ratio is ~1/2.0e4 (Schiff faster).
    #     WRONG SIGN. Rejected with a chi-square statistic of ~7.5e4.
    #   * THE YAML claims ~3.3e4 (Amadori faster).
    #     RIGHT SIGN, wrong magnitude by ~700x. Also rejected (~9.1e2).
    # So the reconciliation is not "pick one file": the YAML has the ordering
    # right and the size wrong, and this table has the ordering wrong.
    #
    # WHY THE NUMBERS BELOW STILL DID NOT MOVE. Converting the fitted rates
    # into this table's kcal/mol units requires assuming the screening lane's
    # prefactor (1e11 s^-1, hardcoded family-agnostically at
    # recommend.py:1216) is the true prefactor for each step. Every decade of
    # error in that assumption moves the derived barrier by 1.71 kcal/mol at
    # 100 C — wider than any confidence interval in the fit. The derived
    # barriers are published in full, with the arithmetic, in the calibration
    # artifact; `amadori_rearrangement` derives to 23.22 kcal/mol against the
    # 23.00 shipped here, which is the first independent confirmation of any
    # value in this table from a measured trajectory. `schiff_condensation`
    # is EXCLUDED from that mapping entirely: it is bimolecular, and there is
    # no dimensionally valid conversion without inventing a standard-state
    # concentration. OPEN OWNER ITEM: apply the derived set, or keep the
    # two-lane split with the ordering now on the record.
    "schiff_condensation":  (15.0,  "Yaylayan 1994: condensation ΔG‡ ≈ 12–20 kcal/mol; midpoint"),
    "schiff_base_hydrolysis":(8.0,  "Schiff base reversion; fast"),
    "amadori_rearrangement":(23.0,  "Martins 2003: 1,2-proton shift ΔG‡ ≈ 20–28; midpoint"),
    "heyns_rearrangement":  (24.0,  "Ketose analogue of Amadori; slightly higher barrier"),

    # ── Sulfur pathways ─────────────────────────────────────────────
    "cysteine_thermolysis": (30.0,  "Wedzicha 1984: thermolysis ΔG‡ ≈ 20–30; upper range"),
    "thiol_addition_trimolecular": (24.0, "H2S-mediated sulfur trapping should remain accessible but no longer tie the upstream carbonyl bottleneck"),
    "thiohemiacetal_formation": (23.3, "Furfural-thiohemiacetal formation is favorable but not faster than the dominant carbonyl cascade"),
    "thiol_dehydration":    (26.8, "Thiohemiacetal dehydration remains feasible but should impose a real selectivity cost relative to direct furfural release"),
    "thiol_addition":       (28.60,  "Re-centred 2026-08-26 to the Hofmann1998-only optimum after the Mottram/Farmer quarantine: the joint fit had pushed the value to 28.85, the boundary of Hofmann's admissible window [28.10, 28.85]; 28.60 minimises Hofmann MALE (0.082 -> 0.051)"),
    "thiol_addition_hexose": (29.65, "UNCONSTRAINED LEGACY FIT (2026-08-26 forensics): tuned solely to the quarantined Farmer1999 benchmark; no surviving literature observable constrains it. Kept for behaviour continuity pending re-anchoring"),
    "thiol_oxidation":      (29.02,  "UNCONSTRAINED LEGACY FIT (2026-08-26 forensics): tuned solely to the quarantined Mottram1994 disulfide row; surviving comparators are synthetic-only. Kept for behaviour continuity pending re-anchoring"),

    # ── Enolisation / dehydration ───────────────────────────────────
    "enolisation_intermediate": (21.0,  "Amadori/Heyns deoxyosone formation is the common gateway into furfural and sulfur branches; keep it competitive instead of falling back to the heuristic default barrier"),
    "1,2-enolisation":      (28.0,  "Literature replication calibration: furfural-forming dehydration should be more competitive in benchmark systems"),
    "2,3-enolisation":      (28.0,  "Nursten 2005: enolisation is rate-limiting in advanced Maillard"),
    "dehydration":          (28.0,  "Coupled with enolisation; same approximate range"),

    # ── Strecker cascade ────────────────────────────────────────────
    "strecker_degradation": (22.0,  "UNCONSTRAINED LEGACY FIT (2026-08-26 forensics): calibrated against the quarantined Farmer1999 pyrazine row; no surviving observable constrains it. Kept pending re-anchoring"),
    "aminoketone_condensation": (29.0,  "UNCONSTRAINED LEGACY FIT (2026-08-26 forensics): tuned solely to the quarantined Farmer1999 pyrazine target; no surviving observable constrains it. Kept pending re-anchoring"),

    # ── Retro-aldol ─────────────────────────────────────────────────
    "retro_aldol":          (32.0,  "Hodge 1953: C-C bond cleavage is high-barrier; softened from 35"),
    "lipid_thiazole":       (20.0,  "Lipid thiazole condensation; comparable to Strecker"),

    # ── DHA / β-elimination ─────────────────────────────────────────
    "beta_elimination":     (18.0,  "β-elimination of Ser/Cys; moderate barrier"),
    "dha_crosslinking":     (18.0,  "DHA crosslinking with Lys; same range as β-elim"),

    # ── Lipid trapping & Synergy ──────────────────────────────────────
    "lipid_condensation":   (14.0,  "Lipid aldehyde Schiff base trapping; fast condensation"),
    "lipid_strecker_synergy": (18.0,  "Lipid-Maillard synergy (alkylthiazole/pyrazine) is highly favourable"),

    # ── Thermal / additive degradation ──────────────────────────────
    "thiamine_degradation": (25.0,  "Thiamine thermolysis; moderate barrier"),
    "additive_degradation": (25.0,  "Generic additive degradation"),
    "glutathione_cleavage": (22.0,  "GSH peptide bond cleavage"),

    # ── Lipid Oxidation (Phase 19) ──────────────────────────────────
    # AUDIT 2026-08-27 (Wave I fix 10, red-team H5) — KINETICALLY DEAD. VALUE
    # UNCHANGED, arithmetic stated. 42.0 kcal/mol evaluated through this
    # module's own `arrhenius_rate_constant` (which uses
    # DEFAULT_REFERENCE_PREEXPONENTIAL = 1e11 for this family, since
    # `lipid_homolysis` has no entry in data/lit/arrhenius_params.yml):
    #     k(100 C) = 2.51e-14 1/s -> t1/2 = 8.8e5 years
    #     k(140 C) = 6.04e-12 1/s -> t1/2 = 3.6e3 years
    #     k(150 C) = 2.02e-11 1/s -> t1/2 = 1.09e3 years
    # Cooking happens in minutes. This step therefore never fires, which means
    # the entire radical chain it is supposed to initiate is switched off at the
    # source — the same failure mode as the eight families that were falling
    # through to DEFAULT_BARRIER = 45.0 before Wave G1 fix 8, except this one is
    # explicit and so looked deliberate. Even with a physically proper O-O
    # homolysis prefactor (1e15-1e16 1/s rather than 1e11) the half-life at
    # 150 C is ~40 days, still dead. The real defect is that 42 kcal/mol is the
    # GAS-PHASE unimolecular O-O bond dissociation enthalpy; in a food matrix
    # hydroperoxide decomposition is METAL-CATALYSED (Fe(II)/Fe(III) redox
    # cycling), with apparent barriers in the ~20-25 kcal/mol range — which is
    # why src/lipid_oxidation.py has to model hexanal through a separate
    # empirical Arrhenius path (Ea = 80 kJ/mol = 19.1 kcal/mol) instead of
    # through this network family at all. Reconciling the two lipid lanes is an
    # OPEN OWNER ITEM; Wave I changed no value.
    "lipid_homolysis":      (42.0,  "O-O bond cleavage in hydroperoxides; high barrier"),
    "beta_scission":        (22.0,  "β-scission of alkoxy radicals; moderate barrier"),
    "radical_crosstalk":    (15.0,  "Radical + H2S collisions; fast. NOTE 2026-08-27: the Radical_Crosstalk FAMILY was retired from the engine (it existed to consume fictitious elemental sulfur and quenched radicals by eating MFT). This key is retained because src/family_barrier_progress.py uses it as a family-coverage key"),

    # ── ESTIMATED TIER (added 2026-08-27, Wave G1 fix 8) ─────────────
    # Eight families were emitted by the engine but had NO entry here and no
    # `_canonical_fast_family` route to one, so every one of them fell through
    # to DEFAULT_BARRIER = 45.0 kcal/mol — a half-life of ~39,000 years at
    # 150 C.  They were, in effect, silently switched off: the acrylamide and
    # CML/CEL safety lanes, the furanone lane, the thiamine/GSH additive lane
    # and the entire radical chain.  Giving them explicit values ACTIVATES
    # those lanes; every value below is an ESTIMATE taken from the closest
    # literature-anchored analogue already in this table, named in its comment,
    # and NONE was chosen to reproduce any previous output.
    "safety_risk_acrylamide": (30.83, "ESTIMATED. Direct transfer of the Knol et al. 2009 lumped Asn+sugar->acrylamide anchor Ea=129 kJ/mol (10.1016/j.foodchem.2009.11.049) already carried in DFT_ANCHOR_METADATA['asparagine_sugar_explicit_water_cluster']; that anchor describes exactly this step"),
    "safety_risk_age":        (23.0,  "ESTIMATED. CML/CEL formation is a carbonyl-amine addition on the lysine epsilon-amine followed by rearrangement; taken as the Amadori analogue `amadori_rearrangement` (23.0, Martins 2003), the same 1,2-addition/proton-shift class"),
    "furanone_formation":     (28.0,  "ESTIMATED. Furanone closure is a cyclisation/dehydration of a deoxyosone; taken as the `dehydration` / `2,3-enolisation` analogue (28.0, Nursten 2005)"),
    "furanone_cyclisation":   (28.0,  "ESTIMATED. 1-deoxyosone -> norfuraneol/DMHF cyclodehydration; same analogue as `furanone_formation` (`dehydration`, 28.0). Deliberately NOT tuned to reproduce the output of the retired one-step MFT shortcut"),
    "thiol_addition_norfuraneol": (26.85, "RETIRED 2026-08-27 (Wave N): no step emits this family any more -- the norfuraneol + H2S -> MFT step was removed from `_furanone_and_mft_route` on isotope evidence (Cerny & Davidek 2003, 10.1021/jf026123f: norfuraneol spiked into a [13C5]ribose/cysteine system yields mainly 13C5-labelled MFT, so norfuraneol is 'unimportant as an intermediate'; Hofmann & Schieberle 1998, 10.1021/jf9705983, independently rank norfuraneol/cysteine as LESS effective MFT precursors). Entry kept for provenance only. ANNOTATION 2026-08-27 (Wave S2c): the Wave H fit recorded below was NEVER a literature constraint, and neither this entry nor any document may describe it as one. Wave S2b established that cys_ribose_140C_Hofmann1998's MFT 342 / FFT 200 ppb are interior points of two invented bands in data/benchmarks/maillard_validation_benchmarks.md section 1.3 -- a repo-internal derivation committed in the same commit as the benchmark, not a value from 10.1021/jf9705983. So 26.85 was fitted through a contradicted ROUTE (Wave N's finding, below) AND against a non-measurement (Wave S2b's finding). The value is left at 26.85 only because NO STEP EMITS THIS FAMILY: reverting a barrier that reaches nothing would move no prediction and would only make the provenance record harder to read. Its history: FITTED 2026-08-27 (Wave H) against cys_ribose_140C_Hofmann1998 ONLY. Norfuraneol + H2S -> MFT (van den Ouweland & Peer 1975, 10.1021/jf60199a045 -- a genuine SYNTHESIS route, not the in-situ one). Was 28.60, inherited from `thiol_addition` and that key's PRE-G1 Hofmann window [28.10, 28.85]. Refit over [23.30, 29.65]: profile min 0.6198 dex, incumbent 0.6987 dex; 26.85 was the conservative edge of the indifference band. The 26.85 was therefore fitted THROUGH a contradicted route and does NOT transfer to `thiol_addition_pentodiulose`. Record: results/validation/sulfur_barrier_refit_hofmann.{json,md}"),
    "thiol_addition_pentodiulose": (28.60, "ESTIMATED. HISTORY, in order, because the middle step must not disappear: (1) 28.60 ESTIMATED 2026-08-27 (Wave N) as the un-fitted `thiol_addition` CLASS value; (2) 26.35 FITTED 2026-08-27 (Wave P item 1) against cys_ribose_140C_Hofmann1998 ONLY; (3) 28.60 REVERTED 2026-08-27 (Wave S2c) back to the Wave N class value, and the fit record results/validation/sulfur_barrier_refit_pentodiulose.{json,md} RETRACTED. WHY THE REVERT: Wave S2b (tasks/audit_remediation.md '## Wave S2b') established that the refit's SOLE fit target, the 342 / 200 ppb in cys_ribose_140C_Hofmann1998, is NOT a measurement from 10.1021/jf9705983 or from any other paper. Both values are interior points of two invented bands in data/benchmarks/maillard_validation_benchmarks.md section 1.3 (MFT `~0.02-0.05` mol %, FFT `~0.01-0.03` mol %), an abstract-reconstructed table committed in c7efbbc -- the SAME commit that created the benchmark file. The arithmetic closes exactly on the file's declared 10 mM basis with MW 114.17: 0.0300 mol % -> 342.5 -> 342 ppb, and the geometric mean of the FFT band, 0.017321 mol % -> 197.7 -> 200 ppb. Confidence ~90%. A constant fitted to the repo's own guess is not calibrated; it is circular, and the same treatment was applied in Wave I when the Methional hydrolysate observability constant was reverted after its only two fit targets turned out to be fabricated. THE SULFUR BRANCH THEREFORE HAS ZERO ABSOLUTE LITERATURE ANCHORS, and this constant is honestly an ESTIMATE again rather than dishonestly a fit. THE MECHANISM IS UNCHANGED BY THE REVERT -- only the number's provenance is. WHAT THE RETIRED FIT SAID, kept verbatim for the record: 1,4-dideoxypento-2,3-diulose + H2S -> MFT + 2 H2O, the isotope-evidenced in-situ MFT step (Cerny & Davidek 2003, 10.1021/jf026123f, proposed; Cerny & Davidek 2004, 10.1021/jf035265m, positionally confirmed with [1-13C]ribose). Wave N shipped it at the un-fitted `thiol_addition` class value 28.60, deliberately NOT the 26.85 that Wave H had fitted through the RETIRED norfuraneol route. This refit is owner-approved and was run LAST in the wave, after the chemistry additions, so it fits the network that actually ships. Range [23.30, 29.65] (thiohemiacetal_formation .. thiol_addition_hexose, both already in this table); ONE free parameter, decision rules identical to the Wave H script. Objective 0.2192 -> 0.0935 dex; MFT 151.87 -> 242.44 ppb (2.25x under -> 1.41x under), FFT 243.72 -> 217.98 ppb (1.22x over -> 1.09x over, co-moving because the two lanes share upstream flux -- it was NOT fitted). THE PROFILE MINIMUM SITS AT THE RANGE FLOOR (argmin 23.30, hit_a_bound = true): the data want a value the table's own sulfur-addition class envelope does not contain, so 26.35 is the conservative edge of the indifference band and NOT the optimum -- the residual is not removable by this barrier. THE CAVEAT THIS NUMBER CARRIES, quoted verbatim from the fit target's own `content_verification_note` (Wave K): \"The paper's abstract reports MFT/FFT yields in mol % (e.g. MFT 1.4 mol % for its best precursor system), not ppb. The 342/200 ppb values here require a mol%->ppb conversion (system volume, precursor moles, molecular weights) that is NOT documented anywhere in this repo, and the full text is paywalled so the values could not be confirmed against the paper's tables. This is the panel's tightest contract (1.45x / 0.09 dex) resting on an unverified derivation.\" So this constant is FITTED AGAINST AN ANCHOR WHOSE UNIT CONVERSION IS UNVERIFIED; if the conversion is wrong the resulting error is LOCALISED HERE, in this one barrier, rather than distributed across the route. It must never be quoted as a measured or literature-derived barrier. END OF THE RETIRED FIT'S OWN TEXT. That Wave K caveat is itself now SUPERSEDED: the question was never whether the mol%->ppb conversion was documented, it was that there is no mol % measurement on the other end of the conversion. Record: results/validation/sulfur_barrier_refit_pentodiulose.{json,md}, RETRACTED 2026-08-27 (Wave S2c) but retained for provenance (it also supersedes the stale sulfur_barrier_refit_hofmann.{json,md}); generator: scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py -- DO NOT RE-RUN it against this benchmark; it needs a real target first (Wave S2b's ILL pack, then a rebuild in native mol %)."),
    "deoxyosone_reduction": (28.0, "ESTIMATED 2026-08-27 (Wave N). 1-deoxy-2,3-pentodiulose + 2[H] -> 1,4-dideoxypento-2,3-diulose + H2O: formal C4 deoxygenation of the 1-deoxyosone via the sugar reductone pool, required by the Cerny & Davidek route (see `thiol_addition_pentodiulose`). No direct kinetic measurement exists; taken as the `dehydration`/`furanone_cyclisation` analogue (28.0, Nursten 2005) since the physical step is a dehydration-plus-reduction of the same substrate that `furanone_cyclisation` (28.0) cyclodehydrates -- setting the two equal deliberately expresses NO prior preference between the norfuraneol branch and the MFT branch at the 1-deoxyosone fork. UNCONSTRAINED; not tuned to reproduce any output"),
    # ── WAVE P TIER (added 2026-08-27) ───────────────────────────────────
    # Every key below is ESTIMATED and every one is written EXPLICITLY rather than
    # left to `_canonical_fast_family`. That is not boilerplate: the fallthrough is
    # a known defect class in this table (Wave G1 fix 8 found eight families silently
    # inheriting DEFAULT_BARRIER = 45.0, a half-life of ~39,000 years at 150 C), and
    # two of these names would land somewhere plausible-but-wrong without an entry --
    # `furanone_reductive_opening` contains "ring"-adjacent language and
    # `mercaptoketone_aldol_addition` contains "addition". Explicit keys make the
    # lookup exact and the value arguable. NONE was chosen to reproduce any output;
    # each names the in-table analogue it was taken from.
    "mercaptoketone_formation": (28.60, "ESTIMATED 2026-08-27 (Wave P item 2). H2S addition to an alpha-dicarbonyl with reduction, giving an alpha-mercaptoketone: pyruvaldehyde + H2S + 2[H] -> 1-mercapto-2-propanone + H2O, and 2,3-pentanedione + H2S + 2[H] -> 2-mercapto-3-pentanone + H2O. Mechanism stated generally in Cerny 2015 (10.1016/b978-1-78242-103-0.00009-6): 'dicarbonyl compounds ... can add hydrogen sulfide (originating from L-cysteine decomposition) and then undergo reduction'; demonstrated on the C4 and C5 dicarbonyls by Mottram, Madruga & Whitfield 1995 (10.1021/jf00049a035). Taken as the `thiol_addition` CLASS value (28.60), i.e. the same un-fitted sulfur-addition analogue Wave N gave `thiol_addition_pentodiulose` -- the two ARE the same mechanistic class. UNCONSTRAINED: no kinetic measurement exists for either substrate, and the C3 case is an analogy to the C4/C5 demonstrations, not itself demonstrated"),
    "mercaptoketone_aldol_addition": (22.0, "ESTIMATED 2026-08-27 (Wave P item 2). Aldol ADDITION of glycolaldehyde to 1-mercapto-2-propanone giving Cerny's named intermediate 4,5-dihydroxy-3-mercapto-2-pentanone (10.1016/b978-1-78242-103-0.00009-6, rendering Hofmann & Schieberle 1998's Figure). A C-C bond FORMATION between two small carbonyls; taken as the reverse-direction analogue of `retro_aldol` is NOT appropriate (32.0 is a cleavage barrier), so the value is taken from `strecker_degradation` (22.0), the other carbonyl-addition-then-rearrangement step in this table operating on the same fragment size class. UNCONSTRAINED"),
    "mercaptoketone_cyclodehydration": (28.0, "ESTIMATED 2026-08-27 (Wave P item 2). 4,5-dihydroxy-3-mercapto-2-pentanone -> 2-methyl-3-furanthiol + 2 H2O: closure of the furan ring by double dehydration. Same transformation class as `furanone_cyclisation` (28.0) and `furan_ring_aromatisation` (28.0), and set equal to both so the model expresses NO preference between closing a furanone ring and closing this thiofuran ring. UNCONSTRAINED"),
    "furanone_reductive_opening": (28.0, "ESTIMATED 2026-08-27 (Wave P item 3). norfuraneol + 2 x 2[H] -> 2,3-pentanedione + H2O: reductive ring opening of the furanone to the C5 alpha-dicarbonyl that Whitfield & Mottram 1999 (10.1021/jf980980v) measure as a 'main non-sulfur compound' of the norfuraneol + cysteine/H2S system. It is the reverse transformation of `furanone_cyclisation` (28.0) on the same molecule, so it is set equal to it -- which deliberately expresses NO preference for keeping the ring closed or opening it. UNCONSTRAINED; not tuned to reproduce any output"),
    "furanone_amino_acid_reduction": (28.0, "ESTIMATED 2026-08-27 (Wave P item 6). hexose 1-deoxyosone + amino acid -> DMHF + Strecker aldehyde + CO2 + NH3 + H2O. This key REPLACES the routing of the hexose DMHF step through `furanone_cyclisation`, and is given the IDENTICAL value 28.0 so that any movement in predicted furaneol is attributable to the stoichiometry change (the `[HH]` pool gate was removed; Blank & Fay 1996, 10.1021/jf950439o, and Kerler et al. 2010, 10.1002/9781444317770.ch3, name the amino acid as the reductant) and NOT to a barrier change. Same `dehydration`/`furanone_formation` analogue as its sibling `furanone_formation` (28.0). UNCONSTRAINED"),
    "fructofuranosyl_dehydration": (28.0, "ESTIMATED 2026-08-27 (Wave P item 5). D-fructose -> HMF + 3 H2O by the ring-retained fructofuranosyl route (Antal, Mok & Richards 1990, 10.1016/0008-6215(90)84096-d; 3-deoxyglucosone excluded as the major fructose HMF precursor by Perez Locas & Yaylayan 2008, 10.1021/jf8010245). Set EQUAL to `dehydration` / `1,2-enolisation` (28.0), the barrier the replaced 3-deoxyosone -> HMF step carried, so the correction is purely TOPOLOGICAL: the model encodes no rate preference for the new route even though Perez Locas & Yaylayan measured fructose converting FASTER than 3-deoxyglucosone. That is the conservative direction and it is deliberately not tuned. UNCONSTRAINED"),
    "additive_thermal_degradation": (25.0, "ESTIMATED. Thiamine/GSH thermal cleavage; taken as the existing `thiamine_degradation` / `additive_degradation` analogues (both 25.0)"),
    "generalized_deamination": (21.0,  "ESTIMATED. Hydrolytic deamination of an alpha-aminoketone; taken as `enolisation_intermediate` (21.0), the directly analogous Amadori->deoxyosone + amine C-N cleavage"),
    "radical_propagation_o2":  (5.0,   "ESTIMATED. R. + O2 is effectively barrierless and diffusion-controlled (<2 kcal/mol in the literature); assigned this table's near-barrierless tier value (`ring_opening`/`mutarotation`, 5.0) rather than 0 so it does not become a numerically dominant rate"),
    "peroxy_h_abstraction":    (15.0,  "ESTIMATED. ROO. abstracting an allylic H is an H-atom transfer to a radical, the same class as `radical_crosstalk` (15.0, 'Radical + H2S collisions; fast'), which is the closest analogue in this table. Literature Ea for bis-allylic abstraction is lower (~6-8 kcal/mol); the conservative in-table analogue is used rather than an out-of-table number"),
    "radical_termination":     (5.0,   "ESTIMATED. Russell-mechanism radical-radical termination is diffusion-controlled; assigned this table's near-barrierless tier (`ring_opening`, 5.0)"),
    "furan_ring_aromatisation": (28.0, "ESTIMATED. Cyclodehydration + two-electron aromatisation closing a furan ring (thiamine -> 5-hydroxy-3-mercapto-2-pentanone -> MFT); taken as the `dehydration` analogue (28.0), the same cyclodehydration class as `furanone_cyclisation`"),
}

# Default barrier when no family pattern matches
DEFAULT_BARRIER: float = 45.0
DEFAULT_REFERENCE_PREEXPONENTIAL: float = 1e11
ARRHENIUS_R_KCAL: float = 0.001987

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})

DONOR_REACTIVITY_MULTIPLIERS: Dict[str, Dict[str, float]] = {
    "schiff_condensation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.94,
    },
    "amadori_rearrangement": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.94,
    },
    "heyns_rearrangement": {
        "phosphorylated": 1.04,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.98,
    },
    "enolisation_intermediate": {
        "phosphorylated": 1.12,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.92,
    },
    "1,2-enolisation": {
        "phosphorylated": 1.12,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.92,
    },
    "2,3-enolisation": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.95,
    },
    "dehydration": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.94,
    },
    "strecker_degradation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.88,
    },
    "aminoketone_condensation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.88,
    },
    "thiol_addition": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.90,
    },
    "thiohemiacetal_formation": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.92,
    },
    "thiol_dehydration": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.92,
    },
}

_DONOR_PRIORITY = {
    "phosphorylated": 4,
    "pentose": 3,
    "fructose": 2,
    "glucose": 1,
}

# Reaction keys for which a literature- or family-anchored barrier is
# carried as part of the computational-gap kinetic priors.  Originally
# these were intended as DFT anchors; after retirement of unreliable
# xTB/DFT-derived values (2026-04-21) the keys are kept for runtime
# wiring but the underlying tiers are literature/family/surrogate.
# `asparagine_sugar_explicit_water_cluster` is excluded here because it
# is wired through the safety pipeline (`src/safety.py`) using the
# Knol 2009 lumped Ea, not through the DFT anchor metadata.
_DFT_ANCHOR_REACTION_KEYS = (
    "hexanal_radical_quench",
    "quinone_cys_michael",
    "aa_ring_open_dicarbonyl",
    "pe_schiff_base",
    "pe_amadori",
    "lysinoalanine_crosslink",
)


def _fallback_dft_anchor_metadata() -> Dict[str, Dict[str, Any]]:
    return {
        "hexanal_radical_quench": {
            "slr_family": "11",
            "current_tier": "no_literature_anchor",
            "target_tier": "wet_lab_anchor",
            "active_key": "hexanal_radical_quench_no_anchor",
            "active_barrier_kcal": None,
            "dft_key": "hexanal_radical_quench_dft",
            "uncertainty_kj": None,
            "promotion_ceiling": "wet_lab_required",
            "honest_label": "No reliable kinetic anchor; xTB-derived value retired 2026-04-21. Runtime falls back to literature_runtime barrier_kj_mol=31.72 default.",
        },
        "quinone_cys_michael": {
            "slr_family": "13",
            "current_tier": "literature_family_surrogate",
            "target_tier": "literature_derived_transfer",
            "active_key": "quinone_cys_michael_thiol_addition_family",
            "active_barrier_kcal": 6.93,
            "dft_key": "quinone_cys_michael_dft",
            "uncertainty_kj": 15.0,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "Family surrogate (Michael thiol addition, ~29 kJ/mol literature); xTB value retired 2026-04-21",
        },
        "aa_ring_open_dicarbonyl": {
            "slr_family": "14",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "aa_ring_open_dicarbonyl_hcw_family14_surrogate",
            "active_barrier_kcal": 7.58,
            "dft_key": "aa_ring_open_dicarbonyl_dft",
            "uncertainty_kj": 20.0,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "HCW Family 14 surrogate for an upstream ascorbic-acid dicarbonyl source, DFT validation pending, ±factor 2-5; bounded calibration",
        },
        "pe_schiff_base": {
            "slr_family": "15",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "pe_schiff_base_lit_derived",
            "active_barrier_kcal": 22.21,
            "dft_key": "pe_schiff_base_dft",
            "uncertainty_kj": 20.9,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "SLR 15 literature Ea=92.9 kJ/mol, ±factor 2-5; bounded calibration",
        },
        "pe_amadori": {
            "slr_family": "15",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "pe_amadori_lit_derived",
            "active_barrier_kcal": 19.81,
            "dft_key": "pe_amadori_dft",
            "uncertainty_kj": 20.9,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "SLR 15 literature Ea=82.9 kJ/mol, ±factor 2-5; bounded calibration",
        },
        "lysinoalanine_crosslink": {
            "slr_family": "12",
            "current_tier": "family_rule_surrogate",
            "target_tier": "selective_dft_anchor",
            "active_key": "lysinoalanine_crosslink_dha_family_surrogate",
            "active_barrier_kcal": 16.0,
            "dft_key": "lysinoalanine_crosslink_dft",
            "uncertainty_kj": 42.0,
            "promotion_ceiling": "ranking_only",
            "honest_label": "DHA-plus-lysine family surrogate, assumes prior dehydroalanine formation, ±factor 5-15; ranking only",
        },
        "asparagine_sugar_explicit_water_cluster": {
            "slr_family": "13_safety",
            "current_tier": "literature_derived_transfer",
            "target_tier": "literature_derived_transfer",
            "active_key": "asparagine_sugar_explicit_water_cluster_knol2009",
            "active_barrier_kcal": 30.83,
            "dft_key": "asparagine_sugar_explicit_water_cluster_dft",
            "uncertainty_kj": 10.0,
            "promotion_ceiling": "literature_anchor",
            "honest_label": "Knol et al. 2009 Ea=129 kJ/mol (DOI 10.1016/j.foodchem.2009.11.049); direct lumped Asn+sugar->acrylamide literature anchor",
        },
    }


def _load_dft_anchor_metadata() -> Dict[str, Dict[str, Any]]:
    fallback = _fallback_dft_anchor_metadata()
    try:
        with open(COMPUTATIONAL_PRIORS_FILE, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return fallback

    entries = {
        str(row.get("reaction_key", "")).strip(): row
        for row in payload.get("dft_kinetic_priors", {}).get("entries", [])
        if str(row.get("reaction_key", "")).strip()
    }
    metadata: Dict[str, Dict[str, Any]] = {}
    for reaction_key in _DFT_ANCHOR_REACTION_KEYS:
        row = entries.get(reaction_key)
        if row is None:
            metadata[reaction_key] = dict(fallback[reaction_key])
            continue
        barrier_kcal = row.get("barrier_kcal_mol")
        active_barrier_kcal = (
            float(barrier_kcal)
            if barrier_kcal is not None
            else fallback[reaction_key]["active_barrier_kcal"]
        )
        uncertainty_raw = row.get("uncertainty_kj", fallback[reaction_key]["uncertainty_kj"])
        uncertainty_kj = float(uncertainty_raw) if uncertainty_raw is not None else None
        metadata[reaction_key] = {
            "slr_family": str(row.get("slr_family", fallback[reaction_key]["slr_family"])),
            "current_tier": str(row.get("current_tier", fallback[reaction_key]["current_tier"])),
            "target_tier": str(row.get("target_tier", fallback[reaction_key]["target_tier"])),
            "active_key": str(row.get("active_arrhenius_key", fallback[reaction_key]["active_key"])),
            "active_barrier_kcal": active_barrier_kcal,
            "dft_key": f"{reaction_key}_dft",
            "uncertainty_kj": uncertainty_kj,
            "promotion_ceiling": str(row.get("promotion_ceiling", fallback[reaction_key]["promotion_ceiling"])),
            "honest_label": str(row.get("honest_label", fallback[reaction_key]["honest_label"])),
        }
    return metadata


DFT_ANCHOR_METADATA = _load_dft_anchor_metadata()


def _normalize_donor_context_token(value: str) -> str:
    return " ".join(str(value).strip().lower().replace("_", " ").replace("-", " ").split())


def infer_carbohydrate_donor_identity(reactant_labels: Optional[Sequence[str]] = None) -> Optional[str]:
    if not reactant_labels:
        return None

    detected: Dict[str, int] = {}
    for raw_label in reactant_labels:
        token = _normalize_donor_context_token(str(raw_label))
        if not token:
            continue
        if any(item in token for item in ["ribose 5 phosphate", "glucose 6 phosphate", "fructose 6 phosphate", "r5p"]):
            detected["phosphorylated"] = max(detected.get("phosphorylated", 0), _DONOR_PRIORITY["phosphorylated"])
            continue
        if any(item in token for item in ["ribose", "xylose", "arabinose"]):
            detected["pentose"] = max(detected.get("pentose", 0), _DONOR_PRIORITY["pentose"])
            continue
        if "fructose" in token:
            detected["fructose"] = max(detected.get("fructose", 0), _DONOR_PRIORITY["fructose"])
            continue
        if "glucose" in token:
            detected["glucose"] = max(detected.get("glucose", 0), _DONOR_PRIORITY["glucose"])

    if not detected:
        return None
    return max(detected.items(), key=lambda item: item[1])[0]


def get_donor_reactivity_multiplier(
    reaction_family: Optional[str],
    *,
    reactant_labels: Optional[Sequence[str]] = None,
    donor_identity: Optional[str] = None,
) -> float:
    donor = str(donor_identity or infer_carbohydrate_donor_identity(reactant_labels) or "").strip().lower()
    if not donor:
        return 1.0

    family_key = _canonical_fast_family(reaction_family)
    if not family_key:
        return 1.0

    family_multipliers = DONOR_REACTIVITY_MULTIPLIERS.get(family_key, {})
    return float(family_multipliers.get(donor, 1.0))


def _normalize_family_key(reaction_family: Optional[str]) -> str:
    if not reaction_family:
        return ""
    return reaction_family.lower().replace(" ", "_").replace("-", "_")


def _canonical_fast_family(reaction_family: Optional[str]) -> str:
    fm = _normalize_family_key(reaction_family)
    if not fm:
        return ""
    if fm in FAST_BARRIERS:
        return fm

    if "enolisation" in fm:
        if "1,2" in str(reaction_family) or "1_2" in fm:
            return "1,2-enolisation"
        if "2,3" in str(reaction_family) or "2_3" in fm:
            return "2,3-enolisation"
        if "dha" in fm:
            return "beta_elimination"
        if "beta" in fm or "elimination" in fm:
            return "beta_elimination"
        return "1,2-enolisation"
    if "schiff" in fm:
        if "hydrolysis" not in fm and "reversion" not in fm:
            return "schiff_condensation"
        return "schiff_base_hydrolysis"
    if "retro" in fm:
        return "retro_aldol"
    if "lipid" in fm and "synergy" in fm:
        return "lipid_strecker_synergy"
    if "lipid" in fm:
        return "lipid_condensation"
    if "synergy" in fm:
        return "lipid_strecker_synergy"
    if "strecker" in fm:
        return "strecker_degradation"
    if "amadori" in fm:
        return "amadori_rearrangement"
    if "heyns" in fm:
        return "heyns_rearrangement"
    if "cysteine" in fm or "thermolysis" in fm:
        return "cysteine_thermolysis"
    if "thiol" in fm and "oxidation" in fm:
        return "thiol_oxidation"
    if "thiol" in fm and "addition" in fm and "hexose" in fm:
        return "thiol_addition_hexose"
    if "thiol" in fm and "addition" in fm:
        return "thiol_addition"
    if "pyrazine" in fm or "aminoketone" in fm:
        return "aminoketone_condensation"
    if "thiazole" in fm:
        return "lipid_thiazole"
    if "beta" in fm:
        return "beta_elimination"
    if "ring" in fm or "mutarotation" in fm:
        return "ring_opening"
    if "homolysis" in fm:
        return "lipid_homolysis"
    if "beta_scission" in fm:
        return "beta_scission"
    if "crosstalk" in fm:
        return "radical_crosstalk"
    return fm


def _arrhenius_yaml_key(family: Optional[str]) -> Optional[str]:
    canonical_family = _canonical_fast_family(family)
    if not canonical_family:
        return None
    yaml_key_map = {
        "schiff_condensation": "schiff_condensation",
        "schiff_base_hydrolysis": None,
        "amadori_rearrangement": "amadori",
        "heyns_rearrangement": "amadori",
        "1,2-enolisation": "enolisation",
        "2,3-enolisation": "enolisation",
        "enolisation_intermediate": "enolisation",
        "dehydration": "dehydration",
        "strecker_degradation": "strecker",
        "aminoketone_condensation": "pyrazine_condensation",
        "cysteine_thermolysis": "cysteine_thermolysis",
        "thiol_addition": "thiol_addition",
        "thiol_addition_hexose": "thiol_addition",
        "retro_aldol": "retro_aldol",
        "beta_elimination": "beta_elimination_dha",
        "thiamine_degradation": "thiamine_degradation",
        # WAVE T4 (2026-08-27) — AN UNREACHABLE LITERATURE ANCHOR, WIRED.
        # `data/lit/arrhenius_params.yml` carries a `thiamine_degradation` A/Ea
        # pair, and NO EMITTED FAMILY REACHED IT. The engine's thiamine steps
        # (`reaction_templates.py:1943-1966`) emit `Additive_Thermal_Degradation`,
        # which canonicalises to `additive_thermal_degradation` — a FAST_BARRIERS
        # key in its own right since Wave G1 fix 8 — so this map returned None and
        # `cantera_export.py:157-163` silently fell back to the heuristic
        # prefactor. Same defect class as the eight DEFAULT_BARRIER fallthroughs,
        # one lane over.
        #
        # SCOPE: the Cantera export lane ONLY (`get_arrhenius_params` ->
        # `src/cantera_export.py` -> `scripts/run_cantera_kinetics.py`). Per the
        # authority statement at the head of this module, FAST_BARRIERS is what
        # `evaluate_benchmark_payload` and every published prediction read, and
        # FAST_BARRIERS is untouched here. NO SHIPPED NUMBER MOVES.
        #
        # TWO CAVEATS, both stated rather than buried:
        #  1. THE FAMILY IS LUMPED. `Additive_Thermal_Degradation` covers the four
        #     thiamine steps AND the glutathione cleavage at `:2052`, so the GSH
        #     step now inherits the thiamine A/Ea. That is not new physics — the
        #     screening lane already gives both a single lumped 25.0 kcal/mol
        #     (see the `additive_thermal_degradation` note above, which names
        #     `thiamine_degradation` as its analogue) — but the Cantera lane now
        #     lumps them too, and there is no `glutathione_cleavage` entry in the
        #     YAML to split them with.
        #  2. THE ANCHOR ITSELF IS UNSOUND, and the fix does not launder it. That
        #     YAML entry is `source_quality: dead_doi` with an audit_flag: the DOI
        #     404s, the source string hedges "or similar lit", and Ea = 100.8
        #     kJ/mol is unanchored. Wiring it does not make it evidence. It makes
        #     it VISIBLE: `cantera_export` propagates `source_quality` verbatim, so
        #     the export now reports "dead_doi" where it previously reported
        #     "heuristic" for a step whose literature entry existed all along.
        #     Re-anchoring against the real thiamine kinetics the intake registry
        #     already names (Voelker 2018 / 2021) is an open [P].
        "additive_thermal_degradation": "thiamine_degradation",
        "ring_opening": "mutarotation",
        "mutarotation": "mutarotation",
        "lipid_thiazole": "pyrazine_condensation",
    }
    return yaml_key_map.get(canonical_family)


_REFINEMENT_PATCH_CACHE: Optional[Dict[str, float]] = None
_REFINEMENT_PATCH_MTIME: Optional[float] = None


def _load_refinement_surrogate_offsets() -> Dict[str, float]:
    global _REFINEMENT_PATCH_CACHE, _REFINEMENT_PATCH_MTIME
    if not REFINEMENT_PATCH_FILE.exists():
        _REFINEMENT_PATCH_CACHE = {}
        _REFINEMENT_PATCH_MTIME = None
        return {}

    mtime = REFINEMENT_PATCH_FILE.stat().st_mtime
    if _REFINEMENT_PATCH_CACHE is not None and _REFINEMENT_PATCH_MTIME == mtime:
        return dict(_REFINEMENT_PATCH_CACHE)

    try:
        with open(REFINEMENT_PATCH_FILE, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        _REFINEMENT_PATCH_CACHE = {}
        _REFINEMENT_PATCH_MTIME = mtime
        return {}

    offsets = payload.get("accepted_offsets", {}) if isinstance(payload, dict) else {}
    if not isinstance(offsets, dict):
        offsets = {}
    normalized = {str(key): float(value) for key, value in offsets.items()}
    _REFINEMENT_PATCH_CACHE = normalized
    _REFINEMENT_PATCH_MTIME = mtime
    return dict(normalized)


def get_barrier(reaction_family: str) -> Tuple[float, float]:
    """Return the FAST-mode (barrier, uncertainty) for a reaction family.
    
    Uncertainty mapping:
    * Heuristic fallback: ±5.0 kcal/mol
    * Literature calibrated: ±2.5 kcal/mol
    """
    default_unc = 5.0
    
    if not reaction_family:
        return DEFAULT_BARRIER, default_unc
        
    fm = _normalize_family_key(reaction_family)
    
    # --- DYNAMIC CALIBRATION OVERRIDES (Phase 1) ---
    import os
    offsets = _load_refinement_surrogate_offsets()
    if "BARRIER_OFFSETS" in os.environ:
        try:
            runtime_offsets = json.loads(os.environ["BARRIER_OFFSETS"])
            if isinstance(runtime_offsets, dict):
                offsets.update({str(key): float(value) for key, value in runtime_offsets.items()})
        except Exception:
            pass
    
    # Check for family-specific offset
    # Map optuna keys (short) to local fm keys
    offset_map = {
        "schiff": "schiff_condensation",
        "amadori": "amadori_rearrangement",
        "enol": "1,2-enolisation",
        "strecker": "strecker_degradation",
        "cys": "cysteine_thermolysis"
    }
    
    active_offset = 0.0
    if fm in offsets:
        active_offset = float(offsets[fm])
    for short_key, full_key in offset_map.items():
        if short_key in offsets and full_key in fm:
            active_offset = offsets[short_key]
            break
    # Check exact match first
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0] + active_offset, 3.5
    fm = _canonical_fast_family(reaction_family)

    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0] + active_offset, 3.5
    return DEFAULT_BARRIER + active_offset, 5.0

# Global cache for arrhenius parameters
_ARRHENIUS_CACHE = None

def get_arrhenius_params(family: str) -> Optional[Tuple[float, float, str, float]]:
    """
    Retrieve literature-calibrated (A_value, Ea_kcal_mol, source_quality, uncertainty) for a reaction family.
    Returns None if no data is found for the family.
    """
    global _ARRHENIUS_CACHE
    if _ARRHENIUS_CACHE is None:
        if ARRHENIUS_FILE.exists():
            with open(ARRHENIUS_FILE, "r") as f:
                data = yaml.safe_load(f)
                _ARRHENIUS_CACHE = data.get("arrhenius_data", {})
                # Warm cache with TST defaults for missing A values
                for k, v in _ARRHENIUS_CACHE.items():
                    if v.get("A_value") is None or (isinstance(v["A_value"], float) and math.isnan(v["A_value"])):
                        v["A_value"] = 6.25e12 # TST @ 150C
                        v["source_quality"] = "estimated_tst"
        else:
            _ARRHENIUS_CACHE = {}
            
    if not family:
        return None
        
    yaml_key = _arrhenius_yaml_key(family)

    if yaml_key and yaml_key in _ARRHENIUS_CACHE:
        entry = _ARRHENIUS_CACHE[yaml_key]
        A = float(entry.get("A_value", 0.0))
        Ea_kj = float(entry.get("Ea_kj_mol", 0.0))
        quality = entry.get("source_quality", "estimated")
        
        # Convert kJ/mol to kcal/mol
        Ea_kcal = Ea_kj / 4.184
        
        # Assign uncertainty based on quality
        quality_unc_map = {
            "literature": 2.0,
            "estimated_tst": 4.0,
            "heuristic": 5.0,
            "estimated": 3.5
        }
        uncertainty = quality_unc_map.get(quality, 3.5)
            
        return A, Ea_kcal, quality, uncertainty
        
    return None


def get_reference_pre_exponential(family: Optional[str] = None) -> float:
    params = get_arrhenius_params(family or "")
    if params is None:
        return DEFAULT_REFERENCE_PREEXPONENTIAL
    return float(params[0])


def arrhenius_rate_constant(
    barrier_kcal: float,
    temperature_kelvin: float,
    family: Optional[str] = None,
    multiplier: float = 1.0,
) -> float:
    if barrier_kcal >= 99.0:
        return 0.0

    pre_exponential = get_reference_pre_exponential(family)
    exponent = -float(barrier_kcal) / (ARRHENIUS_R_KCAL * float(temperature_kelvin))
    return pre_exponential * math.exp(exponent) * max(float(multiplier), 0.0)


def effective_barrier_from_rate_constant(
    rate_constant: float,
    temperature_kelvin: float,
    family: Optional[str] = None,
    *,
    fallback_barrier: float = 99.0,
) -> float:
    if rate_constant <= 0.0:
        return fallback_barrier

    pre_exponential = get_reference_pre_exponential(family)
    if pre_exponential <= 0.0:
        return fallback_barrier

    ratio = rate_constant / pre_exponential
    if ratio <= 0.0:
        return fallback_barrier

    return max(0.0, -ARRHENIUS_R_KCAL * float(temperature_kelvin) * math.log(ratio))


