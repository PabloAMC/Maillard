"""
The incumbent evidence audit (plan sec. 7 rule S1, owner 2026-11-06): one evidence tier per shipped coordinate.

A PURE READER. It reads the shipped parameter objects (``engine.core_parameters`` per Maillard lane, the lipid and
matrix registries), the frozen fit reports, their Laplace covariances and slice profiles, the B38 identifiability audit,
the envelope's priors (``core_prediction_uncertainty.json``), the sulfur reaction evidence table (flux shares), the data
wishlist, the unpinned-constants register and the conflicts register, and assigns every coordinate (log10 k_ref and Ea
SEPARATELY, and every band-carried parameter) one tier. It integrates nothing, fits nothing and moves no constant; the
only arithmetic is Arrhenius transport of a shipped constant to a measurement's temperature.

Tiers (S1 as amended by the coordinator's review, 2026-11-06):

* ``measured``              a direct measurement of that step under comparable conditions;
* ``analogue``              a same-class measurement transferred to this step or these conditions (S1: through a measured
                            law, with its band; flagged ``band_missing`` when no band is declared, a gap for S4);
* ``bounded``               only a measured or physical bound constrains it; ``bounded_sides`` says one- or two-sided (a
                            measurement of a SUBSET of a lumped step's channels is a one-sided LOWER bound when the
                            channels add; a measured TOTAL loss bounds each channel from above);
* ``fitted_identified``     a converged fit whose MARGINAL profile pins it, every other unidentified coordinate profiled
                            over its band. Any one unidentified criterion (band edge, flat or bound-limited profile, a CI
                            spanning the band or reaching an unphysical sign, no uncertainty, frozen since an unidentified
                            fit) overrides a reported stderr;
* ``fitted_unidentified``   everything else that is fitted. ``identification.status`` separates ``unidentified`` (a
                            local criterion fails) from ``locally_pinned_profile_missing`` (every local criterion passes,
                            but no artifact on disk profiles the other unidentified coordinates over their bands, so the
                            rule's profile requirement is unmet);
* ``structural``            a declared structural assumption (a declared zero, an inert arm, a declared barrier of zero,
                            a declared-fast step, an unsourced declared band);
* ``physical_constant``     a property of a substance (a partition constant), not a rate of this chemistry;
* ``conflicting_measurements``  direct measurements of the step disagree beyond their errors and the shipped value is
                            one of them (or a band across them);
* ``withdrawn``             quarantined or withdrawn from use; carried for the record only.

Aligned with sec. 7 as confirmed on branch sec7-symmetric (066684b7). S2 targets fitted-IDENTIFIED coordinates only (a
fitted-unidentified coordinate has no point to contradict and is carried under S4). It needs the measurement to be of
the step (all lumped channels, the engine's rate law and concentration regime), a 95 % interval of (measurement minus
fitted value) that excludes zero, with the reported noise, the within-lab between-paper spread, the between-lab term
(the step's measured spread where >= 2 laboratories measured it, else the lane's: H0's tau_ABS for the thiol lane) and
the coordinate's own interval in quadrature, undefined when a condition the step depends on is unstated; a full-lane
Holm pass; and two readers (an inferred bound or a patent is a bound only). A measurement of a SUBSET of a lumped
step's channels is a one-sided lower bound when the channels add. The contradiction records below
(``DIRECT_MEASUREMENTS``) carry each of those fields, hand-auditable, each citing its paper key, ledger row (or that
none is written yet), dossier and register entry; each names whether, under S4, it is a band source (it constrains the
value) or a rival (it implies a different rate law or species).
"""
from __future__ import annotations

import json
import math
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

import yaml

from src import data_paths

V = data_paths.VALIDATION_DIR
REF = data_paths.DOCS_ROOT / "reference"

TIERS: Tuple[str, ...] = (
    "measured", "analogue", "bounded", "fitted_identified", "fitted_unidentified",
    "structural", "physical_constant", "conflicting_measurements", "withdrawn",
)
FITTED_TIERS = ("fitted_identified", "fitted_unidentified")
R_KJ = 8.314462618e-3
LN10 = math.log(10.0)
#: An Ea band is shown in dex as the change of log10 k it makes between 100 C and 145 C (the span the shipped lanes'
#: reference temperatures cover): dex per kJ/mol = (1/373.15 - 1/418.15) / (R ln 10).
DEX_PER_KJ_100_145 = (1.0 / 373.15 - 1.0 / 418.15) / (R_KJ * LN10)
#: A coordinate within this fraction of its band's width from a declared bound sits on the band edge.
EDGE_FRACTION = 0.02

INPUTS = {
    "b1": V / "kinetic_core_b1_fit_report.json",
    "b3": V / "kinetic_core_b3_fit_report.json",
    "b6": V / "kinetic_core_b6_fit_report.json",
    "b7": V / "kinetic_core_b7_fit_report.json",
    "b18": V / "kinetic_core_b18_fit_report.json",
    "b20": V / "kinetic_core_b20_fit_report.json",
    "b21": V / "kinetic_core_b21_fit_report.json",
    "b39": V / "kinetic_core_b39_fit_report.json",
    "b60": V / "kinetic_core_b60_fit_report.json",
    "b60_laplace": V / "kinetic_core_b60_laplace_covariance.json",
    "b60_profile": V / "kinetic_core_b60_profile.json",
    "b23": V / "kinetic_core_b2_3_fit_report.json",
    "b38": V / "kinetic_core_b38_identifiability.json",
    "envelope": V / "core_prediction_uncertainty.json",
    "sulfur_evidence": V / "sulfur_reaction_evidence.json",
    "wishlist": V / "data_wishlist.json",
    "unpinned": REF / "unpinned_constants.yml",
    "conflicts": data_paths.CONFLICTS_REGISTER,
    "choices": REF / "choices_register.yml",
    "tier_changes": REF / "evidence_tier_changes.yml",
    "frozen_fit_drift": REF / "frozen_fit_drift_allowlist.yml",
}

#: The record's ``identification.fit`` label -> the frozen fit whose objective that coordinate's shipped value rests on
#: (scripts/generators/build_frozen_fit_reproducibility.py audits each; scripts/ci/frozen_fit_drift_gate.py keeps the
#: drifted ones in docs/reference/frozen_fit_drift_allowlist.yml). The fed 3-deoxy triangle the record attributes to
#: B39 ships B41's optimum (parameters_dicarbonyl.SHIPPING_FIT_REPORT), so B41 is the fit it rests on.
FIT_OF_IDENTIFICATION_LABEL: Dict[str, str] = {
    "B1 (variant A)": "B1", "B2.3 (frozen since)": "B2.3", "B60 (free)": "B60", "B3": "B3", "B6": "B6", "B7": "B7",
    "B18": "B18", "B20": "B20", "B21": "B21", "B39": "B41"}

# ---------------------------------------------------------------------------------------------------------------------
# What the record must cover, and how a tier may change (S1: "upgrades need an artifact and a second reader;
# downgrades need a citation"). ``scripts/ci/evidence_tier_gate.py`` reads these; the generator reads the same ones.
# ---------------------------------------------------------------------------------------------------------------------

#: The two coordinates of every shipped rate constant, tiered separately.
COORDINATES: Tuple[str, ...] = ("log10_k_ref", "ea_kj_mol")
#: Tiers a change can raise a coordinate TO, ranked: a move down this ranking, or out of it, is a downgrade.
UPGRADE_RANK: Dict[str, int] = {"measured": 3, "analogue": 2, "bounded": 1, "fitted_identified": 1}
#: The hand-kept register of tier changes since the incumbent audit (the frozen baseline). Each entry: ``coordinate``
#: (a row id), ``from`` (the tier before; null for a coordinate new since the baseline), ``to``, ``date``,
#: ``citation``; an UPGRADE also needs ``artifact`` (a file under results/validation/) and ``second_reader`` (who
#: re-read the evidence). The generator copies each entry onto its row as ``tier_changes``.
TIER_CHANGES = INPUTS["tier_changes"]


def is_upgrade(before: Optional[str], after: str) -> bool:
    """True when ``after`` is a stronger-evidence tier than ``before`` (a coordinate new since the baseline counts
    from no tier).

    A move into measured / analogue / bounded / fitted_identified is an upgrade unless it moves DOWN their ranking
    (measured -> analogue); a lateral move between equal ranks (bounded <-> fitted_identified) counts as an upgrade.
    Every other move (into fitted_unidentified, structural, conflicting_measurements, withdrawn, ...) is a downgrade.
    """
    if after not in UPGRADE_RANK or before == after:
        return False
    return before not in UPGRADE_RANK or UPGRADE_RANK[after] >= UPGRADE_RANK[before]


def tier_changes(path: Any = None) -> List[dict]:
    """The register's entries (``docs/reference/evidence_tier_changes.yml``), in file order."""
    return list((yaml.safe_load((path or TIER_CHANGES).read_text(encoding="utf-8")) or {}).get("entries") or [])


def shipped_parameters() -> Dict[str, Tuple[str, Any]]:
    """Every rate constant the engine reads, keyed by constant: (its first Maillard lane, the parameter object)."""
    from src.kinetic_core import engine

    parameters: Dict[str, Tuple[str, Any]] = {}
    for lane in engine.MAILLARD_LANES:
        for key, p in engine.core_parameters(lane).items():
            parameters.setdefault(key, (lane, p))
    return parameters


def shipped_coordinate_ids(parameters: Optional[Mapping[str, Tuple[str, Any]]] = None) -> List[str]:
    """The row id of every shipped coordinate (``<constant>.log10_k_ref`` and ``<constant>.ea_kj_mol``)."""
    return [f"{key}.{c}" for key in (shipped_parameters() if parameters is None else parameters) for c in COORDINATES]


# ---------------------------------------------------------------------------------------------------------------------
# Vocabularies the record replaces or maps
# ---------------------------------------------------------------------------------------------------------------------

VOCABULARY_MAP: Dict[str, Any] = {
    "evidence_class (engine KineticParameter / SulfurParameter / AcrylamideParameter)": {
        "measured_rate": "measured, or analogue when the measurement is transferred (another matrix, a recast order, a "
                         "declared transfer), or conflicting_measurements when the register records a disagreement",
        "measured_activation_energy": "as measured_rate, for the barrier coordinate",
        "derived_from_fit_data": "fitted_identified | fitted_unidentified (by the identification criteria); structural "
                                 "when the arm is not shipped and the constant is an inert zero; analogue for B21 (its "
                                 "fit rows are the source's own constants for the same step, transferred)",
        "bounded_from_a_timescale_bracket": "bounded (two-sided), or structural when the source is a declared-fast / "
                                            "unconstrained assumption, or analogue when it is a declared transfer",
        "structural_constant": "structural",
        "keep": True,
        "why": "it records what the source printed, which the tier needs; the tier records what that means for the shipped "
               "coordinate. Kept as an input, not as a tier.",
    },
    "sulfur evidence grade (sulfur_reaction_evidence.json)": {
        "measured": "measured | analogue | conflicting_measurements",
        "fitted": "fitted_identified | fitted_unidentified",
        "bounded": "bounded",
        "structural": "structural",
        "derived": "no own coordinate (the step reads another constant)",
        "keep": False,
        "why": "it is a coarser copy of evidence_class for one lane, and its 'unpinned' boolean is "
               "fitted_unidentified AND matters. The sulfur table should read the tier from this record (follow-up PR, "
               "with the unpinned_constants gate reading it too).",
    },
    "provenance_tier (data/lit/computational_priors.json, retention_reference_payloads.json, process_state_calibrations.json)": {
        "direct_measurement": "measured",
        "literature_derived_transfer": "analogue",
        "literature_bounded_provisional": "bounded",
        "repo_literature_synthesis": "analogue (unreviewed synthesis; band usually missing)",
        "direct_model_assumption": "structural",
        "unsourced_withdrawn": "withdrawn",
        "keep": True,
        "why": "it states the CITATION CHAIN of a registry record (citation_gate.py reads it), a different question; and "
               "the kinetic core reads none of those records (src/matrix_prior_registry.py is outside "
               "src/kinetic_core). Keep, with this map, until those records enter an engine lane.",
    },
    "lipid evidence_class (LipidParameter)": {
        "measured_composition": "measured", "hand_fitted_no_uncertainty": "analogue (the source authors' eyeballed value, "
        "transferred)", "declared_assumption": "structural, or analogue when the band spans measurements",
        "structural_constant": "structural", "derived_from_fit_data": "fitted_*", "keep": True,
        "why": "input, as evidence_class",
    },
    "matrix evidence_class (MatrixParameter)": {
        "measured_ratio": "measured (withdrawn when quarantined)", "derived_from_fit_rows": "fitted_*",
        "measured_bound": "bounded (one-sided)", "structural_gate": "structural", "declared_assumption": "structural",
        "reference_input": "structural", "keep": True, "why": "input, as evidence_class",
    },
}

# ---------------------------------------------------------------------------------------------------------------------
# Hand-declared readings. Each names its basis; the reader applies them, never invents them.
# ---------------------------------------------------------------------------------------------------------------------

#: Flags that make a measured constant an ANALOGUE: the measurement is of a related system or is recast.
TRANSFER_FLAGS = (
    "amine_free_glucose_melt", "declared_transfer", "bimolecular_recast", "alanine_to_glycine_declared",
    "refit_prefactor_not_published", "transfer_declared", "selectivity_collapses_in_a_low_moisture_matrix",
)
#: Flags that declare a barrier zero (no measured temperature dependence): a structural assumption.
ZERO_EA_FLAGS = ("ea_zero_by_declaration", "ea_declared_zero", "ea_author_fixed_to_zero", "ea_fixed_to_zero_by_authors",
                 "rate_is_exactly_zero_by_declaration")
#: Flags that declare a step structurally (fast, unconstrained) rather than bound it.
STRUCTURAL_BRACKET_FLAGS = ("declared_fast", "unconstrained_declared_assumption")

#: Coordinates whose direct measurements disagree beyond their errors (register entries or the constant's own flags).
CONFLICTING: Dict[Tuple[str, str], str] = {
    ("k_tdg_fa", "ea_kj_mol"): "flags ea_conflicts_with_knol_2010 and three_laboratories_disagree_b44 (Martins 2005 30 "
                               "kJ/mol vs Knol 2010; B44)",
    ("k_int1_acr", "ea_kj_mol"): "CF-039: 84-105 kJ/mol (Gokmen 2012, Knol 2005) vs 136-182 (Leuven, Claeys, Bassama); "
                                 "the envelope carries it flat across the disagreement",
    ("k_thioether", "log10_k_ref"): "CF-081: the bimolecular recast divides by a site density that one laboratory gives as "
                                    "0.29 (assay) and 8-10 mmol/g (a dose sentence), x30 apart",
}
#: Envelope prior reasons that mean the shipped coordinate is a band across disagreeing laboratories.
DISAGREE_REASON = "two_laboratories_disagree"

#: Waves whose arm is not shipped: their constants run at k = 0 (ship-rule file named).
NOT_SHIPPED_ARMS = {
    "b11_oxygen": "kinetic_core_b11_ship_rule.md (DO NOT SHIP: the consumers ship as declared-inert)",
    "b17_dimer_release": "kinetic_core_b17_ship_rule.md (DO NOT SHIP)",
    "b17a_site_yield": "kinetic_core_b17a_ship_rule.md (DO NOT SHIP)",
    "b25_thiol_addition": "kinetic_core_b25_ship_rule.md (DO NOT SHIP)",
    "b27_dicarbonyl_redox": "kinetic_core_b27_ship_rule.md (NOT FITTED -- DO NOT SHIP)",
    "b49_nf_shet": "kinetic_core_b49_ship_rule.md (DO NOT SHIP)",
    "b22_methionine": "kinetic_core_b22_ship_rule.md / b22b (DO NOT SHIP)",
    "b24_proline": "kinetic_core_b24_ship_rule.md / b24b (DO NOT SHIP)",
}

#: Band-carried envelope coordinates that are not a constant's own coordinate: (tier, basis, bounded_sides).
BAND_ROW_TIERS: Dict[str, Tuple[str, str, Optional[str]]] = {
    "lipid.q10": ("analogue", "CR-02: two measurements in other matrices (Tazi 2009 almond paste 1.6-3.3; Frankel 1993 "
                  "bulk n-6 oils, Q10 2.15-2.27) transferred to hot protein pots; the band spans both", None),
    "lipid.pea_protein_isolate.lipid_mass_fraction": ("structural", "LIPID_CARRIERS: 'literature-typical values already "
                  "present in the codebase', NOT corpus values (declared_assumption); CF-130 (carrier dilution)", None),
    "lipid.pea_protein_isolate.peroxide_value_meq_per_kg": ("structural", "LIPID_CARRIERS declared_assumption (as above)", None),
    "lipid.soy_protein_isolate.lipid_mass_fraction": ("structural", "LIPID_CARRIERS declared_assumption (as above)", None),
    "lipid.soy_protein_isolate.peroxide_value_meq_per_kg": ("structural", "LIPID_CARRIERS declared_assumption (as above)", None),
    "lipid.frankel_pure_hydroperoxide.lipid_mass_fraction": ("structural", "definitional (fed hydroperoxide)", None),
    "lipid.frankel_pure_hydroperoxide.peroxide_value_meq_per_kg": ("structural", "definitional (fed hydroperoxide)", None),
    "trunk.aw_multiplier_scale": ("bounded", "0 = Bell 1995's fixed-molality plateau, 1.2 = the source's 95 % CI "
                                  "(a measured bracket)", "two-sided"),
    "trunk.amadori_ph_exponent_decades_per_unit": ("analogue", "the six Martins 2003 per-step pH ratios span the band, "
                                                   "transferred to the Amadori steps", None),
    "acrylamide.aw_multiplier": ("conflicting_measurements", "CF-104 (open): De Vleeschouwer's powders give a flat "
                                 "formation multiplier, Bassama 2011 plantain x8.1 from a_w 0.972 to 0.430", None),
    "acrylamide.aw_elimination_deficit_scale": ("conflicting_measurements", "CF-104 (open): elimination x14 in Bassama "
                                                "2011 vs De Vleeschouwer 2007's k_E shape", None),
    "acrylamide.ph_exponent_formation_decades_per_unit": ("analogue", "band from the potato-matrix and phosphate slopes "
                                                          "(two measurements, other systems) +/- their SEs", None),
    "acrylamide.ph_exponent_elimination_decades_per_unit": ("analogue", "band from the two measured slopes +/- SEs", None),
    "sulfur.oxygen.k_cys_ox.log10_k": ("structural", "B11 not shipped: the consumer is zero by declaration", None),
    "sulfur.oxygen.k_red_ox.log10_k": ("structural", "B11 not shipped: the consumer is zero by declaration", None),
    "sulfur.oxygen.reservoir_scale": ("structural", "B11 not shipped: the reservoir is inert", None),
    "furanic.partition_ea_offset_kj_mol": ("structural", "a declared corner band on the furanone partition barrier "
                                           "(no measurement)", None),
    "observable.air_water_partition_constant": ("physical_constant", "K_aw of each compound (a substance property); "
                                                "declared +/- 0.5 dex band", None),
    "observable.hs_spme_same_sample_dispersion": ("measured", "Xin 2026 vs Xin 2026b same-sample 10-23x discrepancy "
                                                  "(Declaration D.6 Module 8)", None),
}

#: Fit rows later withdrawn or corrected, per fit (every free coordinate of that fit inherits the flag).
WITHDRAWN_OR_CORRECTED_FIT_ROWS: Tuple[Dict[str, Any], ...] = (
    {"fit": "B3", "row": "spi_extrusion_130C_acrylamide_ppb", "kind": "withdrawn charge",
     "issue": "the precursor charge (SPI 50, Asn 42, Glc 42 mM, pH 6.0, a_w 0.35) is not printed by Ma 2024; withdrawn "
              "from the bundle and the panel row refused 2026-11-06; the frozen B3 report still carries the row "
              "(residual -3.98 dex, sigma 1.0)",
     "cite": "CF-131 (filed on branch m0-data-corrections-batch1, not yet on main); "
             "results/validation/m0_integrity_batch1_other_lanes.md MF-39",
     "size_dex": -3.98},
    {"fit": "B20", "row": "every B20 row (Nguyen 2016 rate table)", "kind": "fitted under a network the source contradicts",
     "issue": "r_flp_decay returns the lysine that Nguyen 2016's network consumes, so the fitted constants compensate "
              "(k_glyc at 120 C is 0.64x Nguyen's printed k3)",
     "cite": "CF-072 (error_confirmed; CF-052 merged into it)", "size_dex": math.log10(0.64)},
    {"fit": "B60", "row": "the two MFT-disulfide share rows (declared hold-out, in-sample since B60)", "kind": "corrected target",
     "issue": "the B2.1-B2.3 generators use a slipped molar mass for C10H10O2S2 (0.0323 -> 0.0326, 0.0429 -> 0.0433, "
              "0.2711 -> 0.2736); immaterial (+0.004 dex against sigma_log 0.4)",
     "cite": "results/validation/m0_integrity_batch1_strecker_sulfur.md (MFTD share check, 'integrity holds')",
     "size_dex": 0.004},
)

#: The direct measurements that bear on a shipped coordinate: contradictions and one-sided bounds. Fields:
#:   relation: same_channel | superset_upper_bound (the measured total loss bounds each channel from above) |
#:             subset_lower_bound (a measured subset of a lumped step's channels bounds it from below) |
#:             lumped_total_vs_channel (a measured total compared with a lumped step that omits some channels);
#:   measured: {kind: rate | upper_bound | lower_bound | ea, value, unit, temperature_c};
#:   band components (one-sigma dex terms, combined in quadrature): noise, within_lab_between_paper, between_lab (the
#:            step's measured spread when >= 2 laboratories measured it; None = use the lane's term), reaction_class;
#:   depends_on_unstated: conditions the step depends on that the source does not state (the interval is then undefined);
#:   inferred_bound: a bound inferred rather than measured, or a patent (S2 (d): a bound only);
#:   s4_route: band (constrains the value) | rival (implies a different rate law or species);
#:   matches: whether the measured system's channels, rate law, concentration range and stated conditions match the
#:            engine's step (S2 (a)).
DIRECT_MEASUREMENTS: Tuple[Dict[str, Any], ...] = (
    {"id": "DM-01", "constant": "k_nf_decay", "coordinate": "log10_k_ref", "lane": "sulfur",
     "measurement": "norfuraneol fed alone, first-order loss",
     "sources": ["bork2022", "bork2023"], "ledger_rows": [], "ledger_note": "proposed in the bork2022 / bork2023 dossiers, not yet written",
     "cite": "CF-112; bork2022_extraction.md Fig. 1D; bork2023_extraction.md Fig. 4B",
     "lab": "tu_berlin_kanzler (TU Berlin food chemistry; one laboratory)",
     "system": "NF 20 mM alone in water, pH 5 initial, unbuffered, 130 C, 0-300 min",
     "channels": "every NF loss in an NF-alone pot; in the engine that pot loses NF only through r_nf_decay (B61's "
                 "hydroxide term is off), so the channels coincide",
     "rate_law": "first order (0-300 min; 4.80e-3 weighted, 5.1e-3 unweighted, 4.4e-3 two-point 0-120 min)",
     "engine_rate_law": "first order, pH-free", "concentration": "20 mM NF (first order: concentration-free if the order holds)",
     "unstated": ["O2 / headspace", "pH drift (unbuffered)"],
     "depends_on_unstated": ["pH drift (NF loss is pH-dependent: B61's hydroxide term, Mikami 2017)",
                             "O2 / headspace (an O2 or metal route, CF-117)"],
     "s4_route": "band", "s4_route_reason": "same law and species: it constrains the value (M1b part (ii) carries it as a band)",
     "relation": "same_channel",
     "measured": {"kind": "rate", "value": 4.8e-3, "unit": "1/min", "temperature_c": 130.0},
     "band": {"noise": math.log10(4.98 / 4.8), "within_lab_between_paper": math.log10(6.1e-3 / 4.8e-3),
              "within_lab_assay": math.log10(1.4), "between_lab": None, "reaction_class": None},
     "matches": {"channels": True, "rate_law": True, "concentration": True, "conditions_stated": False},
     "note": "the engine value at 130 C uses the carbonyl-sink barrier (143.5 kJ/mol, itself unidentified) to transport "
             "k(145 C); the gap is the rate and the barrier together"},
    {"id": "DM-02", "constant": "k_mft_decay", "coordinate": "log10_k_ref", "lane": "sulfur",
     "measurement": "co-heated [2H3]-MFT standard survives 240 min at 110 C",
     "sources": ["us10076126"], "ledger_rows": ["us10076126_fig1_mft_time_course"],
     "ledger_note": "the survival is a Methods statement; the bound is derived in m0_checks_batch2 CK-23",
     "cite": "CF-013 (explained); results/validation/m0_checks_batch2.md CK-23",
     "lab": "the patent's assignee laboratory (US 10,076,126); one source, no standard amount, recovery or LOD printed",
     "system": "xylose + cysteine 0.25 M, 0.2 M phosphate pH 5.5, 110 C, standard heated with the pot",
     "channels": "every MFT loss in the pot (a superset of ch_mft_decay)",
     "rate_law": "survival only: lambda <= ln(1/f_min)/240 min, f_min 1e-6 - 1e-2 declared (0.058 - 0.019 /min)",
     "engine_rate_law": "first order, pot-independent", "concentration": "trace standard in a 250 mM cysteine pot",
     "unstated": ["standard amount, recovery, LOD", "O2 / headspace"],
     "depends_on_unstated": ["the standard's amount and recovery (f_min is declared, not printed)"],
     "inferred_bound": True,
     "s4_route": "rival", "s4_route_reason": "with Grosch 1992 it implies a cysteine-dependent loss law (case M1, C2), "
                                             "not a value of the pot-independent sink",
     "relation": "superset_upper_bound",
     "measured": {"kind": "upper_bound", "value": 0.058, "unit": "1/min", "temperature_c": 110.0},
     "band": {"noise": None, "within_lab_between_paper": None, "between_lab": None, "reaction_class": None},
     "matches": {"channels": True, "rate_law": False, "concentration": False, "conditions_stated": False},
     "note": "the bound is taken at its loose corner (f_min 1e-6, 0.058 /min); at the tight corner (0.019) the gap is "
             "0.70 dex. Grosch 1992's inferred in-pot loss at 100 C (0.044 /min without cysteine, 0.034 with 0.65 mM) agrees with "
             "the engine (0.048); no pot-independent sink at one Ea in Gigl 2021's range meets both (CK-23), which "
             "nominates a loss that free cysteine slows (case M1, combination C2). The rate law is the mismatch."},
    {"id": "DM-03", "constant": "k_aa_frag", "coordinate": "log10_k_ref", "lane": "trunk",
     "measurement": "acetic acid on a plateau 60-240 min",
     "sources": ["davidek2008"], "ledger_rows": [], "ledger_note": "no ledger row written",
     "cite": "CF-069 (open)", "lab": "Nestle (Davidek); one laboratory",
     "system": "xylose + glycine, 0.2 M phosphate pH 8, 120 C, isotope dilution",
     "channels": "every acetic acid loss (a superset of r_aa_frag), read from a NET plateau",
     "rate_law": "a net plateau (26.9 -> 27.6 mol %): loss <= ~4e-4 /min only if formation has stopped",
     "engine_rate_law": "first order", "concentration": "0.2 M sugar / amine (Martins: 0.2 M glucose + glycine)",
     "unstated": ["formation and loss not separated", "pH 8 vs the fit's 6.8", "pentose, not glucose"],
     "depends_on_unstated": ["acetic acid formation over the plateau (a net plateau bounds the loss only if formation "
                             "has stopped)"],
     "inferred_bound": True,
     "s4_route": "band", "s4_route_reason": "an upper bound on the value of a first-order loss",
     "relation": "superset_upper_bound",
     "measured": {"kind": "upper_bound", "value": 4e-4, "unit": "1/min", "temperature_c": 120.0},
     "band": {"noise": None, "within_lab_between_paper": None, "between_lab": None, "reaction_class": None},
     "matches": {"channels": True, "rate_law": False, "concentration": True, "conditions_stated": False},
     "note": "the register's reading: the fitted loss stands in for another deficiency of the Martins fit (a missing late "
             "acid sink, or the formic / acetic split)"},
    {"id": "DM-04", "constant": "k_mgo_mel", "coordinate": "ea_kj_mol", "lane": "trunk",
     "measurement": "apparent Ea of MGO loss with amino acids, 80 -> 100 C",
     "sources": ["kanzler2026", "kanzler2021"], "ledger_rows": [], "ledger_note": "proposed in kanzler2026 dossier sec. 7, not written",
     "cite": "CF-122 (open); kanzler2026_extraction.md secs. 4.7-6", "lab": "tu_berlin_kanzler (one laboratory)",
     "system": "MGO + Ala / Phe / Pro 100 mM each, pH 5 unbuffered, 60-140 C, 30 min",
     "channels": "MGO loss with an amino acid: the amine channels plus MGO alone (r_mgo_mel lumps them but carries no amine)",
     "rate_law": "MGO alone second order (k2 ~1.85e-4 /mM/min; Kanzler 2021 5.9e-4, 0.51 dex apart); excess over MGO alone "
                 "~[Ala]^0.75", "engine_rate_law": "first order, amine-free, pH-free",
     "concentration": "100 mM MGO vs ~0.1 mM in the trunk's pots (x1000)", "unstated": ["pH drift", "O2", "come-up time"],
     "depends_on_unstated": ["pH drift (unbuffered)", "O2"],
     "s4_route": "rival", "s4_route_reason": "it implies a different rate law (second order in MGO, ~[amine]^0.75): a "
                                             "rival MGO-sink structure (CF-122 remedy), not a value. The MGO-alone k2 of "
                                             "Kanzler 2021 and 2026 conflict (0.51 dex, one laboratory)",
     "relation": "lumped_total_vs_channel",
     "measured": {"kind": "ea", "value": 95.0, "high": 121.0, "unit": "kJ/mol", "temperature_c": None},
     "band": {"noise": None, "within_lab_between_paper": 0.51, "between_lab": None, "reaction_class": None},
     "matches": {"channels": False, "rate_law": False, "concentration": False, "conditions_stated": False},
     "note": "the shipped Ea 20.04 sits on FITTED_EA_BOUNDS' floor with B1's 95 % CI -42 to 82 kJ/mol"},
    {"id": "DM-05", "constant": "k_mgo_mel", "coordinate": "log10_k_ref", "lane": "trunk",
     "measurement": "MGO alone, amine-free first-order loss",
     "sources": ["pfeifer2013"], "ledger_rows": [], "ledger_note": "not written",
     "cite": "CF-122 note (2); pfeifer2013_extraction.md secs. 0 and 4",
     "lab": "tu_berlin_kanzler lineage (same corresponding author and method; not a second laboratory)",
     "system": "MGO ~10 mM alone, pH 5 initial, 130 C, aerated, no replicates",
     "channels": "the amine-free channel only: a SUBSET of the lumped r_mgo_mel",
     "rate_law": "first order, k1 ~0.0053 /min (0.0042 - 0.0065)", "engine_rate_law": "first order",
     "concentration": "~10 mM", "unstated": ["pH drift"],
     "relation": "subset_lower_bound",
     "measured": {"kind": "lower_bound", "value": 0.0042, "unit": "1/min", "temperature_c": 130.0},
     "band": {"noise": math.log10(0.0065 / 0.0042) / 2.0, "within_lab_between_paper": None, "between_lab": None,
              "reaction_class": None},
     "matches": {"channels": False, "rate_law": True, "concentration": False, "conditions_stated": True},
     "note": "a one-sided bound, not a contradiction: the lumped loss may exceed its amine-free part"},
    {"id": "DM-06", "constant": "k_glyc", "coordinate": "ea_kj_mol", "lane": "trunk",
     "measurement": "FDNB-available lysine loss in soy, response-surface barrier",
     "sources": ["wolf1977a"], "ledger_rows": [], "ledger_note": "not written",
     "cite": "CF-100 (open)", "lab": "Wolf / Reineccius (one laboratory)",
     "system": "soy protein + glucose 2-4 %, a_w 0.65, pH 7, 80-130 C (moisture 8-12 g/100 g, not printed)",
     "channels": "the initial lysine-loss slope (rival-independent: both CF-072 structures share it)",
     "rate_law": "glucose order 0.7-1.3 (brackets the engine's 1)", "engine_rate_law": "second order (lysine x glucose)",
     "concentration": "low-moisture powder", "unstated": ["moisture"],
     "depends_on_unstated": ["moisture (8 - 12 g/100 g, not printed; the rate is per water basis)"],
     "s4_route": "band", "s4_route_reason": "constrains the barrier's value; CF-100's decider is a soy / pea time course",
     "relation": "same_channel",
     "measured": {"kind": "ea", "value": 137.7, "high": 137.7, "unit": "kJ/mol", "temperature_c": None},
     "band": {"noise": 8.4 * DEX_PER_KJ_100_145, "within_lab_between_paper": None, "between_lab": None,
              "reaction_class": None},
     "matches": {"channels": True, "rate_law": True, "concentration": False, "conditions_stated": False},
     "note": "the engine over-predicts by +0.24 to +0.42 dex at 100 C and +0.5 to +0.7 at 80 C (matched at 125 C); the "
             "barrier is an ANALOGUE (Martins' k_schiff), so S2 does not apply as written"},
    {"id": "DM-07", "constant": "k_ttca_deg", "coordinate": "log10_k_ref", "lane": "sulfur",
     "measurement": "fed TTCA, total first-order loss",
     "sources": ["zhai2021"], "ledger_rows": [], "ledger_note": "see zhai2021 dossier sec. 2",
     "cite": "CF-108 (use the printed k, never the printed Ea); CF-110", "lab": "Jiangnan (Zhai / Zhang)",
     "system": "TTCA ~10 mM in unbuffered water, natural pH 5.5, 100 C",
     "channels": "every TTCA loss (a superset of r_ttca_deg; the reverse exchange is the other)",
     "rate_law": "first order, k 0.0271 /min", "engine_rate_law": "first order, pH-free",
     "concentration": "~10 mM", "unstated": ["in-pot pH (unbuffered)"],
     "relation": "superset_upper_bound",
     "measured": {"kind": "upper_bound", "value": 0.0271, "unit": "1/min", "temperature_c": 100.0},
     "band": {"noise": None, "within_lab_between_paper": None, "between_lab": None, "reaction_class": None},
     "matches": {"channels": True, "rate_law": True, "concentration": True, "conditions_stated": False},
     "note": "a one-sided bound; the pH dependence (x9.6 faster at pH 8, CF-110) is a structural gap, not this coordinate"},
    {"id": "DM-08", "constant": "k_fur_decay", "coordinate": "ea_kj_mol", "lane": "sulfur",
     "measurement": "furfural loss in a cysteine system, measured barrier",
     "sources": ["yaghmur2005"], "ledger_rows": [], "ledger_note": "not written",
     "cite": "kinetic_core_b62_prereg.md (B62 off until a SHIP verdict); parameters_sulfur.FURFURAL_BARRIERS",
     "lab": "Yaghmur et al. (one laboratory)", "system": "furfural + cysteine, pH 5, 40-70 C",
     "channels": "furfural loss with cysteine present (includes channels r_fur_decay does not carry)",
     "rate_law": "first order in furfural", "engine_rate_law": "first order, carbonyl-sink family barrier",
     "concentration": "model solution", "unstated": ["O2"],
     "s4_route": "rival", "s4_route_reason": "furfural's own barriers (B62's switch) are a rival to the shared "
                                             "carbonyl-sink barrier, not a value of it",
     "relation": "lumped_total_vs_channel",
     "measured": {"kind": "ea", "value": 46.5, "high": 46.5, "unit": "kJ/mol", "temperature_c": None},
     "band": {"noise": None, "within_lab_between_paper": None, "between_lab": None, "reaction_class": None},
     "matches": {"channels": False, "rate_law": True, "concentration": False, "conditions_stated": False},
     "note": "measured at 40-70 C, the engine runs it at 100-145 C; Schirle-Keller 1992 measures net furfural FORMATION "
             "at 147 kJ/mol, so the two direct numbers bracket the shipped 143.5 from different channels"},
)

#: The incumbent's structural choices, with the evidence against each and whether a case or a registered choice ever
#: examined it (review item 5).
STRUCTURAL_CHOICES: Tuple[Dict[str, str], ...] = (
    {"id": "SC-01", "choice": "MGO -> melanoidin as a first-order, amine-free, pH-free sink (r_mgo_mel, k_mgo_mel)",
     "where": "src/kinetic_core/network.py r_mgo_mel", "against": "CF-122 (order, amine dependence, Ea), CF-095",
     "case": "none", "choice_register": "none"},
    {"id": "SC-02", "choice": "the bound Amadori compound's decay returns the protein lysine (r_flp_decay)",
     "where": "src/kinetic_core/network.py r_flp_decay", "against": "CF-072 (error_confirmed), CF-096",
     "case": "none", "choice_register": "none"},
    {"id": "SC-03", "choice": "one pot-independent first-order free-thiol sink (ch_mft_decay; k_mft_decay frozen since B2.3)",
     "where": "src/kinetic_core/sulfur.py ch_mft_decay", "against": "CF-013, CF-107; m0_checks_batch2 CK-23 retires it as the sole answer",
     "case": "M1 (thiol loss, oxidant, pH drift)", "choice_register": "none"},
    {"id": "SC-04", "choice": "one lumped formation barrier (64.08 kJ/mol, unidentified since B2.3) shared by every "
                              "single-temperature sulfur constant", "where": "parameters_sulfur lumped formation Ea",
     "against": "B10 ship rule (RE-MERGE: barriers not identified); B55 samples a measured-class distribution instead",
     "case": "none", "choice_register": "CR-07"},
    {"id": "SC-05", "choice": "one barrier per decay family (thiol sink 101.25 at its 102 ceiling; carbonyl sink 143.5 "
                              "unidentified) shared by MFT, FFT, furfural, osone and norfuraneol losses",
     "where": "B60 decay_Ea_kJ_mol", "against": "CF-112; B62 (furfural's own measured barriers, off); B38 P3",
     "case": "none", "choice_register": "none"},
    {"id": "SC-06", "choice": "norfuraneol loss pH-free (r_nf_decay; B61's hydroxide term off)",
     "where": "src/kinetic_core/sulfur.py r_nf_decay", "against": "CF-112 (Bork: engine 3.1-4.3x too slow at pH 5)",
     "case": "M1b (carries Bork as a band)", "choice_register": "none"},
    {"id": "SC-07", "choice": "the cysteine / pentose <-> TTCA exchange pH-flat, K held at Cao 2017's 90 C / pH 5.5 value, "
                              "TTCA loss without a pH law", "where": "r_ttca_cys, r_ttca_deg",
     "against": "CF-076, CF-011, CF-110", "case": "M1 (M1c supply structure)", "choice_register": "CR-13"},
    {"id": "SC-08", "choice": "no non-Strecker dicarbonyl + amine adduct channel (amine and dicarbonyl budgets omit it)",
     "where": "network (absent step)", "against": "CF-095 (>= 1e3x the engine's amine consumption in fed-MGO pots)",
     "case": "none (an A3 candidate)", "choice_register": "none"},
    {"id": "SC-09", "choice": "H2S-consuming steps see all H2S in solution (no vessel gas / liquid partition)",
     "where": "sulfur lane H2S steps", "against": "m0_checks_batch2 CK-34: the fitted furfural constants absorb it (+0.84 dex)",
     "case": "none", "choice_register": "none"},
    {"id": "SC-10", "choice": "the thioether step recast bimolecular on a melanoidin site density of 9 mmol/g",
     "where": "parameters_sulfur MELANOIDIN_SITE_DENSITY_MMOL_PER_G", "against": "CF-081 (0.29 vs 8-10 mmol/g, x30)",
     "case": "none", "choice_register": "none"},
    {"id": "SC-11", "choice": "the pH state's amine-fate basis (acid released per sink event, cysteine carboxyl carried)",
     "where": "src/kinetic_core/ph_state.py AMINE_FATE_BASIS", "against": "CF-084 (Meynier drift <= 0.2 vs engine -0.79)",
     "case": "M1b (pH balance)", "choice_register": "none"},
    {"id": "SC-12", "choice": "first-order fragmentation sinks to FRAG_C for glucose, formic and acetic acid (B1)",
     "where": "src/kinetic_core/parameters.py k_glc_frag / k_fa_frag / k_aa_frag",
     "against": "CF-069 (acetic acid stable at 120 C); k_glc_frag and k_fa_frag unidentified, k_glc_frag on its floor",
     "case": "none", "choice_register": "none"},
    {"id": "SC-13", "choice": "pentose -> deoxyosone steps with cysteine as a catalyst (PENT + Cys -> TDP / DPO + Cys)",
     "where": "src/kinetic_core/sulfur.py r_pent_tdp / r_pent_dpo", "against": "unpinned_constants.yml (de Roos's "
     "measured route to replace it); k_pent_* unidentified since B2.3", "case": "none (planned sub-model 2)",
     "choice_register": "none"},
    {"id": "SC-14", "choice": "first order in amine for the deoxypentosone -> acetylformoin channel (measured n 0.35-0.71)",
     "where": "k_dpo_af", "against": "the constant's own flag first_order_in_amine_but_measured_n_is_0.35_to_0.71",
     "case": "none", "choice_register": "none"},
    {"id": "SC-15", "choice": "one hydroperoxide decomposition rate anchored at 25 C, carried to cooking temperature by a Q10",
     "where": "parameters_lipid Q10_ASSUMPTION", "against": "the anchor authors licensed adjustment, not ~11 decades of "
     "extrapolation; Tazi 2009 and Frankel 1993 disagree", "case": "none", "choice_register": "CR-02"},
    {"id": "SC-16", "choice": "acrylamide formation flat in a_w over 0.34-0.99", "where": "acrylamide_conditions AW_MULTIPLIER",
     "against": "CF-104 (Bassama 2011 x8.1)", "case": "M5 (buffer, fat, water activity)", "choice_register": "none"},
    {"id": "SC-17", "choice": "Zhou 2024's alanine Strecker constants carried to glycine (B18)", "where": "k_go_ak, k_mgo_ak",
     "against": "B18 barriers on their band edges (B38 AT_BOUND)", "case": "M3 (Strecker aldehydes)", "choice_register": "none"},
)


# ---------------------------------------------------------------------------------------------------------------------
# Small readers
# ---------------------------------------------------------------------------------------------------------------------

def _json(name: str) -> Any:
    return json.loads(INPUTS[name].read_text(encoding="utf-8"))


def _yaml(name: str) -> Any:
    return yaml.safe_load(INPUTS[name].read_text(encoding="utf-8")) or {}


def arrhenius(k_ref: float, ea_kj_mol: Optional[float], t_ref_k: float, temperature_c: float) -> float:
    """k at ``temperature_c`` from k at ``t_ref_k`` (no barrier: held, as the engine holds it)."""
    if ea_kj_mol is None:
        return float(k_ref)
    t = float(temperature_c) + 273.15
    return float(k_ref) * math.exp(-(float(ea_kj_mol) / R_KJ) * (1.0 / t - 1.0 / float(t_ref_k)))


def ea_dex(delta_kj_mol: float) -> float:
    """An Ea difference as the log10 k change it makes between 100 and 145 C."""
    return float(delta_kj_mol) * DEX_PER_KJ_100_145


#: The lane's between-laboratory term (sec. 7, S2 (b)): for the thiol (sulfur) lane, H0's tau_ABS; a lane without H0
#: uses the pooled between-laboratory spread of the reaction class, which this record does not yet carry (None).
LANE_BETWEEN_LAB_DEX: Dict[str, Optional[float]] = {"sulfur": 0.964}
LANE_BETWEEN_LAB_SOURCE = ("H0 / EX-059 D2 branch (i): tau_ABS 0.964 dex (results/validation/core_panel_scores.md); "
                           "other lanes: no H0 and no pooled reaction-class spread on record yet")


def comparison_interval_dex(band: Mapping[str, Optional[float]], lane: str,
                            fitted_sigma_dex: Optional[float] = None) -> Tuple[Optional[float], Optional[float], str]:
    """Sec. 7 S2 (b): the 95 % interval of (measurement - fitted value) with the reported noise, the within-lab
    between-paper spread, the between-lab term and the fitted coordinate's own interval in quadrature. Each component is
    read as a one-sigma dex term. The between-lab term is the step's measured spread where >= 2 laboratories measured it
    (``band['between_lab']``), else the lane's term; with neither, the interval is undefined. Returns (sd, 95 % half
    width, note)."""
    terms = {k: abs(float(v)) for k, v in band.items() if v is not None and k != "between_lab"}
    between = band.get("between_lab")
    between_src = "the step's measured between-lab spread"
    if between is None:
        between, between_src = LANE_BETWEEN_LAB_DEX.get(lane), "the lane's between-lab term"
    if between is None:
        return None, None, ("undefined: no between-lab term (one laboratory measured the step, and the lane has no H0 "
                            "or pooled reaction-class spread on record)")
    terms["between_lab"] = abs(float(between))
    if fitted_sigma_dex is not None:
        terms["fitted_coordinate_interval"] = abs(float(fitted_sigma_dex))
    sd = math.sqrt(sum(v * v for v in terms.values()))
    note = ("quadrature of " + ", ".join(f"{k} {v:.3f}" for k, v in terms.items())
            + f" (between-lab from {between_src}); 95 % = 1.96 x {sd:.3f}")
    return sd, 1.96 * sd, note


def local_unidentified_criteria(*, value: Optional[float], lo: Optional[float], hi: Optional[float],
                                ci95_halfwidth: Optional[float], fit_says_identified: Optional[bool],
                                at_bound: Optional[bool] = None, slice_verdict: Optional[str] = None,
                                b38_verdict: Optional[str] = None, ci_lower: Optional[float] = None,
                                is_barrier: bool = False, frozen_note: Optional[str] = None) -> List[str]:
    """Every local criterion that marks a fitted coordinate unidentified (any one suffices)."""
    out: List[str] = []
    if frozen_note:
        out.append(frozen_note)
    if fit_says_identified is False:
        out.append("the fit's own record: not identified")
    if at_bound:
        out.append("on a declared bound (the fit's record)")
    if value is not None and lo is not None and hi is not None and hi > lo:
        if min(abs(value - lo), abs(hi - value)) <= EDGE_FRACTION * (hi - lo) and not at_bound:
            out.append(f"on the band edge ({value:.4g} within {EDGE_FRACTION:.0%} of [{lo:.4g}, {hi:.4g}])")
        if ci95_halfwidth is not None and 2.0 * ci95_halfwidth >= (hi - lo):
            out.append(f"95 % CI (+/- {ci95_halfwidth:.3g}) spans the band [{lo:.4g}, {hi:.4g}]")
    if ci95_halfwidth is None and not frozen_note:
        out.append("no uncertainty reported")
    if is_barrier and ci_lower is not None and ci_lower <= 0.0:
        out.append(f"95 % CI reaches an unphysical sign (lower {ci_lower:.3g} kJ/mol)")
    if slice_verdict in ("flat", "bound_limited"):
        out.append(f"slice profile {slice_verdict}")
    if b38_verdict in ("UNIDENTIFIED", "AT_BOUND"):
        out.append(f"B38 identifiability audit: {b38_verdict}")
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Identification per fit
# ---------------------------------------------------------------------------------------------------------------------

class _Fits:
    """Every fit artifact the shipped coordinates come from, read once."""

    def __init__(self) -> None:
        self.b1 = _json("b1")["variant_A"]["fitted_parameters"]
        self.b3 = _json("b3")["parameter_intervals"]
        self.b3_rows = {r["id"] for r in _json("b3")["rows"]}
        self.b18 = _json("b18")
        self.b20 = _json("b20")
        self.b21 = _json("b21")
        self.b39 = _json("b39")
        self.b60 = _json("b60")["frozen_parameters"]
        lap = _json("b60_laplace")
        self.lap = {}
        for i, c in enumerate(lap["coordinates"]):
            self.lap[(c["key"], c["kind"])] = {
                "value": lap["optimum"][i], "bounds": lap["bounds"][i], "sigma": lap["sigma"][i],
                "ci95": lap["ci95_halfwidth"][i], "identified": lap["identified"][i], "at_bound": lap["at_bound"][i]}
        self.profile = {(r["key"], r["kind"]): r["verdict"] for r in _json("b60_profile")["rows"]}
        self.b23 = _json("b23")["parameter_intervals"]
        self.b38: Dict[Tuple[str, str], dict] = {}
        for fit in _json("b38")["fits"]:
            wave = fit["fit"].split(" ")[0]
            for c in fit["coordinates"]:
                self.b38[(wave, c["key"])] = c
        from src.kinetic_core.parameters import FITTED_BOUNDS_LOG10K, FITTED_EA_BOUNDS
        from src.kinetic_core.parameters_acrylamide import (FITTED_ACRYLAMIDE_BOUNDS_LOG10K,
                                                           FITTED_ACRYLAMIDE_EA_BOUNDS)
        from src.kinetic_core.parameters_sulfur import FITTED_SULFUR_BOUNDS_LOG10K
        self.b1_bounds, self.b1_ea_bounds = FITTED_BOUNDS_LOG10K, FITTED_EA_BOUNDS
        self.b3_bounds, self.b3_ea_bounds = FITTED_ACRYLAMIDE_BOUNDS_LOG10K, FITTED_ACRYLAMIDE_EA_BOUNDS
        self.sulfur_bounds = FITTED_SULFUR_BOUNDS_LOG10K

    # Each returns {fit, artifact, value, bounds, ci95_halfwidth, labels, criteria}
    def sulfur_rate(self, key: str) -> dict:
        lap = self.lap.get((key, "log10k"))
        if lap is not None:
            lo, hi = lap["bounds"]
            ci = lap["ci95"] if lap["sigma"] is not None else None
            verdict = self.profile.get((key, "log10k"))
            crit = local_unidentified_criteria(value=lap["value"], lo=lo, hi=hi, ci95_halfwidth=ci,
                                               fit_says_identified=lap["identified"], at_bound=lap["at_bound"],
                                               slice_verdict=verdict)
            return {"fit": "B60 (free)", "artifact": "kinetic_core_b60_laplace_covariance.json; kinetic_core_b60_profile.json",
                    "value": lap["value"], "bounds": [lo, hi], "ci95_halfwidth": ci,
                    "labels": {"B60 Laplace": "identified" if lap["identified"] else "unidentified",
                               "B60 slice profile": verdict}, "criteria": crit}
        iv = self.b23.get(key) or {}
        lo, hi = self.sulfur_bounds.get(key, (None, None))
        frozen = ("frozen since B2.3 (not free in B60): " +
                  ("B2.3 identified it, but at other coordinates' values; not re-identified at the shipped optimum"
                   if iv.get("identified") else "B2.3 did not identify it"))
        crit = local_unidentified_criteria(value=iv.get("value"), lo=lo, hi=hi, ci95_halfwidth=iv.get("ci95_halfwidth"),
                                           fit_says_identified=iv.get("identified"), frozen_note=frozen)
        return {"fit": "B2.3 (frozen since)", "artifact": "kinetic_core_b2_3_fit_report.json parameter_intervals",
                "value": iv.get("value"), "bounds": [lo, hi], "ci95_halfwidth": iv.get("ci95_halfwidth"),
                "labels": {"B2.3 report": "identified" if iv.get("identified") else "unidentified"}, "criteria": crit}

    def sulfur_family_ea(self, family: str) -> dict:
        lap = self.lap[(family, "Ea")]
        lo, hi = lap["bounds"]
        ci = lap["ci95"] if lap["sigma"] is not None else None
        b23 = self.b23.get(f"Ea_decay_{family}") or {}
        crit = local_unidentified_criteria(value=lap["value"], lo=lo, hi=hi, ci95_halfwidth=ci,
                                           fit_says_identified=lap["identified"], at_bound=lap["at_bound"], is_barrier=True,
                                           ci_lower=(lap["value"] - ci) if ci is not None else None)
        return {"fit": "B60 (free)", "artifact": "kinetic_core_b60_laplace_covariance.json (decay_Ea_kJ_mol)",
                "value": lap["value"], "bounds": [lo, hi], "ci95_halfwidth": ci,
                "labels": {"B60 Laplace": "identified" if lap["identified"] else "unidentified",
                           "B2.3 report": "identified" if b23.get("identified") else "unidentified"}, "criteria": crit}

    def sulfur_lumped_ea(self) -> dict:
        from src.kinetic_core.parameters_sulfur import LUMPED_FORMATION_EA_BOUNDS
        iv = self.b23.get("Ea_lumped_formation") or {}
        lo, hi = LUMPED_FORMATION_EA_BOUNDS
        crit = local_unidentified_criteria(value=iv.get("value"), lo=lo, hi=hi, ci95_halfwidth=iv.get("ci95_halfwidth"),
                                           fit_says_identified=iv.get("identified"), is_barrier=True,
                                           frozen_note="frozen since B2.3 (B10 tried two route barriers: RE-MERGE, not identified)")
        return {"fit": "B2.3 (frozen since)", "artifact": "kinetic_core_b2_3_fit_report.json Ea_lumped_formation; "
                "kinetic_core_b10_ship_rule.md", "value": iv.get("value"), "bounds": [lo, hi], "ci95_halfwidth": None,
                "labels": {"B2.3 report": "unidentified"}, "criteria": crit}

    def b1_coord(self, key: str, coordinate: str) -> dict:
        fp = self.b1[key]
        if coordinate == "log10_k_ref":
            lo, hi = self.b1_bounds[key]
            se = fp.get("log10_k_ref_stderr")
            ci = 1.96 * se if se is not None else None
            crit = local_unidentified_criteria(value=fp["log10_k_ref_100C"], lo=lo, hi=hi, ci95_halfwidth=ci,
                                               fit_says_identified=None if se is not None else False,
                                               at_bound=fp.get("at_bound"))
            label = fp.get("identifiability_k")
            value = fp["log10_k_ref_100C"]
        else:
            lo, hi = self.b1_ea_bounds
            se = fp.get("ea_stderr_kj_mol")
            ci = 1.96 * se if se is not None else None
            cil = (fp.get("ea_ci95_kj_mol") or [None])[0]
            crit = local_unidentified_criteria(value=fp["ea_kj_mol"], lo=lo, hi=hi, ci95_halfwidth=ci,
                                               fit_says_identified=None if se is not None else False, is_barrier=True,
                                               ci_lower=cil)
            label = fp.get("identifiability_ea")
            value = fp["ea_kj_mol"]
        return {"fit": "B1 (variant A)", "artifact": "kinetic_core_b1_fit_report.json variant_A.fitted_parameters",
                "value": value, "bounds": [lo, hi], "ci95_halfwidth": ci, "labels": {"B1 report": label}, "criteria": crit}

    def b3_coord(self, name: str, coordinate: str, key: str) -> dict:
        iv = self.b3[name]
        lo, hi = (self.b3_bounds.get(key, (None, None)) if coordinate == "log10_k_ref" else self.b3_ea_bounds)
        b38 = self.b38.get(("B3", name)) or {}
        crit = local_unidentified_criteria(value=iv["value"], lo=lo, hi=hi, ci95_halfwidth=iv.get("ci95_halfwidth"),
                                           fit_says_identified=iv.get("identified"), at_bound=b38.get("at_bound"),
                                           b38_verdict=b38.get("verdict"), is_barrier=coordinate == "ea_kj_mol",
                                           ci_lower=(iv["value"] - iv["ci95_halfwidth"]) if coordinate == "ea_kj_mol" else None)
        return {"fit": "B3", "artifact": "kinetic_core_b3_fit_report.json parameter_intervals; "
                "kinetic_core_b38_identifiability.json", "value": iv["value"], "bounds": [lo, hi],
                "ci95_halfwidth": iv.get("ci95_halfwidth"),
                "labels": {"B3 report": "identified" if iv.get("identified") else "unidentified", "B38": b38.get("verdict")},
                "criteria": crit}

    def wave_coord(self, wave: str, report: dict, block: str, name: str, *, is_barrier: bool = False) -> dict:
        frozen = report["frozen_parameters"][block]
        value = frozen[name]
        lap = report.get("laplace") or {}
        if isinstance(lap.get("sigma"), list):  # B18: lists in the frozen block's order
            idx = list(frozen).index(name)
            sigma, ident, onb = lap["sigma"][idx], lap["identified"][idx], lap["on_bound"][idx]
        else:
            sigma = (lap.get("sigma") or {}).get(name)
            ident = (lap.get("identified") or {}).get(name)
            onb = (lap.get("on_bound") or {}).get(name)
        if "bounds" in report and name in report["bounds"]:
            lo, hi = report["bounds"][name]
        elif "bands_log10" in report:
            lo, hi = report["bands_log10"]["lower"].get(name), report["bands_log10"]["upper"].get(name)
        else:
            lo = hi = None
        ci = 1.96 * sigma if sigma is not None else None
        b38 = self.b38.get((wave, name)) or {}
        crit = local_unidentified_criteria(value=value, lo=lo, hi=hi, ci95_halfwidth=ci, fit_says_identified=ident,
                                           at_bound=onb, b38_verdict=b38.get("verdict"), is_barrier=is_barrier,
                                           ci_lower=(value - ci) if (is_barrier and ci is not None) else None)
        return {"fit": wave, "artifact": f"kinetic_core_{wave.lower()}_fit_report.json laplace; "
                "kinetic_core_b38_identifiability.json", "value": value, "bounds": [lo, hi], "ci95_halfwidth": ci,
                "labels": {f"{wave} Laplace": "identified" if ident else "unidentified", "B38": b38.get("verdict")},
                "criteria": crit}


# ---------------------------------------------------------------------------------------------------------------------
# Envelope priors -> coordinates
# ---------------------------------------------------------------------------------------------------------------------

def _prior_targets(prior_key: str, parameters: Mapping[str, Any]) -> List[Tuple[str, str]]:
    """The (constant, coordinate) pairs one envelope prior samples; [] for a band-carried parameter of its own."""
    parts = prior_key.split(".")
    wave = parts[0]
    if wave in ("lipid", "observable", "furanic") or prior_key.startswith("sulfur.oxygen") or prior_key.startswith("trunk."):
        return []
    if wave == "acrylamide":
        return [(parts[1], "ea_kj_mol")] if len(parts) == 3 and parts[1] in parameters else []
    if wave == "b3" and parts[1].startswith("Ea_"):
        return {"Ea_int1_mel": [("k_int1_mel", "ea_kj_mol")], "Ea_acr_dp": [("k_acr_dp", "ea_kj_mol")],
                "Ea_competitor_sugar": [(k, "ea_kj_mol") for k in ("k_gln_glc", "k_lys_glc", "k_ala_glc")]}[parts[1]]
    if wave == "b8" and parts[1] == "lumped_formation_Ea_kJ_mol":
        return [(k, "ea_kj_mol") for k, p in parameters.items() if "shares_lumped_formation_Ea" in getattr(p, "flags", ())
                and not any(f.startswith("b2_2_decay_family:") for f in p.flags)]
    if wave == "b8" and parts[1] == "decay_Ea_kJ_mol":
        fam = parts[2]
        return [(k, "ea_kj_mol") for k, p in parameters.items() if f"b2_2_decay_family:{fam}" in getattr(p, "flags", ())]
    if wave == "b8" and parts[1] == "ph_drift":
        return []
    if wave == "b18":
        m = re.match(r"(log10_|ea_)(k_)?([a-z]+_ak)", parts[1])
        if m and "ph_slope" not in parts[1]:
            return [("k_" + m.group(3), "log10_k_ref" if m.group(1) == "log10_" else "ea_kj_mol")]
        return []
    if wave == "b39":
        m = re.match(r"log10_(k_[a-z_]+)_100C", parts[1])
        return [(m.group(1), "log10_k_ref")] if m else []
    if len(parts) >= 3 and parts[1] in parameters:
        return [(parts[1], "ea_kj_mol" if parts[2].startswith("ea") else "log10_k_ref")]
    return []


def _band_from_prior(prior: dict, coordinate: str) -> dict:
    lo, hi = (prior.get("band") or [None, None])[:2]
    sigma = prior.get("sigma")
    dist = prior.get("distribution")
    width = None
    if dist == "normal_log10" and sigma:
        width = 2 * 1.645 * float(sigma)
    elif dist == "normal" and sigma:
        width = 2 * 1.645 * float(sigma) * (DEX_PER_KJ_100_145 if coordinate == "ea_kj_mol" else 1.0)
    elif lo is not None and hi is not None:
        width = (float(hi) - float(lo)) * (DEX_PER_KJ_100_145 if coordinate == "ea_kj_mol" else 1.0)
        if dist in ("log_uniform", "log_uniform_dispersion") and lo and hi and lo > 0:
            width = math.log10(float(hi) / float(lo))
    elif dist == "uniform_band" and sigma:
        width = 2 * 1.645 * float(sigma) * (DEX_PER_KJ_100_145 if coordinate == "ea_kj_mol" else 1.0)
    return {"source": "envelope prior " + prior["key"], "distribution": dist, "lo": lo, "hi": hi, "sigma": sigma,
            "sampled": bool(prior.get("sampled")), "width_dex": width, "reason": prior.get("reason")}


def _prior_label(reason: str) -> Optional[str]:
    r = (reason or "").lower()
    if "unidentified" in r or "frozen in the sulfur fit" in r:
        return "unidentified"
    if "identified" in r or "laplace_covariance" in r:
        return "identified"
    return None


# ---------------------------------------------------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------------------------------------------------

def _t_ref_k(p: Any) -> float:
    from src.kinetic_core.parameters import T_REF_K
    return float(getattr(p, "t_ref_k", None) or T_REF_K)


def _has(flags: Iterable[str], needles: Iterable[str]) -> Optional[str]:
    for f in flags:
        for n in needles:
            if n in f:
                return f
    return None


def _inert_arm(flags: Iterable[str]) -> Optional[str]:
    for f in flags:
        if f in NOT_SHIPPED_ARMS:
            return NOT_SHIPPED_ARMS[f]
    return None


def _conflict_hits(key: str, register: List[dict]) -> List[dict]:
    out = []
    for c in register:
        if re.search(rf"\b{re.escape(key)}\b", json.dumps(c, default=str)):
            out.append({"id": c["id"], "kind": c["kind"], "status": c["status"]})
    return out


def _outputs(key: str, lane: str) -> List[str]:
    from src.kinetic_core import engine
    skip = ("FRAG", "CBX", "OXV", "MEL", "ACID", "LYSP")
    out = set()
    for r in engine._lane_reactions(lane):
        if r.parameter_key == key:
            out |= {s for s in tuple(r.reactants) + tuple(r.products) if not s.startswith(skip)}
    return sorted(out)


def _measured_tier(p: Any, coordinate: str) -> Tuple[str, str]:
    """measured | analogue | conflicting_measurements for a measured constant's coordinate."""
    flags = tuple(p.flags)
    if (p.key, coordinate) in CONFLICTING:
        return "conflicting_measurements", CONFLICTING[(p.key, coordinate)]
    f = _has(flags, TRANSFER_FLAGS)
    if f:
        return "analogue", f"measured in another system and transferred (flag {f}): {p.source_anchor[:160]}"
    return "measured", p.source_anchor[:200]


def _coordinate(p: Any, lane: str, coordinate: str, fits: _Fits, priors: Dict[Tuple[str, str], dict]) -> dict:
    flags = tuple(p.flags)
    k_ref, ea = p.k_ref, p.ea_kj_mol
    row: Dict[str, Any] = {
        "id": f"{p.key}.{coordinate}", "constant": p.key, "lane": lane, "coordinate": coordinate,
        "transformation": p.transformation, "evidence_class": p.evidence_class,
        "reference_temperature_c": round(_t_ref_k(p) - 273.15, 2),
        "value": (math.log10(k_ref) if (k_ref and k_ref > 0) else (None if k_ref is None else 0.0))
        if coordinate == "log10_k_ref" else ea,
        "unit": f"log10({p.unit})" if coordinate == "log10_k_ref" else "kJ/mol",
        "dossiers": [], "bounded_sides": None, "band_missing": False, "identification": None,
    }
    if coordinate == "log10_k_ref" and k_ref == 0:
        row["value"] = None
        row["value_note"] = "k = 0"
    inert = _inert_arm(flags)
    prior = priors.get((p.key, coordinate))
    tier, basis = None, None
    ident = None
    fitted = p.evidence_class == "derived_from_fit_data"
    # --- structural zeros and inert arms
    if inert and (k_ref == 0 or k_ref is None):
        tier, basis = "structural", f"inert zero: the arm is not shipped ({inert})"
    elif p.evidence_class == "structural_constant":
        tier, basis = "structural", "structural constant: " + ", ".join(flags)
    elif prior is not None and DISAGREE_REASON in (prior.get("reason") or ""):
        tier, basis = "conflicting_measurements", f"envelope {prior['key']}: {prior['reason']}"
    elif coordinate == "ea_kj_mol" and ea is None:
        tier, basis = "structural", "no barrier carried: the constant is held at its measured temperature"
    elif coordinate == "ea_kj_mol" and _has(flags, ZERO_EA_FLAGS):
        tier, basis = "structural", f"barrier declared zero (flag {_has(flags, ZERO_EA_FLAGS)})"
    elif p.evidence_class == "bounded_from_a_timescale_bracket":
        if _has(flags, STRUCTURAL_BRACKET_FLAGS):
            tier, basis = "structural", f"declared, not bounded (flag {_has(flags, STRUCTURAL_BRACKET_FLAGS)})"
        elif _has(flags, ("declared_transfer",)):
            tier, basis = "analogue", "a declared transfer from the pentose measurement: " + p.source_anchor[:160]
        else:
            tier, basis = "bounded", "a timescale bracket, no rate: " + p.source_anchor[:160]
            row["bounded_sides"] = "two-sided"
    elif p.evidence_class in ("measured_rate", "measured_activation_energy"):
        tier, basis = _measured_tier(p, coordinate)
    elif fitted and coordinate == "ea_kj_mol":
        tier, basis, ident = _fitted_barrier(p, lane, fits)
    elif fitted:
        tier, basis, ident = _fitted_rate(p, lane, fits)
    if tier is None:
        tier, basis = "unclassified", "no rule fired"
    row.update(tier=tier, basis=basis, identification=ident)
    # --- band
    # S4: "a search bound is never a band". The fit's search span is recorded for ranking only.
    scale = DEX_PER_KJ_100_145 if coordinate == "ea_kj_mol" else 1.0
    bounds = (ident or {}).get("bounds")
    row["search_bounds_span_dex"] = ((bounds[1] - bounds[0]) * scale
                                     if bounds and None not in bounds else None)
    row["band_is_search_bound"] = False
    if prior is not None:
        row["band"] = _band_from_prior(prior, coordinate)
        b = row["band"]
        if (bounds and None not in bounds and b.get("lo") is not None and b.get("hi") is not None
                and abs(float(b["lo"]) - bounds[0]) < 1e-6 and abs(float(b["hi"]) - bounds[1]) < 1e-6):
            row["band_is_search_bound"] = True
            b["note"] = "the envelope samples the fit's SEARCH bounds: under S4 a search bound is never a band"
    else:
        ci = p.k_ref_ci95 if coordinate == "log10_k_ref" else p.ea_ci95_kj_mol
        if ci and (coordinate == "ea_kj_mol" or (k_ref and k_ref > 0)):
            if coordinate == "log10_k_ref":
                width = math.log10((k_ref + ci) / max(k_ref - ci, k_ref * 1e-3))
            else:
                width = 2 * ci * DEX_PER_KJ_100_145
            row["band"] = {"source": "the source's printed 95 % interval", "distribution": "not sampled",
                           "half_width": ci, "sampled": False, "width_dex": width}
        else:
            row["band"] = None
    if tier == "analogue" and row["band"] is None:
        row["band_missing"] = True
    return row


def _label_conflicts(ident: Optional[dict], prior: Optional[dict]) -> List[dict]:
    if not ident or not ident.get("criteria"):
        return []
    out = []
    labels = dict(ident.get("labels") or {})
    if prior is not None:
        labels[f"envelope {prior['key']}"] = _prior_label(prior.get("reason") or "")
    for src, lab in labels.items():
        lab_s = str(lab or "").lower()
        says_identified = lab_s in ("identified", "pinned") or lab_s.startswith("constrained")
        if says_identified:
            out.append({"label_source": src, "label": lab, "fails": ident["criteria"]})
    return out


def _fitted_rate(p: Any, lane: str, fits: _Fits) -> Tuple[str, str, dict]:
    flags = tuple(p.flags)
    key = p.key
    if lane == "sulfur" and key in fits.b60["log10_k_ref_at_145C"]:
        ident = fits.sulfur_rate(key)
    elif key in fits.b1:
        ident = fits.b1_coord(key, "log10_k_ref")
    elif key in fits.b3:
        ident = fits.b3_coord(key, "log10_k_ref", key)
    elif key == "k_dpo_af":
        ident = {"fit": "B7", "artifact": "kinetic_core_b7_fit_report.json frozen_parameters", "value": math.log10(p.k_ref),
                 "bounds": None, "ci95_halfwidth": None,
                 "labels": {"envelope": "the fit report carries no parameter stderr; residual sigma used"},
                 "criteria": ["no parameter uncertainty in the fit report (the envelope uses the residual sigma)"]}
    elif "fitted_wave_b18" in flags:
        ident = fits.wave_coord("B18", fits.b18, "pyrazine", f"log10_{key}_100C")
    elif "fitted_wave_b20" in flags:
        ident = fits.wave_coord("B20", fits.b20, "glycation", f"log10_{key}_100C")
    elif "fitted_wave_b21" in flags:
        name = "log10_k_g_go_aqueous_100C" if key == "k_g_go" else f"log10_{key}_100C"
        ident = fits.wave_coord("B21", fits.b21, "aqueous_glyoxal", name)
        return ("analogue", "B21's fit rows are Hamzalioglu 2026's own first-order constants for this step (lactulosyl-"
                "lysine in milk), carried to fructosyl-glycine in water as a declared transfer (+/- 0.5 dex, "
                "AQUEOUS_TRANSFER_BAND_DECADES)", ident)
    elif "fitted_wave_b39" in flags:
        ident = fits.wave_coord("B39", fits.b39, "fed_3deoxy", f"log10_{key}_100C")
    else:
        return "unclassified", "fitted, but no fit artifact is mapped for it", None
    return _fitted_tier(ident)


def _fitted_tier(ident: dict) -> Tuple[str, str, dict]:
    if ident["criteria"]:
        ident["status"] = "unidentified"
        return "fitted_unidentified", f"{ident['fit']}: " + "; ".join(ident["criteria"]), ident
    ident["status"] = "locally_pinned_profile_missing"
    ident["criteria_profile"] = ("no artifact profiles this coordinate with the other unidentified coordinates profiled over "
                                 "their bands (slice profiles hold them at the optimum; Laplace is local)")
    return ("fitted_unidentified", f"{ident['fit']}: every local criterion passes; the profile requirement is unmet", ident)


def _fitted_barrier(p: Any, lane: str, fits: _Fits) -> Tuple[str, str, Optional[dict]]:
    flags = tuple(p.flags)
    fam = next((f.split(":", 1)[1] for f in flags if f.startswith("b2_2_decay_family:")), None)
    if "measured_Ea_override" in flags:
        from src.kinetic_core.parameters_sulfur import MEASURED_EA_OVERRIDE_ANCHORS
        anchor = MEASURED_EA_OVERRIDE_ANCHORS.get(p.key, "")
        return ("analogue", "a measured barrier of a lumped or related step carried onto this channel (fit cannot move "
                f"it): {str(anchor)[:200]}", None)
    if fam:
        return _fitted_tier(fits.sulfur_family_ea(fam))
    if lane == "sulfur" and ("shares_lumped_formation_Ea" in flags or p.key in ("k_dimer_decay", "k_h2s_loss")):
        return _fitted_tier(fits.sulfur_lumped_ea())
    if p.key in fits.b1:
        return _fitted_tier(fits.b1_coord(p.key, "ea_kj_mol"))
    if p.ea_ci95_kj_mol and _has(flags, TRANSFER_FLAGS):
        return ("analogue", f"the source's measured barrier (+/- {p.ea_ci95_kj_mol:g} kJ/mol) in another system "
                f"(flag {_has(flags, TRANSFER_FLAGS)}); the rate is fitted, the barrier is not: {p.source_anchor[:120]}", None)
    if lane == "acrylamide":
        pol = next((f.split(":", 1)[1] for f in flags if f.startswith("ea_policy:")), "")
        if pol == "measured_ea_e2":
            return ("analogue", "k_acr_cys's measured 51.3 kJ/mol (De Vleeschouwer 2009 Part II) carried to this "
                    "amine's channel (ea_policy measured_ea_e2)", None)
        name = {"k_int1_mel": "Ea_int1_mel", "k_acr_dp": "Ea_acr_dp"}.get(p.key, "Ea_competitor_sugar")
        return _fitted_tier(fits.b3_coord(name, "ea_kj_mol", p.key))
    if "fitted_wave_b18" in flags:
        name = {"k_go_ak": "ea_go_ak_kj_mol", "k_mgo_ak": "ea_mgo_ak_kj_mol"}[p.key]
        return _fitted_tier(fits.wave_coord("B18", fits.b18, "pyrazine", name, is_barrier=True))
    f = _has(flags, ("barrier_from_", "barrier_declared", "barrier_measured_", "barrier_is_"))
    if f:
        return ("analogue", f"a barrier borrowed from another step or system (flag {f}); not fitted", None)
    if "barrier_consistent_with_zero" in flags:
        return ("analogue", "Hamzalioglu 2026 Table 2's glucosone -> glyoxal barrier in milk (4.2 +/- 15.7 kJ/mol), "
                "carried to water", None)
    return "unclassified", "fitted constant's barrier: no rule fired", None


def _engine_value_at(p: Any, coordinate: str, measured: Mapping[str, Any]) -> Optional[float]:
    if measured["kind"] == "ea":
        return p.ea_kj_mol
    return arrhenius(p.k_ref, p.ea_kj_mol, _t_ref_k(p), measured["temperature_c"])


def evaluate_measurement(dm: Mapping[str, Any], p: Any, coordinate_tier: str, lane: str = "",
                         fitted_sigma_dex: Optional[float] = None) -> dict:
    """One contradiction record: the engine at the measurement's conditions, the size, sec. 7 S2 (b)'s 95 % interval,
    and the route: S2 (fitted-identified coordinates only, every condition (a) - (d) met), else an S4 band or rival
    source, or a satisfied one-sided bound."""
    meas = dm["measured"]
    engine_value = _engine_value_at(p, dm["coordinate"], meas)
    sd, half, band_note = comparison_interval_dex(dm["band"], lane or dm.get("lane", ""), fitted_sigma_dex)
    if dm.get("depends_on_unstated"):
        sd = half = None
        band_note = ("undefined: a condition the step depends on is unstated (" + "; ".join(dm["depends_on_unstated"])
                     + "), so the measurement enters under S4, not S2")
    kind = meas["kind"]
    if kind == "ea":
        lo_ea, hi_ea = float(meas["value"]), float(meas.get("high", meas["value"]))
        nearest = lo_ea if abs(engine_value - lo_ea) <= abs(engine_value - hi_ea) else hi_ea
        size = 0.0 if lo_ea <= engine_value <= hi_ea else ea_dex(engine_value - nearest)
        span = f"{lo_ea:g}" if lo_ea == hi_ea else f"{lo_ea:g}-{hi_ea:g}"
        size_note = f"engine Ea {engine_value:.1f} vs measured {span} kJ/mol, as the log10 k change between 100 and 145 C"
        violated = size != 0.0
    else:
        size = math.log10(engine_value / float(meas["value"]))
        if kind == "upper_bound":
            violated = engine_value > float(meas["value"])
            size_note = f"engine {engine_value:.3g} vs upper bound {meas['value']:g} {meas['unit']} at {meas['temperature_c']:g} C"
        elif kind == "lower_bound":
            violated = engine_value < float(meas["value"])
            size_note = f"engine {engine_value:.3g} vs lower bound {meas['value']:g} {meas['unit']} at {meas['temperature_c']:g} C"
        else:
            violated = True
            size_note = f"engine {engine_value:.3g} vs measured {meas['value']:g} {meas['unit']} at {meas['temperature_c']:g} C"
    z = p_two_sided = None
    if sd and violated:
        z = size / sd
        p_two_sided = math.erfc(abs(z) / math.sqrt(2.0))
    beyond = violated and (half is not None and abs(size) > half)
    mismatched = [k for k, v in dm["matches"].items() if not v]
    route = dm.get("s4_route", "band")
    if kind in ("upper_bound", "lower_bound") and not violated:
        status = "consistent: the one-sided bound is satisfied"
    elif not violated:
        status = "consistent: inside the measured range"
    elif coordinate_tier == "fitted_identified":
        blocks = [f"(a) {m} does not match" for m in mismatched if m != "conditions_stated"]
        if half is None:
            blocks.append("(b) the interval is undefined")
        elif not beyond:
            blocks.append("(b) the 95 % interval includes zero")
        if dm.get("inferred_bound"):
            blocks.append("(d) an inferred bound or a patent: a bound only")
        status = ("S2-eligible, pending the full-lane Holm pass (c) and two readers (d)" if not blocks
                  else "S2 blocked: " + "; ".join(blocks))
    elif coordinate_tier == "fitted_unidentified":
        status = (f"not S2 (the coordinate is fitted_unidentified: no point to contradict; carried under S4); an S4 "
                  f"{route} source" + ("" if half is None else
                                       f"; {'beyond' if beyond else 'inside'} the 95 % interval"))
    else:
        status = (f"contradicted, not fitted (the coordinate is {coordinate_tier}); S2 does not apply; an S4 {route} "
                  "source and a CF entry")
    return {**{k: dm[k] for k in ("id", "measurement", "sources", "ledger_rows", "ledger_note", "cite", "lab", "system",
                                  "channels", "rate_law", "engine_rate_law", "concentration", "unstated", "relation",
                                  "matches", "note")},
            "depends_on_unstated": list(dm.get("depends_on_unstated") or []),
            "inferred_bound": bool(dm.get("inferred_bound")), "s4_route": route,
            "s4_route_reason": dm.get("s4_route_reason", ""),
            "measured": dict(meas), "engine_value": engine_value, "size_dex": size, "size_note": size_note,
            "violated": violated, "interval_sd_dex": sd, "interval95_halfwidth_dex": half, "interval_rule": band_note,
            "z": z, "p_two_sided": p_two_sided, "band_components_dex": dict(dm["band"]), "s2_status": status}


def _sensitivity(key: str, lane: str, sulfur_constants: Mapping[str, dict], b38_sens: Mapping[str, float],
                 wishlist: Mapping[str, dict]) -> dict:
    out: Dict[str, Any] = {"species_of_its_steps": _outputs(key, lane) if lane in ("trunk", "sulfur", "acrylamide") else []}
    basis = []
    if lane == "sulfur" and key in sulfur_constants:
        c = sulfur_constants[key]
        out.update(max_pot_flux_share=c["max_flux_share"], max_odorant_path_share=c["max_odorant_path_share"],
                   matters=c["matters"], steps=c["steps"])
        basis.append("sulfur_reaction_evidence.json (flux attribution, eight reference pots)")
    if key in b38_sens:
        out["b38_sensitivity_norm"] = b38_sens[key]
        basis.append("kinetic_core_b38_identifiability.json (Jacobian column norm on the fit rows)")
    if key in wishlist:
        out["unlocks"] = wishlist[key].get("unlocks_observables") or wishlist[key].get("unlocks_from_charges")
        basis.append("data_wishlist.json sec. 1")
    high = bool(out.get("matters")) or (out.get("b38_sensitivity_norm") or 0.0) >= 1.0
    out["class"] = "high" if high else ("low" if basis else "not on record")
    out["basis"] = basis or ["no sensitivity artifact covers it (the step's products are listed only)"]
    return out


def _extra_rows(priors: List[dict], fits: _Fits, register: List[dict]) -> List[dict]:
    """Band-carried envelope parameters, the lipid lane, the matrix layer."""
    from src.kinetic_core import parameters_lipid as pl, parameters_lipid_b28 as pb, parameters_matrix as pm
    rows: List[dict] = []
    lap_names = {"b8.ph_drift.acid_yield_per_sink_event": ("acid_yield_per_sink_event", "acid_yield"),
                 "b8.ph_drift.arp_secondary_ammonium_pKa": ("arp_secondary_ammonium_pKa", "pka")}
    for prior in priors:
        key = prior["key"]
        base = {"id": key, "constant": key, "lane": prior["lane"], "coordinate": key.split(".")[-1],
                "value": prior.get("centre"), "unit": prior.get("unit"), "evidence_class": prior.get("kind"),
                "band": _band_from_prior(prior, "band"), "bounded_sides": None, "band_missing": False,
                "identification": None, "transformation": "", "reference_temperature_c": None}
        if key in BAND_ROW_TIERS:
            tier, basis, sides = BAND_ROW_TIERS[key]
            rows.append({**base, "tier": tier, "basis": basis, "bounded_sides": sides})
        elif key in lap_names:
            name, kind = lap_names[key]
            lap = fits.lap[(name, kind)]
            lo, hi = lap["bounds"]
            ci = lap["ci95"] if lap["sigma"] is not None else None
            crit = local_unidentified_criteria(value=lap["value"], lo=lo, hi=hi, ci95_halfwidth=ci,
                                               fit_says_identified=lap["identified"], at_bound=lap["at_bound"])
            ident = {"fit": "B60 (free)", "artifact": "kinetic_core_b60_laplace_covariance.json (ph_drift)",
                     "value": lap["value"], "bounds": [lo, hi], "ci95_halfwidth": ci,
                     "labels": {"B60 Laplace": "identified" if lap["identified"] else "unidentified"}, "criteria": crit}
            tier, basis, ident = _fitted_tier(ident)
            rows.append({**base, "tier": tier, "basis": basis, "identification": ident,
                         "label_conflicts": _label_conflicts(ident, prior)})
        elif key.startswith("b18.ph_slope"):
            ident = fits.wave_coord("B18", fits.b18, "pyrazine", key.split(".")[1])
            tier, basis, ident = _fitted_tier(ident)
            rows.append({**base, "tier": tier, "basis": basis, "identification": ident,
                         "label_conflicts": _label_conflicts(ident, prior)})
    # lipid lane
    for lp in (pl.K_HEXANAL_SCHROEN, pl.K_LOOH_DECOMP_ANCHOR, pl.PROPANAL_HEXANAL_RATIO):
        tier = "measured" if lp.evidence_class == "measured_composition" else "analogue"
        basis = (lp.source_anchor[:160] if tier == "measured" else
                 "the source authors' hand-fitted value (no standard error), transferred to cooking temperature by the "
                 "Q10 band (CR-02): " + lp.source_anchor[:120])
        rows.append({"id": f"lipid.{lp.key}", "constant": lp.key, "lane": "lipid", "coordinate": "value",
                     "value": lp.value, "unit": lp.unit, "evidence_class": lp.evidence_class, "tier": tier, "basis": basis,
                     "band": ({"source": "Q10 band 1.6-3.3 (envelope lipid.q10)", "sampled": True} if tier == "analogue" else None),
                     "bounded_sides": None, "band_missing": False, "identification": None, "transformation": lp.what,
                     "reference_temperature_c": lp.temperature_of_measurement_c})
    db = pl.DECADIENAL_BREAKAGE
    for coord, value, unit in (("k_loss_per_min_at_ref", db.k_loss_per_min_at_ref, "1/min at 200 C"),
                               ("ea_loss_kj_mol", db.ea_loss_kj_mol, "kJ/mol"),
                               ("share_hexanal", db.share_hexanal, "mol/mol"), ("share_octenal", db.share_octenal, "mol/mol")):
        rows.append({"id": f"lipid.decadienal_breakage.{coord}", "constant": "DECADIENAL_BREAKAGE", "lane": "lipid",
                     "coordinate": coord, "value": value, "unit": unit, "evidence_class": "measured (B48)",
                     "tier": "measured", "basis": db.source[:200] + " (hold-out: " + db.holdout[:80] + ")", "band": None,
                     "bounded_sides": None, "band_missing": False, "identification": None,
                     "transformation": "2,4-decadienal breakage", "reference_temperature_c": db.t_ref_c})
    simplexes = _json("b6")["frozen_parameters"]["branch_model"]["simplexes"]
    for name, shares in simplexes.items():
        ident = {"fit": "B6", "artifact": "kinetic_core_b6_fit_report.json fit_quality (18 values, 12 free, 3 dof)",
                 "value": None, "bounds": [0.0, 1.0], "ci95_halfwidth": None, "labels": {},
                 "criteria": ["no uncertainty reported (B6 records residuals only)"]}
        tier, basis, ident = _fitted_tier(ident)
        rows.append({"id": f"lipid.branch_simplex.{name}", "constant": "B6 branch simplex", "lane": "lipid",
                     "coordinate": name, "value": shares, "unit": "fraction", "evidence_class": "derived_from_fit_data",
                     "tier": tier, "basis": "fitted to Frankel 1989's zero-additive columns: " + basis, "band": None,
                     "bounded_sides": None, "band_missing": False, "identification": ident,
                     "transformation": "hydroperoxide branch shares", "reference_temperature_c": 180.0})
    rows.append({"id": "lipid.oleate_molar_anchor", "constant": "OLEATE_MOLAR_ANCHOR_BAND", "lane": "lipid",
                 "coordinate": "band", "value": pb.OLEATE_MOLAR_ANCHOR_CENTRE, "unit": "mol/mol",
                 "evidence_class": "declared band (B28)", "tier": "bounded",
                 "basis": "B28: a declared two-sided band 0.2-1.0 (parameters_lipid_b28.py)", "band": {"lo": 0.2, "hi": 1.0},
                 "bounded_sides": "two-sided", "band_missing": False, "identification": None, "transformation": "",
                 "reference_temperature_c": None})
    # matrix layer
    seen = set()
    for group in (pm.REVERSIBLE_BINDING, pm.ALPHA_BETA_UNSATURATION_OBSERVATIONS):
        for mp in group:
            if mp.key in seen:
                continue
            seen.add(mp.key)
            quarantined = "QUARANTINED" in mp.declaration_role.upper()
            tier = "withdrawn" if quarantined else "measured"
            basis = (f"quarantined: {mp.declaration_role}" if quarantined else
                     f"{mp.method} in {mp.medium} ({mp.temperature_c} C): {mp.source_anchor[:120]}")
            rows.append({"id": f"matrix.{mp.key}", "constant": mp.key, "lane": "matrix", "coordinate": "value",
                         "value": mp.value, "unit": mp.unit, "evidence_class": mp.evidence_class, "tier": tier,
                         "basis": basis, "band": None, "bounded_sides": None, "band_missing": False, "identification": None,
                         "transformation": mp.quantity, "reference_temperature_c": mp.temperature_c})
    for r in rows:
        r.setdefault("conflicts_named", _conflict_hits(r["constant"], register) if r["lane"] != "matrix" else [])
        r.setdefault("contradicted_by", [])
        r.setdefault("label_conflicts", [])
        r.setdefault("fit_rows_withdrawn_or_corrected", [])
        r.setdefault("sensitivity", {"class": "not on record", "basis": ["no sensitivity artifact covers it"]})
    return rows


def build_record() -> dict:
    """The full record: one row per coordinate, the summaries the S2 and S4 rules need."""
    fits = _Fits()
    envelope = _json("envelope")
    register = _yaml("conflicts").get("entries") or []
    sulfur_constants = _json("sulfur_evidence")["constants"]
    wishlist = {c["key"]: c for c in _json("wishlist")["unidentified_coordinates"]}
    unpinned = {e["key"] for e in _yaml("unpinned").get("entries") or []}
    b38_sens = {}
    for (wave, name), c in fits.b38.items():
        if wave in ("B3", "B18", "B20", "B21"):
            m = re.match(r"(?:log10_)?(k_[a-z0-9_]+?)(?:_aqueous)?_100C$", name)
            k = m.group(1) if m else name
            if not k.startswith("k_"):
                continue
            b38_sens[k] = max(b38_sens.get(k, 0.0), c["sensitivity_norm"])
    for fit in _json("b38")["fits"]:  # B8's norms describe the superseded sulfur vector; kept only where B60 has none
        if fit["fit"].startswith("B8"):
            for c in fit["coordinates"]:
                if c["key"].startswith("k_"):
                    b38_sens.setdefault(c["key"], c["sensitivity_norm"])

    parameters = shipped_parameters()
    pobjs = {k: p for k, (_, p) in parameters.items()}
    priors: Dict[Tuple[str, str], dict] = {}
    unmapped_priors: List[dict] = []
    for prior in envelope["priors"]:
        targets = _prior_targets(prior["key"], pobjs)
        for t in targets:
            # B55's measured-class band supersedes B8's fixed lumped prior for the barrier it samples
            if t not in priors or prior["key"].startswith("b55."):
                priors[t] = prior
        if not targets:
            unmapped_priors.append(prior)

    dm_by_coord: Dict[Tuple[str, str], List[dict]] = {}
    for dm in DIRECT_MEASUREMENTS:
        dm_by_coord.setdefault((dm["constant"], dm["coordinate"]), []).append(dm)
    withdrawn = [dict(w) for w in WITHDRAWN_OR_CORRECTED_FIT_ROWS]
    for w in withdrawn:
        if w["fit"] == "B3":
            w["row_present_in_frozen_report"] = w["row"] in fits.b3_rows
    withdrawn_by_fit: Dict[str, List[str]] = {}
    for w in withdrawn:
        withdrawn_by_fit.setdefault(w["fit"], []).append(f"{w['fit']}: {w['row']} ({w['kind']}; {w['cite'].split(' (')[0]})")

    rows: List[dict] = []
    for key, (lane, p) in parameters.items():
        hits = _conflict_hits(key, register)
        sens = _sensitivity(key, lane, sulfur_constants, b38_sens, wishlist)
        for coordinate in COORDINATES:
            row = _coordinate(p, lane, coordinate, fits, priors)
            row["dossiers"] = _dossiers(p)
            row["conflicts_named"] = hits
            row["sensitivity"] = sens
            row["label_conflicts"] = _label_conflicts(row["identification"], priors.get((key, coordinate)))
            ci = (row["identification"] or {}).get("ci95_halfwidth")
            fitted_sd = None if ci is None else (ci / 1.96) * (DEX_PER_KJ_100_145 if coordinate == "ea_kj_mol" else 1.0)
            row["contradicted_by"] = [evaluate_measurement(dm, p, row["tier"], lane, fitted_sd)
                                      for dm in dm_by_coord.get((key, coordinate), [])]
            fit = (row["identification"] or {}).get("fit", "")
            row["fit_rows_withdrawn_or_corrected"] = [
                w for wave, ws in withdrawn_by_fit.items() for w in ws if fit == wave or fit.startswith(wave + " ")]
            rows.append(row)
    rows += _extra_rows(unmapped_priors, fits, register)
    note_fit_drift(rows, _yaml("frozen_fit_drift").get("entries") or [])
    changes: Dict[str, List[dict]] = {}
    for entry in tier_changes():
        changes.setdefault(str(entry.get("coordinate")), []).append(dict(entry))
    for row in rows:
        row["tier_changes"] = changes.pop(row["id"], [])
    return {"rows": rows, "summary": _summary(rows, wishlist, unpinned, sulfur_constants),
            "vocabulary_map": VOCABULARY_MAP, "structural_choices": list(STRUCTURAL_CHOICES),
            "withdrawn_or_corrected_fit_rows": withdrawn,
            "tier_changes_without_a_row": [e for es in changes.values() for e in es],
            "rules": {"tiers": list(TIERS), "edge_fraction": EDGE_FRACTION, "dex_per_kj_mol_100_145C": DEX_PER_KJ_100_145,
                      "comparison_interval": "sec. 7 S2 (b): 95 % interval of (measurement - fitted value) with the "
                                             "reported noise, the within-lab between-paper spread, the between-lab "
                                             "term (the step's measured spread where >= 2 laboratories measured it, "
                                             "else the lane's) and the fitted coordinate's own interval in quadrature; "
                                             "undefined when a condition the step depends on is unstated (S4, not S2)",
                      "lane_between_lab_dex": LANE_BETWEEN_LAB_DEX, "lane_between_lab_source": LANE_BETWEEN_LAB_SOURCE,
                      "s2_trigger": "sec. 7 S2: a fitted-IDENTIFIED coordinate only; (a) the measurement is of the "
                                    "step (all lumped channels, the engine's rate law and concentration regime); "
                                    "(b) the 95 % interval excludes zero; (c) a full-lane Holm pass with every "
                                    "measurement of the step; (d) two readers, OpenAlex verify, not a quote of a quote "
                                    "nor an inferred bound (patents are bounds only). A contradiction of a "
                                    "fitted-unidentified coordinate is an S4 band or rival source",
                      "s4_band": "a measured bound, a same-class analogue, or a reaction-class distribution; a "
                                 "search bound is never a band",
                      "s4_high_sensitivity": "the step matters in the sulfur flux attribution (> 1 % of a pot's flux or "
                                             "of an odorant-path species), or B38's Jacobian column norm >= 1"},
            "inputs": [str(p) for p in INPUTS.values()]}


def note_fit_drift(rows: List[dict], allowlist: Iterable[Mapping[str, Any]]) -> int:
    """Sec. 7 S1: a coordinate whose fit no longer reproduces at its frozen point says so in its identification (the
    allow-list's MATERIAL entries, which the drift gate holds to today's engine). Text only: no tier moves.
    Returns the number of coordinates noted."""
    material = {e["fit"]: e for e in allowlist if e.get("class") == "MATERIAL"}
    n = 0
    for row in rows:
        ident = row.get("identification")
        entry = material.get(FIT_OF_IDENTIFICATION_LABEL.get((ident or {}).get("fit")))
        if entry is None:
            continue
        cite = entry.get("cf") or f"PR #{entry.get('pr')}"
        note = (f"{entry['fit']}'s objective no longer reproduces at its frozen point (S = r . r {float(entry['S']):.6g} "
                f"today; MATERIAL drift, {cite}; frozen_fit_reproducibility.md): no current objective supports this value")
        ident["fit_reproduces"] = False
        ident["fit_drift"] = {"fit": entry["fit"], "class": "MATERIAL", "S_today": float(entry["S"]),
                              "cf": entry.get("cf"), "pr": entry.get("pr"), "note": note,
                              "record": "results/validation/frozen_fit_reproducibility.json"}
        row["basis"] = f"{note}. {row['basis']}"
        n += 1
    return n


def _dossiers(p: Any) -> List[str]:
    from src import evidence_anchors
    return evidence_anchors.dossier_stems(getattr(p, "dossier_anchor", ""))


def _summary(rows: List[dict], wishlist: Mapping[str, dict], unpinned: set, sulfur_constants: Mapping[str, dict]) -> dict:
    lanes = sorted({r["lane"] for r in rows})
    by_lane = {lane: {t: sum(1 for r in rows if r["lane"] == lane and r["tier"] == t) for t in TIERS + ("unclassified",)}
               for lane in lanes}
    by_lane = {lane: {t: n for t, n in d.items() if n} for lane, d in by_lane.items()}
    by_coord: Dict[str, Dict[str, int]] = {}
    for r in rows:
        by_coord.setdefault(r["coordinate"] if r["coordinate"] in ("log10_k_ref", "ea_kj_mol") else "band / other", {})
        d = by_coord[r["coordinate"] if r["coordinate"] in ("log10_k_ref", "ea_kj_mol") else "band / other"]
        d[r["tier"]] = d.get(r["tier"], 0) + 1
    totals = {t: sum(1 for r in rows if r["tier"] == t) for t in TIERS + ("unclassified",)}
    fitted = [r for r in rows if r["tier"] in FITTED_TIERS]
    status = {}
    for r in fitted:
        s = (r["identification"] or {}).get("status", "unknown")
        status[s] = status.get(s, 0) + 1
    # S2 (sec. 7, final): fitted-IDENTIFIED coordinates contradicted by a direct measurement. Contradictions of
    # fitted-unidentified coordinates are S4 band / rival sources (no point to contradict).
    s2, s4_sources, contradicted_other = [], [], []
    for r in rows:
        for c in r["contradicted_by"]:
            if not c["violated"]:
                continue
            entry = {"coordinate": r["id"], "lane": r["lane"], "tier": r["tier"], "measurement": c["id"] + ": " + c["measurement"],
                     "sources": c["sources"], "size_dex": c["size_dex"], "size_note": c["size_note"],
                     "interval95_halfwidth_dex": c["interval95_halfwidth_dex"], "interval_rule": c["interval_rule"],
                     "z": c["z"], "p_two_sided": c["p_two_sided"], "s4_route": c["s4_route"],
                     "s2_status": c["s2_status"], "outputs": r["sensitivity"].get("species_of_its_steps"),
                     "sensitivity": r["sensitivity"]["class"],
                     "odorant_path_share": r["sensitivity"].get("max_odorant_path_share")}
            if r["tier"] == "fitted_identified":
                s2.append(entry)
            elif r["tier"] == "fitted_unidentified":
                s4_sources.append(entry)
            else:
                contradicted_other.append(entry)
    rank = lambda e: (e["sensitivity"] != "high", -abs(e["size_dex"]))  # noqa: E731
    s2.sort(key=rank)
    s4_sources.sort(key=rank)
    bounds_satisfied = [{"coordinate": r["id"], "measurement": c["id"] + ": " + c["measurement"], "status": c["s2_status"]}
                        for r in rows for c in r["contradicted_by"] if c["s2_status"].startswith("consistent")]
    # S4: fitted-unidentified (local) coordinates with high output sensitivity, ranked by band width
    s4 = []
    for r in rows:
        if r["tier"] != "fitted_unidentified" or (r["identification"] or {}).get("status") != "unidentified":
            continue
        if r["sensitivity"]["class"] != "high":
            continue
        band = r.get("band") or {}
        legit = None if r.get("band_is_search_bound") else band.get("width_dex")
        span = legit if legit is not None else r.get("search_bounds_span_dex")
        sens = r["sensitivity"]
        weight = max(sens.get("max_odorant_path_share") or 0.0, sens.get("max_pot_flux_share") or 0.0,
                     min((sens.get("b38_sensitivity_norm") or 0.0) / 10.0, 1.0))
        s4.append({"impact_score": (span or 0.0) * weight, "coordinate": r["id"], "lane": r["lane"],
                   "why_unidentified": r["identification"]["criteria"],
                   "band_width_dex": legit, "search_bounds_span_dex": r.get("search_bounds_span_dex"),
                   "s4_band_source": ("a search bound (never a band under S4): needs a measured bound, an analogue or a "
                                      "reaction-class distribution, else the output is refused") if legit is None
                   else band.get("source"),
                   "envelope_samples_it": bool(band.get("sampled")),
                   "max_odorant_path_share": r["sensitivity"].get("max_odorant_path_share"),
                   "max_pot_flux_share": r["sensitivity"].get("max_pot_flux_share"),
                   "b38_sensitivity_norm": r["sensitivity"].get("b38_sensitivity_norm"),
                   "outputs": r["sensitivity"].get("species_of_its_steps")})
    #: impact = span (dex: the band, else the search span as a stand-in for ranking only) x the step's largest share of
    #: a pot's flux or of an odorant-path species (or B38's norm / 10, capped at 1, where no flux record exists);
    #: coordinates with no legitimate band come first, then those the envelope holds fixed.
    s4.sort(key=lambda e: (e["band_width_dex"] is not None, e["envelope_samples_it"], -e["impact_score"]))
    label_conflicts = [{"coordinate": r["id"], **lc} for r in rows for lc in r["label_conflicts"]]
    # the two existing unidentified lists, reconciled
    wl = set(wishlist)
    record_unid = {r["constant"] for r in rows if r["lane"] == "sulfur" and r["coordinate"] == "log10_k_ref"
                   and (r["identification"] or {}).get("status") == "unidentified"}
    sulfur_unpinned_matter = {k for k, c in sulfur_constants.items() if c.get("unpinned_and_matters")}
    wl_keys = {k.replace("Ea_decay_", "") for k in wl}
    reconciliation = {
        "data_wishlist_sec1": sorted(wl),
        "unpinned_constants_yml": sorted(unpinned),
        "sulfur_evidence_unpinned_that_matter": sorted(sulfur_unpinned_matter),
        "in_both_lists": sorted(wl & unpinned),
        "why_disjoint": "the wishlist reads B60's slice profiles, so it lists only coordinates B60 FREED (bound-limited or "
                        "flat); the register lists constants frozen since B2.3 whose rate the envelope does not sample. "
                        "Neither covers the other; this record covers both, per coordinate.",
        "unidentified_here_on_neither_list": sorted(record_unid - wl_keys - unpinned - sulfur_unpinned_matter),
        "on_a_list_but_not_unidentified_here": sorted((wl_keys | unpinned) - record_unid -
                                                      {"thiol_sink", "ph_acid_yield_per_sink_event"}),
    }
    unclassified = [{"coordinate": r["id"], "basis": r["basis"]} for r in rows if r["tier"] == "unclassified"]
    analogue_no_band = [r["id"] for r in rows if r["tier"] == "analogue" and r.get("band_missing")]
    return {
        "coordinates": len(rows), "by_tier": {t: n for t, n in totals.items() if n}, "by_lane": by_lane,
        "by_coordinate_kind": by_coord, "fitted_identification_status": status,
        "s2_candidates": s2, "s4_sources_contradicting_unidentified": s4_sources,
        "contradicted_non_fitted": contradicted_other, "one_sided_bounds_satisfied": bounds_satisfied,
        "s4_sensitive_unidentified": s4, "label_conflicts": label_conflicts,
        "unidentified_lists_reconciled": reconciliation, "unclassified": unclassified,
        "analogue_without_a_declared_band": analogue_no_band,
        "structural_choices_without_a_case": [s["id"] for s in STRUCTURAL_CHOICES if s["case"].startswith("none")],
    }


__all__ = ["TIERS", "DIRECT_MEASUREMENTS", "STRUCTURAL_CHOICES", "WITHDRAWN_OR_CORRECTED_FIT_ROWS", "VOCABULARY_MAP",
           "arrhenius", "ea_dex", "comparison_interval_dex", "local_unidentified_criteria", "evaluate_measurement",
           "build_record", "INPUTS"]
