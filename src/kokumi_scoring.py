from __future__ import annotations

from typing import Any, Dict, Mapping, Optional


# ---------------------------------------------------------------------------
# Verified CaSR EC50 constants (echo-only; ZERO prediction impact -- audited).
#
# The values previously echoed here (GSH 0.68 mM, gamma-Glu-Val 0.32 mM, and a
# "+141% CCK secretion" figure) were attributed to "Ohsu et al. (2025)", a paper that does
# not exist; the DOI stored alongside it belongs to Yang et al. (2022). The real CaSR
# kokumi paper is:
#
#   Ohsu, Amino, Nagasaki, Yamanaka et al. (2010), "Involvement of the Calcium-sensing
#   Receptor in Human Taste Perception", J. Biol. Chem. 285:1016.
#
# which reports CaSR EC50 values of 76.5 nM (glutathione), 41.9 nM (gamma-Glu-Val-Gly)
# and 1.34 uM (gamma-Glu-Val) -- three to four orders of magnitude from the numbers that
# were stored. It reports NO CCK-secretion percentages at all, so that figure is
# withdrawn rather than corrected. These constants are reported for reference only: they
# feed no calculation in this module (verified -- the only consumers of the EC50 output
# keys are tests), so correcting them moves no prediction.
# ---------------------------------------------------------------------------
GSH_CASR_EC50_NM = 76.5
GAMMA_GLU_VAL_CASR_EC50_NM = 1340.0
GAMMA_GLU_VAL_GLY_CASR_EC50_NM = 41.9

_EC50_PROVENANCE = (
    "Ohsu, Amino, Nagasaki, Yamanaka et al. (2010), J. Biol. Chem. 285:1016 "
    "(CaSR EC50: GSH 76.5 nM, gamma-Glu-Val-Gly 41.9 nM, gamma-Glu-Val 1.34 uM)"
)
_CCK_WITHDRAWN_NOTE = (
    "withdrawn 2026-08-27: the previously echoed +141% CCK-secretion increase is not "
    "reported in Ohsu et al. (2010) or in any verifiable source; no replacement value"
)

_NM_PER_MM = 1.0e6


def _prior_float(priors: Mapping[str, Any], key: str, default: float) -> float:
    """
    Read a float prior, falling back to `default` ONLY when the key is absent or null.

    The previous idiom was `float(priors.get(key, default) or default)`, which also
    replaced a legitimate stored 0.0 (a real "no uplift over the MSG baseline" datum)
    with the default.
    """
    value = priors.get(key, None)
    if value is None:
        return float(default)
    return float(value)


def build_kokumi_support_profile(
    *,
    glutathione_active: bool,
    gamma_glutamyl_peptide_active: bool,
    peptide_accessibility_factor: float,
    kokumi_priors: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    priors = dict(kokumi_priors or {})
    bounded_accessibility = max(0.0, min(1.0, float(peptide_accessibility_factor)))

    msg_baseline_mouthfulness = _prior_float(priors, "msg_baseline_mouthfulness_score", 3.2)
    gsh_mouthfulness = _prior_float(priors, "gsh_mouthfulness_score", msg_baseline_mouthfulness)
    gamma_glu_val_mouthfulness = _prior_float(
        priors, "gamma_glu_val_mouthfulness_score", msg_baseline_mouthfulness
    )
    gsh_mouthfulness_delta = max(0.0, gsh_mouthfulness - msg_baseline_mouthfulness)
    gamma_glu_val_mouthfulness_delta = max(0.0, gamma_glu_val_mouthfulness - msg_baseline_mouthfulness)

    kokumi_delta_candidates = []
    if glutathione_active:
        kokumi_delta_candidates.append(gsh_mouthfulness_delta)
    if gamma_glutamyl_peptide_active:
        kokumi_delta_candidates.append(gamma_glu_val_mouthfulness_delta * bounded_accessibility)

    selected_kokumi_delta = max(kokumi_delta_candidates or [0.0])
    if glutathione_active and gamma_glutamyl_peptide_active:
        # The combined GSH + gamma-glutamyl-peptide bonus is a SYNERGY term, so it is
        # gated on peptide accessibility exactly like the peptide's own contribution.
        # Ungated, an entirely inaccessible peptide (accessibility 0) still bought a
        # +0.25 * gsh_delta bonus, so declaring a peptide the matrix cannot release
        # RAISED the score above the glutathione-only case. With the gate, accessibility 0
        # reduces exactly to the glutathione-only result and the score is monotone
        # non-decreasing in accessibility.
        selected_kokumi_delta = min(
            gamma_glu_val_mouthfulness_delta,
            selected_kokumi_delta + 0.25 * gsh_mouthfulness_delta * bounded_accessibility,
        )

    kokumi_support_active = bool(glutathione_active or gamma_glutamyl_peptide_active)

    # The gamma-glutamyl-peptide uplift is the reference scale the signal is expressed in.
    # When it is zero the scale is undefined -- and the old `max(delta, 1e-6)` denominator
    # inverted that into the STRONGEST possible answer (any positive numerator over 1e-6
    # clipped to 1.0), so "we have no reference uplift" scored 1.0. Guard it explicitly
    # and report the defined neutral score of 0.0 instead.
    reference_scale_defined = gamma_glu_val_mouthfulness_delta > 0.0
    kokumi_signal_score = 0.0
    if kokumi_support_active and reference_scale_defined:
        kokumi_signal_score = max(
            0.0,
            min(1.0, selected_kokumi_delta / gamma_glu_val_mouthfulness_delta),
        )

    kokumi_signal_mode = "inactive"
    if glutathione_active and gamma_glutamyl_peptide_active:
        kokumi_signal_mode = "combined_glutathione_and_gamma_glutamyl_peptide"
    elif gamma_glutamyl_peptide_active:
        kokumi_signal_mode = "gamma_glutamyl_peptide_supported"
    elif glutathione_active:
        kokumi_signal_mode = "glutathione_reference"

    gsh_ec50_nM = _prior_float(priors, "gsh_ec50_nM", GSH_CASR_EC50_NM)
    gamma_glu_val_ec50_nM = _prior_float(
        priors, "gamma_glu_val_ec50_nM", GAMMA_GLU_VAL_CASR_EC50_NM
    )
    gamma_glu_val_gly_ec50_nM = _prior_float(
        priors, "gamma_glu_val_gly_ec50_nM", GAMMA_GLU_VAL_GLY_CASR_EC50_NM
    )

    cck_increase: Optional[float] = priors.get("gsh_cck_secretion_increase_percent_verified", None)

    return {
        "kokumi_support_active": bool(kokumi_support_active),
        "kokumi_signal_mode": kokumi_signal_mode,
        "kokumi_signal_score": float(kokumi_signal_score),
        "kokumi_reference_scale_defined": bool(reference_scale_defined),
        "selected_kokumi_delta": float(selected_kokumi_delta),
        "peptide_accessibility_factor": float(bounded_accessibility),
        "msg_baseline_mouthfulness_score": float(msg_baseline_mouthfulness),
        "gsh_mouthfulness_scale_increase": float(gsh_mouthfulness_delta),
        "gamma_glu_val_mouthfulness_scale_increase": float(gamma_glu_val_mouthfulness_delta),
        # Echo-only reference constants. The nM keys are the verified primary form; the
        # mM keys are retained for backwards compatibility and derived from them, NOT
        # from the withdrawn 0.68 / 0.32 mM values.
        "gsh_casr_ec50_nM": float(gsh_ec50_nM),
        "gamma_glu_val_casr_ec50_nM": float(gamma_glu_val_ec50_nM),
        "gamma_glu_val_gly_casr_ec50_nM": float(gamma_glu_val_gly_ec50_nM),
        "gsh_casr_ec50_mM": float(gsh_ec50_nM) / _NM_PER_MM,
        "gamma_glu_val_casr_ec50_mM": float(gamma_glu_val_ec50_nM) / _NM_PER_MM,
        "casr_ec50_provenance": _EC50_PROVENANCE,
        "gsh_cck_secretion_increase_percent": cck_increase,
        "gsh_cck_secretion_increase_status": _CCK_WITHDRAWN_NOTE,
    }
