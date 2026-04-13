from __future__ import annotations

from typing import Any, Dict, Mapping


def build_kokumi_support_profile(
    *,
    glutathione_active: bool,
    gamma_glutamyl_peptide_active: bool,
    peptide_accessibility_factor: float,
    kokumi_priors: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    priors = dict(kokumi_priors or {})
    bounded_accessibility = max(0.0, min(1.0, float(peptide_accessibility_factor)))

    msg_baseline_mouthfulness = float(priors.get("msg_baseline_mouthfulness_score", 3.2) or 3.2)
    gsh_mouthfulness = float(priors.get("gsh_mouthfulness_score", msg_baseline_mouthfulness) or msg_baseline_mouthfulness)
    gamma_glu_val_mouthfulness = float(
        priors.get("gamma_glu_val_mouthfulness_score", msg_baseline_mouthfulness) or msg_baseline_mouthfulness
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
        selected_kokumi_delta = min(
            gamma_glu_val_mouthfulness_delta,
            selected_kokumi_delta + 0.25 * gsh_mouthfulness_delta,
        )

    kokumi_support_active = bool(glutathione_active or gamma_glutamyl_peptide_active)
    kokumi_signal_score = 0.0
    if kokumi_support_active:
        kokumi_signal_score = max(
            0.0,
            min(1.0, selected_kokumi_delta / max(gamma_glu_val_mouthfulness_delta, 1.0e-6)),
        )

    kokumi_signal_mode = "inactive"
    if glutathione_active and gamma_glutamyl_peptide_active:
        kokumi_signal_mode = "combined_glutathione_and_gamma_glutamyl_peptide"
    elif gamma_glutamyl_peptide_active:
        kokumi_signal_mode = "gamma_glutamyl_peptide_supported"
    elif glutathione_active:
        kokumi_signal_mode = "glutathione_reference"

    return {
        "kokumi_support_active": bool(kokumi_support_active),
        "kokumi_signal_mode": kokumi_signal_mode,
        "kokumi_signal_score": float(kokumi_signal_score),
        "selected_kokumi_delta": float(selected_kokumi_delta),
        "msg_baseline_mouthfulness_score": float(msg_baseline_mouthfulness),
        "gsh_mouthfulness_scale_increase": float(gsh_mouthfulness_delta),
        "gamma_glu_val_mouthfulness_scale_increase": float(gamma_glu_val_mouthfulness_delta),
        "gsh_casr_ec50_mM": float(priors.get("gsh_ec50_mM", 0.68) or 0.68),
        "gamma_glu_val_casr_ec50_mM": float(priors.get("gamma_glu_val_ec50_mM", 0.32) or 0.32),
        "gsh_cck_secretion_increase_percent": float(priors.get("gsh_cck_secretion_increase_percent", 141.0) or 141.0),
    }