import pytest

from src.kokumi_scoring import build_kokumi_support_profile


def test_build_kokumi_support_profile_is_inactive_without_gsh_or_gamma_glu_peptides():
    profile = build_kokumi_support_profile(
        glutathione_active=False,
        gamma_glutamyl_peptide_active=False,
        peptide_accessibility_factor=0.8,
        kokumi_priors={},
    )

    assert profile["kokumi_support_active"] is False
    assert profile["kokumi_signal_mode"] == "inactive"
    assert profile["kokumi_signal_score"] == pytest.approx(0.0)


def test_build_kokumi_support_profile_reweights_gamma_glu_with_accessibility():
    profile = build_kokumi_support_profile(
        glutathione_active=False,
        gamma_glutamyl_peptide_active=True,
        peptide_accessibility_factor=0.5,
        kokumi_priors={
            "msg_baseline_mouthfulness_score": 3.2,
            "gamma_glu_val_mouthfulness_score": 6.8,
        },
    )

    assert profile["kokumi_support_active"] is True
    assert profile["kokumi_signal_mode"] == "gamma_glutamyl_peptide_supported"
    assert profile["selected_kokumi_delta"] == pytest.approx(1.8)
    assert profile["kokumi_signal_score"] == pytest.approx(0.5)