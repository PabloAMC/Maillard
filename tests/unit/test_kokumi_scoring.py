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


# ---------------------------------------------------------------------------
# Audit remediation (wave G5, module-hunt TIER 4) regression tests.
# ---------------------------------------------------------------------------

_PRIORS = {
    "msg_baseline_mouthfulness_score": 3.2,
    "gsh_mouthfulness_score": 5.4,
    "gamma_glu_val_mouthfulness_score": 6.8,
}


def _score(*, gsh, peptide, accessibility, priors=None):
    return build_kokumi_support_profile(
        glutathione_active=gsh,
        gamma_glutamyl_peptide_active=peptide,
        peptide_accessibility_factor=accessibility,
        kokumi_priors=dict(priors if priors is not None else _PRIORS),
    )["kokumi_signal_score"]


def test_kokumi_signal_is_monotone_non_decreasing_in_accessibility():
    grid = [0.0, 0.05, 0.1, 0.25, 0.4, 0.5, 0.7, 0.9, 1.0]
    combined = [_score(gsh=True, peptide=True, accessibility=a) for a in grid]
    assert combined == sorted(combined)

    peptide_only = [_score(gsh=False, peptide=True, accessibility=a) for a in grid]
    assert peptide_only == sorted(peptide_only)


def test_fully_inaccessible_peptide_cannot_raise_the_kokumi_score():
    """
    The 0.25 * gsh_delta combination bonus is a synergy term and is now gated on
    accessibility. Previously an entirely inaccessible peptide still bought the bonus, so
    declaring a peptide the matrix cannot release scored HIGHER than glutathione alone.
    """
    gsh_only = _score(gsh=True, peptide=False, accessibility=0.0)
    inaccessible_peptide = _score(gsh=True, peptide=True, accessibility=0.0)

    assert inaccessible_peptide == pytest.approx(gsh_only)
    assert inaccessible_peptide <= gsh_only + 1e-12

    # Adding accessibility can only help.
    assert _score(gsh=True, peptide=True, accessibility=1.0) > gsh_only


def test_undefined_reference_scale_collapses_to_a_defined_neutral_score():
    """
    When the gamma-glutamyl reference uplift is zero the signal scale is undefined. The
    old max(delta, 1e-6) denominator inverted that into the maximal score of 1.0.
    """
    profile = build_kokumi_support_profile(
        glutathione_active=True,
        gamma_glutamyl_peptide_active=False,
        peptide_accessibility_factor=1.0,
        kokumi_priors={
            "msg_baseline_mouthfulness_score": 3.2,
            "gsh_mouthfulness_score": 5.4,
            "gamma_glu_val_mouthfulness_score": 3.2,  # no uplift over the MSG baseline
        },
    )

    assert profile["kokumi_reference_scale_defined"] is False
    assert profile["kokumi_signal_score"] == pytest.approx(0.0)
    assert profile["kokumi_support_active"] is True


def test_legitimate_zero_prior_is_not_replaced_by_the_default():
    """The `x or default` idiom silently overwrote a stored 0.0 with the fallback."""
    profile = build_kokumi_support_profile(
        glutathione_active=True,
        gamma_glutamyl_peptide_active=True,
        peptide_accessibility_factor=1.0,
        kokumi_priors={
            "msg_baseline_mouthfulness_score": 0.0,
            "gsh_mouthfulness_score": 0.0,
            "gamma_glu_val_mouthfulness_score": 2.0,
        },
    )

    assert profile["msg_baseline_mouthfulness_score"] == pytest.approx(0.0)
    assert profile["gsh_mouthfulness_scale_increase"] == pytest.approx(0.0)
    assert profile["gamma_glu_val_mouthfulness_scale_increase"] == pytest.approx(2.0)


def test_casr_ec50_echo_matches_verified_ohsu_2010_values():
    """
    Echo-only constants (zero prediction impact, audited). The shipped 0.68 / 0.32 mM
    values were attributed to a non-existent "Ohsu et al. (2025)"; Ohsu et al. (2010),
    J. Biol. Chem. 285:1016 reports 76.5 nM (GSH), 41.9 nM (gamma-Glu-Val-Gly) and
    1.34 uM (gamma-Glu-Val). The CCK figure appears nowhere and is withdrawn.
    """
    profile = build_kokumi_support_profile(
        glutathione_active=True,
        gamma_glutamyl_peptide_active=True,
        peptide_accessibility_factor=0.7,
        kokumi_priors=dict(_PRIORS),
    )

    assert profile["gsh_casr_ec50_nM"] == pytest.approx(76.5)
    assert profile["gamma_glu_val_gly_casr_ec50_nM"] == pytest.approx(41.9)
    assert profile["gamma_glu_val_casr_ec50_nM"] == pytest.approx(1340.0)
    assert profile["gsh_casr_ec50_mM"] == pytest.approx(76.5e-6)
    assert "Ohsu" in profile["casr_ec50_provenance"]

    assert profile["gsh_cck_secretion_increase_percent"] is None
    assert "withdrawn" in profile["gsh_cck_secretion_increase_status"]


def test_unanchored_mM_ec50_priors_do_not_leak_back_into_the_echo():
    """The stale payload numbers must not be able to re-enter through the priors mapping."""
    profile = build_kokumi_support_profile(
        glutathione_active=True,
        gamma_glutamyl_peptide_active=True,
        peptide_accessibility_factor=0.7,
        kokumi_priors={
            **_PRIORS,
            "gsh_ec50_mM": 0.68,
            "gamma_glu_val_ec50_mM": 0.32,
            "gsh_cck_secretion_increase_percent": 141.0,
        },
    )

    assert profile["gsh_casr_ec50_mM"] == pytest.approx(76.5e-6)
    assert profile["gsh_casr_ec50_mM"] != pytest.approx(0.68)
    assert profile["gsh_cck_secretion_increase_percent"] is None