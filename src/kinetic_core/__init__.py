"""
src/kinetic_core -- THE MASS-ACTION KINETIC CORE (Build Wave B1, 2026-08-28).

Module 1 of the approved kinetic-core rebuild: the TRUNK NETWORK plus the
MELANOIDIN MASS SINK. Later modules (sulfur, acrylamide, competition) attach to
the species and the parameter registry declared here.

WHAT THIS PACKAGE IS
--------------------
An explicit mass-action ODE network over the glucose/glycine trunk, carried in
the units its source data are printed in (mmol/L, minutes, Kelvin, kJ/mol),
integrated with a stiff-capable solver, with:

  * a formation AND a consumption term for every intermediate;
  * exact carbon and nitrogen conservation, enforced at import time by
    ``network.validate_balance`` and asserted by the unit tests;
  * a terminal melanoidin pool carried ELEMENTALLY;
  * every rate constant carrying its source anchor, units, conditions and
    ``pH_of_measurement``;
  * no DFT-derived number anywhere, by policy and by construction;
  * no pH term and no a_w term, because the corpus licenses neither.

WHAT THIS PACKAGE IS NOT
------------------------
Nothing in the shipped prediction path imports it. It does not modify the FAST
screening lane, ``benchmark_validation``, any governance module, or anything
under ``data/benchmarks/``. Like ``src/trunk_kinetics``, which it supersedes
for calibration purposes, it is a calibration and prediction lane and only
that. ``src/trunk_kinetics`` is left in place, untouched, because the audit
artefacts of wave S3 refer to it.
"""

from .integrate import CoreRun, flux_budget, integrate  # noqa: F401
from .network import (  # noqa: F401
    REACTIONS,
    STOICHIOMETRY,
    Reaction,
    derivatives,
    describe,
    rate_constants_at,
    reaction_rates,
    validate_balance,
)
from .parameters import (  # noqa: F401
    CROSS_LAB_COMPARATORS,
    FITTED_BOUNDS_LOG10K,
    FITTED_EA_BOUNDS,
    FITTED_KEYS,
    MARTINS_M4,
    NETWORK_PH,
    R_KJ,
    SCHIFF_AMADORI_SPLIT,
    T_REF_K,
    KineticParameter,
    assert_no_dft,
    check_ph_homogeneity,
    fitted_placeholders,
    registry_metadata,
    with_fitted_values,
)
from .species import (  # noqa: F401
    INDEX,
    MEASURED_LABEL_TO_KEY,
    N_SPECIES,
    SPECIES,
    SPECIES_KEYS,
    Species,
    initial_state,
    melanoidin_c_over_n,
    melanoidin_repeat_units,
    state_as_dict,
    total_carbon,
    total_nitrogen,
)


# --- Build Wave B2: the sulfur module ------------------------------------
# Imported lazily-by-name rather than star-imported, so that B1's namespace is
# unchanged for every existing caller and the sulfur symbols are additive.
from .parameters_sulfur import (  # noqa: F401,E402
    MEASURED_SULFUR,
    PROHIBITED_DERIVATIONS,
    THIOL_CHANNELS,
    SulfurParameter,
    T_REF_S_K,
    assert_no_dft_sulfur,
    sulfur_placeholders,
    sulfur_registry_metadata,
    with_fitted_sulfur,
)
from .species_sulfur import (  # noqa: F401,E402
    ODOUR_THRESHOLD_UG_PER_L,
    SULFUR_INDEX,
    SULFUR_STATE,
    odour_activity_values,
    total_sulfur_extended,
)
from .sulfur import (  # noqa: F401,E402
    FULL_REACTIONS,
    OUT_OF_SCOPE,
    SULFUR_REACTIONS,
    SulfurRun,
    branch_shares,
    describe_sulfur,
    integrate_sulfur,
    sulfur_flux_budget,
    validate_sulfur_balance,
)


def operative_parameters(fitted):
    """
    Assemble the full operative parameter set: measured backbone + fitted steps.

    ``fitted`` is ``{key: (k_ref_at_100C, Ea_kJ_per_mol)}`` for the four fitted
    consumption steps. There is no default: a network without them refuses to
    integrate.
    """
    parameters = dict(MARTINS_M4)
    parameters.update(with_fitted_values(fitted))
    return parameters


def operative_sulfur_parameters(b1_fitted, sulfur_log10k, lumped_formation_ea):
    """
    The full operative set for the SULFUR network: B1's trunk (measured +
    fitted), the two measured sulfur constants, and the fitted sulfur steps.

    ``sulfur_log10k`` is ``{key: log10 k_ref at 145 C}``. Formation steps
    receive ``lumped_formation_ea``; consumption steps receive NO activation
    energy at all, which is Module 2's central policy and not an omission.
    """
    parameters = dict(operative_parameters(b1_fitted))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(sulfur_log10k, lumped_formation_ea))
    return parameters


__all__ = [
    "CoreRun",
    "MARTINS_M4",
    "MEASURED_SULFUR",
    "REACTIONS",
    "SPECIES",
    "SULFUR_REACTIONS",
    "SULFUR_STATE",
    "branch_shares",
    "integrate",
    "integrate_sulfur",
    "odour_activity_values",
    "operative_parameters",
    "operative_sulfur_parameters",
]
