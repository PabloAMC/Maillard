"""B63 (2026-10-28): single-temperature sulfur constants take their measured class's median barrier."""
import dataclasses

import numpy as np

from src.kinetic_core import ea_class_prior
from src.kinetic_core import parameters_sulfur as ps


def test_route_medians_are_the_class_priors_computed_medians():
    values = ea_class_prior.class_values()
    for route, klass in (("sugar_trunk", ea_class_prior.CARBOHYDRATE), ("thiol_assembly", ea_class_prior.SULFUR_CLASS)):
        assert ps.ROUTE_CLASS_MEDIAN_EA_KJ_MOL[route] == float(np.median(values[klass])), route


def test_switch_moves_only_route_sharing_barriers(monkeypatch):
    fitted = {"k_pent_dpo": -1.0, "k_nf_mft": -2.0, "k_cys_thermal": -3.0, "k_fur_decay": -0.4}
    off = ps.with_fitted_sulfur(fitted, 64.08, decay_ea={"carbonyl_sink": 143.5, "thiol_sink": 101.3})
    monkeypatch.setattr(ps, "SINGLE_TEMPERATURE_BARRIERS", dataclasses.replace(ps.SINGLE_TEMPERATURE_BARRIERS, enabled=True))
    on = ps.with_fitted_sulfur(fitted, 64.08, decay_ea={"carbonyl_sink": 143.5, "thiol_sink": 101.3})
    assert off["k_pent_dpo"].ea_kj_mol == 64.08 and on["k_pent_dpo"].ea_kj_mol == 103.7   # sugar trunk
    assert off["k_nf_mft"].ea_kj_mol == 64.08 and on["k_nf_mft"].ea_kj_mol == 100.9       # thiol assembly
    assert on["k_cys_thermal"].ea_kj_mol == off["k_cys_thermal"].ea_kj_mol               # measured override
    assert on["k_fur_decay"].ea_kj_mol == off["k_fur_decay"].ea_kj_mol == 143.5           # decay family
    for key in fitted:
        assert on[key].k_ref == off[key].k_ref                                            # k(145 C) held
