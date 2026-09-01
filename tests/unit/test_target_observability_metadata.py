from src.chem_utils import Species
from src.recommend import _canon, _headspace_observability_metadata


def test_headspace_observability_metadata_marks_hmf_as_observable():
    # RENAMED + INVERTED 2026-08-27 (cause: Henry's-law table rebuilt against Sander
    # 5.0.0, owner-approved 2026-08-26). HMF's Kaw was raised 1.0e-10 -> 5.0e-8, which
    # is above the 1.0e-8 `_NON_OBSERVABLE_KAW_THRESHOLD`, so its headspace_class flips
    # from "low_headspace" to "observable". No direct measurement of HMF's Henry
    # constant exists, but every available route lands above the gate (Sander/HSDB 2015
    # -> 2.2e-8; OPERA 2.6 -> 9.4e-8; vapour-pressure/solubility 1.7e-8 to 1.4e-6), so
    # the shipped 1e-10 was 200-1000x low. The assertion is inverted rather than
    # deleted so the classification stays pinned in both directions.
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")
    target_lookup = {
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "desirable", "data": {}},
    }

    metadata = _headspace_observability_metadata(hmf, target_lookup)

    assert metadata["headspace_observable"] is True
    assert metadata["headspace_class"] == "observable"
    assert metadata["henry_kaw_25c"] is not None
    assert metadata["henry_kaw_25c"] >= 1.0e-8


def test_headspace_observability_metadata_still_gates_a_genuinely_nonvolatile_marker():
    # Added 2026-08-27 alongside the HMF inversion above: the low_headspace class is
    # still reachable and must stay tested. CML (Kaw 1.0e-12, zwitterionic at food pH)
    # is the honest remaining exemplar.
    cml = Species("Nε-(Carboxymethyl)lysine (CML)", "NCCCCC(N)C(=O)O")
    target_lookup = {
        _canon(cml.smiles): {
            "name": "Nε-(Carboxymethyl)lysine (CML)",
            "type": "toxic",
            "data": {},
        },
    }

    metadata = _headspace_observability_metadata(cml, target_lookup)

    assert metadata["headspace_observable"] is False
    assert metadata["headspace_class"] == "low_headspace"
    assert metadata["henry_kaw_25c"] is not None
    assert metadata["henry_kaw_25c"] < 1.0e-8
