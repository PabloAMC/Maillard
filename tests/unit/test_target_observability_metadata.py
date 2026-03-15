from src.pathway_extractor import Species
from src.recommend import _canon, _headspace_observability_metadata


def test_headspace_observability_metadata_marks_hmf_as_low_headspace():
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")
    target_lookup = {
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "desirable", "data": {}},
    }

    metadata = _headspace_observability_metadata(hmf, target_lookup)

    assert metadata["headspace_observable"] is False
    assert metadata["headspace_class"] == "low_headspace"
    assert metadata["henry_kaw_25c"] is not None