"""
The Python API: the same payloads the verbs print, as functions.

    from src import api
    api.compare(spec_a, spec_b)          # per-compound ratios, both declarations
    api.predict(spec)                    # one formulation with intervals
    api.explain("2-methyl-3-furanthiol")
    api.score(document)                  # your measurements against the model
    api.calibrate(document, lab="my lab")   # a per-laboratory calibration and its card

A spec is a mapping with precursors (name -> mM), temp_C, time_min, ph, aw and optionally matrix,
buffer, targets; `api.template()` gives one. A calibration is a Calibration object or the path of
one written by `calibrate`. Nothing here differs from the command line; the verbs call these.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import yaml

from src.comparative_cli import SPEC_TEMPLATE, compare_core, predict_core, validate_spec

CalibrationLike = Union[None, str, Path, "Calibration"]


def _calibration(value: CalibrationLike):
    if value is None:
        return None
    from src.kinetic_core.calibration import Calibration

    if isinstance(value, Calibration):
        return value
    return Calibration.load(value)


def template() -> Dict[str, Any]:
    """The compare template as a document: {"a": spec, "b": spec}."""
    return dict(yaml.safe_load(SPEC_TEMPLATE))


def compare(spec_a: Mapping[str, Any], spec_b: Mapping[str, Any], *, targets: Optional[Sequence[str]] = None,
            calibration: CalibrationLike = None) -> Dict[str, Any]:
    a = validate_spec(spec_a, label="a")
    b = validate_spec(spec_b, label="b")
    return compare_core(a, b, targets=targets, calibration=_calibration(calibration))


def predict(spec: Mapping[str, Any], *, targets: Optional[Sequence[str]] = None,
            calibration: CalibrationLike = None) -> Dict[str, Any]:
    return predict_core(validate_spec(spec, label="spec"), targets=targets, calibration=_calibration(calibration))


def explain(compound: str) -> Dict[str, Any]:
    from src.explain_compound import explain as _explain

    return _explain(compound)


def score(document: Mapping[str, Any], *, calibration: CalibrationLike = None) -> Dict[str, Any]:
    from src.kinetic_core.user_scoring import score_document

    return score_document(document, calibration=_calibration(calibration))


def calibrate(document: Mapping[str, Any], lab: str, *, max_coordinates: int = 4,
              coordinates: Optional[Sequence[str]] = None) -> Tuple["Calibration", Dict[str, Any]]:
    from src.kinetic_core.user_fit import calibrate as _calibrate

    return _calibrate(document, lab, max_coordinates=max_coordinates, coordinates=coordinates)


__all__ = ["calibrate", "compare", "explain", "predict", "score", "template"]
