"""
A per-laboratory CALIBRATION: an overlay on the shipped model, made from the user's own data.

Pre-registered in results/validation/calibration_prereg.md (2026-09-08). Two parts:

* **response factors** -- one log10 offset per compound between what the laboratory measures and
  what the shipped model predicts. A property of the measurement chain, so the user's LEVELS may
  set it (the repository's rule: levels validate, they do not fit chemistry).
* **overrides** -- the few fit-report-space coordinates the user's CONTRASTS could identify, each
  with its prior (the shipped value and sigma) and its posterior. Applied through the engine's own
  override mapping (:func:`engine.core_parameters` ``frozen=``), the same mechanism the envelope
  draws through; nothing is forked.

A calibration lives under results/user/<lab>/ and is read only when a caller passes it. The shipped
parameters, the panel scorecard and every tracked artifact never see it (asserted by
tests/unit/test_calibration.py).
"""
from __future__ import annotations

import copy
import dataclasses
import json
import math
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from src import data_paths

USER_RESULTS_DIR: Path = data_paths.RESULTS_ROOT / "user"


@dataclass(frozen=True)
class Coordinate:
    """One address in a lane's fit-report-space vector (:func:`engine.frozen_parameters`).

    ``block`` is the top-level key; ``key`` the entry inside a dict block (``None`` for a scalar block);
    ``field`` the entry inside a nested dict (the B1 variant stores ``{k_ref_100C, ea_kj_mol}``).
    ``log10`` says the calibration works in log10 of a value the block stores linearly.
    """

    lane: str
    block: str
    key: Optional[str]
    field: Optional[str] = None
    log10: bool = False
    prior_key: str = ""     # the envelope's prior name, e.g. "b8.k_nf_mft.log10_k_ref_145C"
    unit: str = ""

    @property
    def name(self) -> str:
        return self.prior_key or ".".join(x for x in (self.block, self.key, self.field) if x)

    def read(self, vector: Mapping[str, Any]) -> float:
        node: Any = vector[self.block]
        if self.key is not None:
            node = node[self.key]
        if self.field is not None:
            node = node[self.field]
        value = float(node)
        return math.log10(value) if self.log10 else value

    def write(self, vector: Dict[str, Any], value: float) -> None:
        stored = 10.0 ** value if self.log10 else float(value)
        if self.key is None:
            vector[self.block] = stored
            return
        block = vector[self.block]
        if self.field is None:
            block[self.key] = stored
        else:
            block[self.key][self.field] = stored


@dataclass(frozen=True)
class Override:
    coordinate: Coordinate
    prior_value: float
    prior_sigma: float
    value: float
    sigma: float
    band: Optional[Tuple[float, float]]

    @property
    def shift(self) -> float:
        return self.value - self.prior_value


@dataclass(frozen=True)
class ResponseFactor:
    compound: str
    log10: float
    sigma: float
    n_rows: int

    @property
    def factor(self) -> float:
        return 10.0 ** self.log10


@dataclass(frozen=True)
class Calibration:
    lab: str
    base_wave: str
    created: str
    matrix: str
    response_factors: Mapping[str, ResponseFactor]
    overrides: Tuple[Override, ...]
    fit_records: Tuple[str, ...]
    validate_records: Tuple[str, ...]
    notes: Tuple[str, ...] = ()
    provenance: Mapping[str, Any] = field(default_factory=dict)

    # -- overlaying the engine ----------------------------------------------------------------
    def lanes(self) -> Tuple[str, ...]:
        return tuple(sorted({o.coordinate.lane for o in self.overrides}))

    def overlay(self, lane: str) -> Optional[Dict[str, Any]]:
        """The lane's full fit-report-space vector with this calibration's overrides written in,
        or ``None`` when the calibration moves nothing on that lane (the engine then reads the
        frozen reports, byte-identically)."""
        from .engine import LIPID, frozen_parameters

        if lane == LIPID:
            return None
        mine = [o for o in self.overrides if o.coordinate.lane == lane]
        if not mine:
            return None
        vector = copy.deepcopy(frozen_parameters(lane))
        for o in mine:
            o.coordinate.write(vector, o.value)
        return vector

    def operative(self, lane: str) -> Optional[Dict[str, Any]]:
        """The operative parameter dict for ``lane`` under this calibration, or ``None``."""
        from .engine import core_parameters

        vector = self.overlay(lane)
        return None if vector is None else core_parameters(lane, frozen=vector)

    def apply_factors(self, run):
        """A copy of ``run`` with the response factors applied to its concentrations, and the
        calibration named in its metadata. A compound without a factor is left as predicted."""
        if not run.answered or not self.response_factors:
            return run
        scaled = dict(run.concentrations_ug_per_l)
        applied: Dict[str, float] = {}
        for compound, value in run.concentrations_ug_per_l.items():
            rf = self.response_factors.get(compound)
            if rf is not None:
                scaled[compound] = float(value) * rf.factor
                applied[compound] = rf.factor
        metadata = dict(run.run_metadata)
        metadata["calibration"] = {"lab": self.lab, "response_factors_applied": applied,
                                   "matrix_of_records": self.matrix}
        # the factor's own uncertainty widens the interval: 1.645 sigma in decades for a 90 % band
        metadata["calibration_extra_decades"] = {
            compound: 1.645 * self.response_factors[compound].sigma for compound in applied
        }
        spec_matrix = str(getattr(run.spec.process, "matrix", "water") or "water")
        if spec_matrix != self.matrix:
            metadata["calibration"]["warning"] = (
                f"the calibration was made on records in matrix {self.matrix!r}; this spec is in "
                f"{spec_matrix!r}. A response factor does not transfer across matrices."
            )
        return dataclasses.replace(run, concentrations_ug_per_l=scaled, run_metadata=metadata)

    # -- serialisation ------------------------------------------------------------------------
    def as_dict(self) -> Dict[str, Any]:
        return {
            "artifact": "maillard_calibration",
            "lab": self.lab,
            "base_wave": self.base_wave,
            "created": self.created,
            "matrix": self.matrix,
            "response_factors": {
                c: {"log10": rf.log10, "factor": rf.factor, "sigma_log10": rf.sigma, "n_rows": rf.n_rows}
                for c, rf in self.response_factors.items()
            },
            "overrides": [
                {
                    "coordinate": dataclasses.asdict(o.coordinate),
                    "prior_value": o.prior_value, "prior_sigma": o.prior_sigma,
                    "value": o.value, "sigma": o.sigma, "shift": o.shift,
                    "band": list(o.band) if o.band else None,
                }
                for o in self.overrides
            ],
            "fit_records": list(self.fit_records),
            "validate_records": list(self.validate_records),
            "notes": list(self.notes),
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "Calibration":
        if d.get("artifact") != "maillard_calibration":
            raise ValueError("not a calibration file (artifact != maillard_calibration)")
        factors = {
            c: ResponseFactor(c, float(v["log10"]), float(v["sigma_log10"]), int(v["n_rows"]))
            for c, v in d.get("response_factors", {}).items()
        }
        overrides = tuple(
            Override(Coordinate(**o["coordinate"]), float(o["prior_value"]), float(o["prior_sigma"]),
                     float(o["value"]), float(o["sigma"]), tuple(o["band"]) if o.get("band") else None)
            for o in d.get("overrides", [])
        )
        return cls(str(d["lab"]), str(d["base_wave"]), str(d["created"]), str(d.get("matrix", "water")),
                   factors, overrides, tuple(d.get("fit_records", [])), tuple(d.get("validate_records", [])),
                   tuple(d.get("notes", [])), dict(d.get("provenance", {})))

    def save(self, directory: Path | str = USER_RESULTS_DIR) -> Path:
        out_dir = Path(directory) / _slug(self.lab)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"calibration_{self.created}.json"
        path.write_text(json.dumps(self.as_dict(), indent=2, sort_keys=False) + "\n", encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: Path | str) -> "Calibration":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def _slug(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in text.strip()) or "lab"


# -- the calibrated prediction ---------------------------------------------------------------
def predict_calibrated(spec, targets: Sequence[str], calibration: Optional[Calibration]):
    """:func:`engine.predict` under a calibration: the lane's overrides through the engine's own
    override mapping, then the response factors on the answer. ``calibration=None`` is byte-identical
    to the plain engine."""
    from .engine import LIPID, declare_envelope, predict

    if calibration is None:
        return predict(spec, targets)
    declaration = declare_envelope(spec, targets)
    parameters = None
    if declaration.is_answerable:
        lanes = declaration.lanes or ((declaration.lane,) if declaration.lane else ())
        lane = next((l for l in lanes if l != LIPID), None)
        if lane is not None:
            parameters = calibration.operative(lane)
    run = predict(spec, targets, parameters=parameters)
    return calibration.apply_factors(run)


# -- the candidate coordinates ---------------------------------------------------------------
def candidate_coordinates(lane: str) -> List[Tuple[Coordinate, float, float, Optional[Tuple[float, float]]]]:
    """The coordinates a calibration may move on ``lane``: the envelope's fitted priors with a
    finite shipped sigma, addressed in fit-report space. Returns (coordinate, prior value, prior
    sigma, band). A coordinate whose shipped sigma is essentially zero (a flat or bound-limited
    profile) is not offered: the user's contrasts would move it without any check."""
    from . import uncertainty
    from .engine import ACRYLAMIDE, B1_VARIANT, SULFUR, TRUNK

    out = []
    for p in uncertainty.core_priors():
        if p.kind not in ("fitted_rate", "fitted_ea") or p.distribution not in ("normal_log10", "normal"):
            continue
        if p.sigma is None or p.centre is None or p.sigma < 1e-3:
            continue
        parts = p.key.split(".")
        coord: Optional[Coordinate] = None
        if parts[0] == "b1" and len(parts) == 3 and p.lane == TRUNK:
            field_name = "k_ref_100C" if parts[2] == "log10_k_ref_100C" else "ea_kj_mol"
            coord = Coordinate(TRUNK, B1_VARIANT, parts[1], field_name, field_name == "k_ref_100C", p.key, p.unit)
        elif parts[0] == "b3" and p.lane == ACRYLAMIDE:
            if len(parts) == 3 and parts[2] == "log10_k_ref_160C":
                coord = Coordinate(ACRYLAMIDE, "log10_k_ref_at_160C", parts[1], None, False, p.key, p.unit)
            elif len(parts) == 2:
                coord = Coordinate(ACRYLAMIDE, "fitted_Ea_kJ_mol", parts[1], None, False, p.key, p.unit)
        elif parts[0] == "b8" and p.lane == SULFUR:
            if len(parts) == 3 and parts[2] == "log10_k_ref_145C":
                coord = Coordinate(SULFUR, "log10_k_ref_at_145C", parts[1], None, False, p.key, p.unit)
            elif parts[1] == "lumped_formation_Ea_kJ_mol":
                coord = Coordinate(SULFUR, "lumped_formation_Ea_kJ_mol", None, None, False, p.key, p.unit)
            elif parts[1] == "decay_Ea_kJ_mol" and len(parts) == 3:
                coord = Coordinate(SULFUR, "decay_Ea_kJ_mol", parts[2], None, False, p.key, p.unit)
        elif parts[0] == "b7" and p.lane == TRUNK:
            coord = Coordinate(TRUNK, "k_dpo_af", None, None, True, p.key, p.unit)
        if coord is None or coord.lane != lane:
            continue
        out.append((coord, float(p.centre), float(p.sigma), tuple(p.band) if p.band else None))
    return out


def shipped_wave() -> str:
    """The wave the engine reads for its sulfur parameters (the newest report it finds)."""
    from .engine import _B2_FIT_REPORT

    return Path(_B2_FIT_REPORT).stem.replace("_fit_report", "")


def today() -> str:
    return date.today().isoformat()
