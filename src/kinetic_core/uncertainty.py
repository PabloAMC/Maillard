"""
src/kinetic_core/uncertainty.py -- THE CORE MONTE-CARLO ENVELOPE (retirement step B2).

Replaces the legacy engine's ``src/uncertainty_propagation.py`` sampler with
one computed ON THE KINETIC CORE. The legacy sampler perturbed barrier
offsets of a model that is being retired; this one perturbs the core's OWN
fitted coordinates, in the space the fit reports declare them in, and the
declared assumptions the core already carries as bands.

WHAT IS SAMPLED, AND FROM WHERE (the design is settled in
``tasks/data_restructure_plan.md``, "B2 core Monte-Carlo envelope")
--------------------------------------------------------------------------
  * B1 trunk: the two IDENTIFIED fitted pairs (``k_mgo_mel``, ``k_aa_frag``)
    from their reported ``log10_k_ref_stderr`` / ``ea_stderr_kj_mol``; the
    other two carry ``None`` stderr in the report and are FIXED.
  * B3 acrylamide: ``k_acr_dp`` / ``Ea_acr_dp`` from ``ci95_halfwidth / 1.96``
    (the only two the report marks ``identified``); the rest FIXED.
  * B7 furanic: the single fitted constant ``k_dpo_af`` from the fit's
    ``sigma_log10`` (the report carries no parameter stderr; the residual
    sigma is the proxy the plan names). ``k_odg_af`` follows it, as in the fit.
  * DECLARED BANDS, sampled UNIFORM / log-uniform over the band exactly as
    declared -- not as a lognormal 90 %: Q10 [2, 3]; each carrier's lipid
    fraction and peroxide value lo/hi; the furanone partition barrier
    +/-50 kJ/mol; the air/water partition constant +/-0.5 dex; the HS-SPME
    same-sample dispersion (10, 23).
  * B8 sulfur -- MOST OF THE PANEL -- carries NO uncertainty in its fit
    report. Every sulfur prior is listed ``sampled: false`` with reason
    ``no_uncertainty_in_fit_report``. Nothing is invented for it; a Laplace
    covariance at the B8 optimum is a later step.

Draws are independent, seeded from ``numpy.random.SeedSequence(seed).spawn(n)``,
so the artifact is reproducible to the byte for a given ``(n, seed)`` and
does not depend on the number of worker processes.
"""

from __future__ import annotations

import json
import math
import multiprocessing
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from src import artifact_io, data_access, data_paths, provenance
from src.report_format import fmt_number as _fmt
from src.benchmark_metadata import (
    benchmark_signal_origin,
    get_benchmark_metadata,
)

from . import engine
from .fit_targets import core_evidence_role, fit_target_of, in_core_fit
from .engine import (
    ACRYLAMIDE,
    LIPID,
    SULFUR,
    TRUNK,
    CoreDraw,
    engine_metadata,
    frozen_parameters,
    predict,
)
from .panel import (
    RATIO_UNIT_FACTORS,
    SHARED_WITH_HOLDOUT_PANEL,
    bundle_targets,
    core_native_value,
    core_spec,
    limiting_precursor_molar,
    load_bundle,
    measured_value,
    panel_bundles,
    quantification_family,
    QUANTIFICATION_EXTRACTION,
)

CI_LEVEL_PCT = 90
#: Below this width an envelope is degenerate and cannot meaningfully contain
#: or miss a measurement; same threshold as the legacy artifact.
MIN_EVALUABLE_CI_WIDTH_LOG10 = 0.01

NO_UNCERTAINTY = "no_uncertainty_in_fit_report"
#: B8 (2026-09-03): a Laplace covariance at the sulfur optimum exists; identified
#: directions are sampled JOINTLY from it, flat ones stay at the optimum.
LAPLACE_SAMPLED = "laplace_covariance_at_b8_optimum"
LAPLACE_FLAT = "unidentified_direction_in_laplace_covariance"
#: The covariance of the SHIPPED sulfur wave (B9 since 2026-09-03): the file beside the
#: fit report the engine reads, named kinetic_core_<wave>_laplace_covariance.json.
LAPLACE_PATH = data_paths.VALIDATION_DIR / engine._B2_FIT_REPORT.name.replace("_fit_report.json", "_laplace_covariance.json")


# ---------------------------------------------------------------------------
# The priors table
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CorePrior:
    """One row of the priors table. Everything the artifact says about a draw."""

    key: str
    lane: str
    #: "fitted_rate" | "fitted_ea" | "fitted_ph_drift" | "declared_band" | "observable"
    kind: str
    #: "normal_log10" | "normal" | "uniform" | "log_uniform" |
    #: "log_uniform_dispersion" | "fixed"
    distribution: str
    centre: Optional[float]
    sigma: Optional[float]
    band: Optional[Tuple[float, float]]
    unit: str
    source: str
    sampled: bool
    reason: str

    def as_dict(self) -> Dict[str, Any]:
        return {
            "key": self.key,
            "lane": self.lane,
            "kind": self.kind,
            "distribution": self.distribution,
            "centre": self.centre,
            "sigma": self.sigma,
            "band": list(self.band) if self.band is not None else None,
            "unit": self.unit,
            "source": self.source,
            "sampled": self.sampled,
            "reason": self.reason,
        }


def _report(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise SystemExit(
            f"{path} not found. The envelope never fits anything; it reads the "
            f"frozen fit reports. Regenerate them first."
        )
    return json.loads(path.read_text())


def _b1_priors() -> List[CorePrior]:
    from .parameters import FITTED_EA_BOUNDS

    report = _report(engine._B1_FIT_REPORT)
    src = data_paths.rel(engine._B1_FIT_REPORT)
    frozen = report["frozen_parameters"][engine.B1_VARIANT]
    fitted = report["variant_A"]["fitted_parameters"]
    out: List[CorePrior] = []
    for key in frozen:
        k_ref = float(frozen[key]["k_ref_100C"])
        ea = float(frozen[key]["ea_kj_mol"])
        row = fitted.get(key) or {}
        k_sigma = row.get("log10_k_ref_stderr")
        ea_sigma = row.get("ea_stderr_kj_mol")
        out.append(
            CorePrior(
                key=f"b1.{key}.log10_k_ref_100C", lane=TRUNK, kind="fitted_rate",
                distribution="normal_log10" if k_sigma is not None else "fixed",
                centre=math.log10(k_ref),
                sigma=float(k_sigma) if k_sigma is not None else None,
                band=None, unit="log10(k at 100 C)",
                source=f"{src}: variant_A.fitted_parameters.{key}.log10_k_ref_stderr",
                sampled=k_sigma is not None,
                reason=(
                    "identified in the fit report (stderr reported)"
                    if k_sigma is not None else
                    "unidentified in the fit report (log10_k_ref_stderr is null)"
                ),
            )
        )
        out.append(
            CorePrior(
                key=f"b1.{key}.ea_kj_mol", lane=TRUNK, kind="fitted_ea",
                distribution="normal" if ea_sigma is not None else "fixed",
                centre=ea,
                sigma=float(ea_sigma) if ea_sigma is not None else None,
                band=tuple(float(b) for b in FITTED_EA_BOUNDS) if ea_sigma is not None else None,
                unit="kJ/mol",
                source=f"{src}: variant_A.fitted_parameters.{key}.ea_stderr_kj_mol",
                sampled=ea_sigma is not None,
                reason=(
                    "identified in the fit report (stderr reported); draws are "
                    "CLIPPED to the fit's own search bounds FITTED_EA_BOUNDS, "
                    "which the fit itself could not leave"
                    if ea_sigma is not None else
                    "unidentified in the fit report (ea_stderr_kj_mol is null)"
                ),
            )
        )
    return out


def _b3_priors() -> List[CorePrior]:
    report = _report(engine._B3_FIT_REPORT)
    src = data_paths.rel(engine._B3_FIT_REPORT)
    frozen = report["frozen_parameters"]
    intervals = report.get("parameter_intervals") or {}
    out: List[CorePrior] = []
    for key, value in frozen["log10_k_ref_at_160C"].items():
        row = intervals.get(key) or {}
        identified = bool(row.get("identified"))
        half = row.get("ci95_halfwidth")
        out.append(
            CorePrior(
                key=f"b3.{key}.log10_k_ref_160C", lane=ACRYLAMIDE, kind="fitted_rate",
                distribution="normal_log10" if identified else "fixed",
                centre=float(value),
                sigma=float(half) / 1.96 if (identified and half is not None) else None,
                band=None, unit="log10(k at 160 C)",
                source=f"{src}: parameter_intervals.{key}.ci95_halfwidth / 1.96",
                sampled=identified,
                reason=(
                    "identified in the fit report (ci95_halfwidth below the "
                    "identified_threshold)"
                    if identified else
                    f"unidentified in the fit report (ci95_halfwidth {half!r} "
                    f"above identified_threshold {row.get('identified_threshold')!r})"
                ),
            )
        )
    for key, value in frozen["fitted_Ea_kJ_mol"].items():
        row = intervals.get(key) or {}
        identified = bool(row.get("identified"))
        half = row.get("ci95_halfwidth")
        out.append(
            CorePrior(
                key=f"b3.{key}", lane=ACRYLAMIDE, kind="fitted_ea",
                distribution="normal" if identified else "fixed",
                centre=float(value),
                sigma=float(half) / 1.96 if (identified and half is not None) else None,
                band=None, unit="kJ/mol",
                source=f"{src}: parameter_intervals.{key}.ci95_halfwidth / 1.96",
                sampled=identified,
                reason=(
                    "identified in the fit report (ci95_halfwidth below the "
                    "identified_threshold)"
                    if identified else
                    f"unidentified in the fit report (ci95_halfwidth {half!r} "
                    f"above identified_threshold {row.get('identified_threshold')!r})"
                ),
            )
        )
    return out


def _b7_priors() -> List[CorePrior]:
    if not engine._B7_FIT_REPORT.exists():
        return []
    report = _report(engine._B7_FIT_REPORT)
    src = data_paths.rel(engine._B7_FIT_REPORT)
    k = float(report["frozen_parameters"]["k_dpo_af"])
    sigma = report.get("fit", {}).get("sigma_log10")
    return [
        CorePrior(
            key="b7.k_dpo_af.log10_k", lane=TRUNK, kind="fitted_rate",
            distribution="normal_log10" if sigma is not None else "fixed",
            centre=math.log10(k),
            sigma=float(sigma) if sigma is not None else None,
            band=None, unit="log10(k, L/(mmol*min))",
            source=f"{src}: fit.sigma_log10 ({report.get('fit', {}).get('sigma_basis', '')})",
            sampled=sigma is not None,
            reason=(
                "the B7 report carries no parameter stderr; its residual "
                "sigma_log10 is used as the log10-k stderr proxy, per the plan. "
                "k_odg_af is DERIVED from k_dpo_af and moves with it, as in the fit"
                if sigma is not None else NO_UNCERTAINTY
            ),
        )
    ]


def _laplace() -> Optional[Dict[str, Any]]:
    """The shipped wave's Laplace covariance artifact (LAPLACE_PATH follows the fit report the
    engine reads), or None when it has not been generated."""
    payload = data_access.load_json(LAPLACE_PATH, missing_ok=True)
    if not payload or not str(payload.get("artifact", "")).endswith("_laplace_covariance"):
        return None
    return dict(payload)


def _laplace_lookup(lap: Optional[Dict[str, Any]]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """``{(block, key): {sigma, identified, bounds, index}}`` for every free coordinate."""
    if lap is None:
        return {}
    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for i, (coord, sigma, ok, bounds) in enumerate(
        zip(lap["coordinates"], lap["sigma"], lap["identified"], lap["bounds"])
    ):
        out[(coord["block"], coord["key"])] = {
            "sigma": float(sigma), "identified": bool(ok), "bounds": tuple(bounds), "index": i,
        }
    return out


def _b8_priors() -> List[CorePrior]:
    report = _report(engine._B2_FIT_REPORT)
    src = data_paths.rel(engine._B2_FIT_REPORT)
    frozen = report["frozen_parameters"]
    lap = _laplace()
    lookup = _laplace_lookup(lap)
    lap_src = data_paths.rel(LAPLACE_PATH) if lap else None

    def prior(key, block, sub, centre, unit, kind, distribution):
        entry = lookup.get((block, sub))
        if entry is None:
            return CorePrior(
                key=key, lane=SULFUR, kind=kind, distribution="fixed", centre=float(centre),
                sigma=None, band=None, unit=unit, source=f"{src}: frozen_parameters.{block}",
                sampled=False, reason=NO_UNCERTAINTY if lap is None else "frozen in the B8 fit (not a free coordinate)",
            )
        if not entry["identified"]:
            return CorePrior(
                key=key, lane=SULFUR, kind=kind, distribution="fixed", centre=float(centre),
                sigma=entry["sigma"], band=entry["bounds"], unit=unit,
                source=f"{lap_src}: sigma above the identification threshold",
                sampled=False, reason=LAPLACE_FLAT,
            )
        return CorePrior(
            key=key, lane=SULFUR, kind=kind, distribution=distribution, centre=float(centre),
            sigma=entry["sigma"], band=entry["bounds"], unit=unit,
            source=f"{lap_src}: Gauss-Newton sigma at the frozen optimum (sampled JOINTLY from the covariance)",
            sampled=True, reason=LAPLACE_SAMPLED,
        )

    out: List[CorePrior] = []
    for key, value in frozen["log10_k_ref_at_145C"].items():
        out.append(prior(f"b8.{key}.log10_k_ref_145C", "log10_k_ref_at_145C", key, value,
                         "log10(k at 145 C)", "fitted_rate", "normal_log10"))
    out.append(prior("b8.lumped_formation_Ea_kJ_mol", "lumped_formation_Ea_kJ_mol", "",
                     frozen["lumped_formation_Ea_kJ_mol"], "kJ/mol", "fitted_ea", "normal"))
    for family, value in (frozen.get("decay_Ea_kJ_mol") or {}).items():
        out.append(prior(f"b8.decay_Ea_kJ_mol.{family}", "decay_Ea_kJ_mol", family, value,
                         "kJ/mol", "fitted_ea", "normal"))
    for key, value in (frozen.get("ph_drift") or {}).items():
        out.append(prior(f"b8.ph_drift.{key}", "ph_drift", key, value,
                         "pKa" if "pKa" in key else "mol acid per sink event", "fitted_ph_drift", "normal"))
    return out


def sulfur_joint_draw(rng: np.random.Generator) -> Optional[Dict[str, Any]]:
    """
    ONE joint draw of the identified sulfur coordinates from the shipped wave's Laplace covariance.

    Returns ``{"maillard_blocks": {...}, "ph_drift": PhDrift | None, "coords": {...}}``
    or None when no covariance artifact exists. The draw is multivariate normal on the
    identified sub-space (Cholesky of the sub-covariance, jittered), added to the
    frozen optimum and clipped to each coordinate's declared bounds. Unidentified
    coordinates stay at the optimum.
    """
    lap = _laplace()
    if lap is None:
        return None
    idx = [i for i, ok in enumerate(lap["identified"]) if ok]
    if not idx:
        return None
    cov = np.asarray(lap["covariance"], dtype=float)[np.ix_(idx, idx)]
    cov = 0.5 * (cov + cov.T) + np.eye(len(idx)) * 1e-12
    chol = np.linalg.cholesky(cov)
    z = chol @ rng.standard_normal(len(idx))
    optimum = np.asarray(lap["optimum"], dtype=float)
    values = optimum.copy()
    for k, i in enumerate(idx):
        lo, hi = lap["bounds"][i]
        values[i] = min(max(optimum[i] + z[k], lo), hi)
    blocks: Dict[str, Any] = {}
    ph: Dict[str, float] = {}
    coords: Dict[str, float] = {}
    for i, coord in enumerate(lap["coordinates"]):
        if i not in idx:
            continue
        block, key = coord["block"], coord["key"]
        coords[f"b8.{block}.{key}".rstrip(".")] = float(values[i])
        if block == "ph_drift":
            ph[key] = float(values[i])
        elif key:
            blocks.setdefault(block, {})[key] = float(values[i])
        else:
            blocks[block] = float(values[i])
    drift = None
    if ph:
        from .ph_state import PhDrift
        frozen_ph = _report(engine._B2_FIT_REPORT)["frozen_parameters"]["ph_drift"]
        drift = PhDrift(
            acid_yield=ph.get("acid_yield_per_sink_event", frozen_ph["acid_yield_per_sink_event"]),
            arp_amine_pka=ph.get("arp_secondary_ammonium_pKa", frozen_ph["arp_secondary_ammonium_pKa"]),
        )
    return {"maillard_blocks": blocks, "ph_drift": drift, "coords": coords}


def _declared_band_priors() -> List[CorePrior]:
    from .parameters_furanic import FURANONE_PARTITION_EA_BAND_KJ_MOL
    from .parameters_lipid import LIPID_CARRIERS, Q10_ASSUMPTION
    from .parameters_matrix import (
        HS_SPME_DISPERSION_SOURCE,
        HS_SPME_SAME_SAMPLE_DISPERSION,
        K_AW_UNCERTAINTY_DECADES,
        K_AW_UNCERTAINTY_SOURCE,
    )

    out: List[CorePrior] = [
        CorePrior(
            key="lipid.q10", lane=LIPID, kind="declared_band", distribution="uniform",
            centre=float(Q10_ASSUMPTION.default), sigma=None,
            band=(float(Q10_ASSUMPTION.lo), float(Q10_ASSUMPTION.hi)), unit="per 10 C",
            source="parameters_lipid.Q10_ASSUMPTION (Schroen & Berton-Carabin 2022)",
            sampled=True, reason="declared corner band, sampled uniform over it",
        )
    ]
    for key, carrier in LIPID_CARRIERS.items():
        degenerate = carrier.lipid_lo == carrier.lipid_hi
        out.append(
            CorePrior(
                key=f"lipid.{key}.lipid_mass_fraction", lane=LIPID, kind="declared_band",
                distribution="fixed" if degenerate else "log_uniform",
                centre=float(carrier.lipid_mass_fraction), sigma=None,
                band=(float(carrier.lipid_lo), float(carrier.lipid_hi)), unit="kg fat / kg",
                source=f"parameters_lipid.LIPID_CARRIERS[{key!r}] ({carrier.evidence_class})",
                sampled=not degenerate,
                reason=(
                    "degenerate band (fed hydroperoxide: the fraction is the definition)"
                    if degenerate else
                    "declared corner band, sampled log-uniform over it as ONE "
                    "scale shared by every carrier in a draw (CoreDraw."
                    "lipid_fraction_scale), clipped to this carrier's own band"
                ),
            )
        )
        degenerate = carrier.pv_lo == carrier.pv_hi
        out.append(
            CorePrior(
                key=f"lipid.{key}.peroxide_value_meq_per_kg", lane=LIPID, kind="declared_band",
                distribution="fixed" if degenerate else "log_uniform",
                centre=float(carrier.peroxide_value_meq_per_kg), sigma=None,
                band=(float(carrier.pv_lo), float(carrier.pv_hi)), unit="meq O2 / kg fat",
                source=f"parameters_lipid.LIPID_CARRIERS[{key!r}] ({carrier.evidence_class})",
                sampled=not degenerate,
                reason=(
                    "degenerate band (fed hydroperoxide: PV is the definition)"
                    if degenerate else
                    "declared corner band, sampled log-uniform over it as ONE "
                    "scale shared by every carrier in a draw (CoreDraw."
                    "peroxide_scale), clipped to this carrier's own band"
                ),
            )
        )
    band = float(FURANONE_PARTITION_EA_BAND_KJ_MOL)
    out.append(
        CorePrior(
            key="furanic.partition_ea_offset_kj_mol", lane=TRUNK, kind="declared_band",
            distribution="uniform", centre=0.0, sigma=None, band=(-band, band), unit="kJ/mol",
            source="parameters_furanic.FURANONE_EA_ASSUMPTION (declared, no furanone Ea exists)",
            sampled=True,
            reason="declared corner band on the furanone PARTITION barrier, sampled uniform",
        )
    )
    kaw = float(K_AW_UNCERTAINTY_DECADES)
    out.append(
        CorePrior(
            key="observable.air_water_partition_constant", lane="observable",
            kind="observable", distribution="log_uniform", centre=0.0, sigma=None,
            band=(-kaw, kaw), unit="dex (multiplier on the ppb)",
            source=K_AW_UNCERTAINTY_SOURCE, sampled=True,
            reason="declared +/-0.5 dex band on K_aw, sampled uniform in log10",
        )
    )
    lo, hi = HS_SPME_SAME_SAMPLE_DISPERSION
    out.append(
        CorePrior(
            key="observable.hs_spme_same_sample_dispersion", lane="observable",
            kind="observable", distribution="log_uniform_dispersion",
            centre=math.sqrt(float(lo) * float(hi)), sigma=None,
            band=(float(lo), float(hi)), unit="x (same-sample dispersion D)",
            source=HS_SPME_DISPERSION_SOURCE, sampled=True,
            reason=(
                "measured 10-23x same-sample dispersion: D is drawn log-uniform "
                "over [10, 23] and the multiplier is D^u with u uniform on "
                "[-1/2, +1/2], i.e. uniform in log over [1/sqrt(D), sqrt(D)]"
            ),
        )
    )
    return out


def core_priors() -> Tuple[CorePrior, ...]:
    """The full priors table, in the order the sampler consumes it."""
    return tuple(
        _b1_priors() + _b3_priors() + _b7_priors() + _b8_priors() + _declared_band_priors()
    )


#: Built from the fit reports at import, as the plan requires.
CORE_PRIORS: Tuple[CorePrior, ...] = core_priors()


def parameter_sources() -> List[Dict[str, Any]]:
    """The fit reports the draws are made from, with their sha256 (``report`` + ``sha256``)."""
    return [
        {"report": s["path"], "sha256": s["sha256"]}
        for s in provenance.input_sources(engine.fit_report_paths())
        if s["sha256"] is not None
    ]


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PanelDraw:
    """One draw: the core's ``CoreDraw`` plus the two observable multipliers."""

    index: int
    core: CoreDraw
    k_aw_multiplier: float
    hs_spme_multiplier: float
    #: the raw sampled coordinates, keyed by prior key, for the record
    coordinates: Mapping[str, float]


def _lipid_scale_band(priors: Sequence[CorePrior], suffix: str) -> Tuple[float, float]:
    """The union, over carriers with a real band, of (lo/centre, hi/centre)."""
    los, his = [], []
    for p in priors:
        if p.kind == "declared_band" and p.key.endswith(suffix) and p.sampled and p.band:
            los.append(p.band[0] / p.centre)
            his.append(p.band[1] / p.centre)
    if not los:
        return 1.0, 1.0
    return min(los), max(his)


def draw_from_rng(
    rng: np.random.Generator, index: int, priors: Sequence[CorePrior] = CORE_PRIORS
) -> PanelDraw:
    """One ``PanelDraw`` from one generator, consuming the priors in table order."""
    b1 = {k: dict(v) for k, v in frozen_parameters(TRUNK)[engine.B1_VARIANT].items()}
    b3_k: Dict[str, float] = {}
    b3_ea: Dict[str, float] = {}
    k_dpo_af: Optional[float] = None
    coords: Dict[str, float] = {}
    q10 = None
    furanone = None
    lipid_u = None
    pv_u = None
    k_aw = 1.0
    hs = 1.0

    lipid_lo, lipid_hi = _lipid_scale_band(priors, ".lipid_mass_fraction")
    pv_lo, pv_hi = _lipid_scale_band(priors, ".peroxide_value_meq_per_kg")

    def _scale(u, lo, hi):
        if u is None:
            return None
        return float(10.0 ** (math.log10(lo) + u * (math.log10(hi) - math.log10(lo))))

    for p in priors:
        if not p.sampled or p.key.startswith("b8."):
            continue  # b8 coordinates are drawn JOINTLY below
        if p.distribution == "normal_log10":
            value = float(rng.normal(p.centre, p.sigma))
        elif p.distribution == "normal":
            value = float(rng.normal(p.centre, p.sigma))
            if p.band is not None:
                value = min(max(value, p.band[0]), p.band[1])
        elif p.distribution == "uniform":
            value = float(rng.uniform(p.band[0], p.band[1]))
        elif p.distribution == "log_uniform":
            if p.kind == "observable":
                value = float(rng.uniform(p.band[0], p.band[1]))  # already in dex
            elif p.key.endswith(".lipid_mass_fraction"):
                # ONE quantile shared by every carrier (see the prior's reason):
                # drawn on the first carrier, reused on the rest.
                if lipid_u is None:
                    lipid_u = float(rng.uniform(0.0, 1.0))
                value = float(p.centre) * _scale(lipid_u, lipid_lo, lipid_hi)
                value = min(max(value, p.band[0]), p.band[1])
            elif p.key.endswith(".peroxide_value_meq_per_kg"):
                if pv_u is None:
                    pv_u = float(rng.uniform(0.0, 1.0))
                value = float(p.centre) * _scale(pv_u, pv_lo, pv_hi)
                value = min(max(value, p.band[0]), p.band[1])
            else:
                raise ValueError(f"{p.key}: no log_uniform routing")
        elif p.distribution == "log_uniform_dispersion":
            log_d = float(rng.uniform(math.log10(p.band[0]), math.log10(p.band[1])))
            u = float(rng.uniform(-0.5, 0.5))
            value = 10.0 ** (log_d * u)
            coords[p.key + ".D"] = 10.0 ** log_d
        else:
            raise ValueError(f"{p.key}: distribution {p.distribution!r} is not samplable")
        coords[p.key] = value

        # --- route the coordinate to where the core reads it ----------------
        if p.key.startswith("b1."):
            _, key, field = p.key.split(".")
            if field == "log10_k_ref_100C":
                b1[key]["k_ref_100C"] = 10.0 ** value
            else:
                b1[key]["ea_kj_mol"] = value
        elif p.key.startswith("b3.") and p.key.endswith(".log10_k_ref_160C"):
            b3_k[p.key.split(".")[1]] = value
        elif p.key.startswith("b3."):
            b3_ea[p.key.split(".")[1]] = value
        elif p.key == "b7.k_dpo_af.log10_k":
            k_dpo_af = 10.0 ** value
        elif p.key == "lipid.q10":
            q10 = value
        elif p.key.endswith((".lipid_mass_fraction", ".peroxide_value_meq_per_kg")):
            pass  # routed through lipid_u / pv_u above
        elif p.key == "furanic.partition_ea_offset_kj_mol":
            furanone = value
        elif p.key == "observable.air_water_partition_constant":
            k_aw = 10.0 ** value
        elif p.key == "observable.hs_spme_same_sample_dispersion":
            hs = value

    maillard: Dict[str, Any] = {engine.B1_VARIANT: b1}
    sulfur = sulfur_joint_draw(rng) if any(p.sampled and p.key.startswith("b8.") for p in priors) else None
    ph_drift = None
    if sulfur is not None:
        frozen_sulfur = frozen_parameters(SULFUR)
        for block, value in sulfur["maillard_blocks"].items():
            if isinstance(value, dict):
                merged = dict(frozen_sulfur[block]); merged.update(value); maillard[block] = merged
            else:
                maillard[block] = value
        for block in ("log10_k_ref_at_145C", "lumped_formation_Ea_kJ_mol", "decay_Ea_kJ_mol"):
            maillard.setdefault(block, frozen_sulfur[block])
        ph_drift = sulfur["ph_drift"]
        coords.update(sulfur["coords"])
    if b3_k or b3_ea:
        frozen_b3 = frozen_parameters(ACRYLAMIDE)
        k_block = dict(frozen_b3["log10_k_ref_at_160C"]); k_block.update(b3_k)
        ea_block = dict(frozen_b3["fitted_Ea_kJ_mol"]); ea_block.update(b3_ea)
        maillard["log10_k_ref_at_160C"] = k_block
        maillard["fitted_Ea_kJ_mol"] = ea_block
    if k_dpo_af is not None:
        maillard["k_dpo_af"] = k_dpo_af

    core = CoreDraw(
        maillard=maillard,
        q10=q10,
        lipid_fraction_scale=_scale(lipid_u, lipid_lo, lipid_hi),
        peroxide_scale=_scale(pv_u, pv_lo, pv_hi),
        furanone_partition_ea_kj_mol=furanone,
        ph_drift=ph_drift,  # the wave's Laplace covariance when present, else the frozen calibration
    )
    coords["lipid.fraction_scale"] = core.lipid_fraction_scale if core.lipid_fraction_scale is not None else float("nan")
    coords["lipid.peroxide_scale"] = core.peroxide_scale if core.peroxide_scale is not None else float("nan")
    return PanelDraw(
        index=index, core=core, k_aw_multiplier=k_aw, hs_spme_multiplier=hs,
        coordinates=coords,
    )


def sample_draws(
    n: int, seed: int = 0, priors: Sequence[CorePrior] = CORE_PRIORS
) -> Tuple[PanelDraw, ...]:
    """``n`` independent draws from ``SeedSequence(seed).spawn(n)``."""
    children = np.random.SeedSequence(int(seed)).spawn(int(n))
    return tuple(
        draw_from_rng(np.random.default_rng(child), i, priors)
        for i, child in enumerate(children)
    )


# ---------------------------------------------------------------------------
# Propagation over the panel
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Job:
    """One integration per draw: a spec and the compounds it answers together."""

    bench_index: int
    lane: Optional[str]
    spec: Any
    compounds: Tuple[str, ...]
    units: Tuple[str, ...]
    #: per compound: multiply by the draw's observable (K_aw, HS-SPME) bands?
    observable: Tuple[bool, ...]
    limiting_molar: Optional[float]


def _run_draw(args: Tuple[PanelDraw, Tuple[_Job, ...]]) -> List[List[Optional[float]]]:
    draw, jobs = args
    out: List[List[Optional[float]]] = []
    for job in jobs:
        run = predict(job.spec, list(job.compounds), draw=draw.core, size_declared_bands=False)
        values: List[Optional[float]] = []
        for compound, unit, observable in zip(job.compounds, job.units, job.observable):
            value = (
                core_native_value(run, compound, unit, job.limiting_molar)
                if run.answered else None
            )
            if value is not None and observable:
                value = float(value) * draw.k_aw_multiplier * draw.hs_spme_multiplier
            values.append(None if value is None else float(value))
        out.append(values)
    return out


def _percentile(values: Sequence[float], pct: float) -> float:
    return float(np.percentile(np.asarray(values, dtype=float), pct, method="linear"))


def propagate_panel(
    bundles: Optional[Sequence[Tuple[Path, str]]] = None,
    *,
    n_samples: int = 200,
    seed: int = 0,
    workers: int = 1,
    priors: Sequence[CorePrior] = CORE_PRIORS,
) -> Dict[str, Any]:
    """
    Re-integrate every answerable (benchmark, compound) of the panel per draw.

    ``bundles`` is ``[(path, panel_tag), ...]``; ``None`` means the union panel
    from :func:`panel_bundles`. ``workers > 1`` fans the draws out over
    processes; the result is identical for any worker count because each
    draw is seeded independently and collected by index.
    """
    if bundles is None:
        bundles, skipped = panel_bundles()
    else:
        bundles, skipped = list(bundles), []

    # ---- the deterministic pass: what is answerable, and its point ----------
    benches: List[Dict[str, Any]] = []
    jobs: List[_Job] = []
    refused: List[Dict[str, Any]] = []
    for bench_index, (path, panel_tag) in enumerate(bundles):
        bench = load_bundle(path)
        benchmark_id = str(bench.get("benchmark_id") or Path(path).stem)
        spec = core_spec(bench, use_buffer=True)
        _, limiting_molar = limiting_precursor_molar(bench)
        family, family_source = quantification_family(bench)
        entry: Dict[str, Any] = {
            "benchmark_id": benchmark_id,
            "bench_file": data_paths.rel(path),
            "panel": panel_tag,
            "execution_path": get_benchmark_metadata(bench).execution_path,
            "protein_type": str(bench.get("protein_type", "") or ""),
            "signal_origin": benchmark_signal_origin(Path(path)),
            # B3: the CORE's evidence role (src/kinetic_core/fit_targets.py).
            "evidence_role": core_evidence_role(benchmark_id, Path(path)),
            "fit_target_of": fit_target_of(benchmark_id),
            "buffer_applied": True,
            "quantification_family": family,
            "quantification_source": family_source,
            "rows": {},          # compound -> row dict (filled below)
            "refused_compounds": [],
        }
        entry["fitted_row"] = entry["evidence_role"] == "fit_recovery"
        by_lane: Dict[Optional[str], List[Tuple[str, str, bool]]] = {}
        for compound, target_spec in bundle_targets(bench).items():
            unit = str(target_spec.get("target_unit", "ppb"))
            measured = measured_value(bench, compound, target_spec)
            point_run = predict(spec, [compound])
            declaration = point_run.declaration
            reason = None
            if not point_run.answered:
                reason = " | ".join(declaration.reasons) or "out of envelope"
            elif unit != "ppb" and unit not in RATIO_UNIT_FACTORS:
                reason = f"target unit {unit!r} has no core conversion"
            elif measured is None:
                reason = "the bundle carries no measured value for this compound"
            if reason is not None:
                item = {
                    "benchmark_id": benchmark_id, "panel": panel_tag,
                    "compound": compound, "reason": reason,
                    "envelope_state": declaration.state, "lane": declaration.lane,
                }
                entry["refused_compounds"].append(item)
                refused.append(item)
                continue
            point = core_native_value(point_run, compound, unit, limiting_molar)
            # The K_aw and HS-SPME bands are facts about a HEADSPACE number. A
            # SIDA / HPLC / LC-MS value of the liquid never passed through an
            # air/water partition or a fibre, so they do not apply to it.
            observable = unit == "ppb" and family != QUANTIFICATION_EXTRACTION
            entry["rows"][compound] = {
                "compound": compound,
                "target_unit": unit,
                "measured": float(measured),
                "measured_ppb": float(measured) if unit == "ppb" else None,
                "predicted_point": None if point is None else float(point),
                "lane": declaration.lane,
                "envelope_state": declaration.state,
                "shared_with": SHARED_WITH_HOLDOUT_PANEL.get((benchmark_id, compound)),
                "observable_multipliers_applied": observable,
                "in_core_fit": in_core_fit(benchmark_id, compound),
                "samples": [],
            }
            by_lane.setdefault(declaration.lane, []).append((compound, unit, observable))
        for lane, pairs in by_lane.items():
            jobs.append(
                _Job(
                    bench_index=bench_index, lane=lane, spec=spec,
                    compounds=tuple(c for c, _, _ in pairs),
                    units=tuple(u for _, u, _ in pairs),
                    observable=tuple(o for _, _, o in pairs),
                    limiting_molar=limiting_molar,
                )
            )
        benches.append(entry)

    # ---- the draws -------------------------------------------------------
    draws = sample_draws(n_samples, seed, priors)
    tasks = [(draw, tuple(jobs)) for draw in draws]
    if workers and workers > 1 and tasks:
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(processes=int(workers)) as pool:
            results = pool.map(_run_draw, tasks, chunksize=1)
    else:
        results = [_run_draw(task) for task in tasks]
    for per_draw in results:
        for job, values in zip(jobs, per_draw):
            rows = benches[job.bench_index]["rows"]
            for compound, value in zip(job.compounds, values):
                if value is not None and math.isfinite(value):
                    rows[compound]["samples"].append(value)

    # ---- percentiles, coverage, the split -------------------------------
    minimum = max(5, int(n_samples) // 4)
    coverage_hits = 0
    coverage_total = 0
    split: Dict[str, Dict[str, Any]] = {
        name: {"hits": 0, "total": 0, "not_evaluable": 0, "ci_widths_log10": []}
        for name in ("external_literature", "internal_synthetic", "fitted_row", "in_core_fit")
    }
    panel_split: Dict[str, Dict[str, Any]] = {}
    out_benches: List[Dict[str, Any]] = []
    for entry in benches:
        compounds: List[Dict[str, Any]] = []
        for compound, row in entry.pop("rows").items():
            samples = row.pop("samples")
            if len(samples) < minimum:
                item = {
                    "benchmark_id": entry["benchmark_id"], "panel": entry["panel"],
                    "compound": compound,
                    "reason": (
                        f"only {len(samples)} finite samples of {n_samples} "
                        f"(minimum {minimum}); no envelope formed"
                    ),
                    "envelope_state": row["envelope_state"], "lane": row["lane"],
                }
                entry["refused_compounds"].append(item)
                refused.append(item)
                continue
            p5, p50, p95 = (
                _percentile(samples, 5.0), _percentile(samples, 50.0), _percentile(samples, 95.0)
            )
            measured = row["measured"]
            inside = bool(p5 <= measured <= p95)
            width = math.log10(p95 / p5) if (p5 > 0.0 and p95 > 0.0) else None
            row.update(
                {
                    "predicted_p5": p5, "predicted_p50": p50, "predicted_p95": p95,
                    "predicted_mean": statistics.fmean(samples),
                    "predicted_std": statistics.pstdev(samples) if len(samples) > 1 else 0.0,
                    "inside_ci": inside,
                    "ci_width_log10": width,
                    "n_finite_samples": len(samples),
                }
            )
            compounds.append(row)
            coverage_total += 1
            if inside:
                coverage_hits += 1
            bucket = (
                split["fitted_row"] if entry["evidence_role"] == "fit_recovery"
                else split[entry["signal_origin"]]
            )
            pbucket = panel_split.setdefault(
                entry["panel"], {"hits": 0, "total": 0, "not_evaluable": 0, "ci_widths_log10": []}
            )
            buckets = [bucket, pbucket]
            if bucket is split["external_literature"] and row["in_core_fit"]:
                buckets.append(split["in_core_fit"])
            for b in buckets:
                if width is None or width <= MIN_EVALUABLE_CI_WIDTH_LOG10:
                    b["not_evaluable"] += 1
                    continue
                b["total"] += 1
                b["ci_widths_log10"].append(width)
                if inside:
                    b["hits"] += 1
        entry["compounds"] = compounds
        entry["matched_compounds"] = len(compounds)
        n_targets = len(compounds) + len(entry["refused_compounds"])
        entry["coverage_rate"] = (len(compounds) / n_targets) if n_targets else 0.0
        out_benches.append(entry)

    for table in (split, panel_split):
        for bucket in table.values():
            widths = sorted(bucket.pop("ci_widths_log10"))
            bucket["median_ci_width_log10"] = _percentile(widths, 50.0) if widths else None
            bucket["ci_coverage_rate"] = bucket["hits"] / bucket["total"] if bucket["total"] else None

    lit = split["external_literature"]
    fitted = split["fitted_row"]
    sampled = [p for p in priors if p.sampled]
    payload: Dict[str, Any] = {
        "artifact": "core_prediction_uncertainty",
        "engine": engine_metadata(),
        "provenance": provenance.provenance_block(
            "core_prediction_uncertainty", generated_by="src/kinetic_core/uncertainty.py",
            inputs=(*engine.fit_report_paths(), LAPLACE_PATH, *(path for path, _tag in bundles)),
        ),
        "summary": {
            "n_samples": int(n_samples),
            "seed": int(seed),
            "benchmark_count": sum(1 for b in out_benches if b["compounds"]),
            "panel_benchmark_count": len(out_benches),
            "matched_compound_count": coverage_total,
            "refused_compound_count": len(refused),
            "ci_coverage_hits": coverage_hits,
            "ci_coverage_rate": (coverage_hits / coverage_total) if coverage_total else None,
            "ci_level_pct": CI_LEVEL_PCT,
            "signal_origin_split": split,
            "panel_split": panel_split,
            "honest_literature_coverage": {
                "hits": lit["hits"],
                "total": lit["total"],
                "rate": lit["ci_coverage_rate"],
                "not_evaluable": lit["not_evaluable"],
                "median_ci_width_log10": lit["median_ci_width_log10"],
                "excluded_fitted_rows": fitted["total"] + fitted["not_evaluable"],
                "excluded_fitted_rows_that_would_have_been_hits": fitted["hits"],
                "definition": (
                    "External-literature rows only, with fitted rows (constants "
                    "back-solved from the benchmark) and internal synthetic "
                    "reproducibility rows removed from BOTH numerator and "
                    "denominator; computed identically across the trust-loop, "
                    "maillard_path hold-out and external-matrix panels. Read it "
                    "together with median_ci_width_log10: a wide interval covering "
                    "a measurement is a weak claim. The SULFUR lane's identified "
                    "coordinates are sampled jointly from the shipped wave's Laplace covariance "
                    "(see priors[]); unidentified ones stay at the optimum."
                ),
            },
            "observable_multiplier_policy": {
                "rule": (
                    "the K_aw (+/-0.5 dex) and HS-SPME same-sample dispersion bands "
                    "multiply a ppb row only when the bundle's quantification_class "
                    "is a headspace method, or is undeclared (then the core's own "
                    "absolute_concentration convention is kept and the row says so); "
                    "extraction-based quantification (SIDA, HPLC, LC-MS/MS) never "
                    "passed through an air/water partition or a fibre."
                ),
                "rows_by_family": {
                    fam: sum(
                        1 for b in out_benches if b["quantification_family"] == fam
                        for _ in b["compounds"]
                    )
                    for fam in ("headspace", "extraction", "undeclared")
                },
                "undeclared_bundles": sorted(
                    b["benchmark_id"] for b in out_benches
                    if b["quantification_family"] == "undeclared" and b["compounds"]
                ),
            },
            "out_of_sample_literature_coverage": {
                "hits": lit["hits"] - split["in_core_fit"]["hits"],
                "total": lit["total"] - split["in_core_fit"]["total"],
                "not_evaluable": lit["not_evaluable"] - split["in_core_fit"]["not_evaluable"],
                "rows_the_core_fit_read": {
                    k: split["in_core_fit"][k] for k in ("hits", "total", "not_evaluable")
                },
                "definition": (
                    "honest_literature_coverage MINUS every row a core fit read "
                    "(src/kinetic_core/fit_targets.py, from the shipped wave's fit-target record; "
                    "since B9 only the step-level bundles). The stronger claim."
                ),
            },
            "sampled_prior_count": len(sampled),
            "fixed_prior_count": len(priors) - len(sampled),
            "unsampled_lanes": sorted(
                {p.lane for p in priors if not p.sampled and p.reason == NO_UNCERTAINTY}
            ),
            "sulfur_laplace": (
                None if _laplace() is None else {
                    "path": data_paths.rel(LAPLACE_PATH),
                    "identified": int(_laplace()["identified_count"]),
                    "free": int(_laplace()["n_free"]),
                    "reduced_chi_square": float(_laplace()["reduced_chi_square"]),
                }
            ),
            "parameter_sources": parameter_sources(),
            "shared_with_holdout_panel": {
                "n": sum(1 for b in out_benches for c in b["compounds"] if c.get("shared_with")),
                "declaration": (
                    "The four Hofmann 1998 ribose/cysteine pH-3 and pH-7 rows are the "
                    "SAME MEASUREMENTS as four rows of the B2.x hold-out panel: NOT "
                    "independent evidence. Each carries `shared_with`."
                ),
            },
        },
        "priors": [p.as_dict() for p in priors],
        "benchmarks": out_benches,
        "refused_compounds": refused,
        "skipped_bundles": skipped,
    }
    return payload


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render_markdown(payload: Mapping[str, Any]) -> str:
    s = payload["summary"]
    lit = s["honest_literature_coverage"]
    out = [
        "# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)",
        "",
        f"n_samples = {s['n_samples']}, seed = {s['seed']}, CI level = {s['ci_level_pct']} %.",
        "",
        f"* benchmarks with an envelope: **{s['benchmark_count']}** of "
        f"{s['panel_benchmark_count']} on the panel; matched rows **{s['matched_compound_count']}**; "
        f"refused rows {s['refused_compound_count']}",
        f"* mixed-population coverage: {s['ci_coverage_hits']}/{s['matched_compound_count']} "
        f"({_fmt(s['ci_coverage_rate'])})",
        f"* **honest literature coverage: {lit['hits']}/{lit['total']} ({_fmt(lit['rate'])})**, "
        f"median CI width {_fmt(lit['median_ci_width_log10'])} log10; "
        f"{lit['not_evaluable']} not evaluable; {lit['excluded_fitted_rows']} fitted rows excluded",
        f"* out-of-sample literature coverage: {s['out_of_sample_literature_coverage']['hits']}/"
        f"{s['out_of_sample_literature_coverage']['total']} "
        f"({s['out_of_sample_literature_coverage']['not_evaluable']} not evaluable); rows the core fit read: "
        f"{s['out_of_sample_literature_coverage']['rows_the_core_fit_read']}",
        f"* sampled priors {s['sampled_prior_count']}, fixed {s['fixed_prior_count']}; "
        f"lanes with NO sampled fit uncertainty: {', '.join(s['unsampled_lanes']) or 'none'}",
        "* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: "
        + ", ".join(f"{k} {v}" for k, v in s["observable_multiplier_policy"]["rows_by_family"].items())
        + (
            "; undeclared bundles: " + ", ".join(s["observable_multiplier_policy"]["undeclared_bundles"])
            if s["observable_multiplier_policy"]["undeclared_bundles"] else ""
        ),
        "",
        "## Per panel",
        "",
        "| panel | hits | total | rate | median width (log10) | not evaluable |",
        "|---|---|---|---|---|---|",
    ]
    for name, b in sorted(s["panel_split"].items()):
        out.append(
            f"| {name} | {b['hits']} | {b['total']} | {_fmt(b['ci_coverage_rate'])} | "
            f"{_fmt(b['median_ci_width_log10'])} | {b['not_evaluable']} |"
        )
    out += [
        "",
        "## Rows",
        "",
        "| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for b in payload["benchmarks"]:
        for c in b["compounds"]:
            out.append(
                f"| {b['benchmark_id']} | {b['panel']} | {c['compound']} | {c['target_unit']} | "
                f"{_fmt(c['measured'])} | {_fmt(c['predicted_point'])} | {_fmt(c['predicted_p5'])} | "
                f"{_fmt(c['predicted_p50'])} | {_fmt(c['predicted_p95'])} | "
                f"{'yes' if c['inside_ci'] else 'no'} | {_fmt(c['ci_width_log10'])} | "
                f"{('yes (' + b['quantification_family'] + ')') if c['observable_multipliers_applied'] else 'no (' + b['quantification_family'] + ')'} | "
                f"{c['lane']} | {b['evidence_role']}"
                f"{' [in core fit]' if c.get('in_core_fit') else ''}"
                f"{' (shared: ' + c['shared_with'] + ')' if c.get('shared_with') else ''} |"
            )
    if payload["refused_compounds"]:
        out += ["", "## Refused rows", "", "| benchmark | panel | compound | reason |", "|---|---|---|---|"]
        for r in payload["refused_compounds"]:
            reason = str(r["reason"]).replace("|", "/")
            out.append(f"| {r['benchmark_id']} | {r['panel']} | {r['compound']} | {reason[:220]} |")
    if payload.get("skipped_bundles"):
        out += ["", "## Bundles kept off the scored panel", ""]
        for r in payload["skipped_bundles"]:
            out.append(f"- {r['benchmark_id']} ({r['panel']}): {r['reason']}")
    out += [
        "",
        "## Priors",
        "",
        "| key | lane | kind | distribution | centre | sigma | band | sampled | reason |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for p in payload["priors"]:
        band = "-" if p["band"] is None else f"[{_fmt(p['band'][0])}, {_fmt(p['band'][1])}]"
        out.append(
            f"| {p['key']} | {p['lane']} | {p['kind']} | {p['distribution']} | {_fmt(p['centre'])} | "
            f"{_fmt(p['sigma'])} | {band} | {'yes' if p['sampled'] else 'no'} | "
            f"{str(p['reason']).replace('|', '/')[:120]} |"
        )
    out.append("")
    return "\n".join(out)


def write_artifact(
    payload: Mapping[str, Any], json_path: Path, md_path: Optional[Path] = None
) -> Tuple[Path, Path]:
    """JSON (sorted keys) + markdown twin through the one artifact writer (src/artifact_io.py)."""
    written, md_written = artifact_io.write_artifact(
        payload, json_path, render=render_markdown, md_path=md_path, sort_keys=True
    )
    return written, Path(md_written)


__all__ = [
    "CI_LEVEL_PCT",
    "CORE_PRIORS",
    "CorePrior",
    "NO_UNCERTAINTY",
    "PanelDraw",
    "core_priors",
    "draw_from_rng",
    "parameter_sources",
    "propagate_panel",
    "render_markdown",
    "sample_draws",
    "write_artifact",
]
