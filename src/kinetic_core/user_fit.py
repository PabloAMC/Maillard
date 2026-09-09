"""
The two-stage fit behind `maillard calibrate` (results/validation/calibration_prereg.md).

Stage one: the laboratory's RESPONSE FACTORS from the fit records' levels.
Stage two: the KINETIC OVERRIDES from CONTRASTS between fit records (log-ratios, from which the
response factor cancels), restricted to the coordinates the contrasts can identify and pulled toward
the shipped values by the shipped uncertainty (a maximum-a-posteriori least squares).
Then stage one again at the fitted kinetics, and the validate records scored before and after.

The engine is reached only through ``predict_fn`` so the machinery is testable on a stand-in model;
the default is the real engine under :func:`calibration.predict_calibrated`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from src.kinetic_core import calibration as C
from src.kinetic_core.user_scoring import systems_of

#: Below this the log10 sigma of a measured level is floored: a laboratory does not know its own
#: number to better than about 12 %.
SIGMA_FLOOR_LOG10 = 0.05
#: A measurement without a stated uncertainty gets this log10 sigma (about 40 %).
SIGMA_DEFAULT_LOG10 = 0.15
#: A candidate coordinate is kept when its pivoted-QR diagonal is above this fraction of the largest.
IDENTIFIABILITY_TOL = 1e-2
#: A coordinate whose shipped sigma moves no contrast residual by more than this is not identifiable.
IDENTIFIABILITY_FLOOR = 1e-3
FOLD_PASS = 3.0


@dataclass(frozen=True)
class Record:
    name: str
    spec: Mapping[str, Any]       # the validated system (precursors, temp_C, time_min, ph, aw, matrix, measured)
    role: Optional[str]           # "fit" | "validate" | None (assigned later)

    @property
    def measured(self) -> Mapping[str, Mapping[str, Any]]:
        return self.spec["measured"]

    def sort_key(self) -> Tuple:
        recipe = tuple(sorted((k, float(v)) for k, v in self.spec["precursors"].items()))
        return (float(self.spec["temp_C"]), float(self.spec["time_min"]), float(self.spec["ph"]), recipe, self.name)


@dataclass(frozen=True)
class Prediction:
    values: Mapping[str, float]   # predicted ug/L per compound answered
    lane: Optional[str]


PredictFn = Callable[[Mapping[str, Any], Sequence[str], Optional[C.Calibration]], Prediction]


def engine_predict(spec: Mapping[str, Any], compounds: Sequence[str], calibration: Optional[C.Calibration]) -> Prediction:
    from src.comparative_cli import spec_to_core

    run = C.predict_calibrated(spec_to_core(spec), list(compounds), calibration)
    if not run.answered:
        return Prediction({}, None)
    return Prediction({c: float(v) for c, v in run.concentrations_ug_per_l.items()}, run.declaration.lane)


@dataclass(frozen=True)
class Contrast:
    compound: str
    a: str        # record names
    b: str
    log10_ratio: float
    sigma: float


def records_from_document(document: Mapping[str, Any]) -> List[Record]:
    return [Record(str(s["name"]), s, s.get("role")) for s in systems_of(document)]


def assign_roles(records: Sequence[Record]) -> Tuple[List[Record], List[Record], str]:
    """Tagged roles are kept. Untagged records: with four or more, every second one in sorted order
    of (temperature, time, pH, recipe) is held out; fewer than four, nothing is."""
    tagged = [r for r in records if r.role in ("fit", "validate")]
    untagged = [r for r in records if r.role not in ("fit", "validate")]
    fit = [r for r in tagged if r.role == "fit"]
    validate = [r for r in tagged if r.role == "validate"]
    if untagged:
        ordered = sorted(untagged, key=Record.sort_key)
        if len(ordered) >= 4 and not validate:
            for i, r in enumerate(ordered):
                (validate if i % 2 == 1 else fit).append(Record(r.name, r.spec, "validate" if i % 2 == 1 else "fit"))
            note = "hold-out chosen before the fit: every second untagged record in sorted order of temperature, time, pH and recipe"
        else:
            fit.extend(Record(r.name, r.spec, "fit") for r in ordered)
            note = ("all untagged records fit (fewer than four): the calibration is UNVALIDATED"
                    if not validate else "untagged records fit; the tagged validate records are the hold-out")
    else:
        note = "roles as tagged in the document"
    if not fit:
        raise ValueError("no record is available to fit: tag at least one `role: fit`")
    return fit, validate, note


def sigma_of(unc_pct: Optional[float]) -> float:
    if unc_pct is None:
        return SIGMA_DEFAULT_LOG10
    return max(SIGMA_FLOOR_LOG10, math.log10(1.0 + float(unc_pct) / 100.0))


def contrasts_of(fit: Sequence[Record]) -> List[Contrast]:
    """Consecutive pairs, per compound, of the fit records that measure it: n records give n - 1
    independent log-ratios, and the laboratory's factor cancels from each."""
    out: List[Contrast] = []
    compounds = sorted({c for r in fit for c in r.measured})
    for compound in compounds:
        have = sorted([r for r in fit if compound in r.measured], key=Record.sort_key)
        for a, b in zip(have, have[1:]):
            ma, mb = a.measured[compound], b.measured[compound]
            if ma["value"] <= 0 or mb["value"] <= 0:
                continue
            sigma = math.hypot(sigma_of(ma.get("uncertainty_pct")), sigma_of(mb.get("uncertainty_pct")))
            out.append(Contrast(compound, a.name, b.name, math.log10(mb["value"] / ma["value"]), sigma))
    return out


class Fitter:
    """Holds the records and the predictor, and logs every record a residual reads (the T3 guard)."""

    def __init__(self, fit: Sequence[Record], validate: Sequence[Record], predict_fn: PredictFn = engine_predict):
        self.fit = list(fit)
        self.validate = list(validate)
        self.predict_fn = predict_fn
        self.reads_during_fit: Set[str] = set()
        self._fitting = False
        self.evaluations = 0

    # -- predictions -------------------------------------------------------------------------
    def predict(self, records: Sequence[Record], calibration: Optional[C.Calibration]) -> Dict[str, Prediction]:
        out = {}
        for r in records:
            if self._fitting:
                self.reads_during_fit.add(r.name)
            self.evaluations += 1
            out[r.name] = self.predict_fn(r.spec, sorted(r.measured), calibration)
        return out

    def lanes(self) -> List[str]:
        preds = self.predict(self.fit, None)
        return sorted({p.lane for p in preds.values() if p.lane})

    # -- stage one ---------------------------------------------------------------------------
    def response_factors(self, calibration: Optional[C.Calibration]) -> Dict[str, C.ResponseFactor]:
        preds = self.predict(self.fit, calibration)
        rows: Dict[str, List[Tuple[float, float]]] = {}
        for r in self.fit:
            p = preds[r.name]
            for compound, m in r.measured.items():
                pv = p.values.get(compound)
                if pv is None or pv <= 0 or m["value"] <= 0:
                    continue
                rows.setdefault(compound, []).append((math.log10(m["value"] / pv), sigma_of(m.get("uncertainty_pct"))))
        out = {}
        for compound, pairs in rows.items():
            w = np.array([1.0 / s**2 for _, s in pairs])
            d = np.array([x for x, _ in pairs])
            mean = float(np.sum(w * d) / np.sum(w))
            se = float(1.0 / math.sqrt(np.sum(w)))
            spread = float(np.sqrt(np.sum(w * (d - mean) ** 2) / np.sum(w))) if len(pairs) > 1 else 0.0
            out[compound] = C.ResponseFactor(compound, mean, max(se, spread / math.sqrt(len(pairs)), SIGMA_FLOOR_LOG10), len(pairs))
        return out

    # -- stage two ---------------------------------------------------------------------------
    def fit_overrides(self, lanes: Sequence[str], *, max_coordinates: int = 4, max_nfev: int = 40,
                      coordinates: Optional[Sequence[str]] = None) -> Tuple[List[C.Override], List[str], Dict[str, Any]]:
        """``coordinates`` restricts the candidates to the named ones (the envelope's prior names, e.g.
        ``b8.k_fft_decay.log10_k_ref_145C``) for a laboratory that knows which step its pot differs in;
        without it every calibratable coordinate of the resolved lanes is a candidate and the contrasts
        choose. Contrasts from one series at one temperature alias a sink rate against formation and
        osone rates, so the unrestricted choice is the best-explaining set, not necessarily the true one;
        the card lists both the chosen and the unchosen."""
        contrasts = contrasts_of(self.fit)
        notes: List[str] = []
        diagnostics: Dict[str, Any] = {"contrasts": len(contrasts), "candidates": 0, "identified": [], "not_identified": []}
        if not contrasts:
            notes.append("no contrasts: every fit record differs from every other in no measured compound, so no kinetic coordinate moves")
            return [], notes, diagnostics
        candidates = [c for lane in lanes for c in C.candidate_coordinates(lane)]
        if coordinates:
            wanted = set(coordinates)
            unknown = wanted - {c[0].name for c in candidates}
            if unknown:
                raise ValueError(f"unknown or non-calibratable coordinate(s) {sorted(unknown)}; the candidates are {[c[0].name for c in candidates]}")
            candidates = [c for c in candidates if c[0].name in wanted]
            notes.append(f"candidates restricted by the caller to {sorted(wanted)}")
        diagnostics["candidates"] = len(candidates)
        if not candidates:
            notes.append("the resolved lane offers no calibratable coordinate (no shipped uncertainty); only the response factors move")
            return [], notes, diagnostics
        x0 = np.array([c[1] for c in candidates])
        sig = np.array([c[2] for c in candidates])
        by_name = {r.name: r for r in self.fit}
        c_sigma = np.array([c.sigma for c in contrasts])
        c_meas = np.array([c.log10_ratio for c in contrasts])
        needed = sorted({c.a for c in contrasts} | {c.b for c in contrasts})

        def contrast_residuals(x: np.ndarray, subset: Sequence[int]) -> np.ndarray:
            overrides = tuple(
                C.Override(candidates[i][0], float(x0[i]), float(sig[i]), float(x[j]), float(sig[i]), candidates[i][3])
                for j, i in enumerate(subset) if abs(float(x[j]) - float(x0[i])) > 0.0
            )
            cal = C.Calibration("_fit", "", "", "", {}, overrides, (), ()) if overrides else None
            preds = self.predict([by_name[n] for n in needed], cal)
            out = np.zeros(len(contrasts))
            for k, c in enumerate(contrasts):
                pa, pb = preds[c.a].values.get(c.compound), preds[c.b].values.get(c.compound)
                if pa is None or pb is None or pa <= 0 or pb <= 0:
                    out[k] = 0.0        # a refused arm contributes nothing rather than a fake residual
                    continue
                out[k] = (c_meas[k] - math.log10(pb / pa)) / c_sigma[k]
            return out

        self._fitting = True
        try:
            r0 = contrast_residuals(x0, range(len(candidates)))
            # Jacobian of the contrast residuals at the shipped values, one forward difference per candidate
            J = np.zeros((len(contrasts), len(candidates)))
            for i in range(len(candidates)):
                h = min(0.5, max(0.02, 0.25 * float(sig[i])))
                xi = x0.copy()
                xi[i] += h
                J[:, i] = (contrast_residuals(xi, range(len(candidates))) - r0) / h
            # identifiability: a column must move a contrast by a detectable amount at all (an absolute
            # floor, so floating-point noise never selects a coordinate), then pivoted QR on the scaled
            # Jacobian keeps the leading pivots
            scaled = J * sig[None, :]
            norms = np.linalg.norm(scaled, axis=0)
            usable = [i for i in range(len(candidates)) if norms[i] > IDENTIFIABILITY_FLOOR]
            selected: List[int] = []
            if usable:
                from scipy.linalg import qr

                _, R, piv = qr(scaled[:, usable], mode="economic", pivoting=True)
                lead = abs(R[0, 0]) if R.size else 0.0
                for k in range(min(R.shape[0], R.shape[1], len(contrasts), max_coordinates)):
                    if lead > IDENTIFIABILITY_FLOOR and abs(R[k, k]) >= IDENTIFIABILITY_TOL * lead:
                        selected.append(int(usable[piv[k]]))
            diagnostics["identified"] = [candidates[i][0].name for i in selected]
            diagnostics["not_identified"] = [candidates[i][0].name for i in range(len(candidates)) if i not in selected]
            if not selected:
                notes.append("the contrasts identify no kinetic coordinate (Jacobian rank test); only the response factors move")
                return [], notes, diagnostics
            from scipy.optimize import least_squares

            lo = np.array([candidates[i][3][0] if candidates[i][3] else x0[i] - 4 * sig[i] for i in selected])
            hi = np.array([candidates[i][3][1] if candidates[i][3] else x0[i] + 4 * sig[i] for i in selected])
            start = np.clip(x0[selected], lo + 1e-9, hi - 1e-9)

            def full_residuals(x: np.ndarray) -> np.ndarray:
                prior = (x - x0[selected]) / sig[selected]
                return np.concatenate([contrast_residuals(x, selected), prior])

            result = least_squares(full_residuals, start, bounds=(lo, hi), x_scale=sig[selected], max_nfev=max_nfev, diff_step=0.05)
            # posterior sigma from the Jacobian at the solution (contrast part + prior part)
            Jsel = result.jac
            info = Jsel.T @ Jsel
            try:
                cov = np.linalg.pinv(info)
                post = np.sqrt(np.clip(np.diag(cov), 0.0, None))
            except np.linalg.LinAlgError:
                post = sig[selected]
            overrides = [
                C.Override(candidates[i][0], float(x0[i]), float(sig[i]), float(result.x[j]),
                           float(min(post[j], sig[i])) if post[j] > 0 else float(sig[i]), candidates[i][3])
                for j, i in enumerate(selected)
            ]
            diagnostics.update({"cost_before": float(0.5 * np.sum(r0**2)), "cost_after": float(result.cost),
                                "nfev": int(result.nfev), "status": int(result.status)})
            return overrides, notes, diagnostics
        finally:
            self._fitting = False

    # -- the hold-out ------------------------------------------------------------------------
    def holdout_scores(self, calibration: Optional[C.Calibration]) -> Dict[str, Any]:
        preds = self.predict(self.validate, calibration)
        rows = []
        for r in self.validate:
            p = preds[r.name]
            for compound, m in r.measured.items():
                pv = p.values.get(compound)
                if pv is None or pv <= 0 or m["value"] <= 0:
                    rows.append({"record": r.name, "compound": compound, "measured_ppb": m["value"], "predicted_ppb": pv, "fold_error": None})
                    continue
                fold = max(pv / m["value"], m["value"] / pv)
                rows.append({"record": r.name, "compound": compound, "measured_ppb": m["value"], "predicted_ppb": pv, "fold_error": fold})
        folds = [x["fold_error"] for x in rows if x["fold_error"] is not None]
        return {"rows": rows, "n": len(folds), "median_fold": float(np.median(folds)) if folds else None,
                "within_3x": sum(1 for f in folds if f <= FOLD_PASS)}


def calibrate(document: Mapping[str, Any], lab: str, *, predict_fn: PredictFn = engine_predict,
              max_coordinates: int = 4, max_nfev: int = 40,
              coordinates: Optional[Sequence[str]] = None) -> Tuple[C.Calibration, Dict[str, Any]]:
    """The whole procedure; returns the calibration and the card payload."""
    from src import provenance

    records = records_from_document(document)
    fit, validate, role_note = assign_roles(records)
    matrices = sorted({str(r.spec.get("matrix") or "water") for r in records})
    fitter = Fitter(fit, validate, predict_fn)
    lanes = fitter.lanes()
    notes = [role_note]
    if len(matrices) > 1:
        notes.append(f"the records span more than one matrix ({', '.join(matrices)}); a calibration is per matrix, the first is recorded")
    before = fitter.holdout_scores(None) if validate else None
    factors_shipped = fitter.response_factors(None)
    overrides, fit_notes, diagnostics = fitter.fit_overrides(lanes, max_coordinates=max_coordinates, max_nfev=max_nfev, coordinates=coordinates)
    notes += fit_notes
    partial = C.Calibration(lab, C.shipped_wave(), C.today(), matrices[0] if matrices else "water", {}, tuple(overrides), tuple(r.name for r in fit), tuple(r.name for r in validate))
    factors = fitter.response_factors(partial if overrides else None)
    leaked = fitter.reads_during_fit & {r.name for r in validate}
    if leaked:
        raise RuntimeError(f"hold-out records were read during the fit: {sorted(leaked)}")
    try:
        prov = provenance.provenance_block("maillard_calibration", generated_by="src/kinetic_core/user_fit.py")
    except Exception:  # noqa: BLE001 - provenance is informative, never blocking
        prov = {}
    cal = C.Calibration(lab, partial.base_wave, partial.created, partial.matrix, factors, tuple(overrides),
                        partial.fit_records, partial.validate_records, tuple(notes), prov)
    after = fitter.holdout_scores(cal) if validate else None
    card = {
        "artifact": "maillard_calibration_card",
        "lab": lab,
        "base_wave": cal.base_wave,
        "created": cal.created,
        "matrix": cal.matrix,
        "lanes": lanes,
        "records": {"fit": list(cal.fit_records), "validate": list(cal.validate_records)},
        "response_factors_at_shipped_kinetics": {c: rf.log10 for c, rf in factors_shipped.items()},
        "response_factors": {c: {"factor": rf.factor, "log10": rf.log10, "sigma_log10": rf.sigma, "n_rows": rf.n_rows} for c, rf in factors.items()},
        "overrides": [o.__dict__ | {"coordinate": o.coordinate.name, "shift": o.shift} for o in overrides],
        "diagnostics": diagnostics,
        "holdout": {"before": before, "after": after},
        "notes": notes,
        "evaluations": fitter.evaluations,
        "reads_during_fit": sorted(fitter.reads_during_fit),
    }
    return cal, card


def render_card(card: Mapping[str, Any]) -> str:
    out = ["=" * 96, f"  CALIBRATION CARD   laboratory: {card['lab']}   base: {card['base_wave']}   matrix: {card['matrix']}", "=" * 96, ""]
    out.append(f"  records: fit {', '.join(card['records']['fit']) or '-'}")
    out.append(f"           validate {', '.join(card['records']['validate']) or '- (unvalidated)'}")
    out.append("")
    out.append("  RESPONSE FACTORS  (from the fit records' levels; a property of the measurement, not the chemistry)")
    if not card["response_factors"]:
        out.append("    none: no fit record's compound was answered by the engine")
    for c, rf in card["response_factors"].items():
        out.append(f"    {c:<32} x{rf['factor']:.3g}   (log10 {rf['log10']:+.2f} +/- {rf['sigma_log10']:.2f}, {rf['n_rows']} rows)")
    out.append("")
    out.append("  KINETIC OVERRIDES  (from contrasts only; pulled toward the shipped value by its shipped sigma)")
    if not card["overrides"]:
        out.append("    none moved")
    for o in card["overrides"]:
        out.append(f"    {o['coordinate']:<40} {o['prior_value']:+.3f} -> {o['value']:+.3f}   (shift {o['shift']:+.2f}; prior sigma {o['prior_sigma']:.2f}, posterior {o['sigma']:.2f})")
    d = card["diagnostics"]
    out.append(f"    contrasts {d.get('contrasts', 0)}, candidates {d.get('candidates', 0)}, identified {len(d.get('identified', []))}")
    if d.get("not_identified"):
        out.append("    not identified by the contrasts (left at the shipped value): " + ", ".join(d["not_identified"][:12]) + (" ..." if len(d["not_identified"]) > 12 else ""))
    out.append("")
    h = card["holdout"]
    out.append("  HOLD-OUT  (validate records, never fitted)")
    if not h["before"]:
        out.append("    none held out: the calibration is unvalidated")
    else:
        b, a = h["before"], h["after"]
        out.append(f"    median fold error   shipped {b['median_fold']:.3g}   calibrated {a['median_fold']:.3g}   (within 3x: {b['within_3x']}/{b['n']} -> {a['within_3x']}/{a['n']})")
        for rb, ra in zip(b["rows"], a["rows"]):
            fb = f"{rb['fold_error']:.3g}x" if rb["fold_error"] else "refused"
            fa = f"{ra['fold_error']:.3g}x" if ra["fold_error"] else "refused"
            out.append(f"      {rb['record']:<28} {rb['compound']:<26} measured {rb['measured_ppb']:.3g}   shipped {fb:>9}   calibrated {fa:>9}")
    out.append("")
    for n in card["notes"]:
        out.append(f"  note: {n}")
    out.append("")
    out.append("  The shipped parameters are untouched. Apply this file with --calibration on compare, predict or score.")
    return "\n".join(out)
