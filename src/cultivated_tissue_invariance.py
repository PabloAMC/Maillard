"""
Does the kinetic layer change what to add to cultivated tissue?

The pre-registration is ``results/validation/cultivated_tissue_invariance_prereg.md``; read its
sections 3 and 6 before this file. Nothing here predicts: the engine is asked through the front
door (``src.comparative_cli.spec_to_core`` -> ``kinetic_core.engine.predict``), exactly as a
user's spec would be, once per arm per draw.

WHAT IS COMPUTED
----------------
For each draw of the composition box (section 6, A7): a beef reference and a cultivated-muscle
composition, each precursor log-uniform in its declared range. Then

* **N**, the naive ranking: candidate precursors ordered by fractional deficit
  ``1 - cultivated / beef`` (A4), ties by absolute deficit; a precursor at or above its beef
  level has nothing to restore and is dropped from that draw.
* **E**, the engine ranking: each restorable precursor set to its beef level in turn, the
  programme integrated, precursors ordered by the change in summed odour-activity over the
  targets that carry a measured water threshold (A3: MFT, FFT, furfural).
* ``top(E) == top(N)``, Kendall tau over the restorable set, the concentration ratios of the
  no-threshold targets, and every refusal the engine issued.

For every draw in which the tops disagree, the parameter envelope (A8) re-runs the two arms
under joint draws of the core's priors and records whether the same precursor still wins.

The verdict is then read off section 3 mechanically. The thresholds are the pre-registration's;
this file does not own them and does not adjust them.

THE COMPOSITION BOX IS A SENSITIVITY DEVICE, NOT DATA
-----------------------------------------------------
Every range below is labelled ``stub``. None is a measurement. They are wide on purpose: if the
ordering is stable across the whole box the values never mattered, and if it flips inside the
box the finding is that composition must be measured before anything is predicted. The box is
echoed into the artifact directory so the record shows exactly what was swept; it is declared
here, in code, next to the other declared assumptions, rather than under ``data/`` where it
would read as curated.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from src import artifact_io, data_paths, provenance
from src.comparative_cli import spec_to_core, validate_spec
from src.kinetic_core import engine
from src.kinetic_core.uncertainty import sample_draws
from src.report_format import fmt_number as _fmt

ARTIFACT = "cultivated_tissue_invariance"
OUTPUT_DIR = data_paths.RESULTS_ROOT / "cultivated_tissue_invariance"
OUTPUT_JSON = OUTPUT_DIR / "cultivated_tissue_invariance.json"
BOX_ECHO = OUTPUT_DIR / "composition_box.yml"
PREREG = data_paths.VALIDATION_DIR / "cultivated_tissue_invariance_prereg.md"

# ---------------------------------------------------------------------------
# Declared before the run (pre-registration sections 2 and 6)
# ---------------------------------------------------------------------------

#: The candidates the engine can take on the sulfur arm (A1, A2). Order is cosmetic.
CANDIDATES: Tuple[str, ...] = ("ribose", "cysteine", "thiamine", "glucose")

#: Declared candidates the engine refuses (A2), with the engine's own reason class.
UNRANKABLE: Mapping[str, str] = {
    "leucine": "UNMAPPED PRECURSOR: not a species in any core lane",
    "IMP": "UNMAPPED PRECURSOR: not a species in any core lane",
    "ribose-5-phosphate": "UNMAPPED PRECURSOR: not a species in any core lane",
}

#: The sulfur arm's targets (A1). Requested by the names the front door resolves.
TARGETS: Tuple[str, ...] = (
    "2-methyl-3-furanthiol",
    "2-furfurylthiol",
    "furfural",
    "hydrogen sulfide",
    "methanethiol",
    "furaneol",
)
#: Those with a measured water threshold in the corpus: the decision metric (A3).
METRIC_TARGETS: Tuple[str, ...] = ("2-methyl-3-furanthiol", "2-furfurylthiol", "furfural")

#: Two programmes, as declared (section 2, A6).
PROGRAMMES: Tuple[Tuple[str, float, float], ...] = (
    ("100C_20min", 100.0, 20.0),
    ("140C_5min", 140.0, 5.0),
)
FIXED_CONDITIONS: Mapping[str, Any] = {
    "ph": 6.0,
    "aw": 0.98,
    "matrix": "water",
    "protein_type": "free",
    "buffer": {
        "kind": "phosphate",
        "phosphate_mol_l": 0.03,
        "source": "declared (A6): the muscle phosphate pool, order of magnitude; identical in every arm",
    },
}

N_DRAWS = 200
SEED = 0
N_ENVELOPE = 50


@dataclass(frozen=True)
class Range:
    lo_mM: float
    hi_mM: float
    label: str  # "stub" or "sourced"
    note: str


#: THE COMPOSITION BOX. Every entry is a stub (see the module docstring). mM in tissue water.
BOX: Mapping[str, Mapping[str, Range]] = {
    "beef": {
        "ribose": Range(0.3, 5.0, "stub", "post-mortem ribose from IMP breakdown in aged beef; order of magnitude only"),
        "cysteine": Range(0.05, 0.5, "stub", "free cysteine in raw muscle; order of magnitude only"),
        "thiamine": Range(0.001, 0.01, "stub", "beef thiamine ~0.05-0.15 mg/100 g; order of magnitude only"),
        "glucose": Range(1.0, 10.0, "stub", "free glucose in post-mortem muscle; order of magnitude only"),
        "leucine": Range(0.3, 1.5, "stub", "free leucine; unrankable by the engine, reported for the naive ranking only"),
        "IMP": Range(1.0, 8.0, "stub", "inosine monophosphate in aged beef; unrankable by the engine"),
        "ribose-5-phosphate": Range(0.01, 0.1, "stub", "unrankable by the engine"),
    },
    "cultivated_muscle": {
        "ribose": Range(0.01, 2.0, "stub", "no ageing step; IMP-to-ribose route depends on post-harvest handling; deliberately wide"),
        "cysteine": Range(0.02, 0.5, "stub", "medium-fed cells; deliberately wide"),
        "thiamine": Range(0.0003, 0.01, "stub", "DMEM carries ~12 uM thiamine; intracellular pool unknown; deliberately wide"),
        "glucose": Range(0.1, 10.0, "stub", "depends on the harvest wash; deliberately wide"),
        "leucine": Range(0.2, 3.0, "stub", "medium is leucine-rich; unrankable by the engine"),
        "IMP": Range(0.1, 3.0, "stub", "unrankable by the engine"),
        "ribose-5-phosphate": Range(0.01, 0.1, "stub", "unrankable by the engine"),
    },
}

#: Section 3, as declared and amended (A5). Read, never edited, here.
THRESHOLDS: Mapping[str, float] = {
    "T1_min_disagreement_fraction": 0.20,
    "T1_min_concentration_share": 0.50,
    "T1_min_envelope_survival": 0.80,
    "T2_min_top_agreement": 0.90,
    "T2_min_mean_tau": 0.75,
    "T3_min_refusal_fraction": 0.50,
}


# ---------------------------------------------------------------------------
# Rankings
# ---------------------------------------------------------------------------


def _log_uniform(rng: np.random.Generator, r: Range) -> float:
    return float(10 ** rng.uniform(math.log10(r.lo_mM), math.log10(r.hi_mM)))


def draw_compositions(rng: np.random.Generator) -> Tuple[Dict[str, float], Dict[str, float]]:
    beef = {p: _log_uniform(rng, r) for p, r in BOX["beef"].items()}
    cult = {p: _log_uniform(rng, r) for p, r in BOX["cultivated_muscle"].items()}
    return beef, cult


def naive_ranking(beef: Mapping[str, float], cult: Mapping[str, float], names: Sequence[str]) -> List[str]:
    """A4: fractional deficit, ties by absolute deficit; nothing to restore -> dropped."""
    restorable = [p for p in names if beef[p] > cult[p]]
    return sorted(restorable, key=lambda p: (-(1.0 - cult[p] / beef[p]), -(beef[p] - cult[p])))


def kendall_tau(a: Sequence[str], b: Sequence[str]) -> Optional[float]:
    """Kendall's tau-a between two orderings of the same set; None below two items."""
    items = list(a)
    if len(items) < 2 or set(items) != set(b):
        return None
    pa = {x: i for i, x in enumerate(a)}
    pb = {x: i for i, x in enumerate(b)}
    conc = disc = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            x, y = items[i], items[j]
            s = (pa[x] - pa[y]) * (pb[x] - pb[y])
            if s > 0:
                conc += 1
            elif s < 0:
                disc += 1
    n = len(items)
    return (conc - disc) / (n * (n - 1) / 2)


# ---------------------------------------------------------------------------
# The engine, through the front door
# ---------------------------------------------------------------------------


def _spec(precursors: Mapping[str, float], temp_c: float, time_min: float, name: str) -> Mapping[str, Any]:
    spec = {
        "name": name,
        "precursors": {p: float(v) for p, v in precursors.items() if p in CANDIDATES},
        "temp_C": float(temp_c),
        "time_min": float(time_min),
        "targets": list(TARGETS),
    }
    spec.update(FIXED_CONDITIONS)
    return validate_spec(spec, label=name)


@dataclass(frozen=True)
class ArmResult:
    answered: bool
    state: str
    metric: Optional[float]                 # summed OAV over METRIC_TARGETS
    ug_per_l: Mapping[str, float]           # every answered target
    oav: Mapping[str, Optional[float]]      # per metric target
    refusals: Tuple[str, ...]               # reasons + refused targets, verbatim


def run_arm(precursors: Mapping[str, float], temp_c: float, time_min: float, name: str,
            draw=None) -> ArmResult:
    core = spec_to_core(_spec(precursors, temp_c, time_min, name))
    if draw is None:
        run = engine.predict(core, list(TARGETS))
    else:
        run = engine.predict(core, list(TARGETS), draw=draw, size_declared_bands=False)
    decl = run.declaration
    refusals = tuple(decl.reasons) + tuple(
        f"refused target: {t}" for t in (run.run_metadata.get("refused_targets") or [])
    )
    if not run.answered:
        return ArmResult(False, decl.state, None, {}, {}, refusals)
    ug: Dict[str, float] = {}
    oav: Dict[str, Optional[float]] = {}
    for row in run.interval_rows():
        ug[row["compound"]] = float(row["predicted_ug_per_l"])
        if row["compound"] in METRIC_TARGETS:
            entry = row["oav"] or {}
            oav[row["compound"]] = entry.get("OAV_point")
    present = [v for v in oav.values() if v is not None]
    metric = float(sum(present)) if present else None
    if metric is None:
        refusals = refusals + ("no metric target carried an odour-activity value",)
    return ArmResult(True, decl.state, metric, ug, oav, refusals)


# ---------------------------------------------------------------------------
# One draw, one programme
# ---------------------------------------------------------------------------


def evaluate_draw(index: int, beef: Mapping[str, float], cult: Mapping[str, float],
                  temp_c: float, time_min: float) -> Dict[str, Any]:
    naive_all = naive_ranking(beef, cult, list(BOX["beef"].keys()))
    naive = naive_ranking(beef, cult, CANDIDATES)
    base = run_arm(cult, temp_c, time_min, f"draw{index}-cultivated")
    record: Dict[str, Any] = {
        "draw": index,
        "beef_mM": {p: beef[p] for p in BOX["beef"]},
        "cultivated_mM": {p: cult[p] for p in BOX["cultivated_muscle"]},
        "restorable": naive,
        "naive_ranking": naive,
        "naive_ranking_all_declared": naive_all,
        "baseline": {
            "answered": base.answered, "state": base.state, "metric": base.metric,
            "ug_per_l": base.ug_per_l, "refusals": list(base.refusals),
        },
        "arms": {},
        "engine_ranking": None,
        "top_agree": None,
        "kendall_tau": None,
        "engine_refused": (not base.answered) or base.metric is None,
    }
    if record["engine_refused"] or len(naive) < 2:
        return record
    deltas: Dict[str, float] = {}
    for p in naive:
        restored = dict(cult)
        restored[p] = beef[p]
        arm = run_arm(restored, temp_c, time_min, f"draw{index}-restore-{p}")
        record["arms"][p] = {
            "answered": arm.answered, "state": arm.state, "metric": arm.metric,
            "delta_metric": (arm.metric - base.metric) if (arm.answered and arm.metric is not None) else None,
            "ug_per_l": arm.ug_per_l, "refusals": list(arm.refusals),
        }
        if arm.answered and arm.metric is not None:
            deltas[p] = arm.metric - base.metric
        else:
            record["engine_refused"] = True
            return record
    ranking = sorted(naive, key=lambda p: -deltas[p])
    record["engine_ranking"] = ranking
    record["delta_metric"] = deltas
    record["top_agree"] = ranking[0] == naive[0]
    record["kendall_tau"] = kendall_tau(naive, ranking)
    return record


def envelope_on_reversal(record: Mapping[str, Any], temp_c: float, time_min: float,
                         n: int, seed: int) -> Dict[str, Any]:
    """A8: does the engine's winner still beat the naive winner under the core's priors?"""
    a, b = record["engine_ranking"][0], record["naive_ranking"][0]
    beef, cult = record["beef_mM"], record["cultivated_mM"]
    wins = refused = 0
    for d in sample_draws(n, seed=seed):
        base = run_arm(cult, temp_c, time_min, "env-base", draw=d.core)
        ra = dict(cult); ra[a] = beef[a]
        rb = dict(cult); rb[b] = beef[b]
        arm_a = run_arm(ra, temp_c, time_min, f"env-{a}", draw=d.core)
        arm_b = run_arm(rb, temp_c, time_min, f"env-{b}", draw=d.core)
        ok = all(x.answered and x.metric is not None for x in (base, arm_a, arm_b))
        if not ok:
            refused += 1
            continue
        if (arm_a.metric - base.metric) > (arm_b.metric - base.metric):
            wins += 1
    evaluated = n - refused
    return {
        "draw": record["draw"], "engine_top": a, "naive_top": b,
        "envelope_draws": n, "refused": refused,
        "survival": (wins / evaluated) if evaluated else None,
    }


# ---------------------------------------------------------------------------
# The sweep and the verdict
# ---------------------------------------------------------------------------


def sweep_programme(label: str, temp_c: float, time_min: float, *, n_draws: int, seed: int,
                    n_envelope: int) -> Dict[str, Any]:
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    t0 = time.time()
    draws = []
    for i in range(n_draws):
        beef, cult = draw_compositions(rng)
        draws.append(evaluate_draw(i, beef, cult, temp_c, time_min))
    evaluated = [d for d in draws if d["top_agree"] is not None]
    refused = [d for d in draws if d["engine_refused"]]
    too_few = [d for d in draws if not d["engine_refused"] and d["top_agree"] is None]
    disagree = [d for d in evaluated if not d["top_agree"]]
    taus = [d["kendall_tau"] for d in evaluated if d["kendall_tau"] is not None]

    pairs: Dict[str, int] = {}
    for d in disagree:
        key = f"{d['engine_ranking'][0]} over {d['naive_ranking'][0]}"
        pairs[key] = pairs.get(key, 0) + 1
    dominant_pair = max(pairs.items(), key=lambda kv: kv[1]) if pairs else None

    envelope = []
    if disagree:
        for k, d in enumerate(disagree):
            envelope.append(envelope_on_reversal(d, temp_c, time_min, n_envelope, seed=seed * 1000 + k + 1))
    survivals = [e["survival"] for e in envelope if e["survival"] is not None]
    # T1 asks about THE dominant reversal, not the average reversal.
    dominant_survivals = [
        e["survival"] for e in envelope
        if e["survival"] is not None and dominant_pair is not None
        and f"{e['engine_top']} over {e['naive_top']}" == dominant_pair[0]
    ]

    # Where does each precursor land, on average, in the two rankings?
    def mean_rank(key: str) -> Dict[str, Optional[float]]:
        out: Dict[str, Optional[float]] = {}
        for p in CANDIDATES:
            pos = [d[key].index(p) + 1 for d in evaluated if p in d[key]]
            out[p] = (sum(pos) / len(pos)) if pos else None
        return out

    # The refusals the engine issued, counted by text.
    refusal_counts: Dict[str, int] = {}
    for d in draws:
        texts = list(d["baseline"]["refusals"])
        for arm in d["arms"].values():
            texts += arm["refusals"]
        for t in set(texts):
            refusal_counts[t] = refusal_counts.get(t, 0) + 1

    summary = {
        "programme": label, "temp_C": temp_c, "time_min": time_min,
        "draws": n_draws,
        "evaluated": len(evaluated),
        "engine_refused": len(refused),
        "fewer_than_two_restorable": len(too_few),
        "refusal_fraction": len(refused) / n_draws,
        "top_agreement_fraction": (1 - len(disagree) / len(evaluated)) if evaluated else None,
        "disagreement_fraction": (len(disagree) / len(evaluated)) if evaluated else None,
        "mean_kendall_tau": (sum(taus) / len(taus)) if taus else None,
        "disagreement_pairs": pairs,
        "dominant_pair": {"pair": dominant_pair[0], "count": dominant_pair[1],
                          "share_of_disagreements": dominant_pair[1] / len(disagree)} if dominant_pair else None,
        "envelope_reversals_tested": len(envelope),
        "mean_envelope_survival": (sum(survivals) / len(survivals)) if survivals else None,
        "dominant_pair_envelope_survival": (sum(dominant_survivals) / len(dominant_survivals)) if dominant_survivals else None,
        "mean_rank_naive": mean_rank("naive_ranking"),
        "mean_rank_engine": mean_rank("engine_ranking"),
        "refusals_by_text": refusal_counts,
        "wall_seconds": round(time.time() - t0, 1),
    }
    return {"summary": summary, "draws": draws, "envelope": envelope}


def verdict(summary: Mapping[str, Any]) -> Dict[str, Any]:
    """Section 3, applied mechanically. Returns the verdict and each clause's truth value."""
    th = THRESHOLDS
    c: Dict[str, Optional[bool]] = {}
    c["T3_refusal"] = summary["refusal_fraction"] >= th["T3_min_refusal_fraction"]
    dis = summary["disagreement_fraction"]
    dom = summary["dominant_pair"]
    surv = summary["dominant_pair_envelope_survival"]
    c["T1_disagreement"] = (dis is not None) and dis >= th["T1_min_disagreement_fraction"]
    c["T1_concentrated"] = (dom is not None) and dom["share_of_disagreements"] >= th["T1_min_concentration_share"]
    c["T1_survives_envelope"] = (surv is not None) and surv >= th["T1_min_envelope_survival"]
    agree = summary["top_agreement_fraction"]
    tau = summary["mean_kendall_tau"]
    c["T2_top_agreement"] = (agree is not None) and agree >= th["T2_min_top_agreement"]
    c["T2_tau"] = (tau is not None) and tau >= th["T2_min_mean_tau"]
    if c["T3_refusal"]:
        result = "T3"
    elif c["T1_disagreement"] and c["T1_concentrated"] and c["T1_survives_envelope"]:
        result = "T1"
    elif c["T2_top_agreement"] and c["T2_tau"]:
        result = "T2"
    else:
        result = "indeterminate"
    return {"verdict": result, "clauses": c, "thresholds": dict(th)}


def build(*, n_draws: int = N_DRAWS, seed: int = SEED, n_envelope: int = N_ENVELOPE) -> Dict[str, Any]:
    programmes = []
    for label, temp_c, time_min in PROGRAMMES:
        result = sweep_programme(label, temp_c, time_min, n_draws=n_draws, seed=seed, n_envelope=n_envelope)
        result["verdict"] = verdict(result["summary"])
        programmes.append(result)
    box_echo = {
        tissue: {p: {"lo_mM": r.lo_mM, "hi_mM": r.hi_mM, "label": r.label, "note": r.note} for p, r in rows.items()}
        for tissue, rows in BOX.items()
    }
    return {
        "provenance": provenance.provenance_block(
            ARTIFACT, generated_by="src/cultivated_tissue_invariance.py", inputs=[PREREG]
        ),
        "pre_registration": data_paths.rel(PREREG),
        "design": {
            "candidates": list(CANDIDATES),
            "unrankable_declared_candidates": dict(UNRANKABLE),
            "targets": list(TARGETS),
            "metric_targets": list(METRIC_TARGETS),
            "fixed_conditions": dict(FIXED_CONDITIONS),
            "n_draws": n_draws, "seed": seed, "n_envelope": n_envelope,
            "composition_box": box_echo,
            "structural_refusals": [
                "trunk arm (methional, pyrazines, furaneol on glucose + glycine): methional refused, pyrazines ~1e-13 ug/L, no water threshold on any target; no decision metric (A1)",
                "cultivated fat: no route from tissue lipid to any lane (A1)",
                "hexanal: needs a lipid carrier; no tissue-fat carrier exists (A1)",
                "leucine, IMP, ribose-5-phosphate: not species in any core lane (A2)",
            ],
        },
        "programmes": programmes,
        "overall": {
            "verdicts": {p["summary"]["programme"]: p["verdict"]["verdict"] for p in programmes},
        },
    }


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------


def _pct(x: Optional[float]) -> str:
    return "—" if x is None else f"{100 * x:.0f} %"


def render_markdown(payload: Mapping[str, Any]) -> str:
    L: List[str] = []
    L.append("# Does the kinetic layer change what to add to cultivated tissue?")
    L.append("")
    L.append(f"*Generated {payload['provenance']['generated_on']} by `{payload['provenance']['generated_by']}`. "
             f"Pre-registration: `{payload['pre_registration']}`. The composition box is a sensitivity device; "
             "every range is a stub, none is a measurement.*")
    L.append("")
    d = payload["design"]
    L.append("## Verdict, by programme")
    L.append("")
    L.append("| programme | verdict | draws evaluated | engine refused | < 2 restorable | top(E)=top(N) | mean τ | disagreements | dominant reversal | its envelope survival |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for p in payload["programmes"]:
        s, v = p["summary"], p["verdict"]
        dom = s["dominant_pair"]
        L.append(
            f"| {s['programme']} | **{v['verdict']}** | {s['evaluated']} / {s['draws']} | {s['engine_refused']} | "
            f"{s['fewer_than_two_restorable']} | "
            f"{_pct(s['top_agreement_fraction'])} | {_fmt(s['mean_kendall_tau']) if s['mean_kendall_tau'] is not None else '—'} | "
            f"{_pct(s['disagreement_fraction'])} | "
            f"{(dom['pair'] + ' (' + _pct(dom['share_of_disagreements']) + ' of disagreements)') if dom else '—'} | "
            f"{_pct(s['dominant_pair_envelope_survival'])} |"
        )
    L.append("")
    L.append("Read against section 3 of the pre-registration: T1 needs ≥ 20 % disagreement, one reversal carrying "
             "≥ 50 % of it, and ≥ 80 % envelope survival; T2 needs ≥ 90 % top agreement and mean τ ≥ 0.75; T3 is "
             "≥ 50 % of draws refused. Anything else is indeterminate and is not a licence to build.")
    L.append("")
    L.append("## Where each precursor lands")
    L.append("")
    for p in payload["programmes"]:
        s = p["summary"]
        L.append(f"**{s['programme']}** (mean position, 1 = restore first; over the draws the engine answered)")
        L.append("")
        L.append("| precursor | naive N | engine E |")
        L.append("|---|---|---|")
        for c in d["candidates"]:
            n, e = s["mean_rank_naive"].get(c), s["mean_rank_engine"].get(c)
            L.append(f"| {c} | {_fmt(n) if n is not None else '—'} | {_fmt(e) if e is not None else '—'} |")
        L.append("")
        if s["disagreement_pairs"]:
            L.append("Disagreements, engine's winner over naive winner: " +
                     "; ".join(f"{k} ×{v}" for k, v in sorted(s["disagreement_pairs"].items(), key=lambda kv: -kv[1])) + ".")
            L.append("")
    L.append("## What the engine refused")
    L.append("")
    L.append("Structural, before any draw:")
    L.append("")
    for r in d["structural_refusals"]:
        L.append(f"- {r}")
    L.append("")
    for p in payload["programmes"]:
        s = p["summary"]
        if s["refusals_by_text"]:
            L.append(f"During the {s['programme']} sweep (draws in which the text appeared at least once):")
            L.append("")
            for t, n in sorted(s["refusals_by_text"].items(), key=lambda kv: -kv[1]):
                L.append(f"- ×{n}: {t[:220]}{'…' if len(t) > 220 else ''}")
            L.append("")
    L.append("## The composition box that was swept")
    L.append("")
    L.append("| tissue | precursor | lo mM | hi mM | label | note |")
    L.append("|---|---|---|---|---|---|")
    for tissue, rows in d["composition_box"].items():
        for c, r in rows.items():
            L.append(f"| {tissue} | {c} | {r['lo_mM']} | {r['hi_mM']} | {r['label']} | {r['note']} |")
    L.append("")
    L.append(f"Fixed in every spec: pH {d['fixed_conditions']['ph']}, a_w {d['fixed_conditions']['aw']}, matrix "
             f"{d['fixed_conditions']['matrix']}, phosphate {d['fixed_conditions']['buffer']['phosphate_mol_l']} mol/L. "
             f"{d['n_draws']} draws per programme, seed {d['seed']}, {d['n_envelope']} envelope draws per reversal.")
    L.append("")
    return "\n".join(L)


def write(payload: Mapping[str, Any]) -> Tuple[Any, Any]:
    import yaml

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    BOX_ECHO.write_text(
        "# Echo of the declared composition box (src/cultivated_tissue_invariance.py). Every range is a stub.\n"
        "# This file is generated; edit the module, not this file.\n"
        + yaml.safe_dump(payload["design"]["composition_box"], sort_keys=False),
        encoding="utf-8",
    )
    return artifact_io.write_artifact(payload, OUTPUT_JSON, render=render_markdown)


__all__ = [
    "BOX", "CANDIDATES", "METRIC_TARGETS", "OUTPUT_JSON", "PROGRAMMES", "TARGETS", "THRESHOLDS",
    "build", "evaluate_draw", "kendall_tau", "naive_ranking", "render_markdown", "verdict", "write",
]
