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

TWO BOXES
---------
``STUB_BOX`` is the box the declared run of 2026-09-13 swept: every range a stub, none a
measurement, wide on purpose. It is frozen here because the artifact of that run echoes it and
the two must stay identical (a test checks). ``BOX`` is the current box: the same shape, with a
range replaced wherever the literature read of 2026-09-14 found a measurement, each such range
naming the extraction dossier it came from. Three labels:

* ``stub``      -- no measurement found; the range is a sensitivity device and nothing else;
* ``secondary`` -- a number read from a review or an abstract, not from the table that measured it;
* ``sourced``   -- read from the measuring paper's own table. The dossier's "Source on disk" line says
                   whether that table was read from the PDF by eye. As of the evening of 2026-09-14 every
                   dossier the box cites has been (koutsidis2008a, koutsidis2008b, bischof2023, kim2024b,
                   muroya2019, joo2022, lombardiboccia2005).

The box is declared here, in code, next to the other declared assumptions, rather than under
``data/`` where it would read as curated. Both boxes are echoed into the artifact directory.
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


LABELS: Tuple[str, ...] = ("stub", "secondary", "sourced")


@dataclass(frozen=True)
class Range:
    lo_mM: float
    hi_mM: float
    label: str  # one of LABELS
    note: str
    #: extraction dossier stems under data/lit/extraction_dossiers/ (empty for a stub)
    dossiers: Tuple[str, ...] = ()


#: THE BOX THE DECLARED RUN SWEPT (2026-09-13). Frozen. Every entry a stub. mM in tissue water.
STUB_BOX: Mapping[str, Mapping[str, Range]] = {
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

#: THE CURRENT BOX (literature read of 2026-09-14). mM in tissue water. A stub here is a stub
#: because the read found NO measurement of that pool in that tissue; the gap map in the artifact
#: lists them. Beef moisture taken as 75 %, the pig construct of kim2024b as 90 % (its text and Fig. 3A).
#: Beef sugar and amino-acid ranges were re-set on 2026-09-14 from the PDFs of Koutsidis 2008a/b and the
#: corrected Table 1 of Bischof 2023 (pre-registration section 9).
BOX: Mapping[str, Mapping[str, Range]] = {
    "beef": {
        "ribose": Range(0.33, 2.2, "sourced",
                        "Koutsidis 2008b Table 2: 0.25 mmol/kg at day 1 to 1.67 at day 21 (n = 16 steers, GC-MS), "
                        "0.33-2.2 mM at 75 % moisture; Koutsidis 2008a Table 1: 0.57-1.08 mmol/kg across 30 steers at "
                        "10 d (0.76-1.44 mM) sits inside; the cited 0.26 mg/g point (Aliani 2013 via Hwang 2026, 2.3 mM) "
                        "sits at the top edge and is no longer a corner",
                        ("koutsidis2008b", "koutsidis2008a")),
        "cysteine": Range(0.002, 0.23, "sourced",
                          "Muroya 2019 Table 1: 1.6 nmol/g at D0 to 107 nmol/g at D14, n = 3 steers, CE-TOFMS (D0 at the "
                          "detection floor); Koutsidis 2008b Table 4: 0.05-0.16 mmol/kg over 21 d, n = 16, GC-MS; "
                          "Koutsidis 2008a Table 3: 0.05-0.17 mmol/kg across 30 steers at 10 d, whose top is the upper "
                          "corner (0.23 mM); the span is ageing plus animal spread",
                          ("muroya2019", "koutsidis2008b", "koutsidis2008a")),
        "thiamine": Range(0.00044, 0.0040, "sourced",
                          "Lombardi-Boccia 2005 Table 2: total thiamine 0.01-0.08 (+/- 0.01) mg/100 g across five raw beef "
                          "cuts by HPLC after acid hydrolysis, lowest mean to highest mean + SD; not detected in any cut after "
                          "cooking; about twofold below the stub",
                          ("lombardiboccia2005",)),
        "glucose": Range(2.4, 15.0, "sourced",
                         "Bischof 2023 Table 1 with the alpha- and beta-glucose rows SUMMED (the first read took one anomer "
                         "row as the total): 4.43 +/- 2.61 to 10.01 +/- 1.42 umol/g wet across two breeds and 28 d, mean "
                         "-/+ SD = 2.4-15 mM; Koutsidis 2008b Table 2 (7.33-10.3 mmol/kg, 9.8-13.7 mM) and 2008a Table 1 "
                         "(6.94-10.6 mmol/kg across 30 steers) sit inside; the cited 1.48 mg/g (11 mM) too",
                         ("bischof2023", "koutsidis2008b", "koutsidis2008a")),
        "leucine": Range(0.29, 3.2, "sourced",
                         "Muroya 2019: 263-827 nmol/g (0.35-1.1 mM); Bischof 2023 (corrected rows): 0.29 +/- 0.07 to "
                         "1.53 +/- 0.52 umol/g (0.29-2.7 mM); Koutsidis 2008b: 0.43-1.75 mmol/kg; Koutsidis 2008a: "
                         "0.78-2.40 mmol/kg across 30 steers (to 3.2 mM); unrankable by the engine",
                         ("muroya2019", "bischof2023", "koutsidis2008b", "koutsidis2008a")),
        "IMP": Range(0.1, 10.0, "sourced",
                     "Muroya 2019: 78 nmol/g pre-rigor to 7574 nmol/g at D1 (0.10-10 mM); Bischof 2023 (corrected rows, "
                     "1.1-4.8 mM), Koutsidis 2008b (3.5-8.4 mM) and 2008a (3.3-5.9 mM) inside it; unrankable by the engine",
                     ("muroya2019", "bischof2023", "koutsidis2008b", "koutsidis2008a")),
        "ribose-5-phosphate": Range(0.005, 0.1, "sourced",
                                    "Muroya 2019: non-detect at D0, 57-70 nmol/g at D1-D14 (to 0.093 mM); Koutsidis 2008b: "
                                    "0.04 mmol/kg flat over 21 d (0.053 mM); lower corner set at 0.005 because the draw is "
                                    "log-uniform; unrankable by the engine",
                                    ("muroya2019", "koutsidis2008b")),
    },
    "cultivated_muscle": {
        "ribose": Range(0.01, 2.0, "stub",
                        "NO MEASUREMENT FOUND in cultured muscle of any species (read of 2026-09-14); range unchanged from the stub box"),
        "cysteine": Range(0.02, 0.5, "stub",
                          "NO MEASUREMENT FOUND: Joo 2022 prints cysteine only as a percent of total amino acids with no absolute "
                          "total and no free/hydrolysed statement; Kim 2024b's free-amino-acid table omits it; range unchanged"),
        "thiamine": Range(0.0003, 0.01, "stub",
                          "NO MEASUREMENT FOUND in cultured muscle; DMEM carries ~12 uM thiamine HCl, what a washed construct "
                          "retains is unknown; range unchanged"),
        "glucose": Range(0.1, 10.0, "stub",
                         "NO MEASUREMENT FOUND in cultured muscle; note Joo 2022 proliferated in glucose-free DMEM; range unchanged"),
        "leucine": Range(0.1, 1.0, "sourced",
                         "Kim 2024b Table 4: 35.8 mg/kg free leucine in a pig gelatin-scaffold construct at 90 % moisture "
                         "(0.30 mM), CONFOUNDED by 49.6 mg/kg in the scaffold-only arm; one point, pig, widened tenfold; "
                         "unrankable by the engine",
                         ("kim2024b",)),
        "IMP": Range(0.0003, 3.0, "sourced",
                     "two primaries four decades apart: Kim 2024b pig construct 0.11 mg/kg (0.00035 mM); Joo 2022 bovine "
                     "2D tissue 1.98 mmol/kg (2.6 mM); the box carries both corners; unrankable by the engine",
                     ("kim2024b", "joo2022")),
        "ribose-5-phosphate": Range(0.01, 0.1, "stub",
                                    "NO MEASUREMENT FOUND in cultured muscle; range unchanged; unrankable by the engine"),
    },
}

#: For every stub in the current box: the measurement that would close it. One experiment covers
#: all of them (see the artifact's gap map).
CLOSES_WITH: Mapping[str, str] = {
    "ribose": "free ribose by GC-MS (oxime-TMS) or enzymatic assay on the washed construct extract; CE-MS panels do not carry it",
    "cysteine": "free cysteine by CE-TOFMS with thiol protection at extraction, as Muroya 2019 did for beef; report cystine alongside",
    "thiamine": "thiamine and its phosphates by HPLC-fluorescence (thiochrome) on the same extract",
    "glucose": "free glucose by enzymatic assay or GC-MS on the same extract; state the harvest wash",
    "ribose-5-phosphate": "on the CE-TOFMS panel with cysteine; it was quantified in beef by that method",
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


def draw_compositions(rng: np.random.Generator, box: Mapping[str, Mapping[str, Range]]
                      ) -> Tuple[Dict[str, float], Dict[str, float]]:
    beef = {p: _log_uniform(rng, r) for p, r in box["beef"].items()}
    cult = {p: _log_uniform(rng, r) for p, r in box["cultivated_muscle"].items()}
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
    naive_all = naive_ranking(beef, cult, list(beef.keys()))
    naive = naive_ranking(beef, cult, CANDIDATES)
    base = run_arm(cult, temp_c, time_min, f"draw{index}-cultivated")
    record: Dict[str, Any] = {
        "draw": index,
        "beef_mM": dict(beef),
        "cultivated_mM": dict(cult),
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
                    n_envelope: int, box: Mapping[str, Mapping[str, Range]] = BOX) -> Dict[str, Any]:
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    t0 = time.time()
    draws = []
    for i in range(n_draws):
        beef, cult = draw_compositions(rng, box)
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


def box_echo(box: Mapping[str, Mapping[str, Range]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    return {
        tissue: {p: {"lo_mM": r.lo_mM, "hi_mM": r.hi_mM, "label": r.label, "note": r.note,
                     "dossiers": list(r.dossiers)} for p, r in rows.items()}
        for tissue, rows in box.items()
    }


def gap_map(box: Mapping[str, Mapping[str, Range]]) -> List[Dict[str, Any]]:
    """Every stub in ``box``: what it is, and the measurement that would close it."""
    out = []
    for tissue, rows in box.items():
        for p, r in rows.items():
            if r.label == "stub":
                out.append({"tissue": tissue, "precursor": p, "lo_mM": r.lo_mM, "hi_mM": r.hi_mM,
                            "engine_rankable": p in CANDIDATES, "note": r.note,
                            "closes_with": CLOSES_WITH.get(p, "")})
    return out


def build(*, n_draws: int = N_DRAWS, seed: int = SEED, n_envelope: int = N_ENVELOPE,
          box: Mapping[str, Mapping[str, Range]] = BOX, run_kind: str = "prediction") -> Dict[str, Any]:
    """``run_kind`` is ``declared`` for the pre-registered stub-box run of 2026-09-13 and
    ``prediction`` for every run on a box with measured ranges (pre-registration section 7)."""
    programmes = []
    for label, temp_c, time_min in PROGRAMMES:
        result = sweep_programme(label, temp_c, time_min, n_draws=n_draws, seed=seed, n_envelope=n_envelope, box=box)
        result["verdict"] = verdict(result["summary"])
        programmes.append(result)
    labels = [r.label for rows in box.values() for r in rows.values()]
    return {
        "provenance": provenance.provenance_block(
            ARTIFACT, generated_by="src/cultivated_tissue_invariance.py", inputs=[PREREG]
        ),
        "pre_registration": data_paths.rel(PREREG),
        "run_kind": run_kind,
        "box_labels": {lab: labels.count(lab) for lab in LABELS},
        "gap_map": gap_map(box),
        "design": {
            "candidates": list(CANDIDATES),
            "unrankable_declared_candidates": dict(UNRANKABLE),
            "targets": list(TARGETS),
            "metric_targets": list(METRIC_TARGETS),
            "fixed_conditions": dict(FIXED_CONDITIONS),
            "n_draws": n_draws, "seed": seed, "n_envelope": n_envelope,
            "composition_box": box_echo(box),
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
    kind = payload.get("run_kind", "declared")
    counts = payload.get("box_labels", {})
    L.append("# Does the kinetic layer change what to add to cultivated tissue?")
    L.append("")
    if kind == "declared":
        L.append(f"*Generated {payload['provenance']['generated_on']} by `{payload['provenance']['generated_by']}`. "
                 f"Pre-registration: `{payload['pre_registration']}`. THE DECLARED RUN: the composition box is a "
                 "sensitivity device; every range is a stub, none is a measurement.*")
    else:
        L.append(f"*Generated {payload['provenance']['generated_on']} by `{payload['provenance']['generated_by']}`. "
                 f"Pre-registration: `{payload['pre_registration']}`. A PREDICTION RUN on the current box "
                 f"({counts.get('sourced', 0)} sourced, {counts.get('secondary', 0)} secondary, {counts.get('stub', 0)} stub "
                 "ranges). It is not a resolution of the declared run's verdict (pre-registration section 7): the "
                 "same statistics are reported against the same thresholds so the two runs can be read side by side.*")
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
    gaps = payload.get("gap_map") or []
    if gaps:
        L.append("## Gap map: what has not been measured")
        L.append("")
        L.append("Every stub in the box this run swept. A stub is a range with no measurement behind it; the engine's "
                 "ranking depends on the cultivated-side values of the four rankable precursors, so a stub there "
                 "is a gap in the answer, not only in the table.")
        L.append("")
        L.append("| tissue | precursor | engine can rank it | swept range mM | what would close it |")
        L.append("|---|---|---|---|---|")
        for g in gaps:
            L.append(f"| {g['tissue']} | {g['precursor']} | {'yes' if g['engine_rankable'] else 'no'} | "
                     f"{g['lo_mM']} to {g['hi_mM']} | {g['closes_with'] or '—'} |")
        L.append("")
        rankable_gaps = [g for g in gaps if g["engine_rankable"] and g["tissue"] == "cultivated_muscle"]
        if rankable_gaps:
            L.append(f"**One experiment closes the rankable gaps.** {len(rankable_gaps)} of the four precursors the engine "
                     "can rank have no published measurement in cultured muscle. Muroya 2019 quantified cysteine, "
                     "ribose 5-phosphate, IMP and leucine in beef on one CE-TOFMS run; the same panel on washed "
                     "cultured bovine myotubes, with free ribose and glucose by GC-MS or enzymatic assay and thiamine "
                     "by thiochrome HPLC on the same extract, beside a beef sample handled identically, turns every "
                     "cultivated stub into a sourced range in one campaign. Three biological replicates, two harvest "
                     "washes (none; PBS), one ageing arm (24 h at 2 °C) to see whether the IMP-to-ribose route runs "
                     "in a construct at all.")
            L.append("")
    L.append("## The composition box that was swept")
    L.append("")
    L.append("| tissue | precursor | lo mM | hi mM | label | dossiers | note |")
    L.append("|---|---|---|---|---|---|---|")
    for tissue, rows in d["composition_box"].items():
        for c, r in rows.items():
            L.append(f"| {tissue} | {c} | {r['lo_mM']} | {r['hi_mM']} | {r['label']} | "
                     f"{', '.join(r.get('dossiers') or []) or '—'} | {r['note']} |")
    L.append("")
    L.append(f"Fixed in every spec: pH {d['fixed_conditions']['ph']}, a_w {d['fixed_conditions']['aw']}, matrix "
             f"{d['fixed_conditions']['matrix']}, phosphate {d['fixed_conditions']['buffer']['phosphate_mol_l']} mol/L. "
             f"{d['n_draws']} draws per programme, seed {d['seed']}, {d['n_envelope']} envelope draws per reversal.")
    L.append("")
    return "\n".join(L)


def write(payload: Mapping[str, Any], stem: Optional[str] = None) -> Tuple[Any, Any]:
    """Write the JSON, its markdown twin and the box echo. ``stem`` names the files; the
    declared run of 2026-09-13 owns the default stem and is not overwritten by a prediction run."""
    import yaml

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT_JSON if stem is None else OUTPUT_DIR / f"{stem}.json"
    echo_path = BOX_ECHO if stem is None else OUTPUT_DIR / f"{stem}_box.yml"
    echo_path.write_text(
        f"# Echo of the composition box this run swept (src/cultivated_tissue_invariance.py; run_kind = "
        f"{payload.get('run_kind', 'declared')}). Labels: stub = no measurement found; secondary = read from a review "
        "or abstract; sourced = read from the measuring paper's table. This file is generated; edit the module, not this file.\n"
        + yaml.safe_dump(payload["design"]["composition_box"], sort_keys=False),
        encoding="utf-8",
    )
    return artifact_io.write_artifact(payload, json_path, render=render_markdown)


__all__ = [
    "BOX", "CANDIDATES", "CLOSES_WITH", "LABELS", "METRIC_TARGETS", "OUTPUT_JSON", "PROGRAMMES", "STUB_BOX",
    "TARGETS", "THRESHOLDS", "box_echo", "build", "evaluate_draw", "gap_map", "kendall_tau", "naive_ranking",
    "render_markdown", "verdict", "write",
]
