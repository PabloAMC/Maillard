#!/usr/bin/env python
"""The front door: one CLI, three verbs, ratios first.

    python scripts/maillard.py compare        two formulations in, per-compound ratios out
    python scripts/maillard.py predict        one formulation, with the model card inline
    python scripts/maillard.py rank-experiments   what to measure next, by value of information

WHY THIS ENTRY POINT EXISTS ALONGSIDE run_pipeline.py
-----------------------------------------------------
2026-08-28 (Wave S5). ``scripts/run_pipeline.py`` is the full forward/inverse designer and is
unchanged. What it does not do is put the model's *measured* strength first: it prints ppb
numbers, and every out-of-sample measurement this repository has made says those numbers are
wrong by 6x to 94x. The directional panel says the ORDINAL claims are right 21 times in 29.

So this verb set leads with ratios, tags every comparison with the panel's measured
reliability on the axis that comparison actually moves, and prints the absolute numbers only
when asked and only next to their caveat. It is orchestration over existing machinery; the
science all lives in src/ and none of it changed to make this work.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.comparative_cli import (  # noqa: E402
    SPEC_TEMPLATE,
    SpecError,
    compare_core,
    compare_systems,
    evaluate_system,
    load_spec_document,
    predict_core,
    predict_system,
    render_compare_core_text,
    render_compare_text,
    render_predict_core_text,
    render_predict_text,
    render_rank_text,
    screening_payload,
    select_system,
    split_comparison_document,
    to_json,
)
from src.config import DEFAULTS  # noqa: E402


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--json", action="store_true", help="emit the machine-readable payload instead of the table")
    parser.add_argument("--top", type=int, default=None, help="show only the N largest rows")
    parser.add_argument("--target-tag", default=DEFAULTS.default_target_tag)
    parser.add_argument("--minimize-tag", default=DEFAULTS.default_minimize_tag)
    parser.add_argument(
        "--lane",
        choices=("core", "fast"),
        default="core",
        help=(
            "core (default): the mass-action kinetic core -- absolute concentrations with "
            "an explicit envelope declaration, and a refusal where it cannot predict. "
            "fast: the ORDINAL SCREENING lane -- rankings only; its absolute ppb are "
            "withheld from every user-facing surface (Wave B5)."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="maillard",
        description=(
            "Comparative interface to the Maillard reaction-network model. Ratios lead because "
            "the model's measured skill is ordinal (directional panel 21/29); absolute "
            "concentrations are unreliable (6x-94x out of sample) and are printed only behind "
            "--absolute, with the caveat attached."
        ),
    )
    verbs = parser.add_subparsers(dest="verb", required=True)

    compare = verbs.add_parser(
        "compare",
        help="two formulations in, per-compound ratios (A/B) out",
        description=(
            "Compare two formulation/process specs. Reports per-compound A/B ratios, the "
            "conservative independent-error bound on each ratio, the dominant (rate-limiting) "
            "pathway, and the directional panel's measured reliability on every axis the two "
            "arms differ along."
        ),
    )
    compare.add_argument("spec", nargs="*", help="one two-arm spec file, or two single-system spec files")
    compare.add_argument("--absolute", action="store_true", help="additionally print absolute ppb, with the caveat")
    compare.add_argument("--template", action="store_true", help="print a starter spec and exit")
    _add_common(compare)

    predict = verbs.add_parser(
        "predict",
        help="one formulation, reported with the model-card caveats inline",
        description=(
            "Score a single formulation. Reports a RANGE, not a point, wherever the "
            "uncertainty panel supports one, with the lane reliability and the dominant "
            "pathway per compound. Every absolute number carries its caveat."
        ),
    )
    predict.add_argument("spec", help="a single-system spec file")
    predict.add_argument("--system", default=None, help="if the file has 'a:'/'b:' arms, which one")
    _add_common(predict)

    rank = verbs.add_parser(
        "rank-experiments",
        help="rank the measurements that would most reduce the model's error",
        description=(
            "Surface the value-of-information loop as a first-class verb. Reads the cached "
            "Monte-Carlo uncertainty panel and ranks candidate wet-lab experiments by how much "
            "each would move the model."
        ),
    )
    rank.add_argument("--top", type=int, default=10)
    rank.add_argument("--matrix", default=None, help="comma-separated matrix families to keep")
    rank.add_argument("--prediction-path", default=None, help="override the uncertainty panel path")
    rank.add_argument("--json", action="store_true")
    rank.add_argument("--write-artifact", action="store_true", help="also write results/validation/experiment_value_ranking.{json,md}")

    return parser


def _load_two_arms(paths, ):
    if len(paths) == 1:
        document = load_spec_document(paths[0])
        return split_comparison_document(document, source=str(paths[0]))
    if len(paths) == 2:
        return (
            select_system(load_spec_document(paths[0]), source=str(paths[0]), arm=None),
            select_system(load_spec_document(paths[1]), source=str(paths[1]), arm=None),
        )
    raise SpecError(
        "compare takes either one two-arm spec file or two single-system spec files. "
        "Run `maillard compare --template` to see the two-arm format."
    )


def run_compare(args: argparse.Namespace) -> int:
    if args.template:
        print(SPEC_TEMPLATE)
        return 0
    if not args.spec:
        raise SpecError("compare needs a spec. Try `maillard compare --template`.")
    spec_a, spec_b = _load_two_arms(args.spec)

    if args.lane == "core":
        if args.absolute:
            print(
                "note: --absolute is redundant on the core lane; the kinetic core reports "
                "absolute concentrations by default, with their envelope declaration.",
                file=sys.stderr,
            )
        payload = compare_core(spec_a, spec_b)
        print(to_json(payload) if args.json else render_compare_core_text(payload))
        return 0

    if args.absolute:
        print(
            "error: --absolute is not available on the screening lane. The FAST lane's\n"
            "  absolute ppb are withheld from every user-facing surface (Wave B5): its\n"
            "  measured skill is ordinal, not quantitative. Re-run with --lane core.",
            file=sys.stderr,
        )
        return 2
    run_a = evaluate_system(spec_a, target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    run_b = evaluate_system(spec_b, target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    payload = screening_payload(compare_systems(run_a, run_b, top_n=args.top))
    print(to_json(payload) if args.json else render_compare_text(payload, show_absolute=False))
    return 0


def run_predict(args: argparse.Namespace) -> int:
    document = load_spec_document(args.spec)
    spec = select_system(document, source=str(args.spec), arm=args.system)

    if args.lane == "core":
        payload = predict_core(spec)
        print(to_json(payload) if args.json else render_predict_core_text(payload))
        return 0

    run = evaluate_system(spec, target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    payload = screening_payload(predict_system(run, top_n=args.top))
    print(to_json(payload) if args.json else render_predict_text(payload))
    return 0


def run_rank(args: argparse.Namespace) -> int:
    from src.experiment_value import (
        PREDICTION_UNCERTAINTY_PATH,
        build_ranking_payload,
        filter_by_matrix,
        load_prediction_payload,
        rank_experiments,
        write_artifact,
    )

    source = Path(args.prediction_path) if args.prediction_path else PREDICTION_UNCERTAINTY_PATH
    if not Path(source).exists():
        print(
            f"error: no uncertainty panel at {source}.\n"
            "  The VoI ranking is computed from the cached Monte-Carlo artifact, which is\n"
            "  gitignored and therefore absent in a fresh clone. Generate it first:\n"
            "    python scripts/generators/generate_prediction_uncertainty.py --n-samples 200 --seed 0",
            file=sys.stderr,
        )
        return 2

    matrix_filter = [item.strip() for item in args.matrix.split(",")] if args.matrix else None
    candidates = rank_experiments(load_prediction_payload(source), top_n=None)
    if matrix_filter:
        candidates = filter_by_matrix(candidates, matrix_filter)
    if args.top is not None:
        candidates = candidates[: max(args.top, 0)]
    payload = build_ranking_payload(candidates, source_path=Path(source), matrix_filter=matrix_filter)
    if args.write_artifact:
        paths = write_artifact(payload)
        print(f"wrote {paths['json']} and {paths['md']}", file=sys.stderr)
    print(to_json(payload) if args.json else render_rank_text(payload))
    return 0


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {"compare": run_compare, "predict": run_predict, "rank-experiments": run_rank}
    try:
        return handlers[args.verb](args)
    except SpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
