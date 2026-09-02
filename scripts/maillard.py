#!/usr/bin/env python
"""The front door: one CLI, four verbs, ratios first.

    python scripts/maillard.py compare        two formulations in, per-compound ratios out
    python scripts/maillard.py predict        one formulation, with the model card inline
    python scripts/maillard.py explain        where a compound comes from, and on what evidence
    python scripts/maillard.py rank-experiments   what to measure next, by value of information

``compare`` and ``predict`` also take ``--report PATH.html``, which writes a
self-contained HTML report (src/report_html.py): rankings, an OAV bar chart on a log
axis with measured intervals, refusal cards, the residual decomposition, and every
declared assumption the run used, each with its band. See docs/USING_THE_TOOL.md.

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
    envelope_error_text,
    load_spec_document,
    predict_core,
    render_compare_core_text,
    render_predict_core_text,
    render_rank_text,
    select_system,
    split_comparison_document,
    to_json,
)

# ---------------------------------------------------------------------------
# Help text, written for a bench scientist (Build Wave V1)
# ---------------------------------------------------------------------------
# The previous help assumed the reader already knew what a "lane", a "target
# tag" and an "envelope declaration" were. A flavour chemist opening this for
# the first time does not, and the cost of that is silent misuse: a defaulted
# unit or an un-read reliability tag is a wrong decision, not a wrong printout.
# Every option below therefore says what the input MEANS PHYSICALLY, and every
# verb carries one worked example that runs as written.

_SPEC_FIELDS = """\
WHAT GOES IN A SPEC FILE (YAML or JSON; the two are read by the same loader)

  precursors:    what you WEIGHED IN, as name -> MILLIMOLAR (mM) in the reacting
                 phase. Not grams, not %w/w. The engine knows a fixed vocabulary
                 of names (glucose, fructose, glycine, ribose/xylose, cysteine,
                 thiamine, asparagine, lysine, ...) plus a few LIPID CARRIERS
                 (pea protein isolate, soy protein isolate, methyl linoleate
                 hydroperoxide) whose declared charge is IGNORED because 'mM of
                 an isolate' has no molar basis. Anything else is refused by
                 name, never guessed at.
  temp_C:        the temperature of the hold, in degrees Celsius. Measured
                 product temperature, not oven set point -- the rate constants
                 are Arrhenius in the real temperature and a 10 C error is a
                 factor of ~2-3.
  time_min:      the length of the hold, in minutes, at temp_C. Come-up and
                 cool-down are NOT modelled; if they matter, split them into a
                 multi-segment program in your own script.
  ph:            pH of the mix. Read the reliability tag before you act on it:
                 the TRUNK and ACRYLAMIDE lanes carry NO pH term at all, so pH
                 there is recorded and ignored, and the panel scores pH claims
                 at or below chance.
  aw:            water activity, 0-1. NO lane has an a_w term; nothing in the
                 fit corpus varies it. It is metadata. It is asked for so that
                 it is on the record, not so that it changes an answer.
  protein_type / matrix:
                 the medium. 'free' or 'water' means an aqueous model system
                 and selects the aqueous odour thresholds. A named food matrix
                 (gelatin_3pct, skim_milk, ...) selects that matrix's MEASURED
                 thresholds where they exist and an explicit
                 'no measured threshold' where they do not. Nothing is ever
                 borrowed from another matrix.
  measured_matrix_ratios:  (optional) compound -> your own MEASURED
                 matrix/water ratio. Supplying it turns the HTML report's
                 residual decomposition on: the model shows how much of YOUR
                 measured shift its named terms explain, and how much is
                 unexplained residual.

RELIABILITY TAGS, IN PLAIN WORDS
  trust        >= 80 % agreement on >= 3 independent claims in the directional
               panel. Act on it.
  caution      >= 60 % agreement, OR too few claims to establish either way.
               Use it to choose what to measure, not what to ship.
  do-not-use   < 60 % agreement, or the axis was never measured. Absence of
               evidence is reported as do-not-use ON PURPOSE.
  NOT RESOLVED a ratio that falls INSIDE the measured same-sample dispersion
               band (~4.8x). It is not a small effect; it is no effect this
               method can see.
  REFUSED      the engine cannot name a species or a lane for your request and
               emits NO NUMBER. This is an answer. It names what is missing.
"""


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit the machine-readable payload instead of the table (everything "
             "the renderer had, plus the full envelope declaration)",
    )
    parser.add_argument(
        "--report",
        metavar="PATH.html",
        default=None,
        help="ALSO write a self-contained HTML report to PATH: rankings, an OAV "
             "bar chart on a log axis with measured intervals, refusal cards, the "
             "residual decomposition and every declared assumption the run used. "
             "One file, no network, prints cleanly.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        help="show only the N largest rows (the table only; the report is complete)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="maillard",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Predict aroma-compound formation in a Maillard / lipid-oxidation\n"
            "system, and -- more usefully -- COMPARE two systems.\n"
            "\n"
            "READ THIS FIRST. This model's measured skill is ORDINAL. Out of sample\n"
            "its absolute concentrations are wrong by 6x-94x; its directional and\n"
            "ranking claims score 24/36 on strictly independent claims, and 18/23\n"
            "once pH and water activity are set aside. So: compare two formulations\n"
            "and read the RATIO. Never quote an absolute number as a specification.\n"
            "Treat any pH or moisture direction as unsupported.\n"
            "\n"
            "Four verbs:\n"
            "  compare           two formulations in, per-compound A/B ratios out\n"
            "  predict           one formulation, with intervals and its caveats inline\n"
            "  explain           where a compound comes from in the model, and on what evidence\n"
            "  rank-experiments  which measurement would most reduce the model's error"
        ),
        epilog=_SPEC_FIELDS,
    )
    verbs = parser.add_subparsers(dest="verb", required=True)

    compare = verbs.add_parser(
        "compare",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="two formulations in, per-compound ratios (A/B) out -- START HERE",
        description=(
            "Compare two formulation/process specs and report the per-compound\n"
            "A/B ratio.\n"
            "\n"
            "WHY THE RATIO AND NOT THE NUMBER. The two biggest error sources in an\n"
            "absolute prediction -- the HS-SPME calibration offset and the air/water\n"
            "partition constant -- are SHARED between two arms run through the same\n"
            "model and cancel exactly in a within-run ratio. What is left is the part\n"
            "that was measured to work.\n"
            "\n"
            "HOW TO READ IT. A ratio inside the measured same-sample dispersion band\n"
            "(~4.8x) is reported NOT RESOLVED, not as a small effect: it is a\n"
            "difference this analytical method cannot see. A ratio against a REFUSAL\n"
            "is not a ratio, so if either arm is out of envelope the comparison\n"
            "declines and names what is missing."
        ),
        epilog=(
            "WORKED EXAMPLE -- does swapping ribose for glucose buy you more meaty thiol?\n"
            "\n"
            "  python scripts/maillard.py compare docs/examples/compare_ribose_vs_glucose.yml \\\n"
            "      --report /tmp/ribose_vs_glucose.html\n"
            "\n"
            "Both arms are 10 mM cysteine + 10 mM sugar, 140 C, 30 min, pH 5.0. The only\n"
            "axis that moves is sugar identity, which is the panel's strongest axis "
            "(9/11).\n"
            "You will see a table of A/B ratios per compound, each tagged resolved or NOT\n"
            "RESOLVED, and -- with --report -- an HTML page with the OAV chart, the "
            "refusals\nand every declared assumption the run used.\n"
            + _SPEC_FIELDS
        ),
    )
    compare.add_argument(
        "spec",
        nargs="*",
        help="one two-arm spec file (top-level 'a:' and 'b:'), or two single-system spec files",
    )
    compare.add_argument(
        "--template", action="store_true", help="print a ready-to-edit two-arm spec and exit"
    )
    _add_common(compare)

    predict = verbs.add_parser(
        "predict",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="one formulation, with intervals and the model-card caveats inline",
        description=(
            "Integrate ONE formulation through the kinetic core and report what\n"
            "comes out.\n"
            "\n"
            "WHAT YOU GET. Per compound: an absolute ug/L (= ppb in water) with its\n"
            "MEASURED reliability interval, the odour activity value (concentration\n"
            "/ odour threshold) with its interval, and the threshold's own source and\n"
            "provenance flag. A compound with no measured threshold in your matrix is\n"
            "reported as having none -- never given a borrowed one.\n"
            "\n"
            "WHAT YOU DO NOT GET. A number for anything the core cannot name. HMF,\n"
            "DMHF, 2-pentylfuran and 1-hexanol are refused BY NAME with the reason,\n"
            "because a plausible-looking float with nothing behind it is worse than a\n"
            "gap.\n"
            "\n"
            "For a DECISION, use `compare`. The absolutes here are order-of-magnitude\n"
            "hypotheses that tell you which experiment to run."
        ),
        epilog=(
            "WORKED EXAMPLE -- what does a cysteine/ribose reaction flavour actually smell of?\n"
            "\n"
            "  python scripts/maillard.py predict docs/examples/compare_ribose_vs_glucose.yml \\\n"
            "      --system a --report /tmp/cys_ribose.html\n"
            "\n"
            "10 mM L-cysteine + 10 mM D-ribose, 140 C, 30 min, pH 5.0, aqueous. The core\n"
            "resolves the SULFUR lane and reports MFT, FFT, the MFT dimer, furfural,\n"
            "2-acetylthiazole and methanethiol. Read the OAV column, not the ug/L column:\n"
            "the dimer is ~15.6x more potent than its own monomer, so mass lost to\n"
            "dimerisation is NOT aroma lost.\n"
            + _SPEC_FIELDS
        ),
    )
    predict.add_argument("spec", help="a single-system spec file")
    predict.add_argument(
        "--system",
        default=None,
        help="if the file has 'a:'/'b:' arms, which one to score (a or b)",
    )
    _add_common(predict)

    explain = verbs.add_parser(
        "explain",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="where a compound comes from in this model, and on what evidence",
        description=(
            "Print every route the model has to a compound, the EVIDENCE CLASS of\n"
            "each step's rate constant, and the literature anchors those steps rest\n"
            "on.\n"
            "\n"
            "The four evidence classes:\n"
            "  measured  a rate constant or activation energy printed in a paper\n"
            "  fitted    estimated here by least squares on declared FIT rows\n"
            "  derived   computed from another constant by a stated relation\n"
            "  pinned    held at a value nothing measured -- including held at zero\n"
            "\n"
            "Asking about a compound the model refuses (HMF, DMHF, 2-pentylfuran,\n"
            "1-hexanol, propanal, 2-nonenal) prints the declared reason rather than\n"
            "an error: it tells you what would have to be measured for the answer to\n"
            "exist."
        ),
        epilog=(
            "WORKED EXAMPLE -- why does the model think 2-methyl-3-furanthiol forms at all?\n"
            "\n"
            "  python scripts/maillard.py explain MFT\n"
            "  python scripts/maillard.py explain hexanal      # the lipid lane\n"
            "  python scripts/maillard.py explain HMF          # a refusal, with its reason\n"
            "\n"
            "Nothing on that page is new data: every line is read from a frozen registry\n"
            "at run time, so it moves the day the model moves.\n"
        ),
    )
    explain.add_argument(
        "compound",
        help="a compound name or alias -- MFT, 2-furfurylthiol, hexanal, acrylamide, "
             "furfural, melanoidins, ...",
    )
    explain.add_argument("--json", action="store_true", help="emit the payload instead of the text")

    score = verbs.add_parser(
        "score",
        help="score YOUR measured concentrations against the core (bring your own data)",
        description=(
            "Score one or more of your own measured systems the way the panel scorecard scores a\n"
            "benchmark: fold error, the 3x band, the reliability interval, and a named refusal for\n"
            "any compound the core cannot represent. Writes a bundle-shaped record per system under\n"
            "results/user/ (never under data/) that a future fit wave can read. This is scoring,\n"
            "not calibration: nothing is refitted."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    score.add_argument("spec", nargs="?", help="a measurement document (YAML or JSON); see --template")
    score.add_argument("--template", action="store_true", help="print an example document and exit")
    score.add_argument("--json", action="store_true", help="emit the payload instead of the table")
    score.add_argument("--out", default=None, help="directory for the records (default results/user/)")
    score.add_argument("--no-write", action="store_true", help="score only; write no record")

    rank = verbs.add_parser(
        "rank-experiments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="rank the measurements that would most reduce the model's error",
        description=(
            "Which experiment should you book next?\n"
            "\n"
            "Reads the cached Monte-Carlo uncertainty panel and ranks candidate\n"
            "wet-lab measurements by value of information. Every ranked row is a\n"
            "place the model is MEASURABLY WRONG, converted into a bookable\n"
            "measurement -- so this is the one claim type that does NOT depend on the\n"
            "model being right. It depends on the model being wrong in a located,\n"
            "quantified way, which it demonstrably is."
        ),
        epilog=(
            "WORKED EXAMPLE -- the three measurements that would move the model most:\n"
            "\n"
            "  python scripts/maillard.py rank-experiments --top 3\n"
            "\n"
            "If it exits saying the panel is absent, generate it first:\n"
            "  python scripts/generators/generate_prediction_uncertainty.py --n-samples 200 --seed 0\n"
        ),
    )
    rank.add_argument("--top", type=int, default=10, help="how many candidates to show")
    rank.add_argument(
        "--matrix",
        default=None,
        help="comma-separated matrix families to keep (e.g. pea_isolate,soy_isolate)",
    )
    rank.add_argument(
        "--prediction-path", default=None, help="override the uncertainty panel path"
    )
    rank.add_argument("--json", action="store_true", help="emit the payload instead of the table")
    rank.add_argument(
        "--write-artifact",
        action="store_true",
        help="also write results/validation/experiment_value_ranking.{json,md}",
    )

    return parser


def _write_report(payload, path: str) -> None:
    """Render the HTML report and say where it went. Never silent."""
    from src.report_html import write_report

    target = write_report(payload, path)
    print(f"wrote {target}", file=sys.stderr)


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
    payload = compare_core(spec_a, spec_b)
    print(to_json(payload) if args.json else render_compare_core_text(payload))
    comparison = payload.get("comparison") or {}
    if not comparison.get("comparable"):
        for label, key in (("A", "declaration_a"), ("B", "declaration_b")):
            declaration = comparison.get(key) or {}
            if declaration.get("state") == "out_of_envelope":
                print(f"arm {label}: {envelope_error_text(declaration)}", file=sys.stderr)
    if args.report:
        _write_report(payload, args.report)
    return 0


def run_predict(args: argparse.Namespace) -> int:
    document = load_spec_document(args.spec)
    spec = select_system(document, source=str(args.spec), arm=args.system)
    payload = predict_core(spec)
    print(to_json(payload) if args.json else render_predict_core_text(payload))
    if not payload.get("answered"):
        print(envelope_error_text(payload.get("declaration") or {}), file=sys.stderr)
    if args.report:
        _write_report(payload, args.report)
    return 0


def run_score(args: argparse.Namespace) -> int:
    from src.kinetic_core.user_scoring import (
        TEMPLATE,
        MeasurementSpecError,
        load_document,
        render_text,
        score_document,
        write_records,
    )

    if args.template:
        print(TEMPLATE)
        return 0
    if not args.spec:
        raise SpecError("score needs a measurement document. Try `maillard score --template`.")
    try:
        payload = score_document(load_document(args.spec))
    except MeasurementSpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    written = write_records(payload, args.out) if not args.no_write else []
    print(to_json(payload) if args.json else render_text(payload))
    for path in written:
        print(f"wrote {path}", file=sys.stderr)
    return 0


def run_explain(args: argparse.Namespace) -> int:
    from src.explain_compound import explain, render_explain_text

    payload = explain(args.compound)
    print(to_json(payload) if args.json else render_explain_text(payload))
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
            f"error: no core uncertainty envelope at {source}.\n"
            "  The VoI ranking is computed from the kinetic core's Monte-Carlo envelope.\n"
            "  Generate it first (about 40 min at n=200):\n"
            "    python scripts/generators/generate_core_prediction_uncertainty.py --n-samples 200 --seed 0 --workers 4",
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
    handlers = {
        "compare": run_compare,
        "predict": run_predict,
        "explain": run_explain,
        "score": run_score,
        "rank-experiments": run_rank,
    }
    try:
        return handlers[args.verb](args)
    except SpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
