#!/usr/bin/env python
"""Generate the model card and splice it into README.md between its markers.

    python scripts/generators/generate_model_card.py
    python scripts/generators/generate_model_card.py --check      # CI: fail if stale
    python scripts/generators/generate_model_card.py --stdout     # print, touch nothing

The card is a validity domain -- claim type x system class -> measured reliability and a
trust / caution / do-not-use verdict -- computed from live artifacts by ``src.model_card``.
It replaces the block between ``<!-- BEGIN GENERATED: model-card -->`` and its END marker in
README.md; nothing outside those markers is touched, and no standalone document is created.

``--check`` is the useful mode for a wave: it regenerates and diffs, so "the README's honest
numbers drifted from the artifacts again" becomes a failing command instead of a discovery
three waves later.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.model_card import (  # noqa: E402
    build_model_card,
    render_model_card_markdown,
    splice_into_readme,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--readme", default=str(ROOT / "README.md"))
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "results" / "validation"),
        help="where the machine-readable model_card.json is written",
    )
    parser.add_argument("--stdout", action="store_true", help="print the markdown and exit")
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit non-zero if README.md is not what this generator would write",
    )
    parser.add_argument(
        "--skip-gates",
        action="store_true",
        help="do not run the three CI gates (faster; the gate line is then omitted)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    card = build_model_card(run_gate_checks=not args.skip_gates)
    markdown = render_model_card_markdown(card)

    if args.stdout:
        print(markdown)
        return 0

    readme_path = Path(args.readme)
    updated, changed = splice_into_readme(markdown, readme_path)

    if args.check:
        if changed:
            print(
                "model card is STALE: README.md does not match what the artifacts say.\n"
                "  Regenerate with: python scripts/generators/generate_model_card.py",
                file=sys.stderr,
            )
            return 1
        print("model card: up to date")
        return 0

    readme_path.write_text(updated, encoding="utf-8")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "model_card.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True, default=str), encoding="utf-8")

    print(f"model card {'updated' if changed else 'unchanged'} in {readme_path}")
    print(f"machine-readable payload: {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
