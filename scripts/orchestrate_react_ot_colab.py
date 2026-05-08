"""End-to-end orchestrator for the React-OT Colab loop.

This wraps three smaller scripts so a scientist can drive the full
loop with one command:

* prepare:  build the upload bundle + print the Colab URL + step-by-step
            instructions.
* finish:   import a downloaded artifacts zip + run the seed coverage
            report.
* full:     run prepare then (optionally, when ``--archive`` is passed)
            also run finish in the same invocation.

React-OT is treated as a geometry-only pathfinder; the coverage report
is the artifact gating any downstream promotion into the Sella DFT seed
pipeline (per the observable-first governance rule).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.open_react_ot_colab import colab_url, github_url  # noqa: E402
from scripts.prepare_react_ot_colab_bundle import (  # noqa: E402
    DEFAULT_OUTPUT as DEFAULT_BUNDLE_PATH,
    ELIGIBLE_TARGETS,
)


PYTHON = sys.executable


def _run(cmd: List[str]) -> None:
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def _instructions(bundle_path: Path) -> str:
    lines = [
        "",
        "============================================================",
        " React-OT Colab loop — next steps",
        "============================================================",
        f" 1. Bundle ready: {bundle_path}",
        f" 2. Open notebook (Colab):  {colab_url()}",
        f"    Or browse on GitHub:    {github_url()}",
        " 3. In Colab: Runtime → Change runtime type → GPU.",
        " 4. Run all cells; when prompted upload the bundle from step 1.",
        " 5. Download the emitted react_ot_colab_artifacts.zip.",
        " 6. Back here, run:",
        "      ./scripts/docker_maillard.sh react-ot-orchestrate \\",
        "          --finish --archive PATH/TO/react_ot_colab_artifacts.zip",
        "============================================================",
        "",
    ]
    return "\n".join(lines)


def cmd_prepare(targets: Optional[List[str]], bundle_path: Path) -> Path:
    cmd = [PYTHON, str(REPO_ROOT / "scripts" / "prepare_react_ot_colab_bundle.py"),
           "--out", str(bundle_path)]
    for t in targets or []:
        cmd.extend(["--target", t])
    _run(cmd)
    print(_instructions(bundle_path))
    return bundle_path


def cmd_finish(archive: Path, out_dir: Optional[Path]) -> None:
    importer_cmd = [
        PYTHON,
        str(REPO_ROOT / "scripts" / "import_react_ot_colab_artifacts.py"),
        str(archive),
    ]
    if out_dir is not None:
        importer_cmd.extend(["--out-dir", str(out_dir)])
    _run(importer_cmd)

    coverage_cmd = [
        PYTHON,
        str(REPO_ROOT / "scripts" / "report_react_ot_seed_coverage.py"),
    ]
    if out_dir is not None:
        coverage_cmd.extend(["--seed-dir", str(out_dir)])
    _run(coverage_cmd)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        action="append",
        choices=ELIGIBLE_TARGETS,
        help="Restrict bundle to one or more targets (default: all).",
    )
    parser.add_argument(
        "--bundle",
        default=str(DEFAULT_BUNDLE_PATH),
        help="Output path for the upload bundle.",
    )
    parser.add_argument(
        "--archive",
        default=None,
        help="Path to react_ot_colab_artifacts.zip emitted by Colab.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Destination directory for imported seeds (default: results/computational_gap_refinement).",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only build the bundle and print the Colab instructions.",
    )
    parser.add_argument(
        "--finish",
        action="store_true",
        help="Skip prepare; only import the archive and run the coverage report.",
    )
    args = parser.parse_args()

    bundle_path = Path(args.bundle)
    out_dir = Path(args.out_dir) if args.out_dir else None
    archive = Path(args.archive) if args.archive else None

    if args.finish:
        if archive is None:
            parser.error("--finish requires --archive PATH")
        cmd_finish(archive, out_dir)
        return 0

    cmd_prepare(args.target, bundle_path)
    if args.prepare_only:
        return 0

    if archive is not None:
        cmd_finish(archive, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
