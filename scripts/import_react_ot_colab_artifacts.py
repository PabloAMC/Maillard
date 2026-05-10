"""Import React-OT Colab artifacts into the local repo workspace.

Consumes the `react_ot_colab_artifacts.zip` archive downloaded from
notebooks/react_ot_colab_gpu.ipynb and copies only the repo-compatible output
files into results/computational_gap_refinement/.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "results" / "computational_gap_refinement"
RECEIPT_NAME = "react_ot_colab_import_receipt.json"


def _allowed_member(name: str) -> bool:
    base = Path(name).name
    if base == "react_ot_pilot_manifest.json":
        return True
    return base.endswith("_react_ot_seed.xyz") or base.endswith("_react_ot_seed.json")


def import_archive(archive_path: Path, output_dir: Path) -> dict:
    if not archive_path.exists():
        raise FileNotFoundError(f"missing archive: {archive_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    imported_files: list[str] = []
    ignored_files: list[str] = []

    with ZipFile(archive_path) as archive:
        for member in archive.namelist():
            if member.endswith("/"):
                continue
            base = Path(member).name
            if not _allowed_member(member):
                ignored_files.append(member)
                continue
            destination = output_dir / base
            destination.write_bytes(archive.read(member))
            imported_files.append(str(destination.relative_to(ROOT)))

    if not any(path.endswith("react_ot_pilot_manifest.json") for path in imported_files):
        raise ValueError("archive missing react_ot_pilot_manifest.json")

    receipt = {
        "status": "imported",
        "source_archive": str(archive_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "imported_files": sorted(imported_files),
        "ignored_files": sorted(ignored_files),
    }
    (output_dir / RECEIPT_NAME).write_text(json.dumps(receipt, indent=2, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", help="Path to react_ot_colab_artifacts.zip")
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Destination directory for imported React-OT seed artifacts.",
    )
    args = parser.parse_args()

    receipt = import_archive(Path(args.archive), Path(args.out_dir))
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())