"""Package CHON-eligible target systems for the React-OT Colab notebook.

This script gathers reactant/product XYZ pairs from
data/geometries/xtb_inputs/<target>/ and builds one uploadable bundle for
Google Colab. The bundle contains:

* reactants.zip  -> <target>-r.xyz entries
* products.zip   -> <target>-p.xyz entries
* manifest.json  -> target list, source paths, and trust posture

The notebook notebooks/react_ot_colab_gpu.ipynb consumes this bundle and runs
React-OT with GPU acceleration as an offline TS-seed generator only.
"""
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "results" / "validation" / "react_ot_colab_bundle.zip"
ELIGIBLE_TARGETS = (
    "lysinoalanine_crosslink",
    "aa_ring_open_dicarbonyl",
    "pe_schiff_base",
    "asparagine_sugar_explicit_water_cluster",
)


def _source_dir(target: str) -> Path:
    return REPO_ROOT / "data" / "geometries" / "xtb_inputs" / target


def _copy_target_pair(target: str, reactants_dir: Path, products_dir: Path) -> dict:
    source_dir = _source_dir(target)
    reactant = source_dir / "reactant.xyz"
    product = source_dir / "product.xyz"
    if not reactant.exists() or not product.exists():
        raise FileNotFoundError(
            f"missing reactant/product xyz for {target} under {source_dir}"
        )

    reactant_name = f"{target}-r.xyz"
    product_name = f"{target}-p.xyz"
    shutil.copyfile(reactant, reactants_dir / reactant_name)
    shutil.copyfile(product, products_dir / product_name)
    return {
        "target": target,
        "reactant": str(reactant.resolve()),
        "product": str(product.resolve()),
        "reactant_upload_name": reactant_name,
        "product_upload_name": product_name,
    }


def _zip_tree(source_dir: Path, destination_zip: Path) -> None:
    with ZipFile(destination_zip, "w", compression=ZIP_DEFLATED) as zip_file:
        for path in sorted(source_dir.rglob("*")):
            if path.is_file():
                zip_file.write(path, arcname=path.relative_to(source_dir.parent))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        action="append",
        choices=ELIGIBLE_TARGETS,
        help="Eligible CHON target to include. Repeat for multiple; defaults to all.",
    )
    parser.add_argument(
        "--out",
        default=str(DEFAULT_OUTPUT),
        help="Output bundle zip for upload to the Colab notebook.",
    )
    args = parser.parse_args()

    targets = args.target or list(ELIGIBLE_TARGETS)
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="react_ot_colab_bundle_") as tmp:
        staging = Path(tmp)
        reactants_dir = staging / "reactants"
        products_dir = staging / "products"
        reactants_dir.mkdir(parents=True, exist_ok=True)
        products_dir.mkdir(parents=True, exist_ok=True)

        entries = [_copy_target_pair(target, reactants_dir, products_dir) for target in targets]
        reactants_zip = staging / "reactants.zip"
        products_zip = staging / "products.zip"
        _zip_tree(reactants_dir, reactants_zip)
        _zip_tree(products_dir, products_zip)

        manifest = {
            "status": "prepared",
            "bundle_version": 1,
            "targets": entries,
            "trust_posture": {
                "role": "ts_initialization_geometry_only",
                "is_runtime_authority": False,
                "energy_use_allowed": False,
            },
        }
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))

        with ZipFile(output_path, "w", compression=ZIP_DEFLATED) as bundle:
            bundle.write(reactants_zip, arcname="reactants.zip")
            bundle.write(products_zip, arcname="products.zip")
            bundle.write(staging / "manifest.json", arcname="manifest.json")

    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())