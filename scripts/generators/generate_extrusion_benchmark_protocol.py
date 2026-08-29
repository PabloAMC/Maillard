#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.doe_generator import (  # noqa: E402
    export_extrusion_benchmark_protocol,
    export_extrusion_sop_lock_register,
    render_extrusion_benchmark_protocol_markdown,
    render_extrusion_sop_lock_register_markdown,
)
from src.extrusion_benchmark_landing import (  # noqa: E402
    export_extrusion_disulfide_follow_on_package,
    export_extrusion_disulfide_follow_on_workbook,
    export_extrusion_external_closure_package,
    export_extrusion_external_closure_workbook,
    render_extrusion_disulfide_follow_on_markdown,
    render_extrusion_disulfide_follow_on_workbook_markdown,
    render_extrusion_external_closure_package_markdown,
    render_extrusion_external_closure_workbook_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Export the extrusion benchmark protocol, its SOP lock register, "
            "the external-closure package and workbook and the disulfide "
            "follow-on package and workbook into results/validation/."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "results" / "validation"),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    payload = export_extrusion_benchmark_protocol(str(output_dir), root=ROOT)
    lock_payload = export_extrusion_sop_lock_register(str(output_dir), root=ROOT)
    closure_payload = export_extrusion_external_closure_package(str(output_dir), root=ROOT)
    closure_workbook = export_extrusion_external_closure_workbook(str(output_dir), root=ROOT)
    follow_on_payload = export_extrusion_disulfide_follow_on_package(str(output_dir), root=ROOT)
    follow_on_workbook = export_extrusion_disulfide_follow_on_workbook(str(output_dir), root=ROOT)
    print(render_extrusion_benchmark_protocol_markdown(payload))
    print(render_extrusion_sop_lock_register_markdown(lock_payload))
    print(render_extrusion_external_closure_package_markdown(closure_payload))
    print(render_extrusion_external_closure_workbook_markdown(closure_workbook))
    print(render_extrusion_disulfide_follow_on_markdown(follow_on_payload))
    print(render_extrusion_disulfide_follow_on_workbook_markdown(follow_on_workbook))
    print(f"Wrote {output_dir / 'extrusion_benchmark_protocol.md'}")
    print(f"Wrote {output_dir / 'extrusion_benchmark_protocol.json'}")
    print(f"Wrote {output_dir / 'extrusion_sop_lock_register.md'}")
    print(f"Wrote {output_dir / 'extrusion_sop_lock_register.json'}")
    print(f"Wrote {output_dir / 'extrusion_external_closure_package.md'}")
    print(f"Wrote {output_dir / 'extrusion_external_closure_package.json'}")
    print(f"Wrote {output_dir / 'extrusion_external_closure_workbook.md'}")
    print(f"Wrote {output_dir / 'extrusion_external_closure_workbook.yaml'}")
    print(f"Wrote {output_dir / 'extrusion_disulfide_follow_on_package.md'}")
    print(f"Wrote {output_dir / 'extrusion_disulfide_follow_on_package.json'}")
    print(f"Wrote {output_dir / 'extrusion_disulfide_follow_on_workbook.md'}")
    print(f"Wrote {output_dir / 'extrusion_disulfide_follow_on_workbook.yaml'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())