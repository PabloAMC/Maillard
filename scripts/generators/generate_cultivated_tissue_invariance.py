#!/usr/bin/env python3
"""Run the cultivated-tissue invariance sweep and write results/cultivated_tissue_invariance/.

The question, the design and the decision rule are in
results/validation/cultivated_tissue_invariance_prereg.md (sections 3 and 6). This script does not
own any threshold; it sweeps the declared composition box through the front door, draws the
parameter envelope on every reversal, and reads the verdict off the pre-registration.

    python scripts/generators/generate_cultivated_tissue_invariance.py                # the declared run
    python scripts/generators/generate_cultivated_tissue_invariance.py --draws 20 --envelope 5   # a smoke run

A smoke run writes to the same paths; do not commit one as the declared run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.cultivated_tissue_invariance import N_DRAWS, N_ENVELOPE, SEED, build, write  # noqa: E402


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--draws", type=int, default=N_DRAWS, help=f"composition draws per programme (declared: {N_DRAWS})")
    parser.add_argument("--envelope", type=int, default=N_ENVELOPE, help=f"parameter draws per reversal (declared: {N_ENVELOPE})")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)
    payload = build(n_draws=args.draws, seed=args.seed, n_envelope=args.envelope)
    json_path, md_path = write(payload)
    for p in payload["programmes"]:
        s, v = p["summary"], p["verdict"]
        print(
            f"{s['programme']}: {v['verdict']} | evaluated {s['evaluated']}/{s['draws']} | refused {s['engine_refused']} | "
            f"top agreement {s['top_agreement_fraction']} | mean tau {s['mean_kendall_tau']} | "
            f"dominant {s['dominant_pair']} | envelope survival {s['mean_envelope_survival']} | {s['wall_seconds']} s"
        )
    print(f"wrote {data_paths.rel(json_path)} and {data_paths.rel(md_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
