#!/usr/bin/env python3
"""Smoke-test pyGSM end-to-end on one or more computational-gap targets.

For each target this:
  1. Reads reactant + product XYZ from the DFT job manifest.
  2. Runs ``GSMRunner.run_de_gsm`` at xTB / GFN2 + ALPB(water).
  3. Probes the resulting TS guess at xTB to check whether it now
     satisfies the saddle gate (1 imaginary frequency, |f| > 50 cm⁻¹).
  4. Persists per-target audit JSON under ``results/validation/``.

The xTB barrier emitted by pyGSM is **never** consumed downstream — it is
recorded for diagnostics only.  This script answers the question: does
pyGSM produce a TS guess that the hardened saddle gate would accept?
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src._native_io import suppress_native_output
from src.dft_refiner import _classify_ts_frequencies
from src.gsm_backend import GSMRunner
from src.xtb_backend import probe_ts_guess_xtb, run_xtb_hessian_thermo

DEFAULT_TARGETS = [
    "hexanal_radical_quench",
    "lysinoalanine_crosslink",
    "aa_ring_open_dicarbonyl",
]
DEFAULT_MANIFEST = (
    "results/computational_gap_refinement/computational_gap_dft_job_manifest.json"
)


def _load_jobs(manifest: Path) -> Dict[str, Dict[str, Any]]:
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    return {str(job.get("target_id")): job for job in payload.get("jobs", [])}


def _run_one(job: Dict[str, Any], scratch_root: Path, *,
             n_nodes: int, max_iters: int, timeout_s: float) -> Dict[str, Any]:
    target_id = str(job.get("target_id"))
    reactant_xyz = (ROOT / str(job["reactant_path"])).read_text(encoding="utf-8")
    product_xyz = (ROOT / str(job["product_path"])).read_text(encoding="utf-8")
    charge = int(job.get("charge", 0) or 0)
    spin = int(job.get("spin", 0) or 0)

    work_dir = scratch_root / target_id
    runner = GSMRunner(
        work_dir=work_dir,
        charge=charge,
        spin=spin,
        n_nodes=n_nodes,
        max_iters=max_iters,
        timeout_s=timeout_s,
    )
    print(f"[{target_id}] Launching pyGSM DE-GSM "
          f"(nodes={n_nodes}, max_iters={max_iters})...", flush=True)
    result = runner.run_de_gsm(reactant_xyz, product_xyz)
    print(f"[{target_id}]   converged={result.converged} "
          f"reason={result.reason}", flush=True)

    record: Dict[str, Any] = {
        "target_id": target_id,
        "gsm": result.to_audit_dict(),
        "passes_saddle_gate": False,
    }

    if result.converged and result.ts_xyz:
        with suppress_native_output():
            probe = probe_ts_guess_xtb(result.ts_xyz, charge, spin)
            hess = run_xtb_hessian_thermo(
                result.ts_xyz, charge, spin,
                electronic_energy_h=probe.get("energy_eh") or 0.0, temp_k=298.15,
            )
        record["xtb_probe"] = {
            "n_imag": probe.get("n_imag"),
            "fmax_ev_ang": probe.get("fmax_ev_ang"),
            "energy_eh": probe.get("energy_eh"),
        }
        freqs = list(hess[0]) if hess else []
        validation = _classify_ts_frequencies(freqs)
        record["ts_validation"] = validation
        record["passes_saddle_gate"] = bool(validation.get("is_true_saddle", False))
        print(f"[{target_id}]   xTB probe: n_imag={probe.get('n_imag')} "
              f"verdict={validation.get('verdict')} "
              f"passes_gate={record['passes_saddle_gate']}", flush=True)

    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST)
    parser.add_argument("--targets", nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--output", default="results/validation/pygsm_smoke_2026_04_19.json")
    parser.add_argument("--scratch", default="results/validation/pygsm_scratch")
    parser.add_argument("--nodes", type=int, default=11)
    parser.add_argument("--max-iters", type=int, default=40,
                        help="Capped low for smoke; full run uses 80.")
    parser.add_argument("--timeout-s", type=float, default=1800.0)
    args = parser.parse_args()

    manifest_path = ROOT / args.manifest
    jobs = _load_jobs(manifest_path)
    scratch_root = ROOT / args.scratch
    scratch_root.mkdir(parents=True, exist_ok=True)

    records: List[Dict[str, Any]] = []
    for target in args.targets:
        job = jobs.get(target)
        if job is None:
            records.append({"target_id": target, "error": "missing in manifest"})
            continue
        try:
            records.append(_run_one(
                job, scratch_root,
                n_nodes=args.nodes,
                max_iters=args.max_iters,
                timeout_s=args.timeout_s,
            ))
        except Exception as exc:
            records.append({
                "target_id": target,
                "error": f"{type(exc).__name__}: {exc}",
            })

    pass_count = sum(1 for r in records if r.get("passes_saddle_gate"))
    payload = {
        "n_targets": len(records),
        "pass_count": pass_count,
        "verdict": ("GO" if pass_count >= max(1, (2 * len(records)) // 3) else "REVIEW"),
        "records": records,
    }
    out_path = ROOT / args.output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nWrote {out_path}: pass_count={pass_count}/{len(records)} "
          f"verdict={payload['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
