from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_GEOMETRY_BENCHMARK_FILE = _repo_root() / "data" / "lit" / "geometry_benchmark_set.json"


@dataclass(frozen=True)
class GeometryBenchmarkEntry:
    benchmark_id: str
    benchmark_kind: str
    chemistry_family: str
    source_tier: str
    source: str
    xyz: str
    recommended_roles: List[str]
    benchmark_value: str = "geometry_sanity"


def load_geometry_benchmark_payload(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_GEOMETRY_BENCHMARK_FILE
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_geometry_benchmark_entries(file_path: Optional[Path | str] = None) -> List[GeometryBenchmarkEntry]:
    payload = load_geometry_benchmark_payload(file_path)
    return [GeometryBenchmarkEntry(**row) for row in payload.get("entries", [])]


def build_geometry_benchmark_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    entries = load_geometry_benchmark_entries(file_path)
    return {
        "summary": {
            "benchmark_id": "mlp_geometry_benchmark_v1",
            "entry_count": len(entries),
            "ground_state_count": sum(1 for entry in entries if entry.benchmark_kind == "ground_state"),
            "ts_seed_count": sum(1 for entry in entries if entry.benchmark_kind == "ts_seed"),
            "source_file": str(Path(file_path) if file_path is not None else DEFAULT_GEOMETRY_BENCHMARK_FILE),
        },
        "entries": [asdict(entry) for entry in entries],
    }


def render_geometry_benchmark_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# MLP Geometry Benchmark",
        "",
        "| Benchmark | Kind | Chemistry Family | Source Tier | Recommended Roles | Source |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("entries", []):
        lines.append(
            f"| {row.get('benchmark_id', 'unknown')} | {row.get('benchmark_kind', 'unknown')} | "
            f"{row.get('chemistry_family', 'unknown')} | {row.get('source_tier', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('recommended_roles', []))} | {row.get('source', 'unknown')} |"
        )
    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark id: {summary.get('benchmark_id', 'unknown')}",
            f"Entries: {int(summary.get('entry_count', 0))}",
            f"Ground-state entries: {int(summary.get('ground_state_count', 0))}",
            f"TS-seed entries: {int(summary.get('ts_seed_count', 0))}",
        ]
    )
    return "\n".join(lines) + "\n"