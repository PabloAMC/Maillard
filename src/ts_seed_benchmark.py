from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src.artifact_io import repo_root


DEFAULT_TS_SEED_BENCHMARK_FILE = repo_root() / "data" / "lit" / "ts_seed_benchmark_set.json"


@dataclass(frozen=True)
class TSSeedBenchmarkEntry:
    benchmark_id: str
    chemistry_family: str
    source_tier: str
    source: str
    reference_xyz: str
    challenged_seed_xyz: str
    recommended_roles: List[str]
    benchmark_value: str = "ts_seed_recovery"


def load_ts_seed_benchmark_payload(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_TS_SEED_BENCHMARK_FILE
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_ts_seed_benchmark_entries(file_path: Optional[Path | str] = None) -> List[TSSeedBenchmarkEntry]:
    payload = load_ts_seed_benchmark_payload(file_path)
    return [TSSeedBenchmarkEntry(**row) for row in payload.get("entries", [])]


def build_ts_seed_benchmark_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    entries = load_ts_seed_benchmark_entries(file_path)
    return {
        "summary": {
            "benchmark_id": "mlp_ts_seed_benchmark_v1",
            "entry_count": len(entries),
            "source_file": str(Path(file_path) if file_path is not None else DEFAULT_TS_SEED_BENCHMARK_FILE),
        },
        "entries": [asdict(entry) for entry in entries],
    }


def render_ts_seed_benchmark_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# TS Seed Recovery Benchmark",
        "",
        "| Benchmark | Chemistry Family | Source Tier | Recommended Roles | Source |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("entries", []):
        lines.append(
            f"| {row.get('benchmark_id', 'unknown')} | {row.get('chemistry_family', 'unknown')} | {row.get('source_tier', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('recommended_roles', [])) or 'none'} | {row.get('source', 'unknown')} |"
        )
    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark id: {summary.get('benchmark_id', 'unknown')}",
            f"Entries: {int(summary.get('entry_count', 0))}",
        ]
    )
    return "\n".join(lines) + "\n"