from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.ts_seed_benchmark import (  # noqa: E402
    build_ts_seed_benchmark_artifact,
    load_ts_seed_benchmark_entries,
    render_ts_seed_benchmark_markdown,
)


def test_ts_seed_benchmark_loads_entries():
    entries = load_ts_seed_benchmark_entries()

    assert len(entries) >= 2
    assert any(entry.chemistry_family == "proton_transfer_rearrangement" for entry in entries)


def test_ts_seed_benchmark_renders_markdown():
    payload = build_ts_seed_benchmark_artifact()
    markdown = render_ts_seed_benchmark_markdown(payload)

    assert payload["summary"]["benchmark_id"] == "mlp_ts_seed_benchmark_v1"
    assert "TS Seed Recovery Benchmark" in markdown
    assert "amadori_ts_seed_recovery" in markdown