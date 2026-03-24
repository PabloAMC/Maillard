import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.geometry_benchmark import (  # noqa: E402
    build_geometry_benchmark_artifact,
    load_geometry_benchmark_entries,
    render_geometry_benchmark_markdown,
)


def test_geometry_benchmark_loads_ground_state_and_ts_seed_entries():
    entries = load_geometry_benchmark_entries()

    assert len(entries) >= 5
    assert any(entry.benchmark_kind == "ground_state" for entry in entries)
    assert any(entry.benchmark_kind == "ts_seed" for entry in entries)
    assert any(entry.chemistry_family == "sulfur_reference" for entry in entries)


def test_geometry_benchmark_renders_markdown():
    payload = build_geometry_benchmark_artifact()
    markdown = render_geometry_benchmark_markdown(payload)

    assert payload["summary"]["benchmark_id"] == "p4_geometry_benchmark_v1"
    assert "P4 Geometry Benchmark" in markdown
    assert "cysteine_sulfur_probe" in markdown