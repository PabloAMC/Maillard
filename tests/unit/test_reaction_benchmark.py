import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.reaction_benchmark import (  # noqa: E402
    build_reaction_benchmark_alias_index,
    build_reaction_benchmark_artifact,
    load_reaction_benchmark_entries,
    render_reaction_benchmark_markdown,
)


def test_reaction_benchmark_loads_entries_and_aliases():
    entries = load_reaction_benchmark_entries()

    assert len(entries) >= 6
    alias_index = build_reaction_benchmark_alias_index(entries)
    assert "strecker" in alias_index
    assert alias_index["strecker"].reaction_family == "strecker_degradation"
    assert "enolisation" in alias_index


def test_reaction_benchmark_artifact_renders_markdown():
    payload = build_reaction_benchmark_artifact()
    markdown = render_reaction_benchmark_markdown(payload)

    assert payload["summary"]["benchmark_id"] == "maillard_reaction_benchmark_v1"
    assert "Reaction Benchmark Set" in markdown
    assert "strecker_degradation" in markdown
    assert "Allowed Roles" in markdown