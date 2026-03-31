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

    assert len(entries) >= 11
    alias_index = build_reaction_benchmark_alias_index(entries)
    assert "strecker" in alias_index
    assert alias_index["strecker"].reaction_family == "strecker_degradation"
    assert "enolisation" in alias_index
    assert "trapping_hexanal" in alias_index
    assert alias_index["trapping_hexanal"].reaction_family == "lipid_strecker_synergy"
    assert "dha_elimination" in alias_index
    assert alias_index["dha_elimination"].reaction_family == "beta_elimination"
    assert "lysinoalanine" in alias_index
    assert alias_index["lysinoalanine"].reaction_family == "beta_elimination"
    assert "thiamine_fragmentation" in alias_index
    assert alias_index["thiamine_fragmentation"].reaction_family == "thiamine_degradation"
    assert "vitamin_b1" in alias_index
    assert alias_index["vitamin_b1"].reaction_family == "thiamine_degradation"
    assert "glutathione_release" in alias_index
    assert alias_index["glutathione_release"].reaction_family == "glutathione_cleavage"
    assert "gsh_release" in alias_index
    assert alias_index["gsh_release"].reaction_family == "glutathione_cleavage"


def test_reaction_benchmark_artifact_renders_markdown():
    payload = build_reaction_benchmark_artifact()
    markdown = render_reaction_benchmark_markdown(payload)

    assert payload["summary"]["benchmark_id"] == "maillard_reaction_benchmark_v1"
    assert "Reaction Benchmark Set" in markdown
    assert "strecker_degradation" in markdown
    assert "lipid_strecker_synergy" in markdown
    assert "beta_elimination" in markdown
    assert "thiamine_degradation" in markdown
    assert "glutathione_cleavage" in markdown
    assert "Allowed Roles" in markdown