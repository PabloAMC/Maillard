import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.refinement_watchlist import build_refinement_watchlist, render_refinement_watchlist_markdown


def test_refinement_watchlist_ranks_benchmark_visible_reaction_families_and_materializes_jobs():
    payload = build_refinement_watchlist([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])

    assert payload["summary"]["candidate_count"] > 0
    assert payload["summary"]["run_now"] >= 1
    assert payload["offline_jobs"]
    assert any("thiol" in row["reaction_family"].lower() or "strecker" in row["reaction_family"].lower() for row in payload["candidates"])

    markdown = render_refinement_watchlist_markdown(payload)
    assert "Refinement Watchlist" in markdown
    assert "Reaction Family" in markdown
    assert "Run-now candidates" in markdown
