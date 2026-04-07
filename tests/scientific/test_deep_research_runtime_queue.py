import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.deep_research_runtime_queue import build_deep_research_runtime_queue, render_deep_research_runtime_queue_markdown


def test_deep_research_runtime_queue_selects_non_benchmark_runtime_batch_from_backlog():
    payload = build_deep_research_runtime_queue(ROOT)
    rows = {row["citation"]: row for row in payload["selected_candidates"]}
    summary = payload["summary"]
    prepared_next_batch = payload["prepared_next_batch"]
    excluded = {row["citation"]: row for row in payload["excluded_candidates"]}

    assert payload["batch_id"] == "2026-04-05-runtime-first-batch-03"
    assert summary["selected_candidate_count"] == len(rows)
    assert summary["selected_candidate_count"] == 0
    assert summary["process_state_calibration_count"] == 0
    assert summary["landed_runtime_citation_count"] == 6
    assert summary["computational_prior_count"] == 0
    assert summary["safety_reference_count"] == 0
    assert summary["excluded_candidate_count"] == 6

    assert rows == {}

    assert excluded["Ordoudi et al. (2014 / PMC12484514)"]["reason"] == "already_landed_in_runtime_registry"
    assert excluded["Hrncirik & Zeelenberg (2014)"]["reason"] == "already_landed_in_runtime_registry"
    assert excluded["Aliani & Farmer (2005)"]["reason"] == "already_landed_in_runtime_registry"
    assert excluded["Blank, Devaud & Grosch (2003)"]["reason"] == "already_landed_in_runtime_registry"
    assert excluded["Glomb & Monnier (1995)"]["reason"] == "already_landed_in_runtime_registry"
    assert excluded["Hidalgo & Zamora (2004)"]["reason"] == "already_landed_in_runtime_registry"

    assert prepared_next_batch["batch_id"] == ""
    assert prepared_next_batch["summary"]["candidate_count"] == 0

    assert "ACS JAFC 3c05991 / PMC10739987" not in rows

    markdown = render_deep_research_runtime_queue_markdown(payload)
    assert "Deep Research Runtime Queue" in markdown
    assert "Out of scope for this batch: benchmark payload promotion." in markdown
    assert "Hrncirik & Zeelenberg (2014)" in markdown
    assert "Blank, Devaud & Grosch (2003)" in markdown
    assert "no curated candidates selected for this lane" in markdown
    assert "Already landed in runtime registries" in markdown
    assert "Prepared Next Batch" not in markdown
    assert "Glomb & Monnier (1995)" in markdown