import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_intake_registry import build_literature_backlog_artifact, render_literature_backlog_markdown


def test_literature_backlog_queues_are_exclusive_and_surface_minimum_primary_experiment():
    payload = build_literature_backlog_artifact(ROOT)
    summary = payload["summary"]
    encoded_rows = {row["id"]: row for row in payload["encoded_reference_rows"]}

    assert summary["encoded_reference_count"] >= 5
    assert summary["queue_conflict_count"] == 0
    assert payload["conflicts"]["encoded_and_ready_ids"] == []

    queue_ids = {row["id"] for row in payload["ready_runtime"] + payload["ready_benchmark"]}
    encoded_ids = set(encoded_rows)
    assert queue_ids.isdisjoint(encoded_ids)
    assert encoded_rows["rizzello_2024_fermentation_cleanup"]["triage_status"] == "ready_calibration"
    assert encoded_rows["zhao_2022_moromi_precursor_release"]["triage_status"] == "ready_calibration"
    assert encoded_rows["ordoudi_2014_hmf_peak_window"]["triage_status"] == "ready_calibration"

    wet_lab_gaps = {row["gap_id"] for row in payload["wet_lab_blocked"]}
    assert "ppi_meaty_positive_matrix_benchmark" in wet_lab_gaps
    assert "spi_meaty_positive_matrix_benchmark" in wet_lab_gaps
    assert "meaty_off_flavour_safety_tradeoff_panel" in wet_lab_gaps

    minimum_primary_experiment = payload["minimum_primary_experiment"]
    assert minimum_primary_experiment["matrices"] == ["PPI 5% buffered slurry", "SPI 5% buffered slurry"]
    assert "HS-SPME-GC-MS" in minimum_primary_experiment["instrumentation"]

    markdown = render_literature_backlog_markdown(payload)
    assert "Literature Backlog" in markdown
    assert "Ready Runtime" in markdown
    assert "Ready Benchmark" in markdown
    assert "Wet-Lab Blocked" in markdown
    assert "Minimum Primary Experiment" in markdown