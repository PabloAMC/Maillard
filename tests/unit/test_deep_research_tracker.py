from pathlib import Path

from scripts.deep_research_tracker import (
    build_audit_payload,
    collapse_occurrences,
    load_registry_entries,
    parse_markdowns,
    summarize_priority_targets,
)


def test_parse_markdowns_filters_to_numbered_reports_and_min_score(tmp_path: Path):
    (tmp_path / "00_valid.md").write_text(
        "1. `Paper A (2024)` — useful dataset. 6/8.\n"
        "2. `Paper B (2024)` — lower score. 5/8.\n",
        encoding="utf-8",
    )
    (tmp_path / "notes.md").write_text(
        "1. `Paper C (2024)` — should be ignored by filename filter. 8/8.\n",
        encoding="utf-8",
    )

    rows = parse_markdowns(tmp_path, min_score=6)

    assert [row["citation"] for row in rows] == ["Paper A (2024)"]


def test_collapse_occurrences_marks_registry_only_and_runtime_bound(tmp_path: Path):
    artifact_path = tmp_path / "data" / "benchmarks" / "example.json"
    artifact_path.parent.mkdir(parents=True)
    artifact_path.write_text("{}", encoding="utf-8")

    registry_path = tmp_path / "data" / "lit" / "benchmark_intake_registry.json"
    registry_path.parent.mkdir(parents=True)
    registry_path.write_text(
        """
        {
          "eligible_references": [
            {
              "id": "bound",
                                    "citation": "Alpha Study (2024)",
              "runtime_artifacts": [{"path": "data/benchmarks/example.json"}]
            },
            {
              "id": "registry-only",
                                    "citation": "Beta Study (2024)",
              "runtime_artifacts": []
            }
          ]
        }
        """,
        encoding="utf-8",
    )

    registry_entries = load_registry_entries(registry_path, root=tmp_path)
    rows = collapse_occurrences(
        [
            {
                "citation": "Alpha Study (2024)",
                "normalized_citation": "alpha study 2024",
                "author_year_key": "alpha 2024",
                "score": "6/8",
                "score_value": 6,
                "file": "00_valid.md",
                "line": 1,
                "description": "bound",
                "raw_line": "raw",
            },
            {
                "citation": "Beta Study (2024)",
                "normalized_citation": "beta study 2024",
                "author_year_key": "beta 2024",
                "score": "7/8",
                "score_value": 7,
                "file": "00_valid.md",
                "line": 2,
                "description": "registry only",
                "raw_line": "raw",
            },
        ],
        registry_entries,
    )

    by_citation = {row["citation"]: row for row in rows}
    assert by_citation["Alpha Study (2024)"]["status"] == "RUNTIME_BOUND"
    assert by_citation["Beta Study (2024)"]["status"] == "REGISTRY_ONLY"


def test_priority_summary_flags_53_targets_as_runtime_bound_in_current_repo():
    payload = build_audit_payload(min_score=6)
    target_summary = summarize_priority_targets(payload["items"])

    assert target_summary["spi_hvp_xylose"]["status"] == "RUNTIME_BOUND"
    assert target_summary["wheat_gluten_hvp_xylose"]["status"] == "RUNTIME_BOUND"
    assert target_summary["acrylamide_fast_kinetics"]["status"] == "RUNTIME_BOUND"
    assert target_summary["cml_cel_commercial_pba_ranges"]["status"] == "RUNTIME_BOUND"
    assert target_summary["furosine_crossover"]["status"] == "RUNTIME_BOUND"