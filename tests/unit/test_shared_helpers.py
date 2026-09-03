"""The shared helpers that replaced per-module copies (2026-09-03, improvement backlog).

provenance, report_format, artifact_io.write_artifact, data_paths.bundle_path, and the two
CLI surfaces that now READ their headline numbers from the artifacts instead of typing them.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import artifact_io, data_paths, provenance
from src.report_format import fmt_number, fmt_rate, md_table


# ---------------------------------------------------------------------------- provenance
def test_git_head_has_the_full_shape_and_dirty_is_a_bool():
    head = provenance.git_head()
    assert set(head) == {"commit", "short", "branch", "dirty"}
    assert isinstance(head["dirty"], bool)
    assert head["short"] == head["commit"][:8]


def test_git_head_outside_a_repository_reports_unknown_rather_than_raising(tmp_path):
    head = provenance.git_head(tmp_path)
    assert head["commit"] == "unknown" and head["branch"] == "unknown" and head["dirty"] is False


def test_sha256_of_a_missing_file_is_none_and_inputs_keep_the_row(tmp_path):
    present = tmp_path / "a.json"
    present.write_text("{}", encoding="utf-8")
    assert len(provenance.sha256_of(present)) == 64
    assert provenance.sha256_of(tmp_path / "missing.json") is None
    rows = provenance.input_sources([present, tmp_path / "missing.json"])
    assert [r["sha256"] is None for r in rows] == [False, True]


def test_provenance_block_carries_date_git_and_inputs(tmp_path):
    block = provenance.provenance_block("x", generated_by="src/x.py", wave="B9", inputs=[tmp_path])
    assert block["artifact"] == "x" and block["generated_by"] == "src/x.py" and block["wave"] == "B9"
    assert len(block["generated_on"]) == 10 and set(block["git"]) == {"commit", "short", "branch", "dirty"}
    assert block["inputs"][0]["sha256"] is None  # a directory is not a file
    assert "generated_on" in provenance.VOLATILE_KEYS and "git" in provenance.VOLATILE_KEYS


# ------------------------------------------------------------------------- report_format
@pytest.mark.parametrize(
    "value, style, expected",
    [(None, "auto", "-"), (True, "auto", "yes"), (False, "auto", "no"), (0.5, "auto", "0.500"),
     (12345.6, "auto", "1.23e+04"), (0.001234, "auto", "0.00123"), (0.0, "auto", "0.000"),
     (2.0, "3g", "2"), (float("nan"), "auto", "-"), (float("inf"), "3g", "-")],
)
def test_fmt_number_rules(value, style, expected):
    assert fmt_number(value, style=style) == expected


def test_fmt_number_dash_is_configurable_and_strings_pass_through():
    assert fmt_number(None, dash="--") == "--"
    assert fmt_number("n/a") == "n/a"


def test_fmt_rate():
    assert fmt_rate(4, 5) == "4/5 (80%)"
    assert fmt_rate(0, 0) == "0/0"


def test_md_table_escapes_pipes_and_has_a_separator_row():
    text = md_table(["a", "b"], [[1, "x|y"]])
    assert text.splitlines() == ["| a | b |", "|---|---|", "| 1 | x\\|y |"]


# ---------------------------------------------------------------------------- artifact_io
def test_write_artifact_writes_json_with_trailing_newline_and_optional_markdown(tmp_path):
    payload = {"b": 1, "a": {"nested": Path("p")}}
    json_path, md_path = artifact_io.write_artifact(payload, tmp_path / "out" / "x.json")
    assert md_path is None and json_path.read_text(encoding="utf-8").endswith("}\n")
    assert json.loads(json_path.read_text())["a"]["nested"] == "p"  # default=str
    json_path, md_path = artifact_io.write_artifact(
        payload, tmp_path / "y.json", render=lambda p: f"# {len(p)} keys\n", sort_keys=True
    )
    assert md_path == tmp_path / "y.md" and md_path.read_text() == "# 2 keys\n"
    assert json_path.read_text().index('"a"') < json_path.read_text().index('"b"')


def test_load_json_mapping_is_strict(tmp_path):
    with pytest.raises(Exception):
        artifact_io.load_json_mapping(tmp_path / "nope.json")


# ------------------------------------------------------------------------------ data_paths
def test_bundle_path_resolves_across_the_three_panel_directories():
    trust = data_paths.bundle_path("thiamine_cys_glucose_120C_Bolton1994")
    holdout = data_paths.bundle_path("mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5")
    matrix = data_paths.bundle_path("external_validation_bi_2020_raw_pea_hexanal")
    assert trust.parent == data_paths.BENCHMARKS_DIR
    assert holdout.parent == data_paths.MAILLARD_PATH_HOLDOUT_DIR
    assert matrix.parent == data_paths.EXTERNAL_VALIDATION_DIR
    with pytest.raises(FileNotFoundError):
        data_paths.bundle_path("no_such_bundle")
    # the trust-loop-only helper is unchanged (frozen generators call it)
    assert data_paths.benchmark_path("x") == data_paths.BENCHMARKS_DIR / "x.json"


# ------------------------------------------------------------- headline numbers are read
def test_the_core_caveat_reads_its_number_from_the_scorecard(monkeypatch, tmp_path):
    from src import comparative_cli

    fake = tmp_path / "scores.json"
    fake.write_text(json.dumps({"summary": {"out_of_sample": {"within_band": 7, "rows": 99}}}))
    monkeypatch.setattr(data_paths, "CORE_PANEL_SCORES", fake)
    assert "7 of 99 out-of-sample rows within 3x" in comparative_cli.core_caveat()
    monkeypatch.setattr(data_paths, "CORE_PANEL_SCORES", tmp_path / "absent.json")
    assert "scorecard not generated" in comparative_cli.core_caveat()


def test_the_cli_description_reads_the_directional_headline(monkeypatch, tmp_path):
    import importlib

    maillard = importlib.import_module("scripts.maillard")
    fake = tmp_path / "directional.json"
    fake.write_text(json.dumps({"summary": {
        "headline": [3, 4],
        "independent": {"excluding_ph_aw": {"agree": 2, "evaluable": 3},
                        "by_category": {"sugar_identity": {"agree": 1, "evaluable": 2}}},
    }}))
    monkeypatch.setattr(data_paths, "CORE_DIRECTIONAL_SCORES", fake)
    text = maillard._headline_lines()
    assert "3/4 on strictly" in text and "2/3 once pH" in text
    assert maillard._axis_note("sugar_identity") == ", scored 1/2 on the directional panel"
    assert maillard._axis_note("moisture_aw") == ""
    monkeypatch.setattr(data_paths, "CORE_DIRECTIONAL_SCORES", tmp_path / "absent.json")
    assert "not generated" in maillard._headline_lines()


# ------------------------------------------------------------ freshness comparison
def test_strip_volatile_removes_clock_and_checkout_keys_recursively():
    payload = {"generated_on": "x", "summary": {"git": {}, "hits": 1, "rows": [{"wall_seconds": 2, "v": 3}]}}
    assert provenance.strip_volatile(payload) == {"summary": {"hits": 1, "rows": [{"v": 3}]}}


def test_payload_differences_ignores_volatile_keys_and_float_noise_but_reports_real_moves():
    tracked = {"generated_on": "2026-01-01", "summary": {"hits": 10, "rate": 0.2}, "rows": [1, 2]}
    live = {"generated_on": "2026-09-03", "summary": {"hits": 10, "rate": 0.2 + 1e-15}, "rows": [1, 2]}
    assert provenance.payload_differences(tracked, live) == []
    live["summary"]["hits"] = 11
    live["rows"] = [1]
    diffs = provenance.payload_differences(tracked, live)
    assert any("summary.hits: 10 != 11" in d for d in diffs) and any("rows: length 2 != 1" in d for d in diffs)
    assert provenance.payload_differences({"a": 1}, {"b": 1}) == ["$.a: missing in live", "$.b: missing in tracked"]


def test_stale_inputs_reports_moved_and_missing_inputs(tmp_path, monkeypatch):
    stable = tmp_path / "stable.json"
    stable.write_text("{}", encoding="utf-8")
    moving = tmp_path / "moving.json"
    moving.write_text("1", encoding="utf-8")
    monkeypatch.setattr(data_paths, "REPO_ROOT", tmp_path)
    payload = {"provenance": {"inputs": [
        {"path": "stable.json", "sha256": provenance.sha256_of(stable)},
        {"path": "moving.json", "sha256": provenance.sha256_of(moving)},
        {"path": "gone.json", "sha256": "abc"},
    ]}}
    assert provenance.stale_inputs(payload) == ["gone.json: recorded abc now None"]
    moving.write_text("2", encoding="utf-8")
    assert [d.split(":")[0] for d in provenance.stale_inputs(payload)] == ["moving.json", "gone.json"]
    assert provenance.stale_inputs({}) == []
