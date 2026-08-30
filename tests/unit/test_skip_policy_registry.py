import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.skip_policy_registry import (  # noqa: E402
    SKIP_LINE_PATTERN,
    build_skip_policy_registry,
)


def test_skip_line_pattern_parses_pytest_rs_output():
    line = "SKIPPED [1] tests/qm/test_optional_backend.py:17: crest binary not found"
    match = SKIP_LINE_PATTERN.match(line)

    assert match is not None
    assert match.group("count") == "1"
    assert match.group("path").endswith("test_optional_backend.py")


def test_skip_policy_registry_classifies_a_real_skip_into_a_lane_and_owner():
    """The registry's actual job: turn pytest -rs lines into classified, owned rows.

    REPLACED 2026-08-27 (Wave J2, red-team finding: no-op tests). The previous body was:

        payload = build_skip_policy_registry(["tests/benchmarks"])
        assert payload["summary"]["skip_count"] == 0
        assert payload["summary"]["lane_counts"] == {}
        assert payload["rows"] == []

    Every one of those assertions was TRUE BY CONSTRUCTION AND VACUOUS. `tests/benchmarks/`
    contained no tests -- only an orphaned `_lane_policy.py` helper for a QM authority lane
    deleted on 2026-04-21 -- so the scan collected nothing and the test asserted that
    nothing is nothing. It spawned a pytest subprocess to do it. The classification logic
    the module exists for (`_classify_lane`, `_classify_dependency`, `_owner_for_class`,
    `_unblock_criteria`, the per-lane counting) was never touched.

    That directory has now been deleted, which would have made the old test's premise
    disappear anyway. This test drives the parse-and-classify path directly with synthetic
    `pytest -rs` output, so it is fast, deterministic, and independent of how many optional
    QM backends happen to be installed on the machine running it -- a dependency the old
    formulation would have acquired the moment it pointed at a lane with real skips.
    """
    import src.skip_policy_registry as registry

    scan_output = "\n".join(
        [
            "SKIPPED [3] tests/qm/test_backend_plumbing.py:43: xTB binary not found in PATH",
            "SKIPPED [1] tests/qm/test_sella_ts.py:16: Sella not installed",
            "1 passed, 4 skipped in 0.10s",
        ]
    )

    original = registry.run_skip_scan
    registry.run_skip_scan = lambda targets=None: {
        "command": "python -m pytest tests/qm -rs -q",
        "returncode": 0,
        "output": scan_output,
    }
    try:
        payload = build_skip_policy_registry(["tests/qm"])
    finally:
        registry.run_skip_scan = original

    summary = payload["summary"]
    assert summary["skip_count"] == 4, "counts must sum the [N] multiplicity, not the lines"
    assert summary["skip_cluster_count"] == 2
    assert summary["quasi_harmonic_helper_active"] is False

    # The classification actually happened: lanes and classes are populated, not empty.
    assert summary["lane_counts"], "no lane was assigned; classification did not run"
    assert sum(summary["lane_counts"].values()) == 4
    assert summary["dependency_class_counts"]
    assert sum(summary["dependency_class_counts"].values()) == 4

    rows = payload["rows"]
    assert [row["path"] for row in rows] == [
        "tests/qm/test_backend_plumbing.py",
        "tests/qm/test_sella_ts.py",
    ], "rows must be sorted by descending count"

    xtb_row = rows[0]
    assert xtb_row["count"] == 3
    assert xtb_row["line"] == 43
    assert xtb_row["reason"] == "xTB binary not found in PATH"
    # Every row must carry the fields that make the registry actionable rather than
    # decorative: who owns the skip and what would unblock it.
    for row in rows:
        assert row["dependency_class"]
        assert row["lane"]
        assert row["owner"]
        assert row["unblock_criteria"]