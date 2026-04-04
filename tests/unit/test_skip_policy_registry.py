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
    line = "SKIPPED [1] tests/benchmarks/test_quasi_harmonic_correction.py:17: Implementation pending for Phase 3.2"
    match = SKIP_LINE_PATTERN.match(line)

    assert match is not None
    assert match.group("count") == "1"
    assert match.group("path").endswith("test_quasi_harmonic_correction.py")


def test_skip_policy_registry_summary_exposes_lane_counts():
    payload = build_skip_policy_registry(["tests/benchmarks"])

    assert payload["summary"]["quasi_harmonic_helper_active"] is True

    irc_rows = [row for row in payload["rows"] if row["path"].endswith("test_irc_validation.py")]
    if payload["summary"]["skip_count"] == 0:
        assert payload["summary"]["lane_counts"] == {}
        assert irc_rows == []
    else:
        assert "optional_dft_authority_lane" in payload["summary"]["lane_counts"]
        assert irc_rows
        assert all(row["dependency_class"] == "missing_external_dataset" for row in irc_rows)