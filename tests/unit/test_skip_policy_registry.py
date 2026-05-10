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


def test_skip_policy_registry_summary_exposes_lane_counts():
    payload = build_skip_policy_registry(["tests/benchmarks"])

    assert payload["summary"]["quasi_harmonic_helper_active"] is False
    assert payload["summary"]["skip_count"] == 0
    assert payload["summary"]["lane_counts"] == {}
    assert payload["rows"] == []