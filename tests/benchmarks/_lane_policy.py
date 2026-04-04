from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

PHASE33_DATASET = ROOT / "data" / "qm" / "phase33_barrier_benchmarks.json"
PHASE35_DATASET = ROOT / "data" / "qm" / "phase35_double_hybrid_benchmarks.json"
IRC_VALIDATION_FIXTURES = ROOT / "data" / "qm" / "irc_validation_cases"
DFT_REFINER_FILE = ROOT / "src" / "dft_refiner.py"


def _has_real_irc_fixtures(fixtures_root: Path) -> bool:
	if not fixtures_root.exists() or not fixtures_root.is_dir():
		return False
	for path in fixtures_root.rglob("*"):
		if not path.is_file():
			continue
		if path.name.lower() == "readme.md":
			continue
		return True
	return False

_DFT_REFINER_TEXT = DFT_REFINER_FILE.read_text(encoding="utf-8") if DFT_REFINER_FILE.exists() else ""

HAS_PHASE33_DATASET = PHASE33_DATASET.exists()
HAS_PHASE35_DATASET = PHASE35_DATASET.exists()
HAS_IRC_IMPLEMENTATION = "def compute_irc(" in _DFT_REFINER_TEXT
HAS_IRC_FIXTURES = _has_real_irc_fixtures(IRC_VALIDATION_FIXTURES)

PHASE33_SKIP_REASON = f"Phase 3.3 benchmark dataset not mounted at {PHASE33_DATASET}"
PHASE35_SKIP_REASON = f"Phase 3.5 double-hybrid dataset not mounted at {PHASE35_DATASET}"
if not HAS_IRC_IMPLEMENTATION:
	IRC_SKIP_REASON = "IRC authority lane unavailable: compute_irc implementation missing"
elif not HAS_IRC_FIXTURES:
	IRC_SKIP_REASON = f"IRC authority lane unavailable: fixture triplets not mounted at {IRC_VALIDATION_FIXTURES}"
else:
	IRC_SKIP_REASON = "IRC authority lane available"