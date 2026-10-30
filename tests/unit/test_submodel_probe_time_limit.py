"""build_submodel_probe.run_pot_limited (2026-10-30): one stiff pot must not stall a fit's pool."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))
import build_submodel_probe as probe  # noqa: E402


def test_a_pot_over_its_limit_raises_and_a_fast_one_returns(monkeypatch):
    monkeypatch.setattr(probe, "run_pot", lambda job: time.sleep(5) or {"states": []})
    t0 = time.time()
    with pytest.raises(probe.PotTimeout):
        probe.run_pot_limited({"pot": "slow"}, seconds=1)
    assert time.time() - t0 < 3
    monkeypatch.setattr(probe, "run_pot", lambda job: {"states": [1]})
    assert probe.run_pot_limited({"pot": "fast"}, seconds=1) == {"states": [1]}
