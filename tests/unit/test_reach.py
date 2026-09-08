"""Reach (roadmap programme 4): the console entry point, the Python API, the spec schema and the page."""
from __future__ import annotations

import json
import subprocess
import sys
import threading
import urllib.request
from pathlib import Path

import pytest
import yaml

from src import api
from src.comparative_cli import SPEC_TEMPLATE, SpecError, validate_spec

ROOT = Path(__file__).resolve().parents[2]


def test_the_cli_module_is_the_front_door_and_the_script_is_a_shim():
    from src import cli

    assert cli.main(["compare", "--template"]) == 0
    shim = (ROOT / "scripts" / "maillard.py").read_text(encoding="utf-8")
    assert "from src.cli import main" in shim and len(shim.splitlines()) < 30


def test_the_console_script_is_declared():
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'maillard = "src.cli:main"' in text


def test_running_the_module_lists_every_verb():
    out = subprocess.run([sys.executable, "-m", "src.cli", "--help"], cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0
    for verb in ("compare", "predict", "explain", "score", "calibrate", "wishlist", "ui"):
        assert verb in out.stdout, verb


def test_the_api_returns_the_verbs_payloads():
    doc = api.template()
    cmp = api.compare(doc["a"], doc["b"])
    assert cmp["artifact"] == "maillard_compare_core"
    pred = api.predict(doc["a"])
    assert pred["artifact"] == "maillard_predict_core" and pred["answered"]
    exp = api.explain("2-methyl-3-furanthiol")
    assert exp["answered"] and exp["routes"]


def test_the_schema_refuses_what_the_hand_validator_used_to_pass_silently():
    doc = api.template()
    bad = dict(doc["a"])
    bad["ph"] = 19
    with pytest.raises(SpecError) as info:
        validate_spec(bad, label="a")
    assert "ph" in str(info.value) and "spec.schema.json" in str(info.value)
    bad = dict(doc["a"])
    bad["precursors"] = {"D-Ribose": "ten"}
    with pytest.raises(SpecError):
        validate_spec(bad, label="a")


def test_the_page_round_trips_the_template():
    from src import ui

    server = ui.make_server(port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        page = urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=30).read().decode("utf-8")
        assert "<textarea" in page and "cysteine_ribose" in page
        body = json.dumps({"verb": "compare", "spec": SPEC_TEMPLATE, "calibration": None}).encode("utf-8")
        req = urllib.request.Request(f"http://127.0.0.1:{port}/run", data=body, headers={"Content-Type": "application/json"})
        html = urllib.request.urlopen(req, timeout=120).read().decode("utf-8")
        assert "<html" in html.lower() and "cysteine_ribose" in html
        bad = json.dumps({"verb": "compare", "spec": "a: 1", "calibration": None}).encode("utf-8")
        req = urllib.request.Request(f"http://127.0.0.1:{port}/run", data=bad, headers={"Content-Type": "application/json"})
        with pytest.raises(urllib.error.HTTPError) as info:
            urllib.request.urlopen(req, timeout=30)
        assert info.value.code == 400
    finally:
        server.shutdown()
        server.server_close()
