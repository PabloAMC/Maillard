"""
The figures' freshness record (2026-09-09, backlog FIG-01).

The three figure builders (`build_reaction_tree.py`, `build_thiol_sink_figures.py`,
`build_story_figures.py`) draw from tracked artifacts, curated data and the engine's own source
files, and nothing compared a PNG with the records behind it, so the figures lagged the
artifacts twice. Each builder now calls :func:`record` when it finishes; the call writes
``results/validation/figure_inputs.json`` with a ``provenance`` block whose ``inputs`` are the
union of every builder's inputs, plus the SHA-256 of every figure written. The artifact
freshness gate (`scripts/ci/artifact_freshness_gate.py`, check 3) already re-hashes every input
of every `results/validation/*.json` that carries a provenance block, so a change to a scorecard,
a ship rule, a species table or a builder script fails the gate until the figures are rebuilt,
and a figure edited or deleted by hand fails it too.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable, Sequence

from src import data_paths, provenance

MANIFEST = data_paths.VALIDATION_DIR / "figure_inputs.json"
FIGURE_DIR = data_paths.REPO_ROOT / "docs" / "assets" / "thiol_sink"


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def record(builder: str, inputs: Iterable[Path | str], figures: Sequence[str]) -> Path:
    """Merge this builder's inputs and figure hashes into the manifest and rewrite it."""
    payload = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else {}
    builders = dict(payload.get("builders") or {})
    fig_paths = [FIGURE_DIR / f for f in figures]
    builders[builder] = {
        "inputs": sorted({data_paths.rel(Path(p)) if Path(p).is_absolute() else str(p) for p in inputs}),
        "figures": {f: _sha(FIGURE_DIR / f) for f in figures if (FIGURE_DIR / f).exists()},
    }
    every_input = sorted({i for b in builders.values() for i in b["inputs"]})
    every_figure = sorted({data_paths.rel(FIGURE_DIR / f) for b in builders.values() for f in b["figures"]})
    out = {
        "artifact": "figure_inputs",
        "what": ("the records behind docs/assets/thiol_sink/*.png: the freshness gate re-hashes every input and every "
                 "figure listed here, so the figures cannot lag the artifacts without failing CI"),
        "builders": builders,
        "provenance": provenance.provenance_block(
            "figure_inputs", generated_by="scripts/generators/figure_manifest.py",
            inputs=[data_paths.REPO_ROOT / i for i in every_input] + [data_paths.REPO_ROOT / f for f in every_figure],
        ),
    }
    MANIFEST.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return MANIFEST


#: The engine source files every tree and map is drawn from.
ENGINE_SOURCES = tuple(
    data_paths.REPO_ROOT / "src" / "kinetic_core" / name
    for name in ("network.py", "sulfur.py", "acrylamide.py", "species.py", "species_sulfur.py", "species_acrylamide.py",
                 "parameters.py", "parameters_sulfur.py", "parameters_acrylamide.py", "parameters_dicarbonyl.py",
                 "parameters_furanic.py", "parameters_pyrazine.py", "parameters_lipid.py", "lipid.py", "matrix_sites.py")
)
