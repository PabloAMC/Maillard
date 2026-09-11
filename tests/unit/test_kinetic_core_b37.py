"""
WAVE B37 (2026-09-11): sixteen sources arrive.

Pre-registration: results/validation/kinetic_core_b37_prereg.md.

This wave fits nothing. It reads sixteen papers, writes their dossiers, corrects one
provenance note that was true when written and false now, and makes the 3-deoxy dicarbonyl
series CHARGEABLE so that a fed-intermediate pot can be expressed at all. The tests below
hold the two things that could go wrong: the new aliases must be inert, and each dossier
must name the file it was read from.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import pytest

from src.kinetic_core import engine

ROOT = Path(__file__).resolve().parents[2]

#: The sixteen papers the owner downloaded on 2026-09-11 against the B36 reading list,
#: each with the PDF its dossier must name.
B37_SOURCES = {
    "mittelmaier2010": "data/articles/mittelmaier2010.pdf",
    "frankel1993": "data/articles/frankel1993.pdf",
    "tazi2009": "data/articles/tazi2009.pdf",
    "yu2012": "data/articles/yu2012.pdf",
    "chen2010": "data/articles/chen2010.pdf",
    "ebert2021": "data/articles/ebert2021.pdf",
    "zhou2026": "data/articles/Zhou2026.pdf",
    "xu2024": "data/articles/Xu2024.pdf",
    "kong2024": "data/articles/Kong2024.pdf",
    "liu2025": "data/articles/Liu2025.pdf",
    "cai2021": "data/articles/cai2021.pdf",
    "zhang2012": "data/articles/zhang2012.pdf",
    "moisio2015": "data/articles/moisio2015.pdf",
    "wang2026b": "data/articles/Wang2026b.pdf",
    "zhai2023c": "data/articles/Zhai2023c.pdf",
    "hernandez2023": "data/articles/hernandez2023.pdf",
}

#: The three species B37 made chargeable. B13 made glyoxal, glucosone, diacetyl and
#: methylglyoxal chargeable and left these as targets only.
B37_NEW_PRECURSORS = {
    "3-deoxyglucosone": "TDG",
    "3-dg": "TDG",
    "3,4-dideoxyglucosone": "DDG",
    "3,4-dideoxyglucosone-3-ene": "DDG",
    "3,4-dge": "DDG",
    "1-deoxyglucosone": "ODG",
    "1-dg": "ODG",
}


@pytest.mark.parametrize("stem", sorted(B37_SOURCES), ids=sorted(B37_SOURCES))
def test_every_b37_source_has_a_dossier_naming_the_pdf_it_was_read_from(stem: str):
    path = ROOT / "data" / "lit" / "extraction_dossiers" / f"{stem}_extraction.md"
    assert path.exists(), f"{stem}: no dossier"
    text = path.read_text(encoding="utf-8")
    assert B37_SOURCES[stem] in text, f"{stem}: dossier does not name its PDF"
    assert (ROOT / B37_SOURCES[stem]).exists(), f"{stem}: the PDF it names is not on disk"


def test_the_three_deoxy_dicarbonyls_are_chargeable():
    """
    THE ONE CODE CHANGE IN B37. Mittelmaier et al. 2011 charge pure 3-deoxyglucosone at
    120 C and follow 3,4-dideoxyglucosone -- the experiment docs/guides/EXPERIMENTS.md asks
    for by name -- and no FormulationSpec could express it before this wave.
    """
    for name, key in B37_NEW_PRECURSORS.items():
        assert engine.PRECURSOR_ALIASES.get(name) == key, name


def test_the_new_aliases_are_inert():
    """
    AND IT MUST CHANGE NOTHING. No benchmark bundle, no directional claim and no fit row
    charges any of the three, so every frozen prediction is bit-for-bit what it was.
    """
    charged = set()
    for pattern in ("data/benchmarks/**/*.json", "data/lit/**/*.json", "docs/validation/*.yml"):
        for path in glob.glob(str(ROOT / pattern), recursive=True):
            if "quarantined" in path:
                continue
            raw = Path(path).read_text(encoding="utf-8", errors="replace").lower()
            for name in B37_NEW_PRECURSORS:
                # a bare mention as a TARGET is fine; a "precursors" block naming it is not
                if f'"{name}"' not in raw and f"'{name}'" not in raw:
                    continue
                try:
                    payload = json.loads(Path(path).read_text(encoding="utf-8"))
                except Exception:
                    continue
                for block in _precursor_blocks(payload):
                    if any(str(k).lower() == name for k in block):
                        charged.add(f"{Path(path).name}:{name}")
    assert not charged, f"B37's new precursor aliases are no longer inert: {sorted(charged)}"


def _precursor_blocks(payload):
    if isinstance(payload, dict):
        block = payload.get("precursors")
        if isinstance(block, dict):
            yield block
        for value in payload.values():
            yield from _precursor_blocks(value)
    elif isinstance(payload, list):
        for value in payload:
            yield from _precursor_blocks(value)


def test_the_hernandez_note_is_corrected_and_keeps_its_prior_claim():
    """
    The B34/B35/B36 practice: the new reading first, the false claim after it, labelled.
    """
    path = ROOT / "data" / "benchmarks" / "resconi_2023_pbma_beef_identity_benchmark.json"
    note = json.loads(path.read_text())["conditions"]["vessel"]["provenance_note"]
    assert note.startswith("WAVE B37")
    assert "THE SOURCE IS ON DISK" in note
    assert "data/articles/hernandez2023.pdf" in note
    assert "SUPERSEDED 2026-09-11 AND RETAINED AS THE AUDIT RECORD" in note
    assert note.index("THE SOURCE IS ON DISK") < note.index("SUPERSEDED")
    # the print gives a skillet cook, and the bundle's conditions are NOT that cook
    assert "NO COOK TIME" in note
