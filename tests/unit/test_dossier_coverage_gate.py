"""scripts/ci/dossier_coverage_gate.py: its two failure modes, on tmp_path fixtures.

(a) a file names a ``*_extraction.md`` that is not in the dossier directory;
(b) an article on disk has no dossier of its stem and no line in the alias table.

Dossier names in this file are built at runtime (``SUFFIX``) so that the gate, which scans
tests/, does not read the fixtures' deliberately missing names as real references.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUFFIX = "_extraction" + ".md"


def _load_gate():
    spec = importlib.util.spec_from_file_location(
        "dossier_coverage_gate", ROOT / "scripts" / "ci" / "dossier_coverage_gate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = _load_gate()


@pytest.fixture()
def dossiers(tmp_path: Path) -> Path:
    d = tmp_path / "extraction_dossiers"
    d.mkdir()
    for stem in ("smith2020", "charlesbernard2005", "renamed2019"):
        (d / f"{stem}{SUFFIX}").write_text("# dossier\n", encoding="utf-8")
    (d / "k1_kinetic_parameters.md").write_text("# synthesis\n", encoding="utf-8")
    return d


def _table(aliases=None, triaged=None, pending=None):
    return {"aliases": dict(aliases or {}), "triaged": dict(triaged or {}), "pending": dict(pending or {})}


# --------------------------------------------------------------------------- (a)


def test_reference_to_a_missing_dossier_fails(tmp_path, dossiers):
    cites = tmp_path / "notes.md"
    cites.write_text(
        f"read in `smith2020{SUFFIX}` and in `ghost2021{SUFFIX}`; see data/lit/extraction_dossiers/Smith2020{SUFFIX}\n",
        encoding="utf-8",
    )
    missing = gate.missing_references([cites], dossiers, tmp_path)
    assert list(missing) == [f"ghost2021{SUFFIX}"]
    assert missing[f"ghost2021{SUFFIX}"] == ["notes.md"]


def test_reference_matching_ignores_case_and_punctuation_and_shorthand(tmp_path, dossiers):
    cites = tmp_path / "code.py"
    cites.write_text(
        f'x = "charles-bernard2005{SUFFIX}"  # hyphenated spelling of an existing dossier\n'
        f'y = "zamora2007/2008/2015{SUFFIX}"  # shorthand, not a name\n'
        f'z = "dossiers/\\nsmith2020{SUFFIX}"  # an escape sequence, then a name\n',
        encoding="utf-8",
    )
    assert gate.missing_references([cites], dossiers, tmp_path) == {}


# --------------------------------------------------------------------------- (b)


def test_article_without_dossier_or_alias_fails(tmp_path, dossiers):
    articles = tmp_path / "articles"
    articles.mkdir()
    for name in ("Smith2020.pdf", "charles-bernard2005.pdf", "orphan2022.pdf", "Renamed.pdf", "patent1.pdf", "queued2026.pdf", ".DS_Store"):
        (articles / name).write_bytes(b"%PDF")
    table = _table(
        aliases={"Renamed": "renamed2019"},
        triaged={"patent1": "presence list only"},
        pending={"queued2026": "downloaded, not yet read"},
    )
    uncovered, absent = gate.uncovered_articles(articles, dossiers, table)
    assert uncovered == ["orphan2022.pdf"]
    assert absent == []
    assert gate.alias_table_problems(table, dossiers) == []


def test_alias_table_defects_fail_without_the_articles(dossiers):
    table = _table(
        aliases={"Renamed": "no_such_dossier", "Synth": "k1_kinetic_parameters"},
        triaged={"patent1": "triaged:  "},
        pending={"Smith2020": "downloaded"},  # a dossier of this name now exists: stale
    )
    problems = gate.alias_table_problems(table, dossiers)
    assert any("no_such_dossier" in p for p in problems)
    assert any("patent1" in p and "no reason" in p for p in problems)
    assert any("Smith2020" in p and "STALE" in p for p in problems)
    assert not any("Synth" in p for p in problems)  # a synthesis dossier is a valid target
    assert len(problems) == 3


def test_the_real_alias_table_is_coherent():
    from src import data_paths

    table = gate.load_alias_table(data_paths.PDF_DOSSIER_ALIASES)
    assert table["aliases"], "the alias table should not be empty"
    assert gate.alias_table_problems(table, data_paths.EXTRACTION_DOSSIERS_DIR) == []
