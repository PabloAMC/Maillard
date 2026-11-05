"""scripts/tools/openalex.py (plan W1.4): pure helpers and the three subcommands on recorded fixtures (no network)."""
from __future__ import annotations

from scripts.tools import openalex as oa


def _work(wid, doi, title, year, authors, insts=(), retracted=False, refs=()):
    return {"id": f"https://openalex.org/{wid}", "doi": f"https://doi.org/{doi}", "display_name": title,
            "publication_year": year, "is_retracted": retracted, "cited_by_count": 3, "referenced_works": list(refs),
            "authorships": [{"author": {"display_name": a}, "institutions": [{"display_name": i, "ror": None,
                             "country_code": "CH"} for i in insts]} for a in authors],
            "primary_location": {"source": {"display_name": "J. Agric. Food Chem."}},
            "open_access": {"oa_url": None}, "abstract_inverted_index": {"cysteine": [0], "ribose": [1]}}


SEED = _work("W1", "10.1021/jf991254a", "Structured fluids as microreactors", 2000,
             ["Sylvie Vauthey", "Christian Milo"], ["Nestec"], refs=["https://openalex.org/W9"])
CITER = _work("W2", "10.1000/new", "A citing paper", 2024, ["Ana Lopez"], ["U. Tennessee"])
HELD = _work("W3", "10.1000/held", "A held paper", 2010, ["B. Held"])
REF = _work("W9", "10.1000/ref", "An older paper", 1995, ["T. Hofmann"])


def _fake_fetch(table):
    def fetch(url):
        for key, value in table.items():
            if key in url:
                return value
        return None
    return fetch


FETCH = _fake_fetch({
    "works/doi:10.1021/jf991254a": SEED,
    "works/doi:10.1000/retracted": _work("W4", "10.1000/retracted", "Withdrawn", 2011, ["C. Smith"], retracted=True),
    "works/doi:10.1000/wrongyear": _work("W5", "10.1000/wrongyear", "Late", 2019, ["D. Ma"]),
    "works/doi:10.1000/wrongauthor": _work("W6", "10.1000/wrongauthor", "Lopez paper", 2019, ["E. Lopez", "F. Wu"]),
    "filter=cites:W1": {"results": [CITER, HELD], "meta": {"next_cursor": None}},
    "filter=openalex:W9": {"results": [REF], "meta": {}},
})


def test_clean_doi_and_parse_citation():
    assert oa.clean_doi("https://doi.org/10.1021/JF991254A.") == "10.1021/jf991254a"
    assert oa.parse_citation("van Boekel (2001)") == {"surname": "van boekel", "year": 2001}
    assert oa.parse_citation("Ma et al. (2024)") == {"surname": "ma", "year": 2024}
    assert oa.parse_citation("Hofmann & Schieberle 2000")["surname"] == "hofmann"
    assert oa.parse_citation("Münch (1997)")["surname"] == "munch"
    assert oa.parse_citation(None) == {"surname": None, "year": None}


def test_check_entry_flags():
    ok = oa.check_entry({"paper_id": "v", "doi": "10.1021/jf991254a", "citation": "Vauthey et al. (2000)"},
                        FETCH("works/doi:10.1021/jf991254a"))
    assert ok["flags"] == [] and ok["openalex_first_author"] == "Sylvie Vauthey"
    assert oa.check_entry({"doi": "10.1000/missing"}, None)["flags"] == ["doi_not_found"]
    assert "retracted" in oa.check_entry({"doi": "10.1000/retracted", "citation": "Smith (2011)"},
                                         FETCH("works/doi:10.1000/retracted"))["flags"]
    assert "year_mismatch" in oa.check_entry({"doi": "10.1000/wrongyear", "citation": "Ma (2014)"},
                                             FETCH("works/doi:10.1000/wrongyear"))["flags"]
    # the 2026-11-05 error this tool exists to catch: a first author invented from a title-only row
    assert "author_not_in_record" in oa.check_entry({"doi": "10.1000/wrongauthor", "citation": "Zhou et al. (2019)"},
                                                    FETCH("works/doi:10.1000/wrongauthor"))["flags"]


def test_run_verify_counts_and_skips_non_dois():
    registry = {"papers": [{"paper_id": "a", "doi": "10.1021/jf991254a", "citation": "Vauthey (2000)"},
                           {"paper_id": "b", "doi": "10.1000/retracted", "citation": "Smith (2011)"},
                           {"paper_id": "c", "doi": "no-doi-slug", "citation": "X (1990)"}]}
    out = oa.run_verify(registry, FETCH)
    assert out["checked"] == 2 and out["clean"] == 1 and out["flag_counts"] == {"retracted": 1}
    md = oa.render_verify_markdown({**out, "date": "2026-11-05", "registry": "data/keys/papers.yml",
                                    "registry_sha256": "0" * 64})
    assert "retracted" in md and "| b |" in md


def test_run_sweep_forward_backward_and_held():
    known = {"10.1000/held"}
    fwd = oa.run_sweep(["10.1021/jf991254a"], FETCH, known, "forward")
    assert fwd["n_works"] == 2 and fwd["n_held"] == 1
    assert {w["doi"] for w in fwd["works"]} == {"10.1000/new", "10.1000/held"}
    both = oa.run_sweep(["10.1021/jf991254a"], FETCH, known, "both")
    assert {w["doi"] for w in both["works"]} == {"10.1000/new", "10.1000/held", "10.1000/ref"}
    ref = next(w for w in both["works"] if w["doi"] == "10.1000/ref")
    assert ref["via"] == [{"seed": "10.1021/jf991254a", "direction": "backward"}]
    missing = oa.run_sweep(["10.1000/missing"], FETCH, known, "forward")
    assert missing["seeds"] == [{"doi": "10.1000/missing", "status": "not_found"}] and missing["n_works"] == 0


def test_run_labs_reads_institutions():
    rows = oa.run_labs(["10.1021/jf991254a", "10.1000/missing"], FETCH)
    assert rows[0]["authorships"][0]["institutions"][0]["name"] == "Nestec"
    assert rows[1] == {"doi": "10.1000/missing", "found": False, "authorships": []}


def test_sweep_refuses_output_under_data(tmp_path):
    import pytest
    with pytest.raises(SystemExit):
        oa.main(["sweep", "--seed", "10.1021/jf991254a", "--out", str(oa.data_paths.DATA_ROOT / "x.json")])
