"""Referential integrity of the paper key registry (data/keys/papers.yml)."""
from __future__ import annotations

import subprocess
import sys
from collections import Counter

from src import data_paths, paper_keys
from scripts.generators.build_paper_registry import gather


def test_every_doi_under_data_has_exactly_one_paper_row_and_vice_versa():
    live = {p["doi"] for p in gather()["papers"]}
    registry = {p.doi for p in paper_keys.all_papers()}
    assert live == registry, (
        f"DOIs under data/ not in the registry: {sorted(live - registry)}; "
        f"registry rows citing nothing: {sorted(registry - live)}"
    )


def test_paper_ids_are_unique_and_every_paper_is_cited_somewhere():
    ids = Counter(p.paper_id for p in paper_keys.all_papers())
    assert not [i for i, c in ids.items() if c > 1]
    assert all(p.record_ids for p in paper_keys.all_papers())


def test_doi_normalisation_handles_wiley_and_markdown():
    assert paper_keys.normalise_doi("**10.1111/joss.12567**") == "10.1111/joss.12567"
    assert paper_keys.normalise_doi("https://doi.org/10.3168/JDS.2019-17495 ✔") == "10.3168/jds.2019-17495"
    assert (
        paper_keys.normalise_doi("10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9")
        == "10.1002/1521-3803(20010601)45:3<150::aid-food150>3.0.co;2-9"
    )
    assert paper_keys.normalise_doi("(see 10.1021/jf0200826)") == "10.1021/jf0200826"
    assert paper_keys.normalise_doi("10.0000/example-pea-matrix-package") is None


def test_registry_is_current():
    completed = subprocess.run(
        [sys.executable, "scripts/generators/build_paper_registry.py", "--check"],
        cwd=data_paths.REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout


# --- linkage: which dossier is a paper's dossier ----------------------------------------
#
# Until 2026-09-11 a dossier was linked only when it printed its DOI in a table row spelled
# exactly ``| DOI | ... |`` with no code span around the value. 89 dossiers spell it some
# other way and the code span swallowed the value on 9 more, so 24 papers whose DOI is
# printed inside a dossier were recorded as having none -- six of them cited from the engine
# itself. These cases hold the widened rule in place.

DOSSIER_IDENTITY_SPELLINGS = {
    # DOI in a code span: the backtick used to be read as part of the DOI.
    "10.1016/j.foodchem.2004.04.006": "martins2005_extraction.md",
    # "| DOI / PII | ... |" -- the label carries a second field.
    "10.1016/j.foodchem.2015.06.110": "nguyen2016_extraction.md",
    # "| **DOI** | * **`...`** |" -- bold label, starred and bolded code-span value.
    "10.1016/j.foodchem.2017.07.131": "hamzalioglu2018_extraction.md",
    "10.1021/acs.jafc.1c06163": "gigl2021_extraction.md",
    # No identity table at all: the DOI is a prose line in the header preamble.
    "10.1016/j.foodres.2026.119010": "Xin2026b_extraction.md",
}


def test_a_dossier_that_prints_a_doi_in_its_prose_is_linked_to_that_paper():
    for doi, dossier in DOSSIER_IDENTITY_SPELLINGS.items():
        key = paper_keys.for_doi(doi)
        assert key is not None, f"{doi} is not in the registry at all"
        assert key.dossier == f"data/lit/extraction_dossiers/{dossier}", (
            f"{doi} should link to {dossier}, got {key.dossier}"
        )


def test_a_doi_merely_mentioned_by_a_dossier_is_not_linked_to_it():
    """The linkage is the dossier's identity block, not any DOI anywhere in its text.

    Matching anywhere in the text would reassign 37 papers. zhang2024b_extraction.md names
    a competing paper's DOI in its opening paragraph while its own identity row reads
    10.1016/j.foodres.2024.114149; blank1996_extraction.md quotes the DOI of Blank 1997,
    whose own dossier is blank1997_extraction.md.
    """
    mentioned = paper_keys.for_doi("10.1021/acs.jafc.4c05736")
    assert mentioned is None or mentioned.dossier != "data/lit/extraction_dossiers/zhang2024b_extraction.md"
    assert paper_keys.for_doi("10.1016/j.foodres.2024.114149").dossier == (
        "data/lit/extraction_dossiers/zhang2024b_extraction.md"
    )
    assert paper_keys.for_doi("10.1021/jf960997i").dossier == (
        "data/lit/extraction_dossiers/blank1997_extraction.md"
    )


# --- no registry entry is a truncated DOI ------------------------------------------------


def test_no_registry_entry_is_a_bare_doi_prefix():
    """A DOI wrapped across two source lines used to be registered as its own paper.

    ``10.1016/j.foodchem`` and ``10.1021/acs.jafc`` were rows in data/keys/papers.yml until
    2026-09-11: the extractor read physical lines, so it stopped at the wrap, and the
    trailing "." was then trimmed as punctuation. A registrant prefix is never a DOI, and
    no DOI in this corpus is a strict prefix of another.
    """
    dois = sorted(p.doi for p in paper_keys.all_papers())
    for fragment in ("10.1016/j.foodchem", "10.1021/acs.jafc"):
        assert fragment not in dois, f"{fragment} is a registrant prefix, not a DOI"
    known = set(dois)
    prefixes = [(d, e) for d in dois for e in known if e != d and e.startswith(d)]
    assert not prefixes, f"registry entries that are a prefix of another entry: {prefixes}"
    assert not [d for d in dois if "`" in d or d.endswith((".", "*", ","))], (
        f"registry entries carrying markdown punctuation: {[d for d in dois if '`' in d]}"
    )


def test_normalise_doi_stops_at_a_markdown_code_span():
    assert paper_keys.normalise_doi("| DOI | `10.1016/j.foodchem.2004.04.006` | p. 257 |") == (
        "10.1016/j.foodchem.2004.04.006"
    )
    assert paper_keys.normalise_doi("**`10.1021/acs.jafc.1c06163`**") == "10.1021/acs.jafc.1c06163"


# --- the extractors themselves, on synthetic input ---------------------------------------


def test_code_doi_lines_rejoins_a_doi_split_across_a_source_wrap():
    from scripts.generators.build_paper_registry import code_doi_lines

    # implicitly concatenated string literal, wrapped after the "."
    literal = [
        '    "CITE IT. Yu, Seow, Ong & Zhou 2018 (Food Chem. 268:2, 10.1016/j.foodchem."',
        '    "2018.06.108; yu2018_extraction.md Table 1, step 5) measure the same "',
    ]
    assert paper_keys.normalise_doi(code_doi_lines(literal)[0]) == "10.1016/j.foodchem.2018.06.108"

    # continued comment, wrapped before the "."
    comment = [
        "# `data/articles/Kocadagli2016.pdf` (LONGER stem) = Food Chem 10.1016/j.foodchem",
        "#                                  .2016.05.150 = glucose/wheat flour = NOT this",
    ]
    assert paper_keys.normalise_doi(code_doi_lines(comment)[0]) == "10.1016/j.foodchem.2016.05.150"

    # wrapped immediately after the registrant slash
    slash = [
        '    "Kocadagli2016.pdf (LONGER stem, 613 kB) is Food Chem 10.1016/"',
        '    "j.foodchem.2016.05.150, glucose/wheat flour -- and its text layer is "',
    ]
    assert paper_keys.normalise_doi(code_doi_lines(slash)[0]) == "10.1016/j.foodchem.2016.05.150"

    # a complete DOI that merely ends a line is left alone
    complete = [
        "# `data/articles/Kocada2016.pdf`   (SHORTER stem) = JAFC 10.1021/acs.jafc.6b01862",
        "#                                  = glucose +/- NaCl caramelization = THIS source.",
    ]
    assert paper_keys.normalise_doi(code_doi_lines(complete)[0]) == "10.1021/acs.jafc.6b01862"


def test_dossier_identity_doi_reads_the_identity_block_and_nothing_below_it():
    from scripts.generators.build_paper_registry import dossier_identity_doi

    assert dossier_identity_doi("| DOI / PII | 10.1016/j.foodchem.2015.06.110 / S0308-8146 |") == (
        "10.1016/j.foodchem.2015.06.110"
    )
    assert dossier_identity_doi("| **DOI** | * **`10.1021/acs.jafc.1c06163`** | p. 15334 |") == (
        "10.1021/acs.jafc.1c06163"
    )
    # A dossier whose identity row says no DOI is printed has none, whatever follows.
    assert dossier_identity_doi(
        "| DOI | **NO DOI IS PRINTED IN THE PDF.** |\n"
        "| compare | 10.1021/jf0200826 |\n"
    ) is None
    # A prose DOI counts in the header preamble, not once the body has started.
    assert dossier_identity_doi("# Xin et al. 2026b\n\nDOI 10.1016/j.foodres.2026.119010.\n") == (
        "10.1016/j.foodres.2026.119010"
    )
    assert dossier_identity_doi("# K3 inventory\n\n## body\n\nDOI 10.1021/acs.jafc.2c08360 and should be renamed.\n") is None
