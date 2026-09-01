"""Referential integrity of the compound key registry (data/keys/compounds.yml).

Phase 3 of the data restructure (2026-09-01). Three things are pinned:

1. Every compound spelling that occurs anywhere under data/ resolves to exactly one
   registry entry -- benchmarks, species YAML, flavor/safety payloads, Henry table,
   decision panel. A new spelling that resolves to nothing fails here, which is the
   point: names no longer float free.
2. The registry is internally consistent: no spelling claimed twice, no two molecules
   with the same InChIKey, every class member points at an existing entry.
3. The registry covers everything the retired hand-written BENCHMARK_NAME_ALIASES table
   covered, so replacing that table with `compound_keys.match_norms` cannot have lost a
   match on the panel.
"""
from __future__ import annotations

import glob
import json
from collections import Counter

from src import compound_keys, data_access, data_paths
from src.text_utils import normalize_compound_name as norm

# The literal that lived in src/benchmark_validation.py until 2026-09-01, verbatim.
RETIRED_BENCHMARK_NAME_ALIASES = {
    "2methyl3furanthiol": {"2methyl3furanthiolmft", "2methyl3furylthiol", "2methylfuran3thiol", "mft"},
    "2furfurylthiol": {"2furfurylthiolfft", "2furylmethanethiol", "fft"},
    "bis2methyl3furyldisulfide": {"bis2methyl3furyldisulfide", "2methyl3furyl2methyl3furyldisulfide"},
    "pyrazine": {
        "25dimethylpyrazine",
        "23dimethylpyrazine",
        "26dimethylpyrazine",
        "2ethyl35dimethylpyrazine",
        "2ethylpyrazine",
        "trimethylpyrazine",
        "tetramethylpyrazine",
        "dimethylpyrazine",
    },
    "3methylbutanal": {"3methylbutanal", "isovaleraldehyde"},
    "2methylbutanal": {"2methylbutanal", "2methylbutyraldehyde"},
    "methional": {"3methylthiopropanal", "methional"},
    "acrylamide": {"acrylamide", "2-propenamide"},
}


def _repository_spellings() -> dict[str, set[str]]:
    seen: dict[str, set[str]] = {}

    def add(name, where):
        if isinstance(name, str) and name.strip():
            seen.setdefault(name.strip(), set()).add(where)

    for f in glob.glob(str(data_paths.BENCHMARKS_DIR / "**" / "*.json"), recursive=True):
        d = json.load(open(f, encoding="utf-8"))
        for key in ("measured_volatiles", "reference_volatiles", "holdout_targets"):
            for name in (d.get(key) or {}):
                add(name, data_paths.rel(f))
    for p in (data_paths.DESIRABLE_TARGETS, data_paths.OFF_FLAVOUR_TARGETS, data_paths.TOXIC_MARKERS):
        for c in (data_access.load_yaml(p) or {}).get("compounds", []):
            add(c.get("name"), data_paths.rel(p))
    for c in (data_access.load_yaml(data_paths.HENRY_CONSTANTS) or {}).get("constants", []):
        add(c.get("name"), "henry_constants")
    for k, v in data_access.load_json(data_paths.FLAVOR_REFERENCE_PAYLOADS).items():
        if isinstance(v, list):
            for r in v:
                if isinstance(r, dict):
                    add(r.get("compound"), f"flavor:{k}")
    for r in data_access.load_json(data_paths.SAFETY_REFERENCE_PAYLOADS).get("entries", []):
        add(r.get("analyte"), "safety:analyte")
    for k, v in data_access.load_json(data_paths.MATRIX_DECISION_PANEL).get("entries", {}).items():
        add(k, "panel")
        add(v.get("display_name"), "panel")
        for a in v.get("aliases", []):
            add(a, "panel")
    return seen


def test_every_repository_spelling_resolves_to_exactly_one_entry():
    unresolved = {
        spelling: sorted(where)
        for spelling, where in _repository_spellings().items()
        if compound_keys.resolve(spelling) is None
    }
    assert not unresolved, (
        "compound spellings under data/ that resolve to no registry entry "
        "(add a seed to scripts/generators/build_compound_registry.py and regenerate):\n"
        + "\n".join(f"  {s!r}: {w}" for s, w in unresolved.items())
    )


def test_no_spelling_is_claimed_by_two_entries():
    claims = Counter()
    for entry in compound_keys.all_compounds():
        for n in entry.norms:
            claims[n] += 1
    assert not [n for n, c in claims.items() if c > 1]


def test_molecules_have_unique_inchikeys_and_classes_resolve():
    keys = Counter(e.inchikey for e in compound_keys.all_compounds() if e.inchikey)
    assert not [k for k, c in keys.items() if c > 1]
    ids = {e.id for e in compound_keys.all_compounds()}
    for entry in compound_keys.all_compounds():
        for member in entry.members:
            assert member in ids, f"{entry.id} lists unknown member {member}"
        if entry.kind == "compound_class":
            assert compound_keys.members_of(entry), f"class {entry.id} has no members"
    molecules_without_structure = sorted(
        e.id for e in compound_keys.all_compounds() if e.kind == "molecule" and not e.inchikey
    )
    # The two nucleotides carry a CAS but no hand-entered structure. Anything else
    # appearing here is a regression.
    assert molecules_without_structure == ["gmp", "imp"]


def test_registry_covers_the_retired_hand_written_alias_table():
    table = compound_keys.benchmark_alias_table()
    missing = {}
    for canonical, aliases in RETIRED_BENCHMARK_NAME_ALIASES.items():
        expected = {norm(a) for a in aliases} | {canonical}
        got = table.get(canonical, frozenset())
        if not expected <= got:
            missing[canonical] = sorted(expected - got)
    assert not missing, f"registry match sets lost aliases the literal table had: {missing}"


def test_match_norms_is_exact_for_molecules_and_expands_classes():
    assert compound_keys.same_compound("Hexanal", "hexanal")
    assert compound_keys.same_compound("Furan-2-aldehyde (FA)", "furfural")
    assert not compound_keys.same_compound("Hexanal", "Nonanal")
    assert "25dimethylpyrazine" in compound_keys.match_norms("pyrazine")
    assert "nonanal" not in compound_keys.match_norms("hexanal")
    assert compound_keys.resolve("no such compound zzz") is None


def test_generator_output_is_current(tmp_path):
    """The committed registry must equal what the generator produces from the seeds."""
    import subprocess
    import sys

    completed = subprocess.run(
        [sys.executable, "scripts/generators/build_compound_registry.py", "--check"],
        cwd=data_paths.REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
