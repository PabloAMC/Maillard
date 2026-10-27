"""Data-side headline guards (moved 2026-09-03, retirement step B5b).

These pins survived the deletion of ``test_honest_headline_guards.py`` because they are
about the CURATED DATA, not about the retired engine: the ``no_verifiable_source`` census
that README.md and the model card quote. Same rule as before: a moved number must move
the prose in the same change.
"""
import json
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"


def _doc_text(path: Path) -> str:
    assert path.exists(), f"{path.name} is missing; the headline claims have no home"
    return path.read_text(encoding="utf-8")


def _assert_quoted(text: str, token: str, where: str, what: str) -> None:
    assert token in text, (
        f"{where} no longer quotes {token!r} for {what}. Either the number changed and this "
        f"guard was not updated with it, or the prose was edited away from the evidence. "
        f"Docs and evidence must not diverge silently."
    )


def _iter_records(node):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from _iter_records(value)
    elif isinstance(node, list):
        for value in node:
            yield from _iter_records(value)


def _contains_number(node) -> bool:
    if isinstance(node, bool):
        return False
    if isinstance(node, (int, float)):
        return True
    if isinstance(node, dict):
        return any(_contains_number(v) for v in node.values())
    if isinstance(node, list):
        return any(_contains_number(v) for v in node)
    return False


def _no_verifiable_source_census():
    """Return {relative_path: (flagged_count, flagged_records_containing_a_number)}.

    Scope and definitions, established by measurement 2026-08-27 and chosen because they
    reproduce README.md's published triple (102 / 80 / 62) exactly:

        flagged  = any object, at any depth, with source_status == "no_verifiable_source"
        numeric  = a flagged record containing at least one numeric value anywhere inside it
        scope    = every parseable .json / .yml / .yaml under data/ plus results/literature/

    No generator emits this census, so before this test the README's counts were
    unreproducible assertions. They are now recomputed on every run.
    """
    import yaml  # noqa: PLC0415 - only this census needs the YAML surface

    census = {}
    # 2026-09-02: the two literature ledgers live under results/literature/ now; they stay
    # in the census so the count is a property of the records, not of the file layout.
    scopes = ("data/**/", "results/literature/")
    for pattern in ("*.json", "*.yml", "*.yaml"):
        for path in sorted(p for scope in scopes for p in ROOT.glob(f"{scope}{pattern}")):
            try:
                text = path.read_text(encoding="utf-8")
                payload = json.loads(text) if path.suffix == ".json" else yaml.safe_load(text)
            except Exception:  # unparseable / non-record file; the citation gate owns those
                continue
            flagged = [
                record
                for record in _iter_records(payload)
                if isinstance(record.get("source_status"), str)
                and record["source_status"].strip() == "no_verifiable_source"
            ]
            if flagged:
                census[str(path.relative_to(ROOT))] = (
                    len(flagged),
                    sum(1 for record in flagged if _contains_number(record)),
                )
    return census


def test_no_record_carries_no_verifiable_source():
    """no_verifiable_source: 0 records. An unsourced record is deleted, never kept with a flag.

    RE-PINNED 2026-10-27: 84/63/63 -> 0/0/0, and the "runtime" third of the pin retired. All
    84 flagged records were DELETED (owner decision): 24 intake candidates, 29 research-backlog
    items, 18 SLR ledger rows, 7 flavour anchors, 3 retention records (the two pH-release
    surrogates read only by the retired `src/headspace.py`, and the Fadel 2015 retention), the
    Henle 2005 CML anchor, and the two flagged computational priors. Every one had been searched
    for (CrossRef sweeps 2026-08-26/27, gap search 2026-10-09) and was fabricated or could not
    be tied to its number. The intake citations survive as strings in
    `closed_unsourced_citations` so the research queue never reopens them. A fall by deletion,
    not verification.

    The retired "runtime" figure counted every flagged record under data/lit/ or
    results/literature/ as consumed at runtime. That was a location proxy, and wrong: the trace
    in `scripts/trace_runtime_data_reads.py` (pinned below) shows the front door never opens
    any of those files. What the model reads is now measured, not inferred.

    History of the earlier pins, kept because it explains how the flags got there. The last
    one read: 84 flagged · 63 carrying numbers · 63 of those "reaching runtime".

    RE-PINNED 2026-10-09: 87/65/65 -> 84/63/63. Six process-state calibration records with no
    source were DELETED from `data/lit/process_state_calibrations.json` (owner decision; their
    intake entries were already `no_verifiable_source` with a repair basis). Three of the six
    carried the flag themselves and two of those carried numbers; the other three were
    unflagged. A fall by deletion, not verification.

    RE-PINNED 2026-09-02: 102/80/80 -> 87/65/65. `data/lit/protein_source_registry.json`
    (15 records, every value self-labelled `mocked_placeholder`) was WITHDRAWN together
    with the code that consumed it; nothing was verified, the numbers simply no longer ship.
    Same day: the two literature ledgers moved to results/literature/ and stay in scope here.

    RE-PINNED 2026-09-01 (cleaning, Phase 1a): 120/98/80 -> 102/80/80. `data/qm/` was
    deleted together with the QM/DFT lane, taking its 18 unverifiable barrier records with
    it. Nothing runtime-reachable changed: those 18 never reached the model, so the runtime
    figure stays at 80 and the total falls by exactly the data/qm share. The history below
    is kept because it explains why the split (data/lit vs. everything else) was pinned.

    Pinned 2026-08-27 against README.md's "On literature provenance" note. The three
    numbers are pinned separately because they fail for different reasons and a maintainer
    needs to know which one moved.

    The 84 -> 102 step earlier the same day is the reason this guard exists at all. It was
    not a regression: `data/qm/` had been hidden from git by `.gitignore`, so every sweep
    this repo ever ran silently excluded it. Tracking it exposed 18 further records. That is
    exactly the failure mode a prose-only number cannot protect against -- the scope of the
    measurement changed and the published figure did not -- so the SPLIT is pinned below,
    not just the total.

    The runtime figure stays at 62 because the 18 data/qm records are read only by a loader
    test and by the skip-heavy Phase 3 authority lane; none reaches the model. 62 is the
    count that actually matters: an unverifiable citation attached to a number the runtime
    consumes is a fabricated parameter, not merely a bad footnote.

    RE-PINNED 2026-08-27 (Wave T3): 102/80/62 -> 120/98/80. This is the SECOND rise in one
    day and, like the first, it is a labelling change, not a data change. Eighteen records
    that were already shipping unverifiable numbers were finally labelled:

      * 15 in `data/lit/protein_source_registry.json` -- the file-level record plus all 14
        protein-source profiles. That file has always described itself as "Mocked values for
        14 protein sources based on Report 06 requirements", but the sentence sat in a JSON
        field nothing read while the numbers underneath it drove `matrix_uncertainty_factor`
        and the meaty-potential score at prediction time (Wave T1 finding T1-01).
      * 2 in `data/lit/retention_reference_payloads.json` -- the two `runtime_surrogate`
        blocks whose `log_slope = 0.235` is exactly ln(1.60)/2, back-solved from an invented
        "~55-65%" band in an in-repo brief (T1-02).
      * 1 in `data/lit/computational_priors.json` -- `ref41_ppi_sulfur_volatile_binding_v1`,
        which cited reference *number* 41 inside an LLM research dump (T1-08).

    NO VALUE WAS ADDED, CHANGED OR INVENTED IN THAT WAVE. All 18 are in data/lit and all
    carry numbers, so the runtime figure moves with the total: 62 -> 80. The honest reading
    is that 80 was always the true runtime figure and 62 was an undercount.

    The pin is zero and stays zero: a record whose source cannot be found is deleted (its
    citation string may go to `closed_unsourced_citations`), or its value becomes a declared
    assumption with a band in the code that uses it. It is never shipped with a flag.
    """
    census = _no_verifiable_source_census()
    assert census == {}, (
        "Records flagged source_status == no_verifiable_source are back: "
        f"{ {k: v[0] for k, v in census.items()} }. Delete them, or source them."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "No record in the repository is flagged `no_verifiable_source`",
                   "README.md", "the zero census")


#: Every curated file the front door opens, traced (scripts/trace_runtime_data_reads.py).
FRONT_DOOR_DATA_READS = [
    "data/keys/compounds.yml",
    "data/schemas/spec.schema.json",
    "data/species/desirable_targets.yml",
    "data/species/off_flavour_targets.yml",
    "data/species/protein_matrices.yml",
]


def test_the_front_door_reads_only_the_pinned_data_files():
    """compare / predict / explain / rank / score / calibrate, traced under an audit hook.

    Pinned 2026-10-27. The rate constants live in src/kinetic_core/ with their evidence class;
    the curated data the front door reads is these five files. A new data dependency of the
    model must show up here, and in README.md, in the same change.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    from trace_runtime_data_reads import trace_front_door  # noqa: PLC0415

    traced = trace_front_door()
    union = sorted({path for paths in traced.values() for path in paths})
    assert union == FRONT_DOOR_DATA_READS, f"front-door data reads moved: {traced}"

    readme = _doc_text(README)
    for path in FRONT_DOOR_DATA_READS:
        _assert_quoted(readme, path, "README.md", "the traced front-door data reads")

