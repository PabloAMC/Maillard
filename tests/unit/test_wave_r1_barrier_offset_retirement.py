"""Wave R1 (2026-08-28) — the barrier-offset auto-acceptance must not come back.

WHAT WAS RETIRED.  `data/lit/refinement_surrogate_patches.json` is read by
`src.barrier_constants.get_barrier()` and its `accepted_offsets` are ADDED to the audited
`FAST_BARRIERS` value before the barrier is returned.  Until 2026-08-28,
`src.refinement_campaign.build_refinement_impact_artifact()` filled that map
automatically with nine offsets of exactly +/-3.0 kcal/mol, accepted BECAUSE they improved
the benchmark panel the model is then scored against.  That is a fit to the evaluation set,
it was never declared to `scripts/ci/fit_target_gate.py`, and every offset sat exactly on
the +/-3.0 search bound -- a bound report, not an optimum.  With the offsets armed the
shipped barriers were `Schiff_Base_Formation` 18.0 (table 15.0),
`Retro_Aldol_Fragmentation` 29.0 (table 32.0) and `Thiol_Addition` 31.6 (table 28.6).  At
150 C, 3.0 kcal/mol is a ~35x rate factor.

WHY A GUARD RATHER THAN A DELETION.  The mechanism is a legitimate DIAGNOSTIC: the
campaign still searches, and still reports what it would have chosen, under
`candidate_offsets_not_applied`.  What must never happen again is that the search's output
reaches the shipped barrier table without being declared.  Two things therefore have to
hold at once and both are asserted here:

  1. the tracked patch file's `accepted_offsets` is empty, and
  2. `get_barrier(family)` returns the FAST_BARRIERS value for every family, so the file
     is not merely empty but demonstrably not overriding anything.

HOW THIS DEFECT WAS ABLE TO HIDE, measured in Wave R1 and the reason (1) is not enough on
its own: at git HEAD the tracked patch file carried an EMPTY `accepted_offsets`, so a clean
checkout shipped the audited table.  The offsets were written into that tracked file by
`scripts/generators/generate_refinement_governance.py`, which runs in the MIDDLE of
`scientific_lane()` in `scripts/docker_maillard.sh` -- ahead of the figure generators, the
envelope report and the scientific-regression pytest lane.  Running the repository's own
documented regeneration sequence therefore ARMED a +/-3.0 kcal/mol override and scored
everything after it against the armed model.  A guard that only reads the file on a clean
tree would have passed throughout.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.barrier_constants as barrier_constants

PATCH_FILE = ROOT / "data" / "lit" / "refinement_surrogate_patches.json"

# The nine offsets that were live before the retirement. Named here so that a reviewer
# reading only this test sees exactly what may not return.
RETIRED_OFFSETS = {
    "retro_aldol": -3.0,
    "retro_aldol_fragmentation": -3.0,
    "lipid_schiff_base": 3.0,
    "schiff": 3.0,
    "schiff_base_formation": 3.0,
    "schiff_condensation": 3.0,
    "thiol_addition": 3.0,
    "thiol_addition_h2": 3.0,
    "thiol_addition_legacy_shortcut": 3.0,
}

# The three families whose SHIPPED barrier was measurably wrong while the offsets were
# armed, with the audited table value each one must now return.
MEASURED_OVERRIDES = {
    "Schiff_Base_Formation": (18.0, 15.0),
    "Retro_Aldol_Fragmentation": (29.0, 32.0),
    "Thiol_Addition": (31.6, 28.6),
}


@pytest.fixture(autouse=True)
def _no_runtime_offsets(monkeypatch):
    """No test in this file may run with a runtime override in the environment.

    `get_barrier()` merges `os.environ["BARRIER_OFFSETS"]` on top of the patch file, so a
    stray variable in the shell would make these assertions meaningless in either
    direction. The module-level patch cache is also cleared, because it is keyed on the
    file's mtime and another test may have pointed it at a tmp_path fixture.
    """
    monkeypatch.delenv("BARRIER_OFFSETS", raising=False)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_CACHE", None)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_MTIME", None)


def test_accepted_offsets_in_the_tracked_patch_file_is_empty():
    payload = json.loads(PATCH_FILE.read_text(encoding="utf-8"))
    accepted = payload.get("accepted_offsets")

    assert accepted == {}, (
        "data/lit/refinement_surrogate_patches.json has a NON-EMPTY accepted_offsets. "
        "That map is added to the audited FAST_BARRIERS value inside get_barrier(), so "
        "every shipped prediction in this working tree is running on barriers the "
        "documented table does not contain. This is the 2026-08-28 defect returning.\n"
        f"  found: {accepted}\n"
        "If a regeneration pass put it there, re-read the retirement note in the patch file before changing anything. "
        "Offsets may only ever be applied through an EXPLICIT, DECLARED fit record that "
        "scripts/ci/fit_target_gate.py can see, with the fitted rows removed from scored "
        "evidence -- never by a search that optimises the evaluation panel."
    )
    assert payload.get("accepted_candidates") == [], (
        "accepted_candidates is non-empty: the campaign accepted something. See above."
    )


def test_the_retirement_note_keeps_the_provenance_of_what_was_retired():
    """The tracked patch file must keep the record of the nine retired offsets.

    Until the 2026-08-30 QM/DFT prune the note was also emitted from
    `src.refinement_campaign.RETIREMENT_NOTE`; that module and the generator that could
    overwrite the file are gone, so the file itself is now the only copy and this test
    guards it directly.
    """
    payload = json.loads(PATCH_FILE.read_text(encoding="utf-8"))
    note = payload.get("retirement_note") or {}

    assert note.get("date") == "2026-08-28"
    assert note.get("retired_offsets_kept_for_provenance") == RETIRED_OFFSETS, (
        "the provenance record of what was retired has been altered or lost"
    )


@pytest.mark.parametrize("family", sorted(barrier_constants.FAST_BARRIERS))
def test_get_barrier_returns_the_documented_table_value_for_every_family(family):
    """No family may ship a barrier that differs from the audited table.

    This is the assertion that (a) the patch file is empty AND (b) nothing else has
    quietly grown an override path. It is parametrised over the whole table rather than a
    sample because the retired offsets reached families through THREE different routes --
    an exact key (`thiol_addition`), a normalised alias (`retro_aldol_fragmentation`) and
    the short-key substring map inside get_barrier (`schiff` -> any family whose
    normalised name contains `schiff_condensation`) -- and a hand-picked sample would have
    to guess which route the next one takes.
    """
    barrier, _ = barrier_constants.get_barrier(family)

    assert barrier == pytest.approx(barrier_constants.FAST_BARRIERS[family][0]), (
        f"get_barrier({family!r}) does not return the FAST_BARRIERS table value. "
        "Something is adding an offset to the audited barrier."
    )


@pytest.mark.parametrize("family,armed,table", [(k, v[0], v[1]) for k, v in MEASURED_OVERRIDES.items()])
def test_the_three_measured_overrides_are_gone(family, armed, table):
    """The public-facing family names whose shipped barrier was measured as wrong.

    These are the exact numbers recorded on 2026-08-28 with the offsets armed. They are
    pinned in both directions: the shipped value must equal the table, and must NOT equal
    what the armed model returned.
    """
    barrier, _ = barrier_constants.get_barrier(family)

    assert barrier == pytest.approx(table)
    assert barrier != pytest.approx(armed), (
        f"{family} is shipping {armed}, the value it had while the +/-3.0 kcal/mol "
        "auto-accepted offsets were armed. At 150 C that is a ~35x rate error."
    )


def test_the_short_key_substring_map_cannot_be_armed_from_an_empty_file():
    """`schiff`/`amadori`/`enol`/`strecker`/`cys` are SUBSTRING keys inside get_barrier.

    A single `schiff` entry reaches every family whose normalised name contains
    `schiff_condensation`, so the blast radius of one line in the patch file is much
    larger than the line suggests. With the file empty the map must be inert; this pins
    that, and names the mechanism so the next reader does not have to rediscover it.
    """
    for family in ("schiff_condensation", "amadori_rearrangement", "strecker_degradation",
                   "cysteine_thermolysis", "1,2-enolisation"):
        barrier, _ = barrier_constants.get_barrier(family)
        assert barrier == pytest.approx(barrier_constants.FAST_BARRIERS[family][0])
