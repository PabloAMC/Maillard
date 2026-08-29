"""
Build Wave V1 -- the usability surfaces: HTML report, network map, `explain`.

WHAT THESE TESTS ARE FOR
------------------------
V1 adds no science. It adds three ways to LOOK at the science, and the failure
mode of a presentation layer is that it quietly stops telling the truth: a
refusal rendered as an empty table, a NOT-RESOLVED ratio styled like a resolved
one, a report that needs the network to render, a network map that changes
between two runs so nobody can review its diff.

So these tests pin BEHAVIOUR AND HONESTY, never a predicted number. Every
assertion below would still hold after a refit; none of them would hold if the
presentation layer started softening an output.

  1. SELF-CONTAINMENT. No report may reference an external host, in any form.
     A lab laptop with the wifi off must render exactly what a connected one
     does.
  2. A REFUSAL IS AN ANSWER. An out-of-envelope request must produce first-class
     refusal cards carrying the engine's own named reasons -- not an empty
     table, not an error page.
  3. NOT RESOLVED IS NOT A SMALL EFFECT. A ratio inside the same-sample
     dispersion band must be rendered as its own visually distinct state.
  4. DETERMINISM. The network map must be byte-identical across two runs, or
     nobody can review a change to it.
  5. EVIDENCE TRAVELS. `explain` must print the evidence class of every route,
     because "where does this number come from" is the question the command
     exists to answer.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import comparative_cli  # noqa: E402
from src import explain_compound  # noqa: E402
from src import report_html  # noqa: E402

CLI = ROOT / "scripts" / "maillard.py"
NETWORK_MAP_GENERATOR = ROOT / "scripts" / "generators" / "generate_network_map.py"

#: A trunk-lane spec: the cheapest thing that produces a real, answered payload.
TRUNK_SPEC = {
    "name": "glucose_glycine_trunk",
    "precursors": {"D-Glucose": 10.0, "Glycine": 10.0},
    "temp_C": 120.0,
    "time_min": 10.0,
    "ph": 6.8,
    "aw": 1.0,
    "protein_type": "free",
}

#: An out-of-envelope spec. `declare_envelope` short-circuits before any
#: integration, so this is fast as well as being the case that matters most.
REFUSED_SPEC = {
    "name": "refused_targets",
    "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
    # WAVE B7 SWAPPED THE REFUSED EXAMPLE. This spec's job is to be
    # OUT OF ENVELOPE, and "HMF" stopped being a refusal when B7 gave 5-HMF a
    # route (two parallel sources ingested whole from Kocadagli & Gokmen 2016's
    # amine-free glucose system, plus Hamzalioglu's measured sink). "HEMF" --
    # homofuraneol -- replaces it and is refused for a SHARPER reason: the
    # compound is real and its route is understood, but it needs a C2 Strecker
    # donor and the core cannot put alanine and a pentose in the same lane.
    "targets": ["2-methyl-3-furanthiol (MFT)", "HEMF", "2-pentylfuran"],
    "temp_C": 140.0,
    "time_min": 30.0,
    "ph": 5.0,
    "aw": 0.98,
    "protein_type": "free",
}

#: Every scheme that would make the page reach off the machine.
_EXTERNAL_REFERENCE = re.compile(
    r"https?://|//cdn\.|\bsrc\s*=\s*[\"']//|@import\s+url\(|"
    r"\bfetch\s*\(|XMLHttpRequest|new\s+WebSocket|integrity\s*=",
    re.IGNORECASE,
)


def _squash(text: str) -> str:
    """Collapse the renderers' soft wrapping so a long anchor can be matched."""
    return " ".join(str(text).split())


@pytest.fixture(scope="module")
def trunk_payload():
    return comparative_cli.predict_core(TRUNK_SPEC)


@pytest.fixture(scope="module")
def refused_payload():
    return comparative_cli.predict_core(REFUSED_SPEC)


# ---------------------------------------------------------------------------
# 1. Self-containment
# ---------------------------------------------------------------------------


def test_predict_report_references_no_external_host(trunk_payload, tmp_path):
    target = report_html.write_report(trunk_payload, tmp_path / "predict.html")
    text = target.read_text(encoding="utf-8")
    match = _EXTERNAL_REFERENCE.search(text)
    assert match is None, (
        f"the report reached off the machine at {match.group(0)!r}. A report that "
        "needs the network is a report that renders differently in a lab."
    )
    assert text.startswith("<!doctype html>")
    assert "<style>" in text, "the CSS must be inline, not linked"


def test_report_carries_the_model_version_and_the_validation_headline(
    trunk_payload, tmp_path
):
    """The footer's job: you can always tell WHICH model wrote a printed page."""
    text = report_html.write_report(
        trunk_payload, tmp_path / "predict.html"
    ).read_text(encoding="utf-8")
    assert "model version" in text
    headline = report_html.model_card_headline()
    if headline["available"]:
        # Not the exact sentence -- it moves with the artifacts, which is the
        # point of reading it rather than writing it.
        assert "Absolute concentrations are unreliable" in text
    else:  # pragma: no cover - only on a checkout without the artifact
        assert "cannot show the generated validation headline" in text


def test_every_absolute_in_the_report_is_shown_with_an_interval(
    trunk_payload, tmp_path
):
    text = report_html.write_report(
        trunk_payload, tmp_path / "predict.html"
    ).read_text(encoding="utf-8")
    assert "The point is not the answer; the interval is." in text
    assert "HS-SPME same-sample dispersion" in text, (
        "the band that sets the floor width of every interval must be declared"
    )
    assert "air/water partition constant" in text


def test_declared_assumptions_are_lane_specific_not_a_blanket_list():
    """A trunk run must not claim to have used the lipid lane's Q10."""
    trunk = report_html.declared_assumptions(
        [{"lane": "trunk", "lanes": ["trunk"]}], [{"matrix": "water"}]
    )
    names = {row["name"] for row in trunk}
    assert "Q10 of hydroperoxide decomposition" not in names
    assert any("network pH" in name for name in names)

    lipid = report_html.declared_assumptions(
        [{"lane": "lipid", "lanes": ["lipid"], "lipid_carriers": ["pea_protein_isolate"]}],
        [{"matrix": "pea_protein_isolate"}],
    )
    lipid_names = {row["name"] for row in lipid}
    assert "Q10 of hydroperoxide decomposition" in lipid_names
    assert any("peroxide value" in name for name in lipid_names)
    for row in lipid:
        assert row["band"], f"{row['name']} was listed with no band"


# ---------------------------------------------------------------------------
# 2. A refusal is an answer
# ---------------------------------------------------------------------------


def test_an_out_of_envelope_request_renders_refusal_cards_with_named_reasons(
    refused_payload, tmp_path
):
    assert refused_payload["answered"] is False
    text = report_html.write_report(
        refused_payload, tmp_path / "refused.html"
    ).read_text(encoding="utf-8")

    assert 'class="card refusal"' in text, "refusals must be first-class cards"
    assert "No number is emitted" in text
    # The engine's own named reasons, verbatim -- not a generic apology.
    assert "Strecker donor" in text and "do not compose" in text
    assert "not in Frankel 1989" in text or "six-product slate" in text
    # And the page must say what CAN be asked instead.
    assert "Targets each lane can report" in text
    assert "Compounds the core deliberately refuses" in text


def test_the_refused_report_contains_no_predicted_number_table(
    refused_payload, tmp_path
):
    """The one thing a refusal page must never grow is a ranking."""
    text = report_html.write_report(
        refused_payload, tmp_path / "refused.html"
    ).read_text(encoding="utf-8")
    assert "<h2>Ranking</h2>" not in text
    assert "Odour activity (log scale" not in text


def test_the_cli_names_the_lane_and_the_missing_species_on_stderr(refused_payload):
    message = comparative_cli.envelope_error_text(refused_payload["declaration"])
    assert "OUT OF ENVELOPE" in message
    assert "Lane resolved: sulfur" in message
    assert "missing species" in message
    assert "HEMF" in message and "2-pentylfuran" in message
    # It reuses the EnvelopeDeclaration's own text rather than paraphrasing it.
    for reason in refused_payload["declaration"]["reasons"]:
        assert reason.split(":")[0] in message


def test_a_not_comparable_compare_renders_the_refusal_rather_than_an_empty_table(
    tmp_path,
):
    payload = {
        "artifact": "maillard_compare_core",
        "a": {"name": "A", "spec": TRUNK_SPEC},
        "b": {"name": "B", "spec": TRUNK_SPEC},
        "comparison": {
            "comparable": False,
            "reason": "at least one arm is out of envelope; a ratio against a refusal is not a ratio.",
            "declaration_a": {
                "state": "out_of_envelope",
                "lane": None,
                "lanes": [],
                "reasons": ["UNMAPPED PRECURSORS 'wheat flour': not a species in any core lane."],
                "warnings": [],
                "unmapped_precursors": ["wheat flour"],
                "unrepresented_targets": [],
                "mapped_precursors": {},
                "mapped_targets": {},
            },
            "declaration_b": {"state": "in_envelope", "lane": "trunk", "lanes": ["trunk"]},
        },
        "caveats": {"core": "core caveat"},
    }
    text = report_html.render_compare_report(payload)
    assert "Not comparable" in text
    assert "a ratio against a refusal is not a ratio" in text
    assert "wheat flour" in text
    assert 'class="card refusal"' in text


# ---------------------------------------------------------------------------
# 3. NOT RESOLVED is a state, not a small effect
# ---------------------------------------------------------------------------


def _compare_payload_with(rows):
    return {
        "artifact": "maillard_compare_core",
        "a": {"name": "A", "spec": TRUNK_SPEC},
        "b": {"name": "B", "spec": TRUNK_SPEC},
        "comparison": {
            "comparable": True,
            "declaration_a": {"state": "in_envelope", "lane": "sulfur", "lanes": ["sulfur"]},
            "declaration_b": {"state": "in_envelope", "lane": "sulfur", "lanes": ["sulfur"]},
            "run_a": {},
            "run_b": {},
            "ratios": {
                "formulation_a": "A",
                "formulation_b": "B",
                "n_compared": len(rows),
                # Q1: this fixture used to compute n_resolved as
                # ``not within_reliability_band``, which is the arithmetic the
                # science layer USED to do -- so it silently kept testing the
                # bug after the layer was fixed. It now mirrors
                # ``matrix_oav.compare_formulations``: an UNDEFINED ratio
                # resolves nothing and is counted separately.
                "n_resolved": sum(
                    1 for r in rows
                    if r["direction"] != "undefined"
                    and not r["within_reliability_band"]
                ),
                "n_undefined": sum(1 for r in rows if r["direction"] == "undefined"),
                "reliability_band_x": 4.79583152331272,
                "rows": rows,
            },
        },
        "caveats": {"core": "core caveat", "ratio": "ratio caveat"},
    }


def test_a_not_resolved_ratio_is_rendered_in_its_own_visually_distinct_state():
    payload = _compare_payload_with(
        [
            {
                "compound": "MFT",
                "ratio_a_over_b": 31.0,
                "direction": "higher_in_A",
                "within_reliability_band": False,
                "note": "resolved above the same-sample dispersion band",
            },
            {
                "compound": "FUR",
                "ratio_a_over_b": 1.4,
                "direction": "higher_in_A",
                "within_reliability_band": True,
                "note": "NOT RESOLVED: inside the measured same-sample dispersion band",
            },
        ]
    )
    text = report_html.render_compare_report(payload)

    # A distinct row CLASS, and a distinct badge -- two independent signals, so
    # a stylesheet change cannot silently collapse the two states into one.
    assert 'class="notresolved"' in text
    assert "not resolved" in text and "resolved" in text
    assert 'class="badge b-unres"' in text
    assert 'class="badge b-ok"' in text
    # And the page must say what NOT RESOLVED means, in words.
    assert "not as a small effect" in text
    assert "1 of 2" in text or "<strong>1</strong> of 2" in text
    # The stylesheet must actually distinguish the row.
    assert "tr.notresolved td" in text


def test_a_zero_arm_is_reported_as_undefined_and_never_as_an_enormous_ratio():
    payload = _compare_payload_with(
        [
            {
                "compound": "MFTD",
                "ratio_a_over_b": float("inf"),
                "direction": "undefined",
                "within_reliability_band": False,
                "note": "B is zero; a ratio is undefined",
            }
        ]
    )
    text = report_html.render_compare_report(payload)
    assert "A only (B is zero)" in text
    # No raw float infinity may reach a table cell. Checked as a CELL, not as a
    # substring: the word "informative" appears legitimately in a hold-out note.
    assert ">inf<" not in text.lower()
    assert "inf&times;" not in text.lower()
    assert "infinity" not in text.lower()


def test_the_compare_report_references_no_external_host():
    payload = _compare_payload_with(
        [
            {
                "compound": "MFT",
                "ratio_a_over_b": 31.0,
                "direction": "higher_in_A",
                "within_reliability_band": False,
                "note": "resolved",
            }
        ]
    )
    assert _EXTERNAL_REFERENCE.search(report_html.render_compare_report(payload)) is None


# ---------------------------------------------------------------------------
# 4. The network map: determinism and honesty
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_the_network_map_is_byte_identical_across_two_runs(tmp_path):
    """
    A generated artifact that changes between runs cannot be reviewed.

    Two full generations, same tree, same spec set, compared as BYTES. This is
    the guard on the whole determinism discipline in the generator: no
    timestamp, sorted iteration, fixed float formatting.
    """
    first, second = tmp_path / "a.html", tmp_path / "b.html"
    for target in (first, second):
        proc = subprocess.run(
            [sys.executable, str(NETWORK_MAP_GENERATOR), "--out", str(target)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=900,
        )
        assert proc.returncode == 0, proc.stderr[-3000:]
    assert first.read_bytes() == second.read_bytes(), (
        "the network map differs between two runs on the same tree; something in "
        "the generator is reading a clock, a set iteration order or an unrounded "
        "float"
    )


@pytest.mark.slow
def test_the_network_map_is_self_contained_and_carries_all_four_lanes(tmp_path):
    target = tmp_path / "map.html"
    proc = subprocess.run(
        [sys.executable, str(NETWORK_MAP_GENERATOR), "--out", str(target)],
        cwd=str(ROOT), capture_output=True, text=True, timeout=900,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    text = target.read_text(encoding="utf-8")

    assert _EXTERNAL_REFERENCE.search(text) is None, "the map reached off the machine"
    for lane in ("trunk", "sulfur", "acrylamide", "lipid"):
        assert f"<h2>{lane} lane</h2>" in text, f"the {lane} lane is missing"
    for evidence in ("measured", "fitted", "derived", "pinned"):
        assert f'data-ev="{evidence}"' in text, f"no edge is classed {evidence}"
    # Flux mode ran, and the widths mean something.
    assert "TIME-INTEGRATED ABSOLUTE" in text
    assert "Dominant steps" in text
    assert "flux, mmol/L" in text
    # Tooltips carry provenance, not just a name.
    assert "validity:" in text and "source:" in text


def test_the_map_layering_is_cycle_safe():
    """
    Glucose <-> fructose is a real cycle. A longest-path layering unrolls it and
    produced a 12 000-pixel-wide lane on the first attempt; this pins the fix.
    """
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    import generate_network_map as gnm

    for lane in ("trunk", "sulfur", "acrylamide"):
        edges = gnm.edge_records(lane)
        species = gnm.lane_species(lane)
        _positions, width, _height = gnm.layout(species, edges)
        assert width < 2600, (
            f"the {lane} lane laid out {width} px wide -- the layering is being "
            "inflated by a cycle again"
        )


# ---------------------------------------------------------------------------
# 5. `explain` -- evidence travels with every route
# ---------------------------------------------------------------------------


def test_explain_prints_an_evidence_class_for_every_formation_route():
    payload = explain_compound.explain("MFT")
    assert payload["answered"] is True
    assert payload["species_key"] == "MFT"
    assert payload["lane"] == "sulfur"
    assert payload["routes"], "MFT must have formation routes in the sulfur lane"
    for route in payload["routes"]:
        assert route["evidence"] in explain_compound.EVIDENCE_ORDER
        assert route["step"]

    text = explain_compound.render_explain_text(payload)
    assert "EVIDENCE CLASSES" in text
    for evidence in explain_compound.EVIDENCE_ORDER:
        assert evidence in text, f"the {evidence} class is not explained in the output"
    assert "FORMATION ROUTES" in text
    assert "TOP LITERATURE ANCHORS" in text
    # The anchors are the registries' own strings, not a paraphrase.
    assert payload["anchors"]
    # Whitespace-normalised: the renderer soft-wraps long anchors, and the point
    # is that the registry's own string survives, not its line breaks.
    assert _squash(payload["anchors"][0]["source_anchor"]) in _squash(text)


def test_explain_marks_the_lipid_rate_as_an_assumption_and_the_branches_as_fitted():
    payload = explain_compound.explain("hexanal")
    assert payload["lane"] == "lipid"
    classes = {route["evidence"] for route in payload["routes"]}
    assert explain_compound.FITTED in classes, "the branch model is a fit"
    assert explain_compound.PINNED in classes, (
        "the lipid RATE is a declared assumption and must never be shown as measured"
    )
    text = _squash(explain_compound.render_explain_text(payload)).upper()
    assert "RATE IS AN ASSUMPTION" in text


def test_explain_answers_a_refused_compound_with_the_declared_reason():
    # B7: "HMF" is answered now; "HEMF" is the sharper refusal that replaces it.
    payload = explain_compound.explain("HEMF")
    assert payload["answered"] is False
    assert payload["state"] == "refused"
    assert "Strecker donor" in payload["reason"]
    text = explain_compound.render_explain_text(payload)
    assert "REFUSED" in text
    assert "A refusal is an answer" in text
    # It must also say what CAN be explained.
    assert "2-methyl-3-furanthiol (MFT)" in text


def test_explain_does_not_crash_on_an_unknown_name_and_says_so_honestly():
    payload = explain_compound.explain("unobtainium sulfoxide")
    assert payload["answered"] is False
    assert payload["state"] == "no_vocabulary_entry"
    assert "gap in the vocabulary, not a claim about the chemistry" in payload["reason"]


def test_the_evidence_class_mapping_never_calls_an_assumption_measured():
    """The one mapping every downstream surface trusts."""
    assert explain_compound.evidence_class_of("measured_rate") == explain_compound.MEASURED
    assert (
        explain_compound.evidence_class_of("measured_activation_energy")
        == explain_compound.MEASURED
    )
    assert (
        explain_compound.evidence_class_of("derived_from_fit_data", ("fitted_here",))
        == explain_compound.FITTED
    )
    assert (
        explain_compound.evidence_class_of("declared_assumption")
        == explain_compound.PINNED
    )
    assert (
        explain_compound.evidence_class_of("bounded_from_a_timescale_bracket")
        == explain_compound.PINNED
    )
    # An unrecognised class must fall to the WEAKEST reading, never the strongest.
    assert explain_compound.evidence_class_of("something_new") == explain_compound.PINNED


@pytest.mark.slow
def test_explain_runs_as_a_subprocess_and_exits_zero():
    proc = subprocess.run(
        [sys.executable, str(CLI), "explain", "MFT"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "FORMATION ROUTES" in proc.stdout
    assert "measured" in proc.stdout and "pinned" in proc.stdout
