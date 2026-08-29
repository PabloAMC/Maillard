"""
src/report_html.py -- THE SELF-CONTAINED HTML REPORT (Build Wave V1, usability).

WHAT THIS IS
------------
One file in, one file out: a ``predict`` or ``compare`` payload from
``src/comparative_cli.py`` becomes a single ``.html`` document that a bench
scientist can open, read, print, and mail to a colleague. There is NO network
dependency of any kind -- the CSS is inline, the JavaScript is inline (and
optional), and every chart is hand-rolled SVG. Opening the file on a laptop
with the wifi off must render exactly what it renders online, because a report
that silently degrades is a report you cannot trust in a lab.

WHAT IT MAY AND MAY NOT DO
--------------------------
May: read a payload, re-key it, call the B4 output layer's PURE arithmetic
(``matrix_oav.oav_table``, ``predict_matrix_shift``, ``decompose_residual``),
read the frozen registries for the declared assumptions it must list, and lay
the result out.

May NOT: integrate anything, fit anything, invent a number, or soften a
refusal. In particular:

  * a REFUSAL IS AN ANSWER and is rendered as a first-class card with its named
    reason, not as an empty table or an error toast;
  * a ratio inside the same-sample dispersion band is rendered NOT RESOLVED and
    is visually distinct from a resolved one -- never as a small effect;
  * an absolute is never drawn without its interval;
  * every declared assumption the run actually used is listed in the footer
    with its band, because the bands are the product.

THIS MODULE DOES NO ARITHMETIC OF ITS OWN (Q1)
----------------------------------------------
It used to do one piece. ``engine.compare`` returned its arms flattened through
``CorePrediction.as_dict()``, which drops the object and therefore drops
``.oav()``, so a compare carried no OAV table -- and this module rebuilt one by
hand from the run dict. The V1 wave recorded that as API feedback and predicted
it would drift. It had: the copy was written in B6 and never learned about B7's
furanone declared-assumption band, so a compare drew NARROWER intervals than a
predict of the identical arm.

Q1 fixed the API instead of the copy. ``engine.compare`` emits ``oav_table_a``
/ ``oav_table_b`` and ``rows_a``/``rows_b`` from the live objects, and this
module reads them. Likewise the three key-spaces a compound wears (display
name, species key, B4 key) are converted ONLY through
``kinetic_core.keyspaces``, never by an inline ``.get`` here -- getting that
hop wrong does not raise, it silently reports "no measured threshold" for a
compound that has one.
"""

from __future__ import annotations

import html as _html
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
MODEL_CARD_PATH = ROOT / "results" / "validation" / "model_card.json"

#: Rendered above the charts, once, in the loudest style the page has.
INTERVAL_RULE = (
    "Every bar on this page is drawn on a LOG axis with its measured interval. "
    "The point is not the answer; the interval is."
)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------


def esc(value: Any) -> str:
    """HTML-escape anything, including None."""
    if value is None:
        return ""
    return _html.escape(str(value), quote=True)


def sig(value: Optional[float], places: int = 3) -> str:
    """Significant-figure formatting that never lies about precision."""
    if value is None:
        return "&mdash;"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return esc(value)
    if number != number:  # NaN
        return "NaN"
    if number in (float("inf"), float("-inf")):
        return "&infin;" if number > 0 else "-&infin;"
    if number == 0.0:
        return "0"
    if abs(number) >= 1e5 or abs(number) < 1e-3:
        mantissa, exponent = f"{number:.{places - 1}e}".split("e")
        return f"{mantissa}&times;10<sup>{int(exponent)}</sup>"
    return f"{number:.{places}g}"


def _git_version() -> str:
    """``git describe`` for the footer, or an honest placeholder."""
    for command in (
        ["git", "describe", "--tags", "--always", "--dirty"],
        ["git", "rev-parse", "--short", "HEAD"],
    ):
        try:
            out = subprocess.run(
                command, cwd=str(ROOT), capture_output=True, text=True, timeout=10
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    return "unknown (not a git checkout)"


def model_card_headline() -> Dict[str, Any]:
    """
    The honest validation one-liner, read from the generated model card.

    Read, never written: if the artifact is missing the footer says so rather
    than substituting a cheerful sentence of its own.
    """
    if not MODEL_CARD_PATH.exists():
        return {
            "available": False,
            "sentences": [],
            "note": (
                f"{MODEL_CARD_PATH.relative_to(ROOT)} is absent from this "
                "checkout, so this report cannot show the generated validation "
                "headline. Regenerate it with "
                "scripts/generators/generate_model_card.py. NO substitute "
                "sentence is shown: an unbacked claim is worse than a gap."
            ),
        }
    try:
        card = json.loads(MODEL_CARD_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover - defensive
        return {"available": False, "sentences": [], "note": f"unreadable: {exc}"}
    sentences = [str(s) for s in (card.get("headline_sentences") or [])]
    return {
        "available": bool(sentences),
        "sentences": sentences,
        "one_liner": (
            "Compare two formulations and read the RATIO. Never quote the "
            "absolute number. Treat any pH or moisture direction as "
            "unsupported."
        ),
        "note": "",
    }


def _markdownish(text: str) -> str:
    """The model card's sentences carry **bold** markers. Honour just that one."""
    escaped = esc(text)
    parts = escaped.split("**")
    out: List[str] = []
    for index, part in enumerate(parts):
        out.append(f"<strong>{part}</strong>" if index % 2 else part)
    return "".join(out)


# ---------------------------------------------------------------------------
# SVG charts -- hand-rolled, no libraries, no external anything
# ---------------------------------------------------------------------------


def _decade_ticks(lo10: int, hi10: int, max_ticks: int = 9) -> List[int]:
    span = hi10 - lo10
    step = max(1, math.ceil(span / max(max_ticks - 1, 1)))
    return list(range(lo10, hi10 + 1, step)) or [lo10]


def _tick_label(exponent: int) -> str:
    if -3 <= exponent <= 4:
        value = 10.0 ** exponent
        return f"{value:g}"
    return f"1e{exponent}"


def svg_log_bars(
    rows: Sequence[Mapping[str, Any]],
    *,
    x_title: str,
    width: int = 900,
    row_height: int = 30,
    label_width: int = 250,
    reference: Optional[float] = None,
    reference_label: str = "",
) -> str:
    """
    Horizontal bars on a base-10 log axis, each with an interval whisker.

    ``rows`` entries: ``label``, ``point``, ``lo``, ``hi``, optional
    ``highlight`` (bool), ``tag`` (short badge string) and ``note``.

    A row whose point is zero or absent is drawn as an explicit ZERO MARKER at
    the left edge with the word "0" -- not dropped, because a compound the model
    says is absent is a claim and deserves a line.
    """
    rows = list(rows)
    if not rows:
        return "<p class='empty'>No compound in this run carries a value to plot.</p>"

    positive = [
        float(v)
        for row in rows
        for v in (row.get("lo"), row.get("point"), row.get("hi"))
        if isinstance(v, (int, float)) and v and float(v) > 0.0
    ]
    if not positive:
        return (
            "<p class='empty'>Every value in this run is zero or unavailable, "
            "so there is no log axis to draw.</p>"
        )

    lo10 = math.floor(math.log10(min(positive)))
    hi10 = math.ceil(math.log10(max(positive)))
    if hi10 <= lo10:
        hi10 = lo10 + 1
    ticks = _decade_ticks(lo10, hi10)
    if ticks[-1] != hi10:
        ticks.append(hi10)

    top_pad, bottom_pad = 26, 44
    plot_x = label_width + 12
    plot_w = width - plot_x - 108
    height = top_pad + row_height * len(rows) + bottom_pad

    def x_of(value: float) -> float:
        clamped = max(min(math.log10(value), hi10), lo10)
        return plot_x + plot_w * (clamped - lo10) / (hi10 - lo10)

    out: List[str] = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" width="100%" '
        f'height="{height}" role="img" aria-label="{esc(x_title)}">',
        f'<rect x="{plot_x}" y="{top_pad - 6}" width="{plot_w}" '
        f'height="{row_height * len(rows) + 8}" class="plotbg"/>',
    ]

    for exponent in ticks:
        x = x_of(10.0 ** exponent)
        out.append(
            f'<line x1="{x:.1f}" y1="{top_pad - 6}" x2="{x:.1f}" '
            f'y2="{top_pad + row_height * len(rows) + 2}" class="grid"/>'
        )
        out.append(
            f'<text x="{x:.1f}" y="{top_pad + row_height * len(rows) + 18}" '
            f'class="tick" text-anchor="middle">{esc(_tick_label(exponent))}</text>'
        )

    if reference is not None and reference > 0:
        x = x_of(float(reference))
        out.append(
            f'<line x1="{x:.1f}" y1="{top_pad - 10}" x2="{x:.1f}" '
            f'y2="{top_pad + row_height * len(rows) + 2}" class="refline"/>'
        )
        if reference_label:
            out.append(
                f'<text x="{x:.1f}" y="{top_pad - 14}" class="reflabel" '
                f'text-anchor="middle">{esc(reference_label)}</text>'
            )

    for index, row in enumerate(rows):
        y = top_pad + index * row_height
        mid = y + row_height / 2 - 2
        point = row.get("point")
        lo, hi = row.get("lo"), row.get("hi")
        highlight = bool(row.get("highlight"))
        klass = "bar hi" if highlight else "bar"

        label = esc(row.get("label", "?"))
        tag = row.get("tag")
        out.append(
            f'<text x="{label_width}" y="{mid + 4}" class="rowlabel'
            f'{" hi" if highlight else ""}" text-anchor="end">{label}</text>'
        )
        if tag:
            out.append(
                f'<text x="{label_width}" y="{mid + 15}" class="rowtag" '
                f'text-anchor="end">{esc(tag)}</text>'
            )

        if not isinstance(point, (int, float)) or float(point) <= 0.0:
            out.append(
                f'<text x="{plot_x + 4}" y="{mid + 4}" class="zero">'
                f'0 &mdash; the model says this compound is not formed</text>'
            )
            continue

        x_point = x_of(float(point))
        out.append(
            f'<rect x="{plot_x}" y="{y + 7}" width="{max(x_point - plot_x, 1):.1f}" '
            f'height="{row_height - 16}" class="{klass}"/>'
        )
        if (
            isinstance(lo, (int, float))
            and isinstance(hi, (int, float))
            and float(lo) > 0.0
            and float(hi) > float(lo)
        ):
            x_lo, x_hi = x_of(float(lo)), x_of(float(hi))
            out.append(
                f'<line x1="{x_lo:.1f}" y1="{mid:.1f}" x2="{x_hi:.1f}" '
                f'y2="{mid:.1f}" class="whisker"/>'
            )
            for cap in (x_lo, x_hi):
                out.append(
                    f'<line x1="{cap:.1f}" y1="{y + 5}" x2="{cap:.1f}" '
                    f'y2="{y + row_height - 9}" class="whisker"/>'
                )
        out.append(
            f'<text x="{plot_x + plot_w + 8}" y="{mid + 4}" class="value">'
            f'{sig(point)}</text>'
        )

    out.append(
        f'<text x="{plot_x + plot_w / 2:.0f}" y="{height - 6}" class="axistitle" '
        f'text-anchor="middle">{esc(x_title)}</text>'
    )
    out.append("</svg>")
    return "".join(out)


def svg_stacked_decades(
    rows: Sequence[Mapping[str, Any]],
    *,
    width: int = 900,
    row_height: int = 34,
    label_width: int = 250,
) -> str:
    """
    Stacked horizontal bars in DECADES, for the residual decomposition.

    ``rows`` entries: ``label`` and ``segments``: a sequence of
    ``(name, decades, css_class)``. Negative segments run leftwards from zero;
    this is deliberate, because a matrix that makes a compound MORE volatile is
    a real, measured outcome and folding it to a magnitude would hide the sign.
    """
    rows = list(rows)
    if not rows:
        return "<p class='empty'>No compound in this run has a decomposable matrix shift.</p>"

    extents = [0.0]
    for row in rows:
        running_pos = running_neg = 0.0
        for _, decades, _ in row["segments"]:
            if decades >= 0:
                running_pos += decades
            else:
                running_neg += decades
        extents.extend([running_pos, running_neg])
    hi = max(max(extents), 0.5)
    lo = min(min(extents), 0.0)
    span = max(hi - lo, 0.5)

    top_pad, bottom_pad = 20, 40
    plot_x = label_width + 12
    plot_w = width - plot_x - 96
    height = top_pad + row_height * len(rows) + bottom_pad

    def x_of(value: float) -> float:
        return plot_x + plot_w * (value - lo) / span

    zero_x = x_of(0.0)
    out: List[str] = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" width="100%" '
        f'height="{height}" role="img" aria-label="residual decomposition in decades">',
        f'<rect x="{plot_x}" y="{top_pad - 4}" width="{plot_w}" '
        f'height="{row_height * len(rows) + 6}" class="plotbg"/>',
        f'<line x1="{zero_x:.1f}" y1="{top_pad - 4}" x2="{zero_x:.1f}" '
        f'y2="{top_pad + row_height * len(rows) + 2}" class="refline"/>',
    ]

    for index, row in enumerate(rows):
        y = top_pad + index * row_height
        mid = y + row_height / 2 - 2
        out.append(
            f'<text x="{label_width}" y="{mid + 4}" class="rowlabel" '
            f'text-anchor="end">{esc(row["label"])}</text>'
        )
        cursor_pos = cursor_neg = 0.0
        for name, decades, css in row["segments"]:
            if abs(decades) < 1e-12:
                continue
            if decades > 0:
                x0, x1 = x_of(cursor_pos), x_of(cursor_pos + decades)
                cursor_pos += decades
            else:
                x1, x0 = x_of(cursor_neg), x_of(cursor_neg + decades)
                cursor_neg += decades
            out.append(
                f'<rect x="{min(x0, x1):.1f}" y="{y + 8}" '
                f'width="{max(abs(x1 - x0), 1.0):.1f}" height="{row_height - 18}" '
                f'class="seg {css}"><title>{esc(name)}: '
                f'{decades:+.3f} decades ({10.0 ** decades:.3g}x)</title></rect>'
            )
        total = cursor_pos + cursor_neg
        out.append(
            f'<text x="{plot_x + plot_w + 8}" y="{mid + 4}" class="value">'
            f'{10.0 ** total:.3g}x</text>'
        )

    out.append(
        f'<text x="{plot_x + plot_w / 2:.0f}" y="{height - 6}" class="axistitle" '
        f'text-anchor="middle">log&#8321;&#8320; shift, decades (right of the '
        f'line = suppressed in matrix)</text>'
    )
    out.append("</svg>")
    return "".join(out)


# ---------------------------------------------------------------------------
# Payload readers
# ---------------------------------------------------------------------------


def _species_key_for(compound: str, declaration: Mapping[str, Any]) -> str:
    """DISPLAY NAME -> SPECIES KEY. See ``kinetic_core.keyspaces``."""
    from src.kinetic_core.keyspaces import species_key_for

    return species_key_for(compound, declaration.get("mapped_targets") or {})


def _oav_entry(oav_table_payload: Mapping[str, Any], species_key: str) -> Dict[str, Any]:
    """
    One arm's OAV entry for a species, looked up in the B4 KEY-SPACE.

    Q1: the species -> B4 hop used to be written out here as a third private
    copy. It lives in ``kinetic_core.keyspaces.b4_key`` now, because getting it
    wrong does not raise -- it returns "no measured threshold" for a compound
    that has one.
    """
    from src.kinetic_core.keyspaces import b4_key

    per_species = dict(oav_table_payload.get("per_species") or {})
    b4 = b4_key(species_key)
    entry = per_species.get(b4) if b4 else None
    return dict(entry) if isinstance(entry, Mapping) else {}


# ---------------------------------------------------------------------------
# Declared assumptions -- the footer's real payload
# ---------------------------------------------------------------------------


def declared_assumptions(
    declarations: Sequence[Mapping[str, Any]],
    metadatas: Sequence[Mapping[str, Any]],
) -> List[Dict[str, str]]:
    """
    Every declared assumption the run ACTUALLY USED, with its band.

    Driven off the run's own lane and metadata, not off a static list, so a
    trunk-only run does not claim to have used the lipid Q10 and a clamped-pH
    run does not claim the unbuffered default.
    """
    lanes = set()
    for declaration in declarations:
        lanes.update(declaration.get("lanes") or [])
        if declaration.get("lane"):
            lanes.add(declaration["lane"])
    carriers = sorted(
        {c for d in declarations for c in (d.get("lipid_carriers") or [])}
    )

    rows: List[Dict[str, str]] = []

    # -- always: the B4 reliability band the absolutes are drawn with --------
    from src.kinetic_core.parameters_matrix import (
        HS_SPME_SAME_SAMPLE_DISPERSION,
        K_AW_UNCERTAINTY_DECADES,
    )

    rows.append(
        {
            "name": "HS-SPME same-sample dispersion",
            "value": f"{HS_SPME_SAME_SAMPLE_DISPERSION[0]:g}x&ndash;"
                     f"{HS_SPME_SAME_SAMPLE_DISPERSION[1]:g}x",
            "band": f"&plusmn;{math.log10(math.sqrt(HS_SPME_SAME_SAMPLE_DISPERSION[1])):.3f} decades",
            "class": "measured",
            "why": (
                "Two papers measuring the SAME samples disagree by this much. "
                "It is a calibration fact, not a fitted error, and it sets the "
                "floor width of every interval on this page -- and the "
                "band inside which a ratio is reported NOT RESOLVED."
            ),
        }
    )
    rows.append(
        {
            "name": "air/water partition constant K_aw",
            "value": "carried, not chosen",
            "band": f"&plusmn;{K_AW_UNCERTAINTY_DECADES:g} decades",
            "class": "declared band",
            "why": (
                "The literature spread on hexanal K_aw alone is 9.5x. The "
                "ruling is to CARRY the band rather than pick a constant, so "
                "it is added in quadrature to every absolute that passed "
                "through a partition step."
            ),
        }
    )

    # -- the trunk's one pinned split ---------------------------------------
    if lanes & {"trunk", "sulfur", "acrylamide"}:
        from src.kinetic_core.parameters import SCHIFF_AMADORI_SPLIT

        rows.append(
            {
                "name": "Schiff / Amadori split",
                "value": f"ratio {float(SCHIFF_AMADORI_SPLIT['ratio_amadori_over_schiff_pseudo_first_order']):.3g} "
                         f"at {float(SCHIFF_AMADORI_SPLIT['amine_loading_mmol_L_for_the_ratio']):g} mM amine",
                "band": "no band &mdash; PINNED",
                "class": str(SCHIFF_AMADORI_SPLIT.get("evidence_class", "derived_from_fit_data")),
                "why": (
                    "Martins refuses to split the condensation from the "
                    "rearrangement, so the pair is a composite and the "
                    "Amadori rate is DERIVED from the condensation by a pinned "
                    "ratio. It is the only derived rate in the trunk."
                ),
            }
        )

    # -- the network pH the trunk/acrylamide lanes actually sit at -----------
    from src.kinetic_core.parameters import AW_OF_MEASUREMENT, NETWORK_PH

    if lanes & {"trunk", "acrylamide"}:
        rows.append(
            {
                "name": "network pH (trunk / acrylamide)",
                "value": f"pH {NETWORK_PH:g}, fixed",
                "band": "no band &mdash; there is NO pH term",
                "class": "structural",
                "why": (
                    "Both lanes are homogeneous at the pH their parameters were "
                    "measured at. A pH you supply is RECORDED AND IGNORED; it "
                    "changes no rate. Five independent sign-crossings in the "
                    "corpus are why no family-level pH term was fitted."
                ),
            }
        )
        rows.append(
            {
                "name": "water activity a_w",
                "value": f"{AW_OF_MEASUREMENT:g} (dilute aqueous)",
                "band": "no band &mdash; there is NO a_w term",
                "class": "structural",
                "why": (
                    "Nothing in the fit corpus varies a_w, so an a_w dependence "
                    "would be invented. Your a_w is metadata on these lanes."
                ),
            }
        )

    # -- the sulfur lane's buffer + drift -----------------------------------
    if "sulfur" in lanes:
        buffers = [dict(m.get("buffer") or {}) for m in metadatas]
        kind = next((b.get("kind") for b in buffers if b), None)
        declared = any(b.get("declared_by_source") for b in buffers)
        rows.append(
            {
                "name": "buffer (sulfur lane pH trajectory)",
                "value": f"{esc(kind or 'none')}"
                         + ("" if declared else " &mdash; DEFAULTED, not declared"),
                "band": "the whole pH trajectory",
                "class": "measured" if declared else "declared assumption",
                "why": (
                    "The pH is an OUTPUT of this lane, not only an input. With "
                    "no buffer declared, the trajectory is computed from water "
                    "autoprotolysis and the charged solutes alone -- so if "
                    "your experiment WAS buffered, every pH-dependent rate in "
                    "this run drifts too far."
                ),
            }
        )
        drift = next((dict(m.get("ph_drift_constants") or {}) for m in metadatas if m.get("ph_drift_constants")), {})
        if drift:
            rows.append(
                {
                    "name": "pH-drift constants (frozen)",
                    "value": "acid yield "
                             f"{float(drift.get('acid_yield_per_sink_event', 0.0)):.4g} / sink event; "
                             f"ARP pKa {float(drift.get('arp_secondary_ammonium_pKa', 0.0)):.4g}",
                    "band": "frozen at the B2 fit; not re-fitted here",
                    "class": "derived_from_fit_data",
                    "why": (
                        "Two constants, fitted once and frozen. The engine "
                        "never constructs its own; it reads the fit report."
                    ),
                }
            )
        from src.kinetic_core.parameters_sulfur import NO_MEASURED_RATE

        for key, reason in NO_MEASURED_RATE.items():
            rows.append(
                {
                    "name": f"{key} &mdash; PINNED AT ZERO",
                    "value": "0.0",
                    "band": "none: the channel is switched off",
                    "class": "pinned (hold-out protected)",
                    "why": reason,
                }
            )

    # -- the lipid lane's three declared assumptions -------------------------
    if "lipid" in lanes:
        from src.kinetic_core.parameters_lipid import (
            K_LOOH_DECOMP_ANCHOR,
            LIPID_CARRIERS,
            Q10_ASSUMPTION,
        )

        rows.append(
            {
                "name": "Q10 of hydroperoxide decomposition",
                "value": f"{Q10_ASSUMPTION.default:.3f} per 10 &deg;C",
                "band": f"{Q10_ASSUMPTION.lo:g}&ndash;{Q10_ASSUMPTION.hi:g}, "
                        f"licensed only over "
                        f"{Q10_ASSUMPTION.licensed_span_c[0]:g}&ndash;"
                        f"{Q10_ASSUMPTION.licensed_span_c[1]:g} &deg;C",
                "class": "declared assumption",
                "why": Q10_ASSUMPTION.warning,
            }
        )
        rows.append(
            {
                "name": "LOOH decomposition rate anchor",
                "value": f"{float(K_LOOH_DECOMP_ANCHOR.value):g} /h at "
                         f"{K_LOOH_DECOMP_ANCHOR.temperature_of_measurement_c:g} &deg;C, "
                         f"pH {K_LOOH_DECOMP_ANCHOR.ph_of_measurement:g}",
                "band": "no standard error is printed anywhere in the source",
                "class": str(getattr(K_LOOH_DECOMP_ANCHOR, "evidence_class", "declared")),
                "why": (
                    "Hand-fitted by visual agreement in a rapeseed O/W emulsion, "
                    "and an explicit LUMP over all secondary products. The "
                    "branch DISTRIBUTION is measured; this RATE is not."
                ),
            }
        )
        for carrier_key in carriers:
            carrier = LIPID_CARRIERS.get(carrier_key)
            if carrier is None:
                continue
            rows.append(
                {
                    "name": f"lipid fraction &mdash; {esc(carrier.display)}",
                    "value": f"{carrier.lipid_mass_fraction:g} mass fraction",
                    "band": f"{carrier.lipid_lo:g}&ndash;{carrier.lipid_hi:g}",
                    "class": str(carrier.evidence_class),
                    "why": carrier.source_anchor,
                }
            )
            rows.append(
                {
                    "name": f"peroxide value &mdash; {esc(carrier.display)}",
                    "value": f"{carrier.peroxide_value_meq_per_kg:g} meq/kg fat",
                    "band": f"{carrier.pv_lo:g}&ndash;{carrier.pv_hi:g} meq/kg",
                    "class": str(carrier.evidence_class),
                    "why": (
                        (carrier.note or "")
                        + " MEASURE THIS AND THE LARGEST BAND IN THIS REPORT "
                        "COLLAPSES: peroxide value is a standard food-analysis "
                        "assay and the model consumes it directly."
                    ),
                }
            )

    # -- the matrix loading, if the run is not aqueous -----------------------
    from src.kinetic_core.parameters_matrix import MATRIX_LOADING

    for metadata in metadatas:
        matrix = str(metadata.get("matrix") or "")
        loading = MATRIX_LOADING.get(matrix)
        if loading is None or matrix in {r["name"] for r in rows}:
            continue
        rows.append(
            {
                "name": f"protein loading &mdash; {esc(matrix)}",
                "value": f"{(loading.protein_g_per_l or 0.0):g} g/L",
                "band": f"{(loading.protein_lo_g_per_l or 0.0):g}&ndash;"
                        f"{(loading.protein_hi_g_per_l or 0.0):g} g/L",
                "class": str(loading.evidence_class),
                "why": (loading.notes or "").strip() or "see parameters_matrix.MATRIX_LOADING",
            }
        )

    return rows


# ---------------------------------------------------------------------------
# Residual decomposition
# ---------------------------------------------------------------------------


def residual_section(
    compounds: Sequence[str],
    declaration: Mapping[str, Any],
    metadata: Mapping[str, Any],
    measured_ratios: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """
    The B4 residual decomposition for this run's matrix, per compound.

    Honest by construction: the UNEXPLAINED share exists only when a MEASURED
    matrix/water ratio is supplied (spec key ``measured_matrix_ratios``). With
    no measurement there is nothing to be unexplained BY, and the section says
    so instead of drawing a zero residual and implying the model closed the gap.
    """
    from src.kinetic_core.engine import resolve_matrix
    from src.kinetic_core.matrix_oav import decompose_residual, predict_matrix_shift
    from src.kinetic_core.parameters_matrix import MATRIX_LOADING

    matrix = resolve_matrix(str(metadata.get("matrix") or "water"))
    ph = metadata.get("ph")
    measured = {str(k): float(v) for k, v in dict(measured_ratios or {}).items()}

    if matrix == "water":
        return {
            "state": "aqueous",
            "matrix": matrix,
            "rows": [],
            "note": (
                "This run is AQUEOUS, so there is no matrix/water shift to "
                "decompose. The residual decomposition answers &lsquo;how much "
                "of the shift between water and a real food matrix do the "
                "model&rsquo;s named terms explain?&rsquo; &mdash; run the same "
                "formulation with a matrix descriptor that has a loading record "
                f"({', '.join(sorted(k for k in MATRIX_LOADING if k != 'water'))}) "
                "to see it."
            ),
        }
    if matrix not in MATRIX_LOADING:
        return {
            "state": "unknown_matrix",
            "matrix": matrix,
            "rows": [],
            "note": (
                f"No protein-loading record exists for the matrix "
                f"{matrix!r}. The B4 layer refuses to decompose a shift it "
                f"cannot size the matrix for; matrices with a record are: "
                f"{', '.join(sorted(MATRIX_LOADING))}."
            ),
        }

    rows: List[Dict[str, Any]] = []
    for compound in compounds:
        b4 = _b4_key_for(_species_key_for(compound, declaration))
        if b4 is None:
            continue
        prediction = predict_matrix_shift(b4, matrix, ph=ph)
        if prediction.state in ("unknown_compound", "unknown_matrix"):
            rows.append(
                {
                    "compound": compound,
                    "state": prediction.state,
                    "segments": [],
                    "warnings": list(prediction.warnings),
                    "measured_ratio": measured.get(compound),
                }
            )
            continue
        terms = dict(prediction.terms)
        segments: List[Tuple[str, float, str]] = []
        for name, css in (
            ("reversible_binding", "seg-bind"),
            ("alpha_beta_unsaturation", "seg-unsat"),
            ("covalent_ceiling", "seg-cov"),
        ):
            term = terms.get(name) or {}
            factor = term.get("factor_x") if isinstance(term, Mapping) else None
            decades = (
                math.log10(float(factor))
                if isinstance(factor, (int, float)) and float(factor) > 0
                else 0.0
            )
            segments.append((name, decades, css))
        row: Dict[str, Any] = {
            "compound": compound,
            "b4_key": b4,
            "state": prediction.state,
            "predicted_ratio": prediction.predicted_ratio,
            "predicted_interval": [prediction.predicted_lo, prediction.predicted_hi],
            "segments": segments,
            "terms": terms,
            "warnings": list(prediction.warnings),
            "measured_ratio": measured.get(compound),
            "flags": [],
        }
        if compound in measured:
            decomposition = decompose_residual(prediction, measured[compound])
            row["decomposition"] = decomposition.as_dict()
            row["flags"] = list(decomposition.flags)
            segments.append(
                ("UNEXPLAINED RESIDUAL", decomposition.residual_decades, "seg-resid")
            )
        rows.append(row)

    return {
        "state": "decomposed",
        "matrix": matrix,
        "rows": rows,
        "note": (
            "Explained shares are the model&rsquo;s three named terms. The "
            "UNEXPLAINED RESIDUAL bar appears only for a compound whose "
            "measured matrix/water ratio you supplied under "
            "<code>measured_matrix_ratios:</code> in the spec &mdash; without a "
            "measurement there is nothing for the model to be unexplained by, "
            "and drawing a zero residual would be a claim, not a gap."
        ),
    }


# ---------------------------------------------------------------------------
# The page
# ---------------------------------------------------------------------------

CSS = """
:root{
  --ink:#1b1a17; --muted:#5f5b53; --faint:#8b857c;
  --rule:#ddd8cf; --paper:#fbfaf7; --card:#ffffff;
  --accent:#7a3e12; --accent2:#2e5d4b;
  --warn:#8a5a00; --warnbg:#fdf4e0; --warnrule:#e8c887;
  --refuse:#8f2b1e; --refusebg:#fdeeeb; --refuserule:#e6b0a6;
  --ok:#2e5d4b; --okbg:#eef5f1;
  --unres:#6b6660; --unresbg:#f0eeea;
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
  font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
.wrap{max-width:1000px;margin:0 auto;padding:28px 22px 64px}
h1{font-size:25px;margin:0 0 4px;letter-spacing:-.01em}
h2{font-size:17px;margin:34px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--rule)}
h3{font-size:14px;margin:18px 0 6px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
p{margin:8px 0}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.9em}
.sub{color:var(--muted);margin:0 0 16px}
.badge{display:inline-block;padding:2px 9px;border-radius:11px;font-size:11.5px;
  font-weight:700;letter-spacing:.05em;text-transform:uppercase;vertical-align:2px}
.b-ok{background:var(--okbg);color:var(--ok);border:1px solid #b9d3c7}
.b-warn{background:var(--warnbg);color:var(--warn);border:1px solid var(--warnrule)}
.b-refuse{background:var(--refusebg);color:var(--refuse);border:1px solid var(--refuserule)}
.b-lane{background:#eee9e0;color:var(--accent);border:1px solid var(--rule)}
.b-unres{background:var(--unresbg);color:var(--unres);border:1px solid #cfcac2}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media (max-width:760px){.grid2{grid-template-columns:1fr}}
.card{background:var(--card);border:1px solid var(--rule);border-radius:7px;padding:14px 16px;margin:12px 0}
.card h4{margin:0 0 6px;font-size:14.5px}
.card.refusal{border-left:5px solid var(--refuse);background:var(--refusebg)}
.card.warn{border-left:5px solid var(--warnrule);background:var(--warnbg)}
.card.assume{border-left:5px solid var(--accent2)}
.kv{width:100%;border-collapse:collapse;font-size:13.5px}
.kv td{padding:3px 8px 3px 0;vertical-align:top;border:0}
.kv td:first-child{color:var(--muted);white-space:nowrap;width:1%}
table.data{width:100%;border-collapse:collapse;font-size:13.5px;margin:6px 0 4px}
table.data th{text-align:left;font-size:11.5px;text-transform:uppercase;letter-spacing:.05em;
  color:var(--muted);border-bottom:2px solid var(--rule);padding:5px 8px 5px 0;font-weight:700}
table.data td{padding:6px 8px 6px 0;border-bottom:1px solid var(--rule);vertical-align:top}
table.data td.num,table.data th.num{text-align:right;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
tr.notresolved td{background:var(--unresbg);color:var(--unres);font-style:italic}
tr.notresolved td.num{text-decoration:line-through;text-decoration-thickness:1px}
tr.dimer td{background:#f6f1ff}
.note{color:var(--muted);font-size:13px}
.tiny{color:var(--faint);font-size:12px}
.rule{color:var(--muted);font-size:13.5px;border-left:3px solid var(--accent);padding-left:10px;margin:10px 0}
svg.chart{display:block;margin:8px 0 4px;max-width:100%}
.chart .plotbg{fill:#fff;stroke:var(--rule)}
.chart .grid{stroke:#eceae5;stroke-width:1}
.chart .refline{stroke:var(--accent);stroke-width:1.2;stroke-dasharray:4 3}
.chart .reflabel{fill:var(--accent);font-size:10.5px}
.chart .bar{fill:#8fa8bd}
.chart .bar.hi{fill:#7a5cae}
.chart .whisker{stroke:#3d4a55;stroke-width:1.4}
.chart .rowlabel{fill:var(--ink);font-size:12.5px}
.chart .rowlabel.hi{fill:#5b3f8c;font-weight:700}
.chart .rowtag{fill:var(--faint);font-size:10px}
.chart .tick{fill:var(--muted);font-size:10.5px}
.chart .value{fill:var(--ink);font-size:11.5px;font-family:ui-monospace,Menlo,monospace}
.chart .zero{fill:var(--faint);font-size:11.5px}
.chart .axistitle{fill:var(--muted);font-size:11.5px}
.chart .seg{stroke:#fff;stroke-width:.8}
.chart .seg-bind{fill:#8fa8bd}
.chart .seg-unsat{fill:#c9a227}
.chart .seg-cov{fill:#cfcac2}
.chart .seg-resid{fill:#b45a4a}
.legend{font-size:12px;color:var(--muted);margin:2px 0 10px}
.legend span{display:inline-block;margin-right:14px}
.swatch{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:-1px;margin-right:4px}
footer{margin-top:44px;border-top:2px solid var(--rule);padding-top:16px}
details{margin:8px 0}
summary{cursor:pointer;color:var(--accent);font-size:13.5px}
@media print{
  body{background:#fff}
  .wrap{max-width:none;padding:0}
  .card,table.data,svg.chart{break-inside:avoid;page-break-inside:avoid}
  h2{break-after:avoid}
  details{display:block}
  details>summary{display:none}
  a[href]:after{content:""}
}
"""


def _page(title: str, body: str) -> str:
    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{esc(title)}</title>"
        f"<style>{CSS}</style>"
        "</head><body><div class=\"wrap\">"
        + body
        + "</div></body></html>\n"
    )


def _declaration_badge(state: str) -> str:
    if state == "in_envelope":
        return '<span class="badge b-ok">in envelope</span>'
    if state == "in_envelope_extrapolated":
        return '<span class="badge b-warn">extrapolated &mdash; answered with warnings</span>'
    return '<span class="badge b-refuse">out of envelope &mdash; refused</span>'


def _spec_table(spec: Mapping[str, Any]) -> str:
    precursors = dict(spec.get("precursors") or {})
    charge = ", ".join(
        f"{esc(name)} {float(value):g} mM" for name, value in sorted(precursors.items())
    )
    rows = [
        ("charge", charge or "&mdash;"),
        (
            "process",
            f"{float(spec.get('temp_C', 0)):g} &deg;C for "
            f"{float(spec.get('time_min', 0)):g} min",
        ),
        ("pH", f"{float(spec.get('ph', 0)):g}"),
        ("a_w", f"{float(spec.get('aw', 0)):g}" if spec.get("aw") is not None else "&mdash;"),
        ("matrix", esc(spec.get("matrix") or spec.get("protein_type") or "water")),
    ]
    body = "".join(f"<tr><td>{label}</td><td>{value}</td></tr>" for label, value in rows)
    return f'<table class="kv">{body}</table>'


def _refusal_cards(declaration: Mapping[str, Any]) -> str:
    """A refusal is an answer. Every named reason gets its own card."""
    parts: List[str] = []
    for reason in declaration.get("reasons") or []:
        parts.append(
            f'<div class="card refusal"><h4>Refused</h4><p>{esc(reason)}</p></div>'
        )
    for entry in declaration.get("unrepresented_targets") or []:
        if isinstance(entry, Mapping):
            compound, reason = entry.get("compound"), entry.get("reason")
        else:  # pragma: no cover - tuple form, defensive
            compound, reason = entry
        parts.append(
            f'<div class="card refusal"><h4>Cannot name: <code>{esc(compound)}</code>'
            f"</h4><p>{esc(reason)}</p></div>"
        )
    unmapped = declaration.get("unmapped_precursors") or []
    if unmapped:
        parts.append(
            '<div class="card refusal"><h4>Unmapped precursors</h4><p>'
            + ", ".join(f"<code>{esc(u)}</code>" for u in unmapped)
            + " &mdash; not a species in any core lane. The core is a NAMED "
            "small-molecule network; an intact protein, an isolate or a flour "
            "is not a precursor it can charge.</p></div>"
        )
    return "".join(parts)


def _warning_cards(declaration: Mapping[str, Any]) -> str:
    warnings = list(declaration.get("warnings") or [])
    if not warnings:
        return ""
    items = "".join(f"<div class=\"card warn\"><p>{esc(w)}</p></div>" for w in warnings)
    return (
        f'<h3>Declared extrapolations ({len(warnings)})</h3>'
        "<p class=\"note\">The engine answered, and these are the reasons it says "
        "the answer is being extrapolated. They are not footnotes.</p>" + items
    )


def _declaration_block(declaration: Mapping[str, Any], *, label: str = "") -> str:
    state = str(declaration.get("state", "unknown"))
    lanes = declaration.get("lanes") or ([declaration.get("lane")] if declaration.get("lane") else [])
    lane_badges = " ".join(
        f'<span class="badge b-lane">{esc(lane)} lane</span>' for lane in lanes if lane
    )
    mapped = dict(declaration.get("mapped_precursors") or {})
    mapped_text = (
        ", ".join(f"<code>{esc(k)}</code> {float(v):g} mM" for k, v in sorted(mapped.items()))
        or "&mdash;"
    )
    heading = f"<h3>{esc(label)} envelope</h3>" if label else "<h3>Envelope</h3>"
    return (
        heading
        + f"<p>{_declaration_badge(state)} {lane_badges}</p>"
        + f'<table class="kv"><tr><td>mapped precursors</td><td>{mapped_text}</td></tr>'
        + f'<tr><td>mapped targets</td><td>'
        + (
            ", ".join(
                f"{esc(c)} &rarr; <code>{esc(k)}</code>"
                for c, k in sorted(dict(declaration.get("mapped_targets") or {}).items())
            )
            or "&mdash;"
        )
        + "</td></tr>"
        + (
            f'<tr><td>lipid carriers</td><td>'
            + ", ".join(f"<code>{esc(c)}</code>" for c in declaration.get("lipid_carriers") or [])
            + "</td></tr>"
            if declaration.get("lipid_carriers")
            else ""
        )
        + "</table>"
        + _refusal_cards(declaration)
        + _warning_cards(declaration)
    )


def _assumptions_block(rows: Sequence[Mapping[str, str]]) -> str:
    if not rows:
        return ""
    cards = "".join(
        '<div class="card assume"><h4>{name} '
        '<span class="badge b-lane">{cls}</span></h4>'
        '<table class="kv"><tr><td>value</td><td>{value}</td></tr>'
        "<tr><td>band</td><td>{band}</td></tr></table>"
        '<p class="note">{why}</p></div>'.format(
            name=row["name"],
            cls=esc(row["class"]),
            value=row["value"],
            band=row["band"],
            why=esc(row["why"]),
        )
        for row in rows
    )
    return (
        f"<h2>Declared assumptions this run used ({len(rows)})</h2>"
        "<p class=\"note\">Not a general disclaimer list &mdash; each of these was "
        "read by THIS run, on the lanes it selected. Every one carries its band, "
        "because the band is the product.</p>" + cards
    )


def _plain(text: str) -> str:
    """Strip the handful of HTML entities the headings carry, for plain-text slots."""
    for entity, replacement in (
        ("&mdash;", "--"),
        ("&ndash;", "-"),
        ("&middot;", "."),
        ("&amp;", "&"),
        ("&times;", "x"),
        ("&deg;", "deg"),
    ):
        text = text.replace(entity, replacement)
    return text


def _footer(title: str) -> str:
    title = _plain(title)
    card = model_card_headline()
    sentences = "".join(
        f'<p class="note">{_markdownish(s)}</p>' for s in card.get("sentences", [])
    )
    if not card.get("available"):
        sentences = f'<p class="note">{esc(card.get("note", ""))}</p>'
    return (
        "<footer>"
        f"<h3>How to read this model</h3>"
        + (
            f'<p class="rule">{esc(card.get("one_liner", ""))}</p>'
            if card.get("one_liner")
            else ""
        )
        + sentences
        + '<table class="kv">'
        f"<tr><td>report</td><td>{esc(title)}</td></tr>"
        f"<tr><td>model version</td><td><code>{esc(_git_version())}</code></td></tr>"
        f"<tr><td>generated</td><td>{esc(datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'))}</td></tr>"
        "<tr><td>generated by</td><td><code>src/report_html.py</code> "
        "(Build Wave V1) via <code>scripts/maillard.py --report</code></td></tr>"
        "<tr><td>self-contained</td><td>no network requests, no external assets, "
        "no scripts &mdash; this file renders identically offline</td></tr>"
        "</table></footer>"
    )


# ---------------------------------------------------------------------------
# PREDICT
# ---------------------------------------------------------------------------


def render_predict_report(payload: Mapping[str, Any]) -> str:
    system = dict(payload.get("system") or {})
    spec = dict(system.get("spec") or {})
    name = str(system.get("name") or "system")
    declaration = dict(payload.get("declaration") or {})
    metadata = dict(payload.get("run_metadata") or {})
    oav_payload = dict(payload.get("oav_table") or {})
    rows = list(payload.get("rows") or [])
    answered = bool(payload.get("answered"))
    title = f"Maillard predict &mdash; {name}"

    body: List[str] = [
        f"<h1>{title}</h1>",
        f'<p class="sub">Kinetic core prediction &middot; '
        f"{esc(metadata.get('thermal_program') or '')} &middot; "
        f"matrix <code>{esc(metadata.get('matrix') or spec.get('matrix') or 'water')}</code></p>",
        '<div class="grid2">',
        f'<div class="card"><h4>Formulation &amp; process</h4>{_spec_table(spec)}</div>',
        f'<div class="card">{_declaration_block(declaration)}</div>',
        "</div>",
    ]

    if not answered:
        body.append(
            '<div class="card refusal"><h4>No number is emitted</h4>'
            "<p>The core declined this request for the reason(s) above. That is "
            "an output, not a failure: the alternative is a plausible-looking "
            "float with nothing behind it, and every documented accuracy defect "
            "in this repository began as a number that should not have "
            "existed.</p></div>"
        )
        body.append(
            f"<h2>What the core can answer instead</h2>{_vocabulary_block()}"
        )
        body.append(_assumptions_block(declared_assumptions([declaration], [metadata])))
        body.append(_footer(title))
        return _page(f"Maillard predict - {name} (refused)", "".join(body))

    # ---- rankings table ---------------------------------------------------
    per_species = dict(oav_payload.get("per_species") or {})
    table_rows: List[str] = []
    chart_rows: List[Dict[str, Any]] = []
    dimer_key = "MFTD"
    for row in rows:
        compound = str(row.get("compound"))
        species_key = str(row.get("species_key") or _species_key_for(compound, declaration))
        entry = _oav_entry(oav_payload, species_key)
        concentration = dict(entry.get("concentration") or {})
        # Q1: the row carries its own interval now. The fallback to the OAV
        # table is kept for payloads written before that, but it is a FALLBACK:
        # the table has no entry at all for a NO_B4_RECORD species, so reading
        # the interval out of it printed a bare point for four of the lipid
        # lane's seven products -- against this module's own rule that an
        # absolute is never drawn without its interval.
        interval = (
            row.get("interval_ug_per_l")
            or concentration.get("interval_ug_per_L")
            or [None, None]
        )
        band_x = row.get("band_x", concentration.get("band_x"))
        oav_summary = dict(row.get("oav") or {})
        oav_interval = entry.get("OAV_interval") or [None, None]
        is_dimer = species_key == dimer_key

        if oav_summary.get("available"):
            oav_cell = (
                f"{sig(oav_summary.get('oav'))}"
                f'<div class="tiny">{sig(oav_interval[0])} &ndash; {sig(oav_interval[1])}</div>'
            )
            threshold_cell = (
                f"{sig(oav_summary.get('threshold_ug_per_l'))}"
                f'<div class="tiny">{esc(oav_summary.get("threshold_source"))}</div>'
            )
            chart_rows.append(
                {
                    "label": compound,
                    "point": oav_summary.get("oav"),
                    "lo": oav_interval[0],
                    "hi": oav_interval[1],
                    "highlight": is_dimer,
                    "tag": "potency-weighted dimer" if is_dimer else None,
                }
            )
        else:
            oav_cell = (
                '<span class="badge b-unres">no threshold</span>'
                f'<div class="tiny">{esc(oav_summary.get("reason"))}</div>'
            )
            threshold_cell = '<span class="tiny">not measured in this matrix</span>'

        table_rows.append(
            '<tr class="{cls}"><td>{compound}<div class="tiny"><code>{key}</code> '
            "&middot; {lane} lane</div></td>"
            '<td class="num">{point}<div class="tiny">{lo} &ndash; {hi}</div></td>'
            '<td class="num">{band}</td><td class="num">{oav}</td><td>{thr}</td></tr>'.format(
                cls="dimer" if is_dimer else "",
                compound=esc(compound),
                key=esc(species_key),
                lane=esc(row.get("lane")),
                point=sig(row.get("predicted_ug_per_l")),
                lo=sig(interval[0]),
                hi=sig(interval[1]),
                band=(
                    sig(band_x) + "&times;"
                    if float(row.get("predicted_ug_per_l") or 0.0) > 0.0
                    else "n/a"
                ),
                oav=oav_cell,
                thr=threshold_cell,
            )
        )

    body.append("<h2>Ranking</h2>")
    body.append(f'<p class="rule">{esc(INTERVAL_RULE)}</p>')
    body.append(
        '<table class="data"><thead><tr>'
        "<th>compound</th>"
        '<th class="num">&micro;g/L (= ppb in water)<br>interval</th>'
        '<th class="num">band</th><th class="num">OAV<br>interval</th>'
        "<th>odour threshold &micro;g/L &amp; source</th>"
        "</tr></thead><tbody>" + "".join(table_rows) + "</tbody></table>"
    )

    if dimer_key in per_species:
        ratio = oav_payload.get("dimer_potency_ratio_over_monomer")
        body.append(
            '<div class="card"><h4>The MFT dimer is plotted at its own potency</h4>'
            f"<p>bis(2-methyl-3-furyl) disulfide is <strong>{sig(ratio)}&times; more "
            "potent</strong> than its own monomer, so it is charted at its "
            "potency-weighted value (its OAV), never at its mass share. "
            "<strong>Mass lost to dimerisation is not aroma lost</strong>, and any "
            "objective scoring the dimerisation channel as a pure loss is wrong "
            "by roughly the threshold ratio.</p>"
            + (
                f'<p class="note">dimer OAV / monomer OAV in this run: '
                f'{sig(oav_payload.get("dimer_over_monomer_OAV"))}</p>'
                if oav_payload.get("dimer_over_monomer_OAV") is not None
                else ""
            )
            + "</div>"
        )

    body.append("<h2>Odour activity (log scale, measured intervals)</h2>")
    body.append(
        svg_log_bars(
            chart_rows,
            x_title="OAV = concentration / odour threshold (log scale)",
            reference=1.0,
            reference_label="OAV = 1 (threshold)",
        )
    )
    body.append(
        '<p class="legend"><span><span class="swatch" style="background:#8fa8bd"></span>'
        "OAV point</span><span><span class=\"swatch\" style=\"background:#7a5cae\"></span>"
        "MFT dimer, at its potency-weighted value</span>"
        "<span>&#9500;&#9472;&#9508; measured interval</span></p>"
    )
    if oav_payload.get("n_without_threshold"):
        body.append(
            f'<p class="note">{int(oav_payload["n_without_threshold"])} compound(s) '
            "carry no measured threshold in this matrix and are therefore absent "
            "from the chart but present in the table. Nothing is borrowed from "
            "another matrix.</p>"
        )

    # ---- residual decomposition -------------------------------------------
    residual = residual_section(
        [str(r.get("compound")) for r in rows],
        declaration,
        metadata,
        measured_ratios=spec.get("measured_matrix_ratios"),
    )
    body.append("<h2>Matrix-shift residual decomposition</h2>")
    body.append(f'<p class="note">{residual["note"]}</p>')
    if residual["rows"]:
        body.append(
            svg_stacked_decades(
                [
                    {"label": r["compound"], "segments": r["segments"]}
                    for r in residual["rows"]
                    if r.get("segments")
                ]
            )
        )
        body.append(
            '<p class="legend">'
            '<span><span class="swatch" style="background:#8fa8bd"></span>reversible binding</span>'
            '<span><span class="swatch" style="background:#c9a227"></span>&alpha;,&beta;-unsaturation</span>'
            '<span><span class="swatch" style="background:#cfcac2"></span>covalent ceiling (INERT, contributes 0)</span>'
            '<span><span class="swatch" style="background:#b45a4a"></span>UNEXPLAINED residual</span>'
            "</p>"
        )
        for row in residual["rows"]:
            if row.get("flags") or row.get("warnings"):
                body.append(
                    f'<div class="card warn"><h4>{esc(row["compound"])}</h4>'
                    + "".join(f"<p>{esc(f)}</p>" for f in row.get("flags", []))
                    + "".join(f'<p class="note">{esc(w)}</p>' for w in row.get("warnings", []))
                    + "</div>"
                )

    body.append(_assumptions_block(declared_assumptions([declaration], [metadata])))
    body.append(
        '<div class="card"><h4>Caveat, in full</h4><p class="note">'
        + esc(dict(payload.get("caveats") or {}).get("core", ""))
        + "</p></div>"
    )
    body.append(_footer(title))
    return _page(f"Maillard predict - {name}", "".join(body))


def _vocabulary_block() -> str:
    """What the engine CAN be asked, printed next to every refusal."""
    from src.kinetic_core.engine import (
        LANE_DEFAULT_TARGETS,
        PRECURSOR_ALIASES,
        UNREPRESENTED_COMPOUNDS,
    )

    lanes = "".join(
        f"<tr><td><code>{esc(lane)}</code></td><td>"
        + ", ".join(esc(t) for t in targets)
        + "</td></tr>"
        for lane, targets in LANE_DEFAULT_TARGETS.items()
    )
    precursors = ", ".join(f"<code>{esc(p)}</code>" for p in sorted(set(PRECURSOR_ALIASES)))
    refused = "".join(
        f"<tr><td><code>{esc(compound)}</code></td><td>{esc(reason)}</td></tr>"
        for compound, reason in sorted(UNREPRESENTED_COMPOUNDS.items())
    )
    return (
        '<div class="card"><h4>Targets each lane can report</h4>'
        f'<table class="data"><thead><tr><th>lane</th><th>targets</th></tr></thead>'
        f"<tbody>{lanes}</tbody></table></div>"
        f'<div class="card"><h4>Precursor names the engine understands</h4>'
        f'<p class="note">{precursors}</p></div>'
        '<details><summary>Compounds the core deliberately refuses, and why</summary>'
        f'<table class="data"><thead><tr><th>compound</th><th>declared reason</th></tr>'
        f"</thead><tbody>{refused}</tbody></table></details>"
    )


# ---------------------------------------------------------------------------
# COMPARE
# ---------------------------------------------------------------------------


def render_compare_report(payload: Mapping[str, Any]) -> str:
    a = dict(payload.get("a") or {})
    b = dict(payload.get("b") or {})
    comparison = dict(payload.get("comparison") or {})
    name_a, name_b = str(a.get("name") or "A"), str(b.get("name") or "B")
    title = f"Maillard compare &mdash; {name_a} vs {name_b}"
    declaration_a = dict(comparison.get("declaration_a") or {})
    declaration_b = dict(comparison.get("declaration_b") or {})

    body: List[str] = [
        f"<h1>{title}</h1>",
        '<p class="sub">Kinetic core comparison &middot; <strong>ratios lead</strong>: '
        "the HS-SPME calibration offset and the air/water partition constant are "
        "shared between the two arms and cancel exactly in a within-run ratio.</p>",
        '<div class="grid2">',
        f'<div class="card"><h4>A &mdash; {esc(name_a)}</h4>'
        f'{_spec_table(dict(a.get("spec") or {}))}'
        f"{_declaration_block(declaration_a, label='A')}</div>",
        f'<div class="card"><h4>B &mdash; {esc(name_b)}</h4>'
        f'{_spec_table(dict(b.get("spec") or {}))}'
        f"{_declaration_block(declaration_b, label='B')}</div>",
        "</div>",
    ]

    if not comparison.get("comparable"):
        body.append(
            '<div class="card refusal"><h4>Not comparable</h4><p>'
            + esc(comparison.get("reason", ""))
            + "</p><p>A ratio against a refusal is not a ratio. The refusals "
            "above name what is missing; fix one of them and the comparison "
            "becomes answerable.</p></div>"
        )
        body.append(f"<h2>What the core can answer instead</h2>{_vocabulary_block()}")
        body.append(
            _assumptions_block(
                declared_assumptions(
                    [declaration_a, declaration_b],
                    [
                        dict((comparison.get("run_a") or {}).get("run_metadata") or {}),
                        dict((comparison.get("run_b") or {}).get("run_metadata") or {}),
                    ],
                )
            )
        )
        body.append(_footer(title))
        return _page(f"Maillard compare - {name_a} vs {name_b} (refused)", "".join(body))

    ratios = dict(comparison.get("ratios") or {})
    rows = list(ratios.get("rows") or [])
    band = ratios.get("reliability_band_x")
    n_resolved = int(ratios.get("n_resolved", 0))
    n_compared = int(ratios.get("n_compared", 0))
    # Q1: published by the science layer, which now excludes an undefined ratio
    # from n_resolved rather than counting it. Read, never re-derived.
    n_undefined = int(ratios.get("n_undefined", 0))

    body.append("<h2>Per-compound ratios (A / B)</h2>")
    body.append(
        f'<p class="rule"><strong>{n_resolved} of {n_compared}</strong> ratios resolve '
        f"above the same-sample dispersion band"
        + (f" ({sig(band)}&times;)" if isinstance(band, (int, float)) else "")
        + ". A ratio INSIDE that band is reported <strong>NOT RESOLVED</strong> &mdash; "
        "not as a small effect. The two are different claims and are styled "
        "differently below."
        + (
            f" {n_undefined} of the {n_compared} is not a ratio at all "
            "(one arm at exactly zero) and resolves nothing."
            if n_undefined else ""
        )
        + "</p>"
    )

    table_rows: List[str] = []
    for row in rows:
        ratio = row.get("ratio_a_over_b")
        unresolved = bool(row.get("within_reliability_band"))
        undefined = (
            not isinstance(ratio, (int, float))
            or ratio != ratio
            or ratio in (float("inf"), float("-inf"))
            or ratio == 0.0
        )
        if not isinstance(ratio, (int, float)) or ratio != ratio:
            ratio_cell = "&mdash;"
        elif ratio in (float("inf"), float("-inf")):
            ratio_cell = "A only (B is zero)"
        elif ratio == 0.0:
            ratio_cell = "B only (A is zero)"
        else:
            ratio_cell = f"{ratio:.3g}&times;"
        # An UNDEFINED ratio is a third state: it is neither resolved nor
        # unresolved, because there is no ratio to resolve. Badging it
        # "resolved" would read as a claim about a quantity that does not exist.
        if undefined:
            state_cell = '<span class="badge b-warn">undefined</span>'
        elif unresolved:
            state_cell = '<span class="badge b-unres">not resolved</span>'
        else:
            state_cell = '<span class="badge b-ok">resolved</span>'
        table_rows.append(
            '<tr class="{cls}"><td>{compound}</td><td class="num">{ratio}</td>'
            "<td>{direction}</td><td>{state}</td><td class=\"note\">{note}</td></tr>".format(
                cls="notresolved" if unresolved else "",
                compound=esc(row.get("compound")),
                ratio=ratio_cell,
                direction=esc(str(row.get("direction", "-")).replace("_", " ")),
                state=state_cell,
                note=esc(row.get("note", "")),
            )
        )
    body.append(
        '<table class="data"><thead><tr><th>compound</th>'
        '<th class="num">A / B</th><th>direction</th><th>state</th><th>note</th>'
        "</tr></thead><tbody>" + "".join(table_rows) + "</tbody></table>"
    )

    resolved_chart = [
        {
            "label": str(row.get("compound")),
            "point": float(row["ratio_a_over_b"]),
            "lo": None,
            "hi": None,
            "highlight": False,
            "tag": None if not row.get("within_reliability_band") else "NOT RESOLVED",
        }
        for row in rows
        if isinstance(row.get("ratio_a_over_b"), (int, float))
        and 0.0 < float(row["ratio_a_over_b"]) < float("inf")
    ]
    if resolved_chart:
        body.append("<h3>Ratios on a log axis</h3>")
        body.append(
            svg_log_bars(
                resolved_chart,
                x_title="A / B ratio, log scale (line = 1.0, no difference)",
                reference=1.0,
                reference_label="A = B",
            )
        )
        if isinstance(band, (int, float)):
            body.append(
                f'<p class="note">Anything between {1.0 / float(band):.3g}&times; and '
                f"{float(band):.3g}&times; is inside the measured same-sample "
                "dispersion band and is <strong>not resolved</strong>, whatever "
                "the bar looks like.</p>"
            )

    # ---- per-arm OAV, AS THE ENGINE EMITTED IT ----------------------------
    # Q1: this block used to rebuild each arm's OAV table here, by hand, from
    # the run dict -- ``engine.compare`` returned the arms flattened through
    # ``as_dict()``, which drops the object and with it ``.oav()``. The copy had
    # already drifted (it never learned B7's furanone band, so a compare drew
    # narrower intervals than a predict of the same arm), which is exactly the
    # failure its own docstring predicted. ``engine.compare`` now emits
    # ``oav_table_a``/``oav_table_b`` from the live objects and this module
    # reads them.
    run_a = dict(comparison.get("run_a") or {})
    run_b = dict(comparison.get("run_b") or {})
    for label, run, declaration, oav_payload in (
        (name_a, run_a, declaration_a, dict(comparison.get("oav_table_a") or {})),
        (name_b, run_b, declaration_b, dict(comparison.get("oav_table_b") or {})),
    ):
        if not run.get("concentrations_ug_per_l"):
            continue
        chart_rows = []
        for species_key in oav_payload.get("ranking_by_OAV") or []:
            entry = dict((oav_payload.get("per_species") or {}).get(species_key) or {})
            interval = entry.get("OAV_interval") or [None, None]
            chart_rows.append(
                {
                    "label": species_key,
                    "point": entry.get("OAV_point"),
                    "lo": interval[0],
                    "hi": interval[1],
                    "highlight": species_key == "MFTD",
                    "tag": "potency-weighted dimer" if species_key == "MFTD" else None,
                }
            )
        if chart_rows:
            body.append(f"<h3>Odour activity &mdash; {esc(label)}</h3>")
            body.append(
                svg_log_bars(
                    chart_rows,
                    x_title="OAV (log scale)",
                    reference=1.0,
                    reference_label="OAV = 1",
                )
            )

    body.append(
        _assumptions_block(
            declared_assumptions(
                [declaration_a, declaration_b],
                [
                    dict(run_a.get("run_metadata") or {}),
                    dict(run_b.get("run_metadata") or {}),
                ],
            )
        )
    )
    caveats = dict(payload.get("caveats") or {})
    body.append(
        '<div class="card"><h4>Caveats, in full</h4>'
        + "".join(
            f'<p class="note">{esc(text)}</p>' for text in caveats.values() if text
        )
        + "</div>"
    )
    body.append(_footer(title))
    return _page(f"Maillard compare - {name_a} vs {name_b}", "".join(body))


# ---------------------------------------------------------------------------


def write_report(payload: Mapping[str, Any], path: Path | str) -> Path:
    """
    Render the right report for this payload and write it.

    Q1 FIXED THE ROUTING PREDICATE. This dispatched on
    ``artifact.startswith("maillard_compare")``, and the FAST lane's artifact is
    named ``"maillard_compare"`` -- a prefix of the core's
    ``"maillard_compare_core"``. So the guard below could never fire for the
    payload it was written to catch: a FAST payload matched the prefix, was
    routed into the CORE renderer, and died on a missing key, while the
    friendly error explaining that ``--report`` is core-only sat unreachable.
    The dispatch is now on the exact core artifact names, so a FAST payload gets
    the error the author intended. No core payload changes route.
    """
    artifact = str(payload.get("artifact") or "")
    if artifact == "maillard_compare_core":
        text = render_compare_report(payload)
    elif artifact == "maillard_predict_core":
        text = render_predict_report(payload)
    else:
        raise ValueError(
            f"no HTML report exists for artifact {artifact!r}. --report is "
            "available on `predict` and `compare` on the core lane."
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return target


__all__ = [
    "declared_assumptions",
    "model_card_headline",
    "render_compare_report",
    "render_predict_report",
    "residual_section",
    "svg_log_bars",
    "svg_stacked_decades",
    "write_report",
]
