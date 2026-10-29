#!/usr/bin/env python
"""
build_sulfur_reaction_evidence.py -- what each sulfur-lane step rests on, and which steps are unclear (strategy M4).

Owner, 2026-10-29: "figure out if we have clarity on the reactions, and flag where not". For every step of the shipped
sulfur network this table joins five records that already exist, so it adds no number of its own:
  * the step's constant: evidence class (measured / fitted / bounded / derived), source, and the extraction dossiers
    its anchor cites (``src/evidence_anchors.py``; every anchor resolves, ``evidence_anchor_gate``);
  * whether the shipped fit identifies it (``kinetic_core_b60_laplace_covariance.json``: identified / on a bound /
    unidentified, per coordinate; "not free" when B60 did not fit it);
  * how much flux it carries: its largest share of a pot's TOTAL integrated reaction flux over the eight reference pots
    of ``sulfur_flux_attribution.json``, and (for context) its largest share of one species' production or consumption;
  * the direct measurements of its species in the measurement ledgers (``data/lit/measurement_ledger``): how many
    data sets measure a reactant or product, how many are time series, how many the model already uses;
  * the conflicts register entries that name its constant.
A step is UNCLEAR when it carries more than 1 % of a pot's total reaction flux and its constant is fitted but not identified, sits on a
bound, is bounded from a timescale bracket, or nothing in the corpus measures any of its species over time; and when
an open conflict names it. The unclear steps, ranked by flux, are what the sub-models must test first.

    python scripts/generators/build_sulfur_reaction_evidence.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, List

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, evidence_anchors, provenance  # noqa: E402
from src.kinetic_core import engine, sulfur  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "sulfur_reaction_evidence.json"
FLUX = data_paths.VALIDATION_DIR / "sulfur_flux_attribution.json"
LAPLACE = data_paths.VALIDATION_DIR / "kinetic_core_b60_laplace_covariance.json"
#: Where a constant B60 did not free was last fitted: B2.3, whose report states per constant whether it was identified.
B23_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b2_3_fit_report.json"
ENVELOPE = data_paths.VALIDATION_DIR / "core_prediction_uncertainty.json"
FLUX_THRESHOLD = 0.01
#: The odorants the lane predicts and their immediate precursors. A step can carry a negligible share of the bulk flux
#: and still decide every prediction (thiol formation and loss), so a step also counts when it carries more than
#: FLUX_THRESHOLD of one of these species' formation or loss in some pot.
ODORANT_PATH = ("MFT", "FFT", "MFTD", "FFTD", "H2S", "NF", "FUR")

GRADE = {"measured_rate": "measured", "measured_activation_energy": "measured", "derived_from_fit_data": "fitted",
         "bounded_from_a_timescale_bracket": "bounded", "structural_constant": "structural"}

#: How each engine species is named in the measurement ledgers' free-text species lists (lower-case patterns).
ALIASES: Dict[str, List[str]] = {
    "PENT": [r"ribose", r"xylose", r"pentose", r"arabinose"], "Cys": [r"\bcys\b", r"cysteine(?!.*amadori)"],
    "THI": [r"thiamin"], "H2S": [r"\bh2s\b", r"hydrogen sulfide", r"\bsulfide\b"], "GSH": [r"glutathione", r"\bgsh\b"],
    "ARP": [r"\barp\b", r"xylose.*amadori", r"amadori.*xyl"], "CAMA": [r"\bcama\b", r"cys.?amadori", r"cysteine amadori"],
    "TTCA": [r"\bttca\b", r"thiazolidine"], "DPO": [r"1-deoxypentosone", r"1-dp\b", r"\bdpo\b", r"1-deoxyosone"],
    "TDP": [r"3-deoxypentosone", r"3-deoxyribosulose", r"3-deoxyxylosone", r"\btdp\b", r"3-dp\b"],
    "DDP": [r"dideoxypentosone"], "NF": [r"norfuraneol", r"4-hydroxy-5-methyl-3\(2h\)-furanone"],
    "PTR": [r"pentanetrione"], "FUR": [r"furfural(?!.*thiol)", r"2-furaldehyde", r"furan-2-aldehyde"],
    "HA": [r"hydroxyacetaldehyde", r"glycolaldehyde"], "MP": [r"mercapto-2-propanone", r"mercaptopropanone"],
    "HMP": [r"5-hydroxy-3-mercapto"], "MESH": [r"methanethiol"], "MFT": [r"\bmft\b", r"2-methyl-3-furanthiol"],
    "FFT": [r"\bfft\b", r"furfurylthiol", r"furanmethanethiol"],
    "MFTD": [r"bis\(2-methyl-3-furyl\)", r"\bmftd\b", r"mft dimer"], "FFTD": [r"difurfuryl", r"bis\(2-furfuryl\)"],
    "MP3P": [r"2-mercapto-3-pentanone"], "MP2P": [r"3-mercapto-2-pentanone"], "ACTZ": [r"acetylthiazole"],
    "Glc": [r"glucose"], "Fru": [r"fructose"], "TDG": [r"3-deoxyglucosone", r"\b3-dg\b"], "MGO": [r"methylglyoxal", r"\bmgo\b"],
    "GO": [r"\bglyoxal\b", r"\bgo\b"], "DA": [r"diacetyl", r"butanedione"], "HMF": [r"\bhmf\b", r"hydroxymethylfurfural"],
    "DMHF": [r"furaneol", r"\bhdmf\b", r"\bdmhf\b"], "AMA": [r"glucose.*amadori", r"\bdfg\b"],
    "FA": [r"formic"], "AA": [r"acetic acid"], "OX": [r"\boxygen\b"],
}


def _ledger_entries() -> List[dict]:
    out = []
    for path in sorted(data_paths.MEASUREMENT_LEDGER_DIR.glob("*.yml")):
        for entry in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("entries") or []:
            entry["_ledger"] = path.stem
            entry["_species_text"] = " | ".join(str(s) for s in entry.get("species") or []).lower()
            out.append(entry)
    return out


def _measures(entry: dict, species_key: str) -> bool:
    return any(re.search(p, entry["_species_text"]) for p in ALIASES.get(species_key, []))


def _flux_shares() -> Dict[str, dict]:
    """reaction -> its largest share of a pot's TOTAL integrated reaction flux, and in which pot; plus its largest share
    of one tracked species' production or consumption (a step that is a species' only source scores 100 % there)."""
    best: Dict[str, dict] = {}
    for pot, record in json.loads(FLUX.read_text(encoding="utf-8"))["pots"].items():
        fluxes = record["reaction_flux_mmol_L"]
        total = sum(fluxes.values())
        for rxn, f in fluxes.items():
            entry = best.setdefault(rxn, {"share": 0.0, "pot": pot, "species_share": 0.0, "species": None,
                                          "odorant_share": 0.0, "odorant": None})
            if total > 0 and f / total > entry["share"]:
                entry.update(share=f / total, pot=pot)
        for species, sides in record["species"].items():
            for side in ("made", "used"):
                for item in sides[side]["top"]:
                    entry = best.setdefault(item["reaction"], {"share": 0.0, "pot": pot, "species_share": 0.0, "species": None,
                                                               "odorant_share": 0.0, "odorant": None})
                    if item["share"] > entry["species_share"]:
                        entry.update(species_share=item["share"], species=f"{species} {side}")
                    if species in ODORANT_PATH and item["share"] > entry["odorant_share"]:
                        entry.update(odorant_share=item["share"], odorant=f"{species} {side} ({pot.split('_')[0]})")
    return best


def _identification() -> Dict[str, Dict[str, str]]:
    """parameter -> {coordinate kind: verdict} from B60's Laplace; for a constant B60 kept frozen, B2.3's verdict, the
    last fit that freed it (key 'B2.3 log10k')."""
    lap = json.loads(LAPLACE.read_text(encoding="utf-8"))
    out: Dict[str, Dict[str, str]] = {}
    for coord, ident, bound in zip(lap["coordinates"], lap["identified"], lap["at_bound"]):
        out.setdefault(coord["key"], {})[coord["kind"]] = "on a bound" if bound else ("identified" if ident else "unidentified")
    for key, interval in json.loads(B23_REPORT.read_text(encoding="utf-8"))["parameter_intervals"].items():
        if key not in out and isinstance(interval, dict) and "identified" in interval:
            out[key] = {"B2.3 log10k": "identified" if interval["identified"] else "unidentified, frozen since"}
    return out


def _envelope_sampling() -> Dict[str, Dict[str, bool]]:
    """parameter -> {coordinate: sampled?} as the shipped envelope records it ('log10k' for the rate, 'Ea' for the
    barrier); prior keys read '<wave>.<parameter>.<coordinate>'."""
    out: Dict[str, Dict[str, bool]] = {}
    for prior in json.loads(ENVELOPE.read_text(encoding="utf-8"))["priors"]:
        parts = prior["key"].split(".")
        if len(parts) >= 3:
            kind = "Ea" if parts[-1].startswith("ea") else "log10k"
            out.setdefault(parts[1], {})[kind] = out.get(parts[1], {}).get(kind, False) or bool(prior["sampled"])
    return out


def _conflicts() -> List[dict]:
    return (yaml.safe_load(data_paths.CONFLICTS_REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []


def _constant_status(steps: List[dict], parameters, ident, sampling) -> Dict[str, dict]:
    """Per constant of a step that matters: how it is pinned (2026-10-29, owner: "How can we prevent this from happening
    in the future?"). Ten constants carrying most of the lane's flux were fitted unidentified in B2.3 and frozen since,
    the status lost on the way and their rates absent from the envelope; this record carries the status forward."""
    out: Dict[str, dict] = {}
    for st in steps:
        key = st["parameter"]
        if not key or key not in parameters:
            continue
        c = out.setdefault(key, {"grade": st["grade"], "identification": ident.get(key) or None,
                                 "envelope_samples": sampling.get(key, {}), "steps": [], "max_flux_share": 0.0,
                                 "max_odorant_path_share": 0.0, "matters": False})
        c["steps"].append(st["reaction"])
        c["max_flux_share"] = max(c["max_flux_share"], st["max_flux_share"])
        c["max_odorant_path_share"] = max(c["max_odorant_path_share"], st["max_odorant_path_share"])
        c["matters"] = c["matters"] or st["matters"]
    for c in out.values():
        verdicts = (c["identification"] or {}).values()
        c["unpinned"] = c["grade"] == "fitted" and (not verdicts or any(v != "identified" for v in verdicts))
        c["unpinned_and_matters"] = c["unpinned"] and c["matters"]
    return out


def build() -> dict:
    parameters = engine.core_parameters(engine.SULFUR)
    shares, ident, ledger, conflicts = _flux_shares(), _identification(), _ledger_entries(), _conflicts()
    sampling = _envelope_sampling()
    steps = []
    for rxn in sulfur.FULL_REACTIONS:
        p = parameters.get(rxn.parameter_key) if rxn.parameter_key else None
        grade = "derived" if p is None else GRADE.get(getattr(p, "evidence_class", ""), "unclassified")
        species = sorted(set(rxn.reactants) | set(rxn.products))
        measured = {s: [e for e in ledger if _measures(e, s)] for s in species}
        series = sorted({s for s, es in measured.items() if any(e.get("kind") == "time_series" for e in es)})
        hits = [c["id"] for c in conflicts if rxn.parameter_key and (
            rxn.parameter_key in json.dumps(c.get("affects") or []) or re.search(rf"\b{re.escape(rxn.parameter_key)}\b", str(c.get("observation", ""))))]
        open_hits = [c["id"] for c in conflicts if c["id"] in hits and c["status"] in ("open", "error_confirmed", "duplicate_confirmed")]
        fit_state = ident.get(rxn.parameter_key or "", {})
        flux = shares.get(rxn.key)
        share = flux["share"] if flux else 0.0
        odorant_share = flux["odorant_share"] if flux else 0.0
        matters = share > FLUX_THRESHOLD or odorant_share > FLUX_THRESHOLD
        reasons = []
        if matters:
            if grade == "fitted" and (not fit_state or any(v != "identified" for v in fit_state.values())):
                reasons.append("fitted, " + (", ".join(f"{k} {v}" for k, v in fit_state.items()) if fit_state
                                             else "no identification on record"))
            if grade == "bounded":
                reasons.append("bounded from a timescale bracket, not measured")
            if not series:
                reasons.append("no species of the step is measured over time in the corpus")
        if open_hits:
            reasons.append("conflicts: " + ", ".join(open_hits))
        steps.append({
            "reaction": rxn.key,
            "equation": " + ".join(rxn.reactants) + " -> " + " + ".join(rxn.products),
            "parameter": rxn.parameter_key,
            "grade": grade,
            "source": (getattr(p, "source_anchor", "") or "")[:160] if p is not None else "",
            "dossiers": evidence_anchors.dossier_stems(getattr(p, "dossier_anchor", "")) if p is not None else [],
            "identification_b60": fit_state or None,
            "max_flux_share": share,
            "max_odorant_path_share": odorant_share,
            "matters": matters,
            "where": flux,
            "ledger": {s: {"data_sets": len(es), "time_series": sum(e.get("kind") == "time_series" for e in es),
                           "used_by_the_model": sum(str(e.get("status", "")).startswith("in_model") for e in es)}
                       for s, es in measured.items()},
            "species_measured_over_time": series,
            "conflicts": hits,
            "unclear": bool(reasons),
            "why_unclear": reasons,
        })
    steps.sort(key=lambda s: (-s["max_flux_share"], -s["max_odorant_path_share"]))
    constants = _constant_status(steps, parameters, ident, sampling)
    unpinned_steps = {st["reaction"] for st in steps if st["parameter"] in constants and constants[st["parameter"]]["unpinned"]}
    odorant_on_unpinned = {}
    for pot, record in json.loads(FLUX.read_text(encoding="utf-8"))["pots"].items():
        odorant_on_unpinned[pot] = {sp: round(sum(i["share"] for i in record["species"][sp]["made"]["top"]
                                                  if i["reaction"] in unpinned_steps), 3)
                                    for sp in ("MFT", "FFT") if sp in record["species"]}
    carrying = [s for s in steps if s["max_flux_share"] > FLUX_THRESHOLD]
    odorant = [s for s in steps if s["max_odorant_path_share"] > FLUX_THRESHOLD]
    matter = [s for s in steps if s["matters"]]
    return {
        "artifact": "sulfur_reaction_evidence",
        "provenance": provenance.provenance_block(
            "sulfur_reaction_evidence", generated_by="scripts/generators/build_sulfur_reaction_evidence.py",
            inputs=[FLUX, LAPLACE, B23_REPORT, ENVELOPE, data_paths.CONFLICTS_REGISTER, ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/parameters_sulfur.py", ROOT / "src/kinetic_core/parameters.py",
                    *sorted(data_paths.MEASUREMENT_LEDGER_DIR.glob("*.yml"))]),
        "summary": {
            "steps": len(steps),
            "carrying_more_than_1_percent": len(carrying),
            "by_grade_among_carrying": {g: sum(s["grade"] == g for s in carrying) for g in sorted({s["grade"] for s in carrying})},
            "unclear_among_carrying": sum(s["unclear"] for s in carrying),
            "on_the_odorant_path": len(odorant),
            "odorant_path_species": list(ODORANT_PATH),
            "steps_that_matter": len(matter),
            # the headline the waves must move (2026-10-29): the share of MFT's and FFT's formation flux, per reference
            # pot, through steps whose constant is unpinned
            "mft_fft_formation_on_unpinned_steps": odorant_on_unpinned,
            "unpinned_constants_that_matter": sorted(k for k, c in constants.items() if c["unpinned_and_matters"]),
            "unpinned_and_rate_not_sampled": sorted(k for k, c in constants.items()
                                                    if c["unpinned_and_matters"] and not c["envelope_samples"].get("log10k")),
            "unclear_among_steps_that_matter": sum(s["unclear"] for s in matter),
            "threshold": FLUX_THRESHOLD,
        },
        "constants": constants,
        "steps": steps,
    }


def markdown(payload: dict) -> str:
    s = payload["summary"]

    def table(rows, heading):
        out = [heading, "",
               "| step | equation | grade | B60 fit | share of pot flux | of an odorant-path species | measured over time | why unclear |",
               "| --- | --- | --- | --- | --- | --- | --- | --- |"]
        for st in rows:
            w = st["where"] or {}
            fit = ", ".join(f"{k} {v}" for k, v in (st["identification_b60"] or {}).items()) or ("not free" if st["grade"] == "fitted" else "")
            out.append(f"| `{st['reaction']}` | {st['equation']} | {st['grade']} | {fit} | {st['max_flux_share']:.1%} "
                       f"({(w.get('pot') or '').split('_')[0]}) | {st['max_odorant_path_share']:.0%} {w.get('odorant') or ''} "
                       f"| {', '.join(st['species_measured_over_time']) or '**none**'} | {'; '.join(st['why_unclear']) or '—'} |")
        return out + [""]

    steps = payload["steps"]
    bulk = [st for st in steps if st["max_flux_share"] > s["threshold"]]
    odorant = sorted([st for st in steps if st["max_odorant_path_share"] > s["threshold"] and st["max_flux_share"] <= s["threshold"]],
                     key=lambda st: -st["max_odorant_path_share"])
    unclear = [st for st in steps if st["matters"] and st["unclear"]]
    lines = [
        "# Sulfur-lane reaction evidence (strategy M4)", "",
        "*Generated by `scripts/generators/build_sulfur_reaction_evidence.py` from the shipped constants, B60's Laplace, the "
        "flux attribution's eight reference pots, the measurement ledgers and the conflicts register. A step MATTERS when it "
        "carries more than 1 % of a pot's total reaction flux, or more than 1 % of the formation or loss of an odorant-path "
        f"species ({', '.join(s['odorant_path_species'])}). A step that matters is UNCLEAR when its constant is fitted but "
        "not identified (by B60's Laplace, or, for a constant B60 kept frozen, by B2.3, the last fit that freed it), sits on a bound, or is bounded from a "
        "bracket; when none of its species is measured over time; or when an open conflict names it.*", "",
        f"{s['steps']} steps; **{s['steps_that_matter']} matter** ({s['carrying_more_than_1_percent']} by bulk flux, "
        f"{s['on_the_odorant_path']} on the odorant path); **{s['unclear_among_steps_that_matter']} of them are unclear.**", "",
        "## The unclear steps, which the sub-models must test first", "",
    ]
    lines += [f"- `{st['reaction']}` ({st['equation']}): {'; '.join(st['why_unclear'])}" for st in unclear] + [""]
    unpinned = [(k, c) for k, c in payload["constants"].items() if c["unpinned_and_matters"]]
    lines += ["## The headline to move: MFT and FFT formation through unpinned steps", "",
              "Share of each thiol's formation flux made by a step whose own constant is unpinned. This counts the LAST step "
              "only: MFT's last step (norfuraneol + H2S) is identified, but norfuraneol itself comes through unpinned steps "
              "(pentose + cysteine -> 1-deoxypentosone -> norfuraneol) in every pot, so read a low MFT number as 'the final "
              "step is pinned', not 'MFT is pinned'.", "",
              "| pot | MFT (last step) | FFT (last step) |", "| --- | --- | --- |"]
    lines += [f"| {pot} | {v.get('MFT', 0):.0%} | {v.get('FFT', 0):.0%} |" for pot, v in s["mft_fft_formation_on_unpinned_steps"].items()] + [""]
    lines += ["## Unpinned constants that matter", "",
              "Fitted but not identified (by B60, or by B2.3 for a constant frozen since), on a step that matters. "
              "`docs/reference/unpinned_constants.yml` must list each one whose rate the envelope does not sample "
              "(`scripts/ci/unpinned_constants_gate.py`).", "",
              "| constant | identification | envelope samples the rate / the barrier | steps |", "| --- | --- | --- | --- |"]
    lines += [f"| `{k}` | {', '.join(f'{a} {b}' for a, b in (c['identification'] or {}).items()) or 'none on record'} | "
              f"{'yes' if c['envelope_samples'].get('log10k') else '**no**'} / {'yes' if c['envelope_samples'].get('Ea') else 'no'} "
              f"| {', '.join(c['steps'])} |" for k, c in unpinned] + [""]
    lines += table(bulk, "## Steps carrying more than 1 % of a pot's total reaction flux")
    lines += table(odorant, "## Odorant-path steps (below 1 % of the bulk flux, above 1 % of an odorant-path species)")
    minor = [st["reaction"] for st in steps if not st["matters"]]
    lines += [f"## Steps that do not matter in the eight reference pots ({len(minor)})", "", ", ".join(f"`{m}`" for m in minor), ""]
    return "\n".join(lines)


def main() -> int:
    payload = build()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUT.with_suffix(".md").write_text(markdown(payload), encoding="utf-8")
    s = payload["summary"]
    print(f"wrote {data_paths.rel(OUT)}: {s['steps']} steps, {s['steps_that_matter']} matter "
          f"({s['carrying_more_than_1_percent']} bulk, {s['on_the_odorant_path']} odorant path), "
          f"{s['unclear_among_steps_that_matter']} unclear")
    return 0


if __name__ == "__main__":
    sys.exit(main())
