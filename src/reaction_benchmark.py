from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_REACTION_BENCHMARK_FILE = _repo_root() / "data" / "lit" / "reaction_benchmark_set.json"


def _normalize_alias(value: str) -> str:
    text = str(value).strip().lower()
    for token in ("-", " ", "/"):
        text = text.replace(token, "_")
    return text


@dataclass(frozen=True)
class ReactionBenchmarkEntry:
    reaction_id: str
    reaction_family: str
    motif_family: str
    benchmark_visible_gap: str
    reference_barrier_kcal: float
    barrier_uncertainty_kcal: float
    reference_method: str
    source: str
    candidate_aliases: List[str]
    allowed_roles: List[str]
    decision_importance: str = "high"

    @property
    def aliases(self) -> List[str]:
        values = [self.reaction_id, self.reaction_family, *self.candidate_aliases]
        normalized: List[str] = []
        seen = set()
        for value in values:
            alias = _normalize_alias(value)
            if alias and alias not in seen:
                normalized.append(alias)
                seen.add(alias)
        return normalized


def load_reaction_benchmark_payload(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_REACTION_BENCHMARK_FILE
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_reaction_benchmark_entries(file_path: Optional[Path | str] = None) -> List[ReactionBenchmarkEntry]:
    payload = load_reaction_benchmark_payload(file_path)
    return [ReactionBenchmarkEntry(**row) for row in payload.get("entries", [])]


def build_reaction_benchmark_alias_index(
    entries: Iterable[ReactionBenchmarkEntry],
) -> Dict[str, ReactionBenchmarkEntry]:
    index: Dict[str, ReactionBenchmarkEntry] = {}
    for entry in entries:
        for alias in entry.aliases:
            index[alias] = entry
    return index


def build_reaction_benchmark_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    entries = load_reaction_benchmark_entries(file_path)
    payload = {
        "summary": {
            "benchmark_id": "maillard_reaction_benchmark_v1",
            "entry_count": len(entries),
            "motif_family_count": len({entry.motif_family for entry in entries}),
            "reaction_families": [entry.reaction_family for entry in entries],
            "source_file": str(Path(file_path) if file_path is not None else DEFAULT_REACTION_BENCHMARK_FILE),
        },
        "entries": [asdict(entry) for entry in entries],
    }
    return payload


def render_reaction_benchmark_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Reaction Benchmark Set",
        "",
        "| Reaction Family | Motif Family | Reference Barrier (kcal/mol) | Uncertainty | Method | Allowed Roles | Source |",
        "| --- | --- | ---: | ---: | --- | --- | --- |",
    ]
    for row in payload.get("entries", []):
        lines.append(
            f"| {row.get('reaction_family', 'unknown')} | {row.get('motif_family', 'unknown')} | "
            f"{float(row.get('reference_barrier_kcal', 0.0)):.1f} | {float(row.get('barrier_uncertainty_kcal', 0.0)):.1f} | "
            f"{row.get('reference_method', 'unknown')} | {', '.join(str(item) for item in row.get('allowed_roles', []))} | "
            f"{row.get('source', 'unknown')} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark id: {summary.get('benchmark_id', 'unknown')}",
            f"Entries: {int(summary.get('entry_count', 0))}",
            f"Motif families: {int(summary.get('motif_family_count', 0))}",
        ]
    )
    return "\n".join(lines) + "\n"