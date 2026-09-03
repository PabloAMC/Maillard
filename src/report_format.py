"""Number and table formatting shared by the markdown reports.

2026-09-03 (improvement backlog, "shared report helpers"). Eight ``_fmt`` copies with three
different dash sentinels, and forty-four files hand-rolling ``|---|---|`` rows. The frozen wave
generators keep their copies (a change to them is a new wave); everything live uses this.
"""
from __future__ import annotations

import math
from typing import Any, Iterable, Optional, Sequence


def fmt_number(value: Optional[float], *, style: str = "auto", dash: str = "-") -> str:
    """Render a number for a report cell.

    ``style="auto"`` (the scorecard rule): three decimals in the readable range, three
    significant figures for |v| >= 1000 or |v| < 0.01. ``style="3g"``: always three
    significant figures. ``None`` and non-finite values render as ``dash``; bools as yes/no.
    """
    if value is None:
        return dash
    if isinstance(value, bool):
        return "yes" if value else "no"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return dash
    if style == "3g":
        return f"{number:.3g}"
    if abs(number) >= 1000 or (abs(number) < 0.01 and number != 0):
        return f"{number:.3g}"
    return f"{number:.3f}"


def fmt_rate(agree: int, evaluable: int) -> str:
    """``"4/5 (80%)"``; ``"0/0"`` when nothing was evaluable."""
    if evaluable <= 0:
        return f"{agree}/{evaluable}"
    return f"{agree}/{evaluable} ({agree / evaluable:.0%})"


def md_table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    """A GitHub-flavoured markdown table; cells are rendered with ``str`` (pipes escaped)."""

    def cell(value: Any) -> str:
        return str(value).replace("|", "\\|")

    lines = ["| " + " | ".join(cell(h) for h in headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(cell(v) for v in row) + " |")
    return "\n".join(lines)


__all__ = ["fmt_number", "fmt_rate", "md_table"]
