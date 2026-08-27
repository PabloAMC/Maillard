from __future__ import annotations

import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from functools import lru_cache
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401


logger = logging.getLogger(__name__)

# 2026-08-27 (Wave I): environment switch that restores the pre-Wave-I strict
# behaviour (raise instead of degrade). CI figure jobs set this so that a
# LaTeX-less runner is still a hard failure there, while the documented local
# conda environment -- which has no `dvipng` and cannot get one, because the
# TeX Live install on the reference machine is root-owned -- can still produce
# a report.
STRICT_LATEX_ENV_VAR = "MAILLARD_STRICT_LATEX"

_TRUTHY = {"1", "true", "yes", "on"}


def strict_latex_requested() -> bool:
    """True when the caller demands a real LaTeX toolchain (no mathtext fallback)."""
    return os.environ.get(STRICT_LATEX_ENV_VAR, "").strip().lower() in _TRUTHY


@lru_cache(maxsize=1)
def latex_contract_failure() -> Optional[str]:
    """Return None when the LaTeX toolchain is usable, else a human-readable reason.

    2026-08-27 (Wave I): this used to be `_assert_latex_contract()`, which raised.
    The probing logic is unchanged; only the reporting changed from "raise" to
    "return a reason", so that `configure_science_plot_style()` can decide
    between raising (strict mode) and degrading to mathtext.
    """
    latex = shutil.which("latex")
    if latex is None:
        return (
            "LaTeX is required for Maillard validation figures, but the 'latex' binary is not available on PATH. "
            "Install BasicTeX or expose /Library/TeX/texbin before generating assets."
        )
    if shutil.which("dvipng") is None:
        return (
            "LaTeX is present but the Matplotlib TeX backend also requires 'dvipng'. "
            "Install dvipng in the local toolchain or in the Docker image before generating assets."
        )

    tex_source = """\\documentclass{article}
\\usepackage{type1cm}
\\usepackage{type1ec}
\\begin{document}
lp
\\end{document}
"""
    with tempfile.TemporaryDirectory(prefix="maillard-latex-check-") as tmpdir:
        tex_path = Path(tmpdir) / "check.tex"
        tex_path.write_text(tex_source, encoding="utf-8")
        completed = subprocess.run(
            [latex, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            return (
                "LaTeX is present but incomplete for Matplotlib science plots. "
                "Install the required TeX packages such as type1cm/cm-super, or provide a complete Docker/local toolchain.\n\n"
                f"latex stdout tail:\n{'\n'.join(completed.stdout.splitlines()[-20:])}\n\n"
                f"latex stderr tail:\n{'\n'.join(completed.stderr.splitlines()[-20:])}"
            )
    return None


def _assert_latex_contract() -> None:
    """Backwards-compatible strict probe: raises when the toolchain is unusable."""
    reason = latex_contract_failure()
    if reason is not None:
        raise RuntimeError(reason)


_FALLBACK_WARNING_EMITTED = False


def _warn_fallback_once(reason: str) -> None:
    global _FALLBACK_WARNING_EMITTED
    if _FALLBACK_WARNING_EMITTED:
        return
    _FALLBACK_WARNING_EMITTED = True
    logger.warning(
        "Falling back to Matplotlib's mathtext renderer (text.usetex=False): the LaTeX "
        "toolchain contract is not satisfied on this machine. Figures will render, but "
        "glyphs will not be typeset by LaTeX. Set %s=1 to make this a hard error instead.\n%s",
        STRICT_LATEX_ENV_VAR,
        reason,
    )


def configure_science_plot_style() -> None:
    """Install the `scienceplots` style, degrading to mathtext when LaTeX is missing.

    2026-08-27 (Wave I): previously this raised `RuntimeError` whenever `latex`,
    `dvipng` or the TeX packages were missing, which made `--report` unusable on
    the documented conda environment. It now warns once and switches to
    matplotlib's built-in mathtext renderer. Set `MAILLARD_STRICT_LATEX=1` to
    restore the raising contract (used by the CI figure jobs).
    """
    texbin = Path("/Library/TeX/texbin")
    if texbin.exists():
        current_path = os.environ.get("PATH", "")
        texbin_str = str(texbin)
        if texbin_str not in current_path.split(":"):
            os.environ["PATH"] = f"{texbin_str}:{current_path}" if current_path else texbin_str

    matplotlib.use("Agg")

    reason = latex_contract_failure()
    if reason is None:
        plt.style.use(["science"])
        return

    if strict_latex_requested():
        raise RuntimeError(reason)

    _warn_fallback_once(reason)
    plt.style.use(["science", "no-latex"])
    # `no-latex` already clears usetex; pin it explicitly so a stale rcParam from
    # an earlier successful configure cannot leak a LaTeX request into this run.
    matplotlib.rcParams["text.usetex"] = False
