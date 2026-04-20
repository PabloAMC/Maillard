from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from functools import lru_cache

import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401


@lru_cache(maxsize=1)
def _assert_latex_contract() -> None:
    latex = shutil.which("latex")
    if latex is None:
        raise RuntimeError(
            "LaTeX is required for Maillard validation figures, but the 'latex' binary is not available on PATH. "
            "Install BasicTeX or expose /Library/TeX/texbin before generating assets."
        )
    if shutil.which("dvipng") is None:
        raise RuntimeError(
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
            raise RuntimeError(
                "LaTeX is present but incomplete for Matplotlib science plots. "
                "Install the required TeX packages such as type1cm/cm-super, or provide a complete Docker/local toolchain.\n\n"
                f"latex stdout tail:\n{'\n'.join(completed.stdout.splitlines()[-20:])}\n\n"
                f"latex stderr tail:\n{'\n'.join(completed.stderr.splitlines()[-20:])}"
            )


def configure_science_plot_style() -> None:
    texbin = Path("/Library/TeX/texbin")
    if texbin.exists():
        current_path = os.environ.get("PATH", "")
        texbin_str = str(texbin)
        if texbin_str not in current_path.split(":"):
            os.environ["PATH"] = f"{texbin_str}:{current_path}" if current_path else texbin_str

    matplotlib.use("Agg")
    _assert_latex_contract()
    plt.style.use(["science"])
