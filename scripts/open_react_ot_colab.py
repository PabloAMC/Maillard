"""Print (and optionally open) the React-OT Colab notebook URL.

The notebook lives at notebooks/react_ot_colab_gpu.ipynb on the
``qm-barriers`` branch. This helper centralises the canonical URL so other
scripts and docs do not drift.

Usage:
    python scripts/open_react_ot_colab.py            # print URL
    python scripts/open_react_ot_colab.py --open     # try to launch a browser
    python scripts/open_react_ot_colab.py --github   # print GitHub URL instead
"""
from __future__ import annotations

import argparse
import sys
import webbrowser

GITHUB_REPO = "PabloAMC/Maillard"
DEFAULT_BRANCH = "qm-barriers"
NOTEBOOK_PATH = "notebooks/react_ot_colab_gpu.ipynb"


def colab_url(branch: str = DEFAULT_BRANCH) -> str:
    return f"https://colab.research.google.com/github/{GITHUB_REPO}/blob/{branch}/{NOTEBOOK_PATH}"


def github_url(branch: str = DEFAULT_BRANCH) -> str:
    return f"https://github.com/{GITHUB_REPO}/blob/{branch}/{NOTEBOOK_PATH}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch", default=DEFAULT_BRANCH)
    parser.add_argument("--github", action="store_true", help="Print GitHub URL instead of Colab.")
    parser.add_argument("--open", dest="do_open", action="store_true", help="Attempt to open in a browser.")
    args = parser.parse_args()

    url = github_url(args.branch) if args.github else colab_url(args.branch)
    print(url)
    if args.do_open:
        opened = webbrowser.open(url)
        if not opened:
            print("warning: webbrowser.open returned False (no display?)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
