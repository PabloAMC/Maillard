#!/usr/bin/env python3
"""
Environment smoke-test for the Maillard Reaction Computational Framework.

Verifies that all required packages are importable and key binaries work.
Run with: python scripts/check_env.py

Exit 0 = all checks pass.
Exit 1 = one or more checks failed.
"""

import sys
import subprocess
import importlib
from pathlib import Path

# ── Colour helpers ──────────────────────────────────────────────────────────
GREEN = "\033[92m"
RED   = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"
BOLD  = "\033[1m"

def ok(msg: str) -> None:
    print(f"  {GREEN}✓{RESET} {msg}")

def fail(msg: str) -> None:
    print(f"  {RED}✗{RESET} {msg}")

def warn(msg: str) -> None:
    print(f"  {YELLOW}⚠{RESET} {msg}")

def section(title: str) -> None:
    print(f"\n{BOLD}{title}{RESET}")
    print("─" * 60)


# ── Check functions ─────────────────────────────────────────────────────────

def check_python_version() -> bool:
    major, minor = sys.version_info[:2]
    if major == 3 and minor >= 10:
        ok(f"Python {major}.{minor} (≥ 3.10 required)")
        return True
    else:
        fail(f"Python {major}.{minor} — need Python ≥ 3.10")
        return False


def check_import(package: str, import_as: str | None = None, version_attr: str = "__version__") -> bool:
    target = import_as or package
    try:
        mod = importlib.import_module(target)
        version = getattr(mod, version_attr, "unknown")
        ok(f"{package} ({version})")
        return True
    except ImportError as e:
        fail(f"{package} — {e}")
        return False


def check_binary(binary: str, args: list[str] | None = None, label: str | None = None) -> bool:
    label = label or binary
    cmd = [binary] + (args or ["--version"])
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        first_line = (result.stdout or result.stderr).strip().split("\n")[0]
        ok(f"{label}: {first_line[:80]}")
        return True
    except FileNotFoundError:
        fail(f"{label} — binary not found in PATH")
        return False
    except subprocess.TimeoutExpired:
        fail(f"{label} — timed out")
        return False


def check_rdkit_smiles_roundtrip() -> bool:
    """Validate that known SMILES round-trip and have correct heavy-atom counts."""
    from rdkit import Chem
    try:
        smiles_pairs = [
            ("2-Methyl-3-furanthiol (MFT)", "Cc1occc1S", None),
            ("2-Furfurylthiol (FFT)", "SCc1ccco1", None),
            ("Methional", "CSCCC=O", None),
            ("L-Cysteine", "N[C@@H](CS)C(=O)O", None),
            ("D-Ribose (Open)", "OC[C@@H](O)[C@@H](O)[C@@H](O)C=O", 5), # 5 carbons
            ("D-Glucose (Open)", "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O", 6), # 6 carbons
        ]
        all_ok = True
        for name, smi, expected_c in smiles_pairs:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                fail(f"  RDKit SMILES parse failed: {name} ({smi})")
                all_ok = False
            else:
                out_smi = Chem.MolToSmiles(mol)
                msg = f"  RDKit round-trip '{name}': {smi} → {out_smi}"
                
                # Check explicit carbon count if requested
                if expected_c is not None:
                    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
                    if c_count != expected_c:
                        fail(f"  RDKit atom count failed for '{name}': Expected {expected_c} carbons, got {c_count}")
                        all_ok = False
                        continue
                    else:
                        msg += f" ({c_count}C confirmed)"
                        
                ok(msg)
        return all_ok
    except Exception as e:
        fail(f"RDKit SMILES round-trip test: {e}")
        return False


def check_yaml_files() -> bool:
    """Verify all required data YAML files exist and are parseable."""
    import yaml
    expected = [
        "data/species/desirable_targets.yml",
        "data/species/off_flavour_targets.yml",
        "data/species/toxic_markers.yml",
        "data/species/precursors.yml",
        "data/reactions/reaction_families.yml",
    ]
    root = Path(__file__).parent.parent
    all_ok = True
    for rel_path in expected:
        p = root / rel_path
        if not p.exists():
            fail(f"Missing: {rel_path}")
            all_ok = False
        else:
            try:
                with open(p) as f:
                    yaml.safe_load(f)
                ok(f"YAML parseable: {rel_path}")
            except yaml.YAMLError as e:
                fail(f"YAML parse error in {rel_path}: {e}")
                all_ok = False
    return all_ok


# ── Main ────────────────────────────────────────────────────────────────────

def main() -> int:
    print(f"\n{BOLD}{'='*60}")
    print(" Maillard Framework — Environment Check")
    print(f"{'='*60}{RESET}")

    failures: list[str] = []

    # ── Python ──────────────────────────────────────────────────────────────
    section("1. Python Runtime")
    if not check_python_version():
        failures.append("python_version")

    # ── Core chemistry packages ──────────────────────────────────────────────
    section("2. Core Chemistry Packages")
    for pkg, import_name in [
        ("rdkit",    "rdkit"),
        ("numpy",    "numpy"),
        ("scipy",    "scipy"),
        ("pandas",   "pandas"),
        ("networkx", "networkx"),
        ("yaml",     "yaml"),
        ("ase",      "ase"),
        ("tqdm",     "tqdm"),
    ]:
        if not check_import(pkg, import_name):
            failures.append(pkg)

    # ── xTB binary ──────────────────────────────────────────────────────────
    section("3. xTB Binary (Tier 1 screening)")
    if not check_binary("xtb", ["--version"], "xtb"):
        warn("xtb not found — install via: conda install -c conda-forge xtb")
        failures.append("xtb_binary")

    # ── PySCF (Tier 2 DFT) ──────────────────────────────────────────────────
    section("4. PySCF (Tier 2 DFT backend)")
    if not check_import("pyscf"):
        warn("Install: pip install pyscf>=2.4.0")
        failures.append("pyscf")

    # ── Skala XC functional ──────────────────────────────────────────────────
    section("5. Microsoft Skala XC Functional (Tier 2)")
    if not check_import("skala", version_attr="__version__"):
        warn("Skala not installed. Install: pip install git+https://github.com/microsoft/skala.git")
        failures.append("skala")

    # ── RMG-Py (Tier 0) ─────────────────────────────────────────────────────
    section("6. RMG-Py (Tier 0 pathway generation)")
    # RMG-Py requires a separate conda environment; check for rmgpy module
    has_rmg = check_import("rmgpy", version_attr="__version__")
    if not has_rmg:
        warn("RMG-Py not installed in current env. Install separately:")
        warn("  conda install -c rmg rmgpy rmg-database")
        warn("  See: https://github.com/ReactionMechanismGenerator/RMG-Py")
        failures.append("rmgpy")

    # ── RDKit SMILES validation ──────────────────────────────────────────────
    section("7. RDKit SMILES Round-trip (target compounds)")
    if "rdkit" not in failures:
        if not check_rdkit_smiles_roundtrip():
            failures.append("rdkit_smiles")

    # ── Data files ───────────────────────────────────────────────────────────
    section("8. Project Data Files")
    if not check_yaml_files():
        failures.append("data_files")

    # ── Summary ─────────────────────────────────────────────────────────────
    print(f"\n{BOLD}{'='*60}{RESET}")
    if not failures:
        print(f"{GREEN}{BOLD}  ✓ All checks passed.{RESET}")
        return 0
    else:
        critical = [f for f in failures if f not in ("rmgpy", "skala")]
        if critical:
            print(f"{RED}{BOLD}  ✗ {len(failures)} check(s) failed:{RESET}")
            for f in failures:
                print(f"      - {f}")
            return 1
        else:
            print(f"{YELLOW}{BOLD}  ⚠ Core checks passed, optional tools missing ({', '.join(failures)}).{RESET}")
            return 0


if __name__ == "__main__":
    sys.exit(main())
