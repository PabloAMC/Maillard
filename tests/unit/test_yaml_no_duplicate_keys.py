"""No tracked YAML file repeats a mapping key.

YAML keeps the last of two equal keys without a warning. On 2026-11-24 the claims panel's ``sources`` map held
``mottram2002`` twice; the later entry, with a wrong DOI, silently replaced the one checked against the PDF.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


_BASE = getattr(yaml, "CSafeLoader", yaml.SafeLoader)      # the C parser: the lever-table contrast files run to 2 MB


class _StrictLoader(_BASE):
    """A safe loader that records every repeated key in a mapping."""

    duplicates: list


def _construct_mapping(loader, node, deep=False):
    seen = {}
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            hash(key)
        except TypeError:
            continue
        line = key_node.start_mark.line + 1
        if key in seen:
            loader.duplicates.append(f"{key!r} at lines {seen[key]} and {line}")
        seen[key] = line
    return _BASE.construct_mapping(loader, node, deep)


_StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def _tracked_yaml():
    try:
        out = subprocess.run(["git", "ls-files", "*.yml", "*.yaml"], cwd=ROOT, capture_output=True, text=True, check=True)
        return [ROOT / f for f in out.stdout.split()]
    except (OSError, subprocess.CalledProcessError):
        return sorted(p for p in ROOT.rglob("*.y*ml") if ".git" not in p.parts)


def _duplicates(path: Path) -> list:
    loader = _StrictLoader(path.read_text(encoding="utf-8"))
    loader.duplicates = []
    try:
        while loader.check_data():
            loader.get_data()
    except yaml.YAMLError:
        return []          # not this test's concern: a file that does not parse fails where it is read
    finally:
        loader.dispose()
    return loader.duplicates


CONTRASTS = "docs/reference/lever_table/contrasts/"      # the keyed lever-table contrasts: ~200 files, up to 2 MB each


def _assert_none_repeat(files):
    assert files, "no YAML files found"
    found = {str(p.relative_to(ROOT)): d for p in files if p.exists() and (d := _duplicates(p))}
    assert not found, "repeated YAML keys (the later one silently wins):\n" + "\n".join(
        f"  {f}: {', '.join(d)}" for f, d in found.items())


def test_no_curated_yaml_repeats_a_key():
    """Every tracked YAML file but the lever-table contrasts (the claims panel, data/, the registers)."""
    _assert_none_repeat([p for p in _tracked_yaml() if not str(p.relative_to(ROOT)).startswith(CONTRASTS)])


@pytest.mark.slow
def test_no_lever_table_contrast_repeats_a_key():
    _assert_none_repeat([p for p in _tracked_yaml() if str(p.relative_to(ROOT)).startswith(CONTRASTS)])


def test_the_loader_catches_a_repeated_key(tmp_path):
    p = tmp_path / "dup.yml"
    p.write_text("sources:\n  a: 1\n  b: 2\n  a: 3\n", encoding="utf-8")
    assert _duplicates(p) == ["'a' at lines 2 and 4"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
