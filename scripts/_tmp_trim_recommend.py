"""Temporary script: remove redundant function bodies from recommend.py (P2.2 refactor)."""
import re

with open("src/recommend.py", "r") as f:
    content = f.read()

start_sentinel = "\n\ndef _henry_entry_for_species(species:"
end_sentinel = "\n    return projected_ppb\n\n\nclass Recommender:"

start_idx = content.find(start_sentinel)
end_idx = content.find(end_sentinel)

if start_idx == -1:
    raise RuntimeError("start sentinel not found")
if end_idx == -1:
    raise RuntimeError("end sentinel not found")

new_content = (
    content[:start_idx]
    + "\n\n# ── Functions above moved to src/projection.py (P2.2). "
    + "Re-exported via import at top. ──\n\n\nclass Recommender:"
    + content[end_idx + len(end_sentinel):]
)

with open("src/recommend.py", "w") as f:
    f.write(new_content)

orig_lines = len(content.splitlines())
new_lines = len(new_content.splitlines())
print(f"Done. {orig_lines} → {new_lines} lines (removed {orig_lines - new_lines} lines).")
