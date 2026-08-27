import json
import re
from pathlib import Path

ROOT = Path(".")
GDR_DIR = ROOT / "data" / "Gemini_Deep_Research" / "raw"
OUTPUT_PATH = ROOT / "scripts" / "get_details_fgh_output.txt"

files = ["Gemini_deep_research_F.md", "Gemini_deep_research_G.md", "Gemini_deep_research_H.md"]

output_lines = []

for fn in files:
    fp = GDR_DIR / fn
    if not fp.exists():
        continue
    output_lines.append(f"\n--- FILE: {fn} ---")
    with open(fp, "r", encoding="utf-8") as f:
        content = f.read()
    
    # Split by lines
    lines = content.splitlines()
    for i, line in enumerate(lines):
        # Match digit followed by dot and optional space, then bibliography content
        # E.g. "1. Solís-Calero" or "1\. Smuda"
        match = re.match(r'^(\d+)[\\.]+\s+(.*)', line.strip())
        if match:
            # Let's check if it has a score in the format [1-8]/8
            # E.g. "5/8." or "6/8"
            score_match = re.search(r'\d/8', line)
            if score_match or any(term in line.lower() for term in ["doi", "matrix", "analyte"]):
                output_lines.append(f"L{i+1}: {line}")
                # Print next lines if they are part of the description
                for offset in range(1, 8):
                    if i + offset < len(lines):
                        next_line = lines[i+offset].strip()
                        if next_line and not re.match(r'^(\d+)[\\.]+\s+', next_line):
                            output_lines.append(f"   + {lines[i+offset]}")
                        else:
                            break

with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
    f.write("\n".join(output_lines))
print(f"Details written to {OUTPUT_PATH}")
