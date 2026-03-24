#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
assets = ROOT / "docs" / "assets"
imgs = [assets / "family_parity.png", assets / "family_benchmark_accuracy.png", assets / "family_coverage.png"]
out = assets / "family_validation_thumbnail.png"

images = [Image.open(p) for p in imgs]
# Normalize heights
max_h = max(im.height for im in images)
resized = []
for im in images:
    if im.height != max_h:
        new_w = int(im.width * (max_h / im.height))
        im = im.resize((new_w, max_h), Image.LANCZOS)
    resized.append(im)

# Add small padding between images
pad = 8
total_w = sum(im.width for im in resized) + pad * (len(resized) - 1)
canvas = Image.new("RGBA", (total_w, max_h), (255, 255, 255, 255))

x = 0
for im in resized:
    canvas.paste(im, (x, 0))
    x += im.width + pad

canvas.save(out, dpi=(150,150))
print(f"Wrote {out}")
