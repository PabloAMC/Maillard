# Engel & Schieberle 2004 — EXTRACTION (cysteamine + 2-13C-fructose, 0.1 M phosphate pH 7.0, 7 mL, laboratory autoclave "raising the temperature within 20 min from 20°C to 145°C"; labelling of 2-acetyl-2-thiazoline / 2-acetylthiazole carbon; a dry-silica arm)

**Source on disk:** `data/articles/engel2004.pdf` (owner's download, 2026-10-28; open access). Read via `pdftotext -q`.

| field | value |
|---|---|
| Venue | Czech J. Food Sci. 22 (Special Issue), S120-S122 (2004) |
| DOI | 10.17221/10632-CJFS (from the journal's landing page; check against the print) |
| Laboratory | Deutsche Forschungsanstalt für Lebensmittelchemie / TU München (Schieberle), the Hofmann 1995-1998 laboratory |

## Why it matters

Wave B64 (`results/validation/kinetic_core_b64_prereg.md`) cannot decide whether the anchor laboratory's "20 min at
145 °C" is a ramp or a ramp plus hold. This paper, from the same laboratory and apparatus class, prints the ramp
wording again, for a 7 mL charge [M]:

> "reacted with 2-13C-fructose (0.7 mmol) in a laboratory autoclave by raising the temperature within 20 min from
> 20°C to 145°C. To mimic dry heating conditions, ..." (Experimental, "Model reactions")

No hold and no cooling are stated; the next sentence turns to the dry arm. A 7 mL charge still taking "within 20 min"
suggests a fixed, controlled ramp rather than full-power heating of the liquid [! inference].

## What the repo takes

- **Heating fact** for `docs/validation/thermal_history_facts.yml` (Munich laboratory, time_convention `ramp_only`
  stated): supports B64's ramp-only reading; it does not settle it, because the 1998 paper's own pots are not this
  paper's pots (register CR-03 stays open).
- Nothing for the kinetic constants: the paper's quantities are isotopologue distributions of cysteamine-derived
  heterocycles, not rates.
