# Koelsch, Downes & Labuza 1991 — EXTRACTION (hexanal formation from soybean oil at 23 °C as a function of oxygen concentration)

**Source on disk:** `data/articles/koelsch1991.pdf` (0.6 MB scan with a text layer; downloaded
2026-09-11 at the 2026-09-11 reading-list row's request, which named it "the nearest lipid-rate
candidate"). Read 2026-09-11 via `pdftotext -layout`. Wave B36. The earlier verdict in
`research_round3_channels.md` (D.2, "DO NOT ORDER") stands and is now backed by the full text.

| field | value |
|---|---|
| Title | "Hexanal Formation via Lipid Oxidation as a Function of Oxygen Concentration: Measurement and Kinetics" |
| Venue | Journal of Food Science 56(3):816–820 (1991) |
| DOI | 10.1111/j.1365-2621.1991.tb05389.x |
| Systems | a freeze-dried model system of soybean oil on a solid support, held at 23 ± 1 % relative humidity and **23 ± 2 °C in the dark**, under four constant oxygen concentrations in a flow-through cell; hexanal followed by headspace GC |
| What is measured | hexanal (ppm w/w in oil) against time at each oxygen level; two derived kinetic models (a cubic model for the initial stage, an extended model for the accelerated stage); rate constants against oxygen concentration (Table 2 of the paper) |

## 1. Why it is not the lipid temperature source

The lipid lane's open ask is the **temperature dependence** of hydroperoxide decomposition to
hexanal in a food matrix between 60 and 140 °C. This paper has **one temperature, 23 °C**. Its
rate constants describe the dependence on oxygen concentration at room temperature, in a dry model
system, and its own text derives temperature effects only by citing the Arrhenius relationship, not
by measurement. Nothing here can be transcribed into a barrier or a Q10, and nothing is.

## 2. What it does establish, for the record

The rate of hexanal formation rises with oxygen concentration with the hyperbolic form the paper
expected from a monomolecular initiation (inverse rate against reciprocal oxygen is linear), and the
break between the initial and accelerated stages moves with oxygen logarithmically. The lipid lane's
oxygen-access declaration (`vessel.atmosphere`) is in the spirit of this result; no number is taken.
