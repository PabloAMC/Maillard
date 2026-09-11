# Zhang, Guo, Liu & Chang 2012 — EXTRACTION (soymilk: raw against stove cooking and two UHT schedules, in ppm, with three grinding temperatures)

**Source on disk:** `data/articles/zhang2012.pdf` (0.34 MB; downloaded 2026-09-11). Read 2026-09-11
via `pdftotext -layout`, Table 1 verified against 400-dpi page crops. Wave B37.

| field | value |
|---|---|
| Title | "Off-Flavor Related Volatiles in Soymilk As Affected by Soybean Variety, Grinding, and Heat-Processing Methods" |
| Venue | Journal of Agricultural and Food Chemistry 2012, 60, 7457–7462 |
| DOI | 10.1021/jf3016199 |
| System | soymilk from two varieties (Prosoy; a black soybean), bean:water 1:10 (w/w) |
| Grinding | **cold (~2 °C), ambient (20 °C), hot (80.5 °C water)**, 3 min at 10 000 rpm, then filtered through muslin |
| Heat processes | **stove: boiled then held 20 min**; **one-phase UHT 140 °C / 5 s** (F0 6.62); **two-phase UHT 120 °C / 80 s + 140 °C / 4 s** (F0 6.35); the UHT line carries a 50 kPa vacuum chamber, applied once and twice respectively |
| Quantification | SPME–GC-FID; **compound-specific standard curves** in a 2 % cow's-milk matrix with **2-methyl-3-heptanone internal standard**; detection limits printed per compound (hexanal ~1 ppb, 1-hexanol 2.5, 1-octen-3-ol 5, (E)-2-nonenal 10, (E,E)-2,4-decadienal 25 ppb); n = 6 |

## 1. Why this paper matters

**The strongest quantification of the sixteen**, and it prints a raw column beside three defined heat
processes on the same material. The repository already scores a pea-protein UHT bundle
(`pea_isolate_uht_140C_Trikusuma2019`) and declares its unheated column under Amendment 37; this is
the soy analogue, at the same UHT temperature, from an independent laboratory.

## 2. Table 1 — hexanal (ppm, SD in parentheses, n = 6)

| variety | grinding | raw | stove (boil + 20 min) | one-phase UHT | two-phase UHT |
|---|---|---:|---:|---:|---:|
| Prosoy | cold | 6.60 (0.16) | 0.27 (0.04) | 0.25 (0.04) | 0.14 (0.01) |
| Prosoy | ambient | 3.23 (0.79) | 0.54 (0.20) | 0.34 (0.01) | 0.00 (0.00) |
| Prosoy | hot | 0.051 (0.004) | 0.006 (0.002) | 0.005 (0.003) | 0.00 (0.00) |
| black | cold | 7.12 (0.47) | 0.63 (0.06) | 0.26 (0.04) | 0.17 (0.01) |
| black | ambient | 7.16 (0.89) | 1.19 (0.34) | 0.52 (0.03) | 0.048 (0.0005) |
| black | hot | 0.16 (0.03) | 0.027 (0.004) | 0.012 (0.005) | 0.006 (0.0003) |

## 3. Table 1 — the other seven compounds, cold-ground Prosoy and cold-ground black (the least-perturbed arms)

| compound | Prosoy cold raw | Prosoy cold stove | Prosoy cold 1-UHT | Prosoy cold 2-UHT | black cold raw | black cold stove | black cold 1-UHT | black cold 2-UHT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1-Hexanol | 0.34 (0.03) | 0.00 | 0.00 | 0.00 | 0.16 (0.02) | 0.00 | 0.00 | 0.00 |
| **2-Pentylfuran** | 0.060 (0.001) | **0.24 (0.01)** | 0.062 (0.001) | 0.064 (0.001) | 0.061 (0.001) | **0.16 (0.07)** | 0.058 (0.001) | 0.058 (0.004) |
| 1-Octen-3-one | 0.39 (0.02) | 0.13 (0.01) | 0.13 (0.01) | 0.11 (0.004) | 0.40 (0.02) | 0.18 (0.02) | 0.17 (0.01) | 0.14 (0.02) |
| 1-Octen-3-ol | 0.47 (0.001) | 0.032 (0.006) | 0.037 (0.003) | 0.019 (0.002) | 0.22 (0.01) | 0.047 (0.004) | 0.034 (0.001) | 0.017 (0.002) |
| (E)-2-Nonenal | 0.032 (0.002) | 0.010 (0.002) | 0.00 | 0.00 | 0.037 (0.004) | 0.037 (0.005) | 0.002 (0.001) | 0.00 |
| (E,E)-2,4-Nonadienal | 0.11 (0.002) | 0.090 (0.004) | 0.078 (0.002) | 0.070 (0.002) | 0.11 (0.003) | 0.10 (0.002) | 0.081 (0.002) | 0.074 (0.001) |
| **(E,E)-2,4-Decadienal** | 0.061 (0.017) | **0.41 (0.11)** | 0.30 (0.03) | 0.25 (0.01) | 0.15 (0.04) | **0.78 (0.02)** | 0.47 (0.04) | 0.56 (0.07) |

Ambient- and hot-ground arms are in the source and in the reading notes; the cold arms are printed
here because they are the ones with a genuinely unheated starting state.

**Only these eight compounds are measured.** Nonanal, pentanal, heptanal, benzaldehyde, octanal,
decanal, every pyrazine and furfural are **not** in this paper.

## 4. The two compounds that RISE, and why they matter

Six of the eight fall with heat. **2-pentylfuran and (E,E)-2,4-decadienal rise on stove cooking** —
2-pentylfuran 0.060 → 0.24 ppm, 2,4-decadienal 0.061 → 0.41 ppm in cold-ground Prosoy — and then fall
again under UHT. The authors say so verbatim of 2-pentylfuran: "in contrast to other odor compounds,
it increased". That is a thermal formation this model's lipid lane claims to make, measured against a
declared starting level, with a boil-plus-20-minutes thermal program: the closest thing in the sixteen
papers to a scoreable lipid formation row.

## 5. Limits

1. **"Raw" is raw soymilk, not raw bean** — it is post-grinding and post-filtration, and grinding is
   itself thermal: hot grinding at 80.5 °C destroys about 99 % of lipoxygenase before any cook. The
   unperturbed baseline is **cold-ground raw**, not "raw" generally.
2. **The UHT arms carry a vacuum step** (50 kPa, once for one-phase and twice for two-phase) that the
   authors say removes volatiles. A fall across a UHT column is heat plus vacuum stripping, not heat
   alone. The stove arm has no vacuum.
3. **Stove "boiling" has no stated temperature**, only "20 min after boiling".
4. **F0 = 6.62 / 6.35 at Z = 10 are microbial lethality values**, not chemical kinetics, and must not
   be repurposed as a barrier.
5. Three cells in the printed table carry a typesetting defect where the value and the significance
   code run together ("0.00 5B1", "0.06 1A2", "0.06 3A2"); they read as 0.005, 0.061 and 0.063 and
   are recorded here with that note rather than silently normalised.

## 6. Verdict

The best-quantified unheated-versus-heated pair among the sixteen, and a candidate both for an
Amendment 37 carried-level declaration and for a scoreable lipid formation row on the stove arm.
Neither is installed in B37; both are named as pre-registered questions for the wave that follows.
