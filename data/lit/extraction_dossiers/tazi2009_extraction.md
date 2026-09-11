# Tazi, Plantevin, Di Falco, Puigserver & Ajandouz 2009 — EXTRACTION (lipoxidation kinetics in almond paste: activation energies and a measured Q10 table over 60–120 °C and three water activities)

**Source on disk:** `data/articles/tazi2009.pdf` (0.35 MB; downloaded 2026-09-11 at this
repository's request). Read 2026-09-11 via `pdftotext -layout`. Wave B37.

| field | value |
|---|---|
| Title | "Effects of light, temperature and water activity on the kinetics of lipoxidation in almond-based products" |
| Venue | Food Chemistry 115 (2009) 958–964 |
| DOI | 10.1016/j.foodchem.2009.01.017 — **read from the printed footer**, not inferred; an earlier draft of this dossier wrote `.018`, reconstructed from the article's PII, and it was wrong |
| Systems | almond paste (an equal mixture of peeled almonds and sugar) and the finished Calisson product — **a real lipid + protein + carbohydrate matrix** |
| Treatments | initial water activity set by storage over saturated salts to **0.38 (LiCl), 0.57 (no salt, unchanged) and 0.72 (KCl)**, then heating at **60 °C (180, 360, 540 min), 80 °C, 100 °C and 120 °C** |
| What is measured | chemiluminescence (CL) and TBARS, as two markers of lipoxidation extent; rate constants, activation energies and Q10 factors derived from them |

## 1. Why this paper matters to this model

The lipid lane bridges its 25 °C anchor to cooking temperatures with `Q10_ASSUMPTION`, a **constant
2–3 with no water-activity term**, and the module's own note says the temperature dependence is
"MEASURED NOWHERE". This paper measures it, in a food matrix, across the model's own temperature
window, and **resolves it by water activity** — which the model carries as an axis on every pot.

## 2. Table 1 — measured Q10, verbatim ("CL and TBARS-based values of Q10 in almond paste in the 60 °C–120 °C temperature interval, depending on the initial water activity value")

| T (°C) | aw 0.38 CL | aw 0.38 TBARS | aw 0.57 CL | aw 0.57 TBARS | aw 0.72 CL | aw 0.72 TBARS |
|---|---:|---:|---:|---:|---:|---:|
| 60–70 | 3.3 | 2.8 | 2.2 | 1.9 | 2.0 | 1.9 |
| 80–90 | 2.9 | 2.5 | 2.1 | 1.8 | 1.8 | 1.8 |
| 100–110 | 2.6 | 2.3 | 1.9 | 1.7 | 1.7 | 1.7 |
| 120–130 | 2.4 | 2.1 | 1.8 | 1.6 | 1.6 | 1.6 |

## 3. Activation energies, verbatim from the Results (r² > 0.956)

> "The activation energy values (r² > 0.956) reached a maximum at aw 0.38 in the case of both CL
> (114 kJ mol⁻¹) and TBARS (100 kJ mol⁻¹). The EaCL values then decreased to 65 kJ mol⁻¹ at aw 0.72;
> whereas the EaTBARS values decreased to 61 kJ mol⁻¹ at aw 0.57, but remained completely unchanged
> (62 kJ mol⁻¹) at aw 0.72."

| marker | aw 0.38 | aw 0.57 | aw 0.72 |
|---|---:|---:|---:|
| CL | **114 kJ/mol** | (between) | **65 kJ/mol** |
| TBARS | **100 kJ/mol** | **61 kJ/mol** | **62 kJ/mol** |

The abstract states the same range as "from 110 kJ mol⁻¹ to 60 kJ mol⁻¹" with Q10 "from 3.3 to 1.6".

## 4. Table 2 — acceleration factors (context)

"Predicted values of the heat factors accelerating lipoxidation in almond paste at various initial
water activity conditions."

| T range (°C) | aw 0.38 CL | aw 0.38 TBARS | aw 0.57 CL | aw 0.57 TBARS | aw 0.72 CL | aw 0.72 TBARS |
|---|---:|---:|---:|---:|---:|---:|
| 60–120 | 531 | 211 | 68 | 30 | 36 | 30 |
| 20–120 | 93028 | 17584 | 2232 | 500 | 724 | 522 |

The 20–120 °C row is **the authors' own extrapolation**, not a measurement — they say so, and give
the regression equations of Q10 against temperature used to make it. The repository does not take
extrapolated rows; the 60–120 row is measured.

## 5. What the repository may take, and the three limits

**Take:** that the temperature dependence of lipoxidation in a real nut matrix at 60–130 °C sits at
Q10 ≈ 1.6–2.4 once the matrix is moist (aw 0.57 and 0.72), rising to 2.4–3.3 only when it is dry
(aw 0.38); and the corresponding Ea of 61–114 kJ/mol, falling as the matrix gets wetter.

**Limits, all three of which must travel with the number:**

1. **CL and TBARS are not hexanal.** They are bulk markers of lipoxidation extent — TBARS in
   particular responds to non-enzymatic browning carbonyls too, which the authors say themselves.
   The model's step is hydroperoxide → hexanal. This paper constrains the lane's overall temperature
   sensitivity; it does not measure the aldehyde-forming branch.
2. **The water-activity range stops at 0.72.** The model's aqueous pots run at aw 0.95–0.99. The
   measured trend is monotone and would extrapolate to a Q10 at or below 1.6 there, but that is an
   extrapolation beyond the data and is recorded as such, not used.
3. **Q10 is not constant, and the paper is the demonstration.** It falls by about 25–30 % from the
   60–70 °C interval to the 120–130 °C interval at every water activity, which is what an Arrhenius
   barrier does. A model carrying a constant Q10 referenced at 25 °C is making an approximation this
   table can price.

## 6. Verdict

The second independent measurement of the lipid lane's missing temperature term, in a matrix rather
than in oil, agreeing with `frankel1993_extraction.md` at the dry end (114 kJ/mol CL against
Frankel's 113.8–122.2 kJ/mol for hexanal from n-6 seed oils at 130–160 °C) and falling well below it
once the matrix holds water. Two laboratories, two markers, two phases — and they do not agree on
one number, which is itself the finding. FIT-class evidence, not acted on in B37.
