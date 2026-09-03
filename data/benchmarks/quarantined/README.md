# Quarantined benchmarks

Files in this directory are **excluded from the validation panel**. The benchmark loader
(`panel_bundles` in `src/kinetic_core/panel.py`) uses a non-recursive
`benchmark_dir.glob("*.json")`, so any JSON placed in a subdirectory of `data/benchmarks/`
is invisible to every panel-wide generator, report and headline count. (The same mechanism
already keeps `data/benchmarks/external_validation/` out of the panel.)

Quarantine is the *first* step: the file is parked here so the original values and tolerance
contracts stay available for forensics and so the decision is reversible if a real source turns
up. It is not the resting place. When source recovery concludes that no source exists, the file
is **deleted** (git history still preserves it) and its section here is kept as the record — see
"Verdicts executed" below, where three of the four files listed have since been deleted.

## Why these files are here

A citation audit on **2026-08-26** could not locate a real literature source for any of them.

| File | Claimed DOI | Audit finding | Current state |
|---|---|---|---|
| `cys_ribose_150C_Mottram1994.json` | `10.1021/jf00049a015` | The DOI resolves to **Zarkadas et al. (1995)**, a beef-bone protein-quality paper containing **no volatile data at all**. No Mottram 1994 cysteine/ribose paper could be located. Nearest genuine candidate: **Mottram & Whitfield (1995), `10.1021/jf00052a027`** — a different system, different conditions, and its values were *not* checked against the numbers below. | **Deleted 2026-08-26** after source recovery confirmed no source exists (see §(a) below). Git history preserves the file. |
| `cys_glucose_150C_Farmer1999.json` | `10.1016/S0308-8146(98)00174-8` | Dead DOI. No Farmer 1999 cysteine/glucose source reporting these volatiles could be located under any DOI. | **Deleted 2026-08-26** after source recovery confirmed no source exists (see §(b) below). Git history preserves the file. |
| `thiamine_cys_ribose_100C_Hofmann1996.json` | `10.1021/jf960062o` | Dead DOI. The file reports **MFT = 3.14 ppb**, i.e. numerically pi — a confabulation tell rather than a measurement. | **Deleted 2026-08-26; replaced by `data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json`**, rebuilt from the verified Bolton et al. 1994 chapter per the recipe in §(c) below. |
| `acrylamide_asparagine_glucose_Parker2012.json` | `10.1021/jf3032779` | Dead DOI, no identifiable source; the real Parker 2012 acrylamide paper is a different system. **Quarantined 2026-08-26** (second pass) — see the dedicated section below. | Quarantined; file retained. |
| `spi_hvp_xylose_120C_PMC9905368.json` | `10.1007/s10068-022-01194-w` | **DOI is LIVE and correct — the paper is real, the content is not.** It uses glucose/fructose at pH 7.5 for 90 min and reports only *relative peak areas*; it never mentions FFT or MFT. The file's xylose/cysteine **absolute ppb** values have no possible source. **Quarantined 2026-08-27** (Wave I) — see the dedicated section below. | Quarantined; file retained. |
| `wheat_gluten_hvp_xylose_120C_PMC9905368.json` | `10.1007/s10068-022-01194-w` | Same source, same finding as the row above. | Quarantined; file retained. |

## Treat the values *and the tolerances* as suspect

Two things in these files are unverified, and the second matters as much as the first:

1. **The `measured_volatiles` values** have no traceable provenance. Do not cite them, do not
   fit against them, do not use them as priors.
2. **The `validation_contract.scale_thresholds` are unusually tight** — these three files
   carried three of the five tightest scale-tolerance contracts in the whole panel
   (`mean_abs_log10_error` of 0.06 / 0.10 / 0.18, versus 0.20–1.25 for most real literature
   benchmarks). A tolerance that narrow is normally *earned* by a well-characterised
   measurement; here it is far more likely a **fitting tell** — a tolerance chosen after the
   fact so the model would pass, on top of values chosen so the model could pass. Any
   "in-panel pass" these benchmarks previously contributed should be read as
   circular, not as evidence.

Note also `scripts/calibrate_barriers.py`, which fits barrier offsets against
`cys_ribose_150C_Mottram1994` and `cys_glucose_150C_Farmer1999`. Its output
(`data/lit/calibration_offsets.json`) is currently **unwired** and must stay that way — see
the warning in that script's module docstring.

## Un-quarantine procedure

Quarantine is reversible, but only by a human, and only in this order:

1. **Verify a real source.** Obtain the actual paper (not a search-result snippet, not a
   model-generated citation) and confirm it reports the system, conditions and analytes this
   benchmark claims.
2. **Fix the DOI** to the verified one, and record the real citation.
3. **Replace the values** with the ones actually reported in that paper. Do not keep the
   existing numbers unless the paper independently reproduces them.
4. **Re-derive the tolerances** from the reported analytical uncertainty of that paper, not
   from what the model currently predicts. If the paper does not support a tight contract,
   the contract must be loose.
5. **Remove** the `source_metadata.quarantined` / `quarantine_reason` / `quarantine_date`
   keys and `git mv` the file back up to `data/benchmarks/`.
6. **Regenerate** the tracked validation artifacts (see `scripts/generators/`) and update any
   test that pins panel counts or headline numbers.

If step 1 fails, the correct outcome is to **delete** the file, not to leave it here
indefinitely.

## Tests that load from this directory

**None, as of 2026-08-26 (second pass).** Earlier, six tests deliberately still pointed at the
three fabricated files to exercise reaction-network / SMIRKS coverage wiring. Because those
files are now deleted, every such test has been rehomed onto a fixture that does not depend on
an unverifiable file:

| Old test | Where its regression value now lives |
|---|---|
| `test_mottram_coverage.py` (deleted) | `tests/scientific/test_pentose_hexose_sulfur_ordering.py` — synthetic cys+ribose payload, coverage assertions only |
| `test_farmer_coverage.py` (deleted) | `tests/scientific/test_pentose_hexose_sulfur_ordering.py` — synthetic cys+glucose payload, coverage assertions only |
| `test_thiamine_fragmentation_benchmarks.py::…hofmann_1996…` | replaced by `…bolton_1994…`, running the rebuilt `data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json` |
| `test_free_aa_quantitative_regression.py` | the two quarantined entries were dropped; only the verified Hofmann1998 panel benchmark remains. Renamed 2026-08-27 (Wave H): the table is now `BENCHMARK_EXPECTED_FOLD_ERRORS` and holds two-sided pins on a quantified under-prediction rather than one-sided calibration ceilings |
| `test_lipid_oxidation_guard.py`, `test_time_propagation.py` | already retargeted (first pass) onto `cys_ribose_140C_Hofmann1998.json` |

`acrylamide_asparagine_glucose_Parker2012.json` is loaded by no test — its two audit-era
consumers were restructured when it was quarantined (see below).

---

## Source recovery findings (2026-08-26)

A follow-up literature search (CrossRef REST API, OpenAlex, PubMed/NCBI, publisher abstract
pages) went after the "nearest real candidate" leads recorded above. Every DOI below was
verified against CrossRef metadata (title + authors + volume/pages) unless flagged otherwise;
**no full text behind an ACS/Wiley paywall was obtained**, so where a number comes from an
indexed abstract rather than a table, that is stated.

**Headline: none of the three benchmarks can be un-quarantined as written. Two should be
deleted; one can be rebuilt from a real source, but only after changing its system,
conditions, analyte list and tolerance.**

### (a) `cys_ribose_150C_Mottram1994` — verdict: **DELETE the 3-analyte ppb triple**

What "Mottram 1994" most plausibly refers to is real, but it is the wrong paper for this
benchmark:

| Real paper | DOI | System | Data type |
|---|---|---|---|
| Mottram & Whitfield, *Aroma Volatiles from Meatlike Maillard Systems*, ACS Symp. Ser. 543 (Thermally Generated Flavors), 180–191 | `10.1021/bk-1994-0543.ch014` | 4-hydroxy-5-methyl-3(2H)-furanone **+ cysteine** (no ribose); pH-profile study | Qualitative — identification of mercaptoketones, thiols and novel disulfides; the abstract reports *"the pH of the reaction mixture had a major effect on the profile"*, no concentrations |
| Mottram & Whitfield, JAFC **43**:984–988 (1995) — the candidate named above | `10.1021/jf00052a027` | cysteine + ribose + **phospholipid**, **low moisture** | Not recoverable: no abstract in CrossRef, OpenAlex or Semantic Scholar (all three return the ACS boilerplate or null); closed access, no repository copy found |
| Mottram & Whitfield, *Maillard–Lipid Interactions in Nonaqueous Systems*, JAFC **43**:1302–1306 (1995) | `10.1021/jf00053a033` | cysteine + ribose + phosphatidylcholine, **nonaqueous** | ditto |
| Mottram & Nóbrega, JAFC (2002) | `10.1021/jf0200826` | cysteine + ribose / ribose-5-P / IMP, buffered and unbuffered, headspace GC-MS | **Comparative** quantities ("much lower", "similar quantities") — the abstract reports no absolute concentrations |
| Cerny & Davídek, JAFC (2003) | `10.1021/jf026123f` | cysteine + ribose/[13C5]ribose, phosphate buffer pH 5, 95 °C, 4 h | **13C-labelling only** — pathway assignment, no yields |
| Hofmann & Schieberle, JAFC **43**:2187–2194 (1995) | `10.1021/jf00056a042` | cysteine + ribose, **145 °C, 20 min** | **AEDA / FD factors only** |

Three separate reasons this benchmark cannot be rescued by re-pointing the DOI:

1. **Conditions don't match and can't be made to match.** The benchmark declares an aqueous
   system at `water_activity: 0.95`; the 1995 Mottram & Whitfield papers are explicitly
   *low-moisture* and *nonaqueous*, and both contain phospholipid, which the benchmark's
   precursor block does not. Repointing the DOI would attach a real citation to conditions
   that citation does not describe — the same failure mode, with better camouflage.
2. **The analyte list is not jointly measurable in one real source.** No located paper reports
   absolute concentrations for 2-methyl-3-furanthiol, bis(2-methyl-3-furyl) disulfide *and*
   furfural from one cysteine/ribose system. The quantitative cysteine+sugar literature of
   this era reports either **FD factors** (Hofmann & Schieberle) or **mol % yields**
   (`10.1021/jf9705983`), never a ppb triple.
3. **FD factors are not convertible to ppb.** An FD factor is a dimensionless stepwise
   dilution encoding concentration *divided by* that compound's odour threshold; the threshold
   is not recoverable from the paper. Any ppb derived from one would be manufactured.

**Recipe, if a cysteine/ribose anchor is wanted at all.** Do not rebuild this as a ppb
benchmark. Two honest options:

- *Ordinal (free, defensible today):* assert `MFT(pentose + cys) >> MFT(hexose + cys)`,
  citing Hofmann & Schieberle, JAFC **46**:235–241 (1998), `10.1021/jf9705983`, whose abstract
  states verbatim that *"pentoses generated much higher amounts of FFT and MFT than hexoses
  when heated in the presence of cysteine"*. That is a real constraint the model can be scored
  against, and it is exactly what the sources support.
- *Quantitative (needs ACS access):* Table 1 of the same paper carries per-precursor
  **mol % yields** for MFT and FFT at controlled T / pH / water content, convertible to ppb
  given the benchmark's precursor loading. This gives **one or two analytes, not three** — no
  furfural, no disulfide — and the yields span ~0.001–1.4 mol % across conditions, so the
  tolerance must be **order-of-magnitude** (a factor of 10 band), not the current
  `mean_abs_log10_error: 0.10`.

### (b) `cys_glucose_150C_Farmer1999` — verdict: **DELETE**

The dead DOI is not a typo. `10.1016/S0308-8146(98)00174-8` fails its check digit; the live
neighbours in that block (`…00174-5`, `…00173-3`, `…00175-7`) are a *dadih* lactose paper and
folate/lupin/jackfruit papers. An ISSN-scoped CrossRef author query returns **zero Linda J.
Farmer papers in Food Chemistry for 1996–2001**.

More decisively, **Farmer never published a cysteine + glucose model system.** Her full
model-system corpus, enumerated in her own 1994 review (Proc. Nutr. Soc. **53**:327–333,
`10.1079/pns19940038`, open access at Cambridge), is cysteine + **ribose** throughout:

| Real paper | DOI | Data type |
|---|---|---|
| Farmer, Mottram & Whitfield, JSFA **49**:347–368 (1989) | `10.1002/jsfa.2740490311` | Identification (GC-MS/GC-IR) + GC-effluent sensory; abstract states no quantities |
| Farmer & Mottram, JSFA **53**:505–525 (1990) | `10.1002/jsfa.2740530409` | ~60 products **compared across four lipids** — comparative, not absolute |
| Farmer & Mottram, JSFA **60**:489–497 (1992) | `10.1002/jsfa.2740600414` | comparative |
| Farmer & Patterson, Food Chem. **40**:201–205 (1991) | `10.1016/0308-8146(91)90103-u` | cooked beef, not a model system |
| Farmer, ACS Symp. Ser. 633:48–58 (1996) | `10.1021/bk-1996-0633.ch005` | review |
| Chevance & Farmer (1999), frankfurters | `10.1021/jf990515d`, `10.1021/jf9905166` | real foods, not model systems |

Four chemistry tells independently mark the numbers as invented, so this is not merely an
unlocatable source:

1. **MFT from glucose is the wrong sugar** — MFT is a pentose/ribose marker
   (`10.1021/jf9705983`, quoted above).
2. **Furfural at 450 ppb as the dominant product of a hexose system is wrong-signed** —
   furfural is the 3-deoxy*pentos*one product; glucose gives HMF and 5-methylfurfural.
3. **"Pyrazine" is not a real analyte** — the literature reports specific alkylpyrazines, and
   cysteine is a poor pyrazine former relative to lysine/alanine.
4. **`water_activity: 0.95` with 10 mM solutes is self-contradictory** — a 10 mM aqueous
   solution has a_w ≈ 1.00.

The genuine glucose/cysteine paper is Hofmann & Schieberle, JAFC **45**:898–906 (1997),
`10.1021/jf960456t` (20 min at 145 °C). Its top odorants are 2-furfurylthiol, 5-acetyl-2,3-
dihydro-1,4-thiazine, 3-mercapto-2-butanone, 3-mercapto-2-pentanone and Furaneol — **MFT does
not appear among them**, contradicting the benchmark directly — and it reports FD factors, not
ppb. There is no replacement source. **Delete the file rather than re-point it.**

### (c) `thiamine_cys_ribose_100C_Hofmann1996` — verdict: **REBUILDABLE from a different paper, at different conditions**

The claimed DOI is confirmed nonexistent, not merely paywalled: CrossRef `/works/10.1021/jf960062o`
and `/works/…/agency` both return *Resource not found* (unregistered with any DOI agency), and
`https://doi.org/10.1021/jf960062o` returns **404** where every other ACS DOI tested returns a
Cloudflare **403** — a clean discriminator between "blocked" and "does not exist". The only
Hofmann paper in JAFC in 1996 is Hofmann, Schieberle & Grosch, *Model Studies on the Oxidative
Stability of Odor-Active Thiols*, **44**:251–255, `10.1021/jf9500703`, which measures thiol
*degradation*, not formation, and so cannot support an MFT-formation benchmark.

Real thiamine literature, verified:

| Paper | DOI | What it gives |
|---|---|---|
| **Bolton, Reineccius, Liardon & Huynh-Ba**, *Role of Cysteine in the Formation of 2-Methyl-3-furanthiol in a Thiamine–Cysteine Model System*, ACS Symp. Ser. **543**:270–278 (1994) | `10.1021/bk-1994-0543.ch022` | **The only quantitative thiamine+cysteine MFT measurement found.** Water + MSG + NaCl + IMP + D-glucose + thiamine·HCl + cysteine·HCl (± D-xylose), **120 °C, 1 h**, sealed vials, 4 variants; MFT quantified by GC-MS-SIM against an internal standard: **389–489 ng per 33.3 g = 11.7–14.7 ppb**. No thiamine → no MFT. Only 7.5–8.0 % of MFT sulfur came from 34S-cysteine, i.e. **thiamine is the dominant S source** |
| Hofmann & Schieberle, JAFC **46**:235–241 (1998) | `10.1021/jf9705983` | mol % yields via SIDA; states verbatim that *"thiamin and norfuraneol/cysteine were less effective precursors of MFT"* |
| Güntert et al., *Thermally Degraded Thiamin*, ACS Symp. Ser. **490**:140–163 (1992) | `10.1021/bk-1992-0490.ch011` | thiamin at pH 1.5 / 7.0 / 9.5 → 38 / 32 / 59 compounds; **qualitative** |
| Cerny & Briffod, JAFC **55**:1552–1556 (2007) | `10.1021/jf062874w` | [13C5]xylose + cysteine + thiamin, pH 4.0–7.0; MFT draws carbon from both — **relative, not ppb** |
| Kerscher & Grosch, JAFC **46**:1954–1958 (1998) | `10.1021/jf970892v` | reality anchor: **MFT 5–28 µg/kg** in real boiled meat (beef 7–28, pork 6–9, lamb 5–11) |
| Dreher, Rouseff & Naim, JAFC **51**:3097–3102 (2003) | `10.1021/jf034023j` | thiamin 0.024 mM in model orange juice, pH 3.8, **35 °C**; MFT peaks at **112 µg/L** — treat as soft: it implies ~4 mol % conversion, which sits awkwardly against Hofmann & Schieberle's "less effective precursor", and is at/above the top calibration standard |

**Un-quarantine recipe (the only one of the three that has one):**

- **Source:** Bolton, T. A.; Reineccius, G. A.; Liardon, R.; Huynh-Ba, T. *ACS Symp. Ser.*
  **1994**, *543*, 270–278. `10.1021/bk-1994-0543.ch022`
- **Conditions — these must change:** 120 °C / 60 min (not 100 °C / 30 min); precursors become
  thiamine·HCl + cysteine·HCl + **D-glucose** (+ IMP, MSG, NaCl), not thiamine + cysteine +
  **ribose**. If the benchmark specifically needs ribose, there is **no** quantitative source
  and the file should be deleted instead of substituted.
- **pH:** not stated in the indexed abstract. Leave `ph` unpinned, or obtain the chapter full
  text before pinning a value. **Do not invent one** — an invented pH is how the file got here.
- **Target:** MFT = **13 ppb** (midpoint of the reported 11.7–14.7).
- **Tolerance:** pass within a **factor of 3** (≈4–40 ppb), i.e. `max_ratio: 3.0`,
  `mean_abs_log10_error: 0.48`. Justification: the published spread across model variants is
  only ±13 %, but the unpinned pH, the ribose↔glucose substitution and the 120 vs 100 °C gap
  each move MFT by roughly 2–3×. A factor-of-2 band would be over-claiming. The existing
  `0.18` contract is not defensible under any reading.
- **Independent corroboration:** 11.7–14.7 ppb sits inside Kerscher & Grosch's 5–28 µg/kg for
  real boiled meat — two independent methods agreeing within ~2×.
- **Caveat to carry into the file's `notes`:** the composition, temperature, time and the
  389–489 ng/33.3 g figure come from the **indexed abstract**, not from the chapter's tables.
  Someone with ACS access should confirm the pH and exact molarities before these are pinned
  as ground truth.

### Also settled by this pass

- `bk-1994-0543` (ACS Symp. Ser. 543, *Thermally Generated Flavors*) contains **both** the real
  "Mottram 1994" (ch. 14) and the Bolton thiamine chapter (ch. 22). If the original synthetic
  citations were seeded from a reading list, this volume is a plausible common ancestor of the
  "Mottram1994" and "Hofmann1996" labels.
- `scripts/calibrate_barriers.py` fits barrier offsets against (a) and (b). Both are now
  recommended for **deletion**, which means that script's output
  (`data/lit/calibration_offsets.json`) is fitted against two sources that do not exist. It is
  currently unwired; the recovery result is that it must **stay** unwired, and the script
  should be retargeted or removed rather than re-enabled.

---

## Verdicts executed (2026-08-26, owner-approved)

The source-recovery verdicts above were approved and applied on the same date. Record:

| File | Verdict | Action taken |
|---|---|---|
| `cys_ribose_150C_Mottram1994.json` | delete | **Deleted 2026-08-26 after source-recovery confirmed no source exists.** `git rm`; history preserves the values. Its lost chemistry regression is replaced by the ordinal pentose ≫ hexose MFT test (below). |
| `cys_glucose_150C_Farmer1999.json` | delete, no replacement | **Deleted 2026-08-26 after source-recovery confirmed no source exists.** `git rm`; history preserves the values. |
| `thiamine_cys_ribose_100C_Hofmann1996.json` | rebuild from Bolton 1994 | **Deleted 2026-08-26; replaced by `thiamine_cys_glucose_120C_Bolton1994`** (see below). |

### Replacement for the deleted cysteine/sugar chemistry regression

`tests/scientific/test_pentose_hexose_sulfur_ordering.py` encodes the one constraint the real
literature does support: **pentose-derived MFT ≫ hexose-derived MFT under matched conditions**,
citing Hofmann & Schieberle, JAFC **46**:235–241 (1998), `10.1021/jf9705983` ("pentoses
generated much higher amounts of FFT and MFT than hexoses when heated in the presence of
cysteine"). It runs the model on two synthetic formulations — cysteine + ribose and cysteine +
glucose at identical T / pH / a_w / time / molarity — and asserts an ordering with a margin,
plus the SMIRKS-coverage assertions inherited from the deleted `test_mottram_coverage.py` /
`test_farmer_coverage.py`. No unverifiable reference number is involved.

### The rebuilt thiamine benchmark

`data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json` follows the un-quarantine recipe in
§(c) verbatim: thiamine·HCl + cysteine·HCl + D-glucose, 120 °C / 60 min, GC-MS-SIM, target
MFT **13 ppb**, tolerance `max_ratio 3.0` / `mean_abs_log10_error 0.48`, DOI
`10.1021/bk-1994-0543.ch022`. Two things the source does **not** give are marked as such in the
file's `source_metadata`: the **pH** (5.0 is a loader-required placeholder, explicitly not
anchored) and the **molarities** (`assumed_not_from_source`).

**It fails, by a lot, and that is the point.** The model predicts **0.080 ppb** against the
13 ppb target — a **163×** under-prediction, `mean |log10 ratio|` 2.21. A sensitivity sweep over
the assumed values shows the failure is not an artefact of the assumptions: ratios stay between
**35× and 745×** across thiamine 1–30 mM, cysteine/glucose 10–50 mM and pH 4–6. This is the
first genuinely failing thiamine row in the panel, and it replaces a row (`3.14 ppb`, tolerance
`0.18`) that passed because both the value and the tolerance were chosen to make it pass.

---

## `acrylamide_asparagine_glucose_Parker2012` — quarantined 2026-08-26 (second pass)

Quarantined on the same evidence pattern as the three files above, after the full citation
sweep (`results/validation/citation_verification_ledger.md`).

| Aspect | Finding |
|---|---|
| Claimed DOI | `10.1021/jf3032779` — **dead**, does not resolve in CrossRef (HTTP 404) |
| Nearest real paper | Parker et al. (2012), `10.1021/jf302415n` — a kinetic model of acrylamide in **commercial French fries during finish frying**. Not a 10 mM asparagine/glucose aqueous model system: different matrix, different precursor basis, different measurement. Re-pointing the DOI would attach a real citation to conditions that citation does not describe. |
| Tolerance | `max_ratio: 1.05`, `mean_abs_log10_error: 0.02` — the **tightest contract in the entire collection**, tighter than all three files quarantined in the first pass. |
| Observed agreement | predicted 1523.98 ppb vs "measured" 1500 ppb, ratio **1.016** — i.e. the model landed inside a ±5 % band on a value with no source. A 1.6 % agreement on an unsourced number is the fitting tell, not a validation. |
| `water_activity` | 0.30 declared for a 10 mM aqueous solution — the same self-contradiction flagged for `cys_glucose_150C_Farmer1999`. |

**Panel consequence.** This was the last surviving **hexose** benchmark in the panel and the
only free-precursor acrylamide anchor; `acrylamide_spi_extrusion_130C_ACSRef3` remains as the
matrix-path acrylamide reference, and `src/safety.py`'s acrylamide kinetics are still covered by
`tests/scientific/test_safety_benchmarks.py`, which asserts directional behaviour rather than
agreement with this file.

**Panel size is unchanged at 16** (−Parker2012, +Bolton1994). The *quality* of the panel is
not unchanged: the swap removes one passing, strict-ready row (Parker's acrylamide, ratio
1.016 against a 1.05 contract) and adds one failing row (Bolton's MFT, ratio 163 against a
3.0 contract). Do **not** carry the old 35/41 and 8/16 headlines forward with a hand
adjustment — other work in flight on this branch (the sulfur-barrier re-centring and the
Sander 5.0.0 Kaw refresh) moves scored rows too. Recompute the headline from a fresh
regeneration rather than predicting it.

---

## `spi_hvp_xylose_120C_PMC9905368` + `wheat_gluten_hvp_xylose_120C_PMC9905368` — quarantined 2026-08-27 (Wave I)

Found by the **cold-start red team** (2026-08-27), independently by both reviewers (forensic
finding F2, scientific finding C1). These two are a *different* failure mode from everything
above, and the difference is the lesson: **their DOI is live, resolves correctly, and the
metadata matches.** Every previous check — the full 225-DOI citation sweep included — passed
them, because those checks verified that the citation *exists*, not that the paper *says what
the file claims*.

### The evidence

The cited source is:

> Food Science and Biotechnology, `10.1007/s10068-022-01194-w` (PMC9905368).

| | What the paper actually did | What these files claim |
|---|---|---|
| Reducing sugar | **D-glucose and D-fructose** | **D-xylose**, 50 mM |
| pH | **7.5** | **6.0** |
| Time | **90 min** | **30 min** |
| Reported quantity | **Relative GC-MS peak areas (%)** — no internal standard, no absolute quantitation | **Absolute `conc_ppb`**, each with a per-analyte `uncertainty_pct` |
| 2-Furfurylthiol (FFT) | **not mentioned anywhere in the paper** | 0.42 / 0.61 ppb ± 12 % |
| 2-Methyl-3-furanthiol (MFT) | **not mentioned anywhere in the paper** | 0.18 / 0.34 ppb ± 12 % |
| Methional | present as a relative peak area | 1.88 / 3.44 ppb ± 10 % |

The decisive point is not that the conditions differ — conditions can be mis-transcribed. It is
that **a paper reporting only relative peak areas cannot be the source of an absolute ppb value
by any route**: there is no unit conversion, no internal standard, and no calibration curve to
convert from. For FFT and MFT there is not even a peak area, because the paper does not report
those compounds at all. **These six numbers have no possible provenance.**

### What they were doing in the repo

Worse than passively sitting in the panel. They were:

1. **Six of the ~15 scored literature rows** in the coverage headline — roughly 40 % of the
   panel's external-literature evidence base.
2. **The sole fit targets** of `scripts/generators/rederive_hydrolysate_observability.py`
   (Wave H, 2026-08-27), which moved the Methional `base_factor` in
   `src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` from 0.0045 to 0.05623 —
   **and then the two Methional rows it had just been fitted to were scored as literature
   coverage.** They were the *only two hits* in the Wave H headline of "2/11 literature rows
   inside the CI". Fitted to a benchmark, then reported as agreement with that benchmark, on
   data that never existed. The honest count for that headline is **0/11**.
3. **Six of 15 rows** in the fit population of
   `scripts/generators/refit_projection_constants.py`.

### Actions taken (2026-08-27)

| Action | Where |
|---|---|
| Files `git mv`'d here with `source_metadata.quarantined` + full evidence block | this directory |
| Methional `base_factor` **reverted** 0.05623 → 0.0045 | `src/recommend.py` |
| Re-derivation record marked **RETRACTED** (json + md), generator refuses to run | `results/validation/hydrolysate_observability_rederivation.*`, `scripts/generators/rederive_hydrolysate_observability.py` |
| Contamination note added; the two files removed from the fit population by the quarantine itself | `scripts/generators/refit_projection_constants.py` |
| Contamination review added, verdict **uncontaminated, result stands** | `scripts/generators/refit_sulfur_barriers_hofmann.py` |
| New blocking CI gate against fit-then-score circularity | `scripts/ci/fit_target_gate.py` |

### Does the sulfur barrier refit survive?

**Yes — verified, not assumed.** `scripts/generators/refit_sulfur_barriers_hofmann.py` fits
against exactly one benchmark, `cys_ribose_140C_Hofmann1998`
(`10.1021/jf9705983` — Hofmann & Schieberle, JAFC **46**:235–241, a real and verified paper).
The two quarantined files were on that script's *forbidden* list before the quarantine and
contributed **zero rows** to its objective. So `thiol_addition_norfuraneol` 28.60 → 26.85 and
every "keep the incumbent" decision in that record stand unchanged.

One caveat that is **not** about the quarantine but should not be lost: Hofmann1998's own
contract (`max_ratio` 1.45 / `mean_abs_log10_error` 0.09) is the third-tightest in the
collection — the same tightness pattern flagged as a fitting tell above. The tolerance has not
been widened and the panel currently *fails* it, so nothing is concealed; but a single
benchmark with a suspiciously tight contract is a thin anchor for an entire branch of the
chemistry, and that remains an open owner item.

### What does NOT survive

- The Methional `base_factor` of 0.05623, and the re-derivation record that produced it.
- Any headline claiming those Methional rows as literature agreement.
- With the xylose HVP lanes gone, **every** entry in
  `_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` is now an unconstrained legacy estimate with no
  literature constraint of any kind. That is stated at the constant.

### The generalisable lesson

The citation gate checks DOI *identity*. It cannot check DOI *content* — no offline gate can.
Every anchor in this repo that has not been read against its source is exposed to exactly this
failure mode, and the population of such anchors is large. Live-DOI-wrong-paper is not a rarer
class than dead-DOI; the 2026-08-26 sweep measured it as **equally large** (47 metadata
mismatches + 45 topic mismatches vs 45 dead).
