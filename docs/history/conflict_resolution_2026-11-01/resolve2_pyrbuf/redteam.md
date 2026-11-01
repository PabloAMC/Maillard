# Red team: resolve2_pyrbuf + Wolf 1977a (2026-11-01)

Scope: `resolve2_pyrbuf/{prereg,results,register_edits.yml,doc_corrections.md,tests.py,shib_abs.py,piloty_engine_probe.*}`,
`wolf1977a_findings.md` (+ scripts), dossier `wolf1977a_extraction.md`. Checked against the PDFs (larralde1995, piloty1979,
wolf1977, wolf1977a: text layers plus the eq-3 page crop), the dossiers, `split_v2_per_change.yml`, `split.yml`,
`contrast_benchmark/README.md`, `F04`/`F15` atlas files and `conflicts_register.yml`. My scripts are in this folder:
`potman_power.py` (T1 power at run-to-run noise) and `shib_incremental.py` (Shibamoto with rate-based pool proxies). Both
ran through docker_maillard.sh, niced, one process.

## Verdicts

### 1. E5: buffer catalysis is linear per step (P 0.8). WEAKENED
- **Split status.** Potman, Kallen and Thornalley are TEST for E5 under v2. They were DEV in v1, except Potman. **Larralde
  is not TEST under any frozen split.** README sec. 10 lists `heldout_larralde1995.yml` as a keyed held-out table whose lab
  joins TEST only "if a new version admits these labs". Lee 2024 (E13) has the same status. Every "TEST" label on
  Larralde or Lee must read "keyed held-out, candidate TEST for v3".
- **The T1 form test stands statistically.** The dossier's Fig. 4 vs Fig. 7 run-to-run agreement is 3-8 % points, not 1.5.
  At those noise levels the rule still calls F-lin non-saturating 94 % (sd 5) and 85 % (sd 8) of the time. It wrongly
  calls F-sat with K 0.1 M non-saturating 0.7 % and 8 % of the time. So "power 1.0" should read "≥ 0.85".
- **What T1 actually rejects.** It rejects K ≤ ~0.3 M. With K = 0.5 M the maximum ratio is 3.0, which sits on the rule's
  threshold.
- **Systematic confounds the prereg did not list:**
  - Phosphate at pH 5.6 is a weak buffer, 1.6 units below pKa2, so buffer capacity scales with concentration. Potman's pH
    drift is unreported but admitted, and the 0.4 / 0.5 M arms plateau after 3 h. The 0.1 M arm is the most exposed to
    drift.
  - Ionic strength co-varies with buffer in Potman, Bell and Guan. Only Kallen held I at 1.0 M (KCl).
  - Neither confound can plausibly give ×3.4 alone, so the F-sat rejection survives. It does mean "linear" itself is not
    shown.
- **The order is ≥ 1, not exactly 1.** Potman's initial-rate order is 1.34. Its 3-h secant gives κ 48-81 /M, which is
  outside the DEV band, while the 5-h secant gives 22-29. Thornalley is convex. Relative κ is also k0-limited: the 5-h k0
  carries ±24 % and the 3-h k0 ±60 % from read error alone. T2's Potman pass depends on the secant window.
- **Missed ceiling inside Potman's own data.** Series A has phosphate at 1.2 M in a 30-40 % dry-matter pot. It gains only
  ×11.6, i.e. κ ≈ 9 /M, where linear κ 25 predicts about ×31: 0.43 dex short. The Potman dossier itself warns that a linear
  term "will over-predict catalysis in concentrated matrices by ≥ 5×". Wolf's pots carry about 1-2 M buffer in the water
  phase. The remedy's "no Michaelis term" ignores this.
- **F-net at the engine level.** `r_schiff` lumps the amine attack and the dehydration. By Kallen's own F-net algebra, a
  lumped step whose attack can become rate-limiting has a hyperbolic (Michaelis-shaped) factor in c. Its half-saturation is
  set by κ and the ceiling rather than by binding. "No Michaelis term" is therefore right only for steps the network
  resolves, and the engine resolves no uncatalysed attack for glucose + amine.
- **T4 Kallen uses the superseded ceilings.** It used ×5 / 3.5 / 4.3. The two-reader re-read (F04 atlas note; README
  "MINOR") gives C3/C2 = 4.43 / 3.06 / 4.06. Re-scored, the residuals are +0.09 / -0.11 / -0.22 dex, which still pass.
  - The predicted pH trend (×3.6 → ×6.8, rising) is not seen: the ceilings run 4.4 / 3.1 / 4.1. The test checks level
    consistency only, and F-sat fits Kallen equally well.
  - CF-076 already uses Kallen's saturation to reject a linear catalysed exchange. CF-NEW-1 must cross-reference it.
- **T3 Larralde: S-P is not rejected.** Larralde's own ratio precision is kH2O/kD2O = 1.0 ± 0.2. So 0.52 ± 0.10 against
  0.38 (0.33-0.43) is 0.9-1.4 σ. CO2 loss in the sealed tube would raise pD and bias the ratio upward. "Bicarbonate 0.23×
  phosphate" is consistent with zero.
- **Bands contradict their own data.**
  - The prereg set the carboxylate band at [0, 1.5] /M on the premise that Guan's acetate and citrate equal water within
    ±15 %. In fact they are +16 % and +34 %, giving κ 1.6 and 3.4.
  - The results' band of [0, 3] /M (B1, CF-NEW-2) excludes Guan's citrate value, 3.4, which is the band's own source.
  - "≤ 0.1× phosphate" holds at 110 C only. At 130 C phosphate κ is 5-16, so the same carboxylate bound is ≤ 0.3×.
- **"κ roughly T-flat 25-110 C" (B1) is unsupported.** It comes from cross-lab, cross-step data. The only within-lab
  temperature contrast, Guan 110 → 130 C, falls 0.4-0.75 dex (28/38 → 5/16).
- **Davidek pH 8.** The catalysed increment itself falls 0.55 → 0.29 /h, so about 40 % of the ratio fade is lost catalysis.
  Only about 60 % comes from k_u rising, so "the fade comes from k_u" is overstated.
- **Pooling and independence.** No pyrophosphate or McIlvaine buffer is pooled anywhere in the band. Guan's null is
  single-lab.
- **Corrected probabilities:**

  | claim | P |
  |---|---|
  | non-saturating with K ≤ 0.3 M in dilute pots | 0.8 |
  | "linear, no Michaelis term" as an engine law for lumped steps and concentrated foods | 0.45 |
  | S-step | 0.5 |
  | B-chain | 0.7 (upheld) |

### 2. Piloty does not discriminate E4 (P 0.85). UPHELD. CF-NEW-3 is real, not a units artefact
- Checked against the PDF: Table 1 gives MGO 312 / amino acids 208 mM total, i.e. 16 mM each across 13 amino acids. Table 4
  gives glycine loss 85 / 86 / 82 %.
- Units are consistent. k_mgo_ak is 3e-8 L/(mmol·min), about 5e-7 M⁻¹s⁻¹. The implied consumption constant is ≥ 3.5e-5
  L/(mmol·min), about 6e-4 M⁻¹s⁻¹. The probe's 5e-4 loss equals k × ∫MGO, with MGO cut 312 → 5 mM by r_mgo_mel.
- The analyser hydrolyses reversible imines, so the measured loss is irreversible consumption. The deficit is real.
- "~1e3×" is a lower bound: it assumes constant MGO. Sedehi's MGO + Gly gap in its own dossier is ~2×10⁴, not ~1e3.
  CF-NEW-3 cites only the glyoxal 1e3 figure.
- `kind: code_mismatch` is wrong: no source is mis-transcribed. Use `conflict` (or `magnitude`). Drop `r_mgo_ak` from
  `affects`, or tag it "not the error", so that nobody "fixes" the Strecker constant.

### 3. Pyrazines. P0 rejection UPHELD (3 labs, not 5). P2b rejection REFUTED as a test. Lee reading WRONG
- **Lab count.** Shibamoto 1976, Wong 1988 and Bemis-Young 1993 all come from Bernhard's UC Davis lab. So "5 labs" is 3
  independent labs: UC Davis, Koehler and Lee. P0 is still rejected.
- **The pool-invariance proxy is mis-specified.** `shib_abs.py` uses √(cumulative amount) as the instantaneous aminoketone
  pool, but cumulative yield grows even at a constant pool.
  - With rate-based proxies, √(Δamount/Δt), the pools change only 1.1-1.4× between 0.5 and 2 h. P2b then predicts -0.09 to
    -0.37 dex against 0.00 observed, which is |z| ≲ 1.8 and no rejection.
  - After 2 h every species shows net loss (yield 1.25 → 0.57 %), so the 3-48 h points that gave |z| 3-6 carry no
    pool information.
- **Koehler** is one arm that also changes the chemistry.
- **Correction:** "P2b lower-order route not tested (P back to ~0.45)".
- **Lee 2024 face value is not a violation.** Under the prereg's own P1x algebra, R = 4r² has no upper bound; Shibamoto's
  R ≈ 400 passes P1. "R ≤ 4" is P0's bound. Face value (15-26) is P1x-compatible (r ≈ 2-2.5). Only MPZ/10 (0.15-0.26)
  breaks CF-087's R ≥ 0.4 floor. The ×10 issue therefore decides whether any bound fails. Rewrite CF-NEW-5 and the results
  bullet.
- **Silent omissions.** The prereg named E13 TEST contrasts shibamoto1977, low2007, vanlancker2012 and amranihemaimi1995.
  None is scored or excused in the results, against CF-087's "no silent omission" rule.
- **The "160 kJ/mol 25-60 C" figure** rests on 30-day pots at 0.005-0.007 % yield, where CF-087's low-level caveat applies.

### 4. Huber / a_w / Wolf 1977. UPHELD
- **Huber.** The oxidation-level argument checks out: C-1 gains an N and C-4 loses an O, net zero, i.e. the 1-DG level. So
  F15 → F07 is right. Note that F15's 1,4-dideoxyosone (C-1 methyl) is still a net reduction product, so F15's definition
  stays.
- **CF-NEW-7.** Its kind `condition_error` fits poorly for a family mis-assignment.
- **Cui and Wolf 1977.** Both agree with their dossiers.

### 5. Wolf 1977a
- **CF-072 undecided (P 0.1): UPHELD.** Table 2 is at 105 C ("processed at 105 °C") and Table 3 at 95 C. The values
  match the PDF.
- **E centre 32.9 kcal/mol (137.7 kJ/mol): UPHELD** on the page image. Table 4 gives SE ≈ 2 kcal/mol, so the engine's Ea
  96.8 sits about 5 SE below it. The two-temperature end-point design with a plateau is a real caveat; keep it.
- **"Weaker glucose dependence": REFUTED as a conflict claim.** Case 1 gives ×2.3 / ×6.6, but the authors' text says
  ×5.3 at 80 C and ×25 at 130 C (orders 0.7 and 1.3). The engine's order 1 sits inside that.
- **"Case 3 does not reproduce (19.9 vs 32.9)": REFUTED.** The printed "=" before "3.3*a_w" is a misprint. Read as "+",
  Case 3's centre is 33.1 kcal/mol and k₁₀₀ 0.0122, matching Cases 1 and 2. Read as "−", it is 19.9. P("+") ≈ 0.9.
- **"Unsourced 5.0 ± 0.7": WEAKENED.**
  - The companion prints it in Table 1, footnote c, "Wolf et al. (1977a)". Its own reference list maps 1977a to J Food
    Sci 42:1540, so the dossier is faithful.
  - 1977a has no such value. The companion mis-letters 1977a/b elsewhere: p. 271 attributes extrusion to 1977b, p. 277 to
    1977a. So this is a mis-citation in the primary, probably 1977b (submitted).
  - Suggested flag text: `secondary_source_misattributed`.

### 6. Register hygiene
- **Ids.** Assign CF-093 onward, since CF-092 is taken by the other batch. Replace "CF-XXX".
- **Overlaps.**
  - CF-NEW-4's remedy ("download Wolf 1977a") is stale. Fold in the 1977a results. Its plateau depth moves with glucose
    (0.43 vs 0.69, confounded), which lowers L-pool (P 0.5 → ~0.35).
  - Reconcile CF-NEW-4 with the CF-072 appends.
- **No proposal moves the engine directly.** B1 writes numeric κ bands and "divide calibration constants" as the Fix. Mark
  them A3-provisional, pending their own prereg.

## REQUIRED before applying
1. Relabel Larralde and Lee as "keyed held-out (not in a frozen split)" everywhere: prereg, results, CF-NEW-1/2/5 and B1.
2. CF-NEW-1 and B1:
   - say "non-saturating to 0.5 M in dilute pots (order 1.0-1.35)", not "linear";
   - add Potman Series A's 1.2 M shortfall (0.43 dex) and the lumped-step F-net caveat;
   - restrict "no Michaelis term" to steps the network resolves;
   - use Kallen ceilings ×3.1-4.4 and re-score T4 (+0.09 / -0.11 / -0.22);
   - cross-reference CF-076;
   - delete "T-flat"; replace "power 1.0" with "≥ 0.85 at run-to-run noise".
3. CF-NEW-2, B1 and A3:
   - Larralde's bicarbonate result becomes "consistent with zero (±0.2 ratio precision)";
   - widen the carboxylate band on Amadori decay to [0, ~4] /M, or state that citrate is 3.4;
   - write "≤ 0.1× at 110 C, ≤ 0.3× at 130 C";
   - drop P(S-step) to 0.5.
4. CF-011 append: delete "Kallen ... puts any catalysis in Cao's pot on the Amadori rearrangement". Kallen's
   no-catalysis-below-pH-5 test used carboxylates, not phosphate. By microscopic reversibility the catalysed dehydration
   also implies a catalysed ring-opening path.
5. Q3:
   - drop "P2b rejected by pool invariance (P 0.75)" from the results, the EX-NEW conclusion and C8;
   - fix "5 labs" → 3;
   - rewrite CF-NEW-5 / Lee: the face value violates only P0;
   - score or excuse the four unreported E13 TEST contrasts.
6. CF-NEW-3: change the kind; add Sedehi's MGO ~2×10⁴; drop or annotate `r_mgo_ak`; call 1e3 a lower bound.
7. Wolf 1977a:
   - drop the glucose-order clause from the new Ea conflict;
   - correct findings item 7 and the ledger flag (Case 3 consistent under "+");
   - reword the 5.0 ± 0.7 flag as a mis-citation in the companion.
8. Merge CF-NEW-4 with the Wolf 1977a results; assign CF ids after CF-092.

## Optional
- Report κ with k0 propagated (±24 % at 5 h). Note Davidek's pH-8 increment loss (about 40 % of the fade).
- Add CF-NEW-6's sibling shibamoto1977 (same lab and method, R 322-631) to its sources.
- Note the Shibamoto 0.5 / 1 / 2 h R identity (365.6 / 365.9 / 365.3) as a curiosity worth a PDF re-read.
