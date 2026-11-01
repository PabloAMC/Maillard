# Conflict triage: which of CF-001 to CF-076 can be resolved, and how (2026-11-01, read-only)

Inputs read: `docs/reference/conflicts_register.yml` (76 entries), `docs/reference/conflict_analyses/CF-011..017.md`,
`docs/guides/STATUS_AND_PLAN.md` sec. 4-5, `docs/guides/EVIDENCE_MAP_2026-10-31.md` (sec. 1.9, 2.11, 2.12),
`docs/history/{thiol_rivals,pyrazine_rivals,cys_exchange,mechanism_reviews}_2026-10-31/`, the measurement ledgers, the
keyed `data/lit/fit_rows/*.yml`, and spot reads of the dossiers named below. Nothing was edited or run.

Codes for cause classes:
- **OP**: observation operator / analytical method.
- **COND**: condition differences.
- **TU**: transcription, units, duplicate data, citation.
- **CHEM**: genuine chemistry.
- **MOD**: model artefact.

Importance (0-3) follows the brief: 3 means it decides a structural choice for Strecker aldehydes, pyrazines,
MFT / FFT / H2S or the hub. P is P(resolvable from public data on disk or in the queue). Hours are work hours.

---

## 0. Register defects found while triaging (fix first; about 1 h in total)

| item | defect | fix |
|---|---|---|
| CF-052 / CF-072 | The same finding (Nguyen 2016 returns no lysine; B20 `r_flp_decay` writes LYSP back) is registered twice. | Merge into CF-072 (the triage's #2). |
| CF-065 | The Zhai source is **Zhai 2023 JAFC Fig. 6**, not Food Chem: furanoid + **(NH4)2S** 10 mM, pH 7, 120 C, g/L axis, semi-quantitative (`data/lit/source_verification/thiol_submodel.yml` ~l. 470-490). The Food Chem dossier has no fed-sulfide figure. "H2S" hides an NH3 co-reactant. | Correct the source and the sulfide. |
| CF-012 | Still cites `li2011` / `hou2017`. The CF-012 analysis shows both are the wrong papers (the real leads are unread). | Replace the two sources, or drop them. |
| CF-011 | Explanation 2 still carries "Zhai 2019 ... falls to 4.9". The CF-011 analysis sec. 1: Zhai 2019 reports no pH after heating; the 7.0 -> 4.9 figure belongs to Zhai 2023 Food Chem. | Rewrite explanation 2. |
| CF-019 | Says zhou2025b re-publishes zhou2024. CF-060 shows it is not a re-plot. The Leahy pair is already settled: leahy1989a's pH 9 rows are the ch. 7 data reprinted (leahy1989a dossier sec. 0). | Update the entry. |
| CF-055 | Status `open`, but both readers and the adjudicator have confirmed it (`thiol_submodel.yml` l. 138). | Set `error_confirmed`. |
| missing entries | The pyrazine factor-2 bug, which the pyrazine README calls "CF-077", is not in the register. Also unregistered: the zhou2023c ARP + Gly row (MPZ not detected, which every pairing rival forbids; pyrazine red team sec. 5); Liu 2023's EI [M-H]+ artefact (thiol review sec. 0); the prereg sec. 13 statement that Zhai 2021's pots are "0.2 M phosphate", when they are water (hub and cysteine reviews). | Register CF-077 to CF-080. |

---

## 1. Per-conflict assessment (all 76)

Status column: "-> ok" means the register's status is still right; "->" followed by a status is the status I recommend.

| CF | (a) importance: mechanism / odorant | (b) status | (c) causes, ranked | (d) how to resolve without a lab; P; hours |
|---|---|---|---|---|
| 001 | 2: MFT pH law (Meynier hold-out) | error_confirmed -> ok; remedy not applied | TU (conditions copied from the Hofmann pot) | Resolved; apply at the next hold-out generation. P 1; 0.5 h |
| 002 | 1: thiol pH law on fed TTCA (S5 unused) | error_confirmed -> ok (the temperature is still unknown) | TU (caption or pasted column) > COND | Compare S5's pH 5.5 / 8 columns with Zhai 2021 Fig. 4b/5 (100 C) as thiol yield per TTCA lost; ask the author. P 0.3; 2 h |
| 003 | 0 | error_confirmed -> ok | TU | Labels only. P 1; 0.5 h |
| 004 | 3: absolute MFT / FFT level in every objective since B8 | error_confirmed -> ok | OP (f' = 1, semi-quantitative) | Resolved. Carry `semi_quantitative` into every successor (the thiol prereg item 6 already sets tau = inf). P 1; 1 h |
| 005 | 2: the ladder counted twice | duplicate_confirmed -> ok | TU (one experiment printed four times) | Apply count_once. P 1; 0.5 h |
| 006 | 1 | duplicate_confirmed -> ok | TU | Ledger markers; re-cite ATHP. P 1; 0.5 h |
| 007 | 2: thiamin vs pentose MFT share (half of the hold-out is fit) | duplicate_confirmed -> ok | TU | Apply. P 1; 0.5 h |
| 008 | 0 | duplicate_confirmed -> ok ((b) and (c) not checked against the PDFs) | TU | PDF check. P 0.95; 1 h |
| 009 | 1: MFT / FFT level in a hold-out | open -> moot once scored shape-only | OP (two injections summed) ~ TU (loose wording) | Look for Table 2 compounds with an RI on one column only; otherwise ask the author. P 0.2 for the cause, 1.0 that it is moot; 1 h |
| 010 | 2 | fixed -> ok | TU | none |
| 011 | 3: the TTCA -> CAMA step; phosphate catalysis; the hub's entry for MFT / FFT | open -> ok (fix the explanation text, sec. 0) | CHEM/COND phosphate catalysis 0.45 > OP 0.25 (Zhai 2019 ELSD low end; **post-quench TTCA re-formation**, new) > MOD 0.15 (barrier on its 40 kJ floor) > COND pH 0.1 | Run probe 1 of the CF-011 analysis (Zhang 2025b with the measured-analogue barrier; cheapest, never run) inside cluster 5; guan2012 / thornalley1999 as phosphate analogues. P 0.3; 4 h |
| 012 | 3: MFT / FFT pH law | explained -> ok as a source conflict; the residual misfit is MOD (engine norfuraneol at pH 3 is 55x low) | COND (unbuffered drift plus ARP feed) 0.70 > CHEM (temperature-dependent pH law) 0.15 > OP 0.15 | Register fix only; the residual goes to the engine (a fed ribose + H2S pH-ladder probe). P 0.6; 1 h |
| 013 | 3: MFT loss law | explained -> ok; residual MOD (EX-050: k_mft_decay at least 3x too fast at <= 110 C) | COND (temperature plus precursor supply) 0.75 > CHEM (lipid, Wang) 0.1 > OP (IS added before heating; JAFC run joins) 0.1 | Done as a source conflict. |
| 014 | 2: Jiangnan thiol level | open -> **explained** (Zhou 2023c is the outside test: two prints agree, the non-thiols agree) | OP (thiol calibration or post-quench loss in 2022) 0.65 > COND (headspace) 0.15 > OP (2023 high) 0.15 | Close now. The exact cause needs both Tables S1 (public SI, owner download). P 0.7 for closing / 0.5 for the cause; 0.5 + 2 h |
| 015 | 2: MFT / FFT ratio; the FFT branch in 3:1 concentrated recipes | open -> partly explained (the shape is explained) | COND (sugar excess, a_w) 0.40 > OP (patent analytics) 0.25 > OP (IS before heating) 0.25 > COND (T, t) 0.1 | Bound the IS-timing inflation with EX-050's lambda(110 C) <= 0.02-0.05 /min and van Seeventer's lambda_FFT / lambda_MFT; check against Hofmann's dry / aqueous ratio. P 0.25; 3 h |
| 016 | 2: furfural + H2S -> FFT (k_fur_fft_hs moves 13x) | open -> ok | OP (purge-and-trap under-recovery) 0.35 > TU (a 10x slip) 0.25 > OP (Hofmann high) 0.25 > COND (run) 0.15 | Purge-stripping bound from furfural's K_aw over a 20-100 mL/min flow band; check whether Schieberle 1998's other new aqueous polar values run low against any other print. The outside ratios already side with Schieberle. P 0.25; 3 h |
| 017 | 2: thiol temperature ladder (a gate) | error_confirmed -> ok | OP (fixed-time slice) 0.75 > CHEM (an FFT step from 120 to 140 C) 0.15 > TU (run join) 0.1 | Retire `kang_switch_on_*` and the T4 gate; add the ordering test. P 1; 1 h |
| 018 | 1: trunk / Strecker constants (hazelnut scored 3x; huang methylbutanal) | open -> ok | TU (same-lab reprints) | Compare cells pair by pair (all dossiers on disk). P 0.9; 2.5 h |
| 019 | 2: pyrazine pH law and the S1 mixing bound (both Leahy entries in_model_fit); norfuraneol decay (mikami / nakamura) | open -> partly settled (sec. 0) | TU | Compare the Leahy cells (discrepancies in sec. 5) and the mikami / nakamura cells; declare copy_of; check that B18 / B21 / EX-052 read one copy only. P 0.9; 1.5 h |
| 020 | 2: Strecker (cremer2000 temperatures for the forced-input layer; weenen2001's 5x feeds Strecker signature 7); 1 for the rest | open -> ok (hofmann2000 is ruled) | TU (paper-side slips) | Re-read the PDFs of cremer2000, weenen2001 and xia2022 first. P 0.7; 2 h |
| 021 | 0 | error_confirmed | TU | Apply. 0.2 h |
| 022 | 1: lipid hold-out | error_confirmed | TU | Apply. 0.5 h |
| 023 | 0 | open | TU (unit) | PDF check. P 0.9; 0.3 h |
| 024 | 0 | open | TU (stem / DOI swap) | Check the first pages. P 0.95; 0.3 h |
| 025 | 0 | open -> fixed by rule (the row was refused) | TU | Close. 0.1 h |
| 026 | 1 | open | TU | PDF. P 0.6; 0.5 h |
| 027 | 0 | open | TU (provenance) | Read the code. P 0.95; 0.5 h |
| 028 | 1 | open -> fixed by rule | TU / COND | Close. 0.1 h |
| 029 | 1: lipid slates (B6 and B28 may both score them) | open | TU | Compare cells and check B6 / B28. P 0.85; 2 h |
| 030 | 1: lipid Ea | open | COND (bulk oil vs moist paste; temperature window) 0.5 > OP (CL / TBARS vs hexanal) 0.3 > CHEM (initiation regime) 0.2 | Tabulate Ea by matrix and marker. P 0.4; 2 h |
| 031 | 2: MFT / FFT in lipid matrices (mechanism 6) | open | COND (phase structure: bulk oil vs interface; lipid oxidation state) 0.45 > OP (fat-phase partition under SPME) 0.3 > CHEM (enal trap vs oil protection) 0.25 | Table of method x phase x oxidation state. Within-study contrasts: Yiltirak, SIDA (oil alone lowers MFT 2-10x, the emulsion restores it); Xu 2011 (fresh lard 21x, oxidised 143x); du2025 (+19-30 %). P 0.5; 3.5 h |
| 032 | 1: hexanal binding to pea | open | TU (Snel header 10x, CF-046; loading basis) 0.4 > OP (PRV vs APCI) 0.35 > COND 0.25 | Settle CF-046 first, then convert everything to L/g protein. P 0.55; 2 h |
| 033 | 1 | open | CHEM (sugar-derived) 0.35 ~ OP (misidentification) 0.35 > COND (contamination) 0.3 | Check compound lists against known sugar products. P 0.35; 2 h |
| 034 | 1 | error_confirmed | TU | Apply. 0.5 h |
| 035 | 1 | error_confirmed | TU | Apply. 0.5 h |
| 036 | 0 | error_confirmed | TU | Apply. 0.2 h |
| 037 | 0 | open | TU | Read the methods. P 0.95; 0.3 h |
| 038 | 0 | error_confirmed | TU | Apply. 0.2 h |
| 039 | 1: acrylamide | open | COND (a_w) 0.5 > CHEM (apparent incl. elimination) 0.35 > OP 0.15 | Literature-level explanation. P 0.4; 3 h |
| 040 | 1: acrylamide counted twice | open | TU | Check the unit conversion. P 0.9; 1 h |
| 041 | 1: cultivated-meat IMP (a ribose source for MFT) | open | TU (units / basis) 0.6 > COND (cell line) 0.4 | Convert both tables to one basis. P 0.7; 1 h |
| 042 | 1 | open -> explained (absolute static-headspace scale) | OP | Record it in M6. P 0.6; 1 h |
| 043 | 1 | open | OP / COND | Needs the Hall & Anderson dossier (not on disk). P 0.2 |
| 044 | 1 | open | TU (dimer basis) | Check the code and provenance. P 0.9; 1 h |
| 045 | 1 | open | TU (basis) | Convert per gram. P 0.95; 0.5 h |
| 046 | 1 | open | TU (unit) | Compare text values with the header. P 0.8; 0.5 h |
| 047 | 0 | open | COND | Carry it as a condition. 0.2 h |
| 048 | 0 | fixed | -- | -- |
| 049 | 0 | error_confirmed | TU | Notes. 0.5 h |
| 050 | 0 | not_an_error | -- | -- |
| 051 | 1 | error_confirmed | TU (phantom constants) | Apply. 1 h |
| 052 | 2 | error_confirmed; a DUPLICATE of CF-072 | MOD | Merge (sec. 0). |
| 053 | 2: amine catalysis (Cao's 24x glycine arm is cited for CF-011 E3 and for B7b) | open | TU (meaning of the multiplier) | Back-solve [Gly] from Gly-ARP formation vs xylose loss across the arms (pseudo-first order in Gly); check sibling Xie papers; ask the authors. P 0.35; 2 h |
| 054 | 1: pH drift at 180 C (a control for the pH model) | open | COND (endpoint hypersensitive after an alkaline start) 0.6 > TU / OP (reading) 0.4 | Charge balance: does a +-20 % acid yield span 4.96-6.05? P 0.5; 2 h |
| 055 | 2: B16 MFT ratio row | open -> error_confirmed | TU (column slip in the text) | Update the status. 0.1 h |
| 056 | 0 | open -> resolvable now | TU | Zhai 2023 Food Chem Table 1 prints MFT 1.237 +- 0.142 (CF-017 analysis), so S5's +-0.061 is the stray value. Close. 0.1 h |
| 057 | 1 | open | TU (dropped decimal) ~ a real value | Check against row neighbours and any printed column total. P 0.5; 0.5 h |
| 058 | 2: fed sulfur-bearing ARP (the MFT supply gate) | open -> identity resolved | MOD (missing species) | RG-ARP is the ribose-glutathione Amadori product (feng2022 dossier title), not CAMA. Keep the rows out until a GSH-ARP species exists. P 0.9; 0.3 h |
| 059 | 0 | fixed | -- | -- |
| 060 | 3: hub dicarbonyl scale (3-DX, 1-DX, GO, MGO, which feed Strecker and pyrazines) | open -> ok | OP (derivatisation / dilution protocol) 0.4 > COND (ARP batch, co-purified catalyst) 0.25 > TU (unit) 0.2 > other 0.1 > CHEM 0.05 | (i) Carbon closure: dicarbonyl C-mol vs the pentose C of ARP lost at each time (the ARP decay is printed and agrees); a print that exceeds it is wrong. (ii) Do the methods alternate? 2022 and 2024 are high, 2023c and 2025b low. (iii) zhou2024b dossier (methods ref. 20). P 0.45; 3 h |
| 061 | 1-2: GO as a C2 donor | open -> explained (no GO standard in 2019; the 2021 rows are below LOD; the remedy is already in the tables) | OP | Close. 0.3 h |
| 062 | 2 | fixed | OP | Trunk-twin flag. 0.3 h |
| 063 | 1 (2 for the missing mass) | open | OP/TU (fit window, unreported intercept) 0.5 > CHEM (unmeasured sink) 0.4 | Refit the digitised points with a free intercept and lag; pass the missing half of the balance to CF-076 as the sink term. P 0.75; 1 h |
| 064 | 1 | open -> explained by the remedy | TU | Close. 0.1 h |
| 065 | 3: the furfural / NF + sulfide steps of the MFT / FFT routes | open -> ok (fix the source first, sec. 0) | COND (speciation: HS- is ~1 % at pH 5 vs ~50 % at pH 7, about 50x; H2S lost to the autoclave headspace at pH 5) 0.45 > CHEM (NH3 from (NH4)2S) 0.25 > OP (Zhai semi-quantitative) 0.2 > TU 0.1 | Speciation-corrected yields with pKa1(T), Henry's-law headspace loss and a 120 vs 145 C barrier band. Check: Hofmann Table 8 (C2 + C3 + H2S at pH 3 / 5 / 7). P 0.45; 3 h |
| 066 | 1 | duplicate_confirmed | TU | Done. |
| 067 | 2: hub C2 / C3 (Hofmann 1999 is the hub's calibrated MGO / glycolaldehyde set) | open | OP (OPD-limited: 2 umol OPD against ~25 umol of dicarbonyls; an IS ratio can hide it) 0.45 > TU (misprint, 1 mol/L) 0.35 > TU (figure vs table run) 0.2 | If OPD limited the result, GO / MGO / DA (relative to the IS) would dip when osones peak (10 min). Test that shape. P 0.5; 2 h |
| 068 | 3: absolute MFT; thiol loss; O2 as a reagent | open | COND (headspace O2 : Cys 2.0 vs 0.26) 0.45 > OP (thiolate in 1 M NaOH in air; storage at -20 C) 0.3 > OP (Garching over-read) 0.1 > TU (aliquot basis) 0.05 | O2 : Cys for every SIDA thiol pot (Hofmann 0.26; Bolton ~2.0, engine ~20x over; Yiltirak 2.0; both patents; Cerny) against the engine residual. Yiltirak's matrix-matched, same-work-up calibration weighs against the work-up branch. P 0.4; 4 h |
| 069 | 2: acid budget -> pH in unbuffered pots | open | MOD (r_aa_frag absorbing another defect) 0.5 > OP (Martins acetic acid late) 0.25 > CHEM 0.15 > COND 0.1 | Arm E on the keyed ph_acid tables. P 0.5; 4 h |
| 070 | 2: Phe consumption in the forced-input Strecker layer | open | TU (one series printed twice) 0.65 > CHEM (a common step) 0.3 | Phe's own thermal loss at 100 C is negligible, so a common step needs equal donor-driven loss. Test Phe degraded against each pot's products. Deciding WHICH pot owns the series is the hard part. P 0.55; 1.5 h |
| 071 | 3: the hub's Cys trap of 3-DX (signature iii) | open | COND (time window, Cys depletion, start pH 7.4 vs 7) 0.4 > CHEM (the trap fades) 0.25 > TU (copied JAFC cells, CF-066) 0.25 > OP 0.1 (a within-paper double ratio is response-factor-free) | Recompute both double ratios at matched times (<= 60 min, 120 C), drop the copied cell, test for a time trend in JAFC. P 0.5; 2 h |
| 072 | 2 (the triage's #2) | open (absorb CF-052) | MOD | Refit B20 under Nguyen's topology against the S1B mass balance (SI on disk), in its own PR. P 0.7; 6 h |
| 073 | 3: Strecker aldehyde share (the forced-input layer's concentration-independent s) | open | OP (Likens-Nickerson stripping drains a reversible PA-amine pool and stops losses) 0.4 > COND (26 vs 100 mM; citrate pH 3.2) 0.35 > CHEM 0.15 > TU 0.1 | Hofmann's PA plateau bounds the loss branch; granvogl2006 (DFA, 5 mM, 7.3 % PA, but at 170 C); pripisnicolau2000 (third lab, pH 3.5). The imine-pool rival predicts Weenen's yield rises with distillation time. P 0.4; 3 h |
| 074 | 2: Strecker O2 acid route (B29's air / argon ratio 3.5) | open | TU (one run printed twice) 0.75 > coincidence 0.25 | 2000b Table 2's ARP-Phe air / argon pair is not a reprint. If its ratio is close to 3.5, the vessel confound is small. P 0.45; 1 h |
| 075 | 1 | open -> explained by the remedy | TU (text) | Close. 0.1 h |
| 076 | 3: cysteine fate (H2S supply for MFT / FFT; the Cys trap) | open | CHEM (catalysed exchange) 0.35 > OP (**post-quench re-formation**: Zhai 2019 makes 14.5 mM TTCA at room temperature before heating and stores at 4 C; ELSD low end) 0.25 > MOD (shared skeleton fails T2 / T3) 0.15 > COND (in-situ pH) 0.15 > CHEM (second sink) 0.1 | Quench / storage feasibility: k_f(RT) from Zhai 2019's t = 0; K(T) from the kallen1971b / pesek1975 / fife1991 dossiers. It predicts concentrated-pot TTCA nearly flat in T at early times and Cao's flat K_app. Zhai 2021 SI (asked). P 0.3; 5 h (the quench sub-test alone: P 0.25, 2.5 h) |

---

## 2. Top 12 by importance x P / hours

This list is restricted to importance >= 2 and to conflicts that still need an investigation. Register-only closures are
listed separately below.

| rank | CF | I | P | h | score | one-line test |
|---|---|---|---|---|---|---|
| 1 | 019 | 2 | 0.9 | 1.5 | 1.20 | Leahy and mikami / nakamura cell compare; one copy per fit |
| 2 | 074 | 2 | 0.45 | 1 | 0.90 | ARP-Phe air / argon ratio vs Glc / Phe 3.5 |
| 3 | 071 | 3 | 0.5 | 2 | 0.75 | Cys double ratio at matched times; copied cells dropped |
| 4 | 070 | 2 | 0.55 | 1.5 | 0.73 | Phe loss vs each pot's products |
| 5 | 020 | 2 | 0.7 | 2 | 0.70 | PDF re-read: cremer2000, weenen2001, xia2022 |
| 6 | 067 | 2 | 0.5 | 2 | 0.50 | OPD-limitation signature (dip at the osone peak) |
| 7 | 060 | 3 | 0.45 | 3 | 0.45 | Carbon closure of the Zhou dicarbonyls; methods alternation |
| 8 | 065 | 3 | 0.45 | 3 | 0.45 | HS- speciation plus headspace correction (after the source fix) |
| 9 | 073 | 3 | 0.4 | 3 | 0.40 | Stripping / imine-pool vs concentration; third lab |
| 10 | 053 | 2 | 0.35 | 2 | 0.35 | Back-solve Cao's glycine from Gly-ARP formation |
| 11 | 068 | 3 | 0.4 | 4 | 0.30 | O2 : Cys regression over SIDA thiol pots |
| 12 | 076 (+011) | 3 | 0.25 | 2.5 | 0.30 | Post-quench TTCA re-formation feasibility |

Next in line: CF-031 (0.29), CF-069 (0.25), CF-072 (0.23), CF-011 (0.23), CF-016 / CF-015 (0.17).

**Register-only closures (about 2 h in total; do them before the clusters):**
- CF-055 -> error_confirmed.
- CF-056 closed (S5's SD is the stray).
- CF-058: identity settled.
- CF-061, CF-064, CF-075 -> explained by their remedies.
- CF-014 -> explained (Zhou 2023c is the outside test).
- CF-009 -> moot under shape-only scoring.
- CF-025 and CF-028 -> fixed by rule.
- The edits listed in sec. 0.

Confirmed errors whose remedies are still unapplied, batched for the next regeneration: CF-001, 003, 004, 005, 007, 017,
021, 022, 034, 035, 051. All their strings are still in `src/` / `scripts/` (grep, 2026-11-01).

---

## 3. Clusters (one agent each)

### Cluster 1 -- Same-lab reprints and text / table slips (bookkeeping that protects single-lab evidence)

- **Conflicts:** CF-019, 070, 074, 020, 067 (figure / table half), 018, 029, 040, 057.
- **Shared hypothesis:** the Ho network, DFA and Wageningen re-print series across papers; two prints of "independent"
  runs are one measurement, and single-lab mechanism evidence (the Strecker O2 route, the pyrazine S1 bound) is thinner
  than counted.
- **Files and data:**
  - dossier pairs (leahy1989 / leahy1989a, mikami2017 / nakamura2020, hofmann2000 / hofmann2000b, cremer2000,
    weenen2001, xia2022, the trunk pairs);
  - `data/lit/fit_rows/strecker_hofmann2000*.yml`, `deoxyosone_hofmann1999.yml`;
  - `data/lit/measurement_ledger/*.yml` (flags `duplicate_of`, `same_lab_republication_of`).
- **Discriminating test, a fingerprint:**
  - Two series are one measurement if at least 3 cells agree to the printed digit (mean AND SD where printed), whereas
    independent runs of these recipes scatter 7-16 % (CF-075).
  - For CF-070, a common chemical step must also close the Phe mass balance in both pots.
  - For CF-074, the ARP-Phe air / argon ratio is the control.
- **"Resolved" means:** duplicate_confirmed with a copy_of map and every fit reading one copy, or not_a_duplicate with
  the cell differences printed. For CF-074, the B29 ratio is kept if the ARP-Phe ratio lies within 0.15 dex of 3.5;
  otherwise it is flagged as a vessel comparison.
- **Systematic add-on (~3 h):** a script that hashes (value, SD) pairs across all keyed tables and ledgers by lab. It
  would catch the reprints nobody has looked for yet (the Ho network has printed one ladder four times).
- **Wait for dossiers:** none; everything is on disk.
- **Time:** about 10 h.

### Cluster 2 -- Hub dicarbonyl scale and the cysteine trap in the Jiangnan / Xie pots

- **Conflicts:** CF-060, 071, 053, 063 (and close CF-061).
- **Shared hypothesis:** single-lab quantification and handling (quinoxaline protocol, dilution, copied cells, fit
  windows), not chemistry, make the 10x scale differences. Scale-free within-paper ratios survive.
- **Files and data:**
  - `data/lit/fit_rows/deoxyosone_submodel.yml`, `route_submodel.yml`;
  - the zhou2022 / 2023c / 2024 / 2025b, zhai2019 / 2021 / 2023jafc and cao2017 dossiers;
  - zhang2025b Figs. 2a-h.
- **Discriminating tests:**
  - carbon closure against the measured ARP loss (CF-060);
  - matched-time double ratios (CF-071);
  - [Gly] back-solved from Gly-ARP kinetics (CF-053);
  - an intercept refit (CF-063).
- **"Resolved" means:** a print that violates carbon closure by more than its read error is declared wrong and scored
  shape-only. CF-071 is explained if the matched-window ratios agree within 0.15 dex. CF-053's arms re-enter only if the
  back-solved [Gly] is below the solubility limit.
- **Wait for:**
  - zhou2024b (the methods reference of Zhou 2024), cui2022, hou2026, weenen1992;
  - Zhou 2024 SI (asked);
  - for the phosphate side, thornalley1999 and guan2012.
- **Time:** about 9 h.

### Cluster 3 -- In-situ pH and speciation (unbuffered pots, sulfide and amine speciation)

- **Conflicts:** CF-065, 054, 069, the residual of 012, and the pH arm of 011 / 076 (a service to cluster 5).
- **Shared hypothesis:** label pH is not reaction pH. Speciation (HS-, free amine, thiazolidine K(pH)) at the in-situ
  pH and temperature explains apparent rate or yield gaps without a new constant.
- **Files and data:**
  - measured trajectories in `data/lit/fit_rows/ph_acid_submodel.yml`: cao2017_fig4c, zhou2023_fig2 (3.22 / 3.42 /
    5.07), martins2005b drift ladder, deng2021;
  - endpoints: Kang 5.5 / 7 / 8 -> 5.1 / 4.9 / 4.7; Zhai 2023 Food Chem 7.0 -> 4.9; chen2000 / chen2002;
  - the engine's charge balance;
  - Hofmann 1998 Tables 3 / 4 / 8 (`hofmann1998_reconciliation.md`).
- **Discriminating test:**
  - Positive control: the pH model reproduces the measured trajectories within 0.3 pH units.
  - Negative control: 0.5 M buffered pots do not move.
  - Then CF-065's yields are corrected for HS- fraction, headspace H2S and temperature, and CF-054's endpoint
    sensitivity is computed.
- **"Resolved" means:** the corrected gap is <= 0.5 dex, with no free constant added. Otherwise the residual is declared
  chemistry (for CF-065, the NH3 co-reactant).
- **Wait for:**
  - kallen1971b and pesek1975 (thiazolidine K(pH));
  - sedehi2013 and noziere2009 (dicarbonyl + ammonium vs pH; relevant to the (NH4)2S in CF-065).
- **Time:** about 9 h.

### Cluster 4 -- Observation operators and vessel O2 in the thiol and Strecker levels

- **Conflicts:** CF-068, 073, 031, 016, 015, plus the cause half of 014.
- **Shared hypothesis:** the level gaps track the method (SIDA after cooling, IS before heating, HS-SPME with f' = 1,
  Likens-Nickerson stripping, purge-and-trap, thiolate work-up) and the vessel O2 : limiting-precursor ratio, not the
  chemistry.
- **Files and data:**
  - one table of every thiol / Strecker level pot with method, IS timing, work-up, fill / headspace, O2 : Cys and
    yield per mol of limiting precursor: `thiol_submodel.yml`, `strecker_hofmann2000*.yml`, the yiltirak2026 dossier
    sec. 5 (the O2 block), bolton1993/1994, the us10076126 / us8920862, weenen2001 and schieberle1996 / 1998 dossiers;
  - evidence map sec. 2.11.
- **Discriminating test:**
  - Regress the log yield residual (against the engine or a within-lab reference) on the operator covariates.
  - Within-study sign checks: Yiltirak oil vs emulsion, the US 8,920,862 lipid arm, Schieberle 1996's method
    comparison, Hofmann air / argon.
- **"Resolved" means:** one covariate explains >= 70 % of the between-paper log variance and every within-study
  contrast agrees in sign. Otherwise the spread becomes the declared between-lab band (thiol prereg sec. 6).
- **Wait for:** grosch1992 (SDE vs direct MFT), mottram1996b, adams2001, nikolantonaki2012, esterbauer1975, doorn2002,
  kim1996b, pripisnicolau2000 (saved as pripis-nicolau2000.pdf), hidalgo2005, zamora2006, ghiron1988;
  charlesbernard2005a is not on disk.
- **Time:** about 14 h.

### Cluster 5 -- Cysteine exchange: catalysis vs sample handling

- **Conflicts:** CF-076, CF-011 (CF-063's sink as an input).
- **Shared hypothesis, two rivals:**
  - (i) catalysed exchange, k0 + k_aa[AA] + k_P[P];
  - (ii) an operator: TTCA re-forms during cooling and storage in the concentrated pots (Zhai 2019 forms 14.5 mM TTCA
    at room temperature in about 10 min; 4 C storage), which makes them look like fast exchange; plus ELSD at the low
    end.
- **Files and data:**
  - `docs/history/cys_exchange_2026-10-31/` (EX-054 battery T1-T7), `thiol_rivals_2026-10-31/cys_fate.py`;
  - `route_submodel.yml`;
  - the zhai2019 / zhai2021 / cao2017 / zhang2025b dossiers;
  - fife1991, nagasawa1981, roberts1987, schubert1939 (dossiered).
- **Discriminating test:** re-score EX-054's battery with a quench / storage operator (cooling ramp plus hold at 4 C or
  room temperature, using k_f(RT) bounded by Zhai 2019's t = 0 and K(T) from the physical-chemistry dossiers).
  - Rival (ii) predicts early TTCA in the concentrated pots that is nearly temperature-independent, and Cao's flat K_app.
  - Rival (i) predicts neither.
  - Then run the cheap CF-011 probe 1 on Zhang 2025b.
- **"Resolved" means:** a structure passes T4-T7 within the pre-registered bands with no dropped test. If only the
  operator does, the conflict is explained as OP and k_ttca_cys is fitted on the dilute pots only. If neither passes,
  the conflict goes to designed experiment E4.
- **Wait for:** kallen1971b, pesek1975, schubert1936, hou2026 (cystine), guan2012; the Zhai 2021 SI (13C5 label
  fraction; asked).
- **Time:** about 8 h.

Kept out of the clusters as engine work: CF-072 (its own PR, as the owner planned) and the sink residuals behind
CF-012 / 013 / 015 / 017 (EX-050).

---

## 4. Patterns

**Primary cause class, all 76:**

| class | count | entries |
|---|---|---|
| TU | 44 | ~18 repo-side (our transcription, code or dossier); 11 paper-side duplicates (005, 006, 007, 008, 018, 019, 029, 040, 066, 070, 074); ~15 paper-side slips |
| OP | 14 | 004, 009, 014, 016, 017, 032, 042, 043, 060, 061, 062, 063, 067, 073 |
| COND | 13 | 011, 012, 013, 015, 030, 031, 039, 047, 054, 065, 068, 071, 076 |
| MOD | 4 | 052 = 072, 058, 069 |
| CHEM | ~1 | 033 at most |

Counting secondary candidates: OP appears in about 24 entries, COND in about 20, CHEM in about 12, and MOD in about 8
(residual engine misfits).

**Of the ten importance-3 conflicts** (004, 011, 012, 013, 060, 065, 068, 071, 073, 076), none has genuine chemistry as
its clear first cause. What chemistry remains is on the engine side (MFT sink, norfuraneol at pH 3, exchange law). The
literature disagrees mostly about how and under what conditions things were measured, not about the chemistry.

**Four systematic issues, each resolving several conflicts at once:**
1. **Scale-free scoring of the Jiangnan / Xie HS-SPME and OPD-internal-standard data**, applied uniformly (CF-062 did the
   deoxyosones). This dissolves the level halves of CF-004, 009, 014, 060, 061 and 065, and makes CF-071 a
   response-factor-free test. Mostly a measurement-model setting, plus a few hours of ledger work.
2. **A same-lab fingerprint scan** (value + SD hashes across ledgers, by lab). It settles the duplicate class (11
   entries plus CF-070 / 074) and finds unregistered reprints before they double-count single-lab mechanism evidence.
   About 3 h.
3. **A vessel and sample-handling covariate block per pot:** headspace : liquid, O2 : limiting precursor, IS timing,
   quench and storage, derivatisation timing. "Unstated pot variable" appears as a candidate in CF-014, 015, 016, 067,
   068 and 076. Many papers do state the vessel, so recording it turns the hand-waving into a regression covariate.
4. **In-situ pH for every unbuffered pot**, from the existing pH-acid submodel calibrated on the measured trajectories.
   It serves CF-011, 012, 054, 065, 071 and 076.

Rough probabilities:
- P(1 + 2 + 4 together move at least 6 open entries to explained / duplicate_confirmed by 11-07): ~0.6.
- P(cluster 4 explains the 100x thiol level spread, CF-068, without designed data): ~0.35.
