# Red team of resolve2_sulfur (2026-11-01)

Scope: `resolve2_sulfur/{prereg,results,corrections}.md`, `register_proposals.yml`, `phi1_*.{py,json}`. Read-only on the
repo. Own computations: `e2c_points.json` (four point evaluations, E2c coupled, no search), `e2c_ea_de.json` (DE, E2c
coupled, Ea fixed 55 / 50, same settings as resolve2), `fig4.png` (van den Ouweland p. 505 Fig. 4, 400 dpi crop).
Primary sources re-read: karpinski1990.pdf (buffer text), vandenouweland1975.pdf (Table I pp. 502-503, Fig. 4, text),
hofmann1998.pdf (intro), whitfield1999.pdf (intro, Table 1, discussion), zhang2022 dossier SI tables, elmore2002
dossier, split_v2_per_change.yml, F28 atlas, heldout_zhang2021_table3.yml header.

---------------------------------------------------------------------------------------------------------------------

## 1. CF-076 / E10: "EXPLAINED (provisional), conditional on release Ea <= ~45": **WEAKENED**

What holds: the arithmetic of the coupling (k_r45 = k25_6 x 10^(-1.1 s_lo) x Arr) and the reported numbers reproduce
(face E2c: k_r45 1.96e-4 => Ea 26.5; checked). The Ea 64 rejection (11.95) is robust. phi1_face.py's monkey-patch does
reach the spawn workers (the main module is re-imported as `__mp_main__`), and res.fun matches the recomputed z.

What does not hold:

1. **The k25 floor is not self-consistent with either structure it is used for.** The floor 3e-5 is Karpinski's AThz
   1.8e-4 deflated x6 = x3.3 "phosphate / ionic-strength excess" x 1.8 (phi at its lower 30 %-RSD edge).
   - E2c has *no* buffer-catalysis term; a structure that says catalysis is absent cannot deflate its own anchor for
     phosphate catalysis. Self-consistent E2c band: ~[1.2e-4, 2.5e-4] (the RSD only) -> required Ea ~26-29.
   - E4c carries its own phosphate term. Evaluated at Karpinski's conditions (pH 6.0, <= 0.1 M phosphate, so
     HPO4 2- <= 7.4 mM, gate 0.91) fcat is 1.5-3.3 over the cP / rho bands and **1.56 at the E4c-coupled best point
     (cP 101, rho 3.24)**, not 6. Self-consistent k25 ~1.2e-4 -> Ea ~28-30.
   - The x3.3 gap is between a *measured* rate (Karpinski, low ionic strength, buffer unstated) and a *derived* one
     (Kallen k_f / K_app at I = 1.0 M). Assigning the whole gap to Karpinski is a choice in the admission's favour.
   - Fair counterpoint (found here): Bustillo 2018 METH-CYS (2-alkyl, 20 C) at pH 4.4 runs 1.3e-3 /min = Karpinski
     PrThz at pH 4.5 (1.3e-3), but at pH 6.0 Bustillo is x5.6 slower (7.8e-4 vs 4.4e-3). So a second lab does say
     Karpinski's pH 6 runs are fast. But that argues for anchoring at pH 4.5, where the rate is uninflated:
     AThz(25 C, pH 4.5) ~ PrThz(4.5) / 24 = 5.4e-5 /min; the feasible k_r45(100 C) 0.9-3.2e-4 then needs only x2-6
     from 25 to 100 C, i.e. **Ea ~7-22 kJ/mol**. Every anchoring route that avoids the stacked deflation gives Ea < 32.

2. **The fitted corner is below every pH-4.5 analogue.** k_r(25 C, pH 4.5) implied by the feasible points:
   E2c coupled 4.2e-6, E4c coupled 1.2e-5 /min. Analogues: Kallen TC 1.4e-5 (x phi 0.55-1.9: 7.7e-6-2.7e-5); AThz via
   PrThz shape 5.4e-5. E2c sits x1.8 below the lowest analogue edge and x13 below the direct one. With E2c's k25 at its
   floor AND s_lo at 0.80-0.85 of a 0.85 ceiling AND Ea at 40-43 of a 40 floor (three active bounds, all in the "slower
   release" direction), the pattern is the data asking for a release slower than any anchor allows. That is a misfit
   signature, not an explanation. P(the slow-release reading is a structural stand-in for something else: an
   unmodelled sink / TTCA quantification / K) ~0.5.

3. **"Rejected at Ea >= 50" is E4c-only and wrong for E2c** (own check). E2c coupled, k25 at the floor, no search
   (nuisance parameters from the face point): Ea 48.5 / s_lo 0.81 -> max|z| **1.83 (feasible)**; Ea 50 / s_lo 0.85 ->
   **2.05**; Ea 55 -> 3.11; Ea 60 -> 4.82. DE at fixed Ea 55 then 50 (`e2c_ea_de.py`) was still running after 65 min at load ~100 when this was
   written; its result lands in `e2c_ea_de.json` / `.log`. The point evaluations alone settle the claim: an
   unoptimised point is an upper bound on the best max|z|, so E2c is feasible at Ea 48.5. So the
   threshold is ~50 at the floor for E2c, not 45; and since the floor is the questionable part (item 1), the honest
   condition is a ridge, not a single Ea.

4. **pKa_T "no longer binds" is overstated, and the comparison mixes temperatures.** The uncoupled positive control
   sits at pKa_T 5.21 (the new floor 5.2); the coupled points at 5.29-5.55, all in the lower third. More important:
   the battery's pKa_T is temperature-independent and the pots run at 70-130 C, so the fitted value is an effective
   at-temperature pKa; Radomski's 5.2-6.1 is a room-temperature charge-balance inference (no measurement; and any
   Cys release in the 20 mM solution, pI 5.1, biases it upward). The engine itself gives the TTCA ring N a
   Delta H_ion of 40 kJ/mol (goldberg2002 dossier: a class value not in GKL), i.e. ~1.4 pK units lower at 100 C.
   Two consequences: "inside Radomski's band" is not a like-for-like check; and on the acid limb a pKa that falls with
   T makes k_r(100 C, pH 4.5) rise faster than Arrhenius, so the *intrinsic* barrier needed is lower still than the
   fitted effective 26-49 (by up to ~20-30 kJ/mol in phosphate pots, less in amine-buffered ones). Confidence that this
   matters by > 10 kJ/mol: ~0.6.

5. **Is 26-32 kJ/mol plausible for thiazolidine release?** Nothing on disk measures it. Nearest numbers: Fife 1991
   IV OH- term ~58 (composite, folds in K); Zhai 2019's own forced-input Arrhenius for the TTCA -> CAMA flux 69 (same
   system, passes through ring opening; cluster-5 verdict sec. 4 already used it against an Ea-floor reading); Kang
   2026 free-Cys loss 55; Cys-ARP formation 100.9. Ring-chain opening is endothermic, which adds to, not subtracts
   from, an apparent release Ea. P(effective release Ea at fixed pH <= 32) ~0.1; <= 45 ~0.25; <= 50 ~0.35.

6. **Controls.** NC1 / NC2 are parameter-free constants scored on 20 random draws, while rivals get DE: they show the
   battery is not vacuous, nothing more. The flexible rivals rejected in cluster 5 (E1c 4.0, E5 5.2, E6 5.1) are the
   meaningful negative controls; cite them. The "positive control" (uncoupled E4c, new pKa band, 1.29) does not test
   the coupling, and "1.29 < cluster-5's 1.45, so maxiter 20 is adequate" confounds the pKa-band change with search
   depth (the band floor moved 5.6 -> 5.2 and the optimum sits on it). A real search check is cluster-5's original
   bands at maxiter 20 (should reach ~1.45). One seed only.

**Corrected probability:** explanation 1 of CF-076 0.35 -> **0.25**; "something unlisted / second sink" up to 0.3.
Decider unchanged and good (one Arrhenius series of release at fixed pH); add that it must be run at pH 4.5 (the
acid limb) and with the pH of the buffer at temperature recorded.

---------------------------------------------------------------------------------------------------------------------

## 2. van den Ouweland 1975: "not a misnaming; P(actually MFT) ~0.15"; CF-092 erratum: **WEAKENED (erratum REFUTED as worded)**

Upheld: Table I peak 3 is named 4-mercapto-2-methylfuran and drawn with SH on C-4 (remote from CH3); footnote a reads
"Structures of these compounds confirmed by comparison of their spectra with synthesized reference compounds"; the
table lists both thiophene regio-isomers (peak 10, 3-SH; peak 13, 4-SH). P(the print intends MFT) <= 0.05: fair.
Hofmann 1998 intro: "van den Ouweland and Peer (1975) were unable to identify MFT": verified.

Missed or misquoted:
1. **Fig. 4 of the same paper (p. 505) draws MFT.** The postulated pathway for X = O, R = H: NF -> 4-mercapto-5-methyl-
   3(2H)-furanone (H2S replaces the enol OH, the carbon next to the methyl) -> reduction / dehydration -> the aromatic
   thiol with SH adjacent to CH3, i.e. 2-methyl-3-furanthiol, plus its dihydro / tetrahydro forms (`fig4.png`). The
   paper's own mechanism predicts the 3-SH regiochemistry; peak 3's 4-SH assignment contradicts it. This is also the
   pathway Hofmann 1998 Fig. 1 prints. So the paper is a legitimate citation for the *proposed* NF -> MFT pathway; it is
   not a citation for MFT *identified*.
2. **Whitfield & Mottram 1999 attribute MFT to van den Ouweland twice, not once.** Discussion (pp. 1629-1630): compounds 6
   (MFT) and 17 "were also reported by van den Ouweland and Peer (1975) to be major products of the reaction between
   HMF and hydrogen sulfide" (paraphrase), in addition to the intro. corrections.md calls it "one intro sentence".
3. **Whitfield's Table 1 column is a spectral-identification reference, not an attribution of who reported MFT** (it
   sits beside "MS + LRI" and LRI 867; row 7 carries MS ions in the same column). "Whitfield's own Table 1 attributes
   MFT to Evers 1976" is a misreading and cannot be used against the intro.
4. **Four later labs make MFT from NF; none reports 2-methyl-4-furanthiol** (Whitfield 1999 lists 69 volatiles; the
   4-isomer is absent). The MS argument ([M-H]+/M+ 0.50 vs 0.28) compares a 1975 trapped-fraction spectrum with NIST:
   weak, instrument-dependent (the dossier itself says so). Peak 3's IR lists no S-H band (~2550), unlike peaks 7, 8,
   10: a small flag on the identification as a whole.

**Corrected probabilities:** peak 3 is in fact MFT (misassigned regiochemistry, by name + drawing or by the reference
synthesis): **~0.3** (keep the dossier's value; 0.15 is not supported once Fig. 4 is read). CF-092 explanation 1
("mis-citation, citation trail only", P 0.85) -> split: "legacy texts overstate a proposed pathway (Fig. 4) as an
identified / synthesised MFT" 0.65; "the paper did make MFT and misassigned it" 0.3; other 0.05. NF-route existence
without the paper: P 0.97 upheld (four labs verified in their dossiers by the writer; not re-read here).

Six files, current text vs erratum:
- isotope_topology_evidence.md row 12 ("reacts with H2S (van den Ouweland & Peer's original synthesis)") and l. 100
  ("cited to van den Ouweland & Peer 1975"): the correction is warranted, but "report ... NOT MFT" must become "propose
  the NF -> MFT pathway (Fig. 4) but list no MFT in Table I (peak 3 named 4-mercapto-2-methylfuran)".
- AUDIT_legacy_lane l. 358-360 ("describes a synthesis of MFT from norfuraneol"), README_legacy_lane l. 391-392
  ("describing a genuine synthesis of MFT"), audit_remediation l. 7488 ("report the same synthesis"): erratum lines are
  warranted; the proposed line "do not report MFT" is accurate for Table I but must add "their Fig. 4 postulates it".
  audit_remediation l. 480 / 1633 only cite the route: the erratum adds nothing there beyond the same note.
- deep_research_backlog.json: the description is wrong for this citation (verified), but **the remedy "fix the
  source 03_thiamine_degradation.md l. 590 and regenerate" is not executable**: no such file exists in
  data/research_corpus (main checkout or worktree) or in git history under that name. Regenerating may silently drop
  the entry. Say so in the remedy.
- whitfield1999_extraction.md flag: as worded it is wrong (items 2-3).
- vandenouweland1975_extraction.md: do not lower 0.3 -> 0.15; add Fig. 4 and Hofmann's reading.

---------------------------------------------------------------------------------------------------------------------

## 3. Wei 2022: structural replacement of "2 H short": **UPHELD with one chemical error**

Checked: TTCA = TC (C4H7NO2S) with H -> C4H9O4 = C8H15NO6S, M 253.27; Cys + Xyl - H2O = same. 1-DX
CH3-CO-CO-CHOH-CH2OH and 3-DX OHC-CO-CH2-CHOH-CH2OH are both C5H8O4 (132.11); DX + Cys - H2O = C8H13NO5S (235.26);
DX + Cys = C8H15NO6S. The Cys-ARP is isobaric with TTCA. The 1-DX ring would be 2,2-disubstituted (no C-1 aldehyde);
3-DX gives a 2-(2-oxo...) ring, M 235. Correct.

Error: corrections.md says 1-DX -> TTCA "needs the C-1 methyl oxidised back: not a Maillard step". 1-DX + H2O =
C5H10O5: the path is redox-neutral overall (C-1 up, C-2/C-3 down), exactly as the dossier's own corrected text says
(l. 173-176: "redox-balanced overall ... internal C-2 / C-3 redox swap"). It is the reverse of the 2,3-enolisation /
beta-elimination: conjugate addition of water at C-1 of the 1-en-2-ol-3-one tautomer gives a pentulose, which must then
isomerise to the aldose. That is the same "return to the pentose level" as for 3-DX. Writing "oxidised" re-introduces
the withdrawn redox argument in another form. Also: the regenerated sugar can epimerise (ribose / arabinose / lyxose
rings, all m/z 254), which the identification cannot separate.

Split 0.6 / 0.15 / 0.25: defensible. Since any DX route runs through the pentose level, "via DX as a direct trap" is
<= 0.1; the categories overlap (a DX that reverts to pentose is counted twice). Optional: 0.6 / 0.1 / 0.3.

---------------------------------------------------------------------------------------------------------------------

## 4. CF-031 / E7: "free-disulfide shift, P 0.65"; "minimal E7 = one capture term, no partition term": **WEAKENED**

Numbers reproduce (water-phase pool: control 2.97 MFT-eq; raw-oil 1/3.4-1/25, median 1/13.5; ENTH 1/1.0-1/3.5;
Elmore pool ~1/17 with SDs of 63 % and 105 % on 807 and 202). But:

1. **Zhang 2022 measured the oil phase, and it holds MFT and the disulfide.** SI: MFT disulfide in oil 288.49
   (coconut RM), 59.74 (coconut EM), 21.74 (palm RM), 20.64 (peanut EM) ng/mL; free MFT in oil 13.79 (olive RM),
   17.76 (coconut RM), 10.18 (peanut EM). Oil : aqueous is 4 : 1, and the dossier records that the oil holds most
   volatiles (flag 3: oil and water headspace not on a common basis). The "pool falls x1/13" is a **water-phase** pool;
   the system pool is unknowable from these data. The disulfide's disappearance from raw-oil water phases is at least
   as well explained by partition into the oil (lipophilic disulfide) as by capture. (Against pure partition, to be
   fair: ENTH water phases keep 0.6-2.3 ng/mL disulfide under the same 4 : 1 oil.)
2. **The pre-registered rule is not met.** Prereg sec. 4: every lab that measured the disulfide shows pool suppression
   within x3. Zhang's ENTH arms (5 of 10) give x1/1.0-1/3.5 vs Elmore's x1/17 (5-17x apart). The pass was obtained by
   restricting to raw oils after seeing the data. By the rule: (a) not established; at best inconclusive.
3. **TEST data were used as structure evidence.** heldout_zhang2021_table3.yml declares itself "TEST data ... engine
   error E7"; corrections.md / CF-031 use its 6 / 6 contrasts as the fourth supporting lab and as the argument against
   partition ("no bulk oil phase"). And **Xu 2011's lard contrast is E12's TEST contrast F28#20@xu2011** (0-based
   index 20 of F28 = the Huazhong lard cell; verify the indexing convention): the proposed E7 capture term is built
   mainly on it, and would move the E12 score on it. The 6 / 6 check is not a development cell reused; it is worse:
   held-out cells viewed before the structure is chosen.
4. **"No partition term (not supported cross-lab)"** is too strong: partition into a bulk oil phase is physics; the
   only question is whether it is needed for the observables on disk. Rephrase as "not needed to explain the sign data
   on disk"; cluster 4's <= 0.1-0.4 dex bound is the right citation.
5. **Pool channel has one lab.** On the free + disulfide channel only Reading (Elmore, one study, huge SDs) stands;
   Hefei's is water-phase only. The ">= 2 labs on the same channel" rule is met only on the free-MFT channel
   (Huazhong + Reading), and the Huazhong cell is E12 TEST.
6. **US 8,920,862 wrong sign.** Carrying it is right; it is not grounds to reject (single run, patent, structured-lipid
   mesophase, a different matrix). But F25#17 was viewed before the structure was proposed, so it can no longer test
   it blind; the reject / keep rule over the other three E7 TEST contrasts (F30#11 kim1996b, F30#15 zamora2015,
   F30#16 hidalgo1993) must be pre-registered before any E7 engine work, with a statement of whether they were viewed.

**Corrected probabilities:** speciation shift explains the Zhang 2022 vs Xu / Reading sign: 0.65 -> **0.4**; disulfide
partition into the oil + increased formation with oil (lower pH 4.3-4.7, more xylose consumed): ~0.35; other 0.25.
P(an E7 capture term supported by >= 2 DEV labs on one channel, without TEST cells) ~0.3.

---------------------------------------------------------------------------------------------------------------------

## 5. Other items in the proposals

- CF-076 `explanations_first_replace` states "TTCA ring-N pKa 5.2-5.5" as part of the explanation: it is a fitted,
  temperature-free effective value, not a pKa; say "fitted effective pKa_T 5.2-5.5 (T-independent in the battery)".
- CF-076 `remedy_append`: "the declared band k_r(100 C, pH 4.5) 2-8e-4 /min stands only with ..." is inconsistent
  with the run's own feasible points (k_r45 0.88-3.2e-4; four of six below 2e-4). Replace by the run's range.
- CF-076 `tests_append` "rejected at Ea 50 and 64" must say "E4c"; E2c at Ea 48.5 is feasible (own check).
- CF-092 `kind: code_mismatch`: no code is involved; it is a citation / documentation error (use the register's
  citation kind if one exists).
- CF-092 explanation 3 ("dihydro products vs none: conditions P 0.6") is reasonable; Whitfield 1999 itself gives the
  same reason (concentrations, 100-fold H2S).
- No proposal moves the engine. But CF-031 `remedy_replace` specifies the E7 term in a form chosen after viewing TEST
  rows (Zhang 2021, F25#17) and an E12 TEST cell (Xu 2011): as written it would steer an engine change without a
  clean pre-registered test.

---------------------------------------------------------------------------------------------------------------------

## REQUIRED changes before applying

1. CF-076: restate the condition as the (k25, Ea) ridge: printed Karpinski rate -> Ea ~26; mid-band -> ~32; x6
   deflated floor -> ~41 (E4c) to ~49 (E2c). Say the x6 deflation is not self-consistent with E2c (no catalysis) and
   exceeds E4c's own fcat (~1.6); P 0.35 -> 0.25.
2. CF-076: scope "rejected at Ea >= 50" to E4c; E2c at the floor is feasible to Ea ~48.5 (inconclusive at 50).
3. CF-076: drop "the pKa is no longer the binding unknown"; add the T-independence of pKa_T and the 40 kJ/mol
   Delta H_ion the engine itself uses; replace "2-8e-4" by the run's 0.9-3.2e-4.
4. CF-092 / corrections B: add van den Ouweland Fig. 4 (postulates NF -> MFT); correct the Whitfield 1999 reading
   (intro AND discussion attribute MFT to them; Table 1's column is a spectral reference); keep P(peak 3 is MFT) at 0.3;
   reword every erratum line from "do not report MFT" to "propose the NF -> MFT pathway (Fig. 4) but identify no MFT
   (Table I lists 4-mercapto-2-methylfuran)"; re-split CF-092's explanations (0.65 / 0.3 / 0.05).
5. CF-092 remedy: the backlog source 03_thiamine_degradation.md does not exist; replace "fix at source and
   regenerate" with an annotation plan that survives regeneration (or record the stale entry).
6. corrections A (Wei): replace "C-1 methyl oxidised back" by "rehydration (reverse of the 2,3-enolisation /
   beta-elimination) back to the pentose level; redox-neutral".
7. CF-031: scope the pool argument to the water phase; record the oil-phase MFT and disulfide values; record that the
   pre-registered within-x3 rule fails on the ENTH arms; P 0.65 -> 0.4.
8. CF-031 / E7: remove Zhang 2021 (CMRC, keyed TEST for E7) from the evidence and from `tests_append`; record that its
   rows and F25#17 were viewed during resolve2 (mark them in the split as seen-by-designer, or exclude them from the
   E7 score); flag that Xu 2011's lard contrast is E12 TEST F28#20 and must not calibrate an E7 term.
9. CF-031 remedy: "no partition term" -> "not needed for the sign data on disk"; require a pre-registration of the E7
   term with its reject rule on the remaining unviewed E7 TEST contrasts before any engine change.

## Optional

- Run cluster-5's original bands at maxiter 20 (true search-depth check) and a second DE seed.
- Profile the printed-rate (k25 1.8e-4) E2c family for its upper feasible Ea (only the optimum 26.5 is known).
- Wei split 0.6 / 0.1 / 0.3; note C-2 epimers of the regenerated pentose are isobaric.
- Cite E1c / E5 / E6 (cluster 5) as the battery's flexible negative controls.
- Ask whether peak 3's reference compound in van den Ouweland was synthesised by a route that fixes the regiochemistry
  (not printed); an author / Unilever archive question only.
