# Proposed corrections to our own dossiers / atlas / docs (NOT applied; for the writers)

## A. wei2022 (Q3): TTCA is C8H15NO6S, M 253

Arithmetic (prereg sec. 3, done before re-reading the dossier): Cys C3H7NO2S 121.16 + xylose C5H10O5 150.13 - H2O
= C8H15NO6S 253.27 = TTCA (Wei SI Fig. S3 base ion m/z 254 [M+H]+; Zhai 2021 SI prints the same). Deoxyosone C5H8O4
132.11 + Cys - H2O = C8H13NO5S 235.26 ([M+H]+ 236, present at ~25 % in Fig. S3). Deoxyosone + Cys WITHOUT water loss
= C8H15NO6S, TTCA's formula. So the "2 H short" argument is wrong and must go; but a STRUCTURAL argument, independent
of the formula, replaces it and is stronger:
- 1-DX (CH3-CO-CO-CHOH-CH2OH) has no C-1 aldehyde: any Cys ring on it is 2,2-disubstituted (2-methyl), never TTCA's
  2-H, 2-(tetrahydroxybutyl) ring. 1-DX -> TTCA needs the C-1 methyl oxidised back: not a Maillard step.
- 3-DX (OHC-CO-CH2-CHOH-CH2OH) + Cys gives a 2-(2-oxo-...) ring (M 235); reaching TTCA needs C-2' reduced AND C-3'
  re-hydroxylated, i.e. water added back across the 2,3-enol: the reverse of the beta-elimination that made 3-DX,
  which regenerates the aldose-level 1,2-enediol (that is, the sugar). Any DX -> TTCA path therefore runs through the
  pentose oxidation level: Wei's TTCA is evidence of pentose-level carbonyl (ARP reversion / xylose / Schiff base), not
  of a deoxyosone trap product.
P (re-stated, replacing the dossier's 0.65 -> 0.55): TTCA mainly via ARP reversion / pentose 0.6; via DX 0.15; the
m/z 254 peak is not (only) TTCA (Cys-ARP has the same formula; standard undescribed) 0.25.

**wei2022_extraction.md sec. 1 table row (e)** -> replace the cell with:
"**Not shown.** TTCA identified by UPLC-Q-TOF (SI Fig. S3: base ion m/z 254, [M+H]+ of C8H15NO6S, M 253 = TTCA's
formula). Scheme 2 makes TTCA from 1-DX and 3-DX. The formula cannot reject it (DX + Cys without water loss is also
C8H15NO6S), but the structure does: 1-DX has no C-1 aldehyde (its Cys ring would be 2-methyl-2-substituted, not
TTCA), and 3-DX -> TTCA needs C-3 re-hydroxylated, the reverse of the dehydration that formed 3-DX, i.e. a return to
the pentose level. No labelling, no fed-DX control. (Corrected 2026-11-01: the earlier text gave TTCA as C8H15NO5S,
M 237, and argued a missing 2 H; both wrong.)" Tags [P] / [Z].

**sec. 1 point 2** -> "2. **The TTCA here most likely does not come from the deoxyosones** (my reading, **~0.6**;
revised 2026-11-01): Scheme 2 is structurally implausible (sec. 1 (e)): 1-DX cannot give TTCA's ring and 3-DX can do so
only by reverting to the pentose level. A mass-balanced route exists: ..." (rest unchanged). This keeps the atlas F03
quote verbatim.

**sec. 6 Flags**, replace "[!] **Scheme 2 is not redox-balanced** (deoxyosone + Cys cannot give TTCA without a
reduction); TTCA from DX is asserted, not traced" with "[!] **Scheme 2 is structurally implausible** (1-DX lacks the
C-1 aldehyde; 3-DX must be re-hydrated at C-3); TTCA from DX is asserted, not traced (no labelling, no fed-DX
control). Formula-wise it balances (DX + Cys - 0 H2O = C8H15NO6S)."

**sec. 5** "it would show whether the product is m/z 238 ([M+H]+ of TTCA) or 236" -> "it shows m/z 254 ([M+H]+ of
TTCA, M 253) with 236 at ~25 % (TTCA - H2O in source, or a co-eluting DX-Cys ring, M 235): see the SI section".

**Atlas F04** (Wei cell, contrast) -> "fed Gly-Xyl Amadori 5 mM + Cys 5 mM, temperature: TTCA appears at 5.0 / 9.6 /
22.4 % of the ARP by 100 min; the source carbonyl (deoxyosone trap vs ARP reversal to xylose) is not shown; TTCA
identity by m/z 254 (M 253, C8H15NO6S; the Cys-ARP has the same formula) against an undescribed standard"

**Atlas F03** (Wei cell, contrast) -> "GX-ARP 5 mM +- cysteine 5 mM, initial pH 7.0: TTCA appears; the authors route it
through the deoxyosones, the dossier through ARP reversion to xylose (1-DX cannot give TTCA's 2-H ring; 3-DX only by
re-hydration at C-3, i.e. back to the pentose level); free xylose 'not observed' (no detection limit)";
`where: "sec. 1 point 2 (dossier reading ~0.6)"`.

**Consequence for "cysteine traps deoxyosones -> TTCA":** none of the TTCA leg ever rested on a measurement; with the
mass corrected it still gets no support from Wei (structure, not formula). The remaining DX-trap evidence is the DX
drop with Cys (net levels; CF-071: not thiol-diagnostic, P 0.8). Engine: no change (TTCA forms only from PENT + Cys,
r_ttca_form; there is no DX -> TTCA reaction to remove). If Wei's TTCA rows are keyed, tag them
`source_unattributed: ARP+Cys` and never as a DX-trap rate.

## B. van den Ouweland 1975 (Q2)

**What the print says (read myself from scratchpad ing1101/vohi-2r.png, Table I p. 502):** peak 3 is named
"4-Mercapto-2-methylfuran" and DRAWN with SH on the ring carbon remote from the methyl (same drawing convention as
peak 8, 4-mercapto-2-methyl-2,3-dihydrothiophene); the same table names and draws the adjacent isomer correctly for
peaks 4 (3-mercapto-2-methyl-4,5-dihydrofuran), 10 and 11 (3-mercapto-2-methylthiophene / -2,3-dihydrothiophene),
and lists BOTH thiophene regio-isomers (peak 10, 3-SH; peak 13, 4-SH). Footnote a: identity checked against a
synthesised reference. MS 114 (100) / 113 (50) / 85 (60): [M-H]+/M+ 0.50 vs MFT's 0.28 (NIST, as used in CF-086).
Name, drawing, the authors' own regio-distinction and the MS agree: **not a misnaming in the print** (P that the
print intends MFT <= 0.05). Whether the authors' identification was right (an MFT peak assigned to a synthesised
4-isomer on packed GC): P ~0.15, not 0.3.

**What later work cites it for:** Hofmann & Schieberle 1998 (Garching) state that van den Ouweland & Peer could not
identify MFT in NF systems (hofmann1998.pdf p. 1, intro; Zeiler-Hilgart 1994 then confirmed the NF / Cys system is a
weak MFT source). Hofmann 1995 cites them for 2-methyl-3-thiophenethiol (Table 1 footnote e, ref 2), and cites MFT to
refs 4 and 6. Whitfield & Mottram 1999's intro says they identified MFT (a mis-summary), but Whitfield's own Table 1
attributes MFT to Evers et al. 1976 and cites van den Ouweland only for 2-methyl-3-thiophenethiol, and its discussion
cites them for the dihydro derivatives. So: the primary print and two later labs agree there is no MFT in van den
Ouweland; one intro sentence (Whitfield 1999) and our legacy docs mis-cite it.

**Does the NF -> MFT route depend on it? No** (P 0.97). Fed-NF MFT on disk, four laboratories: Nestle (Cerny 2003
pot D: unlabelled MFT from fed NF against [13C5]ribose, labelling, 95 C, 7 % of MFT), DFA Garching (Hofmann 1998
Table 4: NF + H2S, 145 C, 0.19 mol %), Reading (Whitfield & Mottram 1999: NF + H2S / Cys, 140 C, 0.12-0.15 mol %;
Whitfield 2001), Givaudan (de Roos 2005: fed NF + Cys, 0.9 mg MFT / mol Cys). E6 is about the NF SHARE, which van den
Ouweland never bounded (no yields): no change to E6, its TEST contrasts (BTBU, CAAS) or k_nf_mft.

**Text to correct (erratum notes, history files are not rewritten):**
- `docs/validation/isotope_topology_evidence.md` row 12 "reacts with H2S (van den Ouweland & Peer's original
  synthesis)" -> "reacts with H2S to MFT (Hofmann 1998, Whitfield 1999, Cerny 2003 pot D; van den Ouweland & Peer
  1975 report 2-methyl-4-furanthiol and dihydro-MFT, NOT MFT)"; line ~100 "cited to van den Ouweland & Peer 1975" ->
  add "(a mis-citation: the paper reports no MFT; see vandenouweland1975_extraction.md)".
- `docs/history/AUDIT_legacy_lane_2026-08.md` (~l. 358-360), `docs/history/README_legacy_lane_2026-09-03.md`
  (~l. 391-392), `docs/history/tasks/audit_remediation.md` (~l. 480, 1633, 7488 "Van den Ouweland & Peer 1975 ...
  report the same synthesis"): append one erratum line each: "Erratum 2026-11-01: van den Ouweland & Peer 1975 do not
  report MFT (Table I: 4-mercapto-2-methylfuran, dihydro-MFT, 3-mercapto-2-methylthiophene); the NF -> MFT synthesis
  is Whitfield et al. 1993 / Mottram & Whitfield 1994 (cited by Hofmann 1998), measured on disk by Hofmann 1998,
  Whitfield 1999, Cerny 2003, de Roos 2005."
- `results/literature/deep_research_backlog.json` entry "van den Ouweland & Peer (1975)": description "MFT in cooked
  pork by SIDA; species-level thiamine-MFT correlation" is wrong for this citation (a generated artifact: fix in its
  source file 03_thiamine_degradation.md line 590 and regenerate; do not hand-edit).
- `data/lit/extraction_dossiers/whitfield1999_extraction.md`: add a flag "[!] the intro says van den Ouweland & Peer
  identified MFT; their Table I does not (they report 2-methyl-4-furanthiol); Whitfield's own Table 1 cites MFT to
  Evers 1976".
- `data/lit/extraction_dossiers/vandenouweland1975_extraction.md` sec. 1 answer 1: "P ~0.3" -> "P ~0.15 (not a
  misnaming: name, drawing and the table's own 3- / 4-isomer distinction agree; a misidentification on packed GC
  remains possible)"; add Hofmann 1998's reading ("unable to identify MFT").

## C. Zhang 2022 / CF-031 (Q4): which MFT species each lab measured

| lab (study) | split role for E7 | MFT species measured | method | lipid effect on what was measured |
|---|---|---|---|---|
| Huazhong (Xu 2011) | DEV | free MFT only (no disulfide in Table 1) | whole-vial HS-SPME 60 C, lard in the vial | x1/21 fresh, x1/143 oxidised (FFT x1/3, 1/6.4) |
| Reading (Elmore 2002) | DEV | free MFT AND bis(2-methyl-3-furyl) disulfide | HS-SPME RT | free 807 -> 58 (linoleate) -> tr; disulfide 202 -> tr -> ND; **pool x1/17 and < 1/100**; FFT 92 -> 87 -> 17 |
| Reading (Whitfield 1988) | DEV | free MFT | Tenax headspace | -0.47 dex (lecithin) |
| Reading (Yiltirak 2026) | DEV | free MFT only (PFBBr in 1 M NaOH; no reductant) | SIDA | buffer + oil x1/9.6; emulsion x0.5-0.9 |
| Hefei (Zhang 2022) | new lab | free MFT AND the disulfide, water phase separated | HS-SPME, RF 1 | free x0.36-4.8 (up); disulfide 2.81 -> ND (RM); **pool x1/3.4-1/25 (RM), x1/1.0-1/3.5 (ENTH)**; FFT x0.93-6.5 |
| CMRC Beijing (Zhang 2021, keyed TEST) | TEST | free MFT only (no disulfide in Table 3) | whole-sample HS-SPME 50 C | MFT x0.60 / 0.36 / 0.43 (lecithin / cephalin / both, Duncan-different from A1), reheated x0.51 / 0.53 / 0.31; FFT x1.15 / 0.70 / 1.34, reheated x1.71 / 1.48 / 1.83 |
| Nestle (US 8,920,862) | TEST (E7 F25#17) | free MFT | SIDA-SPME | structured-lipid mesophase x10 (single run) |

Reading: the ONLY free-MFT rise against lipid in a same-chemistry pot (Xyl or Rib / Cys, 120-140 C) is Zhang 2022's
water phase, and Zhang's control holds 95 % of its MFT pool as the disulfide (2.81 vs 0.14 ng/mL; unbuffered pH 7.2
-> 4.7, air in a 50 mL tube). With oil the disulfide vanishes from the water phase, so free = pool there, and the
"x3.6 rise" is the pool's x1/5.8 fall divided by the control's 0.05 free fraction. Both labs that measured the
disulfide (Reading / Elmore, Hefei / Zhang 2022) show the POOL falling by comparable factors (x1/17; median x1/13 for
raw oils), and in Xu's buffered pH 5 vial (disulfide unreported) the free MFT falls x1/21 with fresh lard: same order.
The sign disagreement is a measurement-species artefact (speciation shift free <-> disulfide), not chemistry.
Residual, not explained by species: the oxidation direction (Xu and Elmore: more oxidised / more unsaturated -> more
suppression; Zhang ENTH -> less). ENTH is lipase-hydrolysed (free fatty acids) and only mildly oxidised (PV <= 4.3 vs
Xu's 184 meq/kg); du2025 (CAAS) also sees free fatty acids RAISE MFT (+0.10 dex): a hydrolysis effect, P ~0.5.

Partition is not needed and not supported as a separate E7 term: Zhang 2021 lowers MFT x0.4-0.6 with 0.07-0.24 %
phospholipid (no bulk oil phase to partition into); Zhang 2022's water-phase FFT RISES under 4:1 oil (pure partition
would cut it x1/20-1/400); cluster 4 bounded a headspace partition operator at <= 0.1-0.4 dex.

**Minimal E7 structure supported by >= 2 labs on the same channel:** a lipid-electrophile capture of the MFT POOL
(free + disulfide), first order in MFT and in lipid-carbonyl load, MFT-selective over FFT (log suppression ratio MFT
/ FFT >= ~3: Xu 1.32 vs 0.48 dex; Elmore 1.24 vs 0.02; Zhang 2022 pool 0.5-1.4 vs FFT +0.03 to -0.64 (raw oils); Zhang 2021 MFT
0.22-0.51 vs FFT -0.26 to +0.15). Four labs (Huazhong, Reading, Hefei DEV-side reading; CMRC TEST consistent, 6 / 6
contrasts MFT falls and falls more than FFT, a post hoc sign check, not a score). NOT supported cross-lab: an
oil-phase partition term; an oxidation-load dose law (2 labs for, 1 confounded against); a "protection" term raising
free MFT (Zhang 2022 water phase + the Nestle mesophase only, and the Nestle row is TEST: it may not build the
structure). Consequence: this structure predicts the TEST contrast F25#17@us8920862 (MFT x10 up) WRONG in sign;
that failure must be carried, not tuned away. Observation operator: every E7 row tagged `mft_species: free | pool`;
compare free rows with free MFT only when the pot's disulfide share is known or bounded; otherwise score the sign of
the pool only where the disulfide was measured.
