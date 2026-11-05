# Canle, Lawley, McManus & More O'Ferrall 1996 — EXTRACTION
DOI: 10.1351/pac199668040813 (not printed in the PDF; Crossref 2026-11-05: "Rate and equilibrium constants for oxazolidine and thiazolidine ring-opening reactions", Canle, Lawley, McManus, O'Ferrall, Pure Appl. Chem. 68(4), 813-818, 1996. **Matches the brief and the PDF.**)
### A short conference paper (6 pp.). Ring-chain (ring ⇌ iminium) equilibria and rates at 25 °C in water for 2-[4-(dimethylamino)phenyl] 1,3-OXAZOLIDINES (N-alkyl / N-phenyl) and for the N-H 1,3-THIAZOLIDINE of 4-(dimethylamino)benzaldehyde, compared with Fife's N-butyl / N-phenyl cinnamaldehyde thiazolidines. For the N-H thiazolidine: apparent ring-closing pKa 5.7, at most ~50 % ring-opened (between pH ~4 and ~6), ring opening k1 = 2.5 × 10⁵ s⁻¹. **No cysteine, no carboxylate, no sugar, no buffer catalysis, 25 °C only.**

**Read 2026-11-05 from `data/articles/canle1996.pdf`** (6 pp., Pure Appl. Chem. pp. 813-818; publisher PDF with a text
layer). Prose from `pdftotext -layout`; **Tables 1-2, Scheme 1, eqn. (3) and Figs. 1-3 read from 150 dpi renders**
(scratchpad `img/can-2.png` … `can-5.png`), because the text layer drops the exponents ("25x105" is 2.5 × 10⁵).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Chemistry, University College Dublin (More O'Ferrall; physical organic). Not in
`lab_ids.yml`; proposed id `ucd_more_oferrall`. More O'Ferrall co-authors `davis1991b_extraction.md` (Huddersfield +
Dublin), so the two share a lab link under merge rule 2. The thiazolidine comparison rows come from Fife 1991 (ref 11,
`fife1991_extraction.md`), a third laboratory.

## 1. What this paper answers for the model

The case (M1c sec. 4.3) asked for "the shape of P1a's K and of ring opening": how the open-chain (iminium / Schiff
base) fraction of a thiazolidine and its ring-opening rate depend on pH.

| model need | what Canle 1996 gives | tag |
|---|---|---|
| open-chain fraction vs pH, a thiazolidine | N-H 2-[4-(dimethylamino)phenyl]thiazolidine (2 × 10⁻⁴ M): iminium absorbance at 400 nm rises from pH 7 to 5 and falls again from pH 4 to 2 (Fig. 2); greatest ring opening "about 50 %" | [F] / [M] |
| ring-chain pKa values | Scheme 1: 5.7 (ring 9 ⇌ open iminium 10, and 9 ⇌ 12) and 3.9 (10 ⇌ diprotonated ring 11, and 12 ⇌ 11) | [M] |
| ring-opening rate | **k1 = 2.5 × 10⁵ s⁻¹** (N-H, benzaldehyde-derived, 25 °C), from the measured ring closure **k−1·Ka = 0.5 s⁻¹** and pKa_app 5.7 | [M] |
| N-substituent effect (Fife's cinnamaldehyde rings) | N-Bu k1 1.4 s⁻¹, N-Ph 5.7 × 10⁻³ s⁻¹; pKa_app 6.2 / 6.3 / 7.1 (N-H / N-Bu / N-Ph) | [C] (printed in Table 2) |
| hydroxide attack on the iminium (= the hydrolysis entry) | 5.5 × 10³ (N-H, B); 1.3 × 10⁴ (N-H, C); 9.1 × 10³ (N-Bu); 1.7 × 10² (N-Ph) M⁻¹ s⁻¹ | [M] / [C] |
| buffer catalysis | **none in this paper** (the brief and queue mention it; it is not here) | [!] |
| anything for a sugar or cysteine ring, any T ≠ 25 °C | **none** | — |

## 2. Conditions

| item | value |
|---|---|
| Substrates | oxazolidines 3: 2-[4-(dimethylamino)phenyl]-N-R-1,3-oxazolidine, R = Me, Et, Bu, iPr, tBu, Ph; thiazolidine 9: N-H 2-[4-(dimethylamino)phenyl]-1,3-thiazolidine (from cysteamine-type amino thiol; no 4-carboxylate) at 2.0 × 10⁻⁴ M (Fig. 2); Fife's p-(dimethylamino)cinnamaldehyde thiazolidines 8 (R = H, Bu, Ph) for comparison |
| Medium | aqueous, **25 °C**; buffers, ionic strength and pH scale **not stated** in this paper |
| Methods | UV-vis (iminium at 400 nm); stopped-flow; ring closure measured by quenching the ring-opened ion (pH 5-7) into strong acid (pH 0-2; Fig. 3 plots k_obs vs [HCl] 0.05-0.4 M) |
| Errors | **none printed** for any constant |

## 3. Tables and figures re-typed

### Table 1 (p. 814). N-substituted 2-[4-(dimethylamino)phenyl]-1,3-oxazolidines, 25 °C [M]

| N-substituent | Me | Et | Bu | iPr | tBu | Ph |
|---|---|---|---|---|---|---|
| pKa_app (ring closing) | 9.3 | 9.1 | 9.3 | 9.1 | 8.4 | 6.9 |
| k1 / s⁻¹ (ring opening) | — | 0.7 | 1.4 | — | 0.35 | 0.025 |

Fig. 1 (p. 814), N-Bu oxazolidine: upward break at pKa_app 9.28 [M], downward break at pH 6.5 (ring protonation) [M].
Derived: hydroxide-catalysed ring closure 4 × 10⁴ M⁻¹ s⁻¹, external OH⁻ on C=N 1.3 × 10³ M⁻¹ s⁻¹, k_intra/k_inter =
25 vs an equilibrium analogue of 50 [M].

### Scheme 1 and Fig. 2 (p. 816). N-H thiazolidine 9 of 4-(dimethylamino)benzaldehyde, 25 °C

- Fig. 2 [F]: A400 ≈ 0.01-0.03 at pH 1-3, maximum ≈ 0.50 near pH 5.3-5.5, ≈ 0.20-0.27 at pH 6.5, ≈ 0.08 at pH
  7.5, ≈ 0.02 at pH ≥ 9 (read error ± 0.02 A, ± 0.2 pH).
- Scheme 1 [M]: pKa 5.7 on the two edges from the neutral ring 9 (to the iminium-thiol 10 and to the
  dimethylamino-protonated ring 12); pKa 3.9 on the two edges into the diprotonated ring 11.
- Maximum degree of ring opening "about 50 %" [M], from the apparent extinction coefficient against the model
  iminium ion 14.

### Eqn. (3) and Fig. 3 (p. 816). Ring closure of the iminium-thiol 10

Ring closure shows normal and inverse [H⁺] dependences in 0.05-0.4 M HCl (Fig. 3, k_obs minimum ≈ 22 at [HCl] ≈ 0.15
M, axis unit not printed [F, !]); the inverse term is pre-equilibrium thiol ionisation (10 → 13) then closure of the
zwitterion 13 → 9 [P]. **k−1·Ka = 0.5 s⁻¹** [M].

### Table 2 (p. 817). Ring opening of N-substituted thiazolidines, 25 °C

| N-sub | aldehyde | pKa_app | k1 / s⁻¹ | k_OH / M⁻¹ s⁻¹ |
|---|---|---|---|---|
| H | B (4-Me2N-benzaldehyde) | 5.7 | **2.5 × 10⁵** | 5.5 × 10³ |
| H | C (4-Me2N-cinnamaldehyde) | 6.2 (Fife, ref 11) | — | 1.3 × 10⁴ |
| Bu | C | 6.3 | 1.4 | 9.1 × 10³ |
| Ph | C | 7.1 | 5.7 × 10⁻³ | 1.7 × 10² |

Text: k1 varies "nearly 10⁸" from H to Ph while pKa_app varies < 30-fold [M/Q]. (2.5 × 10⁵ / 5.7 × 10⁻³ = 4.4 × 10⁷
[Z].)

## 4. Arithmetic (mine) [Z]

**1. Ring closure.** k−1 = 0.5 s⁻¹ / Ka(thiol of 10). With the thiol pKa of the analogous iminium estimated at 9.7
(the authors' ref 7, cited on p. 815), k−1 ≈ 0.5 × 10^9.7 ≈ 2.5 × 10⁹ s⁻¹: ring closure of the thiolate is
near the diffusion / conformational limit, consistent with Davis 1991's > 10⁸ s⁻¹ (`davis1991_extraction.md`).

**2. pH shape of the open fraction.** Above pKa_app (5.7) the open iminium fraction falls ~1 dex per pH unit
(f_open ∝ [H⁺]); Fig. 2's fall from ≈ 0.50 at pH 5.4 to ≈ 0.08 at pH 7.5 is −0.38 dex/pH over that span, flattened by
the ~50 % cap near the maximum. Below pH ~4 the dimethylamino protonation closes the ring again; **that limb is
specific to the 4-Me2N-aryl substituent** and has no TTCA counterpart. So for a thiazolidine with no such group,
the transferable statement is: **the open-chain iminium is favoured by acid and disfavoured by base above the ring-N
pKa**, i.e. a ring-chain K_open that falls ~1 dex per pH unit there.

**3. Why the magnitude does not transfer.** The 4-Me2N-phenyl group stabilises the iminium cation by resonance (the
authors note the thiazolidine ring is ~10⁵ more favoured relative to the open form than the oxazolidine, p. 815).
A polyhydroxyalkyl C-2 has no such stabilisation; Gajda 1990 and Radomski 1989 see no open chain by NMR. The open
fraction of TTCA is therefore orders of magnitude below 50 %; this paper cannot say how far.

## 5. What it can and cannot test in the model

- **Can:** give the class SHAPE of the ring-chain equilibrium (open iminium favoured on the acid side of the ring-N
  pKa) and the class SPEED of the ring-chain step (opening 10⁵ s⁻¹ for an N-H aryl ring at 25 °C; closure ≥ 10⁹
  s⁻¹): on any pot timescale the P1a / ring exchange is at equilibrium, never rate-limiting.
- **Cannot:** give f_oc for TTCA, any K with free cysteine + aldehyde, any hydrolysis (release) rate constant for a
  ring (only the OH⁻ attack on the pre-formed iminium), any temperature dependence, any buffer term.

## 6. Conflicts with existing records

| record | what it says | what Canle says | verdict |
|---|---|---|---|
| `fife1991_extraction.md` (ring opening fast and reversible; K_eq 5.3 × 10⁻⁷ M for ring + H3O⁺ ⇌ iminium-SH, compound II) | pK 6.28 | Table 2 quotes Fife's N-H cinnamaldehyde ring at pKa_app 6.2 | **consistent** (same compound, same quantity) |
| `fife1991_extraction.md` (general-base catalysis of hydrolysis by phosphate) | buffer term measured | no buffer term here | no conflict; this paper simply does not address it |
| `nagasawa1981_extraction.md` (ring-chain step hours-scale at room T, penicillamine + acetaldehyde, acid-catalysed) | slow | 2.5 × 10⁵ s⁻¹ for an N-H aryl ring | not a conflict: aryl iminium stabilisation; Nagasawa's aliphatic, gem-dimethyl ring opens ~10⁹-fold slower. Both say acid favours the open form |
| `davis1991_extraction.md` (penicilloate ring opening 5 × 10⁻⁵ s⁻¹, 30 °C) | slow, aliphatic C-2 | fast, aryl C-2 | same direction as above: the C-2 substituent sets the rate over ~10 orders |
| `radomski1989_extraction.md`, `gajda1990_extraction.md` (no open chain seen) | f_oc ~0 | up to 50 % for the aryl ring | consistent: different class |

## 7. Flags

1. **[!] Not N-aryl, and no buffer catalysis.** The download queue described "N-aryl thiazolidines … with buffer
   catalysis". The thiazolidine studied here is N-H, 2-ARYL; N-aryl (N-phenyl) values are Fife's, quoted; no buffer
   catalysis is reported.
2. **[!] No experimental section, no errors, no buffers, no ionic strength.** A conference summary; every constant
   is a single printed value.
3. Fig. 3's y-axis has no unit (probably s⁻¹).
4. The ring opening k1 is derived (k−1·Ka and pKa_app combined), not measured directly ("too rapid for kinetic
   measurements", p. 816 [Q]).

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source, distant class (the queue said N-aryl): the shape of P1a's K and
of ring opening.

**Band it can supply:**
- **Coordinate 1, P1a's K(pH) shape:** the open-chain iminium fraction is maximal on the ACID side of the ring-N pKa
  and falls ~1 dex per pH unit above it (N-H 2-aryl thiazolidine, 25 °C, pKa_app 5.7; maximum ~50 % at pH ~5.4).
  For TTCA this argues that f_oc, if P1a is non-zero, should DECREASE from label pH 5.5 to 8 and from Cao's pH 4.4 to
  Zhang's 6.5, not increase. It gives no magnitude: the aryl substituent stabilises the iminium by orders of
  magnitude that TTCA lacks, so the 0-0.05 band (Radomski, Karpinski) is not widened or narrowed by it.
- **Coordinate 2, P1a / P2 exchange speed:** ring opening 2.5 × 10⁵ s⁻¹, closure ~10⁹ s⁻¹ (25 °C, aryl): a class
  statement that the ring-chain step is never rate-limiting. For TTCA's aliphatic class the slower analogues
  (Nagasawa; Davis 1991: ~3 × 10⁻³ min⁻¹ at 30 °C) set the band; this paper only confirms the direction.
- **Transfer:** distant (aryl C-2, no carboxylate, no sugar, 25 °C). Shape only.
- **Not decisive:** 25 °C, no cysteine-pentose ring. Nothing the case would call decisive.
