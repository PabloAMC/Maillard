"""
src/m4_reduced_power.py -- case M4 sec. 5.4: the reduced forced-input mini-model for the power check.

Design: docs/reference/cases/M4_glucose_mft_protein_amine.md (draft PR #90, branch m4-case-ready) sec. 5.4, with the
rivals of sec. 3, the decisive units of sec. 4, the measurement model of sec. 5.1 and the bands of sec. 5.2. This is
NEW code in the style of src/m3_strecker_layer.py (forced inputs, rival classes, draw sets, a power harness). It
imports no shipped lane: the five shipped constants it needs (R0's C2 + C3 entry, B21's Amadori -> glucosone and
glucosone -> glyoxal, the trunk's amine-free Fru -> 1-DG and Glc -> 3-DG, and 1-DG's acetic-acid exit) are FROZEN
LITERALS below, asserted equal to the engine by tests/unit/test_m4_reduced_power.py.

SYNTHETIC ONLY. No measured value of any decisive unit is read: D1 / D2 (Hofmann 1998 Tables 2 / 5), D3 (Du 2023,
Zhao 2019b; M2's TEST), the screen's pots and every reserved set are never opened. The only printed numbers used are
band sources quoted in the case (sec. 5.2) and, for job (a), D3's printed "<= 1 %" mixed-channel cell as the case
itself quotes it.

THE POT (one well-mixed aqueous pot, isothermal, label pH held: every design is 0.5 M phosphate)
=============================================================================================
States (mM): glucose G (forced: held at charge apart from the modelled exits), free Cys, the Cys-hexose Amadori pool AP,
1-DG, glucosone GOS, 3-DG TDG, 1-DPO by carbon family (K = C-1 kept, from 1-DG; L = C-1 lost, from glucosone), 3-DPO,
NF by family, MFT by family (K, L, D = R2's DDP route), R2's 1,4-dideoxy C6 (C6DD) and DDP, H2S (forced: made in
proportion to Cys loss, declared yield), formic acid, and per-reaction sink ledgers (carbon and label positions).
Steps (sec. 5.4 item 2):
  G + Cys -> AP                     formation, Li 2011-centred (x0.1 - 3, FM9), Ea 86 - 136, pH slope s_form
  AP -> 1-DG + Cys(rho) | 3-DG + Cys(rho)     decay, Deng 2021-centred (x0.3 - 3), Ea 34 - 107, split phi(pH)
  AP -> GOS + Cys (x O2)            the class's oxidative exit, 0 to B21's k_ama_g + 1.5 dex (FM6)
  G -> 1-DG (amine-free; x cat in R4'), G -> 3-DG          the trunk's shipped amine-free steps
  1-DG -> other (x O2 drain 1 - 4.4) | 1-DPO_K + C1 (R1: C-6 lost) | + amino acid -> C6DD (R2)
  GOS -> glyoxal + C4 (B21 shipped) | formic (C-1) + C5 enediol -> 1-DPO_L / 3-DPO (R5; Glomb ceiling, Ea 40 - 135)
  1-DPO -> NF | -> MFT (direct, x h(H2S)) | other;  NF -> MFT (x h(H2S)) | other;  MFT -> loss   (M2's shares)
  C6DD -> DDP + C1 -> MFT_D (x h(H2S)) | other                                                       (R2)
  Cys -> H2S (y) + sink;  H2S -> loss
R0's C2 + C3 route is carried as a SINK-FREE UPPER BOUND in every rival (k_glc_ha at its shipped value times Hofmann
1998 Table 7's and Table 8's largest measured yields, Table 6's for NF): it favours R0 at every MGO-law corner, so the
MGO law (Ea 47 / 121 / 20, orders 1 / 2) is not a coordinate here (it can only lower R0's bound).

INTEGRATION: implicit (backward) Euler in topological order, every flux computed once and booked on both sides, so
mass, carbon and label positions close to rounding (check < 1e-10). The Cys / AP pair is solved jointly (2 x 2).

LABEL MAPS (sec. 5.4 item 3). MFT's carbons are read as (methyl, C2, C3, C4, C5), NF's as the 1-DPO skeleton.
  family K (1-DG -> 1-DPO; R1, R4', R2's DDP route): (1, 2, 3, 4, 5), C-6 lost
  family L (glucosone -> 1-DPO; R5):               (2, 3, 4, 5, 6), C-1 lost as formic acid
  C-1 / C-6 equilibration on the hexosulose (Tressl 1999 Fig. 14): with probability e the chain is read reversed
  (position p -> 7 - p) before the cleavage; e is integrated on a grid (U[0, 0.5]) in every rival's density.
  C2 + C3 (R0's route): MFT (methyl, C2, C3) = a C3 fragment, glucose (1, 2, 3) or (6, 5, 4) with p 1/2 each, and
  (C4, C5) = a C2 fragment, (2, 1) or (5, 6), from a DIFFERENT glucose molecule, paired at random.
Reads (each a share, 0 - 1): a single positional label's labelled fraction on the molecular ion (M) and on a
methyl-locating fragment ("ring" = [M - CH3]+, every carbon but the methyl; "acyl" = CH3CO+, methyl + ring C2);
[1,2-13C2]'s M+1 and M+2 fractions; the 1:1 [13C6] / [12C6] pool's mixed fraction (M+1 - M+4).

DECLARED CHOICES WHERE THE CASE IS SILENT (each is in the record; none sets a band from a scored pot)
====================================================================================================
D-1  Glucose is a state held at charge apart from AP formation and the two amine-free steps (sec. 5.4 item 1, "held at
     charge"); R4''s amine-free 1-DG supply is the shipped Fru -> 1-DG constant applied to the whole hexose charge (an
     upper bound on the trunk's Glc -> Fru -> 1-DG chain, favouring R4'), times the Cys-catalysis factor.
D-2  AP formation's pH law is 10^(s_form (pH - 7)), s_form on [0, 1] (0 = pH-flat, 1 = the free-amine limit); unsourced,
     required agreement at both ends.
D-3  AP decay's 1-DG : 3-DG split has the 2,3-enolisation odds slope 0.44 dex/pH (Beck / Nedvidek x12.5 over pH
     4.5 -> 7) around a reference 1-DG share phi on [0.2, 0.6] at pH 5.5 (the trunk's 0.59 at pH 5 inside); R5's C5
     enediol uses the same split for 1-DPO : 3-DPO ("the network's pentose enediol split", sec. 5.2).
D-4  Every pool rival carries the class's oxidative exit AP -> GOS + Cys (the B21 step), banded as sec. 5.2 bands it for
     R5; the rivals differ in what glucosone and 1-DG feed: R1 has r_odg_dpo on and glucosone -> C5 off, R5 the
     reverse, the claim both. R4''s pool has no exit but its inert sink (no 1-DG, 3-DG or glucosone from AP).
D-5  1-DG's other exits are centred on the shipped r_odg_aa (1.45 /min at 100 C; Glomb 2010's 50 C half-life
     extrapolated with its 76 kJ/mol gives 1.0, inside x0.3 - 3), pH slope s_dg1 on [0, 0.5] (rising, Davidek 2006a's
     acetic share); r_odg_dpo is pH-flat, its share at 100 C, pH 5 on [5e-4, 0.5] (log odds [-3.3, 0]; the ceiling
     is the unaccounted ~0.5 of Voigt 2009's 7 d exits), barrier 40 - 135.
D-6  R5's beta-cleavage: k_b5 = u 1.4e-3 /min at 50 C (Glomb's total glucosone loss, the ceiling), u log-uniform on
     [1e-3, 1], Ea 40 - 135; it competes with B21's shipped k_g_go on the one glucosone pool.
D-7  Downstream (M2's shares as a factor): 1-DPO total loss k_dpo [0.3, 30] /min at 145 C; NF share a_nf [0.01, 0.3];
     NF -> MFT share b [0.005, 0.1] at h = 1 (Whitfield 1999's 2.6 % inside); the share of MFT through free NF
     s_NF on [0.07, 0.999] (Cerny pot D's 7 % to E6's 99.9 %) sets the direct share a_dir = a_nf b (1 - s) / s;
     k_nf [0.01, 1], k_mft [0.001, 0.05] /min; every downstream step Ea 80 kJ/mol; h = [H2S] / (1 mM + [H2S]).
D-8  Cys -> H2S: k_cd [1e-3, 1e-2] /min at 145 C (Ea 100), yield [0.1, 0.5]; H2S loss [0.01, 0.1] /min (Ea 20).
D-9  R2: k_str 3e-5 - 3e-2 /(mM min) at 100 C (M3's dicarbonyl Strecker class band, kappa_MGO D0 10^-2.5 - 10^0.5 at
     100 mM), Ea 65 - 85 (M3's exits), pH slope [0, 0.5]; the amino acid is Cys + Gly and is not consumed (it is in
     excess); C6DD -> DDP + C1 at 1 /min (145 C), DDP [0.01, 1] /min, DDP -> MFT share [0.01, 0.3] at h = 1.
D-10 3-DG loss [0.1, 10] /min at 145 C; 1-DG's other exits give formic acid at 0.03 mol/mol (Smuda 2010's 27 mmol/mol).
D-11 O2: air factor 1, N2 0 (sec. 5.4); the cysteine-budget form scavenges the headspace O2 (8.56 mM-equivalent at
     Hofmann's half fill) at 4 Cys per O2 before any oxidative step: factor max(0, 1 - Cys0 / (4 O2_0)).
D-12 The claim {R1, R5} is ONE rival: both entries on, the R5 share of 1-DPO formed in the reference cell (145 C, pH 5,
     air, 20 min, arm 1's design) drawn f* ~ U[0.001, 0.999] and the barrier difference dEa ~ U[-60, 60] (FM4); each
     draw's two entry constants are rescaled symmetrically (bisection) to hit f*, so f(T) moves as a logistic in 1/T
     where both entries are minor exits; a rescaled draw outside either entry band is redrawn (counted).
D-13 Measurement (sec. 5.1): every level carries omega 0.419 dex per value plus SIDA CV 10 % / sqrt(n), and H0's
     tau_ABS 0.964 dex as ONE offset per species shared by every value of that species in the laboratory (an engine
     level against an absolute measurement), integrated by Gauss-Hermite; shares carry 1 % / sqrt(n) per value plus
     a 0.5 % calibration floor (the label-precision limit), no omega, no tau. LODs (mM): MFT 1e-6 (~0.1 ug/L), NF 1e-5,
     OPD quinoxalines 1e-4, formic 1e-2, AP and Cys 1e-3; a share is read only where the unlabelled pot's MFT (NF) is
     >= 3 LOD. The d3 standard, when present in labelled pots, adds 10 % / sqrt(n) to the 1:1 mixed share's error. The
     OPD artefact is a glucosone term a AP(t) (a 0.002 post-incubation, 0.05 continuous), blank-subtracted with the
     blank's 10 % / sqrt(n) error.
D-14 R0's C2 + C3 bound is R0's alone. In the entry rivals the same shipped route makes ~1e-14 mM MFT (the case's
     probe C), ~8 dex below the LOD, so it is set to 0 there; carrying R0's generous bound into them would plant a
     spurious mixed share in every low-entry draw.
D-15 The claim's scoring draws take f_ref logit-uniform (mass near the R1 and R5 edges, where the label shares pin f to
     ~1 %) with the prior / proposal ratio carried; every density is the self-normalised importance estimate over the
     draw set's feasible draws, so the target is the band prior restricted to feasible draws. Truth draws come from
     the prior itself.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.special import log_ndtr

# ======================================================================================================== constants

R_GAS = 8.314462618e-3                     # kJ/(mol K)
T100, T145, T50 = 373.15, 418.15, 323.15
DT = 0.1                                   # min

# ---- shipped constants (FROZEN LITERALS; tests/unit/test_m4_reduced_power.py asserts them against the engine)
K_GLC_HA_145 = 2.575871252387084e-09       # /min at 145 C, sulfur lane k_glc_ha (Glc -> HA + MGO + C1)
EA_GLC_HA = 64.0822869313293
LOG10_K_AMA_G_100 = -2.0339128364819725     # B21 FROZEN_B21["log10_k_ama_g_100C"]
EA_AMA_G = 75.9
LOG10_K_G_GO_100 = -0.5104286909031046      # B21 aqueous glucosone -> glyoxal
EA_G_GO = 4.2
K_FRU_ODG_100, EA_FRU_ODG = 7.418950997820022e-06, 99.3   # trunk Fru -> 1-DG (amine-free)
K_GLC_TDG_100, EA_GLC_TDG = 9.39823324272715e-06, 107.2   # trunk Glc -> 3-DG (amine-free)
K_ODG_AA_100 = 1.45                         # trunk 1-DG -> acetic acid + C4 (/min at 100 C)
EA_ODG_AA = 76.0                            # its barrier (asserted to 0.1 kJ/mol against the engine)

# ---- band sources (sec. 5.2) and fixed measured constants
K_FORM_LI_100_PH7 = 0.0011 / 200.0         # /(mM min): Li 2011's free-Cys loss 0.0011 /min at 0.2 M glucose, pH 7, 100 C
K_DEC_DENG_100 = 0.0048                    # /min: Deng 2021, Met-glucose ARP, label pH 5.5, 100 C
S_DEC_DENG = math.log10(0.012 / 0.0048) / 2.0      # dex/pH from Deng's pH 7.5 arm (0.199)
S_SPLIT_ODDS = math.log10(12.5) / 2.5      # 0.44 dex/pH: Beck 1988 / Nedvidek 1992, 1-DG : 3-DG x12.5 over 4.5 -> 7
K_GOS_GLOMB_50 = 1.4e-3                    # /min: Glomb 2010, glucosone t1/2 8 h, 50 C, pH 7.4 + Lys (R5 ceiling)
KAG_CEILING_DEX = 1.5                      # FM6: R5's ceiling is the shipped k_ama_g + 1.5 dex
Y_MP_T7, Y_MFT_T8, Y_NF_T6 = 0.040, 0.0027, 0.0318   # Hofmann 1998 Tables 7 / 8 / 6 largest measured yields
Y_FA_DG1 = 0.03                            # Smuda 2010: formic 27 mmol/mol fed 1-DG
K_H2S_HALF = 1.0                           # mM, D-7's h(H2S)
EA_DOWN, EA_CD, EA_HL = 80.0, 100.0, 20.0
K_C6DD_145 = 1.0
#: mM-equivalent O2 at Hofmann's half fill (100 mL air at 25 C over 100 mL solution): 0.2095 / 24.465 mol/L = 8.56 mM
O2_MM_HALF_FILL = 0.2095 / 24.465 * 1000.0

# ---- measurement (sec. 5.1; D-13)
OMEGA_DEX = 0.419
OMEGA_SENS_DEX = 0.390
TAU_ABS_DEX = 0.964
CV_SIDA = 0.10
CV_DEX = CV_SIDA / math.log(10.0)
SHARE_SD = 0.01
SHARE_FLOOR = 0.005
LOQ_OVER_LOD = 3.0
LOD_MM = {"MFT": 1e-6, "NF": 1e-5, "DPO1": 1e-4, "DPO3": 1e-4, "DDQ": 1e-4, "GOS": 1e-4, "FA": 1e-2,
          "AP": 1e-3, "CYS": 1e-3}
OPD_ARTEFACT = {"post": 0.002, "continuous": 0.05}
D3_EXTRA_SD = 0.10
GH_NODES = 7
THRESHOLD_NATS = 1.92
E_GRID = np.linspace(0.0, 0.5, 11)

RIVALS = ("R0", "R1", "R5", "CLAIM", "R2", "R4P")
POOL_RIVALS = ("R1", "R5", "CLAIM", "R2", "R4P")
RIVAL_LABEL = {"R0": "R0", "R1": "R1", "R5": "R5", "CLAIM": "claim {R1, R5}", "R2": "R2", "R4P": "R4'"}


def arr(ea, t_k: float, t_ref: float):
    return np.exp(-np.asarray(ea, dtype=float) / R_GAS * (1.0 / t_k - 1.0 / t_ref))


# ============================================================================================ coordinates and forms

@dataclass(frozen=True)
class Coord:
    name: str
    lo: float
    hi: float
    rivals: Tuple[str, ...]
    basis: str
    required_agreement: bool = False       # unsourced width: the verdict must hold at both ends


_POOL = POOL_RIVALS
_NUIS = RIVALS                             # nuisance coordinates every rival draws (R0 included: Cys, H2S, sinks)
COORDS: Tuple[Coord, ...] = (
    Coord("log10_mult_form", math.log10(0.1), math.log10(3.0), _POOL, "AP formation x Li 2011 (FM9 low end x0.1)", True),
    Coord("ea_form", 86.0, 136.0, _POOL, "Zhang 2026 ~101 to Li 2011 121, +-15"),
    Coord("s_form", 0.0, 1.0, _POOL, "D-2 formation pH slope (unsourced)", True),
    Coord("log10_mult_dec", math.log10(0.3), math.log10(3.0), _POOL, "AP decay x Deng 2021", True),
    Coord("ea_dec", 34.0, 107.0, _POOL, "Luo 2023 34 - 56 (first reader) to Martins 107"),
    Coord("rho", 0.54, 1.0, _POOL, "Cys returned per decay (Luo 2023 low end, first reader)"),
    Coord("phi_ref", 0.2, 0.6, _POOL, "D-3 1-DG share of AP decay at pH 5.5 (unsourced)", True),
    Coord("log10_kag_rel", -3.0, 0.0, ("R1", "R5", "CLAIM", "R2"), "AP -> GOS, dex below the FM6 ceiling", True),
    Coord("log10_mult_dg1", math.log10(0.3), math.log10(3.0), _NUIS, "1-DG other exits x shipped r_odg_aa (D-5)", True),
    Coord("s_dg1", 0.0, 0.5, _NUIS, "D-5 1-DG exits' pH slope (unsourced)", True),
    Coord("m_o2", 1.0, 4.4, _NUIS, "O2 drain on 1-DG at air (Gobert 3.4, Glomb 4.4 first reader)", True),
    Coord("log10_kcd", -3.0, -2.0, _NUIS, "D-8 Cys -> H2S at 145 C"),
    Coord("y_h2s", 0.1, 0.5, _NUIS, "D-8 H2S yield"),
    Coord("log10_khl", -2.0, -1.0, _NUIS, "D-8 H2S loss"),
    Coord("log10_kdpo", math.log10(0.3), math.log10(30.0), _NUIS, "D-7 1-DPO loss at 145 C"),
    Coord("a_nf", 0.01, 0.3, _NUIS, "D-7 NF share of 1-DPO"),
    Coord("log10_b", math.log10(0.005), math.log10(0.1), _NUIS, "D-7 NF -> MFT share"),
    Coord("s_nf", 0.07, 0.999, _NUIS, "M2's factor: MFT through free NF (Cerny pot D to E6)", True),
    Coord("log10_knf", -2.0, 0.0, _NUIS, "D-7 NF loss at 145 C"),
    Coord("log10_kmft", -3.0, math.log10(0.05), _NUIS, "D-7 MFT loss"),
    Coord("log10_ktdg", -1.0, 1.0, _NUIS, "D-10 3-DG loss at 145 C"),
    Coord("log10_odds_odg", math.log10(5e-4 / (1 - 5e-4)), 0.0, ("R1", "CLAIM", "R4P"), "r_odg_dpo share at 100 C, pH 5 (log odds; 0 - ceiling)", True),
    Coord("ea_odg", 40.0, 135.0, ("R1", "CLAIM", "R4P"), "r_odg_dpo barrier"),
    Coord("log10_u5", -3.0, 0.0, ("R5", "CLAIM"), "R5 beta-cleavage / Glomb ceiling (0 - ceiling)", True),
    Coord("ea_b5", 40.0, 135.0, ("R5", "CLAIM"), "R5 beta-cleavage barrier"),
    Coord("f_ref", 0.001, 0.999, ("CLAIM",), "FM4 R5 share of 1-DPO at 145 C, pH 5, air"),
    Coord("dea", -60.0, 60.0, ("CLAIM",), "FM4 Ea(R5) - Ea(R1)"),
    Coord("log10_kstr", math.log10(3e-5), math.log10(3e-2), ("R2",), "M3's Strecker class (0 - class constant)", True),
    Coord("ea_str", 65.0, 85.0, ("R2",), "M3's exits"),
    Coord("s_str", 0.0, 0.5, ("R2",), "D-9 Strecker pH slope (unsourced)"),
    Coord("log10_kddp", -2.0, 0.0, ("R2",), "D-9 DDP loss"),
    Coord("log10_bddp", math.log10(0.01), math.log10(0.3), ("R2",), "M2's DDP share"),
    Coord("log10_cat", 0.0, 1.0, ("R4P",), "R4' Cys-catalysis factor 1 - 10 (no source)", True),
)
COORD = {c.name: c for c in COORDS}


@dataclass(frozen=True)
class Form:
    """One setting of the declared choices. Layer fields change the predictions (new draw sets, truth AND scoring
    pinned: "the verdict must hold if the world is at that end"); data fields change only the noise model and the
    synthetic data; truth fields change only the generator."""
    name: str = "primary"
    pins: Tuple[Tuple[str, float], ...] = ()
    bands: Tuple[Tuple[str, float, float], ...] = ()
    kag_ph: bool = False
    o2_budget: bool = False
    n2_residual: float = 0.0
    # data / likelihood
    omega: float = OMEGA_DEX
    lod_scale: float = 1.0
    label_limit: bool = False
    d3_collision: bool = False
    opd: str = "post"
    # truth only
    truth_e: float = 0.0

    @property
    def layer_key(self) -> Tuple:
        return (self.pins, self.bands, self.kag_ph, self.o2_budget, self.n2_residual)


PRIMARY = Form()


def coords(rival: str, form: Form = PRIMARY) -> List[Coord]:
    bands = {n: (lo, hi) for n, lo, hi in form.bands}
    pins = dict(form.pins)
    out = []
    for c in COORDS:
        if rival not in c.rivals:
            continue
        lo, hi = bands.get(c.name, (c.lo, c.hi))
        if c.name in pins:
            lo = hi = float(pins[c.name])
        out.append(replace(c, lo=lo, hi=hi))
    return out


def form_touches(form: Form, rival: str) -> bool:
    """Does a layer form change this rival's predictions?"""
    names = {c.name for c in COORDS if rival in c.rivals}
    if any(n in names for n, _ in form.pins) or any(n in names for n, _, _ in form.bands):
        return True
    if rival in ("R1", "R5", "CLAIM", "R2") and form.kag_ph:
        return True
    if form.o2_budget or form.n2_residual:
        return rival != "R0"
    return False


# ================================================================================================== the pots

@dataclass(frozen=True)
class Pot:
    name: str
    arm: str
    t_c: float
    ph: float
    atm: str
    g0: float = 0.0
    cys0: float = 0.0
    ap0: float = 0.0
    dg10: float = 0.0
    gos0: float = 0.0
    h2s0: float = 0.0
    gly: float = 0.0
    t_end: float = 20.0
    times: Tuple[float, ...] = (20.0,)

    @property
    def t_k(self) -> float:
        return self.t_c + 273.15


def _arm1() -> List[Pot]:
    pots = []
    for t_c, phs in ((145.0, (3.0, 5.0, 7.0)), (120.0, (5.0,))):
        for ph in phs:
            for atm in ("air", "N2"):
                pots.append(Pot(f"a1_{int(t_c)}_pH{int(ph)}_{atm}", "arm1", t_c, ph, atm, g0=100.0, cys0=33.0))
    return pots


def _arm2(variant: str) -> List[Pot]:
    pots = []
    atms = ("N2", "air") if variant == "a2air" else ("N2",)
    for dic in ("DG1", "GOS"):
        for ph in (5.0, 7.0):
            for gly in (0.0, 33.0):
                for atm in atms:
                    if variant == "a2air" and atm == "N2":
                        continue
                    kw = dict(dg10=2.0) if dic == "DG1" else dict(gos0=2.0)
                    if variant == "a2cys":
                        kw["cys0"] = 6.0
                    else:
                        kw["h2s0"] = 6.0
                    pots.append(Pot(f"{variant}_{dic}_pH{int(ph)}_gly{int(gly)}_{atm}", variant, 145.0, ph, atm,
                                    gly=gly, times=(5.0, 20.0), **kw))
    return pots


def _arm3(variant: str) -> List[Pot]:
    pots = []
    for t_c in (100.0, 145.0):
        for ph in (3.0, 5.0, 7.0):
            for atm in ("air", "N2"):
                cys_fed = 10.0 if variant == "a3m" else 0.0
                pots.append(Pot(f"{variant}_fed_{int(t_c)}_pH{int(ph)}_{atm}", variant, t_c, ph, atm, ap0=10.0,
                                cys0=cys_fed, t_end=60.0, times=(20.0, 60.0)))
                pots.append(Pot(f"{variant}_sug_{int(t_c)}_pH{int(ph)}_{atm}", variant, t_c, ph, atm, g0=30.0,
                                cys0=10.0, t_end=60.0, times=(20.0, 60.0)))
    return pots


ARM_VARIANTS = ("a1", "a2", "a2air", "a2cys", "a3", "a3m")
POTS: Dict[str, List[Pot]] = {"a1": _arm1(), "a2": _arm2("a2"), "a2air": _arm2("a2air"), "a2cys": _arm2("a2cys"),
                              "a3": _arm3("a3"), "a3m": _arm3("a3m")}
REF_POT = Pot("ref_145_pH5_air", "ref", 145.0, 5.0, "air", g0=100.0, cys0=33.0)


# ============================================================================================ draws

def draw(rival: str, form: Form, n: int, rng: np.random.Generator, sobol: bool = False,
         proposal: bool = False) -> Dict[str, np.ndarray]:
    """Band-prior draws (uniform on every coordinate). With ``proposal`` the claim's f_ref is drawn logit-uniform
    instead (mass at the R1 and R5 edges, where label shares pin it to ~1 %), and ``_logw`` carries log p / q for the
    self-normalised density (D-15). Truth draws always come from the prior."""
    cs = coords(rival, form)
    if sobol:
        from scipy.stats import qmc
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = qmc.Sobol(d=max(len(cs), 1), scramble=True, seed=rng).random(n)
    else:
        u = rng.random((n, max(len(cs), 1)))
    th = {c.name: c.lo + (c.hi - c.lo) * u[:, i] for i, c in enumerate(cs)} | {"_n": np.zeros(n), "_logw": np.zeros(n)}
    if proposal and rival == "CLAIM":
        c = next(c for c in cs if c.name == "f_ref")
        if c.hi > c.lo:
            lo, hi = math.log(c.lo / (1 - c.lo)), math.log(c.hi / (1 - c.hi))
            x = lo + (hi - lo) * u[:, [cc.name for cc in cs].index("f_ref")]
            f = 1.0 / (1.0 + np.exp(-x))
            th["f_ref"] = f
            th["_logw"] = np.log((hi - lo) * f * (1.0 - f) / (c.hi - c.lo))
    return th


def n_of(th: Mapping[str, np.ndarray]) -> int:
    return len(th["_n"])


def take(th: Mapping[str, np.ndarray], idx) -> Dict[str, np.ndarray]:
    return {k: np.asarray(v)[idx] for k, v in th.items()}


def concat(a: Mapping[str, np.ndarray], b: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    return {k: np.concatenate([a[k], b[k]]) for k in a}


# ============================================================================================ rates

def o2_factor(pot: Pot, form: Form) -> float:
    if pot.atm == "N2":
        return float(form.n2_residual)
    if form.o2_budget:
        return max(0.0, 1.0 - pot.cys0 / (4.0 * O2_MM_HALF_FILL))
    return 1.0


def rates(th: Mapping[str, np.ndarray], rival: str, form: Form, pot: Pot) -> Dict[str, np.ndarray]:
    n = n_of(th)
    z = np.zeros(n)
    one = np.ones(n)
    t, ph = pot.t_k, pot.ph
    o2 = o2_factor(pot, form)
    k: Dict[str, np.ndarray] = {}
    # nuisance (every rival)
    k["kcd"] = 10.0 ** th["log10_kcd"] * arr(EA_CD, t, T145)
    k["y_h2s"] = th["y_h2s"]
    k["khl"] = 10.0 ** th["log10_khl"] * arr(EA_HL, t, T145)
    k["kdg1o"] = K_ODG_AA_100 * 10.0 ** th["log10_mult_dg1"] * arr(EA_ODG_AA, t, T100) * 10.0 ** (th["s_dg1"] * (ph - 5.0))
    k["kdg1o2"] = k["kdg1o"] * (th["m_o2"] - 1.0) * o2
    k["kdpo"] = 10.0 ** th["log10_kdpo"] * arr(EA_DOWN, t, T145)
    k["a_nf"] = th["a_nf"]
    b = 10.0 ** th["log10_b"]
    k["b"] = b
    s = th["s_nf"]
    k["a_dir"] = np.minimum(th["a_nf"] * b * (1.0 - s) / s, 1.0 - th["a_nf"])
    k["knf"] = 10.0 ** th["log10_knf"] * arr(EA_DOWN, t, T145)
    k["kmft"] = 10.0 ** th["log10_kmft"] * arr(EA_DOWN, t, T145)
    k["ktdg"] = 10.0 ** th["log10_ktdg"] * arr(EA_DOWN, t, T145)
    k["kgo"] = 10.0 ** LOG10_K_G_GO_100 * arr(EA_G_GO, t, T100) * one
    k["kamf"] = K_FRU_ODG_100 * arr(EA_FRU_ODG, t, T100) * (10.0 ** th["log10_cat"] if rival == "R4P" else one)
    k["kgtdg"] = K_GLC_TDG_100 * arr(EA_GLC_TDG, t, T100) * one
    # R0's C2 + C3 bound (sink-free, largest fed-pot yields): R0 only. In the entry rivals the shipped route is the
    # engine's own ~1e-14 mM (probe C), ~8 dex below the MFT LOD; carrying R0's BOUND there would plant a spurious
    # mixed share in every low-entry draw, so it is 0 (D-14).
    k["kha"] = K_GLC_HA_145 * arr(EA_GLC_HA, t, T145) * (one if rival == "R0" else z)
    k["kc6"] = K_C6DD_145 * arr(EA_DOWN, t, T145) * one
    # the pool
    if rival in POOL_RIVALS:
        k["kf"] = K_FORM_LI_100_PH7 * 10.0 ** th["log10_mult_form"] * arr(th["ea_form"], t, T100) * \
            10.0 ** (th["s_form"] * (ph - 7.0))
        k["kdec"] = K_DEC_DENG_100 * 10.0 ** th["log10_mult_dec"] * arr(th["ea_dec"], t, T100) * \
            10.0 ** (S_DEC_DENG * (ph - 5.5))
        k["rho"] = th["rho"]
        odds = th["phi_ref"] / (1.0 - th["phi_ref"]) * 10.0 ** (S_SPLIT_ODDS * (ph - 5.5))
        k["phi"] = odds / (1.0 + odds)
    else:
        k["kf"], k["kdec"], k["rho"], k["phi"] = z, z, one, 0.5 * one
    if rival in ("R1", "R5", "CLAIM", "R2"):
        kag = 10.0 ** (LOG10_K_AMA_G_100 + KAG_CEILING_DEX + th["log10_kag_rel"]) * arr(EA_AMA_G, t, T100)
        if form.kag_ph:
            kag = kag * 10.0 ** (S_SPLIT_ODDS * (ph - 6.7))
        k["kag"] = kag * o2
    else:
        k["kag"] = z
    k["ap_sink"] = one if rival == "R4P" else z
    # entries
    if rival in ("R1", "CLAIM", "R4P"):
        k["kodg"] = 10.0 ** (th["log10_odds_odg"] + th.get("_dlog_odg", 0.0)) * K_ODG_AA_100 * \
            10.0 ** th["log10_mult_dg1"] * arr(th["ea_odg"], t, T100)
    else:
        k["kodg"] = z
    if rival in ("R5", "CLAIM"):
        ea5 = th["ea_b5"] if rival == "R5" else th["_ea_b5"]
        k["kb5"] = 10.0 ** (th["log10_u5"] + th.get("_dlog_b5", 0.0)) * K_GOS_GLOMB_50 * arr(ea5, t, T50)
    else:
        k["kb5"] = z
    if rival == "R2":
        k["kstr"] = 10.0 ** th["log10_kstr"] * arr(th["ea_str"], t, T100) * 10.0 ** (th["s_str"] * (ph - 5.0))
        k["kddp"] = 10.0 ** th["log10_kddp"] * arr(EA_DOWN, t, T145)
        k["bddp"] = 10.0 ** th["log10_bddp"]
    else:
        k["kstr"], k["kddp"], k["bddp"] = z, one, z
    return k


# ============================================================================================ the integrator

STATE = ("G", "CYS", "AP", "DG1", "GOS", "TDG", "DPOK", "DPOL", "DPO3", "NFK", "NFL", "MFTK", "MFTL", "C6DD", "DDP",
         "MFTD", "H2S", "FA")
HEX = ("G", "AP", "DG1", "GOS", "TDG", "C6DD")
C5_K = ("DPOK", "NFK", "MFTK", "DDP", "MFTD")
C5_L = ("DPOL", "DPO3", "NFL", "MFTL")
LEDGER = ("S6", "S5K", "S1K", "S5L", "S1L", "CYS_BOUND", "CYS_DEG", "H2S_LOSS", "INT_G")


def simulate(th: Mapping[str, np.ndarray], rival: str, form: Form, pot: Pot, dt: float = DT,
             full: bool = False) -> Dict[str, Any]:
    """Integrate one pot for every draw. Returns {(species, t): array} for the reads, plus the final state and the
    ledger when ``full``."""
    k = rates(th, rival, form, pot)
    n = n_of(th)
    x = {s: np.zeros(n) for s in STATE}
    x["G"][:] = pot.g0
    x["CYS"][:] = pot.cys0
    x["AP"][:] = pot.ap0
    x["DG1"][:] = pot.dg10
    x["GOS"][:] = pot.gos0
    x["H2S"][:] = pot.h2s0
    led = {s: np.zeros(n) for s in LEDGER}
    steps = int(round(pot.t_end / dt))
    read_steps = {int(round(t / dt)): t for t in pot.times}
    out: Dict[Any, np.ndarray] = {}
    kf, kdec, kag, rho, phi = k["kf"], k["kdec"], k["kag"], k["rho"], k["phi"]
    ka = kdec + kag
    kret = rho * kdec + kag
    sink = k["ap_sink"]
    kd_tot = k["kdg1o"] + k["kdg1o2"] + k["kodg"]
    gly = pot.gly
    for step in range(1, steps + 1):
        g_old = x["G"]
        c0, a0 = x["CYS"], x["AP"]
        a11 = 1.0 + dt * (kf * g_old + k["kcd"])
        det = a11 * (1.0 + dt * ka) - (dt * kret) * (dt * kf * g_old)
        c = (c0 * (1.0 + dt * ka) + dt * kret * a0) / det
        a = (a11 * a0 + dt * kf * g_old * c0) / det
        f_form = dt * kf * g_old * c
        f_dec = dt * kdec * a
        f_ag = dt * kag * a
        f_cd = dt * k["kcd"] * c
        led["CYS_BOUND"] += f_dec * (1.0 - rho)
        led["CYS_DEG"] += f_cd
        f_amf = dt * k["kamf"] * g_old
        f_gtdg = dt * k["kgtdg"] * g_old
        led["INT_G"] += dt * g_old
        x["G"] = g_old - f_form - f_amf - f_gtdg
        x["CYS"], x["AP"] = c, a
        h2s = (x["H2S"] + k["y_h2s"] * f_cd) / (1.0 + dt * k["khl"])
        led["H2S_LOSS"] += dt * k["khl"] * h2s
        x["H2S"] = h2s
        h = h2s / (K_H2S_HALF + h2s)
        # 1-DG
        kstr_eff = k["kstr"] * (c + gly)
        dg1 = (x["DG1"] + phi * f_dec * (1.0 - sink) + f_amf) / (1.0 + dt * (kd_tot + kstr_eff))
        f_odg = dt * k["kodg"] * dg1
        f_str = dt * kstr_eff * dg1
        f_dg1o = dt * (k["kdg1o"] + k["kdg1o2"]) * dg1
        x["DG1"] = dg1
        led["S6"] += f_dg1o + sink * f_dec
        led["S1K"] += f_odg
        x["FA"] = x["FA"] + Y_FA_DG1 * f_dg1o
        # glucosone
        gos = (x["GOS"] + f_ag) / (1.0 + dt * (k["kgo"] + k["kb5"]))
        f_b5 = dt * k["kb5"] * gos
        led["S6"] += dt * k["kgo"] * gos
        led["S1L"] += f_b5
        x["FA"] = x["FA"] + f_b5
        x["GOS"] = gos
        # 3-DG
        tdg = (x["TDG"] + (1.0 - phi) * f_dec * (1.0 - sink) + f_gtdg) / (1.0 + dt * k["ktdg"])
        led["S6"] += dt * k["ktdg"] * tdg
        x["TDG"] = tdg
        # 1-DPO, 3-DPO
        kdpo = k["kdpo"]
        dpok = (x["DPOK"] + f_odg) / (1.0 + dt * kdpo)
        dpol = (x["DPOL"] + phi * f_b5) / (1.0 + dt * kdpo)
        dpo3 = (x["DPO3"] + (1.0 - phi) * f_b5) / (1.0 + dt * kdpo)
        fk, fl, f3 = dt * kdpo * dpok, dt * kdpo * dpol, dt * kdpo * dpo3
        x["DPOK"], x["DPOL"], x["DPO3"] = dpok, dpol, dpo3
        a_nf, a_dir = k["a_nf"], k["a_dir"] * h
        led["S5K"] += fk * (1.0 - a_nf - a_dir)
        led["S5L"] += fl * (1.0 - a_nf - a_dir) + f3
        # NF
        nfk = (x["NFK"] + a_nf * fk) / (1.0 + dt * k["knf"])
        nfl = (x["NFL"] + a_nf * fl) / (1.0 + dt * k["knf"])
        gk, gl = dt * k["knf"] * nfk, dt * k["knf"] * nfl
        bh = k["b"] * h
        led["S5K"] += gk * (1.0 - bh)
        led["S5L"] += gl * (1.0 - bh)
        x["NFK"], x["NFL"] = nfk, nfl
        # MFT
        mftk = (x["MFTK"] + a_dir * fk + bh * gk) / (1.0 + dt * k["kmft"])
        mftl = (x["MFTL"] + a_dir * fl + bh * gl) / (1.0 + dt * k["kmft"])
        led["S5K"] += dt * k["kmft"] * mftk
        led["S5L"] += dt * k["kmft"] * mftl
        x["MFTK"], x["MFTL"] = mftk, mftl
        # R2
        c6 = (x["C6DD"] + f_str) / (1.0 + dt * k["kc6"])
        f_c6 = dt * k["kc6"] * c6
        led["S1K"] += f_c6
        ddp = (x["DDP"] + f_c6) / (1.0 + dt * k["kddp"])
        f_ddp = dt * k["kddp"] * ddp
        bd = k["bddp"] * h
        led["S5K"] += f_ddp * (1.0 - bd)
        mftd = (x["MFTD"] + bd * f_ddp) / (1.0 + dt * k["kmft"])
        led["S5K"] += dt * k["kmft"] * mftd
        x["C6DD"], x["DDP"], x["MFTD"] = c6, ddp, mftd
        if step in read_steps:
            t = read_steps[step]
            c23 = k["kha"] * led["INT_G"]
            reads = {"MFTK": x["MFTK"], "MFTL": x["MFTL"], "MFTD": x["MFTD"], "MFTC": c23 * Y_MP_T7 * Y_MFT_T8,
                     "NFK": x["NFK"], "NFL": x["NFL"], "NFC": c23 * Y_NF_T6, "DPO1": x["DPOK"] + x["DPOL"],
                     "DPO3": x["DPO3"], "DDQ": x["C6DD"] + x["DDP"], "FA": x["FA"], "AP": x["AP"], "CYS": x["CYS"],
                     "GOS": x["GOS"], "DG1": x["DG1"], "TDG": x["TDG"]}
            for sp, v in reads.items():
                out[(sp, t)] = v.copy()
    if full:
        out["_state"] = {s: v.copy() for s, v in x.items()}
        out["_ledger"] = {s: v.copy() for s, v in led.items()}
        out["_k"] = k
    return out


def balances(pot: Pot, res: Mapping[str, Any]) -> Dict[str, float]:
    """Max absolute residuals (mM, or mM carbon) of the hexose-unit, carbon, Cys, sulfur and label-position books."""
    x, led = res["_state"], res["_ledger"]
    units0 = pot.g0 + pot.ap0 + pot.dg10 + pot.gos0
    hexp = sum(x[s] for s in HEX)
    c5 = sum(x[s] for s in C5_K + C5_L)
    units = hexp + c5 + led["S6"] + led["S5K"] + led["S5L"]
    carbon = 6.0 * hexp + 5.0 * c5 + 6.0 * led["S6"] + 5.0 * (led["S5K"] + led["S5L"]) + led["S1K"] + led["S1L"]
    cys = x["CYS"] + x["AP"] + led["CYS_BOUND"] + led["CYS_DEG"]
    sul = x["H2S"] + led["H2S_LOSS"]
    sul0 = pot.h2s0 + res["_k"]["y_h2s"] * led["CYS_DEG"]
    # label positions at e = 0 and e = 1: every glucose position is held once per hexose unit
    c5k = sum(x[s] for s in C5_K) + led["S5K"]
    c5l = sum(x[s] for s in C5_L) + led["S5L"]
    pos_err = 0.0
    for e in (0, 1):
        kept = {1, 2, 3, 4, 5} if e == 0 else {2, 3, 4, 5, 6}
        lost = {2, 3, 4, 5, 6} if e == 0 else {1, 2, 3, 4, 5}
        for p in range(1, 7):
            tot = hexp + led["S6"] + c5k * (p in kept) + led["S1K"] * (p not in kept) + \
                c5l * (p in lost) + led["S1L"] * (p not in lost)
            pos_err = max(pos_err, float(np.max(np.abs(tot - units0))))
    return {"hexose_units": float(np.max(np.abs(units - units0))),
            "carbon": float(np.max(np.abs(carbon - 6.0 * units0))),
            "cysteine": float(np.max(np.abs(cys - pot.cys0 - pot.ap0))),
            "sulfide": float(np.max(np.abs(sul - sul0))),
            "label_positions": pos_err}


# ============================================================================================ label maps

KEPT_U, KEPT_S = (1, 2, 3, 4, 5), (6, 5, 4, 3, 2)
LOST_U, LOST_S = (2, 3, 4, 5, 6), (5, 4, 3, 2, 1)
C3_PARTS, C2_PARTS = ((1, 2, 3), (6, 5, 4)), ((2, 1), (5, 6))
FRAG_IDX = {"M": (0, 1, 2, 3, 4), "ring": (1, 2, 3, 4), "acyl": (0, 1)}


def _intact_prob(m: Tuple[int, ...], read: str) -> float:
    """P(read) for an intact C5 with glucose positions ``m`` (methyl first). Reads: 'L{k}_{frag}' (single label at k,
    labelled fraction of the fragment), 'L12_M1' / 'L12_M2', 'CAM' (mixed share)."""
    if read == "CAM":
        return 0.0
    if read.startswith("L12_"):
        cnt = sum(1 for p in m if p in (1, 2))
        return float(cnt == int(read[-1]))
    lab, frag = read.split("_")
    kpos = int(lab[1:])
    return float(kpos in [m[i] for i in FRAG_IDX[frag]])


def _c23_prob(read: str) -> float:
    """The same reads for R0's C2 + C3 product (C3 part + C2 part from two molecules, each choice p 1/2)."""
    if read == "CAM":
        return 0.5                                     # 2 q (1 - q), q = 1/2: the binomial mixed share
    tot = 0.0
    for c3 in C3_PARTS:
        for c2 in C2_PARTS:
            m = c3 + c2
            if read.startswith("L12_"):
                cnt = sum(1 for p in m if p in (1, 2))
                tot += 0.25 * float(cnt == int(read[-1]))
            else:
                lab, frag = read.split("_")
                tot += 0.25 * float(int(lab[1:]) in [m[i] for i in FRAG_IDX[frag]])
    return tot


def label_table(read: str) -> Dict[str, float]:
    """Hand-checkable per-family probabilities: K / L unswapped (u) and swapped (s), and C2 + C3."""
    return {"K_u": _intact_prob(KEPT_U, read), "K_s": _intact_prob(KEPT_S, read),
            "L_u": _intact_prob(LOST_U, read), "L_s": _intact_prob(LOST_S, read), "C23": _c23_prob(read)}


def share_coeffs(a_k: np.ndarray, a_l: np.ndarray, a_c: np.ndarray, read: str) -> Tuple[np.ndarray, np.ndarray]:
    """share(e) = A + B e for the pool (K family amount a_k, L family a_l, C2 + C3 a_c)."""
    tb = label_table(read)
    tot = a_k + a_l + a_c
    tot = np.where(tot > 0, tot, 1e-300)
    a = (a_k * tb["K_u"] + a_l * tb["L_u"] + a_c * tb["C23"]) / tot
    b = (a_k * (tb["K_s"] - tb["K_u"]) + a_l * (tb["L_s"] - tb["L_u"])) / tot
    return a, b


# ============================================================================================ the claim's mixture (FM4)

def ref_share(th: Mapping[str, np.ndarray], form: Form, dt: float = DT) -> np.ndarray:
    """The R5 share of 1-DPO formed in the reference cell, from the ledger (S1L = beta-cleavage events)."""
    res = simulate(th, "CLAIM", form, REF_POT, dt=dt, full=True)
    led, k = res["_ledger"], res["_k"]
    r5 = k["phi"] * led["S1L"]
    r1 = led["S1K"] - 0.0                                 # no R2 step in the claim: S1K is r_odg_dpo only
    return r5 / np.maximum(r5 + r1, 1e-300)


def claim_rescale(th: Dict[str, np.ndarray], form: Form, iters: int = 40) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Rescale the two entries symmetrically (log10 shift +mu/2 on R5, -mu/2 on R1) so the reference share equals
    f_ref; Ea(R5) = clip(Ea(R1) + dEa, 40, 135). Returns the draws and a feasibility mask (both rescaled entries
    inside their bands). When a form pins an entry, f emerges instead (no rescale)."""
    th = dict(th)
    th["_ea_b5"] = np.clip(th["ea_odg"] + th["dea"], COORD["ea_b5"].lo, COORD["ea_b5"].hi)
    pins = dict(form.pins)
    n = n_of(th)
    if "log10_odds_odg" in pins or "log10_u5" in pins:
        th["_dlog_odg"] = np.zeros(n)
        th["_dlog_b5"] = np.zeros(n)
        return th, np.ones(n, bool)
    lo, hi = np.full(n, -14.0), np.full(n, 14.0)
    target = th["f_ref"]
    for _ in range(iters):
        mu = 0.5 * (lo + hi)
        th["_dlog_odg"], th["_dlog_b5"] = -mu / 2.0, mu / 2.0
        f = ref_share(th, form)
        up = f < target
        lo, hi = np.where(up, mu, lo), np.where(up, hi, mu)
    mu = 0.5 * (lo + hi)
    th["_dlog_odg"], th["_dlog_b5"] = -mu / 2.0, mu / 2.0
    band = {c.name: (c.lo, c.hi) for c in coords("CLAIM", form)}
    eo = th["log10_odds_odg"] + th["_dlog_odg"]
    eu = th["log10_u5"] + th["_dlog_b5"]
    ok = (eo >= band["log10_odds_odg"][0] - 1e-9) & (eo <= band["log10_odds_odg"][1] + 1e-9) & \
         (eu >= band["log10_u5"][0] - 1e-9) & (eu <= band["log10_u5"][1] + 1e-9)
    return th, ok


def draw_set(rival: str, form: Form, n: int, rng: np.random.Generator, sobol: bool = False,
             proposal: bool = False, max_rounds: int = 40) -> Tuple[Dict[str, np.ndarray], Dict[str, int]]:
    """n feasible draws (the claim's rescale can reject); returns the draws and the counts. The density of a draw
    set is self-normalised over its feasible draws (``lppd``), so the target is the band prior restricted to the
    feasible set whichever proposal drew it."""
    if rival != "CLAIM":
        return draw(rival, form, n, rng, sobol), {"feasible": n, "drawn": n}
    got, drawn = None, 0
    for _ in range(max_rounds):
        m = max(2 * n, 256)
        th = draw(rival, form, m, rng, sobol=False, proposal=proposal)
        th, ok = claim_rescale(th, form)
        drawn += m
        sel = take(th, np.where(ok)[0])
        got = sel if got is None else concat(got, sel)
        if n_of(got) >= n:
            break
    return take(got, np.arange(min(n, n_of(got)))), {"feasible": int(min(n, n_of(got))), "drawn": drawn,
                                                     "accepted_total": int(n_of(got))}


# ============================================================================================ observables

@dataclass(frozen=True)
class Obs:
    group: str          # e.g. 'a1:MFT' (level, species) or 'a1:L1_ring' (share)
    kind: str           # 'level' | 'share'
    species: str        # MFT / NF / ... (levels: the tau offset's key; shares: whose pool)
    pot: str
    t: float
    read: str = ""      # share read
    t_c: float = 145.0


LEVEL_READS = {
    "a1": (("MFT", 20.0), ("NF", 20.0)),
    "a2": (("MFT", 20.0), ("NF", 20.0), ("FA", 20.0), ("DPO1", 5.0), ("DPO3", 5.0), ("DDQ", 5.0)),
    "a3": (("MFT", 60.0), ("AP", 60.0), ("CYS", 60.0), ("GOS", 60.0)),
}
SHARE_READS_145 = ("L1_M", "L1_ring", "L1_acyl", "L6_M", "L6_ring", "L6_acyl", "L2_M", "L2_ring", "L2_acyl",
                   "L12_M1", "L12_M2", "CAM")
SHARE_READS_120 = ("L1_M", "L1_ring", "L1_acyl", "L2_M", "L2_ring", "L2_acyl")
NF_SHARE_READS = ("CAM", "L1_M")


def layout() -> List[Obs]:
    obs: List[Obs] = []
    for av in ARM_VARIANTS:
        base = av[:2]
        for pot in POTS[av]:
            for sp, t in LEVEL_READS[base]:
                if base == "a3":
                    fed = "_fed_" in pot.name
                    if sp == "AP" and not fed or sp == "CYS" and fed or sp == "GOS" and pot.t_c != 145.0:
                        continue
                obs.append(Obs(f"{av}:{sp}", "level", sp, pot.name, t, t_c=pot.t_c))
    for pot in POTS["a1"]:
        reads = SHARE_READS_145 if pot.t_c == 145.0 else SHARE_READS_120
        sfx = "" if pot.t_c == 145.0 else "_120"
        for r in reads:
            obs.append(Obs(f"a1:{r}{sfx}", "share", "MFT", pot.name, 20.0, r, pot.t_c))
        if pot.t_c == 145.0:
            for r in NF_SHARE_READS:
                obs.append(Obs(f"a1:NF_{r}", "share", "NF", pot.name, 20.0, r, pot.t_c))
    return obs


OBS = layout()
LEVEL_IDX = [i for i, o in enumerate(OBS) if o.kind == "level"]
SHARE_IDX = [i for i, o in enumerate(OBS) if o.kind == "share"]
GROUPS: Dict[str, List[int]] = {}
for _i, _o in enumerate(OBS):
    GROUPS.setdefault(_o.group, []).append(_i)
POT_BY_NAME = {p.name: p for av in ARM_VARIANTS for p in POTS[av]}


@dataclass
class Prediction:
    """Every draw's predicted observables: levels (mM, linear), shares as A + B e, the pool amounts behind each share
    (for readability) and AP at the read time (the OPD artefact)."""
    level: np.ndarray       # [N, n_obs] (share columns unused)
    share_a: np.ndarray
    share_b: np.ndarray
    ap_at: np.ndarray       # AP at each obs's pot and time (for the OPD artefact on GOS)
    n: int
    logw: Optional[np.ndarray] = None       # log prior / proposal (D-15); 0 for prior draws
    loglevel: Optional[np.ndarray] = None   # log10 of level, floored at 1e-30

    def __post_init__(self):
        if self.logw is None:
            self.logw = np.zeros(self.n)
        if self.loglevel is None:
            with np.errstate(invalid="ignore"):
                self.loglevel = np.log10(np.maximum(self.level, 1e-30))


def predict(th: Mapping[str, np.ndarray], rival: str, form: Form, arms: Sequence[str] = ARM_VARIANTS,
            dt: float = DT) -> Prediction:
    n = n_of(th)
    nobs = len(OBS)
    level = np.full((n, nobs), np.nan)
    sa = np.full((n, nobs), np.nan)
    sb = np.zeros((n, nobs))
    ap_at = np.zeros((n, nobs))
    sims: Dict[str, Dict[Any, np.ndarray]] = {}
    for av in arms:
        for pot in POTS[av]:
            sims[pot.name] = simulate(th, rival, form, pot, dt=dt)
    for i, o in enumerate(OBS):
        r = sims.get(o.pot)
        if r is None:
            continue
        g = lambda sp: r[(sp, o.t)]  # noqa: E731
        if o.kind == "level":
            if o.species == "MFT":
                v = g("MFTK") + g("MFTL") + g("MFTD") + g("MFTC")
            elif o.species == "NF":
                v = g("NFK") + g("NFL") + g("NFC")
            else:
                v = g(o.species)
            level[:, i] = v
            ap_at[:, i] = g("AP")
        else:
            if o.species == "MFT":
                ak, al, ac = g("MFTK") + g("MFTD"), g("MFTL"), g("MFTC")
            else:
                ak, al, ac = g("NFK"), g("NFL"), g("NFC")
            sa[:, i], sb[:, i] = share_coeffs(ak, al, ac, o.read)
            level[:, i] = ak + al + ac
    return Prediction(level, sa, sb, ap_at, n, logw=np.asarray(th.get("_logw", np.zeros(n)), float))


def take_pred(p: Prediction, idx) -> Prediction:
    idx = np.atleast_1d(idx)
    return Prediction(p.level[idx], p.share_a[idx], p.share_b[idx], p.ap_at[idx], len(idx), logw=p.logw[idx],
                      loglevel=p.loglevel[idx])


# ============================================================================================ measurement model

def lod_of(o: Obs, form: Form) -> float:
    return LOD_MM[o.species] * form.lod_scale


def level_sd(o: Obs, n_heat: int, form: Form) -> float:
    return math.sqrt(form.omega ** 2 + CV_DEX ** 2 / n_heat)


def share_sd(o: Obs, n_heat: int, form: Form) -> float:
    base = SHARE_SD if form.label_limit else SHARE_SD / math.sqrt(n_heat)
    sd2 = base ** 2 + SHARE_FLOOR ** 2
    if form.d3_collision and o.read == "CAM" and o.species == "MFT":
        sd2 += D3_EXTRA_SD ** 2 / n_heat
    return math.sqrt(sd2)


def opd_extra_dex(o: Obs, pred_level: np.ndarray, ap_at: np.ndarray, n_heat: int, form: Form) -> np.ndarray:
    """The blank-subtracted OPD artefact's residual error on a glucosone level, in dex."""
    if o.species != "GOS":
        return np.zeros_like(pred_level)
    a = OPD_ARTEFACT[form.opd]
    rel = CV_SIDA * a * ap_at / np.maximum(pred_level, 1e-30) / math.sqrt(n_heat)
    return np.log10(1.0 + rel)


@dataclass
class Data:
    y_log: np.ndarray       # [n_obs] log10 level (NaN for shares)
    censored: np.ndarray    # [n_obs] bool (levels)
    share: np.ndarray       # [n_obs] share value (NaN where not read)
    read: np.ndarray        # [n_obs] bool: share read (levels: always True)
    n_heat: int


def synthesize(p: Prediction, j: int, n_heat: int, form: Form, rng: np.random.Generator) -> Data:
    """One synthetic laboratory from truth draw j: one tau offset per species, omega + CV per level value, the OPD
    artefact residual on glucosone, shares at the truth's e (form.truth_e) with the share error; censoring at the
    LODs; a share is read only where its unlabelled pool is >= 3 LOD."""
    nobs = len(OBS)
    y = np.full(nobs, np.nan)
    cens = np.zeros(nobs, bool)
    sh = np.full(nobs, np.nan)
    rd = np.ones(nobs, bool)
    off = {sp: rng.normal(0.0, TAU_ABS_DEX) for sp in LOD_MM}
    for i in LEVEL_IDX:
        o = OBS[i]
        v = p.level[j, i]
        sd = math.sqrt(level_sd(o, n_heat, form) ** 2 +
                       float(opd_extra_dex(o, np.array([v]), np.array([p.ap_at[j, i]]), n_heat, form)[0]) ** 2)
        yl = math.log10(max(v, 1e-30)) + off[o.species] + rng.normal(0.0, sd)
        lod = math.log10(lod_of(o, form))
        if yl < lod:
            cens[i] = True
            y[i] = lod
        else:
            y[i] = yl
    # shares: readable iff the same cell's unlabelled level (synthetic) is >= LOQ
    lev_of = {(OBS[i].pot, OBS[i].species): (y[i], cens[i]) for i in LEVEL_IDX if OBS[i].group.startswith("a1:")}
    for i in SHARE_IDX:
        o = OBS[i]
        yl, c = lev_of[(o.pot, o.species)]
        if c or yl < math.log10(LOQ_OVER_LOD * LOD_MM[o.species] * form.lod_scale):
            rd[i] = False
            continue
        s_true = p.share_a[j, i] + p.share_b[j, i] * form.truth_e
        sh[i] = s_true + rng.normal(0.0, share_sd(o, n_heat, form))
    return Data(y, cens, sh, rd, n_heat)


_GH_X, _GH_W = np.polynomial.hermite.hermgauss(GH_NODES)
GH_DELTA = math.sqrt(2.0) * TAU_ABS_DEX * _GH_X
GH_LOGW = np.log(_GH_W / math.sqrt(math.pi))


def group_loglik(p: Prediction, d: Data, form: Form, groups: Sequence[str]) -> Dict[str, np.ndarray]:
    """Per group: levels -> [N, GH_NODES] (log L at each tau node of the group's species, summed over the group's
    values); shares -> [N, len(E_GRID)]."""
    out: Dict[str, np.ndarray] = {}
    for gname in groups:
        idx = GROUPS[gname]
        o0 = OBS[idx[0]]
        if o0.kind == "level":
            ii = np.array(idx)
            pred = p.loglevel[:, ii]                                          # [N, V]
            sd = np.repeat(np.array([level_sd(OBS[i], d.n_heat, form) for i in ii])[None, :], p.n, 0)
            if o0.species == "GOS":
                extra = np.stack([opd_extra_dex(OBS[i], p.level[:, i], p.ap_at[:, i], d.n_heat, form) for i in ii], 1)
                sd = np.sqrt(sd ** 2 + extra ** 2)
            cen = d.censored[ii]
            tot = np.zeros((p.n, GH_NODES))
            if (~cen).any():
                o_ = ~cen
                z = (d.y_log[ii][o_][None, :, None] - pred[:, o_][:, :, None] - GH_DELTA[None, None, :]) / \
                    sd[:, o_][:, :, None]
                tot += (-0.5 * z ** 2).sum(axis=1) - np.log(sd[:, o_]).sum(axis=1)[:, None] - \
                    0.5 * math.log(2 * math.pi) * int(o_.sum())
            if cen.any():
                z = (d.y_log[ii][cen][None, :, None] - pred[:, cen][:, :, None] - GH_DELTA[None, None, :]) / \
                    sd[:, cen][:, :, None]
                tot += log_ndtr(z).sum(axis=1)
            out[gname] = tot                                                  # [N, K]
        else:
            ii = np.array([i for i in idx if d.read[i]])
            if ii.size == 0:
                out[gname] = np.zeros((p.n, len(E_GRID)))
                continue
            s = p.share_a[:, ii][:, :, None] + p.share_b[:, ii][:, :, None] * E_GRID[None, None, :]
            sd = np.array([share_sd(OBS[i], d.n_heat, form) for i in ii])[None, :, None]
            z = (d.share[ii][None, :, None] - s) / sd
            out[gname] = (-0.5 * z ** 2 - np.log(sd) - 0.5 * math.log(2 * math.pi)).sum(axis=1)
    return out


def _lse(x: np.ndarray, axis: int) -> np.ndarray:
    m = np.max(x, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    return np.squeeze(m, axis=axis) + np.log(np.sum(np.exp(x - m), axis=axis))


def design_parts(gl: Mapping[str, np.ndarray], groups: Sequence[str]) -> Dict[str, Any]:
    """A design's per-draw pieces: the share log-lik on the e grid [N, G] (or None) and, per species, the level
    log-lik at each tau node [N, K] (summed over every arm of the design: one laboratory, one offset per species)."""
    sh = [g for g in groups if OBS[GROUPS[g][0]].kind == "share"]
    share = sum(gl[g] for g in sh) if sh else None
    by_sp: Dict[str, np.ndarray] = {}
    for g in groups:
        o = OBS[GROUPS[g][0]]
        if o.kind == "level":
            by_sp[o.species] = by_sp.get(o.species, 0.0) + gl[g]
    return {"share": share, "species": by_sp, "groups": list(groups)}


def design_loglik(gl: Mapping[str, np.ndarray], groups: Sequence[str],
                  parts: Optional[Mapping[str, Any]] = None) -> np.ndarray:
    """Per-draw log L of a design: shares integrated over e (uniform grid), each species' tau offset integrated over
    its Gauss-Hermite nodes jointly across every arm in the design (one laboratory)."""
    parts = parts or design_parts(gl, groups)
    n = next(iter(gl.values())).shape[0]
    tot = np.zeros(n)
    if parts["share"] is not None:
        tot += _lse(parts["share"], 1) - math.log(len(E_GRID))
    for v in parts["species"].values():
        tot += _lse(v + GH_LOGW[None, :], 1)
    return tot


def _sn(llf: np.ndarray, lw: np.ndarray) -> float:
    return float(_lse((llf + lw)[None, :], 1)[0] - _lse(lw[None, :], 1)[0])


def lppd(ll: np.ndarray, logw: Optional[np.ndarray] = None, half_check: bool = False) -> Dict[str, float]:
    """The self-normalised importance estimate log sum w L / sum w over the draw set (w = prior / proposal; 1 for
    band-prior draws), with PSIS k-hat and ESS of w L (predictive_densities.importance_lppd); with ``half_check`` also
    the raw estimate on each half of the draw set (a convergence check on the decision)."""
    from src.kinetic_core import predictive_densities as pdens
    llf = np.where(np.isfinite(ll), ll, -np.inf)
    lw = np.zeros_like(llf) if logw is None else np.asarray(logw, float)
    est = pdens.importance_lppd(llf + lw, np.zeros_like(llf))
    norm = float(_lse(lw[None, :], 1)[0] - math.log(len(lw)))
    out = {"lppd": float(est.lppd) - norm, "lppd_raw": _sn(llf, lw), "k_hat": float(est.k_hat),
           "ess": float(est.ess), "resolved": bool(est.accepted)}
    if half_check:
        h = len(llf) // 2
        out["lppd_half1"] = _sn(llf[:h], lw[:h])
        out["lppd_half2"] = _sn(llf[h:], lw[h:])
    return out


def adequacy(p: Prediction, d: Data, form: Form, parts: Mapping[str, Any]) -> Tuple[float, np.ndarray, np.ndarray]:
    """Posterior-predictive p of the design's read, non-censored values: a diagonal Gaussian with the posterior
    predictive mean and spread over (draw, e, tau offset), plus the measurement sd. The joint posterior over draws and
    the e grid weights the shares; each species' offset enters through its posterior mean and variance per draw.
    Returns (p, predictive means, obs indices)."""
    from scipy.stats import chi2
    n = p.n
    lev_tot = np.zeros(n)
    dmean: Dict[str, np.ndarray] = {}
    dvar: Dict[str, np.ndarray] = {}
    for sp, v in parts["species"].items():
        lv = v + GH_LOGW[None, :]
        z = _lse(lv, 1)
        lev_tot += z
        q = np.exp(lv - z[:, None])
        dmean[sp] = q @ GH_DELTA
        dvar[sp] = np.maximum(q @ GH_DELTA ** 2 - dmean[sp] ** 2, 0.0)
    joint = (lev_tot + p.logw)[:, None] + (parts["share"] if parts["share"] is not None else np.zeros((n, 1)))
    joint = joint - np.max(joint)
    wj = np.exp(joint)
    wj /= wj.sum()
    w = wj.sum(axis=1)
    idx = [i for g in parts["groups"] for i in GROUPS[g]]
    lev = np.array([i for i in idx if OBS[i].kind == "level" and not d.censored[i]], int)
    shi = np.array([i for i in idx if OBS[i].kind == "share" and d.read[i]], int)
    mus, vs, ys = [], [], []
    if lev.size:
        g = p.loglevel[:, lev] + np.stack([dmean[OBS[i].species] for i in lev], 1)
        ovar = np.stack([dvar[OBS[i].species] for i in lev], 1)
        m = w @ g
        v = w @ (g - m) ** 2 + w @ ovar + np.array([level_sd(OBS[i], d.n_heat, form) for i in lev]) ** 2
        mus.append(m), vs.append(v), ys.append(d.y_log[lev])
    if shi.size:
        eg = E_GRID if parts["share"] is not None else np.zeros(1)
        s = p.share_a[:, shi][:, :, None] + p.share_b[:, shi][:, :, None] * eg[None, None, :]   # [N, V, G]
        m = np.einsum("ng,nvg->v", wj, s)
        v = np.einsum("ng,nvg->v", wj, (s - m[None, :, None]) ** 2) + \
            np.array([share_sd(OBS[i], d.n_heat, form) for i in shi]) ** 2
        mus.append(m), vs.append(v), ys.append(d.share[shi])
    if not mus:
        return 1.0, np.array([]), np.array([], int)
    mu, var, yv = np.concatenate(mus), np.concatenate(vs), np.concatenate(ys)
    t = float(np.sum((yv - mu) ** 2 / var))
    return float(chi2.sf(t, len(mu))), mu, np.concatenate([lev, shi])


def means_differ(mu_a: np.ndarray, idx_a: np.ndarray, mu_b: np.ndarray, idx_b: np.ndarray, form: Form) -> bool:
    """The adequacy rule's escape: the rivals' predictive means differ by more than omega on a level (0.05 on a
    share, the declared analogue) somewhere both read."""
    common = np.intersect1d(idx_a, idx_b)
    if common.size == 0:
        return False
    ma = dict(zip(idx_a.tolist(), mu_a.tolist()))
    mb = dict(zip(idx_b.tolist(), mu_b.tolist()))
    for i in common.tolist():
        tol = form.omega if OBS[i].kind == "level" else 0.05
        if abs(ma[i] - mb[i]) > tol:
            return True
    return False


# ============================================================================================ designs

def _g(prefix: str, names: Sequence[str]) -> List[str]:
    return [f"{prefix}:{n}" for n in names]


A1_LEVELS = _g("a1", ("MFT", "NF"))
A1_LABELS_RING = _g("a1", ("L1_M", "L1_ring", "L6_M", "L6_ring", "L2_M", "L2_ring", "CAM", "NF_CAM", "NF_L1_M",
                           "L1_M_120", "L1_ring_120", "L2_M_120", "L2_ring_120"))


def _swap_frag(groups: Sequence[str], frag: Optional[str]) -> List[str]:
    out = []
    for g in groups:
        if g.endswith("_ring") or "_ring_" in g:
            if frag is None:
                continue
            g = g.replace("_ring", f"_{frag}")
        out.append(g)
    return out


A2_SPECIES = ("MFT", "NF", "FA", "DPO1", "DPO3", "DDQ")
A3_SPECIES = ("MFT", "AP", "CYS")         # the decider: MFT per mol Amadori consumed against per mol Cys bound
DESIGNS: Dict[str, List[str]] = {
    "arm1": A1_LEVELS + A1_LABELS_RING,
    "arm1_frag_acyl": A1_LEVELS + _swap_frag(A1_LABELS_RING, "acyl"),
    "arm1_frag_none": A1_LEVELS + _swap_frag(A1_LABELS_RING, None),
    "arm1_L12_for_L2": A1_LEVELS + [g for g in A1_LABELS_RING if not g.startswith("a1:L2_")] + _g("a1", ("L12_M1", "L12_M2"))
    + _g("a1", ("L2_M_120", "L2_ring_120")),
    "arm1_no120": A1_LEVELS + [g for g in A1_LABELS_RING if not g.endswith("_120")],
    "arm1_L1_L6_CAM_only": A1_LEVELS + _g("a1", ("L1_M", "L1_ring", "L6_M", "L6_ring", "CAM", "NF_CAM", "NF_L1_M")),
    "arm1_levels_only": A1_LEVELS,
    "arm2": _g("a2", A2_SPECIES),
    "arm2_plus_air": _g("a2", A2_SPECIES) + _g("a2air", A2_SPECIES),
    "arm2_cys_sulfur": _g("a2cys", A2_SPECIES),
    "arm3": _g("a3", A3_SPECIES),
    "arm3_sulfur_matched": _g("a3m", A3_SPECIES),
    "arm3_plus_opd_glucosone": _g("a3", A3_SPECIES + ("GOS",)),
}
DESIGNS["arm1+arm2"] = DESIGNS["arm1"] + DESIGNS["arm2"]
DESIGNS["arm1+arm3"] = DESIGNS["arm1"] + DESIGNS["arm3"]
DESIGNS["arm2+arm3"] = DESIGNS["arm2"] + DESIGNS["arm3"]
DESIGNS["arm1+arm2+arm3"] = DESIGNS["arm1"] + DESIGNS["arm2"] + DESIGNS["arm3"]
PRIMARY_DESIGNS = ("arm1", "arm2", "arm3", "arm1+arm2", "arm1+arm3", "arm2+arm3", "arm1+arm2+arm3")
DESIGN_ARMS = {"arm1": ("a1",), "arm2": ("a2",), "arm3": ("a3",)}


def arms_of(design: str) -> Tuple[str, ...]:
    av = sorted({g.split(":")[0] for g in DESIGNS[design]})
    return tuple(av)


def arm_split(design: str) -> List[List[str]]:
    """A design's groups by designed arm (for Holm across arms)."""
    out: Dict[str, List[str]] = {}
    for g in DESIGNS[design]:
        out.setdefault(g.split(":")[0][:2], []).append(g)
    return list(out.values())


def heatings(design: str, n_heat: int) -> int:
    """Heatings in a design: unlabelled pots + one labelled pot per label per cell, each x n."""
    pots = set()
    lab = set()
    for g in DESIGNS[design]:
        for i in GROUPS[g]:
            o = OBS[i]
            if o.kind == "level":
                pots.add(o.pot)
            else:
                lab.add((o.pot, "L12" if o.read.startswith("L12") else ("CAM" if o.read == "CAM" else o.read.split("_")[0])))
    return n_heat * (len(pots) + len(lab))


# ============================================================================================ job (a): label arithmetic

def d3_arithmetic(q: float = 0.5, printed_floor: float = 0.01, share_sd_abs: float = 0.01) -> Dict[str, Any]:
    """C2 + C3 in a [13C6] / [12C6] pool of labelled fraction q: two fragments from independent molecules, so the mixed
    share is 2 q (1 - q) times the route's share f. Against D3's printed '<= 1 %' mixed cell (+-1 % absolute, a 1 %
    floor), constant-free: no rate constant enters."""
    mixed = 2.0 * q * (1.0 - q)
    bound = printed_floor + 1.96 * share_sd_abs
    z_full = (mixed - printed_floor) / share_sd_abs
    return {"q": q, "mixed_share_c2c3": mixed, "q_band_0.45_0.55": [2 * 0.45 * 0.55, mixed],
            "z_full_route": z_full, "log_lik_ratio_nats_full_route": 0.5 * z_full ** 2,
            "max_route_share_f": bound / mixed, "r1_plus_f_r3_at_f_0.1": 0.1 * mixed,
            "intact_routes_mixed_share": 0.0}
