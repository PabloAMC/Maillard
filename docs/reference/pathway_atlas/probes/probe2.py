from src.kinetic_core import engine, sulfur
from src.kinetic_core.ph_state import BufferSpec
P = engine.core_parameters("sulfur")
cl = BufferSpec(kind="clamped")
def g(f, *ks): return sum(f.get(k, 0.0) for k in ks)
MFT = ("r_ddp_mft", "r_ddp_mft_hs", "r_nf_mft", "r_ha_mp_mft", "r_hmp_mft")
print("E. Hofmann 1998 T1 analogue: sugar 100 + Cys 33 mM, pH 5, 145 C, 20 min (lit: ribose/glucose MFT 10.4x, FFT 4.3x)")
out = {}
for s in ("PENT", "Glc"):
    f = sulfur.sulfur_flux_budget(P, 418.15, {s: 100.0, "Cys": 33.0}, 20.0, ph=5.0, buffer_spec=cl)
    out[s] = (g(f, *MFT), g(f, "r_fur_fft", "r_fur_fft_hs"))
    print(f"  {s}: MFT made {out[s][0]:.3e} mM, FFT made {out[s][1]:.3e} mM")
print(f"  ratio PENT/Glc: MFT {out['PENT'][0]/max(out['Glc'][0],1e-300):.3g}, FFT {out['PENT'][1]/max(out['Glc'][1],1e-300):.3g}")
for key in ("k_nf_decay", "k_mft_decay", "k_tdp_fur", "k_pent_tdp", "k_pent_caramel"):
    p = P.get(key)
    try:
        print(key, "k(100C)=", p.k_at(373.15), "k(145C)=", p.k_at(418.15), "ph_kind", sulfur.REACTION_PH_FACTOR.get("r"+key[1:], ""))
    except Exception as e:
        print(key, type(p), e)
