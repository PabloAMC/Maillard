# Verification probe (read-only): shipped sulfur lane, clamped pH, flux budgets for three mismatch claims.
from src.kinetic_core import engine, sulfur
from src.kinetic_core.ph_state import BufferSpec
P = engine.core_parameters("sulfur")
cl = BufferSpec(kind="clamped")
def fb(init, T, t, ph):
    return sulfur.sulfur_flux_budget(P, T + 273.15, init, t, ph=ph, buffer_spec=cl)
def g(f, *ks): return sum(f.get(k, 0.0) for k in ks)
print("A. MFT route split, ribose+Cys 30 mM, 95 C, 240 min (Cerny-like)")
for ph in (5.0, 6.0):
    f = fb({"PENT": 30.0, "Cys": 30.0}, 95.0, 240.0, ph)
    tot = g(f, "r_ddp_mft", "r_ddp_mft_hs", "r_nf_mft", "r_ha_mp_mft", "r_hmp_mft")
    print(f"  pH {ph}: MFT made {tot:.3e} mM; NF {g(f,'r_nf_mft')/tot:.3f}  DDP {g(f,'r_ddp_mft','r_ddp_mft_hs')/tot:.3f}  C2+C3 {g(f,'r_ha_mp_mft')/tot:.3f}")
print("B. FFT formation vs pH, ribose+Cys 30 mM, 145 C, 20 min")
for ph in (3.0, 4.5, 5.0, 6.0, 7.0):
    f = fb({"PENT": 30.0, "Cys": 30.0}, 145.0, 20.0, ph)
    print(f"  pH {ph}: FUR made {g(f,'r_tdp_fur','r_glc_fur'):.3e}  FFT made {g(f,'r_fur_fft','r_fur_fft_hs'):.3e} (HS- share {g(f,'r_fur_fft_hs')/max(g(f,'r_fur_fft','r_fur_fft_hs'),1e-300):.3f})  FFT lost thiolate {g(f,'ch_thiolate_loss_fft'):.3e} decay {g(f,'ch_fft_decay'):.3e}")
print("C. MFT route split, glucose+Cys 30 mM, 145 C, 20 min, pH 6 (hexose)")
f = fb({"Glc": 30.0, "Cys": 30.0}, 145.0, 20.0, 6.0)
tot = g(f, "r_ddp_mft", "r_ddp_mft_hs", "r_nf_mft", "r_ha_mp_mft", "r_hmp_mft")
print(f"  MFT made {tot:.3e}; NF {g(f,'r_nf_mft')/tot if tot else 0:.3f} DDP {g(f,'r_ddp_mft','r_ddp_mft_hs')/tot if tot else 0:.3f} C2+C3 {g(f,'r_ha_mp_mft')/tot if tot else 0:.3f}; r_odg_dpo flux {g(f,'r_odg_dpo'):.2e}")
print("D. MFT vs pH, ribose+Cys 30 mM 145 C 20 min")
for ph in (4.5, 5.0, 6.0, 7.0):
    f = fb({"PENT": 30.0, "Cys": 30.0}, 145.0, 20.0, ph)
    print(f"  pH {ph}: MFT made {g(f,'r_ddp_mft','r_ddp_mft_hs','r_nf_mft','r_ha_mp_mft'):.3e}; NF share {g(f,'r_nf_mft')/max(g(f,'r_ddp_mft','r_ddp_mft_hs','r_nf_mft','r_ha_mp_mft'),1e-300):.3f}")
