# Wave B18 ship rule: SHIP

*Rule: SHIP if T1 (six Zhou rates within 0.3 dex, barriers inside band), T2 (no scored panel value moves 0.05 dex) and T5 (log10 k sigma below one decade) hold. Pre-registration `results/validation/kinetic_core_b18_prereg.md`.*

| test | result | pass |
|---|---|---|
| T1 Zhou rates | worst zhou_DMP_rate_100C +0.071 dex; barriers inside band {'ea_go_ak_kj_mol': True, 'ea_mgo_ak_kj_mol': True} | True |
| T2 panel unchanged | 39 predicted numbers compared; worst 0.0000 dex (mp_holdout_glucose_only_autoclave_121C_Steinhagen2021) | True |
| T3 Leahy distribution | pyrazine : 2,5-DMP model 5.94e-08 vs 23.2 (-8.59 dex); methylpyrazine/pyrazine model 3.95e+03 vs 0.744 (+3.73 dex) | False (mixed route ships: False) |
| T4 Leahy total | model 17.4 vs 1.31e+04 ug/L (-2.88 dex) | False |
| T5 identification | log10 k sigma {'log10_k_go_ak_100C': 0.08183260861717218, 'log10_k_mgo_ak_100C': 0.07532045007455257}; barrier sigma {'ea_go_ak_kj_mol': 16.157518353952135, 'ea_mgo_ak_kj_mol': 16.096566410277887} kJ/mol, on bound {'ea_go_ak_kj_mol': True, 'ea_mgo_ak_kj_mol': True} | True |
| T6 pH direction | k(5)/k(9) = 0.0257 (band 1/60 to 1/20) | True |

## Barriers against the hold-out laboratories

- Yu 2018, 2,5-dimethylpyrazine, glucose + glycine pH 10, 70-90 C: model apparent 390.3 kJ/mol vs 99.8 +/- 6.7
- Leahy 1989, PZ, pH 9, 75-95 C: model apparent 459.3 vs whole-cascade 149.8
- Leahy 1989, MPZ, pH 9, 75-95 C: model apparent 378.8 vs whole-cascade 153.1
- Leahy 1989, DMP, pH 9, 75-95 C: model apparent 303.5 vs whole-cascade 177.0

## The glyoxal-sink conditionality

- The B13 glyoxal sink (Kocadagli's dry glass at 180 C, barrier fixed to zero) removes 98 % of a fed 20 mM glyoxal in two hours at 100 C, so the pyrazine growth in the modelled Zhou pot is not linear (Zhou's is, figure-only) and the fitted glyoxal Strecker constant is higher by the shift above than it would be over a constant pool. The methylglyoxal pot loses its dicarbonyl too, through the trunk's own B1 melanoidin sink and the B7 furanone step; those are fitted or measured constants and are not zeroed here. Recorded as a conditionality on every pyrazine answer; the fix is a wave on the dicarbonyl sinks in water, not a move of this one.

- glyoxal_remaining_fraction_at_120min_in_zhou_pot: {'zhou_go_100_GO_remaining_fraction_at_120min': 0.019972110306741266, 'zhou_go_120_GO_remaining_fraction_at_120min': 0.019847375241134658, 'zhou_mgo_100_MGO_remaining_fraction_at_120min': 0.06539707414672412, 'zhou_mgo_120_MGO_remaining_fraction_at_120min': 0.022614623994501926}
- methylglyoxal_remaining_fraction_at_120min_in_zhou_pot: {'zhou_mgo_100_MGO_remaining_fraction_at_120min': 0.06539707414672412, 'zhou_mgo_120_MGO_remaining_fraction_at_120min': 0.022614623994501926}
- pyrazine_rate_0_60_over_0_120_min: 1.749919179776949
- log10_k_go_ak_shift_when_the_dry_glass_glyoxal_sink_is_zeroed: -0.5860313013873828
- log10_k_mgo_ak_shift: 0.01295200166630206
- nosink_cost: 2.415066227158778
- b18_cost: 2.644103944369364

## Leahy 95 C / 2 h, pH 9 (glycine for lysine)

- model ug/L: {'PZ': 1.0354932963337353e-06, 'MPZ': 0.004090716199320298, 'DMP': 17.431728259589768}
- model shares %: {'PZ': 5.938884983119162e-06, 'MPZ': 0.02346156473669305, 'DMP': 99.97653249637833}
- Leahy shares %: {'PZ': 55.8, 'MPZ': 41.5, 'DMP': 2.4}
