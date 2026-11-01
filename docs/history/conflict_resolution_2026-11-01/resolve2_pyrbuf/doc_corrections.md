# resolve2_pyrbuf -- proposed corrections to atlas, dossiers and docs (nothing applied; writer's call)

## A. Pathway atlas

**A1. F15 -> F07 (Huber 1990; CF-NEW-7).** In `F15_redox_reductones.yml`, the McGill (Yaylayan) cell
(`T_C: 250`, `quote: "~5 % via the 1-amino-1,4-dideoxy-2,3-hexodiulose RA[4,5] (B)"`) becomes two cells.

(a) It stays in F15, narrowed to the redox side experiment:
```yaml
  - key: {ph: neutral, temperature: moderate, water: dilute}
    tags: {ph_basis: label, T_C: 100, sugar: none, amine: none, oxygen: not_stated, buffer: none_unbuffered, metals: not_controlled}
    contrast: "pyruvaldehyde alone, 100 C, 8 min: 1,2-propanediol (MGO + 2H) without added reductant"
    status: operating
    quote: "1,2-propanediol 18 %"
    quantity: "1,2-propanediol 18 %, acetic acid 30 %, acetol 7 % of peak area (redox without oxidant; disproportionation)"
    evidence: fed_intermediate
    lab: "McGill (Yaylayan)"
    dossiers: [huyghuesdespointes1996_extraction.md]
    where: "redox side experiments"
    note: "Peak-area shares, not yields. The 1-amino-1,4-dideoxy RA[4,5] share moved to F07 (Huber & Ledl 1990: an elimination product, not a reduction)."
```
Check the water bin against the dossier first (neat or aqueous MGO).

(b) A new F07 cell for the share, plus a Huber cell:
```yaml
  - key: {ph: neutral, temperature: very_high, water: low_aw}
    tags: {ph_basis: label, T_C: 250, sugar: hexose_aldose, amine: other_free_aa, oxygen: not_stated, buffer: none_unbuffered, metals: not_controlled}
    contrast: "six singly labelled glucoses + proline, pyrolysis 20 s: C2 / C4 products by route"
    status: minor
    quote: "~5 % via the 1-amino-1,4-dideoxy-2,3-hexodiulose RA[4,5] (B)"
    quantity: "the 1-amino-1,4-dideoxy-2,3-hexodiulose (C-4 elimination from the 2,3-enediol, amine retained) retro-aldol route ~5 % of C2 / C4 products; 1-DG RA[4,5] ~25 % (author inference from labels)"
    evidence: labelling
    lab: "McGill (Yaylayan)"
    dossiers: [huyghuesdespointes1996_extraction.md, huber1990_extraction.md]
    where: "Scheme 1, Tables 4-6"
    note: "Moved from F15 (2026-11-01): Huber & Ledl 1990 assign the 1-amino-1,4-dideoxy compound to beta-elimination at C-4 (the 1-deoxyosone's oxidation level). Pathway shares are inferences; dry pyrolysis; secondary amine."
  - key: {ph: mild_acid, temperature: moderate, water: dilute}
    tags: {ph_basis: label, T_C: 90, sugar: disaccharide, amine: other_amine, oxygen: not_stated, buffer: phosphate, metals: not_controlled}
    contrast: "maltose vs glucose propylamino ARPs, pH 5, reflux in water / dioxane 1:1, OPD present: a third exit from the 2,3-enediol"
    status: operating
    quote: "formed in higher amounts"
    quantity: "the 1-amino-1,4-dideoxy-2,3-hexodiulose (quinoxaline 7b) from the maltose ARP; from the glucose ARP 'possibly' (not definite); no yields"
    evidence: trapping
    lab: "Ledl, Universitat Stuttgart"
    dossiers: [huber1990_extraction.md]
    where: "Results; Scheme"
    note: "No numbers. A 4-O-glycosyl leaving group favours C-4 elimination over amine release; for free glucose ARPs at most minor. The engine's hexose 2,3-branch is unchanged."
```
Before writing, check the bin vocabulary (`mild_acid` for pH 5; `T_C` ~88 - 90 C, i.e. "moderate") against
families.yml. Update F15's `condition_dependence` and `markers` text where they cite McGill for the 1,4-dideoxyosone
route.

**A2. Atlas summary P2 (README / ENGINE_VS_LITERATURE sec. 1.4 line "P2 ... fading at pH 8")** becomes: "Phosphate
accelerates early condensation, Amadori decomposition and free-sugar degradation (x1.3 - 12). The law is linear in
buffer (no saturation to 0.5 M; Kallen's thiazolidine plateau is a change of rate-determining step). The catalysed
increment peaks broadly at pH 6 - 7, and the RATIO falls at pH 8 because the uncatalysed rate rises."
F03 summary row: replace "fades at pH 8" with "ratio falls at pH 8 (uncatalysed rate rises)".

**A3. F01 / F02 / F03 / F14 condition-dependence text:** add the specificity finding. Carboxylates and bicarbonate act
on free-sugar degradation and browning (Saunders, Rizzi, Davies, Larralde). They do not act on Amadori decay (Guan,
single lab).

**A4. F05:** state that Cui 2026's a_w maxima are 20-min end points of consumed intermediates (a_w set before heating, open
dishes). Add Wolf 1977's plateau (pass 2: no loss) as a cell when keyed, with the CF-NEW-4 pointer.

## B. ENGINE_VS_LITERATURE.md

**B1. Sec. 3, row 5 (E5).** Replace "**Fix:** a buffer-anion factor on the Amadori-decay and enolisation steps (saturating,
fading at pH 8)" with:

> **Fix:** a per-step factor linear in buffer, k = k_u (1 + sum kappa_b c_b). Phosphate kappa is 5 - 40 /M at pH 5 - 7.4
> (DEV: Guan, Davidek, Zyzak, Bell) and roughly T-flat 25 - 110 C. Carboxylate kappa is 0 - 3 /M on Amadori steps, and
> 0.1 - 1x phosphate on free-sugar steps. There is no Michaelis term; the ratio fade at pH 8 comes from k_u. Calibration
> constants are divided by their own buffer factor. TEST rows Potman, Larralde, Kallen and Thornalley pass
> (resolve2_pyrbuf).

**B2. Sec. 3.1, row E5 (designed experiment).** Add a citrate arm "at equal in-situ pH". Name the decider: "citrate = water
on Amadori decay (step-specific) vs citrate 0.3 - 1x phosphate (general)". Drop "saturation" from the question; the
form is settled at the 0.1 - 0.5 M scale.

**B3. Sec. 3.1, row E4.** Add: "measure the Strecker aldehyde (and acid) and the aminoketone, not amino-acid loss.
Amino-acid consumption by fed MGO is dominated by non-Strecker adducts, ~1e3x (Piloty 1979; CF-NEW-3)."

## C. Dossiers (our own text)

**C1. `piloty1979_extraction.md` sec. 5, E4 bullet.** "SUPPORTS the remaining E4 claim ... Weight: moderate (P ~0.6)" and
"would predict Ala loss of ~3 % at pH 4.5 where 32 % is measured" assume that all consumption is the Strecker step.
Replace with:

> Not discriminating (P 0.85). Under the shipped step law, pH-flat consumption only requires the Strecker share f <= 0.04
> (Ala), 0 (Leu), 0.2 - 0.3 (Gly, Phe, Ile). In the engine the step is 1e-5 - 1e-3 of the measured loss under both rivals.
> It bears on E4 only with a pH-independent Strecker share (an assumption; LR ~1.5 for supply). It exposes a missing
> adduct sink (CF-NEW-3).

`scratchpad/ingest_findings.md` batch-2 item 4 ("Piloty supports E4's step-law error") needs the same downgrade.

**C2. `kallen1971b_extraction.md` sec. 1 / answer 3.** After "Phosphate catalysis SATURATES", add: "the plateau is the
uncatalysed amine attack. The buffer-free attack plateau (Fig. 6, 10^3.4) predicts the ceilings within 0.2 dex. This is a
change of rate-determining step, not buffer binding, so it does not license a Michaelis form on other steps (E5)."

**C3. `saunders1966_extraction.md` sec. 5.** Add: "The pH-8 'fade' in the E5 fix is a ratio on the Amadori step (Davidek: the
uncatalysed rate doubles from pH 7 to 8). Browning at pH 7.5 - 8 integrates >= 2 catalysed steps. No contradiction in
kind; the fix text is rewritten (B1)."

**C4. `potman1989_extraction.md` sec. 1.** κ_eff "rises from ~30 M⁻¹ at 0.1 M to ~55 M⁻¹ at 0.5 M ... super-linear" is on
initial-rate extrapolations. On the 5-h secant the per-M increment is flat, 22 - 29 /M (0.1 - 0.5 M); at 3 h it is
48 - 81 /M. Add this as the form result. The super-linearity is not robust to the rate extraction.

**C5. `larralde1995_extraction.md` (and the keyed table's TEXT CHECKS).** Add the derived bicarbonate activity: 0.23x
phosphate per mole, from 0.52 and the rate law in min-1. It holds only if the bicarbonate run is at 100 C (P 0.8) and
its pD held (CO2 loss would raise pD and the rate, which biases toward catalysis).

**C6. `wolf1977_extraction.md` sec. 1, answer 2.** Add: "the shipped lysine-returning structure fails too (pass 2 -36 %,
pseudo-steady free fraction 0.21). The plateau is not a CF-072 decider (CF-NEW-4)."

**C7. `lee2024_extraction.md` / keyed table.** Add: "R ratios within a pot (pool-invariance tests) are independent of the
MPZ x10. Only levels need the dual reading."

**C8. `shibamoto1976_extraction.md`.** Add: "R stays at 366 - 496 while the total yield moves x7.6 (30 min - 48 h). Pool
invariance rejects a lower-order extra MPZ route at |z| 3 - 6. The 15-min point (R 1380, yield 0.013 %) is the
exception." Add a pointer to CF-NEW-6 (Wong 14.6 vs 357 - 522, same lab).

**C9. `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` row 10.** As flagged in the Saunders dossier: no citrate arm
at 100 C.

## D. Owner calls (not ours)

- `lab_ids.yml`: davies2005 -> `leeds_wedzicha`. Recommend merging (the conservative rule); then it is DEV for E5, which
  this analysis assumed nowhere (it used Davies only for form, as unsplit).
- New lab ids: `oklahoma_state_koehler`, `cerebos_edinburgh` (Saunders), `tu_berlin_baltes` (Piloty), `nankai_cui`.
- Download (DOI): Wolf, Thompson & Reineccius 1977a, J. Food Sci., 10.1111/j.1365-2621.1977.tb08420.x (the lysine-loss
  constants and glucose levels behind CF-NEW-4).
