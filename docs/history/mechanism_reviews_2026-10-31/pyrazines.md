# Pyrazine kinetics: review (2026-10-31)

**Current lane** (`parameters_pyrazine.py`, B18):
- Structure: GO/MGO + Gly -> aminoketone, fitted on Zhou 2024 as "pyrazine rate = Strecker/2"; condensation constant declared fast; pH law from Leahy 1989a.
- Deficit: Leahy lysine + glucose, 95 C / 2 h / pH 9: model 21.9 vs 13100 ug/L (-2.78 dex; B21 report).

## 1. Rate-law structure (ranked)

**A. Strecker-limited, fast 2nd-order condensation, competing aminoketone sink — P 0.50**
- Hofmann 2000 Table 4 (MGO + Phe, 0.1 M, pH 7, 98 C, 60 min): Strecker products 1.58 mol %; B18's constant predicts 85x less (1.9 dex).
- Adams 2008 Table 3: aminoacetone goes to 2,5-diacetyl-3-methylpyrrole (1.1 - 3.7 mol %), against 2,5(6)-DMP at 0.02 %; the pyrrole grows with excess MGO.
- So pyrazine yield per aminoketone is pot-dependent: about 1 % with fed dicarbonyl, higher in sugar pots.

**B. Today's lumped yield constant — P 0.15**
- It holds within fed pots: Zhou's 20 mM rate scaled by c^2 matches Adams' 850 mM rate within 10 %.
- It does not transfer to sugar pots.

**C. Condensation- or oxidation-limited (O2) — P 0.15**
- Jousse 2002 calls the condensation (reaction R10) "fast" but prints no number.
- Shibamoto 1977 proposes dehydration of a hydroxy-dihydropyrazine, which needs no oxidant.
- There is no pyrazine O2 data on disk.

**D. Mostly non-Strecker nitrogen (Amadori-retained N, NH3, Namiki C2) — P 0.20**
- Zhou 2024 Table 1: only 54 % of pyrazine N comes from free 15N-Ala.
- Leahy's products rise linearly from t ~ 0 with no lag.

**Order:** second order in aminoketone.
- Amrani-Hemaimi 1995: no amino-acid carbon in methylpyrazine or 2,5-DMP.
- Aldehyde addition is real: 25 - 30 % of a 2,6-DMP methyl, and 20 - 100 % of ethyl / trimethyl substituents, come from the amino acid.
- Cerny 1994: ethyl isomers 3,6 : 3,5 = 10 : 1.

**pH (Leahy 1989a):**
- k9/k7 = 2.7 and k7/k5 = 14, i.e. 0.58 dex/unit below pH 7 and 0.20 above.
- This is a whole-cascade law, yet the lane puts it on the Strecker steps only.
- The dimethylpyrazines are not quantifiable at pH 5.

**Temperature:**
- Whole cascade: Ea 114 - 187 kJ/mol (Leahy ch. 7 Table II).
- Single steps: 100 - 118 (Zhou; Yu 2018 Table 1).
- Model apparent Ea: 300 - 460, the signature of a lagged supply.

## 2. Causes of the 2.8-decade deficit (P = this is the largest contributor)

1. **0.55 — a yield constant used as a Strecker flux.** Worth about 1.9 dex.
2. **0.25 — missing amines and nitrogen routes, worth about 0.3 - 0.8 dex.**
   - Lysine's epsilon-amine: Hwang 1995 finds the Gly + Lys pot gives 4.3x glycine's total; Leahy finds lysine about 6x asparagine.
   - Amadori-retained nitrogen and NH3 are also absent.
3. **0.12 — wrong supply shape and temperature dependence.** The model lags; Leahy is linear from t ~ 0.
4. **0.08 — measurement or units.**
   - Leahy's ppm basis (x1.5).
   - Zhou's SPME calibration matrix is unstated.

## 3. Constants on disk and scale-free signatures

**Constants on disk:** Zhou Table 2; Leahy ch. 7 Tables I - II and ch. 18 (pH 5 / 7 / 9); Yu 2018 Table 1 (time unit unresolved); Adams Table 3; Cerny 1994; Hwang 1995; Hofmann 2000 Tables 3 - 4.

**Not found, on disk or on the web:** an aminoketone-dimerisation constant or a dihydropyrazine-oxidation constant.

**Scale-free signatures:**
- **Mixing test:** MPZ^2 / (PZ x 2,5-DMP) = 11.7 / 14.4 / 17.3 at 95 / 85 / 75 C (my arithmetic, Leahy Table I). One shared condensation constant predicts 4, so methylpyrazine has an extra route.
- Also: Leahy's pH-5 pyrazine share (63 %), Zhou's 15N share (54 %), Amrani's carbon origins, Cerny's isomer ratio, Adams' pyrrole / pyrazine partition.
- **Peptides** (Van Lancker 2012, peak-area ratios): GlyGly / Gly = 34 and AlaLys / (Ala + Lys) = 5.1 for DMP from MGO. Isolates need amine classes: N-terminal, epsilon-Lys, and NH3 from deamidation.

**Next:** add Hofmann's Strecker constant and Adams' aminoketone + dicarbonyl -> pyrrole sink, refit, and test on the mixing ratio and Leahy's total.

## 4. Papers not on disk

Leahy & Reineccius are already on disk.

1. **Koehler & Odell 1970**, JAFC 18:895, "Factors affecting the formation of pyrazine compounds in sugar-amine reactions", 10.1021/jf60171a041.
   - DOI checked on Crossref.
   - Abstract not seen; content unverified.
2. **Shibamoto & Bernhard 1976**, JAFC 24:847, "Effect of time, temperature, and reactant ratio on pyrazine formation in model systems", 10.1021/jf60206a045.
   - DOI checked on Crossref.
   - Abstract not seen; content unverified.
3. **Shibamoto & Bernhard 1977**, Agric. Biol. Chem. 41:143, "Investigation of pyrazine formation pathways in glucose-ammonia model systems", 10.1271/bbb1961.41.143.
   - Checked on Crossref and J-STAGE.
   - Reported to contain alpha-aminocarbonyl intermediates and ring closure by dehydration, studied with O2, metal ions and NaOH.
4. **Koutsidis et al. 2008**, JAFC 56:6105, "Acrylamide and pyrazine formation in model systems containing asparagine", 10.1021/jf703744k.
   - Checked on Crossref and PubMed.
   - Reported to contain kinetic data, and the unsubstituted pyrazine as the only pyrazine from glyoxal + Asn.
5. **Low, Parker & Mottram 2007**, JAFC 55:4087, "Mechanisms of alkylpyrazine formation in a potato model system containing added glycine", 10.1021/jf070044s.
   - Checked on PubMed.
   - Reported to contain [2-13C]glycine methyl groups in alkylpyrazines and two pathways with their relative contributions.
6. **Scalone et al. 2015**, JAFC 63:5364, "Influence of free amino acids, oligopeptides, and polypeptides on the formation of pyrazines in Maillard model systems", 10.1021/acs.jafc.5b01129.
   - Checked on Crossref and PubMed.
   - Reported to show oligopeptides raising pyrazine formation more than free amino acids, with more pyrazine at aw 0.33 and above 120 C.
7. **Van Lancker, Adams & De Kimpe 2010**, JAFC 58:2470, "Formation of pyrazines in Maillard model systems of lysine-containing dipeptides", 10.1021/jf903898t.
   - Checked on Crossref and PubMed.
   - Reported to show more pyrazines from Lys dipeptides than from free amino acids.
8. **Hayashi & Namiki 1980**, Agric. Biol. Chem. 44:2575, "Formation of two-carbon sugar fragment at an early stage of the browning reaction of sugar with amine", 10.1271/bbb1961.44.2575.
   - Checked on Crossref and J-STAGE.
   - Reported to contain an early C2 glyoxal-dialkylimine from glucose + cyclohexylamine.
