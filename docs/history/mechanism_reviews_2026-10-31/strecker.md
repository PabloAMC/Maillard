# Strecker degradation: rate-law review (2026-10-31)

**The core problem: the fed-donor yield is a branching share, not a rate.**
- In Hofmann 2000 Table 3, phenylacetaldehyde (PA) reaches 60 % of its 500-min value by 10 min, and PA + phenylacetic acid (PAA) levels off at about 4 % of the methylglyoxal (MGO) charge.
- The donor is gone in about 15 min.
- So the rows identify the Strecker share s_D = k_D / kappa_D, not the rate constant k_D. The forced-input prereg already handles this.

## 1. Rate law (ranked)

**H1 (P 0.50): competitive mass action on total amino acid, J = k_D f_i g_D(pH) [D][AA_i].** It also carries amine-dependent and amine-free sinks of the donor.
- pH factor g is shallow. From pH 3 to 9, PA + PAA rises ×2.5 for MGO, ×2.1 for glyoxal (GO), ×1.4 for 3-deoxyglucosone (3-DG) and ×1.2 for glucosone (Hofmann Figs 4-5); that is ≤ 0.07 dex per pH unit.
- Oxygen:
  - Under α-dicarbonyls the aldehyde is O2-independent (PA air / argon 0.85-0.93).
  - O2 adds a separate acid flux (PAA ×1.6-5.5, Table 4); it does not divert the aldehyde.
  - Less than 1 % of PA becomes PAA (Table 5).
- Competition: equimolar Leu, Ile, Val or Met halves Phe's products (Hidalgo 2016).

**H2 (P 0.20): Jencks bell.** g = free-amine fraction (pKa about 7.6 at 98 C) × acid catalysis. It matches H1 from pH 3 to 7, so separating the two needs pH × T data.

**H3 (P 0.15): saturating in donor.**
- For: Zamora 2013.
- Against: response is linear in carbonyl in Hidalgo 2013 and Zamora 2015.

**H4 (P 0.05): slope-1 free-amine law.** Refuted by more than 3 dex. It predicts near-zero yield at pH 3, but Hofmann measures half the pH 7 value. Caveat: Hofmann did not print the conditions at pH 3, 5 and 9.

**Partition between aldehyde and acid.**
- PA / PAA runs from 2.5 to 0.19 between pH 3 and 9 for MGO, and from 36 to 1.0 for glucosone.
- The acid is favoured by base, O2 and C2-C3 donors.
- Likely acid route (P 0.5): an α-keto acid formed by transamination, then oxidative decarboxylation.

**Donor classes.**
- **α-Dicarbonyls (aldehyde forms without oxidation):**
  - At pH 7, MGO : GO : 3-DG : glucosone = 1 : 0.83 : 0.76 : 0.73.
  - At pH 3.2, MGO : GO : 3-DG = 35 : 15 : 1 (Weenen 2001).
  - Weenen gets 21 % from MGO versus Hofmann's 1 %, a ×20 gap between labs; Weenen strips the aldehyde out by steam distillation.
- **Lipid enals (aldehyde forms by oxidation):**
  - PAC is ×4-9 higher in air than in N2; without O2 the product is the amine phenylethylamine (PEA).
  - Ea is 27.6 kJ/mol in air versus 78.0 in N2 (Zamora 2013, 2015).
  - PAC falls as pH rises.

## 2. Constants

**Measured, all apparent:**

| source | conditions | value |
|---|---|---|
| Chan 1994b | glucose, pH 6-8, 75-115 C | Ea 80 kJ/mol (3-methylbutanal, 3-MB), 90 (PA) |
| Cremer 2000 | a_w 0.52 | Ea 115-124 kJ/mol; tracks the sugar supply, not the Strecker step |
| Balagiannis 2009 | 130 C | k1 1.36e-2 /min (Ea 137); 3-MB loss 2.7e-3 /min (Ea 78 ± 39); 2-methylbutanal (2-MB) loss ≈ 0 |
| Seville group | pH 3 | PAC 28-81 kJ/mol; PEA 47-65 kJ/mol |

None of these is an α-dicarbonyl Strecker step in water: gap G1 stays open.

**Identifiable from the 111 keyed rows:**
- k_MGO at 98 C, pH 7: about 2e-5 L/(mmol min) from the plateau. The 7e-6 chord is a lower bound.
- s_D, phi_D and omega_D for each donor.
- A weak aldehyde-loss constant k_loss.
- Relative pH shapes.
- Leu / Phe and Ala / Phe ratios.

**Not identifiable:**
- Ea.
- Order in amino acid.
- Absolute k for GO, 3-DG and glucosone.
- O2 dose.
- Whether H1 or H2 holds.

## 3. Scale-free signatures on disk

1. Hofmann:
   - Table 4: air / argon per donor.
   - Table 5: less than 1 % PA → PAA.
   - Figs 4-5: PA peaks at pH 5; PAA rises monotonically.
   - Table 3: plateau.
2. Donor order at pH 7 (Hofmann) versus pH 3.2 (Weenen Table III).
3. Amino-acid identity:
   - Leu / Phe 2.1 for the aldehyde and 1.7 for the acid (Hofmann Table 2). This is side-chain-specific, not set by pKa (P 0.85).
   - Ile / Leu 1.6 (Balagiannis).
   - Cremer's barriers fall within 9 kJ/mol of each other.
   - Ea order Pro < Leu < Phe (Chan).
4. Competition (Hidalgo 2016):
   - Leu, Ile, Val and Met halve Phe's products.
   - Arg lowers PAC linearly.
5. Lipid class:
   - Air / N2 effect and PEA / PAC split (Zamora 2015).
   - Ea in air versus N2 (Zamora 2013).
   - PAC falls with pH.
6. Order in carbonyl: linear in Hidalgo 2013 and Zamora 2015, saturating in Zamora 2013.
7. Amadori product versus glucose: ×4 at pH 7, ×270 at pH 3.2 (Weenen Table II).
8. Losses:
   - 3-MB is lost; 2-MB is not.
   - PA falls 10 % by 500 min.
   - PA + 3-MB forms an aldol product (Chan).

## 4. Papers not on disk

All DOIs were checked against their Crossref records. Content comes only from abstracts actually seen.

1. **Piloty M, Baltes W 1979.** "...I. Reaktivität von Aminosäuren bei Reaktion mit α-Dicarbonylverbindungen." *Z Lebensm Unters Forsch* 168:368-373. 10.1007/BF01186494. Part II: 10.1007/BF01186495.
   - **Abstract unverified.** Only the title is known: amino-acid reactivity toward α-dicarbonyls.
2. **Granvogl M, Bugan S, Schieberle P 2006.** *J Agric Food Chem* 54:1730-1739. 10.1021/jf0525939.
   - Abstract (TUM portal): Ile, Leu, Phe and Val each heated with MGO; aldehydes far exceed amines.
3. **Pripis-Nicolau L et al. 2000.** *J Agric Food Chem* 48:3761-3766. 10.1021/jf991024w.
   - Abstract (AGRIS): diacetyl, pentanedione, GO, MGO, acetoin and acetol with wine amino acids at pH 3.5; temperature and pH studied with cysteine.
4. **Cremer DR, Eichner K 2000.** "The influence of the pH value on the formation of Strecker aldehydes..." *Eur Food Res Technol* 211:247-251. 10.1007/s002170000192.
   - **Abstract unverified.**
5. **Hidalgo FJ, Gallardo E, Zamora R 2005.** *J Agric Food Chem* 53:10254-10259. 10.1021/jf052240+.
   - Abstract (AGRIS): 4-hydroxynonenal + Phe at 37 / 60 / 80 C gives 17 % PA.
6. **Zamora R, Navarro JL, Gallardo E, Hidalgo FJ 2006.** *J Agric Food Chem* 54:6101-6105. 10.1021/jf061239n.
   - Abstract (PubMed): epoxyalkenals convert amino acids into α-keto acids, with a pH 3 optimum. Phenylpyruvate exceeds PA at 37 C and decomposes at 60 C.
7. **Rizzi GP 2006.** *J Agric Food Chem* 54:1893-1897. 10.1021/jf052781z.
   - Abstract (PubMed): polyphenol quinones with Met or Phe at pH 7.17, 22 C, give 0.032-0.42 % molar yield of methional and PA.
8. **Ghiron AF et al. 1988.** *J Agric Food Chem* 36:677-680. 10.1021/jf00082a002.
   - **Abstract unverified.** The title reports 3-DG "involvement in a Strecker degradation".

Also wanted: the Supporting Information of Zamora 2015 (10.1021/acs.jafc.5b02619), which holds the rate series.
