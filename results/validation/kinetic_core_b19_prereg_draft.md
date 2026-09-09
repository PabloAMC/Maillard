# Pre-registration DRAFT: wave B19, amino-acid identity on the sugar path (written 2026-09-08, sources not yet read)

*A draft, not a pre-registration: it states the structure, the rows it will need and the tests, and
names the sources the rows must come from. It becomes a pre-registration when those sources are on
disk as dossiers and the rows are written with their numbers; nothing runs before that.*

## 1. Why

The trunk carries one amine, glycine. Six of the fourteen desirable odorants the engine cannot name
are Strecker aldehydes (methional, 3-methylbutanal, 2-methylbutanal, 2-methylpropanal,
phenylacetaldehyde) or their sulfur children (dimethyl disulfide, dimethyl trisulfide), and the
roasty pyrazines beyond 2,5-dimethylpyrazine need a Strecker aldehyde to add to the ring
(`tasks/roadmap_for_scientists.md` section 5b). B18 built the Strecker step for glycine on the small
dicarbonyls; this wave gives the step an amino-acid identity.

## 2. The structure

- Species: Leu, Ile, Val, Met, Phe, Ala, Pro as trunk reactants beside Gly; their Strecker
  aldehydes 3-methylbutanal, 2-methylbutanal, 2-methylpropanal, methional, phenylacetaldehyde,
  acetaldehyde; for proline the 1-pyrroline that makes 2-acetyl-1-pyrroline; methanethiol, dimethyl
  disulfide and dimethyl trisulfide from methional on the sulfur lane's oxidant pool.
- Steps: rule R07 per amino acid on glyoxal and methylglyoxal (aminoketone + aldehyde + CO2), the
  aldehyde-addition step to the aminoketone pair that makes the ethyl- and trimethyl-pyrazines,
  methional → methanethiol (retro-Michael), methanethiol oxidation to the disulfides with the same
  oxidant pool B17 named as limiting. Each step's constant is measured, fitted on a measured rate or
  within-study ratio, or the class is refused by name.
- The pH term of B18 is shared by all Strecker steps unless a source separates them.

## 3. The rows it needs (to be written from dossiers)

| quantity | candidate source | status |
|---|---|---|
| Strecker aldehyde yield or rate per amino acid at two or more temperatures | READ 2026-09-08 (`huang2017`, `hidalgo2013`, `zamora2015`, `weenen2001`, `parker2013`, `kocadagli2021`) and 2026-09-09 (`chan1994b`, `cremer2000`, `balagiannis2009`, `desclaux2006`, `deng2022`, `pan2025`): **no paper on disk prints a per-amino-acid second-order constant in water, and the four Parker 2013 named as the sources do not either.** Chan & Reineccius 1994 (RSC chapter) prints barriers only, whole-cascade and response-normalised: leucine aldehyde 80.3, phenylacetaldehyde 90.0, 2-acetyl-1-pyrroline 60.2 kJ/mol, no rate constant, glucose 1.25 M with four amino acids at 0.19 M each. Cremer & Eichner 2000: low-moisture (aw 0.52) barriers 124 (Leu), 120 (Ile), 115 (Val), 115 (Ala) ± 6 kJ/mol, rates figure-only, no volume to divide by. Balagiannis 2009 (liver extract): the Strecker step declared diffusion-controlled and not fitted; what carries identity is the yield fraction per amino acid, F(Ile)/F(Leu) = 1.6 ± 16 %. Desclaux 2006: two ARP constants at 100 °C (xylose + glycine, pH 6: ARP → 1-deoxyosone 2.79e-2, → 3-deoxyosone 5.50e-4 per minute), nothing on aldehydes or small dicarbonyls; the tabulated dicarbonyl courses Parker 2013 credits it with are in the Reading thesis, not in these four pages. Pan 2025: methional zero-order in a fruit-sugar pot, 0.268 mM methionine, unit inferred. Deng 2022: methional levels at 120 °C with response factor 1 | read; the rate rows do not exist in the literature on disk |
| amino-acid identity ratios in one pot | Amrani-Hemaimi 1995 Table 2 (40 isotope fractions; stranded since B2) | on disk, never used |
| methional → methanethiol, and the disulfides | to find (see the Scholar questions of 2026-09-08) | none |
| ethyl- and trimethyl-pyrazine formation from an aminoketone + aldehyde | Leahy 1989 distributions (on disk); Yu 2018 barriers (on disk) | on disk |
| 2-acetyl-1-pyrroline from proline + dicarbonyl | to find | none |

## 4. Tests, to be fixed with the rows

T1 the fit rows within 0.3 dex; T2 no scored panel row moves more than 0.05 dex; T3 the panel's
methional and 3-methylbutanal rows (other laboratories) within threefold; T4 Leahy's 95 °C
distribution including the ethyl- and dimethyl-pyrazines; T5 identification. Ship rule to be
declared when the rows are.

## 5. Verdict on pre-registration (2026-09-09)

This wave cannot be pre-registered in the form section 2 describes, because the rows section 3
asks for do not exist in the eleven papers now read, including the four that the literature itself
names as the kinetic sources. What the corpus holds is (i) apparent barriers of three kinds that
disagree by 40 kJ/mol for the same aldehyde (water 80, dry glass 124, liver extract 137 for the
glucose entry), (ii) within-study identity ratios: Balagiannis' yield fractions (Ile : Leu 1.6),
Kocadagli 2021's aldehyde ratios at fivefold loading, Amrani-Hemaimi 1995's isotope fractions
(stranded since B2), and (iii) two methionine time courses with hidden or inferred units.

The honest structure, if the wave is written, is therefore not a per-amino-acid rate fit. It is an
identity layer on the one Strecker step already fitted: glycine's two constants from the pyrazine
step (`FROZEN_B18`) anchor the rate and its barrier; each other amino acid enters as a partition
ratio of the same dicarbonyl pool, FIT on the within-study ratios in (ii) (the owner's rule: ratios
fit), with its barrier DECLARED equal to glycine's unless a source separates them (Cremer's
Leu − Ala = 9 ± 8 kJ/mol says no source does). Its ship rule would test the ratios it fitted (T1),
that no scored panel row moves (T2), and Deng 2022's methional level at 120 °C to an order of
magnitude (T3, the only external level). What it cannot claim is an absolute aldehyde level from a
sugar pot, for the reason the pyrazine step already carries as a caveat: the dicarbonyl supply in
water is not measured. That pre-registration is a separate document to write when the three ratio
sources have been re-read for the ratio rows; it is not written tonight, and nothing runs before it.

