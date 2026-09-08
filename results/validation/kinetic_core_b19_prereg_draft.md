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
| Strecker aldehyde yield or rate per amino acid at two or more temperatures | READ 2026-09-08 (`huang2017`, `hidalgo2013`, `zamora2015`, `weenen2001`, `parker2013`, `kocadagli2021` dossiers): none gives a per-amino-acid second-order constant in water. Huang 2017 has the leucine and isoleucine ladder (90 to 130 °C, k1 for amino-acid loss and k2 for the aldehyde with barriers 83 to 121 kJ/mol) but prints no concentrations, so nothing converts, and its Arrhenius lines miss its own table by fivefold. Hidalgo 2013: one barrier (38 kJ/mol) on the glyoxylic-acid route, a transamination, not the dicarbonyl step. Zamora 2015: the aldehyde-to-amine split on lipid carbonyls, low-moisture slurry. Weenen 2001: yields at one temperature. Still to fetch, named by Parker 2013: Chan & Reineccius 1994 (the OTHER 1994 chapter: 3-methylbutanal and phenylacetaldehyde pseudo-zero-order rates and barriers at pH 6 to 8, 75 to 115 °C), Cremer & Eichner 2000 (barriers 115 to 124), Balagiannis 2009 JAFC (leucine and isoleucine at 120 to 140 °C), Desclaux 2006 (glyoxal, methylglyoxal, butanedione time courses 80 to 120 °C); and Deng 2022 (methional against time at 100, 120, 130 °C, on the download list) | partly read; the rate rows are not yet on disk |
| amino-acid identity ratios in one pot | Amrani-Hemaimi 1995 Table 2 (40 isotope fractions; stranded since B2) | on disk, never used |
| methional → methanethiol, and the disulfides | to find (see the Scholar questions of 2026-09-08) | none |
| ethyl- and trimethyl-pyrazine formation from an aminoketone + aldehyde | Leahy 1989 distributions (on disk); Yu 2018 barriers (on disk) | on disk |
| 2-acetyl-1-pyrroline from proline + dicarbonyl | to find | none |

## 4. Tests, to be fixed with the rows

T1 the fit rows within 0.3 dex; T2 no scored panel row moves more than 0.05 dex; T3 the panel's
methional and 3-methylbutanal rows (other laboratories) within threefold; T4 Leahy's 95 °C
distribution including the ethyl- and dimethyl-pyrazines; T5 identification. Ship rule to be
declared when the rows are.
