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
| Strecker aldehyde yield or rate per amino acid at two or more temperatures | Cremer & Eichner 2000 (cited by Balagiannis 2015 with Ea 115-124 kJ/mol); Hofmann & Schieberle 2000b (ARP-Phe at 100 °C, on disk); Chan & Reineccius 1994 (on disk) | to read / to re-read for rows |
| amino-acid identity ratios in one pot | Amrani-Hemaimi 1995 Table 2 (40 isotope fractions; stranded since B2) | on disk, never used |
| methional → methanethiol, and the disulfides | to find (see the Scholar questions of 2026-09-08) | none |
| ethyl- and trimethyl-pyrazine formation from an aminoketone + aldehyde | Leahy 1989 distributions (on disk); Yu 2018 barriers (on disk) | on disk |
| 2-acetyl-1-pyrroline from proline + dicarbonyl | to find | none |

## 4. Tests, to be fixed with the rows

T1 the fit rows within 0.3 dex; T2 no scored panel row moves more than 0.05 dex; T3 the panel's
methional and 3-methylbutanal rows (other laboratories) within threefold; T4 Leahy's 95 °C
distribution including the ethyl- and dimethyl-pyrazines; T5 identification. Ship rule to be
declared when the rows are.
