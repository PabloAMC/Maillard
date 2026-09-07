# What the field knows about the Maillard reaction, and what it does not

*A short orientation for a reader who knows what the Maillard reaction is and nothing else. This
page is about the published literature, not about this repository's model; the companion guide,
[What this model can and cannot predict](STATE_OF_THE_MODEL.md), is about the model. Updated
2026-09-07; the figures are generated from the papers the repository has read.*

## The chemistry, as the field draws it

![The Maillard reaction as the field draws it, and how well each part has been measured](../assets/thiol_sink/14_field_scheme.png)

The scheme has been stable since Hodge drew it in 1953. A reducing sugar and an amino group
condense to an Amadori compound. That compound decays along two branches: 1,2-enolisation to the
3-deoxyosone, the parent of HMF and furfural, and 2,3-enolisation to the 1-deoxyosone, the parent of
the caramel furanones. Both branches also fragment into small dicarbonyls (glyoxal, methylglyoxal,
diacetyl), which react with amino acids in the Strecker degradation to give the amino-acid-specific
aldehydes and, from those, pyrazines. Everything reactive eventually condenses into melanoidins, the
brown polymers. Cysteine adds a sulfur branch: it releases hydrogen sulfide, which adds to the
furanones and furfural to make the two thiols that smell of cooked meat. Asparagine has its own branch
to acrylamide. Fats oxidise alongside and their aldehydes cross into the Maillard chemistry.

The colours in the figure are the point of this page. Dark blue steps have published rate constants at
several temperatures; green steps a rate or a yield at one temperature; grey steps are known as
mechanisms from isotope-labelling work but have no rate; the dashed red step has no published rate at
all. Read that way, the sugar and acrylamide branches are quantified, the thiol formation branch is
measured once, and the removal of thiols, the step every meaty-aroma prediction ends on, has never
been measured in a cooking pot.

## Where the measurements sit

![Where the quantitative measurements sit](../assets/thiol_sink/15_field_coverage.png)

Each bar is one study and the temperatures it measured at. Three things stand out. Almost every
quantitative study sits between 100 and 145 °C, in water; below 100 °C only milk-powder and storage
studies exist, above 150 °C only dry glasses and roasting matrices. The thiol branch is measured at
single temperatures by different laboratories, so no laboratory has a temperature series of its own
except the five-rung grids published as figures. And the matrix changes the answer: the same
dicarbonyls come out in the opposite order in a sugar glass at 180 °C and in water at 120 °C.

## Settled, quantified, open

**Settled by mechanism.** The intermediates and their connections. Labelled-sugar experiments
(Cerny and Davidek in Lausanne, Schieberle's group in Munich, the Beijing carbon-module work) fixed
which carbons of the sugar end up in which product: the meaty thiol MFT keeps the intact five-carbon
skeleton of ribose and goes through norfuraneol; FFT goes through furfural; the mercaptoketones come
from an intact chain. Nobody disputes the map.

**Quantified.** The sugar side, in water at 100 to 120 °C, by the Wageningen multiresponse models
(Martins and van Boekel; Brands and van Boekel): every step of the glucose-glycine cascade with a
rate and a temperature dependence. Acrylamide, by the Leuven group (De Vleeschouwer and Hendrickx)
and Knol: formation, elimination, and their dependence on temperature, pH and water activity, from
one laboratory over 120 to 200 °C. Thiol formation, once: Hofmann and Schieberle's 1998 feeding
experiments give the yield of each thiol from each fed intermediate at 145 °C, and Whitfield and
Mottram the same for norfuraneol at 140 °C at two pH values.

**Open.** How thiols are removed once formed, at what rate and with what temperature dependence;
the field has identified the products (disulfides, adducts) and measured loss only in a coffee brew
at 80 °C and in storage at 50 °C. What oxygen does: one measurement exists, for the Amadori
compound's Strecker aldehyde (nine times higher under air than argon), and none for the thiols.
Whether constants transfer between matrices: the evidence says they do not, and no study has yet
measured the same step in water and in a glass. Anything below 100 °C for aroma, and anything in a
real food matrix rather than a model pot.

## Ten papers to read first

1. **Hodge 1953**, J. Agric. Food Chem. 1:928. The scheme everyone still draws.
2. **Martins & van Boekel 2005**, Food Chem. 90:257. The glucose-glycine cascade with rate constants at three temperatures; the template for every multiresponse model since.
3. **Hofmann & Schieberle 1998**, JAFC 46:235. Each intermediate fed on its own, the yield of MFT and FFT from each; the quantitative backbone of meaty-aroma chemistry.
4. **Cerny & Davidek 2003**, JAFC 51:2714. Labelled ribose with cysteine: which carbons become which sulfur compound.
5. **Whitfield & Mottram 1999 and 2001**, JAFC 47:1626 and 49:816. Norfuraneol fed with cysteine or H2S at pH 4.5 and 6.5: the pH switch in thiol formation.
6. **De Vleeschouwer et al. 2009**, Food Chem. 114:116. Acrylamide formation and elimination by multiresponse modelling; the best-parameterised branch in the field.
7. **Kocadağlı & Gökmen 2016**, JAFC 64:6333. Dicarbonyls and HMF in a glucose glass at 160 to 200 °C; the dry-side counterpart of the Wageningen work.
8. **Schieberle, Hofmann & Münch 2000**, ACS Symp. Ser. 756 ch. 10. The one published time series of the thiols at 100 °C.
9. **Hofmann & Schieberle 2000**, JAFC 48:4301. The Amadori compound's oxidative Strecker route: the only measured oxygen effect.
10. **Yiltirak et al. 2026**, Food Res. Int. The most recent thiol ladder, with the vessel and buffer stated, from an independent laboratory.

Nine of the ten are on disk in `data/articles` (Hodge 1953 is cited for its scheme); six have an
extraction dossier in `data/lit/extraction_dossiers`, and Cerny 2003, Whitfield 1999 and
De Vleeschouwer 2009 are read into the model's constants directly, without a dossier of their own.
