# Pre-registration: the protein matrix as chemistry (2026-09-08)

## 1. What it is

A spec may state its protein loading (`protein_g_per_l`) and either name a matrix whose reactive
sites a dossier states (`data/species/protein_matrices.yml`) or give its own site densities
(`protein_sites`, mmol per gram: free thiol, disulfide, amine). Two things then happen, and nothing
happens without a stated loading:

1. The sulfur lane's protein-disulfide pool `PROT_SS` is charged (disulfide sites × g/L). The
   thiol-to-protein exchange channel, carried inert since it was declared, runs with a real pool at
   its bracketed ambient rate (no activation energy is measured; the rate is held).
2. Aldehydes and HMF are bound after integration by a pseudo-first-order factor over the thermal
   programme, with the rate brackets and activation-energy bands of the adduct dossiers
   (`src/kinetic_core/matrix_sites.py`): hexanal-class aldehydes and enals to the amine pool
   (Meynier 2004, Anantharamkrishnan 2020b, the Shepelev 2024 barrier), HMF to the thiol and amine
   pools (Hamzalıoğlu 2018). The bracket's corners price an interval width, the way the lipid lane
   prices its Q10.

Nothing is fitted. Every number is a count from a dossier over a molar mass, or a bracket with its
anchor. The only matrix on file is β-lactoglobulin; any other isolate must state its sites, because
no dossier on disk gives the free thiol and disulfide content of pea or soy isolates.

## 2. What counts as success, declared before the run

- **T1, a matrix that is absent changes nothing.** Every prediction without `protein_g_per_l` is
  byte-identical to before; the panel scorecard does not move (the freshness gate).
- **T2, the numbers reproduce from the dossiers.** β-lactoglobulin's site densities equal 15 / 1 / 2
  per 18 362 Da; the binding brackets are the ones printed in the adduct synthesis.
- **T3, the mechanism behaves.** The bound fraction rises with loading and with time, is zero at
  zero loading, and the bracket's lower corner never exceeds its upper. Hexanal in 1 % BLG at 20 °C
  for seven days binds a few percent, inside the synthesis's own "18 to 74 days" half-life bracket.
- **T4, the sulfur channel is honest about its weakness.** FFT in 1 % BLG at 145 °C for 20 min loses
  under 1 % to the protein disulfides at the held ambient rate, and the answer says the rate is held.
- **T5, refusal over defaults.** A named matrix without a loading is an error; a loading with a
  matrix that has no sites on file charges nothing and says why on the answer.

The hold-out the roadmap named (the four external matrix pots and the three isolate pots) cannot be
scored by this layer today: the bundles state the isolate but not its loading or its sites, and no
dossier on disk gives pea or soy site densities. That is recorded here rather than assumed.

## 3. What it will not do

It will not invent a site density for pea or soy, sample the brackets in the Monte-Carlo envelope
(the point interval carries them; the envelope's draw is a follow-up), or produce the electrophile
pool from browning (the thioether channel's partner stays a matrix input).

## 4. Outcome (2026-09-08, first run)

Shipped: `data/species/protein_matrices.yml` (β-lactoglobulin, from the two Anantharamkrishnan
dossiers), `src/kinetic_core/matrix_sites.py`, the `protein_g_per_l` and `protein_sites` spec fields
in the schema, the charge of `PROT_SS`, the post-integration binding factor with its interval, the
matrix block on the answer, `tests/unit/test_matrix_sites.py`.

- **T1 passed**: without a loading nothing is charged and every prediction is unchanged; the panel
  scorecard did not move (the freshness gate).
- **T2 passed**: the site densities are computed from the counts and the molar mass at run time,
  and a test asserts the file's rounded values agree with them; the brackets are the synthesis's.
- **T3 passed**, with the numbers stated: hexanal in 1 % BLG at 20 °C for seven days binds 3.8 %
  (bracket 1.2 to 11.6 %). The same hexanal at 160 °C for 30 minutes binds 0.1 % (0.03 to 0.5 %);
  2,4-decadienal at 100 °C for an hour 0.6 %; HMF at 130 °C for two hours 1.2 % (0.4 to 4.3 %).
- **T4 passed**: FFT in 1 % BLG at 145 °C for 20 minutes loses 0.003 % to the disulfide pool at the
  held ambient rate; the answer prints the pool and the source.
- **T5 passed**: a named matrix without a loading is a spec error at the front door; a loading with
  an unlisted matrix charges nothing and the answer says why.

**What the run says about the chemistry.** At the rates the adduct dossiers measured, covalent
binding of aldehydes and HMF to protein over cooking times is a percent-level effect, and over a week
of storage a few percent. It cannot explain the thousandfold misses on the pea and soy hexanal rows,
which is the adduct synthesis's own verdict ("does not matter at process temperature") reproduced
by the engine. Those rows are a storage and lipoxygenase question (roadmap, programme 5), not a
binding question. The layer's value is that a laboratory that states its isolate's sites gets the
binding declared, bracketed and printed, rather than absent.
