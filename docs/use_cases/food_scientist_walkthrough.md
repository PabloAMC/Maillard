# Formulating a Pea-Protein Burger: A Walkthrough for Food Scientists

This narrative outlines how a food scientist can use the Maillard framework to solve a common R&D challenge: improving the meatiness of a pea-protein isolate base while neutralizing its inherent "beany" off-notes.

---

## 1. The Scenario
You are developing a new plant-based burger patty using a commercial pea protein isolate (`pea_iso`). 
In previous sensory panels, your base recipe proved texturally adequate, but the flavor was described as "too grassy/beany," lacking the savory, roasted umami notes of real beef. 

**Your goals:**
1. Maximize `meaty` and `roasted` notes.
2. Minimize `beany` and `green` lipid oxidation notes (driven primarily by residual hexanal).
3. Identify which precursor (sugar/amino acid) combinations you should trial in the lab next week.

## 2. Setting Up the Pipeline
Instead of arbitrarily testing dozens of amino acids in the lab, you turn to the Maillard framework to screen combinations digitally. 

You know your extrusion/cooking parameters:
- **Temperature:** 130°C
- **Cooking Time:** 20 minutes
- **pH:** 6.5

You decide to test two formulation hypotheses:
**Hypothesis A:** Supplementing with Glucose + Glycine.
**Hypothesis B:** Supplementing with Ribose + Cysteine.

## 3. Running a "Dry Run"
First, check if your target conditions and ingredients are within the model's validated limits using the `--dry-run` flag.

```bash
python scripts/run_pipeline.py \
  --sugars ribose \
  --amino-acids cysteine \
  --lipids hexanal \
  --protein-type pea_iso \
  --temp 130 \
  --ph 6.5 \
  --dry-run
```

*Output:* The model warns you that the `pea_iso` matrix is currently treated with "Moderate" confidence because true quantitative benchmarks for matrix trapping are sparse. You acknowledge that your findings will be *directional priorities* (determining relative ranking), not absolute quantitative assertions (exact ppb values).

## 4. Predicting Hypothesis A (Glucose + Glycine)

Run the full prediction for Hypothesis A.

```bash
python scripts/run_pipeline.py \
  --sugars glucose \
  --amino-acids glycine \
  --lipids hexanal \
  --protein-type pea_iso \
  --temp 130 \
  --ph 6.5 \
  --target meaty \
  --minimize beany
```

**Observation:**
The terminal outputs a table. While some pyrazines form (contributing minor roasted notes), the system struggles to generate potent meaty aromatics. Glycine cannot provide the necessary sulfur backbone for key compounds like 2-methyl-3-furanthiol (MFT). Furthermore, glucose's high barrier to enolization means the reaction remains stunted under just 20 minutes at 130°C. Hexanal remains largely unreacted ("beany" penalty is high).

## 5. Predicting Hypothesis B (Ribose + Cysteine)

Now run the exact same conditions with Hypothesis B.

```bash
python scripts/run_pipeline.py \
  --sugars ribose \
  --amino-acids cysteine \
  --lipids hexanal \
  --protein-type pea_iso \
  --temp 130 \
  --ph 6.5 \
  --target meaty \
  --minimize beany \
  --report \
  --output-dir results/pea_iso_test_run
```

**Observation:**
The result is drastically different. 
- **Meaty notes:** Cysteine's sulfur group acts as the golden key. MFT and 2-Furfurylthiol (FFT) appear in the dominant desirable compounds list. 
- **Trapping Efficiency:** Ribose, a highly reactive pentose, rapidly forms Amadori intermediates. The amino acids effectively "trap" the residual hexanal into non-volatile Schiff bases, violently suppressing the `beany` penalty score.
- **Decision:** Hypothesis B is the clear winner. 

Adding the `--report` flag saved a cleanly formatted Markdown and JSON bundle to `results/pea_iso_test_run`, which you can attach to your experimental proposal ticket.

## 6. Going Further: The Optimizer

Instead of typing combinations manually, you can ask the framework to search for the optimal relative ratios of Ribose to Cysteine to Hexanal.

```bash
python scripts/optimize_formulation.py \
  --sugars ribose \
  --amino-acids cysteine \
  --lipids hexanal \
  --target-tag meaty \
  --minimize-tag beany \
  --protein-type pea_iso 
```

The optimizer sweeps through parameter variants and returns the formulation with the absolute best ratio of meaty aroma generation to beany suppression, constrained by the lysine budget. You now have your exact weighing instructions for your wet-lab experiments.

## 7. Comparing Head-to-Head

Before finalizing your lab order, you want to see exactly how Hypothesis A and B stack up against each other in a single report. You can use the `compare_formulations.py` tool for this.

Note: This tool compares formulations already defined in the system's "grid" or provided via standard names.

```bash
python scripts/compare_formulations.py \
  --names "Hypothesis A (Glucose+Glycine),Hypothesis B (Ribose+Cysteine)" \
  --ph 6.5 \
  --temp 130 \
  --target-tag meaty \
  --minimize-tag beany \
  --output-dir results/comparison_report
```

This generates a side-by-side Markdown table allowing you to compare:
- **Target Scores:** Direct comparison of meaty potential.
- **Risk Penalties:** Visualizing the beany suppression delta.
- **Trapping Efficiency:** Quantitative proof of which sugar is more effective at neutralizing off-notes.

## 8. Campaign Generation & Data Handoff

Finally, you need to share these results with your manager and the procurement team. Instead of sending screenshots of your terminal, you build a **Campaign Package**.

```bash
./scripts/docker_maillard.sh campaign \
  data/campaigns/pea_protein_screening.yml \
  results/pea_burger_campaign
```

This creates a polished, portable directory in `results/pea_burger_campaign/` containing:
- Individual reports for every trial.
- A campaign-level summary `campaign.md` with a leaderboard of the best-performing formulations.
- Full provenance metadata, ensuring that anyone who receives the folder knows exactly which version of the science was used to generate the results.

You zip this folder and attach it to your R&D report. You are now ready for the lab.

