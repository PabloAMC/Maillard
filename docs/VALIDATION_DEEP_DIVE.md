# Scientific Validation Deep Dive: Why "The Numbers" Matter

A common question in computational food science is: *How do we trust the numbers?* We use a three-layered validation approach to ensure our predictions are both **valid** (chemically accurate) and **useful** (industrially actionable).

## 1. Directional Validity (Correlation vs. Absolute Yields)
In complex systems like the Maillard reaction, absolute yields (e.g., exactly 12.44 ppm) are difficult to predict because of minor variations in trace moisture or metals. Instead, we validate for **Directional Abundance**.

- **The Test**: We run `scripts/compare_sim_to_lit.py`.
- **The Success Factor**: If the literature says `Furfural > 2,5-DMP` and our engine predicts `Furfural > 2,5-DMP` with a **Pearson R > 0.85**, the numbers are **Valid**. It means the engine's internal "Kinetics Map" accurately reflects the competitive branching of the chemistry.

## 2. Industrial "Bang for Buck" (Sensitivity Analysis)
The numbers are **Useful** if they identify the primary levers for a formulation. We use **Sensitivity Coefficients** ($S$) to measure this.

$$S = \frac{\% \Delta \text{ Flavor Score}}{\% \Delta \text{ Precursor Conc}}$$

- If $S_{Cysteine} = 5.0$ and $S_{Glucose} = 1.2$, the scientist knows that 1g of Cysteine is 4x more effective than 1g of Glucose for reaching the target. 
- **Validation**: We verify these coefficients against literature "dosage-response" curves.

## 3. Psychophysical Reality (The Human Meter)
A 10x increase in concentration **does not** result in a 10x increase in smell. This is why our "Numbers" for flavor scores use **Stevens' Power Law**.

- **Logic**: Perceived Intensity $\propto Concentration^{0.5}$.
- **Usefulness**: This prevents the tool from suggesting unrealistic "over-supplying" of ingredients. It correctly predicts the "Diminishing Returns" region of the formulation curve, saving ingredient costs.

## 4. Physical Sanity Checks (Regression)
Every code change is validated against a set of **Physical Invariants**:
1. **Arrhenius Check**: Does $T_2 > T_1 \implies Rate_2 > Rate_1$?
2. **pH Branching**: Does $pH > 7 \implies Pyrazines \uparrow$? 
3. **Mass Conservation**: Are we losing atoms? All SMIRKS rules are stoichiometrically balanced.

## 5. Hierarchy of Fidelity
We don't just guess. We use a funnel of increasing accuracy:
- **Tier 0**: Fast flux (Heuristic).
- **Tier 1**: ML (MACE-OFF24) - Near-DFT barriers.
- **Tier 2**: Full DFT - Research-grade validation.

By using the Bayesian Optimizer in the tool, you are effectively performing thousands of these "Sanity Checks" in parallel to find the mathematically most robust solution.
