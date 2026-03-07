# **State-of-the-Art Computational Architecture for Maillard Reaction Modeling (2026)**

## **Executive Summary**

This document outlines the algorithmic and hardware architecture required to model the Maillard reaction network for plant-based meat alternatives. The model balances proven industrial production methodologies with cutting-edge machine learning approaches, achieving Density Functional Theory (DFT) accuracy without the traditional computational bottlenecks. Furthermore, it integrates the complete biochemical context of plant matrices, including competitive pathways, exogenous precursor mechanics, and inverse design strategies.

## **1\. The Core Paradigm: Algorithmic Decoupling**

The guiding principle of this architecture is the strict decoupling of geometric optimization from the calculation of electronic structure. The modern, validated production workflow is defined as follows:

1. **Reaction Generation:** Network mapping (e.g., smirks\_engine.py).  
2. **Initial Guess:** Pre-optimization and conformational search using semi-empirical methods.  
3. **Transition State (TS) Search:** Rigorous geometric optimization using Machine Learning Potentials (MLPs) and eigenvector following.  
4. **Electronic Refinement:** Single-Point energy calculations with high-fidelity DFT.  
5. **Microkinetic Simulation:** Integration of rate constants to predict concentrations over time.

## **2\. Machine Learning Potentials (MLPs) and Geometric Optimization**

To reduce optimization times from hours to seconds, MLPs replace DFT in exploring the potential energy surface.

**SOTA Production Recommendation: MACE (Multi-ACE) with Fine-Tuning**

MACE is currently the most robust equivariant architecture for condensed-phase chemistry. However, applying a pre-trained foundation model directly to the Maillard reaction will produce severe errors.

* **The Limitation:** Standard training datasets lack sufficient representation of sulfur chemistry, sugar fragmentation, zwitterionic intermediates, and complex proton transfers.  
* **The Production Requirement:** The use of MACE is only state-of-the-art *if* a fine-tuning loop is implemented using specific DFT data derived from Maillard precursors (e.g., cysteine, ribose).

## **3\. Transition State Search and Optimization**

**The Problem:** Blind exploration of transition states (TS) is the largest computational resource sink.

**SOTA Production Recommendation: xTB** ![][image1] **Sella \+ MACE**

The most efficient and industrially validated approach requires a high-quality initial guess followed by direct mathematical optimization.

1. **TS Guess via GFN2-xTB:** Despite its absolute energetic errors (5–12 kcal/mol), xTB remains the unparalleled tool for rapid conformational sampling and generating plausible initial TS geometries.  
2. **Optimization via Sella (Eigenvector-Following):** Using the TS guessed by xTB, the Sella optimizer (guided by MACE forces) walks directly toward the saddle point. This method is drastically faster than optimizing entire reaction paths.  
3. **Fallback: NEB.** If the xTB guess is poor and Sella fails, the Nudged Elastic Band (NEB) method is maintained as a robust fallback alternative, as it does not require a highly accurate initial guess.

**Alternative Production Pathway: Generative Diffusion Models**

While xTB provides a first-principles conformational search, generative diffusion models can predict TS geometries directly from 2D molecular graphs, bypassing quantum mechanical sampling entirely.

* **SOTA Diffusion Recommendation: React-TS**  
  Among the available diffusion architectures (e.g., TSDiff, DiffTS), **React-TS** stands as the 2026 state-of-the-art choice. It utilizes SE(3)-equivariant stochastic diffusion to generate 3D saddle points with sub-angstrom accuracy relative to DFT. It should be deployed as the primary first-pass TS generator for well-characterized reaction classes (e.g., standard proton transfers and peptide scissions) to maximize throughput. If React-TS yields a low-confidence geometry due to out-of-distribution complex sugar rearrangements, the pipeline must automatically fall back to the xTB conformational search.

## **4\. Electronic Structure Refinement (Tier 2\)**

**SOTA Recommendation:** ![][image2]**B97M-V / def2-TZVPP via Microsoft Skala**

Once Sella and MACE converge on the geometry, exactly one Single-Point energy calculation is performed using DFT. The range-separated hybrid functional ![][image2]B97M-V is widely considered one of the best for predicting thermochemical barriers and non-covalent interactions, systematically overcoming the inherent errors of MLPs and semi-empirical methods.

* **Infrastructure Integration:** This theoretical standard must be seamlessly integrated with your existing compute backend (src/skala\_refiner.py). The Microsoft Skala AI-accelerated DFT engine should be explicitly configured to evaluate the ![][image2]B97M-V functional. This harmonizes state-of-the-art functional accuracy with your specific high-performance hardware acceleration.

## **5\. The Missing Requirement: Explicit Solvation**

**The Problem:** Maillard chemistry occurs in aqueous matrices with high ionic strength. Implicit solvation models (PCM/SMD) systematically fail when modeling proton transfers, dehydration reactions (crucial in Amadori degradation), and sugar rearrangements.

**SOTA Recommendation: Explicit Solvent Modeling via CREST (QCG)**

The calculation of rate-limiting barriers must include **explicit water clusters (3 to 6** ![][image3] **molecules)** around the reactive center. These act as catalytic bridges for proton exchange. Ignoring explicit water will result in artificially high activation barriers (often deviating by ![][image4] kcal/mol).

* **Automation Tooling:** Manually predicting the lowest-energy hydrogen-bonding network around a reactive center is statistically improbable. The pipeline must utilize the **Conformer-Rotamer Ensemble Sampling Tool (CREST)**, specifically its **Quantum Cluster Growth (QCG)** algorithm. By inputting the bare reactant geometry and specifying the solvent count, CREST computationally automates the geometric placement of water molecules and performs a rapid conformational search prior to the transition state optimization.

## **6\. Thermodynamic Precision: High-Temperature Entropy**

**The Problem:** At extrusion temperatures (120°C to 177°C), flexible organic molecules and hydrogen bonds exhibit low-frequency vibrational modes. Treating these as harmonic oscillators introduces massive entropic errors (easily ![][image5] kcal/mol).

**SOTA Recommendation: Grimme's Quasi-Harmonic Correction**

It is imperative to apply Grimme's free-rotor interpolation (qRRHO module) to the calculated vibrational frequencies. This correction prevents free energy barriers (![][image6]) from artificially inflating due to purely mathematical errors in the standard harmonic model.

## **7\. Network-Level Scaling: ![][image7]\-Machine Learning (![][image7]\-ML)**

To predict the exact DFT energy across the entire network of thousands of generated reactions without incurring massive supercomputing costs:

1. Calculate the base energy using MACE (![][image8]).  
2. Train a parametric correction model on a highly diverse subset of **500 to 1000 DFT-calculated Maillard reactions** (covering sulfur, dehydration, and scission). A subset of 200 reactions is insufficient to capture the variance of this network.  
3. Apply the correction: ![][image9].

## **8\. The Ultimate Goal: Microkinetic Modeling**

Even with perfect DFT barriers, predicting the final aroma yield is impossible based solely on thermodynamics. Stoichiometric competition, concentration changes over time, and temperature ramps must be modeled.

**SOTA Production Recommendation: Integration of Kinetic Solvers (Cantera)**

Thermochemical barriers (![][image6]) must be converted into temperature-dependent rate constants (Arrhenius/Eyring equation) and integrated into an ordinary differential equation (ODE) solver like Cantera. This allows for simulating the concentration profiles of target volatile species versus time, which is the only metric directly comparable with analytical GC-MS data.

## **9\. Evaluated Options and Their Positioning**

To maintain rigor, the following industry-standard methodologies are classified according to their current role:

* **GFN2-xTB for Final Barriers:** *Inadequate.* Its error of 5-12 kcal/mol is too high for exponential rate constants, though it remains state-of-the-art for pre-optimization.  
* **Legacy Functionals (M06-2X, B3LYP-D4):** *Acceptable but suboptimal.* M06-2X is still widely used and respected, and B3LYP-D4 is functional. They are not "invalid" choices, but ![][image2]B97M-V offers systematically superior benchmarked thermochemical performance.  
* **TS Diffusion Models (TSDiff):** *Frontier research.* Highly promising, but published success rates are heavily biased toward pharmaceutical datasets. They require exhaustive validation before replacing xTB/Sella in the production pipeline.

## **10\. Reaction Network Generation and Environmental Conditions**

**The Problem:** Determining which reaction pathways exist before calculating their energies. The Maillard network comprises thousands of competitive reactions. Omitting key pathways invalidates final concentration predictions.

**SOTA Recommendation: Strict Rule-Based Generative Engines**

It is mandatory to use a robust generative engine (e.g., smirks\_engine.py) based on chemoinformatic transformations (SMARTS/SMIRKS). This system must:

* Guarantee strict mass balance at every node in the network.  
* Feature native support for modeling critical environmental variables in food systems, specifically **water activity (![][image10])**, **pH**, and **heme iron catalysis**—factors that drastically modulate selectivity toward pyrazine formation.

## **11\. Chemical Crosstalk and Competitive Pathways in Plant Matrices**

**The Problem:** Unlike animal tissue, plant proteins (pea, soy, gluten) lack the native precursors (ribose, cysteine) to generate desirable meaty volatiles (e.g., 2-Methyl-3-furanthiol \[MFT\], 2-furfurylthiol \[FFT\]). Instead, they generate lipid oxidation off-flavors (e.g., hexanal, nonanal) and compete thermally for available amino acids.

**SOTA Recommendation: Inclusion of Matrix Formulation Metrics**

The computational model cannot treat the Maillard reaction as an isolated system. It must explicitly calculate:

* **Exogenous Sulfur Precursor Mechanics:** To overcome the plant-based sulfur deficit, the generative network (smirks\_engine.py) must explicitly model the thermal degradation of synthetic additives. This strictly includes the thermal fragmentation of **thiamine** (yielding ![][image11] and 2-methylthiophene) and the precise cleavage of **glutathione** (yielding glutamic acid and reactive cysteinyl glycine dipeptides).  
* **Lipid Trapping Efficiency:** Quantification of condensation reactions where added amino acids bind to the electrophilic carbonyls of off-flavor aldehydes (like hexanal), forming stable Schiff bases and neutralizing their sensory impact.  
* **The Lysine Budget and the DHA Pathway:** The dehydroalanine (DHA) pathway occurs under severe heat and shear, where serine and cysteine undergo ![][image12]\-elimination. DHA reacts with lysine to form lysinoalanine (LAL). Because this pathway consumes lysine—a primary Maillard precursor—omitting it will result in a systemic overestimation of aroma generation.

## **12\. Inverse Formulation Design**

**The Problem:** Evaluating individual reactions in a bottom-up manner is prescriptively inefficient for product development. Food scientists need to know the optimal initial matrix, not just the outputs of a given matrix.

**SOTA Recommendation: Automated Computational Search**

Implement inverse design modules (inverse\_design.py) that iterate over formulation matrices (initial precursor concentrations). The mathematical objective of this search is to maximize a reward function based on the production of target aromatic compounds (sulfur heterocycles, pyrazines) while simultaneously minimizing toxic advanced Maillard by-products (e.g., 5-hydroxymethylfurfural \[HMF\], LAL).

## **13\. Topological Validation (Intrinsic Reaction Coordinate \- IRC)**

**The Problem:** After locating a transition state (TS) structure with Sella or NEB, there is no inherent geometric guarantee that the saddle point accurately connects the reactant and product predicted by the SMIRKS engine. A false TS will derail the microkinetic model.

**SOTA Recommendation: Automated IRC Protocol**

The computational framework must enforce a topology validation routine. In high-throughput systems, this is implemented as a "Displace+Optimize" algorithm: the TS geometry is perturbed along the vector of the imaginary frequency in both directions (forward and backward), proceeding to optimize until local minima are reached. The resulting structures must be algorithmically compared (e.g., via RMSD calculation) with the original reactant and product structures to confirm the reactive pathway without manual intervention.

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABMAAAAXCAYAAADpwXTaAAAAdUlEQVR4XmNgGAWjgHpAXl5+L7oY2QBo2D90MbKBnJycDRCXoYuTDYCuO6egoGCOLs4gKytrQg4GGnYLaOg+dMP8yMFAg66BMNAIFhQDSQVAV00EGuSNLk4yABqiCDSsE12cLAA07BO6GNkAaNhhdLFRMNwAADgJIGwPRW62AAAAAElFTkSuQmCC>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA0AAAAYCAYAAAAh8HdUAAAAr0lEQVR4XmNgGAVDCsjLyyehi8GAtLS0DFC+GEVQQUFhFxBnoAgiAaCG/yCMLvgPRQANKCkpyaFokpOTKwMKaCKpARnyQ0tLiw1NDKEJyPmLJMegoqIiiuEUBrC6vcgcFAVAfhS6GNC/HUAxI7gASIGxsTErEv8tEL+E8RUVFcXRDQEp+g4SBOJnIBrox0dQPgj/AtFAjeoomoBWC8AUAdkzoQalIGm0RNEwCgYCAAD69DKk6qv3SQAAAABJRU5ErkJggg==>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACcAAAAYCAYAAAB5j+RNAAACHUlEQVR4Xu1VO0gcURTVCJEYY0hgWdjf7G4WVrYRW8HGToWkEwISJFoIdgqClWAgbWxVtBMrLYSkUbHVJpBSkPQa/BADMYWfc5Y78jyzMyMoSGAPXN6bc8697857s28bGuqo42Hhed4wYiGfz79XLQAYlxC/EVdOHJjciPmhaMeIvltFYpBIJFqRc4H4msvlXpHDOM16GAfVH4C/uPIE+G1qpVKpTbU4IG+euWxQNexeu63bqZoL7hBN31UgohqPAhafC2vMR2xtiGM0ZLPZt6oRVuCf8lHgLlveimouWDeuuYMwAxrutkVmVIsC/JdhNV1Y7XCfb8AxvMPYj+hD9DLwwe5Sy2QyzzQvDMVi8aXV/Kua4q7N7SEm3EBjk7HJNQD/InPwslOqKSLre/a9eSHXg2mB7w3cqmn7eGwSrbogPok3Lq/Ay38274ZqVUD4RYPyBPguS/4k/JYz/6j5eD5XrhasNn2NqlXhGAIAv0EtlUq1CM+cHfe5XC6/cJ6XyaXT6YzPKaCPW50J1Xw0meE+99sT9bApy112eR+FQiFp+s0JBIAz/0ITxgHV+Ou8S3PQf/KyrcFvMlcvYHA9Vrf2jkFYQZwhThBHiFPEBbVkMvkc8z/GUaPnnNeM1sFLfQA/q7wP6tYIa/3gHP519T04cDQdWGiU80ql8hT/CM3qeRTgyF9jF9bQ3BDGEcQ39Twa7KhuhXrqqON/wTWhxM0OmxofRQAAAABJRU5ErkJggg==>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAXCAYAAAB50g0VAAABkElEQVR4Xu2VvUrEQBSFk0VFwTZgkX/TGBCFgKWF+BK+gbCFiJWdhb6CnSiID6AgKgo2ljZ2NqIgi53CChZxYT3XnYTxkp3NxAVB8sFhmXvuzZzNr2HU1Px/RjzPS3mR8H3/GN6u4zjTWJpBEMy5rnuO+jzvHTrY+B7qZuI+gfqN3CN0xPsMJF/gtWGBDduKgNc4W+v4PYQ2uJ8Tx/EYGp7oH2Fpcv83DAh4ZWju18DQHfRo2/YEN6swIOCloRkwB8NndHBcginu6aAKiAfiAlqDn+L3gPqw3yrvU4LBfQx2cEZnuVcGVUCEOYG2pPU49WLPZamtHBjcoWE8UIvcU6EKWAT16vTnYKgpLsEK91SoAiZJMspr2gHRvC2GlrhXhn4BwzB0xXFP5XrpgGjagz6hGe7p0C8g7rOQ6vx+EwE7cu0H9KlBw1sURRb3qoBjfRQFJES9ka1x+2xSjR4Wqe0bE8Yt9FBkVgHHeodeoGehFvQqvrsZ9N5Nvd6V6pIsy5qU/B7iU1fthVlTU/O3fAF6kH6HKfjafQAAAABJRU5ErkJggg==>

[image5]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAXCAYAAAAcP/9qAAABVUlEQVR4Xu2Uu0rEUBCGd1W8gCAWASExl00KDdhaWmjjA9go6Av4DHaKb2AnvoX6AGplYyGCKIpgvSKIgoJ+I2eXs4OJ2RiwyQfDcOafP5Nz9mwajZqakgRB8OT7/jox2Wq1Jsgr1Nq6r3IY8qnDdV2vpymKovmeQgXIIHa5S94jL2n9mzRNh2m4J05YNrVeBhmsa3kMYLgg7jzPG9NiP/Q7uAvGI+I5DMMprRVBBuO9Jl8SZ8QH5SHdlwm/z4GYOIE5reUhg5MkGbHWh6VOAdOOGLmIC1orArufET+xpbVcMGya41vVWgaDai33RwZfqfrP0LhtDItay4LeG/OSo52a4zjj5jnyz8mGhn3inZjV2m9wLx7wvdg1XmLZDF6z610wHSO2uRiO1ooSx/E0z7i1a6zfiFe7JjQpnkuzfTx/gQ1smB0+mnyqezqfzEq+WDU1/84X6gNVachJshYAAAAASUVORK5CYII=>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACYAAAAYCAYAAACWTY9zAAAB9UlEQVR4Xu2VvUsDQRDF/Rb9B5SY5JIYRAKKkMoPEC1EsLOwsBVEED8KK8FCAiJaWIi1WCgWUYS0thZW2ohWWllpoSKEoKhvwpzZvFzuEhBMkR8Me/vmzc7u5pLU1FQAPp+vlbV/JxaLNVmWlWa9IggGg2esuYKTZPx+fwvrpRIKhfrQdAOxE41Gm20d6yZ0/EZ8Ip7kGf7uXHURYFzVwkfOeYGaB61NYVNDiEE83yK28HyM8dr0Y0MH5twVXTgb8h5w3gncbofW3HBOgP6s+UlDrhXdmBcHJ1iGOYFY04Xu2cPIR67eO87ZwNMjHtbRb501R8xibSbzesNSgOFzxcFTB+2FtEJgmkFs23OcZlObOn48AnJJ8eD9meUcA9+Vg3bEWgEOJ/K8Da98Mew6xIeM4XC4lz1ZcOIpGPZYh7Yrhchfci4QCHTq4l+cKxXUJlnLw+3U9ukc9BXd9D7n4vF4I3ITiFG8EsMYR3ArY2RrgP5OWg4Ujrv9nkhj3dw56Uuqz5u6EIlEupBfRO5QPWnMF9iHvu2s/SKFrDG6eJ4PN9Cm+oWpm6DxnHqmOecKCgYQJ6wz8Jxqgzyv04ZNLP3LYd0Te+FywqzHF8Cn+qupC8ZtlbcxvAdBblpKyB8zLSU/lG+GJ6NjvyQxpshfpUqVKn/BDzVMv8f4ke7sAAAAAElFTkSuQmCC>

[image7]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAdCAYAAACqhkzFAAAAx0lEQVR4XmNgGAXIQF5e/qeMjAwnujhZAGhYFRD/B+Kn6HJkAahhYKylpcWGLk8SUFBQKAAa1AzEtVBD76GrIQmADEFmQ/nMSEqIB0DNyUDcDeMDXdsBNfQqsjqiAbLrkMWwiRMEcnJyYUCNU9HFgWKTQQYC5U+iy+EF+FxBsiuBYeUBxAvRxWEA6Lr5UEP3oMthBcTYTrQrgYqsgHgtujg6AKpZBzUUv1qYzaRgdDPgQElJSQ5dMTEYGKb96GaNglEwCoYPAAAREVhS5Qkp0AAAAABJRU5ErkJggg==>

[image8]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAD0AAAAYCAYAAABJA/VsAAAC7UlEQVR4Xu2WX4jMURTHZzO0RSlMs83Mzm9mmjxMHrCFTSklaf2X4kl5ELUv0oa8SCLeyJaUFw8brSJ5kQflwYPwsB6U8u+BIlFCu/YPPmfn3NmzZ3YmtZvaab51uvd+z/fcc8+de+9vYrEmmmiiYZHL5e5HUfTnX8zHznrUK6y9vT1VyzebMUeLHvCOgIYrmoJ6pCh+0W2W5+gfNpqGK/qbL4pxVzab3RfGbMBq65/1CPcZ68TWU+xRvwmNhnCf31LsWdpe7P3/KJocKz3nwQlro4l7flqg0GNSIJNvsjzc69BPp9OLaVrCGO0p/C9rbYzEio+5L8dMnIdoCoXCQs8LyPEc/03ZmHw+v4r2GTYW/PT7sNGofDWvo+8n30OtZb+dqwqIvvvFd3R0zIXrMppf1q/cZ4lLpVJLHL+CpBdkkZb3YIE3JJ52rXPFdeGtloTbjA07TnTbLccGbcxkMkXLVUECfdEWxWIxgf+B54XDPpJ0j+X1FIxFZtOmgFypbl30IetQbo3lAvCdcePKutm8rdLyBVo2oZgCpVJpnhb9xPsCdOKqO8WOrpNfiwVeCRzaW9rW3EQB/i/aDmK9hj9dJ7ZF1hsG5N5gtXXiJkOPYdX3WX07dEOqjjYJD6jmIP5X0k8kEgvk7nO0FtVbADHL0aSlj25A7mHwab67E+raQPdIchN/Mlf+K93nNZOAqD8qPwK/NZE34UewQVmkj4f/IC2P0FLRS5/k96RFf5H+Y6u3ED2aO7S3sU/YV+vDdll9ADHn7FjnqTy+5Cwo/64imkn4X4fkx81YNmzK+4zuvB2jOxE2TcdS9BGrERDX5h88G2cQR7vXk9MGk7bKfQpjXWinHYe+RTKZnI/vkuUY77R65t3i4yUf3FPLyaZ6nb5Po5abEUTlb6J84oZy+mrTf6HtVewn9gMbcXFDGif8+Hc7V75ioherfH8pfLcUJNeF9g3WXZkoNj7XsPidyekaZc5rVttEE000Lv4CAR/+dihdwF4AAAAASUVORK5CYII=>

[image9]: <data:image/png;base64,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>

[image10]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABUAAAAYCAYAAAAVibZIAAABD0lEQVR4XmNgGAWjYOgBZgUFhUo5Obl8dAmyANCgdnl5+V8gtri4ODeIDcT/0dURDYCao0EGAF3JgSR2llJD/wPxcyxiX5DFiAZAb4dAXZmOLA4Vq0QWIxoADd2B7k0lJSU5kJixsTErsjjRAKh5CrqhQP4SmBiQXgqitbS02IDss0DX30JS9xlKHwWKT4CJw2IabqiioqIbNDxhhsLot8h8KSkpESD7L4gNNNAAyF4MMwMMgEHgDDMIiLNBYkD6H4gvIyMjBFMH1NwoD3U5kD0RyO5DMsMHxiYJgCwBupALyv4LsxDIDkZVSQKAeR3GhkUkkH0OoYpEANRsCQ2Wt8AUogakvwDxPXR1o4D6AADNlk4K/v0RTwAAAABJRU5ErkJggg==>

[image11]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACUAAAAYCAYAAAB9ejRwAAACJElEQVR4Xu2VOYgUQRiFV9f7CoRxYK6eCwZGwWAQQTQwdBXNDEXEQNlEEAXBSMFUs2UTQw1EQxMRUxU0MRJEA6P1RjfwwF2/f/irKd/0bI8oazIPHtX93vv/qupzYmKMMf4c9Xr9TJIkFxuNRhK0arW6K870Qeg6/AwXI865vYLjN+J9gFO/NclBrVaboWYB7m+32wXGG/CZ9cNepfkUYVLVDegPzKPhFvXywIKuUPtD9dBT9Rh2RWxRT9UwLLXgPFhd1i3iVu7Gm1c9Bea0Fx9Wz+CL+qZ6HrhKPatl3K4ec+3AO6d6Csy5YVeC4n2+qEvq5cGvhtUulkqlDeoviVBIkyOMB+EUPGBkl4/Mq1Qq67VuFITeEd+Vy+WK5gbg4efwbEwWdD4005pRwUbXUf8q9Alkk1s1myLx5ykZ8pq7N/A8od127wWnk+pnodfrrSb/2utuqZ8C862FVDeg7/EGl0W/Hx2fyKpHm1YtwHvaZrLhgYGmBvR75ulD6jUP4/NOp7NZM/F5DPO4rRdUD5j0Cf7m+7QyK5OlGXhOm8O8PghctQDjUfXsbRtlUfgv2fVsrNHvkNc+ifVCobDJ59sb631g3IRf4Ef4Hn6CP80rFosbOZ53zTzLfLXPhfah+TH0a6qjPW42mzXG074443e40Gq1tmn+n4E//U4mPWXH3W53Df/FtZpZVtg3hl3fYVHHGU/Cu5pZdkS3JKVmxhjjf+MXB5HD3Oq0ybsAAAAASUVORK5CYII=>

[image12]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAwAAAAYCAYAAADOMhxqAAAA8UlEQVR4XtWQPQrCQBCF11osJSj5KxZyCUEUW+8i2FgIgp3Yag7hKfxrPIMW3kGxisY3cbNsJhv7PBiyfPNm92WEqI+iKGoFQbBBTXmvpDAMdzCu6ey6rsQ55R4tmGe+7/dMRgOoucm00LhaGA0sOadG7DhO08IpUoPzvCGklG1E6yuW5OeSkP0AQ0dFoEpQT+7LhMbQdhMN4qIF57SdI2ckDOxRL851fi7wD+rBOeW/cEZS/7ItQJhHgLcCFFnMifVlwDPqLYxd0wLI7Hle17D+lN+C711FSDGw4j6tqg1ZpfIPOK8UzCfO/govjDmrmb6bzzzuTnhi3wAAAABJRU5ErkJggg==>