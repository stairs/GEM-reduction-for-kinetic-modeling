# Kinetic-Parameter Sampling

This directory contains MATLAB scripts to perform Monte-Carlo sampling of enzyme‐kinetic parameters and compute flux control coefficients (FCCs) for:

- **Komkova et al.** GEM-embedded kinetic model
- **Shestov-derived** aerobic glycolysis reference model

---

## Prerequisites

- **MATLAB R2020b** or later  
- **SimBiology Toolbox**  
- **Symbolic Math Toolbox**  

---

## Quick Start

1. Optionally adjust at the top of each `sample_*.m`:
   - `iterations` (e.g. 10 000)  
   - Sampling bounds `F1` and `F2`  
   - Random seed via `seed`  
2. Run one of the sampling scripts:
   ```matlab
   >> sample_komkova
   >> sample_shestov
   ```
3. Locate results in `../out/` (e.g. `mc_komkova.mat`).

---

## Script Descriptions

### `ReadSBML.m`
- Imports SBML with `sbmlimport`  
- Extracts stoichiometric matrix, conserved moieties, steady‐state concentrations & fluxes from CSV  

### `sample_komkova.m` / `sample_shestov.m`
1. Load model & steady‐state (via `ReadSBML`)  
2. Generate symbolic Jacobian terms (via `compSymDeriv`)  
3. Define sampling ranges for kinetic parameters  
4. Monte Carlo loop:
   - Sample parameters
   - Build numeric Jacobian, test stability (max Re(λ) < 0)  
   - Compute FCCs
5. Save FCC arrays (`CJ_rec`), eigenvalues (`MaxRealEigens`), etc.

---

## Configuration Tips

- **Iterations:** Increase for better convergence  
- **Bounds (`F1`, `F2`):** Default covers 0.1×–10× nominal values  
- **Reproducibility:** Set `seed` before sampling  
- **Stability filter:** Only samples with negative eigenvalues retained  

---

## Acknowledgments

Original sampling implementation adapted from [klamt-lab/Models_E.coli_High_ATP_Demand](https://github.com/klamt-lab/Models_E.coli_High_ATP_Demand)  
Komkova _et al._ (2025), Shestov _et al._ (2014)
