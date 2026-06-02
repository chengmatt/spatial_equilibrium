Code repository accompanying:

**Ensuring Consistency in Spatially-Stratified Stock Assessment Models: 
An Analytical Solution for Equilibrium Population Structure with 
Connectivity for Age- and Size-Structured Models**

Cheng MLH, Miller TJ, Goethel DR, Cunningham CJ. *Fisheries Research* (in review)

---

## Overview

This repository contains R code for deriving equilibrium population 
structure in spatially-stratified stock assessment models, including 
both analytical and iterative approaches for age- and size-structured 
populations. The analytical solution is based on a matrix geometric 
series and explicitly accounts for movement among regions, providing 
a consistent and internally coherent approach to equilibrium 
calculations relative to commonly used approximations.

---

## Repository Contents

| File | Description |
|------|-------------|
| `init_naa_demo.R` | Demonstration of analytical and iterative equilibrium solutions for an **age-structured** population across two regions (corresponds to deterministic calculations in the manuscript) |
| `init_nas_demo.R` | Demonstration of analytical and iterative equilibrium solutions for a **size-structured** population across two regions (corresponds to deterministic calculations in the manuscript) |
| `spt_equil_pop_code_demo.R` | Additional standalone demonstration of the analytical equilibrium solution |
| `spatial_plus_group_model_runs.R` | Code for fitting the Alaska sablefish multi-region stock assessment models under the four equilibrium formulations (no_move_all, no_move_plus, iter_bins, matrix) and running deterministic and stochastic population projections (corresponds to the Alaska sablefish case study in the manuscript) |
| `model_results.RDS` | Saved model outputs from the Alaska sablefish case study |

### Output Figures

| File | Description |
|------|-------------|
| `plus_group_demo_plot.png` | Deterministic equilibrium abundance across scenarios — age-structured (Figure 1 in manuscript) |
| `size_struct_demo_plot.png` | Deterministic equilibrium abundance across scenarios — size-structured (Figure S2.2 in manuscript) |
| `movement_scenarios.png` | Illustration of movement scenarios used in deterministic calculations (Figure S2.1 in manuscript) |
| `sablefish_plus_group_plot.png` | Alaska sablefish case study model estimates (Figure 2 in manuscript) |
| `sablefish_plus_group_proj_plot.png` | Alaska sablefish deterministic and stochastic projections (Figure 3 in manuscript) |
| `sablefish_rec_move_plot.png` | Sablefish recruitment and movement estimates (Figure S2.3 in manuscript) |

---

## Dependencies

All code is written in R. The sablefish case study requires the 
[SPoRC](https://github.com/chengmatt/SPoRC) stock assessment package.

---
