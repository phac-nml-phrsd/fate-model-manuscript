# fate-model-manuscript

**Mechanistic modelling of in-sewer viral fate and transport of SARS-CoV-2 to enhance wastewater disease surveillance strategies**

`fate-model-manuscript` is an R-based program developed to simulate the fate and transport of SARS-CoV-2 within urban sewer networks. It supports the estimation of viral losses and transport dynamics using data from the three major Winnipeg wastewater systems — North, South, and West plants.

This repository contains the mechanistic model, simulation workflows, and visualization scripts used in the associated preprint:\
🔗 [Research Square Preprint](https://www.researchsquare.com/article/rs-7774650/v1)

## Overview

This repository contains three main R scripts that should be executed in order.

------------------------------------------------------------------------

### 1. `calc-loss-stochastic.R`

This script estimates **in-sewer viral loss** by solid categories, for each individual conduit (pipe) and flow path.\
The estimation is **stochastic**, and model input parameters are defined as distributions in `parameters.R`.

**Execution requirements:** - Define: - `num.sim` → number of stochastic simulations - `seed` → random seed - Provide **hydraulic flow state data** for each pipe, obtained from the *Winnipeg Hydraulic InfoWorks ICM* model.\
The data file should be placed in: data/dwf_path_states/weekday.csv

**Outputs:**\
Simulation results are saved in the `out/` folder: - `out/stoch.params.rds` – stochastic parameters for each simulation run - `out/sim_df_loss_solid.rds` – simulation results by solid category - `out/sim_df_loss_total.rds` – aggregated simulation results (solid + liquid fractions)

**Estimate run time:** Using demo file for flow states, the approximate run time on normal laptop is 70 seconds for `num.sim = 2`.

### 2. `simu-analysis.R`

This script analyzes the simulation outputs by estimating: - Mean viral loss values - Loss per FSA (neighborhood) - Population-weighted viral loss - Minimum infection rate for detection

**Required input files:**\
- Outputs from `calc-loss-stochastic.R`: - `out/sim_df_loss_solid.rds` - `out/sim_df_loss_total.rds` - `out/stoch.params.rds` - `out/df.flow.rds` (flow states - cleaned version)

**Spatial and demographic data:** - Wastewater treatment plant polygons: `data/iw.polygon.csv` and `data/iw.wwtp.csv` - Demographic data (Statistics Canada 2021 Census): `data/census_English_CSV_data.csv`\
[Census Data Link](https://www12.statcan.gc.ca/census-recensement/2021/dp-pd/prof/details/page.cfm?Lang=E&SearchText=Winnipeg&DGUIDlist=2021A00054611040&GENDERlist=1,2,3&STATISTIClist=1,4&HEADERlist=0) - FSA (neighborhood) shapefiles: `data/shapefiles/postal-code.shp`

**Estimate run time:** Using real simulation results and demographics file, the approximate runtime on normal laptop is 112 seconds.

------------------------------------------------------------------------

### 3. `make-figures.R`

This script generates all the **figures** used in the *fate-model-manuscript*.\
To run this script, first execute `simu-analysis.R`.

**Output:** - All figures are saved in the `figs/` folder.

**Estimate run time:** Using real simulation results and demographics file, the approximate runtime on normal laptop is 95 seconds.

## System Requirements

-   **R**: \>= 4.2 (test on 4.2, 4.3, 4.4)

-   **OS**: macOS, Linux, Windows

-   **RAM/CPU**: depends on dataset size; typical laptop is fine.

## Install from Github

`git clone https://github.com/phac-nml-phrsd/fate-model-manuscript`

## Recreate the exact library

`install.packages("renv")`\
`renv::restore()`
