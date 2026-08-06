---
editor_options:
  markdown:
    wrap: 72
---

# fate-model-manuscript

**Mechanistic modelling of in-sewer viral fate and transport of
SARS-CoV-2 to inform wastewater disease surveillance**

`fate-model-manuscript` is an R-based program developed to simulate the
fate and transport of SARS-CoV-2 within urban sewer networks. It
supports the estimation of viral losses and transport dynamics using
data from the three major Winnipeg wastewater systems --- North, South,
and West plants.

This repository contains the mechanistic model, simulation workflows,
and visualization scripts used in the associated preprint:\
🔗 [Research Square
Preprint](https://www.researchsquare.com/article/rs-7774650/v1)

## Reviewer quick start

From a fresh clone, restore the R environment, run the synthetic model
demo, rebuild the published summary tables, and reproduce every
manuscript figure:

``` sh
git clone https://github.com/phac-nml-phrsd/fate-model-manuscript.git
cd fate-model-manuscript
Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')"
Rscript -e "renv::restore(prompt = FALSE)"
Rscript calc-loss-stochastic.R
Rscript scripts/publish-summarized-manuscript-outputs.R
Rscript make-figures.R --strict
```

The model demo and manuscript reproduction are deliberately separate.
The demo executes the mechanistic simulation on a synthetic network.
Manuscript figures use the included compact summaries of the full
500-simulation run.

------------------------------------------------------------------------

### 1. `calc-loss-stochastic.R`

This script estimates **in-sewer viral loss** by solid categories, for
each individual conduit (pipe) and flow path.\
The estimation is **stochastic**, and model input parameters are defined
as distributions in `parameters.R`.

The repository includes a synthetic 60-conduit reviewer demonstration:

``` r
Rscript calc-loss-stochastic.R
```

Defaults are `NUM_SIM=1`, `SEED=123`, input
`data/demo_dwf_path_states (weekday).csv`, and output
`out/demo-simulation/`. Optional `NUM_SIM`, `SEED`, `FLOW_INPUT`, and
`OUTPUT_DIR` environment variables override these values. The run writes
cleaned flow, geometry, path-index, parameter, solid-loss, and
total-loss RDS files. It completed in approximately 1.4 seconds in the
test environment.

See `data/README_demo.md` for the demo schema and fitted-TSS fallback.
The manuscript results use 500 simulations and are published separately
as compact summaries under `out/manuscript-summary/`.

### 2. `simu-analysis.R`

This is the authors' full-data analysis pipeline for conduit means, FSA
loss, population-equivalent loss, effective surveillance coverage, and
minimum detectable prevalence. It requires restricted Winnipeg hydraulic
products and the full Statistics Canada spatial workflow. Its compact
full-manuscript results are already published under `out/figure-data/`
and `out/manuscript-summary/`.

The public one-simulation demo outputs are isolated under
`out/demo-simulation/` and are not substituted for manuscript results.

------------------------------------------------------------------------

### 3. `make-figures.R`

This script generates every revised manuscript figure from the public
Ottawa workbook and compact files under `out/`. It is safe in a clean R
session and does not use `.RData` or objects from the global
environment.

``` r
Rscript make-figures.R
Rscript make-figures.R --strict  # fail if any expected input is unavailable
```

The legacy `utils/utils_figures.R` and `utils/utils_plot.R` remain the
shared base modules. `--strict` fails immediately if a required public
input is missing; a complete checkout passes this check and creates all
figures.

All figures are saved in `figs/`. The verified clean-session strict run
created 18 figures in approximately 34 seconds in the test environment.

## System Requirements

-   **R**: \>= 4.2 (test on 4.2, 4.3, 4.4)

-   **OS**: macOS, Linux, Windows

-   **RAM/CPU**: depends on dataset size; typical laptop is fine for
    demo simulation and creating manuscript figures.

## Public data layout

-   `data/demo_dwf_path_states (weekday).csv`: synthetic reviewer
    network.
-   `data/demo_tss_lognormal.csv`: public fitted TSS parameters for the
    demo.
-   `data/Ottawa/Ottawa_data.xlsx`: public WWTP--neighbourhood-D
    benchmark (field-work dataset provided by Dr. Banu Ormecy).
-   `out/demo-simulation/`: verified one-simulation demonstration
    outputs.
-   `out/figure-data/`: compact plot-ready CSV/RDS/PDF inputs.
-   `out/manuscript-summary/`: publication-facing CSV summaries,
    dictionary, units, seed, simulation count, provenance, and
    manuscript mapping.

## Hourly-resolved hydraulic sensitivity

The revised manuscript compares the daily-averaged model with an
hourly-resolved implementation. Hydraulic thresholds and path survival
are calculated for each hour before the remaining fractions are averaged
into a 24-hour composite.

The authors' full-data workflow runs in this order:

``` r
Rscript scripts/prepare-hourly-path-states.R
Rscript scripts/calc-loss-hourly-sensitivity.R
Rscript scripts/summarise-hourly-sensitivity.R
```

These three commands require the complete restricted hourly series and
are included for method audit, not as reviewer quick-start steps. The
public configuration uses one stochastic simulation (`num.sim <- 1`) and
the hard-coded example date `20010101`. The manuscript analysis used 500
simulations and 24 hourly hydraulic files per day.

One excerpt from a real 00:00 InfoWorks hourly export is provided under
`data/hourly-data/20010101/dwf/` to document the input schema. The full
hourly exports and matching Winnipeg path/topology template are
restricted, so this single-hour excerpt does not reproduce the reported
24-hour comparison. See `data/hourly-data/README.md` for the expected
columns and directory layout.

### Temporal threshold diagnostic

The supplementary daily-versus-hourly settling and resuspension
diagnostic uses the complete prepared hourly path-state files:

``` r
Rscript scripts/calc-temporal-threshold-diagnostic.R
```

The default `THRESHOLD_MODE=midpoint` uses the midpoint of each particle
class's settling and resuspension ranges. `THRESHOLD_MODE=stochastic`
samples the same uniform ranges used by the main Monte Carlo model.
`NUM_SIM` controls the number of stochastic threshold draws and `SEED`
makes those draws reproducible. Outputs are written to
`out/temporal-threshold-diagnostic/`, with figures under its `figures/`
subdirectory. Compact manuscript-run summary tables are included in the
output directory.

Reviewers can inspect the published diagnostic tables directly under
`out/temporal-threshold-diagnostic/` and
`out/manuscript-summary/threshold_classification_diagnostic.csv`;
rerunning the full diagnostic requires the restricted hourly series.

## Ottawa neighbourhood--WWTP benchmark

The revised manuscript's indirect Ottawa benchmark is reproduced with:

``` r
Rscript scripts/upstream-loss-normalized.R
```

The script reads `data/Ottawa/Ottawa_data.xlsx`, calculates raw,
PMMoV-normalized, total-solids-normalized, and per-capita-flow-adjusted
paired WWTP-to-neighbourhood ratios, and writes processed data, summary
statistics, and the revised figures to `out/ottawa-benchmark/`.

See `data/Ottawa/README.md` for provenance, units, sample-fraction
definitions, exclusions, and the redistribution-status note.

## Summarized manuscript outputs

Compact, non-restricted full-manuscript results are published under
`out/manuscript-summary/`. They support flow-path loss distributions,
FSA loss by sample fraction, population-equivalent loss, effective
surveillance coverage, minimum detectable prevalence,
daily-versus-hourly hydraulic aggregation, threshold classification, and
particle-association sensitivity.

Each table is CSV. `out/manuscript-summary/output_metadata.csv` records
the simulation count, seed, producing script, manuscript figure/table
mapping, and whether the table represents a demonstration or full
manuscript run. `out/manuscript-summary/data_dictionary.csv` documents
columns and units.

The underlying plot-ready spatial objects are retained as compact RDS
files in `out/figure-data/`. Raw InfoWorks hourly states and
multi-gigabyte path-level simulation objects are not distributed.

Rebuild all publication-facing CSV tables from the included compact
figure-data bundle:

``` r
Rscript scripts/publish-summarized-manuscript-outputs.R
```

The author-side scripts `scripts/prepare-public-figure-data.R` and
`scripts/prepare-public-particle-figure-data.R` document how the compact
bundle was derived from the full restricted run. They are provenance
utilities, not required reviewer steps.

To rebuild every manuscript figure from a clean R session:

``` r
Rscript make-figures.R --strict
```
