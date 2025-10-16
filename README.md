# fate-model-manuscript

**Mechanistic modelling of in-sewer viral fate and transport of SARS-CoV-2 to enhance wastewater disease surveillance strategies**

[![R-CMD-check](https://github.com/phac-nml-phrsd/fate-model-manuscript/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phac-nml-phrsd/fate-model-manuscript/actions/workflows/R-CMD-check.yaml)

`fate-model-manuscript` is an R-based program developed to simulate the fate and transport of SARS-CoV-2 within urban sewer networks. It supports the estimation of viral losses and transport dynamics using data from the three major Winnipeg wastewater systems — North, South, and West plants.

This repository contains the mechanistic model, simulation workflows, and visualization scripts used in the associated preprint:\
🔗 [Research Square Preprint](https://www.researchsquare.com/article/rs-7774650/v1)

## Overview

A paragraph describing the problem, your approach, and what users get out of it (mirroring mgcpy’s clear “what this is” tone).\
If relevant, link to a website or docs.

## System Requirements

-   **R**: \>= 4.2 (test on 4.2, 4.3, 4.4)
-   **OS**: macOS, Linux, Windows
-   **RAM/CPU**: depends on dataset size; typical laptop is fine for examples.

## Installation

\`\`\`r \# Stable install.packages("pak") \# optional but fast pak::pak("github::USER/REPO") \# or CRAN if you publish there

# Dev (from source)

pak::pak("github::[USER/REPO\@main](mailto:USER/REPO@main){.email}")
