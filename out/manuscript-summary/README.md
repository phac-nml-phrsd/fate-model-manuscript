# Summarized manuscript outputs

This directory contains compact, non-restricted CSV tables supporting the
revised fate-model manuscript. These are full-manuscript summaries, not the
small public demonstration run, except that the threshold-classification table
uses one deterministic midpoint-threshold realization by design.

Raw Winnipeg InfoWorks hourly states, full flow paths, and multi-gigabyte
simulation objects are not included. The daily-versus-hourly and particle
tables retain only aggregated results required to reproduce manuscript plots.

`output_metadata.csv` records the simulation count, random seed, producing
script, manuscript mapping, and run scope for every published output.
`data_dictionary.csv` documents every column and its units.

Regenerate these tables from the compact figure-data bundle with:

```r
Rscript scripts/publish-summarized-manuscript-outputs.R
```
