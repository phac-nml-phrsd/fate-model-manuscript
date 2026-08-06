# Public stochastic-model demonstration input

`demo_dwf_path_states (weekday).csv` is a synthetic 60-conduit network with
the same schema and list-field encoding expected from the cleaned InfoWorks
path-state export. It contains no Winnipeg network measurements or asset IDs.

`demo_tss_lognormal.csv` contains fitted lognormal parameters for winter TSS
in log(mg/L), allowing the demo to run without the operational TSS workbook.

Regenerate and run the demonstration with:

```r
Rscript scripts/make-demo-dwf_path_states_weekday.R
Rscript calc-loss-stochastic.R
```

Defaults are `NUM_SIM=1`, `SEED=123`, and `out/demo-simulation/`. These are
demonstration results, not the full 500-simulation manuscript output.
