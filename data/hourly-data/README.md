# Hourly hydraulic input

This directory documents the raw hourly InfoWorks ICM input used by the
hourly-resolved sensitivity workflow.

## Public example

One excerpt from a real hourly export is included at:

`20010101/dwf/conduit_states_dwf_20010101_0000.csv`

It contains real pipe records from the hydraulic state at 00:00 on the
example date. The excerpt is provided to document the real column names,
units, and file layout and to exercise the conversion step when the matching
path/topology template is available. It is intentionally a small subset of
the full citywide export.

The example file contains these hydraulic fields:

- pipe identifier in the first unnamed column (`Unnamed: 0` when read by R);
- upstream and downstream node identifiers;
- conduit geometry and physical attributes;
- `depth`, `flow`, and `vol`;
- conduit and wetted cross-sectional quantities;
- `flow_velocity` and `residence_time`;
- `conduit_av`, `conduit_hrt`, and `conduit_hrt_inv`;
- `conduit_roughness`, `conduit_ff`, and `conduit_ss`.

## Full manuscript run

The revised manuscript analysis used 24 hourly files per analysed day. The
remaining hourly records, other InfoWorks exports, and the full Winnipeg
path/topology template are not distributed publicly. Therefore, the single
file in this repository is an input example, not sufficient to reproduce the
reported 24-hour composite loss estimates.

With authorized inputs, use this layout:

```text
data/hourly-data/
  YYYYMMDD/
    dwf/
      conduit_states_dwf_YYYYMMDD_0000.csv
      conduit_states_dwf_YYYYMMDD_0100.csv
      ...
      conduit_states_dwf_YYYYMMDD_2300.csv
```

The matching template must be placed at
`data/dwf_path_states (weekday).csv`. Then run from the repository root:

```sh
Rscript scripts/prepare-hourly-path-states.R
Rscript scripts/calc-loss-hourly-sensitivity.R
Rscript scripts/summarise-hourly-sensitivity.R
```
