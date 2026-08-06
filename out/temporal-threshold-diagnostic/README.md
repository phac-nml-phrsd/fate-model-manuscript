# Temporal threshold diagnostic summaries

These compact CSV files summarize the full midpoint-threshold analysis used
to create Supplementary Figure
`appendix_hydraulics_and_threshold_classification_day3.pdf`.

- `threshold_grid_used.csv`: midpoint settling and resuspension thresholds for
  the eight particle classes.
- `wwtp_particle_threshold_summary.csv`: pipe-day classification and exposure
  summaries by wastewater treatment plant and particle class.
- `wwtp_particle_daily_vs_hourly_exposure.csv`: daily-versus-hourly,
  HRT-weighted exposure by plant, day, and particle class.
- `wwtp_particle_text_summary.csv`: compact averages used to support the
  supplementary interpretation.

Thresholds are in N/m2. Exposure and difference columns are proportions on
the 0–1 scale; columns beginning with `pct_` are percentages. The published
tables use `threshold_mode = midpoint`, which has one simulation because each
class uses its fixed range midpoint.

The full pipe-level table is not distributed because it is approximately
243 MB and is derived from restricted Winnipeg hydraulic-model outputs. It
can be recreated with authorized inputs by setting
`save_pipe_level_table <- TRUE` in the diagnostic script.
