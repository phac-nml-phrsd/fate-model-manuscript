# Reproducible Ottawa neighbourhood-WWTP benchmark used in the revised paper.
# Run from the repository root:
#   Rscript scripts/upstream-loss-normalized.R

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

source("utils/utils_ottawa.R")

input_file <- file.path("data", "Ottawa", "Ottawa_data.xlsx")
out_dir <- file.path("out", "ottawa-benchmark")
fig_dir <- file.path(out_dir, "figures")

if (!file.exists(input_file)) stop("Ottawa input workbook not found: ", input_file)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

df_ww_raw <- load_data(input_file)
df_ww <- add_total_and_total_normalized(df_ww_raw)

box_start <- as.Date("2021-03-22")
box_end <- as.Date("2021-05-18")
ratio_result <- plot_log10_ratio_WWTP_over_D_boxplots(
  df_ww,
  box_start = box_start,
  box_end = box_end
)
loss_summary <- estimate_IQR_loss_from_log10_ratio(ratio_result$data)
figure_ottawa_loss <- combine_boxplot_and_loss_heatmap(df_ww)
figure_ottawa_raw <- plot_combined_biomarker_and_hydraulics_pub(df_ww)

readr::write_csv(df_ww, file.path(out_dir, "ottawa_processed_long.csv"))
readr::write_csv(ratio_result$data, file.path(out_dir, "ottawa_paired_ratios.csv"))
readr::write_csv(loss_summary, file.path(out_dir, "ottawa_apparent_loss_summary.csv"))

ggplot2::ggsave(
  file.path(fig_dir, "figure_Ottawa_loss.png"), figure_ottawa_loss,
  width = 8, height = 8, dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(fig_dir, "figure_Ottawa_loss.pdf"), figure_ottawa_loss,
  width = 8, height = 8, bg = "white"
)
ggplot2::ggsave(
  file.path(fig_dir, "figure_Ottawa_raw.png"), figure_ottawa_raw,
  width = 10, height = 7, dpi = 300, bg = "white"
)

message("Ottawa benchmark complete.")
message("Paired analysis window: ", box_start, " to ", box_end)
message("Summary: ", file.path(out_dir, "ottawa_apparent_loss_summary.csv"))
message("Revised figure: ", file.path(fig_dir, "figure_Ottawa_loss.png"))
