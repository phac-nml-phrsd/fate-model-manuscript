# Figure wrappers for the public Ottawa WWTP--neighbourhood-D benchmark.
# The data transformations and plot components live in utils_ottawa.R.

figure_Ottawa_raw_data <- function(df_ww,
                                   filename = file.path("figs", "figure_Ottawa_raw.png")) {
  plot <- plot_combined_biomarker_and_hydraulics_pub(df_ww)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = 10, height = 7,
                  dpi = 300, units = "in", bg = "white")
  invisible(plot)
}

figure_Ottawa_loss_estimate <- function(df_ww,
                                        filename = file.path("figs", "figure_Ottawa_loss.png")) {
  plot <- combine_boxplot_and_loss_heatmap(df_ww)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = 8, height = 8,
                  dpi = 300, units = "in", bg = "white")
  invisible(plot)
}
