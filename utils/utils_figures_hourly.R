# Figure helpers for the public hourly aggregation sensitivity workflow.

figure_hourly_loss_comparison <- function(summary_file,
                                          filename = file.path("figs", "figure_hourly_loss_comparison.png")) {
  if (!file.exists(summary_file)) {
    stop("Hourly sensitivity summary not found: ", summary_file)
  }
  dat <- readr::read_csv(summary_file, show_col_types = FALSE)
  required <- c("wwtp", "fate_process", "model_version", "median_loss_percent")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop("Hourly summary is missing required columns: ", paste(missing, collapse = ", "))
  }
  plot <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = model_version, y = median_loss_percent / 100, fill = model_version)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.7) +
    ggplot2::facet_grid(wwtp ~ fate_process, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Proportion of median in-sewer loss per path") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = 6.5, height = 4.5,
                  dpi = 300, units = "in", bg = "white")
  invisible(plot)
}
