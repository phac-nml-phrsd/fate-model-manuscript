# Rebuild all figures supported by the inputs available in this public checkout.
# Run from the repository root:
#   Rscript make-figures.R
# Add --strict to fail when any manuscript input is unavailable.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(tidyr)
})

source(file.path("utils", "utils_plot.R"))
source(file.path("utils", "utils_figures.R"))
source(file.path("utils", "utils_figures_core.R"))
source(file.path("utils", "utils_figures_sensitivity.R"))
source(file.path("utils", "utils_ottawa.R"))
source(file.path("utils", "utils_figures_ottawa.R"))
source(file.path("utils", "utils_figures_hourly.R"))
source(file.path("utils", "utils_param.R"))

strict <- "--strict" %in% commandArgs(trailingOnly = TRUE)
fig_dir <- "figs"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
skipped <- character()

skip_figure <- function(label, files) {
  missing <- files[!file.exists(files)]
  skipped <<- c(skipped, paste0(label, ": ", paste(missing, collapse = ", ")))
  message("Skipping ", label, "; missing input(s): ", paste(missing, collapse = ", "))
}

run_figure <- function(label, files, code) {
  if (!all(file.exists(files))) {
    skip_figure(label, files)
    return(invisible(NULL))
  }
  message("Creating ", label, "...")
  force(code)
  invisible(NULL)
}

start <- Sys.time()
message("Start creating figures from explicit public inputs.")

# Ottawa benchmark: fully public and independently reproducible.
ottawa_file <- file.path("data", "Ottawa", "Ottawa_data.xlsx")
run_figure("Ottawa raw-data and attenuation figures", ottawa_file, {
  ottawa_raw <- load_data(ottawa_file)
  ottawa_data <- add_total_and_total_normalized(ottawa_raw)
  figure_Ottawa_raw_data(
    ottawa_data,
    filename = file.path(fig_dir, "figure_Ottawa_raw.png")
  )
  figure_Ottawa_loss_estimate(
    ottawa_data,
    filename = file.path(fig_dir, "figure_Ottawa_loss.png")
  )
})

# Compact public manuscript outputs. No raw InfoWorks states or full path lists
# are needed to rebuild the figures below.
figure_data_dir <- file.path("out", "figure-data")
summary_dir <- file.path("out", "manuscript-summary")
loss_summary_file <- file.path(summary_dir, "flow_path_loss_distributions.csv")
rna_summary_file <- file.path(figure_data_dir, "rna_scenario_loss_histogram.csv")
conduit_file <- file.path(figure_data_dir, "conduit_loss_summary.rds")
map_context_file <- file.path(figure_data_dir, "winnipeg_map_context.rds")
fate_artwork_file <- file.path(figure_data_dir, "figure_fate_processes.pdf")
flow_file <- file.path(figure_data_dir, "flow_hydraulic_figure_data.rds")

run_figure("loss distribution", loss_summary_file, {
  figure_loss_histo_summary(loss_summary_file, file.path(fig_dir, "figure_loss_histo.png"))
})

run_figure("combined fate-process/loss distribution", c(loss_summary_file, fate_artwork_file), {
  figure_loss_histo_fate_processes_summary(
    fate_artwork_file, loss_summary_file,
    file.path(fig_dir, "figure_loss_histo_fate.png")
  )
})

run_figure("conduit loss map", c(conduit_file, map_context_file), {
  map_context <- readRDS(map_context_file)
  figure_combo_lossmap(
    readRDS(conduit_file), map_context$catchments, map_context$wwtps,
    filename = file.path(fig_dir, "figure_combo_lossmap.png")
  )
})

run_figure("travel time, shear stress, and biofilm distributions", flow_file, {
  flow_data <- readRDS(flow_file)
  figure_travel_time(flow_data, file.path(fig_dir, "figure_travel_time.png"))
  figure_shear_stress_pipe(flow_data, file.path(fig_dir, "figure_shear_stress.png"))
  figure_biofilm_pipe(flow_data, file.path(fig_dir, "figure_biofilm.png"))
})

# Revised FSA maps and stochastic boxplots.
pop_loss_file <- file.path(figure_data_dir, "fsa_population_loss.rds")
infection_file <- file.path(figure_data_dir, "fsa_detection_rate.rds")
boxplot_file <- file.path(figure_data_dir, "stochastic_boxplot_metrics.rds")

run_figure("FSA loss map and stochastic effective-coverage boxplots", c(pop_loss_file, boxplot_file), {
  figure_combo_map_loss_fsa_boxplot(
    df_map = readRDS(pop_loss_file),
    df_boxplot = readRDS(boxplot_file),
    filename = file.path(fig_dir, "figure_combo_map_loss_pop_boxplot.png")
  )
})

run_figure("fraction-specific FSA loss map", pop_loss_file, {
  figure_map_loss_fsa_fraction(
    readRDS(pop_loss_file), size.text = 2.5, col.fsa = "white",
    filename = file.path(fig_dir, "figure_map_loss_fsa_fraction.png")
  )
})

run_figure("FSA detection map and stochastic minimum-detection boxplots", c(infection_file, boxplot_file), {
  figure_combo_detection_boxplot(
    df_map = readRDS(infection_file),
    df_boxplot = readRDS(boxplot_file),
    filename = file.path(fig_dir, "figure_combo_detection_boxplot.png")
  )
})

# Particle-class sensitivity figures.
particle_file <- file.path(summary_dir, "particle_association_sensitivity.csv")
particle_labels <- file.path(figure_data_dir, "particle_class_labels.rds")
run_figure("particle-distribution sensitivity figures", c(particle_file, particle_labels), {
  plot_particle_loss_distribution_summary(
    particle_file, particle_labels, "Pipe",
    file.path(fig_dir, "appendix_settling_pipe.pdf")
  )
  plot_particle_loss_distribution_summary(
    particle_file, particle_labels, "Flow path",
    file.path(fig_dir, "appendix_settling_path.pdf")
  )
})

run_figure("RNA-distribution sensitivity figures", c(rna_summary_file, conduit_file, pop_loss_file, infection_file), {
  conduit_data <- readRDS(conduit_file)
  population_data <- readRDS(pop_loss_file)
  infection_data <- readRDS(infection_file)
  figure_sensitivity_loss_histo_summary(
    rna_summary_file, conduit_data,
    file.path(fig_dir, "figure_sensitivity_loss_histo.png")
  )
  figure_sensitivity_map_loss(
    population_data, size.text = 1.5, col.fsa = "white",
    filename = file.path(fig_dir, "figure_sensitivity_loss_map.png")
  )
  figure_sensitivity_bar_plots(
    infection_data,
    filename = file.path(fig_dir, "figure_sensitivity_barplots.png")
  )
})

run_figure("solid-filtration sensitivity figure", character(), {
  figure_solid_filtration(file.path(fig_dir, "figure_solid_filtration.png"))
})

# Hourly aggregation figure from the summary table created by
# scripts/summarise-hourly-sensitivity.R.
hourly_summary <- file.path(figure_data_dir, "hourly_path_level_loss_summary.csv")
run_figure("hourly hydraulic-aggregation comparison", hourly_summary, {
  figure_hourly_loss_comparison(
    hourly_summary,
    filename = file.path(fig_dir, "slow_loss_daily_vs_hourly_true_boxplot.png")
  )
})

if (length(skipped)) {
  message("\nFigures skipped because their explicit inputs are unavailable:")
  message(paste0("- ", skipped, collapse = "\n"))
  if (strict) stop("One or more requested manuscript figures could not be generated.")
}

message("Figure generation complete in ", round(difftime(Sys.time(), start, units = "secs"), 2), " seconds.")
