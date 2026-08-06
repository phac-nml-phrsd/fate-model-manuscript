# Build compact, analysis-ready public inputs for every manuscript figure.
# This is a one-time release-preparation script; it reads the authors' full
# source outputs but never copies raw InfoWorks hourly states or path lists.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(sf)
})

source(file.path("utils", "utils.R"))
source(file.path("utils", "utils_plot.R"))
source(file.path("utils", "utils_param.R"))
source(file.path("utils", "parameters.R"))
source(file.path("utils", "utils_param_stochastic.R"))
source(file.path("utils", "utils_figures_sensitivity.R"))

source_repo <- Sys.getenv(
  "SOURCE_REPO",
  unset = file.path(dirname(normalizePath(".")), "Stochastic-Fate-Model")
)
source_repo <- normalizePath(source_repo, mustWork = TRUE)
output_dir <- file.path("out", "figure-data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_path <- function(...) file.path(source_repo, ...)
require_source <- function(path) if (!file.exists(path)) stop("Required source file not found: ", path)

total_file <- source_path("out", "sim_df_loss_total.rds")
flow_file <- source_path("out", "df.flow.rds")
for (path in c(total_file, flow_file)) require_source(path)

message("Loading full stochastic total-loss output...")
loss <- readRDS(total_file)
loss_dt <- data.table::as.data.table(loss)

# Histogram-ready tables retain exact counts and summary lines but no individual
# paths. Bins match the 0.05 width used by the revised manuscript figures.
summarise_histogram <- function(values, group, group_name) {
  values <- pmin(pmax(values, 0), 1)
  bins <- pmin(floor(values / 0.05), 19L)
  counts <- data.table(wwtp = loss_dt$wwtp, bin = bins)[,
    .(count = .N), by = .(wwtp, bin)
  ][, `:=`(
    bin_left = bin * 0.05,
    bin_right = (bin + 1) * 0.05,
    bin_mid = (bin + 0.5) * 0.05
  )]
  stats <- data.table(wwtp = loss_dt$wwtp, value = values)[,
    .(mean_value = mean(value), median_value = median(value)), by = wwtp
  ]
  counts <- merge(counts, stats, by = "wwtp", all.x = TRUE)
  counts[, (group_name) := group]
  counts
}

loss_distribution <- data.table::rbindlist(list(
  summarise_histogram(1 - loss_dt$remain.set.slow, "Settling / Resuspension", "process"),
  summarise_histogram(1 - loss_dt$remain.bio.path, "Biofilm Adsorption", "process"),
  summarise_histogram(1 - loss_dt$remain.deg.liq.path, "Biodegradation", "process")
), use.names = TRUE)
readr::write_csv(loss_distribution, file.path(output_dir, "loss_distribution_histogram.csv"))

rna_scenario_distribution <- data.table::rbindlist(list(
  summarise_histogram(1 - loss_dt$remain.set.slow, "Skewed (Base)", "scenario"),
  summarise_histogram(1 - loss_dt$remain.set.fast, "Highly Skewed", "scenario"),
  summarise_histogram(1 - loss_dt$remain.set.hmg, "Homogeneous", "scenario")
), use.names = TRUE)
readr::write_csv(rna_scenario_distribution, file.path(output_dir, "rna_scenario_loss_histogram.csv"))

message("Calculating conduit-level means...")
mean_loss <- loss_dt[, .(
  mean.remain.bio.path = mean(remain.bio.path),
  mean.remain.deg.liq.path = mean(remain.deg.liq.path),
  mean.remain.bio.deg.path = mean(remain.bio.deg.path),
  mean.remain.set.hmg = mean(remain.set.hmg),
  mean.remain.set.fast = mean(remain.set.fast),
  mean.remain.set.slow = mean(remain.set.slow),
  mean.loss.total.hmg = mean(loss.total.hmg),
  mean.loss.total.fast = mean(loss.total.fast),
  mean.loss.total.slow = mean(loss.total.slow)
), by = node_id]
mean_loss <- merge(
  unique(loss_dt[, .(node_id, wwtp)]), mean_loss,
  by = "node_id", all.y = TRUE
)

message("Loading compact hydraulic fields and geometry...")
flow <- readRDS(flow_file)
flow_figure <- flow %>%
  dplyr::select(dplyr::any_of(c("wwtp", "path_hrt", "conduit_ss", "conduit_av"))) %>%
  dplyr::distinct()
saveRDS(flow_figure, file.path(output_dir, "flow_hydraulic_figure_data.rds"), compress = "xz")

geometry_lookup <- flow %>%
  dplyr::select(node_id, geometry) %>%
  dplyr::distinct(node_id, .keep_all = TRUE)
conduit_loss <- dplyr::left_join(as.data.frame(mean_loss), geometry_lookup, by = "node_id")
conduit_loss <- conduit_loss %>%
  dplyr::mutate(conduit_id = sprintf("C%05d", dplyr::row_number())) %>%
  dplyr::select(conduit_id, dplyr::everything(), -node_id)
saveRDS(conduit_loss, file.path(output_dir, "conduit_loss_summary.rds"), compress = "xz")
rm(flow_figure, mean_loss)
gc()

# Public map context used by the conduit figure.
iw_polygon <- read.csv(source_path("data", "iw.polygon.csv"))
iw_wwtp <- read.csv(source_path("data", "iw.wwtp.csv"))
saveRDS(list(catchments = iw_polygon, wwtps = iw_wwtp),
        file.path(output_dir, "winnipeg_map_context.rds"), compress = "xz")

message("Loading public census and FSA boundary data...")
fsa_ids <- get_fsa_ids()
demo_file <- source_path("data", "demographics", "census_English_CSV_data.csv")
demo_raw <- data.table::fread(
  demo_file,
  select = c("CENSUS_YEAR", "DGUID", "GEO_NAME", "CHARACTERISTIC_NAME", "C1_COUNT_TOTAL"),
  showProgress = FALSE
)
demographic <- get_demographics(as.data.frame(demo_raw), fsa_ids)
rm(demo_raw)
fsa_polygons <- sf::st_read(
  source_path("data", "shapefiles", "postal-code.shp"),
  quiet = TRUE
) %>%
  dplyr::filter(PRNAME == "Manitoba", CFSAUID %in% fsa_ids)

public_root <- getwd()
setwd(source_repo)
stochastic_parameters <- load_stoch_prms(num.sample = 500, seed = 123)
setwd(public_root)
message("Calculating FSA population-loss and detection-rate products...")
loss_fsa <- calcu_loss_fsa(conduit_loss, demographic, fsa_polygons, stochastic_parameters)
population_loss <- calculate_pop_loss(loss_fsa, stochastic_parameters)
detection_rate <- calcu_inf_rate(population_loss, stochastic_parameters)
saveRDS(population_loss, file.path(output_dir, "fsa_population_loss.rds"), compress = "xz")
saveRDS(detection_rate, file.path(output_dir, "fsa_detection_rate.rds"), compress = "xz")

sim_start <- as.integer(Sys.getenv("SIM_START", unset = "1"))
sim_end <- as.integer(Sys.getenv("SIM_END", unset = "500"))
if (is.na(sim_start) || is.na(sim_end) || sim_start < 1 || sim_end > 500 || sim_start > sim_end) {
  stop("SIM_START and SIM_END must define a range within 1:500.")
}
parts_dir <- file.path(output_dir, "stochastic-metric-parts")
dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
message("Calculating stochastic boxplot metrics for simulations ", sim_start, "--", sim_end, "...")
for (sim_id in seq.int(sim_start, sim_end)) {
  part_file <- file.path(parts_dir, sprintf("sim_%03d.rds", sim_id))
  if (file.exists(part_file)) next
  message("Calculating boxplot metrics for simulation ", sim_id, "...")
  part <- calculate_simulation_boxplot_metrics(
    sim_loss = loss[loss$sim == sim_id, , drop = FALSE],
    df_polygons = geometry_lookup,
    demographic = demographic,
    fsa_polygons = fsa_polygons,
    stochastic_parameters = stochastic_parameters[stochastic_parameters$sim == sim_id, , drop = FALSE]
  )
  saveRDS(part, part_file, compress = "xz")
}

part_files <- file.path(parts_dir, sprintf("sim_%03d.rds", 1:500))
if (all(file.exists(part_files))) {
  boxplot_metrics <- dplyr::bind_rows(lapply(part_files, readRDS))
} else {
  boxplot_metrics <- NULL
  message("Checkpointed requested range; run remaining SIM_START/SIM_END ranges before final assembly.")
}
boxplot_columns <- c(
  "sim", "wwtp",
  "eff.pop.bio.deg.path.wwtp", "eff.pop.set.slow.wwtp", "eff.pop.total.slow.wwtp",
  "eff.pop.bio.deg.path.city", "eff.pop.set.slow.city", "eff.pop.total.slow.city",
  "wwtp.inf.rate.bio.deg.highshed", "wwtp.inf.rate.set.slow.highshed", "wwtp.inf.rate.slow.highshed",
  "city.inf.rate.bio.deg.highshed", "city.inf.rate.set.slow.highshed", "city.inf.rate.slow.highshed"
)
if (!is.null(boxplot_metrics)) {
  boxplot_metrics <- sf::st_drop_geometry(boxplot_metrics) %>%
    dplyr::select(dplyr::all_of(boxplot_columns)) %>%
    dplyr::distinct()
  saveRDS(boxplot_metrics, file.path(output_dir, "stochastic_boxplot_metrics.rds"), compress = "xz")
}

# Small parameter labels required by particle-class plots.
particle_labels <- load_prms()[, c("part.class", "vel.set")]
saveRDS(particle_labels, file.path(output_dir, "particle_class_labels.rds"), compress = "xz")

# Public hourly figure table; no raw hourly states are copied.
hourly_source <- source_path(
  "out", "hourly-sensitivity", "figures",
  "path_level_stats_slow_loss_daily_vs_hourly.csv"
)
require_source(hourly_source)
file.copy(
  hourly_source,
  file.path(output_dir, "hourly_path_level_loss_summary.csv"),
  overwrite = TRUE
)

# Artwork used as panel A of the fate-process/loss composite.
fate_artwork <- source_path("doc", "figs", "figure_fate_processes.pdf")
require_source(fate_artwork)
file.copy(fate_artwork, file.path(output_dir, "figure_fate_processes.pdf"), overwrite = TRUE)

message("Public figure-data bundle complete: ", output_dir)
