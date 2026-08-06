# Publish compact CSV outputs supporting the revised manuscript figures.
# Run after scripts/prepare-public-figure-data.R has created out/figure-data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(tidyr)
})

input_dir <- file.path("out", "figure-data")
output_dir <- file.path("out", "manuscript-summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "loss_distribution_histogram.csv",
  "fsa_population_loss.rds",
  "fsa_detection_rate.rds",
  "hourly_path_level_loss_summary.csv",
  "particle_loss_distribution.csv",
  "rna_scenario_loss_histogram.csv",
  "stochastic_boxplot_metrics.rds"
)
missing <- required[!file.exists(file.path(input_dir, required))]
if (length(missing)) stop("Missing compact figure data: ", paste(missing, collapse = ", "))

# 1. Flow-path loss distributions (pre-binned at the manuscript width of 0.05).
flow_path <- readr::read_csv(
  file.path(input_dir, "loss_distribution_histogram.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    wwtp,
    fate_process = process,
    bin_left_proportion = bin_left,
    bin_right_proportion = bin_right,
    bin_mid_proportion = bin_mid,
    flow_path_count = count,
    mean_loss_proportion = mean_value,
    median_loss_proportion = median_value
  )
write_csv(flow_path, file.path(output_dir, "flow_path_loss_distributions.csv"))

# 2--4. FSA loss, population-equivalent loss and effective coverage.
population <- readRDS(file.path(input_dir, "fsa_population_loss.rds")) %>%
  sf::st_drop_geometry()

fsa_loss <- population %>%
  transmute(
    fsa = CFSAUID,
    wwtp,
    population = C1_COUNT_TOTAL,
    total_loss_proportion = fsa.ave.loss.total.slow,
    solids_loss_proportion = fsa.ave.loss.set.slow,
    liquid_loss_proportion = fsa.ave.loss.bio.deg
  ) %>%
  pivot_longer(
    ends_with("_loss_proportion"),
    names_to = "sample_fraction",
    values_to = "mean_loss_proportion"
  ) %>%
  mutate(sample_fraction = recode(
    sample_fraction,
    total_loss_proportion = "Total",
    solids_loss_proportion = "Solids",
    liquid_loss_proportion = "Liquid"
  ))
write_csv(fsa_loss, file.path(output_dir, "fsa_mean_loss_by_sample_fraction.csv"))

population_equivalent <- population %>%
  transmute(
    fsa = CFSAUID,
    wwtp,
    population = C1_COUNT_TOTAL,
    total_population_equivalent_loss = pop.signal.loss.total.slow.fsa,
    solids_population_equivalent_loss = pop.signal.loss.set.slow.fsa,
    liquid_population_equivalent_loss = pop.signal.loss.bio.deg.path.fsa
  ) %>%
  pivot_longer(
    ends_with("_population_equivalent_loss"),
    names_to = "sample_fraction",
    values_to = "population_equivalent_loss_persons"
  ) %>%
  mutate(sample_fraction = recode(
    sample_fraction,
    total_population_equivalent_loss = "Total",
    solids_population_equivalent_loss = "Solids",
    liquid_population_equivalent_loss = "Liquid"
  ))
write_csv(population_equivalent, file.path(output_dir, "population_equivalent_loss.csv"))

coverage_wwtp <- population %>%
  transmute(
    location = stringr::str_to_title(wwtp),
    Total = eff.pop.total.slow.wwtp,
    Solids = eff.pop.set.slow.wwtp,
    Liquid = eff.pop.bio.deg.path.wwtp
  ) %>%
  distinct()
coverage_city <- population %>%
  transmute(
    location = "City",
    Total = eff.pop.total.slow.city,
    Solids = eff.pop.set.slow.city,
    Liquid = eff.pop.bio.deg.path.city
  ) %>%
  distinct()
effective_coverage <- bind_rows(coverage_wwtp, coverage_city) %>%
  pivot_longer(c(Total, Solids, Liquid), names_to = "sample_fraction",
               values_to = "effective_surveillance_coverage_percent")
write_csv(effective_coverage, file.path(output_dir, "effective_surveillance_coverage.csv"))

stochastic <- readRDS(file.path(input_dir, "stochastic_boxplot_metrics.rds"))
stochastic_coverage <- bind_rows(
  stochastic %>% transmute(
    simulation = sim, location = stringr::str_to_title(wwtp),
    Total = eff.pop.total.slow.wwtp,
    Solids = eff.pop.set.slow.wwtp,
    Liquid = eff.pop.bio.deg.path.wwtp
  ) %>% distinct(),
  stochastic %>% transmute(
    simulation = sim, location = "City",
    Total = eff.pop.total.slow.city,
    Solids = eff.pop.set.slow.city,
    Liquid = eff.pop.bio.deg.path.city
  ) %>% distinct()
) %>%
  pivot_longer(c(Total, Solids, Liquid), names_to = "sample_fraction",
               values_to = "effective_surveillance_coverage_percent")
write_csv(stochastic_coverage,
          file.path(output_dir, "stochastic_effective_surveillance_coverage.csv"))

# 5. Minimum detectable prevalence for FSA, WWTP and city scales.
detection <- readRDS(file.path(input_dir, "fsa_detection_rate.rds")) %>%
  sf::st_drop_geometry()
minimum_prevalence <- detection %>%
  transmute(
    fsa = CFSAUID,
    wwtp,
    fsa_total_high_shedding_percent = inf.rate.slow.highshed,
    fsa_solids_high_shedding_percent = inf.rate.set.slow.highshed,
    fsa_liquid_high_shedding_percent = inf.rate.bio.deg.highshed,
    fsa_no_loss_high_shedding_percent = inf.rate.no.loss.highshed,
    fsa_total_low_shedding_percent = inf.rate.slow.lowshed,
    fsa_no_loss_low_shedding_percent = inf.rate.no.loss.lowshed,
    wwtp_total_high_shedding_percent = wwtp.inf.rate.slow.highshed,
    wwtp_solids_high_shedding_percent = wwtp.inf.rate.set.slow.highshed,
    wwtp_liquid_high_shedding_percent = wwtp.inf.rate.bio.deg.highshed,
    city_total_high_shedding_percent = city.inf.rate.slow.highshed,
    city_solids_high_shedding_percent = city.inf.rate.set.slow.highshed,
    city_liquid_high_shedding_percent = city.inf.rate.bio.deg.highshed
  )
write_csv(minimum_prevalence, file.path(output_dir, "minimum_detectable_prevalence.csv"))

stochastic_prevalence <- bind_rows(
  stochastic %>% transmute(
    simulation = sim, location = stringr::str_to_title(wwtp),
    Total = wwtp.inf.rate.slow.highshed,
    Solids = wwtp.inf.rate.set.slow.highshed,
    Liquid = wwtp.inf.rate.bio.deg.highshed
  ) %>% distinct(),
  stochastic %>% transmute(
    simulation = sim, location = "City",
    Total = city.inf.rate.slow.highshed,
    Solids = city.inf.rate.set.slow.highshed,
    Liquid = city.inf.rate.bio.deg.highshed
  ) %>% distinct()
) %>%
  pivot_longer(c(Total, Solids, Liquid), names_to = "sample_fraction",
               values_to = "minimum_detectable_prevalence_percent")
write_csv(stochastic_prevalence,
          file.path(output_dir, "stochastic_minimum_detectable_prevalence.csv"))

# 6. Daily-versus-hourly hydraulic aggregation summary.
hourly <- read_csv(
  file.path(input_dir, "hourly_path_level_loss_summary.csv"),
  show_col_types = FALSE
)
daily_hourly <- hourly %>%
  group_by(wwtp, fate_process, model_version) %>%
  summarise(
    n_paths = n_distinct(node_id),
    median_of_path_medians_percent = median(median_loss_percent, na.rm = TRUE),
    mean_of_path_medians_percent = mean(median_loss_percent, na.rm = TRUE),
    q25_of_path_medians_percent = quantile(median_loss_percent, 0.25, na.rm = TRUE),
    q75_of_path_medians_percent = quantile(median_loss_percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(daily_hourly, file.path(output_dir, "daily_versus_hourly_comparison.csv"))

# 7. Threshold-classification diagnostic summary.
threshold_file <- file.path(
  "out", "temporal-threshold-diagnostic", "wwtp_particle_threshold_summary.csv"
)
if (!file.exists(threshold_file)) stop("Threshold diagnostic summary not found: ", threshold_file)
file.copy(
  threshold_file,
  file.path(output_dir, "threshold_classification_diagnostic.csv"),
  overwrite = TRUE
)

# 8. Particle-association and RNA partition sensitivity scenarios.
particle <- read_csv(
  file.path(input_dir, "particle_loss_distribution.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    wwtp,
    particle_class = part.class,
    hydraulic_level = level,
    bin_left_proportion = bin_left,
    bin_right_proportion = bin_right,
    bin_mid_proportion = bin_mid,
    observation_count = count
  )
write_csv(particle, file.path(output_dir, "particle_association_sensitivity.csv"))

rna_partition <- read_csv(
  file.path(input_dir, "rna_scenario_loss_histogram.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    wwtp,
    particle_association_scenario = scenario,
    bin_left_proportion = bin_left,
    bin_right_proportion = bin_right,
    bin_mid_proportion = bin_mid,
    flow_path_count = count,
    mean_loss_proportion = mean_value,
    median_loss_proportion = median_value
  )
write_csv(rna_partition, file.path(output_dir, "rna_partition_sensitivity.csv"))

# File-level provenance and manuscript mapping.
metadata <- tibble::tribble(
  ~file, ~simulation_count, ~seed, ~source_script, ~manuscript_mapping, ~run_scope,
  "flow_path_loss_distributions.csv", 500L, 123L, "calc-loss-stochastic.R; scripts/prepare-public-figure-data.R", "Revised loss-distribution and combined fate-process/loss figures", "Full manuscript run",
  "fsa_mean_loss_by_sample_fraction.csv", 500L, 123L, "simu-analysis.R; scripts/prepare-public-figure-data.R", "Revised FSA loss maps, including fraction-specific map", "Full manuscript run",
  "population_equivalent_loss.csv", 500L, 123L, "simu-analysis.R; scripts/prepare-public-figure-data.R", "Population-equivalent loss map and summary", "Full manuscript run",
  "effective_surveillance_coverage.csv", 500L, 123L, "simu-analysis.R; scripts/prepare-public-figure-data.R", "Effective surveillance coverage panels", "Full manuscript run",
  "stochastic_effective_surveillance_coverage.csv", 500L, 123L, "scripts/prepare-public-figure-data.R", "Revised stochastic effective-surveillance-coverage boxplots", "Full manuscript run",
  "minimum_detectable_prevalence.csv", 500L, 123L, "simu-analysis.R; scripts/prepare-public-figure-data.R", "Minimum detectable prevalence maps and summaries", "Full manuscript run",
  "stochastic_minimum_detectable_prevalence.csv", 500L, 123L, "scripts/prepare-public-figure-data.R", "Revised stochastic minimum-detection-rate boxplots", "Full manuscript run",
  "daily_versus_hourly_comparison.csv", 500L, 123L, "scripts/summarise-hourly-sensitivity.R", "Impact of diurnal hydraulic variability / hydraulic aggregation sensitivity", "Full manuscript run; summarized raw restricted hourly hydraulics",
  "threshold_classification_diagnostic.csv", 1L, 123L, "scripts/calc-temporal-threshold-diagnostic.R", "Supplementary threshold-classification diagnostic", "Full manuscript midpoint-threshold diagnostic; one deterministic threshold realization",
  "particle_association_sensitivity.csv", 500L, 123L, "calc-loss-stochastic.R; scripts/prepare-public-particle-figure-data.R", "Particle-class pipe/path sensitivity figures", "Full manuscript run",
  "rna_partition_sensitivity.csv", 500L, 123L, "calc-loss-stochastic.R; scripts/prepare-public-figure-data.R", "Particle-association/RNA-distribution sensitivity figures", "Full manuscript run"
)
write_csv(metadata, file.path(output_dir, "output_metadata.csv"))

# Central column-level data dictionary. Units are explicit and shared metadata
# are repeated per output so each CSV can be interpreted independently.
dictionary_rows <- list()
add_dictionary <- function(file, data, unit_map = character()) {
  dictionary_rows[[length(dictionary_rows) + 1]] <<- tibble(
    file = file,
    column = names(data),
    description = gsub("_", " ", names(data)),
    units = ifelse(names(data) %in% names(unit_map), unit_map[names(data)], "category or identifier")
  )
}
proportion_units <- c(
  bin_left_proportion = "proportion (0-1)", bin_right_proportion = "proportion (0-1)",
  bin_mid_proportion = "proportion (0-1)", mean_loss_proportion = "proportion (0-1)",
  median_loss_proportion = "proportion (0-1)", flow_path_count = "count",
  observation_count = "count", population = "persons",
  population_equivalent_loss_persons = "persons", n_paths = "count",
  median_of_path_medians_percent = "percent", mean_of_path_medians_percent = "percent",
  q25_of_path_medians_percent = "percent", q75_of_path_medians_percent = "percent",
  effective_surveillance_coverage_percent = "percent"
)
add_dictionary("flow_path_loss_distributions.csv", flow_path, proportion_units)
add_dictionary("fsa_mean_loss_by_sample_fraction.csv", fsa_loss, proportion_units)
add_dictionary("population_equivalent_loss.csv", population_equivalent, proportion_units)
add_dictionary("effective_surveillance_coverage.csv", effective_coverage, proportion_units)
add_dictionary("stochastic_effective_surveillance_coverage.csv", stochastic_coverage,
               c(simulation = "simulation index", effective_surveillance_coverage_percent = "percent"))
add_dictionary("minimum_detectable_prevalence.csv", minimum_prevalence,
               setNames(rep("percent", sum(grepl("percent$", names(minimum_prevalence)))),
                        names(minimum_prevalence)[grepl("percent$", names(minimum_prevalence))]))
add_dictionary("stochastic_minimum_detectable_prevalence.csv", stochastic_prevalence,
               c(simulation = "simulation index", minimum_detectable_prevalence_percent = "percent"))
add_dictionary("daily_versus_hourly_comparison.csv", daily_hourly, proportion_units)
threshold <- read_csv(file.path(output_dir, "threshold_classification_diagnostic.csv"), show_col_types = FALSE)
add_dictionary("threshold_classification_diagnostic.csv", threshold,
               setNames(ifelse(grepl("^pct_|^median_F|^q[257]+_F|^median_delta|^q[257]+_delta", names(threshold)), "percent or fraction as named", "category or count"), names(threshold)))
add_dictionary("particle_association_sensitivity.csv", particle, proportion_units)
add_dictionary("rna_partition_sensitivity.csv", rna_partition, proportion_units)
write_csv(bind_rows(dictionary_rows), file.path(output_dir, "data_dictionary.csv"))

message("Published summarized manuscript outputs to ", output_dir)
