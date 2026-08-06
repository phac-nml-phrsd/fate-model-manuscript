# ============================================================
# calc-loss-hourly-sensitivity.R
# ============================================================
# Driver script for hourly sensitivity analysis:
#
#   Version A: daily-averaged hydraulic states + existing fate model
#              Already calculated in the main manuscript pipeline.
#              Loaded here from:
#                out/sim_df_loss_total.rds
#
#   Version B: hourly-resolved hydraulic states.
#              Thresholds/losses are calculated at hourly scale first,
#              then aggregated to the 24-h composite scale.
#
# Public demonstration:
#   num.sim = 1 keeps the example run small. The manuscript results used
#   num.sim = 500 and a complete set of 24 hourly files for each day.
#   The real one-hour input excerpt distributed in data/hourly-data documents
#   the input format but cannot reproduce a 24-hour composite result.
#
# Run from repository root after creating hourly path-state files:
#   Rscript scripts/prepare-hourly-path-states.R
#   Rscript scripts/calc-loss-hourly-sensitivity.R
# ============================================================

suppressMessages({
  library(dplyr)
  library(purrr)
  library(furrr)
  library(future)
  library(readr)
  library(tibble)
})

source("utils/utils.R")
source("utils/utils_param_stochastic.R")
source("utils/utils_simulation_stochastic.R")
source("utils/utils_hourly_sensitivity.R")

# Avoid select() masking by other packages
select <- dplyr::select
all_of <- dplyr::all_of
any_of <- dplyr::any_of

# ---------------- User settings ----------------
num.sim <- 1 # public demonstration; use 500 for the manuscript analysis
seed <- 123

hourly_root <- file.path("data", "hourly-path-states")
out_dir <- file.path("out", "hourly-sensitivity")

# Existing Version A output from manuscript pipeline
version_A_total_file <- file.path("out", "sim_df_loss_total.rds")

# Existing path-ID file from manuscript pipeline
df_long_ids_file <- file.path("out", "df.long.ids.rds")

strict_missing <- FALSE

# The example date is intentionally hard-coded for the public workflow.
# A full manuscript reproduction used the corresponding complete daily set.
run_dates <- "20010101"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Important:
# df.long.ids is very large, so avoid multisession here.
num_cores <- 1
future::plan(sequential)

# ============================================================
# Parameters
# ============================================================

message("Loading stochastic and fixed parameters...")

stoch.params <- load_stoch_prms(
  num.sample = num.sim,
  seed = seed
)

parm.const <- get_const_params()

saveRDS(
  stoch.params,
  file = file.path(out_dir, "stoch.params.hourly_sensitivity.rds")
)

# ============================================================
# Load existing path IDs
# ============================================================

message("Loading existing df.long.ids from manuscript pipeline...")

if (!file.exists(df_long_ids_file)) {
  stop(
    "Could not find df.long.ids file: ",
    df_long_ids_file,
    "\nPlease check the filename in the out directory."
  )
}

df.long.ids <- readRDS(df_long_ids_file)

# ============================================================
# Hourly file index
# ============================================================

message("Indexing hourly path-state files...")

hourly_tbl <- parse_hourly_path_state_files(hourly_root)

available_dates <- sort(unique(as.character(hourly_tbl$date)))

if (is.null(run_dates)) {
  
  dates_to_run <- available_dates
  
} else {
  
  run_dates <- as.character(run_dates)
  
  missing_dates <- setdiff(run_dates, available_dates)
  
  if (length(missing_dates) > 0) {
    stop(
      "The following run_dates were not found in hourly path-state files: ",
      paste(missing_dates, collapse = ", "),
      "\nAvailable dates are: ",
      paste(available_dates, collapse = ", ")
    )
  }
  
  dates_to_run <- run_dates
}

message("Dates found: ", paste(available_dates, collapse = ", "))
message("Dates selected for hourly sensitivity run: ",
        paste(dates_to_run, collapse = ", "))

# ============================================================
# Main loop by selected day(s): Version B only
# ============================================================

all_hourly_total_B <- list()
all_hourly_solid_B <- list()
all_hourly_raw_B <- list()

for (day_i in dates_to_run) {
  
  message("==================================================")
  message("Processing day: ", day_i)
  message("==================================================")
  
  day_files_tbl <- hourly_tbl %>%
    dplyr::mutate(date = as.character(date)) %>%
    dplyr::filter(date == day_i) %>%
    dplyr::arrange(hour_index)
  
  if (nrow(day_files_tbl) == 0) {
    warning("No hourly files found for day: ", day_i)
    next
  }
  
  if (nrow(day_files_tbl) != 24) {
    warning(
      "Day ", day_i, " has ", nrow(day_files_tbl),
      " hourly files, not 24."
    )
  }
  
  # ----------------------------------------------------------
  # Version B: hourly-resolved thresholds and losses
  # ----------------------------------------------------------
  
  message("Running hourly-resolved threshold/loss calculations...")
  
  hourly_results <- purrr::map(
    1:max(stoch.params$sim),
    function(x) {
      
      message("  Simulation ", x, " of ", max(stoch.params$sim))
      
      parms.s <- stoch.params[stoch.params$sim == x, , drop = FALSE]
      total.parms <- cbind(parms.s, parm.const)
      
      simulate_hourly_loss_one_day_sim(
        hourly_file_tbl_day = day_files_tbl,
        total.parms = total.parms,
        df.long.ids = df.long.ids,
        sim = x,
        strict_missing = strict_missing
      )
    }
  )
  
  # ----------------------------------------------------------
  # Accept either:
  #   1) a data frame returned directly, or
  #   2) a list containing df.total / df.solid / df.hourly.raw
  # ----------------------------------------------------------
  
  if (is.data.frame(hourly_results[[1]])) {
    
    sim_df_loss_hourly_raw_B <- dplyr::bind_rows(hourly_results)
    
    saveRDS(
      sim_df_loss_hourly_raw_B,
      file = file.path(out_dir, paste0("sim_df_loss_hourly_raw_B_", day_i, ".rds"))
    )
    
    sim_df_loss_total_B <- aggregate_24h_composite_survival(
      sim_df_loss_hourly_raw_B
    ) %>%
      dplyr::mutate(
        date = day_i,
        version = "B_hourly_resolved"
      )
    
    if (exists("aggregate_24h_composite_solid")) {
      
      sim_df_loss_solid_B <- aggregate_24h_composite_solid(
        sim_df_loss_hourly_raw_B
      ) %>%
        dplyr::mutate(
          date = day_i,
          version = "B_hourly_resolved"
        )
      
    } else {
      
      sim_df_loss_solid_B <- sim_df_loss_hourly_raw_B %>%
        dplyr::mutate(
          date = day_i,
          version = "B_hourly_resolved"
        )
    }
    
  } else {
    
    sim_df_loss_total_B <- dplyr::bind_rows(
      purrr::map(hourly_results, "df.total")
    ) %>%
      dplyr::mutate(
        date = day_i,
        version = "B_hourly_resolved"
      )
    
    sim_df_loss_solid_B <- dplyr::bind_rows(
      purrr::map(hourly_results, "df.solid")
    ) %>%
      dplyr::mutate(
        date = day_i,
        version = "B_hourly_resolved"
      )
    
    if ("df.hourly.raw" %in% names(hourly_results[[1]])) {
      
      sim_df_loss_hourly_raw_B <- dplyr::bind_rows(
        purrr::map(hourly_results, "df.hourly.raw")
      )
      
    } else {
      
      sim_df_loss_hourly_raw_B <- tibble::tibble()
    }
    
    if (nrow(sim_df_loss_hourly_raw_B) > 0) {
      saveRDS(
        sim_df_loss_hourly_raw_B,
        file = file.path(out_dir, paste0("sim_df_loss_hourly_raw_B_", day_i, ".rds"))
      )
    }
  }
  
  # ----------------------------------------------------------
  # Save daily Version B outputs
  # ----------------------------------------------------------
  
  saveRDS(
    sim_df_loss_total_B,
    file = file.path(out_dir, paste0("sim_df_loss_total_B_hourly_", day_i, ".rds"))
  )
  
  saveRDS(
    sim_df_loss_solid_B,
    file = file.path(out_dir, paste0("sim_df_loss_solid_B_hourly_", day_i, ".rds"))
  )
  
  all_hourly_total_B[[day_i]] <- sim_df_loss_total_B
  all_hourly_solid_B[[day_i]] <- sim_df_loss_solid_B
  all_hourly_raw_B[[day_i]] <- sim_df_loss_hourly_raw_B
  
  rm(hourly_results)
  gc()
}

# ============================================================
# Save combined Version B outputs for selected day(s)
# ============================================================

message("Saving combined hourly sensitivity outputs...")

sim_df_loss_total_B_all <- dplyr::bind_rows(all_hourly_total_B)
sim_df_loss_solid_B_all <- dplyr::bind_rows(all_hourly_solid_B)
sim_df_loss_hourly_raw_B_all <- dplyr::bind_rows(all_hourly_raw_B)

saveRDS(
  sim_df_loss_total_B_all,
  file = file.path(out_dir, "sim_df_loss_total.rds")
)

saveRDS(
  sim_df_loss_solid_B_all,
  file = file.path(out_dir, "sim_df_loss_solid.rds")
)

saveRDS(
  sim_df_loss_hourly_raw_B_all,
  file = file.path(out_dir, "sim_df_loss_hourly_raw.rds")
)

readr::write_csv(
  sim_df_loss_total_B_all,
  file.path(out_dir, "sim_df_loss_total.csv")
)

readr::write_csv(
  sim_df_loss_solid_B_all,
  file.path(out_dir, "sim_df_loss_solid.csv")
)

# ============================================================
# Final comparison: Version A daily vs Version B hourly
# Compare only selected day(s), when date exists
# ============================================================

message("Loading existing Version A daily output from manuscript pipeline...")

if (!file.exists(version_A_total_file)) {
  stop(
    "Could not find Version A output file: ",
    version_A_total_file,
    "\nRun the main manuscript pipeline first to create out/sim_df_loss_total.rds."
  )
}

sim_df_loss_total_A_all <- readRDS(version_A_total_file)

if ("date" %in% names(sim_df_loss_total_A_all)) {
  
  sim_df_loss_total_A_all <- sim_df_loss_total_A_all %>%
    dplyr::mutate(
      date = as.character(.data$date),
      version = "A_daily_average"
    )
  
} else {
  
  warning(
    "Version A output does not contain a date column. ",
    "Adding date based on selected run_dates."
  )
  
  if (length(dates_to_run) == 1) {
    
    sim_df_loss_total_A_all <- sim_df_loss_total_A_all %>%
      dplyr::mutate(
        date = dates_to_run[1],
        version = "A_daily_average"
      )
    
  } else {
    
    sim_df_loss_total_A_all <- sim_df_loss_total_A_all %>%
      dplyr::mutate(
        version = "A_daily_average"
      )
  }
}

if ("date" %in% names(sim_df_loss_total_B_all)) {
  
  sim_df_loss_total_B_all <- sim_df_loss_total_B_all %>%
    dplyr::mutate(
      date = as.character(.data$date),
      version = "B_hourly_resolved"
    )
  
} else {
  
  sim_df_loss_total_B_all <- sim_df_loss_total_B_all %>%
    dplyr::mutate(
      date = dates_to_run[1],
      version = "B_hourly_resolved"
    )
}

if ("date" %in% names(sim_df_loss_total_A_all)) {
  
  sim_df_loss_total_A_compare <- sim_df_loss_total_A_all %>%
    dplyr::filter(.data$date %in% dates_to_run)
  
} else {
  
  warning(
    "Version A output does not contain a date column. ",
    "Comparison will use the full Version A dataset."
  )
  
  sim_df_loss_total_A_compare <- sim_df_loss_total_A_all
}

if ("date" %in% names(sim_df_loss_total_B_all)) {
  
  sim_df_loss_total_B_compare <- sim_df_loss_total_B_all %>%
    dplyr::filter(.data$date %in% dates_to_run)
  
} else {
  
  sim_df_loss_total_B_compare <- sim_df_loss_total_B_all
}

message("Comparing Version A daily-averaged losses with Version B hourly-resolved losses...")
message("Comparison dates: ", paste(dates_to_run, collapse = ", "))

comparison_A_vs_B_all <- compare_daily_vs_hourly(
  sim_df_loss_total_A_compare,
  sim_df_loss_total_B_compare
)

saveRDS(
  comparison_A_vs_B_all,
  file = file.path(out_dir, "comparison_A_vs_B_selected_days.rds")
)

readr::write_csv(
  comparison_A_vs_B_all,
  file.path(out_dir, "comparison_A_vs_B_selected_days.csv")
)

saveRDS(
  comparison_A_vs_B_all,
  file = file.path(out_dir, "comparison_A_vs_B_all_days.rds")
)

readr::write_csv(
  comparison_A_vs_B_all,
  file.path(out_dir, "comparison_A_vs_B_all_days.csv")
)


message("Done.")
message("Dates analyzed: ", paste(dates_to_run, collapse = ", "))
message("Hourly Version B total output: ", file.path(out_dir, "sim_df_loss_total.rds"))
message("Hourly Version B solid output: ", file.path(out_dir, "sim_df_loss_solid.rds"))
message("Daily vs hourly comparison: ", file.path(out_dir, "comparison_A_vs_B_selected_days.rds"))
