
#########
#########  TEMPORAL THRESHOLD-EXPOSURE DIAGNOSTIC
#########  Sensitivity analysis for sub-daily settling/resuspension hydraulics
#########
#########  Purpose:
#########    This script compares hourly hydraulic threshold exposure with the
#########    daily-mean shear-stress classification currently used in the fate-model
#########    manuscript pipeline. It is designed as a reviewer-response diagnostic
#########    for whether daily averaging masks sub-daily settling/resuspension.
#########
#########  Main outputs:
#########    out/temporal-threshold-diagnostic/wwtp_particle_threshold_summary.csv
#########    out/temporal-threshold-diagnostic/wwtp_particle_daily_vs_hourly_exposure.csv
#########    out/temporal-threshold-diagnostic/wwtp_particle_text_summary.csv
#########    out/temporal-threshold-diagnostic/figures/
#########      appendix_hydraulics_and_threshold_classification_day3.{png,pdf}
#########
#########  Expected hourly hydraulic input:
#########    The standardized output of scripts/prepare-hourly-path-states.R:
#########      data/hourly-path-states/YYYYMMDD/
#########        dwf_path_states_YYYYMMDD_HHMM.csv
#########    The files must contain the same hydraulic columns used by the main
#########    pipeline, especially:
#########      reference_id or node_id
#########      reference_ids or path_nodes_ids
#########      path_target or wwtp
#########      conduit_ss       # bed shear stress, Pa or N/m2
#########      conduit_hrt      # conduit HRT, hours
#########      depth            # optional, not used in this diagnostic
#########
#########  Usage from repo root:
#########    Rscript scripts/calc-temporal-threshold-diagnostic.R
#########
#########  Optional command line usage:
#########    Rscript scripts/calc-temporal-threshold-diagnostic.R data/hourly-path-states out/temporal-threshold-diagnostic
#########
#########  Optional environment variables:
#########    THRESHOLD_MODE=midpoint or stochastic
#########    NUM_SIM=500
#########    SEED=123
#########
#########    THRESHOLD_MODE controls the threshold grid. `midpoint` uses one
#########    midpoint per particle-class range. `stochastic` draws uniform
#########    thresholds from the same ranges as the main Monte Carlo model.
#########    NUM_SIM is used only in stochastic mode. SEED makes those draws
#########    reproducible.
#########
#########  Recommendation for manuscript appendix:
#########    Use THRESHOLD_MODE=midpoint first. It is transparent and fast.
#########    Use THRESHOLD_MODE=stochastic only if you want to propagate the
#########    sampled threshold distributions from the stochastic model.
#########
#########

suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
  library(fs)
  library(patchwork)
})

ggplot2::theme_set(theme_bw())

# -------------------------------------------------------------------------
# 0. USER SETTINGS
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Edit this path to your folder containing the 24 hourly files for each day.
# If this path is a single CSV, the script reads it as a long-format hourly file.
hourly_input <- ifelse(length(args) >= 1, args[1], "data/hourly-path-states")

# Output directory.
out_dir <- ifelse(length(args) >= 2, args[2], "out/temporal-threshold-diagnostic")
fig_dir <- file.path(out_dir, "figures")

dir_create(out_dir, recurse = TRUE)
dir_create(fig_dir, recurse = TRUE)

# Threshold mode:
#   midpoint   = use midpoint of the class-specific uniform ranges; fast and clear.
#   stochastic = sample ss.crt.set and ss.crt.res from the same ranges used in

#                the stochastic model; can be large because it multiplies hourly
#                hydraulics by particle class and simulation.
threshold_mode <- Sys.getenv("THRESHOLD_MODE", unset = "midpoint")
threshold_mode <- tolower(threshold_mode)

num.sim <- as.integer(Sys.getenv("NUM_SIM", unset = "500"))
seed <- as.integer(Sys.getenv("SEED", unset = "123"))

if (!threshold_mode %in% c("midpoint", "stochastic")) {
  stop("THRESHOLD_MODE must be 'midpoint' or 'stochastic'.")
}
if (is.na(num.sim) || num.sim < 1) {
  stop("NUM_SIM must be a positive integer.")
}
if (is.na(seed)) {
  stop("SEED must be an integer.")
}

# If TRUE, the script saves pipe-day-particle results. For large stochastic runs,
# this file can be large. The summary tables are always saved.
save_pipe_level_table <- FALSE

# A small tolerance avoids classification noise from numerical precision.
eps <- 1e-12

message("\nTemporal threshold-exposure diagnostic")
message("Input: ", hourly_input)
message("Output: ", out_dir)
message("Threshold mode: ", threshold_mode)
message("Number of threshold simulations if stochastic: ", num.sim)

# -------------------------------------------------------------------------
# 1. SOURCE EXISTING PIPELINE FUNCTIONS WHEN AVAILABLE
# -------------------------------------------------------------------------

if (file.exists("utils/utils.R")) {
  source("utils/utils.R")
} else if (file.exists("utils.R")) {
  source("utils.R")
} else {
  warning("Could not find utils/utils.R or utils.R. The script will use local cleaning helpers only.")
}

if (!file.exists("parameters-stochastic.R")) {
  stop("Run from the repository root; parameters-stochastic.R was not found.")
}
source("parameters-stochastic.R")

# -------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

safe_first <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

infer_datetime_from_file <- function(filepath) {
  fname <- path_file(filepath)
  parent <- path_file(path_dir(filepath))
  text <- paste(fname, parent, sep = " ")

  # Preferred pattern for hourly path-state files:
  #   dwf_path_states_YYYYMMDD_HHMM.csv
  # Example:
  #   dwf_path_states_20010101_0100.csv
  x <- stringr::str_match(text, "(\\d{8})[_-](\\d{2})(\\d{2})")
  if (!is.na(x[1, 1])) {
    return(tibble(
      date = as.Date(lubridate::ymd(x[1, 2])),
      hour = as.integer(x[1, 3]),
      minute = as.integer(x[1, 4])
    ))
  }

  tibble(
    date = as.Date(NA),
    hour = NA_integer_,
    minute = NA_integer_
  )
}

infer_date_from_file <- function(filepath) {


  # First try the hourly filename pattern YYYYMMDD_HHMM.
  dt <- infer_datetime_from_file(filepath)
  if (!is.na(dt$date[1])) return(dt$date[1])

  fname <- path_file(filepath)
  parent <- path_file(path_dir(filepath))
  text <- paste(fname, parent, sep = " ")

  # YYYY-MM-DD or YYYY_MM_DD
  x <- str_extract(text, "\\d{4}[-_]\\d{2}[-_]\\d{2}")
  if (!is.na(x)) return(as.Date(str_replace_all(x, "_", "-")))

  # YYYYMMDD
  x <- str_extract(text, "\\d{8}")
  if (!is.na(x)) return(as.Date(ymd(x)))

  as.Date(NA)
}

infer_day_index_from_file <- function(filepath, fallback_index = NA_integer_) {
  fname <- path_file(filepath)
  parent <- path_file(path_dir(filepath))
  text <- paste(fname, parent, sep = " ")

  x <- str_match(tolower(text), "day[^0-9]*([0-9]{1,3})")[, 2]
  if (!is.na(x)) return(as.integer(x))

  if (!is.na(fallback_index)) return(floor((fallback_index - 1) / 24) + 1)

  NA_integer_
}

infer_hour_from_file <- function(filepath, fallback_index = NA_integer_) {

  # First try the hourly filename pattern YYYYMMDD_HHMM.
  # For example, 20010101_0100 gives hour = 1.
  dt <- infer_datetime_from_file(filepath)
  if (!is.na(dt$hour[1])) return(dt$hour[1])

  fname <- path_file(filepath)
  fname_low <- tolower(fname)

  # hour03, hour_03, hr03, h03
  x <- str_match(fname_low, "(?:hour|hr|h)[^0-9]*([0-9]{1,2})")[, 2]
  if (!is.na(x)) return(as.integer(x))

  # _03.csv, -03.csv, space03.csv at the end of the filename
  x <- str_match(fname_low, "(?:_|-|\\s)([0-9]{1,2})(?:\\.[a-z0-9]+)$")[, 2]
  if (!is.na(x)) return(as.integer(x))

  if (!is.na(fallback_index)) return((fallback_index - 1) %% 24)

  NA_integer_
}

infer_minute_from_file <- function(filepath) {
  dt <- infer_datetime_from_file(filepath)
  dt$minute[1]
}

sort_hourly_files_chronologically <- function(files) {

  # Sort by parsed date/hour/minute when filenames contain YYYYMMDD_HHMM.
  # Important: infer_date_from_file() returns a Date, so we convert to character
  # before making the audit tibble. This avoids the as.Date.numeric(origin) error
  # caused by unlisting Date objects.
  file_index <- tibble(file = as.character(files)) %>%
    mutate(
      parsed_date_chr = purrr::map_chr(file, function(f) {
        x <- infer_date_from_file(f)
        if (is.na(x)) NA_character_ else as.character(x)
      }),
      parsed_date = as.Date(parsed_date_chr),
      parsed_hour = purrr::map_int(file, function(f) {
        infer_hour_from_file(f, fallback_index = NA_integer_)
      }),
      parsed_minute = purrr::map_int(file, infer_minute_from_file),
      sort_date = dplyr::coalesce(parsed_date, as.Date("9999-12-31")),
      sort_hour = dplyr::coalesce(parsed_hour, 999L),

      sort_minute = dplyr::coalesce(parsed_minute, 999L)
    ) %>%
    arrange(sort_date, sort_hour, sort_minute, file)

  file_index$file
}

standardize_hydraulic_columns <- function(dat) {

  # Harmonize common names without forcing a specific InfoWorks export format.
  if (!"node_id" %in% names(dat) && "reference_id" %in% names(dat)) {
    dat <- dat %>% rename(node_id = reference_id)
  }

  if (!"path_nodes_ids" %in% names(dat) && "reference_ids" %in% names(dat)) {
    dat <- dat %>% rename(path_nodes_ids = reference_ids)
  }

  if (!"us_node" %in% names(dat) && "us_node_id" %in% names(dat)) {
    dat <- dat %>% rename(us_node = us_node_id)
  }

  if (!"ds_node" %in% names(dat) && "ds_node_id" %in% names(dat)) {
    dat <- dat %>% rename(ds_node = ds_node_id)
  }

  # Infer WWTP/catchment from path_target, matching the existing clean() logic.
  if (!"wwtp" %in% names(dat)) {
    if ("path_target" %in% names(dat)) {
      dat <- dat %>%
        mutate(wwtp = case_when(
          grepl("^NEWPCC", path_target) ~ "north",
          grepl("^S-PL", path_target) ~ "south",
          TRUE ~ "west"
        ))
    } else {
      stop("Cannot identify WWTP/catchment. Add either 'wwtp' or 'path_target' to the hourly hydraulic files.")
    }
  }

  # Normalize case to match current pipeline outputs.
  dat <- dat %>%
    mutate(wwtp = tolower(as.character(wwtp)))

  required_cols <- c("node_id", "wwtp", "conduit_ss", "conduit_hrt")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required hydraulic columns: ", paste(missing_cols, collapse = ", "))
  }

  dat
}

clean_hourly_flow <- function(dat) {

  dat <- standardize_hydraulic_columns(dat)

  # Apply equivalent filters to the main clean() function where possible:
  #   conduit_hrt > 0, path_hrt < 30, non-blank path_target.
  dat <- dat %>%
    filter(!is.na(conduit_hrt), conduit_hrt > 0)

  if ("path_hrt" %in% names(dat)) {
    dat <- dat %>% filter(is.na(path_hrt) | path_hrt < 30)
  }

  if ("path_target" %in% names(dat)) {
    dat <- dat %>% filter(!is.na(path_target), path_target != "")
  }

  dat
}

read_one_hourly_file <- function(filepath, file_index = NA_integer_) {

  message("Reading ", filepath)

  dat <- suppressMessages(readr::read_csv(filepath, show_col_types = FALSE))

  dat <- clean_hourly_flow(dat)


  # For hourly path-state exports, date/hour are expected in the filename:
  #   dwf_path_states_YYYYMMDD_HHMM.csv
  # Use filename-derived values when available because the original file may
  # not include an hour column.
  date_from_file <- infer_date_from_file(filepath)
  hour_from_file <- infer_hour_from_file(filepath, fallback_index = file_index)
  minute_from_file <- infer_minute_from_file(filepath)

  if (!"date" %in% names(dat) || all(is.na(dat$date))) {
    dat$date <- date_from_file
  } else {
    dat$date <- as.Date(dat$date)
    if (!is.na(date_from_file)) {
      dat$date_from_file <- date_from_file
    }
  }

  if (!"day_index" %in% names(dat)) {
    dat$day_index <- infer_day_index_from_file(filepath, fallback_index = file_index)
  }

  if (!"hour" %in% names(dat) || all(is.na(dat$hour))) {
    dat$hour <- hour_from_file
  } else {
    dat$hour <- as.integer(dat$hour)

    # If filename contains a valid hour, keep a copy of the original and use
    # the filename-derived hour. This protects against stale/fixed hour fields.
    if (!is.na(hour_from_file)) {
      dat$hour_original <- dat$hour
      dat$hour <- hour_from_file
    }
  }

  if (!"minute" %in% names(dat) || all(is.na(dat$minute))) {
    dat$minute <- minute_from_file
  }

  dat$source_file <- path_file(filepath)

  dat
}

read_hourly_hydraulics <- function(input_path) {

  if (file.exists(input_path) && !dir_exists(input_path)) {
    message("Reading one long-format hourly file.")
    dat <- suppressMessages(readr::read_csv(input_path, show_col_types = FALSE))
    dat <- clean_hourly_flow(dat)

    if (!"day_index" %in% names(dat)) {
      if ("date" %in% names(dat)) {
        dat <- dat %>% mutate(day_index = as.integer(as.factor(as.Date(date))))
      } else {
        stop("A single long-format file must include either 'day_index' or 'date'.")
      }
    }

    if (!"hour" %in% names(dat)) {
      stop("A single long-format file must include an 'hour' column.")
    }

    if (!"date" %in% names(dat)) {
      dat$date <- as.Date(NA)
    } else {
      dat$date <- as.Date(dat$date)
    }

    dat$source_file <- path_file(input_path)
    return(dat)
  }

  if (!dir_exists(input_path)) {
    stop("Hourly input path does not exist: ", input_path)
  }

  hourly_files <- dir_ls(input_path, regexp = "\\.csv$", recurse = TRUE)
  hourly_files <- hourly_files[
    !grepl(

      "join_report|missing_pipe_ids|report",
      basename(hourly_files),
      ignore.case = TRUE
    )
  ]

  hourly_files <- sort_hourly_files_chronologically(hourly_files)

  if (length(hourly_files) == 0) {
    stop("No CSV files found in hourly input directory: ", input_path)
  }

  message("Found ", length(hourly_files), " hourly CSV files.")

  file_parse_audit <- tibble(
    file_index = seq_along(hourly_files),
    filepath = as.character(hourly_files),
    source_file = path_file(hourly_files)
  ) %>%
    mutate(
      parsed_date_chr = purrr::map_chr(filepath, function(f) {
        x <- infer_date_from_file(f)
        if (is.na(x)) NA_character_ else as.character(x)
      }),
      parsed_date = as.Date(parsed_date_chr),
      parsed_hour = purrr::map_int(filepath, function(f) {
        infer_hour_from_file(f, fallback_index = NA_integer_)
      }),
      parsed_minute = purrr::map_int(filepath, infer_minute_from_file)
    ) %>%
    dplyr::select(file_index, source_file, parsed_date, parsed_hour, parsed_minute, filepath)

  write_csv(file_parse_audit, file.path(out_dir, "hourly_file_parse_audit.csv"))

  dat <- map2_dfr(hourly_files, seq_along(hourly_files), read_one_hourly_file)

  # Recompute day_index globally from parsed dates so that all files from the
  # same date get the same day_index. This is safer than assigning day_index
  # inside each file independently.
  if ("date" %in% names(dat) && any(!is.na(dat$date))) {
    date_lookup <- tibble(
      date = sort(unique(dat$date[!is.na(dat$date)])),
      day_index_from_date = seq_along(sort(unique(dat$date[!is.na(dat$date)])))
    )

    dat <- dat %>%
      dplyr::select(-any_of("day_index")) %>%
      left_join(date_lookup, by = "date") %>%
      rename(day_index = day_index_from_date)
  }

  dat
}

# Class-specific ranges come directly from the main stochastic model.
get_threshold_ranges <- function() {
  if (!exists("get_ss.crt.set") || !exists("get_ss.crt.res")) {
    stop("Threshold functions are missing from parameters-stochastic.R.")
  }
  set_list <- get_ss.crt.set()
  res_list <- get_ss.crt.res()
  tibble(
    part.class = 1:8,
    ss.crt.set.min = map_dbl(set_list, 1),
    ss.crt.set.max = map_dbl(set_list, 2),
    ss.crt.res.min = map_dbl(res_list, 1),
    ss.crt.res.max = map_dbl(res_list, 2)
  )
}

validate_threshold_ranges <- function(ranges) {
  expected <- tibble(
    part.class = 1:8,
    ss.crt.set.min = c(0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.140, 0.170),
    ss.crt.set.max = c(0.004, 0.009, 0.019, 0.040, 0.090, 0.130, 0.160, 0.200),
    ss.crt.res.min = c(0.003, 0.006, 0.020, 0.100, 0.180, 0.300, 0.600, 1.400),
    ss.crt.res.max = c(0.006, 0.010, 0.090, 0.170, 0.290, 0.500, 1.300, 2.000)
  )
  if (!isTRUE(all.equal(as.data.frame(ranges), as.data.frame(expected), tolerance = 0))) {
    stop("Threshold ranges no longer match the parameter definitions used in the manuscript.")
  }
  invisible(TRUE)
}

make_threshold_grid <- function(mode = "midpoint", ns = 500, seed = 123) {

  ranges <- get_threshold_ranges()
  validate_threshold_ranges(ranges)

  if (mode == "midpoint") {
    return(
      ranges %>%
        transmute(
          sim = 1L,
          part.class,
          ss.crt.set = (ss.crt.set.min + ss.crt.set.max) / 2,
          ss.crt.res = (ss.crt.res.min + ss.crt.res.max) / 2,
          threshold_mode = "midpoint"
        )
    )
  }

  if (mode == "stochastic") {
    set.seed(seed)

    return(
      ranges %>%
        tidyr::crossing(sim = 1:ns) %>%
        group_by(part.class) %>%
        mutate(
          ss.crt.set = runif(n(), ss.crt.set.min, ss.crt.set.max),
          ss.crt.res = runif(n(), ss.crt.res.min, ss.crt.res.max),
          threshold_mode = "stochastic"
        ) %>%
        ungroup() %>%
        dplyr::select(sim, part.class, ss.crt.set, ss.crt.res, threshold_mode)
    )
  }

  stop("threshold_mode must be either 'midpoint' or 'stochastic'.")
}

classify_daily_vs_hourly <- function(daily_active, hourly_fraction, process_label) {

  # daily_active is TRUE/FALSE, hourly_fraction is 0 to 1.
  case_when(
    daily_active == FALSE & hourly_fraction <= eps ~
      paste0("consistent_never_", process_label),

    daily_active == FALSE & hourly_fraction > eps ~
      paste0("daily_misses_hourly_", process_label),

    daily_active == TRUE & hourly_fraction >= 1 - eps ~
      paste0("consistent_all_day_", process_label),

    daily_active == TRUE & hourly_fraction < 1 - eps ~
      paste0("daily_assumes_all_day_but_hourly_partial_", process_label),

    TRUE ~ "unclassified"
  )
}

# -------------------------------------------------------------------------
# HELPER LABELS FOR FIGURES
# -------------------------------------------------------------------------

wwtp_labeller <- ggplot2::labeller(
  wwtp = c(
    north = "North",
    south = "South",
    west  = "West"
  )
)

classification_levels <- c(
  "Consistent all-day set/resus",
  "Consistent never set/resus",
  "Daily misses hourly",
  "Daily overcalls all-day"
)

base_axis_theme <- theme(
  axis.title = element_text(face = "bold"),

  axis.text = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold")
)

# -------------------------------------------------------------------------
# 3A. LOAD AND PREPARE HOURLY HYDRAULIC DATA
# -------------------------------------------------------------------------

hyd_hourly_raw <- read_hourly_hydraulics(hourly_input)

message("\nRows before pipe-hour de-duplication: ", nrow(hyd_hourly_raw))

# Keep one hydraulic record per pipe / WWTP / day / hour.
# If there are duplicate records due to path-level repetition, average the
# hydraulic states. This prevents over-counting duplicated pipe rows.
hyd_hourly <- hyd_hourly_raw %>%
  mutate(
    node_id = as.character(node_id),
    wwtp = tolower(as.character(wwtp)),
    day_index = as.integer(day_index),
    hour = as.integer(hour),
    conduit_ss = as.numeric(conduit_ss),
    conduit_hrt = as.numeric(conduit_hrt),
    date = as.Date(date)
  ) %>%
  filter(
    !is.na(node_id),
    !is.na(wwtp),
    !is.na(day_index),
    !is.na(hour),
    !is.na(conduit_ss),
    conduit_ss >= 0,
    !is.na(conduit_hrt),
    conduit_hrt > 0
  ) %>%
  group_by(node_id, wwtp, date, day_index, hour) %>%
  summarise(
    conduit_ss = mean(conduit_ss, na.rm = TRUE),
    conduit_hrt = mean(conduit_hrt, na.rm = TRUE),
    n_duplicate_rows = n(),
    .groups = "drop"
  )

message("Rows after pipe-hour de-duplication: ", nrow(hyd_hourly))
message("Unique pipes: ", n_distinct(hyd_hourly$node_id))
message("WWTPs: ", paste(sort(unique(hyd_hourly$wwtp)), collapse = ", "))
message("Days: ", paste(sort(unique(hyd_hourly$day_index)), collapse = ", "))

# Save a quick audit table.
hourly_audit <- hyd_hourly %>%
  group_by(wwtp, day_index) %>%
  summarise(
    n_pipes = n_distinct(node_id),
    n_hours = n_distinct(hour),
    min_hour = min(hour, na.rm = TRUE),
    max_hour = max(hour, na.rm = TRUE),
    median_tau = median(conduit_ss, na.rm = TRUE),
    q25_tau = quantile(conduit_ss, 0.25, na.rm = TRUE),
    q75_tau = quantile(conduit_ss, 0.75, na.rm = TRUE),
    median_conduit_hrt = median(conduit_hrt, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(hourly_audit, file.path(out_dir, "hourly_hydraulic_input_audit.csv"))

# Stop early if the parsed hourly data does not contain 24 distinct hours
# for each WWTP-day. This prevents interpreting a one-hour-per-day dataset
# as a true sub-daily diagnostic.
missing_hour_groups <- hourly_audit %>%
  filter(n_hours < 24)

if (nrow(missing_hour_groups) > 0) {
  print(missing_hour_groups, n = Inf)
  stop(
    "Hourly parsing check failed: at least one WWTP-day has fewer than 24 hours. ",
    "Check out/temporal-threshold-diagnostic/hourly_file_parse_audit.csv and filename date/hour patterns."
  )
}


# Validate that the hourly input really contains sub-daily information.
# The diagnostic is not meaningful if only one hour is present per day.
expected_hours_per_day <- as.integer(Sys.getenv("EXPECTED_HOURS_PER_DAY", unset = "24"))

bad_hourly_coverage <- hourly_audit %>%
  filter(n_hours < expected_hours_per_day)

if (nrow(bad_hourly_coverage) > 0) {
  print(bad_hourly_coverage, n = 100)
  stop(
    "Hourly input does not contain the expected number of hours per WWTP/day. ",
    "Check out/hourly_file_parse_audit.csv and hourly_hydraulic_input_audit.csv. ",
    "For filenames like dwf_path_states_20010101_0100.csv, the script now parses ",
    "hour from the HHMM part of the filename."
  )
}

# -------------------------------------------------------------------------
# FIGURE: Hour-of-day shear-stress distribution by day and WWTP
# -------------------------------------------------------------------------

# Check how many days and hours are available
hyd_hourly %>%
  group_by(wwtp, day_index) %>%
  summarise(
    n_hours = n_distinct(hour),
    min_hour = min(hour, na.rm = TRUE),
    max_hour = max(hour, na.rm = TRUE),
    n_conduits = n_distinct(node_id),
    .groups = "drop"
  ) %>%
  print(n = 100)

# Prepare plotting data
hyd_hourly_hour_of_day_boxplot <- hyd_hourly %>%
  mutate(
    day_index = factor(day_index),
    hour = as.integer(hour),
    hour_factor = factor(hour, levels = sort(unique(hour))),
    
    # Avoid problems with log10 scale if any shear stress values are zero
    conduit_ss_plot = pmax(conduit_ss, 1e-8)
  )


# -------------------------------------------------------------------------
# 4. DAILY-MEAN HYDRAULIC CLASSIFICATION
# -------------------------------------------------------------------------

hyd_daily <- hyd_hourly %>%
  group_by(node_id, wwtp, date, day_index) %>%
  summarise(
    n_hours = n_distinct(hour),
    tau_daily_mean = mean(conduit_ss, na.rm = TRUE),
    tau_daily_HRT_weighted = weighted.mean(conduit_ss, w = conduit_hrt, na.rm = TRUE),
    HRT_day_sum_hours = sum(conduit_hrt, na.rm = TRUE),
    median_hourly_tau = median(conduit_ss, na.rm = TRUE),
    min_hourly_tau = min(conduit_ss, na.rm = TRUE),
    max_hourly_tau = max(conduit_ss, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(hyd_daily, file.path(out_dir, "daily_mean_hydraulic_states_from_hourly.csv"))

# -------------------------------------------------------------------------
# 5. HOURLY THRESHOLD EXPOSURE
# -------------------------------------------------------------------------

threshold_grid <- make_threshold_grid(
  mode = threshold_mode,
  ns = num.sim,
  seed = seed
)

write_csv(threshold_grid, file.path(out_dir, "threshold_grid_used.csv"))

message("\nCalculating hourly threshold exposure...")

# Cross hourly hydraulics with particle thresholds.

# For midpoint mode, this is N_pipe_hours x 8.
# For stochastic mode, this is N_pipe_hours x 8 x num.sim.
hourly_threshold <- hyd_hourly %>%
  tidyr::crossing(threshold_grid) %>%
  mutate(
    settle_active_hourly = conduit_ss <= ss.crt.set + eps,
    resusp_active_hourly = conduit_ss >= ss.crt.res - eps
  )

pipe_day_particle <- hourly_threshold %>%
  group_by(node_id, wwtp, date, day_index, sim, part.class, threshold_mode) %>%
  summarise(
    n_hours = n_distinct(hour),
    ss.crt.set = safe_first(ss.crt.set),
    ss.crt.res = safe_first(ss.crt.res),

    # Unweighted fraction of hours active.
    F_settle = mean(settle_active_hourly, na.rm = TRUE),
    F_resusp = mean(resusp_active_hourly, na.rm = TRUE),

    # HRT-weighted fraction of time/opportunity active.
    F_settle_HRT = sum(conduit_hrt * settle_active_hourly, na.rm = TRUE) /
      sum(conduit_hrt, na.rm = TRUE),
    F_resusp_HRT = sum(conduit_hrt * resusp_active_hourly, na.rm = TRUE) /
      sum(conduit_hrt, na.rm = TRUE),

    HRT_day_sum_hours = sum(conduit_hrt, na.rm = TRUE),
    tau_hourly_min = min(conduit_ss, na.rm = TRUE),
    tau_hourly_median = median(conduit_ss, na.rm = TRUE),
    tau_hourly_max = max(conduit_ss, na.rm = TRUE),
    tau_hourly_IQR = IQR(conduit_ss, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# 6. COMPARE HOURLY EXPOSURE WITH DAILY-MEAN CLASSIFICATION
# -------------------------------------------------------------------------

pipe_day_particle_compare <- pipe_day_particle %>%
  left_join(hyd_daily,
            by = c("node_id", "wwtp", "date", "day_index"),
            suffix = c("", "_daily")) %>%
  mutate(
    daily_settle_active_mean = tau_daily_mean <= ss.crt.set + eps,
    daily_resusp_active_mean = tau_daily_mean >= ss.crt.res - eps,

    daily_settle_active_HRT = tau_daily_HRT_weighted <= ss.crt.set + eps,
    daily_resusp_active_HRT = tau_daily_HRT_weighted >= ss.crt.res - eps,

    settle_classification_mean = classify_daily_vs_hourly(
      daily_settle_active_mean,
      F_settle,
      "settling"
    ),

    resusp_classification_mean = classify_daily_vs_hourly(
      daily_resusp_active_mean,
      F_resusp,
      "resuspension"
    ),

    settle_classification_HRT = classify_daily_vs_hourly(
      daily_settle_active_HRT,
      F_settle_HRT,
      "settling"
    ),

    resusp_classification_HRT = classify_daily_vs_hourly(
      daily_resusp_active_HRT,
      F_resusp_HRT,
      "resuspension"
    ),

    delta_settle_mean = F_settle - as.numeric(daily_settle_active_mean),
    delta_resusp_mean = F_resusp - as.numeric(daily_resusp_active_mean),

    delta_settle_HRT = F_settle_HRT - as.numeric(daily_settle_active_HRT),
    delta_resusp_HRT = F_resusp_HRT - as.numeric(daily_resusp_active_HRT)
  )


if (save_pipe_level_table) {
  write_csv(
    pipe_day_particle_compare,
    file.path(out_dir, "pipe_day_particle_threshold_exposure.csv")
  )
}

# -------------------------------------------------------------------------
# 7. REVIEWER-FACING SUMMARY TABLES
# -------------------------------------------------------------------------

wwtp_particle_threshold_summary <- pipe_day_particle_compare %>%
  group_by(wwtp, part.class, threshold_mode) %>%
  summarise(
    n_pipe_days = n_distinct(paste(node_id, day_index, sep = "__")),
    n_sims = n_distinct(sim),

    median_F_settle = median(F_settle, na.rm = TRUE),
    q25_F_settle = quantile(F_settle, 0.25, na.rm = TRUE),
    q75_F_settle = quantile(F_settle, 0.75, na.rm = TRUE),

    median_F_settle_HRT = median(F_settle_HRT, na.rm = TRUE),
    q25_F_settle_HRT = quantile(F_settle_HRT, 0.25, na.rm = TRUE),
    q75_F_settle_HRT = quantile(F_settle_HRT, 0.75, na.rm = TRUE),

    median_F_resusp = median(F_resusp, na.rm = TRUE),
    q25_F_resusp = quantile(F_resusp, 0.25, na.rm = TRUE),
    q75_F_resusp = quantile(F_resusp, 0.75, na.rm = TRUE),

    median_F_resusp_HRT = median(F_resusp_HRT, na.rm = TRUE),
    q25_F_resusp_HRT = quantile(F_resusp_HRT, 0.25, na.rm = TRUE),
    q75_F_resusp_HRT = quantile(F_resusp_HRT, 0.75, na.rm = TRUE),

    pct_daily_misses_hourly_settling_HRT =
      mean(settle_classification_HRT == "daily_misses_hourly_settling", na.rm = TRUE) * 100,

    pct_daily_assumes_all_day_but_hourly_partial_settling_HRT =
      mean(settle_classification_HRT == "daily_assumes_all_day_but_hourly_partial_settling", na.rm = TRUE) * 100,

    pct_daily_misses_hourly_resuspension_HRT =
      mean(resusp_classification_HRT == "daily_misses_hourly_resuspension", na.rm = TRUE) * 100,

    pct_daily_assumes_all_day_but_hourly_partial_resuspension_HRT =
      mean(resusp_classification_HRT == "daily_assumes_all_day_but_hourly_partial_resuspension", na.rm = TRUE) * 100,

    median_delta_settle_HRT = median(delta_settle_HRT, na.rm = TRUE),
    q25_delta_settle_HRT = quantile(delta_settle_HRT, 0.25, na.rm = TRUE),
    q75_delta_settle_HRT = quantile(delta_settle_HRT, 0.75, na.rm = TRUE),

    median_delta_resusp_HRT = median(delta_resusp_HRT, na.rm = TRUE),
    q25_delta_resusp_HRT = quantile(delta_resusp_HRT, 0.25, na.rm = TRUE),
    q75_delta_resusp_HRT = quantile(delta_resusp_HRT, 0.75, na.rm = TRUE),

    .groups = "drop"
  )

write_csv(
  wwtp_particle_threshold_summary,
  file.path(out_dir, "wwtp_particle_threshold_summary.csv")
)

# WWTP/day/particle table comparing daily 0/1 exposure with hourly exposure.
# This gives the cleanest numbers for the text because it is weighted by HRT.
wwtp_particle_daily_vs_hourly_exposure <- pipe_day_particle_compare %>%
  group_by(wwtp, day_index, sim, part.class, threshold_mode) %>%
  summarise(
    hourly_settle_exposure_HRT = weighted.mean(F_settle_HRT, w = HRT_day_sum_hours, na.rm = TRUE),
    daily_settle_exposure_HRT = weighted.mean(as.numeric(daily_settle_active_HRT), w = HRT_day_sum_hours, na.rm = TRUE),

    hourly_resusp_exposure_HRT = weighted.mean(F_resusp_HRT, w = HRT_day_sum_hours, na.rm = TRUE),
    daily_resusp_exposure_HRT = weighted.mean(as.numeric(daily_resusp_active_HRT), w = HRT_day_sum_hours, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    difference_settle_hourly_minus_daily = hourly_settle_exposure_HRT - daily_settle_exposure_HRT,
    difference_resusp_hourly_minus_daily = hourly_resusp_exposure_HRT - daily_resusp_exposure_HRT
  )

write_csv(

  wwtp_particle_daily_vs_hourly_exposure,
  file.path(out_dir, "wwtp_particle_daily_vs_hourly_exposure.csv")
)

# Compact table for paper text: average across the five dry-weather days and,
# if stochastic mode, across threshold simulations.
wwtp_particle_text_summary <- wwtp_particle_daily_vs_hourly_exposure %>%
  group_by(wwtp, part.class, threshold_mode) %>%
  summarise(
    mean_hourly_settle_exposure_HRT = mean(hourly_settle_exposure_HRT, na.rm = TRUE),
    mean_daily_settle_exposure_HRT = mean(daily_settle_exposure_HRT, na.rm = TRUE),
    mean_difference_settle_hourly_minus_daily = mean(difference_settle_hourly_minus_daily, na.rm = TRUE),

    mean_hourly_resusp_exposure_HRT = mean(hourly_resusp_exposure_HRT, na.rm = TRUE),
    mean_daily_resusp_exposure_HRT = mean(daily_resusp_exposure_HRT, na.rm = TRUE),
    mean_difference_resusp_hourly_minus_daily = mean(difference_resusp_hourly_minus_daily, na.rm = TRUE),

    .groups = "drop"
  )

write_csv(
  wwtp_particle_text_summary,
  file.path(out_dir, "wwtp_particle_text_summary.csv")
)

# -------------------------------------------------------------------------
# 8. FIGURES
# -------------------------------------------------------------------------

message("\nCreating figures...")

# time-series plot for shear stress showing temporal hydraulic dynamics 
# Boxplot:
# x-axis = hour of day
# fill = day
# each box = distribution across conduits
# facets = WWTP
p_tau_hour_of_day_box <- hyd_hourly_hour_of_day_boxplot %>%
  ggplot(
    aes(
      x = hour_factor,
      y = conduit_ss_plot,
      fill = day_index
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.alpha = 0.08,
    outlier.size = 0.25,
    linewidth = 0.25
  ) +
  facet_wrap(~ wwtp, ncol = 1) +
  scale_y_log10() +
  labs(
    x = "Hour of day",
    y = "Conduit shear stress, conduit_ss, log10 scale",
    fill = "Day",
    title = "Hourly distribution of conduit shear stress by WWTP",
    subtitle = "Each box summarizes conduit-level shear stress for a given hour and day"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(fig_dir, "wwtp_hour_of_day_shear_stress_boxplot_by_day.png"),
  p_tau_hour_of_day_box,
  width = 13,
  height = 8,
  dpi = 300
)

# Figure 1: Hourly settling exposure by WWTP and particle class.
p_settle_box <- pipe_day_particle_compare %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = F_settle_HRT)) +
  geom_boxplot(outlier.alpha = 0.15) +

  facet_wrap(~ wwtp) +
  labs(
    x = "Particle class",
    y = "HRT-weighted fraction of day below settling threshold",
    title = "Hourly settling-threshold exposure by WWTP and particle class"
  )

ggsave(
  file.path(fig_dir, "hourly_settling_exposure_boxplot.png"),
  p_settle_box,
  width = 9,
  height = 5,
  dpi = 300
)

# Figure 2: Hourly resuspension exposure by WWTP and particle class.
p_resusp_box <- pipe_day_particle_compare %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = F_resusp_HRT)) +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_wrap(~ wwtp) +
  labs(
    x = "Particle class",
    y = "HRT-weighted fraction of day above resuspension threshold",
    title = "Hourly resuspension-threshold exposure by WWTP and particle class"
  )

ggsave(
  file.path(fig_dir, "hourly_resuspension_exposure_boxplot.png"),
  p_resusp_box,
  width = 9,
  height = 5,
  dpi = 300
)

# Figure 3: Settling difference, hourly exposure minus daily classification.
p_settle_delta <- wwtp_particle_daily_vs_hourly_exposure %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = difference_settle_hourly_minus_daily)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_wrap(~ wwtp) +
  labs(
    x = "Particle class",
    y = "Hourly exposure minus daily-mean classification",
    title = "Effect of daily averaging on settling-threshold classification"
  )

ggsave(
  file.path(fig_dir, "settling_hourly_minus_daily.png"),
  p_settle_delta,
  width = 9,
  height = 5,
  dpi = 300
)

# Figure 4: Resuspension difference, hourly exposure minus daily classification.
p_resusp_delta <- wwtp_particle_daily_vs_hourly_exposure %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = difference_resusp_hourly_minus_daily)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_wrap(~ wwtp) +
  labs(
    x = "Particle class",
    y = "Hourly exposure minus daily-mean classification",
    title = "Effect of daily averaging on resuspension-threshold classification"
  )

ggsave(
  file.path(fig_dir, "resuspension_hourly_minus_daily.png"),
  p_resusp_delta,
  width = 9,
  height = 5,
  dpi = 300
)

# Figure 5: Settling classification categories.
p_settle_class <- pipe_day_particle_compare %>%
  count(wwtp, part.class, settle_classification_HRT) %>%

  group_by(wwtp, part.class) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    part.class = factor(part.class, levels = 1:8),
    classification_label = dplyr::recode(
      settle_classification_HRT,
      "consistent_all_day_settling" =
        "Consistent all-day set/resus",
      "consistent_never_settling" =
        "Consistent never set/resus",
      "daily_misses_hourly_settling" =
        "Daily misses hourly",
      "daily_assumes_all_day_but_hourly_partial_settling" =
        "Daily overcalls all-day"
    ),
    classification_label = factor(
      classification_label,
      levels = classification_levels
    )
  ) %>%
  ggplot(aes(x = part.class, y = percent, fill = classification_label)) +
  geom_col() +
  facet_wrap(~ wwtp, labeller = wwtp_labeller) +
  labs(
    x = "Particle class",
    y = "Settling Pipe-days (%)",
    fill = "Classification"
  ) +
  scale_fill_discrete(drop = FALSE) +
  theme_bw() +
  base_axis_theme +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "settling_classification_categories.png"),
  p_settle_class,
  width = 11,
  height = 6,
  dpi = 300
)

# Figure 6: Resuspension classification categories.
# -------------------------------------------------------------------------
# HELPER LABELS FOR FIGURES
# -------------------------------------------------------------------------

wwtp_labeller <- ggplot2::labeller(
  wwtp = c(
    north = "North",
    south = "South",
    west  = "West"
  )
)

classification_levels <- c(
  "Consistent all-day set/resus",
  "Consistent never set/resus",
  "Daily misses hourly",
  "Daily overcalls all-day"
)

base_axis_theme <- theme(
  axis.title = element_text(face = "bold"),
  axis.text = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold")
)

# Figure: Settling classification categories.
p_settle_class <- pipe_day_particle_compare %>%
  count(wwtp, part.class, settle_classification_HRT) %>%
  group_by(wwtp, part.class) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    part.class = factor(part.class, levels = 1:8),

    classification_label = dplyr::recode(
      settle_classification_HRT,
      "consistent_all_day_settling" =
        "Consistent all-day set/resus",
      "consistent_never_settling" =
        "Consistent never set/resus",
      "daily_misses_hourly_settling" =
        "Daily misses hourly",
      "daily_assumes_all_day_but_hourly_partial_settling" =
        "Daily overcalls all-day"
    ),
    classification_label = factor(
      classification_label,
      levels = classification_levels
    )
  ) %>%
  ggplot(aes(x = part.class, y = percent, fill = classification_label)) +
  geom_col() +
  facet_wrap(~ wwtp, labeller = wwtp_labeller) +
  labs(
    x = "Particle class",
    y = "Settling pipe-days (%)",
    fill = "Classification"
  ) +
  scale_fill_discrete(drop = FALSE) +
  theme_bw() +
  base_axis_theme +
  theme(
    legend.position = "bottom"
  )

# Figure: Resuspension classification categories.
p_resusp_class <- pipe_day_particle_compare %>%
  count(wwtp, part.class, resusp_classification_HRT) %>%
  group_by(wwtp, part.class) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    part.class = factor(part.class, levels = 1:8),
    classification_label = dplyr::recode(
      resusp_classification_HRT,
      "consistent_all_day_resuspension" =
        "Consistent all-day set/resus",
      "consistent_never_resuspension" =
        "Consistent never set/resus",
      "daily_misses_hourly_resuspension" =
        "Daily misses hourly",
      "daily_assumes_all_day_but_hourly_partial_resuspension" =
        "Daily overcalls all-day"
    ),
    classification_label = factor(
      classification_label,
      levels = classification_levels
    )
  ) %>%
  ggplot(aes(x = part.class, y = percent, fill = classification_label)) +
  geom_col() +
  facet_wrap(~ wwtp, labeller = wwtp_labeller) +
  labs(
    x = "Particle class",
    y = "Resuspension pipe-days (%)",
    fill = "Classification"
  ) +
  scale_fill_discrete(drop = FALSE) +
  theme_bw() +
  base_axis_theme +
  theme(
    legend.position = "bottom"
  )


ggsave(
  file.path(fig_dir, "resuspension_classification_categories.png"),
  p_resusp_class,
  width = 11,
  height = 6,
  dpi = 300
)

# Optional compact heatmap-style tables for visual manuscript supplement.

p_heat_settle <- wwtp_particle_text_summary %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = wwtp, fill = mean_difference_settle_hourly_minus_daily)) +
  geom_tile() +
  geom_text(aes(label = round(mean_difference_settle_hourly_minus_daily, 2)), size = 3) +
  labs(
    x = "Particle class",
    y = "WWTP",
    fill = "Hourly - daily",
    title = "Mean difference in settling exposure: hourly minus daily"
  )

ggsave(
  file.path(fig_dir, "settling_hourly_minus_daily_heatmap.png"),
  p_heat_settle,
  width = 8,
  height = 3.8,
  dpi = 300
)

p_heat_resusp <- wwtp_particle_text_summary %>%
  mutate(part.class = factor(part.class, levels = 1:8)) %>%
  ggplot(aes(x = part.class, y = wwtp, fill = mean_difference_resusp_hourly_minus_daily)) +
  geom_tile() +
  geom_text(aes(label = round(mean_difference_resusp_hourly_minus_daily, 2)), size = 3) +
  labs(
    x = "Particle class",
    y = "WWTP",
    fill = "Hourly - daily",
    title = "Mean difference in resuspension exposure: hourly minus daily"
  )

ggsave(
  file.path(fig_dir, "resuspension_hourly_minus_daily_heatmap.png"),
  p_heat_resusp,
  width = 8,
  height = 3.8,
  dpi = 300
)



# -------------------------------------------------------------------------
# FIGURE: One-day hourly shear-stress distribution by WWTP
# -------------------------------------------------------------------------

# Figure: One-day hourly shear-stress distribution by WWTP
focus_day <- 3

hyd_hourly_one_day_boxplot <- hyd_hourly %>%
  filter(day_index == focus_day) %>%
  mutate(
    hour = as.integer(hour),
    hour_factor = factor(hour, levels = 0:23),
    conduit_ss_plot = pmax(conduit_ss, 1e-8)
  )

p_tau_one_day_box <- hyd_hourly_one_day_boxplot %>%
  ggplot(aes(x = hour_factor, y = conduit_ss_plot)) +
  geom_boxplot(
    outlier.alpha = 0.05,
    outlier.size = 0.2,
    linewidth = 0.25
  ) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = 1),
    linewidth = 0.7,
    color = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.2,
    color = "black"
  ) +
  facet_wrap(~ wwtp, ncol = 1, labeller = wwtp_labeller) +
  scale_y_log10() +
  labs(

    x = "Hour of day",
    y = expression(bold("Conduit shear stress (N " * m^{-2} * ")"))
  )
  theme_bw() +
  base_axis_theme +
  theme(
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(fig_dir, paste0("wwtp_hour_of_day_shear_stress_boxplot_day", focus_day, "_with_median.png")),
  p_tau_one_day_box,
  width = 10,
  height = 8,
  dpi = 300
)


# -------------------------------------------------------------------------
# COMBINED APPENDIX FIGURE:
# Hydraulic dynamics + daily-vs-hourly threshold classification
# -------------------------------------------------------------------------

combined_appendix_figure <-
  p_tau_one_day_box /
  (p_settle_class | p_resusp_class) +
  plot_layout(
    heights = c(1.2, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  file.path(fig_dir, paste0("appendix_hydraulics_and_threshold_classification_day", focus_day, ".png")),
  combined_appendix_figure,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, paste0("appendix_hydraulics_and_threshold_classification_day", focus_day, ".pdf")),
  combined_appendix_figure,
  width = 10,
  height = 8
)




# -------------------------------------------------------------------------
# STACKED BAR FIGURES USING SIMPLE UNWEIGHTED HOURLY FRACTION
# Daily classification based on simple daily mean tau
# Hourly exposure based on fraction of hours active, not HRT-weighted exposure
# -------------------------------------------------------------------------

classification_levels <- c(

  "Consistent all-day set/resus",
  "Consistent never set/resus",
  "Daily misses hourly",
  "Daily overcalls all-day"
)

wwtp_labeller <- ggplot2::labeller(
  wwtp = c(
    north = "North",
    south = "South",
    west  = "West"
  )
)

base_axis_theme <- theme(
  axis.title = element_text(face = "bold"),
  axis.text = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold")
)

# -------------------------------------------------------------------------
# Settling classification categories: simple unweighted version
# -------------------------------------------------------------------------

p_settle_class_simple <- pipe_day_particle_compare %>%
  count(wwtp, part.class, settle_classification_mean) %>%
  group_by(wwtp, part.class) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    part.class = factor(part.class, levels = 1:8),
    classification_label = dplyr::recode(
      settle_classification_mean,
      "consistent_all_day_settling" =
        "Consistent all-day set/resus",
      "consistent_never_settling" =
        "Consistent never set/resus",
      "daily_misses_hourly_settling" =
        "Daily misses hourly",
      "daily_assumes_all_day_but_hourly_partial_settling" =
        "Daily overcalls all-day"
    ),
    classification_label = factor(
      classification_label,
      levels = classification_levels
    )
  ) %>%
  ggplot(aes(x = part.class, y = percent, fill = classification_label)) +
  geom_col() +
  facet_wrap(~ wwtp, labeller = wwtp_labeller) +
  labs(
    x = "Particle class",
    y = "Settling pipe-days (%)",
    fill = "Classification",
    title = "Settling classification based on simple hourly fraction"
  ) +
  scale_fill_discrete(drop = FALSE) +
  theme_bw() +
  base_axis_theme +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "settling_classification_categories_simple_unweighted.png"),
  p_settle_class_simple,
  width = 11,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# Resuspension classification categories: simple unweighted version
# -------------------------------------------------------------------------

p_resusp_class_simple <- pipe_day_particle_compare %>%
  count(wwtp, part.class, resusp_classification_mean) %>%
  group_by(wwtp, part.class) %>%

  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    part.class = factor(part.class, levels = 1:8),
    classification_label = dplyr::recode(
      resusp_classification_mean,
      "consistent_all_day_resuspension" =
        "Consistent all-day set/resus",
      "consistent_never_resuspension" =
        "Consistent never set/resus",
      "daily_misses_hourly_resuspension" =
        "Daily misses hourly",
      "daily_assumes_all_day_but_hourly_partial_resuspension" =
        "Daily overcalls all-day"
    ),
    classification_label = factor(
      classification_label,
      levels = classification_levels
    )
  ) %>%
  ggplot(aes(x = part.class, y = percent, fill = classification_label)) +
  geom_col() +
  facet_wrap(~ wwtp, labeller = wwtp_labeller) +
  labs(
    x = "Particle class",
    y = "Resuspension pipe-days (%)",
    fill = "Classification",
    title = "Resuspension classification based on simple hourly fraction"
  ) +
  scale_fill_discrete(drop = FALSE) +
  theme_bw() +
  base_axis_theme +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "resuspension_classification_categories_simple_unweighted.png"),
  p_resusp_class_simple,
  width = 11,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# Combined appendix-style figure: simple unweighted classification
# -------------------------------------------------------------------------

combined_appendix_figure_simple <-
  p_tau_one_day_box /
  (p_settle_class_simple | p_resusp_class_simple) +
  patchwork::plot_layout(
    heights = c(1.2, 1),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  file.path(fig_dir, paste0("appendix_hydraulics_and_threshold_classification_day", focus_day, "_simple_unweighted.png")),
  combined_appendix_figure_simple,
  width = 10,

  height = 8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, paste0("appendix_hydraulics_and_threshold_classification_day", focus_day, "_simple_unweighted.pdf")),
  combined_appendix_figure_simple,
  width = 10,
  height = 8
)
