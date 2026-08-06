# ============================================================
# utils_hourly_sensitivity.R
# ============================================================
# Hourly sensitivity functions for fate-model-manuscript.
#
# Key principle:
#   Apply hydraulic thresholds at the hourly pipe scale first,
#   calculate path survival for each entry hour, then average
#   the 24 hourly surviving fractions for the composite sample.
#
# Expected workflow:
#   1) Run scripts/prepare-hourly-path-states.R to create:
#        data/hourly-path-states/YYYYMMDD/dwf_path_states_YYYYMMDD_HHMM.csv
#   2) Source this file from a new driver script.
#   3) Run Version A: daily averaged hydraulics with existing simulate_calc_loss().
#   4) Run Version B: hourly-resolved survival and 24-h composite averaging.
#
# Expected hourly path-state columns:
#   reference_id, path_target, reference_ids/path_nodes_ids, conduit_hrt,
#   conduit_ss, depth, conduit_av, and the other hydraulic columns retained
#   by scripts/prepare-hourly-path-states.R. A manuscript-scale day contains
#   24 files named dwf_path_states_YYYYMMDD_HHMM.csv.
# ============================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

parse_hourly_path_state_files <- function(hourly_root = file.path("data", "hourly-path-states")) {
  files <- list.files(
    hourly_root,
    pattern = "^dwf_path_states_\\d{8}_\\d{4}\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No hourly path-state files found in: ", hourly_root)
  }

  tibble(file = sort(files)) %>%
    mutate(
      fname = basename(file),
      date = str_match(fname, "dwf_path_states_(\\d{8})_(\\d{4})\\.csv")[, 2],
      hour = str_match(fname, "dwf_path_states_(\\d{8})_(\\d{4})\\.csv")[, 3],
      hour_index = as.integer(substr(hour, 1, 2))
    ) %>%
    filter(!is.na(date), !is.na(hour)) %>%
    arrange(date, hour_index)
}

read_path_state_csv <- function(file) {
  readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
}

# Daily mean hydraulic state, preserving the original path/topology columns.
# This produces the input for Version A, so it can be passed directly to clean()
# and then to the existing simulate_calc_loss().
make_daily_average_path_states <- function(day_files, id_col = "reference_id") {
  if (length(day_files) == 0) stop("day_files is empty.")

  hourly_df <- purrr::map_dfr(day_files, read_path_state_csv, .id = ".hour_file_id")

  if (!id_col %in% names(hourly_df)) {
    stop("Cannot find id_col = '", id_col, "' in hourly path-state files.")
  }

  hydraulic_cols_candidate <- c(
    "depth", "flow", "vol", "conduit_height", "conduit_width",
    "conduit_volume", "conduit_xs", "flow_xs", "flow_perimeter",
    "flow_area", "flow_volume", "flow_velocity", "residence_time",
    "conduit_av", "conduit_hrt", "conduit_hrt_inv", "conduit_ff", "conduit_ss"
  )

  hydraulic_cols <- intersect(hydraulic_cols_candidate, names(hourly_df))

  if (length(hydraulic_cols) == 0) {
    stop("No hydraulic columns were found to average.")
  }

  original_cols <- names(hourly_df)
  original_cols <- setdiff(original_cols, ".hour_file_id")
  meta_cols <- setdiff(original_cols, hydraulic_cols)

  daily_df <- hourly_df %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      across(all_of(meta_cols), ~ dplyr::first(.x)),
      across(all_of(hydraulic_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::select(all_of(original_cols))

  daily_df
}

get_f_solid_for_wwtp <- function(wwtp, total.parms) {
  w <- tolower(as.character(wwtp))
  dplyr::case_when(
    w %in% c("north", "n") ~ unique(total.parms$f.solid.N)[1],
    w %in% c("south", "s") ~ unique(total.parms$f.solid.S)[1],
    w %in% c("west",  "w") ~ unique(total.parms$f.solid.W)[1],
    TRUE ~ mean(c(
      unique(total.parms$f.solid.N)[1],
      unique(total.parms$f.solid.S)[1],
      unique(total.parms$f.solid.W)[1]
    ), na.rm = TRUE)
  )
}

check_hourly_path_join <- function(df.long.ids, pipe_df, keys = c("node_id", "part.class", "wwtp"), strict_missing = FALSE) {
  missing <- anti_join(
    df.long.ids %>% dplyr::select(all_of(keys)) %>% distinct(),
    pipe_df %>% dplyr::select(all_of(keys)) %>% distinct(),
    by = keys
  )

  if (nrow(missing) > 0) {
    msg <- paste0(
      "Hourly pipe/path join has ", nrow(missing),
      " missing pipe/class/WWTP combinations. ",
      "This usually means a conduit was dropped from one hourly file."
    )
    if (strict_missing) stop(msg) else warning(msg)
  }

  invisible(missing)
}

# ------------------------------------------------------------
# Version B core: one simulation, one entry hour
# ------------------------------------------------------------

simulate_hourly_loss_one_hour <- function(
    df.flow.hour,
    total.parms,
    df.long.ids,
    date,
    hour,
    sim,
    strict_missing = FALSE
) {
  # ---- Solid phase: threshold is applied here at hourly pipe scale ----
  # calc_loss_solid_pipes() uses the same threshold equations as the manuscript:
  #   ss.ratio.set = max(1 - tau/tau_set_crit, 0)
  #   ss.ratio.res = max(tau/tau_res_crit - 1, 0)
  #   net.rate.pipe.mean = max(k_set - k_res, 0)
  solid_pipe <- calc_loss_solid_pipes(df.flow = df.flow.hour, df.prms = total.parms)

  check_hourly_path_join(
    df.long.ids = df.long.ids,
    pipe_df = solid_pipe,
    strict_missing = strict_missing
  )

  solid_join <- df.long.ids %>%
    left_join(
      solid_pipe %>%
        dplyr::select(
          node_id, part.class, wwtp,
          set.rate.pipe.mean, res.rate.pipe, net.rate.pipe.mean,
          remain.set.deg.pipe.mean
        ) %>%
        distinct(),
      by = c("node_id", "part.class", "wwtp")
    ) %>%
    drop_na(remain.set.deg.pipe.mean)

  # Product across conduits in each path for each particle class.
  solid_path_class <- solid_join %>%
    group_by(node_id_entry, path_nodes_ids, part.class, wwtp) %>%
    summarise(
      remain.set.deg.path.mean = prod(remain.set.deg.pipe.mean),
      mean.k.set.path = mean(set.rate.pipe.mean, na.rm = TRUE),
      mean.k.res.path = mean(res.rate.pipe, na.rm = TRUE),
      mean.k.net.path = mean(net.rate.pipe.mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(node_id = node_id_entry)

  rna_dist <- load_RNA_dist(total.parms)

  solid_agg <- solid_path_class %>%
    left_join(rna_dist, by = "part.class") %>%
    mutate(
      remain.sol.slow.class = F.RNA.skw.slow * remain.set.deg.path.mean,
      remain.sol.fast.class = F.RNA.skw.fast * remain.set.deg.path.mean,
      remain.sol.hmg.class  = F.RNA.hmg      * remain.set.deg.path.mean
    ) %>%
    group_by(node_id, path_nodes_ids, wwtp) %>%
    summarise(
      remain.set.slow = sum(remain.sol.slow.class, na.rm = TRUE),
      remain.set.fast = sum(remain.sol.fast.class, na.rm = TRUE),
      remain.set.hmg  = sum(remain.sol.hmg.class,  na.rm = TRUE),
      mean.k.set.path = mean(mean.k.set.path, na.rm = TRUE),
      mean.k.res.path = mean(mean.k.res.path, na.rm = TRUE),
      mean.k.net.path = mean(mean.k.net.path, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- Liquid phase: hourly biofilm + liquid decay survival ----
  liq_pipe <- loss_liq_pipe(df = df.flow.hour, parms = total.parms)
  liq_path <- loss_liq_all_path(df = liq_pipe, df.long = df.long.ids)

  # ---- Combine solid and liquid survival for this entry hour ----
  total_hour <- liq_path %>%
    dplyr::select(node_id, path_nodes_ids, wwtp, remain.bio.path, remain.deg.liq.path, remain.bio.deg.path) %>%
    left_join(solid_agg, by = c("node_id", "path_nodes_ids", "wwtp")) %>%
    mutate(
      f.solid = get_f_solid_for_wwtp(wwtp, total.parms),
      remain.liq = (1 - f.solid) * remain.bio.deg.path,
      remain.total.slow = f.solid * remain.set.slow + remain.liq,
      remain.total.fast = f.solid * remain.set.fast + remain.liq,
      remain.total.hmg  = f.solid * remain.set.hmg  + remain.liq,
      loss.total.slow = 1 - remain.total.slow,
      loss.total.fast = 1 - remain.total.fast,
      loss.total.hmg  = 1 - remain.total.hmg,
      date = date,
      hour = hour,
      sim = sim
    )

  total_hour
}

# ------------------------------------------------------------
# Version B: one day and one simulation
# ------------------------------------------------------------

simulate_hourly_loss_one_day_sim <- function(
    hourly_file_tbl_day,
    total.parms,
    df.long.ids,
    sim,
    strict_missing = FALSE
) {
  purrr::pmap_dfr(
    list(hourly_file_tbl_day$file, hourly_file_tbl_day$date, hourly_file_tbl_day$hour),
    function(file, date, hour) {
      df.hour <- read_path_state_csv(file) %>% clean()

      simulate_hourly_loss_one_hour(
        df.flow.hour = df.hour,
        total.parms = total.parms,
        df.long.ids = df.long.ids,
        date = date,
        hour = hour,
        sim = sim,
        strict_missing = strict_missing
      )
    }
  )
}

# ------------------------------------------------------------
# 24-h composite aggregation
# ------------------------------------------------------------

aggregate_24h_composite_survival <- function(hourly_results) {
  hourly_results %>%
    group_by(date, sim, node_id, path_nodes_ids, wwtp) %>%
    summarise(
      n_hours = n_distinct(hour),

      # This is the key sensitivity definition:
      # average the surviving fractions after hourly threshold/loss calculations.
      remain.total.slow = mean(remain.total.slow, na.rm = TRUE),
      remain.total.fast = mean(remain.total.fast, na.rm = TRUE),
      remain.total.hmg  = mean(remain.total.hmg,  na.rm = TRUE),

      remain.bio.deg.path = mean(remain.bio.deg.path, na.rm = TRUE),
      remain.set.slow = mean(remain.set.slow, na.rm = TRUE),
      remain.set.fast = mean(remain.set.fast, na.rm = TRUE),
      remain.set.hmg  = mean(remain.set.hmg,  na.rm = TRUE),

      mean.k.set.path = mean(mean.k.set.path, na.rm = TRUE),
      mean.k.res.path = mean(mean.k.res.path, na.rm = TRUE),
      mean.k.net.path = mean(mean.k.net.path, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      loss.total.slow = 1 - remain.total.slow,
      loss.total.fast = 1 - remain.total.fast,
      loss.total.hmg  = 1 - remain.total.hmg
    )
}

# ------------------------------------------------------------
# Comparator table: daily-average Version A vs hourly Version B
# ------------------------------------------------------------

compare_daily_vs_hourly <- function(daily_A, hourly_B) {
  daily_A2 <- daily_A %>%
    transmute(
      sim, node_id, path_nodes_ids, wwtp,
      A_loss_total_slow = loss.total.slow,
      A_loss_total_fast = loss.total.fast,
      A_loss_total_hmg  = loss.total.hmg
    )

  hourly_B2 <- hourly_B %>%
    transmute(
      date, sim, node_id, path_nodes_ids, wwtp,
      B_loss_total_slow = loss.total.slow,
      B_loss_total_fast = loss.total.fast,
      B_loss_total_hmg  = loss.total.hmg,
      n_hours
    )

  left_join(
    hourly_B2,
    daily_A2,
    by = c("sim", "node_id", "path_nodes_ids", "wwtp")
  ) %>%
    mutate(
      diff_B_minus_A_slow = B_loss_total_slow - A_loss_total_slow,
      diff_B_minus_A_fast = B_loss_total_fast - A_loss_total_fast,
      diff_B_minus_A_hmg  = B_loss_total_hmg  - A_loss_total_hmg,
      pct_diff_B_minus_A_slow = 100 * diff_B_minus_A_slow / pmax(A_loss_total_slow, .Machine$double.eps),
      pct_diff_B_minus_A_fast = 100 * diff_B_minus_A_fast / pmax(A_loss_total_fast, .Machine$double.eps),
      pct_diff_B_minus_A_hmg  = 100 * diff_B_minus_A_hmg  / pmax(A_loss_total_hmg,  .Machine$double.eps)
    )
}
