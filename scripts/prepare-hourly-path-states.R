# ============================================================
# prepare-hourly-path-states.R
# ============================================================
# Purpose:
#   Convert raw hourly InfoWorks conduit_states files into the
#   same dwf_path_states format used as input by the fate model.
#
# Public example versus manuscript run:
#   This repository includes rows from one real InfoWorks export for
#   2001-01-01 00:00 to document the input schema and demonstrate conversion.
#   The revised
#   manuscript analysis used all 24 hourly exports for each analysed day.
#   One hour alone is not a reproduction of the 24-hour composite analysis.
#
# Input structure expected:
#
#   data/hourly-data/
#     20010101/
#       dwf/
#         conduit_states_dwf_20010101_0000.csv
#         conduit_states_dwf_20010101_0100.csv
#         ...
#     20010102/
#       dwf/
#         conduit_states_dwf_20010102_0000.csv
#         ...
#
# Main template file:
#
#   data/dwf_path_states (weekday).csv
#
# Required template columns:
#   reference_id, reference_ids, path_target, path_hrt, and the conduit/path
#   attributes used by the main fate model. The join key is reference_id.
#
# Required hourly columns:
#   Unnamed: 0 (pipe ID), depth, flow, vol, conduit_height, conduit_width,
#   conduit_shape, conduit_length, conduit_material, conduit_volume,
#   conduit_xs, flow_xs, flow_perimeter, flow_area, flow_volume,
#   flow_velocity, residence_time, conduit_av, conduit_hrt,
#   conduit_hrt_inv, conduit_roughness, conduit_ff, and conduit_ss.
#
# Logic:
#   - Use dwf_path_states as the path/topology template.
#   - Join raw hourly conduit_states by:
#
#       dwf_path_states$reference_id = conduit_states$`Unnamed: 0`
#
#   - Keep only conduits that exist in dwf_path_states.
#   - Report missing conduits for each hourly file.
#   - Drop missing conduits.
#   - Preserve the original dwf_path_states column structure.
#   - Replace only hydraulic-state columns with hourly values.
#
# Output:
#
#   data/hourly-path-states/
#     20010101/
#       dwf_path_states_20010101_0000.csv
#       dwf_path_states_20010101_0100.csv
#       ...
#     20010102/
#       ...
#
#   data/hourly-path-states/hourly_join_report.csv
#   data/hourly-path-states/missing_pipe_ids_by_file.csv
#
# ============================================================


suppressMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(tibble)
})


# ============================================================
# Helper: read either a CSV path or an existing dataframe
# ============================================================

.read_csv_or_df <- function(x) {
  if (is.character(x)) {
    data.table::fread(
      x,
      data.table = FALSE,
      check.names = FALSE
    )
  } else {
    as.data.frame(x, check.names = FALSE)
  }
}


# ============================================================
# Function 1:
# Convert one raw hourly conduit_states file into one
# dwf_path_states-style hourly file
# ============================================================

make_dwf_path_states_hourly <- function(
    dwf_path_states,
    conduit_states_hourly,
    output_file = NULL,
    dwf_id_col = "reference_id",
    drop_missing = TRUE,
    verbose = TRUE,
    return_report = FALSE
) {
  
  # ---- Read files or dataframes ----
  dwf <- .read_csv_or_df(dwf_path_states)
  hourly <- .read_csv_or_df(conduit_states_hourly)
  
  # ---- Clean column names ----
  names(dwf) <- trimws(gsub("^\ufeff", "", names(dwf)))
  names(hourly) <- trimws(gsub("^\ufeff", "", names(hourly)))
  
  # ---- Check DWF reference_id exists ----
  if (!dwf_id_col %in% names(dwf)) {
    stop("Could not find `", dwf_id_col, "` in dwf_path_states.")
  }
  
  # ---- Preserve original DWF column order ----
  dwf_original_cols <- names(dwf)
  
  # ---- Infer hourly ID column by overlap with dwf$reference_id ----
  dwf_ids <- as.character(dwf[[dwf_id_col]])
  
  overlap_count <- sapply(names(hourly), function(col) {
    vals <- as.character(hourly[[col]])
    sum(vals %in% dwf_ids, na.rm = TRUE)
  })
  
  hourly_id_col <- names(which.max(overlap_count))
  best_overlap <- max(overlap_count, na.rm = TRUE)
  
  if (best_overlap == 0) {
    message("Hourly file column names:")
    print(names(hourly))
    
    stop(
      "Could not identify the hourly conduit ID column. ",
      "No column in the hourly file overlaps with dwf_path_states$",
      dwf_id_col, "."
    )
  }
  
  if (verbose) {
    message("Hourly ID column inferred as: `", hourly_id_col, "`")
    message("Number of matching conduit IDs: ", best_overlap)
  }
  
  # ---- Standardize join IDs ----
  dwf <- dwf %>%
    mutate(.join_pipe_id = trimws(as.character(.data[[dwf_id_col]])))
  
  hourly <- hourly %>%
    mutate(.join_pipe_id = trimws(as.character(.data[[hourly_id_col]])))
  
  # ---- If hourly has duplicated pipe IDs, keep first ----
  duplicated_hourly <- hourly %>%
    count(.join_pipe_id) %>%
    filter(n > 1)
  
  if (nrow(duplicated_hourly) > 0) {
    warning(
      "Duplicated pipe IDs found in hourly conduit_states. ",
      "Keeping the first occurrence for each pipe."
    )
    
    hourly <- hourly %>%
      group_by(.join_pipe_id) %>%
      slice(1) %>%
      ungroup()
  }
  
  # ---- Identify missing and extra conduits ----
  dwf_ids_unique <- unique(dwf$.join_pipe_id)
  hourly_ids_unique <- unique(hourly$.join_pipe_id)
  
  missing_in_hourly <- setdiff(dwf_ids_unique, hourly_ids_unique)
  extra_in_hourly <- setdiff(hourly_ids_unique, dwf_ids_unique)
  
  if (verbose) {
    message("--------------------------------------------------")
    message("Hourly conduit-state join report")
    message("--------------------------------------------------")
    message("Rows in dwf_path_states:        ", nrow(dwf))
    message("Rows in hourly conduit_states:  ", nrow(hourly))
    message("Unique conduits in dwf:         ", length(dwf_ids_unique))
    message("Unique conduits in hourly:      ", length(hourly_ids_unique))
    message("Missing from hourly:            ", length(missing_in_hourly))
    message("Extra in hourly:                ", length(extra_in_hourly))
    
    if (length(missing_in_hourly) > 0) {
      message(
        "Missing conduit IDs from hourly file: ",
        paste(missing_in_hourly, collapse = ", ")
      )
    }
    
    message("--------------------------------------------------")
  }
  
  # ---- Drop DWF conduits missing from hourly file ----
  if (drop_missing) {
    dwf_join_base <- dwf %>%
      filter(.join_pipe_id %in% hourly_ids_unique)
  } else {
    dwf_join_base <- dwf
  }
  
  # ---- Hydraulic columns to replace with hourly values ----
  hydraulic_cols_candidate <- c(
    "geometry",
    "depth",
    "flow",
    "vol",
    "conduit_height",
    "conduit_width",
    "conduit_shape",
    "conduit_length",
    "conduit_material",
    "conduit_volume",
    "conduit_xs",
    "flow_xs",
    "flow_perimeter",
    "flow_area",
    "flow_volume",
    "flow_velocity",
    "residence_time",
    "conduit_av",
    "conduit_hrt",
    "conduit_hrt_inv",
    "conduit_roughness",
    "conduit_ff",
    "conduit_ss"
  )
  
  hydraulic_cols_update <- intersect(
    hydraulic_cols_candidate,
    intersect(names(dwf), names(hourly))
  )
  
  if (length(hydraulic_cols_update) == 0) {
    stop(
      "No shared hydraulic columns were found between dwf_path_states ",
      "and hourly conduit_states."
    )
  }
  
  if (verbose) {
    message("Hydraulic columns replaced with hourly values:")
    message(paste(hydraulic_cols_update, collapse = ", "))
  }
  
  # ---- Join hourly hydraulics onto DWF path template ----
  hourly_small <- hourly %>%
    select(.join_pipe_id, all_of(hydraulic_cols_update))
  
  joined <- dwf_join_base %>%
    left_join(
      hourly_small,
      by = ".join_pipe_id",
      suffix = c("", ".hourly")
    )
  
  # ---- Replace hydraulic columns ----
  for (v in hydraulic_cols_update) {
    hourly_v <- paste0(v, ".hourly")
    
    if (hourly_v %in% names(joined)) {
      joined[[v]] <- joined[[hourly_v]]
    }
  }
  
  # ---- Return to original DWF format ----
  joined <- joined %>%
    select(all_of(dwf_original_cols))
  
  stopifnot(identical(names(joined), dwf_original_cols))
  
  # ---- Save if requested ----
  if (!is.null(output_file)) {
    if (!dir.exists(dirname(output_file))) {
      dir.create(dirname(output_file), recursive = TRUE)
    }
    
    data.table::fwrite(joined, output_file)
    
    if (verbose) {
      message("Saved joined hourly path-state file to: ", output_file)
    }
  }
  
  report <- list(
    n_dwf_rows = nrow(dwf),
    n_hourly_rows = nrow(hourly),
    n_unique_dwf_conduits = length(dwf_ids_unique),
    n_unique_hourly_conduits = length(hourly_ids_unique),
    n_missing_in_hourly = length(missing_in_hourly),
    missing_in_hourly = missing_in_hourly,
    n_extra_in_hourly = length(extra_in_hourly),
    extra_in_hourly = extra_in_hourly,
    hourly_id_col_used = hourly_id_col
  )
  
  attr(joined, "join_report") <- report
  
  if (return_report) {
    return(list(data = joined, report = report))
  }
  
  return(joined)
}
# ============================================================
# Function 2:
# Loop over all 6 days and 24 hours per day for January 2001
# ============================================================

make_all_dwf_path_states_hourly <- function(
    dwf_path_states = file.path("data", "dwf_path_states (weekday).csv"),
    hourly_root = file.path("data", "hourly-data"),
    output_root = file.path("data", "hourly-path-states"),
    day_dirs = NULL,
    hourly_subfolder = "dwf",
    file_pattern = "^conduit_states_dwf_\\d{8}_\\d{4}(\\.csv)?$",
    output_prefix = "dwf_path_states",
    drop_missing = TRUE,
    verbose = TRUE
) {
  
  # ---- Check input folders ----
  if (!file.exists(dwf_path_states)) {
    stop("Cannot find dwf_path_states file: ", dwf_path_states)
  }
  
  if (!dir.exists(hourly_root)) {
    stop("Cannot find hourly root folder: ", hourly_root)
  }
  
  # ---- Create output folder ----
  if (!dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE)
  }
  
  # ---- Identify day folders ----
  if (is.null(day_dirs)) {
    day_dirs <- list.dirs(
      hourly_root,
      recursive = FALSE,
      full.names = FALSE
    )
    
    day_dirs <- day_dirs[str_detect(day_dirs, "^\\d{8}$")]
  }
  
  if (length(day_dirs) == 0) {
    stop("No day folders found in: ", hourly_root)
  }
  
  day_dirs <- sort(day_dirs)
  
  if (verbose) {
    message("Found day folders: ", paste(day_dirs, collapse = ", "))
  }
  
  # ---- Read dwf template once ----
  dwf_template <- data.table::fread(
    dwf_path_states,
    data.table = FALSE,
    check.names = FALSE
  )
  
  all_reports <- list()
  all_missing_ids <- list()
  all_output_files <- character()
  
  # ---- Loop over day folders ----
  for (day in day_dirs) {
    
    input_day_dir <- file.path(hourly_root, day, hourly_subfolder)
    
    if (!dir.exists(input_day_dir)) {
      warning("Skipping day because subfolder was not found: ", input_day_dir)
      next
    }
    
    output_day_dir <- file.path(output_root, day)
    
    if (!dir.exists(output_day_dir)) {
      dir.create(output_day_dir, recursive = TRUE)
    }
    
    hourly_files <- list.files(
      input_day_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
    
    hourly_files <- sort(hourly_files)
    
    if (length(hourly_files) == 0) {
      warning("No hourly conduit-state files found in: ", input_day_dir)
      next
    }
    
    if (verbose) {
      message("==================================================")
      message("Processing day: ", day)
      message("Hourly files found: ", length(hourly_files))
      message("==================================================")
    }
    
    # ---- Loop over hourly files ----
    for (f in hourly_files) {
      
      fname <- basename(f)
      
      parsed <- str_match(
        fname,
        "conduit_states_dwf_(\\d{8})_(\\d{4})"
      )
      
      date_string <- parsed[, 2]
      hour_string <- parsed[, 3]
      
      if (is.na(date_string) || is.na(hour_string)) {
        warning("Could not parse date/hour from file name. Skipping: ", fname)
        next
      }
      
      out_file <- file.path(
        output_day_dir,
        paste0(output_prefix, "_", date_string, "_", hour_string, ".csv")
      )
      
      if (verbose) {
        message("Converting: ", fname)
      }
      
      # ---- Convert one raw hourly conduit file ----
      res <- make_dwf_path_states_hourly(
        dwf_path_states = dwf_template,
        conduit_states_hourly = f,
        output_file = out_file,
        drop_missing = drop_missing,
        verbose = FALSE,
        return_report = TRUE
      )
      
      # ---- Save summary report for this file ----
      report_i <- tibble(
        date = date_string,
        hour = hour_string,
        input_file = f,
        output_file = out_file,
        n_dwf_rows = res$report$n_dwf_rows,
        n_hourly_rows = res$report$n_hourly_rows,
        n_unique_dwf_conduits = res$report$n_unique_dwf_conduits,
        n_unique_hourly_conduits = res$report$n_unique_hourly_conduits,
        n_missing_in_hourly = res$report$n_missing_in_hourly,
        n_extra_in_hourly = res$report$n_extra_in_hourly,
        missing_in_hourly_first20 = paste(
          head(res$report$missing_in_hourly, 20),
          collapse = ";"
        ),
        extra_in_hourly_first20 = paste(
          head(res$report$extra_in_hourly, 20),
          collapse = ";"
        )
      )
      
      all_reports[[length(all_reports) + 1]] <- report_i
      all_output_files <- c(all_output_files, out_file)
      
      # ---- Save all missing IDs in long format ----
      if (length(res$report$missing_in_hourly) > 0) {
        missing_i <- tibble(
          date = date_string,
          hour = hour_string,
          input_file = f,
          missing_pipe_id = res$report$missing_in_hourly
        )
        
        all_missing_ids[[length(all_missing_ids) + 1]] <- missing_i
      }
    }
  }
  
  # ---- Combine reports ----
  report_df <- bind_rows(all_reports)
  
  missing_ids_df <- bind_rows(all_missing_ids)
  
  # ---- Save reports ----
  report_file <- file.path(output_root, "hourly_join_report.csv")
  missing_file <- file.path(output_root, "missing_pipe_ids_by_file.csv")
  
  data.table::fwrite(report_df, report_file)
  
  if (nrow(missing_ids_df) > 0) {
    data.table::fwrite(missing_ids_df, missing_file)
  } else {
    data.table::fwrite(
      tibble(
        message = "No missing pipe IDs were found in any hourly file."
      ),
      missing_file
    )
  }
  
  if (verbose) {
    message("==================================================")
    message("Finished hourly conversion.")
    message("Converted files: ", length(all_output_files))
    message("Join report saved to: ", report_file)
    message("Missing pipe report saved to: ", missing_file)
    message("Output root: ", output_root)
    message("==================================================")
  }
  
  return(
    list(
      converted_files = all_output_files,
      report = report_df,
      missing_pipe_report = missing_ids_df,
      report_file = report_file,
      missing_file = missing_file,
      output_root = output_root
    )
  )
}


# ============================================================
# Function 3:
# Quick summary of converted files
# ============================================================

summarise_hourly_join_report <- function(report_df) {
  
  report_df %>%
    summarise(
      total_files_processed = n(),
      total_days = n_distinct(date),
      min_files_per_day = min(table(date)),
      max_files_per_day = max(table(date)),
      median_missing_pipes = median(n_missing_in_hourly),
      max_missing_pipes = max(n_missing_in_hourly),
      median_extra_pipes = median(n_extra_in_hourly),
      max_extra_pipes = max(n_extra_in_hourly)
    )
}


# ============================================================
# Run script
# ============================================================
# You can run from the repo root using:
#
#   Rscript scripts/prepare-hourly-path-states.R
#
# If you only want to source the functions and not run immediately,
# source this file; the block below runs only when invoked with Rscript.
# ============================================================

if (sys.nframe() == 0) {
  
  res_hourly <- make_all_dwf_path_states_hourly(
    dwf_path_states = file.path("data", "dwf_path_states (weekday).csv"),
    hourly_root = file.path("data", "hourly-data"),
    output_root = file.path("data", "hourly-path-states"),
    hourly_subfolder = "dwf",
    drop_missing = TRUE,
    verbose = TRUE
  )
  
  message("\nQuick summary:")
  print(summarise_hourly_join_report(res_hourly$report))
  
  message("\nFiles per day:")
  print(
    res_hourly$report %>%
      count(date, name = "n_files")
  )
}

