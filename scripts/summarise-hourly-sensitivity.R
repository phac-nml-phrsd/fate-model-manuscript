# ============================================================
# summarise-hourly-sensitivity.R
# ============================================================
# Compare daily-averaged and hourly-resolved slow in-sewer loss
# using only:
#
#   Daily:
#     out/sim_df_loss_total.rds
#
#   Hourly:
#     out/hourly-sensitivity/sim_df_loss_total.rds
#
# This fast version:
#   - does NOT create huge daily_dat/hourly_dat/plot_dat objects
#   - does NOT keep path_nodes_ids
#   - selects only needed columns using transmute()
#   - summarizes to path-level statistics before plotting
#
# Because hourly output does not contain separate:
#   remain.bio.path
#   remain.deg.liq.path
#
# the comparison uses:
#   1. Settling / Resuspension
#   2. Biofilm Adsorption + Biodegradation combined
#   3. Total slow loss
#
# Main figure:
#   True boxplots of path-level median loss,
#   comparing daily-averaged vs hourly-resolved.
#
# Figure layout:
#   - WWTP on right-side row facet labels
#   - fate process on top column facet labels
#   - no node_id on x-axis
#
# This script requires completed daily and hourly model outputs. The public
# one-hour input is an input-format example only; it is insufficient for the
# manuscript's 24-hour composite comparison.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(here)
library(readr)

theme_set(theme_bw())

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

daily_file <- here("out", "sim_df_loss_total.rds")

hourly_file <- here(
  "out",
  "hourly-sensitivity",
  "sim_df_loss_total.rds"
)

fig_dir <- here("out", "hourly-sensitivity", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(daily_file)) {
  stop("Daily file not found: ", daily_file)
}

if (!file.exists(hourly_file)) {
  stop("Hourly file not found: ", hourly_file)
}

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

daily_raw <- readRDS(daily_file)
hourly_raw <- readRDS(hourly_file)

# ------------------------------------------------------------
# Colours
# ------------------------------------------------------------
# WWTP colours are retained for consistency with other scripts,
# although this figure uses model-version fill because WWTP is faceted.

col.wwtp <- c(
  North = "#DEEBF7",
  South = "#9ECAE1",
  West  = "#3182BD"
)

col.model_version <- c(
  "Daily-averaged"  = "#DEEBF7",
  "Hourly-resolved" = "#3182BD"
)

# ------------------------------------------------------------
# Helper: standardize WWTP names
# ------------------------------------------------------------

standardize_wwtp <- function(x) {
  case_when(
    str_to_lower(as.character(x)) == "north" ~ "North",
    str_to_lower(as.character(x)) == "south" ~ "South",
    str_to_lower(as.character(x)) == "west"  ~ "West",
    TRUE ~ as.character(x)
  )
}

# ------------------------------------------------------------
# Fast path-level summary function
# ------------------------------------------------------------
# This function:
#   - keeps only required columns
#   - removes path_nodes_ids and other large columns
#   - computes slow loss values
#   - pivots only a reduced object
#   - summarizes to one row per:
#       WWTP x fate process x node_id x model version
# ------------------------------------------------------------

make_path_level_stats <- function(df, model_version_label) {
  
  df_small <- df %>%
    transmute(
      node_id = as.character(node_id),
      wwtp = standardize_wwtp(wwtp),
      
      `Settling / Resuspension` =
        100 * (1 - remain.set.slow),
      
      `Biofilm Adsorption + Biodegradation` =
        100 * (1 - remain.bio.deg.path),
      
      `Total slow loss` =
        100 * loss.total.slow
    ) %>%
    filter(
      wwtp %in% c("North", "South", "West")
    )
  
  df_small %>%
    pivot_longer(
      cols = c(
        `Settling / Resuspension`,
        `Biofilm Adsorption + Biodegradation`,
        `Total slow loss`
      ),
      names_to = "fate_process",
      values_to = "loss_percent"
    ) %>%
    group_by(
      wwtp,
      fate_process,
      node_id
    ) %>%
    summarise(
      n = n(),
      mean_loss_percent = mean(loss_percent, na.rm = TRUE),
      median_loss_percent = median(loss_percent, na.rm = TRUE),
      q05_loss_percent = quantile(loss_percent, 0.05, na.rm = TRUE),
      q25_loss_percent = quantile(loss_percent, 0.25, na.rm = TRUE),
      q75_loss_percent = quantile(loss_percent, 0.75, na.rm = TRUE),
      q95_loss_percent = quantile(loss_percent, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      model_version = model_version_label
    )
}

# ------------------------------------------------------------
# Create small path-level summary data
# ------------------------------------------------------------

path_stats_daily <- make_path_level_stats(
  daily_raw,
  model_version_label = "Daily-averaged"
)

path_stats_hourly <- make_path_level_stats(
  hourly_raw,
  model_version_label = "Hourly-resolved"
)

path_stats <- bind_rows(
  path_stats_daily,
  path_stats_hourly
) %>%
  mutate(
    wwtp = factor(wwtp, levels = c("North", "South", "West")),
    model_version = factor(
      model_version,
      levels = c("Daily-averaged", "Hourly-resolved")
    ),
    fate_process = factor(
      fate_process,
      levels = c(
        "Settling / Resuspension",
        "Biofilm Adsorption + Biodegradation",
        "Total slow loss"
      ),
      labels = c(
        "Settling / Resuspension",
        "Biofilm Adsorption + Biodegradation",
        "Total Loss"
      )
    ),
    node_id = as.character(node_id)
  )

write_csv(
  path_stats,
  file.path(fig_dir, "path_level_stats_slow_loss_daily_vs_hourly.csv")
)

# ------------------------------------------------------------
# Main Figure:
# TRUE boxplot of path-level median loss
#
# Each point underlying the boxplot is one path-level median.
# The boxplot shows the distribution of path medians within
# each WWTP, fate process, and model version.
#
# Layout:
#   - rows: WWTP
#   - columns: fate process
#   - x-axis: model version
#   - no node_id on the x-axis
# ------------------------------------------------------------

p_loss_boxplot <- ggplot(
  path_stats,
  aes(
    x = model_version,
    y = median_loss_percent / 100,
    fill = model_version
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.alpha = 0.20,
    outlier.size = 0.7,
    linewidth = 0.35
  ) +
  facet_grid(
    wwtp ~ fate_process,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = col.model_version
  ) +
  labs(
    x = NULL,
    y = "Proportion of median in-sewer loss per path",
    fill = NULL
    #title = "Slow-loss in-sewer loss: daily-averaged vs hourly-resolved"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(
      fill = "#2b4c7e",
      color = NA
    ),
    strip.text = element_text(
      color = "white",
      face = "bold",
      size = 13
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text = element_text(
      size = 12
    ),
    axis.text.x = element_text(
      face = "bold",
      angle = 20,
      hjust = 1
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    
    legend.position = "none",
    
    plot.margin = margin(12, 18, 10, 10)
  )

ggsave(
  filename = file.path(fig_dir, "slow_loss_daily_vs_hourly_true_boxplot.png"),
  plot = p_loss_boxplot,
  width = 11,
  height = 7.5,
  dpi = 250
)

ggsave(
  filename = file.path(fig_dir, "slow_loss_daily_vs_hourly_true_boxplot.pdf"),
  plot = p_loss_boxplot,
  width = 11,
  height = 7.5
)

# ------------------------------------------------------------
# Optional compact summary by WWTP and fate process
# ------------------------------------------------------------

summary_by_wwtp_process <- path_stats %>%
  group_by(
    wwtp,
    fate_process,
    model_version
  ) %>%
  summarise(
    n_paths = n_distinct(node_id),
    median_of_path_medians =
      median(median_loss_percent, na.rm = TRUE),
    mean_of_path_medians =
      mean(median_loss_percent, na.rm = TRUE),
    q25_of_path_medians =
      quantile(median_loss_percent, 0.25, na.rm = TRUE),
    q75_of_path_medians =
      quantile(median_loss_percent, 0.75, na.rm = TRUE),
    median_of_path_means =
      median(mean_loss_percent, na.rm = TRUE),
    mean_of_path_means =
      mean(mean_loss_percent, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_by_wwtp_process,
  file.path(fig_dir, "summary_slow_loss_by_wwtp_process_fast.csv")
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

message("Done.")
message("Daily file used: ", daily_file)
message("Hourly file used: ", hourly_file)
message("Figure and tables saved to: ", fig_dir)

message("\nPath-level stats rows:")
print(path_stats %>% count(model_version, wwtp, fate_process))

message("\nCompact summary:")
print(summary_by_wwtp_process)
