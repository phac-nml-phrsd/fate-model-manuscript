# Reduce the manuscript's multi-gigabyte particle simulation object to the
# histogram counts needed by the public pipe/path sensitivity figures.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
})

source_repo <- Sys.getenv(
  "SOURCE_REPO",
  unset = file.path(dirname(normalizePath(".")), "Stochastic-Fate-Model")
)
source_file <- file.path(source_repo, "out", "sim_df_loss_solid.rds")
output_dir <- file.path("out", "figure-data")
output_file <- file.path(output_dir, "particle_loss_distribution.csv")

if (!file.exists(source_file)) stop("Source particle output not found: ", source_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading source particle simulations: ", source_file)
raw <- readRDS(source_file)
required <- c("wwtp", "part.class", "remain.set.pipe.raw", "remain.set.path.raw")
missing <- setdiff(required, names(raw))
if (length(missing)) stop("Particle output is missing: ", paste(missing, collapse = ", "))

dt <- data.table::as.data.table(raw)[, ..required]
rm(raw)
gc()

summarise_loss <- function(column, level) {
  values <- pmin(pmax(1 - dt[[column]], 0), 1)
  bins <- pmin(floor(values / 0.05), 19L)
  data.table(wwtp = dt$wwtp, part.class = dt$part.class, bin = bins)[,
    .(count = .N), by = .(wwtp, part.class, bin)
  ][, `:=`(
    level = level,
    bin_left = bin * 0.05,
    bin_right = (bin + 1) * 0.05,
    bin_mid = (bin + 0.5) * 0.05
  )]
}

summary <- rbindlist(list(
  summarise_loss("remain.set.pipe.raw", "Pipe"),
  summarise_loss("remain.set.path.raw", "Flow path")
), use.names = TRUE)

readr::write_csv(summary, output_file)
message("Wrote ", output_file, " (", nrow(summary), " rows).")
