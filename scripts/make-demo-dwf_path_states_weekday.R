# Creates a synthetic (demo) of InfoWorks ICM simulation results
# with the same schema as the real file, and with:
# - per-row path_target (NEWPCC_EFFLUENT_OUT / WEWPCC / S-PL70008977)
# - reference_id and reference_ids set so your utils::rename() can
#   safely rename them to node_id and path_nodes_ids without duplicates.
# - geometry built as a CHAIN of LINESTRINGs: each row's end == next row's start.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sf)
  library(tibble)
  library(stringr)
  library(readr)
  library(jsonlite)
})

set.seed(42)

# ---- helpers ----
mk_ids <- function(n, prefix = "S") {
  paste0(prefix, "-", stringr::str_pad(sample(1e5:9e5, n), 8, pad = "0"))
}

# Build a chained set of points, then turn consecutive pairs into LINESTRINGs
mk_chained_lines <- function(n_edges,
                             x0 = 633361.5,    # start close to your example
                             y0 = 5528403,
                             dx_mean = 0, dx_sd = 6,   # small step in x
                             dy_mean = 10, dy_sd = 8)  # gentle step north-ish in y
{
  x <- numeric(n_edges + 1)
  y <- numeric(n_edges + 1)
  x[1] <- x0; y[1] <- y0
  for (i in seq_len(n_edges)) {
    x[i + 1] <- x[i] + rnorm(1, dx_mean, dx_sd)
    y[i + 1] <- y[i] + rnorm(1, dy_mean, dy_sd)
  }
  geoms <- lapply(seq_len(n_edges), function(i) {
    st_linestring(matrix(c(x[i], x[i + 1], y[i], y[i + 1]), ncol = 2), dim = "XY")
  })
  st_sfc(geoms, crs = 26914)  # NAD83 / UTM zone 14N (as before)
}

mk_paths <- function(nodes) {
  edges <- purrr::map2(nodes[-length(nodes)], nodes[-1], ~ c(.x, .y))
  list(
    nodes = nodes,
    edges = edges,
    reference_ids = paste0(nodes, ".1")
  )
}

# ---- parameters ----
n_pipes <- 60L
shapes  <- c("CIRC", "POLYF")
mats    <- c("PVC", "VC", "POLYF")

# ---- build a single main spine (so each row's path is a suffix of it) ----
spine_us_ids  <- mk_ids(n_pipes + 1)               # node IDs along the spine (no ".1")
spine_ref_ids <- paste0(spine_us_ids, ".1")        # node_id-style with ".1"

# Conduits are consecutive segments of the spine
us_ids <- spine_us_ids[1:n_pipes]
ds_ids <- spine_us_ids[2:(n_pipes + 1)]

# Common tail to resemble your example
tail_nodes <- c(
  "NEWPCC_SURGE.7", "BYPASS_CHAMBER.1", "S-PL00000313.1",
  "SPLIT_NEWPCC_OUTFALL2.1", "S-TE00000316.1", "SPLIT_NEWPCC_OUTFALL1.1",
  "S-MH00014483.1", "S-CO00007374.1", "S-CO00007375.1"
)

# For each conduit i (starting at spine_ref_ids[i]), path is suffix of spine + tail
path_nodes_ids <- lapply(seq_len(n_pipes), function(i) {
  c(spine_ref_ids[i:length(spine_ref_ids)], tail_nodes)
})

# ---- geometry: CHAINED lines like the real file ----
geom <- mk_chained_lines(n_pipes)

# Random hydraulic-ish values (bounded, plausible ranges)
conduit_length <- runif(n_pipes, 10, 80)
diam_m         <- runif(n_pipes, 0.25, 1.0)
depth_m        <- pmin(diam_m, pmax(0.05, rbeta(n_pipes, 2, 5) * diam_m))
vel_mps        <- runif(n_pipes, 0.1, 1.2)
flow_cms       <- vel_mps * (pi * (diam_m/2)^2) * runif(n_pipes, 0.1, 0.8)

# Cross-sections (very rough)
conduit_xs     <- pi * (diam_m/2)^2
flow_frac      <- pmin(1, depth_m / diam_m)
flow_xs        <- conduit_xs * flow_frac
flow_area      <- flow_xs
flow_perim     <- pi * diam_m * sqrt(flow_frac) * 0.5

res_time       <- conduit_length / pmax(vel_mps, 1e-3)
hrt            <- res_time
hrt_inv        <- 1 / pmax(hrt, 1e-6)

conduit_width  <- diam_m
conduit_height <- diam_m
conduit_shape  <- sample(shapes, n_pipes, replace = TRUE)
conduit_mat    <- sample(mats, n_pipes, replace = TRUE)

# Additional columns; keep plausible scales
conduit_av     <- runif(n_pipes, 0.05, 2)
conduit_rough  <- round(runif(n_pipes, 1e-3, 0.02), 6)
conduit_ff     <- runif(n_pipes, 0.05, 5)
conduit_ss     <- runif(n_pipes, 0.01, 0.5)

# Keep a small global path block for schema parity (unused by rename to node paths)
all_nodes <- unique(c(us_ids, ds_ids))
path_nodes_small <- all_nodes[seq_len(min(length(all_nodes), 20))]
paths <- mk_paths(path_nodes_small)

# --- per-row path_target mix (adjust weights as you like) ---
demo_terminals <- c("NEWPCC_EFFLUENT_OUT", "WEWPCC", "S-PL70008977")
demo_probs     <- c(0.75, 0.10, 0.15)
path_target_vec <- sample(demo_terminals, n_pipes, replace = TRUE, prob = demo_probs)

# (Optional) if your pipeline expects this param later, include a neutral default:
skw.fast <- rep(0.0, n_pipes)

# Assemble data frame
iw.flow.sf <- tibble(
  X                   = seq_len(n_pipes),
  us_node_id          = us_ids,
  ds_node_id          = ds_ids,
  geometry            = geom,                       
  depth               = depth_m,
  flow                = flow_cms,
  vol                 = flow_cms * res_time,
  conduit_height      = conduit_height,
  conduit_width       = conduit_width,
  conduit_shape       = conduit_shape,
  conduit_length      = conduit_length,
  conduit_material    = conduit_mat,
  conduit_volume      = conduit_xs * conduit_length,
  conduit_xs          = conduit_xs,
  flow_xs             = flow_xs,
  flow_perimeter      = flow_perim,
  flow_area           = flow_area,
  flow_volume         = flow_cms * res_time,
  flow_velocity       = vel_mps,
  residence_time      = res_time,
  conduit_av          = conduit_av,
  conduit_hrt         = hrt,
  conduit_hrt_inv     = hrt_inv,
  conduit_roughness   = conduit_rough,
  conduit_ff          = conduit_ff,
  conduit_ss          = conduit_ss,
  
  # Leave these so your utils::rename() can map them
  reference_id        = paste0(us_ids, ".1"),  # -> node_id
  path_source         = first(path_nodes_small),
  path_target         = path_target_vec,
  
  # List-cols (your utils can rename reference_ids -> path_nodes_ids)
  nodes               = list(paths$nodes),
  edges               = list(paths$edges),
  reference_ids       = path_nodes_ids,        # -> path_nodes_ids
  
  path_hrt            = runif(n_pipes, 0.0005, 0.02),
  path_weights        = runif(n_pipes, 0.0, 0.2),
  path_av             = runif(n_pipes, 0.01, 2),
  path_ss             = runif(n_pipes, 0.005, 0.5),
  path_ss_mean        = runif(n_pipes, 0.01, 0.4),
  path_length         = runif(n_pipes, 50, 500),
  color               = sample(c("red","blue"), n_pipes, replace = TRUE),
  skw.fast            = skw.fast
) |>
  st_as_sf()

# Write demo CSV (serialize list-cols; geometry as WKT)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

iw.flow.out <- iw.flow.sf %>%
  dplyr::mutate(geometry_wkt = sf::st_as_text(geometry)) %>%  # make WKT in a new column
  sf::st_drop_geometry() %>%                                   # drop the sf sfc column
  dplyr::rename(geometry = geometry_wkt) %>%                   # keep the name "geometry"
  dplyr::mutate(across(where(is.list),
                       ~ vapply(., function(x) jsonlite::toJSON(x, auto_unbox = TRUE),
                                character(1))))

readr::write_csv(iw.flow.out, "data/demo_dwf_path_states (weekday).csv")

message("Wrote demo to data/demo_dwf_path_states (weekday).csv")

