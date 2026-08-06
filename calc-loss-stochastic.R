#########
#########  MIAN SCRIPT FOR STOCHASTIC SIMULATION
#########  LOSS ESTIMATE IN ALL PATHS
#########
#########


suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(snowfall)
  library(sf)
  library(furrr)
  library(future)
  library(tibble)
})

source('utils/utils.R')
source('utils/utils_param_stochastic.R')
source('utils/utils_simulation_stochastic.R')

start <- Sys.time()
#========== SET SEED AND NUMBER OF SIMULATION ============================
num.sim <- as.integer(Sys.getenv("NUM_SIM", unset = "1"))
seed <- as.integer(Sys.getenv("SEED", unset = "123"))
input_file <- Sys.getenv(
  "FLOW_INPUT",
  unset = file.path("data", "demo_dwf_path_states (weekday).csv")
)
out_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path("out", "demo-simulation"))
if (is.na(num.sim) || num.sim < 1) stop("NUM_SIM must be a positive integer.")
if (is.na(seed)) stop("SEED must be an integer.")
if (!file.exists(input_file)) stop("Hydraulic input not found: ", input_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Detect number of available cores, limit to avoid overloading
num_cores <- max(1L, min(num.sim, parallel::detectCores() - 2L))
if (num_cores == 1L) future::plan(future::sequential) else future::plan(future::multisession, workers = num_cores)

#========== LOAD INFOWORKS FLOW DATA FILE ================================
message('Loading data...')
iw.flow <- read.csv(input_file)
message('Data loaded.\n')

# clean
# remove pipe HRT < 0  and flow path HRT > 30 hours
df.flow = clean(data = iw.flow)
saveRDS(df.flow, file = file.path(out_dir, "df.flow.rds"))
message('Data cleaned.')

# extract geometry from raw flow data 
df.polygons = dplyr::select(df.flow, c(node_id,geometry))
saveRDS(df.polygons, file = file.path(out_dir, "df.polygons.rds"))
message('Polygons saved.')

#======== LOAD PARAMETERS ===================================
# stochastic parameters 
stoch.params = load_stoch_prms(num.sim, seed)

# fixed parameters 
parm.const = get_const_params()



#========== IDENTIFY FLOW PATHS ================================
message('\nIdentify flow paths and produce a long-format dataframe')

# get number of solid classes
solid.class = length(unique(stoch.params$part.class))

#(oddly, takes longer with multi-cores...)
df.long.ids = get_path_ids_longformat(df= df.flow,n.cores = 1,
                                        num.class = solid.class)%>%
  dplyr::select(- "node_id")%>%
  rename(node_id = node_id_entry)
saveRDS(df.long.ids, file = file.path(out_dir, "df.long.ids.rds"))
message('Long-format path ids created.')



#======== MONTE CARLO SIMULATION FUNCTION ===================  
system.time({
  # Run all simulations using lapply
  sim_results <- future_map(1:max(stoch.params$sim), function(x) {
    
    parms.s <- stoch.params[stoch.params$sim == x, , drop = FALSE]
    total.parms <- cbind(parms.s, parm.const)
    
    #result <- simulate_calc_loss(df.flow, total.parms, df.long.ids, n.cores = 1, sim = x)
    result <- simulate_calc_loss(readRDS(file.path(out_dir, "df.flow.rds")), total.parms, readRDS(file.path(out_dir, "df.long.ids.rds")), n.cores = 1, sim = parms.s$sim[1])
    
    message(paste0('Simulation completed for sim = ', x))
    
    return(result)  # result is a list with df.total and df.loss.solid
  }, .options = furrr_options(seed = TRUE))  # Ensures reproducibility
  
  # Extract and rbind df.total from all simulations
  # Total loss from all processes: liquid + solid
  sim_df_loss_total <- do.call(rbind, lapply(sim_results, `[[`, "df.total"))
  
  # Extract and rbind df.loss.solid from all simulations
  # Only loss from settling/resuspension in each particle class
  sim_df_loss_solid <- do.call(rbind, lapply(sim_results, `[[`, "df.loss.solid"))
  
})



#=========== Save results as R object 
saveRDS(stoch.params,      file = file.path(out_dir, "stoch.params.rds"))
saveRDS(sim_df_loss_solid, file = file.path(out_dir, "sim_df_loss_solid.rds"))
saveRDS(sim_df_loss_total, file = file.path(out_dir, "sim_df_loss_total.rds"))

writeLines(c(
  "Public demonstration simulation",
  paste0("Input: ", input_file),
  paste0("NUM_SIM: ", num.sim),
  paste0("SEED: ", seed),
  "Scope: synthetic hydraulic network; not the full manuscript run"
), file.path(out_dir, "README.txt"))


end <- Sys.time()

cat("Total runtime:", round(difftime(end, start, units = "secs"), 2), "seconds\n")








