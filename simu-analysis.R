
#########
#########  MIAN SCRIPT FOR RESULT ANALYSIS
#########  STOCHASTIC SIMULATIONS
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

#source('parameters.R')
source('utils/utils.R')
source('utils/utils_plot.R')
source('utils/parameters.R')
source('utils/utils_figures.R')
source('utils/utils_param_stochastic.R')



#=========== View and analyze SIMULATION RESULTS
sim_df_loss_solid = readRDS("~/GitHub/fate-model-manuscript/out/sim_df_loss_solid.rds")
sim_df_loss_total = readRDS("~/GitHub/fate-model-manuscript/out/sim_df_loss_total.rds")
df.flow = readRDS("~/GitHub/fate-model-manuscript/out/df.flow.rds")
stoch.params = readRDS("~/GitHub/fate-model-manuscript/out/stoch.params.rds")
  
#========== LOAD WW POLYGON FILES ================================
message('Loading data...')
iw.polygon = read.csv('data/iw.polygon.csv')
iw.wwtp    = read.csv('data/iw.wwtp.csv')
message('Data loaded.\n')

message('Loading demographic and polygons data...')
#https://www12.statcan.gc.ca/census-recensement/2021/dp-pd/prof/details/page.cfm?Lang=E&SearchText=Winnipeg&DGUIDlist=2021A00054611040&GENDERlist=1,2,3&STATISTIClist=1,4&HEADERlist=0
demo = read.csv("data/census_English_CSV_data.csv")
fsa.names = get_fsa_ids() # fsa for City of Winnipeg
demograph = get_demographics(demo=demo, ids = fsa.names)
wpg.fsa   = get_polygon_fsa(ids = fsa.names)
message('Data loaded.\n')




#======================= MEAN SIMULATION =========

#--- 
# Calculate TOTAL MEAN of Simulations
sim.df.loss.total.mean = get_mean_sim(sim_df_loss_total)
# extract geometry from raw data frame
df.polygons = dplyr::select(df.flow, c(node_id,geometry))
# Merge geometry to the simulation data frame
sim.df.loss.total.mean = left_join(sim.df.loss.total.mean,
                                   df.polygons,
                                   by = "node_id")

# Total Histo PLOTTING from MEAN SIMULATION
#plot_total_loss(sim.df.loss.total.mean, meam.plot = TRUE)



# ================= Loss per FSA ==========================
sim.df.loss.fsa = calcu_loss_fsa(loss = sim.df.loss.total.mean,
                                 demo = demograph,
                                 wpg  = wpg.fsa,
                                 total.parms = stoch.params)
message('Loss per fsa is completed.')



#================== Population Loss  ===================
sim.df.pop.loss = calculate_pop_loss(sim.df.loss.fsa, stoch.params)
message('Population Signal Loss Estimate Completed.')


#================== Minimum infection rate  ===================
sim.df.inf.rate = calcu_inf_rate(sim.df.pop.loss, stoch.params)
message('Minimum Infection Rate For Detection completed.')


message('simu-analysis script is completed.')


