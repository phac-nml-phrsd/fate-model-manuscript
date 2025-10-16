#########
#########  STOCHASTIC SIMULATION PROCESSES
#########  LOSS ESTIMATE 
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
})

source('utils/utils.R')



#========== SIMULATION FUNCTION ============================

simulate_calc_loss <- function(df.flow, total.parms, df.long.ids, n.cores, sim){
  
  #--- SOLID PHASE -------------------------------
  #_________________________________
  # Settling/resuspension and Decay
  
  
  #--- Calculate loss for solid phase in every PIPE
  df.loss.solid.pipe = calc_loss_solid_pipes(df.flow = df.flow, df.prms = total.parms)
  message('Solid losses in pipes completed.')
  
  
  #--- CALCULATE Loss in solid phase for PATH
  df.loss.solid = loss_solid_all_paths(df = df.loss.solid.pipe,
                                       df.long = df.long.ids)
  df.loss.solid$sim = sim #simulation number
  message('Solid losses in paths completed.')
  
  
  #--- CALCULATE SETTLING aggregated across all particle classes
  #---------------------
  # Load different viral distribution on solid particles 
  dist.RNA = load_RNA_dist(total.parms)
  
  df.loss.agg = aggregate_loss_solid(df = df.loss.solid,
                                     df.RNA  = dist.RNA)
  message('Aggregated solid losses completed.')
  
  
  
  
  
  
  #---- LIQUID PHASE --------------------------------------
  #__________________________
  # Biofilm and Decay
  
  # Calculate loss in every Pipe - liquid phase
  df.loss.liq.pipe = loss_liq_pipe(df = df.flow, parms = total.parms) 
  
  message('Liquid losses in pipes completed.')
  
  
  
  # Calculate loss in flow paths - liquid phase
  df.loss.liq.path = loss_liq_all_path(df = df.loss.liq.pipe,
                                       df.long = df.long.ids) 
  message('Liquid losses in paths completed.')
  
  
  
  
  
  
  #---- Merge Losses from Solid and Liquid --------------------
  #__________________________
  #-----------------------------
  df.total = calc_total_loss(df.loss.agg,      #solid
                             df.loss.liq.pipe, #liquid
                             df.loss.liq.path, #liquid
                             total.parms)
  df.total$sim = sim #simulation number
  message('Total Losses completed.')
  
  
  # simulation result - one iteration
  return(list(df.total      = df.total,
              df.loss.solid = df.loss.solid))
}
  




