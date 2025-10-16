####
#### Define Input Parameters to calculate loss
####

library(tidyr)
library(dplyr)
library(ggplot2)

ggplot2::theme_set(theme_bw())



#' Define all model parameters
#' 
load_prms <- function(){
  # ----- REFERENCE https://ens.dk/sites/ens.dk/files/Vindenergi/os-tr-006_hydrography_report.pdf
  # sediment-dependent critical shear stress values 
  # Medium sands, Fine sands, Very fine sands, Silts and Clays
  ss.crt.ref = c(0.0831, 0.1201, 0.1530, 0.1895) # unit N/m2
  vel.set.ref = c(0.000519, 0.002279, 0.00868, 0.02874)# unit m/s
  #---------------------------------------------------------
  
  # Fraction of TSS corresponds to the discreet settling velocity (v.set)
  df = data.frame(
    # RNA adsorption distribution across particle classes
    F.RNA.skw.slow = c(0.53, 0.15, 0.15, 0.11, 0.02, 0.02, 0.01, 0.01),
    F.RNA.skw.fast = c(0.01,0.01,0.02, 0.02, 0.11, 0.15, 0.15, 0.53),
    F.RNA.hmg      = c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
    # Solids concentration distribution
    #F.TSS.high     = c(0.1,0.1,0.076, 0.175, 0.161, 0.216, 0.074, 0.098),#WES model   
    #F.TSS.low      = c(0.1,0.3, 0.093, 0.198, 0.127, 0.132, 0.041, 0.009),#WES model
    F.TSS.high     = 1,  
    F.TSS.low      = 1,
    # Resuspension
    ss.crt.res     = c(0.005,0.01,0.05,0.1,0.18,0.3,0.6,1.4), # critical shear stress for resuspention (N/m2)
    erod.const     = c(0.200,0.294,0.340,0.400,0.450,0.500,0.550,0.648), #Erodability const. (kg/m/d)
    n.ratio        = 1,      #define power of erosion
    bed.mass       = 0.01,   #kg/m load per pipe length: Rinse paper 2018
    CR             = 0.4,    #coverage ratio, between 0-1 (assumption)
    # Settling
    ss.crt.set     = c(0.005,0.01,0.02,0.05,0.1,0.14,0.17,0.2), # critical shear stress for settling (N/m2)
    #vel.set        = c(22,65,136,273,545,1500,2700), # settling velocity (m/day) 
    vel.set        = c(1,10,22,65,136,545,1500,2700), # settling velocity (m/day) 
    # Other parameters
    part.class     = c(1:8),
    f.solid.N      = 0.35, # viral fraction partition to solid (North)
    f.solid.S      = 0.35, # viral fraction partition to solid (South)
    f.solid.W      = 0.35, # viral fraction partition to solid (West)
    kappa.deg      = 0.2,   # decay rate per day for viral genes attached to solid
    # Calculate loss in every Pipe - liquid phase
    kappa.bio = get_kappa_sorption_bio()
  )%>%
    rowwise %>%
    mutate(F.TSS.mean = mean(c(F.TSS.high,F.TSS.low)))
  
  # Checks
  stopifnot(sum(df$F.RNA.skw.slow) == 1)
  stopifnot(sum(df$F.RNA.skw.fast) == 1)
  stopifnot(sum(df$F.RNA.hmg) == 1)
  #stopifnot(sum(df$F.TSS.high) == 1)
  #stopifnot(sum(df$F.TSS.low) == 1)
  
  return(df)
}


#-------------------
#BIOFIOLM

# SORPTION RATE TO SEWER BIOFILM (1/d)
get_kappa_sorption_bio <- function(){
  
  # REFERENCE
  # Zhang et al. 2023 https://www.mdpi.com/2073-4441/15/11/2132
  # Li et al 2023 https://www.nature.com/articles/s44221-023-00033-4#:~:text=Biofilms%20could%20further%20accelerate%20virus,wastewater2%2C16%2C18.
  k1.Zhang = 3.92 *24  #(1/d)
  AV.ratio.Zhang = 70.9 # (1/m)
  k1.Li = 0.05 *24#(1/d)
  AV.ratio.Li = 63 # (1/m) 
  
  kappa.sorp.Zhang = k1.Zhang / AV.ratio.Zhang #(m/day)
  kappa.sorp.Li = k1.Li / AV.ratio.Li #(m/day)
  
  return(kappa.sorp.Li)
}





