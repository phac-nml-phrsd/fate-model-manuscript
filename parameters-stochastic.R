####
#### Define Input Parameters to calculate Stochastic loss
####


suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
})


source('utils/utils_tss_rain.R')
source('utils/utils_param.R')


# ==== DECAY RATE (1/d) ====================================
get_k.deg <- function(){
  
  # decay for N gene in ww at T=4 to T=20 
  # Ref: Delatolla et al 2024 "sewer condition"
  # https://www.sciencedirect.com/science/article/pii/S1438463924001585?via%3Dihub#bib50
  # Ref: Li et al. 2023
  # https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1144026/full
  data.deg = c(0.084, 0.114, 0.054, 0.077, 0.060, 0.0190, 
               0.960, 2.160, 0.0175, 0.0240,0.034, 0.134,
               0.09, 0.261, 0.367, 0.169,  #Li et al. 2023
               0.09, 0.07) # Zhang et al. 2023 https://www.mdpi.com/2073-4441/15/11/2132
  
  return(fit_gamma_dist(data.deg))
}


# ==== TOTAL SUSPENDED SOLIDS (mg/l) ====================================
# Use TSS measurement data from City of Winnipeg
# Consider Dec-Feb of years 2020-2023 as dry season
# Fit Gamma Distribution for each plant
#140 - 540 mg/l) Ulaval Roest 2022
get_TSS <- function(){
  
  tss.dry = get_tss_winter_months()
  
  # FIT LogNormal to each WWTP type
  # and return fitted estimates
  lognormal_fits = tss.dry %>%
    group_by(wwtp) %>%
    group_split() %>%
    setNames(unique(tss.dry$wwtp)) %>%
    lapply(function(d) fit_lognormal(d, unique(d$wwtp)))
  
  population = data.frame(wwtp = c('North','South','West'),
                          pop = c(470000, 150000, 120000))
  tss.dry = left_join(tss.dry, population, by='wwtp')
  
  
  a = tss.dry%>%
    group_by(wwtp)%>%
    mutate(gr.per.capita = load / pop *1e3, #(gram) 
           med.per.capita = summary(gr.per.capita)[3])
  
  # return estimates for each plant
  return(lognormal_fits)
}


# ==== SETTLING VELOCITY (m/day) ==============
# REF: Master's Thesis Roest 2022 
# recently revised to capture finer particles which most of the RNA associated with 
# https://libstore.ugent.be/fulltxt/RUG01/003/062/325/RUG01-003062325_2022_0001_AC.pdf
get_set_velocity <- function(){
  
  set_vel = list(
    vel1 = c(0.5, 1),
    vel2 = c(2,   9),
    vel3 = c(10,  21),
    vel4 = c(22,  64),
    vel5 = c(65,  135),
    vel6 = c(136, 545),
    vel7 = c(546 ,1500.0),
    vel8 = c(1501,2829.6)
  )

  return(set_vel)
}

# ==== SETTLING CRITICAL SHEAR STRESS (N/m2) ==============
#REF:https://ens.dk/sites/ens.dk/files/Vindenergi/os-tr-006_hydrography_report.pdf
#REF: https://pubmed.ncbi.nlm.nih.gov/21496881/
# Medium sands, Fine sands, Very fine sands, Silts and Clays
ss.crt.ref = c(0.0831, 0.1201, 0.1530, 0.1895) # unit N/m2
vel.set.ref = c(0.000519, 0.002279, 0.00868, 0.02874)# unit m/s

# critical shear stress for settling for 
# each particle class
get_ss.crt.set <- function(){

  shear_stress_set = list(
    ss.set1 = c(0.002, 0.004),
    ss.set2 = c(0.005, 0.009),
    ss.set3 = c(0.01, 0.019),
    ss.set4 = c(0.02, 0.04),
    ss.set5 = c(0.05, 0.09),
    ss.set6 = c(0.10, 0.13),
    ss.set7 = c(0.14, 0.16),
    ss.set8 = c(0.17, 0.20)
  )
 
  return(shear_stress_set)
}

# ==== RESUSPENSION CRITICAL SHEAR STRESS (N/m2) ==============
#Ref: Rinas 2018 https://pubmed.ncbi.nlm.nih.gov/30566103/
#Ref: Banasiak 2005 https://pubmed.ncbi.nlm.nih.gov/16309729/
#REf: Vezzaro 2011 https://pubmed.ncbi.nlm.nih.gov/21496881/
# critical shear stress for erosion
# for each particle class
get_ss.crt.res <- function(){
  
  shear_stress_res = list(
    ss.res1 = c(0.003, 0.006),
    ss.res2 = c(0.006, 0.01),
    ss.res3 = c(0.02, 0.09),
    ss.res4 = c(0.10, 0.17),
    ss.res5 = c(0.18, 0.29),
    ss.res6 = c(0.3, 0.5),
    ss.res7 = c(0.6, 1.3),
    ss.res8 = c(1.4, 2.00)
  )
  
  return(shear_stress_res)
}


# ==== WATER-SOLID PARTITIONING DISTRIBUTION CONST. (mL/g) ====
get_kd_partition <- function(){
  
  # Distribution constant for solid-liquid partitioning kd = Cs/Cq
  # Ref:Boehm et al 2023 for 5 and 11 WWTP locations in US
  # https://pubs.rsc.org/en/content/articlehtml/2025/ew/d4ew00225c
  # TABLE S10 https://pmc.ncbi.nlm.nih.gov/articles/instance/8969789/bin/EW-008-D1EW00826A-s001.pdf

  data.kd = c(3500,5600,4400,3000,
              5128.6,7943.3, 776.3)
              #,21877.6, 19498.5,12000)

  #Produce distribution around mean value of limited data
  #create_sample_dist(data.kd)
  return(data.kd)
}


# ==== SORPTION RATE TO SEWER BIOFILM (m/d) ==================================
# Sorption rates are normalized to area of biofilm under the experiment took place 
# Ref: Zhang 2023 https://www.mdpi.com/2073-4441/15/11/2132
# Ref: Li 2023 https://www.nature.com/articles/s44221-023-00033-4
get_k.sorption <-function(){
  
  kappa.Zhang  = 3.92*24          # (1/day)
  CI.zhang     = c(2.78,5.44)*24  # confidential interval
  kappa.Li     = 0.05*24          # (1/day)
  sd.range.Li   = 0.0071*24       # (1/day) standard deviation
  
  AV.ratio.Zhang = 70.9     # (1/m)
  AV.ratio.Li    = 63       # (1/m)
  
  # normalize to A/V ratio (m/day) 
  k.zhang = av_normalize(kappa.Zhang, AV.ratio.Zhang)
  
  k.Li    = av_normalize(kappa.Li,    AV.ratio.Li)
  sd.Li   = av_normalize(sd.range.Li, AV.ratio.Li)
  
  return(c(k.Li, sd.Li))
}


# ==== VIRAL RNA FRACTION ASSOCIATED TO SOLID =======================
# Estimate Solid Fraction using distribution coefficient and TSS
get_solid_fraction <- function(kd,tss){
  
  frac = kd * tss * 1e-6 / (1 + kd * tss * 1e-6)
  return(frac)
}



