#########
#########  Utils: STOCHASTIC PARAMETERS
#########  


suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
})


source('parameters-stochastic.R')



#========= Stochastic Parameters ================================
load_stoch_prms <-function(num.sample, seed){
  
  #------ load distribution of parameters------
  decay      = get_k.deg()
  sorption   = get_k.sorption()
  partition  = get_kd_partition()
  velocity   = get_set_velocity()
  shear.set  = get_ss.crt.set()
  shear.res  = get_ss.crt.res()
  susp.solid = get_TSS()
  
  
  
  # ---- Sample from pre-defined distributions
  
  ## DEBUG
  # set.seed(123)
  # ns = 5e2
  set.seed(seed)
  ns = num.sample   
  
  k.deg   = rgamma(n=ns, shape= decay[1] , rate= decay[2])
  k.bio   = rnorm(n=ns,  mean= sorption[1], sd= sorption[2])
  k.part  = replicate(n=ns, mean(sample(partition, size = 3, replace = TRUE)))
  set.vel = lapply(velocity, function(v) runif(n=ns, min= v[1], max= v[2]))
  ss.set  = lapply(shear.set, function(s) runif(n=ns, min= s[1], max= s[2]))
  ss.res  = lapply(shear.res, function(s) runif(n=ns, min= s[1], max= s[2]))
  tss     = lapply(susp.solid, function(susp) rlnorm(n=ns, meanlog= susp[1], sdlog=susp[2]))
  
  # ----- Generate a grid dataframe of parameters
  # here parameters are constant per solid class
  sim.grid = data.frame(
    sim       = rep(1:ns),
    kappa.deg = k.deg, 
    kappa.bio = k.bio,
    kd        = k.part,
    tss.North = tss$North,
    tss.South = tss$South,
    tss.West  = tss$West
  )%>%
    ## Here we estimate what fraction of viral genes go to solids 
    mutate(f.solid.N = get_solid_fraction(kd, tss.North),
           f.solid.S = get_solid_fraction(kd, tss.South),
           f.solid.W = get_solid_fraction(kd, tss.West)
    )
  
  ## here parameter is different for each partivcle class
  part.class.grid = data.frame(
    sim        = rep(1:ns),
    part.class = rep(c(1:8), each=ns),
    vel.set    = c(set.vel[[1]],
                   set.vel[[2]],
                   set.vel[[3]],
                   set.vel[[4]],
                   set.vel[[5]],
                   set.vel[[6]],
                   set.vel[[7]],
                   set.vel[[8]]),
    ss.crt.set  = c(ss.set[[1]],
                    ss.set[[2]],
                    ss.set[[3]],
                    ss.set[[4]],
                    ss.set[[5]],
                    ss.set[[6]],
                    ss.set[[7]],
                    ss.set[[8]]),
    ss.crt.res  = c(ss.res[[1]],
                    ss.res[[2]],
                    ss.res[[3]],
                    ss.res[[4]],
                    ss.res[[5]],
                    ss.res[[6]],
                    ss.res[[7]],
                    ss.res[[8]])
  )
  
  stoch.params = left_join(sim.grid, part.class.grid, by='sim')
  
  return(stoch.params)
}

#========== CONSTANT PARAMETERS ==========================================
get_const_params <- function(){
  
  # Fraction of TSS corresponds to the discreet settling velocity (v.set)
  df = data.frame(
    # RNA adsorption distribution across particle classes
    F.RNA.skw.slow = c(0.53, 0.15, 0.15, 0.11, 0.02, 0.02, 0.01, 0.01),
    # highly skewed 
    F.RNA.skw.fast = c(0.80,0.06,0.05, 0.05, 0.01, 0.01, 0.01, 0.01),
    F.RNA.hmg      = c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
    # Solids concentration distribution
    #F.TSS.high     = c(0.1,0.1,0.076, 0.175, 0.161, 0.216, 0.074, 0.098),#WES model   
    #F.TSS.low      = c(0.1,0.3, 0.093, 0.198, 0.127, 0.132, 0.041, 0.009),#WES model
    F.TSS.high     = 1,  
    F.TSS.low      = 1,
    # Resuspension
    erod.const     = c(0.200,0.294,0.340,0.400,0.450,0.500,0.550,0.648), #Erodability const. (kg/m/d)
    n.ratio        = 1,      #define power of erosion
    bed.mass       = 0.01,   #kg/m load per pipe length: Rinse paper 2018
    CR             = 0.4    #coverage ratio, between 0-1 (assumption)
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

## CHECK Solid Fraction
if(0){
  grid = expand.grid(tss = tss$North,
                     kd  = k.part)%>%
    mutate(f.sol = kd * tss * 1e-6 / (1 + kd * tss * 1e-6))
  
  g = ggplot(grid, aes(x = f.sol)) +
    geom_histogram(bins = 30, fill = "gray80", color = "black") +
    labs(x = "Viral proportion associated to solids")+
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 14)
    )
  # Save the plot
  ggsave(
    "doc/figs/figure_solid.png",
    plot = g,  
    width = 5,     # in inches (adjust as needed)
    height = 3,    # in inches
    dpi = 300,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}
