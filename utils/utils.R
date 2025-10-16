
######################## Functions used in fate script


suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(patchwork)
})
library(snowfall)
library(sf)


#' Clean InfoWorks flow data 
#'
#' @param data Dataframe. 
#' 
clean <- function(data){
  # remove pipes with zero velocity and blank outfalls 
  # (e.g., flow poured into the river not wwtp)
  # renaming wwtps
  data = data %>%
    filter(conduit_hrt > 0 & path_hrt < 30)%>%
    filter(path_target != "")%>%
    mutate(wwtp = case_when(grepl('^N', path_target) ~ 'north',
                            grepl('^S', path_target) ~ 'south',
                            TRUE ~ 'west'))
  
  # Rename some variables for convenience
  data = data |> 
    rename(
      us_node = us_node_id,
      ds_node = ds_node_id,
      node_id = reference_id,
      path_nodes_ids = reference_ids)
  return(data)
}

#' Helper function - Exponential decay
exp_decay <- function(rate, t) {
  res = exp(-rate * t / 24) #unit per day
  return(res)
}

# Helper function - Settling proportion formula
prop_settling <- function(f.tss, velocity, ratio, depth ) {
  set = f.tss * velocity * ratio / (depth) #unit (1/day)
  return(set)
  }

# Helper function - Resuspension rate formula
get_resusp <- function(flow.area, ratio, erod.cons, bed.mass, n, CR, length){
  res = erod.cons * ratio^n * flow.area * CR/ (bed.mass *length)
  return(res)
}

#' Calculate settling rate in all pipes for one particle class 
#' based on daily averaged flow characteristics
#' 
#' @param i  integer. Particle class number (defined as a row of `df.prms`)
#' @param df.prms dataframe. Parameters defining each particle class
#' @param df.flow dataframe of all pipes characteristics.
#' 
calc_set_pipe <- function(i, df.prms, df.flow) {
  # i = 1 #DEBUG
  message('Settlement calculation for all pipes, particle class #',i)
  
  # Extract the parameters associated 
  # with particle class #i 
  p = df.prms[i,] 
  
  # Perform all calculations on all pipes
  tmp = df.flow |> mutate(
    ss.ratio.set      = ifelse(conduit_ss < p$ss.crt.set, 1 - conduit_ss / p$ss.crt.set, 0),
    ss.ratio.res      = ifelse(conduit_ss > p$ss.crt.res, (conduit_ss / p$ss.crt.res) - 1, 0),
    #
    set.rate.pipe.raw   = prop_settling(1, p$vel.set, ss.ratio.set, depth),  
    set.rate.pipe.high  = prop_settling(p$F.TSS.high, p$vel.set, ss.ratio.set, depth),
    set.rate.pipe.low   = prop_settling(p$F.TSS.low,  p$vel.set, ss.ratio.set, depth),
    set.rate.pipe.mean  = prop_settling(p$F.TSS.mean, p$vel.set, ss.ratio.set, depth),
    # 
    res.rate.pipe       = get_resusp(flow_area, ss.ratio.res, p$erod.const, p$bed.mass, p$n.ratio, p$CR, conduit_length),
    #
    net.rate.pipe.raw   = ifelse(set.rate.pipe.raw  > res.rate.pipe , set.rate.pipe.raw  - res.rate.pipe, 0),
    net.rate.pipe.high  = ifelse(set.rate.pipe.high > res.rate.pipe , set.rate.pipe.high - res.rate.pipe, 0),
    net.rate.pipe.low   = ifelse(set.rate.pipe.low  > res.rate.pipe , set.rate.pipe.low  - res.rate.pipe, 0),
    net.rate.pipe.mean  = ifelse(set.rate.pipe.mean > res.rate.pipe , set.rate.pipe.mean - res.rate.pipe, 0),
    #
    remain.set.pipe.raw  = exp_decay(net.rate.pipe.raw, conduit_hrt),  
    remain.set.pipe.high = exp_decay(net.rate.pipe.high, conduit_hrt),
    remain.set.pipe.low  = exp_decay(net.rate.pipe.low, conduit_hrt), 
    remain.set.pipe.mean = exp_decay(net.rate.pipe.mean, conduit_hrt),
    #
    remain.deg.sol.pipe      = exp_decay(p$kappa.deg, conduit_hrt), 
    remain.set.deg.pipe.raw  = exp_decay(net.rate.pipe.raw + p$kappa.deg , conduit_hrt), 
    remain.set.deg.pipe.high = exp_decay(net.rate.pipe.high + p$kappa.deg , conduit_hrt),
    remain.set.deg.pipe.low  = exp_decay(net.rate.pipe.low + p$kappa.deg , conduit_hrt),
    remain.set.deg.pipe.mean = exp_decay(net.rate.pipe.mean + p$kappa.deg , conduit_hrt),
    #
    part.class = p$part.class
  )  
  
  return(tmp)
}


#' Calculate settling loss in every pipe 
#' using `calc_set_pipe()` function
#'
#' @param df.flow dataframe of all pipes characteristics.
#' @param df.prms dataframe. Parameters defining each particle class
#'
calc_loss_solid_pipes <- function(df.flow, df.prms){
  # Loop over all particle classes
  a = lapply(1:nrow(df.prms), calc_set_pipe, 
             df.prms = df.prms, 
             df.flow = df.flow)
  # Merge all results 
  res = bind_rows(a)
  return(res)
}


correct.format <- function(x){
  
  #i=1
  x = sub('\\[', '', x)
  x = sub('\\]', '', x)
  x = as.character(x)
  x = unlist(strsplit(x,","))
  x = sub("'",'',x)
  x = sub("'",'',x)
  x = sub(" ",'',x)
  
  return(x)
}


#' Helper function to pivot longer 
#' 
longformat  <- function(k, node.id.entry, path.id, wwtp,pipes) {
  res = data.frame(
    node_id_entry  = node.id.entry,
    path_nodes_ids = path.id,
    node_id        = pipes,
    wwtp           = wwtp,
    part.class = k
  )
  return(res)
}
#' Helper function that pivot to long format
#' and joins with the losses values for each pipe
#' @param i row of the entry node / path
#'
path_ident <- function(i, df.ids, num.class) {
  
  if(i%%500==0) message('path identification: ',round(i/nrow(df.ids)*100),'%')
  
  # retrieve the entry pipe and its path
  path.id       = df.ids[i,"path_nodes_ids"]
  node.id.entry = df.ids[i,"node_id"]
  wwtp          = df.ids[i,"wwtp"]
  pipes = correct.format(path.id) 
  #message('length pipes = ', length(pipes))
  
  # pivot to a long format for this path only
  # and for all number of  particle classes
  num.part.class = num.class
  tmplong = lapply(1:num.part.class, longformat,
                   node.id.entry = node.id.entry, 
                   path.id = path.id,
                   wwtp = wwtp,
                   pipes = pipes) |> 
    bind_rows() #|>
  
  return(tmplong)
}

get_path_ids_longformat <- function(df,n.cores, num.class){
  
  #df = df.flow <--- DEBUG
  
  # extract the node and path IDs only:
  df.ids = df |> dplyr::select(node_id, path_nodes_ids, wwtp) |> distinct()
  
  
  # identify flow path ids for every pipe and put it in a long format
  sfInit(parallel = n.cores > 1, cpus = n.cores)
  sfExportAll()
  sfLibrary(dplyr)
  tmps = sfLapply(1:nrow(df.ids), path_ident, df.ids=df.ids, num.class=num.class) 
  sfStop()
  df.long = bind_rows(tmps)
  
  return(df.long)
}


#'apply three different RNA distributions for loss calculation
#' @param type is RNA distribution 
apply_RNA_dist_path <- function(type, df){
  
  if(type == 'skw.slow') F.RNA = 'F.RNA.skw.slow'
  if(type == 'skw.fast') F.RNA = 'F.RNA.skw.fast' # highly skewed to small partiles
  if(type == 'hmg')      F.RNA = 'F.RNA.hmg'
  
  res = df |>
    mutate(RNA.remain.set.path.raw  = remain.set.path.raw,
           RNA.remain.set.path.high = .data[[F.RNA]]*remain.set.path.high,
           RNA.remain.set.path.low  = .data[[F.RNA]]*remain.set.path.low,
           RNA.remain.set.path.mean = .data[[F.RNA]]*remain.set.path.mean,
           RNA.remain.deg.sol.path      = .data[[F.RNA]]*remain.deg.sol.path,
           RNA.remain.set.deg.path.raw  = .data[[F.RNA]]*remain.set.deg.path.raw,
           RNA.remain.set.deg.path.high = .data[[F.RNA]]*remain.set.deg.path.high,
           RNA.remain.set.deg.path.low  = .data[[F.RNA]]*remain.set.deg.path.low,
           RNA.remain.set.deg.path.mean = .data[[F.RNA]]*remain.set.deg.path.mean,
           RNA.type = type)
  return(res)
}  


#' Calculate loss in solid phase for *all* flow paths.
#' @param df dataframe as returned by the function `calc_loss_solid_pipes`
#' 
loss_solid_all_paths <- function(df, df.long) {
  
  # df = df.loss.solid.pipe  # <-- DEBUG
  # df.long = df.long.ids

  # extract the node and path IDs only:
  df.ids = df |> dplyr::select(node_id, path_nodes_ids, wwtp) |> distinct()
  
  # Retrieve the losses values (exponential phrases) 
  # for all pipes and all particle classes
  loss.pipe = df |> dplyr::select(
    node_id, 
    part.class, 
    wwtp,
    remain.set.pipe.raw,
    remain.set.pipe.high,
    remain.set.pipe.low,
    remain.set.pipe.mean,
    remain.deg.sol.pipe,
    remain.set.deg.pipe.raw,
    remain.set.deg.pipe.high,
    remain.set.deg.pipe.low,
    remain.set.deg.pipe.mean) |>
    distinct()
  
  # join pipe losses to the path-identified long data frame `df.long`
  dj = left_join(df.long, 
                 loss.pipe, 
                 by = c('part.class','node_id','wwtp'))
  
  # Remove NAs created by unconsistency in "path_ids" versus "node_ids":
  # some "node_ids" were removed at the beginning because of
  # their unrealistic long traveling time (HRT)
  dj = dj |> drop_na()
  
  # LOSS IN EVERY PATH (a path contains multiple pipes)
  # take a product of all pipes in a path
  loss = dj |>
    group_by(path_nodes_ids, part.class, wwtp) |>
    summarise(remain.set.path.raw  = prod(remain.set.pipe.raw),
              remain.set.path.high = prod(remain.set.pipe.high),
              remain.set.path.low  = prod(remain.set.pipe.low),
              remain.set.path.mean = prod(remain.set.pipe.mean),
              remain.deg.sol.path  = prod(remain.deg.sol.pipe),
              remain.set.deg.path.raw  = prod(remain.set.deg.pipe.raw),
              remain.set.deg.path.high = prod(remain.set.deg.pipe.high),
              remain.set.deg.path.low  = prod(remain.set.deg.pipe.low),
              remain.set.deg.path.mean = prod(remain.set.deg.pipe.mean),
              .groups = 'drop')
  
  loss = left_join(loss, 
                   dplyr::select(df.ids, c(node_id,path_nodes_ids)),
                   by = 'path_nodes_ids')
  loss = left_join(loss,
                  dplyr::select(df, c(remain.set.pipe.raw,
                                      remain.set.pipe.high,
                                      remain.set.pipe.low,
                                      remain.set.pipe.mean,
                                      conduit_length,
                                      node_id,
                                      part.class
                                      )),
                  by = c('node_id','part.class'))
  
  return(loss)
}


# Calculate LOSS of settling aggregated across all particle classes 
aggregate_loss_solid <- function(df, df.RNA){
  
  # df = df.loss.solid <-- DEBUG
  # df.RNA = dist.RNA
  
  # Include RNA distributions to loss data.frame
  df = left_join(df, df.RNA, by = 'part.class')
  
  # using three different assumptions on how virus partition to solids
  #' @param F.RAN.skw.slow skewed viral distribution to slow settling particles
  #' @param F.RNA.skw.fast highly skewed viral distribution to very slow settling particles
  #' @param F.RNA.hmg homogeneous distribution
  
  type = list('skw.slow','skw.fast','hmg')
  tmp  = lapply(type, apply_RNA_dist_path,
                df = df)
  RNA.remain = bind_rows(tmp)
  
  loss = RNA.remain %>%
    group_by(node_id,wwtp,RNA.type) %>%
    summarise(remain.set.agg.raw  = sum(RNA.remain.set.path.raw),
              remain.set.agg.high = sum(RNA.remain.set.path.high),
              remain.set.agg.low  = sum(RNA.remain.set.path.low),
              remain.set.agg.mean = sum(RNA.remain.set.path.mean),
              remain.deg.sol.agg  = sum(RNA.remain.deg.sol.path),
              remain.set.deg.agg.raw  = sum(RNA.remain.set.deg.path.raw),
              remain.set.deg.agg.high = sum(RNA.remain.set.deg.path.high),
              remain.set.deg.agg.low  = sum(RNA.remain.set.deg.path.low),
              remain.set.deg.agg.mean = sum(RNA.remain.set.deg.path.mean),
              .groups = 'drop')%>%
    mutate(reference_id = node_id)
  
  return(loss)
}



loss_liq_pipe <- function(df, parms){
  #df = df.flow <---- DEBUG
  #parms = total.parms
  
  
  k.deg = unique(parms$kappa.deg)
  k.bio = unique(parms$kappa.bio)
  
  # Perform calculations on every pipes
  # HRT converted from hours to day
  tmp = df |> 
    group_by(node_id, path_nodes_ids, wwtp)%>%
    summarise(remain.bio.pipe = exp_decay(k.bio * conduit_av, conduit_hrt),
              remain.deg.liq.pipe = exp_decay(k.deg, conduit_hrt),
              remain.bio.deg.pipe = exp_decay(k.deg + k.bio * conduit_av, conduit_hrt),
              .groups = 'drop')
  return(tmp)
}


loss_liq_all_path <- function(df,df.long){
  
  # df = df.loss.liq.pipe <-- DEBUG
  # df.long = df.long.ids

  # remove part.class column as it is irrelevant to liquid phase
  df.long = df.long %>%
    dplyr::select(-part.class) %>%
    distinct()
  
  # join pipe losses to the path-identified long data frame `df.long`
  dj = left_join(df.long, dplyr::select(df, -c(path_nodes_ids)),
                 by = c('node_id','wwtp'))
  
  # Remove NAs created by dis-consistency in path_ids versus node_ids
  # there are node_ids which removed at the beginning of the code for their
  # unrealistic long traveling time (HRT)
  dj = dj |> drop_na()
  
  # LOSS IN EVERY PATH (several number of pipes)
  # take a product of all pipes in a path
  loss = dj |>
    group_by(path_nodes_ids, wwtp) |>
    summarise(remain.bio.path     = prod(remain.bio.pipe),
              remain.deg.liq.path = prod(remain.deg.liq.pipe),
              remain.bio.deg.path = prod(remain.bio.deg.pipe),
              .groups = 'drop') 
  
  
  loss = left_join(loss, dplyr::select(df, c(node_id,path_nodes_ids)),
                   by = 'path_nodes_ids')
  
  return(loss)
}



# ------ Merge losses 
merge_losses <- function(df.solid, df.liq.pipe, df.liq.path){
  
  # merge to total data frame 
  df.loss.agg.mean = df.solid%>%
    dplyr::select(node_id, wwtp, RNA.type, remain.set.agg.mean)%>%
    pivot_wider(names_from = RNA.type, values_from = remain.set.agg.mean)%>%
    rename(remain.set.fast = skw.fast,
           remain.set.slow = skw.slow,
           remain.set.hmg = hmg)
  
  total = left_join(df.liq.pipe,df.liq.path,
                    by = c('node_id', "path_nodes_ids",'wwtp'))
  total = left_join(total, df.loss.agg.mean,
                    by = c('node_id','wwtp'))
  
  return(total)    
}

# calculate total losses from solid and liquid
calc_total_loss <- function(df.solid,
                            df.liq.pipe,
                            df.liq.path,
                            total.parms){
  
  # prms = df.prms <-- DEBUG
  # df.solid = df.loss.agg
  # df.liq.pipe = df.loss.liq.pipe
  # df.liq.path = df.loss.liq.path
  
  # partitioning fraction 
  f.solid = unique(total.parms$f.solid.N)
  
  
  df = merge_losses(df.solid, df.liq.pipe, df.liq.path)%>%
    #------ Estimate total loss in solid and liquid
    mutate(remain.liq      = (1-f.solid)* remain.bio.deg.path,
           remain.sol.hmg  = f.solid* remain.set.hmg,
           remain.sol.fast = f.solid* remain.set.fast,
           remain.sol.slow = f.solid* remain.set.slow,
           loss.total.hmg  = 1 - (remain.sol.hmg  + remain.liq),
           loss.total.fast = 1 - (remain.sol.fast + remain.liq),
           loss.total.slow = 1 - (remain.sol.slow + remain.liq)) %>%
    dplyr::select(-c(remain.liq,
              remain.sol.hmg,
              remain.sol.fast,
              remain.sol.slow))
  return(df)
}


#' Distributions of RNA adsorption used by the model.
#' # using three different assumptions on how virus partition to solids
#' @param F.RAN.skw.slow skewed viral distribution to slow settling particles
#' @param F.RNA.skw.fast highly skewed viral distribution to very slow settling particles
#' @param F.RNA.hmg homogeneous distribution

load_RNA_dist <- function(df){
  
  df = df%>%
    dplyr::select(F.RNA.skw.slow,
                  F.RNA.skw.fast,
                  F.RNA.hmg,
                  part.class)
  
  return(df)
}







# --- OLD CODE - DELETE WHEN SURE
# Calculate settling loss for one single flow path.
#
loss_solid_one_path <- function(idx, df){
  #idxx = idx[1] <- debug
  
  path_nodes_ids = idx
  path_nodes_ids = correct.format(path_nodes_ids)
  node_id = path_nodes_ids[1] # all path starts with entering point, node_id
  df = ungroup(df)
  tmp = df[df$node_id %in% path_nodes_ids,] %>% 
    group_by(part.class)%>%
    summarise(prod.set.path.raw  = prod(remain.set.pipe.raw),
              prod.set.path.high = prod(remain.set.pipe.high),
              prod.set.path.low  = prod(remain.set.pipe.low),
              prod.set.path.mean = prod(remain.set.pipe.mean),
              prod.deg.sol.path      = prod(remain.deg.sol.pipe),
              prod.set.deg.path.raw  = prod(remain.set.deg.pipe.raw),
              prod.set.deg.path.high = prod(remain.set.deg.pipe.high),
              prod.set.deg.path.low  = prod(remain.set.deg.pipe.low),
              prod.set.deg.path.mean = prod(remain.set.deg.pipe.mean),
              .groups = 'drop')%>%
    mutate(node_id = node_id)
  
  
  return(tmp)
}

# --- OLD CODE - DELETE WHEN SURE
# Calculate loss in solid phase for *all* flow paths.
#
calc_loss_solid_path_OLD <- function(df, df.RNA) {
  
  # df = df.loss.solid.pipe  # <-- DEBUG
  
  # extract the node and path IDs only:
  df.ids = df |> dplyr::select(node_id,path_nodes_ids) |> distinct()
  idx  = as.list(df.ids$path_nodes_ids)

system.time({  
  sfInit(parallel = n.cores > 1, cpus = n.cores)
  sfExportAll()
  sfLibrary(dplyr)
  tmp  = sflapply(idx, loss_solid_one_path, df = df)
  sfStop()
  tmp2 = bind_rows(tmp)
})
  
  tmp2 = left_join(tmp2, df.RNA, by = 'part.class')
  
  
  # using two different assumption on how virus partition to solids
  #' @param F.RAN.skw.slow skewed viral distribution to slow settling particles
  #' @param F.RNA.skw.fast skewed viral distribution to fast settling particles
  #' @param F.RNA.hmg homogeneous distribution
  
  type = list('skw.slow','skw.fast','hmg')
  tmp.loss  = lapply(type, apply_RNA_dist_path,
                product = tmp2)
  loss = bind_rows(tmp.loss)
  
  
  wwtp.ids = df |> dplyr::select(node_id,wwtp) |> distinct()
    
  loss = left_join(loss,wwtp.ids,by = 'node_id')
  
  return(loss)
}

get_mean_sim <- function(df){
  df.mean = df %>%
    group_by(node_id)%>%
    summarize(
      # loss from bio and deg in liquid
      mean.remain.bio.path     = mean(remain.bio.path),
      mean.remain.deg.liq.path = mean(remain.deg.liq.path),
      mean.remain.bio.deg.path = mean(remain.bio.deg.path),
      # aggregated loss from settling and resuspension
      mean.remain.set.hmg = mean(remain.set.hmg),
      mean.remain.set.fast = mean(remain.set.fast),
      mean.remain.set.slow = mean(remain.set.slow),
      # total loss from set + bio + deg
      mean.loss.total.hmg = mean(loss.total.hmg),
      mean.loss.total.fast = mean(loss.total.fast),
      mean.loss.total.slow = mean(loss.total.slow)
    ) 
  
  node_meta <- df %>%
    dplyr::select(node_id, wwtp, path_nodes_ids) %>%
    distinct()
  
  df.mean <- left_join(node_meta, df.mean, by = "node_id")
  
  return(df.mean)
}


library(data.table)
get_median_sim_pipe_path <- function(df) {
  setDT(df)
  
  df.median <- df[, .(
    median.remain.set.path.raw  = median(remain.set.path.raw,  na.rm = TRUE),
    median.remain.set.path.high = median(remain.set.path.high, na.rm = TRUE),
    median.remain.set.path.low  = median(remain.set.path.low,  na.rm = TRUE),
    median.remain.set.pipe.raw  = median(remain.set.pipe.raw,  na.rm = TRUE),
    median.remain.set.pipe.high = median(remain.set.pipe.high, na.rm = TRUE),
    median.remain.set.pipe.low  = median(remain.set.pipe.low,  na.rm = TRUE)
  ), by = .(sim, node_id, part.class)]
  
  df.meta <- unique(df[, .(sim, node_id, conduit_length, wwtp, path_nodes_ids)])
  
  df.out <- merge(df.median, df.meta, by = c("sim", "node_id"), all.x = TRUE)
  return(df.out)
}

#------ Estimate infection rate

# Function to calculate minimum detectable positivity rate (P_min)
calculate_p_min <- function(PLoD, q, S, TT) {
  P_min = 100*(PLoD * q) / (S * TT)
  return(P_min)
}




get_infection_rate <- function(N,S,TT){
  
  #-- parameters 
  #S     <- 1e9         # Shedding rate in copies/person/day
                        # ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC8443535/
  #PLoD  <- 18300        # Practical limit of detection in copies/L
  PLoD =  2000         # ref: https://www.sciencedirect.com/science/article/pii/S0043135422000951
  PLoD.vec = c(915, 6175, 452)# 50% probability of detection with one replicate
  
  # daily wastewater generation
  # ref: https://legacy.winnipeg.ca/waterandwaste/dept/wastewaterflow.stm
  q = 270          # L/ capita /day
  #Q = q*N         # Wastewater flow in L/day
  #Minimum detectable positivity rate
  P_min <- calculate_p_min(PLoD, q, S, TT)
  
  return(P_min)
}

get_loss_fsa_popu_data <- function(loss,demo,wpg){
  
  # --- assign fsa to upper nodes of conduits
  df.total.fsa = assign_fsa(loss,wpg)
  
  # join population to total loss and fsa information
  df.total.fsa.popu = left_join(df.total.fsa, demo, by = 'CFSAUID')
  
  return(df.total.fsa.popu)
}

# --- Estimate minimum Infection rate for detection 
calcu_loss_fsa <- function(loss, demo, wpg,
                           total.parms){
  
  #--- join loss, fsa and population data
  ff = get_loss_fsa_popu_data(loss, demo, wpg)
  
  # partitioning fraction 
  f.solid = mean(total.parms$f.solid.W)
  
  d = ff%>%
    group_by(CFSAUID)%>%
    summarise(
      #--- Average per FSA (separate processes)
      fsa.ave.remain.bio.path     = mean(mean.remain.bio.path),
      fsa.ave.remain.deg.liq.path = mean(mean.remain.deg.liq.path),
      fsa.ave.remain.bio.deg.path = mean((1-f.solid)*mean.remain.bio.deg.path),
      #fsa.ave.remain.bio.deg.path = mean(mean.remain.bio.deg.path),
      fsa.ave.remain.set.hmg  = mean(f.solid*mean.remain.set.hmg),
      fsa.ave.remain.set.fast = mean(f.solid*mean.remain.set.fast),
      fsa.ave.remain.set.slow = mean(f.solid*mean.remain.set.slow),
      #---- average per FSA (remained)
      # total
      fsa.ave.remain.total.hmg  = mean(1-mean.loss.total.hmg),
      fsa.ave.remain.total.fast = mean(1-mean.loss.total.fast),
      fsa.ave.remain.total.slow = mean(1-mean.loss.total.slow),
      #--- average per FSA (lost)
      fsa.ave.loss.total.hmg  = 1- fsa.ave.remain.total.hmg,
      fsa.ave.loss.total.fast = 1- fsa.ave.remain.total.fast,
      fsa.ave.loss.total.slow = 1- fsa.ave.remain.total.slow,
      fsa.ave.loss.bio.deg    = 1- fsa.ave.remain.bio.deg.path,
      fsa.ave.loss.set.slow   = 1- fsa.ave.remain.set.slow
      )%>%
    ungroup()
  
  # Extract unique mapping of CFSAUID to wwtp and pp_city
  fsa.meta = ff %>%
    st_drop_geometry() %>%
    count(CFSAUID, wwtp) %>%
    group_by(CFSAUID) %>%
    #this keeps one WWTP per CFSAUID — the most frequent one.
    slice_max(n, with_ties = FALSE) %>%
    left_join(
      ff %>%
        st_drop_geometry() %>%
        dplyr::select(CFSAUID, pop_city, C1_COUNT_TOTAL, fsa.geometry) %>%
        distinct(),
      by = "CFSAUID"
    )
  
  # Join the metadata back
  d = left_join(d, fsa.meta, by = "CFSAUID")
  
  return(d)
}

# --- calculate Effective Surveillance Coverage 
# or Population equivalent Signal Loss  
calculate_pop_loss <- function(df,total.parms) {
  
  # DEBUG
  #df = sim.df.loss.fsa
  #total.parms = df.prms
  
  # partitioning fraction 
  f.solid = mean(total.parms$f.solid.N)
  
  dff <- df %>%
    group_by(CFSAUID) %>%
    slice(1) %>%  # Keep just the first row for each CFSAUID
    ungroup() %>%
    
    # --- calculate Effective Surveillance Coverage or Population Signal Loss  
    group_by(wwtp) %>%
    mutate(pop.wwtp = sum(C1_COUNT_TOTAL)) %>%
    mutate(
      # effective pop coverage
      eff.pop.total.slow.wwtp       = sum(fsa.ave.remain.total.slow      * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.total.hmg.wwtp        = sum(fsa.ave.remain.total.hmg       * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.total.fast.wwtp       = sum(fsa.ave.remain.total.fast      * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.bio.deg.path.wwtp     = sum(fsa.ave.remain.bio.deg.path   * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.set.slow.wwtp         = sum(fsa.ave.remain.set.slow       * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.set.hmg.wwtp          = sum(fsa.ave.remain.set.hmg        * C1_COUNT_TOTAL) / pop.wwtp * 100,
      eff.pop.set.fast.wwtp         = sum(fsa.ave.remain.set.fast       * C1_COUNT_TOTAL) / pop.wwtp * 100,
      # --- Loss
      pop.signal.loss.total.slow.wwtp       = 100 - sum(fsa.ave.remain.total.slow      * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.total.hmg.wwtp        = 100 - sum(fsa.ave.remain.total.hmg       * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.total.fast.wwtp       = 100 - sum(fsa.ave.remain.total.fast      * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.bio.deg.path.wwtp     = 100 - sum(fsa.ave.remain.bio.deg.path    * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.set.slow.wwtp         = 100 - sum(fsa.ave.remain.set.slow        * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.set.hmg.wwtp          = 100 - sum(fsa.ave.remain.set.hmg         * C1_COUNT_TOTAL) / pop.wwtp * 100,
      pop.signal.loss.set.fast.wwtp         = 100 - sum(fsa.ave.remain.set.fast        * C1_COUNT_TOTAL) / pop.wwtp * 100
    ) %>%
    ungroup() %>%
    
    mutate(
      pop.signal.loss.total.slow.fsa        = C1_COUNT_TOTAL * (1 - fsa.ave.remain.total.slow),
      pop.signal.loss.total.hmg.fsa         = C1_COUNT_TOTAL * (1 - fsa.ave.remain.total.hmg),
      pop.signal.loss.total.fast.fsa        = C1_COUNT_TOTAL * (1 - fsa.ave.remain.total.fast),
      pop.signal.loss.bio.deg.path.fsa      = C1_COUNT_TOTAL * (1 - fsa.ave.remain.bio.deg.path),
      pop.signal.loss.set.slow.fsa          = C1_COUNT_TOTAL * (1 - fsa.ave.remain.set.slow),
      pop.signal.loss.set.hmg.fsa           = C1_COUNT_TOTAL * (1 - fsa.ave.remain.set.hmg),
      pop.signal.loss.set.fast.fsa          = C1_COUNT_TOTAL * (1 - fsa.ave.remain.set.fast),
      
      pop.signal.loss.total.slow.fsa.weighted    = C1_COUNT_TOTAL * (1 - fsa.ave.remain.total.slow)^2,
      pop.signal.loss.total.hmg.fsa.weighted     = C1_COUNT_TOTAL * (1 - fsa.ave.remain.total.hmg)^2,
      pop.signal.loss.bio.deg.path.fsa.weighted  = C1_COUNT_TOTAL * (1 - fsa.ave.remain.bio.deg.path)^2,
      pop.signal.loss.set.slow.fsa.weighted      = C1_COUNT_TOTAL * (1 - fsa.ave.remain.set.slow)^2,
      pop.signal.loss.set.hmg.fsa.weighted       = C1_COUNT_TOTAL * (1 - fsa.ave.remain.set.hmg)^2,
      
      pop.signal.loss.total.slow.city       = 100 - sum(fsa.ave.remain.total.slow     * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.total.hmg.city        = 100 - sum(fsa.ave.remain.total.hmg      * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.total.fast.city       = 100 - sum(fsa.ave.remain.total.fast     * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.bio.deg.path.city     = 100 - sum(fsa.ave.remain.bio.deg.path   * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.set.slow.city         = 100 - sum(fsa.ave.remain.set.slow       * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.set.hmg.city          = 100 - sum(fsa.ave.remain.set.hmg        * C1_COUNT_TOTAL) / pop_city * 100,
      pop.signal.loss.set.fast.city         = 100 - sum(fsa.ave.remain.set.fast       * C1_COUNT_TOTAL) / pop_city * 100,
      
      # effective pop coverage
      eff.pop.total.slow.city       = sum(fsa.ave.remain.total.slow     * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.total.hmg.city        = sum(fsa.ave.remain.total.hmg      * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.total.fast.city       = sum(fsa.ave.remain.total.fast     * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.bio.deg.path.city     = sum(fsa.ave.remain.bio.deg.path   * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.set.slow.city         = sum(fsa.ave.remain.set.slow       * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.set.hmg.city          = sum(fsa.ave.remain.set.hmg        * C1_COUNT_TOTAL) / pop_city * 100,
      eff.pop.set.fast.city         = sum(fsa.ave.remain.set.fast       * C1_COUNT_TOTAL) / pop_city * 100
    )
  
  return(dff)
}



# --- Estimate minimum Infection rate for detection 
calcu_inf_rate <- function(sim.df.pop.loss, 
                           total.parms){
  
  
  
  # partitioning fraction 
  f.solid = mean(total.parms$f.solid.W)
  
  
  d = sim.df.pop.loss%>%
    group_by(CFSAUID)%>%
    mutate(
      #--- infection rate HiGH shedding (1e11 gc/person)
      inf.rate.hmg.highshed     = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.total.hmg),
      inf.rate.fast.highshed    = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.total.fast),
      inf.rate.slow.highshed    = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.total.slow),
      inf.rate.no.loss.highshed = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = 1),
      # --- liquid and solid separately 
      inf.rate.bio.deg.highshed     = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.bio.deg.path),
      inf.rate.set.slow.highshed    = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.set.slow),
      inf.rate.set.hmg.highshed     = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.set.hmg),
      inf.rate.set.fast.highshed    = get_infection_rate(N=C1_COUNT_TOTAL, S=1e9, TT = fsa.ave.remain.set.fast),
      #--- infection rate LOW shedding (1e5 gc/person)
      inf.rate.hmg.lowshed      = get_infection_rate(N=C1_COUNT_TOTAL, S=1e5, TT = fsa.ave.remain.total.hmg),
      inf.rate.fast.lowhshed    = get_infection_rate(N=C1_COUNT_TOTAL, S=1e5, TT = fsa.ave.remain.total.fast),
      inf.rate.slow.lowshed     = get_infection_rate(N=C1_COUNT_TOTAL, S=1e5, TT = fsa.ave.remain.total.slow),
      inf.rate.no.loss.lowshed  = get_infection_rate(N=C1_COUNT_TOTAL, S=1e5, TT = 1)
    )%>%
    ungroup()%>%
    # --- Average limit of detection per WWTP site
    group_by(wwtp)%>%
    mutate(
      wwtp.inf.rate.slow.highshed     = mean(inf.rate.slow.highshed /100)*100,
      wwtp.inf.rate.set.slow.highshed = mean(inf.rate.set.slow.highshed /100)*100,
      wwtp.inf.rate.set.hmg.highshed  = mean(inf.rate.set.hmg.highshed /100)*100,
      wwtp.inf.rate.set.fast.highshed = mean(inf.rate.set.fast.highshed /100)*100,
      wwtp.inf.rate.bio.deg.highshed  = mean(inf.rate.bio.deg.highshed /100)*100
    )%>%
    ungroup()%>%
    # --- Average limit of detection for whole city
    mutate(
      city.inf.rate.slow.highshed     = mean(inf.rate.slow.highshed /100)*100,
      city.inf.rate.set.slow.highshed = mean(inf.rate.set.slow.highshed /100)*100,
      city.inf.rate.set.hmg.highshed  = mean(inf.rate.set.hmg.highshed /100)*100,
      city.inf.rate.set.fast.highshed = mean(inf.rate.set.fast.highshed /100)*100,
      city.inf.rate.bio.deg.highshed  = mean(inf.rate.bio.deg.highshed /100)*100
    )
    
  return(d)
}

#--- FSAs in City of Winnipeg
get_fsa_ids <-function(){
  
  a = c("R2C", "R2G", "R2H", "R2J", "R2K", "R2L", "R2M", "R2N", "R2P", "R2R", "R2V", "R2W", "R2X", "R2Y",
        "R3A", "R3B",  "R3C","R3E", "R3G", "R3H", "R3J", "R3K", "R3L", "R3M", "R3N", "R3P", "R3R", "R3S", "R3T", "R3V", "R3W", "R3X",
        "R3Y")
  
  if(0){
    #Original fsa information for City of Winnipeg
    c(
      "R2C", "R2E", "R2G", "R2H", "R2J", "R2K", "R2L", "R2M", "R2N", "R2P", "R2R", "R2V", "R2W", "R2X", "R2Y",
      "R3A", "R3B", "R3C", "R3E", "R3G", "R3H", "R3J", "R3K", "R3L", "R3M", "R3N", "R3P", "R3R", "R3S", "R3T", "R3V", "R3W", "R3X", "R3Y",
      "R6M", "R6W")
  }
  
  return(a)
}


get_demographics <- function(demo,ids){
  
  df.demo = demo%>%
    filter(CHARACTERISTIC_NAME  == "Population, 2021")%>%
    mutate(CFSAUID = GEO_NAME)%>%
    dplyr::select("CENSUS_YEAR","DGUID","CFSAUID","CHARACTERISTIC_NAME",
           "C1_COUNT_TOTAL")%>%
    filter(CFSAUID %in% ids)%>%
    group_by(CENSUS_YEAR,CHARACTERISTIC_NAME)%>%
    mutate(pop_city = sum(C1_COUNT_TOTAL))%>%
    ungroup()
  
  return(df.demo)
} 

#---- Assign an FSA (CFSAUID) to each conduit line
assign_fsa <- function(loss, wpg){
 
  # Convert loss to an sf object
  df <- st_as_sf(loss, wkt = "geometry", crs = 26914)
  wpg <- st_transform(wpg, 26914)
  
  wpg$fsa.geometry <- st_geometry(wpg)
  
  
  # Extract the first point from each LINESTRING
  us_node_geom <- st_geometry(df) %>%
    lapply(function(geom) {
      st_point(st_coordinates(geom)[1, ])
    }) %>%
    st_sfc(crs = st_crs(df))  # repackage into an sfc with CRS
  
  # Create a new sf object for the first points
  loss <- st_sf(df, us.geometry = us_node_geom)
  st_geometry(loss) <- loss$us.geometry  # Set for spatial join
  
  # Spatial join: assign FSA from wpg to each first point
  loss.fsa <- st_join(loss, dplyr::select(wpg, c(CFSAUID,fsa.geometry)), left = TRUE)
  
  return(loss.fsa)
}


