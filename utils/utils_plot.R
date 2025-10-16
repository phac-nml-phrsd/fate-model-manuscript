# ---- Functions help to plot in main script


suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(ggridges)
})
library(snowfall)
library(sf)



plot_ridge <- function(df, RNA.dist){
  # df = df.loss.solid
  
  g =  df %>%
    filter(type == RNA.dist)%>%
    pivot_longer(cols = c(loss.set.path.mean,
                          loss.set.path.high,
                          loss.set.path.low),
                 names_to = 'name',
                 values_to = 'value')%>%
    ggplot(aes(x = value, y = name, fill=name))+
    #geom_histogram(bins = 10, position = "identity")+
    geom_density_ridges() +
    facet_grid(wwtp ~ part.class, scales = 'free_y')+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.7)))+
    labs(title = 'Viral loss', 
         subtitle = paste('Distribution type:', RNA.dist),
         x = 'Proportion of Lost Viral Genes - Solid Phase')
  # g
  return(g)
}

#-- plot settling in pipe
plot_settling_pipe <- function(df, df_labels,
                               filename){
  
  # Create a named vector where part.class is the key and vel.set is the value
  velocity_labels <- setNames(
    paste0("Class ", df_labels$part.class, "\n(V=", df_labels$vel.set, " m/d)"), 
    df_labels$part.class
  )
  
  # Label mappings for facet strips and legend
  part_labeller <- as_labeller(velocity_labels)
  wwtp_labels <- c("north" = "North", 
                   "south" = "South", 
                   "west"  = "West")
  wwtp_labeller <- as_labeller(wwtp_labels)
  
  a.pipe <- df %>%
    ggplot(aes(x = 1 - remain.set.pipe.raw, fill = wwtp)) +
    geom_histogram(bins = 20, position = "identity") +
    facet_grid(wwtp ~ part.class, scales = 'free_y',
               labeller = labeller(part.class = part_labeller,
                                   wwtp = wwtp_labeller)) +
    scale_fill_discrete(
      name = "WWTP",         # Legend title
      labels = wwtp_labels   # Legend labels
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = rel(0.7))
    ) +
    labs(
      #title = 'Particle Settling/Resuspension In Pipes', 
      #subtitle = 'Impact of settling velocity and critical shear stress on particle classes',
      x = 'Proportion of Solid Particle Lost Per Pipe',
      y = 'Frequency'
    )
  
  plot(a.pipe)
  
  # Save the plot for manuscript 
  # 
  pdf(filename, width = 9, height = 4)
  plot(a.pipe)
  dev.off()
}

plot_settling_path <- function(df, loss.type, df_labels,
                               filename){
  # df = df.loss.solid
  # loss.type = 'settling low'
  
  if(loss.type == 'settling raw')  vname = 'remain.set.path.raw'
  if(loss.type == 'settling mean') vname = 'remain.set.path.mean'
  if(loss.type == 'settling high') vname = 'remain.set.path.high'
  if(loss.type == 'settling low')  vname = 'remain.set.path.low'
  if(loss.type == 'solid deg')     vname = 'remain.deg.sol.path'
  if(loss.type == 'settling and deg mean')   vname = 'remain.set.deg.path.mean'
  
  
  
  # load parameters for settling velocity and critical shear stress 
  df.prms = load_prms()
  #plot_particleclass_distrib()
  message('Parameters loaded.')
  
  # Create a named vector where part.class is the key and vel.set is the value
  velocity_labels <- setNames(
    paste0("Class ", df_labels$part.class, "\n(V=", df_labels$vel.set, " m/d)"), 
    df_labels$part.class
  )
  
  # Label mappings for facet strips and legend
  part_labeller <- as_labeller(velocity_labels)
  wwtp_labels <- c("north" = "North", 
                   "south" = "South", 
                   "west"  = "West")
  wwtp_labeller <- as_labeller(wwtp_labels)
  
  g =  df %>%
    mutate(loss = 1 - .data[[vname]])%>%
    ggplot(aes(x = loss, fill=wwtp))+
    geom_histogram(bins = 20, position = "identity")+
    facet_grid(wwtp ~ part.class, scales = 'free_y',
               labeller = labeller(part.class = part_labeller,
                                   wwtp = wwtp_labeller)) +
    scale_fill_discrete(
      name = "WWTP",         # Legend title
      labels = wwtp_labels   # Legend labels
    ) +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,by=0.25)) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.7)))+
    labs(#title = 'Particle Settling/resuspension Per Path', 
         #subtitle = paste('Impact of settling velocity and critical shear stress on particle classes \n
         #Distribution type:'),#, loss.type),
         x = 'Proportion of Solid Particle Lost Per Flow Path',
         y = 'Frequancy')
  g
  
  # Save the plot for manuscript 
  # 
  pdf(filename, width = 9, height = 4)
  plot(g)
  dev.off()
  return(g)
  
  
}

# Plot aggregated settling loss over all particle classes
plot_aggregated_loss_solid <- function(df, type){
  
  g.hist =  df%>%
    filter(RNA.type == type)%>%
    pivot_longer(cols = c(remain.set.agg.high,
                          remain.set.agg.low,
                          remain.set.agg.mean,
                          remain.deg.sol.agg,
                          remain.set.deg.agg.mean),
                 names_to = 'name',
                 values_to = 'value')%>%
    ggplot(aes(x= 1- value, fill=wwtp))+
    geom_histogram(bins = 20, #color='#e9ecef',
                   # alpha = 0.6, 
                   position = "identity")+
    facet_grid(wwtp~name, scales = 'free_y')+
    theme_bw()+
    theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12, face='bold'),
          axis.title.y = element_text(size=12, face='bold'))+
    labs(title = 'Aggragated Viral loss cross Particle Classes', 
         subtitle = paste('RNA Distribution type:', type),
         x = 'Proportion of Lost Viral Genes - Solid Phase')
  
  return(plot(g.hist))  
  
}

#' Plot the distributions of solids used by the model.
plot_particleclass_distrib <- function() {
  
  prm = load_prms()
  
  rna = prm |> select(starts_with('F.RNA'), part.class) |> 
    pivot_longer(cols = starts_with('F.RNA.skw.slow'))
  
  g.rna = rna |> ggplot(aes(x = part.class, y = value, color = name))+
    geom_line(linewidth = 2) + geom_point(size=3)+
    geom_text(aes(label = value), nudge_y = 0.02, size = 3, alpha = 0.8) + 
    scale_x_continuous(breaks = unique(rna$part.class))+
    theme(panel.grid.minor = element_blank()) +
    scale_color_manual(values = c('hotpink', 'gold3', 'steelblue1')) +
    labs(title = 'RNA distributions', x = 'particle class', y = 'proportion')
  # g.rna
  
  g.tss =  prm |> 
    select(starts_with('F.TSS'), part.class) |> 
    pivot_longer(cols = starts_with('F.TSS')) |> 
    ggplot(aes(x = part.class, y = value, color = name))+
    geom_line(linewidth = 2) + geom_point(size=3)+
    geom_text(aes(label = value), nudge_y = 0.02, size = 3, alpha = 0.8) + 
    scale_x_continuous(breaks = unique(rna$part.class))+
    theme(panel.grid.minor = element_blank()) +
    scale_color_manual(values = c('orange3', 'royalblue', 'seagreen2')) +
    labs(title = 'TSS distributions', x = 'particle class', y = 'proportion')
  # g.tss
  
  pdf('figs/distribution_TSS_RNA.pdf', height=10, width = 12)
  plot(g.tss/g.rna)
  plot(g.tss)
  plot(g.rna)
  dev.off()
}

#--- plotting 
summary_stats_solid_loss <- function(df.loss.solid, df.loss.agg) {
  
  # Summary stats
  dss = df.loss.solid |> 
    group_by(wwtp, part.class) |> 
    summarise(loss.sett.low  = 1 - mean(remain.set.deg.path.low),
              loss.sett.mean = 1 - mean(remain.set.deg.path.mean),
              loss.sett.high = 1 - mean(remain.set.deg.path.high), 
              .groups = 'drop')
  # dss
  
  dssl = dss |> 
    pivot_longer(cols = starts_with('loss'))
  
  dssl$name = factor(x=dssl$name, levels= c('loss.sett.low', 'loss.sett.mean', 'loss.sett.high'))
  
  scfill = scale_fill_manual(values = c('pink1', 'pink3', 'red3'))
  
  g.no.RNA = dssl |> 
    ggplot(aes(x=part.class, y=value, fill = name))+
    geom_col(position = 'dodge', color = 'grey40', width=0.5)+ 
    scale_x_continuous(breaks = unique(dss$part.class))+
    scfill +
    theme(panel.grid = element_blank()) +
    facet_wrap(~wwtp)+
    labs(title = 'Mean losses', y = 'Loss')
  
  # With RNA distribution
  
  dss2 = df.loss.solid |> 
    group_by(wwtp) |> 
    summarise(loss.sett.low  = 1 - mean(remain.set.deg.path.low),
              loss.sett.mean = 1 - mean(remain.set.deg.path.mean),
              loss.sett.high = 1 - mean(remain.set.deg.path.high))
  
  dss2l = dss2 |> pivot_longer(cols = starts_with('loss')) 
  
  dss2l$name = factor(x=dss2l$name, levels= c('loss.sett.low', 'loss.sett.mean', 'loss.sett.high'))
  
  dss2l |> ggplot(aes(x=name, y=value)) + 
    geom_col(aes(fill = name)) + 
    scfill +
    facet_wrap(~wwtp) +
    guides(fill = 'none')+
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 30, hjust=1))+
    labs(title = 'Mean losses', x = 'settling type', y='Loss')
  
  ssa = df.loss.agg |> 
    group_by(wwtp, RNA.type) |>
    summarise(loss.sett.low  = 1 - mean(remain.set.agg.low), 
              loss.sett.mean = 1 - mean(remain.set.agg.mean), 
              loss.sett.high = 1 - mean(remain.set.agg.high), 
              .groups = 'drop'  ) |> 
    pivot_longer(cols = starts_with('loss'))
  
  ssa$name = factor(x=ssa$name, levels= c('loss.sett.low', 'loss.sett.mean', 'loss.sett.high'))
  ssa$RNA.type = factor(x=ssa$RNA.type, levels= c('skw.slow', 'hmg', 'skw.fast'))
  
  g.RNA = ssa |> ggplot(aes(x=name, y=value)) + 
    geom_col(aes(fill = name))+
    scfill+
    facet_grid(wwtp ~ RNA.type)+
    guides(fill = 'none') +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust=1)) +
    labs(title = 'Mean loss', y='Loss', x = 'settling type')
  
  pdf('figs/losses-mean.pdf', width = 20, height = 10)
  plot(g.no.RNA)
  plot(g.RNA)
  dev.off()
  
}

#Plot hydroalic residence time
plot_travel_time <- function(df){
  
  # df = df.flow <-- DEBUG
  stats <- df %>%
    group_by(wwtp) %>%
    summarise(mean_val = mean(path_hrt, na.rm = TRUE),
              median_val = median(path_hrt, na.rm = TRUE))
  
  # Plot histogram with mean and median lines, faceted by 'wwtp'
  g = ggplot(df, aes(x = path_hrt)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(data = stats,aes(xintercept = mean_val), color = "red", linetype = "dashed", size = 1) +
    geom_vline(data = stats,aes(xintercept = median_val), color = "purple", linetype = "dotted", size = 1) +
    facet_wrap(~wwtp, scales = "free_y") +  # Facet by 'wwtp' with independent y scales
    labs(title = "Hydraulic Residence Time Per Flow Path TO WWTP",
         x = "Traveling Time (hour)", y = "Counts of Flow Path") +
    geom_text(data = stats, aes(x = mean_val, y = Inf, label = paste("Mean:", round(mean_val, 2))), 
              color = "red", vjust = 3.5, hjust = -0.4) +
    geom_text(data = stats, aes(x = median_val, y = Inf, label = paste("Median:", round(median_val, 2))), 
              color = "purple", vjust = 5, hjust = -0.4) +
    theme_minimal()
  plot(g)
}

#---- Plot loss of aggregated settling for hmg and slow RNA distribution
# RNA distribution has been applied here
plot_settling_agg <- function(df){
  
  # df = df.total <-- DEBUG
  stats_df <- df %>%
    mutate(loss.set.path.slow = 1 - remain.set.slow,
           loss.set.path.hmg = 1 - remain.set.hmg) %>%
    pivot_longer(cols = c(loss.set.path.slow),
                 names_to = 'name',
                 values_to = 'value') %>%
    group_by(name, wwtp) %>%
    summarise(mean_val = mean(value, na.rm = TRUE),
              median_val = median(value, na.rm = TRUE),
              .groups = 'drop')
  
  # Plot with histogram, mean, and median lines
  g.hist = df %>%
    mutate(loss.set.path.slow = 1 - remain.set.slow,
           loss.set.path.hmg = 1 - remain.set.hmg) %>%
    pivot_longer(cols = c(loss.set.path.slow),
                 names_to = 'name',
                 values_to = 'value') %>%
    ggplot(aes(x = value, fill = wwtp)) +
    geom_histogram(binwidth = 0.03, color = '#e9ecef', alpha = 0.6, position = "identity") +
    
    # Add mean lines (blue, dashed)
    geom_vline(data = stats_df, aes(xintercept = mean_val), 
               color = "blue", linetype = "dashed", size = 1.2) +
    
    # Add median lines (red, solid)
    geom_vline(data = stats_df, aes(xintercept = median_val), 
               color = "red", linetype = "solid", size = 1.2) +
    
    facet_grid(wwtp ~ name, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold')) +
    labs(x = 'Proportion of Lost Viral Genes', 
         y = "Counts of Flow Path",
         caption = "Blue dashed = Mean | Red solid = Median")
  
  plot(g.hist)
  
}

plot_bio_deg<-function(df){
  
  # df = df.loss.liq.path <-- DEBUG
  g.hist =  df%>%
    mutate(loss.bio.path = 1 - remain.bio.path,
           loss.deg.path = 1 - remain.deg.liq.path)%>%
    pivot_longer(cols = c(loss.bio.path,
                          loss.deg.path),
                 names_to = 'name',
                 values_to = 'value')%>%
    ggplot(aes(x=value, fill=wwtp))+
    geom_histogram(binwidth = 0.01, color='#e9ecef', alpha = 0.6, position = "identity")+
    facet_wrap(~name, scales = 'free_y')+
    theme_bw()+
    theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12, face='bold'),
          axis.title.y = element_text(size=12, face='bold'))+
    coord_cartesian(xlim = c(0, 0.3))+
    labs(x = 'Proportion of Lost Viral Genes', y = "Counts of Flow Path") 
  plot(g.hist) 
}

plot_total_loss<-function(df, meam.plot = TRUE){
  # df = sim.df.loss.total.mean
  # meam.plot = TRUE <-- DEBUG
  
  if(meam.plot == TRUE) vname = 'mean.loss.total.slow'
  if(meam.plot == FALSE) vname = 'loss.total.slow'

  g.hist =  df%>%
    pivot_longer(cols = c(.data[[vname]]),
                 names_to = 'name',
                 values_to = 'value')%>%
    ggplot(aes(x=value, fill=wwtp))+
    geom_histogram(color='#e9ecef', alpha = 0.6, position = "identity")+
    #geom_histogram(alpha = 0.8 , binwidth = 0.05) +
    facet_wrap(~name, scales = 'free_y')+
    theme_bw()+
    theme(axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12, face='bold'),
          axis.title.y = element_text(size=12, face='bold'))+
    labs(x = 'Proportion of Lost Viral Genes', y = "Counts of Flow Path") 
  plot(g.hist) 
  
  
  # Summary stats
  stats_df <- df%>%
    pivot_longer(cols = c(.data[[vname]]),
                 names_to = 'name',
                 values_to = 'value')%>%
    group_by(wwtp) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      q1_val = quantile(value, 0.25),
      q3_val = quantile(value, 0.75),
      .groups = 'drop'
    )
}


#----------- Plotting Maps ------------
plot_map_conduit <-function(df,
                    catchment.poly,
                    wwtp.poly,
                    meam.plot=FALSE){
  
  #df = sim.df.loss.total.mean <---DEBUG
  
  if(meam.plot == TRUE)  vname = 'mean.loss.total.slow'
  if(meam.plot == FALSE) vname = 'loss.total.slow'
  
  catchment.sf   = st_as_sf(catchment.poly, wkt = "geometry",crs = 26914)
  wwtp.sf = st_as_sf(wwtp.poly, wkt = 'geometry',crs = 26914,
                     coords = c("long", "lat"))
  
  df_sf = st_as_sf(df, wkt = "geometry",crs = 26914)%>%
    pivot_longer(cols = c(.data[[vname]]),
                 names_to = 'name',
                 values_to = 'value')
  
  g.map = ggplot(df_sf) + 
    theme_bw()+
    geom_sf(data = catchment.sf, aes(fill = location))+
    geom_sf(aes(color = value), linewidth= 1)+
    geom_sf(data = wwtp.sf, fill = 'black', size = 4, shape=22)+
    facet_wrap(~name)+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_brewer()+
    scale_color_gradient(low = 'grey90', high = 'red3')+
    # Rename legend for 'value'
    labs(color = "Viral Loss") 
  plot(g.map)
}
  
  
#--- PLOT Map of average loss per FSA  
plot_map_loss_fsa <-function(df,
                        catchment.poly,
                        wwtp.poly,
                        meam.plot=FALSE){
  
  #Select and prepare your data
  dff <- df %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(
      wwtp,CFSAUID, ave.loss.total.slow
    ) %>%
    distinct()
  
  ggplot(dff) + 
    geom_sf(aes(fill = ave.loss.total.slow), color = "black", linewidth = 0.3) +
    geom_sf_text(aes(label = CFSAUID), size = 3, color = "black") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_fill_gradient(low = 'grey90', high = 'red3') +
    labs(fill = "Viral Loss (fraction)")
}


# --- Plot Map Population Loss per FSA
plot_map_pop_loss_fsa <- function(df){
  
  #Select and prepare your data
  dff <- df %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(
      CFSAUID,wwtp,
      # linear
      pop.signal.loss.total.slow.fsa,
      # weighted
      pop.signal.loss.total.slow.fsa.weighted
    ) %>%
    distinct()
  
  # Unweighted Population Signal Loss
  ggplot(dff) +
    # Linear
    geom_sf(aes(fill = pop.signal.loss.total.slow.fsa), color = NA) +
    # weighted
    #geom_sf(aes(fill = pop.signal.loss.total.slow.fsa.weighted), color = NA) +
    #----color
    #scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Unweighted Loss (%)") +
    scale_fill_gradient(low = "honeydew", high = "forestgreen") +  # continuous fill scale
    #scale_color_gradient(low = 'grey90', high = 'steelblue4') +  # optional color gradient
    labs(title = "Unweighted Number of Population Missing (Total Slow) by FSA", fill = "Population Loss (number)") +
    geom_sf_text(aes(label = CFSAUID), size = 3, color = "black") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
}

# --- Plot Map Population Loss per FSA
plot_map_inf_rate <- function(df, shed.type='low'){
  
  # Prepare data
  dff <- df %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(CFSAUID, wwtp,
           inf.rate.slow.highshed,
           inf.rate.slow.lowshed) %>%
    distinct()
  
  if (shed.type == 'high') {
    
    dff_highshed <- dff %>%
      filter(!is.na(inf.rate.slow.highshed))
    
    # Plot Highshed
    p <- ggplot(dff_highshed) +
      geom_sf(aes(fill = inf.rate.slow.highshed), color = "white", size = 0.1) +
      scale_fill_viridis_c(
        option = "plasma",
        name = "Minimum Infection Rate for Detection\n
        (out of 100,000 population)",
        na.value = "grey90",
        labels = function(x) {
          paste0(
            format(x, nsmall = 3), "%\u2003 (", 
            round(x/100*1e5),")"
          )
        }
      ) +
      labs(
        title = "Infection Rate (Slow) - Highshed Areas"
      ) +
      geom_sf_text(aes(label = CFSAUID), size = 3, color = "black") +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, "lines")
      )
  }
  
  if (shed.type == 'low') {
    
    dff_lowshed <- dff %>%
      filter(!is.na(inf.rate.slow.lowshed))
    
    # Plot Lowshed
    p <- ggplot(dff_lowshed) +
      geom_sf(aes(fill = inf.rate.slow.lowshed), color = "white", size = 0.1) +
      scale_fill_viridis_c(
        option = "plasma",
        name = "Infection Rate (%) (Inverse)",
        na.value = "grey90",
        labels = function(x) {
          paste0(
            x, "% (1 in ", 
            round(x*100*1e5), " persons)"
          )
        }
      ) +
      labs(
        title = "Minimum Infection Rate (Slow) - Lowshed Areas",
        caption = "Source: your data"
      ) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      )
  }
  
  plot(p)
}



# --- Winnipeg polygons based on FSA
get_polygon_fsa <- function(ids){
  
  # read polygons for FSA
  shapefile_path <- "data/shapefiles/postal-code.shp"  # Update with your actual path
  
  # Read the shapefile
  shape_data <- st_read(shapefile_path)
  
  wpg = shape_data%>%
    filter(PRNAME == 'Manitoba')%>%
    filter(CFSAUID %in% ids)  # Select only FSAs in Winnipeg
  
  if(0){
    # CHECK FSA PLOTS
    ggplot(data = wpg) +
      geom_sf(fill = "lightblue", color = "black") +
      geom_sf_text(aes(label = CFSAUID), size = 3, color = "darkblue") +
      labs(title = "Forward Sortation Areas (FSA) in Winnipeg") +
      theme_minimal()
  }
  
  return(wpg)
}


#===== Map total population in City ===== 
plot_pop_fsa <- function(pop,wpg){
  
  wpg = left_join(wpg.fsa, demograph, by = 'CFSAUID')
  
  ggplot(wpg) +
    geom_sf(aes(fill = C1_COUNT_TOTAL), color = "black", linewidth = 0.2) +  # FSA boundaries filled by population
    geom_sf_text(aes(label = CFSAUID), size = 3) +  # Labels for FSA
    scale_fill_gradient(low = "honeydew", high = "forestgreen") +  # continuous fill scale
    scale_color_gradient(low = 'grey90', high = 'steelblue4') +  # optional color gradient
    labs(title = "FSA and Population - City of Winnipeg", fill = "Population") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
}





#---- PLOT POPULATION SIGNAL LOSS PERCENTAGE -----
barplot_pop_signal_loss <-function(df){
  
  # 1. Select and prepare your data
  df_plot <- df %>%
    st_drop_geometry() %>%
    dplyr::select(
      wwtp,
      # liquid
      pop.signal.loss.bio.deg.path.city,
      pop.signal.loss.bio.deg.path.wwtp,
      # solid
      pop.signal.loss.set.slow.city,
      pop.signal.loss.set.slow.wwtp,
      # total
      pop.signal.loss.total.slow.city,
      pop.signal.loss.total.slow.wwtp
    ) %>%
    distinct()
  
  # 2. Pivot to long format
  df_long <- df_plot %>%
    pivot_longer(
      cols = -wwtp,
      names_to = "loss_type",
      values_to = "loss_value"
    ) %>%
    mutate(
      fraction = case_when(
        grepl("bio.deg.path", loss_type) ~ "Liquid",
        grepl("set.slow", loss_type) ~ "Solid",
        grepl("total.slow", loss_type) ~ "Total",
        TRUE ~ "Other"
      ),
      level = case_when(
        grepl("\\.city$", loss_type) ~ "City",
        grepl("\\.wwtp$", loss_type) ~ "WWTP",
        TRUE ~ "Other"
      ),
      location = ifelse(level == "City", "City", wwtp)
    )
  
  # 3. Now plot
  ggplot(df_long %>% filter(level %in% c("WWTP", "City")),
         aes(x = fraction, y = loss_value, fill = level)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +  # narrower bars
    facet_wrap(~ location, scales = "free_x") +
    labs(
      title = "Population Signal Loss by Fraction (WWTPs + City)",
      x = "Fraction Type",
      y = "Population Signal Loss (%)",
      fill = "Location Level"
    ) +
    scale_fill_manual(values = c("WWTP" = "#1f78b4", "City" = "#33a02c")) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),  # bold facet labels
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}




#---- PLOT EFFECTIVE POPULATION COVERAGE -----
barplot_eff_pop_coverage <-function(df){
  
  # 1. Select and prepare your data
  df_plot <- df %>%
    st_drop_geometry() %>%
    dplyr::select(
      wwtp,
      eff.pop.total.slow.city,
      eff.pop.total.slow.wwtp
    ) %>%
    distinct()
  
  # 2. Pivot to long format
  df_long <- df_plot %>%
    pivot_longer(
      cols = -wwtp,
      names_to = "loss_type",
      values_to = "loss_value"
    ) %>%
    mutate(
      fraction = "Total",
      level = case_when(
        grepl("\\.city$", loss_type) ~ "City",
        grepl("\\.wwtp$", loss_type) ~ "WWTP",
        TRUE ~ "Other"
      ),
      location = ifelse(level == "City", "City", wwtp),
      label = paste(location, level, sep = " - ")
    ) %>%
    filter(level == "WWTP")  # ⭐️ remove City
  
  # 3. Plot only WWTP total bars
  ggplot(df_long, aes(x = location, y = loss_value, fill = level)) +
    geom_bar(stat = "identity", width = 0.3) +
    geom_text(
      aes(label = round(loss_value, 1)),
      vjust = -0.5,
      size = 4.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Total Effective Population (WWTP Only)",
      x = "WWTP",
      y = "Effective Population (%)",
      fill = "Level"
    ) +
    scale_fill_manual(values = c("WWTP" = "#33a02c")) +  # 🌿 green color
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}
if(0){
  # 1. Select and prepare your data
  df_plot <- df %>%
    st_drop_geometry() %>%
    dplyr::select(
      wwtp,
      # liquid
      eff.pop.bio.deg.path.city,
      eff.pop.bio.deg.path.wwtp,
      # solid
      eff.pop.set.slow.city,
      eff.pop.set.slow.wwtp,
      # total
      eff.pop.total.slow.city,
      eff.pop.total.slow.wwtp
    ) %>%
    distinct()
  
  # 2. Pivot to long format
  df_long <- df_plot %>%
    pivot_longer(
      cols = -wwtp,
      names_to = "loss_type",
      values_to = "loss_value"
    ) %>%
    mutate(
      fraction = case_when(
        grepl("bio.deg.path", loss_type) ~ "Liquid",
        grepl("set.slow", loss_type) ~ "Solid",
        grepl("total.slow", loss_type) ~ "Total",
        TRUE ~ "Other"
      ),
      level = case_when(
        grepl("\\.city$", loss_type) ~ "City",
        grepl("\\.wwtp$", loss_type) ~ "WWTP",
        TRUE ~ "Other"
      ),
      location = ifelse(level == "City", "City", wwtp)
    )
  
  # 3. Now plot
  ggplot(df_long %>% filter(level %in% c("WWTP", "City")),
         aes(x = fraction, y = loss_value, fill = level)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
    geom_text(
      aes(label = round(loss_value, 1)),
      position = position_dodge(width = 0.6),
      vjust = -0.5,
      size = 4.5
    ) +
    facet_wrap(~ location, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +  # ⭐️ this line
    labs(
      title = "Effective Population by Fraction (WWTPs + City)",
      x = "Sample Fraction",
      y = "Effective Population (%)",
      fill = "Location Level"
    ) +
    scale_fill_manual(values = c("WWTP" = "#1f78b4", "City" = "#33a02c")) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}  
  


#---- Bar PLOT Minimum Detection per sampling site -----
barplot_inf_rate_highshed <-function(df){
  
  # Prepare data
  df_plot <- df %>%
    st_drop_geometry() %>%
    dplyr::select(
      wwtp,
      ave.inf.rate.bio.deg.highshed,    # Liquid
      ave.inf.rate.set.slow.highshed,   # Solid
      ave.inf.rate.slow.highshed        # Total
    ) %>%
    distinct()
  
  df_long <- df_plot %>%
    pivot_longer(
      cols = -wwtp,
      names_to = "loss_type",
      values_to = "loss_value"
    ) %>%
    mutate(
      fraction = case_when(
        grepl("bio\\.deg", loss_type) ~ "Liquid",
        grepl("set\\.slow", loss_type) ~ "Solid",
        grepl("slow\\.highshed$", loss_type) ~ "Total",
        TRUE ~ "Other"
      )
    )
  
  # Shades of blue for each fraction
  fraction_colors <- c("Liquid" = "#9ecae1", "Solid" = "#6baed6", "Total" = "#2171b5")
  
  # Plot
  ggplot(df_long, aes(x = fraction, y = loss_value, fill = fraction)) +
    geom_bar(stat = "identity", width = 0.6, color = "#1f78b4") +
    geom_text(
      aes(label = round(loss_value, 2)),
      vjust = -0.5,
      size = 4.2
    ) +
    facet_wrap(~ wwtp, strip.position = "top") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = fraction_colors) +
    labs(
      title = "Average Infection Rates by Fraction Type per WWTP",
      x = " Sample Fraction Type",
      y = "Average Infection Rate For Detection(%)",
      fill = "Fraction"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    )
}

