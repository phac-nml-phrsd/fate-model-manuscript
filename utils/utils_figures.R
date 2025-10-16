library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)
library(rlang)

# DEFINE COLOUR FOR WWTP - consistent in all plots
col.wwtp <- c(
  North = "#DEEBF7",  # light blue
  South = "#9ECAE1",  # medium blue
  West  = "#3182BD"   # dark blue
)

helper_sim_df <- function(df) {
  #DEBUG
  #df = sim_df_loss_total
  
  res =  df %>%
    mutate(loss.set.path.slow = 1 - remain.set.slow,
           loss.bio.path = 1 - remain.bio.path,
           loss.deg.path = 1 - remain.deg.liq.path)%>%
    pivot_longer(cols = c(loss.set.path.slow,
                          loss.bio.path,
                          loss.deg.path),
                 names_to = 'name',
                 values_to = 'value') %>%
    dplyr::select(sim, node_id, wwtp, name, value) 
  res
}

figure_loss_histo <- function(
    sim_df_loss_total, 
    filename = 'doc/figs/figure_loss_histo.pdf') {
  
  
  # Prepare and label data
  dat <- helper_sim_df(sim_df_loss_total) %>% 
    mutate(
      nameplot = case_when(
        name == 'loss.set.path.slow' ~ 'Settling / Resuspension',
        name == 'loss.bio.path' ~ 'Biofilm Adsorption',
        name == 'loss.deg.path' ~ 'Biodegradation'
      ),
      wwtp = case_when(
        wwtp == "north" ~ "North",
        wwtp == "south" ~ "South",
        wwtp == "west" ~ "West",
        TRUE ~ wwtp
      )
    )
  
  # Summary stats
  stats_df <- dat %>%
    group_by(wwtp, nameplot) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create plot
  g <- dat %>%
    ggplot(aes(x = value)) +
    geom_histogram(
      binwidth = 0.05,
      fill = "#7da0b1",
      color = "white",
      alpha = 0.9
    ) +
    geom_vline(data = stats_df, aes(xintercept = mean_val),
               color = "#1f78b4", linetype = "dashed", size = 1) +
    geom_vline(data = stats_df, aes(xintercept = median_val),
               color = "#ff7f00", linetype = "solid", size = 1) +
    facet_grid(wwtp ~ nameplot, scales = 'free_y') +
    theme_bw(base_size = 14) +  # White background with panel borders
    theme(
      panel.grid.major.x = element_line(color = "gray90"),  # Vertical gridlines
      panel.grid.major.y = element_blank(),  # Remove horizontal
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#2b4c7e", color = NA),
      strip.text = element_text(color = "white", face = "bold", size = 13),
      axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      plot.margin = margin(12, 18, 10, 10)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    labs(
      x = 'Proportion of in-sewer viral biomarker loss',
      y = 'Count of flow paths (500 simulations)'
    )
  
  # Display plot
  g
  # give stats
  stats_df
  
  # Save the plot
  ggsave(
    filename,
    plot = g,  
    width = 10,     # in inches (adjust as needed)
    height = 7,    # in inches
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}



figure_total_loss <- function(df, meam.plot=TRUE) {
  
  # DEBUG
  # df = sim.df.loss.total.mean
  if(meam.plot == TRUE) vname = 'mean.loss.total.slow'
  if(meam.plot == FALSE) vname = 'loss.total.slow'
  
  g =  df %>%
    mutate(wwtp = case_when(
      wwtp == "north" ~ "North",
      wwtp == "south" ~ "South",
      wwtp == "west" ~ "West"))%>%
    pivot_longer(cols = c(.data[[vname]]),
                 names_to = 'name',
                 values_to = 'value')%>%
    ggplot(aes(x=value, fill=wwtp))+
    geom_histogram(color='black',
                   binwidth = 0.05,
                   alpha = 1, position = "identity")+
    scale_fill_manual(values = col.wwtp)+
    facet_wrap(~wwtp, scales = 'free_y', ncol = 1)+
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = rel(0.8)),
          strip.background = element_rect(fill = 'steelblue4'),
          strip.text = element_text(face = 'bold', color = 'white'))+
    guides(fill = 'none')+
    labs(x = 'Proportion of viral biomarkers loss', 
         y = "Counts of flow paths (mean of 500 simulations)") 
  g 
}


figure_map_conduit_loss <-function(df,
                                   catchment.poly,
                                   wwtp.poly,
                                   meam.plot=TRUE){
  
  if(0){
    df             = sim.df.loss.total.mean # or df.loss.inf
    catchment.poly = iw.polygon
    wwtp.poly      = iw.wwtp
  }
  
  if(meam.plot == TRUE)  vname = 'mean.loss.total.slow'
  if(meam.plot == FALSE) vname = 'loss.total.slow'
  
  catchment.sf   = st_as_sf(catchment.poly, wkt = "geometry",crs = 26914)|>
    mutate(wwtp = case_when(
      location == 'N' ~ 'North',
      location == 'SE' ~ 'South',
      location == 'W' ~ 'West'
    ))
  
  wwtp.sf = st_as_sf(wwtp.poly, wkt = 'geometry', crs = 26914) %>%
    mutate(wwtp = case_when(
      location == 'N' ~ 'North',
      location == 'SE' ~ 'South',
      location == 'W' ~ 'West',
      TRUE ~ as.character(location)  # fallback if needed
    ))
  
  
  df_sf = st_as_sf(df, wkt = "geometry",crs = 26914)%>%
    pivot_longer(cols = c(.data[[vname]]),
                 names_to = 'name',
                 values_to = 'value') |>
    rename(loss = value)
  
  
  g.map = ggplot(df_sf) + 
    theme_bw() +
    geom_sf(data = catchment.sf, aes(fill = wwtp)) +
    scale_fill_manual(values = col.wwtp) +
    geom_sf(aes(color = loss), linewidth = 0.8) +
    geom_sf_text(data = wwtp.sf, aes(label = wwtp), 
                 size = 6, nudge_y = 1800) +  # Adjust nudge_y as needed
    geom_sf(data = wwtp.sf, fill = 'black', size = 4, shape = 22) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_color_gradient(low = 'grey96', high = 'red4',
                         labels = scales::percent_format()) +
    labs(color = "In-sewer Viral Loss", 
         fill = "WWTP Catchments")
  
  g.map
}



figure_combo_lossmap <- function(df, 
                                 iw.polygon, 
                                 iw.wwtp, 
                                 filename) {
  #DEBUG
  #df = sim.df.loss.total.mean 
  
  g.map = figure_map_conduit_loss(df, 
                                  catchment.poly = iw.polygon, 
                                  wwtp.poly = iw.wwtp, 
                                  meam.plot = TRUE)
  
  g.histo = figure_total_loss(df, meam.plot = TRUE) 
  
  
  g.combo = g.map + g.histo + 
    plot_layout(ncol = 2, widths = c(2, 1)) + 
    # add panel letters
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = 'collect') & theme(legend.position = 'left')
  
  plot(g.combo)
  
  #pdf(filename, width = 10, height = 5)
  #plot(g.combo)
  #dev.off()
  ggsave(
    filename,
    plot = g.combo,  
    width = 10,     # in inches (adjust as needed)
    height = 5,    # in inches
    dpi = 400,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}



figure_map_loss_fsa <-function(df.fsa, size.text, col.fsa){
  # prepare  data
  dff <- df.fsa %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(
      wwtp,CFSAUID, fsa.ave.loss.total.slow
    ) %>%
    distinct()
  
  g = ggplot(dff) + 
    geom_sf(aes(fill = fsa.ave.loss.total.slow), 
            color = col.fsa, linewidth = 0.1) +
    geom_sf_text(aes(label = CFSAUID), 
                 size = size.text, 
                 color = "black") +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_fill_gradient(low = 'grey95', high = 'red3', 
                        label = scales::percent_format()) +
    labs(fill = "Average in-sewer viral loss \nper neighborhood")
  
  g <- g + theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 8),
    legend.key.width = unit(1, "cm"),
    legend.title = element_text(size = 9))
  g
}

#- map of fsa loss by sample fraction 
figure_map_loss_fsa_fraction <- function(df.fsa, 
                                         size.text,
                                         col.fsa,
                                         filename) {
  # Prepare data: convert to long format
  dff <- df.fsa %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(
      wwtp, CFSAUID, fsa.geometry,
      Total = fsa.ave.loss.total.slow,
      Solid = fsa.ave.loss.set.slow,
      Liquid = fsa.ave.loss.bio.deg
    ) %>%
    distinct() %>%
    pivot_longer(
      cols = c(Total, Solid, Liquid),
      names_to = "LossType",
      values_to = "LossValue"
    ) %>%
    st_as_sf()
  
  # Create faceted map
  g <- ggplot(dff) +
    geom_sf(aes(fill = LossValue), color = col.fsa, linewidth = 0.1) +
    geom_sf_text(aes(label = CFSAUID), size = size.text, color = "black") +
    facet_wrap(~ LossType, ncol = 3) +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.width = unit(2, "cm"),
      legend.title = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold", hjust = 0)  # align left
    ) +
    scale_fill_gradientn(
      # colours = c("grey90", "#FEE0D2", "tomato", "firebrick3", "#67000D"),
      # values  = scales::rescale(c(0, 0.20, 0.50, 0.60, 1)),
      colours = c("grey90", "#FEE0D2","skyblue", "firebrick3", "#67000D"),
      values  = scales::rescale(c(0, 0.20, 0.60, 0.60, 1)),
      limits  = c(0, 1),
      breaks  = c(0, .2, .4, .5, .8, 1),
      labels  = scales::percent
    ) +
    labs(fill = "Average in-sewer viral loss")
  g
  ggsave(
    filename,
    plot = g,  
    width = 9,    
    height = 5,    
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
  return(g)
}




figure_map_pop_loss_fsa <- function(df.fsa, size.text, col.fsa){
  
  #Select and prepare your data
  dff <- df.fsa %>%
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
  g = ggplot(dff) +
    geom_sf(aes(fill = pop.signal.loss.total.slow.fsa), 
            color = col.fsa) +
    scale_fill_gradient(low = "honeydew", high = "forestgreen", labels = scales::scientific) + 
    labs(#title = "Unweighted Number of Population Missing (Total Slow) by FSA", 
      fill = "Population\nequivalent loss") +
    geom_sf_text(aes(label = CFSAUID), 
                 size = size.text, 
                 color = "black") +
    theme_minimal() +
    theme(legend.position = 'bottom',
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()
    )
  
  g <- g + theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 8),
    legend.key.width = unit(1, "cm"),
    legend.title = element_text(size = 9))
  g
}


figure_combo_map_loss_fsa <- function(df.fsa, filename){
  
  size.text = 2
  col.fsa = 'white'
  
  g.pop  = figure_map_pop_loss_fsa(df.fsa, size.text, col.fsa)
  g.pipe = figure_map_loss_fsa(df.fsa, size.text, col.fsa)
  g.bar  = plot_eff_pop_coverage(df.fsa)
  
  # Combine top two maps horizontally
  top_row <- g.pipe + g.pop + plot_layout(ncol = 2)
  
  # Combine with bar plot at the bottom
  g.combo <- top_row / g.bar + 
    plot_layout(heights = c(2, 1)) + 
    plot_annotation(tag_levels = 'A')
  g.combo

  #saving png format
  ggsave(
    filename,
    plot = g.combo,  
    width = 9,    
    height = 10,    
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}



figure_Ceff <- function(df.fsa) {
  
  df_plot <- df.fsa %>%
    st_drop_geometry() %>%
    dplyr::select(
      wwtp,
      eff.pop.total.slow.city,
      eff.pop.total.slow.wwtp
    ) %>%
    distinct()
  
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
    filter(level == "WWTP") 
  
  
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
    scale_fill_manual(values = c("WWTP" = "#33a02c")) +  
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) 
}

figure_map_inf_rate <- function(df.fsa, shed.type='high'){
  
  # Prepare data
  dff <- df.fsa %>%
    st_set_geometry("fsa.geometry") %>%
    dplyr::select(CFSAUID, wwtp,
                  inf.rate.slow.highshed,
                  inf.rate.slow.lowshed) %>%
    distinct()
  
  
  if(shed.type == 'high') dff$rate.per.100k = dff$inf.rate.slow.highshed * 1e3
  if(shed.type == 'low')  dff$rate.per.100k = dff$inf.rate.slow.lowshed * 1e3    
  
  dff_plot <- dff %>%
    filter(!is.na(rate.per.100k))
  
  
  p <- ggplot(dff_plot) +
    #geom_sf(aes(fill = rate.per.100k), 
     #       color = "white", size = 0.1)+
    geom_sf(aes(fill = inf.rate.slow.highshed), 
             color = "white", size = 0.1) +
    scale_fill_gradient(
      name = "Minimum infection rate %\n(per 100,000 people)",
      na.value = "grey90",
      low = 'khaki1', high = 'goldenrod3',
      labels = function(x) {
        paste0(formatC(x, format = "f", digits = 4),
               "  (", formatC(x * 1000, format = "f", digits = 1), ")")
      }
    ) +
    geom_sf_text(aes(label = CFSAUID), 
                 size = 2, color = "grey20") +
    theme_minimal() +
    theme(
      legend.position = 'right',
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      # legend.title = element_text(face = "bold", size = 13),
      # legend.text = element_text(size = 12),
      # legend.key.size = unit(1.5, "lines")
    )
    
  p
}


figure_detection_fractiontype  <-function(df.fsa){
  
  # Prepare data
  df_plot <- df.fsa %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = str_to_title(wwtp),  # Capitalize the first letter
      liquid = wwtp.inf.rate.bio.deg.highshed,
      solid = wwtp.inf.rate.set.slow.highshed,
      total = wwtp.inf.rate.slow.highshed
    )
  
  # Add city row
  df_city <- df.fsa %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = "City",
      liquid = city.inf.rate.bio.deg.highshed,
      solid = city.inf.rate.set.slow.highshed,
      total = city.inf.rate.slow.highshed
    )
  
  # Combine WWTP and City
  df_combined <- bind_rows(df_plot, df_city)%>%
    distinct()
  
  # Pivot longer for plotting
  df_long <- df_combined %>%
    pivot_longer(cols = c(liquid, solid, total),
                 names_to = "fraction", values_to = "loss_value") %>%
    mutate(
      fraction = case_when(
        fraction == "liquid" ~ "Liquid",
        fraction == "solid" ~ "Solid",
        fraction == "total" ~ "Total"
      )
    )
  
  # Plot
  g <- ggplot(df_long, aes(x = wwtp, y = loss_value, fill = wwtp)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_grid(~fraction) +
    scale_fill_manual(values = c("North" = "#DEEBF7",  # light blue
                                 "South" = "#9ECAE1",  # medium blue
                                 "West"  = "#3182BD",   # dark blue 
                                 "City"  = 'goldenrod3')) +
    labs(
      x = "City and catchments of wastewater treatment plants",
      y = "Minimum infection rate\nfor detection (%)"
    ) +
    guides(fill = 'none') +
    theme(
      legend.position = 'none',
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(face = "bold", color = 'white'),
      strip.background = element_rect(fill = 'steelblue4'),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  g
  
  return(g)
}


figure_combo_detection <- function(df.fsa, filename) {
  
  g.map <- figure_map_inf_rate(df.fsa, shed.type = 'high') 
  g.bar = figure_detection_fractiontype(df.fsa)
  
  g.map <- g.map + 
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.margin = margin(t = 10, b = 10)
    )
  
  g.combo <- g.map + g.bar +
    plot_layout(ncol = 1, heights = c(2, 1)) +
    plot_annotation(tag_levels = 'A')
  # Save the plot
  #pdf(filename, width = 7, height = 8)
  #plot(g.combo)
  #dev.off()
  
  ggsave(
    filename,
    plot = g.combo,  
    width = 8,     # in inches (adjust as needed)
    height = 6,    # in inches
    dpi = 400,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}


# --- plot effective surveillance coverage
#---- PLOT EFFECTIVE POPULATION COVERAGE -----
plot_eff_pop_coverage <-function(df){
  
  # Prepare data
  df_plot <- df %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = str_to_title(wwtp),  # Capitalize the first letter
      liquid = eff.pop.bio.deg.path.wwtp,
      solid  = eff.pop.set.slow.wwtp,
      total  = eff.pop.total.slow.wwtp
    )
  
  # Add city row
  df_city <- df %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = "City",
      liquid = eff.pop.bio.deg.path.city,
      solid  = eff.pop.set.slow.city,
      total  = eff.pop.total.slow.city
    )
  
  # Combine WWTP and City
  df_combined <- bind_rows(df_plot, df_city)%>%
    distinct()
  
  # Pivot longer for plotting
  df_long <- df_combined %>%
    pivot_longer(cols = c(liquid, solid, total),
                 names_to = "fraction", values_to = "loss_value") %>%
    mutate(
      fraction = case_when(
        fraction == "liquid" ~ "Liquid",
        fraction == "solid" ~ "Solid",
        fraction == "total" ~ "Total"
      )
    )
  
  # Plot
  g <- ggplot(df_long, aes(x = wwtp, y = loss_value, fill = wwtp)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(
      aes(label = round(loss_value, 1)),
      vjust = -0.5,
      size = 4.5
    ) +
    facet_grid(~fraction) +
    scale_fill_manual(values = c("North" = "#DEEBF7",  # light blue
                                 "South" = "#9ECAE1",  # medium blue
                                 "West"  = "#3182BD",   # dark blue 
                                 "City"  = "#33a02c")) +
    scale_y_continuous(limits = c(0, 101)) +
    labs(
      x = "City and catchments of wastewater treatment plants",
      y = "Effective population coverage (%)"
    ) +
    guides(fill = 'none') +
    theme(
      legend.position = 'none',
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(face = "bold", color = 'white'),
      strip.background = element_rect(fill = 'steelblue4'),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  g
}


#Plot hydroalic residence time
figure_travel_time <- function(df, filename){
  
  df = df%>%
    mutate(
      wwtp = str_to_title(wwtp))  # Capitalize the first letter
      
  # df = df.flow <-- DEBUG
  stats <- df %>%
    group_by(wwtp) %>%
    summarise(mean_val = mean(path_hrt, na.rm = TRUE),
              median_val = median(path_hrt, na.rm = TRUE))
  
  # Plot histogram with mean and median lines, faceted by 'wwtp'
  g = ggplot(df, aes(x = path_hrt, fill = wwtp)) +
    geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
    geom_vline(data = stats,aes(xintercept = mean_val), color = "red", linetype = "dashed", size = 1) +
    geom_vline(data = stats,aes(xintercept = median_val), color = "purple", linetype = "dotted", size = 1) +
    scale_fill_manual(values = c("North" = "#DEEBF7",  # light blue
                                 "South" = "#9ECAE1",  # medium blue
                                 "West"  = "#3182BD")) +
    facet_wrap(~wwtp, scales = "free_y") +  # Facet by 'wwtp' with independent y scales
    labs(#title = "Hydraulic Residence Time Per Flow Path TO WWTP",
         x = "Traveling Time (hour)", y = "Counts of Flow Paths") +
    geom_text(data = stats, aes(x = mean_val, y = Inf, label = paste("Mean:", round(mean_val, 2))), 
              color = "red", vjust = 3.5, hjust = -0.4) +
    geom_text(data = stats, aes(x = median_val, y = Inf, label = paste("Median:", round(median_val, 2))), 
              color = "purple", vjust = 5, hjust = -0.4) +
    theme_minimal()+
    theme(
      strip.text = element_text(face = "bold"))  # Bold facet titles
  plot(g)
  
  ggsave(
    filename,
    plot = g,  
    width = 8,     # in inches (adjust as needed)
    height = 4,    # in inches
    dpi = 400,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}

#Plot hydroalic residence time
figure_shear_stress_pipe <- function(df, filename){
  
  df = df%>%
    mutate(
      wwtp = str_to_title(wwtp))  # Capitalize the first letter
  
  # df = df.flow <-- DEBUG
  stats <- df %>%
    group_by(wwtp) %>%
    summarise(mean_val = mean(conduit_ss, na.rm = TRUE),
              median_val = median(conduit_ss, na.rm = TRUE))
  
  # Plot histogram with mean and median lines, faceted by 'wwtp'
  g = ggplot(df, aes(x = conduit_ss, fill = wwtp)) +
    geom_histogram(color = "black") +
    geom_vline(data = stats,aes(xintercept = mean_val), color = "red", linetype = "dashed", size = 1) +
    geom_vline(data = stats,aes(xintercept = median_val), color = "purple", linetype = "dotted", size = 1) +
    scale_fill_manual(values = c("North" = "#DEEBF7",  # light blue
                                 "South" = "#9ECAE1",  # medium blue
                                 "West"  = "#3182BD")) +
    facet_wrap(~wwtp, scales = "free_y") +  # Facet by 'wwtp' with independent y scales
    labs(#title = "Hydraulic Residence Time Per Flow Path TO WWTP",
      x = expression("Bed Shear Stress (N/"*m^2*")"), y = "Counts of Conduits") +
    geom_text(data = stats, aes(x = mean_val, y = Inf, label = paste("Mean:", round(mean_val, 2))), 
              color = "red", vjust = 3.5, hjust = -0.4) +
    geom_text(data = stats, aes(x = median_val, y = Inf, label = paste("Median:", round(median_val, 2))), 
              color = "purple", vjust = 5, hjust = -0.6) +
    theme_minimal()+
    theme(
      strip.text = element_text(face = "bold"))  # Bold facet titles
  plot(g)
  
  ggsave(
    filename,
    plot = g,  
    width = 8,     # in inches (adjust as needed)
    height = 4,    # in inches
    dpi = 400,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}


#Plot Results from Solid Filtration Experiment Result
# Rav's results
figure_solid_filtration <- function(filename){
  
  # Data retrieved from NML lab's spreadsheet named:
  # 2025-03-10 Filtration PCR Worksheet_DTM.xlsx edits
  # at HC-SC PHAC-ASPC\Project Solid Filtration Experiment (NMLB - SOR) - Documents\Knowledge_Documentation\Filters\Test-Polycarbonate-filters-Vancouver
  df <- data.frame(
    pore_size = c(12, 8, 5, 2, 1, 0.8, 0.6, 0.4, 0.01),
    percent_greater = c(31, 35, 37, 45, 55, 55, 60, 65, 100)
  )
  
  # Plot with regression line
  g = ggplot(df, aes(x = pore_size, y = percent_greater)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    scale_x_log10() +   # log scale for pore size
    labs(
      x = expression("Filter pore size ("*mu*"m)"),
      y = "SARS-CoV-2 detections as % of total"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 12),
      panel.grid.minor = element_line(color = "grey85", size = 0.3, linetype = "dashed"),
      panel.grid.major = element_line(color = "grey80", size = 0.4, linetype = "dashed")
    )
  
  ggsave(
    filename,
    plot = g,  
    width = 8,     
    height = 4,   
    dpi = 300,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}





###### Sensitivity Analysis Figures
figure_sensitivity_loss_histo <- function(
    sim_df_loss_total,
    sim.df.loss.total.mean,
    filename = 'doc/figs/figure_sensitivity_loss_histo.png') {
  
  # settling plot
  g.set = get_histo_plot(sim_df_loss_total,
                         type = 'set')
  
  # total plot
  g.total = get_histo_plot(sim.df.loss.total.mean,
                           type = 'total')
  
  g.combo <- g.set + g.total +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')
  
  # Save combined plot
  ggsave(
    filename,
    plot = g.combo,  
    width = 7.5,    # in inches (adjust as needed)
    height = 8,    # in inches
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}

get_histo_plot <- function(sim_df,
                           type = 'set') {  # 'set' or 'total'
  
  if (type == 'set') {
    # Prepare and label data
    pre_dat <- sim_df %>%
      mutate(
        loss.set.path.slow = 1 - remain.set.slow,
        loss.set.path.fast = 1 - remain.set.fast,
        loss.set.path.hmg  = 1 - remain.set.hmg
      ) %>%
      pivot_longer(
        cols = c(loss.set.path.slow,
                 loss.set.path.fast,
                 loss.set.path.hmg),
        names_to = 'name',
        values_to = 'value'
      ) %>%
      dplyr::select(sim, node_id, wwtp, name, value) %>%
      mutate(
        nameplot = case_when(
          name == 'loss.set.path.fast' ~ 'Highly Skewed',
          name == 'loss.set.path.slow' ~ 'Skewed (Base)',
          name == 'loss.set.path.hmg'  ~ 'Homogeneous'
        ),
        nameplot = factor(nameplot,
                          levels = c("Skewed (Base)", "Highly Skewed", "Homogeneous"))
      )
    
    label.name <- 'Count of flow paths (500 simulations)'
    label.x <- 'Proportion of in-sewer settling/resuspension viral loss'
  }
  
  if (type == 'total') {
    # Prepare and label data
    pre_dat <- sim_df %>%
      pivot_longer(
        cols = c(mean.loss.total.slow,
                 mean.loss.total.fast,
                 mean.loss.total.hmg),
        names_to = 'name',
        values_to = 'value'
      ) %>%
      dplyr::select(node_id, wwtp, name, value) %>%
      mutate(
        nameplot = case_when(
          name == 'mean.loss.total.fast' ~ 'Highly Skewed',
          name == 'mean.loss.total.slow' ~ 'Skewed (Base)',
          name == 'mean.loss.total.hmg'  ~ 'Homogeneous'
        ),
        nameplot = factor(nameplot,
                          levels = c("Skewed (Base)", "Highly Skewed", "Homogeneous"))
      )
    
    label.name <- 'Count of flow paths \n(mean of 500 simulations)'
    label.x <- 'Proportion of in-sewer total viral loss'
  }
  
  dat <- pre_dat %>%
    mutate(
      wwtp = case_when(
        wwtp == "north" ~ "North",
        wwtp == "south" ~ "South",
        wwtp == "west"  ~ "West",
        TRUE ~ wwtp
      )
    )
  
  # Summary stats
  stats_df <- dat %>%
    group_by(wwtp, nameplot) %>%
    summarise(
      mean_val   = mean(value, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create plot
  g <- dat %>%
    ggplot(aes(x = value)) +
    geom_histogram(
      binwidth = 0.05,
      fill = "#7da0b1",
      color = "white",
      alpha = 0.9
    ) +
    geom_vline(data = stats_df, aes(xintercept = mean_val),
               color = "#1f78b4", linetype = "dashed", size = 1) +
    geom_vline(data = stats_df, aes(xintercept = median_val),
               color = "#ff7f00", linetype = "solid", size = 1) +
    facet_grid(wwtp ~ nameplot, scales = 'free_y') +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "#2b4c7e", color = NA),
      strip.text         = element_text(color = "white", face = "bold", size = 10),
      axis.title.x       = element_text(face = "bold", size = 12, margin = margin(t = 10)),
      axis.title.y       = element_text(face = "bold", size = 12, margin = margin(r = 10)),
      axis.text          = element_text(size = 10),
      plot.margin        = margin(12, 18, 10, 10)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_y_continuous(labels = scales::label_scientific()) +
    labs(
      x = label.x,
      y = label.name
    )
  
  return(g)
}


#- map of fsa and population loss
# by RNA distribution scenarios 
# facet wrap of three maps (A and B) 
figure_sensitivity_map_loss <- function(df.fsa,
                                        size.text,
                                        col.fsa,
                                        filename) {
  
  
  # DEBUG
  # df.fsa = sim.df.pop.loss
  
  g.virus = plot_map_loss(df.fsa, type='virus')
  g.popu = plot_map_loss(df.fsa, type='popu')
  
  g.combo <- g.virus + g.popu +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')
  
  # Save combined plot
  ggsave(
    filename,
    plot = g.combo,  
    width = 7,     # in inches (adjust as needed)
    height = 5,    # in inches
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}


figure_sensitivity_bar_plots <- function(df,filename) {
  
  
  # DEBUG
  # df = sim.df.inf.rate
  
  # scenarios of RNA distributions
  g.inf = plot_sensitivity_barplots(df, type='inf_rate')
  g.pop = plot_sensitivity_barplots(df, type='pop_coverage')
  
  g.combo <- g.pop + g.inf +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')
  
  # Save combined plot
  ggsave(
    filename,
    plot = g.combo,  
    width = 7,     # in inches (adjust as needed)
    height = 6,    # in inches
    dpi = 250,     # high-resolution (300–600 recommended for journals)
    units = "in",
    bg = "white",  # set background color
    device = "png"
  )
}


plot_map_loss <- function(df.fsa,
                          type = c("virus", "popu"),
                          size.text = 1.5,
                          col.fsa = "grey40") {
  type <- match.arg(type)
  
  # Prepare data based on type
  dff <- get_prepared_data(df.fsa, type = type) %>%
    dplyr::mutate(
      ViralScenario = factor(
        ViralScenario,
        levels = c("sk", "hsk", "hmg"),  # enforce order
        labels = c("Skewed (base)", "Highly Skewed", "Homogeneous")
      )
    )
  
  # Legend title, colors, labels depend on type
  if (type == "virus") {
    legend_title <- "Average in-sewer viral loss"
    low_col <- "grey96"
      high_col <- "red4"
        label_fun <- scales::percent_format()
        
        # emphasize 20–30%
        scale_fun <- scale_fill_gradientn(
          colours = c(low_col, "pink2", high_col),
          values  = scales::rescale(c(0, 0.3, 0.4, 1)),
          labels  = label_fun
        )
        
  } else {
    legend_title <- "Population equivalent loss"
    low_col <- "honeydew"
      high_col <- "forestgreen"
        label_fun <- scales::comma_format()
        
        # emphasize 4000–8000
        scale_fun <- scale_fill_gradientn(
          colours = c(low_col, "darkolivegreen2", high_col),
          values  = scales::rescale(c(0, 4000, 8000, max(dff$LossValue, na.rm = TRUE))),
          labels  = label_fun
        )
  }
  
  g <- ggplot(dff) +
    geom_sf(aes(fill = LossValue), color = col.fsa, linewidth = 0.1) +
    geom_sf_text(aes(label = CFSAUID), size = size.text, color = "black") +
    facet_wrap(~ ViralScenario, ncol = 3) +
    theme_minimal() +
    theme(
      legend.position = "left",
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.4, "cm"),
      legend.title = element_text(size = 9),
      strip.text  = element_text(size = 10, face = "bold", hjust = 0)
    ) +
    scale_fun +
    labs(fill = legend_title)
  
  return(g)
}



get_prepared_data <- function(df.fsa, type = c("virus", "popu")) {
  type <- match.arg(type)
  
  geometry_col <- "fsa.geometry"
  id_col       <- "CFSAUID"
  wwtp_col     <- "wwtp"
  
  if (type == "virus") {
    dff <- df.fsa %>%
      dplyr::select(
        !!wwtp_col := !!sym(wwtp_col),
        !!id_col   := !!sym(id_col),
        !!geometry_col := !!sym(geometry_col),
        sk  = fsa.ave.loss.total.slow,
        hsk = fsa.ave.loss.total.fast,
        hmg = fsa.ave.loss.total.hmg
      )
  } else {
    dff <- df.fsa %>%
      dplyr::select(
        !!wwtp_col := !!sym(wwtp_col),
        !!id_col   := !!sym(id_col),
        !!geometry_col := !!sym(geometry_col),
        sk  = pop.signal.loss.total.slow.fsa,
        hsk = pop.signal.loss.total.fast.fsa,
        hmg = pop.signal.loss.total.hmg.fsa
      )
  }
  
  dff %>%
    distinct() %>%
    st_set_geometry(geometry_col) %>%
    tidyr::pivot_longer(
      cols = c("sk", "hsk", "hmg"),
      names_to = "ViralScenario",
      values_to = "LossValue"
    ) %>%
    sf::st_as_sf()
}





plot_sensitivity_barplots <- function(df, type) {

  # ---- Define mapping for each type ----
  config <- list(
    inf_rate = list(
      cols_wwtp = c(sk   = "wwtp.inf.rate.set.slow.highshed",
                    hmg  = "wwtp.inf.rate.set.hmg.highshed",
                    hsk  = "wwtp.inf.rate.set.fast.highshed"),
      cols_city = c(sk   = "city.inf.rate.set.slow.highshed",
                    hmg  = "city.inf.rate.set.hmg.highshed",
                    hsk  = "city.inf.rate.set.fast.highshed"),
      y_name    = "Minimum infection rate \nfor detection (%)",
      city_col  = "goldenrod3",
      lm_range = c(0, 0.4),
      dec = 3
    ),
    pop_coverage = list(
      cols_wwtp = c(sk   = "eff.pop.set.slow.wwtp",
                    hmg  = "eff.pop.set.hmg.wwtp",
                    hsk  = "eff.pop.set.fast.wwtp"),
      cols_city = c(sk   = "eff.pop.set.slow.city",
                    hmg  = "eff.pop.set.hmg.city",
                    hsk  = "eff.pop.set.fast.city"),
      y_name    = "Effective population \ncoverage (%)",
      city_col  = "forestgreen",
      lm_range = c(0, 55),
      dec = 1
    )
  )
  
  # ---- Grab config for chosen type ----
  cfg <- config[[type]]
  if (is.null(cfg)) stop("type must be 'inf_rate' or 'pop_coverage'")
  
  # ---- Prepare WWTP rows ----
  df_plot <- df %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = str_to_title(wwtp),
      !!!set_names(syms(cfg$cols_wwtp), names(cfg$cols_wwtp))
    )
  
  # ---- Prepare City row ----
  df_city <- df %>%
    st_drop_geometry() %>%
    distinct() %>%
    transmute(
      wwtp = "City",
      !!!set_names(syms(cfg$cols_city), names(cfg$cols_city))
    )
  
  # ---- Combine and pivot ----
  df_long <- bind_rows(df_plot, df_city) %>%
    distinct() %>%
    pivot_longer(cols = c(sk, hmg, hsk),
                 names_to = "fraction", values_to = "loss_value") %>%
    mutate(
      loss_value = as.numeric(loss_value),
      fraction = factor(fraction,
                        levels = c("sk", "hsk", "hmg"),
                        labels = c("Skewed (base)", "Highly Skewed", "Homogeneous"))
    )
  
  # ---- Plot ----
  ggplot(df_long, aes(x = wwtp, y = loss_value, fill = wwtp)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = round(loss_value, cfg$dec)), vjust = -0.5, size = 3.5) +
    facet_grid(~fraction) +
    scale_fill_manual(values = c("North" = "#DEEBF7",
                                 "South" = "#9ECAE1",
                                 "West"  = "#3182BD",
                                 "City"  = cfg$city_col)) +
    scale_y_continuous(limits = cfg$lm_range) +
    labs(x = "City and catchments of wastewater treatment plants",
         y = cfg$y_name) +
    guides(fill = "none") +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(face = "bold", color = "white", size = 10),
      strip.background = element_rect(fill = "steelblue4"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}
