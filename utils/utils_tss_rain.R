
# Script for helper functions on TSS and Rain data from Winnipeg
#----------------------------------------------------------------
  
suppressMessages({
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(lubridate)
    library(stringr)
    library(patchwork)
    library(snowfall)
    library(sf)
    library(readxl)
    library(grid)
    library(gridExtra)
})



get_tss <- function(){
  
  # --- Read Plant's TSS and Flow Data
  name = paste0('data/City Of Winnipeg Treatment Plant TSS Data.xlsx')
  TSS.North = read_excel(name, sheet = 'North End Plant')
  TSS.South = read_excel(name, sheet = 'South End Plant')
  TSS.West = read_excel(name, sheet = 'West End Plant')
  
  # --- clean and Merge
  # Units are mg/l for TSS and ML for flow
  N = TSS.North%>%
    filter(!is.na(Time))%>%
    rename(flow = "SCADA-NEWPCC.150MAC04",
           TSS = "LIMS-SM.NEWPCC.N_RSEW24.TSS_Report.SOLIDS")%>%
    mutate(wwtp = 'North')
  
  S = TSS.South%>%
    filter(!is.na(Time))%>%
    rename(flow = "SCADA-SEWPCC.125SAC06",
           TSS = "LIMS-SM.SEWPCC.S_RSEWCOMP.TSS_Report.SOLIDS")%>%
    mutate(wwtp = 'South')
  
  W = TSS.West%>%
    filter(!is.na(Time))%>%
    rename(flow = "SCADA-WEWPCC.WDM600FT",
           TSS = "LIMS-SM.WEWPCC.W_RSEW.TSS_Report.SOLIDS")%>%
    mutate(wwtp = 'West')
  
  dat = rbind(N,S)
  dat = rbind(dat,W)%>%
    mutate(flow = as.numeric(flow),
           TSS = as.numeric(TSS),
           Time = as.Date(as.POSIXct(Time), 'GMT'))
  
  # remove zero TSS
  dat = dat%>%
    filter(TSS > 0)
  return(dat)
} 


get_rain <- function(){
  
  # --------- Retrieve Rain Data from prcip.dat
  
  rainData <- read.table('data/precip.dat', header = FALSE, skip = 1)
  # Reformat header and column names
  vars.rain = c("Station","Year","Month","Day","Hour","Minute","rain_mm")
  names(rainData) <- vars.rain
  
  rainData = rainData |> 
    mutate(Time = ymd_hm(paste(Year, Month, Day, Hour, Minute, sep='-')))%>%
    dplyr::select(Time, rain_mm)
  
  return(rainData)
}



plot_tss_rain <- function(){
  
  dat = get_tss()
  rainData = get_rain()
  
  #---------- Plot Winnipeg TSS data
  
  # time series
  g1 = dat %>%
    filter(TSS != 0)%>%
    pivot_longer(-c(wwtp,Time)) %>% 
    ggplot(aes(x=Time,y=value, color=name)) + 
    geom_point(size=0.5) +
    geom_line() +
    facet_wrap(~wwtp, scales = 'free_y')+
    scale_y_log10()+
    scale_color_brewer(palette='Set2')+
    #xaxis + th + xlab('')+
    ggtitle('Influent Flow and TSS - Winnipeg')
  g1
  
  #------------ Plot Precipitation and TSS together
  g2 = rainData%>%
    ggplot()+
    geom_line(aes(x=Time, y=rain_mm), color = 'blue') +
    theme_bw() +
    ylab("Precip (mm)") +
    scale_y_reverse()
  g2
  
  g1 <- ggplot_gtable(ggplot_build(g1))
  g2 <- ggplot_gtable(ggplot_build(g2))
  maxWidth = unit.pmax(g1$widths[2:3], g2$widths[2:3])
  
  g1$widths[2:3] <- maxWidth
  g2$widths[2:3] <- maxWidth
  g = grid.arrange(g2, g1, ncol = 1, heights = c(1, 3))  
  
  return(g)
}


plot_hist_tss <- function(dat, type='dry+wet'){
  
  if(type == 'dry') dat = dat%>%
      filter(weather == 'dry')
    
  # histogram plot
  g.tss = dat%>%
    filter(TSS !=0,
           TSS < 800)%>%
    ggplot(aes(x=TSS))+
    geom_histogram(aes(y= ..density..),
                   color = 1, fill = 'pink1')+
    geom_density(size = 1.5, color = 'brown')+
    facet_wrap(~wwtp)
  
  return(g.tss)
}

#--- filter wet-weather tss values
generate_weather_condition <-function(rain.threshold = 5){
  
  #----- retrieve data
  tss  = get_tss()
  rain = get_rain()
  
  #----- merge
  
  # convert 15-minutes rain values to daily 
  rain2 = rain%>%
    mutate(date = as.Date(as.POSIXct(Time), 'GMT'))%>%
    group_by(date)%>%
    summarise(rain_daily = sum(rain_mm))
  
  tss = tss%>%
    mutate(date = Time)%>%
    dplyr::select(-Time)
  
  # combine tss and rain
  dat = left_join(tss, rain2, by = 'date')
  
  # identify wet weather evernts
  threshold = rain.threshold 
  #tss_threshold
  #lag_time
  #dry_period 
  
  dat = dat%>%
    mutate(weather = ifelse(rain_daily > threshold, 'wet', 'dry'))
  
  return(dat)
}

get_tss_dry <- function(){
  
  dat.tss = generate_weather_condition()
  
  #remove wet weather dates
  tmp = dat.tss%>%
    filter(weather == 'dry')
  
  return(tmp)
}



# Estimate statistical summary for each plant
get_summary_stat <- function(tss){
  
  summary_stat = tss %>%
    group_by(wwtp) %>%
    summarise(
      Q1 = quantile(TSS, 0.25),
      Q3 = quantile(TSS, 0.75),  
      Mean = mean(TSS),
      Median = median(TSS),
      StdDev = sd(TSS))
  
  return(summary_stat)
}

# Helper function for Main script of simulation
# Load tss data and estimate normal distribution
get_tss_dist <- function(){
  
  #----- retrieve data
  tss  = get_tss()
  
  # Only consider data from Dec-Feb of each year
  tss.dry = tss%>%
    filter(format(Time, "%m") %in% c("12", "01", "02"))
  
  # statistical summary for each plant
  quant = get_summary_stat(tss.dry)
  
  return(quant)
}



# Helper function for Main script of simulation
# Load tss data and estimate normal distribution
get_tss_winter_months <- function(){
  
  #----- retrieve data
  tss  = get_tss()
  
  # Only consider data from Dec-Feb of each year
  tss.dry = tss%>%
    filter(format(Time, "%m") %in% c("12", "01", "02"))%>%
    mutate(load = flow*TSS) # in (kg)
  
  return(tss.dry)
}
  