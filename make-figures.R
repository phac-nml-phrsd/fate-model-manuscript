
suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(snowfall)
  library(sf)
  library(ggridges)
})

source('utils/utils_plot.R')
source('utils/utils_figures.R')

start <- Sys.time()
message('Start creating figures.')
# ===== FIGURES MANUSCRIPT ======
# Main text

# Histogram
figure_loss_histo (
  sim_df_loss_total,
  filename = 'figs/figure_loss_histo.png')

# Conduit map and loss distribution
figure_combo_lossmap (
  sim.df.loss.total.mean,
  iw.polygon, 
  iw.wwtp,
  filename = 'figs/figure_combo_lossmap.png')

# FSA map and loss map and bar plot
figure_combo_map_loss_fsa(
  sim.df.pop.loss, 
  filename = 'figs/figure_combo_map_loss_pop.png')


# Infection Limit map and bar plot
figure_combo_detection(sim.df.inf.rate,
  filename = 'figs/figure_combo_detection.png')



#========= Appendix
# PIPES & PATHS
# for every particle classes
# NO RNA dist. applied yet
# impact of settling velocity and critical shear stress
plot_settling_pipe(sim_df_loss_solid,df_labels,
    filename = 'figs/appendix_settling_pipe.pdf')
plot_settling_path(sim_df_loss_solid,loss.type = 'settling raw', df_labels,
    filename = 'figs/appendix_settling_path.pdf')

#figure_travel_time(df.flow,
#                   filename = 'figs/figure_travel_time.png')

#figure_shear_stress_pipe(df.flow,
#                    filename = 'figs/figure_shear_stress.png')

# by fraction
figure_map_loss_fsa_fraction(sim.df.pop.loss, 
                             size.text = 2.5, 
                             col.fsa = "white",
                             filename = 'figs/figure_map_loss_fsa_fraction.png')


figure_solid_filtration(filename = 'figs/figure_solid_filtration.png')


#========================== Sensitivity Analysis (RNA scenarios)
# Histogram of sewer loss
figure_sensitivity_loss_histo(
  sim_df_loss_total,
  sim.df.loss.total.mean,
  filename = 'figs/figure_sensitivity_loss_histo.png')


# map neighborhood loss by RNA scenarios
figure_sensitivity_map_loss(sim.df.pop.loss, 
                            size.text = 1.5, 
                            col.fsa = "white",
                            filename = 'figs/figure_sensitivity_loss_map.png')

# bar plots for population coverage and infection rates
figure_sensitivity_bar_plots(sim.df.inf.rate, 
                             filename = 'figs/figure_sensitivity_barplots.png')


message('Figures are all created.')
end <- Sys.time()

cat("Total runtime:", round(difftime(end, start, units = "secs"), 2), "seconds\n")
