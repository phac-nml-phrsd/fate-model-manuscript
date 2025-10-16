#------------
#------------
# Helper functions for parameters.R
#=====================================

suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(MASS)
})


ggplot2::theme_set(theme_bw())


#' Custom Beta distribution
#' 
#' This is used to conveniently parameterized the shape 
#' of the distributions for suspended solids and RNA adsorption.
#'
#' @param shape1 non-negative parameters of the Beta distribution as expected by `dbeta()`
#' @param shape2 non-negative parameters of the Beta distribution as expected by `dbeta()`
#' @param values Numerical vector representing the values of the parameter.
#'
#' @return Dataframe of values and associated probabilities
#' 
beta_custom <- function(shape1, shape2, values) {
  # Rescale to [0;1] for the standard Beta distribution
  x = (values - min(values)) / max(values)
  # Retrieve Beta distribution with the requested shapes
  b = dbeta(x = x, shape1 = shape1, shape2 = shape2)  
  # save as a dataframe 
  res = data.frame( x = values, y = b / sum(b))
  return(res)
}


if(0){ # DEBUG/TEST code
  values = c(1, 50, 300, 500, 800, 1000)
  d1 = beta_custom(shape1 = 1, shape2 = 3, values)|> mutate(name='A')
  d2 = beta_custom(shape1 = 1, shape2 = 7, values)|> mutate(name='B')
  d = bind_rows(d1, d2)
  ggplot(d, aes(x=x,y=y, color = name))+geom_line() + geom_point()
}


# Function to fit Lognormal distribution and return parameters
# for each wwtp 
fit_lognormal <- function(data, type) {
  cat("\nFitting Lognormal distribution for WWTP Type:", type, "\n")
  
  # Extract TSS values for the given WWTP type
  TSS_values <- data$TSS
  
  # Fit the Lognormal distribution using Maximum Likelihood Estimation (MLE)
  fit <- fitdistr(TSS_values, "lognormal")
  
  # Extract estimated parameters
  meanlog <- fit$estimate["meanlog"]
  sdlog <- fit$estimate["sdlog"]
  
  return(c(meanlog = meanlog, sdlog = sdlog))
}


#----- FITTING a Gamma distribution and return shape and rate 
fit_gamma_dist <- function(data_points, n_samples = 1000) {
  # Fit the Gamma distribution to the data
  fit <- fitdistr(data_points, "gamma")
  
  # Extract estimated parameters
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  
  return(c(shape = shape, rate = rate))
}

#------ PLOT fitted Gamma distribution to data
plot_gamma_dist <- function(data_points, n_samples = 1000) {
  # Fit the Gamma distribution to the data
  fit <- fitdistr(data_points, "gamma")
  
  # Extract estimated parameters
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  
  
  set.seed(123)  # For reproducibility
  # Generate samples from the fitted Gamma distribution
  simulated_data <- rgamma(n_samples, shape = shape, rate = rate)
  
  # Plot the histogram of the simulated Gamma distribution
  hist(data_points, breaks = 30, col = "lightblue", freq = FALSE, 
       main = "Fitted Gamma Distribution", xlab = "Value")
  
  # Overlay density curve
  lines(density(simulated_data), col = "darkblue", lwd = 2)
  
  return(list(shape = shape, rate = rate, simulated_data = simulated_data))
}

#------ PLOT fitted LogNormal distribution to TSS data
plot_fit_lognormal <- function(data, type){
  
  cat("\nFitting Lognormal distribution for WWTP Type:", type, "\n")
  
  # Extract TSS values for the given WWTP type
  TSS_values <- data$TSS
  
  # Fit the Lognormal distribution using Maximum Likelihood Estimation (MLE)
  fit <- fitdistr(TSS_values, "lognormal")
  
  # Extract estimated parameters
  meanlog <- fit$estimate["meanlog"]
  sdlog <- fit$estimate["sdlog"]
  
  # Generate simulated data from the fitted Lognormal distribution
  simulated_data <- rlnorm(1000, meanlog = meanlog, sdlog = sdlog)
  
  # Plot original vs fitted Lognormal distribution
  hist(TSS_values, breaks = 30, col = "gray", freq = FALSE,
       main = paste("TSS - Lognormal Fit (WWTP Type:", type, ")"), xlab = "TSS Value")
  lines(density(simulated_data), col = "blue", lwd = 2)
  
  return(c(meanlog = meanlog, sdlog = sdlog))
}

plot_TSS_logNormal <- function(){
  
  tss.dry = get_tss_winter_months()
  
  # PLOT - Fit Lognormal to each WWTP type 
  # and return plots
  lognormal_fits <- tss.dry %>%
    group_by(wwtp) %>%
    group_split() %>%
    setNames(unique(tss.dry$wwtp)) %>%
    lapply(function(d) plot_fit_lognormal(d, unique(d$wwtp)))
}



# Produce distribution around mean value of limited data 
create_sample_dist <- function(data){
  
  # Set the sample size (number of points in each sample)
  sample_size = 3
  
  # Set the number of resamples to simulate the sampling distribution
  n_resamples = 1000
  
  # Simulate the sampling distribution of the mean
  set.seed(123)  # For reproducibility
 
  sampling_means = replicate(n_resamples, mean(sample(data, size = sample_size, replace = TRUE)))
  
  # View summary statistics of the sampling distribution
  summary(sampling_means)
  
  return(sampling_means)
}

# Normalize biofilm sorption rates to A/V ratio 
av_normalize <- function(k, av){
  
  kappa.av = k / av 
  
  return(kappa.av)
}

