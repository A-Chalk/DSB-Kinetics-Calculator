# title: "Cas_Solver+ (1.7)"
# author: "Alexander Chalk"
# date: "24th October, 2024"

rm(list = ls())
cas_solver <- "1.7"

# Import necessary libraries
library(deSolve)
library(tidyverse)
library(readxl)
library(foreach)
library(openxlsx)
library(doParallel)

# Track script execution time
start_time <- Sys.time()

# Register parallel backend using most available cores
cl <- makeCluster(detectCores() - 3)
registerDoParallel(cl)

# Number of bootstrap datasets (adjust as needed)
num_bootstrap_reps <- 10
export_results <- TRUE  # Set to TRUE to export results

# Adjust WT values to isolate edited cell population
untransfected <- 0  # Input untransfected %

################################################# Raw Data Processing and Bootstrapping #############################################
# Import raw data
Kinetics_data <- read_excel("data/Cas9_Solver_Input.xlsx")
num_reps <- 3  # Specify number of replicates in raw data (should be the same for each measured population)

# Check for the presence of MMEJ or TI (Targeted Integration) datasets
has_MMEJ <- any(grepl("MMEJ", colnames(Kinetics_data)))
has_TI <- any(grepl("TI", colnames(Kinetics_data)))

# Dynamically define columns based on the presence of MMEJ/TI
RNP_columns <- c(paste0("WT", 1:num_reps), paste0("DSB", 1:num_reps), paste0("Indel", 1:num_reps), paste0("LD", 1:num_reps))
RNP_IN_columns <- paste0(RNP_columns, "_IN")

if (has_MMEJ) {
  RNP_columns <- c(RNP_columns, paste0("MMEJ", 1:num_reps))
  RNP_IN_columns <- c(RNP_IN_columns, paste0("MMEJ", 1:num_reps, "_IN"))
}
if (has_TI) {
  RNP_columns <- c(RNP_columns, paste0("TI", 1:num_reps))
  RNP_IN_columns <- c(RNP_IN_columns, paste0("TI", 1:num_reps, "_IN"))
}

# Split the data into 'RNP' and 'RNP + Inhibitors (IN)'
split_data <- function(data, columns) {
  select(data, Hours, all_of(columns))
}
RNP_data <- split_data(Kinetics_data, RNP_columns)
RNP_IN_data <- split_data(Kinetics_data, RNP_IN_columns)

# Adjust column names for RNP_IN_data
colnames(RNP_IN_data) <- sub("_IN", "", colnames(RNP_IN_data))
Timepoints <- RNP_data$Hours

# Group replicates into data frames based on number of replicates
group_replicates <- function(data) {
  lapply(1:num_reps, function(i) {
    df <- data.frame(
      Timepoints,
      WT_observed = data[[paste0("WT", i)]],
      DSB_observed = data[[paste0("DSB", i)]],
      Indel_observed = data[[paste0("Indel", i)]],
      LD_observed = data[[paste0("LD", i)]]
    )
    if (has_MMEJ) df$MMEJ_observed <- data[[paste0("MMEJ", i)]]
    if (has_TI) df$TI_observed <- data[[paste0("TI", i)]]
    df
  })
}
rep_data_RNP <- group_replicates(RNP_data)
rep_data_RNP_IN <- group_replicates(RNP_IN_data)

# Calculate the average of each variable across replicates
calculate_average <- function(rep_data) {
  variables <- names(rep_data[[1]])[-1]  # Exclude Timepoints
  df <- data.frame(Timepoints)
  for (var in variables) {
    df[[var]] <- rowMeans(do.call(cbind, lapply(rep_data, `[[`, var)))
  }
  df
}
average_data_RNP <- calculate_average(rep_data_RNP)
average_data_RNP_IN <- calculate_average(rep_data_RNP_IN)

# Adjust WT observed values based on untransfected population
adjust_WT <- function(rep_data) {
  for (i in 1:num_reps) {
    rep_data[[i]]$WT_observed <- (rep_data[[i]]$WT_observed - untransfected) / (rep_data[[i]]$WT_observed[1] / 100)
  }
  rep_data
}
rep_data_RNP <- adjust_WT(rep_data_RNP)
rep_data_RNP_IN <- adjust_WT(rep_data_RNP_IN)

# Filter RNP_IN data to the first 24 hours (before inhibitors are removed)
rep_data_RNP_IN <- lapply(rep_data_RNP_IN, function(df) filter(df, Timepoints <= 24))

# Generate bootstrap datasets
generate_bootstrap_datasets <- function(rep_data, num_bs_reps) {
  variables <- names(rep_data[[1]])[-1]  # Exclude Timepoints
  lapply(1:num_bs_reps, function(s) {
    bootstrap_dataset <- data.frame(Timepoints = rep_data[[1]]$Timepoints)
    for (var in variables) {
      bootstrap_dataset[[var]] <- sapply(1:nrow(bootstrap_dataset), function(t) {
        sample(sapply(rep_data, function(df) df[t, var]), 1)
      })
    }
    bootstrap_dataset
  })
}
bootstrap_datasets <- list(
  RNP = generate_bootstrap_datasets(rep_data_RNP, num_bootstrap_reps),
  RNP_IN = generate_bootstrap_datasets(rep_data_RNP_IN, num_bootstrap_reps)
)

################################################### ODE Kinetics Model #######################################################

# Define the ODE function for RNP kinetics
my_ode <- function(t, state, parms, has_MMEJ, has_TI) {
  with(as.list(c(state, parms)), {
    D <- 1 - 2^(-(t / delay)) # Alternative delay formula to play with, just make sure only one D is active
    dWTdt <- -(k_dsb * D * WT) + (k_pr * DSB) # Wildtype Rate
    dIndeldt <- k_in * DSB # Indel Rate
    dLDdt <- k_ld * DSB # Large Deletion Rate
    dPRdt <- (k_pr * DSB) - (k_dsb * D * PR) # Precise Repair Rate
    dC_DSBdt <- k_dsb * D * WT # Cumulative DSBs (total DSBs over the time course)
    dC_PRdt <- k_pr * DSB # Cumulative PR
    dDSBdt <- (k_dsb * D * WT) - (k_pr + k_in + k_ld) * DSB # DSB Rate
    dMHdt <- if (has_MMEJ) k_mh * DSB else 0 # MMEJ Rate if MMEJ data exists
    dTIdt <- if (has_TI) k_ti * DSB else 0 # TI Rate if TI data exists
    if (has_MMEJ) dDSBdt <- dDSBdt - k_mh * DSB # Adjust DSB Rate if k_mh exists
    if (has_TI) dDSBdt <- dDSBdt - k_ti * DSB # Adjust DSB Rate if k_ti exists
    
    derivatives <- c(dWTdt, dDSBdt, dIndeldt, dLDdt, dPRdt, dC_DSBdt, dC_PRdt)
    if (has_MMEJ) derivatives <- c(derivatives, dMHdt)
    if (has_TI) derivatives <- c(derivatives, dTIdt)
    list(derivatives)
  })
}

# Define initial values
init_values <- function(observed_data) {
  init <- c(
    WT = observed_data$WT_observed[1], 
    DSB = observed_data$DSB_observed[1], 
    Indel = observed_data$Indel_observed[1], 
    LD = observed_data$LD_observed[1], 
    PR = 0,
    C_DSB = 0,
    C_PR = 0,
    MMEJ = if (has_MMEJ) observed_data$MMEJ_observed[1] else NULL,
    TI = if (has_TI) observed_data$TI_observed[1] else NULL
  )
  init
}

# ODE residuals function to take bootstrap data as an argument with tolerances for the solver
ode_residuals <- function(parameters, observed_data_RNP, observed_data_RNP_IN, has_MMEJ, has_TI) {
  # Shared parameters for both datasets
  parms_shared <- list(k_dsb = parameters["k_dsb"], delay = parameters["delay"])
  
  # Function to build parameters for each dataset
  build_parms <- function(pars, suffix) {
    list(
      k_pr = pars[paste0("k_pr_", suffix)],
      k_in = pars[paste0("k_in_", suffix)],
      k_ld = pars[paste0("k_ld_", suffix)],
      k_mh = if (has_MMEJ) pars[paste0("k_mh_", suffix)] else NULL,
      k_ti = if (has_TI) pars[paste0("k_ti_", suffix)] else NULL
    )
  }
  
  # Build parameters for RNP and RNP_IN
  parms_RNP <- c(parms_shared, build_parms(parameters, "RNP"))
  parms_RNP_IN <- c(parms_shared, build_parms(parameters, "RNP_IN"))
  
  
  # Function to run ODE solver and return simulated data
  run_ode_solver <- function(init, times, parms) {
    tryCatch(
      as.data.frame(ode(y = init, times = times, func = function(t, state, parms) my_ode(t, state, parms, has_MMEJ, has_TI), parms = parms, rtol = 1e-12, atol = 1e-12)),
      error = function(e) return(NULL)
    )
  }
  
  # Run the solver for both datasets
  out_RNP <- run_ode_solver(init_values(observed_data_RNP), observed_data_RNP$Timepoints, parms_RNP)
  out_RNP_IN <- run_ode_solver(init_values(observed_data_RNP_IN), observed_data_RNP_IN$Timepoints, parms_RNP_IN)
  
  # Check if either solver fails and return a penalty
  if (is.null(out_RNP) || is.null(out_RNP_IN)) return(1e10)
  
  # Define column names dynamically
  col_names <- c("Timepoints", "WT", "DSB", "Indel", "LD", "PR", "C_DSB", "C_PR")
  if (has_MMEJ) col_names <- c(col_names, "MMEJ")
  if (has_TI) col_names <- c(col_names, "TI")
  
  colnames(out_RNP) <- col_names
  colnames(out_RNP_IN) <- col_names
  
  # Function to calculate residuals
  calculate_residuals <- function(observed_data, simulated_data) {
    vars <- c("WT", "DSB", "Indel", "LD", if (has_MMEJ) "MMEJ" else NULL, if (has_TI) "TI" else NULL)
    unlist(lapply(vars, function(var) observed_data[[paste0(var, "_observed")]] - simulated_data[[var]]))
  }
  
  # Combine residuals from both datasets
  residuals <- c(calculate_residuals(observed_data_RNP, out_RNP), calculate_residuals(observed_data_RNP_IN, out_RNP_IN))
  return(sum(residuals^2))
}


################################################### Parameter Solver #########################################################

# Initial parameter guess
init_params <- c(k_dsb = 1, delay = 1, k_pr_RNP = 1, k_in_RNP = 0.1, k_ld_RNP = 0.1, k_pr_RNP_IN = 0.1, k_in_RNP_IN = 0.1, k_ld_RNP_IN = 0.1)
if (has_MMEJ) init_params <- c(init_params, k_mh_RNP = 0.1, k_mh_RNP_IN = 0.1)
if (has_TI) init_params <- c(init_params, k_ti_RNP = 0.1, k_ti_RNP_IN = 0.1)

# Lower and upper bounds
lower_bounds <- c(k_dsb = 0, delay = 0, k_pr_RNP = 0, k_in_RNP = 0, k_ld_RNP = 0, k_pr_RNP_IN = 0, k_in_RNP_IN = 0, k_ld_RNP_IN = 0)
upper_bounds <- c(k_dsb = 1, delay = 1, k_pr_RNP = 1, k_in_RNP = 1, k_ld_RNP = 1, k_pr_RNP_IN = 1, k_in_RNP_IN = 1, k_ld_RNP_IN = 1)
if (has_MMEJ) {
  lower_bounds <- c(lower_bounds, k_mh_RNP = 0, k_mh_RNP_IN = 0)
  upper_bounds <- c(upper_bounds, k_mh_RNP = 1, k_mh_RNP_IN = 1)
}
if (has_TI) {
  lower_bounds <- c(lower_bounds, k_ti_RNP = 0, k_ti_RNP_IN = 0)
  upper_bounds <- c(upper_bounds, k_ti_RNP = 1, k_ti_RNP_IN = 1)
}

# Save initial parameters and bounds to the final export
params_data <- list(init_params = init_params, lower_bounds = lower_bounds, upper_bounds = upper_bounds)

# Parallel loop for bootstrap parameter estimation with specified packages
bootstrap_params_estimates <- foreach(i = 1:num_bootstrap_reps, .combine = rbind, .packages = "deSolve") %dopar% {
  # Select the current bootstrap datasets
  observed_data_RNP <- bootstrap_datasets$RNP[[i]]
  observed_data_RNP_IN <- bootstrap_datasets$RNP_IN[[i]]
  
  # Optimizing for the RNP_IN dataset first
  fit_combined <- optim(par = init_params, fn = function(pars) ode_residuals(pars, observed_data_RNP, observed_data_RNP_IN, has_MMEJ, has_TI),
                        method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  fit_combined$par
}

# Separate the parameters estimates
param_names <- c("k_dsb", "delay", "k_pr_RNP", "k_in_RNP", "k_ld_RNP", "k_pr_RNP_IN", "k_in_RNP_IN", "k_ld_RNP_IN")
if (has_MMEJ) param_names <- c(param_names, "k_mh_RNP", "k_mh_RNP_IN")
if (has_TI) param_names <- c(param_names, "k_ti_RNP", "k_ti_RNP_IN")

colnames(bootstrap_params_estimates) <- param_names

# Calculate the mean and standard deviations of the bootstrap parameter estimates
mean_parameters <- colMeans(bootstrap_params_estimates, na.rm = TRUE)
sd_parameters <- apply(bootstrap_params_estimates, 2, sd, na.rm = TRUE)

# Calculate Standard Error and 95% CI
se_parameters <- sd_parameters / sqrt(num_bootstrap_reps)
t_score <- qt(0.975, df = num_bootstrap_reps - 1)
ci_parameters <- rbind(
  pmax(mean_parameters - t_score * se_parameters, 0), # Ensuring CI Lower is not below zero
  mean_parameters + t_score * se_parameters
)
rownames(ci_parameters) <- c("CI Lower", "CI Upper")

# Convert mean_parameters to a named list for the ODE function
mean_parameters_list <- setNames(as.list(mean_parameters), names(mean_parameters))

# Use the average observed data for comparison
observed_data_RNP <- average_data_RNP
observed_data_RNP_IN <- average_data_RNP_IN

# Define a sequence of minute-length time points spanning the entire timeseries
all_timepoints <- seq(from = 0, to = max(Timepoints), by = 1/60)

# Define initial values once and reuse
init_values_RNP <- init_values(average_data_RNP)
init_values_RNP_IN <- init_values(average_data_RNP_IN)

# In the ODE solver call, reuse these initial values
simulate_ode <- function(init_values, parms, has_MMEJ, has_TI, timepoints) {
  as.data.frame(ode(y = init_values, times = timepoints, func = function(t, state, parms) my_ode(t, state, parms, has_MMEJ, has_TI), parms = parms, rtol = 1e-6, atol = 1e-6))
}

column_names <- c("Timepoints", "WT", "DSB", "Indel", "LD", "PR", "C_DSB", "C_PR")
if (has_MMEJ) {
  column_names <- c(column_names, "MMEJ")
}
if (has_TI) {
  column_names <- c(column_names, "TI")
}

# Run the ODE solver for the average data
simulated_values_avg_df_RNP <- simulate_ode(init_values_RNP, list(
  k_dsb = mean_parameters_list$k_dsb,
  delay = mean_parameters_list$delay,
  k_pr = mean_parameters_list$k_pr_RNP,
  k_in = mean_parameters_list$k_in_RNP,
  k_ld = mean_parameters_list$k_ld_RNP,
  k_mh = if (has_MMEJ) mean_parameters_list$k_mh_RNP else NULL,
  k_ti = if (has_TI) mean_parameters_list$k_ti_RNP else NULL
), has_MMEJ, has_TI, all_timepoints)

names(simulated_values_avg_df_RNP) <- c(column_names)

# Run the ODE solver for the first 24 hours using the RNP_IN coefficients
simulated_values_24h_RNP_IN <- simulate_ode(init_values_RNP_IN, list(
  k_dsb = mean_parameters_list$k_dsb,
  delay = mean_parameters_list$delay,
  k_pr = mean_parameters_list$k_pr_RNP_IN,
  k_in = mean_parameters_list$k_in_RNP_IN,
  k_ld = mean_parameters_list$k_ld_RNP_IN,
  k_mh = if (has_MMEJ) mean_parameters_list$k_mh_RNP_IN else NULL,
  k_ti = if (has_TI) mean_parameters_list$k_ti_RNP_IN else NULL
), has_MMEJ, has_TI, seq(from = 0, to = 24, by = 1/60))

names(simulated_values_24h_RNP_IN) <- names(simulated_values_avg_df_RNP)

# Use the final state from the first simulation as the initial state for the second simulation
if (max(Timepoints) > 24) {
  final_state_24h <- setNames(as.numeric(simulated_values_24h_RNP_IN[nrow(simulated_values_24h_RNP_IN), -1]), names(simulated_values_24h_RNP_IN)[-1])
  remaining_timepoints <- seq(from = 24.01, to = max(Timepoints), by = 1/60)
  simulated_values_remaining_RNP <- simulate_ode(final_state_24h, list(
    k_dsb = mean_parameters_list$k_dsb,
    delay = mean_parameters_list$delay,
    k_pr = mean_parameters_list$k_pr_RNP,
    k_in = mean_parameters_list$k_in_RNP,
    k_ld = mean_parameters_list$k_ld_RNP,
    k_mh = if (has_MMEJ) mean_parameters_list$k_mh_RNP else NULL,
    k_ti = if (has_TI) mean_parameters_list$k_ti_RNP else NULL
  ), has_MMEJ, has_TI, remaining_timepoints)
  names(simulated_values_remaining_RNP) <- names(simulated_values_avg_df_RNP)
  simulated_values_avg_df_RNP_IN <- rbind(simulated_values_24h_RNP_IN, simulated_values_remaining_RNP)
} else {
  simulated_values_avg_df_RNP_IN <- simulated_values_24h_RNP_IN
}

adjust_wt_back <- function(simulated_values) {
  # Apply the adjustment logic to the WT column
  simulated_values$WT <- simulated_values$WT * ((100 - untransfected) / 100) + untransfected

  # Return the entire data frame, not just the WT column
  return(simulated_values)
}

# Reassign the entire data frame to ensure the WT column is updated
simulated_values_avg_df_RNP <- adjust_wt_back(simulated_values_avg_df_RNP)
simulated_values_avg_df_RNP_IN <- adjust_wt_back(simulated_values_avg_df_RNP_IN)

# Function to calculate correlations and residuals
calculate_metrics <- function(observed_data, simulated_data) {
  # Merge observed and simulated data by Timepoints
  merged_data <- merge(observed_data, simulated_data, by = "Timepoints")
  
  # Variables to calculate
  variables <- c("WT", "DSB", "Indel", "LD")
  if (has_MMEJ) variables <- c(variables, "MMEJ")
  if (has_TI) variables <- c(variables, "TI")
  
  # Calculate correlations and residuals
  correlations <- sapply(variables, function(var) cor(merged_data[[paste0(var, "_observed")]], merged_data[[var]], use = "complete.obs"))
  residuals <- sapply(variables, function(var) merged_data[[paste0(var, "_observed")]] - merged_data[[var]])
  
  list(correlations = correlations, residuals = as.data.frame(residuals))
}

# Calculate metrics for RNP and RNP_IN
metrics_RNP <- calculate_metrics(observed_data_RNP, simulated_values_avg_df_RNP)
metrics_RNP_IN <- calculate_metrics(observed_data_RNP_IN, simulated_values_avg_df_RNP_IN)

# Extract correlations and residuals
correlation_results_RNP <- metrics_RNP$correlations
correlation_results_RNP_IN <- metrics_RNP_IN$correlations
residuals_RNP <- metrics_RNP$residuals
residuals_RNP_IN <- metrics_RNP_IN$residuals

# Optionally, you can sum squared residuals as a performance measure
sum_squared_residuals_RNP <- sum(residuals_RNP^2)
sum_squared_residuals_RNP_IN <- sum(residuals_RNP_IN^2)


# Stop the cluster after completion
stopCluster(cl)

################################################### Plotting Results #########################################################

# Generate results table based on presence of MMEJ and/or TI data
if (has_MMEJ && has_TI) {
  results_combined <- data.frame(
    Mean = c(mean_parameters["k_dsb"], mean_parameters["delay"],
             mean_parameters["k_pr_RNP"], mean_parameters["k_in_RNP"], mean_parameters["k_ld_RNP"], mean_parameters["k_mh_RNP"], mean_parameters["k_ti_RNP"],
             mean_parameters["k_pr_RNP_IN"], mean_parameters["k_in_RNP_IN"], mean_parameters["k_ld_RNP_IN"], mean_parameters["k_mh_RNP_IN"], mean_parameters["k_ti_RNP_IN"]),
    SD = c(sd_parameters["k_dsb"], sd_parameters["delay"],
           sd_parameters["k_pr_RNP"], sd_parameters["k_in_RNP"], sd_parameters["k_ld_RNP"], sd_parameters["k_mh_RNP"], sd_parameters["k_ti_RNP"],
           sd_parameters["k_pr_RNP_IN"], sd_parameters["k_in_RNP_IN"], sd_parameters["k_ld_RNP_IN"], sd_parameters["k_mh_RNP_IN"], sd_parameters["k_ti_RNP_IN"]),
    CI_Lower = c(ci_parameters["CI Lower", "k_dsb"], ci_parameters["CI Lower", "delay"],
                 ci_parameters["CI Lower", "k_pr_RNP"], ci_parameters["CI Lower", "k_in_RNP"], ci_parameters["CI Lower", "k_ld_RNP"], ci_parameters["CI Lower", "k_mh_RNP"], ci_parameters["CI Lower", "k_ti_RNP"],
                 ci_parameters["CI Lower", "k_pr_RNP_IN"], ci_parameters["CI Lower", "k_in_RNP_IN"], ci_parameters["CI Lower", "k_ld_RNP_IN"], ci_parameters["CI Lower", "k_mh_RNP_IN"], ci_parameters["CI Lower", "k_ti_RNP_IN"]),
    CI_Upper = c(ci_parameters["CI Upper", "k_dsb"], ci_parameters["CI Upper", "delay"],
                 ci_parameters["CI Upper", "k_pr_RNP"], ci_parameters["CI Upper", "k_in_RNP"], ci_parameters["CI Upper", "k_ld_RNP"], ci_parameters["CI Upper", "k_mh_RNP"], ci_parameters["CI Upper", "k_ti_RNP"],
                 ci_parameters["CI Upper", "k_pr_RNP_IN"], ci_parameters["CI Upper", "k_in_RNP_IN"], ci_parameters["CI Upper", "k_ld_RNP_IN"], ci_parameters["CI Upper", "k_mh_RNP_IN"], ci_parameters["CI Upper", "k_ti_RNP_IN"])
  )
} else if (has_MMEJ) {
  results_combined <- data.frame(
    Mean = c(mean_parameters["k_dsb"], mean_parameters["delay"],
             mean_parameters["k_pr_RNP"], mean_parameters["k_in_RNP"], mean_parameters["k_ld_RNP"], mean_parameters["k_mh_RNP"],
             mean_parameters["k_pr_RNP_IN"], mean_parameters["k_in_RNP_IN"], mean_parameters["k_ld_RNP_IN"], mean_parameters["k_mh_RNP_IN"]),
    SD = c(sd_parameters["k_dsb"], sd_parameters["delay"],
           sd_parameters["k_pr_RNP"], sd_parameters["k_in_RNP"], sd_parameters["k_ld_RNP"], sd_parameters["k_mh_RNP"],
           sd_parameters["k_pr_RNP_IN"], sd_parameters["k_in_RNP_IN"], sd_parameters["k_ld_RNP_IN"], sd_parameters["k_mh_RNP_IN"]),
    CI_Lower = c(ci_parameters["CI Lower", "k_dsb"], ci_parameters["CI Lower", "delay"],
                 ci_parameters["CI Lower", "k_pr_RNP"], ci_parameters["CI Lower", "k_in_RNP"], ci_parameters["CI Lower", "k_ld_RNP"], ci_parameters["CI Lower", "k_mh_RNP"],
                 ci_parameters["CI Lower", "k_pr_RNP_IN"], ci_parameters["CI Lower", "k_in_RNP_IN"], ci_parameters["CI Lower", "k_ld_RNP_IN"], ci_parameters["CI Lower", "k_mh_RNP_IN"]),
    CI_Upper = c(ci_parameters["CI Upper", "k_dsb"], ci_parameters["CI Upper", "delay"],
                 ci_parameters["CI Upper", "k_pr_RNP"], ci_parameters["CI Upper", "k_in_RNP"], ci_parameters["CI Upper", "k_ld_RNP"], ci_parameters["CI Upper", "k_mh_RNP"],
                 ci_parameters["CI Upper", "k_pr_RNP_IN"], ci_parameters["CI Upper", "k_in_RNP_IN"], ci_parameters["CI Upper", "k_ld_RNP_IN"], ci_parameters["CI Upper", "k_mh_RNP_IN"])
  )
} else if (has_TI) {
  results_combined <- data.frame(
    Mean = c(mean_parameters["k_dsb"], mean_parameters["delay"],
             mean_parameters["k_pr_RNP"], mean_parameters["k_in_RNP"], mean_parameters["k_ld_RNP"], mean_parameters["k_ti_RNP"],
             mean_parameters["k_pr_RNP_IN"], mean_parameters["k_in_RNP_IN"], mean_parameters["k_ld_RNP_IN"], mean_parameters["k_ti_RNP_IN"]),
    SD = c(sd_parameters["k_dsb"], sd_parameters["delay"],
           sd_parameters["k_pr_RNP"], sd_parameters["k_in_RNP"], sd_parameters["k_ld_RNP"], sd_parameters["k_ti_RNP"],
           sd_parameters["k_pr_RNP_IN"], sd_parameters["k_in_RNP_IN"], sd_parameters["k_ld_RNP_IN"], sd_parameters["k_ti_RNP_IN"]),
    CI_Lower = c(ci_parameters["CI Lower", "k_dsb"], ci_parameters["CI Lower", "delay"],
                 ci_parameters["CI Lower", "k_pr_RNP"], ci_parameters["CI Lower", "k_in_RNP"], ci_parameters["CI Lower", "k_ld_RNP"], ci_parameters["CI Lower", "k_ti_RNP"],
                 ci_parameters["CI Lower", "k_pr_RNP_IN"], ci_parameters["CI Lower", "k_in_RNP_IN"], ci_parameters["CI Lower", "k_ld_RNP_IN"], ci_parameters["CI Lower", "k_ti_RNP_IN"]),
    CI_Upper = c(ci_parameters["CI Upper", "k_dsb"], ci_parameters["CI Upper", "delay"],
                 ci_parameters["CI Upper", "k_pr_RNP"], ci_parameters["CI Upper", "k_in_RNP"], ci_parameters["CI Upper", "k_ld_RNP"], ci_parameters["CI Upper", "k_ti_RNP"],
                 ci_parameters["CI Upper", "k_pr_RNP_IN"], ci_parameters["CI Upper", "k_in_RNP_IN"], ci_parameters["CI Upper", "k_ld_RNP_IN"], ci_parameters["CI Upper", "k_ti_RNP_IN"])
  )
} else {
  results_combined <- data.frame(
    Mean = c(mean_parameters["k_dsb"], mean_parameters["delay"],
             mean_parameters["k_pr_RNP"], mean_parameters["k_in_RNP"], mean_parameters["k_ld_RNP"],
             mean_parameters["k_pr_RNP_IN"], mean_parameters["k_in_RNP_IN"], mean_parameters["k_ld_RNP_IN"]),
    SD = c(sd_parameters["k_dsb"], sd_parameters["delay"],
           sd_parameters["k_pr_RNP"], sd_parameters["k_in_RNP"], sd_parameters["k_ld_RNP"],
           sd_parameters["k_pr_RNP_IN"], sd_parameters["k_in_RNP_IN"], sd_parameters["k_ld_RNP_IN"]),
    CI_Lower = c(ci_parameters["CI Lower", "k_dsb"], ci_parameters["CI Lower", "delay"],
                 ci_parameters["CI Lower", "k_pr_RNP"], ci_parameters["CI Lower", "k_in_RNP"], ci_parameters["CI Lower", "k_ld_RNP"],
                 ci_parameters["CI Lower", "k_pr_RNP_IN"], ci_parameters["CI Lower", "k_in_RNP_IN"], ci_parameters["CI Lower", "k_ld_RNP_IN"]),
    CI_Upper = c(ci_parameters["CI Upper", "k_dsb"], ci_parameters["CI Upper", "delay"],
                 ci_parameters["CI Upper", "k_pr_RNP"], ci_parameters["CI Upper", "k_in_RNP"], ci_parameters["CI Upper", "k_ld_RNP"],
                 ci_parameters["CI Upper", "k_pr_RNP_IN"], ci_parameters["CI Upper", "k_in_RNP_IN"], ci_parameters["CI Upper", "k_ld_RNP_IN"])
  )
}

# Define legend labels and colors
legend_labels <- c("WT + PR", "PR", "DSB", "Indel", "LD", "WT")
legend_colors <- c("black", "purple", "#80B1F9", "#988D4D", "#AFAFB1", "black")
if (has_MMEJ) {
  legend_labels <- c(legend_labels, "MMEJ")
  legend_colors <- c(legend_colors, "orange")
}
if (has_TI) {
  legend_labels <- c(legend_labels, "TI")
  legend_colors <- c(legend_colors, "green")
}

# Layout for Rstudio plots
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))

# Plot the model with smooth lines based on the selected parameters for the average data
plot_model <- function(simulated_values, title, log_scale = FALSE) {
  if (log_scale) {
    plot(simulated_values$Timepoints, simulated_values$WT, type = "n", xlab = "Hours", ylab = "Alleles (%)", xlim = c(0.1, max(Timepoints)), ylim = c(0, 100), main = title, log = 'x')
  } else {
    plot(simulated_values$Timepoints, simulated_values$WT, type = "n", xlab = "Hours", ylab = "Alleles (%)", xlim = c(0.1, max(Timepoints)), ylim = c(0, 100), main = title)
  }
  lines(simulated_values$Timepoints, simulated_values$WT, col = "black", lty = 1, lwd = 1.5)  # WT + PR
  lines(simulated_values$Timepoints, simulated_values$PR, col = "purple", lty = 1, lwd = 1.5) # PR
  lines(simulated_values$Timepoints, simulated_values$DSB, col = '#80B1F9', lty = 1, lwd = 1.5)  # DSB
  lines(simulated_values$Timepoints, simulated_values$Indel, col = '#988D4D', lty = 1, lwd = 1.5)  # Indel
  lines(simulated_values$Timepoints, simulated_values$LD, col = '#AFAFB1', lty = 1, lwd = 1.5)  # LD
  if (has_MMEJ) lines(simulated_values$Timepoints, simulated_values$MMEJ, col = 'orange', lty = 1, lwd = 1.5)  # MMEJ
  if (has_TI) lines(simulated_values$Timepoints, simulated_values$TI, col = 'green', lty = 1, lwd = 1.5)  # TI
  legend("topright", legend = legend_labels, col = legend_colors, lty = c(1, 1, 1, 1, 1, 2), lwd = 1.5)
}

# Plotting RNP model
plot_model(simulated_values_avg_df_RNP, "RNP ODE Model (Average)", log_scale = (max(Timepoints) > 24))
points(average_data_RNP$Timepoints, average_data_RNP$WT_observed, col = "black", pch = 16)
points(average_data_RNP$Timepoints, average_data_RNP$DSB_observed, col = '#80B1F9', pch = 16)
points(average_data_RNP$Timepoints, average_data_RNP$Indel_observed, col = "#988D4D", pch = 16)
points(average_data_RNP$Timepoints, average_data_RNP$LD_observed, col = "#AFAFB1", pch = 16)
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$WT - simulated_values_avg_df_RNP_IN$PR, col = "black", lty = 2, lwd = 1.5)  # WT
if (has_MMEJ) points(average_data_RNP$Timepoints, average_data_RNP$MMEJ_observed, col = "orange", pch = 16)
if (has_TI) points(average_data_RNP$Timepoints, average_data_RNP$TI_observed, col = "green", pch = 16)

# Calculate derivatives for each modeled line based on the mean parameters
calculate_derivatives <- function(simulated_values, all_timepoints) {
  derivatives <- as.matrix(simulated_values[, -1])
  diff(derivatives) / diff(all_timepoints) / 60
}

slopes_avg_RNP <- calculate_derivatives(simulated_values_avg_df_RNP, all_timepoints)

# Plot the derivatives with slopes for the average data
plot(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "WT"], type = "l", col = "black", lty = 1, lwd = 2, xlab = "Hours", ylab = "Change rate", xlim = c(0.1, max(Timepoints)), ylim = c(min(slopes_avg_RNP[, "WT"]), max(slopes_avg_RNP[, "DSB"])), main = "Event rates (Average)", log = if (max(Timepoints) > 24) 'x' else '')
lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "DSB"], col = '#80B1F9', lty = 1, lwd = 2)  # DSB
lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "Indel"], col = '#988D4D', lty = 1, lwd = 2)  # Indel
lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "LD"], col = '#AFAFB1', lty = 1, lwd = 2)  # LD
lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "PR"], col = 'purple', lty = 1, lwd = 2)  # PR
if (has_MMEJ) lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "MMEJ"], col = 'orange', lty = 1, lwd = 2)  # MMEJ
if (has_TI) lines(all_timepoints[-length(all_timepoints)], slopes_avg_RNP[, "TI"], col = 'green', lty = 1, lwd = 2)  # TI
abline(h = 0, col = "black")

# Plot cumulative DSBs and PR, and Intact DNA (WT without PR) for RNP
plot(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$C_DSB, type = "l", xlab = "Time", ylab = "Cumulative DSBs", col = "#80B1F9", xlim = c(0.1, max(Timepoints)), main = "Cumulative Products", log = if (max(Timepoints) > 24) 'x' else '')
lines(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$C_PR, col = "purple")  # Cumulative PR
lines(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$Indel, col = '#988D4D')  # Indel
lines(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$LD, col = '#AFAFB1')  # LD
if (has_MMEJ) lines(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$MMEJ, col = 'orange')  # MMEJ
if (has_TI) lines(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$TI, col = 'green')  # TI
legend("topleft", legend = legend_labels, col = legend_colors, lty = 1, cex = 0.8)

# Plotting RNP + IN model
plot_model(simulated_values_avg_df_RNP_IN, "RNP + IN ODE Model (Average)", log_scale = (max(Timepoints) > 24))
points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$WT_observed, col = "black", pch = 16)
points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$DSB_observed, col = '#80B1F9', pch = 16)
points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$Indel_observed, col = "#988D4D", pch = 16)
points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$LD_observed, col = "#AFAFB1", pch = 16)
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$WT - simulated_values_avg_df_RNP_IN$PR, col = "black", lty = 2, lwd = 1.5)  # WT
if (has_MMEJ) points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$MMEJ_observed, col = "orange", pch = 16)
if (has_TI) points(average_data_RNP_IN$Timepoints, average_data_RNP_IN$TI_observed, col = "green", pch = 16)
abline(v = 24, lty = 2, col = "red")

# Calculate derivatives for RNP + IN model
slopes_avg_RNP_IN <- calculate_derivatives(simulated_values_avg_df_RNP_IN, simulated_values_avg_df_RNP_IN$Timepoints)

# Plot the derivatives for the RNP + IN model
plot(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "WT"], type = "l", col = "black", lty = 1, lwd = 2, 
     xlab = "Hours", ylab = "Change rate", xlim = c(0.1, max(Timepoints)), 
     ylim = c(min(slopes_avg_RNP_IN[, "WT"]), max(slopes_avg_RNP_IN[, "DSB"])), main = "Event rates (RNP + IN)", log = if (max(Timepoints) > 24) 'x' else '')
lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "DSB"], col = '#80B1F9', lty = 1, lwd = 2)  # DSB
lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "Indel"], col = '#988D4D', lty = 1, lwd = 2)  # Indel
lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "LD"], col = '#AFAFB1', lty = 1, lwd = 2)  # LD
lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "PR"], col = 'purple', lty = 1, lwd = 2)  # PR
if (has_MMEJ) lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "MMEJ"], col = 'orange', lty = 1, lwd = 2)  # MMEJ
if (has_TI) lines(simulated_values_avg_df_RNP_IN$Timepoints[-1], slopes_avg_RNP_IN[, "TI"], col = 'green', lty = 1, lwd = 2)  # TI
abline(v = 24, lty = 2, col = "red")
abline(h = 0, col = "black")

# Plot cumulative DSBs and PR, and Intact DNA (WT without PR) for RNP + IN
plot(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$C_DSB, type = "l", xlab = "Time", 
     ylab = "Cumulative DSBs", col = "#80B1F9", xlim = c(0, max(Timepoints)), main = "Cumulative Products (RNP + IN)", 
     log = if (max(Timepoints) > 24) 'x' else '')
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$C_PR, col = "purple")  # Cumulative PR
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$Indel, col = '#988D4D')  # Indel
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$LD, col = '#AFAFB1')  # LD
if (has_MMEJ) lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$MMEJ, col = 'orange')  # MMEJ
if (has_TI) lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP_IN$TI, col = 'green')  # TI
legend("topleft", legend = legend_labels, col = legend_colors, lty = 1, cex = 0.8)

# Plot the difference between PR and WT simulated values for RNP and RNP + IN
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE)) 

plot(simulated_values_avg_df_RNP$Timepoints, simulated_values_avg_df_RNP$PR - simulated_values_avg_df_RNP_IN$PR, type = "l", 
     col = "purple", xlab = "Time", ylab = "Difference", ylim = c(0, max(100 - simulated_values_avg_df_RNP_IN$WT)), 
     main = "Difference in PR and WT Simulated Values")
lines(simulated_values_avg_df_RNP_IN$Timepoints, 100 - simulated_values_avg_df_RNP_IN$WT, col = "black", lty = 2, lwd = 1.5)  # WT + PR
lines(simulated_values_avg_df_RNP_IN$Timepoints, 100 - simulated_values_avg_df_RNP$WT, col = "black", lty = 2, lwd = 1.5)  # WT + PR
lines(simulated_values_avg_df_RNP_IN$Timepoints, simulated_values_avg_df_RNP$WT - simulated_values_avg_df_RNP_IN$WT, col = "black", lty = 2, lwd = 1.5)
abline(v = 24, lty = 2, col = "red")
legend("topright", legend = c("PR Difference", "WT Difference"), col = c("purple", "black"), lty = 1, cex = 0.8)

# Repair likelihood pie charts #
# For RNP Repair Likelihood - Directly from mean_parameters_list
k_dsb_RNP <- mean_parameters_list$k_dsb
k_pr_RNP <- mean_parameters_list$k_pr_RNP
k_in_RNP <- mean_parameters_list$k_in_RNP
k_ld_RNP <- mean_parameters_list$k_ld_RNP
k_mh_RNP <- if (has_MMEJ) mean_parameters_list$k_mh_RNP else 0
k_ti_RNP <- if (has_TI) mean_parameters_list$k_ti_RNP else 0

# Calculate total repair and unresolved DSBs
total_repair_RNP <- k_pr_RNP + k_in_RNP + k_ld_RNP + k_mh_RNP + k_ti_RNP
unresolved_RNP <- if (k_dsb_RNP > total_repair_RNP) k_dsb_RNP - total_repair_RNP else 0

# Create the data for the pie chart (Unresolved before MMEJ and TI)
repair_likelihood_RNP <- data.frame(
  category = c("Precise Repair", "Indels", "Large Deletions", "Unresolved", 
               if (has_MMEJ) "MMEJ" else NULL, if (has_TI) "TI" else NULL),
  value = c(k_pr_RNP, k_in_RNP, k_ld_RNP, unresolved_RNP, 
            if (has_MMEJ) k_mh_RNP else NULL, if (has_TI) k_ti_RNP else NULL)
)

# Calculate percentages
repair_likelihood_RNP$percentage <- round(100 * repair_likelihood_RNP$value / sum(repair_likelihood_RNP$value), 1)

# Plot the pie chart for RNP
pie(repair_likelihood_RNP$value, 
    labels = paste0(repair_likelihood_RNP$category, " (", repair_likelihood_RNP$percentage, "%)"), 
    col = c("purple", "#988D4D", "#AFAFB1", "#80B1F9", "green", "orange"), 
    main = "DSB Repair Outcomes (RNP)")

# For RNP + IN Repair Likelihood
k_pr_RNP_IN <- mean_parameters_list$k_pr_RNP_IN
k_in_RNP_IN <- mean_parameters_list$k_in_RNP_IN
k_ld_RNP_IN <- mean_parameters_list$k_ld_RNP_IN
k_mh_RNP_IN <- if (has_MMEJ) mean_parameters_list$k_mh_RNP_IN else 0
k_ti_RNP_IN <- if (has_TI) mean_parameters_list$k_ti_RNP_IN else 0

total_repair_RNP_IN <- k_pr_RNP_IN + k_in_RNP_IN + k_ld_RNP_IN + k_mh_RNP_IN + k_ti_RNP_IN
unresolved_RNP_IN <- if (mean_parameters_list$k_dsb > total_repair_RNP_IN) mean_parameters_list$k_dsb - total_repair_RNP_IN else 0

repair_likelihood_RNP_IN <- data.frame(
  category = c("Precise Repair", "Indels", "Large Deletions", "Unresolved", 
               if (has_MMEJ) "MMEJ" else NULL, if (has_TI) "TI" else NULL),
  value = c(k_pr_RNP_IN, k_in_RNP_IN, k_ld_RNP_IN, unresolved_RNP_IN, 
            if (has_MMEJ) k_mh_RNP_IN else NULL, if (has_TI) k_ti_RNP_IN else NULL)
)

repair_likelihood_RNP_IN$percentage <- round(100 * repair_likelihood_RNP_IN$value / sum(repair_likelihood_RNP_IN$value), 1)

# Plot the pie chart for RNP + IN
pie(repair_likelihood_RNP_IN$value, 
    labels = paste0(repair_likelihood_RNP_IN$category, " (", repair_likelihood_RNP_IN$percentage, "%)"), 
    col = c("purple", "#988D4D", "#AFAFB1", "#80B1F9", "green", "orange"), 
    main = "DSB Repair Outcomes (RNP + IN)")


# Constants
ln2 <- log(2)

# DSB Generation and Resolution Histograms
# For RNP (using ci_parameters for k_dsb and k_pr_RNP)
dsb_data_RNP <- data.frame(
  Event = c("DSB Generation", "DSB Resolution"),
  Mean = c(ln2 / mean_parameters_list$k_dsb, ln2 / (mean_parameters_list$k_pr_RNP + mean_parameters_list$k_in_RNP + mean_parameters_list$k_ld_RNP)),  # Mean from mean_parameters_list
  CI_Lower = c(ln2 / ci_parameters["CI Upper", "k_dsb"], ln2 / (ci_parameters["CI Upper", "k_pr_RNP"] + ci_parameters["CI Upper", "k_in_RNP"] + ci_parameters["CI Upper", "k_ld_RNP"])),  # Apply formula to CI
  CI_Upper = c(ln2 / ci_parameters["CI Lower", "k_dsb"], ln2 / (ci_parameters["CI Lower", "k_pr_RNP"] + ci_parameters["CI Lower", "k_in_RNP"] + ci_parameters["CI Lower", "k_ld_RNP"]))   # Apply formula to CI
)

# For RNP + IN (using ci_parameters for k_dsb and k_pr_RNP_IN)
dsb_data_RNP_IN <- data.frame(
  Event = c("DSB Generation", "DSB Resolution"),
  Mean = c(ln2 / mean_parameters_list$k_dsb, ln2 / (mean_parameters_list$k_pr_RNP_IN + mean_parameters_list$k_in_RNP_IN + mean_parameters_list$k_ld_RNP_IN)),  # Mean from mean_parameters_list
  CI_Lower = c(ln2 / ci_parameters["CI Upper", "k_dsb"], ln2 / (ci_parameters["CI Upper", "k_pr_RNP_IN"] + ci_parameters["CI Upper", "k_in_RNP_IN"] + ci_parameters["CI Upper", "k_ld_RNP_IN"])),  # Apply formula to CI
  CI_Upper = c(ln2 / ci_parameters["CI Lower", "k_dsb"], ln2 / (ci_parameters["CI Lower", "k_pr_RNP_IN"] + ci_parameters["CI Lower", "k_in_RNP_IN"] + ci_parameters["CI Lower", "k_ld_RNP_IN"]))   # Apply formula to CI
)

# Determine dynamic ylim based on the maximum of Means and CI_Upper
ylim_max_RNP <- max(dsb_data_RNP$Mean, dsb_data_RNP$CI_Upper) * 1.1  # Add 10% padding
ylim_max_RNP_IN <- max(dsb_data_RNP_IN$Mean, dsb_data_RNP_IN$CI_Upper) * 1.1  # Add 10% padding
ylim_max <- max(ylim_max_RNP, ylim_max_RNP_IN)

# Plot side-by-side histograms for RNP and RNP + IN
par(mfrow = c(1, 2))  # Set layout to side by side

# RNP Dot Plot
plot(1:2, dsb_data_RNP$Mean, pch = 16, col = "black", xaxt = "n", ylim = c(0, ylim_max), 
     main = "DSB Events (RNP)", ylab = "Coefficient", xlab = "Event", cex = 1.5)
arrows(1:2, dsb_data_RNP$CI_Lower, 1:2, dsb_data_RNP$CI_Upper, angle = 90, code = 3, length = 0.1)
axis(1, at = 1:2, labels = dsb_data_RNP$Event)

# RNP + IN Dot Plot
plot(1:2, dsb_data_RNP_IN$Mean, pch = 16, col = "black", xaxt = "n", ylim = c(0, ylim_max), 
     main = "DSB Events (RNP + IN)", ylab = "Coefficient", xlab = "Event", cex = 1.5)
arrows(1:2, dsb_data_RNP_IN$CI_Lower, 1:2, dsb_data_RNP_IN$CI_Upper, angle = 90, code = 3, length = 0.1)
axis(1, at = 1:2, labels = dsb_data_RNP_IN$Event)

# Reset layout back to default
par(mfrow = c(1, 1))


# Mark the end time
end_time <- Sys.time()

# Calculate the total execution time
execution_time <- end_time - start_time


################################################### Exporting Results #########################################################

# Check if export_results is TRUE before exporting
if (export_results) {
  
  # Generate a random alphanumeric run code (5 characters)
  run_code <- paste0(sample(0:9, 5, replace = TRUE), collapse = "")
  
  # Create a unique file path
  desktop_path <- paste0("~/Desktop/Cas9_Solver_Results_", run_code, ".xlsx")
  
  # Create a workbook
  wb <- createWorkbook()
  
  # Add a worksheet
  addWorksheet(wb, "Results")
  
  # Define a bold style
  bold_style <- createStyle(textDecoration = "bold")
  
  # Write the date, number of bootstraps, and results to the worksheet
  writeData(wb, "Results", paste("Date:", Sys.Date()), startRow = 1, startCol = 1)
  writeData(wb, "Results", paste("Number of bootstraps:", num_bootstrap_reps), startRow = 2, startCol = 6)
  writeData(wb, "Results", paste("Untransfected (%):", untransfected), startRow = 2, startCol = 3)
  writeData(wb, "Results", paste("Cas_Solver", cas_solver), startRow = 2, startCol = 1)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 1:2, cols = 1:6, gridExpand = TRUE)
  
  # Write the execution time
  writeData(wb, "Results", paste("Execution Time:", round(execution_time, 2)), startRow = 1, startCol = 3)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 3, cols = 1, gridExpand = TRUE)
  
  # Write the run code on line F1
  writeData(wb, "Results", paste("Run code:", run_code), startRow = 1, startCol = 6)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 1, cols = 6, gridExpand = TRUE)
  
  # Write the results_combined table
  writeData(wb, "Results", "Rate Coefficients", startRow = 4, startCol = 1)
  writeData(wb, "Results", results_combined, startRow = 5, startCol = 1, rowNames = TRUE)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 1, gridExpand = TRUE)
  
  # Write the correlation results for RNP with labels
  writeData(wb, "Results", "Correlation Results (RNP)", startRow = nrow(results_combined) + 7, startCol = 1)
  correlation_results_RNP_labeled <- data.frame(
    Alleles = c("WT_observed", "DSB_observed", "Indel_observed", "LD_observed", 
                if (has_MMEJ) "MMEJ_observed" else NULL, if (has_TI) "TI_observed" else NULL),
    Correlation = correlation_results_RNP
  )
  writeData(wb, "Results", correlation_results_RNP_labeled, startRow = nrow(results_combined) + 8, startCol = 1, rowNames = FALSE)
  addStyle(wb, sheet = "Results", style = bold_style, rows = nrow(results_combined) + 7, cols = 1, gridExpand = TRUE)
  
  # Write the correlation results for RNP + IN with labels
  writeData(wb, "Results", "Correlation Results (RNP + IN)", startRow = nrow(results_combined) + 7, startCol = 3)
  correlation_results_RNP_IN_labeled <- data.frame(
    Alleles = c("WT_observed", "DSB_observed", "Indel_observed", "LD_observed", 
                if (has_MMEJ) "MMEJ_observed" else NULL, if (has_TI) "TI_observed" else NULL),
    Correlation = correlation_results_RNP_IN
  )
  writeData(wb, "Results", correlation_results_RNP_IN_labeled, startRow = nrow(results_combined) + 8, startCol = 3, rowNames = FALSE)
  addStyle(wb, sheet = "Results", style = bold_style, rows = nrow(results_combined) + 7, cols = 3, gridExpand = TRUE)
  
  # Write the RNP Repair Likelihood data to the right of the results_combined table
  repair_likelihood_table_RNP <- data.frame(
    Category = c("Precise Repair", "Indels", "Large Deletions", "Unresolved", 
                 if (has_MMEJ) "MMEJ" else NULL, if (has_TI) "TI" else NULL),
    Value = c(repair_likelihood_RNP$percentage[1], repair_likelihood_RNP$percentage[2], repair_likelihood_RNP$percentage[3], repair_likelihood_RNP$percentage[4], 
              if (has_MMEJ) repair_likelihood_RNP$percentage[5] else NULL, if (has_TI) repair_likelihood_RNP$percentage[6] else NULL)
  )
  writeData(wb, "Results", "Repair Likelihood (RNP)", startRow = 4, startCol = 7)
  writeData(wb, "Results", repair_likelihood_table_RNP, startRow = 5, startCol = 7)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 7, gridExpand = TRUE)
  
  # Write the RNP + IN Repair Likelihood data to the right of the RNP table
  repair_likelihood_table_RNP_IN <- data.frame(
    Category = c("Precise Repair", "Indels", "Large Deletions", "Unresolved", 
                 if (has_MMEJ) "MMEJ" else NULL, if (has_TI) "TI" else NULL),
    Value = c(repair_likelihood_RNP_IN$percentage[1], repair_likelihood_RNP_IN$percentage[2], repair_likelihood_RNP_IN$percentage[3], repair_likelihood_RNP_IN$percentage[4], 
              if (has_MMEJ) repair_likelihood_RNP_IN$percentage[5] else NULL, if (has_TI) repair_likelihood_RNP_IN$percentage[6] else NULL)
  )
  writeData(wb, "Results", "Repair Likelihood (RNP + IN)", startRow = 4, startCol = 9)
  writeData(wb, "Results", repair_likelihood_table_RNP_IN, startRow = 5, startCol = 9)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 9, gridExpand = TRUE)
  
  # Write the DSB Generation/Resolution table for RNP to the right of the Repair Likelihood tables
  dsb_generation_table_RNP <- data.frame(
    Event = c("DSB Generation", "DSB Resolution"),
    Value = c(ln2 / mean_parameters_list$k_dsb, ln2 / (mean_parameters_list$k_pr_RNP + mean_parameters_list$k_in_RNP + mean_parameters_list$k_ld_RNP))  # Consistent with plotting
  )
  writeData(wb, "Results", "DSB Half-Life (RNP)", startRow = 4, startCol = 12)
  writeData(wb, "Results", dsb_generation_table_RNP, startRow = 5, startCol = 12)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 12, gridExpand = TRUE)
  
  # Write the DSB Generation/Resolution table for RNP + IN to the right of the RNP table
  dsb_generation_table_RNP_IN <- data.frame(
    Event = c("DSB Generation", "DSB Resolution"),
    Value = c(ln2 / mean_parameters_list$k_dsb, ln2 / (mean_parameters_list$k_pr_RNP_IN + mean_parameters_list$k_in_RNP_IN + mean_parameters_list$k_ld_RNP_IN))  # Consistent with plotting
  )
  writeData(wb, "Results", "DSB Half-Life (RNP + IN)", startRow = 4, startCol = 14)
  writeData(wb, "Results", dsb_generation_table_RNP_IN, startRow = 5, startCol = 14)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 14, gridExpand = TRUE)
  
  # Use the actual Timepoints from simulated_values_avg_df_RNP
  timepoints_RNP <- simulated_values_avg_df_RNP$Timepoints
  
  # Dynamically select the columns based on what's available in the simulated data
  available_columns_RNP <- c("WT", "DSB", "Indel", "LD", "PR", 
                             if (has_MMEJ) "MMEJ" else NULL, 
                             if (has_TI) "TI" else NULL,
                             "C_DSB", "C_PR")
  
  # Dynamically select columns from simulated_values_avg_df_RNP
  allele_freq_RNP <- simulated_values_avg_df_RNP[, available_columns_RNP, drop = FALSE]
  
  # Repeat for RNP + IN
  available_columns_RNP_IN <- c("WT", "DSB", "Indel", "LD", "PR",
                                if (has_MMEJ) "MMEJ" else NULL, 
                                if (has_TI) "TI" else NULL,
                                "C_DSB", "C_PR")
  allele_freq_RNP_IN <- simulated_values_avg_df_RNP_IN[, available_columns_RNP_IN, drop = FALSE]
  
  # Align the length of Timepoints and Slopes
  allele_freq_table_RNP <- data.frame(Hours = timepoints_RNP[-length(timepoints_RNP)], 
                                      allele_freq_RNP[-length(timepoints_RNP), ], 
                                      Slopes = slopes_avg_RNP)
  
  allele_freq_table_RNP_IN <- data.frame(Hours = timepoints_RNP[-length(timepoints_RNP)], 
                                         allele_freq_RNP_IN[-length(timepoints_RNP), ], 
                                         Slopes = slopes_avg_RNP_IN)
  
  # Add a column of labels
  labels <- names(params_data$init_params)
  
  # Add labels to column F (left side of the table) starting at row 16
  writeData(wb, "Results", labels, startRow = 5, startCol = 17, colNames = FALSE)
  
  # Add params_data table starting at G15 in the export file
  writeData(wb, "Results", "Initial Parameter Guesses", startRow = 4, startCol = 17)
  writeData(wb, "Results", params_data$init_params, startRow = 5, startCol = 18, rowNames = FALSE)
  
  writeData(wb, "Results", "Lower Bounds", startRow = 4, startCol = 19)
  writeData(wb, "Results", params_data$lower_bounds, startRow = 5, startCol = 19, rowNames = FALSE)
  
  writeData(wb, "Results", "Upper Bounds", startRow = 4, startCol = 20)
  writeData(wb, "Results", params_data$upper_bounds, startRow = 5, startCol = 20, rowNames = FALSE)
  
  # Apply bold style to headers
  addStyle(wb, sheet = "Results", style = bold_style, rows = 4, cols = 17:20, gridExpand = TRUE)
  
  # Write the allele frequencies below the existing tables
  writeData(wb, "Results", "Allele Frequencies (RNP)", startRow = 22, startCol = 1)
  writeData(wb, "Results", allele_freq_table_RNP, startRow = 23, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 22, cols = 1, gridExpand = TRUE)
  
  writeData(wb, "Results", "Allele Frequencies (RNP + IN)", startRow = 22, startCol = ncol(allele_freq_table_RNP) + 2)
  writeData(wb, "Results", allele_freq_table_RNP_IN, startRow = 23, startCol = ncol(allele_freq_table_RNP) + 2, colNames = TRUE)
  addStyle(wb, sheet = "Results", style = bold_style, rows = 22, cols = ncol(allele_freq_table_RNP) + 2, gridExpand = TRUE)
  
  # Save the workbook to an Excel file on your desktop
  saveWorkbook(wb, file = desktop_path, overwrite = TRUE)
  
  print(paste("Results successfully exported to", desktop_path, "with Run code:", run_code))
} else {
  print("Export cancelled.")
}

# Mark the end time and calculate the total execution time
print(paste("Total execution time:", execution_time))


################################################### Table of Results #########################################################

# Printing results
print(untransfected)

print(residuals_RNP)
print(paste("Sum of squared residuals for RNP:", sum_squared_residuals_RNP))
print(residuals_RNP_IN)
print(paste("Sum of squared residuals for RNP_IN:", sum_squared_residuals_RNP_IN))

print(correlation_results_RNP)
print(correlation_results_RNP_IN)
print(t(results_combined))

