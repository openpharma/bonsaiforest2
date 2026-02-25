#--------------------------------------------------------------------------------------------------
#------ GENERATE AND STORE SIMULATION DATASETS FOR ALL CONTINUOUS SCENARIOS
#------ 
#------ This script:
#------ 1. Loads the calibrated scenario parameters
#------ 2. Generates n_simulations datasets for each scenario
#------ 3. Stores each scenario's datasets in a single RDS file for efficient loading
#--------------------------------------------------------------------------------------------------

library(benchtm)
library(dplyr)
library(tidyverse)

set.seed(0)

# ==================== CONFIGURATION ====================

# Number of simulations per scenario
n_simulations <- 1000  # 1000 simulations per scenario, 500 patients each (250 per arm)


# ==================== 1. LOAD SCENARIO PARAMETERS ====================

# Source the SimulationTruth.R to get calibrated selected_scen_par and X with categorical variables
message("Loading scenario parameters...")
source('SimulationTruth.R')

# Note: This will load:
# - selected_scen_par: calibrated scenario parameters
# - X: population data with categorical variables created
# - simulation_truth: saved as truth.RData

# ==================== 2. DEFINE SIMULATION FUNCTION ====================

#' Generate a single simulated dataset for a given scenario
#'
#' @param scenario_no The scenario number (1-12)
#' @param selected_scen_par Data frame with calibrated scenario parameters
#' @param n_patients Sample size (default 500)
#'
#' @return A data frame with columns: trt, X1, X2, X3, X4, X8, X11cat, X14cat, X17cat, Y
simulate_scenario <- function(scenario_no, selected_scen_par, n_patients = 500) {
  
  scen <- selected_scen_par[scenario_no, ]
  
  # Simulate covariates (via resampling from large population)
  df <- benchtm:::X_large_pop[sample(1:50000, n_patients), ]
  
  # Create categorical variables at specified cutpoints
  df$X11cat <- cut(df$X11, c(-Inf, 0.4615385, Inf), labels = c("a", "b"))
  df$X14cat <- cut(df$X14, c(-Inf, 0.2478632, 0.3333833, Inf), labels = c("a", "b", "c"))
  df$X17cat <- cut(df$X17, c(-Inf, 0.5888801, Inf), labels = c("a", "b"))
  
  # Simulate treatment assignment (balanced)
  trt <- rep(c(0, 1), n_patients / 2)
  
  # Simulate outcome Y
  df <- generate_y(
    X = df, 
    trt = trt,
    prog = scen$prog, 
    pred = scen$pred, 
    b0 = scen$b0, 
    b1 = scen$b1, 
    type = scen$type
  )
  
  # Select and return relevant variables
  df |> select(trt, X1, X2, X3, X4, X8, X11cat, X14cat, X17cat, Y)
}


# ==================== 3. GENERATE AND STORE DATASETS FOR EACH SCENARIO ====================

message("Generating simulation datasets for all scenarios...")

for (scenario_no in 1:nrow(selected_scen_par)) {
  
  # Get scenario info for messaging
  scen_info <- selected_scen_par[scenario_no, ]
  message(sprintf("  Scenario %d (%s): Generating %d simulations...", 
                  scenario_no, 
                  as.character(scen_info$scenario),
                  n_simulations))
  
  # Generate all simulations for this scenario as a list
  scenario_data <- replicate(
    n_simulations,
    simulate_scenario(scenario_no, selected_scen_par),
    simplify = FALSE
  )
  
  # Keep as list and name by sim_id (1, 2, 3, ..., 1000)
  names(scenario_data) <- as.character(1:n_simulations)
  
  # Add metadata to each data frame in the list
  scenario_data <- lapply(names(scenario_data), function(sim_id) {
    df <- scenario_data[[sim_id]]
    df %>%
      mutate(
        scenario_no = scenario_no,
        scenario = as.character(scen_info$scenario),
        sim_id = as.integer(sim_id),
        .before = trt
      )
  })
  
  # Name the list elements by sim_id
  names(scenario_data) <- as.character(1:n_simulations)
  
  # Save to RDS file
  output_file <- sprintf("Scenarios/scenario%d.rds", scenario_no)
  saveRDS(scenario_data, file = output_file)
  message(sprintf("    ✓ Saved to %s", output_file))
}

message("✓ All scenario datasets generated and stored successfully!")
message(sprintf("Generated 3 scenarios with %d simulations each", n_simulations))
