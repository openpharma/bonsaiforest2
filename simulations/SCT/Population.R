# -----------------------------------------------------------------
# RUN NAIVE POPULATION ANALYSIS FOR SCT TRIAL
#
# This script runs the naive population (overall) analysis for the
# SCT continuous endpoint (MFM32 change from baseline).
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---
ENDPOINT_ID <- "continuous"
RESULTS_DIR <- "Results"

# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(checkmate)
library(dplyr)
library(broom)
library(MASS)
library(parallel)
library(survival)

# Source helper functions from parent directory  
source('../functions.R')

# Check if required functions are loaded
if (!exists("naivepop") || !exists("fun_analysis") || !exists("compute_results")) {
  stop("Error: Required functions not found in ../functions.R")
}

# --- 2. DEFINE ENDPOINT PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

endpoint_params <- list(
  folder = "SCT",
  resp = "y",
  status = NULL,
  resptype = "continuous"
)

# --- 3. LOAD SCENARIO DATA ---
scenarios_to_run <- as.character(1:6)  # All 6 scenarios for SCT
scenarios_list <- list()

message(paste("Loading scenarios from folder:", endpoint_params$folder))

for (scen in scenarios_to_run) {
  scen_file <- file.path("Scenarios", paste0("SCT_Scenario_", scen, ".rds"))
  
  if (!file.exists(scen_file)) {
    stop(paste("File not found:", scen_file, "- Did you run Scenarios_generation.Rmd?"))
  }
  
  message(paste("   Loading:", scen_file))
  scenarios_list[[scen]] <- readRDS(scen_file)
}

# --- 4. RUN ANALYSIS ---
message("--- Starting Naive Population Analysis ---")

# 1. Create the analysis function using the constructor
population_analysis <- fun_analysis(population_method_endpoint)

# 2. Define the cache file (the final results file)
cache_file <- file.path(RESULTS_DIR, paste0("SCT_population.rds"))

# 3. Run the computation
population_results <- compute_results(
  scenarios = scenarios_list,
  analyze = population_analysis,
  cache = cache_file
)

# --- 5. DONE ---
message("--- Population Analysis Complete ---")
message(paste("Results saved to:", cache_file))
