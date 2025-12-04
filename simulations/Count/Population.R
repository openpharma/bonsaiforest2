# -----------------------------------------------------------------
# RUN NAIVE POPULATION ANALYSIS FOR A SINGLE ENDPOINT
#
# This script:
# 1. Sources all functions from `functions.R` (including `naivepop`).
# 2. Runs the naive population (overall) analysis for the one
#    endpoint defined in the `ENDPOINT_ID` variable.
# 3. Saves the results to the `Results/` folder using the
#    standardized (log-scale) format.
#
# MODIFIED: `population_method_endpoint` is now fully dynamic
#           and no longer hard-coded for survival.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---
# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "count"

# Set the main results directory
RESULTS_DIR <- "Results"
# ------------------------


# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(checkmate)
library(dplyr)
library(broom)
library(MASS) # For glm.nb
library(parallel) # For compute_results
library(survival) # For Surv

# Source all helper functions (from the same directory)
# Make sure this file contains the `naivepop` function you provided
source('functions.R')

# Check if required functions are loaded
if (!exists("naivepop") || !exists("fun_analysis") || !exists("compute_results")) {
  stop("Error: `naivepop`, `fun_analysis`, or `compute_results` not found.
       Please ensure they are in `functions.R`.")
}


# --- 2. DEFINE ENDPOINT PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

# This switch sets all the dynamic parameters for this endpoint
endpoint_params <- switch(ENDPOINT_ID,
                          "tte" = list(
                            folder = "TTE",
                            resp = "tt_pfs", status = "ev_pfs", resptype = "survival"
                          ),
                          "binary" = list(
                            folder = "Binary", # <-- Corrected capitalization
                            resp = "y", status = NULL, resptype = "binary"
                          ),
                          "count" = list(
                            folder = "Count",
                            resp = "y", status = NULL, resptype = "count"
                          ),
                          "continuous" = list(
                            folder = "Continuous", # <-- Corrected capitalization
                            resp = "y", status = NULL, resptype = "continuous"
                          ),
                          stop("Invalid ENDPOINT_ID. Must be one of: 'tte', 'binary', 'count', 'continuous'")
)

# --- 3. DEFINE THE CORRECTED ANALYSIS METHOD ---

# --- 4. LOAD SCENARIO DATA ---
scenarios_to_run <- as.character(1:6)
scenarios_list <- list()

message(paste("Loading 6 scenarios from folder:", endpoint_params$folder))

for (scen in scenarios_to_run) {
  scen_file <- file.path(endpoint_params$folder, "Scenarios", paste0("scenario", scen, ".rds"))

  if (!file.exists(scen_file)) {
    stop(paste("File not found:", scen_file, "- Did you run the simulation generation script?"))
  }

  message(paste("   Loading:", scen_file))
  scenarios_list[[scen]] <- readRDS(scen_file)
}


# --- 5. RUN ANALYSIS ---
message("--- Starting Naive Population Analysis ---")

# 1. Create the analysis function using the constructor
population_analysis <- fun_analysis(population_method_endpoint)

# 2. Define the cache file (the final results file)
#    This will be e.g. "Results/binary_population.rds"
cache_file <- file.path(endpoint_params$folder,RESULTS_DIR, paste0(ENDPOINT_ID, "_population.rds"))

# 3. Run the computation
population_results <- compute_results(
  scenarios = scenarios_list,
  analyze = population_analysis,
  cache = cache_file
)

message("--- Analysis Complete ---")
message(paste("Results saved to:", cache_file))
