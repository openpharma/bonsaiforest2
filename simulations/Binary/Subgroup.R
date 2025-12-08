# -----------------------------------------------------------------
# RUN NAIVE SUBGROUP ANALYSIS FOR A SINGLE ENDPOINT
#
# This script:
# 1. Sources all functions from `functions.R`.
# 2. Runs the naive subgroup analysis for the one endpoint
#    defined in the `ENDPOINT_ID` variable.
# 3. Saves the results to the `Results/` folder.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---
# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "binary"

# Set the main results directory
RESULTS_DIR <- "Results"

# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(checkmate)
library(dplyr)
library(broom)
library(MASS) # For glm.nb (used in naive)
library(parallel) # For compute_results
library(survival) # For Surv (used in naive)

# Source all helper functions (e.g., naive, generate_stacked_data, etc.)
source("functions.R")

# --- 2. DEFINE ENDPOINT PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

# This switch sets all the dynamic parameters for this endpoint
endpoint_params <- switch(ENDPOINT_ID,
                          "tte" = list(
                            folder = "TTE",
                            resp = "tt_pfs",
                            status = "ev_pfs",
                            resptype = "survival"
                          ),
                          "binary" = list(
                            folder = "Binary",
                            resp = "y",
                            status = NULL,
                            resptype = "binary"
                          ),
                          "count" = list(
                            folder = "Count",
                            resp = "y",
                            status = NULL,
                            resptype = "count"
                          ),
                          "continuous" = list(
                            folder = "Continuous",
                            resp = "y",
                            status = NULL,
                            resptype = "continuous"
                          ),
                          # Default error case
                          stop("Invalid ENDPOINT_ID. Must be one of: 'tte', 'binary', 'count', 'continuous'")
)


# --- 3. LOAD SCENARIO DATA ---
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

# --- 4. RUN ANALYSIS ---
message("--- Starting Naive Subgroup Analysis ---")

# 1. Create the analysis function using the constructor
#    'subgroup_method_endpoint' is the method we defined above
subgroup_analysis <- fun_analysis(subgroup_method_endpoint)

# 2. Define the cache file (the final results file)
#    This will be e.g. "Results/tte_subgroup.rds"
cache_file <- file.path(endpoint_params$folder,RESULTS_DIR, paste0(ENDPOINT_ID, "_subgroup.rds"))

# 3. Run the computation
#    This will apply `subgroup_analysis` to the `scenarios_list`
#    and save the result to `cache_file`.
subgroup_results <- compute_results(
  scenarios = scenarios_list,
  analyze = subgroup_analysis,
  cache = cache_file
)

message("--- Analysis Complete ---")
message(paste("Results saved to:", cache_file))
