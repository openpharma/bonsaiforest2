#--------------------------------------------------------------------------------------------------
#------ GENERATE AND STORE SIMULATION DATASETS FOR TTE SCENARIOS (1, 2, 3, 4)
#------
#------ This script:
#------ 1. Generates n_simulations datasets for each selected scenario
#------ 2. Stores each scenario's datasets in a single RDS file for efficient loading
#------ 3. Uses the Wolbers et al. (2025) framework with renumbered scenarios
#--------------------------------------------------------------------------------------------------

library(checkmate)
library(tidyverse)

set.seed(0)
RNGkind('Mersenne-Twister')

# ==================== CONFIGURATION ====================

# Number of simulations per scenario
n_simulations <- 1000

# Scenarios to generate (1, 2, 3, 4 - renumbered from original 1, 2, 4, 5)
scenarios_to_use <- c(1, 2, 3, 4)

# Load functions from parent directory
source('../functions.R')

# ==================== TTE PARAMETERS ====================

tte_params <- list(
  n_patients = 1000,
  n_events = 247,
  sigma_aft = 0.85,
  recr_duration = 3,
  rate_cens = 0.02,
  add_uncensored_pfs = FALSE,
  inflation_factor = 1
)

# ==================== CREATE OUTPUT DIRECTORY ====================

output_dir <- "Scenarios"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
}

# ==================== GENERATE SCENARIOS ====================

message("Starting TTE scenario generation (Scenarios 1, 2, 3, 4)...")
message(paste("Number of simulations per scenario:", n_simulations))

for (scenario in scenarios_to_use) {

  start_time <- Sys.time()
  message(paste("\nGenerating Scenario", scenario, "..."))

  # Generate n_simulations datasets for this scenario
  scenario_data <- simul_scenario(
    scenario = as.character(scenario),
    n_datasets = n_simulations,
    n_patients = tte_params$n_patients,
    n_events = tte_params$n_events,
    sigma_aft = tte_params$sigma_aft,
    recr_duration = tte_params$recr_duration,
    rate_cens = tte_params$rate_cens,
    add_uncensored_pfs = tte_params$add_uncensored_pfs,
    inflation_factor = tte_params$inflation_factor
  )

  # Save to RDS file
  output_file <- file.path(output_dir, paste0("scenario", scenario, ".rds"))
  saveRDS(scenario_data, file = output_file)

  elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
  file_size_mb <- file.size(output_file) / (1024^2)

  message(paste("    Scenario", scenario, "complete"))
  message(paste("    File:", basename(output_file)))
  message(paste("    Size:", round(file_size_mb, 1), "MB"))
  message(paste("    Time:", round(elapsed_time, 1), "minutes"))
}

message("\nAll TTE scenarios generated successfully!")
message(paste("Scenarios:", paste(scenarios_to_use, collapse = ", ")))
message(paste("Location:", file.path(getwd(), output_dir)))
