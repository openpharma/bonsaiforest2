#--------------------------------------------------------------------------------------------------
#------ SIMULATION TRUTH FOR TTE OUTCOMES (Scenarios 1, 2, 3, 4)
#------ 
#------ This script:
#------ 1. Uses the renumbered scenarios 1, 2, 3, 4 (originally 1, 2, 4, 5)
#------ 2. Computes empirical ground truth treatment effects from large simulation
#------ 3. Saves truth.RData with true effects for all subgroups
#--------------------------------------------------------------------------------------------------

library(survival)
library(tidyverse)
library(checkmate)
library(MASS)
library(broom)

set.seed(0)
RNGkind('Mersenne-Twister')

# ==================== CONFIGURATION ====================

# Load functions from parent directory
source('../functions.R')

# Scenarios to generate (1, 2, 3, 4 - renumbered from original 1, 2, 4, 5)
scenarios_to_use <- c(1, 2, 3, 4)

# Output directory
output_dir <- "Scenarios"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
}

# ==================== DEFINE SUBGROUPS ====================

# All 25 subgroups from 10 variables
subgroup_names <- c(
  "x_1.a", "x_1.b", "x_2.a", "x_2.b", "x_3.a", "x_3.b", "x_4.a",
  "x_4.b", "x_4.c", "x_5.a", "x_5.b", "x_5.c", "x_5.d", "x_6.a",
  "x_6.b", "x_7.a", "x_7.b", "x_8.a", "x_8.b", "x_8.c", "x_9.a",
  "x_9.b", "x_10.a", "x_10.b", "x_10.c"
)

# Formulas for truth calculation
base_model_tte <- Surv(tt_pfs, ev_pfs) ~ arm
subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10

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

# ==================== GENERATE TRUTH ====================

message("Generating empirical truth for TTE scenarios...")

truth_effects <- list()

for (scenario in scenarios_to_use) {
  
  message(paste("Processing Scenario", scenario, "..."))
  
  # Generate a large simulation dataset to compute truth
  # Using n_patients_truth = 10000 for better empirical estimation
  large_sim <- simul_scenario(
    scenario = as.character(scenario),
    n_datasets = 1,
    n_patients = 10000,
    n_events = 2470,  # Scaled proportionally for larger n
    sigma_aft = tte_params$sigma_aft,
    recr_duration = tte_params$recr_duration,
    rate_cens = tte_params$rate_cens,
    add_uncensored_pfs = tte_params$add_uncensored_pfs,
    inflation_factor = tte_params$inflation_factor
  )
  
  # Extract the single dataset
  truth_data <- large_sim[[1]]
  
  # Fit base model to get population effect
  pop_model <- coxph(base_model_tte, data = truth_data)
  pop_log_hr <- coef(pop_model)[["armTreated"]]
  
  message(paste("  Population log(HR):", round(pop_log_hr, 4)))
  
  # Compute subgroup-specific effects
  subgroup_effects <- rep(NA, length(subgroup_names))
  names(subgroup_effects) <- subgroup_names
  
  for (i in 1:length(subgroup_names)) {
    tryCatch({
      # Fit model for this specific subgroup
      sg_fit <- coxph(base_model_tte, data = truth_data, subset = truth_data$.subgroup_id == i)
      subgroup_effects[i] <- coef(sg_fit)[["armTreated"]]
    }, error = function(e) {
      warning(paste("Could not estimate effect for subgroup", subgroup_names[i]))
    })
  }
  
  truth_effects[[paste0("scenario", scenario)]] <- subgroup_effects
  
  message(paste("  Completed scenario", scenario))
}

# ==================== SAVE TRUTH ====================

message(paste("Saving truth to", file.path(output_dir, "truth.RData")))
save(truth_effects, file = file.path(output_dir, "truth.RData"))

message("\n✓ Truth generation complete!")
message(paste("Scenarios processed:", paste(scenarios_to_use, collapse = ", ")))
