# -----------------------------------------------------------------
# TEST OVAT APPROACH WITH A SINGLE DATASET
# Quick test to verify the model works before running full simulation
# -----------------------------------------------------------------

# --- 1. SETUP ---
library(bonsaiforest2)
library(brms)
library(survival)
library(dplyr)

setwd("/home/pedreram/bonsaiforest2/simulations")
source("functions.R")

# --- 2. LOAD ONE DATASET ---
cat("Loading one test dataset...\n")
scenario1 <- readRDS("TTE/Scenarios/scenario1.rds")
test_data <- scenario1[[1]]  # Just the first dataset from scenario 1

cat("Dataset loaded with", nrow(test_data), "rows\n")
cat("Variables:", names(test_data), "\n")

# --- 3. TEST WITH ONE COVARIATE ---
test_covariate <- "x_4"  # Test with just x_1

cat("\n--- Testing model with covariate:", test_covariate, "---\n")

# Prepare the model formula
prepared_run <- bonsaiforest2::prepare_formula_model(
  data = test_data,
  response_formula = Surv(tt_pfs, ev_pfs) ~ arm,
  response_type = "survival",
  unshrunk_terms_formula = paste("~", test_covariate, "+ arm"),
  shrunk_predictive_formula = paste("~ (0 + arm ||", test_covariate, ")")
)

cat("\nFormula created successfully!\n")
cat("Formula:\n")
print(prepared_run$formula)

# --- 4. FIT THE MODEL (quick test with few iterations) ---
cat("\n--- Fitting model (quick test with 500 iterations) ---\n")

fit_test <- fit_brms_model(
  prepared_model = prepared_run,
  sigma_ref = 1,
  chains = 2,        # Reduced for testing
  iter = 500,        # Reduced for testing
  warmup = 250,      # Reduced for testing
  cores = 1,
  backend = "cmdstanr",
  refresh = 50       # Show progress
)

cat("\n--- Model fitted successfully! ---\n")

# --- 5. SUMMARIZE RESULTS ---
cat("\n--- Getting subgroup effects summary ---\n")

summary_test <- bonsaiforest2::summary_subgroup_effects(fit_test
)

cat("\n--- Results ---\n")
print(summary_test$estimates)

cat("\n=== TEST COMPLETED SUCCESSFULLY ===\n")
cat("\nIf this worked, you can run the full simulation with OVAT_1.R\n")
