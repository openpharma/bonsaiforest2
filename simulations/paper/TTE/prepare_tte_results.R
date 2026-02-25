#!/usr/bin/env Rscript
# Prepare TTE simulation results - merge estimates with truth (PAPER folder)
# This script loads all TTE result files and merges them with truth values
# All metric calculations are done in Paper_Results_Analysis.Rmd

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

cat("=== TTE Results Preparation (PAPER) - Merging Only ===\n\n")

# Set paths for paper simulations
tte_base_path <- "/home/pedreram/bonsaiforest2/simulations/paper/TTE"
tte_results_path <- file.path(tte_base_path, "Results")
tte_output_file <- file.path(tte_base_path, "tte_results_prepared.rds")

# Load truth data
tte_truth_file <- file.path(tte_base_path, "Scenarios", "truth.RData")
load(tte_truth_file)

cat("Loaded truth data from:", tte_truth_file, "\n")

# Extract subgroup-level truth (log hazard ratios)
tte_truth_subgroup <- simul_parameter$true_subgroup_ahr %>%
  mutate(scenario_no = as.character(row_number())) %>%
  pivot_longer(
    cols = -scenario_no,
    names_to = "subgroup_raw",
    values_to = "truth_ahr"
  ) %>%
  mutate(
    join_key = gsub("x_|\\.", "", subgroup_raw) %>% tolower(),
    truth_log = log(truth_ahr)
  ) %>%
  dplyr::select(scenario_no, join_key, truth_log, truth_ahr)

# Extract population-level truth
tte_truth_population <- simul_parameter$true_overall_results %>%
  mutate(
    scenario_no = as.character(1:6),
    truth_pop_ahr = AHR,
    truth_pop_log = log(AHR)
  ) %>%
  dplyr::select(scenario_no, truth_pop_ahr, truth_pop_log)

cat("Truth subgroup rows:", nrow(tte_truth_subgroup), "\n")
cat("Unique scenarios:", unique(tte_truth_subgroup$scenario_no), "\n\n")

# Load all TTE result files
tte_all_files <- list.files(tte_results_path, pattern = "\\.rds$", full.names = TRUE)

cat("Found", length(tte_all_files), "result files\n\n")

# Function to load and standardize TTE results
load_and_standardize_tte <- function(file_path) {
  
  df <- readRDS(file_path)
  estimator_name <- str_remove(basename(file_path), "^tte_") %>%
                    str_remove("\\.rds$")
  
  # Check if it's a naive estimator (population/subgroup) or bonsaiforest2
  if (estimator_name %in% c("population", "subgroup")) {
    # Naive estimators
    df_clean <- df %>%
      mutate(
        scenario_id = as.integer(scenario_no),
        replication_id = as.integer(simul_no),
        join_key = gsub("S_|Overall", "", subgroup) %>% 
                   str_replace_all("[\\._\\s]", "") %>% 
                   tolower(),
        estimator = estimator_name,
        estimate_log = estimate,
        ci_lower = lower_ci,
        ci_upper = upper_ci
      ) %>%
      filter(join_key != "")
      
  } else {
    # bonsaiforest2 estimators (Global/OVAT)
    df_clean <- df %>%
      filter(Subgroup != "Overall") %>%
      mutate(
        scenario_id = as.integer(scenario_id),
        replication_id = as.integer(replication_id),
        join_key = gsub("x_|: | ", "", Subgroup) %>% 
                   str_replace_all("[\\._\\s]", "") %>% 
                   tolower(),
        estimator = paste(model_type, prior_name, sep = "_"),
        estimate_log = log(Median),
        ci_lower = log(CI_Lower),
        ci_upper = log(CI_Upper)
      )
  }
  
  df_clean %>%
    dplyr::select(scenario_id, replication_id, estimator, join_key, 
           estimate_log, ci_lower, ci_upper) %>%
    mutate(scenario_no = as.character(scenario_id))
}

# Load and process all files
cat("Processing result files:\n")
tte_all_results <- map_dfr(tte_all_files, function(fp) {
  estimator_name <- str_remove(basename(fp), "^tte_") %>% str_remove("\\.rds$")
  cat("  -", estimator_name, "\n")
  load_and_standardize_tte(fp)
})

cat("\nLoaded", nrow(tte_all_results), "total result rows\n")
cat("Estimators found:", n_distinct(tte_all_results$estimator), "\n\n")

# Merge results with truth
tte_results_merged <- tte_all_results %>%
  left_join(tte_truth_subgroup, by = c("scenario_no", "join_key")) %>%
  left_join(tte_truth_population, by = "scenario_no")

# Filter to only rows with truth values and keep scenarios 1, 2, 4, 5
tte_results_merged <- tte_results_merged %>%
  filter(!is.na(truth_log)) %>%
  filter(scenario_no %in% c("1", "2", "4", "5"))

cat("Total estimates with truth:", nrow(tte_results_merged), "\n\n")

cat("Final TTE merged dataframe:\n")
cat("  Rows:", nrow(tte_results_merged), "\n")
cat("  Columns:", ncol(tte_results_merged), "\n")
cat("  Column names:", paste(names(tte_results_merged), collapse = ", "), "\n\n")

# Save the merged dataframe
saveRDS(tte_results_merged, tte_output_file)
cat("✓ Saved TTE merged results to:", tte_output_file, "\n\n")

# Display sample
cat("Sample of TTE results:\n")
print(head(tte_results_merged, 5))
