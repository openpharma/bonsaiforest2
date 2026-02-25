#!/usr/bin/env Rscript
# Prepare Continuous simulation results - merge estimates with truth (PAPER folder)
# This script loads all Continuous result files and merges them with truth values
# All metric calculations are done in Paper_Results_Analysis.Rmd

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(rlang)

cat("=== Continuous Results Preparation (PAPER) - Merging Only ===\n\n")

# Set paths
cont_base_path <- "/home/pedreram/bonsaiforest2/simulations/paper/Continuous"
cont_results_path <- file.path(cont_base_path, "Results")
cont_output_file <- file.path(cont_base_path, "continuous_results_prepared.rds")

# Load truth data
cont_truth_file <- file.path(cont_base_path, "Scenarios", "truth.RData")
load(cont_truth_file)

cat("Loaded truth data from:", cont_truth_file, "\n")

# Extract subgroup-level truth (treatment effect)
cont_truth_subgroup <- simulation_truth_filtered %>%
  mutate(scenario_no = as.character(scenario_no)) %>%
  dplyr::select(scenario_no, scenario, subgroup_var, level, trt_effect) %>%
  rename(truth_trt_effect = trt_effect)

cat("Truth subgroups:", nrow(cont_truth_subgroup), "rows\n")
cat("Unique scenarios:", unique(cont_truth_subgroup$scenario_no), "\n")
cat("Scenarios with subgroups:\n")
print(cont_truth_subgroup %>% distinct(scenario_no, scenario))
cat("\n")

# Load all Continuous result files
cont_all_files <- list.files(cont_results_path, pattern = "\\.rds$", full.names = TRUE)

cat("Found", length(cont_all_files), "result files\n\n")

# Function to load and standardize Continuous results
load_and_standardize_continuous <- function(file_path) {
  
  df <- readRDS(file_path)
  file_name <- basename(file_path)
  
  # Create estimator name
  estimator_name <- str_remove(file_name, "^continuous_") %>%
                    str_remove("\\.rds$")
  
  # Normalize scenario column (some files use scenario_no, others use scenario_id)
  if (!"scenario_no" %in% names(df) && "scenario_id" %in% names(df)) {
    df <- df %>% mutate(scenario_no = scenario_id)
  }
  
  # Check if it's a naive estimator (naivepop/subgroup) or bonsaiforest2
  if (estimator_name %in% c("naivepop", "subgroup")) {
    # Naive estimators - old format
    if (estimator_name == "naivepop") {
      df_clean <- df %>%
        mutate(
          scenario_id = as.integer(scenario_no),
          replication_id = as.integer(sim_id),
          estimator = "population",
          estimate = trt_effect
        )
    } else {
      df_clean <- df %>%
        mutate(
          scenario_id = as.integer(scenario_no),
          replication_id = as.integer(sim_id),
          estimator = "subgroup",
          estimate = trt_effect,
          subgroup_var = subgroup_var,
          level = subgroup_level
        )
    }
  } else {
    # bonsaiforest2 estimators - new format with Median column
    # Handle different column naming conventions
    repl_col <- ifelse("replication_id" %in% names(df), "replication_id", "simulation_id")
    
    df_clean <- df %>%
      filter(Subgroup != "Overall") %>%
      mutate(
        scenario_id = as.integer(scenario_no),
        replication_id = as.integer(!!sym(repl_col)),
        estimator = paste(model_type, prior_name, sep = "_"),
        estimate = Median,  # Use Median column for estimate
        # Parse subgroup column (e.g., "X1: N" -> subgroup_var="X1", level="N")
        subgroup_var = str_extract(Subgroup, "^[^:]+"),
        level = str_trim(str_remove(Subgroup, "^[^:]+:\\s*"))
      )
  }
  
  df_clean %>%
    mutate(scenario_no = as.character(scenario_id)) %>%
    dplyr::select(any_of(c("scenario_id", "replication_id", "scenario_no", "estimator", 
                           "estimate", "subgroup_var", "level")))
}

# Load and process all files
cat("Processing result files:\n")
cont_all_results <- map_dfr(cont_all_files, function(fp) {
  file_name <- basename(fp)
  cat("  -", file_name, "\n")
  load_and_standardize_continuous(fp)
})

cat("\nLoaded", nrow(cont_all_results), "total result rows\n")
cat("Estimators found:", n_distinct(cont_all_results$estimator), "\n\n")

# Merge results with truth
# For population estimator: expand to all subgroups
# For subgroup estimator: match to corresponding subgroup

cont_results_merged <- cont_all_results %>%
  mutate(scenario_no = as.character(scenario_no)) %>%
  {
    pop_results <- filter(., is.na(subgroup_var) | estimator == "population")
    subgroup_results <- filter(., !is.na(subgroup_var) & estimator != "population")
    
    cat("Population estimator results:", nrow(pop_results), "rows\n")
    cat("Subgroup estimator results:", nrow(subgroup_results), "rows\n\n")
    
    # Population: expand to all subgroups
    if (nrow(pop_results) > 0) {
      pop_joined <- pop_results %>%
        dplyr::select(-subgroup_var, -level) %>%
        left_join(
          cont_truth_subgroup %>% dplyr::select(scenario_no, subgroup_var, level, truth_trt_effect),
          by = "scenario_no",
          relationship = "many-to-many"
        )
    } else {
      pop_joined <- tibble()
    }
    
    # Subgroup: match to truth
    if (nrow(subgroup_results) > 0) {
      subgroup_joined <- subgroup_results %>%
        left_join(cont_truth_subgroup,
                  by = c("scenario_no", "subgroup_var", "level"))
    } else {
      subgroup_joined <- tibble()
    }
    
    bind_rows(pop_joined, subgroup_joined)
  }

# Filter to only rows with truth values
cont_results_merged <- cont_results_merged %>%
  filter(!is.na(truth_trt_effect))

cat("Total estimate-truth pairs:", nrow(cont_results_merged), "\n\n")

cat("Final Continuous merged dataframe:\n")
cat("  Rows:", nrow(cont_results_merged), "\n")
cat("  Columns:", ncol(cont_results_merged), "\n")
cat("  Column names:", paste(names(cont_results_merged), collapse = ", "), "\n\n")

# Save the merged dataframe
saveRDS(cont_results_merged, cont_output_file)
cat("✓ Saved Continuous merged results to:", cont_output_file, "\n\n")

# Display sample
cat("Sample of Continuous results:\n")
print(head(cont_results_merged, 5))
