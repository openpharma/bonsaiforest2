#!/usr/bin/env Rscript

#--------------------------------------------------------------------------------------------------
#------ FILTER AND COMBINE RESULTS FOR SELECTED SCENARIOS (1, 2, 4, 5 -> 1, 2, 3, 4)
#------
#------ This script:
#------ 1. Reads results from the two scenario groupings (scenarios_1_3 and scenarios_4_6)
#------ 2. Filters to keep only scenarios 1, 2, 4, 5 from original numbering
#------ 3. Renumbers them as 1, 2, 3, 4
#------ 4. Combines into a single RDS file per model type
#--------------------------------------------------------------------------------------------------

library(tidyverse)
library(checkmate)

# Set working directory to Results folder
setwd("Results")

message("Starting result filtering and combining...")
message("Scenarios: Original (1, 2, 4, 5) -> Renumbered (1, 2, 3, 4)")

# ==================== GLOBAL MODEL RESULTS ====================

process_global_model <- function(base_name) {
  message(paste("\nProcessing:", base_name))
  
  # Read both scenario files
  file_1_3 <- paste0("tte_", base_name, "_scenarios_1_3.rds")
  file_4_6 <- paste0("tte_", base_name, "_scenarios_4_6.rds")
  
  df_1_3 <- readRDS(file_1_3)
  df_4_6 <- readRDS(file_4_6)
  
  # Filter to keep only scenarios 1, 2 from first file, 4, 5 from second file
  df_1_2 <- df_1_3 %>% filter(scenario_id %in% c(1, 2))
  df_4_5 <- df_4_6 %>% filter(scenario_id %in% c(4, 5))
  
  # Renumber scenarios 4, 5 to 3, 4
  df_4_5 <- df_4_5 %>%
    mutate(scenario_id = case_when(
      scenario_id == 4 ~ 3L,
      scenario_id == 5 ~ 4L,
      TRUE ~ scenario_id
    ))
  
  # Combine
  df_combined <- bind_rows(df_1_2, df_4_5)
  
  # Save as single file
  output_file <- paste0("tte_", base_name, ".rds")
  saveRDS(df_combined, file = output_file)
  
  message(paste("  Rows in original scenarios_1_3:", nrow(df_1_3)))
  message(paste("  Rows in original scenarios_4_6:", nrow(df_4_6)))
  message(paste("  Rows in filtered/combined output:", nrow(df_combined)))
  message(paste("  Saved to:", output_file))
}

# Process Global model variants
process_global_model("global_HN_global_phi_1")
process_global_model("global_HN_global_phi_delta_plan")
process_global_model("global_HN_global_phi_delta_plan_half")

# ==================== RHS MODEL RESULTS ====================

process_rhs_model <- function(base_name) {
  message(paste("\nProcessing:", base_name))
  
  # Read both scenario files
  file_1_3 <- paste0("tte_", base_name, "_scenarios_1_3.rds")
  file_4_6 <- paste0("tte_", base_name, "_scenarios_4_6.rds")
  
  df_1_3 <- readRDS(file_1_3)
  df_4_6 <- readRDS(file_4_6)
  
  # Filter to keep only scenarios 1, 2 from first file, 4, 5 from second file
  df_1_2 <- df_1_3 %>% filter(scenario_id %in% c(1, 2))
  df_4_5 <- df_4_6 %>% filter(scenario_id %in% c(4, 5))
  
  # Renumber scenarios 4, 5 to 3, 4
  df_4_5 <- df_4_5 %>%
    mutate(scenario_id = case_when(
      scenario_id == 4 ~ 3L,
      scenario_id == 5 ~ 4L,
      TRUE ~ scenario_id
    ))
  
  # Combine
  df_combined <- bind_rows(df_1_2, df_4_5)
  
  # Save as single file
  output_file <- paste0("tte_", base_name, ".rds")
  saveRDS(df_combined, file = output_file)
  
  message(paste("  Rows in original scenarios_1_3:", nrow(df_1_3)))
  message(paste("  Rows in original scenarios_4_6:", nrow(df_4_6)))
  message(paste("  Rows in filtered/combined output:", nrow(df_combined)))
  message(paste("  Saved to:", output_file))
}

# Process RHS model variants
process_rhs_model("RHS_theta0_1_s_2")
process_rhs_model("RHS_theta0_delta_plan_s_2delta_plan")
process_rhs_model("RHS_theta0_delta_plan_half_s_delta_plan")

# ==================== OVAT MODEL RESULTS ====================

process_ovat_models_by_prior <- function() {
  message("\nProcessing OVAT models by prior (combining all scenarios per prior)")
  
  # Map original scenario numbers to new numbering
  scenario_map <- list(
    "1" = 1L,
    "2" = 2L,
    "4" = 3L,
    "5" = 4L
  )
  
  # List of priors to process
  priors <- list(
    "HN_phi_1" = "1",
    "HN_phi_delta_plan" = "delta_plan",
    "HN_phi_delta_plan_half" = "delta_plan_half"
  )
  
  # For each prior, combine all scenarios
  for (prior_name in names(priors)) {
    message(paste("\nProcessing prior:", prior_name))
    
    prior_suffix <- priors[[prior_name]]
    
    # Read OVAT files for scenarios 1, 2, 4, 5 with this prior
    df_combined <- NULL
    
    for (orig_scenario in c(1, 2, 4, 5)) {
      file_path <- paste0("tte_ovat_", orig_scenario, "_oneway_HN_phi_", prior_suffix, ".rds")
      
      if (file.exists(file_path)) {
        df <- readRDS(file_path)
        
        # Renumber scenario
        new_scenario <- scenario_map[[as.character(orig_scenario)]]
        df <- df %>%
          mutate(scenario_id = new_scenario)
        
        df_combined <- bind_rows(df_combined, df)
        message(paste("  Added original scenario", orig_scenario, "->", new_scenario))
      } else {
        warning(paste("File not found:", file_path))
      }
    }
    
    # Save as single file for this prior
    output_file <- paste0("tte_ovat_oneway_", prior_name, ".rds")
    saveRDS(df_combined, file = output_file)
    
    message(paste("  Total rows for this prior:", nrow(df_combined)))
    message(paste("  Saved to:", output_file))
  }
}

process_ovat_models_by_prior()

# ==================== CLEAN UP ====================

message("\n" %+% "Removing old scenario-split files...")

# Remove old scenarios_1_3 and scenarios_4_6 files
old_files <- list.files(pattern = "scenarios_(1_3|4_6)\\.rds$")
file.remove(old_files)
message(paste("Removed", length(old_files), "old scenario-split files"))

# Remove old individual OVAT files organized by scenario and prior
# (they've been combined by prior instead)
old_ovat_files <- list.files(pattern = "ovat_[1245]_oneway_HN_phi_(1|delta_plan|delta_plan_half)\\.rds$")
file.remove(old_ovat_files)
message(paste("Removed", length(old_ovat_files), "old OVAT scenario/prior split files"))

# ==================== FINAL SUMMARY ====================

message("\n=== PROCESSING COMPLETE ===")
message("Result files structure:")
message("  Global Models (3 files, scenarios 1-4): 1 file per prior variant")
message("  RHS Models (3 files, scenarios 1-4): 1 file per prior variant")
message("  OVAT Models (3 files, scenarios 1-4): 1 file per prior variant")
message("  Baseline Models (2 files): population and subgroup")
message("\nNew files created:")
new_files <- list.files(pattern = "\\.rds$")
for (f in sort(new_files)) {
  size_mb <- file.size(f) / (1024^2)
  message(sprintf("  %-50s (%.1f MB)", f, size_mb))
}
message(paste("\nTotal files:", length(new_files)))
