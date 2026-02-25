# =========================================================================
# === NAIVE SUBGROUP ANALYSIS FOR CONTINUOUS OUTCOMES
# ===
# === This script:
# === 1. Loads all scenario datasets
# === 2. Performs subgroup analyses on each simulation
# === 3. Tests treatment effects within each subgroup
# === 4. Saves results to Results/ folder
# =========================================================================

# --- 0. CONFIGURATION ---
ENDPOINT_ID <- "continuous"
RESULTS_DIR <- "Results"

RNGkind('Mersenne-Twister')
set.seed(0)

# --- 1. LOAD LIBRARIES ---
message("--- Loading Libraries ---")
library(checkmate)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(parallel)

# --- 2. DEFINE ENDPOINT PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

# Determine the script's directory to make paths work from anywhere
# Get script directory using commandArgs (works with Rscript)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep('--file=', script_args, value = TRUE)
if (length(script_file) > 0) {
  script_dir <- dirname(sub('--file=', '', script_file))
} else {
  # Fallback: try sys.frame
  script_dir <- tryCatch({dirname(normalizePath(sys.frame(1)$ofile))}, error = function(e) {getwd()})
}
if (is.na(script_dir) || script_dir == "." || script_dir == "") {
  script_dir <- getwd()
}
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- getwd()
}

endpoint_params <- list(
  folder = script_dir,
  resp = "Y",
  resptype = "continuous"
)

# --- 3. SET UP RESULTS DIRECTORY ---
results_dir_path <- file.path(endpoint_params$folder, RESULTS_DIR)
dir.create(results_dir_path, recursive = TRUE, showWarnings = FALSE)

base_file_name <- paste(ENDPOINT_ID, "subgroup", sep = "_")
results_file <- file.path(results_dir_path, paste0(base_file_name, ".rds"))
log_file <- file.path(results_dir_path, paste0(base_file_name, ".log"))

if (file.exists(log_file)) file.remove(log_file)

cat(sprintf("Results will be saved to: %s\n", results_file))

# --- 4. LOAD SCENARIO DATA ---
message("Loading all scenario datasets...")
scenarios_dir <- file.path(endpoint_params$folder, "Scenarios")

all_scenarios <- list()
for (scen_no in 1:3) {
  scenario_file <- file.path(scenarios_dir, sprintf("scenario%d.rds", scen_no))
  if (file.exists(scenario_file)) {
    all_scenarios[[scen_no]] <- readRDS(scenario_file)
  } else {
    warning(sprintf("Scenario %d not found\n", scen_no))
  }
}

cat(sprintf("Loaded %d scenarios\n", length(all_scenarios)))

# --- 5. DEFINE SUBGROUP ANALYSIS FUNCTION ---

#' Perform naive subgroup analysis for a single simulation dataset
#'
#' @param df A data frame with columns: trt, subgroup variables, Y
#' @param subgroup_vars Names of subgrouping variables to test
#'
#' @return A data frame with subgroup analyses results

perform_subgroup_analysis <- function(df, subgroup_vars) {
  
  results_list <- list()
  
  # Test treatment effect within each subgroup variable
  for (subgroup_var in subgroup_vars) {
    
    # Get unique levels of this subgroup variable
    levels <- unique(df[[subgroup_var]])
    levels <- levels[!is.na(levels)]
    
    # Analyze treatment effect within each level
    var_results <- map_df(levels, function(level) {
      
      # Subset data to this subgroup level
      df_sub <- df %>% filter(!is.na(.data[[subgroup_var]]), .data[[subgroup_var]] == level)
      
      if (nrow(df_sub) < 2) {
        return(tibble(
          subgroup_var = subgroup_var,
          subgroup_level = as.character(level),
          n = nrow(df_sub),
          n_trt0 = 0,
          n_trt1 = 0,
          trt_effect = NA_real_,
          trt_se = NA_real_,
          trt_t_stat = NA_real_,
          trt_p_value = NA_real_,
          mean_y_trt0 = NA_real_,
          mean_y_trt1 = NA_real_,
          sd_y_trt0 = NA_real_,
          sd_y_trt1 = NA_real_
        ))
      }
      
      # Fit linear model: Y ~ trt
      fit <- lm(Y ~ trt, data = df_sub)
      tidy_res <- broom::tidy(fit) %>% filter(term == "trt")
      
      # Get counts and means by treatment group
      summary_stats <- df_sub %>%
        group_by(trt) %>%
        summarise(
          n_trt = n(),
          mean_y = mean(Y),
          sd_y = sd(Y),
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = trt,
          values_from = c(n_trt, mean_y, sd_y),
          names_glue = "{.value}_trt{trt}"
        )
      
      # Combine results
      tibble(
        subgroup_var = subgroup_var,
        subgroup_level = as.character(level),
        n = nrow(df_sub),
        n_trt0 = summary_stats$n_trt_trt0,
        n_trt1 = summary_stats$n_trt_trt1,
        trt_effect = tidy_res$estimate[1],
        trt_se = tidy_res$std.error[1],
        trt_t_stat = tidy_res$statistic[1],
        trt_p_value = tidy_res$p.value[1],
        mean_y_trt0 = summary_stats$mean_y_trt0,
        mean_y_trt1 = summary_stats$mean_y_trt1,
        sd_y_trt0 = summary_stats$sd_y_trt0,
        sd_y_trt1 = summary_stats$sd_y_trt1
      )
    })
    
    results_list[[subgroup_var]] <- var_results
  }
  
  bind_rows(results_list)
}

# --- 6. RUN ANALYSIS FOR ALL SCENARIOS ---

message("--- Starting Naive Subgroup Analysis ---")

# Subgrouping variables to analyze
subgroup_vars <- c("X1", "X2", "X3", "X4", "X8", "X11cat", "X14cat", "X17cat")

all_results <- list()

for (scen_no in seq_along(all_scenarios)) {
  
  # Combine the list of data.frames into a single stacked data.frame
  scenario_data <- bind_rows(all_scenarios[[scen_no]])
  
  cat(sprintf("\n--- Scenario %d (%s) ---\n", 
              scen_no, unique(scenario_data$scenario)[1]))
  
  start_t <- Sys.time()
  
  # Perform subgroup analysis for each simulation
  scenario_results <- scenario_data %>%
    group_by(sim_id) %>%
    nest() %>%
    mutate(
      subgroup_results = map(data, ~perform_subgroup_analysis(.x, subgroup_vars))
    ) %>%
    unnest(subgroup_results) %>%
    select(-data) %>%
    mutate(
      scenario_no = scen_no,
      scenario = unique(scenario_data$scenario)[1],
      .before = sim_id
    )
  
  end_t <- Sys.time()
  elapsed <- difftime(end_t, start_t, units = "secs")
  
  cat(sprintf("Completed in %0.1f seconds\n", elapsed))
  cat(sprintf("Generated %d subgroup estimates\n", nrow(scenario_results)))
  
  all_results[[scen_no]] <- scenario_results
}

# --- 7. COMBINE AND SAVE ---

combined_results <- bind_rows(all_results)

saveRDS(combined_results, file = results_file)
cat(sprintf("\n✓ All results saved to %s\n", results_file))
cat(sprintf("✓ Processed %d scenarios\n", n_distinct(combined_results$scenario_no)))
cat(sprintf("✓ Generated %d total subgroup estimates\n", nrow(combined_results)))

# Summary statistics
cat("\n--- Summary ---\n")
cat(sprintf("Subgroup variables tested: %s\n", paste(subgroup_vars, collapse = ", ")))
cat(sprintf("Scenarios: %d (1-12)\n", n_distinct(combined_results$scenario_no)))
cat(sprintf("Simulations per scenario: %d\n", n_distinct(combined_results$sim_id)))
cat(sprintf("Total subgroup estimates: %d\n", nrow(combined_results)))
