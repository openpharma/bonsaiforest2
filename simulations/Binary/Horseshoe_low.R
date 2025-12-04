# -----------------------------------------------------------------
# RUN GLOBAL MODEL ANALYSIS (bonsaiforest2) FOR A SINGLE ENDPOINT
#
# This script:
# 1. Sources all functions from `../functions.R`.
# 2. Runs the "Global" model analysis using `bonsaiforest2` for
#    the one endpoint defined in `ENDPOINT_ID`.
# 3. Allows custom prior definitions.
# 4. Saves the results to the `Results/` folder.
#
# MODIFICATION: Removed the survival "monkey-patching" (knot-fixing)
# logic. This script WILL NOT work correctly for TTE analysis.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---

# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "binary"

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"

# DEFINE YOUR PRIORS HERE
PRIOR_SPECIFICATIONS <- list(
  horseshoe_low = list(prior = "horseshoe(scale_global = 0.025)", stanvars = NULL)
)
# ------------------------


# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(parallel)
library(parallelly)
library(dplyr)
library(tibble)
library(purrr)
library(checkmate)
library(cmdstanr)     # Dependency for bonsaiforest2
library(RhpcBLASctl)
library(stringr)
library(tidyr)
library(bonsaiforest2) # The main analysis package
library(survival)      # For Surv()

# Set threads to 1 for each parallel process
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

# Source all helper functions from one folder back
source("functions.R")

# --- 2. DEFINE ENDPOINT-SPECIFIC PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))
endpoint_params <- switch(ENDPOINT_ID,
                          "tte" = list(
                            folder = "TTE",
                            resp_formula = "Surv(tt_pfs, ev_pfs) ~ arm",
                            resptype = "survival"
                          ),
                          "binary" = list(
                            folder = "Binary",
                            resp_formula = "y ~ arm",
                            resptype = "binary"
                          ),
                          "count" = list(
                            folder = "Count",
                            resp_formula = "y ~ arm",
                            resptype = "count"
                          ),
                          "continuous" = list(
                            folder = "Continuous",
                            resp_formula = "y ~ arm",
                            resptype = "continuous"
                          ),
                          stop("Invalid ENDPOINT_ID. Must be one of: 'tte', 'binary', 'count', 'continuous'")
)

# --- 3. LOAD AND PROCESS INPUT DATA ---
cat(sprintf("--- Starting FULL GLOBAL MODEL RUN (%s) ---\n", ENDPOINT_ID))
cat("Loading all scenario data...\n")

scenarios_to_run <- as.character(1:6)
scenarios_list <- list()

for (scen in scenarios_to_run) {
  scen_file <- file.path(endpoint_params$folder, "Scenarios", paste0("scenario", scen, ".rds"))

  if (!file.exists(scen_file)) {
    stop(paste("File not found:", scen_file, "- Did you run the simulation generation script?"))
  }

  message(paste("   Loading:", scen_file))
  scenarios_list[[scen]] <- readRDS(scen_file)
}

# Process all datasets from all scenarios into a single flat list
all_datasets_processed <- lapply(seq_along(scenarios_list), function(i) {
  scenario_data <- scenarios_list[[i]]
  setNames(scenario_data, paste(i, 1:length(scenario_data), sep = "_"))
})

flat_named_list <- unlist(all_datasets_processed, recursive = FALSE)

cat(sprintf("Successfully loaded and processed %d total datasets from %d scenarios.\n",
            length(flat_named_list), length(scenarios_list)))

# Clean up memory
rm(scenarios_list, all_datasets_processed)
gc()


# --- 4. DEFINE PARAMETERS AND CREATE TASK LIST ---
RNGkind("L'Ecuyer-CMRG") # Recommended for parallel reproducibility
set.seed(42)

covariate_set <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")
num_cores <- availableCores()

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("Global"), # Only running "Global" model
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)

# --- 5. SET UP LOGGING AND RESULTS FILES ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "global", prior_names_str, sep = "_")

log_file <- file.path(endpoint_params$folder,RESULTS_DIR, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path(endpoint_params$folder,RESULTS_DIR, paste0(base_file_name, ".rds"))

total_tasks <- nrow(task_grid)
cat(sprintf("\nStarting %d total tasks using %d cores...\n", total_tasks, num_cores))
cat(sprintf("Results will be saved to: %s\n", results_file))
cat(sprintf("Progress will be logged to: %s\n\n", log_file))


# --- 6. RUN IN PARALLEL & COMBINE ---

run_single_task <- function(i) {
  current_task <- task_grid[i, ]

  run_global_model_task(
    sim_id = current_task$sim_id,
    model_type = current_task$model_type,
    prior_name = current_task$prior_name,
    all_data = flat_named_list,
    subgr_vars = covariate_set,
    log_path = log_file,
    prior_spec = PRIOR_SPECIFICATIONS[[current_task$prior_name]],
    endpoint_params = endpoint_params
  )
}

# --- Pre-compile step: Run task 1 serially ---
cat("--- Starting pre-compile step for task 1 ---\n")
result_task_1 <- list(run_single_task(1))
cat("--- Pre-compile task 1 finished. ---\n")

# --- Parallel step: Run tasks 2 to N ---
if (total_tasks > 1) {
  cat(sprintf("--- Starting parallel run for remaining %d tasks ---\n", total_tasks - 1))
  results_tasks_2_to_N <- mclapply(
    X = 2:total_tasks,
    FUN = run_single_task,
    mc.cores = num_cores,
    mc.preschedule = FALSE
  )
} else {
  cat("--- Only 1 task in total, parallel step skipped ---\n")
  results_tasks_2_to_N <- list()
}

# --- Combine pre-compile and parallel results ---
results_list <- c(result_task_1, results_tasks_2_to_N)
cat("\nAll tasks finished. Combining results...\n")

final_results_df <- bind_rows(results_list)

# --- 7. FINAL FORMATTING AND SAVING ---
if (nrow(final_results_df) > 0) {

  rows_per_task <- sapply(results_list, function(res) {
    if (is.data.frame(res)) nrow(res) else 1
  })

  task_indices <- rep(1:nrow(task_grid), times = rows_per_task)

  final_results_df <- final_results_df %>%
    mutate(
      sim_id_full = task_grid$sim_id[task_indices],
      scenario_id = as.integer(stringr::str_extract(sim_id_full, "^\\d+")),
      replication_id = as.integer(stringr::str_extract(sim_id_full, "\\d+$")),
      .before = 1
    ) %>%
    select(-sim_id_full)
}

# --- Save the final file ---
dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  object = final_results_df,
  file = results_file
)
cat(sprintf("Final adapted results saved to %s\n", results_file))
