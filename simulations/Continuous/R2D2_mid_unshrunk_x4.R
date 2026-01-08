# -----------------------------------------------------------------
# RUN GLOBAL MODEL ANALYSIS (bonsaiforest2) FOR CONTINUOUS ENDPOINT
# WITH UNSHRUNK PREDICTIVE VARIABLE x_4
#
# This script:
# 1. Sources all functions from `../functions.R`.
# 2. Runs the "Global" model analysis using `bonsaiforest2` with
#    R2D2 mid shrinkage applied to all predictive variables EXCEPT x_4.
# 3. Variable x_4 is specified as an unshrunk predictive variable
#    with a normal prior.
# 4. Saves the results to the `Results/` folder.
#
# Note: Uses unique cmdstanr directories per prior to prevent
# compilation conflicts when running parallel simulations.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---

# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "continuous"

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"

# Set unique seed
RNGkind('Mersenne-Twister')
set.seed(0)

# DEFINE YOUR PRIORS HERE
# This configuration tests the new functionality of unshrunk predictive variables
# x_4 is excluded from shrinkage and given a normal prior
PRIOR_SPECIFICATIONS <- list(
  r2d2_mid_unshrunk_x4 = list(
    prior = "R2D2(mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 0.3)", 
    unshrunk_prior = "normal(0, 5)",
    stanvars = NULL
  )
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

# Source all helper functions from simulations folder
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
# All covariates - x_4 will be split into unshrunk predictive
covariate_set <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")

# Variables that will receive shrinkage (all except x_4 for predictive)
shrunk_vars <- c("x_1", "x_2", "x_3", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")

# x_4 will be unshrunk for predictive effect
unshrunk_var <- "x_4"

num_cores <- 96

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("Global"), # Only running "Global" model
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)

# --- 5. SET UP LOGGING AND RESULTS FILES ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "global", prior_names_str, sep = "_")

log_file <- file.path(endpoint_params$folder, RESULTS_DIR, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path(endpoint_params$folder, RESULTS_DIR, paste0(base_file_name, ".rds"))

total_tasks <- nrow(task_grid)
cat(sprintf("\nStarting %d total tasks using %d cores...\n", total_tasks, num_cores))
cat(sprintf("Results will be saved to: %s\n", results_file))
cat(sprintf("Progress will be logged to: %s\n\n", log_file))


# --- 6. DEFINE CUSTOM TASK FUNCTION WITH UNSHRUNK x_4 ---

run_single_task_unshrunk <- function(i) {
  current_task <- task_grid[i, ]
  
  # --- 1. Thread Safety ---
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(threads = 1)

  # --- 2. Logging Setup ---
  task_description <- sprintf("SimID: %s | Model: %s | Prior: %s | Endpoint: %s",
                              current_task$sim_id, current_task$model_type, 
                              current_task$prior_name, endpoint_params$resptype)
  start_time <- Sys.time()

  cat(sprintf("[%s] STARTING: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S"), task_description),
      file = log_file, append = TRUE)

  # --- 3. Execute the Core Task ---
  results <- tryCatch({
    df <- flat_named_list[[current_task$sim_id]]
    prior_spec <- PRIOR_SPECIFICATIONS[[current_task$prior_name]]

    # --- Model Fitting Logic with x_4 as unshrunk predictive ---
    # All variables as prognostic (unshrunk)
    prognostic_str <- paste(covariate_set, collapse = " + ")
    
    # Shrunk predictive: all except x_4
    shrunk_predictive_str <- paste(paste0("arm:", shrunk_vars), collapse = " + ")
    
    # Unshrunk predictive: only x_4
    unshrunk_predictive_str <- paste0("arm:", unshrunk_var)

    # Set unique cmdstanr output directory per task to avoid conflicts
    job_temp_base <- Sys.getenv("JOB_TEMP_DIR", tempdir())
    safe_sim_id <- gsub("[^[:alnum:]_]", "_", current_task$sim_id)
    safe_prior_name <- gsub("[^[:alnum:]_]", "_", current_task$prior_name)

    # Add random suffix to cmdstan dir to prevent compilation collisions
    random_suffix <- paste0(sample(c(0:9, letters), 6, replace = TRUE), collapse = "")
    cmdstan_output_dir <- file.path(
      job_temp_base,
      sprintf("cmdstan_%s_%s_%s", safe_prior_name, safe_sim_id, random_suffix)
    )
    dir.create(cmdstan_output_dir, recursive = TRUE, showWarnings = FALSE)

    fit <- run_brms_analysis(
      data = df,
      response_formula_str = endpoint_params$resp_formula,
      response_type = endpoint_params$resptype,
      unshrunk_prognostic_formula_str = paste("~", prognostic_str),
      shrunk_predictive_formula_str = paste("~", shrunk_predictive_str),
      unshrunk_predictive_formula_str = paste("~", unshrunk_predictive_str),
      predictive_effect_priors = list(
        shrunk = prior_spec$prior,
        unshrunk = prior_spec$unshrunk_prior
      ),
      stanvars = prior_spec$stanvars,
      chains = 4, iter = 2000, warmup = 1000, cores = 1,
      backend = "cmdstanr",
      output_dir = cmdstan_output_dir
    )

    # Summarize results
    output <- summary_subgroup_effects(
      fit,
      df,
      "arm",
      endpoint_params$resptype,
      subgroup_vars = covariate_set  # Still report all subgroups
    )

    # Clean up cmdstan output directory to save disk space
    unlink(cmdstan_output_dir, recursive = TRUE)

    output$estimates

  }, error = function(e) {
    clean_message <- gsub("\\n", " ", e$message)
    cat(sprintf("[%s] ERROR in %s: %s\n", Sys.time(), task_description, clean_message),
        file = log_file, append = TRUE)
    # Clean up cmdstan directory even on error
    if (exists("cmdstan_output_dir") && dir.exists(cmdstan_output_dir)) {
      unlink(cmdstan_output_dir, recursive = TRUE)
    }
    return(tibble(error = clean_message))
  })

  # --- 4. Final Logging & Formatting ---
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("[%s] FINISHED: %s | Duration: %.2f mins\n", 
              format(end_time, "%Y-%m-%d %H:%M:%S"), task_description, duration),
      file = log_file, append = TRUE)

  if (is.data.frame(results) && nrow(results) > 0) {
    results %>%
      mutate(
        model_type = current_task$model_type,
        prior_name = current_task$prior_name,
        .before = 1
      )
  } else {
    tibble(
      model_type = current_task$model_type,
      prior_name = current_task$prior_name,
      error = ifelse(is.data.frame(results), "No results returned", results$error[1])
    )
  }
}


# --- 7. RUN IN PARALLEL & COMBINE ---

# --- Pre-compile step: Run task 1 serially ---
cat("--- Starting pre-compile step for task 1 ---\n")
result_task_1 <- list(run_single_task_unshrunk(1))
cat("--- Pre-compile task 1 finished. ---\n")

# --- Parallel step: Run tasks 2 to N ---
if (total_tasks > 1) {
  cat(sprintf("--- Starting parallel run for remaining %d tasks ---\n", total_tasks - 1))
  results_tasks_2_to_N <- mclapply(
    X = 2:total_tasks,
    FUN = run_single_task_unshrunk,
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

# --- 8. FINAL FORMATTING AND SAVING ---
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
cat(sprintf("Final results saved to %s\n", results_file))
