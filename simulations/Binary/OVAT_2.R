# -----------------------------------------------------------------
# RUN OVAT MODEL ANALYSIS (bonsaiforest2 - NO UPDATE)
#
# This script:
# 1. Sources all functions from `functions.R`.
# 2. Runs the "OVAT" model analysis (fitting from scratch) for
#    the one endpoint defined in `ENDPOINT_ID`.
# 3. Allows custom prior definitions.
# 4. Saves the results to the `Results/` folder.
#
# !! WARNING: This script is much slower than the `update()`
#    version as it re-fits every model.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---

# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "binary"

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"

# ------------------------


# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(parallel)
library(parallelly)
library(dplyr)
library(tibble)
library(purrr)
library(checkmate)
library(cmdstanr)
library(RhpcBLASctl)
library(stringr)
library(tidyr)
library(bonsaiforest2) # Assumed to be loaded
library(brms)        # For set_prior/stanvar
library(survival)      # For Surv()

# Set threads for parallel safety
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

# DEFINE YOUR PRIORS HERE
cat("Defining custom full hierarchical prior (N(mu, sigma^2))...\n")
stanvars_full_hierarchical <- stanvar(
  scode = "  real mu_pred;\n  real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  stanvar(
    scode = "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); // N(0, 16) -> SD=4\n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); // Half-Normal(1)\n",
    block = "model"
  )
prior_full_hierarchical <- set_prior("normal(mu_pred, sigma_pred)")

PRIOR_SPECIFICATIONS <- list(
  hierarchical = list(prior = prior_full_hierarchical, stanvars = stanvars_full_hierarchical)
)

# Reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# Source all helper functions (from the same directory)
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
cat(sprintf("--- Starting FULL OVAT MODEL RUN (%s) ---\n", ENDPOINT_ID))
cat("Loading all scenario data...\n")

subgr_vars <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")
scenarios_to_run <- as.character(1:6)
scenarios_list <- list()

for (scen in scenarios_to_run) {
  scen_file <- file.path(endpoint_params$folder, "Scenarios", paste0("scenario", scen, ".rds"))
  if (!file.exists(scen_file)) {
    stop(paste("File not found:", scen_file))
  }
  message(paste("    Loading:", scen_file))
  scenarios_list[[scen]] <- readRDS(scen_file)
}

# Process all datasets from all scenarios into a single flat list
all_datasets_processed <- lapply(seq_along(scenarios_list), function(i) {
  scenario_data <- scenarios_list[[i]]
  setNames(scenario_data, paste(i, 1:length(scenario_data), sep = "_"))
})

flat_named_list <- unlist(all_datasets_processed, recursive = FALSE)
flat_named_list <-flat_named_list[1001:2000]

cat(sprintf("Successfully loaded %d datasets.\n", length(flat_named_list)))
rm(scenarios_list, all_datasets_processed)
gc()


# --- 4. DEFINE TASK GRID & FILE NAMES ---
task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("OVAT"),
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)

total_tasks <- nrow(task_grid)
num_cores <- min(total_tasks, parallelly::availableCores()) # Use all tasks if fewer than cores
if (num_cores < 1) num_cores <- 1

# Create dynamic file names
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "ovat_2", prior_names_str, "no_update", sep = "_")

# *** (FIXED FILE PATHS) ***
# Ensure the log and results files go into the correct sub-directory
log_file <- file.path(endpoint_params$folder, RESULTS_DIR, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path(endpoint_params$folder, RESULTS_DIR, paste0(base_file_name, ".rds"))
# *** (END FIX) ***

cat(sprintf("\nStarting %d total tasks using %d cores...\n", total_tasks, num_cores))
cat(sprintf("Progress will be logged to: %s\n\n", log_file))


# --- 5. SURVIVAL KNOT-PATCHING (TTE ONLY) ---
# This ensures that the 10 OVAT models *within* a single task
# don't all recompile. This is still critical for TTE.

if (ENDPOINT_ID == "tte") {
  cat("Applying monkey-patch to fix knots for TTE...\n")

  time_var_name <- stringr::str_trim(stringr::str_match(endpoint_params$resp_formula, "Surv\\((.*?),(.*?)\\)")[, 2])
  FIXED_KNOTS <- .calculate_bhaz_knots(flat_named_list[[1]][[time_var_name]])
  original_knots_fun <- .calculate_bhaz_knots
  fake_knots_fun <- function(time_data) { return(FIXED_KNOTS) }

  try(unlockBinding(".calculate_bhaz_knots", .GlobalEnv), silent = TRUE)
  assign(".calculate_bhaz_knots", fake_knots_fun, .GlobalEnv)

  on.exit({
    cat("Restoring original .calculate_bhaz_knots function...\n")
    try(unlockBinding(".calculate_bhaz_knots", .GlobalEnv), silent = TRUE)
    assign(".calculate_bhaz_knots", original_knots_fun, .GlobalEnv)
    try(lockBinding(".calculate_bhaz_knots", .GlobalEnv), silent = TRUE)
  }, add = TRUE)

  cat("Knots fixed. Patch applied.\n")
}


# --- 6. Define the Parallel "Fit from Scratch" Function ---
# This replaces both the serial compile and parallel update steps

run_parallel_ovat_fit_task <- function(i) {
  # Set threads *again* just to be safe in the new process
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(threads = 1)

  # --- Get task info from grid ---
  current_task <- task_grid[i, ]
  sim_id <- current_task$sim_id
  df <- flat_named_list[[sim_id]]
  prior_spec <- PRIOR_SPECIFICATIONS[[current_task$prior_name]]

  # --- Logging Setup ---
  task_description <- sprintf("SimID: %s | Model: %s | Prior: %s | Endpoint: %s",
                              sim_id, current_task$model_type, current_task$prior_name, ENDPOINT_ID)
  start_time <- Sys.time()
  cat(sprintf("[%s] STARTING (Task %d): %s\n", format(start_time, "%Y-%m-%d %H:%M:%S"), i, task_description),
      file = log_file, append = TRUE)

  result_for_task <- tryCatch({

    # --- Loop over all covariates and FIT FROM SCRATCH ---
    ovat_results_for_task <- map(subgr_vars, function(covariate) {

      # --- A. Prepare data for THIS covariate ---
      prognostic_str <- covariate
      predictive_str <- paste0("arm:", covariate)

      prepared_run <- prepare_formula_model(
        data = df,
        response_formula_str = endpoint_params$resp_formula, # Dynamic
        response_type = endpoint_params$resptype,         # Dynamic
        unshrunk_prognostic_formula_str = paste("~", prognostic_str),
        shrunk_predictive_formula_str = paste("~", predictive_str)
      )

      # --- B. Fit from scratch (slow) ---
      fitted_model_cov <- fit_brms_model(
        formula = prepared_run$formula,
        data = prepared_run$data,
        response_type = endpoint_params$resptype, # Dynamic
        predictive_effect_priors = list(shrunk = prior_spec$prior),
        stanvars = prior_spec$stanvars,
        chains = 4, iter = 2000, warmup = 1000, cores = 1, backend = "cmdstanr"
      )

      # --- C. Get results for this covariate ---
      summary_obj <- summary_subgroup_effects(
        fitted_model_cov, df, "arm",
        endpoint_params$resptype, # Dynamic
        subgroup_vars = 'auto'
      )

      # Return the estimates
      summary_obj$estimates
    })

    # --- D. Combine all covariate summaries into ONE data frame ---
    bind_rows(ovat_results_for_task)

  }, error = function(e) {
    clean_message <- gsub("\\n", " ", e$message)
    cat(sprintf("[%s] ERROR in (Task %d) %s: %s\n", Sys.time(), i, task_description, clean_message),
        file = log_file, append = TRUE)
    return(tibble(error = clean_message)) # Return a tibble with error
  })

  # --- Logging & Final Formatting ---
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("[%s] FINISHED (Task %d): %s | Duration: %.2f mins\n", format(end_time, "%Y-%m-%d %H:%M:%S"), i, task_description, duration),
      file = log_file, append = TRUE)

  # Add identifying columns
  if (is.data.frame(result_for_task) && nrow(result_for_task) > 0) {
    result_for_task %>%
      mutate(
        model_type = current_task$model_type,
        prior_name = current_task$prior_name,
        .before = 1
      )
  } else {
    tibble(
      model_type = current_task$model_type,
      prior_name = current_task$prior_name,
      error = ifelse(is.data.frame(result_for_task), "No results returned", "Unknown error")
    )
  }
}

# --- 7. Run All Tasks (Serial Compile + Parallel Run) ---
cat(sprintf("\n--- Starting SERIAL COMPILE run for task 1... ---\n"))
cat(sprintf("(This will compile all 10 OVAT models for endpoint '%s'...)\n", ENDPOINT_ID))

# Run task 1 serially to compile all models.
# We put its result in a list to match the mclapply output.
results_list_serial <- list(
  run_parallel_ovat_fit_task(1)
)
cat("--- Serial compile task finished. ---\n")


# Initialize an empty list for parallel results
results_list_parallel <- list()

# Only run in parallel if there are tasks left
if (total_tasks > 1) {

  # Adjust number of cores if we used one for the serial task (if total_tasks < num_cores)
  # This is a bit safer but num_cores is usually fine.
  parallel_cores <- min((total_tasks - 1), num_cores)

  cat(sprintf("\n--- Starting PARALLEL run for remaining %d tasks on %d cores... ---\n",
              (total_tasks - 1), parallel_cores))

  # Run the rest in parallel (tasks 2 through total_tasks)
  results_list_parallel <- mclapply(
    2:total_tasks,
    FUN = run_parallel_ovat_fit_task,
    mc.cores = parallel_cores,
    mc.preschedule = FALSE
  )

} else {
  cat("--- Only 1 total task. No parallel run needed. ---\n")
}


# Combine the single serial result with all the parallel results
results_list <- c(results_list_serial, results_list_parallel)


# --- 8. Combine All Results & Save ---
cat("\n--- All tasks finished. Combining results... ---\n")

# Bind all data frames together
final_results_df <- bind_rows(results_list)

# Add back the simulation and replication identifiers.
if (nrow(final_results_df) > 0) {

  # This logic needs to be robust to errors
  rows_per_task <- sapply(results_list, function(res) {
    # Check if it's a valid data frame, otherwise it's an error row (count as 1)
    if (is.data.frame(res) && nrow(res) > 0) {
      nrow(res)
    } else {
      1 # Count 1 for error rows or empty tibbles
    }
  })

  # Check for length mismatch (can happen if sapply failed, though unlikely)
  if(length(rows_per_task) == total_tasks) {
    task_indices <- rep(1:total_tasks, times = rows_per_task)

    final_results_df <- final_results_df %>%
      mutate(
        sim_id_full = task_grid$sim_id[task_indices],
        scenario_id = as.integer(stringr::str_extract(sim_id_full, "^\\d+")),
        replication_id = as.integer(stringr::str_extract(sim_id_full, "\\d+$")),
        .before = 1
      ) %>%
      select(-sim_id_full) # Clean up temporary column
  } else {
    cat("WARNING: Could not re-map task IDs. `rows_per_task` length did not match `total_tasks`.\n")
  }
}

# --- Save the final file ---
dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  object = final_results_df,
  file = results_file
)
cat(sprintf("Final adapted results saved to %s\n", results_file))
cat("\n--- Loop finished. ---\n")
