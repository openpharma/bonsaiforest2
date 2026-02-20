# -----------------------------------------------------------------
# RUN OVAT MODEL ANALYSIS - PIPE-PIPE (||) SYNTAX
#
# This script:
# 1. Sources functions directly from the R folder (library not updated yet)
# 2. Runs the "OVAT" model analysis with PIPE-PIPE (||) syntax
# 3. Continuous endpoint only
# 4. Uses hierarchical priors
# 5. Saves the results to the `Results/` folder.
#
# !! WARNING: This script is much slower than the `update()`
#    version as it re-fits every model.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"
ENDPOINT_ID <- "continuous"

# Set to TRUE to use hierarchical prior, FALSE to use fixed prior
USE_HIERARCHICAL_PRIOR <- FALSE

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
library(brms)          # For set_prior/stanvar
library(survival)      # For Surv()

# Set threads for parallel safety
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

# --- DEFINE PRIORS ---
cat("Defining priors for PIPE-PIPE (||) syntax...\n")

if (USE_HIERARCHICAL_PRIOR) {
  cat("Using HIERARCHICAL prior: N(mu_pred, sigma_pred)\n")
  stanvars_prior <- stanvar(
    scode = "  real mu_pred;\n  real<lower=0> sigma_pred;\n",
    block = "parameters"
  ) +
    stanvar(
      scode = "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); // N(0, 16) -> SD=4\n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); // Half-Normal(1)\n",
      block = "model"
    )
  prior_spec <- "normal(mu_pred, sigma_pred)"

  PRIOR_SPECIFICATIONS <- list(
    hierarchical_pipe_pipe = list(prior = prior_spec, stanvars = stanvars_prior)
  )
} else {
  cat("Using FIXED prior: N(0, tau) with tau = nu * se_trt_ref\n")
  se_trt_ref <- 0.5
  nu <- qnorm(1 - 0.1) + qnorm(1 - 0.025)
  prior_tau <- nu * se_trt_ref

  prior_spec <- paste0("normal(0, ", prior_tau, ")")

  PRIOR_SPECIFICATIONS <- list(
    fixed_pipe_pipe = list(prior = prior_spec, stanvars = NULL)
  )
}

# Reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# --- SOURCE FUNCTIONS DIRECTLY FROM R FOLDER ---
message("--- Sourcing functions from ../R/ folder ---")
source("../R/prepare_formula_model.R")
source("../R/fit_brms_model.R")
source("../R/estimate_subgroup_effects.R")

# Also load helper functions from simulations folder
source("functions.R")

# --- 2. DEFINE ENDPOINT-SPECIFIC PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

endpoint_params <- list(
  folder = "Continuous",
  resp_formula = "y ~ arm",
  resptype = "continuous"
)

# --- 3. LOAD AND PROCESS INPUT DATA ---
cat(sprintf("--- Starting FULL OVAT MODEL RUN with PIPE-PIPE (||) syntax (%s) ---\n", ENDPOINT_ID))
cat("Loading all scenario data...\n")

subgr_vars <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")
scenarios_to_run <- as.character(1:6)
scenarios_list <- list()

for (scen in scenarios_to_run) {
  scen_file <- file.path("Continuous/Scenarios", paste0("scenario", scen, ".rds"))
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
flat_named_list <- flat_named_list[1:1000]  # Limit to 1000 for faster testing

cat(sprintf("Successfully loaded %d datasets.\n", length(flat_named_list)))
rm(scenarios_list, all_datasets_processed)
gc()


# --- 4. DEFINE TASK GRID & FILE NAMES ---
task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("OVAT_PIPE_PIPE"),
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)

total_tasks <- nrow(task_grid)
num_cores <- min(total_tasks, parallelly::availableCores())
if (num_cores < 1) num_cores <- 1

# Create dynamic file names
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "ovat_1_pipe_pipe", prior_names_str, "trial", sep = "_")

# Ensure the log and results files go into the correct sub-directory
log_file <- file.path("Continuous", RESULTS_DIR, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path("Continuous", RESULTS_DIR, paste0(base_file_name, ".rds"))

cat(sprintf("\nStarting %d total tasks using %d cores...\n", total_tasks, num_cores))
cat(sprintf("Progress will be logged to: %s\n\n", log_file))


# --- 5. Define the Parallel OVAT Function with PIPE-PIPE Syntax ---
run_ovat_model_task <- function(sim_id,
                                model_type,
                                prior_name,
                                all_data,
                                subgr_vars,
                                log_path,
                                prior_spec,
                                endpoint_params) {

  # Set threads
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(threads = 1)

  # --- Get task info ---
  df <- all_data[[sim_id]]

  # --- Logging Setup ---
  task_description <- sprintf("SimID: %s | Model: %s | Prior: %s",
                              sim_id, model_type, prior_name)
  start_time <- Sys.time()
  cat(sprintf("[%s] STARTING: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S"), task_description),
      file = log_path, append = TRUE)

  result_for_task <- tryCatch({

    # Initialize list for results
    ovat_results_list <- list()

    # --- Loop over all covariates ---
    for (covariate in subgr_vars) {

      tryCatch({

        # --- A. Prepare data with PIPE-PIPE (||) syntax ---
        prognostic_str <- covariate
        # Use PIPE-PIPE syntax: (arm || covariate)
        predictive_str <- paste0("(arm || ", covariate, ")")

        prepared_run <- prepare_formula_model(
          data = df,
          response_formula_str = endpoint_params$resp_formula,
          response_type = endpoint_params$resptype,
          unshrunk_prognostic_formula_str = paste("~", prognostic_str),
          shrunk_predictive_formula_str = paste("~", predictive_str)
        )

        # --- B. Fit from scratch ---
        fitted_model_cov <- fit_brms_model(
          prepared_model = prepared_run,
          predictive_effect_priors = list(shrunk = prior_spec$prior),
          stanvars = prior_spec$stanvars,
          chains = 4, iter = 2000, warmup = 1000, cores = 1, backend = "cmdstanr"
        )

        # --- C. Get results ---
        summary_obj <- estimate_subgroup_effects(
          brms_fit = fitted_model_cov,
          original_data = df,
          trt_var = "arm",
          response_type = endpoint_params$resptype,
          subgroup_vars = "auto"
        )

        # Extract the estimates tibble and add covariate identifier
        result_df <- as.data.frame(summary_obj$estimates)
        result_df$covariate <- covariate
        ovat_results_list[[covariate]] <- result_df

      }, error = function(e) {
        stop(sprintf("Error in covariate '%s': %s", covariate, e$message))
      })

    } # End for loop

    # --- D. Combine all results ---
    combined_results <- bind_rows(ovat_results_list)

    # Reorder columns to put covariate first (or adjust as needed)
    if ("covariate" %in% names(combined_results)) {
      combined_results <- combined_results %>%
        select(covariate, everything())
    }

    combined_results

  }, error = function(e) {
    clean_message <- gsub("\\n", " ", e$message)
    cat(sprintf("[%s] ERROR in %s: %s\n", Sys.time(), task_description, clean_message),
        file = log_path, append = TRUE)
    return(tibble(error = clean_message))
  })

  # --- Logging & Final Formatting ---
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_status <- if (any(names(result_for_task) == "error")) "ERROR" else "FINISHED"

  cat(sprintf("[%s] %s: %s | Duration: %.2f mins\n",
              format(end_time, "%Y-%m-%d %H:%M:%S"),
              log_status,
              task_description,
              duration),
      file = log_path, append = TRUE)

  if (is.data.frame(result_for_task) && nrow(result_for_task) > 0 && !"error" %in% names(result_for_task)) {
    result_for_task %>%
      mutate(
        model_type = model_type,
        prior_name = prior_name,
        .before = 1
      )
  } else {
    tibble(
      model_type = model_type,
      prior_name = prior_name,
      error = ifelse(is.data.frame(result_for_task) && "error" %in% names(result_for_task),
                     result_for_task$error[1],
                     "Unknown error")
    )
  }
}

# --- 6. Run All Tasks (Serial Compile + Parallel Run) ---
cat(sprintf("\n--- Starting SERIAL COMPILE run for task 1... ---\n"))
cat(sprintf("(This will compile all 10 OVAT models for endpoint '%s' with PIPE-PIPE syntax...)\n", ENDPOINT_ID))

# Run task 1 serially to compile all models
results_list_serial <- list(
  run_ovat_model_task(
    sim_id = task_grid$sim_id[1],
    model_type = task_grid$model_type[1],
    prior_name = task_grid$prior_name[1],
    all_data = flat_named_list,
    subgr_vars = subgr_vars,
    log_path = log_file,
    prior_spec = PRIOR_SPECIFICATIONS[[task_grid$prior_name[1]]],
    endpoint_params = endpoint_params
  )
)
cat("--- Serial compile task finished. ---\n")

# Initialize an empty list for parallel results
results_list_parallel <- list()

# Only run in parallel if there are tasks left
if (total_tasks > 1) {

  parallel_cores <- min((total_tasks - 1), num_cores)

  cat(sprintf("\n--- Starting PARALLEL run for remaining %d tasks on %d cores... ---\n",
              (total_tasks - 1), parallel_cores))

  # Run the rest in parallel
  results_list_parallel <- mclapply(
    2:total_tasks,
    FUN = function(task_idx) {
      run_ovat_model_task(
        sim_id = task_grid$sim_id[task_idx],
        model_type = task_grid$model_type[task_idx],
        prior_name = task_grid$prior_name[task_idx],
        all_data = flat_named_list,
        subgr_vars = subgr_vars,
        log_path = log_file,
        prior_spec = PRIOR_SPECIFICATIONS[[task_grid$prior_name[task_idx]]],
        endpoint_params = endpoint_params
      )
    },
    mc.cores = parallel_cores,
    mc.preschedule = FALSE
  )

} else {
  cat("--- Only 1 total task. No parallel run needed. ---\n")
}

# Combine the single serial result with all the parallel results
results_list <- c(results_list_serial, results_list_parallel)


# --- 7. Combine All Results & Save ---
cat("\n--- All tasks finished. Combining results... ---\n")

# Bind all data frames together
final_results_df <- bind_rows(results_list)

# Add back the simulation and replication identifiers
if (nrow(final_results_df) > 0) {

  rows_per_task <- sapply(results_list, function(res) {
    if (is.data.frame(res) && nrow(res) > 0) {
      nrow(res)
    } else {
      1
    }
  })

  if (length(rows_per_task) == total_tasks) {
    task_indices <- rep(1:total_tasks, times = rows_per_task)

    final_results_df <- final_results_df %>%
      mutate(
        sim_id_full = task_grid$sim_id[task_indices],
        scenario_id = as.integer(stringr::str_extract(sim_id_full, "^\\d+")),
        replication_id = as.integer(stringr::str_extract(sim_id_full, "\\d+$")),
        .before = 1
      ) %>%
      select(-sim_id_full)
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
cat(sprintf("Final results saved to %s\n", results_file))
cat("\n--- SIMULATION FINISHED ---\n")
