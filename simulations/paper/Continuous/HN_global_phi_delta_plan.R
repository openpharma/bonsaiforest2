# --- 0. CONFIGURATION ---

# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "continuous"

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"

# Define delta_plan = 0.35
delta_plan <- 0.35

# Define the prior: Half-Normal with phi = delta_plan
phi <- delta_plan

message(paste("Defining Half-Normal prior with phi = delta_plan =", phi, "..."))

stanvars_hn_phi_delta <- brms::stanvar(
  scode = "real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = paste0("  // Half-Normal prior on sigma_pred\n  target += normal_lpdf(sigma_pred | 0, ", phi, ") - normal_lccdf(0 | 0, ", phi, "); \n"),
    block = "model"
  )

prior_hn_phi_delta <- brms::set_prior("normal(0, sigma_pred)")

PRIOR_SPECIFICATIONS <- list(
  HN_global_phi_delta_plan = list(prior = prior_hn_phi_delta, stanvars = stanvars_hn_phi_delta)
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
                            resp_formula = "Surv(tt_pfs, ev_pfs) ~ trt",
                            resptype = "survival"
                          ),
                          "binary" = list(
                            folder = "binary",
                            resp_formula = "y ~ trt",
                            resptype = "binary"
                          ),
                          "count" = list(
                            folder = "Count",
                            resp_formula = "Y ~ trt",
                            resptype = "count"
                          ),
                          "continuous" = list(
                            folder = "Continuous",
                            resp_formula = "Y ~ trt",
                            resptype = "continuous"
                          ),
                          stop("Invalid ENDPOINT_ID. Must be one of: 'tte', 'binary', 'count', 'continuous'")
)

# --- 3. LOAD AND PROCESS INPUT DATA ---
cat(sprintf("--- Starting GLOBAL SHRINKAGE MODEL RUN (%s) - SCENARIOS 1-3 ---\n", ENDPOINT_ID))
cat("Loading scenario data for scenarios 1-3...\n")

scenarios_to_run <- as.character(1:3)
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
RNGkind('Mersenne-Twister') # Recommended for parallel reproducibility
set.seed(0)

covariate_set <- c("X1", "X2", "X3", "X4", "X8", "X11cat", "X14cat", "X17cat")
num_cores <- 96

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("Global"), # Global shrinkage model
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


# --- 6. RUN IN PARALLEL & COMBINE ---

run_single_task <- function(i) {
  current_task <- task_grid[i, ]
  sim_id <- current_task$sim_id
  df <- flat_named_list[[sim_id]]

  start_t <- Sys.time()
  task_desc <- sprintf("Sim %s | Model: %s | Prior: %s", sim_id, current_task$model_type, current_task$prior_name)

  cat(sprintf("[%s] START: %s\n", format(start_t, "%H:%M:%S"), task_desc), file = log_file, append = TRUE)

  out <- tryCatch({

    # Set up cmdstan directory
    job_temp_base <- Sys.getenv("JOB_TEMP_DIR", tempdir())
    safe_sim_id <- gsub("[^[:alnum:]_]", "_", sim_id)
    safe_prior_name <- gsub("[^[:alnum:]_]", "_", current_task$prior_name)
    random_suffix <- paste0(sample(c(0:9, letters), 6, replace = TRUE), collapse = "")
    cmdstan_output_dir <- file.path(
      job_temp_base,
      sprintf("cmdstan_%s_%s_%s", safe_prior_name, safe_sim_id, random_suffix)
    )
    dir.create(cmdstan_output_dir, recursive = TRUE, showWarnings = FALSE)

    # Build formulas
    prognostic_str <- paste(covariate_set, collapse = " + ")
    predictive_str <- paste(paste0("trt:", covariate_set), collapse = " + ")

    # Run model
    fit <- run_brms_analysis(
      data = df,
      response_formula = endpoint_params$resp_formula,
      response_type = endpoint_params$resptype,
      unshrunk_terms_formula = paste("~", prognostic_str),
      shrunk_predictive_formula = paste("~ 0 +", predictive_str),
      shrunk_predictive_prior = PRIOR_SPECIFICATIONS[[current_task$prior_name]]$prior,
      stanvars = PRIOR_SPECIFICATIONS[[current_task$prior_name]]$stanvars,
      chains = 4, iter = 2000, warmup = 1000, cores = 1,
      backend = "cmdstanr",
      output_dir = cmdstan_output_dir
    )

    # Get estimates
    output <- summary_subgroup_effects(fit)
    est_df <- output$estimates

    # Extract convergence diagnostics
    rhat_vals <- brms::rhat(fit)
    max_rhat <- ifelse(!is.null(rhat_vals), max(rhat_vals, na.rm = TRUE), NA)
    rhat_gt_1_05 <- ifelse(!is.null(rhat_vals), sum(rhat_vals > 1.05, na.rm = TRUE), NA)

    # Get divergent transitions
    n_divergent <- sum(fit$sampler_diagnostics[, , "divergent__"])

    # Get neff ratio
    neff_ratio <- mean(brms::neff_ratio(fit), na.rm = TRUE)

    # Add metadata
    est_df <- est_df %>%
      mutate(
        model_type = current_task$model_type,
        prior_name = current_task$prior_name,
        max_rhat = max_rhat,
        n_rhats_gt_1_05 = rhat_gt_1_05,
        n_divergent = n_divergent,
        mean_neff_ratio = neff_ratio
      )

    # Clean up model directory
    unlink(cmdstan_output_dir, recursive = TRUE)

    est_df

  }, error = function(e) {
    cat(sprintf("[ERROR] %s: %s\n", task_desc, as.character(e)), file = log_file, append = TRUE)
    tibble()
  })

  end_t <- Sys.time()
  cat(sprintf("[%s] DONE: %s (%.1f sec)\n", format(end_t, "%H:%M:%S"), task_desc,
              as.numeric(difftime(end_t, start_t, units = "secs"))), file = log_file, append = TRUE)

  out
}

# Pre-compile task 1
cat("Pre-compiling first task...\n")
results_pre1 <- run_single_task(1)

# Run remaining tasks in parallel
if (total_tasks > 1) {
  cat(sprintf("Running remaining %d tasks in parallel...\n", total_tasks - 1))
  results_remaining <- mclapply(2:total_tasks, run_single_task, mc.cores = num_cores)
} else {
  results_remaining <- list()
}

# Combine results
results_combined <- bind_rows(results_pre1, results_remaining)

# Save results to disk
cat(sprintf("Saving %d estimate rows to: %s\n", nrow(results_combined), results_file))
saveRDS(results_combined, file = results_file, compress = FALSE)

cat("Model run completed successfully!\n")
