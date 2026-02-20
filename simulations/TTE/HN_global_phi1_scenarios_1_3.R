# -----------------------------------------------------------------
# RUN GLOBAL SHRINKAGE MODEL ANALYSIS (bonsaiforest2) FOR A SINGLE ENDPOINT
# Prior: Half-Normal with φ = 1
# Scenarios: 1-3
#
# This script:
# 1. Sources all functions from `functions.R`.
# 2. Runs the global shrinkage model analysis using `bonsaiforest2` for
#    the one endpoint defined in `ENDPOINT_ID`.
# 3. Uses Half-Normal(φ=1) prior on fixed effects.
# 4. Processes SCENARIOS 1-3 only.
# 5. Saves the results to the `Results/` folder.
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---

# SET THIS VARIABLE to the endpoint you want to analyze
# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "tte"

# SET THE MAIN RESULTS DIRECTORY
RESULTS_DIR <- "Results"

# Define delta_plan = abs(log(0.70))
delta_plan <- abs(log(0.70))

# Define the prior: Half-Normal with phi = 1
phi <- 1

message(paste("Defining Half-Normal prior with phi =", phi, "..."))

stanvars_hn_phi_1 <- brms::stanvar(
  scode = "real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = paste0("  // Half-Normal prior on sigma_pred\n  target += normal_lpdf(sigma_pred | 0, ", phi, ") - normal_lccdf(0 | 0, ", phi, "); \n"),
    block = "model"
  )

prior_hn_phi_1 <- brms::set_prior("normal(0, sigma_pred)")

PRIOR_SPECIFICATIONS <- list(
  HN_global_phi_1 = list(prior = prior_hn_phi_1, stanvars = stanvars_hn_phi_1)
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
                            folder = "binary",
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
                            resp_formula = y ~ arm,
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

covariate_set <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")
num_cores <- 96

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("Global"), # Global shrinkage model
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)

# --- 5. SET UP LOGGING AND RESULTS FILES ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "global", prior_names_str, "scenarios_1_3", sep = "_")

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
    predictive_str <- paste(paste0("arm:", covariate_set), collapse = " + ")
    
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
    
    # Get ESS estimates
    neff_ratio <- brms::neff_ratio(fit)
    mean_neff_ratio <- ifelse(!is.null(neff_ratio), mean(neff_ratio, na.rm = TRUE), NA)
    
    # Add diagnostics to each row
    est_df$max_rhat <- max_rhat
    est_df$n_rhats_gt_1_05 <- rhat_gt_1_05
    est_df$n_divergent <- n_divergent
    est_df$mean_neff_ratio <- mean_neff_ratio
    
    # Clean up
    unlink(cmdstan_output_dir, recursive = TRUE)
    
    est_df
    
  }, error = function(e) {
    msg <- gsub("\\n", " ", e$message)
    cat(sprintf("ERROR %s: %s\n", task_desc, msg), file = log_file, append = TRUE)
    if (exists("cmdstan_output_dir") && dir.exists(cmdstan_output_dir)) {
      unlink(cmdstan_output_dir, recursive = TRUE)
    }
    tibble(error = msg)
  })

  # Add identifiers
  if (is.data.frame(out) && !"error" %in% names(out)) {
    out$model_type <- current_task$model_type
    out$prior_name <- current_task$prior_name
    out <- out %>% dplyr::select(model_type, prior_name, everything())
  }

  # Log completion
  dur <- difftime(Sys.time(), start_t, units = "mins")
  cat(sprintf("[%s] DONE: %s (%.2f m)\n", format(Sys.time(), "%H:%M:%S"), task_desc, dur),
      file = log_file, append = TRUE)

  return(out)
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
