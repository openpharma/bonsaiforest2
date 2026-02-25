# =========================================================================
# === ONE-WAY SHRINKAGE MODEL (OVAT) SCRIPT - CONTINUOUS OUTCOMES
# === PARALLEL SIMULATION SCRIPT (OVAT + mclapply + UPDATE METHOD)
# ===
# === STRATEGY:
# === 1. Process only Scenario 1 data
# === 2. Run Sim_ID 1 Serially to compile models for all subgrouping variables
# === 3. Store these compiled models in memory
# === 4. Run remaining tasks in parallel using update() (FAST)
# =========================================================================

# --- 0. CONFIGURATION ---

# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "continuous" # <-- Set to "continuous" for continuous modeling
RESULTS_DIR <- "Results"
SCENARIO_ID <- 1  # Process only scenario 1

# Set unique seed based on script for reproducibility
RNGkind('Mersenne-Twister')
set.seed(0)

# Define delta_plan = abs(log(0.70))
delta_plan <- abs(log(0.70))

# --- Source functions (assuming existence in parent directory) ---
# source('../functions.R')

# --- Define priors for one-way shrinkage models ---
stanvars_hn_phi_1 <- brms::stanvar(
  scode = "real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = "  // Half-Normal prior on sigma_pred\n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); \n",
    block = "model"
  )

stanvars_hn_phi_delta <- brms::stanvar(
  scode = "real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = "  // Half-Normal prior on sigma_pred\n  target += normal_lpdf(sigma_pred | 0, ", delta_plan, ") - normal_lccdf(0 | 0, ", delta_plan, "); \n",
    block = "model"
  )

stanvars_hn_phi_delta_half <- brms::stanvar(
  scode = "real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = "  // Half-Normal prior on sigma_pred\n  target += normal_lpdf(sigma_pred | 0, ", delta_plan / 2, ") - normal_lccdf(0 | 0, ", delta_plan / 2, "); \n",
    block = "model"
  )

PRIOR_SPECIFICATIONS <- list(
  HN_phi_1 = list(prior = brms::set_prior("normal(0, sigma_pred)"), stanvars = stanvars_hn_phi_1),
  HN_phi_delta_plan = list(prior = brms::set_prior("normal(0, sigma_pred)"), stanvars = stanvars_hn_phi_delta),
  HN_phi_delta_plan_half = list(prior = brms::set_prior("normal(0, sigma_pred)"), stanvars = stanvars_hn_phi_delta_half)
)

# --- 1. LOAD LIBRARIES AND SET UP ---
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
library(brms)
library(bonsaiforest2) # <-- Ensure this is loaded before processing

# Set threads to 1 for each parallel process
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

# --- 2. DEFINE ENDPOINT-SPECIFIC PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

# Determine the script's directory to make paths work from anywhere
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- getwd()
}

endpoint_params <- list(
  folder = script_dir,
  resp_formula = "Y ~ trt",
  resptype = "continuous"
)

# --- 3. LOAD AND PROCESS INPUT DATA ---
cat(sprintf("Loading scenario %d data only...\n", SCENARIO_ID))
scenarios_dir <- file.path(endpoint_params$folder, "Scenarios")
if (!dir.exists(scenarios_dir)) stop("Data directory not found")

# Load only the specified scenario
scenario_file <- file.path(scenarios_dir, sprintf("scenario%d.rds", SCENARIO_ID))
if (!file.exists(scenario_file)) stop(sprintf("Scenario %d data not found", SCENARIO_ID))

scenario_data <- readRDS(scenario_file)

# Convert to list format by sim_id
all_datasets_processed <- scenario_data %>%
  group_by(sim_id) %>%
  nest() %>%
  pull(data, name = sim_id) %>%
  as.list()

flat_named_list <- all_datasets_processed
rm(scenario_data, all_datasets_processed); gc()

cat(sprintf("Successfully loaded %d datasets from scenario %d.\n", length(flat_named_list), SCENARIO_ID))

# --- 4. DEFINE PARAMETERS AND GRID ---
# Subgrouping variables to test
subgr_vars <- c("X1", "X2", "X3", "X4", "X8", "X11cat", "X14cat", "X17cat")
num_cores <- as.integer(parallelly::availableCores() - 2)
num_cores <- max(1, num_cores)

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("OVAT"),
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)
total_tasks <- nrow(task_grid)

# --- 5. LOGGING SETUP ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "ovat_1_oneway", prior_names_str, sep = "_")
results_dir_path <- file.path(endpoint_params$folder, RESULTS_DIR)
dir.create(results_dir_path, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(results_dir_path, paste0(base_file_name, "_scenario", SCENARIO_ID, ".log"))
results_file <- file.path(results_dir_path, paste0(base_file_name, "_scenario", SCENARIO_ID, ".rds"))

if (file.exists(log_file)) file.remove(log_file)

cat(sprintf("Total tasks: %d\n", total_tasks))
cat(sprintf("Results: %s\nLog: %s\n", results_file, log_file))
cat(sprintf("Using %d cores for parallel processing\n", num_cores))

# =========================================================================
# === PART A: PRE-COMPILE MODELS (Task 1, Serially)
# =========================================================================

cat("\n=== PART A: Pre-compiling base models (Task 1) ===\n")

# Load first dataset
df <- flat_named_list[[1]]

GLOBAL_COMPILED_MODELS <- list()

for (covariate in subgr_vars) {
  cat(sprintf("  Pre-compiling model for %s...\n", covariate))
  
  # Convert factor to formula-safe name
  covar_formula <- covariate
  
  # Ensure the covariate exists in the data
  if (!covar_formula %in% names(df)) {
    warning(sprintf("Covariate %s not found in data. Skipping.\n", covar_formula))
    next
  }

  prepared_run <- bonsaiforest2::prepare_formula_model(
    data = df,
    response_formula = endpoint_params$resp_formula,
    response_type = endpoint_params$resptype,
    unshrunk_terms_formula = paste("~", covar_formula, "+ trt"),
    shrunk_predictive_formula = paste("~ (0 + trt ||", covar_formula, ")")
  )

  base_prior <- PRIOR_SPECIFICATIONS[[1]]$prior
  base_stanvars <- PRIOR_SPECIFICATIONS[[1]]$stanvars

  GLOBAL_COMPILED_MODELS[[covariate]] <- brms::brm(
    formula = prepared_run$formula,
    data = prepared_run$data,
    family = brms::brmsfamily(endpoint_params$resptype),
    prior = base_prior,
    stanvars = base_stanvars,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    refresh = 0,
    verbose = FALSE,
    backend = "cmdstanr"
  )
}

results_list <- list()

# =========================================================================
# === PART B: RUN FUNCTION FOR update() TASKS
# =========================================================================

run_update_task <- function(i) {
  current_task <- task_grid[i, ]
  sim_id <- current_task$sim_id
  df <- flat_named_list[[sim_id]]
  
  start_t <- Sys.time()
  task_desc <- sprintf("Sim %s | Model: %s | Prior: %s", sim_id, current_task$model_type, current_task$prior_name)
  
  cat(sprintf("[%s] START: %s\n", format(start_t, "%H:%M:%S"), task_desc), file = log_file, append = TRUE)
  
  out <- tryCatch({

    cov_results <- list()

    for (covariate in subgr_vars) {
      
      # Skip if model wasn't compiled (e.g., due to missing covariate)
      if (!covariate %in% names(GLOBAL_COMPILED_MODELS)) {
        next
      }

      # --- KEY STEP: Retrieve the pre-compiled model ---
      base_model <- GLOBAL_COMPILED_MODELS[[covariate]]

      # Prepare NEW data with random effects notation
      prepared_run <- bonsaiforest2::prepare_formula_model(
        data = df,
        response_formula = endpoint_params$resp_formula,
        response_type = endpoint_params$resptype,
        unshrunk_terms_formula = paste("~", covariate, "+ trt"),
        shrunk_predictive_formula = paste("~ (0 + trt ||", covariate, ")")
      )

      # Get the prior for this task
      task_prior <- PRIOR_SPECIFICATIONS[[current_task$prior_name]]$prior
      task_stanvars <- PRIOR_SPECIFICATIONS[[current_task$prior_name]]$stanvars

      # --- KEY STEP: UPDATE() ---
      updated_fit <- update(
        base_model,
        newdata = prepared_run$data,
        prior = task_prior,
        stanvars = task_stanvars,
        recompile = FALSE,
        refresh = 0 # Silence Stan output
      )

      # Summarize
      summ <- bonsaiforest2::summary_subgroup_effects(updated_fit)
      est_df <- summ$estimates
      
      # Extract convergence diagnostics
      rhat_vals <- brms::rhat(updated_fit)
      max_rhat <- ifelse(!is.null(rhat_vals), max(rhat_vals, na.rm = TRUE), NA)
      rhat_gt_1_05 <- ifelse(!is.null(rhat_vals), sum(rhat_vals > 1.05, na.rm = TRUE), NA)
      
      # Get divergent transitions
      n_divergent <- sum(updated_fit$sampler_diagnostics[, , "divergent__"])
      
      # Get ESS estimates
      neff_ratio <- brms::neff_ratio(updated_fit)
      mean_neff_ratio <- ifelse(!is.null(neff_ratio), mean(neff_ratio, na.rm = TRUE), NA)
      
      # Add diagnostics to each row
      est_df$covariate_tested <- covariate
      est_df$max_rhat <- max_rhat
      est_df$n_rhats_gt_1_05 <- rhat_gt_1_05
      est_df$n_divergent <- n_divergent
      est_df$mean_neff_ratio <- mean_neff_ratio
      
      cov_results[[covariate]] <- est_df
    }

    bind_rows(cov_results)

  }, error = function(e) {
    msg <- gsub("\\n", " ", e$message)
    cat(sprintf("ERROR %s: %s\n", task_desc, msg), file = log_file, append = TRUE)
    tibble(error = msg)
  })

  # Add identifiers
  out <- out %>%
    mutate(
      sim_id = sim_id,
      scenario_id = SCENARIO_ID,
      prior_spec = current_task$prior_name,
      .before = 1
    )

  end_t <- Sys.time()
  elapsed <- difftime(end_t, start_t, units = "secs")
  cat(sprintf("  END (%0.1f sec)\n", elapsed), file = log_file, append = TRUE)

  return(out)
}

# =========================================================================
# === PART C: PARALLEL EXECUTION (Tasks 2 through N)
# =========================================================================

cat(sprintf("\n=== PART C: Running %d tasks in parallel across %d cores ===\n", 
            total_tasks, num_cores))

results_list <- mclapply(
  1:total_tasks,
  run_update_task,
  mc.cores = num_cores,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE
)

# --- Combine results ---
all_results <- bind_rows(results_list)

# --- Save results ---
saveRDS(all_results, file = results_file)
cat(sprintf("\n✓ Results saved to %s\n", results_file))
cat(sprintf("✓ Log saved to %s\n", log_file))
cat(sprintf("✓ Processed %d tasks successfully\n", nrow(all_results)))
