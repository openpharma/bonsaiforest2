# =========================================================================
# === PARALLEL SIMULATION SCRIPT (OVAT + mclapply)
# ===
# === STRATEGY:
# === 1. Run Sim_ID 1 Serially to compile models for x_1...x_10
# === 2. Store these compiled models in memory
# === 3. Run remaining tasks in parallel using update() (FAST)
# =========================================================================

# --- 0. CONFIGURATION ---

# Options: "tte", "binary", "count", "continuous"
ENDPOINT_ID <- "continuous" # <-- Set to "tte" for survival modeling
RESULTS_DIR <- "Results"

# Set unique seed based on script name for independent parallel execution
script_name <- "Continuous_OVAT_5"
RNGkind('Mersenne-Twister')
set.seed(0)

# --- Source functions ---
source('functions.R') # Assuming this contains fit_brms_model etc.

# --- Define the "Example 2" hierarchical prior ---
message("Defining custom full hierarchical prior (N(mu, sigma^2))...")
stanvars_full_hierarchical <- brms::stanvar(
  scode = "  real mu_pred;\n  real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); \n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); \n",
    block = "model"
  )
prior_full_hierarchical <- brms::set_prior("normal(mu_pred, sigma_pred)")

PRIOR_SPECIFICATIONS <- list(
  full_hierarchical = list(prior = prior_full_hierarchical, stanvars = stanvars_full_hierarchical)
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
library(survival)
library(bonsaiforest2) # <-- Ensure this is loaded before processing

# Set threads to 1 for each parallel process
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

# --- 2. DEFINE ENDPOINT-SPECIFIC PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))
endpoint_params <- switch(
  ENDPOINT_ID,
  "tte" = list(folder = "TTE", resp_formula = "Surv(tt_pfs, ev_pfs) ~ arm", resptype = "survival"),
  "binary" = list(folder = "binary", resp_formula = "y ~ arm", resptype = "binary"),
  "count" = list(folder = "Count", resp_formula = "y ~ arm", resptype = "count"),
  "continuous" = list(folder = "continuous", resp_formula = "y ~ arm", resptype = "continuous"),
  stop("Invalid ENDPOINT_ID")
)

# --- 3. LOAD AND PROCESS INPUT DATA ---
cat("Loading all scenario data...\n")
scenarios_dir <- file.path(endpoint_params$folder, "Scenarios")
if (!dir.exists(scenarios_dir)) stop("Data directory not found")

scenario_files <- grep("scenario\\d\\.rds$", dir(scenarios_dir, full.names = TRUE), value = TRUE)
scenarios <- lapply(scenario_files, readRDS)
all_datasets_processed <- lapply(seq_along(scenarios), function(i) {
  setNames(scenarios[[i]], paste(i, 1:length(scenarios[[i]]), sep = "_"))
})
flat_named_list <- unlist(all_datasets_processed, recursive = FALSE)

flat_named_list <- flat_named_list[4001:5000]
rm(scenarios, all_datasets_processed); gc()

cat(sprintf("Successfully loaded %d datasets.\n", length(flat_named_list)))

# ... (Lines 1 to 110: Configuration, Library Loading, Parameter Definition) ...

# --- 4. Apply knot-fixing for TTE models ---

if (ENDPOINT_ID == "tte") {
  cat("\nApplying knot-fixing for fast brms::update()...\n")

  # Extract the time variable name from the response formula
  time_var_name <- stringr::str_trim(stringr::str_match(endpoint_params$resp_formula, "Surv\\((.*?),(.*?)\\)")[, 2])

  # 1. Calculate the fixed knots based on the FIRST dataset
  # Use ::: to forcibly access the unexported function for calculation
  FIXED_KNOTS <- bonsaiforest2:::.calculate_bhaz_knots(flat_named_list[[1]][[time_var_name]])

  # 2. Define the fake function that always returns the fixed knots
  fake_knots_fun <- function(time_data) { return(FIXED_KNOTS) }

  # 3. Patch the function within the bonsaiforest2 namespace

  # Get the package environment
  bonsai_env <- as.environment("package:bonsaiforest2")

  # Use ::: to get the original function before overwriting it
  # We must store it so we can restore it later
  original_knots_fun <- bonsaiforest2:::.calculate_bhaz_knots

  try({
    # Unlock binding in the package environment
    # Use the name without ':::'. R knows where to find it in the environment.
    unlockBinding(".calculate_bhaz_knots", bonsai_env)

    # Assign the fake function to the package environment
    assign(".calculate_bhaz_knots", fake_knots_fun, bonsai_env)

    # Lock binding again to protect it
    lockBinding(".calculate_bhaz_knots", bonsai_env)
  }, silent = FALSE) # Set to silent=FALSE temporarily for debugging patches, if needed

  # Ensure the original function is restored when this *entire script* exits
  on.exit({
    cat("Restoring original bonsaiforest2::.calculate_bhaz_knots function...\n")
    try({
      unlockBinding(".calculate_bhaz_knots", bonsai_env)
      # Restore the original function stored earlier
      assign(".calculate_bhaz_knots", original_knots_fun, bonsai_env)
      lockBinding(".calculate_bhaz_knots", bonsai_env)
    }, silent = TRUE)
  }, add = TRUE)

  cat("Knots fixed. Patch applied to bonsaiforest2 namespace.\n")
}



# --- 5. DEFINE PARAMETERS AND GRID (renumbered from 4) ---
subgr_vars <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")
num_cores <- 96

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  model_type = c("OVAT"),
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)
total_tasks <- nrow(task_grid)

# --- 6. LOGGING SETUP (renumbered from 5) ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "ovat_5", prior_names_str, sep = "_")
results_dir_path <- file.path(endpoint_params$folder, RESULTS_DIR)
dir.create(results_dir_path, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_path, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path(results_dir_path, paste0(base_file_name, ".rds"))

cat(sprintf("Results: %s\nLog: %s\n", results_file, log_file))


# =========================================================================
# === PART A: SERIAL PRE-COMPILATION (TASK 1)
# =========================================================================
# ... (rest of Part A remains the same)

GLOBAL_COMPILED_MODELS <- list()
results_list <- list()

cat("\n--- STARTING PRE-COMPILATION (Task 1 Serially) ---\n")
sim_id_1 <- task_grid$sim_id[1]
df_1 <- flat_named_list[[sim_id_1]]
prior_spec_1 <- PRIOR_SPECIFICATIONS[[task_grid$prior_name[1]]]

results_task1_list <- list()

for (covariate in subgr_vars) {
  cat(sprintf("  Compiling model for covariate: %s ... ", covariate))

  # 1. Prepare Data
  prepared_run <- bonsaiforest2::prepare_formula_model(
    data = df_1,
    response_formula_str = endpoint_params$resp_formula,
    response_type = endpoint_params$resptype,
    unshrunk_prognostic_formula_str = paste("~", covariate),
    shrunk_predictive_formula_str = paste("~ arm:", covariate)
  )

  # 2. Fit (Compiles here)
  fit_obj <- fit_brms_model(
    prepared_model = prepared_run,
    predictive_effect_priors = list(shrunk = prior_spec_1$prior),
    stanvars = prior_spec_1$stanvars,
    chains = 4, iter = 2000, warmup = 1000, cores = 1, backend = "cmdstanr"
  )

  # 3. Save the FITTED object to the global list
  GLOBAL_COMPILED_MODELS[[covariate]] <- fit_obj

  # 4. Get Summary for Task 1 results
  summ <- bonsaiforest2::summary_subgroup_effects(
    fit_obj, df_1, "arm", endpoint_params$resptype, subgroup_vars = 'auto'
  )
  results_task1_list[[covariate]] <- summ$estimates
  cat("Done.\n")
}

# Combine results for Task 1
result_task_1_df <- bind_rows(results_task1_list) %>%
  mutate(model_type = "OVAT", prior_name = task_grid$prior_name[1], .before = 1)

# Store in main results list
results_list[[1]] <- result_task_1_df

cat("--- Pre-compilation complete. Models stored in memory. ---\n")


# =========================================================================
# === PART B: PARALLEL UPDATE WORKER
# =========================================================================
# The worker function can remain the same since the patch is now active
# in the bonsaiforest2 namespace when the workers are forked.
run_update_task <- function(i, task_grid, flat_named_list, GLOBAL_COMPILED_MODELS, 
                            endpoint_params, subgr_vars, log_file) {
  # 1. Setup Environment
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(threads = 1)

  library(bonsaiforest2)
  library(purrr)
  library(dplyr)
  if(file.exists("functions.R")) source("functions.R")

  # 2. Get Task Info
  current_task <- task_grid[i, ]
  sim_id <- current_task$sim_id
  df <- flat_named_list[[sim_id]]

  task_desc <- sprintf("SimID: %s", sim_id)
  start_t <- Sys.time()
  # Optional: Less verbose logging for speed
  # cat(sprintf("[%s] START: %s\n", format(start_t, "%H:%M:%S"), task_desc), file = log_file, append = TRUE)

  out <- tryCatch({

    cov_results <- list()

    for (covariate in subgr_vars) {

      # --- KEY STEP: Retrieve the pre-compiled model ---
      base_model <- GLOBAL_COMPILED_MODELS[[covariate]]

      # Prepare NEW data
      prepared_run <- bonsaiforest2::prepare_formula_model(
        data = df,
        response_formula_str = endpoint_params$resp_formula,
        response_type = endpoint_params$resptype,
        unshrunk_prognostic_formula_str = paste("~", covariate),
        shrunk_predictive_formula_str = paste("~ arm:", covariate)
      )

      # --- KEY STEP: UPDATE() ---
      # The patch ensures recompile = FALSE is effectively maintained
      updated_fit <- update(
        base_model,
        newdata = prepared_run$data,
        recompile = FALSE,
        refresh = 0 # Silence Stan output
      )

      # Summarize
      summ <- bonsaiforest2::summary_subgroup_effects(
        updated_fit, df, "arm", endpoint_params$resptype, subgroup_vars = 'auto'
      )
      cov_results[[covariate]] <- summ$estimates
    }

    bind_rows(cov_results)

  }, error = function(e) {
    msg <- gsub("\\n", " ", e$message)
    cat(sprintf("ERROR %s: %s\n", task_desc, msg), file = log_file, append = TRUE)
    tibble(error = msg)
  })

  # Add identifiers
  if (is.data.frame(out) && !"error" %in% names(out)) {
    out$model_type <- current_task$model_type
    out$prior_name <- current_task$prior_name
    out <- out %>% select(model_type, prior_name, everything())
  }

  # Log completion
  dur <- difftime(Sys.time(), start_t, units = "mins")
  cat(sprintf("[%s] DONE: %s (%.2f m)\n", format(Sys.time(), "%H:%M:%S"), task_desc, dur),
      file = log_file, append = TRUE)

  return(out)
}


# =========================================================================
# === PART C: RUN PARALLEL (Task 2 to N)
# =========================================================================

if (total_tasks > 1) {
  cat(sprintf("\n--- Starting parallel run for remaining %d tasks on %d cores ---\n", total_tasks - 1, num_cores))

  results_tasks_2_to_N <- mclapply(
    X = 2:total_tasks,
    FUN = run_update_task,
    task_grid = task_grid,
    flat_named_list = flat_named_list,
    GLOBAL_COMPILED_MODELS = GLOBAL_COMPILED_MODELS,
    endpoint_params = endpoint_params,
    subgr_vars = subgr_vars,
    log_file = log_file,
    mc.cores = num_cores,
    mc.preschedule = FALSE
  )

  results_list <- c(results_list, results_tasks_2_to_N)
} else {
  cat("No parallel tasks needed.\n")
}

# =========================================================================
# === PART D: COMBINE & SAVE
# ========================================================================

cat("\nCombining results...\n")
final_results_df <- bind_rows(results_list)

if (nrow(final_results_df) > 0) {
  rows_per_task <- sapply(results_list, function(res) if(is.data.frame(res)) nrow(res) else 1)
  task_indices <- rep(1:total_tasks, times = rows_per_task)

  final_results_df <- final_results_df %>%
    mutate(
      sim_id_full = task_grid$sim_id[task_indices],
      scenario_id = as.integer(stringr::str_extract(sim_id_full, "^\\d+")),
      replication_id = as.integer(stringr::str_extract(sim_id_full, "\\d+$")),
      .before = 1
    ) %>%
    select(-sim_id_full)
}

saveRDS(final_results_df, file = results_file)
cat(sprintf("Final results saved to %s\n", results_file))

