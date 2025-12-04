# -----------------------------------------------------------------
# RUN GLOBAL MODEL ANALYSIS (bonsaiforest2) FOR SUNFISH TRIAL
# Prior: r2d2_strong
# -----------------------------------------------------------------

# --- 0. CONFIGURATION ---
ENDPOINT_ID <- "continuous"
RESULTS_DIR <- "Results"

PRIOR_SPECIFICATIONS <- list(
  r2d2_strong = list(prior = "R2D2(mean_R2 = 0.5, prec_R2 = 4, cons_D2 = 2)", stanvars = NULL)
)

# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(parallel); library(parallelly); library(dplyr); library(tibble)
library(purrr); library(checkmate); library(cmdstanr); library(RhpcBLASctl)
library(stringr); library(tidyr); library(bonsaiforest2); library(survival)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
data.table::setDTthreads(threads = 1)

source("../functions.R")

endpoint_params <- list(folder = "SUNFISH", resp_formula = "y ~ arm", resptype = "continuous")

# --- 3. LOAD AND PROCESS INPUT DATA ---
scenarios_to_run <- as.character(1:6)  # All 6 scenarios for SUNFISH
scenarios_list <- list()

for (scen in scenarios_to_run) {
  scen_file <- file.path("Scenarios", paste0("SUNFISH_Scenario_", scen, ".rds"))
  if (!file.exists(scen_file)) stop(paste("File not found:", scen_file))
  scenarios_list[[scen]] <- readRDS(scen_file)
}

all_datasets_processed <- lapply(seq_along(scenarios_list), function(i) {
  setNames(scenarios_list[[i]], paste(i, 1:length(scenarios_list[[i]]), sep = "_"))
})
flat_named_list <- unlist(all_datasets_processed, recursive = FALSE)
rm(scenarios_list, all_datasets_processed); gc()

# --- 4. DEFINE PARAMETERS AND CREATE TASK LIST ---
RNGkind("L'Ecuyer-CMRG"); set.seed(42)
covariate_set <- c("x_1", "x_2", "x_3", "x_4", "x_5")
num_cores <- availableCores()

task_grid <- expand.grid(sim_id = names(flat_named_list), model_type = c("Global"),
  prior_name = names(PRIOR_SPECIFICATIONS), stringsAsFactors = FALSE)

prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste("SUNFISH", "global", prior_names_str, sep = "_")
log_file <- file.path(RESULTS_DIR, paste0(base_file_name, ".log"))
if (file.exists(log_file)) file.remove(log_file)
results_file <- file.path(RESULTS_DIR, paste0(base_file_name, ".rds"))

# --- 6. RUN IN PARALLEL & COMBINE ---
run_single_task <- function(i) {
  current_task <- task_grid[i, ]
  run_global_model_task(sim_id = current_task$sim_id, model_type = current_task$model_type,
    prior_name = current_task$prior_name, all_data = flat_named_list, subgr_vars = covariate_set,
    log_path = log_file, prior_spec = PRIOR_SPECIFICATIONS[[current_task$prior_name]],
    endpoint_params = endpoint_params)
}

result_task_1 <- list(run_single_task(1))

if (nrow(task_grid) > 1) {
  results_tasks_2_to_N <- mclapply(X = 2:nrow(task_grid), FUN = run_single_task,
    mc.cores = num_cores, mc.preschedule = FALSE)
} else {
  results_tasks_2_to_N <- list()
}

results_list <- c(result_task_1, results_tasks_2_to_N)
final_results_df <- bind_rows(results_list)

# --- 7. FINAL FORMATTING AND SAVING ---
if (nrow(final_results_df) > 0) {
  rows_per_task <- sapply(results_list, function(res) if (is.data.frame(res)) nrow(res) else 1)
  task_indices <- rep(1:nrow(task_grid), times = rows_per_task)
  final_results_df <- final_results_df %>%
    mutate(sim_id_full = task_grid$sim_id[task_indices],
      scenario_id = as.integer(stringr::str_extract(sim_id_full, "^\\d+")),
      replication_id = as.integer(stringr::str_extract(sim_id_full, "\\d+$")), .before = 1) %>%
    select(-sim_id_full)
}

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(object = final_results_df, file = results_file)

