#!/usr/bin/env Rscript
# =========================================================================
# === ONE-WAY SHRINKAGE MODEL (OVAT) SCRIPT - CONTINUOUS ENDPOINT
# === SCENARIOS 1-3 (CONSOLIDATED)
# === PARALLEL SIMULATION SCRIPT (OVAT + mclapply)
# === PRIOR: HN_phi_sigma_plan (Half-Normal with phi = sigma_plan)
# =========================================================================

# --- 0. CONFIGURATION ---
ENDPOINT_ID <- "continuous"
RESULTS_DIR <- "Results"
RNGkind('Mersenne-Twister')
set.seed(0)
source('functions.R')

sigma_plan <- 1.2
phi <- sigma_plan

PRIOR_SPECIFICATIONS <- list(
  HN_phi_sigma_plan = list(
    name = "Half-Normal(sigma_plan)",
    prior = brms::set_prior(paste0("normal(0, ", phi, ")"), class = "sd")
  )
)

# --- 1. LOAD LIBRARIES AND SET UP ---
message("--- Loading Libraries and Functions ---")
library(parallel); library(parallelly); library(dplyr); library(tibble); library(purrr)
library(checkmate); library(cmdstanr); library(RhpcBLASctl); library(stringr); library(tidyr); library(brms); library(survival); library(bonsaiforest2)

RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1); data.table::setDTthreads(threads = 1)

# --- 2. DEFINE ENDPOINT-SPECIFIC PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))
endpoint_params <- switch(ENDPOINT_ID,
  "tte" = list(folder = "TTE", resp_formula = "Surv(tt_pfs, ev_pfs) ~ trt", resptype = "survival"),
  "binary" = list(folder = "Binary", resp_formula = "y ~ trt", resptype = "binary"),
  "count" = list(folder = "Count", resp_formula = "y ~ trt", resptype = "count"),
  "continuous" = list(folder = "Continuous", resp_formula = "Y ~ trt", resptype = "continuous"),
  stop("Invalid ENDPOINT_ID")
)

# --- 3. LOAD AND PROCESS INPUT DATA (SCENARIOS 1-3) ---
cat("Loading scenario data for scenarios 1-3...\n")
scenarios_dir <- file.path(endpoint_params$folder, "Scenarios")
if (!dir.exists(scenarios_dir)) stop("Data directory not found")

scenarios_to_run <- as.character(1:3)
scenarios <- list()
for (scen in scenarios_to_run) {
  scen_file <- file.path(scenarios_dir, paste0("scenario", scen, ".rds"))
  if (!file.exists(scen_file)) stop(paste("Scenario", scen, "not found"))
  scenarios[[scen]] <- readRDS(scen_file)
}

all_datasets_processed <- lapply(seq_along(scenarios), function(i) {
  setNames(scenarios[[i]], paste(i, 1:length(scenarios[[i]]), sep = "_"))
})
flat_named_list <- unlist(all_datasets_processed, recursive = FALSE)
rm(scenarios, all_datasets_processed); gc()

cat(sprintf("Successfully loaded %d datasets from scenarios 1-3.\n", length(flat_named_list)))

# --- 4. DEFINE PARAMETERS AND GRID ---
subgr_vars <- c("X1", "X2", "X3", "X4", "X8", "X11cat", "X14cat", "X17cat")
num_cores <- 96

task_grid <- expand.grid(
  sim_id = names(flat_named_list),
  covariate = subgr_vars,
  model_type = c("OVAT"),
  prior_name = names(PRIOR_SPECIFICATIONS),
  stringsAsFactors = FALSE
)
total_tasks <- nrow(task_grid)

# --- 5. LOGGING SETUP ---
prior_names_str <- paste(names(PRIOR_SPECIFICATIONS), collapse = "_")
base_file_name <- paste(ENDPOINT_ID, "ovat_oneway", prior_names_str, sep = "_")
results_dir_path <- file.path(endpoint_params$folder, RESULTS_DIR)
dir.create(results_dir_path, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(results_dir_path, paste0(base_file_name, ".log"))
results_file <- file.path(results_dir_path, paste0(base_file_name, ".rds"))

if (file.exists(log_file)) file.remove(log_file)

cat(sprintf("Total tasks: %d\nResults: %s\nLog: %s\n", total_tasks, results_file, log_file))

# --- 6. RUN FUNCTION FOR ONE-WAY MODELS ---

run_ovat_task <- function(i) {
  current_task <- task_grid[i, ]
  sim_id <- current_task$sim_id
  covariate <- current_task$covariate
  df <- flat_named_list[[sim_id]]

  start_t <- Sys.time()
  task_desc <- sprintf("Sim %s | Cov: %s | Prior: %s", sim_id, covariate, current_task$prior_name)

  cat(sprintf("[%s] START: %s\n", format(start_t, "%H:%M:%S"), task_desc), file = log_file, append = TRUE)

  out <- tryCatch({

    task_prior <- PRIOR_SPECIFICATIONS[[current_task$prior_name]]$prior

    fit <- bonsaiforest2::run_brms_analysis(
      data = df,
      response_type = endpoint_params$resptype,
      response_formula = as.formula(endpoint_params$resp_formula),
      unshrunk_terms_formula = as.formula(paste("~", covariate)),
      shrunk_prognostic_formula = NULL,
      shrunk_predictive_formula = as.formula(paste("~ (0 + trt ||", covariate, ")")),
      intercept_prior = "normal(0, 10)",
      unshrunk_prior = NULL,
      shrunk_predictive_prior = task_prior,
      chains = 4, iter = 2000, warmup = 1000, cores = 4, refresh = 0, backend = "cmdstanr"
    )

    summ <- bonsaiforest2::summary_subgroup_effects(fit)
    est_df <- summ$estimates

    rhat_vals <- brms::rhat(fit)
    max_rhat <- ifelse(!is.null(rhat_vals), max(rhat_vals, na.rm = TRUE), NA)
    rhat_gt_1_05 <- ifelse(!is.null(rhat_vals), sum(rhat_vals > 1.05, na.rm = TRUE), NA)

    n_divergent <- sum(fit$sampler_diagnostics[, , "divergent__"])
    neff_ratio <- brms::neff_ratio(fit)
    mean_neff_ratio <- ifelse(!is.null(neff_ratio), mean(neff_ratio, na.rm = TRUE), NA)

    est_df$covariate_tested <- covariate
    est_df$max_rhat <- max_rhat
    est_df$n_rhats_gt_1_05 <- rhat_gt_1_05
    est_df$n_divergent <- n_divergent
    est_df$mean_neff_ratio <- mean_neff_ratio

    est_df

  }, error = function(e) {
    msg <- gsub("\\n", " ", e$message)
    cat(sprintf("ERROR %s: %s\n", task_desc, msg), file = log_file, append = TRUE)
    tibble(error = msg)
  })

  if (is.data.frame(out) && !"error" %in% names(out)) {
    out$model_type <- current_task$model_type
    out$prior_name <- current_task$prior_name
    out$task_idx <- i
    out <- out %>% dplyr::select(task_idx, model_type, prior_name, everything())
  }

  dur <- difftime(Sys.time(), start_t, units = "mins")
  cat(sprintf("[%s] DONE: %s (%.2f m)\n", format(Sys.time(), "%H:%M:%S"), task_desc, dur), file = log_file, append = TRUE)
  return(out)
}

# --- 7. RUN PARALLEL BY COVARIATE ---

all_results_by_cov <- list()

for (cov_idx in seq_along(subgr_vars)) {
  current_cov <- subgr_vars[cov_idx]
  cov_task_indices <- which(task_grid$covariate == current_cov)
  n_cov_tasks <- length(cov_task_indices)
  
  if (n_cov_tasks == 0) next
  
  cat(sprintf("\n--- Processing Covariate %s (%d/%d) | %d tasks ---\n", 
              current_cov, cov_idx, length(subgr_vars), n_cov_tasks))
  
  cat(sprintf("Pre-compiling %s with first dataset...\n", current_cov))
  first_task_idx <- cov_task_indices[1]
  precompile_result <- run_ovat_task(first_task_idx)
  
  if (n_cov_tasks > 1) {
    cat(sprintf("Running remaining %d tasks in parallel...\n", n_cov_tasks - 1))
    remaining_task_indices <- cov_task_indices[-1]
    
    results_parallel <- mclapply(
      X = remaining_task_indices,
      FUN = run_ovat_task,
      mc.cores = num_cores,
      mc.preschedule = FALSE
    )
    
    cov_results <- c(list(precompile_result), results_parallel)
  } else {
    cov_results <- list(precompile_result)
  }
  
  all_results_by_cov[[paste(current_cov, sep = "_")]] <- cov_results
}

# --- 8. COMBINE & SAVE ---

results_list <- unlist(all_results_by_cov, recursive = FALSE)

cat("\nCombining results...\n")
final_results_df <- bind_rows(results_list)

if (nrow(final_results_df) > 0) {
  final_results_df <- final_results_df %>%
    left_join(
      task_grid %>% 
        mutate(task_idx = row_number()) %>%
        select(task_idx, sim_id),
      by = "task_idx"
    ) %>%
    mutate(
      scenario_id = as.integer(stringr::str_extract(sim_id, "^\\d+")),
      replication_id = as.integer(stringr::str_extract(sim_id, "\\d+$")),
      .before = 1
    ) %>%
    dplyr::select(-c(task_idx, sim_id))
}

saveRDS(final_results_df, file = results_file)
cat(sprintf("Final results saved to %s\n", results_file))
