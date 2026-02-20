# =========================================================================
# === NAIVE POPULATION ANALYSIS FOR CONTINUOUS OUTCOMES
# ===
# === This script:
# === 1. Runs the naive population (overall) analysis for continuous outcomes
# === 2. Processes all scenarios
# === 3. Saves results to the Results/ folder
# =========================================================================

# --- 0. CONFIGURATION ---
ENDPOINT_ID <- "continuous"
RESULTS_DIR <- "Results"

RNGkind('Mersenne-Twister')
set.seed(0)

# --- 1. LOAD LIBRARIES AND FUNCTIONS ---
message("--- Loading Libraries and Functions ---")
library(checkmate)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(parallel)

# --- 2. DEFINE ENDPOINT PARAMETERS ---
message(paste("--- Configuring for Endpoint:", ENDPOINT_ID, "---"))

# Determine the script's directory to make paths work from anywhere
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- getwd()
}

endpoint_params <- list(
  folder = script_dir,
  resp = "Y",
  resptype = "continuous"
)

# --- 3. SET UP RESULTS DIRECTORY ---
results_dir_path <- file.path(endpoint_params$folder, RESULTS_DIR)
dir.create(results_dir_path, recursive = TRUE, showWarnings = FALSE)

base_file_name <- paste(ENDPOINT_ID, "naivepop", sep = "_")
results_file <- file.path(results_dir_path, paste0(base_file_name, ".rds"))
log_file <- file.path(results_dir_path, paste0(base_file_name, ".log"))

if (file.exists(log_file)) file.remove(log_file)

cat(sprintf("Results will be saved to: %s\n", results_file))

# --- 4. LOAD SCENARIO DATA ---
message("Loading all scenario datasets...")
scenarios_dir <- file.path(endpoint_params$folder, "Scenarios")

all_scenarios <- list()
for (scen_no in 1:12) {
  scenario_file <- file.path(scenarios_dir, sprintf("scenario%d.rds", scen_no))
  if (file.exists(scenario_file)) {
    all_scenarios[[scen_no]] <- readRDS(scenario_file)
  } else {
    warning(sprintf("Scenario %d not found\n", scen_no))
  }
}

cat(sprintf("Loaded %d scenarios\n", length(all_scenarios)))

# --- 5. DEFINE NAIVE POPULATION ANALYSIS FUNCTION ---

run_naive_population_analysis <- function(scenario_data, scenario_no) {
  
  cat(sprintf("Processing scenario %d: %s simulations\n", 
              scenario_no, n_distinct(scenario_data$sim_id)))
  
  results <- scenario_data %>%
    group_by(sim_id) %>%
    nest() %>%
    mutate(
      fit = map(data, ~lm(Y ~ trt, data = .x)),
      tidy_result = map(fit, ~broom::tidy(.x) %>% 
                        filter(term == "trt") %>%
                        select(term, estimate, std.error, statistic, p.value)),
      glance_result = map(fit, broom::glance)
    ) %>%
    unnest(c(tidy_result, glance_result), keep_empty = TRUE) %>%
    select(sim_id, estimate, std.error, statistic, p.value, r.squared, adj.r.squared) %>%
    rename(
      trt_effect = estimate,
      trt_se = std.error,
      trt_t_stat = statistic,
      trt_p_value = p.value
    ) %>%
    mutate(
      scenario_no = scenario_no,
      scenario = unique(scenario_data$scenario)
    ) %>%
    select(scenario_no, scenario, sim_id, everything())
  
  return(results)
}

# --- 6. RUN ANALYSIS FOR ALL SCENARIOS ---

all_results <- list()

for (scen_no in seq_along(all_scenarios)) {
  
  cat(sprintf("\n--- Scenario %d ---\n", scen_no))
  
  start_t <- Sys.time()
  result <- run_naive_population_analysis(all_scenarios[[scen_no]], scen_no)
  end_t <- Sys.time()
  
  elapsed <- difftime(end_t, start_t, units = "secs")
  cat(sprintf("✓ Completed in %0.1f seconds\n", elapsed))
  
  all_results[[scen_no]] <- result
}

# --- 7. COMBINE AND SAVE ---

combined_results <- bind_rows(all_results)

saveRDS(combined_results, file = results_file)
cat(sprintf("\n✓ All results saved to %s\n", results_file))
cat(sprintf("Processed %d scenarios with %d total simulations\n", 
            n_distinct(combined_results$scenario_no),
            nrow(combined_results)))
