
## --- 1. LIBRARIES & SETUP ---

library(survival)
library(tidyverse) # For piping and dplyr
library(checkmate) # For assertions
library(MASS)      # For mvrnorm (in simul_covariates)
library(broom)     # For tidy()
library(bonsaiforest)

# This script assumes ALL simulation functions are loaded:
source('functions.R')

RNGkind('Mersenne-Twister')
set.seed(0) # For reproducibility


## --- 2. PARAMETERS AND CONTAINERS ---

# Select endpoints to run
all_endpoints <- c("tte")
all_scenarios <- as.character(1:4)

# TTE Model formulas
base_model_tte <- Surv(tt_pfs, ev_pfs) ~ arm
subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10

# Subgroup Definitions
# Standard subgroups (x_1 to x_10) - 25 total subgroups
# Used for TTE endpoint
subgroup_names <- c(
  "x_1.a", "x_1.b", "x_2.a", "x_2.b", "x_3.a", "x_3.b", "x_4.a",
  "x_4.b", "x_4.c", "x_5.a", "x_5.b", "x_5.c", "x_5.d", "x_6.a",
  "x_6.b", "x_7.a", "x_7.b", "x_8.a", "x_8.b", "x_8.c", "x_9.a",
  "x_9.b", "x_10.a", "x_10.b", "x_10.c"
)

# For non-TTE endpoints, use all subgroups (x_1 to x_10) to match data generation
# Same as subgroup_names, but kept as separate variable for clarity
subgroup_names_short <- subgroup_names

# Simulation Settings
n_repetitions <- 10       # Number of large simulations to average over
inflation_factor <- 1000  # Multiplier for N (1000 * 1000 = 1 Million patients)
t_quantile <- 0.95        # For TTE AHR calculation

# Helper to parse "x_1.a" into variable "x_1" and level "a"
get_subgroup_filter <- function(sg_name) {
  parts <- strsplit(sg_name, "\\.")[[1]]
  list(var = parts[1], lvl = parts[2])
}

# Helper to fit marginal models for non-TTE endpoints
fit_marginal_model <- function(data, endpoint) {
  if (endpoint == "binary") {
    # Log-Odds Ratio
    fit <- glm(y ~ arm, data = data, family = binomial(link = "logit"))
    return(coef(fit)["arm1"])
  } else if (endpoint == "count") {
    # Log-Rate Ratio - use Negative Binomial to match data generation and estimators
    fit <- tryCatch({
      MASS::glm.nb(y ~ arm, data = data)
    }, error = function(e) {
      # Fallback to Poisson if glm.nb fails (rare, but possible with small subgroups)
      warning("glm.nb failed, using Poisson GLM as fallback")
      glm(y ~ arm, data = data, family = poisson(link = "log"))
    })
    return(coef(fit)["arm1"])
  } else if (endpoint == "continuous") {
    # Mean Difference
    fit <- lm(y ~ arm, data = data)
    return(coef(fit)["arm1"])
  }
}



for (ep in all_endpoints) {

  # ========================================================
  # BRANCH A: TIME-TO-EVENT (Existing Logic)
  # ========================================================
  if (ep == "tte") {

    message(paste("🚀 Calculating TTE Truth (Empirical) for", n_repetitions, "reps..."))

    # Initialize result frames
    true_overall_results <- init_data_frame(all_scenarios, c("HR", "AHR", "median_C", "median_I"))
    true_subgroup_hr <- init_data_frame(all_scenarios, subgroup_names)
    true_subgroup_ahr <- init_data_frame(all_scenarios, subgroup_names)

    for (scenario in all_scenarios) {
      message(paste("   ...Scenario", scenario))

      # Matrices to store repetitions
      all_overall_results <- matrix(nrow = 4, ncol = n_repetitions)
      all_subgroups_log_hr <- matrix(nrow = length(subgroup_names), ncol = n_repetitions)
      all_subgroups_log_ahr <- matrix(nrow = length(subgroup_names), ncol = n_repetitions)

      for (j in seq_len(n_repetitions)) {
        # Simulate TTE Data
        sim_data <- simul_scenario(
          scenario = scenario, n_datasets = 1, inflation_factor = inflation_factor, add_uncensored_pfs = TRUE
        )[[1]]
        sim_data$arm <- factor(sim_data$arm)

        # 1. Overall
        all_overall_results[1, j] <- coxph(base_model_tte, data = sim_data)$coef
        all_overall_results[2, j] <- log(ahr_from_km("tt_pfs", "arm", sim_data, "ev_pfs", t_quantile))
        all_overall_results[3, j] <- median(sim_data$tt_pfs_uncens[sim_data$arm == "0"])
        all_overall_results[4, j] <- median(sim_data$tt_pfs_uncens[sim_data$arm == "1"])

        # 2. Subgroups (Using existing stacking logic)
        stacked_data <- generate_stacked_data(base_model_tte, subgroup_model, sim_data, resptype = "survival")

        naive_estimates <- stacked_data %>%
          group_by(subgroup) %>%
          do(tidy(coxph(Surv(time, status) ~ arm, data = .))) %>%
          filter(term == "arm1")

        all_subgroups_log_hr[, j] <- naive_estimates$estimate[match(subgroup_names, naive_estimates$subgroup)]

        # AHR Calculation
        for (k in seq_along(subgroup_names)) {
          ds_k <- subset(stacked_data, subgroup == subgroup_names[k])
          if (nrow(ds_k) > 0 && length(unique(ds_k$arm)) == 2) {
            all_subgroups_log_ahr[k, j] <- log(ahr_from_km("time", "arm", ds_k, "status", t_quantile))
          }
        }

      }

      # Averaging results over repetitions
      true_overall_results[scenario, "HR"] <- exp(mean(all_overall_results[1, ], na.rm = TRUE))
      true_overall_results[scenario, "AHR"] <- exp(mean(all_overall_results[2, ], na.rm = TRUE))
      true_overall_results[scenario, "median_C"] <- mean(all_overall_results[3, ], na.rm = TRUE)
      true_overall_results[scenario, "median_I"] <- mean(all_overall_results[4, ], na.rm = TRUE)

      true_subgroup_hr[scenario, ] <- exp(apply(all_subgroups_log_hr, 1, mean, na.rm = TRUE))
      true_subgroup_ahr[scenario, ] <- exp(apply(all_subgroups_log_ahr, 1, mean, na.rm = TRUE))

    }

    # Save TTE Results
    simul_parameter <- list(
      true_overall_results = true_overall_results,
      true_subgroup_hr = true_subgroup_hr,
      true_subgroup_ahr = true_subgroup_ahr
    )
    save_file_tte <- file.path("TTE", "Scenarios", "truth.RData")
    save(simul_parameter, file = save_file_tte)
    message(paste("Saved TTE truth to:", save_file_tte))

}
}


