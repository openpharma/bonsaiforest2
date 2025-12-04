#' SUNFISH Trial Simulation Functions
#'
#' Adapted from the original bonsaiforest2 simulation framework to replicate
#' the SUNFISH Part 2 trial based on the published manuscript and supplementary appendix.
#'
#' TRIAL DESIGN:
#' - N = 180 (120 Risdiplam : 60 Placebo, 2:1 ratio)
#' - Primary endpoint: Change from baseline in MFM32 total score at Month 12
#' - Treatment difference: 1.55 points (Risdiplam 1.36 vs Placebo -0.19)
#' - Standard deviation: 6.0
#'
#' COVARIATES (from Table 1, Figure S3, and Forest Plot):
#' - x_1: Age Group (2-5y, 6-11y, 12-17y, 18-25y)
#' - x_2: SMA Type (Type2, Type3)
#' - x_3: SMN2 Copy Number (2_copies, 3_copies, 4_copies)
#' - x_4: Disease Severity based on Baseline MFM32 quartiles (Severe_leQ1, Moderate, Mild_gtQ3)
#' - x_5: Region (Europe, NorthAmerica, China, Japan, RoW) - from Forest Plot subgroup analysis
#' - baseline_mfm: Baseline MFM32 score (continuous, Mean=46.11, SD=11.46)
#'
#' SCENARIOS:
#'   1. Homogeneous effect (1.55 points, primary analysis result)
#'   2. Age-driven heterogeneity (based on Figure S3 subgroup analysis:
#'      2-5y: 3.14, 6-11y: 1.58, 12-17y: 1.04, 18-25y: -0.65)

# Source the original functions for shared utilities
source('../functions.R')

#' Generate Stacked Data for Subgroup Analysis
#'
#' Creates a long-format dataset where each patient appears once for each
#' subgroup they belong to, with a 'subgroup' indicator.
#'
#' @param base_model Formula for the base model
#' @param subgr_model Formula with subgrouping variables
#' @param data Original data frame
#' @param resptype Type of response ("survival", "binary", "count", "continuous")
#'
#' @return A stacked data frame with subgroup indicators
#'
generate_stacked_data <- function(base_model, subgr_model, data, resptype) {
  
  # Extract variable names from formulas
  base_vars <- all.vars(base_model)
  subgr_vars <- all.vars(subgr_model)
  
  # Get the treatment variable name to preserve factor levels
  trt_var <- if (resptype == "survival") {
    base_vars[length(base_vars)]
  } else {
    base_vars[2]
  }
  
  # Store the original arm factor levels
  arm_levels <- levels(data[[trt_var]])
  
  # Initialize list to hold stacked data for each subgroup
  stacked_list <- list()
  
  # For each subgrouping variable
  for (var in subgr_vars) {
    # Get the levels of this variable
    var_levels <- levels(data[[var]])
    
    # For each level, create a subset of data for that subgroup
    for (level in var_levels) {
      subgroup_data <- data[data[[var]] == level, ]
      
      # Skip subgroups with no observations
      if (nrow(subgroup_data) == 0) next
      
      # Add subgroup identifier
      subgroup_data$subgroup <- paste0(var, ".", level)
      
      # Rename variables for the naive function
      if (resptype == "survival") {
        # For survival: keep original variable names (time, status, arm)
        resp_var <- base_vars[1]  # Assuming first var is response
        subgroup_data$time <- subgroup_data[[resp_var]]
        subgroup_data$arm <- factor(subgroup_data[[trt_var]], levels = arm_levels)
      } else {
        # For binary/count/continuous: standardize to y and arm
        resp_var <- base_vars[1]
        subgroup_data$y <- subgroup_data[[resp_var]]
        subgroup_data$arm <- factor(subgroup_data[[trt_var]], levels = arm_levels)
      }
      
      stacked_list[[paste0(var, ".", level)]] <- subgroup_data
    }
  }
  
  # Combine all subgroups into one stacked data frame
  stacked_data <- do.call(rbind, stacked_list)
  rownames(stacked_data) <- NULL
  
  return(stacked_data)
}

#' Generate Covariates for SUNFISH Part 2 Trial
#'
#' Generates the specific subgroup variables used in the SUNFISH trial:
#' Age Group, SMA Type, SMN2 Copy Number, and Disease Severity.
#'
#' @param n Sample size (Total 180 for SUNFISH Part 2)
#' @param arm_factor Logical. If TRUE, returns arm as a factor
#'
#' @return A data.frame with specific SUNFISH covariates
#'
#' @details
#' Specific SUNFISH covariates based on Table 1, Figure S3, and Forest Plot:
#' - x_1: Age Group (2-5y [31%], 6-11y [32%], 12-17y [25%], 18-25y [12%])
#' - x_2: SMA Type (Type2 [71%], Type3 [29%])
#' - x_3: SMN2 Copy Number (2_copies [3%], 3_copies [87%], 4_copies [10%])
#' - x_4: Disease Severity based on Baseline MFM32 quartiles (Severe_leQ1, Moderate, Mild_gtQ3)
#' - x_5: Region (Europe [67%], NorthAmerica [13%], China [9.5%], Japan [8.5%], RoW [2%])
#' - baseline_mfm: Baseline MFM32 score (continuous, Mean=46.11, SD=11.46)
#'
#' Correlation structure (from SAP simulation documentation):
#' - Biological factors (Age, SMA Type, SMN2, Baseline MFM32) have moderate correlation (0.2)
#' - Reflects known relationships: younger age correlates with Type 2 and lower SMN2 copies
#' - SMN2 copy number correlates with baseline motor function
#'
#' Treatment assignment: 2:1 ratio (Risdiplam=1, Placebo=0)
simul_covariates_sunfish <- function(n = 180, arm_factor = TRUE, max_attempts = 1000) {
  
  # Use rejection sampling to ensure at least 1 patient per subgroup per arm
  for (attempt in 1:max_attempts) {
    
    # 1. Generate underlying continuous normals with biological correlation structure
    # Based on SAP documentation: moderate correlation (0.2) between biological factors
    # This reflects:
    # - Age and SMA Type relationship (Type 2 tends to be diagnosed earlier)
    # - SMN2 copies and disease severity (fewer copies = more severe)
    # - Baseline MFM32 and age (younger patients may have different functional baselines)
    # - Region and demographic factors (geographic enrollment patterns)
    
    sigma <- matrix(0.2, nrow = 5, ncol = 5)  # Expanded to 5x5 for Region
    diag(sigma) <- 1
    
    # Adjust specific correlations based on known biological relationships:
    # - Age and SMN2 copy number: slightly stronger correlation (0.3)
    #   (Type 2 SMA, which is more common in younger patients, often has fewer SMN2 copies)
    sigma[1, 3] <- sigma[3, 1] <- 0.3
    
    # - SMN2 copy number and baseline MFM32: moderate positive correlation (0.35)
    #   (More SMN2 copies generally associated with better motor function)
    sigma[3, 4] <- sigma[4, 3] <- 0.35
    
    # - Age and baseline MFM32: weak negative correlation (-0.15)
    #   (Older patients may have experienced more disease progression over time)
    sigma[1, 4] <- sigma[4, 1] <- -0.15
    
    # Region (x_5) maintains base correlation (0.2) with other factors
    # (Geographic enrollment patterns may relate to age, type, etc.)
    
    z <- MASS::mvrnorm(n, mu = rep(0, 5), Sigma = sigma)  # Now 5 dimensions
    
    x <- data.frame(id = 1:n)
  
  # --- x_1: Age Group (Stratification Factor) ---
  # Distribution from Table 1: 31% (2-5), 32% (6-11), 25% (12-17), 12% (18-25)
  # Cumulative probs: 0.31, 0.63, 0.88, 1.0
  x$x_1 <- cut(pnorm(z[,1]), 
               breaks = c(0, 0.31, 0.63, 0.88, 1.0), 
               labels = c("2-5y", "6-11y", "12-17y", "18-25y"))
  
  # --- x_2: SMA Type ---
  # Distribution from Table 1: ~71% Type 2, ~29% Type 3
  # Correlated with age through z[,2] which has correlation with z[,1]
  x$x_2 <- cut(pnorm(z[,2]), 
               breaks = c(0, 0.71, 1.0), 
               labels = c("Type2", "Type3"))
  
  # --- x_3: SMN2 Copy Number ---
  # Distribution from Table 1: ~3% (2 copies), 87% (3 copies), 10% (4 copies)
  # Cumulative probs: 0.03, 0.90, 1.0
  # Note: Stronger correlation with age (sigma[1,3]=0.3) and baseline MFM32 (sigma[3,4]=0.35)
  # Reflects known biology: fewer copies → younger onset, more severe disease
  x$x_3 <- cut(pnorm(z[,3]), 
               breaks = c(0, 0.03, 0.90, 1.0), 
               labels = c("2_copies", "3_copies", "4_copies"))
  
  # --- Baseline MFM32 Score (Continuous) ---
  # From Article text: Mean 46.11, SD 11.46
  # Correlated with SMN2 copy number (more copies → higher baseline function)
  # and age (older → potentially lower baseline due to progression)
  # We clamp values between 0 and 96 (MFM32 range)
  base_score <- 46.11 + (z[,4] * 11.46)
  x$baseline_mfm <- pmin(pmax(base_score, 0), 96)
  
  # --- x_4: Disease Severity (Categorical Subgroup) ---
  # Based on Appendix p22: Q1=37.5, Q3=54.17
  # Categories: <=Q1 (Severe), >Q1-Q3 (Moderate), >Q3 (Mild)
  # This is deterministically derived from baseline_mfm, so inherits all correlations
  x$x_4 <- cut(x$baseline_mfm,
               breaks = c(-Inf, 37.5, 54.17, Inf),
               labels = c("Severe_leQ1", "Moderate", "Mild_gtQ3"))
  
  # --- x_5: Region (NEW - from Forest Plot Subgroup Analysis) ---
  # Distribution based on N counts in Forest Plot (Risdiplam arm N=115):
  # Europe: 77/115 ≈ 67%
  # North America: 15/115 ≈ 13%
  # China: 11/115 ≈ 9.5%
  # Japan: 10/115 ≈ 8.5%
  # RoW (Rest of World): 2/115 ≈ 2%
  # Cumulative probabilities: 0.67, 0.80, 0.895, 0.98, 1.0
  x$x_5 <- cut(pnorm(z[,5]),
               breaks = c(0, 0.67, 0.80, 0.895, 0.98, 1.0),
               labels = c("Europe", "NorthAmerica", "China", "Japan", "RoW"))
  
    # 2. Treatment Assignment (2:1 Ratio - Risdiplam:Placebo)
    # Risdiplam = 1, Placebo = 0
    n_active <- round(n * (2/3))
    n_placebo <- n - n_active
    # Simple randomization (stratification by age was used in the trial, 
    # but simple random shuffling is sufficient for simulation at this scale)
    trt_arm <- sample(rep(c(1, 0), c(n_active, n_placebo)))
    
    x$arm <- if (arm_factor) factor(trt_arm) else trt_arm
    
    # 3. Check if all subgroup × arm combinations have at least 1 patient
    # This ensures that naive subgroup analysis can fit models for all subgroups
    subgroup_vars <- c("x_1", "x_2", "x_3", "x_4", "x_5")
    
    all_present <- TRUE
    for (var in subgroup_vars) {
      # Check if each level of each subgroup variable has both arms
      cross_tab <- table(x[[var]], x$arm)
      if (any(cross_tab == 0)) {
        all_present <- FALSE
        break
      }
    }
    
    # If all combinations present, return the data
    if (all_present) {
      return(x)
    }
    
    # Otherwise, continue to next attempt
    if (attempt %% 100 == 0) {
      message(sprintf("Attempt %d: Regenerating covariates to ensure all subgroup×arm combinations...", attempt))
    }
  }
  
  # If we exhaust all attempts, issue a warning and return the last attempt
  warning(sprintf("Could not generate data with all subgroup×arm combinations after %d attempts. Returning last attempt.", max_attempts))
  return(x)
}

#' Get Model Parameters for SUNFISH Trial
#'
#' Returns calibrated parameters for the MFM32 change from baseline endpoint.
#'
#' @return List with endpoint, N, intercept, base_effect, and sd
#'
#' @details
#' Calibration based on SUNFISH Part 2 results:
#' - Intercept: -0.19 (Mean change in Placebo, Table 2)
#' - Base effect: 1.55 (Treatment difference, Table 2)
#' - SD: 5.0 (Reduced from 6.0 to decrease variance in estimates)
.get_model_parameters_sunfish <- function() {
  list(
    endpoint = "continuous",
    N = 180,             # Total N from SUNFISH Part 2
    intercept = -0.19,   # Placebo mean change (Table 2)
    base_effect = 1.55,  # Treatment difference (Table 2)
    sd = 5.0             # Reduced SD for lower variance
  )
}

#' Get Scenario Coefficients for SUNFISH Simulations
#'
#' Defines the true coefficient values for each of 6 simulation scenarios.
#' Adapts the original 6-scenario framework to SUNFISH trial context.
#'
#' @param scenario Character, "1" through "6"
#' @param model_params List from .get_model_parameters_sunfish()
#'
#' @return Named vector of coefficients
#'
#' @details
#' Scenario 1 (Homogeneous): 
#'   - Consistent 1.55 point improvement across all subgroups
#'   - Replicates the primary efficacy analysis
#'
#' Scenario 2 (Heterogeneous - Region-Driven):
#'   - Based on Forest Plot subgroup analysis showing dramatic regional differences
#'   - Europe/North America: Strong positive effects (~2.40)
#'   - Asian regions (Japan/China): Negative effects (qualitative interaction)
#'   - Tests ability to detect geographic heterogeneity and sign reversal
#'
#' Scenario 3 (Null/Crossover):
#'   - Overall null effect (mean = 0)
#'   - SMA Type 2 improves, Type 3 worsens (crossover)
#'
#' Scenario 4 (Mild Random Heterogeneity):
#'   - Small random noise (SD=0.5) added to treatment effects
#'   - Tests robustness to mild heterogeneity
#'
#' Scenario 5 (Large Random Heterogeneity):
#'   - Large random noise (SD=2.0) added to treatment effects
#'   - Tests robustness to strong heterogeneity
#'
#' Scenario 6 (Interaction):
#'   - Specific Age × SMA Type interaction
#'   - Treatment works differently for young Type 2 vs Type 3 patients
#'
#' Prognostic effects (all scenarios):
#'   - SMA Type 3 naturally progresses slightly better than Type 2 (0.5 points)
#'   - SMN2 copy number has strong prognostic effect
#'
#' Note: Uses dummy coding where arm1 represents reference group
.get_scenario_coefs_sunfish <- function(scenario, model_params) {
  
  # Set seed for reproducibility of random scenarios
  RNGkind('Mersenne-Twister')
  set.seed(as.integer(scenario) * 123)
  
  base_effect <- model_params$base_effect  # 1.55
  intercept <- model_params$intercept      # -0.19
  
  # Base prognostic factors (treatment-independent)
  constant_coefs <- c(
    "(Intercept)" = intercept,
    "x_2Type3" = 0.5,          # Type 3 vs Type 2 prognostic effect
    "x_32_copies" = -1.2,      # 2 SMN2 copies: worse prognosis
    "x_34_copies" = 0.3        # 4 SMN2 copies: better prognosis
  )
  
  # Helper: names for all subgroup interaction terms (including Region)
  # These match design matrix column names
  all_interaction_names <- c(
    "x_16-11y_arm", "x_112-17y_arm", "x_118-25y_arm",     # Age (ref: 2-5y)
    "x_2Type3_arm",                                         # Type (ref: Type2)
    "x_33_copies_arm", "x_34_copies_arm",                  # SMN2 (ref: 2_copies)
    "x_4Moderate_arm", "x_4Mild_gtQ3_arm",                 # Severity (ref: Severe)
    "x_5NorthAmerica_arm", "x_5China_arm", "x_5Japan_arm", "x_5RoW_arm"  # Region (ref: Europe)
  )
  
  coefs <- switch(
    scenario,
    
    # --- Scenario 1: Homogeneous (Primary Trial Result) ---
    "1" = c(
      constant_coefs,
      "arm1" = base_effect  # 1.55 across all subgroups
    ),
    
    # --- Scenario 2: Heterogeneous (Region-Driven, from Forest Plot) ---
    # Dramatic regional heterogeneity with qualitative interaction (sign reversal)
    # Reference: Europe (strongest Western region effect)
    "2" = c(
      constant_coefs,
      # Reference effect (Europe): +2.40 points
      "arm1" = 2.40,
      
      # Regional treatment interactions (relative to Europe)
      # North America (Target: 2.44) -> 2.44 - 2.40 = +0.04
      "x_5NorthAmerica_arm" = 0.04,
      
      # Japan (Target: -3.58) -> -3.58 - 2.40 = -5.98 (QUALITATIVE INTERACTION)
      "x_5Japan_arm" = -5.98,
      
      # China (Target: -1.85) -> -1.85 - 2.40 = -4.25 (QUALITATIVE INTERACTION)
      "x_5China_arm" = -4.25,
      
      # RoW/Rest of World (Target: 3.30) -> 3.30 - 2.40 = +0.90
      "x_5RoW_arm" = 0.90
    ),
    
    # --- Scenario 3: Null/Crossover (Type-Based) ---
    # Overall null, but Type 2 benefits and Type 3 worsens
    "3" = c(
      constant_coefs,
      "arm1" = 2.0,              # Type 2 (reference): +2.0
      "x_2Type3_arm" = -4.0      # Type 3: 2.0 - 4.0 = -2.0 (crossover)
    ),
    
    # --- Scenario 4: Mild Random Heterogeneity ---
    # Add small random noise to all subgroup interactions
    "4" = {
      random_effects <- rnorm(length(all_interaction_names), mean = 0, sd = 0.5)
      names(random_effects) <- all_interaction_names
      c(constant_coefs, "arm1" = base_effect, random_effects)
    },
    
    # --- Scenario 5: Large Random Heterogeneity ---
    # Add large random noise to all subgroup interactions
    "5" = {
      random_effects <- rnorm(length(all_interaction_names), mean = 0, sd = 2.0)
      names(random_effects) <- all_interaction_names
      c(constant_coefs, "arm1" = base_effect, random_effects)
    },
    
    # --- Scenario 6: Interaction (Age × Type) ---
    # Specific interaction: treatment works best in young Type 2 patients
    "6" = c(
      constant_coefs,
      "arm1" = 2.5,                        # Base (2-5y, Type 2)
      "x_16-11y_arm" = -0.5,               # Older age reduces effect
      "x_112-17y_arm" = -1.0,
      "x_118-25y_arm" = -2.0,
      "x_2Type3_arm" = -0.8,               # Type 3 reduces effect
      # Specific interaction: being both old AND Type 3 is particularly bad
      "x_118-25y_x_2Type3_arm" = -1.5     # Additional penalty for 18-25y Type 3
    )
  )
  
  return(coefs)
}

#' Simulate a Single SUNFISH Dataset
#'
#' Internal helper function to generate one dataset with MFM32 outcomes.
#'
#' @param n Sample size
#' @param model_params List from .get_model_parameters_sunfish()
#' @param coefs Named vector of coefficients
#' @param scenario Character scenario ID ("1" to "6") - needed for Scenario 6 interaction
#'
#' @return Data frame with id, arm, covariates (x_1 to x_5, baseline_mfm), and y (MFM32 change)
.simul_sunfish_single <- function(n, model_params, coefs, scenario = "1") {
  
  # 1. Generate SUNFISH Covariates (N=180, 2:1 Ratio, Age/SMA/SMN2/Severity/Region)
  covariates <- simul_covariates_sunfish(n = n, arm_factor = TRUE)
  
  # 2. Build Design Matrix
  # Scenario 6 requires Age × Type interaction with treatment
  # Other scenarios use standard main effects + arm interactions
  if (scenario == "6") {
    # Add specific 2-way interaction for Age × Type × Arm
    subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_1:x_2
  } else {
    # Standard: all main effects (including Region x_5)
    subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5
  }
  
  # Main effects design matrix (arm + covariates)
  design_main <- model.matrix(update(subgroup_model, ~ arm + .), data = covariates)
  
  # Build interaction matrix (arm * all terms in subgroup_model)
  # This creates columns like: x_16-11y_arm, x_2Type3_arm, x_118-25y_x_2Type3_arm, etc.
  subgroup_vars <- all.vars(subgroup_model)
  design_ia <- NULL
  
  for (j in subgroup_vars) {
    ia_j <- model.matrix(
      as.formula(paste("~", j, "-1")), 
      data = covariates
    ) * as.numeric(covariates$arm == "1")  # 1 for treatment, 0 for control
    design_ia <- cbind(design_ia, ia_j)
  }
  colnames(design_ia) <- paste(colnames(design_ia), "arm", sep = "_")
  
  # For Scenario 6, also add the Age:Type interaction with arm
  if (scenario == "6") {
    # Create Age × Type interaction columns
    age_type_ia <- model.matrix(~ x_1:x_2 - 1, data = covariates) * 
                   as.numeric(covariates$arm == "1")
    colnames(age_type_ia) <- paste(colnames(age_type_ia), "arm", sep = "_")
    design_ia <- cbind(design_ia, age_type_ia)
  }
  
  design_matrix <- cbind(design_main, design_ia)
  
  # 3. Apply Coefficients
  reg_coef <- rep(0, ncol(design_matrix))
  names(reg_coef) <- colnames(design_matrix)
  
  # Map our defined coefs to the matrix columns (safety check for name matching)
  valid_coefs <- intersect(names(coefs), names(reg_coef))
  reg_coef[valid_coefs] <- coefs[valid_coefs]
  
  # 4. Generate Outcome (Continuous MFM32 Change from Baseline)
  lp <- design_matrix %*% reg_coef
  y <- rnorm(n, mean = lp, sd = model_params$sd)
  
  # 5. Return data frame (keeping baseline_mfm as a covariate)
  data.frame(covariates, y = y)
}

#' Simulate SUNFISH Study Data for Multiple Datasets
#'
#' Main wrapper function to generate replicated datasets for SUNFISH scenarios.
#'
#' @param scenario Character, "1" through "6" (see .get_scenario_coefs_sunfish for details)
#' @param n_datasets Number of replicate datasets to generate (default 1000)
#'
#' @return List of data frames, each containing one simulated dataset
#'
#' @examples
#' \dontrun{
#' # Generate 1000 datasets for scenario 1 (homogeneous effect)
#' datasets_s1 <- simul_sunfish_data(scenario = "1", n_datasets = 1000)
#' 
#' # Generate 1000 datasets for scenario 2 (age heterogeneity)
#' datasets_s2 <- simul_sunfish_data(scenario = "2", n_datasets = 1000)
#' 
#' # Generate for all 6 scenarios
#' for (s in 1:6) {
#'   data <- simul_sunfish_data(scenario = as.character(s), n_datasets = 1000)
#'   saveRDS(data, paste0("Scenario_", s, ".rds"))
#' }
#' }
simul_sunfish_data <- function(scenario = c("1", "2", "3", "4", "5", "6"), n_datasets = 1000) {
  
  # Match arguments
  scenario <- match.arg(scenario)
  assert_count(n_datasets)
  
  message(paste("Generating SUNFISH Data: Scenario", scenario))
  
  # Get SUNFISH parameters
  params <- .get_model_parameters_sunfish()
  
  # Get coefficients based on scenario
  coefs <- .get_scenario_coefs_sunfish(scenario, params)
  
  # Generate the list of n_datasets
  replicate(
    n_datasets,
    .simul_sunfish_single(
      n = params$N,
      model_params = params,
      coefs = coefs,
      scenario = scenario
    ),
    simplify = FALSE
  )
}
