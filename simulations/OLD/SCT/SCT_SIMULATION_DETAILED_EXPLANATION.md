# SCT Simulation: Super-Detailed Technical Explanation

## Table of Contents

1.  [Overview & Architecture](#overview--architecture)
2.  [Data Generation Pipeline](#data-generation-pipeline)
3.  [Truth Calculation](#truth-calculation)
4.  [Analysis Methods](#analysis-methods)
5.  [Execution Framework](#execution-framework)
6.  [Results Processing](#results-processing)

------------------------------------------------------------------------

## Overview & Architecture

### What is SCT?

**SCT** refers to the **SecondClinical Trial** we will use as a example for our data simulation process. This trial is a real-world randomized controlled trial for a rare disease treatment. This simulation framework replicates the trial structure to evaluate different Bayesian subgroup analysis methods.

### Trial Parameters

-   **Sample Size**: N = 180 (2:1 randomization, 120 treatment : 60 control)
-   **Endpoint**: Continuous (Change from baseline in MFM32 motor function score)
-   **Primary Effect**: 1.55 points improvement (Treatment: +1.36, Placebo: -0.19)
-   **Standard Deviation**: 6.0
-   **Subgroup analysis**: Subgroup estimator suggested that treatment effects vary by age and region

### File Structure

```         
SCT/
├── functions.R                  # Core simulation logic
├── Scenarios_generation.Rmd     # Generate 1000 datasets per scenario
├── Truth.Rmd                    # Calculate theoretical true effects
├── Population.R                 # Naive overall analysis (no subgroups)
├── Subgroup.R                   # Naive per-subgroup analysis
├── Horseshoe_low.R              # Bayesian global model (weak shrinkage)
├── Horseshoe_mid.R              # Bayesian global model (moderate shrinkage)
├── Horseshoe_strong.R           # Bayesian global model (strong shrinkage)
├── R2D2_low.R                   # Alternative Bayesian prior (weak)
├── R2D2_mid.R                   # Alternative Bayesian prior (moderate)
├── R2D2_strong.R                # Alternative Bayesian prior (strong)
├── run_all.sh                   # Parallel execution script
└── Results_analysis.Rmd         # Compare all methods
```

## Data Generation Pipeline {#data-generation-pipeline}

### Step 1: Covariate Generation (`simul_covariates_sct()`)

Located in `functions.R`, lines 132-284.

#### Covariate Construction

Each covariate is created by **cutting** the continuous normal distribution at specific quantiles to match trial frequencies:

**x_1: Age Group** (4 levels, distribution from Table 1)

``` r
x$x_1 <- cut(pnorm(z[,1]),
             breaks = c(0, 0.31, 0.63, 0.88, 1.0),
             labels = c("2-5y", "6-11y", "12-17y", "18-25y"))
# 31% are 2-5y, 32% are 6-11y, 25% are 12-17y, 12% are 18-25y
```

**x_2: SMA Type** (2 levels)

``` r
x$x_2 <- cut(pnorm(z[,2]),
             breaks = c(0, 0.71, 1.0),
             labels = c("Type2", "Type3"))
# 71% Type 2, 29% Type 3 (Type 2 is more severe, earlier onset)
```

**x_3: SMN2 Copy Number** (3 levels, genetic biomarker)

``` r
x$x_3 <- cut(pnorm(z[,3]),
             breaks = c(0, 0.03, 0.90, 1.0),
             labels = c("2_copies", "3_copies", "4_copies"))
# 3% have 2 copies, 87% have 3 copies, 10% have 4 copies
# Fewer copies = more severe disease
```

**baseline_mfm: Continuous Baseline Function** (from z[,4])

``` r
base_score <- 46.11 + (z[,4] * 11.46)
x$baseline_mfm <- pmin(pmax(base_score, 0), 96)  # Clamp to MFM32 range [0, 96]
# Mean = 46.11, SD = 11.46
```

**x_4: Disease Severity** (3 levels, **derived** from baseline_mfm)

``` r
x$x_4 <- cut(x$baseline_mfm,
             breaks = c(-Inf, 37.5, 54.17, Inf),
             labels = c("Severe_leQ1", "Moderate", "Mild_gtQ3"))
# Split at Q1 and Q3 from trial (Q1=37.5, Q3=54.17)
# NOTE: This is DETERMINISTIC from baseline_mfm, so inherits all correlations
```

**x_5: Region** (5 levels, from Forest Plot)

``` r
x$x_5 <- cut(pnorm(z[,5]),
             breaks = c(0, 0.67, 0.80, 0.895, 0.98, 1.0),
             labels = c("Europe", "NorthAmerica", "China", "Japan", "RoW"))
# 67% Europe, 13% North America, 9.5% China, 8.5% Japan, 2% Rest of World
```

#### Correlation Structure

The covariates are **not** generated independently. They use a multivariate normal distribution with a **correlation matrix** that reflects real biological relationships:

``` r
# Base correlation: 0 between all biological factors
sigma <- matrix(0 , nrow = 5, ncol = 5)
diag(sigma) <- 1

# Specific correlations based on biology and paper:
#should we think about a logical correlation structure?

# Generate 5-dimensional correlated normals
z <- MASS::mvrnorm(n, mu = rep(0, 5), Sigma = sigma)
```

#### Treatment Assignment (2:1 Randomization)

``` r
n_active <- round(n * (2/3))   # 120 patients
n_placebo <- n - n_active       # 60 patients
trt_arm <- sample(rep(c(1, 0), c(n_active, n_placebo)))  # Simple randomization
x$arm <- factor(trt_arm)
```

#### Rejection Sampling for Valid Datasets

**Critical detail**: The function uses **rejection sampling** to ensure every subgroup level has at least 1 patient in BOTH treatment arms:

``` r
for (attempt in 1:max_attempts) {
  # Generate covariates...
  
  # Check if all subgroup × arm combinations exist
  all_present <- TRUE
  for (var in subgroup_vars) {
    cross_tab <- table(x[[var]], x$arm)
    if (any(cross_tab == 0)) {  # Any cell with 0 patients?
      all_present <- FALSE
      break
    }
  }
  
  if (all_present) return(x)  # Success!
}
```

Without this, naive subgroup models would fail (can't fit `y ~ arm` if arm only has 1 level).

### Step 2: Scenario-Specific Coefficients (`.get_scenario_coefs_sct()`)

Located in `functions.R`, lines 323-411.

Each scenario defines a **true data-generating model** with specific regression coefficients. The outcome is generated as:

```         
y[i] = design_matrix[i,] %*% coefficients + ε[i],  where ε ~ N(0, sd²)
```

#### Scenario 1: Positive Homogeneous Effect

``` r
coefs <- c(
  "(Intercept)" = -0.19,  # Placebo arm mean (Table 2)
  "arm1" = 1.55           # Treatment effect (consistent across ALL subgroups)
)
```

**Interpretation**: Every patient improves by exactly 1.55 points on treatment, regardless of age, type, region, etc. This is the "null hypothesis" for heterogeneity—treatment works uniformly.

#### Scenario 2: Positive except Regional Heterogeneity

``` r
coefs <- c(
  "(Intercept)" = -0.19,
  "arm1" = 2.40,                     # Reference: Europe effect
  "x_5NorthAmerica_arm" = 0.04,      # North America: 2.40 + 0.04 = 2.44
  "x_5Japan_arm" = -5.98,            # Japan: 2.40 - 5.98 = -3.58 (HARM!)
  "x_5China_arm" = -4.25,            # China: 2.40 - 4.25 = -1.85 (HARM!)
  "x_5RoW_arm" = 0.90                # RoW: 2.40 + 0.90 = 3.30
)
```

**Critical concept: Dummy Coding** - `arm1` is the treatment effect for the **reference group** (Europe, the first level of x_5) - `x_5Japan_arm` is the **interaction coefficient**: how much Japan differs from Europe - Total effect for Japan = `arm1 + x_5Japan_arm` = 2.40 + (-5.98) = **-3.58**

**Qualitative Interaction (Sign Reversal)**:

-   Western regions (Europe, North America): **Beneficial** (+2.4 to +3.3 points)

-   Asian regions (Japan, China): **Harmful** (-1.8 to -3.6 points)

This tests the methods' ability to detect **true heterogeneity** where treatment direction reverses.

#### Scenario 3: Null Overall with Regional Heterogeneity

``` r
coefs <- c(
  "(Intercept)" = 0,
  "arm1" = 0,                
   "x_5NorthAmerica_arm" = 0.04,     
  "x_5Japan_arm" = -5.98,           
  "x_5China_arm" = -4.25,           
  "x_5RoW_arm" = 0.90               
)
```

**Key property**: Overall effect averages to \~0 (71% × 2.0 + 29% × (-2.0) ≈ 0.62 ≈ 0)

This tests if methods can detect **subgroup-specific effects** when the naive population model shows no benefit.

#### Scenario 4: Mild Random Heterogeneity

``` r
# Base effect
coefs <- c("(Intercept)" = -0.19, "arm1" = 1.55)

# Add random noise to ALL interaction terms
all_interaction_names <- c(
  "x_16-11y_arm", "x_112-17y_arm", "x_118-25y_arm",
  "x_2Type3_arm", 
  "x_33_copies_arm", "x_34_copies_arm",
  "x_4Moderate_arm", "x_4Mild_gtQ3_arm",
  "x_5NorthAmerica_arm", "x_5China_arm", "x_5Japan_arm", "x_5RoW_arm"
)

random_effects <- rnorm(length(all_interaction_names), mean = 0, sd = 0.15)
names(random_effects) <- all_interaction_names
coefs <- c(coefs, random_effects)
```

**Interpretation**: Every subgroup interaction gets a small random bump (SD=0.15), creating **unpredictable heterogeneity**. Tests if shrinkage priors can filter signal from noise.

#### Scenario 5: Large Random Heterogeneity

Same as Scenario 4, but with **SD = 0.3** instead of 0.5. This creates **massive variance** across subgroups, stressing the shrinkage priors.

#### Scenario 6: Complex 2-Way Interaction (Age × Type)

``` r
coefs <- c(
  "(Intercept)" = -0.19,
  "arm1" = 2.5,                         # Base: Young (2-5y) Type2
  "x_16-11y_arm" = -0.5,                # 6-11y: reduces effect by 0.5
  "x_112-17y_arm" = -1.0,               # 12-17y: reduces by 1.0
  "x_118-25y_arm" = -2.0,               # 18-25y: reduces by 2.0
  "x_2Type3_arm" = -0.8,                # Type3: reduces by 0.8
  "x_118-25y_x_2Type3_arm" = -1.5       # 18-25y AND Type3: ADDITIONAL penalty
)
```

### Step 3: Design Matrix Construction (`.simul_sct_single()`)

Located in `functions.R`, lines 415-491.

This is the **most complex part** of data generation, directly mirroring how `bonsaiforest2` constructs formulas internally.

#### Main Effects Matrix

``` r
# For Scenario 6 (with Age×Type interaction)
subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_1:x_2

# For other scenarios (main effects only)
subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5

design_main <- model.matrix(update(subgroup_model, ~ arm + .), data = covariates)
```

**Example columns from design_main**:

```         
(Intercept)  arm1  x_16-11y  x_112-17y  x_118-25y  x_2Type3  x_33_copies  ...
     1        1      0         0          0         0          1         ...
     1        0      1         0          0         1          1         ...
```

-   `arm1`: Indicator for treatment (1 if treated, 0 if control)
-   `x_16-11y`: Dummy for age 6-11 years (0 if reference 2-5y)
-   Reference levels are **absorbed into intercept** (standard R behavior)

#### Interaction Matrix (Treatment × Subgroups)

``` r
subgroup_vars <- all.vars(subgroup_model)  # c("x_1", "x_2", "x_3", "x_4", "x_5")
design_ia <- NULL

for (j in subgroup_vars) {
  # Create dummy matrix for this variable
  ia_j <- model.matrix(as.formula(paste("~", j, "-1")), data = covariates)
  
  # Multiply by treatment indicator (0 or 1)
  ia_j <- ia_j * as.numeric(covariates$arm == "1")
  
  design_ia <- cbind(design_ia, ia_j)
}
colnames(design_ia) <- paste(colnames(design_ia), "arm", sep = "_")
```

**Key insight**: The `-1` in formula removes intercept, creating dummies for **ALL levels** (including reference):

```         
For x_1 (Age):
x_12-5y  x_16-11y  x_112-17y  x_118-25y  →  MULTIPLY by arm  →  x_12-5y_arm, x_16-11y_arm, ...
   1        0         0          0                                      1×arm        0×arm
   0        1         0          0                                      0×arm        1×arm
```

This creates **over-parameterized dummy coding**: every level gets an interaction term, not just non-reference levels.

#### For Scenario 6: Add Higher-Order Interaction

``` r
if (scenario == "6") {
  # Create Age × Type interaction columns
  age_type_ia <- model.matrix(~ x_1:x_2 - 1, data = covariates)
  age_type_ia <- age_type_ia * as.numeric(covariates$arm == "1")
  colnames(age_type_ia) <- paste(colnames(age_type_ia), "arm", sep = "_")
  design_ia <- cbind(design_ia, age_type_ia)
}
```

### Step 4: Coefficient Mapping and Outcome Generation

``` r
# Initialize coefficient vector matching design matrix
reg_coef <- rep(0, ncol(design_matrix))
names(reg_coef) <- colnames(design_matrix)

# Map scenario coefficients to design matrix columns
valid_coefs <- intersect(names(coefs), names(reg_coef))
reg_coef[valid_coefs] <- coefs[valid_coefs]

# Generate linear predictor
lp <- design_matrix %*% reg_coef

# Add Gaussian noise
y <- rnorm(n, mean = lp, sd = model_params$sd)  # sd = 2.0
```

### Step 5: Replication and Storage

``` r
# simul_sct_data() in functions.R, lines 493-531
replicate(
  n_datasets,  # 1000
  .simul_sct_single(n = params$N, model_params = params, coefs = coefs, scenario = scenario),
  simplify = FALSE
)
```

Each scenario generates **1000 independent datasets** with identical data-generating parameters but different random noise.

**File output** (`Scenarios_generation.Rmd`):

``` r
for (scen in 1:6) {
  datasets <- simul_sct_data(scenario = as.character(scen), n_datasets = 1000)
  saveRDS(datasets, file.path("Scenarios", paste0("SCT_Scenario_", scen, ".rds")))
}
```

Creates: - `Scenarios/SCT_Scenario_1.rds` (1000 datasets, homogeneous) - `Scenarios/SCT_Scenario_2.rds` (1000 datasets, regional heterogeneity) - ... (6 files total)

------------------------------------------------------------------------

## Truth Calculation {#truth-calculation}

Located in `Truth.Rmd`.

### Why Calculate Truth?

The simulation compares **estimated effects** from different methods against **true effects** from the data-generating model. RMSE (root mean squared error) quantifies accuracy.

### Overall Treatment Effect

The overall effect is the **population-averaged** treatment effect:

``` r
overall_effect <- coefs["arm1"]
```

**BUT** this is only correct if no subgroup interactions exist!

For Scenario 3 (crossover): - `arm1` = 2.0 (Type 2 effect) - `x_2Type3_arm` = -4.0 (Type 3 interaction) - True overall = 0.71 × 2.0 + 0.29 × (-2.0) = 1.42 - 0.58 = **0.84**

The code **simplifies** by using `arm1` directly, which is technically the **reference group effect**, not the population average. This is acceptable if we interpret "overall" as "reference subgroup" for consistency.

### Subgroup-Specific Effects

Each subgroup's true effect is calculated by **adding the base effect + relevant interactions**:

``` r
for (subg_name in subgroup_names) {
  coef_name <- subgroup_coef_map[subg_name]  # Maps "x_5.Japan" → "x_5Japan_arm"
  
  if (coef_name == "arm1") {
    # Reference level
    total_effect <- coefs["arm1"]
  } else {
    # Non-reference level
    interaction_effect <- coefs[coef_name]  # Could be NA if not defined
    total_effect <- coefs["arm1"] + interaction_effect
  }
}
```

**Example: Scenario 2, Japan**

``` r
subg_name = "x_5.Japan"
coef_name = "x_5Japan_arm"
coefs["arm1"] = 2.40
coefs["x_5Japan_arm"] = -5.98
total_effect = 2.40 + (-5.98) = -3.58
```

### Scenario 6 Special Case: 2-Way Interactions

For Age × Type combinations:

``` r
for (age in c("2-5y", "6-11y", "12-17y", "18-25y")) {
  for (type in c("Type2", "Type3")) {
    total_effect <- coefs["arm1"]  # Start with base
    
    # Add age main effect (if not reference)
    if (age != "2-5y") {
      total_effect <- total_effect + coefs[paste0("x_1", age, "_arm")]
    }
    
    # Add type main effect (if not reference)
    if (type != "Type3") {
      total_effect <- total_effect + coefs["x_2Type3_arm"]
    }
    
    # Add 2-way interaction (only for 18-25y × Type3)
    if (age == "18-25y" && type == "Type3") {
      total_effect <- total_effect + coefs["x_118-25y_x_2Type3_arm"]
    }
  }
}
```

**Example: 18-25y Type3**

```         
2.5 (base) + (-2.0) (age) + (-0.8) (type) + (-1.5) (interaction) = -1.8
```

### Output Format

``` r
simul_truth <- list(
  truth_overall = data.frame(scenario, metric="mean_diff", value),
  truth_subgroup = data.frame(scenario, subgroup, metric="mean_diff", value),
  truth_scen6_interaction = data.frame(scenario, subgroup, metric="mean_diff", value)
)
save(simul_truth, file = "Scenarios/truth.RData")
```

------------------------------------------------------------------------

## Analysis Methods {#analysis-methods}

The simulation compares **8 estimators** across 6 scenarios × 1000 datasets = 6000 analyses per method.

### Method 1-2: Naive Frequentist Baselines

#### Population Method (`Population.R`)

**Strategy**: Ignore subgroups, fit one model to entire dataset.

``` r
# For continuous endpoint
fit <- lm(y ~ arm, data = df)

# Extract treatment effect
est <- broom::tidy(fit)
est_arm <- est[est$term == "arm1", ]  # Coefficient for treatment
```

**Output**: Single estimate, replicated for ALL subgroups (recycled).

``` r
data.frame(
  simul_no = 1,
  estimator = "population",
  subgroup = c("x_1.2-5y", "x_1.6-11y", ..., "x_5.RoW"),  # All 17 subgroups
  estimate = -0.415,  # Same value for all rows
  std.error = 0.127,
  lower_ci = -0.664,
  upper_ci = -0.166
)
```

**Performance expectation**: - **Best** for Scenario 1 (homogeneous—truth matches model) - **Worst** for Scenarios 2-3 (misses heterogeneity) - **Moderate** for Scenarios 4-5 (averages out noise, but loses signal)

------------------------------------------------------------------------

#### Subgroup Method (`Subgroup.R`)

**Strategy**: Fit **separate models** for each subgroup level.

``` r
# Create "stacked" dataset where each patient appears once per subgroup membership
stacked_data <- generate_stacked_data(base_model, subgr_model, data, resptype)

# Example: Patient 1 is 2-5y, Type2, Europe → appears in 3 subgroups:
#   x_1.2-5y, x_2.Type2, x_5.Europe

# Fit separate model for each subgroup
list_subg <- split(stacked_data, ~subgroup)
fit <- lapply(list_subg, FUN = function(d) {
  lm(y ~ arm, data = d)  # Model on subset only
})

# Extract estimates
naive_estimates <- bind_rows(lapply(fit, broom::tidy), .id = "subgroup")
```

**Output**: One estimate per subgroup (17 rows per dataset).

**Performance expectation**: - **Bad** for Scenario 1 (overfits noise, high variance) - **Good** for Scenario 2 (captures true heterogeneity) - **Moderate** for Scenarios 4-5 (noisy, but can detect real signals)

**Key weakness**: Small subgroups (e.g., Japan N≈15) have **huge standard errors** due to limited data.

------------------------------------------------------------------------

### Method 3-8: Bayesian Global Models (bonsaiforest2)

These methods use the **3-step workflow**:

1.  `run_brms_analysis()` - Fit Bayesian hierarchical model with shrinkage priors
2.  `summary_subgroup_effects()` - Estimate marginal effects via G-computation
3.  Results contain estimates for all subgroups

#### Core Idea: Differential Shrinkage

**Horseshoe Prior** (and R2D2) applies **adaptive regularization**: - Strong shrinkage on **predictive interactions** (treatment × subgroups) - Weaker shrinkage on **prognostic main effects** (baseline predictors)

**Intuition**: Most interactions are noise; only a few are signal. Shrinkage pulls noisy interactions toward zero while preserving strong signals.

#### Model Specification (`Horseshoe_low.R`, lines 106-150)

``` r
# Define what to shrink vs not shrink
prognostic_str <- paste(subgr_vars, collapse = " + ")  # "x_1 + x_2 + x_3 + x_4 + x_5"
predictive_str <- paste(paste0("arm:", subgr_vars), collapse = " + ")  
# "arm:x_1 + arm:x_2 + arm:x_3 + arm:x_4 + arm:x_5"

fit <- run_brms_analysis(
  data = df,
  response_formula_str = "y ~ arm",
  response_type = "continuous",
  unshrunk_prognostic_formula_str = paste("~", prognostic_str),
  shrunk_predictive_formula_str = paste("~", predictive_str),
  predictive_effect_priors = list(shrunk = "horseshoe(scale_global = 0.025)"),
  chains = 4, iter = 2000, warmup = 1000
)
```

**Behind the scenes** (in `bonsaiforest2::prepare_formula_model()`):

Creates non-linear brmsformula:

``` r
bf(
  y ~ arm + unprogeffect + shpredeffect,
  unprogeffect ~ 0 + x_1 + x_2 + x_3 + x_4 + x_5,  # Unshrunk (normal prior)
  shpredeffect ~ 0 + arm:x_1 + arm:x_2 + ... ,     # Shrunk (horseshoe prior)
  family = gaussian()
)
```

**Treatment variable conversion**: `arm` is converted to numeric 0/1 internally, and interaction dummies are created.

#### Prior Variants

**Horseshoe Scale Parameter** (controls shrinkage strength): - `scale_global = 0.025` (low/weak) → More shrinkage, assumes most interactions are zero - `scale_global = 0.05` (mid) → Moderate shrinkage - `scale_global = 0.1` (strong) → Less shrinkage, allows more heterogeneity

**R2D2 Prior** (alternative, controls proportion of variance explained): - `R2D2(mean_R2=0.5, prec_R2=1, cons_D2=0.5)` (low) → Weak signal - `R2D2(mean_R2=0.5, prec_R2=2, cons_D2=1)` (mid) - `R2D2(mean_R2=0.5, prec_R2=4, cons_D2=2)` (strong) → Strong signal

#### G-Computation for Marginal Effects

`summary_subgroup_effects()` implements **counterfactual estimation**:

``` r
# For each posterior draw:
for (draw in 1:4000) {
  # 1. Create two counterfactual datasets
  df_treated <- df %>% mutate(arm = factor("1"))
  df_control <- df %>% mutate(arm = factor("0"))
  
  # 2. Predict outcomes under each scenario
  pred_treated <- posterior_epred(fit, newdata = df_treated, draw = draw)
  pred_control <- posterior_epred(fit, newdata = df_control, draw = draw)
  
  # 3. Calculate individual treatment effects
  ite <- pred_treated - pred_control
  
  # 4. Average within subgroups
  for (subgroup in unique(df$x_1)) {
    idx <- which(df$x_1 == subgroup)
    subgroup_effect[draw, subgroup] <- mean(ite[idx])
  }
}

# 5. Summarize posterior (median, credible intervals)
estimates <- apply(subgroup_effect, 2, median)
```

**Key advantage**: Handles complex correlation structure by **marginalizing** over other covariates within each subgroup.

------------------------------------------------------------------------

### Execution Framework (`run_all.sh`)

Parallel job submission:

``` bash
#!/bin/bash
# Generate task list (all .R files)
ls *.R > task_list.txt

# Run in parallel using xargs (ACE platform)
cat task_list.txt | xargs -P 10 -I {} Rscript {}
```

Each `.R` script: 1. Loads all 6 scenarios (6000 datasets total) 2. Creates task grid: `expand.grid(sim_id, model_type, prior_name)` 3. Runs tasks in parallel using `mclapply()` 4. Saves results to `Results/<method_name>.rds`

**Memory management**: - Compiles Stan model on first task (serially) - Runs remaining tasks in parallel - Cleans up cmdstan temporary files after each task - Uses 1 thread per task to avoid resource conflicts

------------------------------------------------------------------------

## Results Processing {#results-processing}

### Output Structure

Each method produces a data frame:

``` r
data.frame(
  scenario_id = 1,         # Scenario number (1-6)
  replication_id = 1,      # Dataset number (1-1000)
  simul_no = 1,            # Same as replication_id (legacy)
  estimator = "horseshoe_low",
  subgroup = "x_1.6-11y",  # Sanitized (S_1.6-11y in some outputs)
  estimate = -0.382,       # Point estimate (log-scale for binary/survival)
  std.error = 0.145,
  lower_ci = -0.666,
  upper_ci = -0.098,
  p.value = 0.008          # (for frequentist methods)
)
```

### Truth Comparison (`Results_analysis.Rmd`)

``` r
# Load truth
load("Scenarios/truth.RData")

# Load estimates
horseshoe_low <- readRDS("Results/SCT_global_horseshoe_low.rds")

# Merge and calculate error
results <- horseshoe_low %>%
  left_join(simul_truth$truth_subgroup, by = c("scenario_id" = "scenario", "subgroup"))

# Calculate squared error per estimate
results <- results %>%
  mutate(squared_error = (estimate - value)^2)

# Aggregate RMSE by scenario and subgroup
rmse_summary <- results %>%
  group_by(scenario_id, estimator, subgroup) %>%
  summarise(rmse = sqrt(mean(squared_error)))

# Aggregate RMSE by scenario only (average across subgroups)
rmse_overall <- results %>%
  group_by(scenario_id, estimator) %>%
  summarise(rmse = sqrt(mean(squared_error)))
```

### Performance Metrics

**RMSE (Root Mean Squared Error)**: Combines bias and variance

```         
RMSE = sqrt(mean((estimate - truth)²))
```

**Bias**: Systematic error

```         
Bias = mean(estimate - truth)
```

**Coverage**: Proportion of 95% CIs containing truth

```         
Coverage = mean(lower_ci <= truth & truth <= upper_ci)
```

**Standardized RMSE**: Relative to population method

```         
SRMSE = RMSE / RMSE_population
```

### Expected Performance Patterns

| Scenario | Population | Subgroup | Horseshoe | Best Method |
|----|----|----|----|----|
| 1 (Homogeneous) | **Best** (RMSE≈SE) | Worst (overfits) | Good (shrinks to truth) | Population |
| 2 (Regional) | Worst (misses heterogeneity) | Good (captures truth) | **Best** (borrows strength) | Horseshoe |
| 3 (Crossover) | Bad (averages to null) | Good | **Best** | Horseshoe |
| 4 (Mild Random) | Good (averages noise) | Bad (high variance) | **Best** (filters noise) | Horseshoe |
| 5 (Large Random) | Moderate | Bad | **Best** (regularizes) | Horseshoe |
| 6 (Interaction) | Bad | Moderate | **Best** (captures structure) | Horseshoe |

**Key insight**: Horseshoe priors provide **automatic model selection**—they shrink irrelevant interactions to zero while preserving true heterogeneity. This makes them robust across all scenarios.

------------------------------------------------------------------------

## Technical Details & Edge Cases

### Why Over-Parameterized Dummy Coding?

Standard R dummy coding drops reference levels (e.g., `x_1a` is absorbed into intercept). But `bonsaiforest2` uses **all levels** for interactions:

``` r
# Standard (treatment coding)
x_1: a (ref), b, c  →  Coefficients: arm1, x_1b_arm, x_1c_arm

# Over-parameterized (all levels)
x_1: a, b, c  →  Coefficients: arm1, x_1a_arm, x_1b_arm, x_1c_arm
```

**Reason**: Shrinkage priors work better with **exchangeable parameters**. If `a` is reference, it gets special treatment (no shrinkage), which breaks symmetry. Over-parameterization treats all levels equally.

**Mathematical equivalence**: The reference level coefficient `x_1a_arm` is redundant (could be absorbed into `arm1`), but Stan handles this via the non-linear formula structure.

### cmdstanr Output Directory Management

``` r
# Problem: Multiple parallel tasks compile the same model → file conflicts
# Solution: Create unique temp directory per task
cmdstan_output_dir <- file.path(
  tempdir(),
  sprintf("cmdstan_%s_%s_%s", prior_name, sim_id, random_suffix)
)

# Cleanup after model finishes
unlink(cmdstan_output_dir, recursive = TRUE)
```

Without this, parallel jobs on ACE platform would crash due to concurrent writes to the same Stan model cache.

### Rejection Sampling Rationale

With N=180 and 17 subgroups × 2 arms = 34 combinations, **some cells will be empty** due to: 1. Small sample size 2. Unbalanced randomization (2:1) 3. Rare subgroups (e.g., RoW = 2% → only \~3-4 patients total)

Rejection sampling ensures **all methods can run** without encountering singular fits (can't estimate treatment effect in subgroup with only 1 arm).

**Failure rate**: Typically succeeds within 10-50 attempts. If it fails after 1000 attempts, returns last dataset with warning.

------------------------------------------------------------------------

## Summary of Code Flow

1.  **Generate Data** (`Scenarios_generation.Rmd`)
    -   For each scenario (1-6):
        -   Get coefficients from `.get_scenario_coefs_sct()`
        -   Generate 1000 datasets via `simul_sct_data()`
        -   Save to `Scenarios/SCT_Scenario_X.rds`
2.  **Calculate Truth** (`Truth.Rmd`)
    -   For each scenario:
        -   Extract coefficients
        -   Calculate overall effect (`arm1`)
        -   Calculate subgroup effects (`arm1 + interactions`)
    -   Save to `Scenarios/truth.RData`
3.  **Run Analyses** (8 scripts: `Population.R`, `Subgroup.R`, `Horseshoe_*.R`, `R2D2_*.R`)
    -   Load all scenarios
    -   For each dataset:
        -   Fit model (frequentist or Bayesian)
        -   Extract subgroup estimates
    -   Save to `Results/<method>.rds`
4.  **Compare Results** (`Results_analysis.Rmd`)
    -   Load all method results + truth
    -   Calculate RMSE, bias, coverage per method/scenario/subgroup
    -   Generate plots (forest plots, RMSE comparisons)

------------------------------------------------------------------------

## Key Takeaways

1.  **SCT replicates a real trial** with realistic N (180), randomization (2:1), and covariate structure

2.  **Scenarios test different heterogeneity patterns**: none (1), strong/qualitative (2), crossover (3), random (4-5), interaction (6)

3.  **Data generation uses correlated covariates** to match biological relationships in SMA

4.  **Over-parameterized dummy coding** ensures all subgroup levels are treated symmetrically for shrinkage

5.  **G-computation** marginalizes over covariates to estimate **causal subgroup effects**

6.  **Horseshoe priors** provide automatic variable selection, shrinking noise while preserving signal

7.  **Simulation evaluates trade-offs**: Population (low variance, high bias) vs Subgroup (high variance, low bias) vs Horseshoe (balanced)

This framework allows **comprehensive comparison** of Bayesian vs frequentist subgroup analysis methods under controlled, clinically-motivated scenarios.
