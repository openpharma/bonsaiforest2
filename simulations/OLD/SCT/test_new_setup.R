# Test script to verify new SCT simulation setup
# 4 subgroup variables, 0.6 correlation, specific percentages, SD=6.0

library(MASS)
source('functions.R')

# Set seed for reproducibility
RNGkind('Mersenne-Twister')
set.seed(42)

cat("Testing new SCT simulation setup\n")
cat("=================================\n\n")

# Test 1: Generate covariates
cat("1. Testing covariate generation...\n")
covariates <- simul_covariates_sct(n = 180, arm_factor = TRUE)

cat("   Sample size:", nrow(covariates), "\n")
cat("   Variables:", paste(names(covariates), collapse = ", "), "\n\n")

# Check distributions
cat("2. Checking distributions:\n")
cat("   x_1 (4 levels, 25% each):\n")
print(prop.table(table(covariates$x_1)))

cat("\n   x_2 (3 levels: 50%, 25%, 25%):\n")
print(prop.table(table(covariates$x_2)))

cat("\n   x_3 (2 levels: 33%, 66%):\n")
print(prop.table(table(covariates$x_3)))

cat("\n   x_4 (3 levels: 20%, 40%, 40%):\n")
print(prop.table(table(covariates$x_4)))

cat("\n   Treatment arm (2:1 ratio):\n")
print(table(covariates$arm))

# Test 2: Check model parameters
cat("\n3. Testing model parameters...\n")
params <- .get_model_parameters_sct()
cat("   SD:", params$sd, "(should be 6.0)\n")
cat("   Base effect:", params$base_effect, "\n")
cat("   Intercept:", params$intercept, "\n")

# Test 3: Generate a single dataset for each scenario
cat("\n4. Testing data generation for each scenario:\n")
for (scen in 1:6) {
  cat("   Scenario", scen, "... ")
  
  tryCatch({
    dat <- .simul_sct_single(
      n = 180,
      model_params = params,
      coefs = .get_scenario_coefs_sct(as.character(scen), params),
      scenario = as.character(scen)
    )
    
    cat("✓ (N=", nrow(dat), ", mean y=", round(mean(dat$y), 2), 
        ", sd y=", round(sd(dat$y), 2), ")\n", sep = "")
  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })
}

# Test 4: Generate multiple datasets
cat("\n5. Testing multiple dataset generation (scenario 1, n=10):\n")
datasets <- simul_sct_data(scenario = "1", n_datasets = 10)
cat("   Generated", length(datasets), "datasets\n")
cat("   First dataset dimensions:", nrow(datasets[[1]]), "x", ncol(datasets[[1]]), "\n")

# Test 5: Check for all subgroup x arm combinations
cat("\n6. Checking subgroup × arm combinations (first dataset):\n")
dat <- datasets[[1]]
for (var in c("x_1", "x_2", "x_3", "x_4")) {
  cross_tab <- table(dat[[var]], dat$arm)
  has_zeros <- any(cross_tab == 0)
  cat("   ", var, ":", ifelse(has_zeros, "✗ Has empty cells", "✓ All combinations present"), "\n")
}

cat("\n✓ All tests completed!\n")
