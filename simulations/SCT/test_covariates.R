#!/usr/bin/env Rscript
# Quick test script to verify SCT covariate generation and correlations

library(MASS)

# Source the functions
source("functions.R")

# Generate a test dataset
set.seed(123)
test_data <- simul_covariates_sct(n = 180, arm_factor = TRUE)

cat("=== SCT Covariate Generation Test ===\n\n")

# 1. Check sample size and structure
cat("1. Dataset dimensions:\n")
cat(sprintf("   - N = %d (expected 180)\n", nrow(test_data)))
cat(sprintf("   - Variables = %d\n", ncol(test_data)))
cat("\n")

# 2. Check treatment allocation
cat("2. Treatment allocation (2:1 ratio):\n")
print(table(test_data$arm))
cat(sprintf("   - Ratio = %.2f:1\n", sum(test_data$arm == 1) / sum(test_data$arm == 0)))
cat("\n")

# 3. Check covariate distributions
cat("3. Covariate distributions:\n\n")

cat("   Age Group (x_1):\n")
print(prop.table(table(test_data$x_1)))
cat("   Expected: 2-5y=0.31, 6-11y=0.32, 12-17y=0.25, 18-25y=0.12\n\n")

cat("   SMA Type (x_2):\n")
print(prop.table(table(test_data$x_2)))
cat("   Expected: Type2=0.71, Type3=0.29\n\n")

cat("   SMN2 Copy Number (x_3):\n")
print(prop.table(table(test_data$x_3)))
cat("   Expected: 2_copies=0.03, 3_copies=0.87, 4_copies=0.10\n\n")

cat("   Disease Severity (x_4):\n")
print(prop.table(table(test_data$x_4)))
cat("   Expected: roughly 25%, 50%, 25% based on quartiles\n\n")

# 4. Check baseline MFM32
cat("4. Baseline MFM32 Score:\n")
cat(sprintf("   - Mean = %.2f (expected 46.11)\n", mean(test_data$baseline_mfm)))
cat(sprintf("   - SD = %.2f (expected 11.46)\n", sd(test_data$baseline_mfm)))
cat(sprintf("   - Range = [%.2f, %.2f] (should be within [0, 96])\n", 
            min(test_data$baseline_mfm), max(test_data$baseline_mfm)))
cat("\n")

# 5. Check relationships
cat("5. Key Relationships (approximate due to categorization):\n\n")

# Convert factors to numeric for correlation
test_data$age_num <- as.numeric(test_data$x_1)
test_data$type_num <- as.numeric(test_data$x_2)
test_data$smn2_num <- as.numeric(test_data$x_3)
test_data$severity_num <- as.numeric(test_data$x_4)

cat("   Correlation: Age Group × SMN2 Copy Number\n")
cat(sprintf("     r = %.3f (expected ~0.3 in latent space)\n", 
            cor(test_data$age_num, test_data$smn2_num)))

cat("   Correlation: SMN2 Copy Number × Baseline MFM32\n")
cat(sprintf("     r = %.3f (expected ~0.35 in latent space)\n", 
            cor(test_data$smn2_num, test_data$baseline_mfm)))

cat("   Correlation: Age Group × Baseline MFM32\n")
cat(sprintf("     r = %.3f (expected ~-0.15 in latent space)\n", 
            cor(test_data$age_num, test_data$baseline_mfm)))

cat("\n   Note: Correlations are attenuated by categorization\n")
cat("\n")

# 6. Generate a scenario 2 dataset to verify coefficients
cat("6. Testing Scenario 2 coefficient structure:\n")

params <- .get_model_parameters_sct()
coefs <- .get_scenario_coefs_sct("2", params)

cat(sprintf("   - Number of non-zero coefficients: %d\n", sum(coefs != 0)))
cat("   - Key coefficients:\n")
cat(sprintf("     * arm1 (2-5y effect): %.2f\n", coefs["arm1"]))
cat(sprintf("     * x_16-11y_arm (6-11y modifier): %.2f\n", coefs["x_16-11y_arm"]))
cat(sprintf("     * x_32_copies_arm (2 SMN2 copies interaction): %.2f\n", 
            coefs["x_32_copies_arm"]))

cat("\n=== Test Complete ===\n")
