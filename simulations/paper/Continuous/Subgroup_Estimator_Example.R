#--------------------------------------------------------------------------------------------------
#------ COMPUTE SUBGROUP TREATMENT EFFECT ESTIMATORS FOR SCENARIO 3
#------
#------ This script:
#------ 1. Loads the first simulation from scenario 3
#------ 2. For each subgroup, fits a linear model on subgroup data only
#------ 3. Extracts the treatment effect estimate for each subgroup
#------ 4. Compares estimates against the true subgroup effects
#------ 5. Creates a comparison plot
#--------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(ggplot2)

# ==================== 1. LOAD DATA ====================

# Set working directory to script location
script_dir <- "~/bonsaiforest2/simulations/Continuous"
setwd(path.expand(script_dir))

# Load ground truth
load("Scenarios/truth.RData")

# Load scenario 3 data
scenario3_data <- readRDS("Scenarios/scenario3.rds")

# Convert to data frame and filter to first simulation
sim_data <- as.data.frame(scenario3_data)
sim_data <- sim_data[sim_data$sim_id == 1, ]

# Get scenario 3 truth data
truth_scenario3 <- simulation_truth[simulation_truth$scenario_no == 3, ]

# ==================== 2. COMPUTE SUBGROUP ESTIMATORS ====================

# Define all subgrouping variables (all already categorical in the data)
subgroup_definitions <- list(
  X1 = unique(na.omit(sim_data$X1)),
  X2 = unique(na.omit(sim_data$X2)),
  X3 = unique(na.omit(sim_data$X3)),
  X4 = unique(na.omit(sim_data$X4)),
  X8 = unique(na.omit(sim_data$X8)),
  X11cat = unique(na.omit(sim_data$X11cat)),
  X14cat = unique(na.omit(sim_data$X14cat)),
  X17cat = unique(na.omit(sim_data$X17cat))
)

# Compute treatment effect estimates for each subgroup
subgroup_estimates <- NULL

for (var in names(subgroup_definitions)) {
  for (level in subgroup_definitions[[var]]) {
    # Filter data to subgroup using base R
    subgroup_data <- sim_data[sim_data[[var]] == level, ]
    
    # Fit linear model: Y ~ trt
    if (nrow(subgroup_data) > 1) {
      model <- lm(Y ~ trt, data = subgroup_data)
      trt_estimate <- coef(model)["trt"]
      n_subgroup <- nrow(subgroup_data)
      
      # Get confidence interval
      ci <- confint(model, "trt", level = 0.95)
      ci_lower <- ci[1]
      ci_upper <- ci[2]
      
      # Get standard error
      se <- summary(model)$coefficients["trt", "Std. Error"]
      
      # Store estimate
      subgroup_estimates <- rbind(subgroup_estimates,
                                  data.frame(
                                    subgroup_var = var,
                                    level = level,
                                    subgroup = paste(var, level, sep = "."),
                                    trt_estimate = trt_estimate,
                                    se = se,
                                    ci_lower = ci_lower,
                                    ci_upper = ci_upper,
                                    n_subgroup = n_subgroup
                                  ))
    }
  }
}

# ==================== 3. MERGE WITH TRUTH ====================

# Merge estimates with truth data
# Truth is computed for each unique level of each variable
truth_for_merge <- truth_scenario3[, c("subgroup_var", "level", "trt_effect")]

# Convert level to character for matching
subgroup_estimates$level_char <- as.character(subgroup_estimates$level)
truth_for_merge$level_char <- as.character(truth_for_merge$level)

# Merge on variable and level
comparison <- merge(subgroup_estimates,
                   truth_for_merge,
                   by.x = c("subgroup_var", "level_char"),
                   by.y = c("subgroup_var", "level_char"),
                   all.x = TRUE)
names(comparison)[names(comparison) == "trt_effect"] <- "trt_truth"

# Add a numeric order for sorting subgroups
comparison$var_order <- as.numeric(factor(comparison$subgroup_var,
                                          levels = names(subgroup_definitions)))

# Sort by variable and level
comparison <- comparison[order(comparison$var_order, as.character(comparison$level_char)), ]

# Create a numeric y-position for each subgroup in the plot
comparison$y_pos <- seq_along(comparison$subgroup)

# ==================== 4. DISPLAY RESULTS ====================

cat("\n========== Scenario 3: Subgroup Treatment Effect Estimates vs Truth (All Covariates) ==========\n")
cat("Covariates: X1-X4, X8, X11cat, X14cat, X17cat\n")
display_cols <- c("subgroup", "trt_estimate", "ci_lower", "ci_upper", "trt_truth", "n_subgroup")
print(comparison[, display_cols])

# ==================== 5. FOREST PLOT ====================

# Create a single forest plot showing both estimates and truth
# Prepare estimate data
estimate_data <- comparison[, c("subgroup", "y_pos", "subgroup_var", "trt_estimate", "ci_lower", "ci_upper")]
estimate_data$type <- "Estimate"
names(estimate_data)[4] <- "value"

# Prepare truth data
truth_data <- comparison[, c("subgroup", "y_pos", "subgroup_var", "trt_truth")]
truth_data$ci_lower <- NA
truth_data$ci_upper <- NA
truth_data$type <- "Truth"
names(truth_data)[4] <- "value"

# Combine for plotting
plot_data <- rbind(estimate_data, truth_data)

# Create forest plot with both estimates and truth
forest <- ggplot(plot_data, aes(x = value, y = reorder(subgroup, y_pos), color = type, shape = type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(data = estimate_data,
                aes(xmin = ci_lower, xmax = ci_upper, color = type),
                height = 0.2, position = position_dodge(width = 0.5), size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
  scale_color_manual(values = c("Estimate" = "darkblue", "Truth" = "darkred")) +
  scale_shape_manual(values = c("Estimate" = 21, "Truth" = 23)) +
  theme_bw() +
  labs(
    title = "Scenario 3: Forest Plot - Subgroup Treatment Effects (Estimates vs Truth)",
    subtitle = "All 8 covariates - First simulation (n=500)",
    x = "Treatment Effect",
    y = "Subgroup",
    color = "Type",
    shape = "Type"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("forest_plot_scenario3.png", plot = forest, width = 14, height = 10, dpi = 300)
message("✓ Forest plot saved as forest_plot_scenario3.png")

# Display plot
print(forest)
