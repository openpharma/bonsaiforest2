#--------------------------------------------------------------------------------------------------
#------ SIMULATION TRUTH AND SCENARIO DATA GENERATION FOR CONTINUOUS OUTCOMES
#------
#------ Select continuous scenarios from Sun et al (https://doi.org/10.1002/bimj.202100337) 
#------ with interaction strengths b1_rel = 0, 1, 2
#------ Re-scale parameter b0 to increase power for overall comparison from 50% to 90%
#------ Select only subgrouping variables X1, X2, X3, X4, X8, X11 (dichotomized at median), 
#------ X14 (cut at Q1 and median), X17 (dichotomized at median) 
#------ (8 subgrouping variables defining 20 subgroups in total)
#--------------------------------------------------------------------------------------------------

# devtools::install_github("Sophie-Sun/benchtm")
library(benchtm)
library(dplyr)
library(tidyverse)
library(ggplot2)

set.seed(0)

# ==================== 1. SELECT AND CALIBRATE SCENARIOS ====================

#------ Select continuous parameter scenarios with interaction strengths b1_rel = 0, 1, or 2
selected_scen_par_power50 <- subset(scen_param, (type == "continuous") & (b1_rel %in% c(0,1,2)))

selected_scen_par_power50 <- 
  selected_scen_par_power50 |>
  mutate(
    scenario_no = 1:12,
    scenario = fct_inorder(paste(rep(1:4, each = 3), case_match(b1_rel, 0 ~ "homo", 1~"hetero_mild", 2 ~"hetero_strong"), sep = "_"))
  ) |>
  relocate(scenario_no, scenario)

#------ Re-scale b0 of selected_scen_par to ensure 90% (rather than 50%) power for group comparison
X <- benchtm:::X_large_pop
selected_scen_par <- selected_scen_par_power50

for (i in (1:nrow(selected_scen_par))){
  # generate initial dataset with 50% power
  scen <- selected_scen_par[i,]
  df_tmp <- generate_y(
    X = X, 
    trt = rep(c(0,1), each = nrow(X)/2),
    prog = scen$prog, pred = scen$pred, 
    b0 = scen$b0, b1 = scen$b1, type = scen$type, 
    include_truth = TRUE)
  
  # get treatment effect delta and re-calibrate b0 to reach 90% power
  delta <- mean(df_tmp$trt_effect)
  new_delta <- delta * abs((qnorm(0.025) + qnorm(0.1)) / (qnorm(0.025) + qnorm(0.5))) # delta for 80% power
  
  selected_scen_par$b0[i] <- selected_scen_par$b0[i] + (new_delta - delta)
}

# ==================== 2. DETERMINE SIMULATION TRUTH FOR SELECTED SUBGROUPS ====================

# Create categorical versions of continuous variables at specified cutpoints
X$X11cat <- cut(X$X11, c(-Inf, 0.4615385, Inf), labels = c("a", "b")) # cut at median
X$X14cat <- cut(X$X14, c(-Inf, 0.2478632, 0.3333833, Inf), labels = c("a", "b", "c")) # cut at Q1 and median
X$X17cat <- cut(X$X17, c(-Inf, 0.5888801, Inf), labels = c("a", "b")) # cut at median 

simulation_truth <- NULL 

for (i in (1:nrow(selected_scen_par))){
  scen <- selected_scen_par[i,]
  
  df <- generate_y(
    X = X, 
    trt = rep(c(0,1), each = nrow(X)/2), 
    prog = scen$prog, pred = scen$pred, 
    b0 = scen$b0, b1 = scen$b1, type = scen$type, 
    include_truth = TRUE)
  
  pop_trt_effect <- mean(df$trt_effect) # overall treatment effect
  sd_approx <- sigma(lm(Y ~ trt, data = df)) # approximate sd

  true_subgroup_effects_scenario <- df |>
    select(X1, X2, X3, X4, X8, X11cat, X14cat, X17cat, trt_effect) |>
    pivot_longer(cols = starts_with("X"), names_to = "subgroup_var", values_to = "level") |>
    mutate(subgroup_var = fct_inorder(subgroup_var)) |>
    group_by(subgroup_var, level) |>
    summarise(trt_effect = mean(trt_effect, na.rm = TRUE), .groups = "drop") |>
    mutate(
      subgroup = fct_inorder(paste(subgroup_var, level, sep = ".")),
      scenario_no = scen$scenario_no,
      scenario = scen$scenario,
      pop_trt_effect = pop_trt_effect,
      sd_approx = sd_approx) |>
    relocate(c(scenario_no, scenario, subgroup))
  
  simulation_truth <- rbind(simulation_truth, true_subgroup_effects_scenario)
}

# ==================== 3. PLOT SIMULATION TRUTH ====================

ggplot(simulation_truth, aes(x = subgroup, y = trt_effect, color = subgroup_var)) +
  geom_point(size = 3) + 
  facet_wrap(~ scenario, ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "True Treatment Effect by Subgroup and Scenario",
    x = "Subgroup",
    y = "True Treatment Effect",
    color = "Variable"
  )

# Save simulation truth
save(simulation_truth, file = "Scenarios/truth.RData")
message("✓ Simulation truth saved to Scenarios/truth.RData")
