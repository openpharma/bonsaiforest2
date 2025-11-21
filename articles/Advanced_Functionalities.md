# Advanced Functionalities

This vignette demonstrates advanced functionalities for Bayesian
shrinkage estimation, including specifying an **offset for count data**,
customizing **prior distributions**, and using **stratification** to
handle heterogeneity in nuisance parameters.

## 1 Handling Exposure in Count Outcomes with Offsets

For count data (like disease exacerbations or event rates), the observed
count often depends on the **exposure time** or **follow-up time**.
Using an **offset** variable is the standard statistical method to
properly account for this time variation by modeling the event **rate**
(count per unit time) instead of the raw count.

In `bonsaiforest2`, you can include the offset directly in the
`response_formula_str` using the standard `brms` syntax:
`outcome + offset(log_time) ~ predictors`.

### 1.1 Example 1: Count Outcome with Offset (Disease Exacerbations)

This scenario models exacerbation counts using a Negative Binomial
distribution and explicitly accounts for the patient’s exposure time.

*Scenario*: Modeling exacerbation counts. Adjust for `baseline_severity`
(unshrunk) and many exploratory `biomarkers` (shrunk). No interaction
terms (overall treatment effect only).

``` r
# Data Simulation 
set.seed(789)
library(bonsaiforest2)
n_patients <- 150
biomarker_data <- as.data.frame(matrix(rnorm(n_patients * 10), ncol = 10))
names(biomarker_data) <- paste0("biomarker_", 1:10)
count_data <- data.frame(
  exacerbation_count = rnbinom(n_patients, size = 1.5, mu = 3),
  medication = sample(0:1, n_patients, replace = TRUE),
  baseline_severity = rnorm(n_patients, 10, 2),
  log_exposure_time = log(runif(n_patients, 0.5, 1.5)) # Example offset
)
count_data <- cbind(count_data, biomarker_data)
count_data$medication <- factor(count_data$medication, levels = c(0, 1))

# Create formula string for all biomarkers (Prognostic effects only here)
shrunk_prog_str <- paste("~", paste(names(biomarker_data), collapse = " + "))
```

``` r
# Model Fitting with Offset
count_model_fit <- run_brms_analysis(
  data = count_data,
  # Include offset(log_exposure_time) directly in the response formula
  response_formula_str = "exacerbation_count + offset(log_exposure_time) ~ medication",
  response_type = "count",
  unshrunk_prognostic_formula_str = "~ 1 + baseline_severity",
  shrunk_prognostic_formula_str = shrunk_prog_str,
  chains = 1, iter = 200, warmup = 100, cores = 1, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> Treatment 'medication' added to unshrunk prognostic terms by default.
#> 
#> Step 2: Fitting the brms model...
#> Using default priors for unspecified effects:
#>   - shrunk prognostic (b): horseshoe(1)
#>   - unshrunk prognostic (b): normal(0, 5)
#>   - prognostic intercept: normal(0, 5)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10 
#> Chain 1 finished in 2.7 seconds.
#> Loading required namespace: rstan
#> 
#> Analysis complete.
```

## 2 Customizing Prior Distributions

The choice of prior is central to Bayesian shrinkage. `bonsaiforest2`
provides sensible defaults, but it allows for full customization using
the `prognostic_effect_priors` and `predictive_effect_priors` arguments.

### 2.1 Prior Specification Mechanics

You specify priors as named lists containing strings, which are then
passed to
[`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html).

| Prior Group                            | Argument                   | Sub-List Names                    | Default Prior (Type)       |
|:---------------------------------------|:---------------------------|:----------------------------------|:---------------------------|
| **Predictive** (Interaction/Treatment) | `predictive_effect_priors` | `shrunk`, `unshrunk`              | `horseshoe(1)` (Shrinkage) |
| **Prognostic** (Covariates)            | `prognostic_effect_priors` | `shrunk`, `unshrunk`, `intercept` | `horseshoe(1)` (Shrinkage) |

#### 2.1.1 Prior Recommendations for Non-Shrunk Terms

For non-shrunk (fixed) effects, choosing a **weakly informative prior**
is critical to stabilize the model without unduly influencing the
results. A robust strategy is to scale the priors based on the
**standard deviation (SD)** of the outcome variable.

\\\sigma\_{outcome} = \text{SD}(\text{Outcome})\\

| Parameter               | Recommended Prior Scale                           | Justification                                                                                               |
|:------------------------|:--------------------------------------------------|:------------------------------------------------------------------------------------------------------------|
| **Intercept**           | \\\text{Normal}(0, 10 \times \sigma\_{outcome})\\ | Wide enough to cover the range of outcomes if all predictors are zero.                                      |
| **Unshrunk Prognostic** | \\\text{Normal}(0, \sigma\_{outcome})\\           | Assumes a typical coefficient’s impact is roughly comparable to the magnitude of the outcome’s variability. |

### 2.2 Practical Examples of Prior Setting

To provide a functional example, we will generate a **synthetic
dataset** that mirrors the structure of the **Tirzepatide (SURPASS-2)
trial** case study described in Wang et al. (2024). This study used a
**continuous endpoint** (change in HbA1c), which fits perfectly with
your code snippet calculating the standard deviation (`sd`).

#### 2.2.1 Dataset Generation and Model Preparation (Synthetic SURPASS-2)

This dataset mimics a clinical trial comparing a Treatment (Tirzepatide)
vs. Control (Semaglutide), looking at subgroups defined by **Sex**,
**Race**, and **Age Group**.

``` r
# Prerequisites: Ensure packages used by your function are loaded
library(survival) # For Surv()
library(brms)     # For brmsformula objects
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.0). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: 'brms'
#> The following object is masked from 'package:survival':
#> 
#>     kidney
#> The following object is masked from 'package:stats':
#> 
#>     ar
library(bonsaiforest2)

# 1. Create Sample Data (matching your @examples section)
set.seed(123)
n <- 100
sim_data <- data.frame(
  time = round(runif(n, 1, 100)),
  status = sample(0:1, n, replace = TRUE),
  trt = sample(0:1, n, replace = TRUE),
  age = rnorm(n, 50, 10),
  region = sample(c("A", "B"), n, replace = TRUE),
  subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
)

# Ensure variables are factors where appropriate
sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
sim_data$region <- as.factor(sim_data$region)
sim_data$subgroup <- as.factor(sim_data$subgroup)

# 2. Run prepare_formula_model
# This transforms the data (creating dummy variables for interactions)
# and builds the non-linear brms formula.
prepared_model <- prepare_formula_model(
  data = sim_data,
  
  # Main response variable and treatment identifier
  response_formula_str = "Surv(time, status) ~ trt",
  
  # Predictive (Interaction) effects to be shrunk
  # Note: The function will automatically generate dummy columns for 'subgroup'
  shrunk_predictive_formula_str = "~ trt:subgroup",
  
  # Prognostic (Main) effects: Unshrunk (Fixed)
  unshrunk_prognostic_formula_str = "~ age",
  
  # Prognostic (Main) effects: Shrunk (Regularized)
  shrunk_prognostic_formula_str = "~ region",
  
  # Model type
  response_type = "survival",
  
  # Stratify the baseline hazard by region
  stratification_formula_str = "~ region"
)
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup

# 3. Inspect the Results

# A. The generated brms formula object
# You should see parts for 'unprogeffect', 'shprogeffect', etc.
print(prepared_model$formula)
#> time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unprogeffect + shprogeffect + shpredeffect 
#> unprogeffect ~ age + trt + subgroup + 0
#> shprogeffect ~ region + 0
#> shpredeffect ~ subgroup_S1_x_trt + subgroup_S2_x_trt + subgroup_S3_x_trt + 0

# B. The processed data
# You should see new columns like 'subgroup_S1_x_trt', 'subgroup_S2_x_trt', etc.
head(prepared_model$data)
#>   time status trt      age region subgroup subgroup_S1_x_trt subgroup_S2_x_trt
#> 1   29      0   1 57.87739      B       S2                 0                 1
#> 2   79      1   1 57.69042      B       S1                 1                 0
#> 3   41      1   1 53.32203      B       S3                 0                 0
#> 4   88      0   0 39.91623      B       S3                 0                 0
#> 5   94      1   1 48.80547      A       S3                 0                 0
#> 6    6      1   1 47.19605      A       S2                 0                 1
#>   subgroup_S3_x_trt
#> 1                 0
#> 2                 0
#> 3                 1
#> 4                 0
#> 5                 1
#> 6                 0
```

#### 2.2.2 Model Preparation

#### 2.2.3 Example 2: Using Recommended Weakly Informative Priors

This is the easiest example we can have. Just use the recommended priors
without any tunning.

``` r
fit_ex3 <- fit_brms_model(
  prepared_model = prepared_model
)
#> Using default priors for unspecified effects:
#>   - shrunk prognostic (b): horseshoe(1)
#>   - unshrunk prognostic (b): normal(0, 5)
#>   - shrunk predictive (b): horseshoe(1)
#> Fitting brms model...
#> Compiling Stan program...
#> Trying to compile a simple C file
#> Running /opt/R/4.5.2/lib/R/bin/R CMD SHLIB foo.c
#> using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
#> gcc -std=gnu2x -I"/opt/R/4.5.2/lib/R/include" -DNDEBUG   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
#> In file included from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Core:19,
#>                  from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Dense:1,
#>                  from /home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22,
#>                  from <command-line>:
#> /home/runner/work/_temp/Library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: cmath: No such file or directory
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> compilation terminated.
#> make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:202: foo.o] Error 1
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000317 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.17 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 5.153 seconds (Warm-up)
#> Chain 1:                1.854 seconds (Sampling)
#> Chain 1:                7.007 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.59 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 5.481 seconds (Warm-up)
#> Chain 2:                1.634 seconds (Sampling)
#> Chain 2:                7.115 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 6.1e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.61 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 4.991 seconds (Warm-up)
#> Chain 3:                1.728 seconds (Sampling)
#> Chain 3:                6.719 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.6 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 5.216 seconds (Warm-up)
#> Chain 4:                2.025 seconds (Sampling)
#> Chain 4:                7.241 seconds (Total)
#> Chain 4:
#> Warning: There were 4 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

#### 2.2.4 Example 3: Using the R2D2 Shrinkage Prior

The R2D2 prior is useful when you want to control the *global* shrinkage
via the coefficient of determination (\\R^2\\) rather than a scale
parameter \\\tau\\. This is often more interpretable for stakeholders.

``` r
# Use a custom R2D2 string for the shrunk predictive effects (interactions)
# mean_R2 = 0.5 implies we expect the subgroups to explain 50% of the variance
r2d2_prior_string <- "R2D2(mean_R2 = 0.5, prec_R2 = 1)" 

fit_ex4 <- fit_brms_model(
  prepared_model = prepared_model,
  predictive_effect_priors = list(shrunk = r2d2_prior_string)
)
#> Using default priors for unspecified effects:
#>   - shrunk prognostic (b): horseshoe(1)
#>   - unshrunk prognostic (b): normal(0, 5)
#> Fitting brms model...
#> Compiling Stan program...
#> Trying to compile a simple C file
#> Running /opt/R/4.5.2/lib/R/bin/R CMD SHLIB foo.c
#> using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
#> gcc -std=gnu2x -I"/opt/R/4.5.2/lib/R/include" -DNDEBUG   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
#> In file included from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Core:19,
#>                  from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Dense:1,
#>                  from /home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22,
#>                  from <command-line>:
#> /home/runner/work/_temp/Library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: cmath: No such file or directory
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> compilation terminated.
#> make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:202: foo.o] Error 1
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000102 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.02 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 5.058 seconds (Warm-up)
#> Chain 1:                2.417 seconds (Sampling)
#> Chain 1:                7.475 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.54 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 5.265 seconds (Warm-up)
#> Chain 2:                2.158 seconds (Sampling)
#> Chain 2:                7.423 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 5.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.53 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 4.621 seconds (Warm-up)
#> Chain 3:                1.557 seconds (Sampling)
#> Chain 3:                6.178 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 5.8e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.58 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 5.041 seconds (Warm-up)
#> Chain 4:                1.624 seconds (Sampling)
#> Chain 4:                6.665 seconds (Total)
#> Chain 4:
#> Warning: There were 6 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

#### 2.2.5 Example 4: Custom Hierarchical Prior (Advanced)

This example demonstrates injecting raw Stan code. This is necessary if
you want to implement a hierarchical structure that `brms` does not
support natively, such as linking specific coefficients via a custom
global parameter \\\tau\_{pred}\\.

``` r
library(brms)

# 1. Define a new hyperparameter 'tau_pred' for the Stan code
stanvars_full_hierarchical <- brms::stanvar(
  scode = "  real mu_pred;\n  real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  brms::stanvar(
    scode = "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); \n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); \n",
    block = "model"
  )

# 
prior_full_hierarchical <- brms::set_prior("normal(mu_pred, sigma_pred)")

# 3. Pass both objects to the function
fit_ex5 <- fit_brms_model(
  prepared_model = prepared_model,
  predictive_effect_priors = list(shrunk = prior_full_hierarchical),
  stanvars = stanvars_full_hierarchical
)
#> Re-targeting 'brmsprior' object for nlpar: shpredeffect class: b
#> Using default priors for unspecified effects:
#>   - shrunk prognostic (b): horseshoe(1)
#>   - unshrunk prognostic (b): normal(0, 5)
#> Fitting brms model...
#> Compiling Stan program...
#> Trying to compile a simple C file
#> Running /opt/R/4.5.2/lib/R/bin/R CMD SHLIB foo.c
#> using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
#> gcc -std=gnu2x -I"/opt/R/4.5.2/lib/R/include" -DNDEBUG   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
#> In file included from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Core:19,
#>                  from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Dense:1,
#>                  from /home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22,
#>                  from <command-line>:
#> /home/runner/work/_temp/Library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: cmath: No such file or directory
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> compilation terminated.
#> make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:202: foo.o] Error 1
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.92 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 8.413 seconds (Warm-up)
#> Chain 1:                8.748 seconds (Sampling)
#> Chain 1:                17.161 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.5 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 7.685 seconds (Warm-up)
#> Chain 2:                0.121 seconds (Sampling)
#> Chain 2:                7.806 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.5 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 7.733 seconds (Warm-up)
#> Chain 3:                3.844 seconds (Sampling)
#> Chain 3:                11.577 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 5.1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.51 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 7.623 seconds (Warm-up)
#> Chain 4:                4.567 seconds (Sampling)
#> Chain 4:                12.19 seconds (Total)
#> Chain 4:
#> Warning: There were 1310 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.56, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
```

\#Stratification for Nuisance Parameters

In many trials, parameters like the observation error variance
(\\\sigma^2\\ for continuous outcomes) or the baseline hazard function
(\\h_0(t)\\ for survival outcomes) are known to vary by site, country,
or other factors. Stratification models this known heterogeneity by
fitting these nuisance parameters separately for each level of a
grouping variable.

Use the `stratification_formula_str` argument to define the grouping
factor(s).

### 2.3 Example 5: Stratified Continuous Model

*Scenario*: We model SBP change where the observation **variance**
\\\sigma^2\\ differs by \\\text{clinic\\site}\\. This is modeled by
stratifying the \\\sigma\\ parameter in the Normal distribution by the
site variable.

``` r
set.seed(42)
n_patients <- 180
sigma_by_site <- c(SiteA = 6, SiteB = 12, SiteC = 18)
sample_data_strat_cont <- data.frame(
    id = 1:n_patients,
    clinic_site = factor(sample(c("SiteA", "SiteB", "SiteC"), n_patients, replace = TRUE))
)
sample_data_strat_cont$trt <- factor(sample(0:1, n_patients, replace = TRUE, prob = c(0.5, 0.5)))
sample_data_strat_cont$baseline_sbp <- rnorm(n_patients, mean = 145, sd = 8)
sample_data_strat_cont$age_group <- factor(sample(c("Under60", "Over60"), n_patients, replace = TRUE))
noise <- rnorm(n_patients, mean = 0, sd = sigma_by_site[sample_data_strat_cont$clinic_site])
sample_data_strat_cont$sbp_change <- -5 - (as.integer(sample_data_strat_cont$trt)-1) * 8 -
                                     0.2 * (sample_data_strat_cont$baseline_sbp - 145) +
                                     (as.integer(sample_data_strat_cont$trt)-1) * ifelse(sample_data_strat_cont$age_group == "Over60", -5, 0) +
                                     noise
```

``` r
# Model Fitting with Stratified Sigma
fit_continuous_stratified <- run_brms_analysis(
  data = sample_data_strat_cont,
  response_formula_str = "sbp_change ~ trt",
  response_type = "continuous",
  unshrunk_prognostic_formula_str = "~ baseline_sbp",
  shrunk_predictive_formula_str = "~ age_group:trt",
  stratification_formula_str = "~ clinic_site",
  chains = 1, iter = 200, warmup = 100, cores = 1, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> Applying stratification: estimating sigma by 'clinic_site'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: age_group
#> 
#> Step 2: Fitting the brms model...
#> Using default priors for unspecified effects:
#>   - unshrunk prognostic (b): normal(0, 5)
#>   - prognostic intercept: normal(0, 5)
#>   - shrunk predictive (b): horseshoe(1)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10 
#> Chain 1 finished in 2.7 seconds.
#> Warning: 82 of 100 (82.0%) transitions hit the maximum treedepth limit of 10.
#> See https://mc-stan.org/misc/warnings for details.
#> 
#> Analysis complete.
```

``` r
strat_continuous_summary <- summary_subgroup_effects(
  brms_fit = fit_continuous_stratified,
  original_data = sample_data_strat_cont,
  trt_var = "trt",
  response_type = "continuous" # Auto-detects age_group interaction
)
#> --- Calculating specific subgroup effects... ---
#> `subgroup_vars` set to 'auto'. Detecting from model interaction terms...
#> ...detected subgroup variable(s): age_group
#> Step 1: Creating counterfactual datasets...
#> ...setting interaction dummy variables for the 'all treatment' scenario.
#> Step 2: Generating posterior predictions...
#> ... (predicting expected outcomes)...
#> Step 3: Calculating marginal effects...
#> ... processing age_group
#> Done.

print(strat_continuous_summary$estimates)
#> # A tibble: 2 × 4
#>   Subgroup           Median CI_Lower CI_Upper
#>   <chr>               <dbl>    <dbl>    <dbl>
#> 1 age_group: Over60  -10.8     -14.5    -7.17
#> 2 age_group: Under60  -7.95    -12.0    -5.00
plot(strat_continuous_summary, title = "Stratified Continuous: Subgroup Effects")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Advanced_Functionalities_files/figure-html/unnamed-chunk-10-1.png)

### 2.4 Example 6: Stratified Survival Model

*Scenario*: We model a time-to-event outcome where the **baseline
hazard** \\h_0(t)\\ differs by \\\text{country}\\. This is particularly
important for multi-regional trials where regional differences in
standard of care affect baseline risk.

``` r
set.seed(123)
n_patients <- 200
lambda_by_country <- c(US = 0.01, EU = 0.03)
surv_data_strat <- data.frame(
    id = 1:n_patients,
    country = factor(sample(c("US", "EU"), n_patients, replace = TRUE)),
    trt = factor(sample(0:1, n_patients, replace = TRUE)),
    age = rnorm(n_patients, 65, 10),
    biomarker = factor(sample(c("Low", "High"), n_patients, replace = TRUE))
)
lp <- (as.integer(surv_data_strat$trt)-1) * -0.6 + (surv_data_strat$age - 65) * 0.03 +
      (as.integer(surv_data_strat$trt)-1) * (as.integer(surv_data_strat$biomarker) - 1) * -0.5
u <- runif(n_patients)
lambda_vec <- lambda_by_country[surv_data_strat$country]
gamma <- 1.5 # Weibull shape parameter
true_event_time <- (-log(u) / (lambda_vec * exp(lp)))^(1/gamma)
censoring_time <- 60 # Administrative censoring time
surv_data_strat$event_status <- ifelse(true_event_time <= censoring_time, 1, 0)
surv_data_strat$event_time <- pmin(true_event_time, censoring_time)
```

``` r
# Model Fitting with Stratified Baseline Hazard
fit_surv_stratified <- run_brms_analysis(
  data = surv_data_strat,
  response_formula_str = "Surv(event_time, event_status) ~ trt",
  response_type = "survival",
  unshrunk_prognostic_formula_str = "~ age",
  shrunk_predictive_formula_str = "~ biomarker:trt",
  stratification_formula_str = "~ country",
  chains = 1, iter = 200, warmup = 100, cores = 1, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'country'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: biomarker
#> 
#> Step 2: Fitting the brms model...
#> Using default priors for unspecified effects:
#>   - unshrunk prognostic (b): normal(0, 5)
#>   - shrunk predictive (b): horseshoe(1)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10 
#> Chain 1 finished in 7.8 seconds.
#> Warning: 100 of 100 (100.0%) transitions hit the maximum treedepth limit of 10.
#> See https://mc-stan.org/misc/warnings for details.
#> 
#> Analysis complete.
```

``` r
strat_surv_summary <- summary_subgroup_effects(
  brms_fit = fit_surv_stratified,
  original_data = surv_data_strat,
  trt_var = "trt",
  response_type = "survival" # Auto-detects biomarker interaction
)
#> --- Calculating specific subgroup effects... ---
#> `subgroup_vars` set to 'auto'. Detecting from model interaction terms...
#> ...detected subgroup variable(s): biomarker
#> Step 1: Creating counterfactual datasets...
#> ...setting interaction dummy variables for the 'all treatment' scenario.
#> Step 2: Generating posterior predictions...
#> ... (reconstructing baseline hazard and getting linear predictors)...
#> ... model is stratified by 'country'. Reconstructing hazard for each stratum.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Step 3: Calculating marginal effects...
#> ... processing biomarker
#> Done.

print(strat_surv_summary$estimates)
#> # A tibble: 2 × 4
#>   Subgroup        Median CI_Lower CI_Upper
#>   <chr>            <dbl>    <dbl>    <dbl>
#> 1 biomarker: High  0.568    0.419    0.760
#> 2 biomarker: Low   0.364    0.261    0.522
plot(strat_surv_summary, title = "Stratified Survival: Subgroup Effects (AHR)")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Advanced_Functionalities_files/figure-html/ex6-summary-plot-1.png)
