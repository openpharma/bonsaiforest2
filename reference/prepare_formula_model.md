# Prepare a Multi-Part brms Formula and Corresponding Data

This function serves as a powerful pre-processor for building complex
Bayesian models with the \`brms\` package. It automates the construction
of a multi-part, non-linear formula by classifying covariates into four
distinct categories: unshrunk prognostic, shrunk prognostic, unshrunk
predictive, and shrunk predictive.

## Usage

``` r
prepare_formula_model(
  data,
  response_formula_str,
  shrunk_predictive_formula_str = NULL,
  unshrunk_prognostic_formula_str = NULL,
  unshrunk_predictive_formula_str = NULL,
  shrunk_prognostic_formula_str = NULL,
  response_type = c("binary", "count", "continuous", "survival"),
  stratification_formula_str = NULL
)
```

## Arguments

- data:

  A data.frame containing all the necessary variables.

- response_formula_str:

  A character string for the response part, e.g., "outcome ~ trt", for
  count models "n_events + offset(log(days)) ~ trt" or for survival
  models "Surv(time,status) ~ trt".

- shrunk_predictive_formula_str:

  Predictive terms to be shrunk ('shpredeffect'). These are typically
  interactions with the treatment variable.

- unshrunk_prognostic_formula_str:

  Prognostic terms not to be shrunk ('unprogeffect'). These are main
  effects assumed to be important.

- unshrunk_predictive_formula_str:

  Predictive terms not to be shrunk ('unpredeffect').

- shrunk_prognostic_formula_str:

  Prognostic terms to be shrunk ('shprogeffect'). These are main effects
  where regularization is desired.

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival".

- stratification_formula_str:

  A formula string specifying the stratification variable, e.g., "~
  strata_var".

## Value

A list with three elements: \`formula\` (a \`brmsformula\` object),
\`data\` (the modified data.frame), and \`response_type\` (the character
string of the response type).

## Details

This classification allows for applying differential shrinkage
(regularization) to different parts of the model. The function also
prepares the corresponding data using R's contrast coding system for
proper factor handling and interactions.

## Key Features

- **Multi-Part Formula Construction:** It generates a \`brmsformula\`
  object with up to four distinct linear components (\`unprogeffect\`,
  \`shprogeffect\`, \`unpredeffect\`, \`shpredeffect\`), which are
  combined in a non-linear model. This allows for assigning different
  priors to each component.

- **Automated Interaction Handling:** Predictive terms (e.g., \`~
  trt:subgroup\`) are automatically handled using R's contrast coding.
  All variables involved in interactions are converted to factors with
  treatment contrasts (reference coding), enabling proper model fitting
  and prediction on new data without manual dummy variable creation.

- **Hierarchical Integrity:** When a predictive term like
  \`trt:subgroup\` is specified, the function automatically ensures that
  the corresponding prognostic (main) effect \`subgroup\` is included in
  the model, adhering to the principle of model hierarchy.

- **Intelligent Defaults:** The main treatment variable is automatically
  added as an unshrunk prognostic term if not specified elsewhere. It
  also provides warnings and resolves overlaps if a term is accidentally
  specified in multiple categories.

## Survival Model Details (Robust Spline Knot Calculation)

For survival models (\`response_type = "survival"\`), the function
explicitly models the baseline hazard using B-splines via
\`brms::bhaz()\`. It implements a highly robust method for calculating
spline knots to ensure stability:

1.  **Boundary Definition:** It first establishes boundary knots by
    taking the range of unique event times and adding a small buffer.
    This creates a "safe" interval for placing internal knots.

2.  **Internal Knot Placement:** It calculates internal knots using
    quantiles of the event times that fall strictly within the boundary
    knots. This prevents knots from being placed at the exact edges of
    the data, which can cause numerical instability.

3.  **Fallback for Sparse Data:** If there are too few unique event
    times to calculate quantile-based knots, it gracefully falls back to
    placing evenly spaced knots within the defined boundaries.

## Data Transformation

It is critical to note that this function returns a modified
\`data.frame\`. The treatment variable and any variables involved in
interactions are converted to factors with treatment contrasts
(reference coding). This contrast-based approach enables proper model
fitting and easy prediction on new data, as the same contrasts are
automatically applied. The returned \`data\` object should be used in
the subsequent call to \`brms::brm()\`, not the original data.

## Stratification

The \`stratification_formula_str\` argument allows for estimating
certain parameters separately for different groups. Its behavior depends
on the \`response_type\`:

- `survival`: Estimates a separate baseline hazard function (\`bhaz\`)
  for each level of the stratification variable.

- `continuous`: Models the residual standard deviation (\`sigma\`) as
  varying across levels of the stratification variable.

- `count`: Models the overdispersion parameter (\`shape\`) as varying
  across levels of the stratification variable.

## Examples

``` r
if (require("brms") && require("survival")) {
  # 1. Create Sample Data
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
  sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
  sim_data$region <- as.factor(sim_data$region)
  sim_data$subgroup <- as.factor(sim_data$subgroup)

  # 2. Run the function
  prepared_model <- prepare_formula_model(
    data = sim_data,
    response_formula_str = "Surv(time, status) ~ trt",
    shrunk_predictive_formula_str = "~ trt:subgroup",
    unshrunk_prognostic_formula_str = "~ age",
    shrunk_prognostic_formula_str = "~ region",
    response_type = "survival",
    stratification_formula_str = "~ region"
  )

  # 3. View the results
  print(prepared_model$formula)
  print(head(prepared_model$data))
}
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup
#> time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unprogeffect + shprogeffect + shpredeffect 
#> unprogeffect ~ age + trt + subgroup + 0
#> shprogeffect ~ region + 0
#> shpredeffect ~ trt_subgroupS1 + trt_subgroupS2 + trt_subgroupS3 + 0
#>   time status trt      age region subgroup trt_subgroupS1 trt_subgroupS2
#> 1   29      0   1 57.87739      B       S2              0              1
#> 2   79      1   1 57.69042      B       S1              1              0
#> 3   41      1   1 53.32203      B       S3              0              0
#> 4   88      0   0 39.91623      B       S3              0              0
#> 5   94      1   1 48.80547      A       S3              0              0
#> 6    6      1   1 47.19605      A       S2              0              1
#>   trt_subgroupS3
#> 1              0
#> 2              0
#> 3              1
#> 4              0
#> 5              1
#> 6              0
```
