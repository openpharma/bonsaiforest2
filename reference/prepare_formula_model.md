# Prepare a Multi-Part brms Formula and Corresponding Data

This function serves as a pre-processor for building complex Bayesian
models with the `brms` package. It automates the construction of a
multi-part, non-linear formula by classifying covariates into three
distinct categories: unshrunk terms, shrunk prognostic, and shrunk
predictive.

## Usage

``` r
prepare_formula_model(
  data,
  response_formula,
  unshrunk_terms_formula = NULL,
  shrunk_prognostic_formula = NULL,
  shrunk_predictive_formula = NULL,
  response_type = c("binary", "count", "continuous", "survival"),
  stratification_formula = NULL
)
```

## Arguments

- data:

  A data frame. Dataset containing all necessary variables for model
  fitting. Must include the response variable, treatment variable, and
  all covariates specified in the formula arguments. This data will be
  modified (contrast coding applied) and returned for use in model
  fitting.

- response_formula:

  A formula object. The response specification defining the outcome
  variable and treatment. Examples: `outcome ~ trt` for continuous
  outcomes, `n_events ~ trt + offset(log(days))` for count outcomes, or
  `Surv(time, status) ~ trt` for survival models. The treatment variable
  (right-hand side) will be automatically extracted and converted to a
  numeric binary variable (0/1).

- unshrunk_terms_formula:

  A formula object or `NULL`. Formula specifying unshrunk terms for the
  `unshrunktermeffect` component. May include main effects and treatment
  interactions without regularization. Supports three syntaxes that can
  be combined: colon notation (e.g., `~ age + sex + trt:biomarker`),
  star notation (e.g., `~ age + trt*biomarker`), or random effects
  notation (e.g., `~ age + (trt || biomarker)`).

- shrunk_prognostic_formula:

  A formula object or `NULL`. Formula specifying prognostic main effects
  to be regularized in the `shprogeffect` component. These are
  covariates where shrinkage/regularization is desired. Should use
  `~ 0 + ...` syntax to apply one-hot encoding for symmetric
  regularization across all levels without privileging a reference
  group.

- shrunk_predictive_formula:

  A formula object or `NULL`. Formula specifying predictive terms
  (treatment interactions) to be regularized in the `shpredeffect`
  component. Supports three syntaxes that can be combined: colon
  notation (e.g., `~ 0 + trt:region`), star notation (e.g.,
  `~ 0 + trt*region`), or random effects notation (e.g.,
  `~ (0 + trt || region)`). Should use `~ 0 + ...` syntax for one-hot
  encoding.

- response_type:

  A character string. Type of outcome variable, one of `"binary"`,
  `"count"`, `"continuous"`, or `"survival"`. This determines the
  appropriate likelihood function and link function for the model.

- stratification_formula:

  A formula object or `NULL`. Formula specifying stratification variable
  (e.g., `~ strata_var`) for modeling baseline hazard (survival models)
  or distributional parameters (other models). For survival models,
  estimates separate baseline hazard functions (`bhaz`) for each
  stratum. For continuous models, models varying residual standard
  deviation (`sigma`). For count models, models varying overdispersion
  (`shape`).

## Value

`list` with six named elements:

- `formula`:

  `brmsformula` object containing the complete multi-part model
  specification

- `data`:

  `data.frame` with modified contrast coding (must be used for model
  fitting)

- `response_type`:

  `character(1)` indicating the outcome type (`"binary"`, `"count"`,
  `"continuous"`, or `"survival"`)

- `trt_var`:

  `character(1)` name of treatment variable extracted from response
  formula

- `stan_variable_names`:

  `list` of character vectors showing Stan parameter names for each
  model component

- `has_intercept`:

  `logical(1)` indicating whether the unshrunk formula includes an
  intercept (`TRUE`) or was specified with `~ 0 + ...` (`FALSE`)

## Details

This classification allows for applying differential shrinkage
(regularization) to different parts of the model. The function also
prepares the corresponding data using R's contrast coding system for
proper factor handling and interactions.

## Key Features

- **Multi-Part Formula Construction:** It generates a `brmsformula`
  object with up to three distinct linear components
  (`unshrunktermeffect`, `shprogeffect`, `shpredeffect`), which are
  combined in a non-linear model. This allows for assigning different
  priors to each component.

- **Automated Interaction Handling:** Predictive terms support three
  syntaxes that can be used together or separately: colon notation
  (e.g., `~ trt:subgroup`), star notation (e.g., `~ trt*subgroup`), or
  random effects notation (e.g., `~ (trt || subgroup)`). You can mix
  these syntaxes in the same model (e.g.,
  `~ trt:var1 + trt*var2 + (trt || var3)`). All variables involved in
  interactions are converted to factors with appropriate contrasts
  (reference coding for unshrunk, one-hot for shrunk), enabling proper
  model fitting and prediction on new data without manual dummy variable
  creation.

- **Hierarchical Integrity Check:** When a predictive term like
  `trt:subgroup` is specified, the function checks whether the
  corresponding prognostic effect `subgroup` is included in the model.
  If missing, it issues a warning message to alert you about the
  violation of the marginality principle, allowing you to decide whether
  to add the main effect. Note: Star notation (`*`) automatically
  includes main effects, so these variables are excluded from the
  marginality check.

- **Intelligent Defaults:** The main treatment variable is automatically
  added as an unshrunk prognostic term if not specified elsewhere. It
  also provides warnings and resolves overlaps if a term is accidentally
  specified in multiple categories.

## Survival Model Details (Robust Spline Knot Calculation)

For survival models (`response_type = "survival"`), the function
explicitly models the baseline hazard using B-splines via
`brms::bhaz()`. It implements a highly robust method for calculating
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

## Data Transformation and Contrast Coding

This function returns a modified `data.frame` that must be used in
subsequent calls to
[`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md),
not the original data. The transformations are:

- **Treatment variable:** Converted to numeric binary (0/1) to avoid
  multi-level factor interactions. The first level (alphabetically or
  numerically) becomes 0, the second becomes 1.

- **Factor covariates:** Automatic contrast coding based on shrinkage
  type:

  - **Shrunk terms:** One-hot encoding (all factor levels represented).
    For proper regularization, specify these using `~ 0 + var` to
    explicitly remove the intercept and treat all subgroups
    symmetrically without privileging a specific reference group
    (ensures the exchangeability assumption).

  - **Unshrunk terms:** Dummy encoding (reference level dropped). Use
    standard formula syntax `~ var` with intercept.

- **Custom contrasts:** If you have pre-specified contrasts via
  `contrasts(data$var) <- <matrix>` before calling this function, they
  will be preserved. Otherwise, defaults are applied based on shrinkage
  category.

## Stratification

The `stratification_formula_str` argument allows for estimating certain
parameters separately for different groups. Its behavior depends on the
`response_type`:

- `survival`: Estimates a separate baseline hazard function (`bhaz`) for
  each level of the stratification variable.

- `continuous`: Models the residual standard deviation (`sigma`) as
  varying across levels of the stratification variable.

- `count`: Models the overdispersion parameter (`shape`) as varying
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

  # 2a. Example with colon interaction syntax
  prepared_model <- prepare_formula_model(
    data = sim_data,
    response_formula = Surv(time, status) ~ trt,
    unshrunk_terms_formula = ~ age + subgroup,
    shrunk_predictive_formula = ~ trt:subgroup,
    response_type = "survival",
    stratification_formula = ~ region
  )

  # 2b. Alternatively, using pipe-pipe (||) syntax
  # prepared_model <- prepare_formula_model(
  #   data = sim_data,
  #   response_formula = Surv(time, status) ~ trt,
  #   unshrunk_terms_formula = ~ age,
  #   shrunk_predictive_formula = ~ (0 + trt || subgroup),
  #   response_type = "survival",
  #   stratification_formula = ~ region
  # )

  # 3. View the results
  print(prepared_model$formula)
  print(head(prepared_model$data))
}
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Note: Marginality principle not followed - interaction term 'subgroup_onehot' is used without its main effect. Consider adding 'subgroup_onehot' to prognostic terms for proper model hierarchy.
#> Warning: Formula 'shpredeffect' contains an intercept. For proper regularization/interpretation, consider removing it by adding '~ 0 + ...' or '~ -1 + ...' to your input formula.
#> time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unshrunktermeffect + shpredeffect 
#> unshrunktermeffect ~ 0 + age + subgroup + trt
#> shpredeffect ~ trt:subgroup_onehot
#>   time status trt      age region subgroup subgroup_onehot
#> 1   29      0   1 57.87739      B       S2              S2
#> 2   79      1   1 57.69042      B       S1              S1
#> 3   41      1   1 53.32203      B       S3              S3
#> 4   88      0   0 39.91623      B       S3              S3
#> 5   94      1   1 48.80547      A       S3              S3
#> 6    6      1   1 47.19605      A       S2              S2
```
