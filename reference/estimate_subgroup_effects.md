# Estimate Marginal Subgroup Treatment Effects

The function uses a counterfactual, marginal approach based on the
posterior predictive distribution. It averages over all other covariates
to provide robust estimates of subgroup-specific effects.

## Usage

``` r
estimate_subgroup_effects(
  brms_fit,
  trt_var = NULL,
  data = NULL,
  subgroup_vars = "auto",
  response_type = NULL,
  ndraws = NULL
)
```

## Arguments

- brms_fit:

  \`brmsfit\`. Fitted model object from \`fit_brms_model()\` or
  \`run_brms_analysis()\`.

- trt_var:

  \`character(1)\` or \`NULL\`. Treatment variable name. If \`NULL\`,
  extracted from model attributes (set by \`fit_brms_model()\`).

- data:

  \`data.frame\` or \`NULL\`. Dataset used for model fitting. If
  \`NULL\`, extracted from model attributes (set by
  \`fit_brms_model()\`).

- subgroup_vars:

  \`character\` or \`"auto"\`. Subgroup variable names. If \`"auto"\`,
  automatically detects treatment interaction terms (colon syntax) and
  random effect grouping factors (pipe syntax) from all formula
  components.

- response_type:

  \`character(1)\` or \`NULL\`. Outcome type: \`"binary"\`, \`"count"\`,
  \`"continuous"\`, or \`"survival"\`. If \`NULL\`, extracted from model
  attributes (set by \`fit_brms_model()\`).

- ndraws:

  \`integer(1)\` or \`NULL\`. Number of posterior draws to use. If
  \`NULL\` (default), all available draws are used.

## Value

\`list\` with two named elements:

- \`estimates\`:

  \`tibble\` where each row represents a subgroup (or "Overall" effect),
  with columns for \`Subgroup\`, \`Median\`, \`CI_Lower\`, and
  \`CI_Upper\` from posterior distribution

- \`draws\`:

  \`data.frame\` containing posterior draws for each subgroup

## Details

This post-processing function estimates marginal treatment effects for
specified subgroups from a fitted \`brmsfit\` object.
