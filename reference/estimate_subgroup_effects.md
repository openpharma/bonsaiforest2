# Estimate Marginal Subgroup Treatment Effects

The function uses a counterfactual, marginal approach based on the
posterior predictive distribution. It averages over all other covariates
to provide robust estimates of subgroup-specific effects.

## Usage

``` r
estimate_subgroup_effects(
  brms_fit,
  original_data,
  trt_var,
  subgroup_vars = "auto",
  response_type = c("continuous", "binary", "count", "survival"),
  ndraws = NULL
)
```

## Arguments

- brms_fit:

  A fitted \`brmsfit\` object from \`fit_brms_model()\` or
  \`run_brms_analysis()\`.

- original_data:

  A \`data.frame\` containing the data. While this parameter is called
  "original_data" for backward compatibility, the function internally
  uses the processed data from \`brms_fit\$data\` to ensure consistency
  with the model's factor coding and contrasts. This parameter is mainly
  used for validation purposes.

- trt_var:

  A character string specifying the name of the treatment variable.

- subgroup_vars:

  A character vector of subgroup variable names found in
  \`original_data\`. If set to \`"auto"\` (the default), the function
  attempts to automatically identify subgroup variables from the model's
  formula.

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival".

- ndraws:

  An integer specifying the number of posterior draws to use. If
  \`NULL\` (default), all available draws are used.

## Value

A \`data.frame\` (tibble) where each row corresponds to a subgroup (or
the "Overall" effect), providing the estimated marginal effect and
posterior summaries (e.g., \`mean\`, \`sd\`, \`q2.5\`, \`q97.5\`).

## Details

This post-processing function estimates marginal treatment effects for
specified subgroups from a fitted \`brmsfit\` object.
