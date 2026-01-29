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

  A `brmsfit` object. Fitted model object from
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)
  or
  [`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md).
  Must contain the necessary attributes for extracting treatment
  variable and response type information.

- trt_var:

  A character string or `NULL`. Treatment variable name. If `NULL`,
  automatically extracted from model attributes (set by
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
  Must be a binary variable coded as 0/1 in the dataset.

- data:

  A data frame or `NULL`. Dataset used for model fitting. If `NULL`,
  automatically extracted from model attributes (set by
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
  This dataset is used for generating counterfactual predictions.

- subgroup_vars:

  A character vector or `"auto"`. Subgroup variable names for which to
  estimate treatment effects. If `"auto"` (default), automatically
  detects treatment interaction terms (colon syntax) and random effect
  grouping factors (pipe syntax) from all formula components
  (`unshrunktermeffect`, `shprogeffect`, `shpredeffect`).

- response_type:

  A character string or `NULL`. Outcome type, one of `"binary"`,
  `"count"`, `"continuous"`, or `"survival"`. If `NULL`, automatically
  extracted from model attributes (set by
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
  This determines the appropriate scale for effect estimation.

- ndraws:

  An integer or `NULL`. Number of posterior draws to use for estimation.
  If `NULL` (default), all available posterior draws are used. Reducing
  this can speed up computation at the cost of precision.

## Value

`list` with two named elements:

- `estimates`:

  `tibble` where each row represents a subgroup (or "Overall" effect),
  with columns for `Subgroup`, `Median`, `CI_Lower`, and `CI_Upper` from
  posterior distribution

- `draws`:

  `data.frame` containing posterior draws for each subgroup

## Details

This post-processing function estimates marginal treatment effects for
specified subgroups from a fitted `brmsfit` object.
