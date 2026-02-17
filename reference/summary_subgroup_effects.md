# Create a Summary of Marginal Subgroup Treatment Effects

This function orchestrates calls to `estimate_subgroup_effects` to
generate a summary table for specific subgroups.

## Usage

``` r
summary_subgroup_effects(
  brms_fit,
  trt_var = NULL,
  response_type = NULL,
  subgroup_vars = "auto"
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

  A character string or `NULL`. Treatment variable name (coded 0/1). If
  `NULL`, automatically extracted from model attributes (set by
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
  This should match the treatment variable used in model fitting.

- response_type:

  A character string or `NULL`. Outcome type, one of `"binary"`,
  `"count"`, `"continuous"`, or `"survival"`. If `NULL`, automatically
  extracted from model attributes (set by
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
  This determines the appropriate scale for summarizing effects.

- subgroup_vars:

  A character vector or `"auto"`. Subgrouping variable names for which
  to generate summaries. Defaults to `"auto"`, which automatically
  detects treatment interactions from all formula components
  (`unshrunktermeffect`, `shprogeffect`, `shpredeffect`). Cannot be
  `NULL`.

## Value

`subgroup_summary`. S3 object containing:

- `estimates`:

  `tibble` where each row represents a subgroup, with columns for
  `Subgroup`, `Median`, `CI_Lower`, and `CI_Upper` from posterior
  distribution

- `response_type`:

  `character(1)` indicating outcome type

- `ci_level`:

  `numeric(1)` credible interval level (default 0.95)

- `trt_var`:

  `character(1)` treatment variable name
