# Create a Summary of Marginal Subgroup Treatment Effects

This function orchestrates calls to \`estimate_subgroup_effects\` to
generate a summary table for specific subgroups.

## Usage

``` r
summary_subgroup_effects(
  brms_fit,
  original_data,
  trt_var,
  response_type = c("binary", "count", "continuous", "survival"),
  subgroup_vars = "auto"
)
```

## Arguments

- brms_fit:

  A fitted \`brmsfit\` object.

- original_data:

  The original \`data.frame\` used for model fitting.

- trt_var:

  The name of the treatment variable (coded 0/1).

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival".

- subgroup_vars:

  A character vector of subgrouping variable names. Defaults to "auto".
  Cannot be NULL.

## Value

An object of class \`subgroup_summary\`.
