# Create a Summary of Marginal Subgroup Treatment Effects

This function orchestrates calls to \`estimate_subgroup_effects\` to
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

  \`brmsfit\`. Fitted model object.

- trt_var:

  \`character(1)\` or \`NULL\`. Treatment variable name (coded 0/1). If
  \`NULL\`, extracted from model attributes (set by
  \`fit_brms_model()\`).

- response_type:

  \`character(1)\` or \`NULL\`. Outcome type: \`"binary"\`, \`"count"\`,
  \`"continuous"\`, or \`"survival"\`. If \`NULL\`, extracted from model
  attributes.

- subgroup_vars:

  \`character\` or \`"auto"\`. Subgrouping variable names. Defaults to
  \`"auto"\`, which automatically detects treatment interactions from
  all formula components (\`unshrunktermeffect\`, \`shprogeffect\`,
  \`shpredeffect\`). Cannot be \`NULL\`.

## Value

\`subgroup_summary\`. S3 object containing:

- \`estimates\`:

  \`tibble\` with subgroup-specific treatment effect estimates

- \`response_type\`:

  \`character(1)\` indicating outcome type

- \`ci_level\`:

  \`numeric(1)\` credible interval level (default 0.95)

- \`trt_var\`:

  \`character(1)\` treatment variable name
