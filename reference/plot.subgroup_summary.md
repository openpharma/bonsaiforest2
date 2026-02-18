# Plot Marginal Subgroup Treatment Effects

Creates a forest plot from a `subgroup_summary` object or list of
objects.

## Usage

``` r
# S3 method for class 'subgroup_summary'
plot(x, x_lab = NULL, title = NULL, ...)
```

## Arguments

- x:

  `subgroup_summary` or `list`. Object created by
  [`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md),
  or a named list of such objects for comparison plots.

- x_lab:

  `character(1)` or `NULL`. Custom label for x-axis.

- title:

  `character(1)` or `NULL`. Custom title for plot.

- ...:

  Additional arguments (currently unused).

## Value

`ggplot`. Forest plot visualization of subgroup treatment effects.

## Details

This function creates forest plots for subgroup treatment effects. It
supports two modes:

**Single Model Plot**: When `x` is a single `subgroup_summary` object,
creates a traditional forest plot with estimates and confidence
intervals displayed as text.

**Multiple Model Comparison**: When `x` is a named list of
`subgroup_summary` objects, creates a comparative plot where different
models are distinguished by colors and shapes. This is useful for
comparing one-way models vs global models, different prior
specifications, or sensitivity analyses.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single model plot
summary1 <- summary_subgroup_effects(brms_fit = model1)
plot(summary1, title = "Subgroup Effects")

# Multiple model comparison
summary2 <- summary_subgroup_effects(brms_fit = model2)
comparison <- list(
  "Model 1" = summary1,
  "Model 2" = summary2
)
plot(comparison, title = "Model Comparison")
} # }
```
