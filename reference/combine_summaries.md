# Combine Multiple Subgroup Summaries for Comparison Plotting

Combines estimates from multiple subgroup summary objects into a single
summary object that can be plotted to compare models.

## Usage

``` r
combine_summaries(summary_list)
```

## Arguments

- summary_list:

  Named list of `subgroup_summary` objects.

## Value

`subgroup_summary`. Combined object with a "Model" column in estimates.

## Examples

``` r
if (FALSE) { # \dontrun{
summary1 <- summary_subgroup_effects(brms_fit = model1)
summary2 <- summary_subgroup_effects(brms_fit = model2)

combined <- combine_summaries(list(
  "One-way model" = summary1,
  "Global" = summary2
))
plot(combined)
} # }
```
