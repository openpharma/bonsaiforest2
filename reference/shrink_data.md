# Simulated data for shrinkage analysis

A simulated dataset with 500 observations designed for demonstrating
subgroup shrinkage methods across multiple outcome types.

## Usage

``` r
shrink_data
```

## Format

A data frame with 500 rows and 11 columns:

- id:

  Subject identifier (1-500)

- trt:

  Binary treatment indicator (0 or 1)

- x1:

  Categorical covariate with 2 levels (a, b)

- x2:

  Categorical covariate with 3 levels (a, b, c)

- x3:

  Categorical covariate with 4 levels (a, b, c, d)

- y:

  Continuous outcome

- response:

  Binary response (0 or 1)

- tt_event:

  Time-to-event outcome

- event_yn:

  Event indicator (0 = censored, 1 = event)

- fup_duration:

  Follow-up duration in months

- count:

  Poisson count outcome
