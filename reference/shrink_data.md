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

  Time-to-event outcome (months)

- event_yn:

  Event indicator (0 = censored, 1 = event)

- fup_duration:

  Follow-up duration in months

- count:

  Poisson count outcome

## Generation process

The dataset was generated with `set.seed(16)` using the following steps:

**Covariates**

Treatment (`trt`) is assigned using permuted blocks of size 4, yielding
a balanced binary indicator (0/1) across the 500 subjects. Three latent
standard-normal variables (`z1`, `z2`, `z3`) are drawn independently,
then `z2` is replaced by `0.4 * z1 + sqrt(1 - 0.4^2) * z2` to induce a
correlation of 0.4 between `z1` and `z2`. The three subgroup covariates
are obtained by discretising these latent variables:

- `x1`: cut of `z1` at 0.5 → binary levels `a` (≈69\\ and `b` (≈31\\

- `x2`: cut of `z2` at −0.5 and 0.5 → three levels `a`, `b`, `c`.
  Because `z2` is correlated with `z1`, `x2` is also moderately
  correlated with `x1`.

- `x3`: cut of `z3` at −0.66, 0, and 0.66 → four approximately
  equal-sized levels `a`, `b`, `c`, `d`.

**Continuous outcome (`y`) and binary response (`response`)**

\$\$y = 5 + 0.5 \cdot \mathbf{1}\[x_1 = a\] + 0.5 \cdot \mathbf{1}\[x_3
= b\] + 0.1 \cdot \mathbf{1}\[\text{trt} = 1\] + 0.4 \cdot
\mathbf{1}\[x_1 = a, \text{trt} = 1\] + \varepsilon, \quad \varepsilon
\sim N(0, 1)\$\$

This gives a population-level treatment effect of 0.1, which is
amplified to 0.5 within the `x1 = a` subgroup (treatment-by-`x1`
interaction of 0.4). `x3 = b` acts as a prognostic covariate only. The
binary outcome `response` is the indicator `y > 6`.

**Time-to-event outcome (`tt_event`, `event_yn`)**

Event times are drawn from an exponential distribution with rate
\$\$\lambda = \exp\\\big(-3.5 + \log(0.5)\cdot\mathbf{1}\[x_1=a\] +
\log(1.5)\cdot\mathbf{1}\[x_3=b\] +
\log(0.9)\cdot\mathbf{1}\[\text{trt}=1\] +
\log(0.5)\cdot\mathbf{1}\[x_1=a,\\\text{trt}=1\]\big)\$\$

so the hazard ratios for treatment are 0.9 in the `x1 = b` subgroup and
`0.9 × 0.5 = 0.45` in the `x1 = a` subgroup. Censoring combines
administrative censoring (uniform accrual over 24 months, 60-month study
duration) and loss to follow-up (exponential with an annual rate of 3\\

**Poisson count outcome (`count`, `fup_duration`)**

Individual follow-up is capped at 24 months or loss to follow-up,
whichever comes first. Counts are drawn from a Poisson distribution with
mean equal to the individual rate times follow-up duration, where the
rate is \$\$\mu = \exp\\\big(-3 + \log(0.8)\cdot\mathbf{1}\[x_1=a\] +
\log(0.6)\cdot\mathbf{1}\[\text{trt}=1\]\big)\$\$

There is no treatment-by-subgroup interaction in this outcome; the rate
ratio for treatment is 0.6 across all subgroups.

The complete generation code is in `data-raw/shrink_data.R`.
