## code to prepare `shrink_data` dataset goes here

library(dplyr)

set.seed(16)

n <- 500

# simulate covariates
df <- data.frame(
  id = 1:500,
  trt = as.vector(replicate(n/4, sample(rep(c(0, 1), each = 2)))),
  z1 = rnorm(n),
  z2 = rnorm(n),
  z3 = rnorm(n)
) |>
  mutate(
    z2 = 0.4*z1+sqrt(1-0.4^2)*z2,

    x1 = cut(z1, c(-Inf, 0.5, Inf), labels = c("a","b")),
    x2 = cut(z2, c(-Inf, -0.5, 0.5, Inf), labels = c("a","b", "c")),
    x3 = cut(z3, c(-Inf, -0.66, 0, 0.66, Inf), labels = c("a","b", "c", "d"))
  )

# add continuous and binary outcomes
df <- df |>
  mutate(
    y = 5+ 0.5*(x1=="a") + 0.5*(x3=="b") +
      0.1*(trt==1) + 0.4*((x1=="a")&(trt==1)) + rnorm(n),
    response = as.numeric(y>6)
  )

# add TTE outcome
df <- df |>
  mutate(
    time = rexp(
      n, rate = exp( -3.5 + log(0.5)*(x1=="a") + log(1.5)*(x3=="b") +
                       log(0.9)*(trt==1) + log(0.5)*((x1=="a")&(trt==1)))),
    # censor assuming 24 month uniform accrual, 60 month study duration,
    # and 3% drop-out per 12 months
    tt_entry = seq(0, 24, length = n),
    tt_admin_censor = 60-tt_entry,
    tt_dropout = rexp(n, rate=-log(0.97) / 12),
    tt_cens = pmin(tt_admin_censor, tt_dropout),
    tt_event = pmin(time, tt_cens),
    event_yn = as.numeric(time<tt_cens)
  )

# add Poisson count/rate outcome
df <- df |>
  mutate(
    fup_duration = pmin(24, tt_dropout),
    poi_rate = exp(-3 + log(0.8)*(x1=="a") + log(0.6)*(trt==1)), # Base lambda
    count = rpois(n, lambda = poi_rate * fup_duration)
  )

# remove unnecessary variables
shrink_data <- df |>
  dplyr::select(-c(z1,z2,z3, time, tt_entry, tt_admin_censor, tt_dropout, tt_cens, poi_rate))

usethis::use_data(shrink_data, overwrite = TRUE)

