#' Simulated data for shrinkage analysis
#'
#' A simulated dataset with 500 observations designed for demonstrating
#' subgroup shrinkage methods across multiple outcome types.
#'
#' @format A data frame with 500 rows and 11 columns:
#' \describe{
#'   \item{id}{Subject identifier (1-500)}
#'   \item{trt}{Binary treatment indicator (0 or 1)}
#'   \item{x1}{Categorical covariate with 2 levels (a, b)}
#'   \item{x2}{Categorical covariate with 3 levels (a, b, c)}
#'   \item{x3}{Categorical covariate with 4 levels (a, b, c, d)}
#'   \item{y}{Continuous outcome}
#'   \item{response}{Binary response (0 or 1)}
#'   \item{tt_event}{Time-to-event outcome (months)}
#'   \item{event_yn}{Event indicator (0 = censored, 1 = event)}
#'   \item{fup_duration}{Follow-up duration in months}
#'   \item{count}{Poisson count outcome}
#' }
#'
#' @section Generation process:
#'
#' The dataset was generated with \code{set.seed(16)} using the following steps:
#'
#' \strong{Covariates}
#'
#' Treatment (\code{trt}) is assigned using permuted blocks of size 4, yielding
#' a balanced binary indicator (0/1) across the 500 subjects. Three latent
#' standard-normal variables (\code{z1}, \code{z2}, \code{z3}) are drawn
#' independently, then \code{z2} is replaced by
#' \code{0.4 * z1 + sqrt(1 - 0.4^2) * z2} to induce a correlation of 0.4
#' between \code{z1} and \code{z2}. The three subgroup covariates are obtained
#' by discretising these latent variables:
#' \itemize{
#'   \item \code{x1}: cut of \code{z1} at 0.5 → binary levels \code{a} (≈69\%)
#'     and \code{b} (≈31\%).
#'   \item \code{x2}: cut of \code{z2} at −0.5 and 0.5 → three levels
#'     \code{a}, \code{b}, \code{c}. Because \code{z2} is correlated with
#'     \code{z1}, \code{x2} is also moderately correlated with \code{x1}.
#'   \item \code{x3}: cut of \code{z3} at −0.66, 0, and 0.66 → four
#'     approximately equal-sized levels \code{a}, \code{b}, \code{c}, \code{d}.
#' }
#'
#' \strong{Continuous outcome (\code{y}) and binary response (\code{response})}
#'
#' \deqn{y = 5 + 0.5 \cdot \mathbf{1}[x_1 = a] + 0.5 \cdot \mathbf{1}[x_3 = b]
#'   + 0.1 \cdot \mathbf{1}[\text{trt} = 1]
#'   + 0.4 \cdot \mathbf{1}[x_1 = a, \text{trt} = 1] + \varepsilon,
#'   \quad \varepsilon \sim N(0, 1)}
#'
#' This gives a population-level treatment effect of 0.1, which is amplified to
#' 0.5 within the \code{x1 = a} subgroup (treatment-by-\code{x1} interaction
#' of 0.4). \code{x3 = b} acts as a prognostic covariate only. The binary
#' outcome \code{response} is the indicator \code{y > 6}.
#'
#' \strong{Time-to-event outcome (\code{tt_event}, \code{event_yn})}
#'
#' Event times are drawn from an exponential distribution with rate
#' \deqn{\lambda = \exp\!\big(-3.5
#'   + \log(0.5)\cdot\mathbf{1}[x_1=a]
#'   + \log(1.5)\cdot\mathbf{1}[x_3=b]
#'   + \log(0.9)\cdot\mathbf{1}[\text{trt}=1]
#'   + \log(0.5)\cdot\mathbf{1}[x_1=a,\,\text{trt}=1]\big)}
#'
#' so the hazard ratios for treatment are 0.9 in the \code{x1 = b} subgroup
#' and \code{0.9 × 0.5 = 0.45} in the \code{x1 = a} subgroup. Censoring
#' combines administrative censoring (uniform accrual over 24 months, 60-month
#' study duration) and loss to follow-up (exponential with an annual rate of
#' 3\%). The observed time \code{tt_event = min(event time, censoring time)}.
#'
#' \strong{Poisson count outcome (\code{count}, \code{fup_duration})}
#'
#' Individual follow-up is capped at 24 months or loss to follow-up, whichever
#' comes first. Counts are drawn from a Poisson distribution with mean equal to
#' the individual rate times follow-up duration, where the rate is
#' \deqn{\mu = \exp\!\big(-3 + \log(0.8)\cdot\mathbf{1}[x_1=a]
#'   + \log(0.6)\cdot\mathbf{1}[\text{trt}=1]\big)}
#'
#' There is no treatment-by-subgroup interaction in this outcome; the rate
#' ratio for treatment is 0.6 across all subgroups.
#'
#' The complete generation code is in \code{data-raw/shrink_data.R}.
"shrink_data"
