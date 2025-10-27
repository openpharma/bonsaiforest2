#' Fit a Bayesian Hierarchical Model using brms
#'
#' This function serves as a wrapper to fit a Bayesian model using the `brms`
#' package. It uses a formula and data prepared by `prepare_formula_model` and
#' simplifies the process of assigning complex priors to the different non-linear
#' components of the model.
#'
#' @section Prior Specification:
#' Priors for prognostic and predictive effects should be provided as **named lists**.
#' This makes the code explicit and prevents errors. For example:
#' `prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 5)")`.
#' The function accepts prior definitions as character strings (e.g., `"normal(0, 1)"`)
#' or as full `brmsprior` objects, which is useful for defining complex
#' hierarchical or parameter-specific priors.
#'
#' @param formula A `brmsformula` object from `prepare_formula_model`.
#' @param data A data.frame from `prepare_formula_model`.
#' @param response_type The type of the outcome variable. One of "binary", "count",
#'   "continuous", or "survival".
#' @param predictive_effect_priors A **named list** with elements `shrunk` and/or `unshrunk`
#'   containing the priors for predictive effects.
#' @param prognostic_effect_priors A **named list** with elements `shrunk` and/or `unshrunk`
#'   containing the priors for prognostic effects.
#' @param stanvars An object created by `brms::stanvar()` to add custom Stan code.
#' @param ... Additional arguments passed directly to `brms::brm()`.
#'
#' @return A fitted `brmsfit` object.
#'
#' @importFrom checkmate assert_class assert_data_frame assert_choice assert_list assert_string
#' @importFrom brms bernoulli negbinomial cox
#' @export
#'
#' @examples
#' if (require("brms") && require("survival")) {
#'   # 1. Create Sample Data
#'   set.seed(123)
#'   n <- 100
#'   sim_data <- data.frame(
#'     time = round(runif(n, 1, 100)),
#'     status = sample(0:1, n, replace = TRUE),
#'     trt = sample(0:1, n, replace = TRUE),
#'     age = rnorm(n, 50, 10),
#'     region = sample(c("A", "B"), n, replace = TRUE),
#'     subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
#'   )
#'   sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
#'   sim_data$region <- as.factor(sim_data$region)
#'   sim_data$subgroup <- as.factor(sim_data$subgroup)
#'
#'   # 2. Prepare the formula and data
#'   prepared_model <- prepare_formula_model(
#'     data = sim_data,
#'     response_formula_str = "Surv(time, status) ~ trt",
#'     shrunk_predictive_formula_str = "~ trt:subgroup",
#'     unshrunk_prognostic_formula_str = "~ age",
#'     shrunk_prognostic_formula_str = "~ region",
#'     response_type = "survival",
#'     stratification_formula_str = "~ region"
#'   )
#'
#'   # 3. Fit the model
#'   # We use \donttest{} because fitting a model can be slow
#'   # and is not suitable for automated CRAN checks.
#'   \donttest{
#'   fit <- fit_brms_model(
#'     formula = prepared_model$formula,
#'     data = prepared_model$data,
#'     response_type = "survival",
#'     prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 5)"),
#'     predictive_effect_priors = list(shrunk = "horseshoe(1)"),
#'     chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
#'   )
#'
#'   print(fit)
#'   }
#' }
fit_brms_model <- function(formula, data, response_type,
                           predictive_effect_priors = list(),
                           prognostic_effect_priors = list(),
                           stanvars = NULL,
                           ...) {

  # --- 1. Argument Validation (with checkmate) ---
  checkmate::assert_class(formula, "brmsformula")
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_choice(response_type,
                           choices = c("binary", "count", "continuous", "survival")
  )
  checkmate::assert_list(predictive_effect_priors, names = "named", null.ok = TRUE)
  checkmate::assert_list(prognostic_effect_priors, names = "named", null.ok = TRUE)
  checkmate::assert_class(stanvars, "stanvars", null.ok = TRUE) # Changed "stanvar" to "stanvars"

  # --- 2. Determine brms Family ---
  model_family <- switch(
    response_type,
    continuous = gaussian(),
    binary     = bernoulli(link = "logit"),
    count      = negbinomial(),
    survival   = cox()
  )
  if (is.null(model_family)) {
    stop("Invalid 'response_type' specified.", call. = FALSE)
  }

  # --- 3. Construct the Prior List  ---

  prior_config <- list(
    list(nlpar = "shprogeffect", user_prior = prognostic_effect_priors$shrunk,   default = "horseshoe(1)"),
    list(nlpar = "unprogeffect", user_prior = prognostic_effect_priors$unshrunk, default = "normal(0, 10)"),
    list(nlpar = "shpredeffect", user_prior = predictive_effect_priors$shrunk,   default = "horseshoe(1)"),
    list(nlpar = "unpredeffect", user_prior = predictive_effect_priors$unshrunk, default = "normal(0, 10)")
  )

  defined_nlpars <- names(formula$pforms)

  # Use lapply to build the list of prior objects, only for nlpars that exist
  priors_list <- lapply(prior_config, function(conf) {
    if (conf$nlpar %in% defined_nlpars) {
      .process_and_retarget_prior(
        user_prior = conf$user_prior,
        target_nlpar = conf$nlpar,
        default_str = conf$default
      )
    } else {
      NULL
    }
  })

  # Filter out NULLs (for nlpars that weren't in the formula)
  priors_list <- priors_list[!sapply(priors_list, is.null)]

  # Combine all prior objects into a single brmsprior object
  final_priors <- if (length(priors_list) > 0) Reduce("+", priors_list) else brms::empty_prior()


  # --- 4. Call brm() ---
  message("Fitting brms model...")
  model_fit <- brms::brm(
    formula = formula,
    data = data,
    family = model_family,
    prior = final_priors,
    stanvars = stanvars,
    ...
  )

  # --- 5. Return the fitted model ---
  return(model_fit)
}


#' Process and Re-target a brms Prior
#'
#' This internal helper function processes a user-provided prior. It handles three cases:
#' 1. If the prior is NULL, it applies a default prior string.
#' 2. If the prior is a character string, it converts it to a `brmsprior` object.
#' 3. If the prior is already a `brmsprior` object (potentially complex, for
#'    hierarchical models), it intelligently re-targets any priors defined for the
#'    general coefficient class 'b' to the specific non-linear parameter (`nlpar`)
#'    required by the model structure.
#'
#' @param user_prior A user-provided prior. Can be `NULL`, a character string,
#'   or a `brmsprior` object.
#' @param target_nlpar A character string specifying the non-linear parameter
#'   (e.g., "shprogeffect") to which the prior should apply.
#' @param default_str A default prior string (e.g., "normal(0, 10)") to use if
#'   `user_prior` is NULL.
#'
#' @return A `brmsprior` object correctly targeted to the specified `nlpar`.
#' @noRd
.process_and_retarget_prior <- function(user_prior, target_nlpar, default_str) {

  # --- Assertions for helper function ---
  checkmate::assert_string(target_nlpar, min.chars = 1)
  checkmate::assert_string(default_str, min.chars = 1)

  if (is.null(user_prior)) {
    # Use the default if nothing is provided
    return(brms::prior_string(default_str, nlpar = target_nlpar))
  }

  if (is.character(user_prior)) {
    # User provided a simple string
    return(brms::prior_string(user_prior, nlpar = target_nlpar))
  }

  if (inherits(user_prior, "brmsprior")) {
    # User provided a brmsprior object.
    # We need to re-target any priors on class 'b' to our target_nlpar.
    message(paste("Re-targeting 'brmsprior' object for nlpar:", target_nlpar))
    modified_prior <- user_prior

    # Find rows with class 'b' and no nlpar set, and retarget them.
    b_class_rows <- modified_prior$class == 'b' & nchar(modified_prior$nlpar) == 0
    if(any(b_class_rows)) {
      modified_prior$nlpar[b_class_rows] <- target_nlpar
      modified_prior$class[b_class_rows] <- ''
    }
    return(modified_prior)
  }

  stop(paste("Prior for", target_nlpar, "must be NULL, a string, or a brmsprior object."), call. = FALSE)
}
