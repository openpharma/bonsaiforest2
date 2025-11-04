#' Fit a Bayesian Hierarchical Model using brms
#'
#' This function serves as a wrapper to fit a Bayesian model using the `brms`
#' package. It uses a formula and data prepared by `prepare_formula_model` and
#' simplifies the process of assigning complex priors to the different non-linear
#' components of the model.
#'
#' @section Prior Specification:
#' Priors for prognostic and predictive effects should be provided as named lists.
#' This makes the code explicit and prevents errors. For example:
#' `prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 5)", intercept = "normal(0, 10)")`.
#' The function accepts prior definitions as character strings (e.g., `"normal(0, 1)"`)
#' or as full `brmsprior` objects, which is useful for defining complex
#' hierarchical or parameter-specific priors.
#'
#' @param prepared_model A list object returned from `prepare_formula_model()`.
#'   It must contain the elements `formula`, `data`, and `response_type`.
#' @param predictive_effect_priors A named list with elements `shrunk` and/or `unshrunk`
#'   containing the priors for predictive effects.
#' @param prognostic_effect_priors A named list with elements `shrunk`, `unshrunk`,
#'   and/or `intercept` containing the priors for prognostic effects.
#' @param stanvars An object created by `brms::stanvar()` to add custom Stan code.
#' @param ... Additional arguments passed directly to `brms::brm()`.
#'
#' @return A fitted `brmsfit` object.
#'
#' @importFrom checkmate assert_class assert_data_frame assert_choice assert_list assert_string
#' @importFrom brms bernoulli negbinomial cox set_prior empty_prior
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
#'   \donttest{
#'   fit <- fit_brms_model(
#'     prepared_model = prepared_model,
#'     # Note: intercept prior is not needed for survival models
#'     prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 2)"),
#'     predictive_effect_priors = list(shrunk = "horseshoe(1)"),
#'     chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
#'   )
#'
#'   print(fit)
#'  }
#' }
fit_brms_model <- function(prepared_model,
                           predictive_effect_priors = list(),
                           prognostic_effect_priors = list(),
                           stanvars = NULL,
                           ...) {

  # --- 1. Argument Validation (Validate-First) ---

  # 1a. Validate the container and its structure
  checkmate::assert_list(prepared_model, names = "named", .var.name = "prepared_model")
  checkmate::assert_names(
    names(prepared_model),
    must.include = c("formula", "data", "response_type"),
    .var.name = "prepared_model"
  )

  # 1b. Validate the contents of the container
  checkmate::assert_class(prepared_model$formula, "brmsformula")
  checkmate::assert_data_frame(prepared_model$data, min.rows = 1)
  checkmate::assert_choice(
    prepared_model$response_type,
    choices = c("binary", "count", "continuous", "survival")
  )
  checkmate::assert_list(predictive_effect_priors, names = "named", null.ok = TRUE)
  checkmate::assert_list(prognostic_effect_priors, names = "named", null.ok = TRUE)
  checkmate::assert_class(stanvars, "stanvars", null.ok = TRUE)

  # --- Unpack  ---
  formula <- prepared_model$formula
  data <- prepared_model$data
  response_type <- prepared_model$response_type

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

  # Define all possible prior components
  # KEY CHANGE: The Intercept is class 'b' with coef 'Intercept'
  prior_config <- list(
    # Shrunk Prognostic (b) - No intercept by definition ( ~ 0 + ...)
    list(nlpar = "shprogeffect", class = "b", coef = NULL,
         user_prior = prognostic_effect_priors$shrunk,
         default = "horseshoe(1)",
         label = "shrunk prognostic (b)"),

    # Unshrunk Prognostic (b) - NON-intercepts
    list(nlpar = "unprogeffect", class = "b", coef = NULL,
         user_prior = prognostic_effect_priors$unshrunk,
         default = "normal(0, 5)",
         label = "unshrunk prognostic (b)"),

    # Unshrunk Prognostic (Intercept)
    list(nlpar = "unprogeffect", class = "b", coef = "Intercept", # <-- CHANGED
         user_prior = prognostic_effect_priors$intercept,
         default = "normal(0, 5)",
         label = "prognostic intercept"),

    # Shrunk Predictive (b) - No intercept by definition
    list(nlpar = "shpredeffect", class = "b", coef = NULL,
         user_prior = predictive_effect_priors$shrunk,
         default = "horseshoe(1)",
         label = "shrunk predictive (b)"),

    # Unshrunk Predictive (b) - No intercept by definition
    list(nlpar = "unpredeffect", class = "b", coef = NULL,
         user_prior = predictive_effect_priors$unshrunk,
         default = "normal(0, 10)",
         label = "unshrunk predictive (b)")
  )

  defined_nlpars <- names(formula$pforms)
  default_messages <- c("Using default priors for unspecified effects:")
  priors_list <- list()

  for (conf in prior_config) {
    # Only add prior if the nlpar exists in the formula
    if (conf$nlpar %in% defined_nlpars) {

      # Special case: Intercept prior only relevant if nlpar is unprogeffect
      # AND the response type is NOT survival (which has no intercept)
      # CHANGED: Check conf$coef now, not conf$class
      if (!is.null(conf$coef) && conf$coef == "Intercept" && response_type == "survival") {
        next # Skip intercept prior for survival models
      }

      processed <- .process_and_retarget_prior(
        user_prior = conf$user_prior,
        target_nlpar = conf$nlpar,
        default_str = conf$default,
        target_class = conf$class,
        target_coef = conf$coef  # <-- PASSING NEW ARG
      )

      priors_list <- c(priors_list, list(processed$prior))

      if (processed$default_used) {
        default_messages <- c(default_messages, paste0("  - ", conf$label, ": ", conf$default))
      }
    }
  }

  # Print messages if any defaults were used
  if (length(default_messages) > 1) {
    message(paste(default_messages, collapse = "\n"))
  }

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
#' 1. If the prior is NULL, it applies a default prior string using `brms::set_prior`.
#' 2. If the prior is a character string, it converts it to a `brmsprior` object
#'    using `brms::set_prior`.
#' 3. If the prior is already a `brmsprior` object, it intelligently re-targets
#'    any priors defined without an `nlpar` to the specific `target_nlpar`,
#'    `target_class`, and `target_coef` required.
#'
#' @param user_prior A user-provided prior. Can be `NULL`, a character string,
#'  or a `brmsprior` object.
#' @param target_nlpar A character string specifying the non-linear parameter.
#' @param default_str A default prior string to use if `user_prior` is NULL.
#' @param target_class A character string specifying the class (e.g., "b").
#' @param target_coef A character string specifying the coef (e.g., "Intercept"), or NULL.
#'
#' @return A list containing `prior` (a `brmsprior` object) and
#'  `default_used` (a boolean).
#' @noRd
.process_and_retarget_prior <- function(user_prior, target_nlpar, default_str, target_class = NULL, target_coef = NULL) {

  # --- Assertions for helper function ---
  checkmate::assert_string(target_nlpar, min.chars = 1)
  checkmate::assert_string(default_str, min.chars = 1)
  checkmate::assert_string(target_class, null.ok = TRUE)
  checkmate::assert_string(target_coef, null.ok = TRUE)

  default_used <- FALSE

  if (is.null(user_prior)) {
    # Use the default if nothing is provided
    prior_to_use <- default_str
    default_used <- TRUE
  } else {
    prior_to_use <- user_prior
  }

  if (is.character(prior_to_use)) {
    # --- THIS BLOCK IS REVISED ---
    # Build a list of arguments, excluding NULLs
    args <- list(
      prior = prior_to_use,
      nlpar = target_nlpar
    )

    # Only add class and coef to the call list if they are not NULL
    if (!is.null(target_class)) {
      args$class <- target_class
    }
    if (!is.null(target_coef)) {
      args$coef <- target_coef
    }

    # Call set_prior dynamically with the constructed argument list
    return(list(
      prior = do.call(brms::set_prior, args),
      default_used = default_used
    ))
  }

  if (inherits(prior_to_use, "brmsprior")) {

    msg <- paste("Re-targeting 'brmsprior' object for nlpar:", target_nlpar)
    if (!is.null(target_class)) msg <- paste(msg, "class:", target_class)
    if (!is.null(target_coef)) msg <- paste(msg, "coef:", target_coef)
    message(msg)

    modified_prior <- prior_to_use

    # Find rows with no nlpar set. These are the ones to retarget.
    rows_to_target <- nchar(modified_prior$nlpar) == 0

    if (any(rows_to_target)) {
      modified_prior$nlpar[rows_to_target] <- target_nlpar

      # Only assign class/coef if they are not NULL
      if (!is.null(target_class)) {
        modified_prior$class[rows_to_target] <- target_class
      }
      if (!is.null(target_coef)) {
        modified_prior$coef[rows_to_target] <- target_coef
      }
    }
    return(list(prior = modified_prior, default_used = FALSE))
    # --- END REVISION ---
  }

  stop(paste("Prior for", target_nlpar, "must be NULL, a string, or a brmsprior object."), call. = FALSE)
}
