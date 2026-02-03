#' Fit a Bayesian Hierarchical Model using brms
#'
#' This function serves as a wrapper to fit a Bayesian model using the `brms`
#' package. It uses a formula and data prepared by `prepare_formula_model` and
#' simplifies the process of assigning complex priors to the different non-linear
#' components of the model.
#'
#' @section Prior Specification:
#' Priors are specified separately for each model component using dedicated parameters.
#' This provides explicit control and prevents errors. The function accepts prior
#' definitions as character strings (e.g., `"normal(0, 2.5 * sigma_ref)"`) or as
#' `brmsprior` objects for complex hierarchical or parameter-specific priors.
#'
#' **Available prior parameters:**
#' - `intercept_prior`: Prior for the intercept in the unshrunk component
#' - `unshrunk_prior`: Prior for all unshrunk terms (main effects you trust)
#' - `shrunk_prognostic_prior`: Prior for shrunk prognostic effects (regularized)
#' - `shrunk_predictive_prior`: Prior for shrunk predictive effects (regularized)
#'
#' **Using sigma_ref in priors:**
#' You can reference `sigma_ref` in prior strings, and it will be automatically
#' substituted. For example: `"normal(0, 2.5 * sigma_ref)"` becomes
#' `"normal(0, 8)"` when sigma_ref = 3.2.
#'
#' **Default Priors by Response Type and Model Structure:**
#'
#' If you don't specify priors, the function uses sensible defaults that adapt to
#' the model structure and response type:
#'
#' **For models with fixed effects (GLOBAL) (colon `:` syntax):**
#'
#' *Continuous outcomes:*
#' - `intercept_prior`: `normal(outcome_mean, 2.5 * sigma_ref)` (weakly informative, centered at outcome mean)
#' - `unshrunk_prior`: `normal(0, 2.5 * sigma_ref)` (weakly informative)
#' - `shrunk_prognostic_prior`: `horseshoe(1)` (strong regularization)
#' - `shrunk_predictive_prior`: `horseshoe(1)` (strong regularization)
#' - `sigma`: `student_t(3, 0, sigma_ref)` (half-t prior for residual SD)
#'
#' *Binary outcomes:*
#' - `intercept_prior`: `normal(0, 1.5)` (weakly informative on logit scale)
#' - `unshrunk_prior`: `normal(0, 1.5)`
#' - `shrunk_prognostic_prior`: `horseshoe(1)`
#' - `shrunk_predictive_prior`: `horseshoe(1)`
#'
#' *Count outcomes:*
#' - `intercept_prior`: `normal(0, 2)` (weakly informative on log scale)
#' - `unshrunk_prior`: `normal(0, 2)`
#' - `shrunk_prognostic_prior`: `horseshoe(1)`
#' - `shrunk_predictive_prior`: `horseshoe(1)`
#' - `shape`: Uses brms default prior for negative binomial shape parameter
#'
#' *Survival outcomes:*
#' - No intercept (Cox models don't have intercepts)
#' - `unshrunk_prior`: `normal(0, 1.5)` (weakly informative on log-hazard scale)
#' - `shrunk_prognostic_prior`: `horseshoe(1)`
#' - `shrunk_predictive_prior`: `horseshoe(1)`
#'
#' **For models with random effects (One-variable-at-a-time (OVAT)) (pipe-pipe `||` syntax):**
#'
#' Random effects in shrunk predictive terms (e.g., `~ (1 + trt || subgroup)`)
#' automatically receive `normal(0, sigma_ref)` priors on the standard deviation
#' scale for each coefficient and group combination, regardless of the
#' shrunk_predictive_prior setting. For example:
#' - `prior(normal(0, sigma_ref), class = "sd", coef = "Intercept", group = "subgroup")`
#' - `prior(normal(0, sigma_ref), class = "sd", coef = "trt", group = "subgroup")`
#'
#' **Example:**
#' ```
#' fit_brms_model(
#'   prepared_model = prepared,
#'   sigma_ref = 1,
#'   unshrunk_prior = "normal(0, 2.5 * sigma_ref)",
#'   shrunk_prognostic_prior = "horseshoe(scale_global = 0.5)",
#'   shrunk_predictive_prior = "horseshoe(scale_global = 0.1)"
#' )
#' ```
#'
#' @param prepared_model A list object. Object returned from `prepare_formula_model()` containing
#'   the elements: `formula` (a `brmsformula` object), `data` (a modified `data.frame` with
#'   appropriate contrast coding), `response_type` (a character string), and `trt_var` (a
#'   character string). The `trt_var` element is stored as an attribute on the fitted model
#'   for use by downstream functions like `estimate_subgroup_effects()`.
#' @param sigma_ref A numeric scalar. Reference scale for prior specification. Default is 1.
#'   For continuous outcomes, recommended values are: (1) Assumed standard deviation from trial
#'   protocol (preferred), or (2) `sd(outcome_variable)` if protocol value unavailable.
#'   For binary, survival and count outcomes, the default value of 1 is typically appropriate.
#'   Can be referenced in prior strings using the placeholder `sigma_ref`
#'   (e.g., `"normal(0, 2.5 * sigma_ref)"`).
#' @param intercept_prior A character string, `brmsprior` object, or `NULL`. Prior specification
#'   for the model intercept in the `unshrunktermeffect` component. Supports `sigma_ref`
#'   placeholder substitution. Not used for survival models (Cox models have no intercept).
#'   Example: `"normal(0, 10)"`.
#' @param unshrunk_prior A character string, `brmsprior` object, or `NULL`. Prior specification
#'   for unshrunk terms in the `unshrunktermeffect` component (excludes intercept). These are
#'   main effects and interactions for which no regularization is desired.
#' @param shrunk_prognostic_prior A character string, `brmsprior` object, or `NULL`. Prior
#'   specification for regularized prognostic effects in the `shprogeffect` component. These
#'   are main effects for which strong shrinkage/regularization is desired.
#' @param shrunk_predictive_prior A character string, `brmsprior` object, or `NULL`. Prior
#'   specification for regularized predictive effects (treatment interactions) in the
#'   `shpredeffect` component. These are treatment interactions for which strong
#'   shrinkage/regularization is desired.
#' @param stanvars A `stanvars` object or `NULL`. Custom Stan code created via `brms::stanvar()`
#'   for implementing hierarchical priors or other advanced Stan functionality. See the brms
#'   documentation for details on creating `stanvars` objects.
#' @param ... Additional arguments passed to `brms::brm()`. Common arguments include `chains`
#'   (number of MCMC chains), `iter` (number of iterations per chain), `cores` (number of
#'   CPU cores to use), `backend` (Stan backend), and `refresh` (progress update frequency).
#'
#' @return `brmsfit`. Fitted Bayesian model object with attributes:
#'   `response_type` (`character`), `model_data` (`data.frame`), and
#'   `trt_var` (`character`) for downstream analysis functions.
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
#'   # 2. Prepare the formula and data using colon syntax
#'   prepared_model <- prepare_formula_model(
#'     data = sim_data,
#'     response_formula = Surv(time, status) ~ trt,
#'     shrunk_predictive_formula = ~ trt:subgroup,
#'     unshrunk_terms_formula = ~ age,
#'     shrunk_prognostic_formula = ~ region,
#'     response_type = "survival",
#'     stratification_formula = ~ region
#'   )
#'
#'   # 3. Fit the model
#'   \dontrun{
#'   # For survival models, typically use sigma_ref = 1
#'   # For continuous outcomes, use protocol sigma (preferred) or sd(outcome) (fallback)
#'   # Example: sigma_ref <- 12.5  # from protocol
#'   # OR: sigma_ref <- sd(sim_data$outcome)  # if protocol value unavailable
#'
#'   fit <- fit_brms_model(
#'     prepared_model = prepared_model,
#'     sigma_ref = 1,
#'     unshrunk_prior = "normal(0, 2 * sigma_ref)",
#'     shrunk_prognostic_prior = "horseshoe(scale_global = sigma_ref)",
#'     shrunk_predictive_prior = "horseshoe(scale_global = sigma_ref)",
#'     chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
#'   )
#'
#'   print(fit)
#'   }
#' }
fit_brms_model <- function(prepared_model,
                           sigma_ref = 1,
                           intercept_prior = NULL,
                           unshrunk_prior = NULL,
                           shrunk_prognostic_prior = NULL,
                           shrunk_predictive_prior = NULL,
                           stanvars = NULL,
                           ...) {

  # --- 1. Argument Validation ---

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

  # Extract trt_var from prepared_model
  trt_var <- prepared_model$trt_var
  if (!is.null(trt_var)) {
    checkmate::assert_string(trt_var, min.chars = 1)
    checkmate::assert_subset(trt_var, names(prepared_model$data))
    message("Using trt_var from prepared_model: ", trt_var)
  }

  # Validate sigma_ref
  checkmate::assert_number(sigma_ref, lower = 0, finite = TRUE,
                          .var.name = "sigma_ref")

  # Warn if using default sigma_ref = 1 for continuous outcomes
  if (prepared_model$response_type == "continuous" && sigma_ref == 1) {
    message("Note: Using default sigma_ref = 1 for continuous outcome. ",
            "Consider setting sigma_ref to the outcome standard deviation ",
            "from the trial protocol (preferred) or sd(outcome) for better prior specification.")
  }

  checkmate::assert_class(stanvars, "stanvars", null.ok = TRUE)

  # --- Unpack  ---
  formula <- prepared_model$formula
  data <- prepared_model$data
  response_type <- prepared_model$response_type
  has_intercept <- prepared_model$has_intercept  # Extract intercept information
  has_random_effects <- prepared_model$has_random_effects  # Extract random effects flag

  message("Using sigma_ref = ", round(sigma_ref, 3))

  # Calculate outcome mean for continuous models (used for intercept prior centering)
  outcome_mean <- NULL

  if (response_type == "continuous") {
    response_vars <- all.vars(formula$formula[[2]])
    response_var <- response_vars[1]
    if (response_var %in% names(data)) {
      outcome_mean <- mean(data[[response_var]], na.rm = TRUE)
      message("Computed outcome mean: ", round(outcome_mean, 3))
    }
  }

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

  # Process prior specifications with sigma_ref substitution
  # Enables user-friendly prior specification syntax: "normal(0, 2.5 * sigma_ref)"
  # The literal string "sigma_ref" is replaced with its numeric value before model compilation
  intercept_prior_processed <- .process_sigma_ref(intercept_prior, sigma_ref)
  unshrunk_prior_processed <- .process_sigma_ref(unshrunk_prior, sigma_ref)
  shrunk_prognostic_prior_processed <- .process_sigma_ref(shrunk_prognostic_prior, sigma_ref)
  shrunk_predictive_prior_processed <- .process_sigma_ref(shrunk_predictive_prior, sigma_ref)

  # Set default priors based on response type and model structure
  # For shrunk terms: use horseshoe(1) for fixed effects, normal(0, sigma_ref) for random effects
  if (response_type == "continuous") {
    intercept_mean <- if (!is.null(outcome_mean)) outcome_mean else 0
    default_intercept <- paste0("normal(", intercept_mean, ", ", 2.5 * sigma_ref, ")")
    default_unshrunk <- paste0("normal(0, ", 2.5 * sigma_ref, ")")
    default_shrunk_prog <- "horseshoe(1)"  # Regularizing horseshoe for fixed effects
    default_shrunk_pred <- "horseshoe(1)"  # For fixed effects; random effects handled separately with normal(0, sigma_ref)
  } else if (response_type == "binary") {
    # Binary outcomes: priors on logit scale
    default_intercept <- "normal(0, 1.5)"  # Weakly informative on logit scale
    default_unshrunk <- "normal(0, 1.5)"
    default_shrunk_prog <- "horseshoe(1)"  # Regularizing horseshoe for fixed effects
    default_shrunk_pred <- "horseshoe(1)"  # For fixed effects; random effects handled separately
  } else if (response_type == "count") {
    # Count outcomes: priors on log scale
    default_intercept <- "normal(0, 2)"  # Weakly informative on log scale
    default_unshrunk <- "normal(0, 2)"
    default_shrunk_prog <- "horseshoe(1)"  # Regularizing horseshoe for fixed effects
    default_shrunk_pred <- "horseshoe(1)"  # For fixed effects; random effects handled separately
  } else {
    # Survival outcomes: priors on log-hazard scale
    default_intercept <- "normal(0, 1.5)"  # Not used (Cox has no intercept)
    default_unshrunk <- "normal(0, 1.5)"
    default_shrunk_prog <- "horseshoe(1)"  # Regularizing horseshoe for fixed effects
    default_shrunk_pred <- "horseshoe(1)"  # For fixed effects; random effects handled separately
  }

  # Define prior configuration for three-component model architecture:
  # - unshrunktermeffect: Consolidated component for all unshrunk terms (no separate prog/pred)
  # - shprogeffect: Shrunk prognostic effects with strong regularization (e.g., horseshoe)
  # - shpredeffect: Shrunk predictive effects (treatment interactions) with strong regularization
  # This architecture enables differential prior specification across model components
  prior_config <- list(
    # Intercept for unshrunk terms
    list(nlpar = "unshrunktermeffect", class = "b", coef = "Intercept",
         user_prior = intercept_prior_processed,
         default = default_intercept,
         label = "intercept"),

    # Unshrunk terms (all non-intercept coefficients in unshrunktermeffect)
    # This includes both prognostic and predictive terms that are not regularized
    list(nlpar = "unshrunktermeffect", class = "b", coef = NULL,
         user_prior = unshrunk_prior_processed,
         default = default_unshrunk,
         label = "unshrunk terms"),

    # Shrunk Prognostic - No intercept (user specifies formulas with ~ 0 + ...)
    # Applied to prognostic biomarkers/covariates with strong regularization
    # Prognostic effects are main effects only (no treatment interactions)
    list(nlpar = "shprogeffect", class = "b", coef = NULL,
         user_prior = shrunk_prognostic_prior_processed,
         default = default_shrunk_prog,
         label = "shrunk prognostic"),

    # Shrunk Predictive - No intercept (user specifies formulas with ~ 0 + ...)
    # Applied to treatment interactions (effect modifiers) with strong regularization
    # Use 'sd' class for random effects, 'b' class for fixed effects
    list(nlpar = "shpredeffect",
         class = if (has_random_effects) "sd" else "b",
         coef = NULL,
         user_prior = shrunk_predictive_prior_processed,
         default = default_shrunk_pred,
         label = "shrunk predictive",
         is_random_effect = has_random_effects)
  )
  
  # Add mixed notation config ONLY if has_random_effects is TRUE
  # This handles the case: (1+trt||subvar1) + trt:subvar2
  # where trt:subvar2 are fixed effects that need horseshoe priors
  if (has_random_effects) {
    prior_config <- c(prior_config, list(
      list(nlpar = "shpredeffect",
           class = "b",
           coef = NULL,
           user_prior = shrunk_predictive_prior_processed,
           default = default_shrunk_pred,
           label = "shrunk predictive (fixed)",
           is_mixed_notation_fixed = TRUE)
    ))
  }

  # Extract random effects structure if present in shrunk predictive component
  random_effects_structure <- list()
  if (has_random_effects) {
    random_effects_structure <- .extract_random_effects_structure(formula, nlpar = "shpredeffect")

    if (length(random_effects_structure) > 0) {
      message("Detected random effects structure in shrunk predictive component:")
      for (re in random_effects_structure) {
        message("  - Group: ", re$group, ", Coefficients: ", paste(re$coefs, collapse = ", "))
      }
    }
  }

  defined_nlpars <- names(formula$pforms)
  default_messages <- c("Using default priors for unspecified effects:")
  priors_list <- list()

  for (conf in prior_config) {
    # Only add prior if the nlpar exists in the formula
    if (conf$nlpar %in% defined_nlpars) {

      # Skip intercept prior for survival models (Cox models have no intercept)
      # Intercept only exists in unshrunktermeffect for non-survival response types
      if (!is.null(conf$coef) && conf$coef == "Intercept" && response_type == "survival") {
        next
      }

      # Skip intercept prior if user specified ~ 0 + ... (no intercept)
      if (!is.null(conf$coef) && conf$coef == "Intercept" && !has_intercept) {
        next
      }

      # Special handling for random effects in shrunk predictive
      if (isTRUE(conf$is_random_effect) && length(random_effects_structure) > 0) {
        # For each random effect structure, create specific priors for each coefficient
        # Random effects always use normal(0, sigma_ref) as the default
        for (re in random_effects_structure) {
          for (coef_name in re$coefs) {
            # Create a prior for this specific coefficient in this specific group
            # Use user's prior if provided, otherwise use normal(0, sigma_ref)
            # If user provided a prior, it should be interpreted as applying to random effects
            # Otherwise, use the standard normal(0, sigma_ref) for random effect variances

            re_default_str <- paste0("normal(0, ", sigma_ref, ")")

            processed <- .process_and_retarget_prior(
              user_prior = conf$user_prior,  # Use user's prior if they provided one
              target_nlpar = conf$nlpar,
              default_str = re_default_str,  # Default for random effects is normal(0, sigma_ref)
              target_class = "sd",
              target_coef = coef_name,
              target_group = re$group
            )

            priors_list <- c(priors_list, list(processed$prior))

            if (processed$default_used) {
              default_messages <- c(default_messages,
                                  paste0("  - ", conf$label, " (sd, ", re$group, ", ", coef_name, "): ", re_default_str))
            }
          }
        }
      } else if (isTRUE(conf$is_mixed_notation_fixed)) {
        # For mixed notation: apply horseshoe(1) to fixed effects (b class)
        # even when random effects are present
        # This handles the case: (1+trt||subvar1) + trt:subvar2
        # where trt:subvar2 are fixed effects that need horseshoe priors
        # ONLY process this if we actually found random effects structure
        # AND if the user didn't explicitly set the prior to NULL (which means no fixed effects)
        if (length(random_effects_structure) > 0 && !is.null(conf$user_prior)) {
          processed <- .process_and_retarget_prior(
            user_prior = conf$user_prior,
            target_nlpar = conf$nlpar,
            default_str = conf$default,  # horseshoe(1) for fixed effects
            target_class = conf$class,  # "b"
            target_coef = conf$coef
          )

          priors_list <- c(priors_list, list(processed$prior))

          if (processed$default_used) {
            default_messages <- c(default_messages, paste0("  - ", conf$label, ": ", conf$default))
          }
        }
      } else {
        # Standard (non-random effect, non-mixed notation) prior processing
        # Process user-provided prior or use default, then target it to the correct nlpar
        # This handles conversion from strings to brmsprior objects and re-targeting
        processed <- .process_and_retarget_prior(
          user_prior = conf$user_prior,
          target_nlpar = conf$nlpar,
          default_str = conf$default,
          target_class = conf$class,
          target_coef = conf$coef
        )

        priors_list <- c(priors_list, list(processed$prior))

        if (processed$default_used) {
          default_messages <- c(default_messages, paste0("  - ", conf$label, ": ", conf$default))
        }
      }
    }
  }

  # Print messages if any defaults were used
  if (length(default_messages) > 1) {
    message(paste(default_messages, collapse = "\n"))
  }

  # Add sigma prior only for continuous models without stratification
  # Count models use negative binomial with shape parameter, not sigma
  # Stratified continuous models have distributional formulas for sigma, so no global sigma prior
  if (response_type == "continuous" && !("sigma" %in% names(formula$pforms))) {
    # Use half-t prior for residual SD (weakly informative, heavy tails)
    sigma_prior_str <- paste0("student_t(3, 0, ", sigma_ref, ")")
    sigma_prior <- brms::set_prior(sigma_prior_str, class = "sigma")
    priors_list <- c(priors_list, list(sigma_prior))
    message("Adding sigma prior: ", sigma_prior_str)
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

  # --- 5. Attach metadata as attributes ---
  # Store response_type, trt_var, and data for downstream functions
  attr(model_fit, "response_type") <- response_type
  attr(model_fit, "model_data") <- data
  if (!is.null(trt_var)) {
    attr(model_fit, "trt_var") <- trt_var
  }

  # --- 6. Return the fitted model ---
  return(model_fit)
}


#' Process Prior String with sigma_ref Substitution
#'
#' This internal helper function processes a prior string or brmsprior object
#' and substitutes the literal string \"sigma_ref\" with the actual numeric value.
#' This allows users to write priors like \"normal(0, 2.5 * sigma_ref)\" which
#' will be evaluated to the actual prior string.
#'
#' @param prior A prior specification. Can be NULL, a character string, or a brmsprior object.
#' @param sigma_ref A numeric value to substitute for \"sigma_ref\" in the prior string.
#'
#' @return The prior with sigma_ref substituted, or NULL if prior is NULL.
#' @noRd
.process_sigma_ref <- function(prior, sigma_ref) {
  if (is.null(prior)) {
    return(NULL)
  }

  if (is.character(prior)) {
    # Replace "sigma_ref" with its numeric value using string substitution
    # Example: "normal(0, 2.5 * sigma_ref)" becomes "normal(0, 2.5 * 3.2)"
    prior_substituted <- gsub("sigma_ref", as.character(sigma_ref), prior)
    return(prior_substituted)
  }

  # If it's a brmsprior object, process each prior string in it
  if (inherits(prior, "brmsprior")) {
    # Process each prior string in the brmsprior object
    for (i in seq_len(nrow(prior))) {
      if (nchar(prior$prior[i]) > 0 && grepl("sigma_ref", prior$prior[i])) {
        # Replace sigma_ref with the actual value
        prior$prior[i] <- gsub("sigma_ref", as.character(sigma_ref), prior$prior[i])
      }
    }
    return(prior)
  }

  return(prior)
}


#' Process and Re-target a brms Prior to a Specific nlpar
#'
#' This internal helper function processes a user-provided prior and ensures it targets
#' the correct formula component (nlpar). It handles three cases:
#' 1. NULL prior: Uses the provided default prior string
#' 2. Character string: Converts to a `brmsprior` object targeting the specified nlpar
#' 3. Existing `brmsprior` object: Intelligently re-targets priors that don't have an
#'    nlpar specified, while preserving coefficient-specific priors that do.
#'
#' This allows users to provide generic priors that get automatically targeted to the
#' correct model component (unshrunktermeffect, shprogeffect, or shpredeffect).
#'
#' @param user_prior A user-provided prior. Can be `NULL`, a character string,
#'  or a `brmsprior` object.
#' @param target_nlpar A character string specifying the non-linear parameter.
#' @param default_str A default prior string to use if `user_prior` is NULL.
#' @param target_class A character string specifying the class (e.g., "b", "sd"), or NULL.
#' @param target_coef A character string specifying the coef (e.g., "Intercept"), or NULL.
#' @param target_group A character string specifying the group for random effects (e.g., "subvar"), or NULL.
#'
#' @return A list containing `prior` (a `brmsprior` object) and
#'  `default_used` (a boolean).
#' @noRd
.process_and_retarget_prior <- function(user_prior, target_nlpar, default_str,
                                        target_class = NULL, target_coef = NULL,
                                        target_group = NULL) {

  # --- Assertions for helper function ---
  checkmate::assert_string(target_nlpar, min.chars = 1)
  checkmate::assert_string(default_str, min.chars = 1)
  checkmate::assert_string(target_class, null.ok = TRUE)
  checkmate::assert_string(target_coef, null.ok = TRUE)
  checkmate::assert_string(target_group, null.ok = TRUE)

  default_used <- FALSE

  if (is.null(user_prior)) {
    # Use the default if nothing is provided
    prior_to_use <- default_str
    default_used <- TRUE
  } else {
    prior_to_use <- user_prior
  }

  if (is.character(prior_to_use)) {
    # Build a list of arguments, excluding NULLs
    args <- list(
      prior = prior_to_use,
      nlpar = target_nlpar
    )

    # Only add class, coef, and group to the call list if they are not NULL
    if (!is.null(target_class)) {
      args$class <- target_class
    }
    if (!is.null(target_coef)) {
      args$coef <- target_coef
    }
    if (!is.null(target_group)) {
      args$group <- target_group
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
    if (!is.null(target_group)) msg <- paste(msg, "group:", target_group)
    message(msg)

    modified_prior <- prior_to_use

    # Find rows with no nlpar set. These are the ones that need the nlpar added.
    rows_to_target <- nchar(modified_prior$nlpar) == 0

    if (any(rows_to_target)) {
      # Add nlpar to all rows that don't have it
      modified_prior$nlpar[rows_to_target] <- target_nlpar

      # For class/coef/group: only set them if the row doesn't already have a coef specified
      # This allows users to pass brmsprior objects with specific coefficient priors
      # that will be preserved

      # Find rows that are class-level priors (no coef specified)
      class_level_rows <- rows_to_target & (nchar(modified_prior$coef) == 0)

      if (any(class_level_rows)) {
        # Only assign class/coef/group to class-level priors, not coefficient-specific ones
        if (!is.null(target_class)) {
          modified_prior$class[class_level_rows] <- target_class
        }
        if (!is.null(target_coef)) {
          modified_prior$coef[class_level_rows] <- target_coef
        }
        if (!is.null(target_group)) {
          modified_prior$group[class_level_rows] <- target_group
        }
      }
    }
    return(list(prior = modified_prior, default_used = FALSE))
  }

  stop(paste("Prior for", target_nlpar, "must be NULL, a string, or a brmsprior object."), call. = FALSE)
}


#' Extract Random Effects Structure from brms Formula
#'
#' This internal helper function parses a brmsformula object to extract
#' random effects specifications from the shrunk predictive component.
#' It identifies patterns like (1+trt||subvar) and returns the coefficient
#' names and grouping variables.
#'
#' @param formula A brmsformula object
#' @param nlpar Character string specifying the non-linear parameter to check (default: "shpredeffect")
#'
#' @return A list with elements for each random effect:
#'   \item{group}{The grouping variable name (e.g., "subvar")}
#'   \item{coefs}{A character vector of coefficient names (e.g., c("Intercept", "trt"))}
#' @noRd
.extract_random_effects_structure <- function(formula, nlpar = "shpredeffect") {
  checkmate::assert_class(formula, "brmsformula")
  checkmate::assert_string(nlpar)

  # Check if the nlpar exists in the formula
  if (!nlpar %in% names(formula$pforms)) {
    return(list())
  }

  # Get the formula for the specified nlpar
  nlpar_formula <- formula$pforms[[nlpar]]

  # The nlpar_formula IS the formula object, not a list with $formula
  # Convert formula to string for pattern matching
  if (inherits(nlpar_formula, "formula")) {
    formula_str <- deparse(nlpar_formula, width.cutoff = 500L)
  } else if (!is.null(nlpar_formula$formula)) {
    formula_str <- deparse(nlpar_formula$formula, width.cutoff = 500L)
  } else {
    return(list())
  }

  formula_str <- paste(formula_str, collapse = " ")

  # Pattern to match random effects: (var1 + var2 || group) or (1 + var || group)
  # This regex captures the content before || and the group name
  pattern <- "\\(\\s*([^|]+?)\\s*\\|\\|\\s*([^)]+?)\\s*\\)"

  matches <- stringr::str_match_all(formula_str, pattern)[[1]]

  if (nrow(matches) == 0) {
    return(list())
  }

  result <- list()

  for (i in seq_len(nrow(matches))) {
    lhs <- stringr::str_squish(matches[i, 2])  # e.g., "1 + trt" or "0 + trt"
    group <- stringr::str_squish(matches[i, 3])  # e.g., "subvar"

    # Remove leading 0 + or -1 + if present
    lhs <- stringr::str_remove(lhs, "^(0|\\-1)\\s*\\+\\s*")
    lhs <- stringr::str_squish(lhs)

    # Split by + to get individual terms
    terms <- stringr::str_squish(stringr::str_split(lhs, "\\+")[[1]])

    # Convert "1" to "Intercept" for brms prior specification
    coefs <- sapply(terms, function(x) {
      if (x == "1") "Intercept" else x
    }, USE.NAMES = FALSE)

    result[[length(result) + 1]] <- list(
      group = group,
      coefs = coefs
    )
  }

  return(result)
}
