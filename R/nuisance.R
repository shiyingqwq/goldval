fit_goldval_nuisance <- function(dat,
                                 outcome_formula = gold_outcome ~ proxy_outcome * qlogis_pred,
                                 verification_formula = verified ~ proxy_outcome + qlogis_pred,
                                 fit_outcome = TRUE,
                                 fit_verification = TRUE) {
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  chart <- dat[dat$verified == 1L, , drop = FALSE]

  outcome <- if (fit_outcome) {
    fit_outcome_model(chart, outcome_formula)
  } else {
    empty_model_fit("not_fit")
  }
  verification <- if (fit_verification) {
    fit_verification_model(dat, verification_formula)
  } else {
    empty_model_fit("not_fit")
  }

  list(
    outcome_model = outcome$fit,
    outcome_model_status = outcome$status,
    outcome_model_warning_n = outcome$warning_n,
    outcome_model_warnings = outcome$warning_messages,
    verification_model = verification$fit,
    verification_model_status = verification$status,
    verification_model_warning_n = verification$warning_n,
    verification_model_warnings = verification$warning_messages,
    q_hat = predict_outcome_probability(outcome$fit, dat),
    pi_hat = predict_verification_probability(verification$fit, dat)
  )
}

empty_model_fit <- function(status) {
  list(fit = NULL, warning_n = 0L, warning_messages = character(0), status = status)
}

fit_outcome_model <- function(chart, formula) {
  fit_glm_with_status(formula, chart)
}

fit_verification_model <- function(dat, formula) {
  fit_glm_with_status(formula, dat)
}

fit_glm_with_status <- function(formula, data) {
  warning_n <- 0L
  warning_messages <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(formula, data = data, family = stats::binomial()),
      warning = function(w) {
        warning_n <<- warning_n + 1L
        warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(fit = NULL, warning_n = warning_n, warning_messages = warning_messages, status = "failed"))
  }
  status <- if (warning_n > 0L) "warning" else "ok"
  list(fit = fit, warning_n = warning_n, warning_messages = warning_messages, status = status)
}

predict_outcome_probability <- function(fit, dat) {
  if (is.null(fit)) return(rep(NA_real_, nrow(dat)))
  pred <- suppressWarnings(stats::predict(fit, newdata = dat, type = "response"))
  pmin(pmax(as.numeric(pred), 1e-6), 1 - 1e-6)
}

predict_verification_probability <- function(fit, dat) {
  if (is.null(fit)) return(rep(NA_real_, nrow(dat)))
  pi_hat <- as.numeric(stats::predict(fit, newdata = dat, type = "response"))
  pmin(pmax(pi_hat, 1e-4), 1 - 1e-4)
}
