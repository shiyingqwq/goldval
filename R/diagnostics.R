#' Diagnose gold-standard verification structure
#'
#' @param object A `goldval_object`.
#' @param verification_formula Formula for the verification model.
#' @param risk_deciles Number of prediction-risk groups.
#'
#' @return A `goldval_diagnostics` object.
#'
#' @details Diagnostic flags use heuristic software thresholds. They are not
#' validated minimum sample-size requirements. The proxy-positive verification
#' odds ratio is descriptive and is not interpreted causally.
#' @export
diagnose_goldval <- function(object,
                             verification_formula = verified ~ proxy_outcome + qlogis_pred,
                             risk_deciles = 10) {
  if (!inherits(object, "goldval_object")) {
    stop("`object` must be a goldval_object.", call. = FALSE)
  }
  if (!is.numeric(risk_deciles) || length(risk_deciles) != 1 || risk_deciles < 2) {
    stop("`risk_deciles` must be a single integer >= 2.", call. = FALSE)
  }

  dat <- object$data
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  dat$risk_decile <- risk_decile(dat$pred, as.integer(risk_deciles))

  fit_info <- fit_verification_model(dat, verification_formula)
  pi_hat <- predict_verification_probability(fit_info$fit, dat)
  weights <- if (all(is.na(pi_hat))) rep(NA_real_, nrow(dat)) else 1 / pi_hat

  risk_table <- stats::aggregate(
    cbind(
      n = rep(1, nrow(dat)),
      reviewed_n = dat$verified,
      gold_event_n = ifelse(dat$verified == 1L, dat$gold_outcome, 0),
      proxy_positive_n = dat$proxy_outcome
    ) ~ risk_decile,
    data = dat,
    FUN = sum
  )
  risk_table$review_rate <- risk_table$reviewed_n / risk_table$n
  risk_table$proxy_positive_fraction <- risk_table$proxy_positive_n / risk_table$n
  risk_table$mean_pred <- as.numeric(tapply(dat$pred, dat$risk_decile, mean))[match(risk_table$risk_decile, sort(unique(dat$risk_decile)))]

  enrichment <- verification_enrichment(dat)
  positivity <- data.frame(
    model_status = fit_info$status,
    min_pi = safe_min(pi_hat),
    p01_pi = safe_quantile(pi_hat, 0.01),
    p05_pi = safe_quantile(pi_hat, 0.05),
    median_pi = safe_median(pi_hat),
    p95_weight = safe_quantile(weights, 0.95),
    p99_weight = safe_quantile(weights, 0.99),
    max_weight = safe_max(weights)
  )

  flags <- diagnostic_flags(object, risk_table, positivity, fit_info$warning_n)
  out <- list(
    summary = data.frame(
      n = object$n,
      n_verified = object$n_verified,
      gold_events = sum(dat$gold_outcome[dat$verified == 1L]),
      proxy_prevalence = mean(dat$proxy_outcome),
      gold_prevalence = mean(dat$gold_outcome[dat$verified == 1L]),
      verification_rate = mean(dat$verified)
    ),
    enrichment = enrichment,
    positivity = positivity,
    risk_deciles = risk_table,
    verification_model = fit_info$fit,
    verification_model_status = fit_info$status,
    verification_model_warning_n = fit_info$warning_n,
    verification_model_warnings = fit_info$warning_messages,
    flags = flags
  )
  class(out) <- "goldval_diagnostics"
  out
}

qlogis_clip <- function(p, eps = 1e-6) {
  stats::qlogis(pmin(pmax(p, eps), 1 - eps))
}

risk_decile <- function(pred, n_groups) {
  ranks <- rank(pred, ties.method = "first")
  as.integer(ceiling(ranks / length(pred) * n_groups))
}

fit_verification_model <- function(dat, formula) {
  warning_n <- 0L
  warning_messages <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(formula, data = dat, family = stats::binomial()),
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
  list(fit = fit, warning_n = warning_n, warning_messages = warning_messages, status = "ok")
}

predict_verification_probability <- function(fit, dat) {
  if (is.null(fit)) return(rep(NA_real_, nrow(dat)))
  pi_hat <- as.numeric(stats::predict(fit, newdata = dat, type = "response"))
  pmin(pmax(pi_hat, 1e-4), 1 - 1e-4)
}

safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, prob, na.rm = TRUE))
}

verification_enrichment <- function(dat) {
  tab <- table(
    proxy_positive = factor(dat$proxy_outcome, levels = c(0, 1)),
    verified = factor(dat$verified, levels = c(0, 1))
  )
  a <- tab["1", "1"] + 0.5
  b <- tab["1", "0"] + 0.5
  c <- tab["0", "1"] + 0.5
  d <- tab["0", "0"] + 0.5
  data.frame(
    proxy_positive_review_or = as.numeric((a * d) / (b * c)),
    proxy_positive_review_rate = mean(dat$verified[dat$proxy_outcome == 1L]),
    proxy_negative_review_rate = mean(dat$verified[dat$proxy_outcome == 0L])
  )
}

diagnostic_flags <- function(object, risk_table, positivity, warning_n) {
  flags <- character(0)
  if (object$n_verified < 100) flags <- c(flags, "LOW_GOLD_N")
  if (sum(object$data$gold_outcome[object$data$verified == 1L]) < 20) flags <- c(flags, "LOW_GOLD_EVENTS")
  if (identical(positivity$model_status, "failed")) flags <- c(flags, "PI_MODEL_FAILED")
  if (!is.na(positivity$p05_pi) && positivity$p05_pi < 0.02) flags <- c(flags, "POSITIVITY_WEAK")
  if (!is.na(positivity$p99_weight) && positivity$p99_weight > 50) flags <- c(flags, "WEIGHT_UNSTABLE")
  if (any(risk_table$reviewed_n < 5)) flags <- c(flags, "TAIL_SPARSE")
  if (warning_n > 0) flags <- c(flags, "PI_MODEL_WARNING")
  if (length(flags) == 0) "OK" else unique(flags)
}

#' @export
print.goldval_diagnostics <- function(x, ...) {
  cat("Goldval verification diagnostics\n")
  cat("Gold sample size:", x$summary$n_verified, "\n")
  cat("Gold prevalence:", round(x$summary$gold_prevalence, 3), "\n")
  cat("Verification enrichment OR(proxy positive):", round(x$enrichment$proxy_positive_review_or, 3), "\n")
  cat("Lowest estimated verification probability:", round(x$positivity$min_pi, 3), "\n")
  cat("1st percentile verification probability:", round(x$positivity$p01_pi, 3), "\n")
  cat("5th percentile verification probability:", round(x$positivity$p05_pi, 3), "\n")
  cat("Max inverse probability weight:", round(x$positivity$max_weight, 1), "\n")
  if (x$verification_model_warning_n > 0) {
    cat("Verification model warnings:", paste(utils::head(x$verification_model_warnings, 3), collapse = " | "), "\n")
  }
  cat("Flags:", paste(x$flags, collapse = ", "), "\n")
  invisible(x)
}
