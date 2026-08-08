#' Estimate validation performance
#'
#' @param object A `goldval_object`.
#' @param method Character vector of methods. Supported methods are `"naive"`,
#'   `"gold_only"`, `"OR"`, and `"AIPW"`.
#' @param metrics Character vector of metrics. Supported metrics are `"auc"`,
#'   `"brier"`, `"cal_intercept"`, and `"cal_slope"`.
#' @param outcome_formula Formula for the gold-outcome regression model used by
#'   corrected methods.
#' @param verification_formula Formula for the verification model used by AIPW.
#'
#' @return A `goldval_performance` object.
#' @export
performance <- function(object,
                        method = c("naive", "gold_only", "OR", "AIPW"),
                        metrics = c("auc", "brier", "cal_intercept", "cal_slope"),
                        outcome_formula = gold_outcome ~ proxy_outcome * qlogis_pred,
                        verification_formula = verified ~ proxy_outcome + qlogis_pred) {
  if (!inherits(object, "goldval_object")) {
    stop("`object` must be a goldval_object.", call. = FALSE)
  }
  method <- match.arg(method, several.ok = TRUE)
  metrics <- match.arg(metrics, c("auc", "brier", "cal_intercept", "cal_slope"), several.ok = TRUE)

  dat <- object$data
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  nuisance <- if (any(method %in% c("OR", "AIPW"))) {
    fit_goldval_nuisance(dat, outcome_formula, verification_formula)
  } else {
    empty_performance_nuisance(dat)
  }

  estimates <- lapply(method, function(m) {
    y <- performance_outcome(dat, m, nuisance)
    metric_values <- compute_performance_metrics(dat, y, metrics, m, nuisance)
    data.frame(
      method = m,
      metric = names(metric_values),
      estimate = as.numeric(metric_values),
      n_used = attr(metric_values, "n_used"),
      stringsAsFactors = FALSE
    )
  })
  estimates <- do.call(rbind, estimates)
  row.names(estimates) <- NULL

  out <- list(
    estimates = estimates,
    methods = method,
    metrics = metrics,
    nuisance = nuisance,
    call = match.call()
  )
  class(out) <- "goldval_performance"
  out
}

empty_performance_nuisance <- function(dat) {
  list(
    outcome_model = NULL,
    outcome_model_status = "not_fit",
    outcome_model_warning_n = 0L,
    outcome_model_warnings = character(0),
    verification_model = NULL,
    verification_model_status = "not_fit",
    verification_model_warning_n = 0L,
    verification_model_warnings = character(0),
    q_hat = rep(NA_real_, nrow(dat)),
    pi_hat = rep(NA_real_, nrow(dat))
  )
}

performance_outcome <- function(dat, method, nuisance) {
  if (method == "naive") {
    return(dat$proxy_outcome)
  }
  if (method == "gold_only") {
    y <- rep(NA_real_, nrow(dat))
    y[dat$verified == 1L] <- dat$gold_outcome[dat$verified == 1L]
    return(y)
  }
  if (method == "OR") {
    return(nuisance$q_hat)
  }
  if (method == "AIPW") {
    if (all(is.na(nuisance$pi_hat))) return(rep(NA_real_, nrow(dat)))
    residual <- rep(0, nrow(dat))
    reviewed <- dat$verified == 1L
    residual[reviewed] <- dat$gold_outcome[reviewed] - nuisance$q_hat[reviewed]
    return(nuisance$q_hat + dat$verified / nuisance$pi_hat * residual)
  }
  stop("Unknown method: ", method, call. = FALSE)
}

compute_performance_metrics <- function(dat, y, metrics, method, nuisance) {
  pred <- dat$pred
  ok <- !is.na(pred) & !is.na(y)
  pred_ok <- pred[ok]
  y_ok <- y[ok]
  values <- vapply(metrics, function(metric) {
    if (metric == "auc") return(weighted_auc(pred_ok, y_ok))
    if (metric == "brier") return(brier_estimate(dat, y, method, nuisance))
    if (metric %in% c("cal_intercept", "cal_slope")) {
      cal <- weak_calibration(pred_ok, y_ok)
      return(cal[[metric]])
    }
    NA_real_
  }, numeric(1))
  attr(values, "n_used") <- sum(ok)
  values
}

brier_estimate <- function(dat, y, method, nuisance) {
  if (method %in% c("naive", "gold_only")) {
    ok <- !is.na(y)
    return(mean((y[ok] - dat$pred[ok])^2))
  }
  q <- nuisance$q_hat
  conditional_loss <- q * (1 - dat$pred)^2 + (1 - q) * dat$pred^2
  if (method == "OR") {
    return(mean(conditional_loss))
  }
  if (method == "AIPW") {
    if (all(is.na(nuisance$pi_hat))) return(NA_real_)
    reviewed <- dat$verified == 1L
    residual <- rep(0, nrow(dat))
    observed_loss <- (dat$gold_outcome[reviewed] - dat$pred[reviewed])^2
    residual[reviewed] <- observed_loss - conditional_loss[reviewed]
    return(mean(conditional_loss + dat$verified / nuisance$pi_hat * residual))
  }
  NA_real_
}

weighted_auc <- function(pred, y) {
  if (length(pred) < 2L) return(NA_real_)
  y <- pmin(pmax(y, 0), 1)
  case_weight <- sum(y)
  control_weight <- sum(1 - y)
  if (case_weight <= 0 || control_weight <= 0) return(NA_real_)

  ord <- order(pred)
  pred <- pred[ord]
  y <- y[ord]
  control <- 1 - y
  r <- rle(pred)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  contribution <- numeric(length(pred))
  for (i in seq_along(starts)) {
    idx <- starts[i]:ends[i]
    controls_below <- if (starts[i] == 1L) 0 else sum(control[seq_len(starts[i] - 1L)])
    contribution[idx] <- controls_below + 0.5 * sum(control[idx])
  }
  sum(y * contribution) / (case_weight * control_weight)
}

weak_calibration <- function(pred, y) {
  if (length(pred) < 10L || stats::var(y) == 0) {
    return(c(cal_intercept = NA_real_, cal_slope = NA_real_))
  }
  lp <- qlogis_clip(pred)
  slope_fit <- try(stats::glm(y ~ lp, family = stats::binomial()), silent = TRUE)
  intercept_fit <- try(stats::glm(y ~ 1, offset = lp, family = stats::binomial()), silent = TRUE)
  c(
    cal_intercept = if (inherits(intercept_fit, "try-error")) NA_real_ else unname(stats::coef(intercept_fit)[["(Intercept)"]]),
    cal_slope = if (inherits(slope_fit, "try-error")) NA_real_ else unname(stats::coef(slope_fit)[["lp"]])
  )
}

#' @export
print.goldval_performance <- function(x, ...) {
  cat("goldval performance estimates\n\n")
  print(x$estimates, row.names = FALSE)
  invisible(x)
}
