#' Estimate calibration summaries and flexible calibration curves
#'
#' @param object A `goldval_object`.
#' @param method Character vector of methods. Supported methods are `"naive"`,
#'   `"gold_only"`, `"OR"`, and `"AIPW"`.
#' @param type Calibration outputs to return: `"intercept"`, `"slope"`, and/or
#'   `"curve"`.
#' @param outcome_formula Formula for the gold-outcome regression model.
#' @param verification_formula Formula for the verification model.
#' @param curve_df Degrees of freedom for the natural spline calibration curve.
#' @param grid Optional prediction-risk grid.
#'
#' @return A `goldval_calibration` object.
#' @export
calibration <- function(object,
                        method = c("naive", "gold_only", "OR", "AIPW"),
                        type = c("intercept", "slope", "curve"),
                        outcome_formula = gold_outcome ~ proxy_outcome * qlogis_pred,
                        verification_formula = verified ~ proxy_outcome + qlogis_pred,
                        curve_df = 4,
                        grid = NULL) {
  if (!inherits(object, "goldval_object")) {
    stop("`object` must be a goldval_object.", call. = FALSE)
  }
  method <- match.arg(method, several.ok = TRUE)
  type <- match.arg(type, c("intercept", "slope", "curve"), several.ok = TRUE)
  if (!is.numeric(curve_df) || length(curve_df) != 1L || curve_df < 1) {
    stop("`curve_df` must be a single positive number.", call. = FALSE)
  }

  dat <- object$data
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  nuisance <- if (any(method %in% c("OR", "AIPW"))) {
    fit_goldval_nuisance(dat, outcome_formula, verification_formula)
  } else {
    empty_performance_nuisance(dat)
  }
  if (is.null(grid)) {
    grid <- seq(
      stats::quantile(dat$pred, 0.02),
      stats::quantile(dat$pred, 0.98),
      length.out = 100L
    )
  }

  parameters <- if (any(type %in% c("intercept", "slope"))) {
    do.call(rbind, lapply(method, function(m) {
      estimate_weak_calibration(dat, m, nuisance)
    }))
  } else {
    data.frame()
  }
  row.names(parameters) <- NULL

  curves <- if ("curve" %in% type) {
    do.call(rbind, lapply(method, function(m) {
      estimate_calibration_curve(dat, m, nuisance, grid, curve_df)
    }))
  } else {
    data.frame()
  }
  row.names(curves) <- NULL

  out <- list(
    parameters = parameters,
    curves = curves,
    nuisance = nuisance,
    settings = list(curve_df = curve_df, grid = grid),
    call = match.call()
  )
  class(out) <- "goldval_calibration"
  out
}

estimate_weak_calibration <- function(dat, method, nuisance) {
  objective <- function(theta) {
    score <- calibration_score(dat, method, nuisance, theta)
    sum(score^2)
  }
  fit <- try(stats::optim(c(0, 1), objective, method = "BFGS"), silent = TRUE)
  if (inherits(fit, "try-error") || fit$convergence != 0) {
    return(data.frame(method = method, cal_intercept = NA_real_, cal_slope = NA_real_, n_used = n_used_for_method(dat, method)))
  }
  data.frame(
    method = method,
    cal_intercept = fit$par[[1]],
    cal_slope = fit$par[[2]],
    n_used = n_used_for_method(dat, method)
  )
}

calibration_score <- function(dat, method, nuisance, theta) {
  mu <- stats::plogis(theta[[1]] + theta[[2]] * qlogis_clip(dat$pred))
  x <- cbind(1, qlogis_clip(dat$pred))

  if (method == "naive") {
    return(colMeans(x * (dat$proxy_outcome - mu)))
  }
  if (method == "gold_only") {
    reviewed <- dat$verified == 1L
    return(colMeans(x[reviewed, , drop = FALSE] * (dat$gold_outcome[reviewed] - mu[reviewed])))
  }
  if (method == "OR") {
    if (all(is.na(nuisance$q_hat))) return(c(NA_real_, NA_real_))
    return(colMeans(x * (nuisance$q_hat - mu)))
  }
  if (method == "AIPW") {
    if (all(is.na(nuisance$q_hat)) || all(is.na(nuisance$pi_hat))) return(c(NA_real_, NA_real_))
    reviewed <- dat$verified == 1L
    residual <- rep(0, nrow(dat))
    residual[reviewed] <- dat$gold_outcome[reviewed] - nuisance$q_hat[reviewed]
    aipw_residual <- nuisance$q_hat - mu + dat$verified / nuisance$pi_hat * residual
    return(colMeans(x * aipw_residual))
  }
  c(NA_real_, NA_real_)
}

estimate_calibration_curve <- function(dat, method, nuisance, grid, curve_df) {
  y <- calibration_curve_outcome(dat, method, nuisance)
  ok <- !is.na(y) & !is.na(dat$pred)
  if (sum(ok) < curve_df + 2L) {
    return(data.frame(method = method, pred = grid, calibrated_risk = NA_real_))
  }
  curve <- fit_spline_curve(y[ok], dat$pred[ok], grid, curve_df)
  data.frame(method = method, pred = grid, calibrated_risk = curve)
}

calibration_curve_outcome <- function(dat, method, nuisance) {
  if (method == "naive") return(dat$proxy_outcome)
  if (method == "gold_only") {
    y <- rep(NA_real_, nrow(dat))
    y[dat$verified == 1L] <- dat$gold_outcome[dat$verified == 1L]
    return(y)
  }
  if (method == "OR") return(nuisance$q_hat)
  if (method == "AIPW") {
    if (all(is.na(nuisance$q_hat)) || all(is.na(nuisance$pi_hat))) return(rep(NA_real_, nrow(dat)))
    reviewed <- dat$verified == 1L
    residual <- rep(0, nrow(dat))
    residual[reviewed] <- dat$gold_outcome[reviewed] - nuisance$q_hat[reviewed]
    return(nuisance$q_hat + dat$verified / nuisance$pi_hat * residual)
  }
  rep(NA_real_, nrow(dat))
}

fit_spline_curve <- function(y, pred, grid, curve_df) {
  dat <- data.frame(y = y, lp = qlogis_clip(pred))
  fit <- stats::lm(y ~ splines::ns(lp, df = curve_df), data = dat)
  out <- stats::predict(fit, newdata = data.frame(lp = qlogis_clip(grid)))
  pmin(pmax(as.numeric(out), 0), 1)
}

n_used_for_method <- function(dat, method) {
  if (method == "gold_only") return(sum(dat$verified == 1L))
  nrow(dat)
}

#' @export
print.goldval_calibration <- function(x, ...) {
  cat("goldval calibration estimates\n\n")
  if (nrow(x$parameters) > 0) print(x$parameters, row.names = FALSE)
  if (nrow(x$curves) > 0) {
    cat("\nCalibration curve grid points:", length(unique(x$curves$pred)), "\n")
  }
  invisible(x)
}

#' @export
plot.goldval_calibration <- function(x, ...) {
  if (nrow(x$curves) == 0) {
    stop("No calibration curve data available.", call. = FALSE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  methods <- unique(x$curves$method)
  cols <- seq_along(methods)
  graphics::plot(
    x$curves$pred,
    x$curves$calibrated_risk,
    type = "n",
    xlab = "Predicted risk",
    ylab = "Estimated observed risk",
    xlim = range(x$curves$pred, na.rm = TRUE),
    ylim = c(0, 1)
  )
  graphics::abline(0, 1, col = "gray60", lty = 2)
  for (i in seq_along(methods)) {
    dat <- x$curves[x$curves$method == methods[[i]], , drop = FALSE]
    graphics::lines(dat$pred, dat$calibrated_risk, col = cols[[i]], lwd = 2)
  }
  graphics::legend("topleft", legend = methods, col = cols, lwd = 2, bty = "n")
  invisible(x)
}
