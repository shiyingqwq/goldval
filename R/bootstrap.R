#' Bootstrap uncertainty for goldval analyses
#'
#' @param object A `goldval_object`.
#' @param analysis Analysis type. v0.1 supports `"performance"` only.
#' @param method Character vector of validation methods.
#' @param metrics Character vector of scalar metrics. Use `"weak"` for weak
#'   calibration intercept and slope.
#' @param B Number of bootstrap replicates.
#' @param seed Optional random seed.
#' @param outcome_formula Formula for the gold-outcome regression model.
#' @param verification_formula Formula for the verification model.
#' @param min_success_rate Minimum replicate success rate required to report a
#'   percentile interval.
#' @param conf_level Confidence level for percentile intervals.
#'
#' @return A `goldval_bootstrap` object.
#' @export
bootstrap_goldval <- function(object,
                              analysis = "performance",
                              method = c("naive", "gold_only", "OR", "AIPW"),
                              metrics = c("auc", "brier", "weak"),
                              B = 1000,
                              seed = NULL,
                              outcome_formula = gold_outcome ~ proxy_outcome * qlogis_pred,
                              verification_formula = verified ~ proxy_outcome + qlogis_pred,
                              min_success_rate = 0.90,
                              conf_level = 0.95) {
  if (!inherits(object, "goldval_object")) {
    stop("`object` must be a goldval_object.", call. = FALSE)
  }
  analysis <- match.arg(analysis, "performance")
  method <- match.arg(method, c("naive", "gold_only", "OR", "AIPW"), several.ok = TRUE)
  metrics <- match.arg(metrics, c("auc", "brier", "weak"), several.ok = TRUE)
  B <- validate_bootstrap_B(B)
  min_success_rate <- validate_probability_scalar(min_success_rate, "min_success_rate")
  conf_level <- validate_probability_scalar(conf_level, "conf_level")

  original <- estimate_bootstrap_target(object, analysis, method, metrics, outcome_formula, verification_formula)
  boot <- run_bootstrap_replicates(
    object,
    analysis = analysis,
    method = method,
    metrics = metrics,
    B = B,
    seed = seed,
    outcome_formula = outcome_formula,
    verification_formula = verification_formula
  )
  intervals <- bootstrap_intervals(original, boot$replicates, B, min_success_rate, conf_level)
  failure_summary <- bootstrap_failure_summary(boot$replicates, B)

  out <- list(
    estimates = original,
    intervals = intervals,
    replicates = boot$replicates,
    replicate_diagnostics = boot$replicate_diagnostics,
    failure_summary = failure_summary,
    settings = list(
      analysis = analysis,
      method = method,
      metrics = metrics,
      B = B,
      seed = seed,
      min_success_rate = min_success_rate,
      conf_level = conf_level,
      bootstrap_unit = "patient_row"
    ),
    call = match.call()
  )
  class(out) <- "goldval_bootstrap"
  out
}

estimate_bootstrap_target <- function(object, analysis, method, metrics, outcome_formula, verification_formula) {
  if (analysis != "performance") stop("Unsupported analysis.", call. = FALSE)
  dat <- object$data
  estimate_bootstrap_target_data(dat, analysis, method, metrics, outcome_formula, verification_formula)
}

estimate_bootstrap_target_data <- function(dat, analysis, method, metrics, outcome_formula, verification_formula) {
  if (analysis != "performance") stop("Unsupported analysis.", call. = FALSE)
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  nuisance <- if (any(method %in% c("OR", "AIPW"))) {
    fit_goldval_nuisance(dat, outcome_formula, verification_formula)
  } else {
    empty_performance_nuisance(dat)
  }
  estimate_bootstrap_target_with_nuisance(dat, method, metrics, nuisance)
}

estimate_bootstrap_target_with_nuisance <- function(dat, method, metrics, nuisance) {
  out <- data.frame()
  performance_metrics <- intersect(metrics, c("auc", "brier"))
  if (length(performance_metrics) > 0L) {
    perf <- do.call(rbind, lapply(method, function(m) {
      estimate_performance_bootstrap_rows(dat, m, performance_metrics, nuisance)
    }))
    out <- rbind(out, perf)
  }
  if ("weak" %in% metrics) {
    weak <- do.call(rbind, lapply(method, function(m) {
      weak_parameters_to_long(estimate_weak_calibration(dat, m, nuisance))
    }))
    out <- rbind(out, weak)
  }
  row.names(out) <- NULL
  out
}

estimate_performance_bootstrap_rows <- function(dat, method, metrics, nuisance) {
  y <- performance_outcome(dat, method, nuisance)
  metric_values <- compute_performance_metrics(dat, y, metrics, method, nuisance)
  estimates <- as.numeric(metric_values)
  data.frame(
    method = method,
    metric = names(metric_values),
    estimate = estimates,
    n_used = attr(metric_values, "n_used"),
    estimate_status = estimate_status_for_performance(method, estimates, nuisance),
    stringsAsFactors = FALSE
  )
}

estimate_status_for_performance <- function(method, estimates, nuisance) {
  status <- rep("ok", length(estimates))
  if (method %in% c("OR", "AIPW") && nuisance$outcome_model_status == "failed") {
    status[] <- "q_failed"
  }
  if (method == "AIPW" && nuisance$verification_model_status == "failed") {
    status[] <- "pi_failed"
  }
  status[!is.finite(estimates) & status == "ok"] <- "estimation_failed"
  status
}

weak_parameters_to_long <- function(parameters) {
  estimate <- c(rbind(parameters$weak_calibration_intercept, parameters$weak_calibration_slope))
  status <- rep(parameters$weak_status, each = 2L)
  status[!is.finite(estimate) & status == "ok"] <- "estimation_failed"
  data.frame(
    method = rep(parameters$method, each = 2L),
    metric = rep(c("weak_intercept", "weak_slope"), times = nrow(parameters)),
    estimate = estimate,
    n_used = rep(parameters$n_used, each = 2L),
    estimate_status = status,
    stringsAsFactors = FALSE
  )
}

run_bootstrap_replicates <- function(object,
                                     analysis,
                                     method,
                                     metrics,
                                     B,
                                     seed,
                                     outcome_formula,
                                     verification_formula) {
  dat <- object$data
  rng_state <- NULL
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rng_state <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
      if (is.null(rng_state)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", rng_state, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  estimates <- vector("list", B)
  diagnostics <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
    boot_dat <- dat[idx, , drop = FALSE]
    one <- run_one_bootstrap_replicate(
      boot_dat,
      covariates = object$covariates,
      replicate = b,
      analysis = analysis,
      method = method,
      metrics = metrics,
      outcome_formula = outcome_formula,
      verification_formula = verification_formula
    )
    estimates[[b]] <- one$estimates
    diagnostics[[b]] <- one$diagnostics
  }
  list(
    replicates = do.call(rbind, estimates),
    replicate_diagnostics = do.call(rbind, diagnostics)
  )
}

run_one_bootstrap_replicate <- function(boot_dat,
                                        covariates,
                                        replicate,
                                        analysis,
                                        method,
                                        metrics,
                                        outcome_formula,
                                        verification_formula) {
  boot_obj <- try(
    goldval(
      pred = "pred",
      proxy_outcome = "proxy_outcome",
      gold_outcome = "gold_outcome",
      verified = "verified",
      data = boot_dat,
      covariates = covariates
    ),
    silent = TRUE
  )
  if (inherits(boot_obj, "try-error")) {
    return(list(
      estimates = failed_estimate_rows(replicate, method, metrics, "estimation_failed"),
      diagnostics = replicate_diagnostic_row(replicate, "estimation_failed", NULL)
    ))
  }

  dat <- boot_obj$data
  dat$qlogis_pred <- qlogis_clip(dat$pred)
  nuisance <- try({
    if (any(method %in% c("OR", "AIPW"))) {
      fit_goldval_nuisance(dat, outcome_formula, verification_formula)
    } else {
      empty_performance_nuisance(dat)
    }
  }, silent = TRUE)

  if (inherits(nuisance, "try-error")) {
    return(list(
      estimates = failed_estimate_rows(replicate, method, metrics, "estimation_failed"),
      diagnostics = replicate_diagnostic_row(replicate, "estimation_failed", NULL)
    ))
  }

  estimates <- try(
    estimate_bootstrap_target_with_nuisance(dat, method, metrics, nuisance),
    silent = TRUE
  )
  if (inherits(estimates, "try-error")) {
    return(list(
      estimates = failed_estimate_rows(replicate, method, metrics, "estimation_failed"),
      diagnostics = replicate_diagnostic_row(replicate, "estimation_failed", nuisance)
    ))
  }

  estimates$replicate <- replicate
  estimates <- estimates[, c("replicate", "method", "metric", "estimate", "n_used", "estimate_status")]
  list(
    estimates = estimates,
    diagnostics = replicate_diagnostic_row(replicate, overall_replicate_status(nuisance), nuisance)
  )
}

overall_replicate_status <- function(nuisance) {
  if (nuisance$outcome_model_status == "failed") return("q_failed")
  if (nuisance$verification_model_status == "failed") return("pi_failed")
  "success"
}

replicate_diagnostic_row <- function(replicate, overall_status, nuisance) {
  if (is.null(nuisance)) {
    return(data.frame(
      replicate = replicate,
      overall_status = overall_status,
      q_status = "not_fit",
      pi_status = "not_fit",
      q_warning_n = 0L,
      pi_warning_n = 0L,
      q_warning_messages = "",
      pi_warning_messages = "",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    replicate = replicate,
    overall_status = overall_status,
    q_status = nuisance$outcome_model_status,
    pi_status = nuisance$verification_model_status,
    q_warning_n = nuisance$outcome_model_warning_n,
    pi_warning_n = nuisance$verification_model_warning_n,
    q_warning_messages = paste(nuisance$outcome_model_warnings, collapse = " | "),
    pi_warning_messages = paste(nuisance$verification_model_warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

failed_estimate_rows <- function(replicate, method, metrics, status) {
  metric_names <- c(intersect(metrics, c("auc", "brier")), if ("weak" %in% metrics) c("weak_intercept", "weak_slope"))
  grid <- expand.grid(
    method = method,
    metric = metric_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data.frame(
    replicate = replicate,
    method = grid$method,
    metric = grid$metric,
    estimate = NA_real_,
    n_used = NA_integer_,
    estimate_status = status,
    stringsAsFactors = FALSE
  )
}

bootstrap_intervals <- function(original, replicates, B, min_success_rate, conf_level) {
  alpha <- 1 - conf_level
  split_keys <- unique(original[c("method", "metric")])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    original_row <- original[original$method == key$method & original$metric == key$metric, , drop = FALSE][1L, ]
    point <- original_row$estimate
    if (!identical(original_row$estimate_status, "ok") || !is.finite(point)) {
      return(data.frame(
        method = key$method,
        metric = key$metric,
        estimate = point,
        conf_low = NA_real_,
        conf_high = NA_real_,
        success_n = 0L,
        success_rate = 0,
        interval_status = "ORIGINAL_ESTIMATE_FAILED",
        stringsAsFactors = FALSE
      ))
    }

    vals <- replicates$estimate[
      replicates$method == key$method &
        replicates$metric == key$metric &
        replicates$estimate_status == "ok"
    ]
    vals <- vals[is.finite(vals)]
    success_n <- length(vals)
    success_rate <- success_n / B
    if (success_rate < min_success_rate || success_n < 2L) {
      return(data.frame(
        method = key$method,
        metric = key$metric,
        estimate = point,
        conf_low = NA_real_,
        conf_high = NA_real_,
        success_n = success_n,
        success_rate = success_rate,
        interval_status = "BOOTSTRAP_UNSTABLE",
        stringsAsFactors = FALSE
      ))
    }
    qs <- stats::quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, na.rm = TRUE)
    data.frame(
      method = key$method,
      metric = key$metric,
      estimate = point,
      conf_low = qs[[1]],
      conf_high = qs[[2]],
      success_n = success_n,
      success_rate = success_rate,
      interval_status = "ok",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

bootstrap_failure_summary <- function(replicates, B) {
  statuses <- sort(unique(replicates$estimate_status))
  split_keys <- unique(replicates[c("method", "metric")])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    x <- replicates[replicates$method == key$method & replicates$metric == key$metric, , drop = FALSE]
    counts <- table(factor(x$estimate_status, levels = statuses))
    data.frame(
      method = key$method,
      metric = key$metric,
      estimate_status = names(counts),
      n = as.integer(counts),
      fraction = as.integer(counts) / B,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

validate_bootstrap_B <- function(B) {
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || !is.finite(B) || B != as.integer(B) || B < 1L) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }
  as.integer(B)
}

validate_probability_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", name, "` must be a single number in (0, 1).", call. = FALSE)
  }
  x
}

#' @export
print.goldval_bootstrap <- function(x, ...) {
  cat("goldval bootstrap intervals\n\n")
  print(x$intervals, row.names = FALSE)
  invisible(x)
}
