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
  replicates <- run_bootstrap_replicates(
    object,
    analysis = analysis,
    method = method,
    metrics = metrics,
    B = B,
    seed = seed,
    outcome_formula = outcome_formula,
    verification_formula = verification_formula
  )
  intervals <- bootstrap_intervals(original, replicates, B, min_success_rate, conf_level)
  failure_summary <- bootstrap_failure_summary(replicates, B)

  out <- list(
    estimates = original,
    intervals = intervals,
    replicates = replicates,
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
  metric_map <- bootstrap_metric_map(metrics)
  out <- data.frame()
  if (length(metric_map$performance_metrics) > 0L) {
    perf <- performance(
      object,
      method = method,
      metrics = metric_map$performance_metrics,
      outcome_formula = outcome_formula,
      verification_formula = verification_formula
    )
    out <- perf$estimates
  }
  if ("weak" %in% metrics) {
    cal <- calibration(
      object,
      method = method,
      type = "weak",
      outcome_formula = outcome_formula,
      verification_formula = verification_formula
    )
    weak <- weak_parameters_to_long(cal$parameters)
    out <- rbind(out, weak)
  }
  out
}

bootstrap_metric_map <- function(metrics) {
  performance_metrics <- intersect(metrics, c("auc", "brier"))
  if (length(performance_metrics) == 0L) performance_metrics <- character(0)
  list(performance_metrics = performance_metrics)
}

weak_parameters_to_long <- function(parameters) {
  data.frame(
    method = rep(parameters$method, each = 2L),
    metric = rep(c("weak_intercept", "weak_slope"), times = nrow(parameters)),
    estimate = c(rbind(parameters$weak_calibration_intercept, parameters$weak_calibration_slope)),
    n_used = rep(parameters$n_used, each = 2L),
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

  out <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
    boot_dat <- dat[idx, , drop = FALSE]
    out[[b]] <- run_one_bootstrap_replicate(
      boot_dat,
      replicate = b,
      analysis = analysis,
      method = method,
      metrics = metrics,
      outcome_formula = outcome_formula,
      verification_formula = verification_formula
    )
  }
  do.call(rbind, out)
}

run_one_bootstrap_replicate <- function(boot_dat,
                                        replicate,
                                        analysis,
                                        method,
                                        metrics,
                                        outcome_formula,
                                        verification_formula) {
  boot_obj <- try(
    goldval(
      pred = boot_dat$pred,
      proxy_outcome = boot_dat$proxy_outcome,
      gold_outcome = boot_dat$gold_outcome,
      verified = boot_dat$verified
    ),
    silent = TRUE
  )
  if (inherits(boot_obj, "try-error")) {
    return(failed_replicate_rows(replicate, method, metrics, "estimation_failed"))
  }
  estimates <- try(
    estimate_bootstrap_target(boot_obj, analysis, method, metrics, outcome_formula, verification_formula),
    silent = TRUE
  )
  if (inherits(estimates, "try-error")) {
    return(failed_replicate_rows(replicate, method, metrics, "estimation_failed"))
  }

  replicate_status <- infer_replicate_status(boot_obj, outcome_formula, verification_formula)
  estimates$replicate <- replicate
  estimates$replicate_status <- status_for_estimate_rows(estimates$method, replicate_status)
  estimates[, c("replicate", "method", "metric", "estimate", "n_used", "replicate_status")]
}

infer_replicate_status <- function(object, outcome_formula, verification_formula) {
  dat <- object$data
  nuisance <- fit_goldval_nuisance(dat, outcome_formula, verification_formula)
  if (nuisance$outcome_model_status == "failed") return("q_failed")
  if (nuisance$verification_model_status == "failed") return("pi_failed")
  "success"
}

status_for_estimate_rows <- function(method, replicate_status) {
  if (replicate_status == "q_failed") {
    return(ifelse(method %in% c("OR", "AIPW"), "q_failed", "success"))
  }
  if (replicate_status == "pi_failed") {
    return(ifelse(method == "AIPW", "pi_failed", "success"))
  }
  rep(replicate_status, length(method))
}

failed_replicate_rows <- function(replicate, method, metrics, status) {
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
    replicate_status = status,
    stringsAsFactors = FALSE
  )
}

bootstrap_intervals <- function(original, replicates, B, min_success_rate, conf_level) {
  alpha <- 1 - conf_level
  split_keys <- unique(original[c("method", "metric")])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    vals <- replicates$estimate[replicates$method == key$method & replicates$metric == key$metric & replicates$replicate_status == "success"]
    vals <- vals[is.finite(vals)]
    success_n <- length(vals)
    success_rate <- success_n / B
    point <- original$estimate[original$method == key$method & original$metric == key$metric][[1]]
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
  statuses <- c("success", "q_failed", "pi_failed", "estimation_failed")
  counts <- table(factor(unique(replicates[c("replicate", "replicate_status")])$replicate_status, levels = statuses))
  data.frame(
    replicate_status = names(counts),
    n = as.integer(counts),
    fraction = as.integer(counts) / B,
    stringsAsFactors = FALSE
  )
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
