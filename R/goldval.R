#' Create a goldval validation object
#'
#' @param pred Numeric vector of predicted probabilities, or a column name when
#'   `data` is supplied.
#' @param proxy_outcome Binary proxy outcome observed for every patient, or a
#'   column name when `data` is supplied.
#' @param gold_outcome Binary gold-standard outcome; missing for unverified
#'   patients, or a column name when `data` is supplied.
#' @param verified Logical or binary indicator for gold-standard review, or a
#'   column name when `data` is supplied.
#' @param data Optional data frame. If supplied, the four core inputs must be
#'   character column names or direct vectors. Expressions are not evaluated in
#'   `data`.
#' @param covariates Optional character vector naming covariates in `data`.
#'
#' @return A `goldval_object`.
#'
#' @details Binary outcomes must be coded as `0/1` or `FALSE/TRUE`. Factors such
#' as `"Yes"`/`"No"` are intentionally rejected and should be recoded explicitly.
#' Exact 0 and 1 predictions are retained in the object; clipping is used only
#' by logit-based diagnostic models.
#' @export
goldval <- function(pred,
                    proxy_outcome,
                    gold_outcome,
                    verified,
                    data = NULL,
                    covariates = NULL) {
  if (!is.null(data) && !is.data.frame(data)) {
    stop("`data` must be a data frame or NULL.", call. = FALSE)
  }

  pred <- resolve_input(pred, data, "pred")
  proxy_outcome <- resolve_input(proxy_outcome, data, "proxy_outcome")
  gold_outcome <- resolve_input(gold_outcome, data, "gold_outcome")
  verified <- resolve_input(verified, data, "verified")

  n <- length(pred)
  assert_same_length(n, proxy_outcome, gold_outcome, verified)
  assert_numeric_probability(pred, "pred")
  proxy_outcome <- normalize_binary(proxy_outcome, "proxy_outcome", allow_na = FALSE)
  verified <- normalize_binary(verified, "verified", allow_na = FALSE)
  gold_outcome <- normalize_binary(gold_outcome, "gold_outcome", allow_na = TRUE)

  if (!any(verified == 1L) || !any(verified == 0L)) {
    stop("`verified` must include both reviewed and unreviewed patients.", call. = FALSE)
  }
  if (any(is.na(gold_outcome[verified == 1L]))) {
    stop("`gold_outcome` must be observed for all verified patients.", call. = FALSE)
  }
  if (any(!is.na(gold_outcome[verified == 0L]))) {
    warning("`gold_outcome` is observed for some unverified patients; these values will be ignored.")
    gold_outcome[verified == 0L] <- NA_integer_
  }

  covariate_data <- NULL
  if (!is.null(covariates)) {
    if (is.null(data)) {
      stop("`covariates` requires `data`.", call. = FALSE)
    }
    missing_covariates <- setdiff(covariates, names(data))
    if (length(missing_covariates) > 0) {
      stop("Missing covariates in `data`: ", paste(missing_covariates, collapse = ", "), call. = FALSE)
    }
    covariate_data <- data[covariates]
  }

  object_data <- data.frame(
    pred = as.numeric(pred),
    proxy_outcome = proxy_outcome,
    gold_outcome = gold_outcome,
    verified = verified
  )
  if (!is.null(covariate_data)) {
    object_data <- cbind(object_data, covariate_data)
  }

  structure(
    list(
      data = object_data,
      n = nrow(object_data),
      n_verified = sum(object_data$verified == 1L),
      covariates = covariates,
      settings = list(prediction_clip_eps = 1e-6),
      call = match.call(),
      fitted = list(),
      diagnostics = list(),
      results = list()
    ),
    class = "goldval_object"
  )
}

resolve_input <- function(x, data, arg_name) {
  if (!is.null(data) && is.character(x) && length(x) == 1L) {
    if (!x %in% names(data)) {
      stop("Column `", x, "` supplied to `", arg_name, "` was not found in `data`.", call. = FALSE)
    }
    return(data[[x]])
  }
  x
}

assert_same_length <- function(n, ...) {
  lengths <- vapply(list(...), length, integer(1))
  if (any(lengths != n)) {
    stop("All input vectors must have the same length.", call. = FALSE)
  }
}

assert_numeric_probability <- function(x, name) {
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`", name, "` must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(x) | x < 0 | x > 1)) {
    stop("`", name, "` must contain probabilities in [0, 1].", call. = FALSE)
  }
}

normalize_binary <- function(x, name, allow_na) {
  if (is.logical(x)) {
    out <- as.integer(x)
  } else if (is.factor(x)) {
    stop("`", name, "` must be coded as 0/1 or FALSE/TRUE. Factors such as Yes/No must be recoded explicitly.", call. = FALSE)
  } else if (is.numeric(x) || is.integer(x)) {
    if (any(!is.na(x) & !is.finite(x))) {
      stop("`", name, "` must contain only finite binary values.", call. = FALSE)
    }
    bad_numeric <- !is.na(x) & !(x %in% c(0, 1))
    if (any(bad_numeric)) {
      stop("`", name, "` must be binary 0/1 or TRUE/FALSE.", call. = FALSE)
    }
    out <- as.integer(x)
  } else {
    stop("`", name, "` must be coded as 0/1 or FALSE/TRUE.", call. = FALSE)
  }

  if (!allow_na && anyNA(out)) {
    stop("`", name, "` must not contain missing values.", call. = FALSE)
  }
  bad <- !is.na(out) & !(out %in% c(0L, 1L))
  if (any(bad)) {
    stop("`", name, "` must be binary 0/1 or TRUE/FALSE.", call. = FALSE)
  }
  out
}

#' @export
print.goldval_object <- function(x, ...) {
  dat <- x$data
  cat("goldval clinical prediction validation object\n\n")
  cat("Patients:        ", x$n, "\n", sep = "")
  cat("Gold reviewed:   ", x$n_verified, " (", round(100 * x$n_verified / x$n, 1), "%)\n", sep = "")
  cat("Proxy positive:  ", round(100 * mean(dat$proxy_outcome == 1L), 1), "%\n", sep = "")
  cat("Gold events:     ", sum(dat$gold_outcome[dat$verified == 1L]), "\n", sep = "")
  cat("Prediction range:", paste(round(range(dat$pred), 3), collapse = " - "), "\n")
  cat("Covariates:      ", length(x$covariates), "\n", sep = "")
  cat("Diagnostics fit: ", if (length(x$diagnostics)) "yes" else "no", "\n", sep = "")
  cat("Results fit:     ", if (length(x$results)) "yes" else "no", "\n", sep = "")
  invisible(x)
}
