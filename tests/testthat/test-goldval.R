test_that("goldval validates input and stores a stable object", {
  pred <- c(0.1, 0.2, 0.8, 0.9)
  proxy <- c(0, 0, 1, 1)
  gold <- c(0, NA, 1, NA)
  verified <- c(1, 0, 1, 0)

  obj <- goldval(pred, proxy, gold, verified)
  expect_s3_class(obj, "goldval_object")
  expect_equal(obj$n, 4)
  expect_equal(obj$n_verified, 2)
  expect_named(obj$data, c("pred", "proxy_outcome", "gold_outcome", "verified"))
})

test_that("goldval supports data plus string column names", {
  dat <- data.frame(
    p = c(0.1, 0.2, 0.8, 0.9),
    ystar = c(0, 0, 1, 1),
    y = c(0, NA, 1, NA),
    r = c(1, 0, 1, 0),
    age = 1:4
  )
  obj <- goldval("p", "ystar", "y", "r", data = dat, covariates = "age")
  expect_s3_class(obj, "goldval_object")
  expect_equal(obj$covariates, "age")
  expect_equal(obj$settings$prediction_clip_eps, 1e-6)
})

test_that("goldval rejects malformed inputs", {
  expect_error(goldval(c(0.1, 1.2), c(0, 1), c(0, NA), c(1, 0)), "probabilities")
  expect_error(goldval(c(0.1, NA), c(0, 1), c(0, NA), c(1, 0)), "missing")
  expect_error(goldval(c(-0.1, 0.2), c(0, 1), c(0, NA), c(1, 0)), "probabilities")
  expect_error(goldval(c(0.1, 0.2), c(0, 2), c(0, NA), c(1, 0)), "binary")
  expect_error(goldval(c(0.1, 0.2), c(0, 1.2), c(0, NA), c(1, 0)), "binary")
  expect_error(goldval(c(0.1, 0.2), c("0", "1"), c(0, NA), c(1, 0)), "coded as 0/1")
  expect_error(goldval(c(0.1, 0.2), factor(c("No", "Yes")), c(0, NA), c(1, 0)), "recoded explicitly")
  expect_error(goldval(c(0.1, 0.2), c(0, 1), c(NA, NA), c(1, 0)), "observed")
  expect_error(goldval(c(0.1, 0.2), c(0, 1), c(0, 1), c(1, 1)), "both reviewed and unreviewed")
  expect_error(goldval(c(0.1, 0.2), c(0, 1, 1), c(0, NA), c(1, 0)), "same length")
  dat <- data.frame(p = c(0.1, 0.2), ystar = c(0, 1), y = c(0, NA), r = c(1, 0))
  expect_error(goldval("missing", "ystar", "y", "r", data = dat), "not found")
  expect_error(goldval("p", "ystar", "y", "r", data = dat, covariates = "age"), "Missing covariates")
})

test_that("gold observed among unverified is warned and ignored", {
  expect_warning(obj <- goldval(c(0.1, 0.2, 0.8), c(0, 1, 1), c(0, 1, 1), c(1, 0, 1)), "unverified")
  expect_true(is.na(obj$data$gold_outcome[2]))
})

test_that("logical outcomes are accepted", {
  obj <- goldval(c(0.1, 0.2, 0.8, 0.9), c(FALSE, FALSE, TRUE, TRUE), c(FALSE, NA, TRUE, NA), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(obj$data$proxy_outcome, c(0L, 0L, 1L, 1L))
  expect_equal(obj$data$verified, c(1L, 0L, 1L, 0L))
})

test_that("diagnose_goldval detects proxy enrichment and positivity", {
  set.seed(1)
  n <- 800
  pred <- stats::plogis(stats::rnorm(n))
  true <- stats::rbinom(n, 1, pred)
  proxy <- true
  proxy[true == 1] <- stats::rbinom(sum(true == 1), 1, 0.75)
  proxy[true == 0] <- stats::rbinom(sum(true == 0), 1, 0.02)
  review_prob <- stats::plogis(-3.4 + log(4) * proxy + 0.8 * as.numeric(scale(qlogis(pred))))
  verified <- stats::rbinom(n, 1, review_prob)
  gold <- ifelse(verified == 1, true, NA)

  obj <- goldval(pred, proxy, gold, verified)
  dx <- diagnose_goldval(obj)

  expect_s3_class(dx, "goldval_diagnostics")
  expect_gt(dx$enrichment$proxy_positive_review_or, 2)
  expect_true(all(c("min_pi", "p05_pi", "median_pi", "max_weight") %in% names(dx$positivity)))
  expect_equal(nrow(dx$risk_deciles), 10)
})

test_that("diagnose_goldval flags weak positivity and sparse tails", {
  set.seed(3)
  n <- 500
  pred <- seq(0.001, 0.999, length.out = n)
  true <- stats::rbinom(n, 1, pred)
  proxy <- true
  verified <- as.integer(pred > 0.9)
  gold <- ifelse(verified == 1, true, NA)
  obj <- goldval(pred, proxy, gold, verified)
  dx <- diagnose_goldval(obj)
  expect_true(any(dx$flags %in% c("POSITIVITY_WEAK", "WEIGHT_UNSTABLE", "TAIL_SPARSE", "PI_MODEL_WARNING")))
})

test_that("diagnose_goldval reports low gold sample and low events", {
  pred <- rep(seq(0.1, 0.9, length.out = 10), 20)
  proxy <- rep(c(0, 1), 100)
  verified <- rep(0L, 200)
  verified[1:20] <- 1L
  gold <- rep(NA_integer_, 200)
  gold[1:20] <- 0L
  gold[20] <- 1L
  obj <- goldval(pred, proxy, gold, verified)
  dx <- diagnose_goldval(obj)
  expect_true("LOW_GOLD_N" %in% dx$flags)
  expect_true("LOW_GOLD_EVENTS" %in% dx$flags)
})

test_that("failed verification model is explicit", {
  pred <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  proxy <- c(0, 0, 1, 1, 0, 1)
  verified <- c(1, 0, 1, 0, 1, 0)
  gold <- ifelse(verified == 1, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)
  dx <- diagnose_goldval(obj, verification_formula = verified ~ missing_column)
  expect_equal(dx$verification_model_status, "failed")
  expect_true("PI_MODEL_FAILED" %in% dx$flags)
  expect_true(all(is.na(dx$positivity$min_pi)))
})

test_that("random verification has approximately flat risk-decile review rates", {
  set.seed(2)
  n <- 1000
  pred <- stats::plogis(stats::rnorm(n))
  true <- stats::rbinom(n, 1, pred)
  proxy <- true
  verified <- stats::rbinom(n, 1, 0.1)
  gold <- ifelse(verified == 1, true, NA)

  obj <- goldval(pred, proxy, gold, verified)
  dx <- diagnose_goldval(obj, verification_formula = verified ~ 1)

  expect_lt(abs(mean(dx$risk_deciles$review_rate) - 0.1), 0.04)
  expect_lt(diff(range(dx$risk_deciles$review_rate)), 0.16)
  expect_lt(dx$positivity$max_weight, 20)
})
