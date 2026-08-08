make_bootstrap_fixture <- function(n = 250, seed = 100) {
  set.seed(seed)
  pred <- stats::plogis(stats::rnorm(n))
  true <- stats::rbinom(n, 1, pred)
  proxy <- true
  proxy[true == 1L] <- stats::rbinom(sum(true == 1L), 1, 0.85)
  proxy[true == 0L] <- stats::rbinom(sum(true == 0L), 1, 0.05)
  review_prob <- stats::plogis(-1.5 + log(3) * proxy + 0.4 * as.numeric(scale(qlogis(pred))))
  verified <- stats::rbinom(n, 1, review_prob)
  gold <- ifelse(verified == 1L, true, NA)
  goldval(pred, proxy, gold, verified)
}

make_bootstrap_covariate_fixture <- function(n = 300, seed = 101) {
  set.seed(seed)
  age <- stats::rnorm(n)
  site <- stats::rbinom(n, 1, 0.45)
  pred <- stats::plogis(-0.1 + 0.8 * age + 0.5 * site + stats::rnorm(n, sd = 0.3))
  true <- stats::rbinom(n, 1, pred)
  proxy <- ifelse(true == 1L, stats::rbinom(n, 1, 0.85), stats::rbinom(n, 1, 0.08))
  verified <- stats::rbinom(n, 1, stats::plogis(-1 + 0.6 * proxy + 0.5 * age + 0.4 * site))
  gold <- ifelse(verified == 1L, true, NA)
  dat <- data.frame(pred, proxy, gold, verified, age, site)
  goldval("pred", "proxy", "gold", "verified", data = dat, covariates = c("age", "site"))
}

test_that("bootstrap_goldval is reproducible with fixed seed", {
  obj <- make_bootstrap_fixture()
  b1 <- bootstrap_goldval(obj, method = c("naive", "OR"), metrics = c("auc", "brier"), B = 20, seed = 123)
  b2 <- bootstrap_goldval(obj, method = c("naive", "OR"), metrics = c("auc", "brier"), B = 20, seed = 123)

  expect_s3_class(b1, "goldval_bootstrap")
  expect_identical(b1$replicates, b2$replicates)
  expect_identical(b1$replicate_diagnostics, b2$replicate_diagnostics)
  expect_identical(b1$intervals, b2$intervals)
})

test_that("bootstrap_goldval changes distribution with different seeds", {
  obj <- make_bootstrap_fixture()
  b1 <- bootstrap_goldval(obj, method = "naive", metrics = "auc", B = 20, seed = 123)
  b2 <- bootstrap_goldval(obj, method = "naive", metrics = "auc", B = 20, seed = 124)

  expect_false(identical(b1$replicates$estimate, b2$replicates$estimate))
})

test_that("bootstrap_goldval supports weak calibration only", {
  obj <- make_bootstrap_fixture(n = 400)
  b <- bootstrap_goldval(obj, method = c("OR", "AIPW"), metrics = "weak", B = 10, seed = 125)

  expect_true(all(b$intervals$metric %in% c("weak_intercept", "weak_slope")))
  expect_equal(nrow(b$replicates), 10 * 2 * 2)
  expect_true("estimate_status" %in% names(b$replicates))
})

test_that("bootstrap failure does not crash and unstable intervals are flagged", {
  obj <- make_bootstrap_fixture()
  b <- bootstrap_goldval(
    obj,
    method = c("OR", "AIPW"),
    metrics = c("auc", "weak"),
    B = 10,
    seed = 126,
    outcome_formula = gold_outcome ~ missing_column,
    min_success_rate = 0.9
  )

  expect_true(any(b$replicate_diagnostics$overall_status == "q_failed"))
  expect_true(any(b$replicates$estimate_status == "q_failed"))
  expect_true(all(is.na(b$intervals$conf_low)))
  expect_true(all(b$intervals$interval_status %in% c("BOOTSTRAP_UNSTABLE", "ORIGINAL_ESTIMATE_FAILED")))
  expect_true(all(b$failure_summary$fraction >= 0 & b$failure_summary$fraction <= 1))
})

test_that("bootstrap validates arguments", {
  obj <- make_bootstrap_fixture()
  expect_error(bootstrap_goldval(obj, B = 1.5), "positive integer")
  expect_error(bootstrap_goldval(obj, min_success_rate = 1), "in \\(0, 1\\)")
  expect_error(bootstrap_goldval(obj, conf_level = 0), "in \\(0, 1\\)")
})

test_that("bootstrap preserves covariates used by nuisance formulas", {
  obj <- make_bootstrap_covariate_fixture()
  b <- bootstrap_goldval(
    obj,
    method = c("OR", "AIPW"),
    metrics = "auc",
    B = 10,
    seed = 127,
    outcome_formula = gold_outcome ~ proxy_outcome + qlogis_pred + age,
    verification_formula = verified ~ proxy_outcome + qlogis_pred + site
  )

  expect_true(all(c("q_status", "pi_status") %in% names(b$replicate_diagnostics)))
  expect_false(all(b$replicates$estimate_status %in% c("q_failed", "pi_failed", "estimation_failed")))
})

test_that("bootstrap records nuisance warning bookkeeping separately", {
  nuisance <- list(
    outcome_model_status = "warning",
    outcome_model_warning_n = 1L,
    outcome_model_warnings = "q warning",
    verification_model_status = "warning",
    verification_model_warning_n = 2L,
    verification_model_warnings = c("pi warning 1", "pi warning 2")
  )
  diag <- replicate_diagnostic_row(3, "success", nuisance)

  expect_equal(diag$overall_status, "success")
  expect_equal(diag$q_status, "warning")
  expect_equal(diag$pi_status, "warning")
  expect_equal(diag$q_warning_n, 1L)
  expect_equal(diag$pi_warning_n, 2L)
  expect_match(diag$pi_warning_messages, "pi warning 1")
})

test_that("original estimate failure blocks interval generation", {
  n <- 120
  pred <- stats::plogis(seq(-2, 2, length.out = n))
  proxy <- rep(c(0L, 1L), length.out = n)
  verified <- rep(c(1L, 0L), length.out = n)
  gold <- ifelse(verified == 1L, 1L, NA)
  obj <- goldval(pred, proxy, gold, verified)

  b <- bootstrap_goldval(obj, method = "gold_only", metrics = "auc", B = 5, seed = 129)

  expect_true(is.na(b$estimates$estimate))
  expect_equal(b$estimates$estimate_status, "estimation_failed")
  expect_equal(b$intervals$interval_status, "ORIGINAL_ESTIMATE_FAILED")
})

test_that("bootstrap restores ambient RNG state when seed is supplied", {
  obj <- make_bootstrap_fixture()
  set.seed(999)
  before <- .Random.seed
  invisible(bootstrap_goldval(obj, method = "naive", metrics = "auc", B = 5, seed = 130))

  expect_identical(.Random.seed, before)
})
