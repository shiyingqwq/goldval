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

test_that("bootstrap_goldval is reproducible with fixed seed", {
  obj <- make_bootstrap_fixture()
  b1 <- bootstrap_goldval(obj, method = c("naive", "OR"), metrics = c("auc", "brier"), B = 20, seed = 123)
  b2 <- bootstrap_goldval(obj, method = c("naive", "OR"), metrics = c("auc", "brier"), B = 20, seed = 123)

  expect_s3_class(b1, "goldval_bootstrap")
  expect_identical(b1$replicates, b2$replicates)
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

  expect_true(any(b$replicates$replicate_status == "q_failed"))
  expect_true(all(is.na(b$intervals$conf_low)))
  expect_true(all(b$intervals$interval_status == "BOOTSTRAP_UNSTABLE"))
})

test_that("bootstrap validates arguments", {
  obj <- make_bootstrap_fixture()
  expect_error(bootstrap_goldval(obj, B = 1.5), "positive integer")
  expect_error(bootstrap_goldval(obj, min_success_rate = 1), "in \\(0, 1\\)")
  expect_error(bootstrap_goldval(obj, conf_level = 0), "in \\(0, 1\\)")
})
