test_that("calibration estimates perfect calibration approximately", {
  set.seed(20)
  n <- 3000
  pred <- stats::plogis(stats::rnorm(n))
  y <- stats::rbinom(n, 1, pred)
  verified <- stats::rbinom(n, 1, 0.3)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, y, gold, verified)

  cal <- calibration(obj, method = "naive", type = "weak")
  expect_s3_class(cal, "goldval_calibration")
  expect_lt(abs(cal$parameters$weak_calibration_intercept), 0.15)
  expect_lt(abs(cal$parameters$weak_calibration_slope - 1), 0.2)
})

test_that("calibration detects systematic intercept and slope", {
  set.seed(21)
  n <- 5000
  pred <- stats::plogis(stats::rnorm(n))
  lp <- qlogis(pred)
  y_prob <- stats::plogis(-0.5 + 0.7 * lp)
  y <- stats::rbinom(n, 1, y_prob)
  verified <- stats::rbinom(n, 1, 0.3)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, y, gold, verified)

  cal <- calibration(obj, method = "naive", type = "weak")
  expect_lt(abs(cal$parameters$weak_calibration_intercept + 0.5), 0.12)
  expect_lt(abs(cal$parameters$weak_calibration_slope - 0.7), 0.12)
})

test_that("naive calibration curve agrees with direct spline fit", {
  set.seed(22)
  n <- 400
  pred <- stats::plogis(stats::rnorm(n))
  y <- stats::rbinom(n, 1, pred)
  verified <- stats::rbinom(n, 1, 0.3)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, y, gold, verified)
  grid <- seq(stats::quantile(pred, 0.02), stats::quantile(pred, 0.98), length.out = 25)

  cal <- calibration(obj, method = "naive", type = "curve", grid = grid, curve_df = 4)
  direct <- fit_spline_curve(y, pred, grid, 4)
  expect_equal(cal$curves$calibrated_risk, direct, tolerance = 1e-10)
})

test_that("gold-only curve ignores unverified gold outcomes", {
  set.seed(23)
  n <- 200
  pred <- stats::plogis(stats::rnorm(n))
  y <- stats::rbinom(n, 1, pred)
  verified <- rep(0L, n)
  verified[seq(1, n, by = 2)] <- 1L
  gold1 <- ifelse(verified == 1L, y, NA)
  gold2 <- gold1
  gold2[verified == 0L] <- 1L - y[verified == 0L]
  obj1 <- goldval(pred, y, gold1, verified)
  obj2 <- suppressWarnings(goldval(pred, y, gold2, verified))
  grid <- seq(0.1, 0.9, length.out = 20)

  c1 <- calibration(obj1, method = "gold_only", type = "curve", grid = grid)
  c2 <- calibration(obj2, method = "gold_only", type = "curve", grid = grid)
  expect_equal(c1$curves$calibrated_risk, c2$curves$calibrated_risk)
})

test_that("OR and AIPW weak calibration match naive with perfect proxy and random verification", {
  set.seed(24)
  n <- 3000
  pred <- stats::plogis(stats::rnorm(n))
  y <- stats::rbinom(n, 1, pred)
  verified <- stats::rbinom(n, 1, 0.5)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, y, gold, verified)

  cal <- calibration(obj, method = c("naive", "OR", "AIPW"), type = "weak")
  params <- cal$parameters
  naive <- params[params$method == "naive", ]
  or <- params[params$method == "OR", ]
  aipw <- params[params$method == "AIPW", ]
  expect_lt(abs(or$weak_calibration_intercept - naive$weak_calibration_intercept), 0.15)
  expect_lt(abs(or$weak_calibration_slope - naive$weak_calibration_slope), 0.15)
  expect_lt(abs(aipw$weak_calibration_intercept - naive$weak_calibration_intercept), 0.15)
  expect_lt(abs(aipw$weak_calibration_slope - naive$weak_calibration_slope), 0.15)
})

test_that("OR and AIPW recover known weak calibration under proxy error", {
  set.seed(26)
  n <- 8000
  pred <- stats::plogis(stats::rnorm(n))
  lp <- qlogis(pred)
  true_prob <- stats::plogis(-0.5 + 0.7 * lp)
  y <- stats::rbinom(n, 1, true_prob)
  proxy <- y
  proxy[y == 1L] <- stats::rbinom(sum(y == 1L), 1, 0.80)
  proxy[y == 0L] <- stats::rbinom(sum(y == 0L), 1, 0.05)
  verified <- stats::rbinom(n, 1, 0.5)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, proxy, gold, verified)

  cal <- calibration(obj, method = c("naive", "OR", "AIPW"), type = "weak")
  params <- cal$parameters
  naive <- params[params$method == "naive", ]
  or <- params[params$method == "OR", ]
  aipw <- params[params$method == "AIPW", ]
  expect_gt(abs(naive$weak_calibration_intercept + 0.5), 0.15)
  expect_lt(abs(or$weak_calibration_intercept + 0.5), 0.15)
  expect_lt(abs(or$weak_calibration_slope - 0.7), 0.15)
  expect_lt(abs(aipw$weak_calibration_intercept + 0.5), 0.15)
  expect_lt(abs(aipw$weak_calibration_slope - 0.7), 0.15)
})

test_that("failed nuisance model yields NA corrected calibration", {
  pred <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  proxy <- c(0, 0, 1, 1, 0, 1)
  verified <- c(1, 0, 1, 0, 1, 0)
  gold <- ifelse(verified == 1L, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)
  cal <- calibration(obj, method = c("OR", "AIPW"), outcome_formula = gold_outcome ~ missing_column)
  expect_equal(cal$nuisance$outcome_model_status, "failed")
  expect_true(all(is.na(cal$parameters$weak_calibration_intercept)))
  expect_true(all(is.na(cal$curves$calibrated_risk)))
})

test_that("calibration validates curve_df and grid", {
  pred <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  proxy <- c(0, 0, 1, 1, 0, 1)
  verified <- c(1, 0, 1, 0, 1, 0)
  gold <- ifelse(verified == 1L, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)
  expect_error(calibration(obj, curve_df = 3.7), "positive integer")
  expect_error(calibration(obj, grid = c(0.1, NA)), "finite probabilities")
  expect_error(calibration(obj, grid = c(-0.1, 0.5)), "finite probabilities")
  cal <- calibration(obj, method = "naive", type = "curve", grid = c(0.8, 0.2, 0.2))
  expect_equal(cal$settings$grid, c(0.2, 0.8))
})

test_that("calibration curve fails clearly with insufficient unique predictions", {
  pred <- rep(0.2, 20)
  proxy <- rep(c(0, 1), 10)
  verified <- rep(c(1, 0), 10)
  gold <- ifelse(verified == 1L, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)
  cal <- calibration(obj, method = "naive", type = "curve", curve_df = 4, grid = c(0.1, 0.2, 0.3))
  expect_true(all(is.na(cal$curves$calibrated_risk)))
  expect_true(all(cal$curves$curve_status == "insufficient_unique_predictions"))
})

test_that("plot.goldval_calibration runs", {
  set.seed(25)
  n <- 200
  pred <- stats::plogis(stats::rnorm(n))
  y <- stats::rbinom(n, 1, pred)
  verified <- stats::rbinom(n, 1, 0.3)
  gold <- ifelse(verified == 1L, y, NA)
  obj <- goldval(pred, y, gold, verified)
  cal <- calibration(obj, method = "naive", type = "curve")
  grDevices::png(tempfile(fileext = ".png"))
  expect_silent(plot(cal))
  grDevices::dev.off()
})
