test_that("weighted_auc agrees with hand-calculated perfect and reversed rankings", {
  expect_equal(weighted_auc(c(0.1, 0.2, 0.8, 0.9), c(0, 0, 1, 1)), 1)
  expect_equal(weighted_auc(c(0.9, 0.8, 0.2, 0.1), c(0, 0, 1, 1)), 0)
  expect_equal(weighted_auc(c(0.5, 0.5), c(0, 1)), 0.5)
})

test_that("performance returns stable estimates for naive and gold-only methods", {
  pred <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  proxy <- c(0, 0, 1, 1, 0, 1)
  verified <- c(1, 0, 1, 0, 1, 0)
  gold <- ifelse(verified == 1, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)

  perf <- performance(obj, method = c("naive", "gold_only"), metrics = c("auc", "brier"))
  expect_s3_class(perf, "goldval_performance")
  expect_equal(nrow(perf$estimates), 4)
  expect_equal(perf$estimates$estimate[perf$estimates$method == "naive" & perf$estimates$metric == "auc"], 1)
  expect_equal(perf$estimates$n_used[perf$estimates$method == "gold_only" & perf$estimates$metric == "auc"], 3)
})

test_that("OR and AIPW performance recover oracle direction in a simple fixture", {
  set.seed(10)
  n <- 600
  pred <- stats::plogis(stats::rnorm(n))
  true <- stats::rbinom(n, 1, pred)
  proxy <- true
  proxy[true == 1L] <- stats::rbinom(sum(true == 1L), 1, 0.75)
  proxy[true == 0L] <- stats::rbinom(sum(true == 0L), 1, 0.02)
  review_prob <- stats::plogis(-2.8 + log(4) * proxy + 0.7 * as.numeric(scale(qlogis(pred))))
  verified <- stats::rbinom(n, 1, review_prob)
  gold <- ifelse(verified == 1L, true, NA)
  obj <- goldval(pred, proxy, gold, verified)

  perf <- performance(obj, method = c("naive", "OR", "AIPW"), metrics = c("auc", "brier"))
  aucs <- subset(perf$estimates, metric == "auc")
  briers <- subset(perf$estimates, metric == "brier")

  expect_true(all(is.finite(aucs$estimate)))
  expect_true(all(is.finite(briers$estimate)))
  expect_true(all(aucs$estimate >= 0 & aucs$estimate <= 1))
  expect_true(all(briers$estimate >= 0))
  expect_equal(perf$nuisance$outcome_model_status, "ok")
  expect_equal(perf$nuisance$verification_model_status, "ok")
})

test_that("performance reports failed nuisance model without crashing", {
  pred <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  proxy <- c(0, 0, 1, 1, 0, 1)
  verified <- c(1, 0, 1, 0, 1, 0)
  gold <- ifelse(verified == 1, proxy, NA)
  obj <- goldval(pred, proxy, gold, verified)
  perf <- performance(obj, method = "AIPW", metrics = "auc", verification_formula = verified ~ missing_column)
  expect_equal(perf$nuisance$verification_model_status, "failed")
  expect_true(is.na(perf$estimates$estimate))
})
