# Uncertainty Plan for goldval v0.1

Status: design protocol. This document must be approved before bootstrap code is
implemented.

## Goal

The goal of uncertainty v0.1 is not to introduce new inference theory. The goal
is to give software users a transparent, reproducible bootstrap procedure that
propagates the main sources of uncertainty in validation with:

- full-cohort predicted risks;
- full-cohort proxy outcomes;
- partial gold-standard review;
- outcome-regression and AIPW correction models.

The package should describe this as practical bootstrap uncertainty, not as a
novel analytic standard-error estimator.

## Estimands

v0.1 scalar estimands:

- AUC;
- Brier score;
- weak calibration intercept;
- weak calibration slope.

v0.1 flexible calibration estimand:

- pointwise flexible calibration curve values on a fixed prediction grid.

However, implementation order is staged:

1. scalar performance uncertainty first;
2. flexible calibration pointwise uncertainty later;
3. no simultaneous calibration bands in v0.1.

## Bootstrap Unit

The bootstrap unit is the validation cohort row:

```text
(pred_i, proxy_outcome_i, verified_i, gold_outcome_i, covariates_i)
```

Each bootstrap replicate samples rows with replacement from the full validation
cohort. This preserves the observed joint data structure and resamples both
reviewed and unreviewed patients.

## Required Bootstrap Workflow

Each bootstrap replicate must refit nuisance models:

```text
bootstrap full validation rows
        |
        v
recreate goldval object
        |
        v
refit q model: P(Y = 1 | proxy outcome, prediction, covariates)
        |
        v
refit pi model: P(verified = 1 | proxy outcome, prediction, covariates)
        |
        v
recalculate naive / gold-only / OR / AIPW estimators
```

The bootstrap must not fix `q_hat` or `pi_hat` from the original sample and only
resample final outputs. That would understate uncertainty by ignoring nuisance
model fitting.

## Supported v0.1 Inference

Initial implementation:

```r
bootstrap_goldval(
  object,
  statistic = c("performance"),
  method = c("naive", "gold_only", "OR", "AIPW"),
  metrics = c("auc", "brier", "weak"),
  B = 500,
  seed = NULL
)
```

Supported output:

- point estimate from original data;
- bootstrap replicate estimates;
- percentile 95% interval;
- bootstrap success/failure counts;
- nuisance warning/failure counts.

The first implementation should support `performance()` only. Flexible
calibration curve intervals are intentionally deferred.

## Unsupported in v0.1

Do not implement in the first uncertainty layer:

- BCa intervals;
- analytic influence-function standard errors;
- simultaneous calibration bands;
- cross-fitting;
- optimal review design;
- double-robustness claims;
- selective-review design optimization.

These are possible later extensions, but they would move the project back toward
statistical methods research rather than software MVP delivery.

## Reproducibility Rules

Bootstrap functions must support a `seed` argument.

With fixed data, fixed `B`, fixed methods/metrics, and fixed `seed`, repeated
calls should return identical bootstrap replicate estimates.

Required test:

```r
set.seed-independent call 1:
  bootstrap_goldval(object, B = 50, seed = 123)

set.seed-independent call 2:
  bootstrap_goldval(object, B = 50, seed = 123)

expect_identical(result1$replicates, result2$replicates)
```

The function should avoid relying on the ambient `.Random.seed` when `seed` is
provided.

## Failure Handling

A single failed bootstrap replicate must not abort the entire bootstrap run.

For each replicate:

- catch fitting errors;
- record failed nuisance models;
- return `NA` estimates for failed method/metric combinations;
- continue to the next replicate.

The returned object should include:

```r
list(
  estimates = data.frame(),
  intervals = data.frame(),
  replicates = data.frame(),
  diagnostics = list(
    B_requested,
    B_completed,
    B_failed,
    nuisance_warning_n,
    nuisance_failure_n
  )
)
```

## Reporting Language

Allowed wording:

- percentile bootstrap interval;
- bootstrap uncertainty propagation;
- nuisance models refitted within bootstrap;
- AIPW-based correction.

Avoid wording:

- novel confidence interval;
- doubly robust inference;
- valid simultaneous calibration band;
- analytic standard error.

## First Implementation Gate

The first bootstrap implementation is accepted only if:

- scalar `performance()` bootstrap works for naive, gold-only, OR, and AIPW;
- `B = 50` reproducibility test passes exactly with fixed seed;
- failed nuisance models are recorded rather than silently ignored;
- `testthat` passes;
- development `R CMD check --no-manual --no-build-vignettes` remains OK.

Only after this gate passes should pointwise flexible calibration curve
uncertainty be implemented.
