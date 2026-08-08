# bootstrap_goldval() API Contract

Status: v0.1 implementation contract.

This file freezes the user-facing bootstrap behavior before implementation.

## Scope

v0.1 supports scalar performance uncertainty only.

Supported:

- `analysis = "performance"`
- methods: `naive`, `gold_only`, `OR`, `AIPW`
- metrics: `auc`, `brier`, `weak`

The `weak` metric returns weak calibration intercept and slope.

Not supported in v0.1:

- flexible calibration curve uncertainty;
- pointwise curve intervals;
- simultaneous calibration bands;
- cluster/site/hospital bootstrap;
- BCa intervals;
- analytic standard errors;
- cross-fitting;
- double-robustness claims.

## Function

```r
bootstrap_goldval(
  object,
  analysis = "performance",
  method = c("naive", "gold_only", "OR", "AIPW"),
  metrics = c("auc", "brier", "weak"),
  B = 1000,
  seed = NULL,
  outcome_formula = gold_outcome ~ proxy_outcome * qlogis_pred,
  verification_formula = verified ~ proxy_outcome + qlogis_pred,
  min_success_rate = 0.90,
  conf_level = 0.95
)
```

## Bootstrap Unit

The bootstrap unit is the patient row:

```text
(pred_i, proxy_outcome_i, verified_i, gold_outcome_i, covariates_i)
```

Each replicate samples `seq_len(n)` with replacement.

No cluster bootstrap is supported in v0.1.

## Replicate Workflow

For each bootstrap replicate:

```text
sample validation rows with replacement
        |
        v
create bootstrap goldval object
        |
        v
refit q model, if OR or AIPW requested
        |
        v
refit pi model, if AIPW requested
        |
        v
recompute scalar performance estimates
```

The implementation must not reuse original-sample `q_hat` or `pi_hat`.

## Return Class

```r
class(x)
#> "goldval_bootstrap"
```

Return structure:

```r
list(
  estimates = data.frame(),
  intervals = data.frame(),
  replicates = data.frame(),
  replicate_diagnostics = data.frame(),
  failure_summary = data.frame(),
  settings = list(),
  call = match.call()
)
```

## Replicate and Estimate Status

Each bootstrap row-sample has one replicate-level diagnostic status stored in
`replicate_diagnostics$overall_status`:

- `success`
- `q_failed`
- `pi_failed`
- `estimation_failed`
- `partial_failure`

`success` means every requested method/metric estimate in that bootstrap
row-sample has `estimate_status = "ok"`. If nuisance fitting succeeds but at
least one requested method/metric estimate fails, `overall_status` is
`partial_failure`. If all requested method/metric estimates fail after nuisance
fitting succeeds, `overall_status` is `estimation_failed`.

The method/metric layer has a separate `replicates$estimate_status`.

Allowed estimate statuses:

- `ok`
- `q_failed`
- `pi_failed`
- `nuisance_failed`
- `optimization_failed`
- `insufficient_events`
- `estimation_failed`

A non-finite estimate must not be marked as `ok`.

`replicate_diagnostics` also records:

- `q_status`
- `pi_status`
- `q_warning_n`
- `pi_warning_n`
- `q_warning_messages`
- `pi_warning_messages`

For OR-only analyses, the verification model is not fitted and `pi_status` is
`not_fit`.

Replicate failure must not abort the entire bootstrap.

No silent fallback to prevalence or intercept-only nuisance models is allowed.

## CI Rule

Intervals are percentile intervals using bootstrap estimates with
`estimate_status = "ok"` and finite estimates.

If the original-sample point estimate is non-finite or has
`estimate_status != "ok"`, do not report a CI for that method/metric. Return
`NA` interval bounds and set:

```text
interval_status = "ORIGINAL_ESTIMATE_FAILED"
```

If success rate for a method/metric is below `min_success_rate`, do not report
a CI for that method/metric. Return `NA` interval bounds and set:

```text
interval_status = "BOOTSTRAP_UNSTABLE"
```

Otherwise:

```text
interval_status = "ok"
```

## Reproducibility

If `seed` is supplied, repeated calls with the same object, arguments, `B`, and
seed must produce identical bootstrap replicate estimates.

Different seeds should generally produce different bootstrap distributions.

The implementation should not rely on ambient `.Random.seed` when `seed` is
supplied.

## Claim Boundary

Allowed:

- percentile bootstrap interval;
- nuisance models refit within bootstrap;
- AIPW-based correction.

Forbidden:

- doubly robust uncertainty estimator;
- valid simultaneous calibration band;
- novel analytic confidence interval.
