goldval
================

`goldval` is an R workflow for clinical prediction model validation when
fixed predicted risks are available for the full validation cohort, a
proxy binary outcome is observed for everyone, and a reference-standard
outcome is available only for a verified subset. Verification may be
selective rather than random. The package is positioned as biomedical
software/workflow infrastructure and does not claim new AUC, Brier, or
calibration-slope estimators.

## Installation

``` r
remotes::install_github("shiyingqwq/goldval")
```

Requires R \>= 4.1.0. The current release is v0.1.0.

## Core functions

- `goldval()`: constructs the validation object from predicted risk,
  proxy outcome, partial reference outcome, verification indicator, and
  optional covariates.
- `diagnose_goldval()`: verification-process diagnostics, including
  enrichment by proxy status and predicted risk, fitted verification
  probabilities, inverse-verification weights, verified-sample
  representation across predicted-risk deciles, and operational warning
  flags (weak positivity, unstable weights, sparse tails).
- `performance()`: AUC, Brier score, and weak-calibration
  intercept/slope for naive, gold-only, outcome-regression (OR), and
  augmented inverse-probability-weighted (AIPW) analyses.
- `calibration()`: weak-calibration intercept/slope from method-specific
  estimating equations and flexible calibration curves (natural spline
  with four degrees of freedom on logit predicted risk), with plotting
  methods.
- `bootstrap_goldval()`: patient-row percentile bootstrap intervals for
  scalar AIPW analyses (AUC, Brier score, weak-calibration
  intercept/slope), refitting nuisance models in every replicate and
  reporting replicate-level status and failure information.

## Quick start

``` r
library(goldval)

# Direct vectors
x <- goldval(
  pred = validation_data$predicted_risk,
  proxy_outcome = validation_data$ehr_outcome,
  gold_outcome = validation_data$chart_review_outcome,
  verified = validation_data$reviewed
)

# Or column names with a data frame
x <- goldval(
  data = validation_data,
  pred = "predicted_risk",
  proxy_outcome = "ehr_outcome",
  gold_outcome = "chart_review_outcome",
  verified = "reviewed"
)

diagnose_goldval(x)

performance(
  x,
  method = c("naive", "gold_only", "OR", "AIPW"),
  metrics = c("auc", "brier")
)

calibration(
  x,
  method = c("naive", "gold_only", "OR", "AIPW"),
  type = c("weak", "curve")
)
```

## Vignette

The `goldval-workflow` vignette demonstrates the end-to-end workflow on
a fully synthetic example: object construction, verification
diagnostics, naive/gold-only/OR/AIPW performance, weak and flexible
calibration, and scalar bootstrap uncertainty.

## Scope

- Binary reference outcomes and fixed predicted risks.
- Naive, gold-only, outcome-regression, and AIPW analyses.
- AUC, Brier score, weak calibration, and flexible calibration curves.
- Patient-row percentile bootstrap uncertainty for supported scalar AIPW
  analyses.

Explicit non-goals:

- no novel AUC, Brier, weak-calibration, or flexible-calibration
  estimator claims;
- no simultaneous calibration bands or flexible-curve confidence
  intervals;
- no DCA, NRI/IDI, survival, fairness, SHAP, or additional metrics in
  the current version;
- bootstrap percentile intervals are validated only for selected AIPW
  scalar targets and evaluated settings.

## License

MIT + file LICENSE (see [LICENSE](LICENSE)).
