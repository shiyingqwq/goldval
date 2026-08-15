# goldval 0.1.0

* First release corresponding to the version evaluated in the manuscript.
* Added `performance()` for naive, gold-only, OR, and AIPW analyses of AUC,
  Brier score, and weak-calibration scalar summaries.
* Added `calibration()` for weak-calibration estimating equations and flexible
  calibration curves (natural spline on logit predicted risk, four degrees of
  freedom).
* Added `bootstrap_goldval()` for patient-row percentile bootstrap intervals
  with nuisance-model refitting and replicate-level status reporting.
* Added verification diagnostics with operational warning flags.
* Updated documentation, vignettes, and package metadata.

# goldval 0.0.0.9000

* Initialized package skeleton.
* Added `goldval()` object constructor.
* Added `diagnose_goldval()` verification and positivity diagnostics.
* Added initial `performance()` infrastructure for naive, gold-only, OR, and AIPW estimates.
* Added initial `calibration()` infrastructure for weak calibration and flexible curves.
* Added initial `bootstrap_goldval()` protocol implementation for scalar performance uncertainty.
* Added initial testthat coverage for input validation and diagnostics.
