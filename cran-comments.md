## Resubmission

This is a resubmission of version 1.0.0, rejected by the incoming pre-tests
with the following NOTE on r-devel-linux-x86_64-debian-gcc:

    Running R code in 'testthat.R' had CPU time 4.6 times elapsed time

The `dr_ml` tests call `grf`, which by default uses every available core.
The tests now pass `num_threads = 2` to every `grf`-backed call, keeping
`R CMD check` within CRAN's two-core limit. Verified locally with
`_R_CHECK_LIMIT_CORES_=TRUE`. The package code is unchanged; the version was
bumped to 1.0.1 for the resubmission.

The incoming-feasibility NOTE "Possibly misspelled words in DESCRIPTION:
Caetano, Callaway, Sant'Anna" lists the surnames of the authors of the paper
the package implements; they are spelled correctly.

## Test environments

- Local Ubuntu 24.04.4 LTS, R 4.6.1, `rcmdcheck::rcmdcheck(args = c("--as-cran"))`:
    - 0 errors, 0 warnings, 1 note (see below).
- Local macOS (Apple Silicon), R 4.6.0, `R CMD check --as-cran` with
  `_R_CHECK_LIMIT_CORES_=TRUE`: 0 errors, 0 warnings, 1 note (new submission).
- GitHub Actions:
    - Windows-latest (R release)
    - Windows-latest (R devel)
    - macOS-latest (R release)
    - Ubuntu-latest (R release)
    - Ubuntu-latest (R devel)
    - All checks passed without issues.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE — New submission (expected;
  this is the package's first release).

## Downstream dependencies

`badcontrols` has no reverse dependencies (this is a new package), checked
via `revdeplite::revdeplite()`.

## Additional comments

* This is the first CRAN release of `badcontrols`, which implements
  imputation, doubly robust, and machine-learning estimators for
  difference-in-differences designs where a time-varying covariate is
  affected by the treatment ("bad controls"), based on Caetano, Callaway,
  Payne, and Sant'Anna (2026) <doi:10.48550/arXiv.2608.03881>.
* Two vignettes are included: a conceptual overview of the identification
  strategies, and a coding walkthrough on a real NLSY79 application
  dataset bundled with the package.
