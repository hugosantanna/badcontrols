## Submission summary

This is a minor release (1.0.1 -> 1.1.0). It requires `ptetools`
(>= 1.1.0), now on CRAN. Main changes (see NEWS.md):

* Support for sampling weights (`weightsname`)
* New `clustervars` argument for clustered standard errors.
* Covariate formulas can transform variables (e.g. `~ factor(state)`),
  which previously failed or was silently ignored.
* Warnings from individual group-time cells (e.g. covariates dropped for
  lack of overlap) are collected and reported once per call.

Results without weights are unchanged.

## Test environments

- Local Ubuntu 26.04.1 LTS, R 4.5.2
- GitHub Actions (all Status: OK, 0 errors | 0 warnings | 0 notes):
    - windows-latest (R release)
    - windows-latest (R devel)
    - macos-latest (R release)
    - ubuntu-latest (R release)
    - ubuntu-latest (R devel)

## R CMD check results

`rcmdcheck::rcmdcheck(args = "--as-cran")`:

0 errors | 0 warnings | 0 notes

The incoming-feasibility check has previously flagged "Caetano, Callaway,
Sant'Anna" as possibly misspelled words in DESCRIPTION; these are the
correctly spelled surnames of the authors of the paper the package
implements.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.
