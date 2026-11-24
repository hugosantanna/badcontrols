## Submission summary

This is a minor release (1.0.1 -> 1.1.0). It requires `ptetools`
(>= 1.1.0), now on CRAN. Main changes (see NEWS.md):

* Sampling weights (`weightsname`) now work: previously they only entered
  the aggregation across groups, so group-time estimates and their
  standard errors were unweighted.
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

(Locally, there are two NOTEs that come from the check machine, not the
package: HTML validation and math rendering in the manual were skipped
because HTML Tidy and V8 are not installed, and the current time could not
be verified.)

The incoming-feasibility check has previously flagged "Caetano, Callaway,
Sant'Anna" as possibly misspelled words in DESCRIPTION; these are the
correctly spelled surnames of the authors of the paper the package
implements.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.
