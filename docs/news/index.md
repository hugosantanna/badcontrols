# Changelog

## badcontrols 1.1.0

- Sampling weights (`weightsname`) now work. They were passed on to
  ptetools, which only used them to aggregate across groups, so
  `ATT(g,t)` and its standard errors were those of the unweighted
  estimator. They now enter both imputation estimators (each regression,
  the average over treated units and the influence function, including
  the logit for a binary bad control) and `dr_ml` (every nuisance
  function, with `grf`’s `sample.weights` for forests, and the doubly
  robust score). Integer weights give the same estimates as duplicating
  units; without a bad control, the weighted imputation estimator
  matches `did::att_gt(est_method = "reg")` with weights, standard
  errors included. Results without weights are unchanged.
- New `clustervars` argument to cluster standard errors at a level above
  the unit, as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html):
  one variable, constant within each unit. It applies to analytical and
  bootstrap standard errors (through ptetools). With
  `est_method = "dr_ml"`, cross-fitting folds are also assigned by
  cluster, so no cluster has units in both the training and held-out
  folds.
- Covariate formulas can now transform variables, e.g.
  `xformula = ~ factor(state) + I(z^2)` or
  `bad_control_cov_formula = ~ log(w)`. Transformed terms in `xformula`
  failed with “object not found”, and in the other covariate formulas
  the transformation was silently ignored. `est_method = "dr_ml"` with
  `nuisance_method = "ml"` now also accepts factor and character
  covariates (they used to fail inside `grf`). A transformed
  `bad_control_formula` (e.g. `~ log(x)`) now gives an error instead of
  silently using the untransformed variable.
- When a covariate does not vary among comparison units (given the
  others) but does among treated units, such as a dummy for a region
  with treated units only, it is now dropped with a warning that the
  estimates extrapolate for those treated units, as in ptetools. Such
  covariates used to be dropped silently by
  [`lm()`](https://rdrr.io/r/stats/lm.html), or make the influence
  function fail. Duplicated covariates and empty factor levels are
  dropped silently.
- Warnings are now given once per call, listing the group-time cells,
  instead of once per cell (and again for every empirical bootstrap
  draw): `dr_ml`’s fallbacks to the imputation estimator (too few
  treated units, or propensity scores above `overlap_threshold`, now
  reported as the number of units rather than their ids), and change
  covariates that are dropped because their change is the same for every
  unit.
- `boot_type = "empirical"` no longer warns “aggte_fun not specified”
  (from ptetools).
- Requires ptetools (\>= 1.1.0).

## badcontrols 1.0.1

CRAN release: 2026-09-12

- Tests now pin `num_threads = 2` for the `grf`-backed estimators so
  that `R CMD check` stays within CRAN’s two-core limit. No user-facing
  changes.

## badcontrols 1.0.0

- First CRAN release.
- Added `num_threads` argument to
  [`didbc()`](https://hsantanna.org/badcontrols/reference/didbc.md)/[`dr_ml_attgt()`](https://hsantanna.org/badcontrols/reference/dr_ml_attgt.md)
  to control `grf`’s thread usage under `nuisance_method = "ml"`.
- Added conceptual and coding vignettes, the latter using the bundled
  `nlsy_job_displacement` application data.

## badcontrols 0.1.1

- Switched the panel treatment effects backend from `pte` to `ptetools`
  (`pte::pte2()` calls now use
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html)).

## badcontrols 0.1.0

- Initial release.
