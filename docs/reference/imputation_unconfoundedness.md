# Imputation estimator under Covariate Unconfoundedness for the bad control

Implements the two-step imputation estimator and its influence function
(`eq:psi-ra-final` in the supplementary appendix) under Covariate
Unconfoundedness: the bad control's untreated potential evolution is
modeled from `(bc_pre, W, Z)` among the comparison group. Called by
[`imputation_attgt`](https://hsantanna.org/badcontrols/reference/imputation_attgt.md)
when `bad_control_identification_strategy = "unconfoundedness"`.

## Usage

``` r
imputation_unconfoundedness(cd, bad_control_binary = FALSE)
```

## Arguments

- cd:

  one (g,t) cell's data, from `bc_cell_data()`

- bad_control_binary:

  Logical; whether the bad control is binary (see
  [`imputation_attgt`](https://hsantanna.org/badcontrols/reference/imputation_attgt.md))
