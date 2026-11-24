# Imputation estimator under parallel trends for the bad control

Implements the two-step imputation estimator and its influence function
under parallel trends for the bad control itself, rather than Covariate
Unconfoundedness: `ass:bad-control-parallel-trends` and
`cor:att-under-bad-control-parallel-trends-and-linearity` in
`app:bad-control-parallel-trends` of the supplementary appendix, and
`dev/bad_control_parallel_trends_influence_function.md` for the
influence function derivation. Step 1 regresses the bad control's own
change, rather than its post-period level, on `(W, Z)` among the
comparison group; Step 2 (the outcome regression) is unchanged from
[`imputation_unconfoundedness`](https://hsantanna.org/badcontrols/reference/imputation_unconfoundedness.md).
Called by
[`imputation_attgt`](https://hsantanna.org/badcontrols/reference/imputation_attgt.md)
when `bad_control_identification_strategy = "did"`.

## Usage

``` r
imputation_did(cd)
```

## Arguments

- cd:

  one (g,t) cell's data, from `bc_cell_data()`
