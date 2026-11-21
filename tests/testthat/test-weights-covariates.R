# =============================================================================
# Title: Tests for sampling weights, covariate formulas and cell warnings
# Description: Sampling weights in didbc() (equivalence with duplicated
#   units, invariance to rescaling, agreement with did), transformed and
#   factor covariates, the warnings collected once per call, and
#   cross-fitting folds by cluster.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-11-20
# =============================================================================

# --- Shared data -----------------------------------------------------------

# Small simulated panel with integer sampling weights (constant within
# unit), and the same data with each unit duplicated `wt` times under new
# ids. Weighting by wt should give the same estimates as the duplicated
# data, and -- clustering the duplicated data on the original unit -- the
# same standard errors.
make_weight_data <- function() {
  set.seed(3)
  sim <- simulate_bad_controls(
    n = 400, T_max = 3, groups = 3, dgp = "dgp1", beta_drift = 0
  )
  d <- sim$data
  d$X2 <- as.integer(d$X > stats::median(d$X))
  d$wt <- stats::ave(sample(1:3, nrow(d), TRUE), d$id,
                     FUN = function(x) x[1])
  dup <- d[rep(seq_len(nrow(d)), d$wt), ]
  dup$copy <- stats::ave(dup$period, dup$id, dup$period, FUN = seq_along)
  dup$orig_id <- dup$id
  dup$id <- dup$id * 10 + dup$copy
  dup <- dup[order(dup$id, dup$period), ]
  list(d = d, dup = dup)
}

run_didbc <- function(data, ...) {
  suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = data, bstrap = FALSE, ...
  ))
}

# --- Sampling weights ------------------------------------------------------

test_that("weights equal duplicating units, for every estimator", {
  wd <- make_weight_data()
  configs <- list(
    unconfoundedness = list(bad_control_formula = ~X, xformula = ~Z,
                            bad_control_cov_formula = ~W),
    did = list(bad_control_formula = ~X, xformula = ~Z,
               bad_control_cov_formula = ~W,
               bad_control_identification_strategy = "did"),
    binary = list(bad_control_formula = ~X2, xformula = ~Z,
                  bad_control_cov_formula = ~W),
    dr_ml = list(bad_control_formula = ~X, xformula = ~Z,
                 bad_control_cov_formula = ~W, est_method = "dr_ml",
                 nuisance_method = "parametric")
  )
  for (cf in configs) {
    set.seed(9)
    weighted <- do.call(run_didbc, c(list(wd$d, weightsname = "wt"), cf))
    set.seed(9)
    duplicated <- do.call(run_didbc,
                          c(list(wd$dup, clustervars = "orig_id"), cf))
    expect_equal(weighted$att_gt$att, duplicated$att_gt$att, tolerance = 1e-8)
    expect_equal(weighted$att_gt$se, duplicated$att_gt$se, tolerance = 1e-8)
    expect_equal(weighted$overall_att$overall.att,
                 duplicated$overall_att$overall.att, tolerance = 1e-8)
  }
})

test_that("weights matter and rescaling them changes nothing", {
  wd <- make_weight_data()
  args <- list(bad_control_formula = ~X, xformula = ~Z,
               bad_control_cov_formula = ~W)
  weighted <- do.call(run_didbc, c(list(wd$d, weightsname = "wt"), args))
  unweighted <- do.call(run_didbc, c(list(wd$d), args))
  d2 <- wd$d
  d2$wt <- d2$wt * 7.3
  rescaled <- do.call(run_didbc, c(list(d2, weightsname = "wt"), args))
  expect_false(isTRUE(all.equal(weighted$att_gt$att, unweighted$att_gt$att)))
  expect_equal(rescaled$att_gt$att, weighted$att_gt$att)
  expect_equal(rescaled$att_gt$se, weighted$att_gt$se)
})

test_that("without a bad control, weighted imputation matches did (reg)", {
  skip_if_not_installed("did")
  wd <- make_weight_data()
  res <- run_didbc(wd$d, xformula = ~Z, weightsname = "wt")
  ref <- suppressWarnings(did::att_gt(
    yname = "Y", tname = "period", idname = "id", gname = "G",
    xformla = ~Z, data = wd$d, est_method = "reg", weightsname = "wt",
    control_group = "notyettreated", bstrap = FALSE, cband = FALSE
  ))
  expect_equal(res$att_gt$att, ref$att, tolerance = 1e-8)
  expect_equal(res$att_gt$se, ref$se, tolerance = 1e-8)
  agg <- did::aggte(ref, type = "simple", bstrap = FALSE)
  expect_equal(res$overall_att$overall.att, agg$overall.att, tolerance = 1e-8)
  expect_equal(res$overall_att$overall.se, agg$overall.se, tolerance = 1e-8)
})

test_that("dr_ml with forests runs with weights, clusters and factors", {
  skip_if_not_installed("grf")
  wd <- make_weight_data()
  d <- wd$d
  d$race <- c("a", "b", "c")[(d$id %% 3) + 1]
  d$cl <- d$id %% 30
  res <- run_didbc(
    d, bad_control_formula = ~X, xformula = ~Z + race,
    bad_control_cov_formula = ~W, est_method = "dr_ml",
    weightsname = "wt", clustervars = "cl", num_threads = 2
  )
  expect_true(is.finite(res$overall_att$overall.att))
  expect_true(is.finite(res$overall_att$overall.se))
})

# --- Covariate formulas ----------------------------------------------------

test_that("covariate formulas can transform variables", {
  wd <- make_weight_data()
  d <- wd$d
  d$Z2 <- d$Z^2
  d$logW <- log(abs(d$W) + 1)
  transformed <- run_didbc(
    d, bad_control_formula = ~X, xformula = ~ Z + I(Z^2),
    bad_control_cov_formula = ~ log(abs(W) + 1)
  )
  precomputed <- run_didbc(
    d, bad_control_formula = ~X, xformula = ~ Z + Z2,
    bad_control_cov_formula = ~logW
  )
  expect_equal(transformed$att_gt$att, precomputed$att_gt$att)
  expect_equal(transformed$att_gt$se, precomputed$att_gt$se)
})

test_that("a transformed bad control gives a clear error", {
  wd <- make_weight_data()
  expect_error(
    run_didbc(wd$d, bad_control_formula = ~ log(abs(X) + 1), xformula = ~Z),
    "not a transformation"
  )
})

# --- Warnings collected once per call ----------------------------------------

# didbc() with warnings captured rather than suppressed
didbc_warnings <- function(data, ...) {
  msgs <- character(0)
  withCallingHandlers(
    didbc(yname = "Y", gname = "G", tname = "period", idname = "id",
          data = data, bstrap = FALSE, ...),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  msgs
}

test_that("a covariate only treated units have is dropped with one warning", {
  wd <- make_weight_data()
  d <- wd$d
  # a dummy that only (some) treated units have
  d$region <- as.integer(d$G > 0 & d$id %% 4 == 0)
  msgs <- didbc_warnings(d, bad_control_formula = ~X,
                         xformula = ~ Z + region,
                         bad_control_cov_formula = ~W)
  dropped <- grep("Covariates dropped", msgs, value = TRUE)
  expect_length(dropped, 1)
  expect_match(dropped, "region")
})

test_that("constant change covariates and dr_ml fallbacks warn once", {
  wd <- make_weight_data()
  msgs <- didbc_warnings(wd$d, bad_control_formula = ~X, xformula = ~Z,
                         d_covs_formula = ~Z)
  expect_length(grep("change is the same for every unit", msgs), 1)

  msgs <- didbc_warnings(wd$d, bad_control_formula = ~X, xformula = ~Z,
                         est_method = "dr_ml",
                         nuisance_method = "parametric",
                         overlap_threshold = 0.3)
  expect_length(grep("overlap_threshold", msgs), 1)
})

# --- Cross-fitting folds -----------------------------------------------------

test_that("cross-fitting folds keep clusters together", {
  set.seed(1)
  D <- rep(c(1, 0), times = c(40, 160))
  cl <- rep(1:50, each = 4)
  folds <- bc_folds(D, n_folds = 5, cl = cl)
  expect_true(all(tapply(folds, cl, function(f) length(unique(f))) == 1))
  expect_setequal(unique(folds), 1:5)
  # treated units (clusters 1-10) are spread over every fold
  expect_setequal(unique(folds[D == 1]), 1:5)
})
