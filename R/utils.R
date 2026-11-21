# =============================================================================
# Title: Shared internal helpers
# Description: Helpers shared by the imputation estimators (R/imputation.R)
#   and dr_ml_attgt() (R/dr_ml.R): building one (g,t) cell's outcomes,
#   weights and covariate matrices, and fitting weighted regressions among
#   comparison units with collinear columns dropped.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-07-27
# =============================================================================

# --- Cell data -------------------------------------------------------------

# Design matrix of a one-sided formula on `data`, without an intercept
# column (the regressions add their own). Built from the original variables
# in one step, so terms that transform them -- factor(), I(), log() -- work,
# as in ptetools' gt_covariates(). NULL gives a matrix with no columns.
#
# Arguments: formula (one-sided or NULL), data (data.frame).
# Returns: numeric matrix with nrow(data) rows.
bc_design <- function(formula, data) {
  if (is.null(formula)) {
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0))
  }
  # na.pass keeps one row per unit, aligned with the other columns; missing
  # values are caught in bc_cell_data()
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  mm <- stats::model.matrix(attr(mf, "terms"), data = mf)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

# Everything an estimator needs from one (g,t) cell, one row per unit, in
# the order of the pre-period rows (which is the order compute.pte()
# expects for the influence function):
#
#   id, D, DeltaY, w      unit id, treatment indicator, outcome change and
#                         sampling weight (normalized to mean 1 in the cell,
#                         as DRDID and ptetools do; all 1 without weights)
#   has_bc, bc_var        whether there is a bad control, and its name
#   bc_pre, bc_post       the bad control in the base and current periods
#                         (X_{t*-1}, X_{t*} in the paper)
#   X, dX                 general covariates: base-period levels
#                         (xformula) and changes (d_covs_formula)
#   W, dW                 auxiliary covariates for the bad control: levels
#                         (bad_control_cov_formula) and changes
#                         (bad_control_d_cov_formula)
#   g, tp                 the cell's group and time period (internal codes)
#
# With panel data the base period's covariates are used, as in ptetools and
# did. Change covariates whose change is the same for every unit in the cell
# (e.g. time-invariant variables) are dropped with a cell warning: they are
# collinear with the intercept and carry no information.
#
# Arguments: gt_data (from ptetools::two_by_two_subset()) and the
#   covariate formulas of didbc().
# Returns: a list with the elements above.
bc_cell_data <- function(gt_data, xformula = ~1, bad_control_formula = NULL,
                         d_covs_formula = NULL,
                         bad_control_cov_formula = NULL,
                         bad_control_d_cov_formula = NULL) {
  pre <- gt_data[gt_data$name == "pre", , drop = FALSE]
  post <- gt_data[gt_data$name == "post", , drop = FALSE]
  post <- post[match(pre$id, post$id), , drop = FALSE]
  if (anyNA(post$id)) {
    stop("Each unit needs an observation in both periods of a (g,t) cell; ",
         "didbc() requires a balanced panel.", call. = FALSE)
  }

  w <- pre$.w
  if (is.null(w)) w <- rep(1, nrow(pre))
  w <- w / mean(w)

  if (is.null(xformula)) xformula <- ~1
  X <- bc_design(xformula, pre)
  W <- bc_design(bad_control_cov_formula, pre)

  # changes between the current and base periods, dropping those that are
  # the same for every unit
  bc_change <- function(formula) {
    dm <- bc_design(formula, post) - bc_design(formula, pre)
    if (ncol(dm) == 0) return(dm)
    colnames(dm) <- paste0("d_", colnames(dm))
    constant <- apply(dm, 2, function(x) {
      isTRUE(all.equal(min(x), max(x)))
    })
    if (any(constant)) {
      ptetools::pte_cell_warning(
        paste0("Change covariates dropped because their change is the same ",
               "for every unit (e.g. time-invariant variables)"),
        paste(colnames(dm)[constant], collapse = ", ")
      )
    }
    dm[, !constant, drop = FALSE]
  }
  dX <- bc_change(d_covs_formula)
  dW <- bc_change(bad_control_d_cov_formula)

  has_bc <- !is.null(bad_control_formula)
  bc_var <- NULL
  bc_pre <- bc_post <- NULL
  if (has_bc) {
    # validated by didbc(): exactly one variable, untransformed
    bc_var <- all.vars(bad_control_formula)[1]
    bc_pre <- pre[[bc_var]]
    bc_post <- post[[bc_var]]
  }

  covs <- list(X, dX, W, dW, bc_pre, bc_post)
  if (any(vapply(covs, anyNA, logical(1)))) {
    stop("Missing values in the covariates or the bad control; ",
         "drop those units before calling didbc().", call. = FALSE)
  }

  list(
    id = pre$id, D = pre$D, DeltaY = post$Y - pre$Y, w = w,
    has_bc = has_bc, bc_var = bc_var, bc_pre = bc_pre, bc_post = bc_post,
    X = X, dX = dX, W = W, dW = dW,
    g = unique(pre$G[pre$D == 1]), tp = unique(post$period)
  )
}

# --- Weighted regressions among comparison units ---------------------------

# Weighted least squares (or, with binary = TRUE, a weighted logit) of y on
# an intercept and the columns of X, fit on the rows fit_idx. Columns that
# are collinear in the fitting sample are dropped first (R's pivoted QR, as
# lm() does). When the fitting sample is the comparison units, a dropped
# column that still varies among treated units means the prediction for
# those units extrapolates; that is reported once per call through a
# ptetools cell warning (the rule of ptetools::unsupported_covariates(), the
# same check as ptetools::pte_attgt()). Columns that are redundant for
# treated units too -- duplicates, factor levels nobody in the cell has --
# are dropped silently.
#
# Sampling weights enter as the regression weights, so the fit estimates
# the population regression; with all weights equal the fit is the
# unweighted one.
#
# Arguments: X (matrix, all units of the cell), y (response for fit_idx),
#   w (weights, all units), fit_idx (rows to fit on), D (treatment
#   indicator, all units; used for the warning), binary (logit instead of
#   OLS), warn (check for unsupported covariates; only meaningful when
#   fit_idx are the comparison units).
# Returns: list with Z (the kept design columns, intercept first, all
#   units), coef, fitted (fitted values / probabilities, all units), resid
#   (y minus fitted, fit_idx rows) and info (Fisher-information weight of
#   each unit: p(1-p) for the logit, 1 for OLS).
bc_fit <- function(X, y, w, fit_idx, D, binary = FALSE, warn = TRUE) {
  Z <- cbind("(Intercept)" = 1, X)
  Zf <- Z[fit_idx, , drop = FALSE]
  qr0 <- qr(Zf)
  keep <- sort(qr0$pivot[seq_len(qr0$rank)])
  if (warn && length(keep) < ncol(Z)) {
    bc_warn_dropped(Z, D, keep)
  }
  Z <- Z[, keep, drop = FALSE]
  Zf <- Zf[, keep, drop = FALSE]

  if (binary) {
    # quasibinomial: same estimates as binomial, without warnings about
    # non-integer successes when the weights are not integers
    fit <- stats::glm.fit(Zf, y, weights = w[fit_idx],
                          family = stats::quasibinomial())
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    fitted <- stats::plogis(drop(Z %*% coef))
    info <- fitted * (1 - fitted)
  } else {
    fit <- stats::lm.wfit(Zf, y, w[fit_idx])
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    fitted <- drop(Z %*% coef)
    info <- rep(1, nrow(Z))
  }
  names(coef) <- colnames(Z)
  list(Z = Z, coef = coef, fitted = fitted, resid = y - fitted[fit_idx],
       info = info, keep = keep, binary = binary)
}

# Predictions from a bc_fit() fit for new rows of the same covariates
# (without the intercept column).
bc_fit_predict <- function(fit, newX) {
  eta <- drop(cbind(1, newX)[, fit$keep, drop = FALSE] %*% fit$coef)
  if (fit$binary) stats::plogis(eta) else eta
}

# Cell warning for covariates dropped for collinearity among comparison
# units that matter for treated units' predictions, with the (unweighted)
# share of treated units affected, worded as in ptetools::pte_attgt().
#
# Arguments: Z (design matrix with intercept, all units), D, keep (kept
#   columns of Z).
bc_warn_dropped <- function(Z, D, keep) {
  flags <- ptetools::unsupported_covariates(Z, D, keep)
  if (ncol(flags) == 0) return(invisible(NULL))
  shares <- colMeans(flags)
  ptetools::pte_cell_warning(
    paste0("Covariates dropped because they do not vary among comparison ",
           "units (given the other covariates) but do among treated units; ",
           "the estimates extrapolate for those treated units (see ",
           "ptetools::check_overlap())"),
    paste0(colnames(flags), " (", round(100 * shares, 1),
           "% of treated units)", collapse = ", ")
  )
}

# Weighted mean of the rows of a matrix.
bc_wmean_rows <- function(M, w) {
  colSums(M * w) / sum(w)
}
