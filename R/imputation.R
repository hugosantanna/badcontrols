# =============================================================================
# Title: Imputation estimators for ATT(g,t) with a bad control
# Description: imputation_attgt() and the two identification strategies it
#   dispatches to (imputation_unconfoundedness(), imputation_did()): impute
#   the bad control's untreated potential value for treated units, then an
#   outcome regression among comparison units, with the influence function
#   that accounts for both estimation steps. Sampling weights enter every
#   regression, mean and influence function.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-07-25
# =============================================================================

#' Imputation ATT(g,t) estimator for bad controls
#'
#' @param gt_data data.frame from ptetools::two_by_two_subset
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1}).
#'   Terms can transform variables, e.g. \code{~ factor(state) + I(z^2)}.
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment).
#'   \code{NULL} (the default) means no bad control at all.
#' @param d_covs_formula One-sided formula for general exogenous covariates,
#'   entered as their change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param bad_control_cov_formula One-sided formula for auxiliary covariates
#'   (W in the paper) used to model the bad control's counterfactual
#'   evolution, entered as their pre-treatment-period level. Any variable
#'   name can be used here, including the outcome itself (to reproduce the
#'   old "lagged outcome as W" behavior). Default \code{NULL} (unused).
#' @param bad_control_d_cov_formula Like \code{bad_control_cov_formula}, but
#'   entered as a change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param bad_control_binary Logical; whether the bad control is binary
#'   (detected automatically by \code{didbc()}). If \code{TRUE}, Step 1 (the
#'   bad-control evolution model) is fit by logistic regression instead of
#'   OLS, and the influence function is adjusted accordingly. The continuous
#'   case (\code{FALSE}, the default) is unaffected. Ignored when
#'   \code{bad_control_identification_strategy = "did"} (see below).
#' @param bad_control_identification_strategy Which assumption identifies the
#'   bad control's untreated potential evolution: \code{"unconfoundedness"}
#'   (default), i.e. Covariate Unconfoundedness for the bad control given
#'   \code{(bc_pre, W, Z)}, or \code{"did"}, i.e. parallel trends for the bad
#'   control itself given \code{(W, Z)}
#'   (`app:bad-control-parallel-trends` in the supplementary appendix).
#'   \code{"did"} requires \code{bad_control_formula} to be non-\code{NULL}
#'   and does not distinguish a binary bad control.
#' @param ... unused
#' @details
#' This function only builds the cell's data (one row per unit, with its
#' sampling weight and covariate matrices; see \code{bc_cell_data()}), then
#' dispatches to \code{\link{imputation_unconfoundedness}} or
#' \code{\link{imputation_did}} depending on
#' \code{bad_control_identification_strategy}.
#'
#' With sampling weights (\code{weightsname} in \code{didbc()}), both steps
#' are weighted regressions, the ATT is the weighted mean over treated units,
#' and the influence function carries each unit's weight, so standard errors
#' are for the weighted estimator.
#' @keywords internal
imputation_attgt <- function(gt_data,
                             xformula = ~1,
                             bad_control_formula = NULL,
                             d_covs_formula = NULL,
                             bad_control_cov_formula = NULL,
                             bad_control_d_cov_formula = NULL,
                             bad_control_binary = FALSE,
                             bad_control_identification_strategy = c("unconfoundedness", "did"),
                             ...) {
  bad_control_identification_strategy <- match.arg(bad_control_identification_strategy)

  cd <- bc_cell_data(
    gt_data,
    xformula = xformula,
    bad_control_formula = bad_control_formula,
    d_covs_formula = d_covs_formula,
    bad_control_cov_formula = bad_control_cov_formula,
    bad_control_d_cov_formula = bad_control_d_cov_formula
  )

  if (bad_control_identification_strategy == "did") {
    imputation_did(cd)
  } else {
    imputation_unconfoundedness(cd, bad_control_binary = bad_control_binary)
  }
}

#' Imputation estimator under Covariate Unconfoundedness for the bad control
#'
#' Implements the two-step imputation estimator and its influence function
#' (`eq:psi-ra-final` in the supplementary appendix) under Covariate
#' Unconfoundedness: the bad control's untreated potential evolution is
#' modeled from \code{(bc_pre, W, Z)} among the comparison group. Called by
#' \code{\link{imputation_attgt}} when
#' \code{bad_control_identification_strategy = "unconfoundedness"}.
#'
#' @param cd one (g,t) cell's data, from \code{bc_cell_data()}
#' @param bad_control_binary Logical; whether the bad control is binary (see
#'   \code{\link{imputation_attgt}})
#' @keywords internal
imputation_unconfoundedness <- function(cd, bad_control_binary = FALSE) {
  comp <- which(cd$D == 0)

  # Step 1: impute the untreated potential bad control X_{t*}(0) for
  # treated units, from a regression of X_{t*} on (X_{t*-1}, W, Z) among
  # comparison units. Binary bad control: logit (Lambda(S'gamma) in
  # eq:ra-x-evol's generalization); continuous: OLS. With no bad control
  # there is nothing to impute and no Step 1.
  step1 <- NULL
  if (cd$has_bc) {
    S <- cbind(cd$bc_pre, cd$W, cd$dW, cd$dX, cd$X)
    colnames(S)[1] <- paste0(cd$bc_var, "_pre")
    step1 <- bc_fit(S, cd$bc_post[comp], cd$w, comp, cd$D,
                    binary = bad_control_binary)
    bc_imp <- ifelse(cd$D == 1, step1$fitted, cd$bc_post)
  }

  imputation_step2(cd, bc_imp = if (cd$has_bc) bc_imp, step1 = step1)
}

#' Imputation estimator under parallel trends for the bad control
#'
#' Implements the two-step imputation estimator and its influence function
#' under parallel trends for the bad control itself, rather than Covariate
#' Unconfoundedness: \code{ass:bad-control-parallel-trends} and
#' \code{cor:att-under-bad-control-parallel-trends-and-linearity} in
#' \code{app:bad-control-parallel-trends} of the supplementary appendix, and
#' \code{dev/bad_control_parallel_trends_influence_function.md} for the
#' influence function derivation. Step 1 regresses the bad control's own
#' change, rather than its post-period level, on \code{(W, Z)} among the
#' comparison group; Step 2 (the outcome regression) is unchanged from
#' \code{\link{imputation_unconfoundedness}}. Called by
#' \code{\link{imputation_attgt}} when
#' \code{bad_control_identification_strategy = "did"}.
#'
#' @inheritParams imputation_unconfoundedness
#' @keywords internal
imputation_did <- function(cd) {
  if (!cd$has_bc) {
    stop(
      "bad_control_identification_strategy = \"did\" requires a bad ",
      "control (bad_control_formula cannot be NULL); parallel trends for ",
      "the bad control is an identification strategy for the bad ",
      "control's own evolution, so it has no meaningful degenerate case ",
      "without one."
    )
  }
  comp <- which(cd$D == 0)

  # Step 1': E[DeltaX_t* | W, Z] = W'delta_1 + Z'delta_2
  # (ass:deltaX-linear-cef), identified among the comparison group under
  # ass:bad-control-parallel-trends. bc_pre is deliberately NOT a regressor
  # here (unlike imputation_unconfoundedness's Step 1), matching the paper's
  # conditioning set. The imputed X_t*(0) is bc_pre plus the predicted
  # change.
  WZ <- cbind(cd$W, cd$dW, cd$dX, cd$X)
  delta_bc <- cd$bc_post - cd$bc_pre
  step1 <- bc_fit(WZ, delta_bc[comp], cd$w, comp, cd$D)
  bc_imp <- ifelse(cd$D == 1, cd$bc_pre + step1$fitted, cd$bc_post)

  # Step 2 and the influence function are those of
  # imputation_unconfoundedness(), with (S, gamma, v, Sigma_S) ->
  # ((W,Z), (delta_1, delta_2), eps, Sigma_WZ) -- see
  # dev/bad_control_parallel_trends_influence_function.md's "Sample analog"
  # section. The derivative of the imputed X_t*(0) with respect to delta is
  # (W,Z) itself, so the Fisher-information weight is 1 (OLS).
  imputation_step2(cd, bc_imp = bc_imp, step1 = step1)
}

# Step 2 of both imputation estimators, the ATT and its influence function.
#
# Step 2 regresses DeltaY on (imputed X_{t*}(0), X_{t*-1}, Z) among
# comparison units (for them the "imputed" value is the observed X_{t*});
# mu-hat = nu_0-hat is its prediction for every unit, and
#   ATT-hat = sum_{D=1} w (DeltaY - mu-hat) / sum_{D=1} w     (eq:att-ra).
#
# Influence function (eq:psi-ra-final in the supplementary appendix, with
# sampling weights w normalized to mean 1 and pi-hat = mean(w D)):
#   treated:     w (DeltaY - mu-hat - ATT) / pi-hat
#   comparison: -w (R' kappa_r u + S' kappa_s v) / (1 - pi-hat)
# where u, v are the Step 2 and Step 1 residuals, R and S their regressors,
#   kappa_r = Sigma_R^{-1} E_w[R~ | D = 1],
#   kappa_s = Sigma_S^{-1} beta_1 E_w[info * S | D = 1],
# Sigma_R, Sigma_S the weighted second-moment matrices among comparison
# units (Sigma_S with the logit's Fisher-information weight p(1-p) for a
# binary bad control, 1 for OLS), R~ the Step 2 regressors at the imputed
# value, and beta_1 Step 2's coefficient on the imputed bad control. The
# comparison-group terms account for estimating Step 2's beta and Step 1's
# gamma (the "OLS asymptotic linear representations" lemma); with weights
# the estimating equations are weighted, so each unit's contribution is
# multiplied by its weight.
#
# Arguments: cd (from bc_cell_data()), bc_imp (imputed bad control, NULL
#   without one), step1 (the Step 1 fit from bc_fit(), NULL without a bad
#   control).
# Returns: a ptetools::attgt_if object.
imputation_step2 <- function(cd, bc_imp, step1) {
  D <- cd$D
  w <- cd$w
  comp <- which(D == 0)
  trt <- which(D == 1)
  n <- length(D)

  R <- cbind(cd$dX, cd$X)
  bc_imp_name <- NULL
  if (cd$has_bc) {
    bc_imp_name <- paste0(cd$bc_var, "_post")
    R <- cbind(bc_imp, cd$bc_pre, R)
    colnames(R)[1:2] <- c(bc_imp_name, paste0(cd$bc_var, "_pre"))
  }
  step2 <- bc_fit(R, cd$DeltaY[comp], w, comp, D)
  mu_hat <- step2$fitted

  att <- sum(w[trt] * (cd$DeltaY[trt] - mu_hat[trt])) / sum(w[trt])

  pi_hat <- mean(w * D)
  w_c <- w[comp]

  # Step 2's contribution. sigma_r is normalized by the comparison units'
  # total weight (n0 without weights), hence the 1 / (1 - pi-hat) below.
  r_mat <- step2$Z[comp, , drop = FALSE]
  sigma_r <- crossprod(r_mat, r_mat * w_c) / sum(w_c)
  kappa_r <- solve(sigma_r, bc_wmean_rows(step2$Z[trt, , drop = FALSE], w[trt]))
  correction <- drop(r_mat %*% kappa_r) * step2$resid

  # Step 1's contribution, through Step 2's coefficient on the imputed bad
  # control (0 if that column was dropped for collinearity)
  if (!is.null(step1)) {
    beta1_hat <- step2$coef[bc_imp_name]
    if (is.na(beta1_hat)) beta1_hat <- 0
    s_mat <- step1$Z[comp, , drop = FALSE]
    sigma_s <- crossprod(s_mat, s_mat * step1$info[comp] * w_c) / sum(w_c)
    # LOGIT-SPECIFIC: the information weight at treated units' own S_i enters
    # through the chain rule when linearizing Lambda(S'gamma) -- see the
    # theory note for eq:nu-linearize's generalization. 1 for OLS.
    s_bar <- bc_wmean_rows(step1$Z[trt, , drop = FALSE] * step1$info[trt],
                           w[trt])
    kappa_s <- solve(sigma_s, beta1_hat * s_bar)
    correction <- correction + drop(s_mat %*% kappa_s) * step1$resid
  }

  inf_func <- rep(0, n)
  inf_func[trt] <- w[trt] * (cd$DeltaY[trt] - mu_hat[trt] - att) / pi_hat
  inf_func[comp] <- -w_c * correction / (1 - pi_hat)

  ptetools::attgt_if(attgt = att, inf_func = inf_func)
}
