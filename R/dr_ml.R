# =============================================================================
# Title: Doubly robust (DR/ML) estimator for ATT(g,t) with a bad control
# Description: dr_ml_attgt() -- the cross-fitted doubly robust estimator of
#   Algorithm 1 in the paper, with nuisance functions from random forests
#   (grf) or linear/logit working models -- and its nuisance-fitting and
#   cross-fitting helpers.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-07-25
# =============================================================================

#' @title Doubly Robust Estimator for Bad Controls
#'
#' @description Implements the semiparametric doubly robust estimator from
#'  Caetano, Callaway, Payne, and Sant'Anna, with nuisance functions
#'  estimated either via machine learning or via parametric (linear/logit)
#'  working models, and K-fold cross-fitting. Estimates ATT(g,t) in the
#'  presence of post-treatment covariates (bad controls).
#'
#' @details Implements Algorithm 1 from the paper: four nuisance functions,
#'  with K-fold cross-fitting:
#'  \itemize{
#'    \item m_0(X_t*, X_t*-1, Z): outcome regression
#'      `E[DeltaY | X_t*, X_t*-1, Z, D=0]`
#'    \item nu_0(X_t*-1, W, Z): `E[m_0(X_t*, X_t*-1, Z) | X_t*-1, W, Z, D=0]`,
#'      a second-stage regression of m_0's fitted values on (X_t*-1, W, Z)
#'      among untreated units
#'    \item p_2(X_t*-1, W, Z): the propensity score P(D=1 | X_t*-1, W, Z)
#'    \item omega_0(X_t*, X_t*-1, Z): `E[p_2/(1-p_2) | X_t*, X_t*-1, Z, D=0]`,
#'      a regression of the fitted propensity odds ratio on (X_t*, X_t*-1, Z)
#'      among untreated units
#'  }
#'
#'  \code{nuisance_method = "ml"} (the default) estimates all four via
#'  cross-fitted random forests (\code{grf} package). \code{nuisance_method
#'  = "parametric"} estimates m_0/nu_0/omega_0 via OLS and p_2 via logistic
#'  regression; this assumes all four working models are correctly
#'  specified, in which case the doubly robust score's Neyman orthogonality
#'  means no further correction to the influence function is needed beyond
#'  what is already here for either choice of \code{nuisance_method} --
#'  parametric estimators converge faster than the ML rate conditions this
#'  orthogonality argument requires. Cross-fitting is retained under
#'  \code{"parametric"} for a single, shared code path, even though it is
#'  not required for validity there.
#'
#'  Unlike the imputation estimator, no counterfactual imputation of the bad
#'  control is needed: X_t* is observed directly for untreated units (X_t* =
#'  X_t*(0)), and nu_0/omega_0 target m_0's/p_2's fitted values through a
#'  second regression rather than plugging in a predicted X_t*(0).
#'
#'  nu_0 and omega_0 are themselves nested conditional expectations of
#'  m_0/p_2 (e.g. nu_0 = E\[m_0 | X_t*-1, W, Z, D=0\]), so regressing m_0's/
#'  p_2's fitted values on the coarser feature set using the same
#'  training-fold observations that fit m_0/p_2 would be a generated-
#'  regressor problem: the pseudo-outcome's estimation error is correlated
#'  with the very sample the nested regression is fit on, not just small in
#'  L^2. Under \code{nuisance_method = "ml"}, this is avoided using each
#'  forest's out-of-bag (OOB) predictions -- \code{grf}'s \code{predict(fit)}
#'  with no \code{newdata}, which averages only over trees that did not
#'  include a given row in their subsample -- as the pseudo-outcome for
#'  nu_0/omega_0's training data, instead of the in-sample
#'  \code{predict(fit, newdata = <training data>)}. This costs nothing extra
#'  to fit (no additional forests, unlike an explicit training-fold split)
#'  and keeps the full training fold available to nu_0/omega_0, at the cost
#'  of a less clean-cut independence argument than literal disjoint
#'  subsamples would give. The main m_0/p_2 fits used directly in the doubly
#'  robust score are unaffected -- they are already evaluated out-of-sample
#'  via the outer K-fold. Under \code{nuisance_method = "parametric"}, no
#'  OOB equivalent applies (\code{lm}/\code{glm} have no such notion); that
#'  path is unchanged, per the paragraph above.
#'
#'  The doubly robust score is consistent if either (m_0, nu_0) or (p_2,
#'  omega_0) are correctly specified.
#'
#'  With no bad control at all (\code{bad_control_formula = NULL}), nu_0 and
#'  omega_0 are not estimated at all: nu_0 = m_0 and omega_0 = p_2/(1-p_2)
#'  exactly (m_0 no longer depends on X_t*, so there is nothing left for
#'  nu_0 to marginalize over; p_2 is already a function of Z alone, so
#'  omega_0's further conditioning on Z is a no-op). The doubly robust score
#'  then reduces exactly to the classical Sant'Anna and Zhao (2020) AIPW-DiD
#'  estimator.
#'
#'  There are two triggers that fall back to \code{imputation_attgt()} for a
#'  whole (g,t) cell, rather than dropping p_2/omega_0 alone -- doing that
#'  would leave a moment that is not Neyman orthogonal, which is exactly why
#'  \code{imputation_attgt()} needs its own first-stage correction terms
#'  that this function's nuisances don't have:
#'  \enumerate{
#'    \item \strong{Small treated group.} If the cell has fewer treated
#'      units than the number of p_2 covariates plus \code{min_group_size}.
#'      \code{didbc()} warns once, listing the cells and treated counts. Fitting
#'      p_2 needs a nontrivial treated sample to identify at all, relative
#'      to how many covariates it has to fit (a logit can hit perfect
#'      separation, a \code{probability_forest} is essentially memorizing a
#'      handful of points), which can happen even when the population-level
#'      covariate distributions overlap fine -- a distinct problem from
#'      overlap below.
#'    \item \strong{Overlap.} Before cross-fitting, a preliminary propensity
#'      model (not used in the final estimation) is fit once on the whole
#'      cell: if any unit's fitted propensity exceeds
#'      \code{overlap_threshold}. \code{didbc()} warns once, listing the
#'      cells and the number of units above the threshold.
#'  }
#'
#'  \strong{Sampling weights.} With \code{weightsname} in \code{didbc()},
#'  every nuisance function is fit with the sampling weights (as
#'  \code{sample.weights} in \code{grf}, or as regression weights), and the
#'  doubly robust score, the share of treated units and the influence
#'  function are weighted.
#'
#'  \strong{Clustering.} With \code{clustervars}, the cross-fitting folds
#'  are assigned by cluster, so that no cluster has units in both the
#'  training and held-out folds, as is usual for cross-fitting with
#'  dependent data (Chernozhukov et al. 2018). Cluster-robust standard
#'  errors themselves come from \code{ptetools}.
#'
#'  \strong{Covariates.} Covariate formulas can transform variables (e.g.
#'  \code{~ factor(state) + I(z^2)}); factors are expanded into dummies
#'  for \code{grf}. With \code{nuisance_method = "parametric"}, covariates
#'  that are collinear among the units a model is fit on are dropped, with a
#'  warning when the outcome model would then extrapolate for treated
#'  units.
#'
#' @param gt_data data.frame from \code{ptetools::two_by_two_subset} with
#'   columns id, D, period, name (pre/post), Y, plus covariate columns
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment).
#'   \code{NULL} means no bad control at all; see Details. Multiple bad
#'   controls are not yet supported.
#' @param d_covs_formula One-sided formula for general exogenous covariates,
#'   entered as their change (post minus pre) rather than a level. Used
#'   only in the outcome regression. Default \code{NULL} (unused).
#' @param bad_control_cov_formula One-sided formula for auxiliary covariates
#'   (W in the paper) used to model the bad control's counterfactual
#'   evolution and the propensity score, entered as their pre-treatment
#'   level. Any variable name can be used, including the outcome itself
#'   (e.g. \code{~Y} reproduces the old "lagged outcome as W" behavior).
#'   Default \code{NULL} (unused, ignored if \code{bad_control_formula} is
#'   \code{NULL}).
#' @param bad_control_d_cov_formula Like \code{bad_control_cov_formula}, but
#'   entered as a change; used only in the outcome regression, not the
#'   propensity score. Default \code{NULL} (unused).
#' @param nuisance_method \code{"ml"} (default) estimates all four nuisance
#'   functions via cross-fitted random forests; \code{"parametric"}
#'   estimates them via OLS/logistic regression, assuming correct
#'   specification. See Details.
#' @param bad_control_binary Logical; whether the bad control is binary
#'   (detected automatically by \code{didbc()}). Unused by
#'   \code{dr_ml_attgt()} itself (the bad control never appears as a
#'   modeled response here), but passed through to \code{imputation_attgt()}
#'   in case of a fallback on overlap; see Details.
#' @param overlap_threshold If any unit's fitted propensity (from a
#'   preliminary, non-cross-fit fit) exceeds this value, estimation falls
#'   back to \code{imputation_attgt()} for that (g,t) cell. Default 0.99.
#' @param min_group_size If a (g,t) cell has fewer treated units than the
#'   number of p_2 covariates plus \code{min_group_size}, estimation falls
#'   back to \code{imputation_attgt()} for that cell (a separate trigger
#'   from \code{overlap_threshold}; see Details). Default 5, matching
#'   \code{did}'s own convention for an analogous check.
#' @param n_folds number of cross-fitting folds (default 5)
#' @param num_threads Number of threads for \code{grf}'s forests under
#'   \code{nuisance_method = "ml"}. Default \code{NULL} uses \code{grf}'s own
#'   auto-detection (all available cores). Set to \code{1} to pin a single
#'   thread per forest -- e.g. for Monte Carlo work, where many reps are run
#'   in outer parallelism (one core each) rather than letting each individual
#'   forest fit claim every core on the machine. Unused under
#'   \code{nuisance_method = "parametric"} (\code{lm}/\code{glm} have no
#'   thread concept).
#' @param clustervars Name of the cluster variable in \code{gt_data}, or
#'   \code{NULL} (the default). If given, cross-fitting folds are assigned
#'   by cluster; see Details.
#' @param ... additional arguments (unused)
#'
#' @return \code{attgt_if} object with ATT estimate and influence function
#'
#' @examples
#' \dontrun{
#' # dr_ml_attgt() is normally called internally as didbc()'s attgt_fun
#' # (via ptetools::pte()), not invoked directly. This is the standard way
#' # to reach it end-to-end; slow due to cross-fitted random forests plus
#' # the multiplier bootstrap, hence \dontrun.
#' sim <- simulate_bad_controls(n = 500)
#' res <- didbc(
#'   yname = "Y", gname = "G", tname = "period", idname = "id",
#'   data = sim$data,
#'   bad_control_formula = ~X,
#'   xformula = ~Z,
#'   est_method = "dr_ml",
#'   nuisance_method = "parametric"
#' )
#' summary(res)
#' }
#'
#' @references
#' Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2026).
#'   "Difference-in-Differences with Bad Controls." arXiv preprint
#'   arXiv:2608.03881. \doi{10.48550/arXiv.2608.03881}
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#'   Newey, W., and Robins, J. (2018). "Double/Debiased Machine Learning for
#'   Treatment and Structural Parameters." \emph{The Econometrics Journal}.
#'
#' Sant'Anna, P.H.C. and Zhao, J. (2020). "Doubly Robust
#'   Difference-in-Differences Estimators." \emph{Journal of Econometrics}.
#'
#' Wager, S. and Athey, S. (2018). "Estimation and Inference of
#'   Heterogeneous Treatment Effects using Random Forests." \emph{Journal
#'   of the American Statistical Association}.
#'
#' @export
dr_ml_attgt <- function(gt_data,
                        xformula = ~1,
                        bad_control_formula = NULL,
                        d_covs_formula = NULL,
                        bad_control_cov_formula = NULL,
                        bad_control_d_cov_formula = NULL,
                        nuisance_method = c("ml", "parametric"),
                        bad_control_binary = FALSE,
                        overlap_threshold = 0.99,
                        min_group_size = 5,
                        n_folds = 5,
                        num_threads = NULL,
                        clustervars = NULL,
                        ...) {

  nuisance_method <- match.arg(nuisance_method)

  if (nuisance_method == "ml" && !requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' required for nuisance_method='ml'. ",
         "Install with: install.packages('grf')")
  }

  # One row per unit, with sampling weights (mean 1) and covariate matrices;
  # transformed and factor covariates are already expanded to numeric
  # columns, which grf needs.
  cd <- bc_cell_data(
    gt_data,
    xformula = xformula,
    bad_control_formula = bad_control_formula,
    d_covs_formula = d_covs_formula,
    bad_control_cov_formula = bad_control_cov_formula,
    bad_control_d_cov_formula = bad_control_d_cov_formula
  )
  D <- cd$D
  w <- cd$w
  DeltaY <- cd$DeltaY
  has_bc <- cd$has_bc
  n <- length(D)
  n1 <- sum(D)
  n0 <- n - n1
  comparison_idx <- which(D == 0)
  treated_idx <- which(D == 1)

  # grf's forests are left unweighted when all weights are equal, so results
  # without sampling weights are exactly those of an unweighted fit
  sample_weights <- if (isTRUE(all.equal(min(w), max(w)))) NULL else w

  # Fall back to the imputation estimator for the whole cell; see Details
  # for why. Reported once per call through a ptetools cell warning.
  fallback <- function(summary, detail) {
    ptetools::pte_cell_warning(summary, detail)
    imputation_unconfoundedness(cd, bad_control_binary = bad_control_binary)
  }

  # m_0 and omega_0 depend on (X_t*, X_t*-1, Z); nu_0 and p_2 depend on
  # (X_t*-1, W, Z). Z/W are generalized to include user-supplied exogenous
  # covariates, in levels or as changes -- the same regressor sets used for
  # the analogous steps in imputation_attgt(). With no bad control, both
  # collapse to just (Z), and nu_0/omega_0 are skipped entirely below.
  if (has_bc) {
    m_features <- cbind(cd$bc_post, cd$bc_pre, cd$dX, cd$X)
    colnames(m_features)[1:2] <- paste0(cd$bc_var, c("_post", "_pre"))
    p_features <- cbind(cd$bc_pre, cd$W, cd$dW, cd$dX, cd$X)
    colnames(p_features)[1] <- paste0(cd$bc_var, "_pre")
  } else {
    m_features <- cbind(cd$dX, cd$X)
    p_features <- m_features
  }

  # Linear working models: drop collinear columns once, on the whole cell,
  # warning when m_0 (fit on comparison units) would extrapolate for treated
  # units, as the imputation estimator does. Forests need neither.
  if (nuisance_method == "parametric") {
    m_check <- bc_fit(m_features, DeltaY[comparison_idx], w, comparison_idx,
                      D)
    m_features <- m_features[, m_check$keep[-1] - 1, drop = FALSE]
    p_check <- qr(cbind(1, p_features))
    p_keep <- sort(p_check$pivot[seq_len(p_check$rank)])
    p_features <- p_features[, p_keep[-1] - 1, drop = FALSE]
  }

  # Small treated group: even under perfect overlap, fitting p_2 (and by
  # extension omega_0) needs a nontrivial treated sample to identify at all
  # relative to how many covariates it has to fit -- with too few treated
  # units, a logit can hit perfect separation and a probability_forest is
  # essentially memorizing a handful of points. Unlike the overlap check
  # below, this doesn't need an unlucky covariate configuration to show up;
  # it can happen even when the population-level covariate distributions
  # overlap fine. Falls back to the imputation estimator, since its
  # Step 1/Step 2 are fit entirely on the comparison group and only ever
  # average over treated units, never fit anything on them.
  min_treated <- ncol(p_features) + min_group_size
  if (n1 < min_treated) {
    return(fallback(
      paste0("dr_ml fell back to the imputation estimator: fewer treated ",
             "units than propensity-score covariates plus min_group_size"),
      paste0(n1, " treated, ", min_treated, " needed")
    ))
  }

  # Overlap check (not cross-fit -- a preliminary diagnostic fit, never used
  # in the final estimation). Uses the same nuisance_method as the real
  # estimation (forest or logit), since the point is to ask whether *this*
  # estimation is going to run into an overlap problem, not a generic proxy.
  # If violated, fall back to the imputation estimator for the whole cell
  # rather than dropping p_2/omega_0 alone -- see Details for why.
  overlap_fit <- fit_prop_nuisance(p_features, D, nuisance_method,
                                   num_threads, sample_weights)
  overlap_p <- predict_prop_nuisance(overlap_fit, p_features, nuisance_method)
  overlap_violation <- overlap_p > overlap_threshold
  if (any(overlap_violation)) {
    return(fallback(
      paste0("dr_ml fell back to the imputation estimator: estimated ",
             "propensity scores above overlap_threshold = ",
             overlap_threshold),
      paste0(sum(overlap_violation), " of ", n, " units")
    ))
  }

  # Cross-fitting folds. With clustervars, whole clusters are assigned to
  # folds, so that no cluster has units in both the training and the
  # held-out folds (the usual practice with dependent data; Chernozhukov et
  # al. 2018). Folds are balanced separately among treated and comparison
  # units (clusters with and without treated units).
  cl <- NULL
  if (!is.null(clustervars)) {
    cl <- gt_data[gt_data$name == "pre", clustervars]
    cl <- cl[match(cd$id, gt_data$id[gt_data$name == "pre"])]
    if (length(unique(cl[treated_idx])) < 2) {
      return(fallback(
        paste0("dr_ml fell back to the imputation estimator: cross-fitting ",
               "by cluster needs treated units in at least two clusters"),
        paste0(length(unique(cl[treated_idx])), " cluster with treated units")
      ))
    }
  }
  # drawn after the overlap check's forest, which takes its seed from R's
  # RNG, so that results for a given seed are those of earlier versions
  fold_ids <- bc_folds(D, n_folds, cl)

  m_hat <- rep(NA_real_, n)
  nu_hat <- rep(NA_real_, n)
  p_hat <- rep(NA_real_, n)
  omega_hat <- rep(NA_real_, n)

  # sampling weights of a subset of units (NULL without weights)
  sw <- function(idx) if (is.null(sample_weights)) NULL else w[idx]

  for (k in 1:n_folds) {
    train_idx <- which(fold_ids != k)
    train_comp_idx <- intersect(train_idx, comparison_idx)
    eval_idx <- which(fold_ids == k)
    if (length(eval_idx) == 0) next

    m_features_train_comp <- m_features[train_comp_idx, , drop = FALSE]
    p_features_train_comp <- p_features[train_comp_idx, , drop = FALSE]
    m_features_eval <- m_features[eval_idx, , drop = FALSE]
    p_features_eval <- p_features[eval_idx, , drop = FALSE]

    # First stage (Algorithm 1, step 2a): m_0 on untreated units, p_2 on all
    # units.
    m_fit <- fit_reg_nuisance(m_features_train_comp, DeltaY[train_comp_idx],
                              nuisance_method, num_threads,
                              sw(train_comp_idx))
    p_fit <- fit_prop_nuisance(
      p_features[train_idx, , drop = FALSE], D[train_idx], nuisance_method,
      num_threads, sw(train_idx)
    )

    # Second stage (Algorithm 1, step 2b): nu_0 regresses m_0's fitted
    # values on (X_t*-1, W, Z); omega_0 regresses p_2's fitted odds ratio on
    # (X_t*, X_t*-1, Z) -- both among untreated units in the training fold.
    # With no bad control, nu_0 = m_0 and omega_0 = p_2/(1-p_2) exactly
    # (m_0 no longer depends on X_t*, so there is nothing left for nu_0 to
    # marginalize over; p_2 is already a function of Z alone, so omega_0's
    # further conditioning on Z is a no-op), so there is nothing to fit
    # here.
    if (has_bc && nuisance_method == "ml") {
      # nu_0/omega_0 are nested conditional expectations of m_0/p_2:
      # regressing m_fit's/p_fit's fitted values on the coarser feature set
      # using the SAME sample that fit m_fit/p_fit is a generated-regressor
      # problem. Use grf's out-of-bag predictions (predict(fit), no
      # newdata) instead of predict(fit, newdata = <training data>) as the
      # pseudo-outcome -- OOB predictions for a row only average over
      # trees that didn't include that row in their subsample, so they're
      # not self-referential the way an in-sample newdata prediction is.
      m_fitted_train <- predict_reg_nuisance_oob(m_fit)
      nu_fit <- fit_reg_nuisance(p_features_train_comp, m_fitted_train,
                                 nuisance_method, num_threads,
                                 sw(train_comp_idx))

      # p_fit's OOB predictions come back for all of train_idx (treated +
      # control, in train_idx's order); subset down to the untreated units
      # (train_comp_idx) that omega_fit is trained on.
      p_fitted_all_oob <- predict_prop_nuisance_oob(p_fit)
      p_fitted_train <- p_fitted_all_oob[match(train_comp_idx, train_idx)]
      odds_fitted_train <- p_fitted_train / (1 - p_fitted_train)
      omega_fit <- fit_reg_nuisance(m_features_train_comp, odds_fitted_train,
                                    nuisance_method, num_threads,
                                    sw(train_comp_idx))
    } else if (has_bc) {
      # nuisance_method == "parametric": no OOB equivalent for linear
      # working models, per the Details section -- root-n first stages
      # don't need it anyway.
      m_fitted_train <- predict_reg_nuisance(m_fit, m_features_train_comp,
                                             nuisance_method)
      nu_fit <- fit_reg_nuisance(p_features_train_comp, m_fitted_train,
                                 nuisance_method, num_threads,
                                 sw(train_comp_idx))

      p_fitted_train <- predict_prop_nuisance(p_fit, p_features_train_comp,
                                              nuisance_method)
      odds_fitted_train <- p_fitted_train / (1 - p_fitted_train)
      omega_fit <- fit_reg_nuisance(m_features_train_comp, odds_fitted_train,
                                    nuisance_method, num_threads,
                                    sw(train_comp_idx))
    }

    # Evaluate all four nuisance functions on the held-out fold (step 3).
    m_hat[eval_idx] <- predict_reg_nuisance(m_fit, m_features_eval,
                                            nuisance_method)
    p_hat[eval_idx] <- predict_prop_nuisance(p_fit, p_features_eval,
                                             nuisance_method)
    nu_hat[eval_idx] <- if (has_bc) {
      predict_reg_nuisance(nu_fit, p_features_eval, nuisance_method)
    } else {
      m_hat[eval_idx]
    }
    omega_hat[eval_idx] <- if (has_bc) {
      predict_reg_nuisance(omega_fit, m_features_eval, nuisance_method)
    } else {
      p_hat[eval_idx] / (1 - p_hat[eval_idx])
    }
  }

  # Doubly robust score (eqn:eif1 / eq:att-dr-sample in the paper). With
  # sampling weights, each unit's score is multiplied by its weight (mean 1
  # in the cell) and pi-hat is the weighted share of treated units.
  pi_hat <- mean(w * D)
  score <- (D / pi_hat) * DeltaY - (D / pi_hat) * nu_hat -
    ((1 - D) / pi_hat) * (m_hat - nu_hat) * (p_hat / (1 - p_hat)) -
    ((1 - D) / pi_hat) * (DeltaY - m_hat) * omega_hat
  phi1 <- w * score

  att_dr <- mean(phi1)

  # Influence function (Algorithm 1, step 5): corrects for pi-hat estimation
  # uncertainty on top of the doubly robust score's own first-order terms.
  # Without weights this is phi1 - att - (att / pi) (D - pi).
  inf_func <- phi1 - att_dr * w * D / pi_hat

  extra_returns <- list(
    group = cd$g, time_period = cd$tp,
    n_control = n0, n_treated = n1,
    est_method = "dr_ml", nuisance_method = nuisance_method,
    bad_control_var = cd$bc_var, n_folds = n_folds
  )

  ptetools::attgt_if(attgt = att_dr, inf_func = inf_func,
                     extra_gt_returns = extra_returns)
}

# Cross-fitting fold of each unit. Without clusters, treated and comparison
# units are each spread evenly over the folds at random (the draws are the
# same as in earlier versions, so results for a given seed are unchanged).
# With clusters (a vector aligned with D), whole clusters are assigned to
# folds, clusters with treated units and those without each spread evenly.
#
# Arguments: D (treatment indicator), n_folds, cl (cluster of each unit, or
#   NULL).
# Returns: integer vector of folds, one per unit.
bc_folds <- function(D, n_folds, cl = NULL) {
  fold_ids <- rep(NA_integer_, length(D))
  if (is.null(cl)) {
    treated_idx <- which(D == 1)
    comparison_idx <- which(D == 0)
    fold_ids[treated_idx] <- sample(rep(1:n_folds,
                                        length.out = length(treated_idx)))
    fold_ids[comparison_idx] <- sample(rep(1:n_folds,
                                           length.out = length(comparison_idx)))
    return(fold_ids)
  }
  clusters <- unique(cl)
  has_treated <- clusters %in% unique(cl[D == 1])
  cluster_fold <- integer(length(clusters))
  for (grp in c(TRUE, FALSE)) {
    idx <- which(has_treated == grp)
    cluster_fold[idx] <- sample(rep(1:n_folds, length.out = length(idx)))
  }
  cluster_fold[match(cl, clusters)]
}

# Fit a regression nuisance function (m_0, nu_0, or omega_0): a random
# forest under nuisance_method = "ml", or OLS under
# nuisance_method = "parametric". `weights` are sampling weights (NULL for
# none).
fit_reg_nuisance <- function(x_mat, y, nuisance_method, num_threads = NULL,
                             weights = NULL) {
  if (nuisance_method == "ml") {
    return(grf::regression_forest(X = x_mat, Y = y, sample.weights = weights,
                                  num.threads = num_threads))
  }
  if (is.null(weights)) weights <- rep(1, length(y))
  bc_fit(x_mat, y, weights, seq_along(y), D = NULL, warn = FALSE)
}

predict_reg_nuisance <- function(fit, newdata_mat, nuisance_method) {
  if (nuisance_method == "ml") {
    return(predict(fit, newdata = newdata_mat)$predictions)
  }
  bc_fit_predict(fit, newdata_mat)
}

# Fit the propensity score p_2: a probability forest under
# nuisance_method = "ml", or logistic regression under
# nuisance_method = "parametric". `weights` are sampling weights (NULL for
# none).
fit_prop_nuisance <- function(x_mat, d, nuisance_method, num_threads = NULL,
                              weights = NULL) {
  if (nuisance_method == "ml") {
    return(grf::probability_forest(X = x_mat, Y = as.factor(d),
                                   sample.weights = weights,
                                   num.threads = num_threads))
  }
  if (is.null(weights)) weights <- rep(1, length(d))
  bc_fit(x_mat, d, weights, seq_along(d), D = NULL, binary = TRUE,
         warn = FALSE)
}

predict_prop_nuisance <- function(fit, newdata_mat, nuisance_method) {
  if (nuisance_method == "ml") {
    return(predict(fit, newdata = newdata_mat)$predictions[, 2])
  }
  bc_fit_predict(fit, newdata_mat)
}

# Out-of-bag predictions for a grf forest's own training data: calling
# predict() with no `newdata` returns, for each training row, an average
# over only the trees that did not include that row in their subsample --
# an out-of-sample-like prediction obtained for free from the already-fit
# forest. Used as the pseudo-outcome for nu_0/omega_0 (see the "ml" branch
# of dr_ml_attgt()'s main loop) instead of predict(fit, newdata =
# <training data>), which would be self-referential. ml-only: linear
# working models (the "parametric" nuisance_method) have no OOB equivalent.
predict_reg_nuisance_oob <- function(fit) {
  predict(fit)$predictions
}

predict_prop_nuisance_oob <- function(fit) {
  predict(fit)$predictions[, 2]
}
