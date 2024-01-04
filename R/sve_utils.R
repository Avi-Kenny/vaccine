#' Stochastic risk evaluation for immune correlates analyses
#'
#' @param ve_dat A rectangular data object (e.g., \code{data.frame}) containing
#'  data on study units from a vaccine efficacy trial. See the arguments below
#'  for details on the variables that ought to be included in \code{ve_dat}.
#' @param Y An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to the outcome of interest (an infection/disease endpoint);
#'  currently, this is limited to binary endpoints. Future extensions include
#'  support for time-to-event endpoints subject to right-censoring.
#' @param C An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of whether a study unit was lost to follow-up
#'  during the VE trial. The default of \code{NA_character_} assumes no loss to
#'  follow-up.
#' @param S An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to the measurements of an immunologic marker of interest's
#'  measured activity post-vaccination.
#' @param An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of vaccination status or of assignment to the
#'  vaccine versus placebo arm in the VE trial.
#' @param W A vector \code{character} naming those variables in \code{ve_dat}
#'  that correspond to measured baseline covariates that are to comprise the
#'  adjustment set when computing estimates of nuisance parameters.
#' @param B A vector \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of inclusion in a designed two-phase sample
#'  (e.g., case-control or case-cohort design). This is used to compute inverse
#'  probability of sampling weights to recover population-level inference.
#' @param delta_grid A vector \code{numeric} of values by which measurements in
#'  the immune marker specified by \code{S} should be mean-shifted to obtain
#'  counterfactual risk estimates.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator. Both estimators are semiparametric-efficient and compatible with
#'  machine learning-based estimation of nuisance parameters.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted", which correspond to how the auxiliary covariate is included
#'  in the logistic tilting regression model in the TML updating step.
#' @param weighting Whether to weight each parameter estimate by the inverse of
#'  its variance (in order to improve stability of the resultant MSM fit) or to
#'  simply weight all parameter estimates equally. The default is the option
#'  \code{"identity"}, weighting all estimates identically.
#' @param ci_level A \code{numeric} indicating the desired coverage level of a
#'  Wald-style confidence interval to be computed for the parameter estimate.
#' @param ci_type Whether to construct a simultaneous confidence band covering
#'  all parameter estimates at once or marginal confidence intervals covering
#'  each parameter estimate separately. The default is to construct marginal
#'  confidence intervals for each parameter estimate rather than a simultaneous
#'  confidence band.
#' @param gps_bound An atomic \code{numeric} specifying the lower limit of the
#'  generalized propensity score estimates to be tolerated (with a default of
#'  0.01). Estimates below this are truncated to max(1/n, \code{gps_bound}).
#' @param gps_est ...
#' @param or_est ...
#' @param tps_est ...
#' @param eif_est ...
#'
#' @importFrom data.table setDT
#' @importFrom txshift msm_vimshift
#' @importFrom sl3 Lrnr_glm_fast Lrnr_density_semiparametric
est_mtprisk <- function(ve_dat,
                        Y,
                        C = NA_character_,
                        S,
                        A,
                        W,
                        B,
                        delta_grid = seq(-0.5, 0.5, by = 0.25),
                        estimator = c("tmle", "onestep"),
                        fluctuation = c("weighted", "standard"),
                        weighting = c("identity", "variance"),
                        ci_level = 0.95,
                        ci_type = c("marginal", "simultaneous"),
                        gps_bound = 0.01,
                        gps_est = sl3::Lrnr_density_semiparametric$new(
                          mean_learner = sl3::Lrnr_glm_fast$new()
                        ),
                        or_est = sl3::Lrnr_glm_fast$new(),
                        tps_est = sl3::Lrnr_glm_fast$new(),
                        eif_est = c("hal", "glm")
                       ) {
  # set defaults for estimation options
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  weighting <- match.arg(weighting)
  ci_type <- match.arg(ci_type)
  eif_est <- match.arg(eif_est)

  # convert input VE trial rectangular data object to DT
  data.table::setDT(ve_dat)

  # check whether censoring variable exists in input data
  if (is.na(C)) {
    C_cens <- rep(1, ve_dat[get(A) == 1, .N])
  } else {
    C_cens <- ve_dat[get(A) == 1, get(C)]
  }

  # check input tps_est: if vector, then assume to be known weights; if an sl3
  # learner, then should be used for estimation of IP sampling weights (below)
  if (is.vector(tps_est)) {
    # compute stochastic interventional risk with MSM summary, in this case
    # using *known* weights
    cop_risk_msm <- txshift::msm_vimshift(
      Y = ve_dat[get(A) == 1, get(Y)],
      C_cens = C_cens,
      A = ve_dat[get(A) == 1, get(S)],
      W = ve_dat[get(A) == 1, ..W],
      C_samp = ve_dat[get(A) == 1, get(B)],
      V = c("W", "Y"),
      delta_grid = delta_grid,
      msm_form = list(type = "linear", knot = NA),
      estimator = estimator,
      fluctuation = fluctuation,
      weighting = weighting,
      ci_level = ci_level,
      ci_type = ci_type,
      # NOTE: arguments passed through to txshift::txshift()
      gps_bound = gps_bound,
      samp_fit_args = list(
        fit_type = "external"
      ),
      g_exp_fit_args = list(
        fit_type = "sl",
        sl_learners_density = gps_est
      ),
      Q_fit_args = list(
        fit_type = "sl",
        sl_learners = or_est
      ),
      eif_reg_type = eif_est,
      # NOTE: need to pass in sampling probabilities, not weights
      samp_fit_ext = 1 / tps_est
    )
  } else {
    # compute stochastic interventional risk with MSM summary, in this case
    # using *estimated* IP sampling weights
    cop_risk_msm <- txshift::msm_vimshift(
      Y = ve_dat[get(A) == 1, get(Y)],
      C_cens = C_cens,
      A = ve_dat[get(A) == 1, get(S)],
      W = ve_dat[get(A) == 1, ..W],
      C_samp = ve_dat[get(A) == 1, get(B)],
      V = c("W", "Y"),
      delta_grid = delta_grid,
      msm_form = list(type = "linear", knot = NA),
      estimator = estimator,
      fluctuation = fluctuation,
      weighting = weighting,
      ci_level = ci_level,
      ci_type = ci_type,
      # NOTE: arguments passed through to txshift::txshift()
      gps_bound = gps_bound,
      samp_fit_args = list(
        fit_type = "sl",
        sl_learners = tps_est
      ),
      g_exp_fit_args = list(
        fit_type = "sl",
        sl_learners_density = gps_est
      ),
      Q_fit_args = list(
        fit_type = "sl",
        sl_learners = or_est
      ),
      eif_reg_type = eif_est
    )
  }
}


#' Stochastic vaccine efficacy evaluation for immune correlates analyses
#'
#' @param ve_dat A rectangular data object (e.g., \code{data.frame}) containing
#'  data on study units from a vaccine efficacy trial. See the arguments below
#'  for details on the variables that ought to be included in \code{ve_dat}.
#' @param Y An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to the outcome of interest (an infection/disease endpoint);
#'  currently, this is limited to binary endpoints. Future extensions include
#'  support for time-to-event endpoints subject to right-censoring.
#' @param C An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of whether a study unit was lost to follow-up
#'  during the VE trial. The default of \code{NA_character_} assumes no loss to
#'  follow-up.
#' @param S An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to the measurements of an immunologic marker of interest's
#'  measured activity post-vaccination.
#' @param An atomic \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of vaccination status or of assignment to the
#'  vaccine versus placebo arm in the VE trial.
#' @param W A vector \code{character} naming those variables in \code{ve_dat}
#'  that correspond to measured baseline covariates that are to comprise the
#'  adjustment set when computing estimates of nuisance parameters.
#' @param B A vector \code{character} naming a variable in \code{ve_dat} that
#'  corresponds to an indicator of inclusion in a designed two-phase sample
#'  (e.g., case-control or case-cohort design). This is used to compute inverse
#'  probability of sampling weights to recover population-level inference.
#' @param delta_grid A vector \code{numeric} of values by which measurements in
#'  the immune marker specified by \code{S} should be mean-shifted to obtain
#'  counterfactual risk estimates.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator. Both estimators are semiparametric-efficient and compatible with
#'  machine learning-based estimation of nuisance parameters.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted", which correspond to how the auxiliary covariate is included
#'  in the logistic tilting regression model in the TML updating step.
#' @param weighting Whether to weight each parameter estimate by the inverse of
#'  its variance (in order to improve stability of the resultant MSM fit) or to
#'  simply weight all parameter estimates equally. The default is the option
#'  \code{"identity"}, weighting all estimates identically.
#' @param ci_level A \code{numeric} indicating the desired coverage level of a
#'  Wald-style confidence interval to be computed for the parameter estimate.
#' @param ci_type Whether to construct a simultaneous confidence band covering
#'  all parameter estimates at once or marginal confidence intervals covering
#'  each parameter estimate separately. The default is to construct marginal
#'  confidence intervals for each parameter estimate rather than a simultaneous
#'  confidence band.
#' @param gps_bound Atomic (\code{numeric}) specifying the lower limit of the
#'  generalized propensity score estimates to be tolerated (with a default of
#'  0.01). Any estimates below this are truncated to max(1/n, this bound).
#' @param gps_est ...
#' @param or_est ...
#' @param tps_est ...
#' @param eif_est ...
#'
#' @importFrom data.table as.data.table copy setnames setDT
#' @importFrom stats as.formula cov glm model.matrix pnorm qnorm gaussian
#' @importFrom matrixStats colVars
#' @importFrom assertthat assert_that
#' @importFrom txshift msm_vimshift
#' @importFrom sl3 Lrnr_glm_fast Lrnr_density_semiparametric
est_sve <- function(ve_dat,
                    Y,
                    C = NA_character_,
                    S,
                    A,
                    W,
                    B,
                    delta_grid = seq(-0.5, 0.5, by = 0.25),
                    estimator = c("tmle", "onestep"),
                    fluctuation = c("weighted", "standard"),
                    weighting = c("identity", "variance"),
                    ci_level = 0.95,
                    ci_type = c("marginal", "simultaneous"),
                    gps_bound = 0.01,
                    gps_est = sl3::Lrnr_density_semiparametric$new(
                      mean_learner = sl3::Lrnr_glm_fast$new()
                    ),
                    or_est = sl3::Lrnr_glm_fast$new(),
                    tps_est = sl3::Lrnr_glm_fast$new(),
                    eif_est = c("hal", "glm")
                   ) {
  # set defaults for estimation options
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  weighting <- match.arg(weighting)
  ci_type <- match.arg(ci_type)
  eif_est <- match.arg(eif_est)

  # convert input VE trial rectangular data object to DT
  data.table::setDT(ve_dat)

  # compute risk estimate for numerator of SVE
  cop_risk_est <- est_mtprisk(
    ve_dat = ve_dat,
    Y = Y,
    C = C,
    S = S,
    A = A,
    W = W,
    B = B,
    delta_grid = delta_grid,
    estimator = estimator,
    fluctuation = fluctuation,
    weighting = weighting,
    ci_level = ci_level,
    ci_type = ci_type,
    gps_bound = gps_bound,
    gps_est = gps_est,
    or_est = or_est,
    tps_est = tps_est,
    eif_est = eif_est
  )

  # get P(Y(0) = 1) = P(Y = 1 | A = 0), by randomization
  placebo_risk_est <- ve_dat[get(A) == 0, mean(get(Y))]
  placebo_risk_eif <- rep(0, ve_dat[, .N])
  placebo_risk_eif[ve_dat[, which(get(A) == 0)]] <-
    ve_dat[get(A) == 0, get(Y)] - placebo_risk_est

  # extract useful components from txshift_msm object
  cop_vax_risk_eif <- replicate(length(delta_grid), rep(0, ve_dat[, .N]))
  cop_vax_risk_eif[ve_dat[, which(get(A) == 1)], ] <- cop_risk_est$.eif_mat

  # apply the delta method to get the EIFs for SVE for each delta
  # NOTE: letting psi_1: interventional vaccinee risk, psi_0: control risk
  #       EIF(SVE) = [EIF(psi_1) / psi_1] - [EIF(psi_0) / psi_0]
  cop_sve_eifs <- lapply(seq_along(delta_grid), function(delta_idx) {
    (cop_vax_risk_eif[, delta_idx] / cop_risk_est$param_est$psi[delta_idx]) -
      (placebo_risk_eif / placebo_risk_est)
  })
  cop_sve_eifs <- do.call(cbind, cop_sve_eifs)

  # compute standard error based on SVE EIFs
  cop_sve_stderr <- sqrt(matrixStats::colVars(cop_sve_eifs) / ve_dat[, .N])

  # compute SVE estimate and Wald-type confidence intervals on log-scale
  # NOTE: signs for CIs flipped due to working on the log-scale
  wald_ci_mult <- abs(stats::qnorm(p = (1 - cop_risk_est$.ci_level) / 2))
  stderr_ci_mult <- cop_sve_stderr * wald_ci_mult
  cop_sve_psi_logscale <- log(cop_risk_est$param_est$psi / placebo_risk_est)
  cop_sve_lwrci_logscale <- cop_sve_psi_logscale + stderr_ci_mult
  cop_sve_uprci_logscale <- cop_sve_psi_logscale - stderr_ci_mult

  # summarize output
  cop_sve_est <- data.table::copy(cop_risk_est$param_est)
  cop_sve_est[, ci_lwr := 1 - exp(cop_sve_lwrci_logscale)]
  cop_sve_est[, psi := 1 - exp(cop_sve_psi_logscale)]
  cop_sve_est[, ci_upr := 1 - exp(cop_sve_uprci_logscale)]

  # construct MSM point estimate on logRR scale
  if (weighting == "variance") {
    weights_eif <- as.numeric(1 / diag(stats::cov(cop_sve_eifs)))
  } else if (weighting == "identity") {
    weights_eif <- rep(1, ncol(cop_sve_eifs))
  }
  x_mat <- stats::model.matrix(
    stats::as.formula(paste0("psi_vec ~ delta")),
    data = data.frame(psi_vec = cop_sve_psi_logscale, delta = delta_grid)
  )
  omega <- diag(weights_eif)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% cop_sve_psi_logscale)

  # get EIFs, variance, SE for MSM parameter (estimates)
  msm_eif <- tcrossprod(cop_sve_eifs, s_mat)
  msm_var <- diag(stats::cov(msm_eif))
  msm_stderr <- sqrt(msm_var / ve_dat[, .N])

  browser()
  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- msm_stderr %*% t(c(-wald_ci_mult, wald_ci_mult)) + msm_param
  stat_msm_param <- msm_param / msm_stderr
  pval_msm_param <- 2 * stats::pnorm(abs(stat_msm_param), lower.tail = FALSE)

  # create summary table for MSM estimates
  msm_results <- data.table::as.data.table(list(
    param = names(msm_stderr),
    ci_lwr = exp(ci_msm_param[, 1]),
    param_est = exp(msm_param),
    ci_upr = exp(ci_msm_param[, 2]),
    param_se = msm_stderr,
    p_value = pval_msm_param
  ))

  # create and rename MSM data for downstream ggplot2 compatibility
  msm_fit <- stats::glm(
    formula = stats::as.formula("y ~ x"), weights = weights_eif,
    family = stats::gaussian(), start = c(0, 0),
    data = data.frame(y = cop_sve_psi_logscale, x = delta_grid)
  )
  # CHECK: GLM coefficients should match exactly with manually calculations
  assertthat::assert_that(all.equal(unname(stats::coef(msm_fit)), msm_param))

  # predict from fitted GLM to create data for plotting MSM VE curves
  delta_grid_seq <- data.frame(
    x = seq(min(delta_grid), max(delta_grid), length.out = 1000L)
  )
  rr_logscale_msm_pred <- stats::predict(msm_fit, delta_grid_seq)
  ve_msm_pred <- 1 - exp(rr_logscale_msm_pred)
  msm_plot_data <- data.table::as.data.table(
    list(psi = ve_msm_pred, delta = delta_grid_seq)
  )

  # overwrite txshift_msm class from txshift::msm_vimshift with SVE results
  cop_sve_results <- cop_risk_est
  cop_sve_results$param_est <- cop_sve_est
  cop_sve_results$msm_est <- msm_results
  cop_sve_results$.msm_data <- msm_plot_data
  cop_sve_results$msm_fit <- msm_fit
  return(cop_sve_results)
}

#' Plot of working marginal structural model for stochastic vaccine efficacy
#'
#' @details Creates a visualization of the intervention-specific counterfactual
#'  means on the log(relative risk) scale as well as the working marginal
#'  structural model summarizing the trend across posited values of the
#'  intervention.
#'
#' @param x Object of class \code{txshift_msm} as produced by a call to
#'  \code{\link{msm_vimshift}}.
#' @param ... Additional arguments passed to \code{plot} as necessary.
#'
#' @import ggplot2
#' @importFrom latex2exp TeX
plot_sve <- function(x, ...) {
  # build geom for MSM in plot
  # NOTE: error bars for marginal CIs but band for simultaneous CIs
  if (x$.ci_type == "marginal") {
    geom_ci <- geom_errorbar(
      data = x$param_est,
      aes_string(
        ymin = "ci_lwr",
        ymax = "ci_upr"
      ),
      linetype = "dotted",
      width = 0.05
    )
  } else if (x$.ci_type == "simultaneous") {
    geom_ci <- geom_ribbon(
      data = x$param_est,
      aes_string(
        ymin = "ci_lwr",
        ymax = "ci_upr"
      ),
      fill = "grey",
      alpha = 0.3
    )
  }

  # create plot
  p_msm <- ggplot(x$param_est, aes_string(x = "delta", y = "psi")) +
    geom_point(size = 3, alpha = 0.75) +
    geom_ci +
    geom_line(
      data = x$.msm_data, aes_string(x = "delta", y = "psi"),
      linetype = "dashed", size = 2, color = "black"
    ) +
    labs(
      x = latex2exp::TeX("Shift in exposure $\\delta$"),
      y = latex2exp::TeX("Counterfactual mean $EY_{A + \\delta(W)}$"),
      title = "Estimated mean counterfactual outcome under shifted exposure",
      subtitle = paste(
        "with", x$.ci_type, "confidence intervals and",
        x$.msm_type, "working MSM summarization"
      )
    ) +
    theme_bw()

  # output plot
  return(p_msm)
}
