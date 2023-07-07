#' Estimate controlled effect curves
#'
#' @description Estimate controlled risk (CR) curves and/or controlled vaccine
#'     efficacy (CVE) curves. See references for definitions of these curves.
#' @param dat A data object returned by load_data
#' @param type One of c("Cox", "NP"). This specifies whether to estimate the
#'     curve(s) using a marginalized Cox proportional hazards model or using a
#'     monotone-constrained nonparametric estimator.
#' @param t_0 Time point of interest
#' @param cr Boolean. If TRUE, the controlled risk (CR) curve is computed and
#'     returned.
#' @param cve Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
#'     computed and returned.
#' @param s_out A numeric vector of s-values (on the biomarker scale) for which
#'     cve(s) and/or cr(s) are computed. Defaults to a grid of 101 points
#'     between the min and max biomarker values.
#' @param ci_type One of c("transformed", "truncated", "regular", "none"). If
#'     ci_type="transformed", confidence intervals are computed on the logit(CR)
#'     and/or log(1-CVE) scale to ensure that confidence limits lie within [0,1]
#'     for CR and/or lie within (-inf,1] for CVE. If ci_type="truncated",
#'     confidence limits are constructed on the CR and/or CVE scale but
#'     truncated to lie within [0,1] for CR and/or lie within (-inf,1] for CVE.
#'     If ci_type="regular", confidence limits are not transformed or truncated.
#'     If ci_type="none", confidence intervals are not computed.
#' @param placebo_risk_method One of c("KM", "Cox"). Method for estimating
#'     overall risk in the placebo group. "KM" computes a Kaplan-Meier estimate
#'     and "Cox" computes an estimate based on a marginalized Cox model survival
#'     curve. Only relevant if cve=TRUE.
#' @param return_extras Boolean; if TRUE, objects useful for debugging are
#'     returned.
#' @param params_cox A list of options returned by
#'     \code{\link{params_ce_cox}} that are relevant if type="Cox".
#' @param params_np A list of options returned by \code{\link{params_ce_np}}
#'     that are relevant if type="NP".
#' @return A list of the form \code{list(cr=list(...), cve=list(...))}
#'     containing CR and/or CVE estimates. Each of the inner lists contains the
#'     following: \itemize{
#'         \item{\code{s}: a vector of marker values corresponding to s_out}
#'         \item{\code{est}: a vector of point estimates}
#'         \item{\code{ci_lower}: a vector of confidence interval lower limits}
#'         \item{\code{ci_upper}: a vector of confidence interval upper limits}
#' }
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578)
#' }
#' @references Gilbert P, Fong Y, Kenny A, and Carone, M (2022). A Controlled
#'     Effects Approach to Assessing Immune Correlates of Protection.
#'     <doi:10.1093/biostatistics/kxac24>
#' @export
est_ce <- function(
    dat, type="Cox", t_0, cr=TRUE, cve=FALSE,
    s_out=seq(from=min(dat$v$s,na.rm=TRUE), to=max(dat$v$s,na.rm=TRUE), l=101),
    ci_type="transformed", placebo_risk_method="KM", return_extras=FALSE,
    params_cox=params_ce_cox(), params_np=params_ce_np()
) {

  # !!!!! Input validation

  # !!!!! Move common functionality (e.g. data prep) into here

  if (type=="Cox") {
    p <- params_cox
    ests <- est_cox(
      dat=dat, t_0=t_0, cr=cr, cve=cve, s_out=s_out, ci_type=ci_type,
      placebo_risk_method=placebo_risk_method, return_extras=return_extras,
      spline_df=p$spline_df, spline_knots=p$spline_knots, edge_ind=p$edge_ind
    )
  }

  if (type=="NP") {
    p <- params_np
    params_np[["dir"]] <- NULL
    params_np[["edge_corr"]] <- NULL
    params_np[["grid_size"]] <- NULL
    ests <- est_np(
      dat=dat, t_0=t_0, cr=cr, cve=cve, s_out=s_out, ci_type=ci_type,
      placebo_risk_method=placebo_risk_method, return_extras=return_extras, dir=p$dir,
      edge_corr=p$edge_corr, params=params_np, grid_size=p$grid_size,
      cf_folds=1
    )
  }

  class(ests) <- "vaccine_est"

  return(ests)

}
