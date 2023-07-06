#' Set parameters controlling Cox model estimation of controlled effect curves
#'
#' @description This should be used in conjunction with \code{\link{est_ce}} to
#'     set parameters controlling Cox model estimation of controlled effect
#'     curves; see examples.
#' @param spline_df An integer; if the marker is modeled flexibly within the Cox
#'     model linear predictor as a natural cubic spline, this option controls
#'     the degrees of freedom in the spline; knots are chosen to be equally
#'     spaced across the range of the marker.
#' @param spline_knots A numeric vector; as an alternative to specifying
#'     \code{spline_df}, the exact locations of the knots in the spline
#'     (including boundary knots) can be specified with this option.
#' @param edge_ind Boolean. If TRUE, an indicator variable corresponding to the
#'     lower limit of the marker will be included in the Cox model linear
#'     predictor.
#' @return A list of options.
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_cox <- est_ce(
#'   dat = dat,
#'   type = "Cox",
#'   t_0 = 578,
#'   params_cox = params_ce_cox(spline_df=4)
#' )
#' }
#' @export
params_ce_cox <- function(spline_df=NA, spline_knots=NA, edge_ind=F) {

  return(list(spline_df=spline_df, spline_knots=spline_knots,
              edge_ind=edge_ind))

}



#' Set parameters controlling nonparametric estimation of controlled effect
#'     curves
#'
#' @description This should be used in conjunction with \code{\link{est_ce}} to
#'     set parameters controlling nonparametric estimation of controlled effect
#'     curves; see examples.
#' @param dir One of c("decr", "incr"); controls the direction of monotonicity.
#'     If dir="decr", it is assumed that CR decreases as a function of the
#'     marker. If dir="incr", it is assumed that CR increases as a function of
#'     the marker.
#' @param edge_corr Boolean. If TRUE, the "edge correction" is performed to
#'     adjust for bias near the marker lower limit (see references). It is
#'     recommended that the edge correction is only performed if there are at
#'     least (roughly) 10 events corresponding to the marker lower limit.
#' @param grid_size A list with keys \code{y}, \code{s}, and \code{x}; controls
#'     the rounding of data values. Decreasing the grid size values results in
#'     shorter computation times, and increasing the values results in more
#'     precise estimates. If grid_size$s=101, this means that a grid of 101
#'     equally-spaced points (defining 100 intervals) will be created from
#'     min(S) to max(S), and each S value will be rounded to the nearest grid
#'     point. For grid_size$y, a grid will be created from 0 to t_0, and then
#'     extended to max(Y). For grid_size$x, a separate grid is created for each
#'     covariate column (binary/categorical covariates are ignored).
#' @param surv_type One of c("Cox", "survSL", "survML-G", "survML-L"); controls
#'     the method to use to estimate the conditional survival and conditional
#'     censoring functions. If type="Cox", a survival function based on a Cox
#'     proportional hazard model will be used. If type="survSL", the Super
#'     Learner method of Westling 2023 (via the \code{survSuperLearner} package)
#'     is used. If type="survML-G", the global survival stacking method of
#'     Wolock 2022 (via the \code{survML} package) is used. If type="survML-L",
#'     the local survival stacking method of Polley 2011 (via the \code{survML}
#'     package) is used.
#' @param density_type One of c("binning", "parametric"); controls the method to
#'     use to estimate the density ratio f(S|X)/f(S).
#' @param density_bins An integer; if density_type="binning", the number of bins
#'     to use. If density_bins=0, the number of bins will be selected via
#'     cross-validation.
#' @param deriv_type One of c("m-spline", "linear"); controls the method to use
#' to estimate the derivative of the CR curve. If deriv_type="linear", a linear
#' spline is constructed based on the midpoints of the jump points of the
#' estimated function (plus the estimated function evaluated at the endpoints),
#' which is then numerically differentiated. deriv_type="m-spline" is similar to
#' deriv_type="linear" but smooths the set of points (using the method of
#' Fritsch and Carlson 1980) before differentiating.
#' @return A list of options.
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_np <- est_ce(
#'   dat = dat,
#'   type = "NP",
#'   t_0 = 578,
#'   params_np = params_ce_np(edge_corr=T, surv_type="survML-L")
#' )
#' }
#' @export
params_ce_np <- function(
    dir = "decr",
    edge_corr = F,
    grid_size = list(y=101,s=101,x=5),
    surv_type = "survML-G",
    density_type = "binning",
    density_bins = 15,
    deriv_type = "m-spline"
) {
  return(list(dir=dir, edge_corr=edge_corr, grid_size=grid_size,
              surv_type=surv_type, density_type=density_type,
              density_bins=density_bins, deriv_type=deriv_type))
}
