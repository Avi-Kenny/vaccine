#' Estimate CVE/CR nonparametrically
#'
#' @description Estimate controlled vaccine efficacy (CVE) and/or controlled
#'     risk (CR) using the monotone-constrainted nonparametric method of Kenny
#'     et al. 2023.
#' @param dat A data object returned by load_data
#' @param t_0 Time point of interest
#' @param cve Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
#'     computed.
#' @param cr Boolean. If TRUE, the controlled risk (CR) curve is computed.
#' @param s_out A numeric vector of s-values (on the biomarker scale) for which
#'     cve(s) and/or cr(s) are computed. Defaults to a grid of 101 points
#'     between the min and max biomarker values.
#' @param params A list of key value pairs that can be used to select or tune
#'     various nuisance estimators. See examples. These include: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @param verbose A Boolean. If set to TRUE, intermediate output will be
#'     displayed.
#' @return A list containing the following: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @examples
#' print("to do")
#' #list(cond_density_type="binning", cond_density_nbins=10)
#' @export
#' @references Kenny, A., Gilbert P., and Carone, M. (2023). Nonparametric
#'     inference for the controlled risk and controlled vaccine efficacy curves.
#' @references Gilbert P., Fong Y., Kenny A., and Carone, M. (2022). A
#'     Controlled Effects Approach to Assessing Immune Correlates of Protection.
#' @export
est_np <- function(
  dat, cve=T, cr=T, s_out=seq(from=min(dat$s), to=max(dat$s), l=101),
  params=list(), verbose=F
) {

  dat_orig <- dat # !!!!! Maybe change this later
  if (class(dat_orig)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  # Fix s_out if needed
  if (any(is.na(dat$s))) {
    if (missing(s_out)) {
      s_out <- seq(from=min(dat$s, na.rm=T), to=max(dat$s, na.rm=T), l=101)
    } else {
      stop("NA values not allowed in s_out.")
    }
  }

  # Set params
  .default_params <- list(
    # !!!!! Rename all of these and condense formatting
    Q_n_type="Super Learner",
    g_n_type="binning",
    deriv_type="linear",
    ecdf_type="linear (mid)",
    gamma_type="Super Learner",
    q_n_type="new",
    omega_n_type="estimated",
    boot_reps=1000,
    ci_type="logit",
    cf_folds=1,
    m=5,
    edge_corr="none",
    lod_shift="none",
    n_bins=5,
    convex_type="GCM",
    f_sIx_n_bins=15
  )
  for (i in c(1:length(.default_params))) {
    p_name <- names(.default_params)[i]
    if (is.null(params[[p_name]])) {
      params[[p_name]] <- .default_params[[i]]
    }
  }
  p <- params
  rm(params)

  # !!!!! CONTINUE

  return(999)

}
