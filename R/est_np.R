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
#' @param ci_type TO DO
#' @param cf_folds An integer representing the number of cross-fitting folds to
#'     use. If cf_folds=1, cross-fitting will not be done.
#' @param edge_corr Boolean. If TRUE, the edge correction is performed. It is
#'     not recommended that the edge correction is performed unless there are
#'     at least 10 events corresponding to the marker lower limit
#' @param params A list of key value pairs that can be used to select or tune
#'     various nuisance estimators. See examples. These include: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @param grid_size A list containing the following three keys: \itemize{
#'     \item{\code{y}: grid size for time values}
#'     \item{\code{s}: grid size for marker values}
#'     \item{\code{x}: grid size for covariate values}
#' }
#'     This controls rounding of covariates, which results in substantially
#'     shorter computation times. If grid_size$s=100, this means that a grid of
#'     equally-spaced points will be created from min(dat$s) to max(dat$s), and
#'     each dat$s value will be rounded to the nearest grid point. For
#'     grid_size$y, a grid will be created from 0 to max(c(dat$y,t_0)). For
#'     grid_size$x, a separate grid is created for each covariate column
#'     (binary and categorical covariates are ignored).
#' @param return_extras Boolean. If set to TRUE, the following quantities (most
#'     of which are mainly useful for debugging) are returned: \itemize{
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
#' @note
#'   - This method assumes that risk decreases as the biomarker increases. If it
#'     is assumed that risk increases as the biomarker increases, reverse the
#'     scale of the biomarker.
#' @references Kenny, A., Gilbert P., and Carone, M. (2023). Nonparametric
#'     inference for the controlled risk and controlled vaccine efficacy curves.
#' @references Gilbert P., Fong Y., Kenny A., and Carone, M. (2022). A
#'     Controlled Effects Approach to Assessing Immune Correlates of Protection.
#' @export
est_np <- function(
  dat, t_0, cve=T, cr=T, s_out=seq(from=min(dat$s), to=max(dat$s), l=101),
  ci_type="logit", cf_folds=1, edge_corr=F, params=list(),
  grid_size=list(y=101, s=101, x=5), return_extras=F, verbose=F
) {

  if (class(dat)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  # !!!!! Validate other inputs; import error handling function from SimEngine

  # Alias variables
  dat_orig <- dat$df_vc # !!!!! Maybe change this later
  .v <- verbose

  # Fix s_out if needed
  if (F) { # !!!!! if clause is temporary
    if (any(is.na(dat$s))) {
      if (missing(s_out)) {
        s_out <- seq(from=min(dat$s, na.rm=T), to=max(dat$s, na.rm=T), l=101)
      } else {
        stop("NA values not allowed in s_out.")
      }
    }
  }

  # Set params
  .default_params <- list(
    surv_type = "Super Learner",
    density_type = "binning",
    density_bins = 15,
    deriv_type = "m-spline",
    gamma_type = "Super Learner",
    q_n_type = "standard",
    convex_type = "GCM"
  )
  for (i in c(1:length(.default_params))) {
    p_name <- names(.default_params)[i]
    if (is.null(params[[p_name]])) { params[[p_name]] <- .default_params[[i]] }
  }
  p <- params
  p$ci_type <- ci_type
  p$cf_folds <- cf_folds
  p$edge_corr <- edge_corr

  # Rescale S to lie in [0,1] and round values
  s_min <- min(dat_orig$s, na.rm=T)
  s_max <- max(dat_orig$s, na.rm=T)
  s_shift <- -1 * s_min
  s_scale <- 1/(s_max-s_min)
  dat_orig$s <- (dat_orig$s+s_shift)*s_scale
  grid <- create_grid(dat_orig, grid_size, t_0)
  dat_orig <- round_dat(dat_orig, grid, grid_size)

  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr) { s_min2 <- min(dat_orig$s[dat_orig$s!=0], na.rm=T) }

  # Rescale/round s_out and remove s_out points outside [0,1]
  s_out_orig <- s_out
  s_out <- (s_out+s_shift)*s_scale
  s_out <- sapply(s_out, function(s) { grid$s[which.min(abs(grid$s-s))] })
  na_head <- sum(s_out<0)
  na_tail <- sum(s_out>1)
  if (na_head>0) { s_out <- s_out[-c(1:na_head)] }
  len_p <- length(s_out)
  if (na_tail>0) { s_out <- s_out[-c((len_p-na_tail+1):len_p)] }

  # Create phase-two data object
  dat <- ss(dat_orig, which(dat_orig$z==1))

  # Prepare precomputation values for conditional survival estimator
  x_distinct <- dplyr::distinct(dat_orig$x)
  x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
  vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
  vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
  vals <- list(
    t = vals_pre$t,
    x = subset(vals_pre, select=-c(t,x_index,s)),
    s = vals_pre$s
  )

  # Fit conditional survival estimator
  srvSL <- construct_Q_n(p$surv_type, dat, vals)
  Q_n <- srvSL$srv
  Qc_n <- srvSL$cens

  # Compute various nuisance functions
  omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
  f_sIx_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins, z1=F)
  f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
  g_n <- construct_g_n(f_sIx_n, f_s_n)
  Phi_n <- construct_Phi_n(ss(dat, which(dat$s!=0))) # !!!!! Make sure Phi_n(1)=1
  n_orig <- attr(dat_orig, "n_orig")
  p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))
  eta_n <- construct_eta_n(dat, Q_n, p_n, t_0)
  r_tilde_Mn <- construct_r_tilde_Mn(dat_orig, Q_n, t_0)
  Gamma_tilde_n <- construct_Gamma_tilde_n(dat, r_tilde_Mn, p_n)
  f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)

  q_n <- construct_q_n(type=p$q_n_type, dat, omega_n, g_n, r_tilde_Mn,
                       Gamma_tilde_n, f_n_srv)

  # !!!!! Profiling q_n
  if (F) {

    # p$q_n_type <- "zero"
    p$q_n_type <- "standard"
    q_n <- construct_q_n(type=p$q_n_type, dat, omega_n, g_n, r_tilde_Mn,
                         Gamma_tilde_n, f_n_srv)
    dat_orig_df <- as_df(dat_orig)
    u <- 0.5; dim_x <- 2;

    # Profile this line
    q_n_do <- as.numeric(apply(dat_orig_df, 1, function(r) {
      q_n(as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]], u)
    }))

  }

  Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n,
                                     q_n, r_tilde_Mn, Gamma_tilde_n)

  # !!!!! Profiling Gamma_os_n
  if (F) {

    # p$q_n_type <- "zero"
    p$q_n_type <- "standard"
    q_n <- construct_q_n(type=p$q_n_type, dat, omega_n, g_n, r_tilde_Mn,
                         Gamma_tilde_n, f_n_srv)
    Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n,
                                       q_n, r_tilde_Mn, Gamma_tilde_n)

    microbenchmark({
      Gamma_os_n(0.6)
    }, times=10L)

  }

  # Compute edge-corrected estimator and standard error
  if (p$edge_corr) {

    g_sn <- construct_g_sn(dat, f_n_srv, g_n, p_n)
    r_Mn_edge_est <- r_Mn_edge(dat_orig, g_sn, g_n, p_n, Q_n, omega_n, t_0)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
                                                     r_Mn_edge_est, p_n, t_0)
    dat_orig_df <- as_df(dat_orig)
    dim_x <- attr(dat_orig, "dim_x")
    sigma2_edge_est <- mean(apply(dat_orig_df, 1, function(r) {
      (infl_fn_r_Mn_edge(r[["z"]], r[["weights"]], r[["s"]],
                         as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]]))^2
    }))

  }

  # Compute GCM (or least squares line) and extract its derivative
  gcm_x_vals <- sapply(sort(unique(dat$s)), Phi_n)
  indices_to_keep <- !base::duplicated(gcm_x_vals)
  gcm_x_vals <- gcm_x_vals[indices_to_keep]
  gcm_y_vals <- -1 * sapply(sort(unique(dat$s))[indices_to_keep], Gamma_os_n)
  if (!any(gcm_x_vals==0)) {
    gcm_x_vals <- c(0, gcm_x_vals)
    gcm_y_vals <- c(0, gcm_y_vals)
  }
  if (p$convex_type=="GCM") {
    gcm <- fdrtool::gcmlcm(x=gcm_x_vals, y=gcm_y_vals, type="gcm")
    dGCM <- approxfun(
      x = gcm$x.knots[-length(gcm$x.knots)],
      y = gcm$slope.knots,
      method = "constant",
      rule = 2,
      f = 0
    )
  } else if (p$convex_type=="CLS") {
    gcm <- function(x) { 1 } # Ignored
    fit <- cvx.lse.reg(t=gcm_x_vals, z=gcm_y_vals)
    pred_x <- round(seq(0,1,0.001),3)
    pred_y <- predict(fit, newdata=pred_x)
    dGCM <- function(u) {
      width <- 0.05
      u1 <- u - width/2; u2 <- u + width/2;
      if (u1<0) { u2 <- u2 - u1; u1 <- 0; }
      if (u2>1) { u1 <- u1 - u2 + 1; u2 <- 1; }
      u1 <- round(u1,3); u2 <- round(u2,3);
      ind1 <- which(pred_x==u1); ind2 <- which(pred_x==u2);
      y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
      return((y2-y1)/width)
    }
  }

  # Construct Grenander-based r_Mn estimator (and truncate to lie within [0,1])
  r_Mn_Gr <- function(u) { min(max(-1*dGCM(Phi_n(u)),0),1) }

  # Compute variance component nuisance estimators
  f_sIx_z1_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins,
                                  z1=T)
  f_s_z1_n <- construct_f_s_n(dat_orig, f_sIx_z1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, type="Super Learner", omega_n,
                               grid)
  g_zn <- construct_g_zn(dat_orig, type="Super Learner", f_sIx_n, f_sIx_z1_n)

  # Create either regular or edge-corrected r_Mn estimator
  if (p$edge_corr) {

    r_Mn <- function(u) {
      if(u==0 || u<s_min2) {
        return(r_Mn_edge_est)
      } else {
        return(min(r_Mn_edge_est, r_Mn_Gr(u)))
      }
    }

  } else {

    r_Mn <- r_Mn_Gr

  }

  # Generate estimates for each point
  ests <- sapply(s_out, r_Mn)

  # Construct variance scale factor
  deriv_r_Mn <- construct_deriv_r_Mn(type=p$deriv_type, r_Mn, grid)
  tau_n <- construct_tau_n(deriv_r_Mn, gamma_n, f_sIx_n, g_zn, dat_orig)

  # Generate confidence limits
  if (p$ci_type=="none") {

    ci_lo <- rep(NA, length(ests))
    ci_hi <- rep(NA, length(ests))

  } else {

    # Generate variance scale factor for each point
    tau_ns <- sapply(s_out, tau_n)

    # Construct CIs
    # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
    qnt <- 1.00
    if (p$ci_type=="regular") {
      ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
      ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
    } else if (p$ci_type=="logit") {
      ci_lo <- expit(
        logit(ests) - (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
      )
      ci_hi <- expit(
        logit(ests) + (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
      )
    }

    # CI edge correction
    if (p$edge_corr) {
      ci_lo2 <- ests[1] - 1.96*sqrt(sigma2_edge_est/n_orig)
      ci_hi2 <- ests[1] + 1.96*sqrt(sigma2_edge_est/n_orig)
      ci_lo <- In(r_Mn_edge_est<=ests)*pmin(ci_lo,ci_lo2) +
        In(r_Mn_edge_est>ests)*ci_lo
      ci_lo[1] <- ci_lo2
      ci_hi <- In(r_Mn_edge_est<=ests)*pmin(ci_hi,ci_hi2) +
        In(r_Mn_edge_est>ests)*ci_hi
      ci_hi[1] <- ci_hi2
    }

  }

  # Create results object
  res <- list(
    s = s_out_orig,
    est = c(rep(NA,na_head), ests, rep(NA,na_tail)),
    ci_lo = c(rep(NA,na_head), ci_lo, rep(NA,na_tail)),
    ci_hi = c(rep(NA,na_head), ci_hi, rep(NA,na_tail))
  )

  if (T) {
    res$extras <- list(
      "Gamma_0.2" = Gamma_os_n(0.2),
      "Gamma_0.5" = Gamma_os_n(0.5),
      "Gamma_tilde_0.2" = Gamma_tilde_n(0.2),
      "Gamma_tilde_0.5" = Gamma_tilde_n(0.5),
      "r_tilde_0.2" = r_tilde_Mn(0.2),
      "r_tilde_0.5" = r_tilde_Mn(0.5),
      "eta_0.2" = eta_n(u=0.2,x=c(0,0)),
      "eta_0.5" = eta_n(u=0.5,x=c(1,1)),
      "g_n_0.2" = g_n(s=0.2,x=c(0,0)),
      "g_n_0.5" = g_n(s=0.5,x=c(1,1)),
      "omega_1" = omega_n(x=c(0,0),s=0.2,y=100,delta=0),
      "omega_2" = omega_n(x=c(0,0),s=0.4,y=150,delta=1),
      "omega_3" = omega_n(x=c(1,1),s=0.6,y=100,delta=0),
      "omega_4" = omega_n(x=c(1,1),s=0.8,y=150,delta=1),
      "p_n" = p_n
    )
  } # DEBUG (3 of 3): !!!!! figure out differences w vaccine package

  if (return_extras) {

    res$extras <- list(
      res$n <- n_orig,
      res$tau_ns <- c(rep(NA,na_head), tau_ns, rep(NA,na_tail))
      # !!!!! continue
    )

  }

  return(res)

}
