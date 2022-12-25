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
#' @param grid_size A list containing the following three keys:
#'     \itemize{
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
  dat, t_0, cve=T, cr=T, s_out=seq(from=min(dat$s), to=max(dat$s), l=101),
  ci_type="logit", cf_folds=1, edge_corr=F, params=list(),
  grid_size=list(y=101, s=101, x=5), verbose=F
) {

  if (class(dat_orig)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  # Alias variables
  dat_orig <- dat$df_vc # !!!!! Maybe change this later
  .v <- verbose
  # if (.v) { message("") } # !!!!!

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
    surv_type = "Super Learner",
    density_type = "binning",
    density_bins = 15,
    deriv_type = "linear",
    gamma_type = "Super Learner",
    q_n_type = "new"
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
  dat_orig <- round_dat(dat_orig, grid)

  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr=="min") { s_min2 <- min(dat_orig$s[dat_orig$s!=0], na.rm=T) }

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
  if (type=="Westling") {
    srvSL <- construct_Q_n(dat, vals, type=p$Q_n_type, print_coeffs=T)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
  } else if (type=="Wolock") {
    # !!!!! TO DO
  } else if (type=="Cox") {
    # !!!!! TO DO
  }

  omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
  f_sIx_n <- construct_f_sIx_n(dat, type, k=0, z1=F)





  # !!!!! CONTINUE
  {
    f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    dat2 <- ss(dat, which(dat$s!=0))
    Phi_n <- construct_Phi_n(dat2, type=p$ecdf_type)
    n_orig <- length(dat_orig$z)
    p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))
    eta_n <- construct_eta_n(dat, Q_n, p_n, vals=NA)
    r_tilde_Mn <- construct_r_tilde_Mn(dat_orig, vals=vlist$S_grid, Q_n)
    Gamma_tilde_n <- construct_Gamma_tilde_n(dat, r_tilde_Mn, p_n, vals=NA)
    f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
    q_n <- construct_q_n(type=p$q_n_type, dat, dat_orig, omega_n=omega_n, g_n=g_n,
                         p_n=p_n, r_tilde_Mn=r_tilde_Mn, Gamma_tilde_n=Gamma_tilde_n,
                         Q_n=Q_n, Qc_n=Qc_n, f_n_srv=f_n_srv)
    Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n,
                                       q_n, r_tilde_Mn, Gamma_tilde_n,
                                       vals=vlist$S_grid)

    # Construct one-step edge estimator
    if (p$edge_corr!="none") {
      g_sn <- construct_g_sn(dat, f_n_srv, g_n, p_n)
      r_Mn_edge_est <- r_Mn_edge(dat_orig, dat, g_sn, g_n, p_n, Q_n, omega_n)
      infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
                                                       r_Mn_edge_est, p_n)
      sigma2_edge_est <- (1/n_orig) * sum((
        infl_fn_r_Mn_edge(dat_orig$z, dat_orig$weights, dat_orig$s, dat_orig$x,
                          dat_orig$y, dat_orig$delta)
      )^2)
    }

    # Compute GCM (or least squares line) and extract its derivative
    grid <- sort(unique(dat$s))
    x_vals <- Phi_n(grid)
    indices_to_keep <- !base::duplicated(x_vals)
    x_vals <- x_vals[indices_to_keep]
    if (dir=="incr") {
      y_vals <- Gamma_os_n(grid[indices_to_keep])
    } else {
      y_vals <- -1 * Gamma_os_n(grid[indices_to_keep])
    }
    if (!any(x_vals==0)) {
      x_vals <- c(0, x_vals)
      y_vals <- c(0, y_vals)
    }
    if (p$convex_type=="GCM") {
      gcm <- gcmlcm(x=x_vals, y=y_vals, type="gcm")
      dGCM <- approxfun(
        x = gcm$x.knots[-length(gcm$x.knots)],
        y = gcm$slope.knots,
        method = "constant",
        rule = 2,
        f = 0
      )
    } else if (p$convex_type=="CLS") {
      gcm <- function(x) { 1 } # Ignored
      fit <- cvx.lse.reg(t=x_vals, z=y_vals)
      pred_x <- round(seq(0,1,0.001),3)
      pred_y <- predict(fit, newdata=pred_x)
      dGCM <- Vectorize(function(u) {
        width <- 0.05
        u1 <- u - width/2; u2 <- u + width/2;
        if (u1<0) { u2 <- u2 - u1; u1 <- 0; }
        if (u2>1) { u1 <- u1 - u2 + 1; u2 <- 1; }
        u1 <- round(u1,3); u2 <- round(u2,3);
        ind1 <- which(pred_x==u1); ind2 <- which(pred_x==u2);
        y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
        return((y2-y1)/width)
      })
    }

    # Construct Grenander-based r_Mn
    if (dir=="incr") {
      r_Mn_Gr <- Vectorize(function(u) { min(max(dGCM(Phi_n(u)),0),1) })
    } else {
      r_Mn_Gr <- Vectorize(function(u) { min(max(-1*dGCM(Phi_n(u)),0),1) })
    }

    # Compute variance component functions
    f_sIx_z1_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                    k=p$f_sIx_n_bins, z1=TRUE,
                                    s_scale=s_scale, s_shift=s_shift)
    f_s_z1_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_z1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, type=p$gamma_type,
                                 vals=vlist$S_grid, omega_n=omega_n,
                                 f_sIx_n=f_sIx_n, f_s_n=f_s_n,
                                 f_s_z1_n=f_s_z1_n)

    g_zn <- construct_g_zn(dat_orig, vals=NA, type="Super Learner", f_sIx_n,
                           f_sIx_z1_n)

    # Edge correction
    if (p$edge_corr=="none") {

      r_Mn <- r_Mn_Gr

    } else if (p$edge_corr=="point") {

      r_Mn <- function(u) {
        if(u==0) { r_Mn_edge_est } else { r_Mn_Gr(u) }
      }

    } else if (p$edge_corr=="min") {

      r_Mn <- Vectorize(function(u) {
        if(u==0 || u<s_min2) {
          r_Mn_edge_est
        } else {
          if (dir=="incr") {
            max(r_Mn_edge_est, r_Mn_Gr(u))
          } else {
            min(r_Mn_edge_est, r_Mn_Gr(u))
          }
        }
      })

      # gren_s_out <- sapply(c(1:length(s_out)), function(i) {
      #   if (dir=="incr") {
      #     as.numeric(r_Mn(s_out[i])>r_Mn_edge_est)
      #   } else {
      #     as.numeric(r_Mn(s_out[i])<r_Mn_edge_est)
      #   }
      # })

    }

    # Generate estimates for each point
    ests <- r_Mn(s_out)

    # Construct variance scale factor
    deriv_r_Mn <- construct_deriv_r_Mn(r_Mn, type=p$deriv_type, dir=dir)

    tau_n <- construct_tau_n(deriv_r_Mn, gamma_n, f_sIx_n, f_s_n, g_zn,
                             dat_orig)

    # Generate confidence limits
    if (p$ci_type=="none") {

      ci_lo <- rep(0, length(ests))
      ci_hi <- rep(0, length(ests))

    } else {

      # Generate variance scale factor for each point
      tau_ns <- tau_n(s_out)

      # Construct CIs
      # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
      # The Normal approximation would use qnorm(0.975, sd=0.52) instead
      qnt <- 1.00
      n_orig <- length(dat_orig$z)
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
      } else if (p$ci_type=="trunc") {
        ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
        ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
        ci_lo %<>% pmax(0) %>% pmin(1)
        ci_hi %<>% pmax(0) %>% pmin(1)
      }

      # Edge correction
      if (p$edge_corr=="point") {
        ci_lo[1] <- ests[1] - 1.96*sqrt(sigma2_edge_est/n_orig)
        ci_hi[1] <- ests[1] + 1.96*sqrt(sigma2_edge_est/n_orig)
      } else if (p$edge_corr=="min") {
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

    res <- list(
      s = s_out_orig,
      est = c(rep(NA,na_head), ests, rep(NA,na_tail)),
      ci_lo = c(rep(NA,na_head), ci_lo, rep(NA,na_tail)),
      ci_hi = c(rep(NA,na_head), ci_hi, rep(NA,na_tail))
    )

    res$tau_ns <- c(rep(NA,na_head), tau_ns, rep(NA,na_tail))
    res$n <- n_orig

  }

  return(999)

}
