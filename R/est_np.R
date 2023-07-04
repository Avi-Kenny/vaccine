#' Estimate CVE/CR nonparametrically
#'
#' @description See docs for est_ce and params_ce_cox
#' @noRd
est_np <- function(
    dat, t_0, cr=T, cve=F,
    s_out=seq(from=min(dat$v$s, na.rm=T), to=max(dat$v$s, na.rm=T), l=101),
    ci_type="transformed", placebo_risk_method="KM", return_extras=F,
    dir="decr", edge_corr=F, params=list(), grid_size=list(y=101, s=101, x=5),
    cf_folds=1
) {

  if (!methods::is(dat,"vaccine_dat")) {
    stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by lo",
                "ad_data()."))
  }

  if (!(attr(dat, "groups") %in% c("vaccine", "both"))) {
    stop("Vaccine group data not detected.")
  }

  if (!(dir %in% c("decr", "incr"))) {
    stop("`dir` must equal one of c('decr','incr').")
  } else {
    dir_factor <- ifelse(dir=="decr", -1, 1)
  }

  # !!!!! Validate other inputs; import error handling function from SimEngine

  # Alias variables
  dat_copy <- dat # Consider keeping this and changing old `dat` variable
  dat_orig <- dat$v

  if (any(is.na(s_out))) { stop("NA values not allowed in s_out.") }

  # Set params
  .default_params <- list(
    surv_type = "Super Learner",
    density_type = "binning",
    density_bins = 15,
    deriv_type = "m-spline",
    gamma_type = "Super Learner",
    q_n_type = "zero", # standard
    convex_type = "GCM",
    mono_cis = T
  )
  for (i in c(1:length(.default_params))) {
    p_name <- names(.default_params)[i]
    if (is.null(params[[p_name]])) { params[[p_name]] <- .default_params[[i]] }
  }
  p <- params
  p$ci_type <- ci_type
  p$cf_folds <- cf_folds
  p$edge_corr <- edge_corr

  # Rescale S to lie in [0,1] and create rounded data object
  s_min <- min(dat_orig$s, na.rm=T)
  s_max <- max(dat_orig$s, na.rm=T)
  s_shift <- -1 * s_min
  s_scale <- 1/(s_max-s_min)
  dat_orig$s <- (dat_orig$s+s_shift)*s_scale
  grid <- create_grid(dat_orig, grid_size, t_0)
  dat_orig_rounded <- round_dat(dat_orig, grid, grid_size)

  # Rescale/round s_out and remove s_out points outside [0,1]
  s_out_orig <- s_out
  s_out <- (s_out+s_shift)*s_scale
  na_head <- sum(s_out<0)
  na_tail <- sum(s_out>1)
  if (na_head>0) { s_out <- s_out[-c(1:na_head)] }
  len_p <- length(s_out)
  if (na_tail>0) { s_out <- s_out[-c((len_p-na_tail+1):len_p)] }
  s_out <- sapply(s_out, function(s) { grid$s[which.min(abs(grid$s-s))] })

  # Prepare precomputation values for conditional survival estimator
  x_distinct <- dplyr::distinct(dat_orig_rounded$x)
  x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
  vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
  vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
  vals <- list(
    t = vals_pre$t,
    x = subset(vals_pre, select=names(dat_orig_rounded$x)),
    s = vals_pre$s
  )

  # Create phase-two data object (unrounded)
  dat <- ss(dat_orig_rounded, which(dat_orig_rounded$z==1)) # !!!!! uses the `dat` name again; change???

  # Fit conditional survival estimator
  srvSL <- construct_Q_n(p$surv_type, dat, vals)
  Q_n <- srvSL$srv
  Qc_n <- srvSL$cens

  # Use rounded data objects moving forward
  dat_orig <- dat_orig_rounded

  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr) { s_min2 <- min(dat_orig$s[dat_orig$s!=0], na.rm=T) }

  # Compute various nuisance functions
  omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
  f_sIx_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins, z1=F)
  f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
  g_n <- construct_g_n(f_sIx_n, f_s_n)
  Phi_n <- construct_Phi_n(dat_orig, dat)
  r_tilde_Mn <- construct_r_tilde_Mn(dat_orig, Q_n, t_0)
  f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)
  q_n <- construct_q_n(type=p$q_n_type, dat, omega_n, g_n, r_tilde_Mn,
                       f_n_srv)

  Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n,
                                     q_n, r_tilde_Mn)

  # Compute edge-corrected estimator and standard error
  n_orig <- attr(dat_orig, "n_orig")
  if (p$edge_corr) {
    p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))
    g_sn <- construct_g_sn(dat, f_n_srv, g_n, p_n)
    r_Mn_edge_est <- r_Mn_edge(dat_orig, g_sn, g_n, p_n, Q_n, omega_n, t_0)
    r_Mn_edge_est <- min(max(r_Mn_edge_est, 0), 1)
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
  gcm_y_vals <- dir_factor *
    sapply(sort(unique(dat$s))[indices_to_keep], Gamma_os_n)
  if (!any(gcm_x_vals==0)) {
    gcm_x_vals <- c(0, gcm_x_vals)
    gcm_y_vals <- c(0, gcm_y_vals)
  }
  if (p$convex_type=="GCM") {
    gcm <- fdrtool::gcmlcm(x=gcm_x_vals, y=gcm_y_vals, type="gcm")
    dGCM <- stats::approxfun(
      x = gcm$x.knots[-length(gcm$x.knots)],
      y = gcm$slope.knots,
      method = "constant",
      rule = 2,
      f = 0
    )
  } else if (p$convex_type=="CLS") {
    # !!!!! This section is experimental; theory not yet developed
    gcm <- function(x) { 1 } # Ignored
    fit <- simest::cvx.lse.reg(t=gcm_x_vals, z=gcm_y_vals)
    pred_x <- round(seq(0,1,0.001),3)
    pred_y <- stats::predict(fit, newdata=pred_x)
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
  r_Mn_Gr <- function(u) { min(max(dir_factor*dGCM(Phi_n(u)), 0), 1) }

  # Compute variance component nuisance estimators
  f_sIx_z1_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins,
                                  z1=T)
  f_s_z1_n <- construct_f_s_n(dat_orig, f_sIx_z1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, type="Super Learner", omega_n,
                               grid)
  g_zn <- construct_g_zn(dat_orig, type="Super Learner", f_sIx_n, f_sIx_z1_n)

  # Create edge-corrected r_Mn estimator
  if (p$edge_corr) {
    r_Mn <- function(u) {
      if(u==0 || u<s_min2) {
        return(r_Mn_edge_est)
      } else {
        if (dir=="decr") {
          return(min(r_Mn_edge_est, r_Mn_Gr(u)))
        } else {
          return(max(r_Mn_edge_est, r_Mn_Gr(u)))
        }
      }
    }
  } else {
    r_Mn <- r_Mn_Gr
  }

  # Generate estimates for each point
  ests_cr <- sapply(s_out, r_Mn)

  # Construct variance scale factor
  deriv_r_Mn <- construct_deriv_r_Mn(type=p$deriv_type, r_Mn, dir, grid)
  tau_n <- construct_tau_n(deriv_r_Mn, gamma_n, f_sIx_n, g_zn, dat_orig)

  # Generate confidence limits
  if (p$ci_type=="none") {

    ci_lo_cr <- rep(NA, length(ests_cr))
    ci_hi_cr <- rep(NA, length(ests_cr))

  } else {

    # Generate variance scale factor for each point
    tau_ns <- sapply(s_out, tau_n)

    # Construct CIs
    # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
    qnt <- 1.00
    if (p$ci_type=="regular") {
      ci_lo_cr <- ests_cr - (qnt*tau_ns)/(n_orig^(1/3))
      ci_hi_cr <- ests_cr + (qnt*tau_ns)/(n_orig^(1/3))
    } else if (p$ci_type=="transformed") {
      ci_lo_cr <- expit(
        logit(ests_cr) - (qnt*tau_ns*deriv_logit(ests_cr))/(n_orig^(1/3))
      )
      ci_hi_cr <- expit(
        logit(ests_cr) + (qnt*tau_ns*deriv_logit(ests_cr))/(n_orig^(1/3))
      )
    }

    # CI edge correction
    if (p$edge_corr) {
      ci_lo_cr2 <- ests_cr[1] - 1.96*sqrt(sigma2_edge_est/n_orig)
      ci_hi_cr2 <- ests_cr[1] + 1.96*sqrt(sigma2_edge_est/n_orig)
      if (dir=="decr") {
        ci_lo_cr <- In(r_Mn_edge_est<=ests_cr)*pmin(ci_lo_cr,ci_lo_cr2) +
          In(r_Mn_edge_est>ests_cr)*ci_lo_cr
        ci_lo_cr[1] <- ci_lo_cr2
        ci_hi_cr <- In(r_Mn_edge_est<=ests_cr)*pmin(ci_hi_cr,ci_hi_cr2) +
          In(r_Mn_edge_est>ests_cr)*ci_hi_cr
        ci_hi_cr[1] <- ci_hi_cr2
      } else {
        ci_lo_cr <- In(r_Mn_edge_est>=ests_cr)*pmax(ci_lo_cr,ci_lo_cr2) +
          In(r_Mn_edge_est<ests_cr)*ci_lo_cr
        ci_lo_cr[1] <- ci_lo_cr2
        ci_hi_cr <- In(r_Mn_edge_est>=ests_cr)*pmax(ci_hi_cr,ci_hi_cr2) +
          In(r_Mn_edge_est<ests_cr)*ci_hi_cr
        ci_hi_cr[1] <- ci_hi_cr2
      }
    }

  }

  # Monotone CI correction
  if (p$mono_cis) {
    val <- ci_lo_cr[1]
    for (i in c(2:length(ci_lo_cr))) {
      if (dir=="decr") {
        if (!is.na(ci_lo_cr[i]) && !is.na(val) && ci_lo_cr[i]>val) {
          ci_lo_cr[i] <- val
        }
      } else {
        if (!is.na(ci_lo_cr[i]) && !is.na(val) && ci_lo_cr[i]<val) {
          ci_lo_cr[i] <- val
        }
      }
      val <- ci_lo_cr[i]
    }
    val <- ci_hi_cr[1]
    for (i in c(2:length(ci_hi_cr))) {
      if (dir=="decr") {
        if (!is.na(ci_hi_cr[i]) && !is.na(val) && ci_hi_cr[i]>val) {
          ci_hi_cr[i] <- val
        }
      } else {
        if (!is.na(ci_hi_cr[i]) && !is.na(val) && ci_hi_cr[i]<val) {
          ci_hi_cr[i] <- val
        }
      }
      val <- ci_hi_cr[i]
    }
  }

  # Create results object
  res <- list(
    "cr" = list(
      s = s_out_orig,
      est = c(rep(NA,na_head), ests_cr, rep(NA,na_tail)),
      ci_lo = c(rep(NA,na_head), ci_lo_cr, rep(NA,na_tail)),
      ci_hi = c(rep(NA,na_head), ci_hi_cr, rep(NA,na_tail))
    ),
    "cve" = list(s=s_out)
  )

  # Compute CVE
  if (cve) {
    if (attr(dat_copy, "groups")!="both") {
      stop("Placebo group data not detected.")
    }
    ov <- overall(dat=dat_copy, t_0=t_0, method=placebo_risk_method, ve=F)
    risk_p <- ov[ov$group=="placebo","est"]
    se_p <- ov[ov$group=="placebo","se"]
    res$cve$est <- 1 - res$cr$est/risk_p
    res$cve$ci_lo <- 1 - res$cr$ci_hi/risk_p
    res$cve$ci_hi <- 1 - res$cr$ci_lo/risk_p
  }

  # !!!!! TO DO: perform finite sample variance correction for CVE

  if (return_extras) {

    res$extras <- list(
      res$n <- n_orig,
      res$tau_ns <- c(rep(NA,na_head), tau_ns, rep(NA,na_tail))
      # !!!!! continue
    )

  }

  return(res)

}
