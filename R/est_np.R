#' Estimate CVE/CR nonparametrically
#'
#' @description See docs for est_ce and params_ce_cox
#' @noRd
est_np <- function(
    dat, t_0, cr, cve, s_out, ci_type, placebo_risk_method, return_p_value,
    return_extras, params, cf_folds, p_val_only=FALSE
) {

  # Set params
  .default_params <- params_ce_np()
  .default_params$gamma_type <- "Super Learner" # !!!!! Move to params or hard-code
  .default_params$q_n_type <- "zero" # !!!!! Temp; implement "standard"
  # .default_params$q_n_type <- "standard" # !!!!! Temp; implement "standard"
  .default_params$mono_cis <- T # !!!!! Move to params or hard-code
  for (i in c(1:length(.default_params))) {
    p_name <- names(.default_params)[i]
    if (is.null(params[[p_name]])) { params[[p_name]] <- .default_params[[i]] }
  }
  p <- params

  if (!(p$dir %in% c("decr", "incr"))) {
    stop("`dir` must equal one of c('decr','incr').")
  } else {
    dir_factor <- ifelse(p$dir=="decr", -1, 1)
  }

  # Create filtered data objects, alias variables
  dat_v <- dat[dat$a==1,]
  dim_x <- attr(dat, "dim_x")
  n_vacc <- attr(dat, "n_vacc")

  if (any(is.na(s_out))) { stop("NA values not allowed in s_out.") }

  # Rescale S to lie in [0,1] and create grid for rounding
  s_min <- min(dat_v$s, na.rm=T)
  s_max <- max(dat_v$s, na.rm=T)
  s_shift <- -1 * s_min
  s_scale <- 1/(s_max-s_min)
  dat_v$s <- (dat_v$s+s_shift)*s_scale
  grid <- create_grid(dat_v, p$grid_size, t_0) # !!!!! feed in dat instead ?????

  # Create additional filtered datasets
  dat_v_rd <- round_dat(dat_v, grid, p$grid_size) # !!!!! see (1) above; also, make grid_size a param, like in est_med
  dat_v2_rd <- dat_v_rd[dat_v_rd$z==1,]
  datx_v_rd <- dat_v_rd[, c(1:dim_x), drop=F]
  class(datx_v_rd) <- "data.frame"

  # Rescale/round s_out and remove s_out points outside [0,1]
  s_out_orig <- s_out
  s_out <- (s_out+s_shift)*s_scale
  na_head <- sum(s_out<0)
  na_tail <- sum(s_out>1)
  if (na_head>0) { s_out <- s_out[-c(1:na_head)] }
  len_p <- length(s_out)
  if (na_tail>0) { s_out <- s_out[-c((len_p-na_tail+1):len_p)] }
  s_out <- sapply(s_out, function(s) { grid$s[which.min(abs(grid$s-s))] })

  # Precomputation values for conditional survival/censoring estimators
  if (attr(dat, "covariates_ph2")) {
    datx_v2_rd <- dat_v2_rd[, c(1:dim_x), drop=F]
    class(datx_v2_rd) <- "data.frame"
    x_distinct <- dplyr::distinct(datx_v2_rd)
  } else {
    x_distinct <- dplyr::distinct(datx_v_rd)
  }
  x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
  vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
  vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
  vals <- list(
    t = vals_pre$t,
    x = subset(vals_pre, select=names(datx_v_rd)),
    s = vals_pre$s
  )

  # Fit conditional survival estimator
  srvSL <- construct_Q_n(p$surv_type, dat_v2_rd, vals)
  Q_n <- srvSL$srv
  Qc_n <- srvSL$cens

  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr) { s_min2 <- min(dat_v_rd$s[dat_v_rd$s!=0], na.rm=T) }

  # Compute various nuisance functions
  omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
  f_sIx_n <- construct_f_sIx_n(dat_v2_rd, type=p$density_type, k=p$density_bins,
                               z1=F)
  f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)
  if (!p_val_only) {
    f_s_n <- construct_f_s_n(dat_v_rd, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    Phi_n <- construct_Phi_n(dat_v2_rd)
    r_tilde_Mn <- construct_r_tilde_Mn(dat_v_rd, Q_n, t_0)
    q_n <- construct_q_n(type=p$q_n_type, dat_v2_rd, omega_n, g_n, r_tilde_Mn,
                         f_n_srv)
    Gamma_os_n <- construct_Gamma_os_n(dat_v_rd, omega_n, g_n, q_n, r_tilde_Mn)
  }

  if (return_p_value) {

    # Helper function to compute P values
    compute_p_val <- function(alt_type, beta_n, var_n) {
      if (alt_type=="incr") {
        p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n), lower.tail=FALSE)
      } else if (alt_type=="decr") {
        p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n))
      } else if (alt_type=="two-tailed") {
        p_val <- pchisq(beta_n^2/var_n, df=1, lower.tail=FALSE)
      }
      return(p_val)
    }

    # Construct additional nuisance functions
    chk(1) # !!!!!
    etastar_n <- construct_etastar_n(Q_n, t_0, vals)
    chk(2) # !!!!!
    q_tilde_n <- construct_q_tilde_n(type=p$q_n_type, f_n_srv, f_sIx_n,
                                     omega_n)
    chk(3) # !!!!!
    Theta_os_n <- construct_Theta_os_n(dat_v_rd, omega_n, f_sIx_n, q_tilde_n,
                                       etastar_n)
    chk(4) # !!!!!
    infl_fn_Theta <- construct_infl_fn_Theta(omega_n, f_sIx_n, q_tilde_n,
                                             etastar_n, Theta_os_n)
    chk(5) # !!!!!
    infl_fn_beta_n <- construct_infl_fn_beta_n(infl_fn_Theta)
    chk(6) # !!!!!

    # Functions needed for edge-corrected test
    # p_n <- (1/n_vacc) * sum(dat$weights * In(dat$s!=0))
    # f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
    # g_n <- construct_g_n(f_sIx_n, f_s_n)
    # g_sn <- construct_g_sn(dat, f_n_srv, g_n, p_n)
    # r_Mn_edge_est <- r_Mn_edge(dat_orig, dat, g_sn, g_n, p_n, Q_n, omega_n)
    # infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
    #                                                  r_Mn_edge_est, p_n)
    # infl_fn_beta_en <- construct_infl_fn_beta_en(infl_fn_Theta,
    #                                              infl_fn_r_Mn_edge)

    if (p$edge_corr) {

      stop("Edge-corrected test not yet implemented.")

      # # Construct pieces needed for beta_n
      # u_mc <- round(seq(0.001,1,0.001),3)
      # m <- length(u_mc)
      # lambda_1 <- mean(u_mc) # ~1/2
      # lambda_2 <- mean((u_mc)^2) # ~1/3
      # lambda_3 <- mean((u_mc)^3) # ~1/4
      # beta_n <- mean((
      #   (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
      #     (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
      # ) * Theta_os_n(u_mc))
      # infl_fn_beta_n <- construct_infl_fn_beta_n(infl_fn_Theta)
      #
      # # Construct pieces needed for beta_en
      # beta_en <- Theta_os_n(1) - r_Mn_edge_est
      #
      # # Calculate variance components
      # sigma2_bn <- 0
      # sigma2_ben <- 0
      # cov_n <- 0
      # for (i in c(1:n_vacc)) {
      #   z_i <- dat_orig$z[i]
      #   s_i <- dat_orig$s[i]
      #   x_i <- as.numeric(dat_orig$x[i,])
      #   y_i <- dat_orig$y[i]
      #   delta_i <- dat_orig$delta[i]
      #   weight_i <- dat_orig$weight[i]
      #
      #   if_bn <- infl_fn_beta_n(s_i, y_i, delta_i, weight_i, x_i)
      #   if_ben <- infl_fn_beta_en(z_i, x_i, y_i, delta_i, s_i, weight_i)
      #
      #   sigma2_bn <- sigma2_bn + if_bn^2
      #   sigma2_ben <- sigma2_ben + if_ben^2
      #   cov_n <- cov_n + if_bn*if_ben
      # }
      # sigma_bn <- sqrt(sigma2_bn/n_vacc)
      # sigma_ben <- sqrt(sigma2_ben/n_vacc)
      # cov_n <- cov_n/n_vacc
      # rho_n <- cov_n/(sigma_bn*sigma_ben)
      #
      # # Calculate combined test statistic
      # beta_star_n <- sqrt(n_vacc/(2+2*rho_n)) *
      #   (beta_n/sigma_bn + beta_en/sigma_ben)
      #
      # res[[length(res)+1]] <- list(
      #   type = "combined 2",
      #   p_val = compute_p_val(alt_type, beta_star_n, 1),
      #   beta_n = beta_star_n,
      #   var_n = 1
      # )

    } else {

      # Compute test statistic and variance estimate
      chk(7) # !!!!!
      u_mc <- round(seq(0.01,1,0.01),2)
      m <- length(u_mc)
      lambda_1 <- mean(u_mc) # ~1/2
      lambda_2 <- mean((u_mc)^2) # ~1/3
      lambda_3 <- mean((u_mc)^3) # ~1/4
      chk(8) # !!!!!

      beta_n <- mean((
        (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
          (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
      ) * Theta_os_n(u_mc))
      chk(9) # !!!!!

      # var_n <- 0
      # for (i in c(1:n_vacc)) {
      #   var_n <- var_n + (
      #     infl_fn_beta_n(dat_orig$s[i], dat_orig$y[i], dat_orig$delta[i],
      #                    dat_orig$weight[i], as.numeric(dat_orig$x[i,]))
      #   )^2
      # }
      # var_n <- var_n/n_vacc^2
      var_n <- (1/n_vacc^2) * sum(apply(dat_v_rd, 1, function(r) {
        (infl_fn_beta_n(r[["s"]], r[["y"]], r[["delta"]], r[["weights"]],
                        as.numeric(r[1:dim_x])))^2
      }))
      chk(10) # !!!!!

      test_res <- list(
        p_val = compute_p_val(alt_type=p$dir, beta_n, var_n),
        beta_n = beta_n,
        var_n = var_n
      )
      chk(11) # !!!!!
      chk(12.1, paste("P-val:", test_res$p_val)) # !!!!!
      chk(12.2, paste("beta_n:", test_res$beta_n)) # !!!!!
      chk(12.3, paste("var_n:", test_res$var_n)) # !!!!!
      if (p_val_only) { return(list(p=test_res$p_val)) }

    }

    if (F) {
      res$extras <- list(
        Theta_1.0 = Theta_os_n(1),
        r_Mn_0.0 = r_Mn_edge_est,
        # var_edge = var_edge,
        # sd_edge = sqrt(var_edge),
        beta_n = beta_n,
        beta_en = beta_en,
        sigma_bn = sigma_bn,
        sigma_ben = sigma_ben,
        rho_n = rho_n
      )

      # res$extras <- list(
      #   Theta_0.1 = Theta_os_n(0.1),
      #   Theta_0.4 = Theta_os_n(0.4),
      #   Theta_0.8 = Theta_os_n(0.8),
      #   etastar_0.1 = mean(sapply(c(1:n_vacc), function(i) {
      #     etastar_n(u=0.1, x=as.numeric(dat_orig$x[i,]))
      #   })),
      #   etastar_0.4 = mean(sapply(c(1:n_vacc), function(i) {
      #     etastar_n(u=0.4, x=as.numeric(dat_orig$x[i,]))
      #   })),
      #   etastar_0.8 = mean(sapply(c(1:n_vacc), function(i) {
      #     etastar_n(u=0.8, x=as.numeric(dat_orig$x[i,]))
      #   }))
      # )
    } # DEBUG: return debugging components

  }

  # Compute edge-corrected estimator and standard error
  if (p$edge_corr) {
    if (attr(dat, "covariates_ph2")) {
      # stop("edge correction not yet available for covariates_ph2==T.")
      warning("edge correction not yet available for covariates_ph2==T.") # !!!!!
    }
    p_n <- (1/n_vacc) * sum(dat_v2_rd$weights * In(dat_v2_rd$s!=0))
    g_sn <- construct_g_sn(dat_v2_rd, f_n_srv, g_n, p_n)
    r_Mn_edge_est <- r_Mn_edge(dat_v_rd, g_sn, g_n, p_n, Q_n, omega_n, t_0)
    r_Mn_edge_est <- min(max(r_Mn_edge_est, 0), 1)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
                                                     r_Mn_edge_est, p_n, t_0)
    sigma2_edge_est <- mean(apply(dat_v_rd, 1, function(r) {
      (infl_fn_r_Mn_edge(r[["z"]], r[["weights"]], r[["s"]],
                         as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]]))^2
    }))
  }

  # Compute GCM (or least squares line) and extract its derivative
  gcm_x_vals <- sapply(sort(unique(dat_v2_rd$s)), Phi_n)
  indices_to_keep <- !base::duplicated(gcm_x_vals)
  gcm_x_vals <- gcm_x_vals[indices_to_keep]
  gcm_y_vals <- dir_factor *
    sapply(sort(unique(dat_v2_rd$s))[indices_to_keep], Gamma_os_n)

  if (any(is.nan(gcm_y_vals)) || any(is.na(gcm_y_vals))) {
    stop("Gamma_os_n produced NAN or NA values.")
  }

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
  if (ci_type!="none") {
    f_sIx_z1_n <- construct_f_sIx_n(dat_v2_rd, type=p$density_type,
                                    k=p$density_bins, z1=T)
    f_s_z1_n <- construct_f_s_n(dat_v_rd, f_sIx_z1_n)
    gamma_n <- construct_gamma_n(dat_v_rd, type="Super Learner", omega_n, grid)
    g_zn <- construct_g_zn(dat_v_rd, type="Super Learner", f_sIx_n, f_sIx_z1_n)
  }

  # Create edge-corrected r_Mn estimator
  if (p$edge_corr) {
    r_Mn <- function(u) {
      if(u==0 || u<s_min2) {
        return(r_Mn_edge_est)
      } else {
        if (p$dir=="decr") {
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

  # Generate confidence limits
  if (ci_type=="none") {

    ci_lo_cr <- rep(NA, length(ests_cr))
    ci_up_cr <- rep(NA, length(ests_cr))

  } else {

    # Construct variance scale factor function
    deriv_r_Mn <- construct_deriv_r_Mn(type=p$deriv_type, r_Mn, p$dir, grid)
    tau_n <- construct_tau_n(dat_v_rd, deriv_r_Mn, gamma_n, f_sIx_n, g_zn)

    # Generate variance scale factor for each point
    tau_ns <- sapply(s_out, tau_n)

    # Construct CIs
    # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
    qnt <- 1.00
    if (ci_type=="regular") {
      ci_lo_cr <- ests_cr - (qnt*tau_ns)/(n_vacc^(1/3))
      ci_up_cr <- ests_cr + (qnt*tau_ns)/(n_vacc^(1/3))
    } else if (ci_type=="truncated") {
      ci_lo_cr <- pmax(ests_cr - (qnt*tau_ns)/(n_vacc^(1/3)), 0)
      ci_up_cr <- pmin(ests_cr + (qnt*tau_ns)/(n_vacc^(1/3)), 1)
    } else if (ci_type=="transformed") {
      ci_lo_cr <- expit(
        logit(ests_cr) - qnt*deriv_logit(ests_cr)*(tau_ns/(n_vacc^(1/3)))
      )
      ci_up_cr <- expit(
        logit(ests_cr) + qnt*deriv_logit(ests_cr)*(tau_ns/(n_vacc^(1/3)))
      )
    } else if (ci_type=="transformed 2") {
      ci_lo_cr <- expit2(
        logit2(ests_cr) - qnt*deriv_logit2(ests_cr)*(tau_ns/(n_vacc^(1/3)))
      )
      ci_up_cr <- expit2(
        logit2(ests_cr) + qnt*deriv_logit2(ests_cr)*(tau_ns/(n_vacc^(1/3)))
      )
    }

    # CI edge correction
    if (p$edge_corr) {

      se_edge_est <- sqrt(sigma2_edge_est/n_vacc)
      if (ci_type=="regular") {
        ci_lo_cr2 <- ests_cr[1] - 1.96*se_edge_est
        ci_up_cr2 <- ests_cr[1] + 1.96*se_edge_est
      } else if (ci_type=="truncated") {
        ci_lo_cr2 <- max(ests_cr[1] - 1.96*se_edge_est, 0)
        ci_up_cr2 <- min(ests_cr[1] - 1.96*se_edge_est, 1)
      } else if (ci_type=="transformed") {
        ci_lo_cr2 <- expit(
          logit(ests_cr[1]) - 1.96*deriv_logit(ests_cr[1])*se_edge_est
        )
        ci_up_cr2 <- expit(
          logit(ests_cr[1]) + 1.96*deriv_logit(ests_cr[1])*se_edge_est
        )
      } else if (ci_type=="transformed 2") {
        ci_lo_cr2 <- expit2(
          logit2(ests_cr[1]) - 1.96*deriv_logit2(ests_cr[1])*se_edge_est
        )
        ci_up_cr2 <- expit2(
          logit2(ests_cr[1]) + 1.96*deriv_logit2(ests_cr[1])*se_edge_est
        )
      }

      if (p$dir=="decr") {
        ind_edge <- In(r_Mn_edge_est<=ests_cr)
        ci_lo_cr <- ind_edge*pmin(ci_lo_cr,ci_lo_cr2) + (1-ind_edge)*ci_lo_cr
        ci_up_cr <- ind_edge*pmin(ci_up_cr,ci_up_cr2) + (1-ind_edge)*ci_up_cr
        ci_lo_cr[1] <- ci_lo_cr2 # !!!!!
        ci_up_cr[1] <- ci_up_cr2 # !!!!!
      } else {
        ind_edge <- In(r_Mn_edge_est>=ests_cr)
        ci_lo_cr <- ind_edge*pmax(ci_lo_cr,ci_lo_cr2) + (1-ind_edge)*ci_lo_cr
        ci_up_cr <- ind_edge*pmax(ci_up_cr,ci_up_cr2) + (1-ind_edge)*ci_up_cr
        ci_lo_cr[1] <- ci_lo_cr2 # !!!!!
        ci_up_cr[1] <- ci_up_cr2 # !!!!!
      }

    }

    # Monotone CI correction
    if (p$mono_cis) {
      new_lims <- monotonize_cis(
        ci_lo = ci_lo_cr,
        ci_up = ci_up_cr,
        dir = p$dir,
        type = "regular"
      )
      ci_lo_cr <- new_lims$ci_lo
      ci_up_cr <- new_lims$ci_up
    }

  }

  # Create results object
  res <- list()
  if (cr) {
    res$cr <- list(
      s = s_out_orig,
      est = c(rep(NA,na_head), ests_cr, rep(NA,na_tail)),
      ci_lower = c(rep(NA,na_head), ci_lo_cr, rep(NA,na_tail)),
      ci_upper = c(rep(NA,na_head), ci_up_cr, rep(NA,na_tail))
    )
  }

  # Compute CVE
  if (cve) {

    res$cve <- list(s=s_out_orig)
    if (attr(dat_v, "groups")!="both") { stop("Placebo group not detected.") }
    ov <- est_overall(dat=dat, t_0=t_0, method=placebo_risk_method, ve=F)
    risk_p <- ov[ov$group=="placebo","est"]
    se_p <- ov[ov$group=="placebo","se"] # This equals sd_p/n_vacc
    res$cve$est <- 1 - res$cr$est/risk_p

    if (ci_type=="none") {

      na_vec <- rep(NA, length(res$cve$s_out))
      res$cve$se <- na_vec
      res$cve$ci_lower <- na_vec
      res$cve$ci_upper <- na_vec

    } else {

      # Variance of Chernoff RV
      var_z <- 0.52^2

      # !!!!! Need to do edge_corr

      # Finite sample corrected SE estimate
      res$cve$se <- sqrt(
        ( res$cr$est^2/risk_p^4 ) * se_p^2 +
          (qnt*tau_ns/(risk_p*n_vacc^(1/3)))^2 * var_z
      )

      if (ci_type=="regular") {
        res$cve$ci_lower <- res$cve$est - 1.96*res$cve$se
        res$cve$ci_upper <- res$cve$est + 1.96*res$cve$se
      } else if (ci_type=="truncated") {
        res$cve$ci_lower <- pmin(res$cve$est - 1.96*res$cve$se, 1)
        res$cve$ci_upper <- pmin(res$cve$est + 1.96*res$cve$se, 1)
      } else if (ci_type=="transformed") {
        res$cve$ci_lower <- 1 - exp(
          log(1-res$cve$est) + 1.96*(1/(1-res$cve$est))*res$cve$se
        )
        res$cve$ci_upper <- 1 - exp(
          log(1-res$cve$est) - 1.96*(1/(1-res$cve$est))*res$cve$se
        )
      } else if (ci_type=="transformed 2") {
        res$cve$ci_lower <- 1 - exp2(
          log2(1-res$cve$est) + 1.96*deriv_log2(1-res$cve$est)*res$cve$se
        )
        res$cve$ci_upper <- 1 - exp2(
          log2(1-res$cve$est) - 1.96*deriv_log2(1-res$cve$est)*res$cve$se
        )
      }

      # CI edge correction
      if (p$edge_corr) {

        # Finite sample corrected SE estimate
        se_edge_est_cve <- sqrt(
          (ests_cr[1]^2/risk_p^4)*se_p^2 + (sigma2_edge_est/n_vacc)/(risk_p^2)
        )

        if (ci_type=="regular") {
          ci_lo_cve2 <- res$cve$est[1] - 1.96*se_edge_est_cve
          ci_up_cve2 <- res$cve$est[1] + 1.96*se_edge_est_cve
        } else if (ci_type=="truncated") {
          ci_lo_cve2 <- min(res$cve$est[1] - 1.96*se_edge_est_cve, 1)
          ci_up_cve2 <- min(res$cve$est[1] - 1.96*se_edge_est_cve, 1)
        } else if (ci_type=="transformed") {
          ci_lo_cve2 <- 1 - exp(
            log(1-res$cve$est[1]) + 1.96*(1/(1-res$cve$est[1]))*se_edge_est_cve
          )
          ci_up_cve2 <- 1 - exp(
            log(1-res$cve$est[1]) - 1.96*(1/(1-res$cve$est[1]))*se_edge_est_cve
          )
        } else if (ci_type=="transformed 2") {
          ci_lo_cve2 <- 1 - exp2(
            log2(1-res$cve$est[1])
            + 1.96*deriv_log2(1-res$cve$est[1])*se_edge_est_cve
          )
          ci_up_cve2 <- 1 - exp2(
            log2(1-res$cve$est[1])
            - 1.96*deriv_log2(1-res$cve$est[1])*se_edge_est_cve
          )
        }

        if (p$dir=="decr") {
          res$cve$ci_lower <- ind_edge*pmax(res$cve$ci_lower,ci_lo_cve2) +
            (1-ind_edge)*res$cve$ci_lower
          res$cve$ci_upper <- ind_edge*pmax(res$cve$ci_upper,ci_up_cve2) +
            (1-ind_edge)*res$cve$ci_upper
          res$cve$ci_lower[1] <- ci_lo_cve2 # !!!!!
          res$cve$ci_upper[1] <- ci_up_cve2 # !!!!!
        } else {
          res$cve$ci_lower <- ind_edge*pmin(res$cve$ci_lower,ci_lo_cve2) +
            (1-ind_edge)*res$cve$ci_lower
          res$cve$ci_upper <- ind_edge*pmin(res$cve$ci_upper,ci_up_cve2) +
            (1-ind_edge)*res$cve$ci_upper
          res$cve$ci_lower[1] <- ci_lo_cve2 # !!!!!
          res$cve$ci_upper[1] <- ci_up_cve2 # !!!!!
        }

      }

      if (p$mono_cis) {
        dir_opposite <- ifelse(p$dir=="decr", "incr", "decr")
        new_lims <- monotonize_cis(
          ci_lo = res$cve$ci_lower,
          ci_up = res$cve$ci_upper,
          dir = dir_opposite,
          type = "regular"
        )
        res$cve$ci_lower <- new_lims$ci_lo
        res$cve$ci_upper <- new_lims$ci_up
      }

    }

  }

  if (return_extras) {

    if (attr(dat, "covariates_ph2")) {
      ind_sample <- sample(c(1:nrow(dat_v2_rd)), size=20)
    } else {
      ind_sample <- sample(c(1:nrow(dat_v_rd)), size=20)
    }
    s1 <- min(grid$s)
    s3 <- max(grid$s)
    s2 <- grid$s[which.min(abs((s3-s1)/2-grid$s))]
    Q_n_df <- data.frame(
      "ind" = double(),
      "t" = double(),
      "s" = double(),
      "est" = double()
    )
    Qc_n_df <- Q_n_df
    for (ind in ind_sample) {
      for (t in grid$y) {
        if (attr(dat, "covariates_ph2")) {
          x_val <- as.numeric(dat_v2_rd[ind,c(1:dim_x)])
        } else {
          x_val <- as.numeric(dat_v_rd[ind,c(1:dim_x)])
        }
        Q_val_s1 <- Q_n(t=t, x=x_val, s=s1)
        Q_val_s2 <- Q_n(t=t, x=x_val, s=s2)
        Q_val_s3 <- Q_n(t=t, x=x_val, s=s3)
        Qc_val_s1 <- Qc_n(t=t, x=x_val, s=s1)
        Qc_val_s2 <- Qc_n(t=t, x=x_val, s=s2)
        Qc_val_s3 <- Qc_n(t=t, x=x_val, s=s3)
        Q_n_df[nrow(Q_n_df)+1,] <- c(ind, t, s1, Q_val_s1)
        Q_n_df[nrow(Q_n_df)+1,] <- c(ind, t, s2, Q_val_s2)
        Q_n_df[nrow(Q_n_df)+1,] <- c(ind, t, s3, Q_val_s3)
        Qc_n_df[nrow(Qc_n_df)+1,] <- c(ind, t, s1, Qc_val_s1)
        Qc_n_df[nrow(Qc_n_df)+1,] <- c(ind, t, s2, Qc_val_s2)
        Qc_n_df[nrow(Qc_n_df)+1,] <- c(ind, t, s3, Qc_val_s3)
      }
    }

    if (return_p_value) { res$p <- test_res$p_val }

    res$extras <- list(
      r_Mn = data.frame(
        s = s_out_orig,
        est = c(rep(NA,na_head), sapply(s_out, r_Mn), rep(NA,na_tail))
      ),
      deriv_r_Mn = data.frame(
        s = s_out_orig,
        est = c(rep(NA,na_head), sapply(s_out, deriv_r_Mn), rep(NA,na_tail))
      ),
      Gamma_os_n = data.frame(
        s = s_out_orig,
        est = c(rep(NA,na_head), sapply(s_out, Gamma_os_n), rep(NA,na_tail))
      ),
      f_s_n = data.frame(
        s = s_out_orig,
        est = c(rep(NA,na_head), sapply(s_out, f_s_n), rep(NA,na_tail))
      ),
      Q_n = Q_n_df,
      Qc_n = Qc_n_df
    )

    # !!!!! TEMP
    if (p$edge_corr) {
      res$extras$r_Mn_edge_est <- r_Mn_edge_est
      res$extras$r_Mn_0 <- r_Mn(0)
      res$extras$r_Mn_Gr_0 <- r_Mn_Gr(0)
      res$extras$sigma2_edge_est <- sigma2_edge_est
      res$extras$risk_p <- risk_p
      res$extras$se_p <- se_p
    }

  }

  return(res)

}
