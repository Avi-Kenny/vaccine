#' Estimate CVE/CR using Cox model
#'
#' @description See docs for est_ce and params_ce_cox
#' @noRd
est_cox <- function(
    dat, t_0, cr, cve, s_out, ci_type, placebo_risk_method, return_extras,
    spline_df, spline_knots, edge_ind
) {

  if (!methods::is(dat,"vaccine_dat")) {
    stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by lo",
                "ad_data()."))
  }

  if (!(attr(dat, "groups") %in% c("vaccine", "both"))) {
    stop("Vaccine group data not detected.")
  }

  # !!!!! Validate other inputs; import error handling function from SimEngine

  if (any(is.na(s_out))) { stop("NA values not allowed in s_out.") }

  # Create filtered data objects
  dat_v <- dat[dat$a==1,]
  dat_p <- dat[dat$a==0,]

  # Create spline basis and function
  dat_v_spl <- data.frame("s1"=dat_v$s)
  if ((!is.na(spline_df) && spline_df!=1) || !is.na(spline_knots[[1]])) {

    if (!is.na(spline_df) && !is.na(spline_knots)) {
      stop("Specify either `spline_df` or `spline_knots`, but not both.")
    }

    if (!is.na(spline_knots) && length(spline_knots)<3) {
      stop("`spline_knots` must be a numeric vector of at least length 3.")
    }

    # Create spline basis
    if (!is.na(spline_df)) {
      spl_bnd_knots <- as.numeric(stats::quantile(dat_v$s, c(0.05,0.95), na.rm=T))
      spl_knots <- seq(spl_bnd_knots[1], spl_bnd_knots[2],
                       length.out=spline_df+1)[c(2:spline_df)]
    } else {
      spl_bnd_knots <- spline_knots[c(1,length(spline_knots))]
      spl_knots <- spline_knots[c(2:(length(spline_knots)-1))]
    }
    spl_basis <- splines::ns(
      x = dat_v$s,
      knots = spl_knots,
      intercept = F,
      Boundary.knots = spl_bnd_knots
    )

    for (i in c(1:(dim(spl_basis)[2]))) {
      dat_v_spl[[paste0("s",i)]] <- spl_basis[,i]
    }
    dim_s <- dim(spl_basis)[2]

    if (edge_ind) {

      dim_s <- dim_s + 1
      min_s <- min(dat_v$s, na.rm=T)
      dat_v_spl[[paste0("s",dim_s)]] <- In(dat_v$s==min_s) # !!!!! Should this be dim_s+1 ?????
      s_to_spl <- function(s) {
        spl <- as.numeric(splines::ns(
          x = s,
          knots = spl_knots,
          intercept = F,
          Boundary.knots = spl_bnd_knots
        ))
        return(c(spl, In(s==min_s)))
      }

    } else {

      s_to_spl <- function(s) {
        as.numeric(splines::ns(
          x = s,
          knots = spl_knots,
          intercept = F,
          Boundary.knots = spl_bnd_knots
        ))
      }

    }

  } else {

    if (edge_ind) {

      dim_s <- 2
      min_s <- min(dat_v$s, na.rm=T)
      dat_v_spl[[paste0("s",dim_s)]] <- In(dat_v$s==min_s)
      s_to_spl <- function(s) { c(s, In(s==min_s)) }

    } else {

      dim_s <- 1
      s_to_spl <- function(s) { s }

    }

  }

  # Create phase-two data objects (unrounded)
  dat_v2 <- dat_v[dat_v$z==1,]
  dat_v2_spl <- dat_v_spl[dat_v$z==1,, drop=F]

  # Alias random variables
  WT <- dat_v2$weights
  ST <- dat_v2$strata
  n_vacc <- attr(dat, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")
  X <- dat_v2[,c(1:dim_x), drop=F]
  class(X) <- "data.frame"
  SP <- dat_v2_spl
  V_ <- t(as.matrix(cbind(X,SP)))
  Y_ <- dat_v2$y
  D_ <- dat_v2$delta
  dim_v <- dim(V_)[1]

  # Create set of event times
  i_ev <- which(D_==1)

  # Fit an IPS-weighted Cox model
  model <- survival::coxph(
    formula = stats::formula(paste0("survival::Surv(y,delta)~",
                             paste(names(X),collapse="+"), "+",
                             paste(names(SP),collapse="+"))),
    data = cbind(y=Y_, delta=D_, X, SP),
    weights = WT
  )
  coeffs <- model$coefficients
  beta_n <- as.numeric(coeffs)

  if (any(is.na(coeffs))) {

    if (any(is.na(coeffs[which(substr(names(coeffs),1,1)=="x")]))) {
      na_coeffs <- names(coeffs)[which(is.na(coeffs))]
      stop(paste0("The following covariate coefficients were NA.: ",
                  paste(na_coeffs, collapse=","),
                  " Try removing these coefficients and rerunning."))
      # !!!!! Automatically omit from the model, as is done for spline coeffs
    }

    if (any(is.na(coeffs[which(substr(names(coeffs),1,1)=="s")]))) {
      warning(paste0("Some spline coefficients were NA; these coefficients hav",
                     "e been automatically omitted from the model"))
      na_inds <- as.numeric(which(is.na(coeffs)))
      s_na_inds <- which(is.na(coeffs[which(substr(names(coeffs),1,1)=="s")]))
      s_na_inds <- as.numeric(s_na_inds)
      beta_n <- beta_n[-na_inds]
      SP <- SP[,-s_na_inds]
      names(SP) <- paste0("s", c(1:length(SP)))
      V_ <- t(as.matrix(cbind(X,SP)))
      dim_v <- dim(V_)[1]
      s_to_spl2 <- s_to_spl
      s_to_spl <- function(s) { s_to_spl2(s)[-s_na_inds] }

    }

  }

  LIN <- as.numeric(t(beta_n)%*%V_)

  # Intermediate functions
  {

    S_0n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            (1/n_vacc) * sum(WT*In(Y_>=x)*exp(LIN))
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()

    S_1n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            (1/n_vacc) * as.numeric(V_ %*% (WT*In(Y_>=x)*exp(LIN)))
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()

    S_2n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            res <- matrix(NA, nrow=dim_v, ncol=dim_v)
            for (i in c(1:dim_v)) {
              for (j in c(1:dim_v)) {
                if (!is.na(res[j,i])) {
                  res[i,j] <- res[j,i]
                } else {
                  res[i,j] <- (1/n_vacc)*sum(WT*In(Y_>=x)*V_[i,]*V_[j,]*exp(LIN))
                }
              }
            }
            return(res)
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()

    m_n <- function(x) { S_1n(x) / S_0n(x) }

  }

  # Estimated information matrix (for an individual)
  I_tilde <- Reduce("+", lapply(i_ev, function(i) {
    WT[i] * ( (S_2n(Y_[i])/S_0n(Y_[i])) - m_n(Y_[i]) %*% t(m_n(Y_[i])) )
  }))
  I_tilde <- (1/n_vacc)*I_tilde
  I_tilde_inv <- solve(I_tilde)

  # Score function (Cox model)
  l_n <- function(v_i,d_i,y_i) {
    d_i*(v_i-m_n(y_i)) - (1/n_vacc)*Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*exp(sum(v_i*beta_n))*In(Y_[j]<=y_i) * (v_i-m_n(Y_[j]))) /
        S_0n(Y_[j])
    }))
  }

  # Influence function: beta_hat (regular Cox model)
  l_tilde <- (function() {
    .cache <- new.env()
    function(v_i,d_i,y_i) {
      key <- paste(c(v_i,d_i,y_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,d_i,y_i) {
          I_tilde_inv %*% l_n(v_i,d_i,y_i)
        })(v_i,d_i,y_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Create p_n and p1_n vectors
  n_strata <- max(dat_v$strata)
  p_n <- c()
  p1_n <- c()
  for (c in c(1:n_strata)) {
    p_n[c] <- mean(c==dat_v$strata)
    p1_n[c] <- mean(c==dat_v$strata & dat_v$z==1)
  }

  # Influence function: beta_hat (est. weights)
  infl_fn_beta <- (function() {
    .cache <- new.env()

    exp_terms <- list()
    for (i in c(1:max(dat_v$strata))) {
      js <- which(ST==i)
      if (length(js)>0) {
        exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
          l_tilde(V_[,j],D_[j],Y_[j])
        }))
      } else {
        exp_terms[[i]] <- as.matrix(rep(0,dim_v)) # Avoid NAs in small samples
      }
    }

    function(v_i,z_i,d_i,y_i,wt_i,st_i) {
      key <- paste(c(v_i,z_i,d_i,y_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,z_i,d_i,y_i,wt_i,st_i) {

          if (z_i) {
            piece_1 <- wt_i*l_tilde(v_i,d_i,y_i)
          } else {
            piece_1 <- as.matrix(rep(0,dim_v))
          }

          if (p1_n[st_i]!=0) {
            piece_2 <- (1 - (z_i*p_n[st_i])/p1_n[st_i]) * exp_terms[[st_i]]
          } else {
            piece_2 <- as.matrix(rep(0,dim_v)) # Avoid NAs in small samples
          }

          return(as.numeric(piece_1+piece_2))

        })(v_i,z_i,d_i,y_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance constant: mu_n
  mu_n <- -1 * (1/n_vacc) * as.numeric(Reduce("+", lapply(i_ev, function(j) {
    (WT[j] * In(Y_[j]<=t_0) * m_n(Y_[j])) / S_0n(Y_[j])
  })))

  # Nuisance function: v_1n
  v_1n <- (function() {
    .cache <- new.env()
    function(st_i,z_i,y_j) {
      key <- paste(c(st_i,z_i,y_j), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,z_i,y_j) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/n_vacc) * (p1_n[st_i]-p_n[st_i]*z_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  In(Y_[k]>=y_j) * exp(sum(beta_n*V_[,k]))
                }))
              )
            )
          } else { return(0) }
        })(st_i,z_i,y_j)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance function: v_2n
  v_2n <- (function() {
    .cache <- new.env()
    function(st_i,z_i) {
      key <- paste(c(st_i,z_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,z_i) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/n_vacc) * (p1_n[st_i]-p_n[st_i]*z_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  ( D_[k] * In(Y_[k]<=t_0) ) / S_0n(Y_[k])
                }))
              )
            )
          } else { return(0) }
        })(st_i,z_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Breslow estimator
  Lambda_n <- function(t) {
    (1/n_vacc) * sum(unlist(lapply(i_ev, function(i) {
      WT[i] * ( In(Y_[i]<=t) / S_0n(Y_[i]) )
    })))
  }

  # Survival estimator (at a point)
  Q_n <- (function() {
    Lambda_n_t_0 <- Lambda_n(t_0)
    function(z) { exp(-exp(sum(z*beta_n))*Lambda_n_t_0) }
  })()

  # Influence function: Breslow estimator (est. weights)
  infl_fn_Lambda <- (function() {
    .cache <- new.env()
    function(v_i,z_i,d_i,y_i,wt_i,st_i) {
      key <- paste(c(v_i,z_i,d_i,y_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,z_i,d_i,y_i,wt_i,st_i) {
          pc_4 <- (1/n_vacc) * sum(unlist(lapply(i_ev, function(j) {
            ( WT[j] * In(Y_[j]<=t_0) * v_1n(st_i,z_i,Y_[j]) ) / (S_0n(Y_[j]))^2
          })))
          pc_5 <- sum(mu_n*infl_fn_beta(v_i,z_i,d_i,y_i,wt_i,st_i))
          pc_2 <- v_2n(st_i,z_i)

          if (z_i==1) {
            pc_1 <- ( wt_i * d_i * In(y_i<=t_0) ) / S_0n(y_i)
            pc_3 <- (1/n_vacc) * sum(unlist(lapply(i_ev, function(j) {
              (WT[j]*In(Y_[j]<=t_0)*wt_i*In(y_i>=Y_[j])*exp(sum(beta_n*v_i))) /
                (S_0n(Y_[j]))^2
            })))
            return(pc_1+pc_2+pc_5-pc_3-pc_4)
          } else {
            return(pc_2+pc_5-pc_4)
          }
        })(v_i,z_i,d_i,y_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Compute marginalized risk
  res_cox <- list()
  Lambda_n_t_0 <- Lambda_n(t_0)
  res_cox$est_marg <- unlist(lapply(s_out, function(s) {
    (1/n_vacc) * sum((apply(dat_v, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      s_spl <- s_to_spl(s)
      exp(-1*exp(sum(beta_n*c(x_i,s_spl)))*Lambda_n_t_0)
    })))
  }))

  # Compute variance estimate
  if (ci_type!="none") {

    res_cox$var_est_marg <- unlist(lapply(s_out, function(s) {

      # Precalculate pieces dependent on s
      s_spl <- s_to_spl(s)
      K_n <- (1/n_vacc) * Reduce("+", apply2(dat_v, 1, function(r) {
        x_i <- as.numeric(r[1:dim_x])
        Q <- Q_n(c(x_i,s_spl))
        explin <- exp(sum(c(x_i,s_spl)*beta_n))
        K_n1 <- Q
        K_n2 <- Q * explin
        K_n3 <- Q * explin * c(x_i,s_spl)
        return(c(K_n1,K_n2,K_n3))
      }, simplify=F))
      K_n1 <- K_n[1]
      K_n2 <- K_n[2]
      K_n3 <- K_n[3:length(K_n)]

      (1/n_vacc^2) * sum((apply(dat_v, 1, function(r) {

        x_i <- as.numeric(r[1:dim_x])
        if (!is.na(r[["s"]])) {
          s_i <- s_to_spl(r[["s"]])
        } else {
          s_i <- NA
        }
        z_i <- r[["z"]]
        d_i <- r[["delta"]]
        y_i <- r[["y"]]
        wt_i <- r[["weights"]]
        st_i <- r[["strata"]]

        pc_1 <- Q_n(c(x_i,s_spl))
        pc_2 <- Lambda_n_t_0 * sum(
          K_n3 * infl_fn_beta(c(x_i,s_i),z_i,d_i,y_i,wt_i,st_i)
        )
        pc_3 <- K_n2 * infl_fn_Lambda(c(x_i,s_i),z_i,d_i,y_i,wt_i,st_i)
        pc_4 <- K_n1

        return((pc_1-pc_2-pc_3-pc_4)^2)

      })))

    }))

  }

  # Extract estimates and SEs
  ests_cr <- 1-res_cox$est_marg

  # Generate confidence limits
  if (ci_type=="none") {

    ci_lo_cr <- ci_up_cr <- ses_cr <- rep(NA, length(ests_cr))

  } else {

    ses_cr <- sqrt(res_cox$var_est_marg)

    if (ci_type=="regular") {
      ci_lo_cr <- ests_cr - 1.96*ses_cr
      ci_up_cr <- ests_cr + 1.96*ses_cr
    } else if (ci_type=="truncated") {
      ci_lo_cr <- pmin(pmax(ests_cr - 1.96*ses_cr, 0), 1)
      ci_up_cr <- pmin(pmax(ests_cr + 1.96*ses_cr, 0), 1)
    } else if (ci_type=="transformed") {
      ci_lo_cr <- expit(logit(ests_cr) - 1.96*deriv_logit(ests_cr)*ses_cr)
      ci_up_cr <- expit(logit(ests_cr) + 1.96*deriv_logit(ests_cr)*ses_cr)
    } else if (ci_type=="transformed 2") {
      ci_lo_cr <- expit2(logit2(ests_cr) - 1.96*deriv_logit2(ests_cr)*ses_cr)
      ci_up_cr <- expit2(logit2(ests_cr) + 1.96*deriv_logit2(ests_cr)*ses_cr)
    }

  }

  # Create results object
  res <- list()
  if (cr) {
    res$cr <- list(s=s_out, est=ests_cr, se=ses_cr, ci_lower=ci_lo_cr,
                   ci_upper=ci_up_cr)
  }

  # Compute CVE
  if (cve) {
    res$cve <- list(s=s_out)
    if (attr(dat, "groups")!="both") { stop("Placebo group not detected.") }
    ov <- est_overall(dat=dat, t_0=t_0, method=placebo_risk_method, ve=F)
    risk_p <- ov[ov$group=="placebo","est"]
    se_p <- ov[ov$group=="placebo","se"] # This equals sd_p/n_orig
    res$cve$est <- 1 - res$cr$est/risk_p

    if (ci_type=="none") {

      na_vec <- rep(NA, length(res$cve$s_out))
      res$cve$se <- na_vec
      res$cve$ci_lower <- na_vec
      res$cve$ci_upper <- na_vec

    } else {

      res$cve$se <- sqrt(ses_cr^2/risk_p^2 + (res$cr$est^2*se_p^2)/risk_p^4)
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

      # # !!!!! OLD !!!!!
      # res$cve$se <- NA
      # res$cve$ci_lower <- 1 - res$cr$ci_up/risk_p # !!!!! Add placebo group correction
      # res$cve$ci_upper <- 1 - res$cr$ci_lo/risk_p # !!!!! Add placebo group correction
    }

  }

  # Return extras
  if (return_extras) { res$extras <- list(model=model) }

  return(res)

}
