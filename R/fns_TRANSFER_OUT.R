if (F) {

#' Helper function for debugging; prints timestamps
#'
#' @param num Number
#' @param msg Message
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  print(str)
}



#' Probability of sampling
#'
#' @param sampling One of c("iid", "two-phase (6%)", "two-phase (72%)",
#'     "two-phase (70% random)", "two-phase (6% random)")
#' @param delta Component of dataset returned by generate_data()
#' @param y Component of dataset returned by generate_data()
#' @param x Component of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
#' @note
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
Pi <- function(sampling, delta, y, x) {

  if (sampling=="iid") {
    probs <- rep(1, length(delta))
  } else if (sampling=="two-phase (70% random)") {
    probs <- rep(0.7, length(delta))
  } else if (sampling=="two-phase (50% random)") {
    probs <- rep(0.5, length(delta))
  } else if (sampling=="two-phase (25% random)") {
    probs <- rep(0.25, length(delta))
  } else if (sampling=="two-phase (6% random)") {
    probs <- rep(0.06, length(delta))
  } else if (sampling=="cycle") {
    probs <- rep(c(0.8,0.6,0.4,0.2), length.out=length(delta))
  } else {
    ev <- In(delta==1 & y<=C$t_0)
    if (sampling=="two-phase (6%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-3.85)
    } else if (sampling=="two-phase (72%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-0.1)
    } else if (sampling=="two-phase (50%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-1)
    } else if (sampling=="two-phase (25%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-2.2)
    } else if (sampling=="x1") {
      probs <- (0.2 + 0.6*x$x1)
    } else if (sampling=="x2") {
      probs <- (0.2 + 0.6*x$x2)
    }

  }

  return(probs)

}



#' Return IP weights
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param scale One of c("none", "stabilized")
#' @param type One of c("true", "estimated")
#' @param return_strata Whether the discrete two-phase sampling strata
#'     membership variable should be returned
#' @return A sum-to-one vector of weights
#' @note
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
wts <- function(dat_orig, scale="stabilized", type="true", return_strata=F) {

  sampling <- attr(dat_orig,"sampling")
  Pi_0 <- Pi(sampling, dat_orig$delta, dat_orig$y, dat_orig$x)
  strata1 <- In(factor(Pi_0))

  if (type=="true") {

    weights <- dat_orig$z / Pi_0

  } else if (type=="estimated") {

    Pi_vals <- c()
    for (i in c(1:max(strata1))) {
      Pi_vals[i] <- sum(In(strata1==i)*dat_orig$z) / sum(In(strata1==i))
      if (Pi_vals[i]==0) {
        # Hack to avoid NA values in small sample sizes
        warning(paste0("stratum ", i, " had no one sampled."))
        Pi_vals[i] <- 1
      }
    }
    weights <- dat_orig$z / Pi_vals[strata1]

  }

  if (scale=="stabilized") {
    sm <- sum(weights) / length(dat_orig$z)
    weights <- weights / sm
  }

  if (!return_strata) {
    return(weights)
  } else {
    return(list(weights=weights, strata=strata1))
  }

}



#' Construct q_tilde_n nuisance estimator function
#'
#' @param x x
#' @return q_tilde_n nuisance estimator function
construct_q_tilde_n <- function(type="standard", f_n_srv, f_sIx_n, omega_n) {

  if (type=="standard") {

    seq_01 <- grid$s[2:length(grid$s)]
    # seq_01 <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))

    fnc <- function(x, y, delta, u) {

      denom <- C$appx$s * sum(sapply(seq_01, function(s) {
        f_n_srv(y, delta, x, s) * f_sIx_n(s,x)
      }))

      if (denom==0) {
        return (0)
      } else {
        if (u==0) {
          return(0)
        } else {
          seq_0u <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
          num <- C$appx$s * sum(sapply(seq_0u, function(s) {
            omega_n(x,s,y,delta) * f_n_srv(y, delta, x, s)
          }))
          return(num/denom)
        }
      }

    }

  }

  if (type=="zero") {

    fnc <- function(x, y, delta, u) { 0 }

  }

  return(memoise2(fnc))

}






#' Construct nuisance estimator eta*_n
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param vals List of values to pre-compute function on
#' @return Estimator function of nuisance eta*_0
construct_etastar_n <- function(Q_n, vals=NA, tmp) {

  fnc <- function(u,x) {
    u <- round(u,-log10(C$appx$s))
    if (u==0) {
      return(0)
    } else {
      s_seq <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
      integral <- C$appx$s * sum(sapply(s_seq, function(s) {
        Q_n(C$t_0, x, s)
      }))
      return(u-integral)
    }
  }

  return(memoise2(fnc))

}



#' lambda estimator
#'
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @return Value of lambda
lambda <- function(dat, k, G) {

  n_orig <- attr(dat, "n_orig")
  lambda <- (1/n_orig) * sum( dat$weights * (G(dat$s))^k )
  return(lambda)

}



#' Construct Theta_os_n primitive one-step estimator
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param vals List of values to pre-compute function on
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_sIx_n Conditional density estimator returned by construct_f_sIx_n
#' @param etastar_n A nuisance estimator returned by construct_etastar_n()
#' @return Gamma_os_n estimator
#' @note This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Theta_os_n <- function(dat, dat_orig, omega_n=NA, f_sIx_n=NA,
                                 q_tilde_n=NA, etastar_n=NA, vals=NA) {

  n_orig <- attr(dat_orig, "n_orig")
  piece_1 <- (dat$weights*omega_n(dat$x,dat$s,dat$y,dat$delta)) /
    f_sIx_n(dat$s,dat$x)
  piece_2 <- (1-dat_orig$weights)

  # Remove large intermediate objects
  rm(omega_n,f_sIx_n)

  fnc <- function(u) {
    (1/(n_orig)) * sum(piece_1*In(dat$s<=u)) +
      (1/n_orig) * sum(
        piece_2 * q_tilde_n(dat_orig$x, dat_orig$y, dat_orig$delta, u) +
          etastar_n(rep(u,length(dat_orig$z)),dat_orig$x)
      )
  }

  # !!!!! DEBUG
  if (F) {

    pc1 <- c()
    pc2 <- c()
    sm <- c()
    for (u in round(seq(0,1,0.1),2)) {
      pc1_ <- (1/(n_orig)) * sum(piece_1*In(dat$s<=u))
      pc2_ <- (1/n_orig) * sum(
        piece_2 * q_tilde_n(dat_orig$x, dat_orig$y, dat_orig$delta, u) +
          etastar_n(rep(u,length(dat_orig$z)),dat_orig$x)
      )
      pc1 <- c(pc1, pc1_)
      pc2 <- c(pc2, pc2_)
      sm <- c(sm, pc1_+pc2_)
    }

    grid <- round(seq(0,1,0.1),2)
    df_plot <- data.frame(
      x = rep(grid, 3),
      y = c(pc1,pc2,sm),
      which = rep(c("omega piece", "etastar piece", "Theta_n (sum)"), each=11)
    )
    ggplot(df_plot, aes(x=x, y=y, color=factor(which))) +
      geom_line() +
      labs(color="Piece")

  }

  return(memoise2(fnc))

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_infl_fn_Theta <- function(omega_n, f_sIx_n, q_tilde_n, etastar_n,
                                    Theta_os_n) {

  fnc <- function(u,x,y,delta,s,wt) {
    if (wt==0) {
      piece_1 <- 0
      piece_2 <- 0
    } else {
      piece_1 <- In(s<=u)
      piece_2 <- omega_n(x,s,y,delta)/f_sIx_n(s,x)
    }
    wt*piece_1*piece_2 +
      (1-wt) * q_tilde_n(x,y,delta,u) +
      etastar_n(u,x) -
      Theta_os_n(round(u,-log10(C$appx$s)))
  }

  return(memoise2(fnc))

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_infl_fn_beta_n <- function(infl_fn_Theta) {

  u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
  m <- length(u_mc)
  lambda_1 <- mean(u_mc) # ~1/2
  lambda_2 <- mean((u_mc)^2) # ~1/3
  lambda_3 <- mean((u_mc)^3) # ~1/4

  fnc <- function(s, y, delta, weight, x) {

    s_m <- rep(s,m)
    y_m <- rep(y,m)
    delta_m <- rep(delta,m)
    weight_m <- rep(weight,m)
    x_m <- as.data.frame(matrix(rep(x,m), ncol=length(x), byrow=T))

    return((1/m) * sum(
      infl_fn_Theta(u=u_mc, x_m, y_m, delta_m, s_m, weight_m) * (
        (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
          (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
      ))
    )

  }

  return(fnc)

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_Gamma_cf_k <- function(dat_train, dat_test, vals=NA, omega_n, g_n,
                                 r_tilde_Mn, eta_n) {

  # !!!!! Needs to be updated

  n_test <- length(dat_test$s)
  dat_test$weights <- wts(dat_test) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d1 <- ss(dat_test, which(dat_test$z==1))
  weights_1 <- d1$weights

  n_train <- length(dat_train$s)
  dat_train$weights <- wts(dat_train) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d2 <- ss(dat_train, which(dat_train$z==1))
  weights_2 <- d2$weights

  fnc <- function(s) {

    piece_1 <- (1/n_test) * sum(weights_1 * (
      In(d1$s<=s) *
        ((omega_n(d1$x,d1$s,d1$y,d1$delta) / g_n(d1$s,d1$x)) +
           r_tilde_Mn(d1$s)
        ) +
        eta_n(s,d1$x)
    ))

    piece_2 <- (1/n_train) * sum(weights_2 * In(d2$s<=s)*r_tilde_Mn(d2$s))

    return(piece_1-piece_2)

  }

  return(memoise2(fnc))

}



#' Construct cross-fitted Gamma_0 estimator
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param params The same params object passed to est_curve
#' @param vlist A list of dataframes returned by create_val_list(); Q_n REQUIRED
#'     but others can be NA
#' @return A cross-fitted Gamma_0 estimator
construct_Gamma_cf <- function(dat_orig, params, vlist) {

  # !!!!! Update this

  # Prep for cross-fitting
  Gamma_cf_k <- list()
  n_orig <- attr(dat_orig, "n_orig")
  rows <- c(1:n_orig)
  folds <- sample(cut(rows, breaks=params$cf_folds, labels=FALSE))

  # Loop through folds
  for (k in 1:params$cf_folds) {

    # Split data
    dat_orig_train <- ss(dat_orig, -which(folds==k))
    dat_orig_test <- ss(dat_orig, which(folds==k))
    dat_train <- ss(dat_orig_train, which(dat_orig_train$z==1))
    dat_test <- ss(dat_orig_test, which(dat_orig_test$z==1))

    # Construct component functions
    Phi_n <- construct_Phi_n(dat_train, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat_train, which="inverse",
                                 type=params$ecdf_type)
    srvSL <- construct_Q_n(dat_train, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    r_tilde_Mn <- construct_r_tilde_Mn(dat_orig_train, vlist$S_grid, Q_n)
    # eta_n <- construct_eta_n(dat_train, vlist$SX_grid, Q_n) # ARCHIVED
    f_sIx_n <- construct_f_sIx_n(dat_train, vlist$SX_grid, type=params$g_n_type,
                                 k=15, s_scale=s_scale, s_shift=s_shift)
    f_s_n <- construct_f_s_n(dat_orig_train, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n,
                                 type=params$omega_n_type)

    # Construct K functions
    Gamma_cf_k[[k]] <- construct_Gamma_cf_k(
      dat_train, dat_test, vlist$S_grid, omega_n, g_n, r_tilde_Mn, eta_n
    )

    # Remove objects
    rm(Phi_n,Phi_n_inv,Q_n,Qc_n,r_tilde_Mn,eta_n,f_sIx_n,f_s_n,g_n,omega_n)

  }

  # Construct cross-fitted Gamma_os_n
  return(Vectorize(function(s) {
    mean(sapply(c(1:params$cf_folds), function(k) {
      Gamma_cf_k[[k]](s)
    }))
  }))

}





#' Construct influence function corresponding to test statistic beta_en
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @return Influence function
construct_infl_fn_beta_en <- function(infl_fn_Theta, infl_fn_r_Mn_edge) {

  return(function(z, x, y, delta, s, weight) {
    infl_fn_Theta(u=1, x, y, delta, s, weight) -
      infl_fn_r_Mn_edge(z, weight, s, x, y, delta)
  })

}





# Estimate the variance of the Cox model based marginalized survival estimator
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param t The end time of interest
#' @param points The A-values of interest
#' @param se_beta Estimates related to parameter vector beta_n
#' @param se_bshz Estimates related to Breslow baseline cuml. hazard estimator
#' @param se_surv Estimates related to survival at a point
#' @param se_marg Estimates related to marginalized survival
#' @param return_extras For debugging
#' @param parallel Parallelize internal lapply calls
#' @param verbose Print status updates and progress bars
#' @return A list containing the selected return values
cox_var <- function(dat_orig, dat, t, points, se_beta=F, se_bshz=F,
                    se_surv=F, se_marg=F, return_extras=F, parallel=F,
                    verbose=F) {
  # Setup
  if (verbose) {
    print(paste("Check 0 (start):", Sys.time()))
    pbapply::pboptions(type="txt", char="#", txt.width=40, style=3)
  } else {
    pbapply::pboptions(type="none")
  }
  if (parallel) {
    n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    if (verbose) {
      print(cl)
    }
  } else {
    cl <- NULL
  }

  # Fit a Cox model
  # Note: scaling the weights affects the SEs but not the estimates; thus, this
  #       is only needed for debugging
  model <- coxph(
    formula = formula(paste0("Surv(y,delta)~",
                             paste(names(dat$x),collapse="+"),"+s")),
    data = cbind(y=dat$y, delta=dat$delta,
                 dat$x, s=dat$s),
    weights = dat$weights * (length(dat$weights)/sum(dat$weights))
  )
  beta_n <- as.numeric(model$coefficients)

  if (verbose) { print(paste("Check 1 (Cox model fit):", Sys.time())) }

  # !!!!! Stabilize weights

  # Alias random variables
  WT <- dat$weights
  ST <- dat$strata
  N <- length(dat_orig$s)
  n <- length(dat$s)
  Z_ <- t(as.matrix(cbind(dat$x,s=dat$s)))
  T_ <- dat$y
  Ds_ <- dat$delta
  lin <- as.numeric(t(beta_n)%*%Z_)
  d <- dim(Z_)[1]

  # Intermediate functions
  {
    S_2n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            res <- matrix(NA, nrow=d, ncol=d)
            for (i in c(1:d)) {
              for (j in c(1:d)) {
                if (!is.na(res[j,i])) {
                  res[i,j] <- res[j,i]
                } else {
                  res[i,j] <- (1/N)*sum(WT*In(T_>=x)*Z_[i,]*Z_[j,]*exp(lin))
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

    S_0n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            (1/N) * sum(WT*In(T_>=x)*exp(lin))
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
            (1/N) * as.numeric(Z_ %*% (WT*In(T_>=x)*exp(lin)))
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()

    m_n <- function(x) {
      S_1n(x) / S_0n(x)
    }

    h <- function(x) {
      (S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x))
    }

  }

  # Create set of event times
  i_ev <- which(Ds_==1)
  t_ev <- T_[which(Ds_==1)]

  # Estimated information matrix (for an individual)
  I_tilde <- Reduce("+", lapply(i_ev, function(i) {
    WT[i] * h(T_[i])
  }))
  I_tilde <- (1/N)*I_tilde
  I_tilde_inv <- solve(I_tilde)

  # Score function (Cox model)
  l_star <- function(z_i,ds_i,t_i) {
    ds_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*exp(sum(z_i*beta_n))*In(T_[j]<=t_i)*(z_i-m_n(T_[j]))) /
        S_0n(T_[j])
    }))
  }

  # Influence function: beta_hat (regular Cox model)
  l_tilde <- (function() {
    .cache <- new.env()
    function(z_i,ds_i,t_i) {
      key <- paste(c(z_i,ds_i,t_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,ds_i,t_i) {
          I_tilde_inv %*% l_star(z_i,ds_i,t_i)
        })(z_i,ds_i,t_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Create p_n and p1_n vectors
  n_strata <- max(dat_orig$strata)
  p_n <- c()
  p1_n <- c()
  for (c in c(1:n_strata)) {
    p_n[c] <- mean(c==dat_orig$strata)
    p1_n[c] <- mean(c==dat_orig$strata & dat_orig$z==1)
  }

  # Influence function: beta_hat (est. weights)
  lstar_tilde <- (function() {
    .cache <- new.env()

    exp_terms <- list()
    for (i in c(1:max(dat_orig$strata))) {
      js <- which(ST==i)
      if (length(js)>0) {
        exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
          l_tilde(Z_[,j],Ds_[j],T_[j])
        }))
      } else {
        exp_terms[[i]] <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
      }
    }

    function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          if (d_i) {
            piece_1 <- wt_i*l_tilde(z_i,ds_i,t_i)
          } else {
            piece_1 <- as.matrix(rep(0,d))
          }

          if (p1_n[st_i]!=0) {
            piece_2 <- (1 - (d_i*p_n[st_i])/p1_n[st_i]) * exp_terms[[st_i]]
          } else {
            piece_2 <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
          }

          return(as.numeric(piece_1+piece_2))
        })(z_i,d_i,ds_i,t_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance constant: mu_n
  mu_n <- (1/N) * as.numeric(Reduce("+", lapply(t_ev, function(t_j) {
    (In(t_j<=t) * S_1n(t_j)) / (S_0n(t_j))^2
  })))

  # Nuisance function: v_n
  v_n <- (function() {
    .cache <- new.env()
    function(st_i,d_i,t_j) {
      key <- paste(c(st_i,d_i,t_j), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,d_i,t_j) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  In(T_[k]>=t_j) * exp(sum(beta_n*Z_[,k]))
                }))
              )
            )
          } else { return(0) }
        })(st_i,d_i,t_j)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance function: vstar_n
  vstar_n <- function(st_i,d_i,t) {
    k_set <- which(ST==st_i & Ds_==1)
    if (length(k_set)>0) {
      return(
        (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
          unlist(lapply(k_set, function(k) {
            In(T_[k]<=t) / S_0n(T_[k])
          }))
        )
      )
    } else { return(0) }
  }

  # Breslow estimator
  Lambda_n <- Vectorize((function() {
    .cache <- new.env()
    function(x) {
      key <- as.character(x)
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(x) {
          (1/N) * sum(unlist(lapply(i_ev, function(i) {
            (WT[i] * In(T_[i]<=x)) / S_0n(T_[i])
          })))
        })(x)
        .cache[[key]] <- val
      }
      return(val)
    }
  })())

  # Survival estimator (at a point)
  Q_n <- (function() {
    .cache <- new.env()
    function(z) {
      key <- paste(z, collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z) {
          exp(-exp(sum(z*beta_n))*Lambda_n(t))
        })(z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: Breslow estimator (est. weights)
  infl_fn_Bres <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          pc_2 <- vstar_n(st_i,d_i,t)
          pc_4 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
            ( In(t_j<=t) * v_n(st_i,d_i,t_j) ) / (S_0n(t_j))^2
          })))
          pc_5 <- sum(mu_n*lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))

          if (d_i==1) {
            pc_1 <- ( wt_i * ds_i * In(t_i<=t) ) / S_0n(t_i)
            pc_3 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
              (In(t_j<=t)*wt_i*In(t_i>=t_j)*exp(sum(beta_n*z_i))) /
                (S_0n(t_j))^2
            })))
            return(pc_1+pc_2-pc_3-pc_4-pc_5)
          } else {
            return(pc_2-pc_4-pc_5)
          }
        })(z_i,d_i,ds_i,t_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: survival at a point (est. weights)
  omega_n <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i,z), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
          explin <- exp(sum(z*beta_n))
          piece_1 <- Lambda_n(t) * explin *
            (t(z)%*%lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))[1]
          piece_2 <- explin * infl_fn_Bres(z_i,d_i,ds_i,t_i,wt_i,st_i)
          return(-1*Q_n(z)*(piece_1+piece_2))
        })(z_i,d_i,ds_i,t_i,wt_i,st_i,z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: marginalized survival (est. weights)
  infl_fn_marg <- (function() {
    .cache <- new.env()
    x_ <- as.list(as.data.frame(t(dat_orig$x)))
    function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
      key <- paste(c(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
          piece_1 <- Q_n(c(x_i,s))
          piece_2 <- (1/N) * sum(unlist(lapply(c(1:N), function(j) {
            x_j <- x_[[j]]
            omgn <- omega_n(c(x_i,s_i),d_i,ds_i,t_i,wt_i,st_i,c(x_j,s))
            return(omgn-Q_n(c(x_j,s)))
          })))
          return(piece_1+piece_2)
        })(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Construct results object
  res <- list(model=model, beta_n=beta_n)

  if (verbose) { print(paste("Check 2 (functions declared):", Sys.time())) }

  # Variance estimate: beta_hat
  if (se_beta) {

    res$est_beta <- beta_n

    if (verbose) { print(paste("Check 3a (var est START: beta):", Sys.time())) }
    res$var_est_beta <- (1/N^2) * Reduce("+", pblapply(c(1:N), function(i) {
      (lstar_tilde(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl))
    if (verbose) { print(paste("Check 4a (var est END: beta):", Sys.time())) }

  }

  # Variance estimate: Breslow estimator
  if (se_bshz) {

    res$est_bshz <- Lambda_n(t)

    if (verbose) { print(paste("Check 3b (var est START: bshz):", Sys.time())) }
    res$var_est_bshz <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (infl_fn_Bres(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4b (var est END: bshz):", Sys.time())) }

  }

  # Variance estimate: survival at a point
  if (se_surv) {

    z_0 <- c(0.3,1,0.5)

    res$est_surv <- Q_n(z_0)

    if (verbose) { print(paste("Check 3c (var est START: surv):", Sys.time())) }
    res$var_est_surv <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (omega_n(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i],
        z = z_0
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4c (var est END: surv):", Sys.time())) }

  }

  # Variance estimate: marginalized survival
  if (se_marg) {

    # !!!! basehaz vs. Lambda_hat
    bh <- basehaz(model, centered=FALSE)
    index <- max(which((bh$time<C$t_0)==T))
    est_bshz <- bh$hazard[index]
    N <- sum(dat$weights)
    res$est_marg <- unlist(lapply(points, function(s) {
      (1/N) * sum(unlist(lapply(c(1:N), function(i) {
        exp(-1*exp(sum(beta_n*c(as.numeric(dat_orig$x[i,]),s)))*est_bshz)
      })))
    }))

    # !!!!! Pre-calculate omega_n values

    if (verbose) { print(paste("Check 3d (var est START: marg):", Sys.time())) }
    res$var_est_marg <- unlist(lapply(points, function(s) {
      if (verbose) { print(paste0("Check 4d (point=",s,"): ", Sys.time())) }
      (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
        (infl_fn_marg(
          x_i = as.numeric(dat_orig$x[i,]),
          s_i = dat_orig$s[i],
          d_i = dat_orig$z[i],
          ds_i = dat_orig$delta[i],
          t_i = dat_orig$y[i],
          wt_i = dat_orig$weights[i],
          st_i = dat_orig$strata[i],
          s = s
        ))^2
      }, cl=cl)))
    }))
    if (verbose) { print(paste("Check 5d (var est END: marg):", Sys.time())) }

  }

  if (return_extras) {
    res$S_0n <- S_0n
    res$S_1n <- S_1n
    res$S_2n <- S_2n
    res$I_tilde <- I_tilde
    res$I_tilde_inv <- I_tilde_inv
    res$l_tilde <- l_tilde
    res$omega_n <- omega_n
    res$Lambda_n <- Lambda_n
  }

  return(res)

}

}
