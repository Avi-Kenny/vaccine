if (F) {

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
Pi <- function(sampling, delta, y, x, t_0) {

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
    ev <- In(delta==1 & y<=t_0)
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
construct_etastar_n <- function(Q_n, t_0) {

  fnc <- function(u,x) {
    u <- round(u,-log10(C$appx$s))
    if (u==0) {
      return(0)
    } else {
      s_seq <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
      integral <- C$appx$s * sum(sapply(s_seq, function(s) {
        Q_n(t_0, x, s)
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

}
