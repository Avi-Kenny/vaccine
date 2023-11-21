#' Construct influence function corresponding to r_Mn_edge
#'
#' @param Q_n Conditional survival function estimator returned by
#'     `construct_Q_n`
#' @param g_sn Propensity score estimator returned by `construct_g_sn`
#' @param omega_n A nuisance influence function returned by `construct_omega_n`
#' @param r_Mn_edge_est Estimate returned by one-step estimator `r_Mn_edge`
#' @return Influence function estimator
#' @noRd
construct_infl_fn_r_Mn_edge <- function(Q_n, g_sn, omega_n, g_n,
                                        r_Mn_edge_est, p_n, t_0) {

  return(memoise2(function(z, weight, s, x, y, delta) {

    if (z==0) {
      s <- 0
      pi_ <- 1
    } else {
      pi_ <- 1/weight
    }

    return(
      1 - Q_n(t_0,x,s=0) +
        ((z*In(s==0)+g_sn(x,y,delta)*(pi_-z)) * omega_n(x,s=0,y,delta)) /
        (pi_*(1-p_n)*g_n(s=0,x)) - r_Mn_edge_est
    )

  }))

}



#' Construct influence function corresponding to Kaplan-Meier estimator
#'
#' @param dat_combined Combined dataset; !!!!! TEMP
#' @return Influence function estimator
#' @note Used by est_med()
#' @noRd
construct_infl_fn_risk_p <- function(dat_p_rd, Q_noS_n, omega_noS_n, t_0,
                                     p_plac) {

  dim_x <- attr(dat_p_rd, "dim_x")

  mean_Q_n <- mean(apply(dat_p_rd, 1, function(r) {
    Q_noS_n(t_0,as.numeric(r[1:dim_x]))
  }))

  infl_fn <- function(a,delta,y,x) {
    if (a==1) {
      return(0)
    } else {
      return((1/p_plac)*(omega_noS_n(x,y,delta)-Q_noS_n(t_0,x)+mean_Q_n))
    }
  }

  return(infl_fn)

}



#' Construct influence function corresponding to Kaplan-Meier estimator
#'
#' @param dat_combined Combined dataset; !!!!! TEMP
#' @return Influence function estimator
#' @note Used by est_med()
#' @noRd
construct_infl_fn_risk_v <- function(dat_v_rd, Q_n, g_n, omega_n, q_tilde_n,
                                     t_0, p_vacc) {

  dat_v2_rd <- dat_v_rd[dat_v_rd$z==1,]
  n_vacc <- attr(dat_v_rd, "n_vacc")
  dim_x <- attr(dat_v_rd, "dim_x")

  # Helper function
  C_n <- memoise2(function(s) {
    mean(apply(dat_v_rd, 1, function(r) {
      x <- as.numeric(r[1:dim_x])
      return(Q_n(t_0,x,s)*g_n(s,x))
    }))
  })

  mean_Q_n_g_n <- (1/n_vacc) * sum(apply(dat_v2_rd, 1, function(r) {
    r[["weights"]] * C_n(r[["s"]]) # Replace with sapply
  }))

  infl_fn <- function(a,z,weight,s,x,y,delta) {

    if (a==0) {

      return(0)

    } else {

      if (z==0) {
        s <- 0
        pi_ <- 1
      } else {
        pi_ <- 1/weight
      }

      val <- (z/pi_)*(omega_n(x,s,y,delta)-Q_n(t_0,x,s)) +
        (1-z/pi_)*q_tilde_n(x,y,delta) + mean_Q_n_g_n

      return((1/p_vacc)*val)

    }
  }

  return(memoise2(infl_fn))

}



#' Construct influence function corresponding to Kaplan-Meier estimator
#'
#' @param dat_combined Combined dataset; !!!!! TEMP
#' @return Influence function estimator
#' @note Used by est_med()
#' @noRd
construct_infl_fn_risk_v_v2 <- function(dat_v_rd, Q_noS_n_v, omega_noS_n_v, t_0,
                                     p_vacc) {

  dim_x <- attr(dat_v_rd, "dim_x")

  mean_Q_n <- mean(apply(dat_v_rd, 1, function(r) {
    Q_noS_n_v(t_0,as.numeric(r[1:dim_x]))
  }))

  infl_fn <- function(a,delta,y,x) {
    if (a==0) {
      return(0)
    } else {
      return((1/p_vacc)*(omega_noS_n_v(x,y,delta)-Q_noS_n_v(t_0,x)+mean_Q_n))
    }
  }

  return(infl_fn)

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

  return(construct_superfunc(fnc, vec=c(1,2,1,1,1,1)))

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
