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
construct_infl_fn_risk_p <- function(dat_p_copy, Q_noS_n, omega_noS_n,
                                     dim_x, t_0, prob) {

  # !!!!! Temporary
  df_p <- cbind(dat_p_copy$x,
                y = dat_p_copy$y,
                delta = dat_p_copy$delta)

  mean_Q_n <- mean(apply(df_p, 1, function(r) {
    Q_noS_n(t_0,as.numeric(r[1:dim_x]))
  }))

  infl_fn <- function(a,delta,y,x) {
    if (a==1) {
      return(0)
    } else {
      return((1/prob)*(omega_noS_n(x,y,delta)-Q_noS_n(t_0,x)+mean_Q_n))
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
construct_infl_fn_risk_v <- function(dat_orig, Q_n, g_n, omega_n,
                                     q_tilde_n, t_0, prob) {

  n_orig <- attr(dat_orig, "n_orig")
  dim_x <- attr(dat_orig, "dim_x")
  dat_orig_df <- as_df(dat_orig)
  dat2 <- dat_orig_df[dat_orig_df$z==1,]

  # Helper function
  C_n <- memoise2(function(s) {
    mean(apply(dat_orig_df, 1, function(r) {
      x <- as.numeric(r[1:dim_x])
      return(Q_n(t_0,x,s)*g_n(s,x))
    }))
  })

  mean_Q_n_g_n <- (1/n_orig) * sum(apply(dat2, 1, function(r) {
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

      return((1/prob)*val)
    }
  }

  return(memoise2(infl_fn))

}
