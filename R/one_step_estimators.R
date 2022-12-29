#' Construct Gamma_os_n primitive one-step estimator (based on EIF)
#'
#' @noRd
construct_Gamma_os_n <- function(dat, dat_orig, omega_n, g_n, eta_n, p_n, q_n,
                                 r_tilde_Mn, Gamma_tilde_n) {

  n_orig <- attr(dat_orig, "n_orig")
  dim_x <- attr(dat_orig, "dim_x")
  piece_1 <- In(dat$s!=0)

  dat_df <- as_df(dat)
  dat_orig_df <- as_df(dat_orig)
  piece_2 <- as.numeric(apply(dat_df, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    s <- r[["s"]]
    y <- r[["y"]]
    delta <- r[["delta"]]
    return((omega_n(x,s,y,delta)/g_n(s,x))+r_tilde_Mn(s))
  }))
  piece_3 <- (1-dat_orig$weights)

  # Remove large intermediate objects
  rm(omega_n,g_n,r_tilde_Mn)

  fnc <- function(u) {

    q_n_do <- as.numeric(apply(dat_orig_df, 1, function(r) {
      q_n(as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]], u)
    }))
    eta_n_xo <- apply(dat_orig$x, 1, function(x) { eta_n(u,as.numeric(x)) })

    return(
      (1/(n_orig*p_n)) * sum(dat$weights * (
        piece_1*In(dat$s<=u)*piece_2 - piece_1*Gamma_tilde_n(u)
      )) +
        (1/n_orig) * sum(piece_3*(q_n_do/p_n)+eta_n_xo)
    )

  }

  return(memoise2(fnc))

}



#' Compute one-step estimator of counterfactual survival at S=0
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param g_sn Propensity score estimator returned by construct_g_sn()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @return Value of one-step estimator
#' @noRd
r_Mn_edge <- function(dat_orig, g_sn, g_n, p_n, Q_n, omega_n, t_0) {

  n_orig <- attr(dat_orig, "n_orig")
  dim_x <- attr(dat_orig, "dim_x")
  dat_orig_df <- as_df(dat_orig)

  v <- (1/n_orig) * sum(apply(dat_orig_df, 1, function(r) {

    y <- r[["y"]]
    delta <- r[["delta"]]
    z <- r[["z"]]
    x <- as.numeric(r[1:dim_x])
    if (z==0) {
      s <- 0
      pi_ <- 1
    } else {
      s <- r[["s"]]
      pi_ <- 1/r[["weights"]]
    }

    return(
      ((z*In(s==0)+g_sn(x,y,delta)*(pi_-z))*omega_n(x,s=0,y,delta)) /
        (pi_*(1-p_n)*g_n(s=0,x)) - Q_n(t_0,x,s=0)
    )

  }))

  return(1+v)

}



