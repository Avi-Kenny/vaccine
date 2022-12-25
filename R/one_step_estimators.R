#' Construct Gamma_os_n primitive one-step estimator (based on EIF)
#'
construct_Gamma_os_n <- function(dat, dat_orig, omega_n, g_n, eta_n, p_n, q_n,
                                 r_tilde_Mn, Gamma_tilde_n) {

  n_orig <- attr(dat_orig, "n_orig")
  dim_x <- attr(dat_orig, "dim_x")
  piece_1 <- In(dat$s!=0)
  dat_df <- cbind(dat$x, s=dat$s, y=dat$y, delta=dat$delta) # !!!!! Possibly use this elsewhere
  dat_orig_df <- cbind(dat_orig$x, s=dat_orig$s, y=dat_orig$y,
                       delta=dat_orig$delta) # !!!!! Possibly use this elsewhere
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



