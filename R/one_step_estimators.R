#' Construct Gamma_os_n primitive one-step estimator (based on EIF)
#'
#' @noRd
construct_Gamma_os_n <- function(dat, dat_orig, omega_n, g_n, q_n, r_tilde_Mn) {

  n_orig <- attr(dat_orig, "n_orig")
  dim_x <- attr(dat_orig, "dim_x")

  dat_df <- as_df(dat)
  dat_orig_df <- as_df(dat_orig)
  piece_1 <- as.numeric(apply(dat_df, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    s <- r[["s"]]
    y <- r[["y"]]
    delta <- r[["delta"]]
    g_n_val <- g_n(s,x)
    if (is.nan(g_n_val)) {
      stop(paste0("One or more g_n values were NAN; density might be zero. Try ",
                  "using a density estimator that guarantees positive density ",
                  "estimates."))
    }
    return((omega_n(x,s,y,delta)/g_n_val)+r_tilde_Mn(s))
  }))
  piece_2 <- (1-dat_orig$weights)

  # Remove large intermediate objects
  rm(omega_n,g_n,r_tilde_Mn)

  fnc <- function(u) {

    q_n_do <- as.numeric(apply(dat_orig_df, 1, function(r) {
      q_n(as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]], u)
    }))

    return((1/n_orig) * sum(dat$weights * (In(dat$s<=u)*piece_1))
           +(1/n_orig) * sum(piece_2*q_n_do))

  }

  return(memoise2(fnc))

}



#' Compute one-step estimator of counterfactual risk at S=0
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



#' Compute one-step estimator of overall risk (vaccine group)
#'
#' @param TODO TO DO
#' @return Value of one-step estimator
#' @noRd
risk_overall_np_v <- function(dat_orig, g_n, Q_n, omega_n, f_n_srv, q_tilde_n,
                              t_0) {

  # !!!!! Port to est_overall

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

    val <- (z/pi_)*(omega_n(x,s,y,delta)-Q_n(t_0,x,s)) + # !!!!! New
      (1-z/pi_)*q_tilde_n(x,y,delta) # !!!!! New
    return(val)

  }))

  return(1+v)

}



#' Compute one-step estimator of overall risk (placebo group)
#'
#' @param TODO TO DO
#' @return Value of one-step estimator
#' @noRd
risk_overall_np_p <- function(dat_orig_rounded_p, dim_x, Q_noS_n, omega_noS_n, t_0) {

  # !!!!! Temporary
  df_p <- cbind(dat_orig_rounded_p$x,
                y = dat_orig_rounded_p$y,
                delta = dat_orig_rounded_p$delta)

  # !!!!! Port to est_overall

  n_orig <- nrow(df_p)
  # dat_orig_df <- as_df(dat_orig)

  v <- (1/n_orig) * sum(apply(df_p, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    omega_noS_n(x,r[["y"]],r[["delta"]])-Q_noS_n(t_0,x)
  }))

  return(1+v)

}
