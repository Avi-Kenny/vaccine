#' Construct Gamma_os_n primitive one-step estimator (based on EIF)
#'
#' @noRd
construct_Gamma_os_n <- function(dat_v, omega_n, g_n, q_n, r_tilde_Mn) {

  dat_v2 <- dat_v[dat_v$z==1,]
  n_vacc <- attr(dat_v, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")

  piece_1 <- as.numeric(apply(dat_v2, 1, function(r) {
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
  piece_2 <- 1 - dat_v$weights

  # Remove large intermediate objects
  rm(omega_n,g_n,r_tilde_Mn)

  fnc <- function(u) {

    q_n_do <- as.numeric(apply(dat_v, 1, function(r) {
      q_n(as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]], u)
    }))

    return((1/n_vacc) * sum(dat_v2$weights * (In(dat_v2$s<=u)*piece_1))
           + (1/n_vacc) * sum(piece_2*q_n_do))

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
r_Mn_edge <- function(dat_v_rd, g_sn, g_n, p_n, Q_n, omega_n, t_0) {

  n_vacc <- attr(dat_v_rd, "n_vacc")
  dim_x <- attr(dat_v_rd, "dim_x")

  v <- (1/n_vacc) * sum(apply(dat_v_rd, 1, function(r) {

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
risk_overall_np_v <- function(dat_v_rd, g_n, Q_n, omega_n, f_n_srv, q_tilde_n,
                              t_0) {

  # !!!!! Port to est_overall

  n_vacc <- attr(dat_v_rd, "n_vacc")
  dim_x <- attr(dat_v_rd, "dim_x")

  v <- (1/n_vacc) * sum(apply(dat_v_rd, 1, function(r) {

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

    val <- (z/pi_)*(omega_n(x,s,y,delta)-Q_n(t_0,x,s)) +
      (1-z/pi_)*q_tilde_n(x,y,delta)
    return(val)

  }))

  return(1+v)

}



#' Compute one-step estimator of overall risk (placebo group)
#'
#' @param TODO TO DO
#' @return Value of one-step estimator
#' @noRd
risk_overall_np_p <- function(dat_p_rd, Q_noS_n, omega_noS_n, t_0) {

  # !!!!! Port to est_overall

  n_plac <- attr(dat_p_rd, "n_plac")
  dim_x <- attr(dat_p_rd, "dim_x")

  v <- (1/n_plac) * sum(apply(dat_p_rd, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    omega_noS_n(x,r[["y"]],r[["delta"]])-Q_noS_n(t_0,x)
  }))

  return(1+v)

}



#' Compute one-step estimator of overall risk (placebo group)
#'
#' @param TODO TO DO
#' @return Value of one-step estimator
#' @noRd
risk_overall_np_v_v2 <- function(dat_v_rd, Q_noS_n_v, omega_noS_n_v, t_0) {

  # !!!!! Port to est_overall

  n_vacc <- attr(dat_v_rd, "n_vacc")
  dim_x <- attr(dat_v_rd, "dim_x")

  v <- (1/n_vacc) * sum(apply(dat_v_rd, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    omega_noS_n_v(x,r[["y"]],r[["delta"]])-Q_noS_n_v(t_0,x)
  }))

  return(1+v)

}



#' Construct Theta_os_n primitive one-step estimator
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_sIx_n Conditional density estimator returned by construct_f_sIx_n
#' @param etastar_n A nuisance estimator returned by construct_etastar_n()
#' @return Gamma_os_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
#' @noRd
construct_Theta_os_n <- function(dat_v, omega_n, f_sIx_n, q_tilde_n,
                                 etastar_n) {

  dat_v2 <- dat_v[dat_v$z==1,]
  n_vacc <- attr(dat_v, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")

  piece_1 <- as.numeric(apply(dat_v2, 1, function(r) {
    x <- as.numeric(r[1:dim_x])
    s <- r[["s"]]
    y <- r[["y"]]
    delta <- r[["delta"]]
    f_sIx_n_val <- f_sIx_n(s,x)
    if (is.nan(f_sIx_n_val)) {
      stop(paste0("One or more f_sIx_n_val values were NAN; density might be z",
                  "ero. Try using a density estimator that guarantees positive",
                  " density estimates."))
    }
    return((dat_v2$weights*omega_n(x,s,y,delta)) / f_sIx_n_val)
  }))

  # Remove large intermediate objects
  rm(omega_n,f_sIx_n)

  fnc <- function(u) {
    (1/n_vacc) * sum(piece_1*In(dat_v2$s<=u)) +
      (1/n_vacc) * sum(as.numeric(apply(dat_v, 1, function(r) {
        weights <- r[["weights"]]
        x <- as.numeric(r[1:dim_x])
        y <- r[["y"]]
        delta <- r[["delta"]]
        return((1-weights) * q_tilde_n(x,y,delta,u) + etastar_n(u,x))
      })))
  }

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

  } # !!!!! DEBUG

  return(memoise2(fnc))

}



