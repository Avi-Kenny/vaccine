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
construct_km_infl_fn <- function(dat_combined, which) {

  num_v <- sum(dat_combined$a==1)
  num_p <- sum(dat_combined$a==0)
  PA1 <- num_v/(num_v+num_p)

  if (which=="placebo") {
    dat2 <- dat_combined[dat_combined$a==0,]
    prob <- 1-PA1
  } else if (which=="vaccine") {
    dat2 <- dat_combined[dat_combined$a==1,]
    prob <- PA1
  } else (
    stop("`which` must equal one of c('vaccine','placebo')")
  )

  Q_n <- construct_km(dat2)
  F_n <- stats::ecdf(dat2$y)
  v <- sort(unique(dat2$y))

  infl_fn <- function(A_i,Delta_i,Y_i,t) {
    if (which=="placebo" && A_i==1) {
      return(0)
    } else if (which=="vaccine" && A_i==0) {
      return(0)
    } else {
      v <- v[v<=min(t,Y_i)]
      if (Delta_i==1 && Y_i<=t) {
        term_1 <- 1 / (1-F_n(Y_i))
      } else {
        term_1 <- 0
      }
      term_2 <- 0
      for (v_j in v) {
        num <- sum(dat2$delta*In(dat2$y==v_j))
        if (num!=0) {
          den <- sum(In(v_j<=dat2$y))
          term_2 <- term_2 + log(1-(num/den)) / (1-F_n(v_j))
        }
      }
      return(-1*(1/prob)*Q_n(t)*(term_1+term_2))
    }
  }

  return(memoise2(infl_fn))

}
