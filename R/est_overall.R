#' Estimate overall risk and vaccine efficacy
#'
#' @description Estimate overall risk and vaccine efficacy.
#' @param dat A data object returned by load_data
#' @param t_0 Time point of interest
#' @param method One of c("KM", "Cox"), corresponding to either a Kaplan-Meier
#'     estimator ("KM") or a marginalized Cox proportional hazards model
#'     ("Cox").
#' @param risk Boolean. If TRUE, the controlled risk (CR) curve is computed.
#' @param ve Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
#'     computed.
#' @return A dataframe containing estimates
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' est_overall(dat=dat, t_0=578, method="KM")
#' @export
est_overall <- function(dat, t_0, method="Cox", risk=TRUE, ve=TRUE) { # ci_type="transformed"

  if (!(method %in% c("Cox", "KM")) || length(method)>1) {
    stop("`method` must equal either 'Cox' or 'KM'.")
  }

  if (!methods::is(dat,"vaccine_dat")) {
    stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by lo",
                "ad_data()."))
  }

  .groups <- attr(dat, "groups")
  if (.groups=="vaccine" && ve) {
    warning(paste0("Vaccine efficacy cannot be calculated because `dat` only c",
                   "ontains data from the vaccine group."))
  }
  if (.groups=="placebo" && ve) {
    warning(paste0("Vaccine efficacy cannot be calculated because `dat` only c",
                   "ontains data from the placebo group."))
  }

  if (.groups=="both") { .groups <- c("vaccine", "placebo") }

  # Container for results
  res <- data.frame(
    "stat" = character(),
    "group" = character(),
    "est" = double(),
    "se" = double(),
    "ci_lower" = double(),
    "ci_upper" = double()
  )

  for (grp in .groups) {

    if (grp=="vaccine") {
      dat_grp <- dat[dat$a==1,]
    } else {
      dat_grp <- dat[dat$a==0,]
    }

    if (method=="Cox") {

      # Alias random variables
      N <- length(dat_grp$s)
      dim_x <- attr(dat, "dim_x")
      X <- dat_grp[,c(1:dim_x), drop=F]
      V_ <- t(as.matrix(X))
      Y_ <- dat_grp$y
      D_ <- dat_grp$delta
      dim_v <- dim(V_)[1]

      # Create set of event times
      i_ev <- which(D_==1)

      # Fit an unweighted Cox model
      model <- survival::coxph(
        formula = stats::formula(paste0("survival::Surv(y,delta)~",
                                        paste(names(X),collapse="+"))),
        data = cbind(y=Y_, delta=D_, X)
      )
      coeffs <- model$coefficients
      beta_n <- as.numeric(coeffs)

      if (any(is.na(coeffs))) {
        na_coeffs <- names(coeffs)[which(is.na(coeffs))]
        stop(paste0("The following covariate coefficients were NA.: ",
                    paste(na_coeffs, collapse=","),
                    " Try removing these coefficients and rerunning."))
        # !!!!! Automatically omit from the model, as is done for spline coeffs
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
                (1/N) * sum(In(Y_>=x)*exp(LIN))
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
                (1/N) * as.numeric(V_ %*% (In(Y_>=x)*exp(LIN)))
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
                      res[i,j] <- (1/N)*sum(In(Y_>=x)*V_[i,]*V_[j,]*exp(LIN))
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

        m_n <- (function() {
          .cache <- new.env()
          function(x) {
            val <- .cache[[as.character(x)]]
            if (is.null(val)) {
              val <- S_1n(x) / S_0n(x)
              .cache[[as.character(x)]] <- val
            }
            return(val)
          }
        })()

      }

      # Estimated information matrix (for an individual)
      I_tilde <- Reduce("+", lapply(i_ev, function(i) {
        (S_2n(Y_[i])/S_0n(Y_[i])) - m_n(Y_[i]) %*% t(m_n(Y_[i]))
      }))
      I_tilde <- (1/N)*I_tilde
      I_tilde_inv <- solve(I_tilde)

      # Score function (Cox model)
      l_n <- (function() {
        .cache <- new.env()
        function(v_i,d_i,y_i) {
          key <- paste(c(v_i,d_i,y_i), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            explin_i <- exp(sum(v_i*beta_n))
            val <- d_i*(v_i-m_n(y_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
              (explin_i*In(Y_[j]<=y_i) * (v_i-m_n(Y_[j]))) / S_0n(Y_[j])
            }))
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # Influence function: beta_hat (regular Cox model)
      infl_fn_beta <- function(v_i,d_i,y_i) { I_tilde_inv %*% l_n(v_i,d_i,y_i) }

      # Nuisance constant: mu_n
      mu_n <- -1 * (1/N) * as.numeric(Reduce("+", lapply(i_ev, function(j) {
        (In(Y_[j]<=t_0) * m_n(Y_[j])) / S_0n(Y_[j])
      })))

      # Breslow estimator
      Lambda_n <- function(t) {
        (1/N) * sum(unlist(lapply(i_ev, function(i) {
          In(Y_[i]<=t) / S_0n(Y_[i])
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
        function(v_i,d_i,y_i) {
          key <- paste(c(v_i,d_i,y_i), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            val <- (function(v_i,d_i,y_i) {
              pc_5 <- sum(mu_n*infl_fn_beta(v_i,d_i,y_i))
              pc_1 <- ( d_i * In(y_i<=t_0) ) / S_0n(y_i)
              explin_i <- exp(sum(beta_n*v_i))
              pc_3 <- (1/N) * sum(unlist(lapply(i_ev, function(j) {
                (In(Y_[j]<=t_0)*In(y_i>=Y_[j])*explin_i) /
                  (S_0n(Y_[j]))^2
              })))
              return(pc_1+pc_5-pc_3)
            })(v_i,d_i,y_i)
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # Compute marginalized risk
      res_cox <- list()
      Lambda_n_t_0 <- Lambda_n(t_0)
      res_cox$est_marg <- (1/N) * sum((apply(dat_grp, 1, function(r) {
        x_i <- as.numeric(r[1:dim_x])
        exp(-1*exp(sum(beta_n*x_i))*Lambda_n_t_0)
      })))

      # Compute variance estimate
      res_cox$var_est_marg <- (function() {

        # Precalculate pieces dependent on s
        K_n <- (1/N) * Reduce("+", apply2(dat_grp, 1, function(r) {
          x_i <- as.numeric(r[1:dim_x])
          Q <- Q_n(x_i)
          explin <- exp(sum(x_i*beta_n))
          K_n1 <- Q
          K_n2 <- Q * explin
          K_n3 <- Q * explin * x_i
          return(c(K_n1,K_n2,K_n3))
        }, simplify=F))
        K_n1 <- K_n[1]
        K_n2 <- K_n[2]
        K_n3 <- K_n[3:length(K_n)]

        (1/N^2) * sum((apply(dat_grp, 1, function(r) {

          x_i <- as.numeric(r[1:dim_x])
          z_i <- r[["z"]]
          d_i <- r[["delta"]]
          y_i <- r[["y"]]

          pc_1 <- Q_n(x_i)
          pc_2 <- Lambda_n_t_0 * sum(K_n3*infl_fn_beta(x_i,d_i,y_i))
          pc_3 <- K_n2 * infl_fn_Lambda(x_i,d_i,y_i)
          pc_4 <- K_n1

          return((pc_1-pc_2-pc_3-pc_4)^2)

        })))

      })()

      # Extract estimates and SEs
      est <- 1-res_cox$est_marg
      se <- sqrt(res_cox$var_est_marg)

      # Generate confidence limits
      # if (ci_type=="none") {
        # ci_lo <- NA
        # ci_up <- NA
      # } else if (ci_type=="regular") {
        # ci_lo <- pmin(pmax(est - 1.96*se, 0), 1)
        # ci_up <- pmin(pmax(est + 1.96*se, 0), 1)
      # } else if (ci_type=="transformed") {
        ci_lo <- expit(logit(est) - 1.96*deriv_logit(est)*se)
        ci_up <- expit(logit(est) + 1.96*deriv_logit(est)*se)
      # }

    }

    if (method=="KM") {

      srv_p <- survival::survfit(survival::Surv(y,delta)~1, data=dat_grp)
      est <- 1 - srv_p$surv[which.min(abs(srv_p$time-t_0))]
      ci_lo <- 1 - srv_p$upper[which.min(abs(srv_p$time-t_0))]
      ci_up <- 1 - srv_p$lower[which.min(abs(srv_p$time-t_0))]
      se <- srv_p$std.err[which.min(abs(srv_p$time-t_0))]

    }

    # Update results dataframe
    res[nrow(res)+1,] <- list(stat="risk", group=grp, est=est, se=se,
                              ci_lower=ci_lo, ci_upper=ci_up)

  }

  if (ve) {

    est_v <- res[res$group=="vaccine","est"]
    est_p <- res[res$group=="placebo","est"]
    est_ve <- 1-(est_v/est_p)

    # CIs on the log(1-VE) scale
    se_v <- res[res$group=="vaccine","se"]
    se_p <- res[res$group=="placebo","se"]
    ci_lo <- 1 - exp(
      log(est_v/est_p) + 1.96*sqrt(se_v^2/est_v^2 + se_p^2/est_p^2)
    )
    ci_up <- 1 - exp(
      log(est_v/est_p) - 1.96*sqrt(se_v^2/est_v^2 + se_p^2/est_p^2)
    )

    # Standard error on the VE scale
    se_ve <- sqrt( se_v^2/est_p^2 + (se_p^2*est_v^2)/est_p^4 )

    # Update results dataframe
    res[nrow(res)+1,] <- list(stat="ve", group="both", est=est_ve, se=se_ve,
                              ci_lower=ci_lo, ci_upper=ci_up)

  }

  class(res) <- c(class(res), "vaccine_overall")

  return(res)

}
