#' Estimate CVE/CR using Cox model
#'
#' @description Estimate controlled vaccine efficacy (CVE) and/or controlled
#'     risk (CR) using a marginalized Cox proportional hazards model
#' @param dat A data object returned by load_data
#' @param t_0 Time point of interest
#' @param cve Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
#'     computed.
#' @param cr Boolean. If TRUE, the controlled risk (CR) curve is computed.
#' @param s_out A numeric vector of s-values (on the biomarker scale) for which
#'     cve(s) and/or cr(s) are computed. Defaults to a grid of 101 points
#'     between the min and max biomarker values.
#' @param ci_type TO DO
#' @param grid_size A list containing the following three keys: \itemize{
#'     \item{\code{y}: grid size for time values}
#'     \item{\code{s}: grid size for marker values}
#'     \item{\code{x}: grid size for covariate values}
#' }
#'     This controls rounding of data values. Decreasing the grid size values
#'     results in shorter computation times, and increasing the grid size values
#'     results in more precise estimates. If grid_size$s=101, this means that a
#'     grid of 101 equally-spaced points (defining 100 intervals) will be
#'     created from min(S) to max(S), and each S value will be rounded to the
#'     nearest grid point. For grid_size$y, a grid will be created from 0 to
#'     t_0, and then extended to max(Y). For grid_size$x, a separate grid is
#'     created for each covariate column (binary/categorical covariates are
#'     ignored).
#' @param return_extras Boolean. If set to TRUE, the following quantities (most
#'     of which are mainly useful for debugging) are returned: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @param verbose A Boolean. If set to TRUE, intermediate output will be
#'     displayed.
#' @return A list containing the following: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @examples
#' print("to do")
#' @export
#' @references Kenny, A., Gilbert P., and Carone, M. (2023). Inference for
#'     controlled risk and controlled vaccine efficacy curves using a
#'     marginalized Cox proportional hazards model
#' @references Gilbert P., Fong Y., Kenny A., and Carone, M. (2022). A
#'     Controlled Effects Approach to Assessing Immune Correlates of Protection.
#' @export
est_cox <- function(
    dat, t_0, cve=T, cr=T, s_out=seq(from=min(dat$v$s), to=max(dat$v$s), l=101),
    ci_type="logit", grid_size=list(y=101, s=101, x=5), return_extras=F,
    verbose=F, spline_df=NA, edge_ind=F # !!!!!
) {

  # !!!!! Testing
  if (F) {
    # dat=list(v=readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_200.rds"));
    dat=list(v=readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_400.rds"));
    # dat=list(v=readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_800.rds"));
    class(dat)="dat_vaccine";
    attr(dat$v, "n_orig") <- length(dat$v$z)
    attr(dat$v, "dim_x") <- 2
    t_0=200; cve=T; cr=T; s_out=round(seq(0,1,0.02),2); ci_type="logit";
    grid_size=list(y=101, s=101, x=5); return_extras=F; verbose=F; spline_df=4;
    source("R/misc_functions.R"); s_out=round(seq(0,1,0.02),2); edge_ind=F;
  }

  # !!!!! Implement this eventually
  # # Setup
  # if (verbose) {
  #   print(paste("Check 0 (start):", Sys.time()))
  #   pbapply::pboptions(type="txt", char="#", txt.width=40, style=3)
  # } else {
  #   pbapply::pboptions(type="none")
  # }
  # if (parallel) {
  #   n_cores <- parallel::detectCores() - 1
  #   cl <- parallel::makeCluster(n_cores)
  #   if (verbose) {
  #     print(cl)
  #   }
  # } else {
  #   cl <- NULL
  # }


  if (class(dat)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  # !!!!! Validate other inputs; import error handling function from SimEngine

  # Alias variables
  .v <- verbose

  # Fix s_out if needed
  if (missing(s_out)) {
    s_out <- seq(from=min(dat$v$s, na.rm=T),
                 to=max(dat$v$s, na.rm=T),
                 l=101)
  } else {
    if (any(is.na(s_out))) { stop("NA values not allowed in s_out.") }
  }

  # Create spline basis and function
  dat$v$spl <- data.frame("s1"=dat$v$s)
  if (!is.na(spline_df) && spline_df!=1) {

    # Create spline basis
    spl_basis <- splines::ns(
      x = dat$v$s,
      df = spline_df,
      intercept = F,
      Boundary.knots = quantile(dat$v$s, c(0.05,0.95), na.rm=T)
    )
    spl_knots <- as.numeric(attr(spl_basis, "knots"))
    spl_bnd_knots <- as.numeric(attr(spl_basis, "Boundary.knots"))

    for (i in c(1:(dim(spl_basis)[2]))) {
      dat$v$spl[[paste0("s",i)]] <- spl_basis[,i]
    }
    dim_s <- dim(spl_basis)[2]

    if (edge_ind) {

      dim_s <- dim_s + 1
      min_s <- min(dat$v$s, na.rm=T)
      dat$v$spl[,paste0("s",dim_s)] <- In(dat$v$s==min_s)
      s_to_spl <- function(s) {
        spl <- as.numeric(splines::ns(
          x = s,
          df = spline_df,
          knots = spl_knots,
          intercept = F,
          Boundary.knots = spl_bnd_knots
        ))
        return(c(spl, In(s==min_s)))
      }

    } else {

      s_to_spl <- function(s) {
        as.numeric(splines::ns(
          x = s,
          df = spline_df,
          knots = spl_knots,
          intercept = F,
          Boundary.knots = spl_bnd_knots
        ))
      }

    }

  } else {

    if (edge_ind) {

      dim_s <- 2
      min_s <- min(dat$v$s, na.rm=T)
      dat$v$spl[,paste0("s",dim_s)] <- In(dat$v$s==min_s)
      s_to_spl <- function(s) { c(s, In(s==min_s)) }

    } else {

      dim_s <- 1
      s_to_spl <- function(s) { s }

    }

  }

  # Create phase-two data object (unrounded)
  dat_v_ph2 <- ss(dat$v, which(dat$v$z==1))

  # Alias random variables
  WT <- dat_v_ph2$weights
  ST <- dat_v_ph2$strata
  N <- length(dat$v$s)
  n <- length(dat_v_ph2$s)
  X <- dat_v_ph2$x
  SP <- dat_v_ph2$spl
  V_ <- t(as.matrix(cbind(X,SP)))
  Y_ <- dat_v_ph2$y
  D_ <- dat_v_ph2$delta

  # Get dimensions
  dim_v <- dim(V_)[1]
  dim_x <- attr(dat$v, "dim_x")

  # Create set of event times
  i_ev <- which(D_==1)

  # Fit an IPS-weighted Cox model
  model <- survival::coxph(
    formula = formula(paste0("survival::Surv(y,delta)~",
                             paste(names(X),collapse="+"), "+",
                             paste(names(SP),collapse="+"))),
    data = cbind(y=Y_, delta=D_, X, SP),
    weights = WT
  )
  coeffs <- model$coefficients
  beta_n <- as.numeric(coeffs)

  if (any(is.na(coeffs))) {

    if (any(is.na(coeffs[which(substr(names(coeffs),1,1)=="x")]))) {
      print(summary(model))
      stop(paste0("Some covariate coefficients were NA. Try removing these coe",
                  "fficients."))
      # !!!!! Automatically omit from the model, as is done for spline coeffs
    }

    if (any(is.na(coeffs[which(substr(names(coeffs),1,1)=="s")]))) {
      print(summary(model))
      warning(paste0("Some spline coefficients were NA; these coefficients hav",
                     "e been automatically omitted from the model"))
      na_inds <- as.numeric(which(is.na(coeffs)))
      s_na_inds <- which(is.na(coeffs[which(substr(names(coeffs),1,1)=="s")]))
      s_na_inds <- as.numeric(s_na_inds)
      beta_n <- beta_n[-na_inds]
      SP <- SP[,-s_na_inds]
      names(SP) <- paste0("s", c(1:length(SP)))
      V_ <- t(as.matrix(cbind(X,SP)))
      dim_v <- dim(V_)[1]
      s_to_spl2 <- s_to_spl
      s_to_spl <- function(s) { s_to_spl2(s)[-s_na_inds] }

    }

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
            (1/N) * sum(WT*In(Y_>=x)*exp(LIN))
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
            (1/N) * as.numeric(V_ %*% (WT*In(Y_>=x)*exp(LIN)))
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
                  res[i,j] <- (1/N)*sum(WT*In(Y_>=x)*V_[i,]*V_[j,]*exp(LIN))
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

    m_n <- function(x) { S_1n(x) / S_0n(x) }

  }

  # Estimated information matrix (for an individual)
  I_tilde <- Reduce("+", lapply(i_ev, function(i) {
    WT[i] * ( (S_2n(Y_[i])/S_0n(Y_[i])) - m_n(Y_[i]) %*% t(m_n(Y_[i])) )
  }))
  I_tilde <- (1/N)*I_tilde
  I_tilde_inv <- solve(I_tilde)

  # Score function (Cox model)
  l_n <- function(v_i,d_i,y_i) {
    d_i*(v_i-m_n(y_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*exp(sum(v_i*beta_n))*In(Y_[j]<=y_i) * (v_i-m_n(Y_[j]))) /
        S_0n(Y_[j])
    }))
  }

  # Influence function: beta_hat (regular Cox model)
  l_tilde <- (function() {
    .cache <- new.env()
    function(v_i,d_i,y_i) {
      key <- paste(c(v_i,d_i,y_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,d_i,y_i) {
          I_tilde_inv %*% l_n(v_i,d_i,y_i)
        })(v_i,d_i,y_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Create p_n and p1_n vectors
  n_strata <- max(dat$v$strata)
  p_n <- c()
  p1_n <- c()
  for (c in c(1:n_strata)) {
    p_n[c] <- mean(c==dat$v$strata)
    p1_n[c] <- mean(c==dat$v$strata & dat$v$z==1)
  }

  # Influence function: beta_hat (est. weights)
  infl_fn_beta <- (function() {
    .cache <- new.env()

    exp_terms <- list()
    for (i in c(1:max(dat$v$strata))) {
      js <- which(ST==i)
      if (length(js)>0) {
        exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
          l_tilde(V_[,j],D_[j],Y_[j])
        }))
      } else {
        exp_terms[[i]] <- as.matrix(rep(0,dim_v)) # Avoid NAs in small samples
      }
    }

    function(v_i,z_i,d_i,y_i,wt_i,st_i) {
      key <- paste(c(v_i,z_i,d_i,y_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,z_i,d_i,y_i,wt_i,st_i) {

          if (z_i) {
            piece_1 <- wt_i*l_tilde(v_i,d_i,y_i)
          } else {
            piece_1 <- as.matrix(rep(0,dim_v))
          }

          if (p1_n[st_i]!=0) {
            piece_2 <- (1 - (z_i*p_n[st_i])/p1_n[st_i]) * exp_terms[[st_i]]
          } else {
            piece_2 <- as.matrix(rep(0,dim_v)) # Avoid NAs in small samples
          }

          return(as.numeric(piece_1+piece_2))

        })(v_i,z_i,d_i,y_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance constant: mu_n
  mu_n <- -1 * (1/N) * as.numeric(Reduce("+", lapply(i_ev, function(j) {
    (WT[j] * In(Y_[j]<=t_0) * m_n(Y_[j])) / S_0n(Y_[j])
  })))

  # Nuisance function: v_1n
  v_1n <- (function() {
    .cache <- new.env()
    function(st_i,z_i,y_j) {
      key <- paste(c(st_i,z_i,y_j), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,z_i,y_j) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/N) * (p1_n[st_i]-p_n[st_i]*z_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  In(Y_[k]>=y_j) * exp(sum(beta_n*V_[,k]))
                }))
              )
            )
          } else { return(0) }
        })(st_i,z_i,y_j)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Nuisance function: v_2n
  v_2n <- (function() {
    .cache <- new.env()
    function(st_i,z_i) {
      key <- paste(c(st_i,z_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,z_i) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/N) * (p1_n[st_i]-p_n[st_i]*z_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  ( D_[k] * In(Y_[k]<=t_0) ) / S_0n(Y_[k])
                }))
              )
            )
          } else { return(0) }
        })(st_i,z_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Breslow estimator
  Lambda_n <- function(t) {
    (1/N) * sum(unlist(lapply(i_ev, function(i) {
      WT[i] * ( In(Y_[i]<=t) / S_0n(Y_[i]) )
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
    function(v_i,z_i,d_i,y_i,wt_i,st_i) {
      key <- paste(c(v_i,z_i,d_i,y_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(v_i,z_i,d_i,y_i,wt_i,st_i) {
          pc_4 <- (1/N) * sum(unlist(lapply(i_ev, function(j) {
            ( WT[j] * In(Y_[j]<=t_0) * v_1n(st_i,z_i,Y_[j]) ) / (S_0n(Y_[j]))^2
          })))
          pc_5 <- sum(mu_n*infl_fn_beta(v_i,z_i,d_i,y_i,wt_i,st_i))
          pc_2 <- v_2n(st_i,z_i)

          if (z_i==1) {
            pc_1 <- ( wt_i * d_i * In(y_i<=t_0) ) / S_0n(y_i)
            pc_3 <- (1/N) * sum(unlist(lapply(i_ev, function(j) {
              (WT[j]*In(Y_[j]<=t_0)*wt_i*In(y_i>=Y_[j])*exp(sum(beta_n*v_i))) /
                (S_0n(Y_[j]))^2
            })))
            return(pc_1+pc_2+pc_5-pc_3-pc_4)
          } else {
            return(pc_2+pc_5-pc_4)
          }
        })(v_i,z_i,d_i,y_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # if (verbose) { print(paste("Check 2 (functions declared):", Sys.time())) }

  # Compute marginalized risk
  res_cox <- list()
  Lambda_n_t_0 <- Lambda_n(t_0)
  res_cox$est_marg <- unlist(lapply(s_out, function(s) {
    (1/N) * sum((apply(dat$v$x, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      s_spl <- s_to_spl(s)
      exp(-1*exp(sum(beta_n*c(x_i,s_spl)))*Lambda_n_t_0)
    })))
  }))

  # if (verbose) { print(paste("Check 3d (var est START: marg):", Sys.time())) }

  # Compute variance estimate
  dat_v_df <- as_df(dat$v, strata=T)
  res_cox$var_est_marg <- unlist(lapply(s_out, function(s) {

    # Precalculate pieces dependent on s
    s_spl <- s_to_spl(s)
    K_n <- (1/N) * Reduce("+", apply2(dat_v_df, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      Q <- Q_n(c(x_i,s_spl))
      explin <- exp(sum(c(x_i,s_spl)*beta_n))
      K_n1 <- Q
      K_n2 <- Q * explin
      K_n3 <- Q * explin * c(x_i,s_spl)
      return(c(K_n1,K_n2,K_n3))
    }, simplify=F))
    K_n1 <- K_n[1]
    K_n2 <- K_n[2]
    K_n3 <- K_n[3:length(K_n)]

    (1/N^2) * sum((apply(dat_v_df, 1, function(r) {

      x_i <- as.numeric(r[1:dim_x])
      if (!is.na(r[["s"]])) {
        s_i <- s_to_spl(r[["s"]])
      } else {
        s_i <- NA
      }
      z_i <- r[["z"]]
      d_i <- r[["delta"]]
      y_i <- r[["y"]]
      wt_i <- r[["weights"]]
      st_i <- r[["strata"]]

      pc_1 <- Q_n(c(x_i,s_spl))
      pc_2 <- Lambda_n_t_0 * sum(
        K_n3 * infl_fn_beta(c(x_i,s_i),z_i,d_i,y_i,wt_i,st_i)
      )
      pc_3 <- K_n2 * infl_fn_Lambda(c(x_i,s_i),z_i,d_i,y_i,wt_i,st_i)
      pc_4 <- K_n1

      return((pc_1-pc_2-pc_3-pc_4)^2)

    })))

  }))

  # if (verbose) { print(paste("Check 5d (var est END: marg):", Sys.time())) }

  # Extract estimates and SEs
  ests <- 1-res_cox$est_marg
  ses <- sqrt(res_cox$var_est_marg)

  # Generate confidence limits
  if (ci_type=="none") {
    ci_lo <- rep(NA, length(ests))
    ci_hi <- rep(NA, length(ests))
  } else if (ci_type=="regular") {
    ci_lo <- (ests - 1.96*ses) %>% pmax(0) %>% pmin(1)
    ci_hi <- (ests + 1.96*ses) %>% pmax(0) %>% pmin(1)
  } else if (ci_type=="logit") {
    ci_lo <- expit(logit(ests) - 1.96*deriv_logit(ests)*ses)
    ci_hi <- expit(logit(ests) + 1.96*deriv_logit(ests)*ses)
  }

  # Create results object
  res <- list(s=s_out, est=ests, ci_lo=ci_lo, ci_hi=ci_hi)

  # Return extras
  if (return_extras) { res$model <- model }

  return(res)

}
