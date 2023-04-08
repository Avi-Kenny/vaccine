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
    verbose=F
) {

  # !!!!! Testing
  if (F) {
    dat=list(v=readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_200.rds"));
    class(dat)="dat_vaccine";
    t_0=200; cve=T; cr=T; s_out=round(seq(0,1,0.02),2); ci_type="logit";
    grid_size=list(y=101, s=101, x=5); return_extras=F; verbose=F;
    source("R/misc_functions.R");
    s_out=round(seq(0,1,0.2),2) # !!!!!
  }

  if (class(dat)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  # !!!!! Validate other inputs; import error handling function from SimEngine

  # Alias variables
  dat_orig <- dat$v # !!!!! Maybe change this later
  .v <- verbose

  # Fix s_out if needed
  if (missing(s_out)) {
    s_out <- seq(from=min(dat_orig$s, na.rm=T),
                 to=max(dat_orig$s, na.rm=T),
                 l=101)
  } else {
    if (any(is.na(s_out))) { stop("NA values not allowed in s_out.") }
  }

  # Set params
  p <- list("ci_type"=ci_type)

  # Rescale S to lie in [0,1] and create rounded data object
  # !!!!! Maybe don't need to round or rescale
  s_min <- min(dat_orig$s, na.rm=T)
  s_max <- max(dat_orig$s, na.rm=T)
  s_shift <- -1 * s_min
  s_scale <- 1/(s_max-s_min)
  dat_orig$s <- (dat_orig$s+s_shift)*s_scale
  grid <- create_grid(dat_orig, grid_size, t_0)
  dat_orig_rounded <- round_dat(dat_orig, grid, grid_size)

  # Rescale/round s_out and remove s_out points outside [0,1]
  # !!!!! Maybe don't need to round or rescale
  s_out_orig <- s_out
  s_out <- (s_out+s_shift)*s_scale
  s_out <- sapply(s_out, function(s) { grid$s[which.min(abs(grid$s-s))] })
  na_head <- sum(s_out<0)
  na_tail <- sum(s_out>1)
  if (na_head>0) { s_out <- s_out[-c(1:na_head)] }
  len_p <- length(s_out)
  if (na_tail>0) { s_out <- s_out[-c((len_p-na_tail+1):len_p)] }

  # Create phase-two data object (unrounded)
  dat <- ss(dat_orig, which(dat_orig$z==1))

  # Fit Cox model and compute variance
  {

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

      # Fit an IPS-weighted Cox model
      model <- survival::coxph(
        formula = formula(paste0("survival::Surv(y,delta)~",
                                 paste(names(dat$x),collapse="+"),"+s")),
        data = cbind(y=dat$y, delta=dat$delta, dat$x, s=dat$s),
        weights = dat$weights
        # Note: scaling the weights affects the SEs but not the estimates; thus, this is only needed for debugging
        # weights = dat$weights * (length(dat$weights)/sum(dat$weights))
      )
      beta_n <- as.numeric(model$coefficients)

      # Alias random variables
      WT <- dat$weights
      ST <- dat$strata
      N <- length(dat_orig$s)
      n <- length(dat$s)
      Z_ <- t(as.matrix(cbind(dat$x,s=dat$s)))
      T_ <- dat$y
      Ds_ <- dat$delta
      lin <- as.numeric(t(beta_n)%*%Z_)
      d <- dim(Z_)[1]

      # Intermediate functions
      {
        S_2n <- (function() {
          .cache <- new.env()
          function(x) {
            val <- .cache[[as.character(x)]]
            if (is.null(val)) {
              val <- (function(x) {
                res <- matrix(NA, nrow=d, ncol=d)
                for (i in c(1:d)) {
                  for (j in c(1:d)) {
                    if (!is.na(res[j,i])) {
                      res[i,j] <- res[j,i]
                    } else {
                      res[i,j] <- (1/N)*sum(WT*In(T_>=x)*Z_[i,]*Z_[j,]*exp(lin))
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

        S_0n <- (function() {
          .cache <- new.env()
          function(x) {
            val <- .cache[[as.character(x)]]
            if (is.null(val)) {
              val <- (function(x) {
                (1/N) * sum(WT*In(T_>=x)*exp(lin))
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
                (1/N) * as.numeric(Z_ %*% (WT*In(T_>=x)*exp(lin)))
              })(x)
              .cache[[as.character(x)]] <- val
            }
            return(val)
          }
        })()

        m_n <- function(x) { S_1n(x) / S_0n(x) }

        h <- function(x) { (S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x)) }

      }

      # Create set of event times
      i_ev <- which(Ds_==1)
      t_ev <- T_[which(Ds_==1)]

      # Estimated information matrix (for an individual)
      I_tilde <- Reduce("+", lapply(i_ev, function(i) {
        WT[i] * h(T_[i])
      }))
      I_tilde <- (1/N)*I_tilde
      I_tilde_inv <- solve(I_tilde)

      # Score function (Cox model)
      l_star <- function(z_i,ds_i,t_i) {
        ds_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
          (WT[j]*exp(sum(z_i*beta_n))*In(T_[j]<=t_i)*(z_i-m_n(T_[j]))) /
            S_0n(T_[j])
        }))
      }

      # Influence function: beta_hat (regular Cox model)
      l_tilde <- (function() {
        .cache <- new.env()
        function(z_i,ds_i,t_i) {
          key <- paste(c(z_i,ds_i,t_i), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            val <- (function(z_i,ds_i,t_i) {
              I_tilde_inv %*% l_star(z_i,ds_i,t_i)
            })(z_i,ds_i,t_i)
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # Create p_n and p1_n vectors
      n_strata <- max(dat_orig$strata)
      p_n <- c()
      p1_n <- c()
      for (c in c(1:n_strata)) {
        p_n[c] <- mean(c==dat_orig$strata)
        p1_n[c] <- mean(c==dat_orig$strata & dat_orig$z==1)
      }

      # Influence function: beta_hat (est. weights)
      infl_fn_beta <- (function() {
        .cache <- new.env()

        exp_terms <- list()
        for (i in c(1:max(dat_orig$strata))) {
          js <- which(ST==i)
          if (length(js)>0) {
            exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
              l_tilde(Z_[,j],Ds_[j],T_[j])
            }))
          } else {
            exp_terms[[i]] <- as.matrix(rep(0,d)) # Avoid NAs in small samples
          }
        }

        function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
              if (d_i) {
                piece_1 <- wt_i*l_tilde(z_i,ds_i,t_i)
              } else {
                piece_1 <- as.matrix(rep(0,d))
              }

              if (p1_n[st_i]!=0) {
                piece_2 <- (1 - (d_i*p_n[st_i])/p1_n[st_i]) * exp_terms[[st_i]]
              } else {
                piece_2 <- as.matrix(rep(0,d)) # Avoid NAs in small samples
              }

              return(as.numeric(piece_1+piece_2))
            })(z_i,d_i,ds_i,t_i,wt_i,st_i)
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # Nuisance constant: mu_n
      mu_n <- (1/N) * as.numeric(Reduce("+", lapply(t_ev, function(t_j) {
        (In(t_j<=t_0) * S_1n(t_j)) / (S_0n(t_j))^2
      })))

      # Nuisance function: v_n
      v_n <- (function() {
        .cache <- new.env()
        function(st_i,d_i,t_j) {
          key <- paste(c(st_i,d_i,t_j), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            val <- (function(st_i,d_i,t_j) {
              k_set <- which(ST==st_i)
              if (length(k_set)>0) {
                return(
                  (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
                    unlist(lapply(k_set, function(k) {
                      In(T_[k]>=t_j) * exp(sum(beta_n*Z_[,k]))
                    }))
                  )
                )
              } else { return(0) }
            })(st_i,d_i,t_j)
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # Breslow estimator
      Lambda_n <- function(t) {
        (1/N) * sum(unlist(lapply(i_ev, function(i) {
          In(T_[i]<=t) / S_0n(T_[i])
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
        function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
          val <- .cache[[key]]
          if (is.null(val)) {
            val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
              pc_4 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
                ( In(t_j<=t_0) * v_n(st_i,d_i,t_j) ) / (S_0n(t_j))^2
              })))
              pc_5 <- sum(mu_n*infl_fn_beta(z_i,d_i,ds_i,t_i,wt_i,st_i))

              if (d_i==1) {
                pc_1 <- ( ds_i * In(t_i<=t_0) ) / S_0n(t_i)
                pc_3 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
                  (In(t_j<=t_0)*wt_i*In(t_i>=t_j)*exp(sum(beta_n*z_i))) /
                    (S_0n(t_j))^2
                })))
                return(pc_1-pc_3-pc_4-pc_5)
              } else {
                return(-1*(pc_4+pc_5))
              }
            })(z_i,d_i,ds_i,t_i,wt_i,st_i)
            .cache[[key]] <- val
          }
          return(val)
        }
      })()

      # memoise2 <- function(fnc) {
      #
      #   htab <- new.env()
      #   ..new_fnc <- function() {
      #     ..e <- parent.env(environment())
      #     ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
      #     key <- rlang::hash(..mc)
      #     val <- ..e$htab[[key]]
      #     if (is.null(val)) {
      #       val <- do.call(..e$fnc, ..mc)
      #       ..e$htab[[key]] <- val
      #     }
      #     return(val)
      #   }
      #
      #   # Set formals and set up environment
      #   formals(..new_fnc) <- formals(fnc)
      #   f_env <- new.env(parent=environment(fnc))
      #   f_env$arg_names <- names(formals(fnc))
      #   f_env$htab <- htab
      #   f_env$fnc <- fnc
      #   environment(..new_fnc) <- f_env
      #
      #   return(..new_fnc)
      #
      # }

      # # Influence function: survival at a point (est. weights)
      # omega_n <- (function() {
      #   .cache <- new.env()
      #   function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
      #     # ..count <<- ..count+1 # !!!!!
      #     # ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame()) # !!!!!
      #     # key <- rlang::hash(..mc) # !!!!!
      #     key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i,z), collapse=" ")
      #     val <- .cache[[key]]
      #     if (is.null(val)) {
      #       val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
      #         explin <- exp(sum(z*beta_n))
      #         piece_1 <- Lambda_n(t_0) * explin *
      #           (t(z)%*%infl_fn_beta(z_i,d_i,ds_i,t_i,wt_i,st_i))[1]
      #         piece_2 <- explin * infl_fn_Lambda(z_i,d_i,ds_i,t_i,wt_i,st_i)
      #         return(-1*Q_n(z)*(piece_1+piece_2))
      #       })(z_i,d_i,ds_i,t_i,wt_i,st_i,z)
      #       .cache[[key]] <- val
      #     }
      #     return(val)
      #   }
      # })()

      # # Influence function: marginalized survival (est. weights)
      # infl_fn_marg <- (function() {
      #   .cache <- new.env()
      #   x_ <- as.list(as.data.frame(t(dat_orig$x)))
      #   function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
      #     key <- paste(c(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s), collapse=" ")
      #     val <- .cache[[key]]
      #     if (is.null(val)) {
      #       val <- (function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
      #         piece_1 <- Q_n(c(x_i,s))
      #         piece_2 <- (1/N) * sum(unlist(lapply(c(1:N), function(j) {
      #           x_j <- x_[[j]]
      #           omgn <- omega_n(c(x_i,s_i),d_i,ds_i,t_i,wt_i,st_i,c(x_j,s))
      #           return(omgn-Q_n(c(x_j,s)))
      #         })))
      #         return(piece_1+piece_2)
      #       })(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s)
      #       .cache[[key]] <- val
      #     }
      #     return(val)
      #   }
      # })()

      # if (verbose) { print(paste("Check 2 (functions declared):", Sys.time())) }

      # !!!! basehaz vs. Lambda_hat
      res_cox <- list()
      # res_cox <- list(model=model, beta_n=beta_n)
      dim_x <- attr(dat_orig, "dim_x")
      bh <- survival::basehaz(model, centered=FALSE)
      index <- max(which((bh$time<t_0)==T))
      est_bshz <- bh$hazard[index]
      N <- sum(dat$weights)
      res_cox$est_marg <- unlist(lapply(s_out, function(s) {
        (1/N) * sum((apply(dat_orig$x, 1, function(r) {
          x_i <- as.numeric(r[1:dim_x])
          exp(-1*exp(sum(beta_n*c(as.numeric(x_i),s)))*est_bshz)
        })))
      }))

      # if (verbose) { print(paste("Check 3d (var est START: marg):", Sys.time())) }
      dat_orig_df <- as_df(dat_orig, strata=T)

      # Pre-calculate Lambda_n(t_0)
      Lambda_n_t_0 <- Lambda_n(t_0)

      # Calculate variance
      # !!!!! Not using parallelization or progress bar for now
      res_cox$var_est_marg <- unlist(lapply(s_out, function(s) {
      # if (verbose) { print(paste0("Check 4d (point=",s,"): ", Sys.time())) }

        # Precalculate pieces dependent on s
        K_n <- (1/N) * Reduce("+", apply(dat_orig_df, 1, function(r) {
          x_i <- as.numeric(r[1:dim_x])
          Q <- Q_n(c(x_i,s))
          explin <- exp(sum(c(x_i,s)*beta_n))
          K_n1 <- Q
          K_n2 <- Q * explin
          K_n3 <- Q * explin * c(x_i,s)
          return(c(K_n1,K_n2,K_n3))
        }, simplify=F))
        K_n1 <- K_n[1]
        K_n2 <- K_n[2]
        K_n3 <- K_n[3:length(K_n)]

        (1/N^2) * sum((apply(dat_orig_df, 1, function(r) {

          x_i <- as.numeric(r[1:dim_x])
          s_i <- r[["s"]]
          z_i <- r[["z"]] # d_i
          dl_i <- r[["delta"]] # ds_i
          y_i <- r[["y"]] # t_i
          wt_i <- r[["weights"]]
          st_i <- r[["strata"]]

          pc_1 <- Q_n(c(x_i,s))
          pc_2 <- Lambda_n_t_0 * sum(
            K_n3 * infl_fn_beta(c(x_i,s_i),z_i,dl_i,y_i,wt_i,st_i)
          )
          pc_3 <- K_n2 * infl_fn_Lambda(c(x_i,s_i),z_i,dl_i,y_i,wt_i,st_i)
          pc_4 <- K_n1

          return((pc_1-pc_2-pc_3-pc_4)^2)

        })))

      }))
      # if (verbose) { print(paste("Check 5d (var est END: marg):", Sys.time())) }

  }

  # Extract estimates and SEs
  ests <- 1-res_cox$est_marg
  ses <- sqrt(res_cox$var_est_marg)

  # Generate confidence limits
  if (p$ci_type=="none") {
    ci_lo <- rep(NA, length(ests))
    ci_hi <- rep(NA, length(ests))
  } else if (p$ci_type=="regular") {
    ci_lo <- (ests - 1.96*ses) %>% pmax(0) %>% pmin(1)
    ci_hi <- (ests + 1.96*ses) %>% pmax(0) %>% pmin(1)
  } else if (p$ci_type=="logit") {
    ci_lo <- expit(logit(ests) - 1.96*deriv_logit(ests)*ses)
    ci_hi <- expit(logit(ests) + 1.96*deriv_logit(ests)*ses)
  }

  # Create results object
  res <- list(
    s = s_out_orig,
    est = c(rep(NA,na_head), ests, rep(NA,na_tail)),
    ci_lo = c(rep(NA,na_head), ci_lo, rep(NA,na_tail)),
    ci_hi = c(rep(NA,na_head), ci_hi, rep(NA,na_tail))
  )

  return(res)

}
