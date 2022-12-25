# !!!!! TO DO: construct_Q_n



#' Construct estimator of nuisance influence function omega_n
#'
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param Qc_n Conditional censoring survival function estimator returned by
#'     construct_Q_n
#' @param type Defaults to "estimated". Override with "true" for debugging. Note
#'     that type="true" only works for surv_true="Cox PH" and assumes that Q_0
#'     and Qc_0 (i.e. the true functions) are passed in.
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(Q_n, Qc_n, t_0, grid) {

  # Construct cumulative hazard estimator
  H_n <- function(t,x,s) { -1 * log(Q_n(t,x,s)) }

  # Construct memoised integral function
  x_vals <- memoise2(function(x,s) {
    diff(sapply(grid$y, function(t) { H_n(t,x,s) }))
  })
  y_vals <- memoise2(function(x,s) {
    sapply(grid$y[2:length(grid$y)], function(t) {
      ( Q_n(t,x,s) * Qc_n(t,x,s) )^-1
    })
  })
  omega_integral <- function(k,x,s) {
    if (k<grid$y[2]) {
      integral <- 0
    } else {
      index <- max(which(k>=grid$y))
      integral <- sum(x_vals(x,s)[1:index]*y_vals(x,s)[1:index])
    }
  }
  # omega_integral <- memoise2(omega_integral) # !!!!! Maybe try un-memoising this

  fnc <- function(x,s,y,delta) {
    Q_n(t_0,x,s) * (
      (delta * In(y<=t_0)) / (Q_n(y,x,s) * Qc_n(y,x,s)) -
        omega_integral(min(y,t_0),x,s)
    )
  }

  return(memoise2(fnc))

}



#' Construct estimator of conditional density of S given X
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param type One of c("parametric", "binning")
#' @param k Number of bins for the binning estimator (if k=0, then the number of
#'     bins will be selected via cross-validation); ignored for the parametric
#'     estimator
#' @param z1 Compute the density conditional on Z=1
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of S is [0,1]
construct_f_sIx_n <- function(dat, type, k=0, z1=F) {

  if (z1) { dat$weights <- rep(1, length(dat$weights)) }

  if (type=="parametric") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

    # Set up weighted likelihood
    wlik <- function(prm) {
      -1 * sum(sapply(c(1:length(dat$s)), function(i) {
        dat$weights[i] * log(pmax(dens_s(s=dat$s[i], x=as.numeric(dat$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer
    opt <- solnp(
      pars = c(rep(0, length(dat$x)+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("construct_f_sIx_n: solnp() did not converge")
    }
    prm <- opt$pars

    # Remove large intermediate objects
    rm(dat,dens_s,opt)

    fnc <- function(s, x) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

  }

  if (type=="parametric (edge)") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( dtruncnorm(s, a=C$appx$s, b=1, mean=mu, sd=sigma) )
    }

    # Estimate p_0
    n_orig <- sum(dat$weights)
    p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))

    # Filter out observations with s==0
    dat <- ss(dat, which(dat$s!=0))

    # Set up weighted likelihood
    wlik <- function(prm) {
      -1 * sum(sapply(c(1:length(dat$s)), function(i) {
        dat$weights[i] *
          log(pmax(dens_s(s=dat$s[i], x=as.numeric(dat$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer
    opt <- solnp(
      pars = c(rep(0, length(dat$x)+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("construct_f_sIx_n: solnp() did not converge")
    }
    prm <- opt$pars

    # Remove large intermediate objects
    rm(dat,dens_s,opt)

    fnc <- function(s, x) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      pc_1 <- In(s==0) * ( (1-p_n) / C$appx$s )
      pc_2 <- In(s!=0) * p_n * dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma)
      return(pc_1+pc_2)
    }

  }

  if (type=="parametric (edge) 2") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[round(length(x)+2)]) # !!!!! Try using exp instead of expit
      return( dtruncnorm(s, a=C$appx$s, b=1, mean=mu, sd=sigma) )
    }

    # Probability P(S=s|X=x) for s==0 or s!=0
    prob_s <- function(s, x, prm) {
      s <- ifelse(s==0, 0, 1)
      lin <- sum(c(1,x)*prm)
      return(expit(lin)^(1-s)*(1-expit(lin))^s)
    }

    # Set up weighted likelihood (edge)
    wlik_1 <- function(prm) {
      -1 * sum(sapply(c(1:length(dat$s)), function(i) {
        dat$weights[i] *
          log(pmax(prob_s(s=dat$s[i], x=as.numeric(dat$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer (edge)
    opt_1 <- solnp(
      pars = rep(0.01, length(dat$x)+1),
      fun = wlik_1
    )
    if (opt_1$convergence!=0) {
      warning("construct_f_sIx_n: opt_1 did not converge")
    }
    prm_1 <- opt_1$pars
    print("prm_1") # !!!!!
    print(prm_1) # !!!!!

    # Filter out observations with s==0
    dat_1 <- ss(dat, which(dat$s!=0))

    # Set up weighted likelihood (Normal)
    wlik_2 <- function(prm) {
      -1 * sum(sapply(c(1:length(dat_1$s)), function(i) {
        dat_1$weights[i] *
          log(pmax(dens_s(s=dat_1$s[i], x=as.numeric(dat_1$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer (Normal)
    opt_2 <- solnp(
      pars = c(rep(0.01, length(dat_1$x)+1), 0.15),
      fun = wlik_2
    )
    if (opt_2$convergence!=0) {
      warning("construct_f_sIx_n: solnp() did not converge")
    }
    prm_2 <- opt_2$pars
    print("prm_2") # !!!!!
    print(prm_2) # !!!!!

    # Remove large intermediate objects
    rm(dat,dat_1,opt_1,opt_2)

    fnc <- function(s, x) {
      if (s==0) {
        return(prob_s(s=0, x, prm_1) / C$appx$s)
      } else {
        return(prob_s(s=1, x, prm_1) * dens_s(s, x, prm_2))
      }
    }

  }

  if (type=="binning") {

    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # prm[1] through prm[k-1] are the hazard components for the bins 1 to k-1
    # prm[k] and prm[k+1] are the coefficients for W1 and W2
    create_dens <- function(k, dat, remove_dat=F) {

      # Cut points
      alphas <- seq(0, 1, length.out=k+1)

      # Density for a single observation
      dens_s <- memoise2(function(s, x, prm) {
        bin <- ifelse(s==1, k, which.min(s>=alphas)-1)
        hz <- expit(
          prm[c(1:(ifelse(bin==k,k-1,bin)))] +
            as.numeric(prm[c(k:(k+length(x)-1))] %*% x)
        )
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        return(k*p1*p2)
      })

      # Set up weighted likelihood
      wlik <- function(prm) {
        -1 * sum(dat$weights * log(pmax(dens_s(s=dat$s, x=dat$x, prm),1e-8)))
      }

      # Run optimizer
      opt <- solnp(
        pars = rep(0.001,k+length(dat$x)-1),
        fun = wlik
      )
      if (opt$convergence!=0) {
        warning("construct_f_sIx_n: solnp() did not converge")
      }
      prm <- opt$pars

      # Remove large intermediate objects
      rm(dens_s,opt)

      if (remove_dat) { rm(dat) }

      fnc <- function(s, x) {

        bin <- ifelse(s==1, k, which.min(s>=alphas)-1)
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(prm[j] + prm[c(k:(k+length(x)-1))] %*% x)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))

        return(k*p1*p2)

      }

    }

    # Select k via cross-validation
    if (k==0) {

      # Prep
      n_folds <- 5
      folds <- sample(cut(c(1:length(dat$s)), breaks=n_folds, labels=FALSE))
      ks <- c(5,10,15,20) # !!!!! Make this an argument
      best <- list(k=999, max_log_lik=999)

      # Cross-validation
      for (k in ks) {

        sum_log_lik <- 0
        for (i in c(1:n_folds)) {
          dat_train <- list(
            s = dat$s[-which(folds==i)],
            weights = dat$weights[-which(folds==i)],
            x = dat$x[-which(folds==i),]
          )
          dat_test <- list(
            s = dat$s[which(folds==i)],
            weights = dat$weights[which(folds==i)],
            x = dat$x[which(folds==i),]
          )
          dens <- create_dens(k, dat_train)
          sum_log_lik <- sum_log_lik +
            sum(log(sapply(c(1:length(dat_test$s)), function(j) {
              dens(dat_test$s[j], as.numeric(dat_test$x[j,]))
            })))
        }

        if (sum_log_lik>best$max_log_lik || best$max_log_lik==999) {
          best$k <- k
          best$max_log_lik <- sum_log_lik
        }

      }

      k <- best$k
      print(paste("Number of bins selected (f_sIx_n):", k))

    }

    fnc <- create_dens(k, dat, remove_dat=T)

  }

  return(memoise2(fnc))

}
