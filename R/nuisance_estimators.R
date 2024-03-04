# !!!!! TO DO: update docs throughout

#' Construct conditional survival estimator Q_n
#'
#' @param type One of c("true", "Cox", "Random Forest", "Super Learner")
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param vals List of values to pre-compute function on; REQUIRED FOR SUPERLEARNER
#' @param return_model Logical; if TRUE, return the model object instead of the
#'     function
#' @return Conditional density estimator function
#'
#' @noRd
construct_Q_n <- function(type, dat_v, vals, return_model=F) {

  dim_x <- attr(dat_v,"dim_x")

  if (type=="Cox") {

    if (F) { weights <- NA } # Prevents CRAN error
    x_names <- names(dat_v)[1:dim_x]
    dat_v$delta2 <- 1 - dat_v$delta

    model_srv <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta)~",
                               paste(x_names,collapse="+"),"+s")),
      data = dat_v,
      weights = weights
    )
    coeffs_srv <- as.numeric(model_srv$coefficients)
    bh_srv <- survival::basehaz(model_srv, centered=F)

    model_cens <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta2)~",
                               paste(x_names,collapse="+"),"+s")),
      data = dat_v,
      weights = weights
    )
    coeffs_cens <- as.numeric(model_cens$coefficients)
    bh_cens <- survival::basehaz(model_cens, centered=F)

    fnc_srv <- function(t, x, s) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_srv$hazard[which.min(abs(bh_srv$time-t))]
        return(exp(-1*Lambda_t*exp(sum(coeffs_srv*as.numeric(c(x,s))))))
      }
    }

    fnc_cens <- function(t, x, s) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_cens$hazard[which.min(abs(bh_cens$time-t))]
        return(exp(-1*Lambda_t*exp(sum(coeffs_cens*as.numeric(c(x,s))))))
      }
    }

    rm(model_srv)
    rm(model_cens)

  } else {
    do.call("library", list("SuperLearner"))
  }

  if (type=="survSL") {

    # Excluding "survSL.rfsrc" for now. survSL.pchSL gives errors.
    methods <- c("survSL.coxph", "survSL.expreg", "survSL.km",
                 "survSL.loglogreg", "survSL.pchreg", "survSL.weibreg")

    # Prevents CRAN Note
    if (F) {
      # x <- earth::earth
      # x <- glmnet::glmnet
      x <- e1071::svm
    }

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)

    # Temporary (until survSuperLearner is on CRAN)
    survSuperLearner <- function() {}
    rm(survSuperLearner)
    tryCatch(
      expr = { do.call("library", list("survSuperLearner")) },
      error = function(e) {
        stop(paste0(
          "To use surv_type='survSL', you must install the `survSuperLearner` ",
          "package from github, using:\n\ndevtools::install_github(repo='tedwe",
          "stling/survSuperLearner')"))
      }
    )

    srv <- suppressWarnings(survSuperLearner(
      time = dat_v$y,
      event = dat_v$delta,
      X = dat_v[,c(1:dim_x,which(names(dat_v)=="s"))],
      newX = newX,
      new.times = new.times,
      event.SL.library = methods,
      cens.SL.library = methods,
      obsWeights = dat_v$weights,
      control = list(initWeightAlg=methods[1], max.SL.iter=10)
    ))

    srv_pred <- srv$event.SL.predict
    cens_pred <- srv$cens.SL.predict
    rm(srv)

    fnc_srv <- function(t, x, s) {
      row <- find_index(c(x,s), newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      row <- find_index(c(x,s), newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      # Note: the max() function is to prevent unreasonably small censoring
      #       probabilities from destabilizing estimates
      return(max(cens_pred[row,col], 0.001))
    }

  }

  if (type %in% c("survML-G","survML-L")) {

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)

    survML_args <- list(
      time = dat_v$y,
      event = dat_v$delta,
      X = dat_v[,c(1:dim_x,which(names(dat_v)=="s"))],
      newX = newX,
      newtimes = new.times,
      bin_size = 0.05,
      time_basis = "continuous",
      SL_control = list(
        # SL.library = rep(c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),2), # Note: rep() is to avoid a SuperLearner bug
        SL.library = rep(c("SL.mean", "SL.glm", "SL.gam"),2), # Note: rep() is to avoid a SuperLearner bug
        V = 5,
        obsWeights = dat_v$weights
      )
    )
    if (type=="survML-G") {

      fit <- suppressWarnings(do.call(survML::stackG, survML_args))
      srv_pred <- fit$S_T_preds
      cens_pred <- fit$S_C_preds

    } else if (type=="survML-L") {

      survML_args2 <- survML_args
      survML_args2$event <- round(1 - survML_args2$event)
      fit_s <- suppressWarnings(do.call(survML::stackL, survML_args))
      fit_c <- suppressWarnings(do.call(survML::stackL, survML_args2))
      srv_pred <- fit_s$S_T_preds
      cens_pred <- fit_c$S_T_preds

    }

    fnc_srv <- function(t, x, s) {
      row <- find_index(c(x,s), newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      row <- find_index(c(x,s), newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      # Note: the max() function is to prevent unreasonably small censoring
      #       probabilities from destabilizing estimates
      return(max(cens_pred[row,col], 0.001))
    }

  }

  if (type=="Cox") {
    sfnc_srv <- fnc_srv
    sfnc_cens <- fnc_cens
  } else {
    sfnc_srv <- memoise2(fnc_srv)
    sfnc_cens <- memoise2(fnc_cens)
  }
  rm("vals", envir=environment(get("fnc_srv",envir=environment(sfnc_srv))))

  return(list(srv=sfnc_srv, cens=sfnc_cens))

}



#' Construct conditional survival estimator Q_n
#'
#' @param type One of c("true", "Cox", "Random Forest", "Super Learner")
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param vals List of values to pre-compute function on; REQUIRED FOR SUPERLEARNER
#' @param return_model Logical; if TRUE, return the model object instead of the
#'     function
#' @return Conditional density estimator function
#'
#' @noRd
construct_Q_noS_n <- function(type, dat, vals, return_model=F) {

  dim_x <- attr(dat,"dim_x")

  # if (type=="Cox") {
  if (type %in% c("Cox", "survML-G", "survML-L")) { # !!!!! Temporary to avoid R-CMD-CHECK failure

    x_names <- names(dat)[1:dim_x]
    dat$delta2 <- 1 - dat$delta

    model_srv <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta)~",
                                      paste(x_names,collapse="+"))),
      data = dat
    )
    coeffs_srv <- as.numeric(model_srv$coefficients)
    bh_srv <- survival::basehaz(model_srv, centered=F)

    model_cens <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta2)~",
                                      paste(x_names,collapse="+"))),
      data = dat
    )
    coeffs_cens <- as.numeric(model_cens$coefficients)
    bh_cens <- survival::basehaz(model_cens, centered=F)

    fnc_srv <- function(t, x) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_srv$hazard[which.min(abs(bh_srv$time-t))]
        return(exp(-1*Lambda_t*exp(sum(coeffs_srv*as.numeric(x)))))
      }
    }

    fnc_cens <- function(t, x) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_cens$hazard[which.min(abs(bh_cens$time-t))]
        return(exp(-1*Lambda_t*exp(sum(coeffs_cens*as.numeric(x)))))
      }
    }

    rm(model_srv)
    rm(model_cens)

  } else {
    do.call("library", list("SuperLearner"))
  }

  if (type=="survSL") {

    # Excluding "survSL.rfsrc" for now. survSL.pchSL gives errors.
    methods <- c("survSL.coxph", "survSL.expreg", "survSL.km",
                 "survSL.loglogreg", "survSL.pchreg", "survSL.weibreg")

    newX <- vals$x[which(vals$t==0),, drop=F]
    new.times <- unique(vals$t)

    # Temporary (until survSuperLearner is on CRAN)
    survSuperLearner <- function() {}
    rm(survSuperLearner)
    tryCatch(
      expr = { do.call("library", list("survSuperLearner")) },
      error = function(e) {
        stop(paste0(
          "To use surv_type='survSL', you must install the `survSuperLearner` ",
          "package from github, using:\n\ndevtools::install_github(repo='tedwe",
          "stling/survSuperLearner')"))
      }
    )

    srv <- suppressWarnings(survSuperLearner(
      time = dat$y,
      event = dat$delta,
      X = dat[,c(1:dim_x), drop=F],
      newX = newX,
      new.times = new.times,
      event.SL.library = methods,
      cens.SL.library = methods,
      control = list(initWeightAlg=methods[1], max.SL.iter=10)
    ))

    srv_pred <- srv$event.SL.predict
    cens_pred <- srv$cens.SL.predict
    rm(srv)

    fnc_srv <- function(t, x) {
      row <- find_index(x, newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x) {
      row <- find_index(x, newX)
      col <- which.min(abs(t-new.times))
      if (length(col)!=1) { stop("length(col)!=1") }
      # Note: the max() function is to prevent unreasonably small censoring
      #       probabilities from destabilizing estimates
      return(max(cens_pred[row,col], 0.001))
    }

  }

  if (type=="Cox") {
    sfnc_srv <- fnc_srv
    sfnc_cens <- fnc_cens
  } else {
    sfnc_srv <- memoise2(fnc_srv)
    sfnc_cens <- memoise2(fnc_cens)
  }
  rm("vals", envir=environment(get("fnc_srv",envir=environment(sfnc_srv))))

  return(list(srv=sfnc_srv, cens=sfnc_cens))

}



#' Construct estimator of nuisance influence function omega_n
#'
#' @param Q_n Conditional survival function estimator returned by
#'     `construct_Q_n`
#' @param Qc_n Conditional censoring survival function estimator returned by
#'     `construct_Q_n`
#' @param type Defaults to "estimated". Override with "true" for debugging. Note
#'     that type="true" only works for surv_true="Cox" and assumes that Q_0
#'     and Qc_0 (i.e. the true functions) are passed in.
#' @return Estimator function of nuisance omega_0
#' @noRd
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
    index <- which.min(abs(k-grid$y)) - 1
    if (index==0) {
      return(0)
    } else {
      return(sum(x_vals(x,s)[1:index]*y_vals(x,s)[1:index]))
    }
  }
  # omega_integral <- memoise2(omega_integral) # !!!!! Does this make a difference (for profiling)?

  fnc <- function(x,s,y,delta) {
    Q_n(t_0,x,s) * (
      (delta * In(y<=t_0)) / (Q_n(y,x,s) * Qc_n(y,x,s)) -
        omega_integral(min(y,t_0),x,s)
    )
  }

  return(memoise2(fnc))

}



#' Construct estimator of nuisance influence function omega_n
#'
#' @param Q_noS_n Conditional survival function estimator returned by
#'     `construct_Q_noS_n` (no biomarker)
#' @param Qc_noS_n Conditional censoring survival function estimator returned by
#'     `construct_Q_noS_n` (no biomarker)
#' @param type Defaults to "estimated". Override with "true" for debugging. Note
#'     that type="true" only works for surv_true="Cox" and assumes that Q_0
#'     and Qc_0 (i.e. the true functions) are passed in.
#' @return Estimator function of nuisance omega_0
#' @noRd
construct_omega_noS_n <- function(Q_noS_n, Qc_noS_n, t_0, grid) {

  Q_n <- Q_noS_n
  Qc_n <- Qc_noS_n

  # Construct cumulative hazard estimator
  H_n <- function(t,x) { -1 * log(Q_n(t,x)) }

  # Construct memoised integral function
  x_vals <- memoise2(function(x) {
    diff(sapply(grid$y, function(t) { H_n(t,x) }))
  })
  y_vals <- memoise2(function(x) {
    sapply(grid$y[2:length(grid$y)], function(t) {
      ( Q_n(t,x) * Qc_n(t,x) )^-1
    })
  })
  omega_integral <- function(k,x) {
    index <- which.min(abs(k-grid$y)) - 1
    if (index==0) {
      return(0)
    } else {
      return(sum(x_vals(x)[1:index]*y_vals(x)[1:index]))
    }
  }

  fnc <- function(x,y,delta) {
    Q_n(t_0,x) * (
      (delta * In(y<=t_0)) / (Q_n(y,x) * Qc_n(y,x)) -
        omega_integral(min(y,t_0),x)
    )
  }

  return(memoise2(fnc))

}



#' Construct estimator of conditional density of S given X
#'
#' @param dat_v2 Dataset returned by `load_data`, subset to ph2 vacc
#' @param type One of c("parametric", "binning")
#' @param k Number of bins for the binning estimator (if k=0, then the number of
#'     bins will be selected via cross-validation); ignored for the parametric
#'     estimator
#' @param z1 Compute the density conditional on Z=1
#' @return Conditional density estimator function
#' @note
#'   - Assumes support of S is [0,1]
#' @noRd
construct_f_sIx_n <- function(dat_v2, type, k=0, z1=F) {

  dim_x <- attr(dat_v2, "dim_x")
  n_vacc2 <- attr(dat_v2, "n_vacc2")

  if (z1) { dat_v2$weights <- rep(1, length(dat_v2$weights)) }

  if (type=="parametric") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( truncnorm::dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

    # Set up weighted likelihood
    wlik <- function(prm) {
      -1 * sum(dat_v2$weights * apply(dat_v2, 1, function(r) {
        log(pmax(dens_s(s=r[["s"]], x=as.numeric(r[1:dim_x]), prm), 1e-8))
      }))
    }

    # Run optimizer
    opt <- Rsolnp::solnp(
      pars = c(rep(0, dim_x+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp() did not converge")
    }
    prm <- opt$pars

    # Remove large intermediate objects
    rm(dat_v2,dens_s,opt)

    fnc <- function(s, x) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( truncnorm::dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

  }

  if (type=="parametric (edge) 2") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[round(length(x)+2)]) # !!!!! Try using exp instead of expit
      return(truncnorm::dtruncnorm(s, a=0.01, b=1, mean=mu, sd=sigma))
    }

    # Probability P(S=s|X=x) for s==0 or s!=0
    prob_s <- function(s, x, prm) {
      s <- ifelse(s==0, 0, 1)
      lin <- sum(c(1,x)*prm)
      return(expit(lin)^(1-s)*(1-expit(lin))^s)
    }

    # Set up weighted likelihood (edge)
    wlik_1 <- function(prm) {
      -1 * sum(dat_v2$weights * apply(dat_v2, 1, function(r) {
        log(pmax(prob_s(s=r[["s"]], x=as.numeric(r[1:dim_x]), prm), 1e-8))
      }))
    }

    # Run optimizer (edge)
    opt_1 <- Rsolnp::solnp(
      pars = rep(0.01, dim_x+1),
      fun = wlik_1
    )
    if (opt_1$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp did not converge")
    }
    prm_1 <- opt_1$pars

    # Filter out observations with s==0
    dat_1 <- dat_v2[dat_v2$s!=0,]

    # Set up weighted likelihood (Normal)
    wlik_2 <- function(prm) {
      -1 * sum(dat_1$weights * apply(dat_1, 1, function(r) {
        log(pmax(dens_s(s=r[["s"]], x=as.numeric(r[c(1:dim_x)]), prm), 1e-8))
      }))
    }

    # Run optimizer (Normal)
    opt_2 <- Rsolnp::solnp(
      pars = c(rep(0.01, dim_x+1), 0.15),
      fun = wlik_2
    )
    if (opt_2$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp() did not converge")
    }
    prm_2 <- opt_2$pars

    # Remove large intermediate objects
    rm(dat_v2,dat_1,opt_1,opt_2)

    fnc <- function(s, x) {
      if (s==0) {
        return(prob_s(s=0, x, prm_1) / 0.01)
      } else {
        return(prob_s(s=1, x, prm_1) * dens_s(s, x, prm_2))
      }
    }

  }

  if (type=="binning") {

    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # prm[1] through prm[k-1] are the hazard components for the bins 1 to k-1
    # prm[k] and prm[k+1] are the coefficients for W1 and W2
    create_dens <- function(k, dat_v2, remove_dat=F) {

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
        -1 * sum(dat_v2$weights * log(pmax(apply(dat_v2, 1, function(r) {
          dens_s(s=r[["s"]], x=as.numeric(r[1:dim_x]), prm)
        }),1e-8)))
      }

      # Run optimizer
      if(.Platform$OS.type=="unix") { sink("/dev/null") } else { sink("NUL") }
      opt <- Rsolnp::solnp(
        pars = rep(0.001,k+dim_x-1),
        fun = wlik
      )
      sink()
      if (opt$convergence!=0) {
        warning("f_sIx_n: Rsolnp::solnp() did not converge")
      }
      prm <- opt$pars

      # Remove large intermediate objects
      rm(dens_s,opt)

      if (remove_dat) { rm(dat_v2) }

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
      folds <- sample(cut(c(1:n_vacc2), breaks=n_folds, labels=FALSE))
      ks <- c(4,8,12,16) # !!!!! Make this an argument
      best <- list(k=999, max_log_lik=999)

      # Cross-validation
      for (k in ks) {

        sum_log_lik <- 0
        for (i in c(1:n_folds)) {
          dat_train <- dat_v2[which(folds!=i),]
          dat_test <- dat_v2[which(folds==i),]
          dens <- create_dens(k, dat_train)
          sum_log_lik <- sum_log_lik + sum(log(apply(dat_test, 1, function(r) {
            dens(r[["s"]], as.numeric(r[c(1:dim_x)]))
          })))
        }

        if (sum_log_lik>best$max_log_lik || best$max_log_lik==999) {
          best$k <- k
          best$max_log_lik <- sum_log_lik
        }

      }

      k <- best$k
      message(paste("Number of bins selected (f_sIx_n):", k)) # !!!!! make this configurable with verbose option

    }

    fnc <- create_dens(k, dat_v2, remove_dat=T)

  }

  return(memoise2(fnc))

}



#' Construct estimator of marginal density of S
#'
#' @param dat_orig Dataset returned by `generate_data`
#' @param f_sIx_n A conditional density estimator returned by
#'     `construct_f_sIx_n`
#' @return Marginal density estimator function
#' @noRd
construct_f_s_n <- function(dat_v, f_sIx_n) {

  n_vacc <- attr(dat_v, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")

  if (attr(dat_v, "covariates_ph2")) {
    dat_v2 <- dat_v[dat_v$z==1,]
    datx_v2 <- dat_v2[, c(1:dim_x), drop=F]
    memoise2(function(s) {
      (1/n_vacc) * sum(dat_v2$weights * apply(datx_v2, 1, function(x) {
        f_sIx_n(s,as.numeric(x))
      }))
    })
  } else {
    datx_v <- dat_v[, c(1:dim_x), drop=F]
    memoise2(function(s) {
      (1/n_vacc) * sum(apply(datx_v, 1, function(x) {
        f_sIx_n(s,as.numeric(x))
      }))
    })
  }

}



#' Construct density ratio estimator g_n
#'
#' @param f_sIx_n Conditional density estimator returned by construct_f_sIx_n
#' @param f_s_n Marginal density estimator returned by construct_f_s_n
#' @return Density ratio estimator function
#' @noRd
construct_g_n <- function(f_sIx_n, f_s_n) {

  memoise2(function(s,x) { f_sIx_n(s,x) / f_s_n(s) })

}



#' Construct Phi_n
#'
#' @param dat_orig Dataset returned by `generate_data`
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param type One of c("step", "linear (mid)")
#' @return IPS-weighted CDF estimator function
#' @noRd
construct_Phi_n <- function (dat_v2, dat, type="linear (mid)") {

  # !!!!! type currently ignored

  n_vacc <- attr(dat_v2, "n_vacc")

  fnc <- memoise2(function(x) {
    (1/n_vacc) * sum(dat_v2$weights*In(dat_v2$s<=x))
  })

  return(memoise2(fnc))

}



#' Construct g-computation estimator function of theta_0
#'
#' @param dat_orig Dataset returned by `generate_data`
#' @param Q_n Conditional survival function estimator returned by
#'     `construct_Q_n`
#' @return G-computation estimator of theta_0
#' @noRd
construct_r_tilde_Mn <- function(dat_v, Q_n, t_0) {

  n_vacc <- attr(dat_v, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")
  if (attr(dat_v, "covariates_ph2")) {
    dat_v2 <- dat_v[dat_v$z==1,]
    datx_v2 <- dat_v2[, c(1:dim_x), drop=F]
    memoise2(function(s) {
      1 - (1/n_vacc) * sum(dat_v2$weights * apply(datx_v2, 1, function(x) {
        Q_n(t_0, as.numeric(x), s)
      }))
    })
  } else {
    datx_v <- dat_v[, c(1:dim_x), drop=F]
    memoise2(function(s) {
      1 - (1/n_vacc) * sum(apply(datx_v, 1, function(x) {
        Q_n(t_0, as.numeric(x), s)
      }))
    })
  }

}



#' Construct nuisance estimator of conditional density f(Y,Delta|X,S)
#'
#' @noRd
construct_f_n_srv <- function(Q_n, Qc_n, grid) {

  # Helper function to calculate derivatives
  construct_surv_deriv <- function(Q) {

    max_index <- length(grid$y)

    fnc <- function(t, x, s) {
      t_index <- which.min(abs(grid$y-t))
      if (t_index==1) {
        t1 <- 0
        t2 <- grid$y[2]
      } else if (t_index==max_index) {
        t1 <- grid$y[round(max_index-1)]
        t2 <- grid$y[max_index]
      } else {
        t1 <- grid$y[round(t_index-1)]
        t2 <- grid$y[round(t_index+1)]
      }
      return((Q(t2,x,s)-Q(t1,x,s))/(t2-t1))
    }

    return(memoise2(fnc))

  }

  # Calculate derivative estimators
  Q_n_deriv <- construct_surv_deriv(Q_n)
  Qc_n_deriv <- construct_surv_deriv(Qc_n)

  fnc <- function(y, delta, x, s) {
    if (delta==1) {
      return(-1*Qc_n(y,x,s)*Q_n_deriv(y,x,s))
    } else {
      return(-1*Q_n(y,x,s)*Qc_n_deriv(y,x,s))
    }
  }

  return(memoise2(fnc))

}



#' Construct q_n nuisance estimator function
#'
#' @noRd
construct_q_n <- function(type="standard", dat_v2, omega_n, g_n, r_tilde_Mn,
                          f_n_srv) {

  if (type=="standard") {

    # !!!!! Can this function be used elsewhere?
    q_n_star_inner <- memoise2(function(x,y,delta,s) {
      omega_n(x,s,y,delta)/g_n(s,x)+r_tilde_Mn(s)
    })

    # q_n_star <- memoise2(function(y,delta,x,s,u) {
    q_n_star <- function(y,delta,x,s,u) {
      if (s>u) {
        return(0)
      } else {
        return(q_n_star_inner(x,y,delta,s))
      }
    # })
    }

    # Helper functions
    # !!!!! Can these vectors/functions be used elsewhere?
    f_n_srv_s <- memoise2(function(y,delta,x) {
      sapply(dat_v2$s, function(s) { f_n_srv(y,delta,x,s) })
    })
    q_n_star_s <- memoise2(function(y,delta,x,u) {
      sapply(dat_v2$s, function(s) { q_n_star(y,delta,x,s,u) })
    })
    g_n_s <- memoise2(function(x) {
      sapply(dat_v2$s, function(s) { g_n(s,x) })
    })

    fnc <- function(x,y,delta,u) {

      f_n_srv_s <- f_n_srv_s(y,delta,x)
      q_n_star_s <- q_n_star_s(y,delta,x,u)
      g_n_s <- g_n_s(x)

      denom <- sum(dat_v2$weights*f_n_srv_s*g_n_s)
      if (denom==0) {
        return (0)
      } else {
        num <- sum(dat_v2$weights*q_n_star_s*f_n_srv_s*g_n_s)
        return(num/denom)
      }

    }

    return(memoise2(fnc))

  }

  if (type=="zero") {

    return(function(x, y, delta, u) { 0 })

  }

}



#' Construct estimator of nuisance g_sn
#'
#' @noRd
construct_g_sn <- function(dat_v2, f_n_srv, g_n, p_n) {

  n_vacc <- attr(dat_v2, "n_vacc")

  return(memoise2(function(x, y, delta) {
    num <- f_n_srv(y, delta, x, 0) * g_n(s=0,x) * (1-p_n)
    den <- (1/n_vacc) * sum(dat_v2$weights * sapply(dat_v2$s, function(s) {
      f_n_srv(y,delta,x,s) * g_n(s,x)
    }))
    if (den==0) { return(0) } else { return(num/den) }
  }))

}



#' Construct gamma_n nuisance estimator function
#'
#' @param dat_orig Dataset returned by `generate_data`
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param vals List of values to pre-compute function on
#' @param type Type of regression; one of c("cubic", "kernel", "kernel2")
#' @param omega_n A nuisance influence function returned by `construct_omega_n`
#' @param f_sIx_n A conditional density estimator returned by
#'     `construct_f_sIx_n`
#' @param f_sIx_n A conditional density estimator returned by
#'     `construct_f_sIx_n` among the observations for which z==1
#' @return gamma_n nuisance estimator function
#' @noRd
construct_gamma_n <- function(dat_v, type="Super Learner", omega_n,
                              grid) {

  dim_x <- attr(dat_v, "dim_x")
  dat_v2 <- dat_v[dat_v$z==1,]
  datx_v <- dat_v[, c(1:dim_x), drop=F]
  class(datx_v) <- "data.frame"

  # Construct pseudo-outcomes
  omega_n_d <- apply(dat_v2, 1, function(r) {
    omega_n(as.numeric(r[1:dim_x]), r[["s"]], r[["y"]], r[["delta"]])
  })
  dat_v2$po <- (dat_v2$weights*as.numeric(omega_n_d))^2

  # Filter out infinite values
  if (sum(!is.finite(dat_v2$po))!=0) {
    dat_v2 <- dat_v2[which(is.finite(dat_v2$po)),]
    warning(paste("gamma_n:", sum(!is.finite(dat_v2$po)),
                  "non-finite pseudo-outcome values"))
  }

  # Setup
  if (attr(dat_v, "covariates_ph2")) {
    datx_v2 <- dat_v2[, c(1:dim_x), drop=F]
    class(datx_v2) <- "data.frame"
    x_distinct <- dplyr::distinct(datx_v2)
  } else {
    x_distinct <- dplyr::distinct(datx_v)
  }
  x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
  newX <- expand.grid(x_index=x_distinct$x_index, s=grid$s)
  newX <- dplyr::inner_join(x_distinct, newX, by="x_index")
  newX$x_index <- NULL

  # Run regression
  if (type=="Super Learner") {

    # Fit SuperLearner regression
    do.call("library", list("SuperLearner"))
    # SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.loess",
    #                 "SL.nnet", "SL.ksvm", "SL.rpartPrune", "SL.svm")
    SL.library <- c("SL.mean", "SL.mean", "SL.gam", "SL.gam", "SL.ranger") # Changed on 2024-02-13; SL.mean written twice to avoid SuperLearner bug

    model_sl <- suppressWarnings(SuperLearner::SuperLearner(
      Y = dat_v2$po,
      X = dat_v2[,c(1:dim_x,which(names(dat_v2)=="s"))],
      newX = newX,
      family = "gaussian",
      SL.library = SL.library,
      verbose = F
    ))
    pred <- as.numeric(model_sl$SL.predict)
    if (sum(pred<0)!=0) {
      warning(paste("gamma_n:", sum(pred<0), "negative predicted values."))
    }

    # Construct regression function
    reg <- function(x,s) { pred[find_index(c(x,s), newX)] }

  }

  # Remove large intermediate objects
  rm(dat_v,dat_v2,datx_v,omega_n,model_sl)

  return(memoise2(reg))

}



#' Construct estimator of g_z0(x,s) = P(Z=1|X=x,S=s)
#'
#' @noRd
construct_g_zn <- function(dat_v, type="Super Learner", f_sIx_n,
                           f_sIx_z1_n) {

  # Prevents CRAN Note
  if (F) {
    x <- ranger::ranger
    x <- gam::gam
  }

  # Set library
  if (type=="Super Learner") {
    do.call("library", list("SuperLearner"))
    # SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
    #                 "SL.glmnet")
    # SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.nnet",
    #                 "SL.glmnet")
    SL.library <- c("SL.mean", "SL.mean", "SL.gam", "SL.gam", "SL.ranger") # Changed 2024-02-13; SL.mean written twice to avoid SuperLearner bug
  } else if (type=="logistic") {
    SL.library <- c("SL.glm")
  }

  # Create data objects
  dim_x <- attr(dat_v, "dim_x")
  dat_v2 <- dat_v[dat_v$z==1,]
  datx_v <- dat_v[, c(1:dim_x), drop=F]
  datx_v2 <- dat_v2[, c(1:dim_x), drop=F]
  class(datx_v) <- "data.frame"
  class(datx_v2) <- "data.frame"

  # Fit SuperLearner regression
  if (attr(dat_v, "covariates_ph2")) {
    newX <- dplyr::distinct(datx_v2)
    model_sl <- suppressWarnings(SuperLearner::SuperLearner(
      Y = dat_v2$z,
      X = datx_v2,
      newX = newX,
      family = "binomial",
      SL.library = SL.library,
      verbose = F
    ))
  } else {
    newX <- dplyr::distinct(datx_v)
    model_sl <- suppressWarnings(SuperLearner::SuperLearner(
      Y = dat_v$z,
      X = datx_v,
      newX = newX,
      family = "binomial",
      SL.library = SL.library,
      verbose = F
    ))
  }

  pred <- as.numeric(model_sl$SL.predict)

  # Construct regression function
  reg <- function(x) { pred[find_index(x, newX)] }

  # Remove large intermediate objects
  rm(dat_v,datx_v,model_sl)

  return(memoise2(function(x,s) { (f_sIx_z1_n(s,x)/f_sIx_n(s,x))*reg(x) }))

}



#' Construct derivative estimator r_Mn'_n
#'
#'
#' @param r_Mn An estimator of r_M0
#' @param type One of c("m-spline", "linear", "line")
#' @param dir One of c("decr", "incr")
#' @param grid TO DO
#' @noRd
construct_deriv_r_Mn <- function(type="m-spline", r_Mn, dir, grid) {

  if (type=="line") {

    fnc <- function(s) { r_Mn(1)-r_Mn(0) }

  } else {

    if (r_Mn(0)==r_Mn(1)) {
      fnc <- function(s) { 0 }
      warning("Estimated function is flat; variance estimation not possible.")
    }

    # Estimate entire function on grid
    r_Mns <- sapply(grid$s, r_Mn)

    # Compute set of midpoints of jump points (plus endpoints)
    points_x <- grid$s[1]
    points_y <- r_Mns[1]
    for (i in 2:length(grid$s)) {
      if (r_Mns[i]-r_Mns[round(i-1)]!=0) {
        points_x <- c(points_x, (grid$s[i]+grid$s[round(i-1)])/2)
        points_y <- c(points_y, mean(c(r_Mns[i],r_Mns[round(i-1)])))
      }
    }
    points_x <- c(points_x, grid$s[length(grid$s)])
    points_y <- c(points_y, r_Mns[length(grid$s)])

    if (type=="linear") {
      fnc_pre <- stats::approxfun(x=points_x, y=points_y, method="linear", rule=2)
    }

    if (type=="m-spline") {
      fnc_pre <- stats::splinefun(x=points_x, y=points_y, method="monoH.FC")
    }

    # Construct numerical derivative
    fnc <- function(s) {

      width <- 0.1 # !!!!! Changed from 0.2
      x1 <- s - width/2
      x2 <- s + width/2
      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }

      if (dir=="decr") {
        return(min((fnc_pre(x2)-fnc_pre(x1))/width,0))
      } else {
        return(max((fnc_pre(x2)-fnc_pre(x1))/width,0))
      }

    }

  }

  return(fnc)

}



#' Construct tau_n Chernoff scale factor function
#'
#' @param deriv_r_Mn A derivative estimator returned by `construct_deriv_r_Mn`
#' @param gamma_n Nuisance function estimator returned by `construct_gamma_n`
#' @param f_s_n Density estimator returned by `construct_f_s_n`
#' @return Chernoff scale factor estimator function
#' @noRd
construct_tau_n <- function(dat_v, deriv_r_Mn, gamma_n, f_sIx_n, g_zn) {

  n_vacc <- attr(dat_v, "n_vacc")
  dim_x <- attr(dat_v, "dim_x")
  dat_v2 <- dat_v[dat_v$z==1,]

  if (attr(dat_v, "covariates_ph2")) {
    dat_apply <- dat_v2
  } else {
    dat_apply <- dat_v
  }

  return(function(u) {
    abs(4*deriv_r_Mn(u) * (1/n_vacc) * sum(apply(dat_apply, 1, function(r) {
      x <- as.numeric(r[1:dim_x])
      return((gamma_n(x,u)*g_zn(x,u))/f_sIx_n(u,x))
    })))^(1/3)
  })

}



#' Construct nuisance estimator eta*_n
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @return Estimator function of nuisance eta*_0
#' @noRd
construct_etastar_n <- function(Q_n, t_0, grid) {

  int_values <- memoise2(function(x) {
    sapply(grid$s, function(s) { Q_n(t_0, x, s) })
  })

  fnc <- function(u,x) {
    indices <- which(grid$s<=u)
    # s_seq <- vals$s[indices]
    if (u==0 || length(indices)==0) {
      return(0)
    } else {
      # integral <- u * mean(sapply(s_seq, function(s) { Q_n(t_0, x, s) }))
      integral <- u * mean(int_values(x)[indices])
      return(u-integral)
    }
  }

  return(memoise2(fnc))

}



#' Construct q_tilde_n nuisance estimator function
#'
#' @param x x
#' @return q_tilde_n nuisance estimator function
#' @noRd
construct_q_tilde_n <- function(type="new", f_n_srv, f_sIx_n, omega_n) {

  if (type=="new") {

    # !!!!! TO DO

    # seq_01 <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    #
    # fnc <- function(x, y, delta, u) {
    #
    #   denom <- C$appx$s * sum(sapply(seq_01, function(s) {
    #     f_n_srv(y, delta, x, s) * f_sIx_n(s,x)
    #   }))
    #
    #   if (denom==0) {
    #     return (0)
    #   } else {
    #     if (u==0) {
    #       return(0)
    #     } else {
    #       seq_0u <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
    #       num <- C$appx$s * sum(sapply(seq_0u, function(s) {
    #         omega_n(x,s,y,delta) * f_n_srv(y, delta, x, s)
    #       }))
    #       return(num/denom)
    #     }
    #   }
    #
    # }

  }

  if (type=="zero") {

    fnc <- function(x, y, delta, u) { 0 }
    return(fnc) # !!!!! TEMP

  }

  return(memoise2(fnc))

}



