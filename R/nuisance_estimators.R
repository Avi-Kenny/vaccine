# !!!!! TO DO: update docs throughout

#' Construct conditional survival estimator Q_n
#'
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param vals List of values to pre-compute function on; REQUIRED FOR SUPERLEARNER
#' @param type One of c("true", "Cox", "Random Forest", "Super Learner")
#' @param return_model Logical; if TRUE, return the model object instead of the
#'     function
#' @param print_coeffs Logical; if TRUE, print the algorithm coefficients for
#'     the conditional survival and censoring functions (applicable only if
#'     type=="Super Learner")
#' @return Conditional density estimator function
#'
#' @noRd
construct_Q_n <- function(type, dat, vals, return_model=F, print_coeffs=F) {

  if (type=="Cox") {

    model_srv <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta)~",
                               paste(names(dat$x),collapse="+"),"+s")),
      data = cbind(y=dat$y, delta=dat$delta, dat$x, s=dat$s),
      weights = dat$weights
    )
    coeffs_srv <- as.numeric(model_srv$coefficients)
    bh_srv <- survival::basehaz(model_srv, centered=F)

    model_cens <- survival::coxph(
      formula = stats::formula(paste0("survival::Surv(y,delta)~",
                               paste(names(dat$x),collapse="+"),"+s")),
      data = cbind(y=dat$y, delta=1-dat$delta, dat$x, s=dat$s),
      weights = dat$weights
    )
    coeffs_cens <- as.numeric(model_cens$coefficients)
    bh_cens <- survival::basehaz(model_cens, centered=F)

    fnc_srv <- function(t, x, s) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_srv$hazard[which.min(abs(bh_srv$time-t))]
        # Lambda_t <- bh_srv$hazard[max(which((bh_srv$time<t)==T))] # This version sometimes throws an error
        return(exp(-1*Lambda_t*exp(sum(coeffs_srv*as.numeric(c(x,s))))))
      }
    }

    fnc_cens <- function(t, x, s) {
      if (t==0) {
        return(1)
      } else {
        Lambda_t <- bh_cens$hazard[which.min(abs(bh_cens$time-t))]
        # Lambda_t <- bh_cens$hazard[max(which((bh_cens$time<t)==T))] # This version sometimes throws an error
        return(exp(-1*Lambda_t*exp(sum(coeffs_cens*as.numeric(c(x,s))))))
      }
    }

    rm(model_srv)
    rm(model_cens)

  }

  if (type=="Super Learner") { # "Super Learner (Westling)"

    # Excluding "survSL.rfsrc" for now. survSL.pchSL gives errors.
    methods <- c("survSL.coxph", "survSL.expreg", "survSL.km",
                 "survSL.loglogreg", "survSL.pchreg", "survSL.weibreg")
    # methods <- c("survSL.km", "survSL.pchreg", "survSL.rfsrc") # !!!!!

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)

    srv <- survSuperLearner::survSuperLearner(
      time = dat$y,
      event = dat$delta,
      X = cbind(dat$x, s=dat$s),
      newX = newX,
      new.times = new.times,
      event.SL.library = methods,
      cens.SL.library = methods,
      obsWeights = dat$weights,
      control = list(initWeightAlg=methods[1], max.SL.iter=10)
    )

    if (print_coeffs) {
      cat("\n------------------------------\n")
      cat("SuperLearner algorithm weights\n")
      cat("------------------------------\n\n")
      cat("event.coef\n")
      cat("----------\n")
      print(sort(srv$event.coef, decreasing=T))
      cat("\ncens.coef\n")
      cat("---------\n")
      print(sort(srv$cens.coef, decreasing=T))
      cat("\n------------------------------\n")
    }

    srv_pred <- srv$event.SL.predict
    cens_pred <- srv$cens.SL.predict
    rm(srv)

    # !!!!! Later consolidate these via a wrapper/constructor function
    fnc_srv <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (methods::is(newX[["s"]][1],"factor")) {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (methods::is(newX[["s"]][1],"factor")) {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(cens_pred[row,col])
    }

  }

  # if (type=="survML") {
  if (substr(type, 1, 6)=="survML") { # "Super Learner (Wolock)"

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)

    survML_args <- list(
      time = dat$y,
      event = dat$delta,
      X = cbind(dat$x, s=dat$s),
      newX = newX,
      newtimes = new.times,
      bin_size = 0.1, # !!!!! bin_size = 0.02
      time_basis = "continuous",
      SL_control = list(
        # SL.library = rep(c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),2), # Note: rep() is to avoid a SuperLearner bug
        SL.library = rep(c("SL.mean", "SL.glm", "SL.gam"),2), # Note: rep() is to avoid a SuperLearner bug
        V = 5,
        obsWeights = dat$weights
      )
    )
    if (type=="survML-G") {
      fit <- do.call(survML::stackG, survML_args)
      srv_pred <- fit$S_T_preds
      cens_pred <- fit$S_C_preds
    } else if (type=="survML-L") {
      survML_args2 <- survML_args
      survML_args2$event <- round(1 - survML_args2$event)
      fit_s <- do.call(survML::stackL, survML_args)
      fit_c <- do.call(survML::stackL, survML_args2)
      srv_pred <- fit_s$S_T_preds
      cens_pred <- fit_c$S_T_preds
    }

    if (print_coeffs) {
      cat("\n-------------------------------------\n")
      cat("survML SuperLearner algorithm weights\n")
      cat("\n-------------------------------------\n")
      cat("P_Delta coeffs\n")
      cat("--------------\n")
      fit$fits$P_Delta$reg.object$coef
      cat("S_Y_1 coeffs\n")
      cat("------------\n")
      fit$fits$S_Y_1$reg.object$coef
      cat("S_Y_0 coeffs\n")
      cat("------------\n")
      fit$fits$S_Y_0$reg.object$coef
      cat("\n-------------------------------------\n")
    }

    fnc_srv <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (methods::is(newX[["s"]][1],"factor")) {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (methods::is(newX[["s"]][1],"factor")) {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(cens_pred[row,col])
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
#' @param dat Subsample of dataset returned by `ss` for which z==1
#' @param type One of c("parametric", "binning")
#' @param k Number of bins for the binning estimator (if k=0, then the number of
#'     bins will be selected via cross-validation); ignored for the parametric
#'     estimator
#' @param z1 Compute the density conditional on Z=1
#' @return Conditional density estimator function
#' @note
#'   - Assumes support of S is [0,1]
#' @noRd
construct_f_sIx_n <- function(dat, type, k=0, z1=F) {

  if (z1) { dat$weights <- rep(1, length(dat$weights)) }

  if (type=="parametric") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( truncnorm::dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

    # Set up weighted likelihood
    wlik <- function(prm) {
      -1 * sum(sapply(c(1:length(dat$s)), function(i) {
        dat$weights[i] *
          log(pmax(dens_s(s=dat$s[i], x=as.numeric(dat$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer
    opt <- Rsolnp::solnp(
      pars = c(rep(0, length(dat$x)+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp() did not converge")
    }
    prm <- opt$pars

    # Remove large intermediate objects
    rm(dat,dens_s,opt)

    fnc <- function(s, x) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( truncnorm::dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma) )
    }

  }

  if (type=="parametric (edge)") {

    # Density for a single observation
    dens_s <- function(s, x, prm) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( truncnorm::dtruncnorm(s, a=0.01, b=1, mean=mu, sd=sigma) )
    }

    # Estimate p_0
    n_orig <- attr(dat, "n_orig")
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
    opt <- Rsolnp::solnp(
      pars = c(rep(0, length(dat$x)+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp() did not converge")
    }
    prm <- opt$pars

    # Remove large intermediate objects
    rm(dat,dens_s,opt)

    fnc <- function(s, x) {
      mu <- sum( c(1,x) * c(expit(prm[1]),prm[2:(length(x)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      pc_1 <- In(s==0) * ( (1-p_n) / 0.01 )
      pc_2 <- In(s!=0) * p_n *
        truncnorm::dtruncnorm(s, a=0, b=1, mean=mu, sd=sigma)
      return(pc_1+pc_2)
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
      -1 * sum(sapply(c(1:length(dat$s)), function(i) {
        dat$weights[i] *
          log(pmax(prob_s(s=dat$s[i], x=as.numeric(dat$x[i,]), prm),1e-8))
      }))
    }

    # Run optimizer (edge)
    opt_1 <- Rsolnp::solnp(
      pars = rep(0.01, length(dat$x)+1),
      fun = wlik_1
    )
    if (opt_1$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp did not converge")
    }
    prm_1 <- opt_1$pars
    # print("prm_1") # !!!!!
    # print(prm_1) # !!!!!

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
    opt_2 <- Rsolnp::solnp(
      pars = c(rep(0.01, length(dat_1$x)+1), 0.15),
      fun = wlik_2
    )
    if (opt_2$convergence!=0) {
      warning("f_sIx_n: Rsolnp::solnp() did not converge")
    }
    prm_2 <- opt_2$pars

    # Remove large intermediate objects
    rm(dat,dat_1,opt_1,opt_2)

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
      dat_df <- as_df(dat)
      dim_x <- attr(dat, "dim_x")
      wlik <- function(prm) {
        -1 * sum(dat$weights * log(pmax(apply(dat_df, 1, function(r) {
          dens_s(s=r[["s"]], x=as.numeric(r[1:dim_x]), prm)
        }),1e-8)))
      }

      # Run optimizer
      opt <- Rsolnp::solnp(
        pars = rep(0.001,k+length(dat$x)-1),
        fun = wlik
      )
      if (opt$convergence!=0) {
        warning("f_sIx_n: Rsolnp::solnp() did not converge")
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
      print(paste("Number of bins selected (f_sIx_n):", k)) # !!!!! make this a message, configurable with verbose option

    }

    fnc <- create_dens(k, dat, remove_dat=T)

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
construct_f_s_n <- function(dat_orig, f_sIx_n) {

  n_orig <- attr(dat_orig, "n_orig")
  memoise2(function(s) {
    (1/n_orig) * sum(apply(dat_orig$x, 1, function(x) {
      f_sIx_n(s,as.numeric(x))
    }))
  })

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
construct_Phi_n <- function (dat_orig, dat, type="linear (mid)") {

  # !!!!! type currently ignored

  n_orig <- attr(dat_orig, "n_orig")

  fn <- memoise::memoise(function(x) {
    (1/n_orig) * sum(dat$weights*In(dat$s<=x))
  })

  return(fn)

  # # n_orig <- attr(dat, "n_orig")
  # n_orig <- sum(dat$weights) # !!!!!
  # df <- data.frame(s=dat$s, weights=dat$weights)
  # df <- dplyr::arrange(df, s)
  # vals_x <- unique(df$s)
  # vals_y <- c()
  #
  # for (j in 1:length(vals_x)) {
  #   indices <- which(df$s==vals_x[j])
  #   weights_j <- df$weights[indices]
  #   new_y_val <- (1/n_orig) * sum(weights_j)
  #   vals_y <- c(vals_y, new_y_val)
  # }
  # vals_y <- cumsum(vals_y)
  #
  # if (type=="step") {
  #   method <- "constant"
  # } else if (type=="linear (mid)") {
  #   vals_x <- c(vals_x[1], vals_x[1:(length(vals_x)-1)]+diff(vals_x)/2,
  #               vals_x[length(vals_x)])
  #   vals_y <- c(0, vals_y[1:(length(vals_y)-1)], 1)
  #   method <- "linear"
  # }
  #
  # return(stats::approxfun(vals_x, vals_y, method=method, yleft=0, yright=1, f=0,
  #                  ties="ordered"))

}



#' Construct nuisance estimator eta_n
#'
#' @noRd
construct_eta_n <- function(dat, Q_n, t_0) {

  n_orig <- attr(dat, "n_orig")

  return(memoise2(function(u,x) {
    (1/n_orig) * sum(
      dat$weights * sapply(dat$s, function(s) { In(s<=u) * (1-Q_n(t_0,x,s)) })
    )
  }))

}



#' Construct g-computation estimator function of theta_0
#'
#' @param dat_orig Dataset returned by `generate_data`
#' @param Q_n Conditional survival function estimator returned by
#'     `construct_Q_n`
#' @return G-computation estimator of theta_0
#' @noRd
construct_r_tilde_Mn <- function(dat_orig, Q_n, t_0) {

  n_orig <- attr(dat_orig, "n_orig")
  memoise2(function(s) {
    1 - (1/n_orig) * sum(apply(dat_orig$x, 1, function(x) {
      Q_n(t_0, as.numeric(x), s)
    }))
  })

}



#' Construct nuisance estimator Gamma_tilde_n
#'
#' @noRd
construct_Gamma_tilde_n <- function(dat, r_tilde_Mn) {

  n_orig <- attr(dat, "n_orig")
  piece_1 <- dat$weights * sapply(dat$s, r_tilde_Mn)

  return(memoise2(function(u) { (1/n_orig) * sum(In(dat$s<=u)*piece_1) }))

}



#' Construct nuisance estimator of conditional density f(Y,Delta|X,S)
#'
#' @noRd
construct_f_n_srv <- function(Q_n, Qc_n, grid) {

  # Helper function to calculate derivatives
  # !!!!! Move outside?
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
construct_q_n <- function(type="standard", dat, omega_n, g_n, r_tilde_Mn,
                          f_n_srv) {

  if (type=="standard") {

    # !!!!! Can this function be used elsewhere?
    q_n_star_inner <- memoise2(function(x,y,delta,s) {
      omega_n(x,s,y,delta)/g_n(s,x)+r_tilde_Mn(s)
    })

    q_n_star <- memoise2(function(y,delta,x,s,u) {
      if (s>u) {
        return(0)
      } else {
        return(q_n_star_inner(x,y,delta,s))
      }
    })

    n <- length(dat$s)

    # Helper functions
    # !!!!! Can these vectors/functions be used elsewhere?
    f_n_srv_s <- memoise2(function(y,delta,x) {
      sapply(dat$s, function(s) { f_n_srv(y,delta,x,s) })
    })
    q_n_star_s <- memoise2(function(y,delta,x,u) {
      sapply(dat$s, function(s) { q_n_star(y,delta,x,s,u) })
    })
    g_n_s <- memoise2(function(x) {
      sapply(dat$s, function(s) { g_n(s,x) })
    })

    fnc <- function(x, y, delta, u) {

      f_n_srv_s <- f_n_srv_s(y,delta,x)
      q_n_star_s <- q_n_star_s(y,delta,x,u)
      g_n_s <- g_n_s(x)

      denom <- sum(dat$weights*f_n_srv_s*g_n_s)
      if (denom==0) {
        return (0)
      } else {
        num <- sum(dat$weights*q_n_star_s*f_n_srv_s*g_n_s)
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
construct_g_sn <- function(dat, f_n_srv, g_n, p_n) {

  n_orig <- attr(dat, "n_orig")
  dat_df <- as_df(dat)

  return(memoise2(function(x, y, delta) {
    num <- f_n_srv(y, delta, x, 0) * g_n(s=0,x) * (1-p_n)
    den <- (1/n_orig) * sum(
      dat$weights * apply(dat_df, 1, function(r) {
        f_n_srv(y, delta, x, r[["s"]]) * g_n(r[["s"]],x) # !!!!! Can this be replaced with sapply ?????
      })
    )
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
construct_gamma_n <- function(dat_orig, dat, type="Super Learner", omega_n,
                              grid) {

  # Construct pseudo-outcomes
  dat_df <- as_df(dat)
  dim_x <- attr(dat_orig, "dim_x")
  omega_n_d <- apply(dat_df, 1, function(r) {
    omega_n(as.numeric(r[1:dim_x]), r[["s"]], r[["y"]], r[["delta"]])
  })
  dat_df$po <- (dat$weights*as.numeric(omega_n_d))^2

  # Filter out infinite values
  if (sum(!is.finite(dat_df$po))!=0) {
    # dat_df <- dplyr::filter(dat_df, is.finite(po))
    dat_df <- dat_df[which(is.finite(dat_df$po)),] # !!!!! New code (to avoid R CMD CHECK error)
    warning(paste("gamma_n:", sum(!is.finite(dat_df$po)),
                  "non-finite pseudo-outcome values"))
  }

  # Setup
  x_distinct <- dplyr::distinct(dat_orig$x)
  x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
  newX <- expand.grid(x_index=x_distinct$x_index, s=grid$s)
  newX <- dplyr::inner_join(x_distinct, newX, by="x_index")
  newX$x_index <- NULL

  # Run regression
  if (type=="Super Learner") {

    # Fit SuperLearner regression
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.loess",
                    "SL.nnet", "SL.ksvm", "SL.caret", "SL.rpartPrune",
                    "SL.svm")

    model_sl <- SuperLearner::SuperLearner(
      Y = dat_df$po,
      X = cbind(dat_df[,c(1:dim_x), drop=F],s=dat_df$s),
      newX = newX,
      family = "gaussian",
      SL.library = SL.library,
      verbose = F
    )
    pred <- as.numeric(model_sl$SL.predict)
    if (sum(pred<0)!=0) {
      warning(paste("gamma_n:", sum(pred<0), "negative predicted values"))
    }
    # print(coef(model_sl))
    rm(model_sl)

    newX$index <- c(1:nrow(newX))
    reg <- function(x,s) {
      # Dynamically filter to select index
      cond <- paste0("round(s,5)==",round(s,5))
      for (i in c(1:length(x))) {
        cond <- paste0(cond," & round(x",i,",5)==",round(x[i],5))
      }
      index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
      if (length(index)!=1) {
        stop(paste0("Error in gamma_n; ", "x=(", paste(x,collapse=","), "), s=",
                    s))
      }
      return(pred[index])
    }

  }

  # Remove large intermediate objects
  rm(dat_orig,dat,omega_n)

  return(memoise2(reg))

}



#' Construct estimator of g_z0(x,s) = P(Z=1|X=x,S=s)
#'
#' @noRd
construct_g_zn <- function(dat_orig, type="Super Learner", f_sIx_n,
                           f_sIx_z1_n) {

  # Set library
  if (type=="Super Learner") {
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
                    "SL.svm", "SL.glmnet")
  } else if (type=="logistic") {
    SL.library <- c("SL.glm")
  }

  # Create data objects
  newX <- dplyr::distinct(dat_orig$x)

  # Fit SuperLearner regression
  model_sl <- SuperLearner::SuperLearner(
    Y = dat_orig$z,
    X = dat_orig$x,
    newX = newX,
    family = "binomial",
    SL.library = SL.library,
    verbose = F
  )

  pred <- as.numeric(model_sl$SL.predict)
  rm(model_sl)

  # Construct function
  newX$index <- c(1:nrow(newX))
  reg <- function(x) {
    # Dynamically filter to select index
    cond <- "1==1"
    for (i in c(1:length(x))) {
      cond <- paste0(cond," & round(x",i,",5)==",round(x[i],5))
    }
    index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
    if (length(index)!=1) {
      stop(paste0("Error in g_zn; ",
                  "x=(", paste(x,collapse=","), ")"))
    }
    return(pred[index])
  }

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
construct_tau_n <- function(deriv_r_Mn, gamma_n, f_sIx_n, g_zn,
                            dat_orig) {

  n_orig <- attr(dat_orig, "n_orig")
  x <- dat_orig$x
  return(Vectorize(function(u) {
    abs(4*deriv_r_Mn(u) * (1/n_orig) * sum(apply(dat_orig$x, 1, function(x) {
      x <- as.numeric(x)
      return((gamma_n(x,u)*g_zn(x,u))/f_sIx_n(u,x))
    })))^(1/3)
  }))

}



