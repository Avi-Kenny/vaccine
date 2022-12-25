if (F) {

#' Helper function for debugging; prints timestamps
#'
#' @param num Number
#' @param msg Message
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  print(str)
}



#' Probability of sampling
#'
#' @param sampling One of c("iid", "two-phase (6%)", "two-phase (72%)",
#'     "two-phase (70% random)", "two-phase (6% random)")
#' @param delta Component of dataset returned by generate_data()
#' @param y Component of dataset returned by generate_data()
#' @param x Component of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
Pi <- function(sampling, delta, y, x) {

  if (sampling=="iid") {
    probs <- rep(1, length(delta))
  } else if (sampling=="two-phase (70% random)") {
    probs <- rep(0.7, length(delta))
  } else if (sampling=="two-phase (50% random)") {
    probs <- rep(0.5, length(delta))
  } else if (sampling=="two-phase (25% random)") {
    probs <- rep(0.25, length(delta))
  } else if (sampling=="two-phase (6% random)") {
    probs <- rep(0.06, length(delta))
  } else if (sampling=="cycle") {
    probs <- rep(c(0.8,0.6,0.4,0.2), length.out=length(delta))
  } else {
    ev <- In(delta==1 & y<=C$t_0)
    if (sampling=="two-phase (6%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-3.85)
    } else if (sampling=="two-phase (72%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-0.1)
    } else if (sampling=="two-phase (50%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-1)
    } else if (sampling=="two-phase (25%)") {
      probs <- ev + (1-ev)*expit(x$x1+x$x2-2.2)
    } else if (sampling=="x1") {
      probs <- (0.2 + 0.6*x$x1)
    } else if (sampling=="x2") {
      probs <- (0.2 + 0.6*x$x2)
    }

  }

  return(probs)

}



#' Return IP weights
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param scale One of c("none", "stabilized")
#' @param type One of c("true", "estimated")
#' @param return_strata Whether the discrete two-phase sampling strata
#'     membership variable should be returned
#' @return A sum-to-one vector of weights
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
wts <- function(dat_orig, scale="stabilized", type="true", return_strata=F) {

  sampling <- attr(dat_orig,"sampling")
  Pi_0 <- Pi(sampling, dat_orig$delta, dat_orig$y, dat_orig$x)
  strata1 <- In(factor(Pi_0))

  if (type=="true") {

    weights <- dat_orig$z / Pi_0

  } else if (type=="estimated") {

    Pi_vals <- c()
    for (i in c(1:max(strata1))) {
      Pi_vals[i] <- sum(In(strata1==i)*dat_orig$z) / sum(In(strata1==i))
      if (Pi_vals[i]==0) {
        # Hack to avoid NA values in small sample sizes
        warning(paste0("stratum ", i, " had no one sampled."))
        Pi_vals[i] <- 1
      }
    }
    weights <- dat_orig$z / Pi_vals[strata1]

  }

  if (scale=="stabilized") {
    sm <- sum(weights) / length(dat_orig$z)
    weights <- weights / sm
  }

  if (!return_strata) {
    return(weights)
  } else {
    return(list(weights=weights, strata=strata1))
  }

}



#' Construct Phi_n and Phi_n^{-1}
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param which One of c("ecdf", "inverse")
#' @param type One of c("true", "step", "linear (top)", "linear (mid)")
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
construct_Phi_n <- function (dat, which="ecdf", type="step") {

  if (type!="true") {

    n_orig <- sum(dat$weights)
    df <- data.frame(s=dat$s, weights=dat$weights)
    df %<>% arrange(s)
    vals_x <- unique(df$s)
    vals_y <- c()

    for (j in 1:length(vals_x)) {
      indices <- which(df$s==vals_x[j])
      weights_j <- df$weights[indices]
      new_y_val <- (1/n_orig) * sum(weights_j)
      vals_y <- c(vals_y, new_y_val)
    }
    vals_y <- cumsum(vals_y)

    if (type=="step") {
      method <- "constant"
    } else if (type=="linear (top)") {
      method <- "linear"
    } else if (type=="linear (mid)") {
      vals_x <- c(vals_x[1], vals_x[1:(length(vals_x)-1)]+diff(vals_x)/2,
                  vals_x[length(vals_x)])
      vals_y <- c(0, vals_y[1:(length(vals_y)-1)], 1)
      method <- "linear"
    }

    if (which=="ecdf") {
      rval <- approxfun(vals_x, vals_y, method=method, yleft=0,
                        yright=1, f=0, ties="ordered")
    } else if (which=="inverse") {
      rval <- approxfun(vals_y, vals_x, method=method, yleft=min(vals_x),
                        yright=max(vals_x), f=1, ties="ordered")
    }
    return(rval)

  } else if (type=="true") {

    if (L$distr_S=="Unif(0,1)") {
      return(function(s) {s})
    } else if (L$distr_S=="N(0.5,0.01)") {
      if (which=="ecdf") {
        return(function(s) { ptruncnorm(s, a=0, b=1, mean=0.5, sd=0.1) })
      } else {
        return(function(s) { qtruncnorm(s, a=0, b=1, mean=0.5, sd=0.1) })
      }
    } else if (L$distr_S=="N(0.5,0.04)") {
      if (which=="ecdf") {
        return(function(s) { ptruncnorm(s, a=0, b=1, mean=0.5, sd=0.2) })
      } else {
        return(function(s) { qtruncnorm(s, a=0, b=1, mean=0.5, sd=0.2) })
      }
    } else if (L$distr_S=="N(0.3+0.4x2,0.09)") {
      stop("ecdf not programmed for N(0.3+0.4x2,0.09)")
    }

  }

}



#' Construct conditional survival estimator Q_n
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param vals List of values to pre-compute function on; REQUIRED FOR SUPERLEARNER
#' @param type One of c("true", "Cox PH", "Random Forest", "Super Learner")
#' @param return_model Logical; if TRUE, return the model object instead of the
#'     function
#' @param print_coeffs Logical; if TRUE, print the algorithm coefficients for
#'     the conditional survival and censoring functions (applicable only if
#'     type=="Super Learner")
#' @return Conditional density estimator function
construct_Q_n <- function(dat, vals, type, return_model=F, print_coeffs=F) {

  # if (type %in% c("Cox PH2", "Cox PH3")) {
  #
  #   # ..p <<- list() # !!!!!
  #   model_srv <- coxph(
  #     formula = formula(paste0("Surv(y,delta)~",
  #                              paste(names(dat$x),collapse="+"),"+s")),
  #     data = cbind(y=dat$y, delta=dat$delta, dat$x, s=dat$s),
  #     weights = dat$weights
  #   )
  #   coeffs_srv <- as.numeric(model_srv$coefficients)
  #   # ..p$coeffs_srv <<- coeffs_srv # !!!!!
  #   bh_srv <- survival::basehaz(model_srv, centered=FALSE)
  #
  #   model_cens <- coxph(
  #     formula = formula(paste0("Surv(y,delta)~",
  #                              paste(names(dat$x),collapse="+"),"+s")),
  #     data = cbind(y=dat$y, delta=1-dat$delta, dat$x, s=dat$s),
  #     weights = dat$weights
  #   )
  #   coeffs_cens <- as.numeric(model_cens$coefficients)
  #   # ..p$coeffs_cens <<- coeffs_cens # !!!!!
  #   bh_cens <- survival::basehaz(model_cens, centered=FALSE)
  #
  #   fnc_srv <- function(t, x, s) {
  #     if (t==0) {
  #       return(1)
  #     } else {
  #       # Lambda_t <- bh_srv$hazard[which.min(abs(bh_srv$time-t))]
  #       Lambda_t <- bh_srv$hazard[max(which((bh_srv$time<t)==T))]
  #       return(exp(-1*Lambda_t*exp(sum(coeffs_srv*as.numeric(c(x,s))))))
  #     }
  #   }
  #
  #   fnc_cens <- function(t, x, s) {
  #     if (t==0) {
  #       return(1)
  #     } else {
  #       Lambda_t <- bh_cens$hazard[max(which((bh_cens$time<t)==T))]
  #       # Lambda_t <- bh_cens$hazard[which.min(abs(bh_cens$time-t))]
  #       return(exp(-1*Lambda_t*exp(sum(coeffs_cens*as.numeric(c(x,s))))))
  #     }
  #   }
  #
  #   rm(model_srv)
  #   rm(model_cens)
  #
  # }

  if (type %in% c("Cox PH", "Super Learner")) {
    if (type=="Cox PH") {
      methods <- c("survSL.coxph")
    } else if (type=="Super Learner") {
      # Excluding "survSL.rfsrc" for now. survSL.pchSL gives errors.
      methods <- c("survSL.coxph", "survSL.expreg", "survSL.km",
                   "survSL.loglogreg", "survSL.pchreg", "survSL.weibreg")
      # methods <- c("survSL.km", "survSL.pchreg", "survSL.rfsrc") # !!!!!
    }

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)
    srv <- survSuperLearner(
      time = dat$y,
      event = dat$delta,
      X = cbind(dat$x, s=dat$s),
      newX = newX,
      new.times = new.times,
      event.SL.library = methods,
      cens.SL.library = methods,
      obsWeights = dat$weights,
      control = list(
        initWeightAlg = methods[1],
        max.SL.iter = 10
      )
    )

    if (print_coeffs && type=="Super Learner") {
      cat("\n------------------------------\n")
      cat("SuperLearner algorithm weights\n")
      cat("------------------------------\n\n")
      cat("event.coef\n")
      cat("----------\n")
      print(sort(srv$event.coef, decr=T))
      cat("\ncens.coef\n")
      cat("---------\n")
      print(sort(srv$cens.coef, decr=T))
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
      if (class(newX[["s"]][1])=="factor") {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (class(newX[["s"]][1])=="factor") {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(cens_pred[row,col])
    }

  }

  if (type=="true") {

    surv_true <- L$surv_true
    alpha_3 <- L$alpha_3
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    lmbd2 <- L$sc_params$lmbd2
    v2 <- L$sc_params$v2
    if (L$surv_true=="exp") {
      v <- 1   # Overwrite
      v2 <- 1 # Overwrite
    }

    fnc_srv <- function(t, x, s) {
      if (L$surv_true=="Non PH") {
        # !!!!! Blank for now
      } else {
        if (L$surv_true=="Cox PH") {
          if (L$dir=="decr") {
            lin <- C$alpha_1*x[1] + C$alpha_2*x[2] + alpha_3*s
          } else {
            lin <- C$alpha_1*x[1] + C$alpha_2*x[2] + alpha_3*(1-s)
          }
        } else if (L$surv_true=="Complex") {
          if (L$dir=="decr") {
            lin <- alpha_3*expit(20*s-10)
          } else {
            lin <- alpha_3*expit(20*(1-s)-10)
          }
        } else if (L$surv_true=="Step") {
          if (L$dir=="decr") {
            lin <- C$alpha_1*x[1] + C$alpha_2*x[2] + (alpha_3/2)*In(s>0)
          } else {
            lin <- C$alpha_1*x[1] + C$alpha_2*x[2] + (alpha_3/2)*In(s<=0)
          }
        } else if (L$surv_true=="exp") {
          lin <- 0
        }
      }
      return(exp(-1*lmbd*(t^v)*exp(lin)))
    }

    fnc_cens <- function(t, x, s) {
      if (L$surv_true %in% c("Cox PH", "Complex", "Non PH")) {
        lin <- C$alpha_1*x[1] + C$alpha_2*x[2]
      } else if (L$surv_true=="exp") {
        lin <- 0
      }
      return(exp(-1*lmbd2*(t^v2)*exp(lin)))
    }

  }

  if (type=="constant") {
    fnc_srv <- function(t, x, s) { 0.5 }
    fnc_cens <- function(t, x, s) { 0.5 }
  }

  if (type=="Random Forest") {

    # Note: not running this through Super Learner for now

    # Setup
    fml <- "Surv(y,delta)~s"
    for (i in 1:length(dat$x)) { fml <- paste0(fml, "+x",i) }
    fml <- formula(fml)
    df <- cbind("y"=dat$y, "delta"=dat$delta, "s"=dat$s,
                dat$x, "weights"=dat$weights)

    # Survival function
    model_srv <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=20,
                       splitrule="logrank", nsplit=0, case.wt=df$weights,
                       samptype="swor")
    newX <- cbind(s=vals$s, vals$x)[which(vals$t==0),]
    pred_srv <- predict(model_srv, newdata=newX)
    fnc_srv <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-pred_srv$time.interest))
      if (length(row)!=1) { stop("Error in construct_Q_n (B)") }
      if (length(col)!=1) { stop("Error in construct_Q_n (C)") }
      return(pred_srv$survival[row,col])
    }

    # Censoring function
    df$delta <- round(1-df$delta)
    model_cens <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=20,
                        splitrule="logrank", nsplit=0, case.wt=df$weights,
                        samptype="swor")
    pred_cens <- predict(model_cens, newdata=newX)
    fnc_cens <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-pred_cens$time.interest))
      if (length(row)!=1) { stop("Error in construct_Q_n (B)") }
      if (length(col)!=1) { stop("Error in construct_Q_n (C)") }
      return(pred_cens$survival[row,col])
    }

  }

  if (type=="survML") {

    newX <- cbind(vals$x, s=vals$s)[which(vals$t==0),]
    new.times <- unique(vals$t)

    fit <- survMLc(
      time = dat$y,
      event = dat$delta,
      X = cbind(dat$x, s=dat$s),
      newX = newX,
      newtimes = new.times,
      bin_size = 0.025,
      time_basis = "continuous",
      time_grid_approx = sort(unique(dat$y)),
      surv_form = "exp",
      SL.library = rep(c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),2),
      V = 5,
      obsWeights = dat$weights
    )

    srv_pred <- fit$S_T_preds
    cens_pred <- fit$S_C_preds

    if (print_coeffs && type=="Super Learner") {
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
      if (class(newX[["s"]][1])=="factor") {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(srv_pred[row,col])
    }

    fnc_cens <- function(t, x, s) {
      r <- list()
      for (i in 1:length(x)) {
        r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
      }
      if (class(newX[["s"]][1])=="factor") {
        r[[length(x)+1]] <- which(s==newX[["s"]])
      } else {
        r[[length(x)+1]] <- which(abs(s-newX[["s"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",x=(",
                    paste(x,collapse=","),"),s=",s,""))
      }
      return(cens_pred[row,col])
    }

  }

  if (type=="Cox PH3") {
    sfnc_srv <- fnc_srv
    sfnc_cens <- fnc_cens
  } else {
    sfnc_srv <- memoise2(fnc_srv)
    sfnc_cens <- memoise2(fnc_cens)
  }
  rm("vals", envir=environment(get("fnc_srv",envir=environment(sfnc_srv))))

  return(list(srv=sfnc_srv, cens=sfnc_cens))

}



#' Construct g-computation estimator function of theta_0
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param vals List of values to pre-compute function on
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @return G-computation estimator of theta_0
construct_r_tilde_Mn <- function(dat_orig, vals=NA, Q_n) {

  fnc <- function(s) {
    1 - mean(Q_n(
      rep(C$t_0, nrow(dat_orig$x)),
      dat_orig$x,
      rep(s, nrow(dat_orig$x))
    ))
  }

  return(memoise2(fnc))

}



#' Construct derivative estimator r_Mn'_n
#'
#'
#' @param r_Mn An estimator of r_M0
#' @param type One of c("gcomp", "linear", "spline")
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @param dat_orig Only used for type="true"
construct_deriv_r_Mn <- function(r_Mn, type, dir="decr", dat_orig=NA) {

  # Estimate entire function on grid
  if (type!="true") {
    grid <- round(seq(0,1,0.01),2)
    r_Mns <- r_Mn(grid)
  }

  if (type=="true") {

    fnc <- function(s) {
      s_idx <- which.min(abs(s-C$points))
      ch_y <- diff(attr(dat_orig, "r_M0"))[s_idx]
      if (is.na(ch_y)) { ch_y <- diff(attr(dat_orig,"r_M0"))[round(s_idx-1)] }
      ch_x <- C$points[2] - C$points[1] # Assumes equal spacing of points
      return(ch_y/ch_x)
    }

  }

  if (type=="linear") {

    grid_width <- grid[2] - grid[1]
    points_x <- c(grid[1])
    points_y <- c(r_Mns[1])
    for (i in 2:length(grid)) {
      if (r_Mns[i]-r_Mns[i-1]!=0) {
        points_x <- c(points_x, grid[i]-(grid_width/2))
        points_y <- c(points_y, mean(c(r_Mns[i],r_Mns[i-1])))
      }
    }
    points_x <- c(points_x, grid[length(grid)])
    points_y <- c(points_y, r_Mns[length(grid)])
    points_sl <- c()
    for (i in 2:length(points_x)) {
      slope <- (points_y[i]-points_y[i-1]) /
        (points_x[i]-points_x[i-1])
      points_sl <- c(points_sl, slope)
    }

    fnc <- function(s) {
      if (s==0) {
        index <- 1
      } else {
        index <- which(round(s,6)<=round(points_x,6))[1]-1
      }
      if (dir=="incr") {
        return(max(points_sl[index],0))
      } else {
        return(min(points_sl[index],0))
      }
    }

  }

  if (type=="line") {

    r_Mn_left <- r_Mn(0) # 0.1
    r_Mn_right <- r_Mn(1) # 0.9
    fnc <- function(s) { (r_Mn_right-r_Mn_left)/1 } # 0.8

  }

  if (type=="spline") {

    # Identify jump points of step function
    jump_points <- c(0)
    for (i in 2:length(grid)) {
      if (r_Mns[i]!=r_Mns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    jump_points <- c(jump_points,grid[length(grid)])

    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)

    # Fit cubic smoothing spline
    r_Mn_smoothed <- smooth.spline(x=midpoints, y=r_Mn(midpoints))

    # Construct derivative function
    fnc <- function(s) {

      width <- 0.2
      x1 <- s - width/2
      x2 <- s + width/2
      if (x1<0) { x2<-width; x1<-0 }
      if (x2>1) { x1<-1-width; x2<-1; }

      y1 <- predict(r_Mn_smoothed, x=x1)$y
      y2 <- predict(r_Mn_smoothed, x=x2)$y

      if (dir=="incr") {
        return(max((y2-y1)/width,0))
      } else {
        return(min((y2-y1)/width,0))
      }

    }

  }

  if (type=="m-spline") {

    # Identify jump points of step function
    jump_points <- c(0)
    for (i in 2:length(grid)) {
      if (r_Mns[i]!=r_Mns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    jump_points <- c(jump_points,grid[length(grid)])

    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)

    if (length(midpoints)>=2) {
      # Fit monotone cubic smoothing spline
      r_Mn_smoothed <- splinefun(
        x = midpoints,
        y = r_Mn(midpoints),
        method = "monoH.FC"
      )
    } else {
      # Fit a straight line instead if there are <2 midpoints
      r_Mn_smoothed <- function(u) {
        slope <- r_Mn(1) - r_Mn(0)
        intercept <- r_Mn(0)
        return(intercept + slope*u)
      }
    }

    # Construct derivative function
    fnc <- function(s) {

      width <- 0.2
      x1 <- s - width/2
      x2 <- s + width/2

      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }

      y1 <- r_Mn_smoothed(x1)
      y2 <- r_Mn_smoothed(x2)

      if (dir=="incr") {
        return(max((y2-y1)/width,0))
      } else {
        return(min((y2-y1)/width,0))
      }

    }

  }

  if (type=="gcomp") {

    fnc <- function(s) {

      # Set derivative appx x-coordinates
      width <- 0.2
      p1 <- s - width/2
      p2 <- s + width/2
      if (p1<0) { p2 <- width; p1 <- 0; }
      if (p2>1) { p1<-1-width; p2<-1; }
      c(p1,p2)

      if (dir=="incr") {
        return(max((r_Mn(p2)-r_Mn(p1))/width,0))
      } else {
        return(min((r_Mn(p2)-r_Mn(p1))/width,0))
      }

    }

  }

  return(Vectorize(fnc))

}



#' Construct tau_n Chernoff scale factor function
#'
#' @param deriv_r_Mn A derivative estimator returned by construct_deriv_r_Mn()
#' @param gamma_n Nuisance function estimator returned by construct_gamma_n()
#' @param f_s_n Density estimator returned by construct_f_s_n()
#' @return Chernoff scale factor estimator function
construct_tau_n <- function(deriv_r_Mn=NA, gamma_n=NA, f_sIx_n=NA, f_s_n=NA,
                            g_zn=NA, dat_orig=NA) {

  n_orig <- length(dat_orig$s)
  x <- dat_orig$x
  # return(Vectorize(function(u) {
  #   abs(
  #     ((4*deriv_r_Mn(u))/(n_orig*f_s_n(u))) *
  #       sum((gamma_n(x,u)*g_zn(x,u))/g_n(u,x))
  #   )^(1/3)
  # }))
  return(Vectorize(function(u) {
    # abs(
    #   ((4*deriv_r_Mn(u))/f_s_n(u)) * (1/n_orig)*sum(
    #     (gamma_n(x,rep(u,n_orig))*g_zn(x,rep(u,n_orig))) / g_n(rep(u,n_orig),x)
    #   ))^(1/3)
    abs(
      4*deriv_r_Mn(u) * (1/n_orig)*sum(
        (gamma_n(x,rep(u,n_orig))*g_zn(x,rep(u,n_orig))) /
          f_sIx_n(rep(u,n_orig),x)
      ))^(1/3)
  }))

}





#' Construct gamma_n nuisance estimator function
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param vals List of values to pre-compute function on
#' @param type Type of regression; one of c("cubic", "kernel", "kernel2")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_sIx_n A conditional density estimator returned by
#'     construct_f_sIx_n()
#' @param f_sIx_n A conditional density estimator returned by
#'     construct_f_sIx_n() among the observations for which z==1
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, dat, type="Super Learner", vals=NA,
                              omega_n=NA, f_sIx_n=NA, f_s_n=NA,
                              f_s_z1_n=NA) {

  # Construct pseudo-outcomes
  po <- (dat$weights*omega_n(dat$x,dat$s,dat$y,dat$delta))^2

  # Create dataframe for regression
  df <- data.frame(s=dat$s, po=po)
  if (sum(!is.finite(df$po))!=0) {
    df %<>% filter(is.finite(po))
    warning(paste("construct_gamma_n:", sum(!is.finite(df$po)),
                  "non-finite po values"))
  }

  # Setup
  x_grid <- round(seq(0,1,C$appx$s), -log10(C$appx$s))
  X <- cbind(dat$x, s=dat$s)
  X_reduced <- distinct(dat_orig$x)
  X_reduced <- cbind("x_index"=c(1:nrow(X_reduced)), X_reduced)
  newX <- expand.grid(
    x_index = X_reduced$x_index,
    s = x_grid
  )
  newX <- inner_join(X_reduced, newX, by="x_index")
  newX$x_index <- NULL

  # Run regression
  if (type=="Super Learner") {

    # Fit SuperLearner regression
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.loess",
                    "SL.nnet", "SL.ksvm", "SL.caret", "SL.rpartPrune",
                    "SL.svm")

    model_sl <- SuperLearner(Y=df$po, X=X, newX=newX, family="gaussian",
                             SL.library=SL.library, verbose=F)
    pred <- as.numeric(model_sl$SL.predict)
    if (sum(pred<0)!=0) {
      warning(paste("construct_gamma_n:", sum(pred<0),
                    "negative predicted values"))
    }
    if (F) {
      print(coef(model_sl))
    } # DEBUG: Print algorithm coeffs

    rm(model_sl)

    newX$index <- c(1:nrow(newX))
    reg <- function(x,s) {
      # Dynamically filter to select index
      cond <- paste0("s==",s,"")
      for (i in c(1:length(x))) { cond <- paste0(cond," & x",i,"==",x[i]) }
      index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
      if (length(index)!=1) {
        stop(paste0("Error in construct_gamma_n; ", "x=(",
                    paste(x,collapse=","), "), s=",s))
      }

      # Return prediction
      return(pred[index])
    }

  }

  # Remove large intermediate objects
  rm(dat_orig,dat,omega_n,f_sIx_n)

  fnc <- reg
  return(memoise2(fnc))

}



#' Construct nuisance estimator of conditional density f(Y,Delta|X,S)
#'
#' @param Q_n xxx
#' @param Qc_n xxx
#' @param width xxx
#' @param vals xxx
#' @return nuisance estimator function
construct_f_n_srv <- function(Q_n=NA, Qc_n=NA, width=40, vals=NA) {

  # Helper function to calculate derivatives
  surv_deriv <- function(Q_n) {
    fnc <- function(t, x, s) {
      t1 <- t - width/2
      t2 <- t + width/2
      if (t1<0) { t2<-width; t1<-0; }
      if (t2>C$t_0) { t1<-C$t_0-width; t2<-C$t_0; }
      t1 <- round(t1,-log10(C$appx$t_0))
      t2 <- round(t2,-log10(C$appx$t_0))
      ch_y <- Q_n(t2,x,s) - Q_n(t1,x,s)
      return(ch_y/width)
    }
    return(memoise2(fnc))
  }

  # Calculate derivative estimators
  Q_n_deriv <- surv_deriv(Q_n)
  Qc_n_deriv <- surv_deriv(Qc_n)

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
#' @param type One of c("new", "zero")
#' @return q_n nuisance estimator function
construct_q_n <- function(type="new", dat, dat_orig,
                          omega_n=NA, g_n=NA, p_n=NA, r_tilde_Mn=NA,
                          Gamma_tilde_n=NA, f_sIx_n=NA, Q_n=NA, Qc_n=NA,
                          f_n_srv=NA, vals=NA) {

  if (type=="new") {

    q_n_star_inner <- function(y, delta, x, s) {
      In(s!=0)*omega_n(x,s,y,delta)/g_n(s,x) + r_tilde_Mn(s)
    }
    q_n_star_inner <- memoise2(q_n_star_inner)

    q_n_star <- function(y, delta, x, s, u) {
      In(s<=u)*q_n_star_inner(y, delta, x, s) -
        In(s!=0)*Gamma_tilde_n(u)
    }
    q_n_star <- memoise2(q_n_star)

    n <- length(dat$s)

    fnc <- function(x, y, delta, u) {

      y_ <- rep(y,n)
      delta_ <- rep(delta,n)
      x_ <- as.data.frame(matrix(rep(x,n), ncol=length(x), byrow=T))
      u_ <- rep(u,n)
      f_n_srv_ <- f_n_srv(y_, delta_, x_, dat$s)
      q_n_star_ <- q_n_star(y, delta_, x_, dat$s, u_)
      g_n_ <- g_n(dat$s,x_)

      denom <- sum(dat$weights*f_n_srv_*g_n_)
      if (denom==0) {
        return (0)
      } else {
        num <- sum(dat$weights*q_n_star_*f_n_srv_*g_n_)
        return(num/denom)
      }

    }

  }

  if (type=="zero") {

    fnc <- function(x, y, delta, u) { 0 }

  }

  return(memoise2(fnc))

}



#' Construct q_tilde_n nuisance estimator function
#'
#' @param x x
#' @return q_tilde_n nuisance estimator function
construct_q_tilde_n <- function(type="new", f_n_srv=NA, f_sIx_n=NA,
                                omega_n=NA, vals=NA) {

  if (type=="new") {

    seq_01 <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))

    fnc <- function(x, y, delta, u) {

      denom <- C$appx$s * sum(sapply(seq_01, function(s) {
        f_n_srv(y, delta, x, s) * f_sIx_n(s,x)
      }))

      if (denom==0) {
        return (0)
      } else {
        if (u==0) {
          return(0)
        } else {
          seq_0u <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
          num <- C$appx$s * sum(sapply(seq_0u, function(s) {
            omega_n(x,s,y,delta) * f_n_srv(y, delta, x, s)
          }))
          return(num/denom)
        }
      }

    }

  }

  if (type=="zero") {

    fnc <- function(x, y, delta, u) { 0 }

  }

  return(memoise2(fnc))

}






#' Construct estimator of marginal density f_S
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param vals List of values to pre-compute function on
#' @param f_sIx_n A conditional density estimator returned by
#'     construct_f_sIx_n()
#' @return Marginal density estimator function
construct_f_s_n <- function(dat_orig, vals=NA, f_sIx_n) {

  fnc <- function(s) {
    mean(f_sIx_n(rep(s,nrow(dat_orig$x)),dat_orig$x))
  }

  return(memoise2(fnc))

}



#' Construct density ratio estimator g_n
#'
#' @param f_sIx_n Conditional density estimator returned by construct_f_sIx_n
#' @param f_s_n Marginal density estimator returned by construct_f_s_n
#' @return Density ratio estimator function
construct_g_n <- function(f_sIx_n, f_s_n, vals=NA) {

  fnc <- function(s,x) { f_sIx_n(s,x) / f_s_n(s) }

  return(memoise2(fnc))

}



#' Construct nuisance estimator eta*_n
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param vals List of values to pre-compute function on
#' @return Estimator function of nuisance eta*_0
construct_etastar_n <- function(Q_n, vals=NA, tmp) {

  fnc <- function(u,x) {
    u <- round(u,-log10(C$appx$s))
    if (u==0) {
      return(0)
    } else {
      s_seq <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
      integral <- C$appx$s * sum(sapply(s_seq, function(s) {
        Q_n(C$t_0, x, s)
      }))
      return(u-integral)
    }
  }

  return(memoise2(fnc))

}



#' lambda estimator
#'
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @return Value of lambda
lambda <- function(dat, k, G) {

  n_orig <- sum(dat$weights)
  lambda <- (1/n_orig) * sum( dat$weights * (G(dat$s))^k )
  return(lambda)

}



#' Construct Theta_os_n primitive one-step estimator
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param vals List of values to pre-compute function on
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_sIx_n Conditional density estimator returned by construct_f_sIx_n
#' @param etastar_n A nuisance estimator returned by construct_etastar_n()
#' @return Gamma_os_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Theta_os_n <- function(dat, dat_orig, omega_n=NA, f_sIx_n=NA,
                                 q_tilde_n=NA, etastar_n=NA, vals=NA) {

  n_orig <- length(dat_orig$z)
  piece_1 <- (dat$weights*omega_n(dat$x,dat$s,dat$y,dat$delta)) /
    f_sIx_n(dat$s,dat$x)
  piece_2 <- (1-dat_orig$weights)

  # Remove large intermediate objects
  rm(omega_n,f_sIx_n)

  fnc <- function(u) {
    (1/(n_orig)) * sum(piece_1*In(dat$s<=u)) +
      (1/n_orig) * sum(
        piece_2 * q_tilde_n(dat_orig$x, dat_orig$y, dat_orig$delta, u) +
          etastar_n(rep(u,length(dat_orig$z)),dat_orig$x)
      )
  }

  # !!!!! DEBUG
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

  }

  return(memoise2(fnc))

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_infl_fn_Theta <- function(omega_n, f_sIx_n, q_tilde_n, etastar_n,
                                    Theta_os_n) {

  fnc <- function(u,x,y,delta,s,wt) {
    if (wt==0) {
      piece_1 <- 0
      piece_2 <- 0
    } else {
      piece_1 <- In(s<=u)
      piece_2 <- omega_n(x,s,y,delta)/f_sIx_n(s,x)
    }
    wt*piece_1*piece_2 +
      (1-wt) * q_tilde_n(x,y,delta,u) +
      etastar_n(u,x) -
      Theta_os_n(round(u,-log10(C$appx$s)))
  }

  return(memoise2(fnc))

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_infl_fn_beta_n <- function(infl_fn_Theta) {

  u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
  m <- length(u_mc)
  lambda_1 <- mean(u_mc) # ~1/2
  lambda_2 <- mean((u_mc)^2) # ~1/3
  lambda_3 <- mean((u_mc)^3) # ~1/4

  fnc <- function(s, y, delta, weight, x) {

    s_m <- rep(s,m)
    y_m <- rep(y,m)
    delta_m <- rep(delta,m)
    weight_m <- rep(weight,m)
    x_m <- as.data.frame(matrix(rep(x,m), ncol=length(x), byrow=T))

    return((1/m) * sum(
      infl_fn_Theta(u=u_mc, x_m, y_m, delta_m, s_m, weight_m) * (
        (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
          (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
      ))
    )

  }

  return(fnc)

}



#' !!!!! document
#'
#' @param x x
#' @return x
construct_Gamma_cf_k <- function(dat_train, dat_test, vals=NA, omega_n, g_n,
                                 r_tilde_Mn, eta_n) {

  # !!!!! Needs to be updated

  n_test <- length(dat_test$s)
  dat_test$weights <- wts(dat_test) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d1 <- ss(dat_test, which(dat_test$z==1))
  weights_1 <- d1$weights

  n_train <- length(dat_train$s)
  dat_train$weights <- wts(dat_train) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d2 <- ss(dat_train, which(dat_train$z==1))
  weights_2 <- d2$weights

  fnc <- function(s) {

    piece_1 <- (1/n_test) * sum(weights_1 * (
      In(d1$s<=s) *
        ((omega_n(d1$x,d1$s,d1$y,d1$delta) / g_n(d1$s,d1$x)) +
           r_tilde_Mn(d1$s)
        ) +
        eta_n(s,d1$x)
    ))

    piece_2 <- (1/n_train) * sum(weights_2 * In(d2$s<=s)*r_tilde_Mn(d2$s))

    return(piece_1-piece_2)

  }

  return(memoise2(fnc))

}



#' Construct cross-fitted Gamma_0 estimator
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param params The same params object passed to est_curve
#' @param vlist A list of dataframes returned by create_val_list(); Q_n REQUIRED
#'     but others can be NA
#' @return A cross-fitted Gamma_0 estimator
construct_Gamma_cf <- function(dat_orig, params, vlist) {

  # !!!!! Update this

  # Prep for cross-fitting
  Gamma_cf_k <- list()
  n_orig <- length(dat_orig$z)
  rows <- c(1:n_orig)
  folds <- sample(cut(rows, breaks=params$cf_folds, labels=FALSE))

  # Loop through folds
  for (k in 1:params$cf_folds) {

    # Split data
    dat_orig_train <- ss(dat_orig, -which(folds==k))
    dat_orig_test <- ss(dat_orig, which(folds==k))
    dat_train <- ss(dat_orig_train, which(dat_orig_train$z==1))
    dat_test <- ss(dat_orig_test, which(dat_orig_test$z==1))

    # Construct component functions
    Phi_n <- construct_Phi_n(dat_train, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat_train, which="inverse",
                                 type=params$ecdf_type)
    srvSL <- construct_Q_n(dat_train, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    r_tilde_Mn <- construct_r_tilde_Mn(dat_orig_train, vlist$S_grid, Q_n)
    # eta_n <- construct_eta_n(dat_train, vlist$SX_grid, Q_n) # ARCHIVED
    f_sIx_n <- construct_f_sIx_n(dat_train, vlist$SX_grid, type=params$g_n_type,
                                 k=15, s_scale=s_scale, s_shift=s_shift)
    f_s_n <- construct_f_s_n(dat_orig_train, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n,
                                 type=params$omega_n_type)

    # Construct K functions
    Gamma_cf_k[[k]] <- construct_Gamma_cf_k(
      dat_train, dat_test, vlist$S_grid, omega_n, g_n, r_tilde_Mn, eta_n
    )

    # Remove objects
    rm(Phi_n,Phi_n_inv,Q_n,Qc_n,r_tilde_Mn,eta_n,f_sIx_n,f_s_n,g_n,omega_n)

  }

  # Construct cross-fitted Gamma_os_n
  return(Vectorize(function(s) {
    mean(sapply(c(1:params$cf_folds), function(k) {
      Gamma_cf_k[[k]](s)
    }))
  }))

}



#' #' Construct propensity score estimator of g_s0 = P(S=0|X=x)
#' #'
#' #' @param dat Subsample of dataset returned by ss() for which z==1
#' #' @param vals List of values to pre-compute function on; REQUIRED FOR SUPERLEARNER
#' #' @param type One of c("true", "logistic", "Super Learner", "generalized"). If
#' #'     type=="true", the only valid value is zero. If type=="generalized", the
#' #'     arguments `f_sIx_n` and `cutoffs` must also be supplied.
#' #' @return Propensity score estimator of g_s0
#' #' @notes For all types except for "generalized", this function constructs the
#' #'     probability P(S=0|X=x). The type "generalized" constructs the probability
#' #'     P(S=s|X=x) for a generic value a that has positive mass.
#' construct_g_sn <- function(dat, vals=NA, type, f_sIx_n=NA, cutoffs=NA) {
#'
#'   # Construct indicator I{S=0}
#'   ind_S0 <- In(dat$s==0)
#'
#'   if (type=="true") {
#'
#'     # Note: These are only valid for val==0
#'     if (L$edge=="expit") {
#'       fnc <- function(x, val) {
#'         if (val==0) { expit(x[1]+x[2]-3.3) }
#'       }
#'     } else if (L$edge=="expit2") {
#'       fnc <- function(x, val) {
#'         if (val==0) { expit(x[1]+x[2]-1) }
#'       }
#'     } else if (L$edge=="Complex") {
#'       fnc <- function(x, val) {
#'         if (val==0) { 0.84*x[2]*pmax(0,1-4*abs(x[1]-0.5)) }
#'       }
#'     } else if (L$edge=="none") {
#'       fnc <- function(x, val) {
#'         if (val==0) { 0 }
#'       }
#'     }
#'
#'   } else if (type=="logistic") {
#'
#'     fml <- "ind_S0~1"
#'     for (i in 1:length(dat$x)) {
#'       fml <- paste0(fml, "+x",i)
#'     }
#'     fml <- formula(fml)
#'     df <- cbind(ind_S0, dat$x)
#'     suppressWarnings({
#'       model <- glm(
#'         fml,
#'         data = df,
#'         family = "binomial",
#'         weights = dat$weights
#'       )
#'     })
#'     coeffs <- model$coefficients
#'
#'     # !!!!! val is unused
#'     fnc <- function(x, val) { as.numeric(expit(coeffs %*% c(1,x))) }
#'
#'   } else if (type=="SL") {
#'
#'     # sl <- SuperLearner(
#'     #   Y = ind_S0,
#'     #   X = dat$x,
#'     #   newX = vals,
#'     #   family = binomial(),
#'     #   SL.library = "SL.earth", # SL.glm SL.gbm SL.ranger SL.earth
#'     #   # SL.library = c("SL.earth", "SL.gam", "SL.ranger"), # SL.glm SL.gbm SL.ranger SL.earth
#'     #   obsWeights = dat$weights,
#'     #   control = list(saveFitLibrary=FALSE)
#'     # )
#'     # assign("sl", sl, envir=.GlobalEnv) # ?????
#'     #
#'     # fnc <- function(x, val) {
#'     #
#'     #   r <- list()
#'     #   for (i in 1:length(x)) {
#'     #     r[[i]] <- which(abs(x[i]-newX[[paste0("x",i)]])<1e-8)
#'     #   }
#'     #   index <- Reduce(intersect, r)
#'     #   return(sl$SL.predict[index])
#'     # }
#'
#'   } else if (type=="generalized") {
#'
#'     fnc <- function(x, val) {
#'
#'       # val should be an index of the bin
#'       bin_start <- cutoffs[as.numeric(val)]
#'       bin_stop <- cutoffs[as.numeric(val)+1]
#'       grid <- seq(bin_start, bin_stop, length.out=11)[1:10]
#'       avg_height <- mean(sapply(grid, function(s) { f_sIx_n(s, x=x) }))
#'
#'       return((bin_stop-bin_start)*avg_height)
#'
#'     }
#'
#'   }
#'
#'   return(memoise2(fnc))
#'
#' }


#' Construct estimator of nuisance g_sn
#'
#' @param x TO DO
#' @return TO DO
construct_g_sn <- function(dat, f_n_srv, g_n, p_n, vals=NA) {

  n_orig <- sum(dat$weights)

  fnc <- function(x, y, delta) {

    num <- f_n_srv(y, delta, x, 0) * g_n(0,x) * (1-p_n)
    den <- (1/n_orig) * sum(unlist(lapply(c(1:length(dat$s)), function(i) {
      dat$weights[i] * f_n_srv(y, delta, x, dat$s[i]) * g_n(dat$s[i],x)
    })))

    if (den==0) {
      warning("Denominator of zero in g_sn()")
    }

    return(num/den)

  }

  return(memoise2(fnc))

}



#' Construct estimator of g_z0(x,s) = P(Z=1|X=x,S=s)
#'
#' @param x !!!!!
construct_g_zn <- function(dat_orig, vals=NA, type="Super Learner",
                           f_sIx_n=NA, f_sIx_z1_n=NA) {

  # Set library
  if (type=="Super Learner") {
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
                    "SL.svm", "SL.glmnet")
  } else if (type=="logistic") {
    SL.library <- c("SL.glm")
  }

  # Create data objects
  X <- dat_orig$x
  newX <- distinct(dat_orig$x)

  # Fit SuperLearner regression
  model_sl <- SuperLearner(Y=dat_orig$z, X=X, newX=newX, family="binomial",
                           SL.library=SL.library, verbose=F)

  pred <- as.numeric(model_sl$SL.predict)
  rm(model_sl)

  # Construct function
  newX$index <- c(1:nrow(newX))
  reg <- function(x) {

    # Dynamically filter to select index
    cond <- "1==1"
    for (i in c(1:length(x))) {
      cond <- paste0(cond," & round(x",i,",8)==",round(x[i],8))
    }
    index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
    if (length(index)!=1) {
      stop(paste0("Error in construct_g_zn; ",
                  "x=(", paste(x,collapse=","), ")"))
    }

    # Return prediction
    return(pred[index])

  }

  if (F) {
    d0 <- ss(dat_orig, which(dat_orig$x$x2==0))
    d1 <- ss(dat_orig, which(dat_orig$x$x2==0))
    grid <- round(seq(0,1,0.1),1)
    preds0 <- sapply(grid, function(x1) { reg(c(x1,0)) })
    preds1 <- sapply(grid, function(x1) { reg(c(x1,1)) })
    ggplot(
      data.frame(x=d0$x$x1, y=d0$z),
      aes(x=x, y=y)
    ) +
      geom_jitter(alpha=0.3, width=0.04, height=0.04) +
      geom_line(data=data.frame(x=grid, y=preds0), color="forestgreen")
    ggplot(
      data.frame(x=d1$x$x1, y=d1$z),
      aes(x=x, y=y)
    ) +
      geom_jitter(alpha=0.3, width=0.04, height=0.04) +
      geom_line(data=data.frame(x=grid, y=preds1), color="forestgreen")
  } # DEBUG

  fnc <- function(x,s) {
    (f_sIx_z1_n(s,x) / f_sIx_n(s,x)) * reg(x)
  }

  return(memoise2(fnc))

}



#' #' Compute one-step estimator of counterfactual survival at S=0
#' #'
#' #' @param dat Subsample of dataset returned by ss() for which z==1
#' #' @param g_sn Propensity score estimator returned by construct_g_sn()
#' #' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' #' @param omega_n A nuisance influence function returned by construct_omega_n()
#' #' @param val Value of S
#' #' @return Value of one-step estiamtor
#' r_Mn_edge <- function(dat, g_sn, Q_n, omega_n, val=0) {
#'
#'   n_orig <- sum(dat$weights)
#'   n_dat <- nrow(dat$x)
#'
#'   return(
#'     1 - (1/n_orig) * sum(dat$weights * (
#'       Q_n(rep(C$t_0,n_dat),dat$x,s=rep(val,n_dat)) - (
#'         (In(dat$s==val)/g_sn(dat$x, rep(val,n_dat))) *
#'           omega_n(dat$x,s=rep(val,n_dat),dat$y,dat$delta)
#'       )
#'     ))
#'   )
#'
#' }



#' Compute one-step estimator of counterfactual survival at S=0
#'
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param g_sn Propensity score estimator returned by construct_g_sn()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param val Value of S
#' @return Value of one-step estimator
r_Mn_edge <- function(dat_orig, dat, g_sn, g_n, p_n, Q_n, omega_n, val=0) {

  n_orig <- sum(dat$weights)

  v <- 0
  for (i in c(1:n_orig)) {

    z_i <- dat_orig$z[i]
    weight_i <- dat_orig$weights[i]
    if (z_i==0) {
      s_i <- 0
      pi_i <- 1
    } else {
      s_i <- dat_orig$s[i]
      pi_i <- 1/weight_i
    }
    x_i <- as.numeric(dat_orig$x[i,])
    y_i <- dat_orig$y[i]
    delta_i <- dat_orig$delta[i]

    v <- v + (
      ( z_i*In(s_i==val) + g_sn(x_i, y_i, delta_i)*(pi_i-z_i) ) *
        omega_n(x_i,s=val,y_i,delta_i)
    ) / (pi_i * (1-p_n) * g_n(val,x_i)) - Q_n(C$t_0,x_i,s=val)

  }
  v <- v / n_orig

  return(1+v)

}



#' #' Construct influence function corresponding to r_Mn_edge
#' #'
#' #' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' #' @param g_sn Propensity score estimator returned by construct_g_sn()
#' #' @param omega_n A nuisance influence function returned by construct_omega_n()
#' #' @param r_Mn_edge_est Estimate returned by one-step estimator r_Mn_edge()
#' #' @param val Value of S
#' #' @return Value of one-step estimator
#' construct_infl_fn_r_Mn_edge <- function(Q_n, g_sn, omega_n, r_Mn_edge_est,
#'                                         val=0, vals=NA) {
#'
#'   fnc <- function(weight, s, x, y, delta) {
#'     if (weight==0) {
#'       return(0)
#'     } else {
#'       return(weight * (
#'         1 - Q_n(C$t_0,x,s=val) + (
#'           (In(s==val)/g_sn(x,val)) * omega_n(x,s=val,y,delta)
#'         ) -
#'           r_Mn_edge_est
#'       ))
#'     }
#'   }
#'
#'   return(memoise2(fnc))
#'
#' }



#' Construct influence function corresponding to r_Mn_edge
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param g_sn Propensity score estimator returned by construct_g_sn()
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param r_Mn_edge_est Estimate returned by one-step estimator r_Mn_edge()
#' @param val Value of S
#' @return Value of one-step estimator
construct_infl_fn_r_Mn_edge <- function(Q_n, g_sn, omega_n, g_n,
                                        r_Mn_edge_est, p_n, val=0, vals=NA) {

  fnc <- function(z, weight, s, x, y, delta) {

    if (z==0) {
      s <- 0
      pi_ <- 1
    } else {
      pi_ <- 1/weight
    }

    return(
      1 - Q_n(C$t_0, x, s=val) + (
        ( z*In(s==val) + g_sn(x, y, delta)*(pi_-z) ) *
          omega_n(x,s=val,y,delta)
      ) / (pi_ * (1-p_n) * g_n(val,x)) -
        r_Mn_edge_est
    )

  }

  return(memoise2(fnc))

}



#' Construct influence function corresponding to test statistic beta_en
#'
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @return Influence function
construct_infl_fn_beta_en <- function(infl_fn_Theta, infl_fn_r_Mn_edge) {

  return(function(z, x, y, delta, s, weight) {
    infl_fn_Theta(u=1, x, y, delta, s, weight) -
      infl_fn_r_Mn_edge(z, weight, s, x, y, delta)
  })

}



#' Construct nuisance estimator alpha*_n
#'
#' @param dat !!!!!
#' @param gcomp !!!!!
#' @param p_n !!!!!
#' @param vals !!!!!
#' @return Nuisance estimator function
construct_Gamma_tilde_n <- function(dat, r_tilde_Mn, p_n, vals=NA) {

  n_orig <- sum(dat$weights)
  piece_1 <- dat$weights * In(dat$s!=0) * r_tilde_Mn(dat$s)

  fnc <- function(u) {
    (1/(n_orig*p_n)) * sum(In(dat$s<=u)*piece_1)
  }

  return(memoise2(fnc))

}



#' Construct nuisance estimator eta_n
#'
#' @param dat !!!!!
#' @param Q_n !!!!!
#' @param p_n !!!!!
#' @param vals !!!!!
#' @return Nuisance estimator function
construct_eta_n <- function(dat, Q_n, p_n, vals=NA) {

  n_orig <- sum(dat$weights)
  piece_1 <- dat$weights * In(dat$s!=0)

  fnc <- function(u,x) {
    x_long <- as.data.frame(
      matrix(rep(x,length(dat$s)), ncol=length(x), byrow=T)
    )
    (1/(n_orig*p_n)) * sum(
      piece_1 * In(dat$s<=u) *
        (1-Q_n(rep(C$t_0,length(dat$s)),x_long,dat$s))
    )
  }

  return(memoise2(fnc))

}



#' Construct Gamma_os_n primitive one-step estimator (based on EIF)
#'
#' @param x !!!!!
construct_Gamma_os_n <- function(dat, dat_orig, omega_n, g_n, eta_n, p_n, q_n,
                                 r_tilde_Mn, Gamma_tilde_n, vals=NA) {

  n_orig <- length(dat_orig$z)
  piece_1 <- In(dat$s!=0)
  piece_2 <- (omega_n(dat$x,dat$s,dat$y,dat$delta) /
                g_n(dat$s,dat$x)) + r_tilde_Mn(dat$s)
  piece_3 <- (1-dat_orig$weights)

  # Remove large intermediate objects
  rm(omega_n,g_n,r_tilde_Mn)

  fnc <- function(u) {

    if (F) {
      eta_n(u,c(0,0))
      return(
        (1/n_orig) * sum(
          eta_n(rep(u,n_orig),dat_orig$x)
        )
      )
    } # DEBUG

    (1/(n_orig*p_n)) * sum(dat$weights * (
      piece_1*In(dat$s<=u)*piece_2 - piece_1*Gamma_tilde_n(u)
    )) +
      (1/n_orig) * sum(
        piece_3 * (
          q_n(dat_orig$x,dat_orig$y,dat_orig$delta,u)/p_n
        ) +
          eta_n(rep(u,n_orig),dat_orig$x)
      )

  }

  return(memoise2(fnc))

}



# Estimate the variance of the Cox model based marginalized survival estimator
#'
#' @param dat_orig Dataset returned by generate_data()
#' @param dat Subsample of dataset returned by ss() for which z==1
#' @param t The end time of interest
#' @param points The A-values of interest
#' @param se_beta Estimates related to parameter vector beta_n
#' @param se_bshz Estimates related to Breslow baseline cuml. hazard estimator
#' @param se_surv Estimates related to survival at a point
#' @param se_marg Estimates related to marginalized survival
#' @param return_extras For debugging
#' @param parallel Parallelize internal lapply calls
#' @param verbose Print status updates and progress bars
#' @return A list containing the selected return values
cox_var <- function(dat_orig, dat, t, points, se_beta=F, se_bshz=F,
                    se_surv=F, se_marg=F, return_extras=F, parallel=F,
                    verbose=F) {
  # Setup
  if (verbose) {
    print(paste("Check 0 (start):", Sys.time()))
    pbapply::pboptions(type="txt", char="#", txt.width=40, style=3)
  } else {
    pbapply::pboptions(type="none")
  }
  if (parallel) {
    n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    if (verbose) {
      print(cl)
    }
  } else {
    cl <- NULL
  }

  # Fit a Cox model
  # Note: scaling the weights affects the SEs but not the estimates; thus, this
  #       is only needed for debugging
  model <- coxph(
    formula = formula(paste0("Surv(y,delta)~",
                             paste(names(dat$x),collapse="+"),"+s")),
    data = cbind(y=dat$y, delta=dat$delta,
                 dat$x, s=dat$s),
    weights = dat$weights * (length(dat$weights)/sum(dat$weights))
  )
  beta_n <- as.numeric(model$coefficients)

  if (verbose) { print(paste("Check 1 (Cox model fit):", Sys.time())) }

  # !!!!! Stabilize weights

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

    m_n <- function(x) {
      S_1n(x) / S_0n(x)
    }

    h <- function(x) {
      (S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x))
    }

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
  lstar_tilde <- (function() {
    .cache <- new.env()

    exp_terms <- list()
    for (i in c(1:max(dat_orig$strata))) {
      js <- which(ST==i)
      if (length(js)>0) {
        exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
          l_tilde(Z_[,j],Ds_[j],T_[j])
        }))
      } else {
        exp_terms[[i]] <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
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
            piece_2 <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
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
    (In(t_j<=t) * S_1n(t_j)) / (S_0n(t_j))^2
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

  # Nuisance function: vstar_n
  vstar_n <- function(st_i,d_i,t) {
    k_set <- which(ST==st_i & Ds_==1)
    if (length(k_set)>0) {
      return(
        (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
          unlist(lapply(k_set, function(k) {
            In(T_[k]<=t) / S_0n(T_[k])
          }))
        )
      )
    } else { return(0) }
  }

  # Breslow estimator
  Lambda_n <- Vectorize((function() {
    .cache <- new.env()
    function(x) {
      key <- as.character(x)
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(x) {
          (1/N) * sum(unlist(lapply(i_ev, function(i) {
            (WT[i] * In(T_[i]<=x)) / S_0n(T_[i])
          })))
        })(x)
        .cache[[key]] <- val
      }
      return(val)
    }
  })())

  # Survival estimator (at a point)
  Q_n <- (function() {
    .cache <- new.env()
    function(z) {
      key <- paste(z, collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z) {
          exp(-exp(sum(z*beta_n))*Lambda_n(t))
        })(z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: Breslow estimator (est. weights)
  infl_fn_Bres <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          pc_2 <- vstar_n(st_i,d_i,t)
          pc_4 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
            ( In(t_j<=t) * v_n(st_i,d_i,t_j) ) / (S_0n(t_j))^2
          })))
          pc_5 <- sum(mu_n*lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))

          if (d_i==1) {
            pc_1 <- ( wt_i * ds_i * In(t_i<=t) ) / S_0n(t_i)
            pc_3 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
              (In(t_j<=t)*wt_i*In(t_i>=t_j)*exp(sum(beta_n*z_i))) /
                (S_0n(t_j))^2
            })))
            return(pc_1+pc_2-pc_3-pc_4-pc_5)
          } else {
            return(pc_2-pc_4-pc_5)
          }
        })(z_i,d_i,ds_i,t_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: survival at a point (est. weights)
  omega_n <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i,z), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
          explin <- exp(sum(z*beta_n))
          piece_1 <- Lambda_n(t) * explin *
            (t(z)%*%lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))[1]
          piece_2 <- explin * infl_fn_Bres(z_i,d_i,ds_i,t_i,wt_i,st_i)
          return(-1*Q_n(z)*(piece_1+piece_2))
        })(z_i,d_i,ds_i,t_i,wt_i,st_i,z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Influence function: marginalized survival (est. weights)
  infl_fn_marg <- (function() {
    .cache <- new.env()
    x_ <- as.list(as.data.frame(t(dat_orig$x)))
    function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
      key <- paste(c(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s) {
          piece_1 <- Q_n(c(x_i,s))
          piece_2 <- (1/N) * sum(unlist(lapply(c(1:N), function(j) {
            x_j <- x_[[j]]
            omgn <- omega_n(c(x_i,s_i),d_i,ds_i,t_i,wt_i,st_i,c(x_j,s))
            return(omgn-Q_n(c(x_j,s)))
          })))
          return(piece_1+piece_2)
        })(x_i,s_i,d_i,ds_i,t_i,wt_i,st_i,s)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()

  # Construct results object
  res <- list(model=model, beta_n=beta_n)

  if (verbose) { print(paste("Check 2 (functions declared):", Sys.time())) }

  # Variance estimate: beta_hat
  if (se_beta) {

    res$est_beta <- beta_n

    if (verbose) { print(paste("Check 3a (var est START: beta):", Sys.time())) }
    res$var_est_beta <- (1/N^2) * Reduce("+", pblapply(c(1:N), function(i) {
      (lstar_tilde(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl))
    if (verbose) { print(paste("Check 4a (var est END: beta):", Sys.time())) }

  }

  # Variance estimate: Breslow estimator
  if (se_bshz) {

    res$est_bshz <- Lambda_n(t)

    if (verbose) { print(paste("Check 3b (var est START: bshz):", Sys.time())) }
    res$var_est_bshz <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (infl_fn_Bres(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4b (var est END: bshz):", Sys.time())) }

  }

  # Variance estimate: survival at a point
  if (se_surv) {

    z_0 <- c(0.3,1,0.5)

    res$est_surv <- Q_n(z_0)

    if (verbose) { print(paste("Check 3c (var est START: surv):", Sys.time())) }
    res$var_est_surv <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (omega_n(
        z_i = as.numeric(c(dat_orig$x[i,], dat_orig$s[i])),
        d_i = dat_orig$z[i],
        ds_i = dat_orig$delta[i],
        t_i = dat_orig$y[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i],
        z = z_0
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4c (var est END: surv):", Sys.time())) }

  }

  # Variance estimate: marginalized survival
  if (se_marg) {

    # !!!! basehaz vs. Lambda_hat
    bh <- basehaz(model, centered=FALSE)
    index <- max(which((bh$time<C$t_0)==T))
    est_bshz <- bh$hazard[index]
    N <- sum(dat$weights)
    res$est_marg <- unlist(lapply(points, function(s) {
      (1/N) * sum(unlist(lapply(c(1:N), function(i) {
        exp(-1*exp(sum(beta_n*c(as.numeric(dat_orig$x[i,]),s)))*est_bshz)
      })))
    }))

    # !!!!! Pre-calculate omega_n values

    if (verbose) { print(paste("Check 3d (var est START: marg):", Sys.time())) }
    res$var_est_marg <- unlist(lapply(points, function(s) {
      if (verbose) { print(paste0("Check 4d (point=",s,"): ", Sys.time())) }
      (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
        (infl_fn_marg(
          x_i = as.numeric(dat_orig$x[i,]),
          s_i = dat_orig$s[i],
          d_i = dat_orig$z[i],
          ds_i = dat_orig$delta[i],
          t_i = dat_orig$y[i],
          wt_i = dat_orig$weights[i],
          st_i = dat_orig$strata[i],
          s = s
        ))^2
      }, cl=cl)))
    }))
    if (verbose) { print(paste("Check 5d (var est END: marg):", Sys.time())) }

  }

  if (return_extras) {
    res$S_0n <- S_0n
    res$S_1n <- S_1n
    res$S_2n <- S_2n
    res$I_tilde <- I_tilde
    res$I_tilde_inv <- I_tilde_inv
    res$l_tilde <- l_tilde
    res$omega_n <- omega_n
    res$Lambda_n <- Lambda_n
  }

  return(res)

}

}
