
# Misc profiling
if (F) {

  # Old code
  microbenchmark({
    K_n1 <- (1/N) * sum((apply(dat_orig_df, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      return(Q_n(c(x_i,s)))
    })))
    K_n2 <- (1/N) * sum((apply(dat_orig_df, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      return(Q_n(c(x_i,s)) * exp(sum(c(x_i,s)*beta_n)))
    })))
    K_n3 <- (1/N) * Reduce("+", apply(dat_orig_df, 1, function(r) {
      x_i <- as.numeric(r[1:dim_x])
      return(Q_n(c(x_i,s)) * exp(sum(c(x_i,s)*beta_n)) * c(x_i,s))
    }, simplify=F))
  }, times=100L)

}

# Profiling est_cox
if (F) {

  df_v <- readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_200.rds")
  # df_v <- readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_400.rds")
  # df_v <- readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig_800.rds")
  attr(df_v, "n_orig") <- length(df_v$z)
  attr(df_v, "dim_x") <- 2 # !!!!!
  dat=list(v=df_v);
  class(dat) <- "dat_vaccine"
  library(vaccine)
  res_cox <- est_cox(dat=dat, t_0=200, cve=F)
  # s_out <- seq(min(dat$v$s,na.rm=T),max(dat$v$s,na.rm=T), l=11)
  # res_cox <- est_cox(dat=dat, t_0=200, cve=F, s_out=s_out)
  ..count


  n <- 16000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  beta1 <- runif(1)
  beta2 <- runif(1)
  beta3 <- runif(1)
  lambda_t0 <- 2

  Q_n <- function(x1,x2,x3,s) {
    exp(-1*lambda_t0*exp(beta1*x1+beta2*x2+beta3*x3))
  }

  system.time({
    mean(sapply(seq(0,1,0.01), function(s) {
      mean(sapply(c(1:n), function(i) {
        Q_n(x1[i],x2[i],x3[i],s)
      }))
    }))
  })

  system.time({
    n <- 400
    ..count <- 0
    x <- mean(sapply(c(1:101), function(s) {
      y <- rnorm(n)
      z <- rnorm(n)
      mean(sapply(y, function(i) {
        mean(sapply(z, function(j) {
          ..count <<- ..count+1
          i+j
        }))
      }))
    }))
    ..count
  })

}

# Profiling q_n
if (F) {

  microbenchmark({

  for (file in c("influence_functions", "load_data", "misc_functions",
                 "nuisance_estimators", "one_step_estimators")) {
    source(paste0("R/",file,".R"))
  }
  # Constants
  C2 <- list(n=800, alpha_1=0.5, alpha_2=0.7, alpha_3=-1, t_0=200,
             sc_params=list(lmbd=2e-4,v=1.5,lmbd2=5e-5,v2=1.5),
             appx=list(t_0=1,x_tol=25,s=0.01))

  # Generate dataset
  if (T) {
    # Covariates and exposure
    x <- data.frame(
      x1 = sample(round(seq(0,1,0.1),1), size=C2$n, replace=T),
      x2 = rbinom(C2$n, size=1, prob=0.5)
    )
    s <- round(runif(C2$n),2)

    # Survival times
    U <- runif(C2$n)
    H_0_inv <- function(t) { ((1/C2$sc_params$lmbd)*t)^(1/C2$sc_params$v) }
    lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2 + C2$alpha_3*s # Cox model
    # lin <- C2$alpha_3*expit(20*s-10) + C2$alpha_2*x$x1*x$x2 # Complex model
    t <- H_0_inv(-1*log(U)*exp(-1*lin))

    # Censoring times
    U <- runif(C2$n)
    H_0_inv2 <- function(t) { ((1/C2$sc_params$lmbd2)*t)^(1/C2$sc_params$v2) }
    lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))

    # Survival variables
    y <- pmin(t,c)
    delta <- as.integer(y==t)

    # Data structure
    expit <- function(x) { 1 / (1+exp(-x)) }
    Pi <- function(delta, y, x) {
      ev <- In(delta==1 & y<=C2$t_0)
      return(ev + (1-ev)*expit(x$x1+x$x2-1))
    }
    Pi_0 <- Pi(delta,y,x)
    z <- rbinom(C2$n, size=1, prob=Pi_0)
    weights <- z / Pi_0
    dat_orig <- list(x=x, s=ifelse(z==1,s,NA), z=z, y=y, delta=delta,
                     weights=weights)
    dat_orig$s <- round(dat_orig$s, -log10(C2$appx$s))
    dat_orig$y <- round(dat_orig$y, -log10(C2$appx$t_0))

    class(dat_orig) <- "dat_vaccine"
    attr(dat_orig, "n_orig") <- C2$n
    attr(dat_orig, "dim_x") <- 2
  }
  dat=list(v=dat_orig);

  # Run start of code
  if (T) {

    t_0=200; cve=T; cr=T; s_out=seq(0,1,0.1); ci_type="logit";
    edge_corr=F; grid_size=list(y=101,s=101,x=5); verbose=F; cf_folds=1;
    params <- list(surv_type="Cox", density_type="parametric", q_n_type="zero",
                   deriv_type="linear")
    dat_orig <- dat$v
    # Set params
    .default_params <- list(
      surv_type = "Super Learner",
      density_type = "binning",
      density_bins = 15,
      deriv_type = "m-spline",
      gamma_type = "Super Learner",
      q_n_type = "standard",
      convex_type = "GCM"
    )
    for (i in c(1:length(.default_params))) {
      p_name <- names(.default_params)[i]
      if (is.null(params[[p_name]])) { params[[p_name]] <- .default_params[[i]] }
    }
    p <- params
    p$ci_type <- ci_type
    p$cf_folds <- cf_folds
    p$edge_corr <- edge_corr

    # Rescale S to lie in [0,1] and round values
    s_min <- min(dat_orig$s, na.rm=T)
    s_max <- max(dat_orig$s, na.rm=T)
    s_shift <- -1 * s_min
    s_scale <- 1/(s_max-s_min)
    dat_orig$s <- (dat_orig$s+s_shift)*s_scale
    grid <- create_grid(dat_orig, grid_size, t_0)
    dat_orig <- round_dat(dat_orig, grid)

    # Obtain minimum value (excluding edge point mass)
    if (p$edge_corr) { s_min2 <- min(dat_orig$s[dat_orig$s!=0], na.rm=T) }

    # Rescale/round s_out and remove s_out points outside [0,1]
    s_out_orig <- s_out
    s_out <- (s_out+s_shift)*s_scale
    s_out <- sapply(s_out, function(s) { grid$s[which.min(abs(grid$s-s))] })
    na_head <- sum(s_out<0)
    na_tail <- sum(s_out>1)
    if (na_head>0) { s_out <- s_out[-c(1:na_head)] }
    len_p <- length(s_out)
    if (na_tail>0) { s_out <- s_out[-c((len_p-na_tail+1):len_p)] }

    # Create phase-two data object
    dat <- ss(dat_orig, which(dat_orig$z==1))

    # Prepare precomputation values for conditional survival estimator
    x_distinct <- dplyr::distinct(dat_orig$x)
    x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
    vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
    vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
    vals <- list(
      t = vals_pre$t,
      x = subset(vals_pre, select=-c(t,x_index,s)),
      s = vals_pre$s
    )

    # Fit conditional survival estimator
    srvSL <- construct_Q_n(p$surv_type, dat, vals)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens

    # Compute various nuisance functions
    omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
    f_sIx_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins, z1=F)
    f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    Phi_n <- construct_Phi_n(ss(dat, which(dat$s!=0))) # !!!!! Make sure Phi_n(1)=1
    n_orig <- attr(dat_orig, "n_orig")
    p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))
    eta_n <- construct_eta_n(dat, Q_n, p_n, t_0)
    r_tilde_Mn <- construct_r_tilde_Mn(dat_orig, Q_n, t_0)
    Gamma_tilde_n <- construct_Gamma_tilde_n(dat, r_tilde_Mn, p_n)
    f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)

  }

  # Construct and profile q_n
  p$q_n_type <- "standard"
  q_n <- construct_q_n(type=p$q_n_type, dat, omega_n, g_n, r_tilde_Mn,
                       Gamma_tilde_n, f_n_srv)
  dat_orig_df <- as_df(dat_orig)
  u <- 0.5; dim_x <- 2;

  # Profile this line
  q_n_do <- as.numeric(apply(dat_orig_df, 1, function(r) {
    q_n(as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]], u)
  }))

  }, times=1L)

}

# Testing est_np
if (F) {

  for (file in c("influence_functions", "load_data", "misc_functions",
                 "nuisance_estimators", "one_step_estimators")) {
    source(paste0("R/",file,".R"))
  }
  dat=list(v=readRDS("C:/Users/avike/OneDrive/Desktop/dat_orig/dat_orig.rds"));

  t_0=200; cve=T; cr=T; s_out=seq(0,1,0.1); ci_type="logit";
  edge_corr=F; grid_size=list(y=101,s=101,x=5); verbose=F; cf_folds=1;
  params <- list(surv_type="Cox", density_type="parametric", q_n_type="zero",
                 deriv_type="linear")

  # Run est_np.R through construct_Q_n()

  vals_df <- data.frame(t=vals$t, x1=vals$x$x1, x2=vals$x$x2, s=vals$s)

  # Run through construct_Gamma_os_n()

  dat_df <- as_df(dat)

  microbenchmark({
    aa <- apply(dat_df, 1, function(r) {
      omega_n(as.numeric(r[1:2]), r[["s"]], r[["y"]], r[["delta"]])
    })
  }, times=100L)

  microbenchmark({
    aa <- apply(dat_df, 1, function(r) {
      g_n(r[["s"]], as.numeric(r[1:2]))
    })
  }, times=100L)

  microbenchmark({
    aa <- sapply(dat$s, r_tilde_Mn)
  }, times=100L)

  microbenchmark({
    aa <- sapply(dat$s, Phi_n)
  }, times=100L)


}

# Testing Cox
if (F) {

  # Generate data
  if (T) {

    # Constants
    C2 <- list(
      # points = round(seq(0,1,0.02),2), # round(seq(0,1,0.1),2)
      n = 50000,
      alpha_1 = 0.5,
      alpha_2 = 0.7,
      alpha_3 = -1,
      t_0 = 200,
      sc_params = list(lmbd=2e-4, v=1.5, lmbd2=5e-5, v2=1.5),
      appx = list(t_0=1, x_tol=25, s=0.01)
    )

    # Generate dataset
    {
      # Covariates and exposure
      x <- data.frame(
        x1 = sample(round(seq(0,1,0.1),1), size=C2$n, replace=T),
        x2 = rbinom(C2$n, size=1, prob=0.5)
      )
      s <- round(runif(C2$n),2)

      # Survival times
      U <- runif(C2$n)
      H_0_inv <- function(t) { ((1/C2$sc_params$lmbd)*t)^(1/C2$sc_params$v) }
      lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2 + C2$alpha_3*s # Cox model
      # lin <- C2$alpha_3*expit(20*s-10) + C2$alpha_2*x$x1*x$x2 # Complex model
      t <- H_0_inv(-1*log(U)*exp(-1*lin))

      # Censoring times
      U <- runif(C2$n)
      H_0_inv2 <- function(t) { ((1/C2$sc_params$lmbd2)*t)^(1/C2$sc_params$v2) }
      lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2
      c <- H_0_inv2(-1*log(U)*exp(-1*lin))

      # Survival variables
      y <- pmin(t,c)
      delta <- as.integer(y==t)

      # Data structure
      expit <- function(x) { 1 / (1+exp(-x)) }
      Pi <- function(delta, y, x) {
        ev <- In(delta==1 & y<=C2$t_0)
        return(ev + (1-ev)*expit(x$x1+x$x2-1))
      }
      Pi_0 <- Pi(delta,y,x)
      z <- rbinom(C2$n, size=1, prob=Pi_0)
      weights <- z / Pi_0
      dat_orig <- list(x=x, s=ifelse(z==1,s,NA), z=z, y=y, delta=delta,
                       weights=weights)
      dat_orig$s <- round(dat_orig$s, -log10(C2$appx$s))
      dat_orig$y <- round(dat_orig$y, -log10(C2$appx$t_0))
    }

  }

  dat <- ss(dat_orig, which(dat_orig$z==1))

  ..p <<- list() # !!!!!
  model_srv <- survival::coxph(
    formula = formula(paste0("survival::Surv(y,delta)~",
                             paste(names(dat$x),collapse="+"),"+s")),
    data = cbind(y=dat$y, delta=dat$delta, dat$x, s=dat$s),
    weights = dat$weights
  )
  coeffs_srv <- as.numeric(model_srv$coefficients)
  coeffs_srv
  c(alpha_1=0.5, alpha_2=0.7, alpha_3=-1)
  # model_srv

  ..p$coeffs_srv <<- coeffs_srv # !!!!!
  bh_srv <- survival::basehaz(model_srv, centered=FALSE)

  model_cens <- survival::coxph(
    formula = formula(paste0("survival::Surv(y,delta)~",
                             paste(names(dat$x),collapse="+"),"+s")),
    data = cbind(y=dat$y, delta=1-dat$delta, dat$x, s=dat$s),
    weights = dat$weights
  )
  coeffs_cens <- as.numeric(model_cens$coefficients)
  ..p$coeffs_cens <<- coeffs_cens # !!!!!
  bh_cens <- survival::basehaz(model_cens, centered=FALSE)

  fnc_srv <- function(t, x, s) {
    if (t==0) {
      return(1)
    } else {
      # Lambda_t <- bh_srv$hazard[which.min(abs(bh_srv$time-t))]
      Lambda_t <- bh_srv$hazard[max(which((bh_srv$time<t)==T))]
      return(exp(-1*Lambda_t*exp(sum(coeffs_srv*as.numeric(c(x,s))))))
    }
  }

  fnc_cens <- function(t, x, s) {
    if (t==0) {
      return(1)
    } else {
      Lambda_t <- bh_cens$hazard[max(which((bh_cens$time<t)==T))]
      return(exp(-1*Lambda_t*exp(sum(coeffs_cens*as.numeric(c(x,s))))))
    }
  }


}

# Comparing two derivative estimators
if (F) {

  # Linear
  # points_x: 0.000 0.195 0.295 0.695 0.895 1.000
  # points_y: 0.00 0.10 0.40 0.75 0.95 1.00
  # points_sl: 0.513 3.000 0.875 1.000 0.476

  # M-Spline
  # jump_points: 0.000 0.195 0.295 0.695 0.895 1.000
  # midpoints: 0.098 0.245 0.495 0.795 0.948

  # !!!!!
  u <- runif(5)
  ecd <- ecdf(u)
  r_Mn <- Vectorize(function(x) {
    return(1-ecd(x))
    # p1 <- 0.2*In(x>=0.2)
    # p2 <- 0.4*In(x>=0.3)
    # p3 <- 0.3*In(x>=0.7)
    # p4 <- 0.1*In(x>=0.9)
    # return(1-(p1+p2+p3+p4))
  })

  if (T) {
    # Estimate entire function on grid
    r_Mns <- r_Mn(grid$s)

    # grid_width <- grid$s[2] - grid$s[1]
    points_x <- grid$s[1]
    points_y <- r_Mns[1]
    for (i in 2:length(grid$s)) {
      if (r_Mns[i]-r_Mns[round(i-1)]!=0) {
        points_x <- c(points_x, (grid$s[i]+grid$s[round(i-1)])/2)
        # points_x <- c(points_x, grid$s[i]-(grid_width/2))
        points_y <- c(points_y, mean(c(r_Mns[i],r_Mns[round(i-1)])))
      }
    }
    points_x <- c(points_x, grid$s[length(grid$s)])
    points_y <- c(points_y, r_Mns[length(grid$s)])
    fnc_pre <- approxfun(x=points_x, y=points_y, method="linear", rule=2)
    fnc_pre_lin <- fnc_pre # !!!!!
    fnc_pre <- splinefun(x=points_x, y=points_y, method="monoH.FC")
    fnc_pre_spl <- fnc_pre # !!!!!
    fnc_lin <- function(s) {
      width <- 0.1 # !!!!! Changed from 0.2
      x1 <- s - width/2
      x2 <- s + width/2
      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }
      return(min((fnc_pre_lin(x2)-fnc_pre_lin(x1))/width,0))
    }
    fnc_spl <- function(s) {
      width <- 0.1 # !!!!! Changed from 0.2
      x1 <- s - width/2
      x2 <- s + width/2
      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }
      return(min((fnc_pre_spl(x2)-fnc_pre_spl(x1))/width,0))
    }
  }

  df_plot <- data.frame(
    x = c(grid$s, grid$s, grid$s),
    y = c(sapply(grid$s,r_Mn),
          sapply(grid$s,fnc_pre_lin),
          sapply(grid$s,fnc_pre_spl)),
    grp = c(rep("r_Mn",length(r_Mns)),
            rep("fnc_pre_lin",length(r_Mns)),
            rep("fnc_pre_spl",length(r_Mns)))
  )
  ggplot(df_plot, aes(x=x, y=y, color=grp)) + geom_line()

  df_plot <- data.frame(
    x = c(grid$s, grid$s),
    y = c(sapply(grid$s,fnc_lin),
          sapply(grid$s,fnc_spl)),
    grp = c(rep("fnc_lin",length(r_Mns)),
            rep("fnc_spl",length(r_Mns)))
  )
  ggplot(df_plot, aes(x=x, y=y, color=grp)) + geom_line()

}

# Reconciling
if (F) {

  # 640
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

  # 694
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

}

# Benchmarking ways to apply over vectors
if (F) {

  ss <- rnorm(20000)

  microbenchmark({
    yy <- apply(as.array(ss), 1, function(s) { s^2 })
  }, times=1000L)

  microbenchmark({
    yy <- sapply(ss, function(s) { s^2 })
  }, times=1000L)

  microbenchmark({
    yy <- unlist(lapply(ss, function(s) { s^2 }))
  }, times=1000L)

}



# Benchmarking ways to apply over rows
if (F) {

  x1 <- rnorm(20000)
  x2 <- runif(20000)
  x3 <- rexp(20000)
  xx <- data.frame(x1=x1, x2=x2, x3=x3)

  microbenchmark({
    yy <- apply(xx, 1, function(x) {
      x <- as.numeric(x)
      x[1]*x[2]*x[3]
    })
  }, times=10L)

  microbenchmark({
    yy <- sapply(c(1:nrow(xx)), function(i) {
      x <- as.numeric(xx[i,])
      x[1]*x[2]*x[3]
    })
  }, times=10L)

  microbenchmark({
    yy <- unlist(lapply(c(1:nrow(xx)), function(i) {
      x <- as.numeric(xx[i,])
      x[1]*x[2]*x[3]
    }))
  }, times=10L)

}



# New superfunc
if (F) {

  # Old superfunc

  # New superfunc
  # Note: memoising will not work if the function itself returns NULL
  construct_superfunc_new <- function(fnc) {

    htab <- new.env()
    ..new_fnc <- function() {
      ..e <- parent.env(environment())
      ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
      keylist <- lapply(..e$arg_names, function(arg_name) {
        as.numeric(..mc[[arg_name]])
      })
      key <- paste(keylist, collapse=";")
      val <- ..e$htab[[key]]
      if (is.null(val)) {
        val <- do.call(..e$fnc, ..mc)
        ..e$htab[[key]] <- val
      }
      res <- val
      return(res)
    }

    # Set formals and set up environment
    formals(..new_fnc) <- formals(fnc)
    f_env <- new.env(parent=environment(fnc))
    f_env$arg_names <- names(formals(fnc))
    f_env$htab <- htab
    f_env$fnc <- fnc
    environment(..new_fnc) <- f_env

    return(..new_fnc)

  }

  memoise2 <- function(fnc) {

    htab <- new.env()
    ..new_fnc <- function() {
      ..e <- parent.env(environment())
      ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
      key <- rlang::hash(..mc)
      val <- ..e$htab[[key]]
      if (is.null(val)) {
        val <- do.call(..e$fnc, ..mc)
        ..e$htab[[key]] <- val
      }
      return(val)
    }

    # Set formals and set up environment
    formals(..new_fnc) <- formals(fnc)
    f_env <- new.env(parent=environment(fnc))
    f_env$arg_names <- names(formals(fnc))
    f_env$htab <- htab
    f_env$fnc <- fnc
    environment(..new_fnc) <- f_env

    return(..new_fnc)

  }

  # microbenchmark
  fnc <- function(a,b,c) { a * (b[1]+b[2]) * c }
  fnc_old <- construct_superfunc(fnc, aux=NA, vec=c(1,1,1), vals=NA)
  fnc_new <- construct_superfunc_new(fnc)
  fnc_new2 <- construct_superfunc_new2(fnc)
  fnc_mem <- memoise(fnc)
  a <- round(runif(4),5)
  m1 <- microbenchmark({ fnc(a[1], c(a[2],a[3]), a[4]) }, times=10000L)
  m2 <- microbenchmark({ fnc_old(a[1], c(a[2],a[3]), a[4]) }, times=10000L)
  m3 <- microbenchmark({ fnc_new(a[1], c(a[2],a[3]), a[4]) }, times=10000L)
  m4 <- microbenchmark({ fnc_new2(a[1], c(a[2],a[3]), a[4]) }, times=10000L)
  m5 <- microbenchmark({ fnc_mem(a[1], c(a[2],a[3]), a[4]) }, times=10000L)
  print(paste("Orig:", round(mean(m1$time)/1000), "ms"))
  print(paste("Old :", round(mean(m2$time)/1000), "ms"))
  print(paste("New :", round(mean(m3$time)/1000), "ms"))
  print(paste("New2:", round(mean(m4$time)/1000), "ms"))
  print(paste("Mem :", round(mean(m5$time)/1000), "ms"))

}



# Generate fake data
if (F) {

  # Constants
  C2 <- list(
    # points = round(seq(0,1,0.02),2), # round(seq(0,1,0.1),2)
    n = 500,
    alpha_1 = 0.5,
    alpha_2 = 0.7,
    alpha_3 = -1,
    t_0 = 200,
    sc_params = list(lmbd=2e-4, v=1.5, lmbd2=5e-5, v2=1.5),
    appx = list(t_0=1, x_tol=25, s=0.01)
  )

  # Generate dataset
  {
    # Covariates and exposure
    x <- data.frame(
      x1 = sample(round(seq(0,1,0.1),1), size=C2$n, replace=T),
      x2 = rbinom(C2$n, size=1, prob=0.5)
    )
    s <- round(runif(C2$n),2)

    # Survival times
    U <- runif(C2$n)
    H_0_inv <- function(t) { ((1/C2$sc_params$lmbd)*t)^(1/C2$sc_params$v) }
    lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2 + C2$alpha_3*s # Cox model
    # lin <- C2$alpha_3*expit(20*s-10) + C2$alpha_2*x$x1*x$x2 # Complex model
    t <- H_0_inv(-1*log(U)*exp(-1*lin))

    # Censoring times
    U <- runif(C2$n)
    H_0_inv2 <- function(t) { ((1/C2$sc_params$lmbd2)*t)^(1/C2$sc_params$v2) }
    lin <- C2$alpha_1*x$x1 + C2$alpha_2*x$x2
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))

    # Survival variables
    y <- pmin(t,c)
    delta <- as.integer(y==t)

    # Data structure
    expit <- function(x) { 1 / (1+exp(-x)) }
    Pi <- function(delta, y, x) {
      ev <- In(delta==1 & y<=C2$t_0)
      return(ev + (1-ev)*expit(x$x1+x$x2-1))
    }
    Pi_0 <- Pi(delta,y,x)
    z <- rbinom(C2$n, size=1, prob=Pi_0)
    weights <- z / Pi_0
    dat_orig <- list(x=x, s=ifelse(z==1,s,NA), z=z, y=y, delta=delta,
                     weights=weights)
    dat_orig$s <- round(dat_orig$s, -log10(C2$appx$s))
    dat_orig$y <- round(dat_orig$y, -log10(C2$appx$t_0))

    class(dat_orig) <- "dat_vaccine"
    attr(dat_orig, "n_orig") <- C2$n
    attr(dat_orig, "dim_x") <- 2
  }

}
