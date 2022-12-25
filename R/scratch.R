
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
    s <- runif(C2$n)

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
    dat_orig <- list(x=x, s=s, z=rep(1,C2$n), y=y, delta=delta)
    dat_orig$s <- round(dat_orig$s, -log10(C2$appx$s))
    dat_orig$y <- round(dat_orig$y, -log10(C2$appx$t_0))

  }

}
