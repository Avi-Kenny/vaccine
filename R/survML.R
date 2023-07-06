#################################.
##### compute_exponential.R #####
#################################.

#' Compute the plug-in exponential survival function estimator
#'
#' @param cdf_uncens Matrix of predictions of the cdf of the uncensored times (F_Y_1) on chosen time grid
#' @param cdf_cens Predictions of the cdf of the censored times (F_Y_0) on chosen time grid
#' @param p_uncens Prediction of the probability of being uncensored
#' @param newtimes Times at which to make the prediction
#' @param time_grid Grid of time points over which to discretize the integral
#' @param denom_method Method of computing the denominator
#'
#' @return A vector of estimates of the survival function over \code{time_grid}
#'
#' @noRd
compute_exponential <- function(cdf_uncens,
                                cdf_cens = NA,
                                cdf_marg = NA,
                                entry_uncens = NA,
                                entry_cens = NA,
                                entry_marg = NA,
                                p_uncens,
                                newtimes,
                                time_grid,
                                denom_method = "stratified",
                                truncation = TRUE){

  estimate_S_T <- function(t){
    curr_length <- sum(time_grid <= t)

    # get S_Y estimates up to t
    S_Y_1_curr <- cdf_uncens[1:curr_length]

    dF_Y_1_pred <- c(S_Y_1_curr[1], diff(S_Y_1_curr))

    S_Y_1_pred_left <- c(1, 1-S_Y_1_curr[-length(S_Y_1_curr)])# probability of being "at risk" at time t
    ### CHECK TO MAKE SURE THIS IS CORRECT WITH THE DISCRETIZATION OF TIME

    if (!truncation){ # truncation
      if (denom_method != "stratified"){# marginal
        S_Y_curr <- cdf_marg[1:curr_length]
        S_Y_pred_left <- c(1, 1-S_Y_curr[-length(S_Y_curr)])
        low <- S_Y_pred_left
      } else{
        S_Y_0_curr <- cdf_cens[1:curr_length]
        S_Y_0_pred_left <- c(1, 1-S_Y_0_curr[-length(S_Y_0_curr)])# probability of being "at risk" at time t
        low_right <- S_Y_0_pred_left * (1 - p_uncens)
        low_left <- S_Y_1_pred_left * p_uncens
      }


      # product form
      if (denom_method == "stratified"){
        #print(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))
        S_T_est <- exp(-sum(p_uncens * dF_Y_1_pred/(low_left + low_right)))
      } else{
        S_T_est <- exp(-sum(p_uncens * dF_Y_1_pred/low))
      }
    } else{ # if there is truncation
      if (denom_method != "stratified"){
        S_Y_curr <- cdf_marg[1:curr_length]
        F_W_curr <- entry_marg[1:curr_length]
        S_Y_pred_left <- c(1, 1-S_Y_curr[-length(S_Y_curr)])
        low <- S_Y_pred_left * F_W_curr
      } else{ # marginal
        S_Y_0_curr <- cdf_cens[1:curr_length]
        S_Y_0_pred_left <- c(1, 1-S_Y_0_curr[-length(S_Y_0_curr)])# probability of being "at risk" at time t
        F_W_1_curr <- entry_uncens[1:curr_length]
        F_W_0_curr <-entry_cens[1:curr_length]
        low_right <- F_W_0_curr * S_Y_0_pred_left * (1 - p_uncens)
        low_left <- F_W_1_curr * S_Y_1_pred_left * p_uncens
      }


      # product form
      if (denom_method == "stratified"){
        #print(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))
        S_T_est <- exp(-sum(p_uncens * dF_Y_1_pred/(low_left + low_right)))
      } else{
        S_T_est <- exp(-sum(p_uncens * dF_Y_1_pred/low))
      }
    }

    if (curr_length == 0){
      S_T_est <- 1
    }
    return(S_T_est)
  }

  S_T_ests <- apply(X = as.matrix(newtimes),
                    MARGIN = 1,
                    FUN = estimate_S_T)

  return(S_T_ests)
}



#############################.
##### compute_prodint.R #####
#############################.

#' Compute the plug-in product-limit survival function estimator
#'
#' @param cdf_uncens Matrix of predictions of the cdf of the uncensored times (F_Y_1) on chosen time grid
#' @param cdf_cens Predictions of the cdf of the censored times (F_Y_0) on chosen time grid
#' @param p_uncens Prediction of the probability of being uncensored
#' @param newtimes Times at which to make the prediction
#' @param time_grid Grid of time points over which to discretize the product integral
#' @param denom_method Method of computing the denominator
#'
#' @return A vector of estimates of the survival function over \code{time_grid}
#'
#' @noRd
compute_prodint <- function(cdf_uncens,
                            cdf_cens = NA,
                            cdf_marg = NA,
                            entry_uncens = NA,
                            entry_cens = NA,
                            entry_marg = NA,
                            p_uncens,
                            newtimes,
                            time_grid,
                            denom_method = "stratified",
                            truncation = TRUE){

  estimate_S_T <- function(t){
    curr_length <- sum(time_grid <= t)

    # get S_Y estimates up to t
    S_Y_1_curr <- cdf_uncens[1:curr_length]


    dF_Y_1_pred <- c(S_Y_1_curr[1], diff(S_Y_1_curr))

    S_Y_1_pred_left <- c(1, 1-S_Y_1_curr[-length(S_Y_1_curr)])# probability of being "at risk" at time t
    ### CHECK TO MAKE SURE THIS IS CORRECT WITH THE DISCRETIZATION OF TIME

    if (!truncation){ # truncation
      if (denom_method != "stratified"){# marginal
        S_Y_curr <- cdf_marg[1:curr_length]
        S_Y_pred_left <- c(1, 1-S_Y_curr[-length(S_Y_curr)])
        low <- S_Y_pred_left
      } else{
        S_Y_0_curr <- cdf_cens[1:curr_length]
        S_Y_0_pred_left <- c(1, 1-S_Y_0_curr[-length(S_Y_0_curr)])# probability of being "at risk" at time t
        low_right <- S_Y_0_pred_left * (1 - p_uncens)
        low_left <- S_Y_1_pred_left * p_uncens
      }

      # product form
      if (denom_method == "stratified"){
        S_T_est <- prod(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))
      } else{
        S_T_est <- prod(1 - p_uncens * dF_Y_1_pred/low)
      }
    } else{ # if there is truncation
      if (denom_method != "stratified"){
        S_Y_curr <- cdf_marg[1:curr_length]
        F_W_curr <- entry_marg[1:curr_length]
        S_Y_pred_left <- c(1, 1-S_Y_curr[-length(S_Y_curr)])
        low <- S_Y_pred_left * F_W_curr
      } else{ # marginal
        S_Y_0_curr <- cdf_cens[1:curr_length]
        S_Y_0_pred_left <- c(1, 1-S_Y_0_curr[-length(S_Y_0_curr)])# probability of being "at risk" at time t
        F_W_1_curr <- entry_uncens[1:curr_length]
        F_W_0_curr <-entry_cens[1:curr_length]
        low_right <- F_W_0_curr * S_Y_0_pred_left * (1 - p_uncens)
        low_left <- F_W_1_curr * S_Y_1_pred_left * p_uncens
      }


      # product form
      if (denom_method == "stratified"){
        #print(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))
        S_T_est <- prod(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))
      } else{
        S_T_est <- prod(1 - p_uncens * dF_Y_1_pred/low)
      }
    }

    if (curr_length == 0){
      S_T_est <- 1
    }

    return(S_T_est)
  }

  S_T_ests <- apply(X = as.matrix(newtimes),
                    MARGIN = 1,
                    FUN = estimate_S_T)

  return(S_T_ests)
}



############################.
##### f_w_algorithms.R #####
############################.

#' Stacked binary regression using SuperLearner
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_w_stack_SuperLearner}
#' @noRd
f_w_stack_SuperLearner <- function(time,
                                   event,
                                   entry,
                                   X,
                                   censored,
                                   bin_size,
                                   SL_control,
                                   time_basis,
                                   direction){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
      obsWeights <- SL_control$obsWeights[!as.logical(event)]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      entry <- entry[as.logical(event)]
      obsWeights <- SL_control$obsWeights[as.logical(event)]
    }
  } else{
    time <- time
    entry <- entry
    X <- X
    obsWeights <- SL_control$obsWeights
  }

  cv_folds <- split(sample(1:length(time)), rep(1:SL_control$V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)


  if (!is.null(bin_size)){
    #time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid <- sort(unique(stats::quantile(time, probs = seq(0, 1, by = bin_size))))
    time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
    #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)

  }

  ids <- seq(1:length(time))

  if (!is.null(obsWeights)){
    stackX <- as.matrix(data.frame(X,
                                   obsWeights = obsWeights,
                                   ids = ids))
  } else{
    stackX <- as.matrix(data.frame(X,
                                   ids = ids))
  }

  stacked <- stack_entry(time = time,
                         entry = entry,
                         X = stackX,
                         time_grid = time_grid,
                         time_basis = "continuous")

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }

  long_obsWeights <- stacked$obsWeights
  stacked_ids <- stacked$ids
  stacked$obsWeights <- NULL
  stacked$ids <- NULL
  .Y <- stacked[,ncol(stacked)]
  .X <- stacked[,-ncol(stacked)]

  get_validRows <- function(fold_sample_ids){
    validRows <- which(stacked_ids %in% fold_sample_ids)
    return(validRows)
  }

  validRows <- lapply(cv_folds, get_validRows)

  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  if (is.null(SL_control$stratifyCV)){
    SL_control$stratifyCV <- FALSE
  }

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL_control$SL.library,
                                    family = stats::binomial(),
                                    method = SL_control$method,
                                    verbose = FALSE,
                                    obsWeights = long_obsWeights,
                                    cvControl = list(V = SL_control$V,
                                                     validRows = validRows,
                                                     stratifyCV = SL_control$stratifyCV))

  fit <- list(reg.object = fit, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_SuperLearner")
  return(fit)
}

#' Prediction function for stacked SuperLearner
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_SuperLearner <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
  } else if (fit$time_basis == "dummy"){
    get_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1)
  }

  return(predictions)
}



#######################.
##### f_w_stack.R #####
#######################.

#' Wrapper for various f_w stacked algorithms
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param entry Study entry variable, if applicable. Defaults to \code{NULL},
#' indicating that there is no truncation.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param censored Logical, indicates whether to run regression on censored
#' observations (\code{event == 0}) vs. uncensored (\code{event == 1}).
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s.
#' @param time_basis How to treat time for training the binary
#' classifier. Options are \code{"continuous"} and \code{"dummy"}, meaning
#' an indicator variable is included for each time in the time grid.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return An fitted pooled binary regression for the truncation distribution
#'
#' @noRd
f_w_stack <- function(time,
                      event,
                      X,
                      entry,
                      censored,
                      bin_size,
                      learner = "SuperLearner",
                      SL_control,
                      xgb_control,
                      time_basis){

  if (learner == "SuperLearner"){
    fit <- f_w_stack_SuperLearner(time = time,
                                  event = event,
                                  X = X,
                                  censored = censored,
                                  bin_size = bin_size,
                                  SL_control = SL_control,
                                  time_basis = time_basis,
                                  entry = entry)

  }

  return(fit)
}



############################.
##### f_y_algorithms.R #####
############################.

#' Stacked binary regression with SuperLearner, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param time_basis How to treat time
#' @param SL_control SuperLearner control parameters
#'
#' @return An object of class \code{f_y_stack_SuperLearner}
#' @noRd
f_y_stack_SuperLearner <- function(time,
                                   event,
                                   X,
                                   censored,
                                   bin_size,
                                   isotonize = TRUE,
                                   SL_control,
                                   time_basis){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      X <- X[!as.logical(event),]
      obsWeights <- SL_control$obsWeights[!as.logical(event)]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      obsWeights <- SL_control$obsWeights[as.logical(event)]
    }
  } else{
    time <- time
    X <- X
    obsWeights <- SL_control$obsWeights
  }

  cv_folds <- split(sample(1:length(time)), rep(1:SL_control$V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  dat <- data.frame(X, time)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- sort(unique(stats::quantile(dat$time, probs = seq(0, 1, by = bin_size))))
    time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
    #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  ids <- seq(1:length(time))

  if (!is.null(obsWeights)){
    stackX <- as.matrix(data.frame(X,
                                   obsWeights = obsWeights,
                                   ids = ids))
  } else{
    stackX <- as.matrix(data.frame(X,
                                   ids = ids))
  }

  stacked <- stack_cdf(time = time,
                       X = stackX,
                       time_grid = time_grid,
                       time_basis = "continuous")

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }

  long_obsWeights <- stacked$obsWeights
  stacked_ids <- stacked$ids
  stacked$obsWeights <- NULL
  stacked$ids <- NULL
  .Y <- stacked[,ncol(stacked)]
  .X <- stacked[,-ncol(stacked)]

  get_validRows <- function(fold_sample_ids){
    validRows <- which(stacked_ids %in% fold_sample_ids)
    return(validRows)
  }

  validRows <- lapply(cv_folds, get_validRows)

  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  if (is.null(SL_control$stratifyCV)){
    SL_control$stratifyCV <- FALSE
  }

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL_control$SL.library,
                                    family = stats::binomial(),
                                    method = SL_control$method,
                                    verbose = FALSE,
                                    obsWeights = long_obsWeights,
                                    cvControl = list(V = SL_control$V,
                                                     validRows = validRows,
                                                     stratifyCV = SL_control$stratifyCV))

  fit <- list(reg.object = fit,
              time_grid = time_grid,
              isotonize = isotonize,
              time_basis = time_basis)
  class(fit) <- c("f_y_stack_SuperLearner")
  return(fit)
}

#' Prediction function for stacked SuperLearner CDF
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stack_SuperLearner <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
  } else{
    get_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1)
  }

  if (fit$isotonize){
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  return(iso.cdf.ests)
}



#######################.
##### f_y_stack.R #####
#######################.

#' Wrapper for various f_y stacked algorithms
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param censored Logical, indicates whether to run regression on censored
#' observations (\code{event == 0}) vs. uncensored (\code{event == 1}).
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s.
#' @param time_basis How to treat time for training the binary
#' classifier. Options are \code{"continuous"} and \code{"dummy"}, meaning
#' an indicator variable is included for each time in the time grid.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return An fitted pooled binary regression for the CDF
#'
#' @noRd
f_y_stack <- function(time,
                      event,
                      X,
                      censored,
                      bin_size,
                      learner = "SuperLearner",
                      SL_control,
                      xgb_control,
                      time_basis){

  if (learner == "SuperLearner"){
    fit <- f_y_stack_SuperLearner(time = time,
                                  event = event,
                                  X = X,
                                  censored = censored,
                                  bin_size = bin_size,
                                  time_basis = time_basis,
                                  SL_control = SL_control)
  }



  return(fit)
}



#####################.
##### p_delta.R #####
#####################.

#' Wrapper for various p_delta algorithms
#'
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return An fitted binary regression for (complement of)
#' probability of censoring
#'
#' @noRd
p_delta <- function(event,
                    X,
                    learner = "SuperLearner",
                    SL_control,
                    xgb_control){
  if (learner == "SuperLearner"){
    fit <- p_delta_SuperLearner(event = event,
                                X = X,
                                SL_control = SL_control)
  }

  return(fit)
}



################################.
##### p_delta_algorithms.R #####
################################.

#' Binary SuperLearner
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param SL_control Super Learner control parameters
#'
#' @return An object of class \code{p_delta_SuperLearner}
#' @noRd
p_delta_SuperLearner <- function(event,
                                 X,
                                 SL_control){

  X <- as.data.frame(X)

  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  if (is.null(SL_control$stratifyCV)){
    SL_control$stratifyCV <- FALSE
  }

  opt_fit <- SuperLearner::SuperLearner(Y = event,
                                        X = X,
                                        family = stats::binomial(),
                                        SL.library = SL_control$SL.library,
                                        method = SL_control$method,
                                        verbose = FALSE,
                                        cvControl = list(V = SL_control$V,
                                                         stratifyCV = SL_control$stratifyCV),
                                        obsWeights = SL_control$obsWeights)

  fit <- list(reg.object = opt_fit)
  class(fit) <- c("p_delta_SuperLearner")
  return(fit)
}

#' Prediction function for p delta SuperLearner
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_SuperLearner <- function(fit,
                                         newX){
  X <- as.data.frame(newX)
  preds <- stats::predict(fit$reg.object, newdata = newX)$pred
  return(preds)
}



####################.
##### stackG.R #####
####################.

#' Estimate a conditional survival function using global survival stacking
#'
#' @noRd
stackG <- function(time,
                   event = rep(1, length(time)),
                   entry = NULL,
                   X,
                   newX = NULL,
                   newtimes = NULL,
                   direction = "prospective",
                   bin_size = NULL,
                   time_basis,
                   time_grid_approx = sort(unique(time)),
                   surv_form = "PI",
                   learner = "SuperLearner",
                   SL_control = list(SL.library = c("SL.mean"),
                                     V = 10,
                                     method = "method.NNLS",
                                     stratifyCV = FALSE),
                   xgb_control = list(tuning_params = list(ntrees = 500,
                                                           max_depth = 2,
                                                           eta = 0.01,
                                                           subsample = 1),
                                      V = 10,
                                      objective = "binary:logistic",
                                      eval_metric = "logloss"),
                   tau = NULL){
  P_Delta_opt <- NULL
  S_Y_opt <- NULL
  F_Y_1_opt <- NULL
  F_Y_0_opt <- NULL
  G_W_1_opt <- NULL
  G_W_0_opt <- NULL
  F_W_opt <- NULL

  tau <- NULL

  if (is.null(newX)){
    newX <- X
  }

  if (is.null(newtimes)){
    newtimes <- time_grid_approx
  }

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    entry <- tau - entry
    event <- rep(1, length(time))
    newtimes <- tau - newtimes
    time_grid_approx <- sort(tau - time_grid_approx)
    P_Delta_opt_preds <- rep(1, nrow(newX))
    F_Y_0_opt_preds <- matrix(0, nrow = nrow(newX), ncol = length(time_grid_approx))
    G_W_0_opt_preds <- matrix(0, nrow = nrow(newX), ncol = length(time_grid_approx))
  }

  # if there is a censoring probability to estimate, i.e. if there is censoring
  if (sum(event == 0) != 0){
    P_Delta_opt <- p_delta(event = event,
                           X = X,
                           learner = learner,
                           SL_control = SL_control,
                           xgb_control = xgb_control)
    P_Delta_opt_preds <- stats::predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms

    F_Y_0_opt <- f_y_stack(time = time,
                           event = event,
                           X = X,
                           censored = TRUE,
                           bin_size = bin_size,
                           learner = learner,
                           SL_control = SL_control,
                           xgb_control = xgb_control,
                           time_basis = time_basis)
    F_Y_0_opt_preds <- stats::predict(F_Y_0_opt,
                                      newX = newX,
                                      newtimes = time_grid_approx)
  }

  F_Y_1_opt <- f_y_stack(time = time,
                         event = event,
                         X = X,
                         censored = FALSE,
                         bin_size = bin_size,
                         learner = learner,
                         SL_control = SL_control,
                         xgb_control = xgb_control,
                         time_basis = time_basis)

  if (!is.null(entry)){ # if a truncation variable is given
    G_W_1_opt <- f_w_stack(time = time,
                           event = event,
                           X = X,
                           censored = FALSE,
                           bin_size = bin_size,
                           learner = learner,
                           SL_control = SL_control,
                           xgb_control = xgb_control,
                           entry = entry,
                           time_basis = time_basis)
    G_W_1_opt_preds <- stats::predict(G_W_1_opt,
                                      newX = newX,
                                      newtimes = time_grid_approx)
    if (sum(event == 0) != 0){ # if there's censoring
      G_W_0_opt <- f_w_stack(time = time,
                             event = event,
                             X = X,
                             censored = TRUE,
                             bin_size = bin_size,
                             learner = learner,
                             SL_control = SL_control,
                             xgb_control = xgb_control,
                             entry = entry,
                             time_basis = time_basis)
      G_W_0_opt_preds <- stats::predict(G_W_0_opt,
                                        newX = newX,
                                        newtimes = time_grid_approx)
    } else{
      G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
    }
  } else{
    G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
    G_W_1_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
  }

  F_Y_1_opt_preds <- stats::predict(F_Y_1_opt,
                                    newX = newX,
                                    newtimes = time_grid_approx)

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    F_Y_1_curr <- F_Y_1_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]
    F_Y_0_curr <- F_Y_0_opt_preds[i,]
    G_W_0_curr <- G_W_0_opt_preds[i,]
    G_W_1_curr <- G_W_1_opt_preds[i,]
    if (surv_form == "PI"){
      S_T_ests <-compute_prodint(cdf_uncens = F_Y_1_curr,
                                 cdf_cens = F_Y_0_curr,
                                 entry_uncens = G_W_1_curr,
                                 entry_cens = G_W_0_curr,
                                 p_uncens = pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
      S_C_ests <-compute_prodint(cdf_uncens = F_Y_0_curr,
                                 cdf_cens = F_Y_1_curr,
                                 entry_uncens = G_W_0_curr,
                                 entry_cens = G_W_1_curr,
                                 p_uncens = 1 - pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
    } else if (surv_form == "exp"){
      S_T_ests <-compute_exponential(cdf_uncens = F_Y_1_curr,
                                     cdf_cens = F_Y_0_curr,
                                     entry_uncens = G_W_1_curr,
                                     entry_cens = G_W_0_curr,
                                     p_uncens = pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
      S_C_ests <-compute_exponential(cdf_uncens = F_Y_0_curr,
                                     cdf_cens = F_Y_1_curr,
                                     entry_uncens = G_W_0_curr,
                                     entry_cens = G_W_1_curr,
                                     p_uncens = 1 - pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
    }
    return(list(S_T_ests = S_T_ests, S_C_ests = S_C_ests))
  }



  preds <- t(matrix(unlist(apply(X = as.matrix(seq(1, nrow(newX))),
                                 MARGIN = 1,
                                 FUN = estimate_S_T)), nrow = 2*length(newtimes)))

  S_T_preds <- preds[,1:length(newtimes)]
  S_C_preds <- preds[,(length(newtimes) + 1):(2*length(newtimes))]

  if (direction == "retrospective"){
    S_T_preds <- 1 - S_T_preds
    S_C_preds <- NULL
  }

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = S_C_preds,
              time_grid_approx = time_grid_approx,
              direction = direction,
              tau = tau,
              surv_form = surv_form,
              time_basis = time_basis,
              learner = learner,
              SL_control = SL_control,
              xgb_control = xgb_control,
              fits = list(P_Delta = P_Delta_opt,
                          F_Y_1 = F_Y_1_opt,
                          F_Y_0 = F_Y_0_opt,
                          G_W_1 = G_W_1_opt,
                          G_W_0 = G_W_0_opt))
  class(res) <- "stackG"
  return(res)
}


#' @noRd
#' @export
predict.stackG <- function(object,
                           newX,
                           newtimes,
                           surv_form = "PI"){

  if (object$direction == "retrospective"){
    newtimes <- object$tau - newtimes
  }

  if (!is.null(object$fits$P_Delta)){
    P_Delta_opt_preds <- stats::predict(object$fits$P_Delta, newX = newX)
  } else{
    P_Delta_opt_preds <- rep(1, nrow(newX))
  }
  if (!is.null(object$fits$G_W_1)){
    G_W_1_opt_preds <- stats::predict(object$fits$G_W_1,
                                      newX = newX,
                                      newtimes = object$time_grid_approx)
  } else{
    G_W_1_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(object$time_grid_approx))
  }
  if (!is.null(object$fits$G_W_0)){
    G_W_0_opt_preds <- stats::predict(object$fits$G_W_0,
                                      newX = newX,
                                      newtimes = object$time_grid_approx)
  } else{
    G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(object$time_grid_approx))
  }
  if (!is.null(object$fits$F_Y_0)){
    F_Y_0_opt_preds <- stats::predict(object$fits$F_Y_0,
                                      newX = newX,
                                      newtimes = object$time_grid_approx)
  } else{
    F_Y_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(object$time_grid_approx))
  }

  F_Y_1_opt_preds <- stats::predict(object$fits$F_Y_1,
                                    newX = newX,
                                    newtimes = object$time_grid_approx)

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    F_Y_1_curr <- F_Y_1_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]
    F_Y_0_curr <- F_Y_0_opt_preds[i,]
    G_W_0_curr <- G_W_0_opt_preds[i,]
    G_W_1_curr <- G_W_1_opt_preds[i,]
    if (surv_form == "PI"){
      S_T_ests <-compute_prodint(cdf_uncens = F_Y_1_curr,
                                 cdf_cens = F_Y_0_curr,
                                 entry_uncens = G_W_1_curr,
                                 entry_cens = G_W_0_curr,
                                 p_uncens = pi_curr,
                                 newtimes = newtimes,
                                 time_grid = object$time_grid_approx)
      S_C_ests <-compute_prodint(cdf_uncens = F_Y_0_curr,
                                 cdf_cens = F_Y_1_curr,
                                 entry_uncens = G_W_0_curr,
                                 entry_cens = G_W_1_curr,
                                 p_uncens = 1 - pi_curr,
                                 newtimes = newtimes,
                                 time_grid = object$time_grid_approx)
    } else if (surv_form == "exp"){
      S_T_ests <-compute_exponential(cdf_uncens = F_Y_1_curr,
                                     cdf_cens = F_Y_0_curr,
                                     entry_uncens = G_W_1_curr,
                                     entry_cens = G_W_0_curr,
                                     p_uncens = pi_curr,
                                     newtimes = newtimes,
                                     time_grid = object$time_grid_approx)
      S_C_ests <-compute_exponential(cdf_uncens = F_Y_0_curr,
                                     cdf_cens = F_Y_1_curr,
                                     entry_uncens = G_W_0_curr,
                                     entry_cens = G_W_1_curr,
                                     p_uncens = 1 - pi_curr,
                                     newtimes = newtimes,
                                     time_grid = object$time_grid_approx)
    }

    return(list(S_T_ests = S_T_ests, S_C_ests = S_C_ests))
  }

  preds <- t(matrix(unlist(apply(X = as.matrix(seq(1, nrow(newX))),
                                 MARGIN = 1,
                                 FUN = estimate_S_T)), nrow = 2*length(newtimes)))

  S_T_preds <- preds[,1:length(newtimes)]
  S_C_preds <- preds[,(length(newtimes) + 1):(2*length(newtimes))]

  if (object$direction == "retrospective"){
    S_T_preds <- 1 - S_T_preds
    S_C_preds <- NULL
  }

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = S_C_preds,
              surv_form = surv_form)
  return(res)

}



####################.
##### stackL.R #####
####################.

#' Estimate a conditional survival function via local survival stacking
#'
#' @noRd
stackL <- function(time,
                   event = rep(1, length(time)),
                   entry = NULL,
                   X,
                   newX,
                   newtimes,
                   direction = "prospective",
                   bin_size = NULL,
                   time_basis = "continuous",
                   SL_control = list(SL.library = c("SL.mean"),
                                     V = 10,
                                     method = "method.NNLS",
                                     stratifyCV = FALSE),
                   tau = NULL){

  if (is.null(newX)){
    newX <- X
  }

  if (is.null(newtimes)){
    newtimes <- time
  }

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    newtimes <- tau - newtimes
    entry <- tau - entry
    event <- rep(1, length(time))
  }

  X <- as.matrix(X)
  time <- as.matrix(time)
  event <- as.matrix(event)
  dat <- data.frame(X, time, event)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- sort(unique(stats::quantile(dat$time[dat$event == 1], probs = seq(0, 1, by = bin_size))))
    time_grid[1] <- 0
  } else{
    time_grid <- sort(unique(dat$time[dat$event == 1]))
    time_grid <- c(0, time_grid)
  }

  # this truncated time grid does not include the first time, since our discretization
  # convention pushes events with t= < time < t + 1 to time t
  trunc_time_grid <- time_grid#[-length(trunc_time_grid)]

  if (!is.null(SL_control$obsWeights)){
    stackX <- as.matrix(data.frame(X, obsWeights = SL_control$obsWeights))
  } else{
    stackX <- X
  }

  # create stacked dataset
  stacked <- stack_haz(time = time,
                       event = event,
                       X = stackX,
                       time_grid = time_grid,
                       entry = entry,
                       time_basis = "continuous")
  #print(stacked$event_indicators)

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }

  long_obsWeights <- stacked$obsWeights
  stacked$obsWeights <- NULL
  .Y <- stacked[,ncol(stacked)]
  .X <- data.frame(stacked[,-ncol(stacked)])
  # fit Super Learner
  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  if (is.null(SL_control$stratifyCV)){
    SL_control$stratifyCV <- FALSE
  }

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL_control$SL.library,
                                    family = stats::binomial(),
                                    method = SL_control$method,
                                    verbose = FALSE,
                                    obsWeights = long_obsWeights,
                                    cvControl = list(V = SL_control$V,
                                                     stratifyCV = SL_control$stratifyCV))

  # create function to get discrete hazard predictions
  if (time_basis == "continuous"){
    get_hazard_preds <- function(index){
      new_stacked <- data.frame(t = trunc_time_grid[index], newX)
      preds <- stats::predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  } else if (time_basis == "dummy"){
    get_hazard_preds <- function(index){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  }

  # don't estimate hazard at t =0
  #hazard_preds <- apply(X = matrix(time_grid), FUN = get_hazard_preds, MARGIN = 1)
  hazard_preds <- apply(X = matrix(1:length(trunc_time_grid)),
                        FUN = get_hazard_preds,
                        MARGIN = 1)

  get_surv_preds <- function(t){
    if (sum(trunc_time_grid <= t) != 0){ # if you don't fall before the first time in the grid
      final_index <- max(which(trunc_time_grid <= t))
      haz <- as.matrix(hazard_preds[,1:final_index])
      anti_haz <- 1 - haz
      surv <- apply(anti_haz, MARGIN = 1, prod)
    } else{
      surv <- rep(1, nrow(hazard_preds))
    }
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)

  if (direction == "retrospective"){
    surv_preds <- 1 - surv_preds
  }

  res <- list(S_T_preds = surv_preds,
              direction = direction,
              time_basis = time_basis,
              time_grid = time_grid,
              tau = tau,
              fit = fit)
  class(res) <- "stackL"
  return(res)
}

#' @noRd
#' @export
predict.stackL <- function(object,
                           newX,
                           newtimes){

  trunc_time_grid <- object$time_grid

  # create function to get discrete hazard predictions
  if (object$time_basis == "continuous"){
    get_hazard_preds <- function(index){
      new_stacked <- data.frame(t = trunc_time_grid[index], newX)
      preds <- stats::predict(object$fit, newdata=new_stacked)$pred
      return(preds)
    }
  } else if (object$time_basis == "dummy"){
    get_hazard_preds <- function(index){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(object$fit, newdata=new_stacked)$pred
      return(preds)
    }
  }

  # don't estimate hazard at t =0
  #hazard_preds <- apply(X = matrix(time_grid), FUN = get_hazard_preds, MARGIN = 1)
  hazard_preds <- apply(X = matrix(1:length(trunc_time_grid)),
                        FUN = get_hazard_preds,
                        MARGIN = 1)

  get_surv_preds <- function(t){
    if (sum(trunc_time_grid <= t) != 0){ # if you don't fall before the first time in the grid
      final_index <- max(which(trunc_time_grid <= t))
      haz <- as.matrix(hazard_preds[,1:final_index])
      anti_haz <- 1 - haz
      surv <- apply(anti_haz, MARGIN = 1, prod)
    } else{
      surv <- rep(1, nrow(hazard_preds))
    }
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)


  if (object$direction == "retrospective"){
    surv_preds <- 1 - surv_preds
  }

  return(list(S_T_preds = surv_preds))
}



#########################.
##### stack_utils.R #####
#########################.

#' Stack a dataset, using the hazard
#'
#' @return A stacked dataset
#' @noRd
stack_haz <- function(time, event, X, time_grid, entry, time_basis){
  trunc_time_grid <- time_grid#[-length(time_grid)]
  time_grid <- c(time_grid, max(time) + 1)
  dat <- data.frame(X, event = event, time = time)

  if (time_basis == "continuous"){
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  } else if (time_basis == "dummy"){
    ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
  }
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){ # can change this to not do anything in last time bin
    if (is.null(entry)){ # no entry (truncation variable) given
      #risk_set <- dat[dat$time > time_grid[i],]
      risk_set <- dat[dat$time >= time_grid[i],]
    } else{ # entry given
      risk_set <- dat[dat$time >= time_grid[i] & entry < time_grid[i+1],]
      #risk_set <- dat[dat$time >= time_grid[i] & entry <= time_grid[i],]
    }
    risk_set_covariates <- risk_set[,1:ncol(X)]
    #event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
    event_indicators <- matrix(ifelse(risk_set$time < time_grid[i + 1 ] & risk_set$event == 1, 1, 0))

    if (time_basis == "continuous"){
      #t <- rep(time_grid[i + 1], nrow(event_indicators))
      t <- rep(time_grid[i], nrow(event_indicators))
    } else if (time_basis == "dummy"){
      t <- matrix(0, ncol = length(trunc_time_grid),
                  nrow = nrow(risk_set))
      t[,i] <- 1
    }
    newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
    stacked <- rbind(stacked, newdata)
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
    if (time_basis == "dummy"){
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
    }
  }
  stacked <- stacked[-1,]
  stacked <- data.frame(stacked)
  return(stacked)
}

#' Stack a dataset, using the entry time conditional on event time
#'
#' @return A stacked dataset
#' @noRd
stack_entry <- function(time, entry, X, time_grid, time_basis){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, time = time, entry = entry)
  if (time_basis == "continuous"){
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(trunc_time_grid))){
      #for (i in 1:(length(time_grid))){# can change this to not do anything in last time bin
      # if (i > 1){
      #   risk_set <- dat[dat$time > time_grid[i-1],]
      # } else{
      #   risk_set <- dat
      # }
      risk_set <- dat[dat$time > time_grid[i],]# maybe this should be >= i+1? Need to think more carefully about this. obv in the limit, doesn't matter
      risk_set_covariates <- risk_set[,1:ncol(X)]
      #event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i], 1, 0))
      event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
      #t <- rep(time_grid[i], nrow(risk_set_covariates))
      t <- rep(time_grid[i + 1], nrow(risk_set_covariates))
      newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  } else if (time_basis ==  "dummy"){
    ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(trunc_time_grid))){
      risk_set <- dat[dat$time > time_grid[i],]# maybe this should be >= i+1? Need to think more carefully about this. obv in the limit, doesn't matter
      risk_set_covariates <- risk_set[,1:ncol(X)]
      event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
      dummies[,i] <- 1
      newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(stacked)[1:(length(time_grid))] <- risk_set_names
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  }
  stacked <- data.frame(stacked)
  return(stacked)
}

#' Stack a dataset for CDF estimation
#'
#' @return A stacked dataset
#' @noRd
stack_cdf <- function(time, X, time_grid, time_basis){

  if (time_basis == "continuous"){
    # we will treat time as continuous
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(time_grid))){ # can change this to not do anything in last time bin
      event_indicators <- matrix(ifelse(time <= time_grid[i], 1, 0))
      t <- time_grid[i]
      newdata <- as.matrix(cbind(t, X, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  } else if (time_basis == "dummy"){
    dat <- data.frame(X, time)
    ncol_stacked <- ncol(X) + length(time_grid) + 1 # covariates, risk set dummies, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(time_grid))){
      risk_set <- dat
      risk_set_covariates <- risk_set[,1:ncol(X)]
      event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i], 1, 0))
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(risk_set))
      dummies[,i] <- 1
      newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }

    stacked <- stacked[-1,]
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(stacked)[1:(length(time_grid))] <- risk_set_names
  }

  stacked <- data.frame(stacked)
  return(stacked)
}
