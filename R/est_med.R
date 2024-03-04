#' Estimate mediation effects
#'
#' @description Estimate mediation effects, including the natural direct effect
#'     (NDE), the natural indirect effect (NIE), and the proportion mediated
#'     (PM). See references for definitions of these objects.
#' @param dat A data object returned by load_data
#' @param type One of c("NP", "Cox"). This specifies whether to estimate the
#'     effects using a marginalized Cox proportional hazards model or using a
#'     nonparametric estimator.
#' @param t_0 Time point of interest
#' @param nde Boolean. If TRUE, the natural direct effect is computed and
#'     returned.
#' @param nie Boolean. If TRUE, the natural indirect effect is computed and
#'     returned.
#' @param pm Boolean. If TRUE, the proportion mediated is computed and returned.
#' @param scale One of c("RR", "VE"). This determines whether NDE and NIE
#'     estimates and CIs are computed on the risk ratio (RR) scale or the
#'     vaccine efficacy (VE) scale. The latter equals one minus the former.
#' @param params_np A list of options returned by \code{\link{params_med_np}}
#'     that are relevant if type="NP".
#' @return A dataframe containing the following columns: \itemize{
#'     \item{\code{effect}: one of c("NDE", "NIE", "PM")}
#'     \item{\code{est}: point estimate}
#'     \item{\code{se}: standard error of point estimate}
#'     \item{\code{ci_lower}: a confidence interval lower limit}
#'     \item{\code{ci_upper}: a confidence interval upper limit}
#' }
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_np <- est_med(dat=dat, type="NP", t_0=578)
#' }
#' @references Fay MP and Follmann DA (2023). Mediation Analyses for the Effect
#'     of Antibodies in Vaccination <doi:10.48550/arXiv.2208.06465>
#' @export
est_med <- function(
    dat, type="NP", t_0, nde=TRUE, nie=TRUE, pm=TRUE, scale="RR",
    # ci_type="transformed", return_extras=FALSE,
    # params_cox=params_med_cox(),
    params_np=params_med_np()
) {

  # !!!!! Need to refactor with est_ce functions

  if (!(scale %in% c("RR", "VE"))) {
    stop("`scale` must equal one of c('RR', 'VE').")
  }

  if (!methods::is(dat,"vaccine_dat")) {
    stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by ",
                "load_data()."))
  }

  if (!(attr(dat, "groups") %in% c("vaccine", "both"))) {
    stop("Vaccine group data not detected.")
  }

  if ((nde||pm) && !(attr(dat, "groups") %in% c("placebo", "both"))) {
    stop("Placebo group data not detected.")
  }

  if (!(nde||nie||pm)) {
    stop("At least one of `nde`, `nie`, or `pm` must be TRUE.")
  }

  if (type=="Cox") { stop("Coming soon!") }

  if (type=="NP") {

    # Set params
    .default_params <- params_med_np()
    for (i in c(1:length(.default_params))) {
      p_name <- names(.default_params)[i]
      if (is.null(params_np[[p_name]])) { params_np[[p_name]] <- .default_params[[i]] }
    }
    p <- params_np

    # Create filtered data objects
    dat_v <- dat[dat$a==1,]
    dat_p <- dat[dat$a==0,]

    # Alias variables
    dim_x <- attr(dat, "dim_x")
    n_vacc <- attr(dat, "n_vacc")
    n_plac <- attr(dat, "n_plac")
    n <- attr(dat, "n")
    p_vacc <- n_vacc/n
    p_plac <- n_plac/n

    # Rescale S to lie in [0,1] and create grid for rounding
    s_min <- min(dat_v$s, na.rm=T)
    s_max <- max(dat_v$s, na.rm=T)
    s_shift <- -1 * s_min
    s_scale <- 1/(s_max-s_min)
    dat_v$s <- (dat_v$s+s_shift)*s_scale
    dat$s <- ifelse(dat$a==1, (dat$s+s_shift)*s_scale, 0)
    grid <- create_grid(dat_v, p$grid_size, t_0) # !!!!! feed in dat instead (1)

    # Create rounded filtered data objects
    dat_v_rd <- round_dat(dat_v, grid, p$grid_size) # !!!!! see (1) above
    dat_p_rd <- round_dat(dat_p, grid, p$grid_size) # !!!!! see (1) above
    datx_v_rd <- dat_v_rd[, c(1:dim_x), drop=F]
    datx_p_rd <- dat_p_rd[, c(1:dim_x), drop=F]
    class(datx_v_rd) <- "data.frame"
    class(datx_p_rd) <- "data.frame"
    dat_v2_rd <- dat_v_rd[dat_v_rd$z==1,]
    dat_rd <- round_dat(dat, grid, p$grid_size) # !!!!! see (1) above

    # Precomputation values for conditional survival/censoring estimators
    x_distinct <- dplyr::distinct(datx_v_rd)
    x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
    vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
    vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
    vals <- list(
      t = vals_pre$t,
      x = subset(vals_pre, select=names(datx_v_rd)),
      s = vals_pre$s
    )

    # Fit conditional survival estimator (vaccine group)
    srvSL <- construct_Q_n(p$surv_type, dat_v2_rd, vals)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens

    # Compute various nuisance functions
    omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
    f_sIx_n <- construct_f_sIx_n(dat_v2_rd, type=p$density_type, k=p$density_bins, z1=F)
    f_s_n <- construct_f_s_n(dat_v_rd, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)

    # Compute edge-corrected estimator and standard error
    p_n <- (1/n_vacc) * sum(dat_v2_rd$weights * In(dat_v2_rd$s!=0))
    g_sn <- construct_g_sn(dat_v2_rd, f_n_srv, g_n, p_n)
    r_Mn_edge_est <- r_Mn_edge(dat_v_rd, g_sn, g_n, p_n, Q_n, omega_n, t_0)
    r_Mn_edge_est <- min(max(r_Mn_edge_est, 0), 1)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
                                                     r_Mn_edge_est, p_n, t_0)

    # Create results object
    res <- data.frame(
      "effect" = character(),
      "est" = double(),
      "se" = double(),
      "ci_lower" = double(),
      "ci_upper" = double()
    )

    # Create edge influence function (relative to the whole dataset)
    infl_fn_edge_2 <- function(a,z,weights,s,x,y,delta) {
      if (a==0) {
        return(0)
      } else {
        return((1/p_vacc)*infl_fn_r_Mn_edge(z,weights,s,x,y,delta))
      }
    }

    # Calculate NDE
    if (nde) {

      # Precomputation values for conditional survival/censoring estimators
      x_distinct_p <- dplyr::distinct(datx_p_rd)
      x_distinct_p <- cbind("x_index"=c(1:nrow(x_distinct_p)), x_distinct_p)
      vals_p_pre <- expand.grid(t=grid$y, x_index=x_distinct_p$x_index)
      vals_p_pre <- dplyr::inner_join(vals_p_pre, x_distinct_p, by="x_index")
      vals_p <- list(
        t = vals_p_pre$t,
        x = subset(vals_p_pre, select=names(datx_p_rd))
      )

      # Fit conditional survival estimator (placebo group)
      srvSL_p <- construct_Q_noS_n(type=p$surv_type, dat_p_rd, vals_p)
      Q_noS_n <- srvSL_p$srv
      Qc_noS_n <- srvSL_p$cens

      # Construct other nuisance estimators
      omega_noS_n <- construct_omega_noS_n(Q_noS_n, Qc_noS_n, t_0, grid)
      risk_p_est <- risk_overall_np_p(dat_p_rd, Q_noS_n, omega_noS_n, t_0)
      infl_fn_risk_p <- construct_infl_fn_risk_p(dat_p_rd, Q_noS_n, omega_noS_n,
                                                 t_0, p_plac)

      # !!!!! NEW CODE
      if (T) {

        # Precomputation values for conditional survival/censoring estimators
        x_distinct_v <- dplyr::distinct(datx_v_rd)
        x_distinct_v <- cbind("x_index"=c(1:nrow(x_distinct_v)), x_distinct_v)
        vals_v_pre <- expand.grid(t=grid$y, x_index=x_distinct_v$x_index)
        vals_v_pre <- dplyr::inner_join(vals_v_pre, x_distinct_v, by="x_index")
        vals_v <- list(
          t = vals_v_pre$t,
          x = subset(vals_v_pre, select=names(datx_v_rd))
        )

        # Fit conditional survival estimator (placebo group)
        srvSL_v <- construct_Q_noS_n(type=p$surv_type, dat_v_rd, vals_v)
        Q_noS_n_v <- srvSL_v$srv
        Qc_noS_n_v <- srvSL_v$cens

        # Construct other nuisance estimators
        omega_noS_n_v <- construct_omega_noS_n(Q_noS_n_v, Qc_noS_n_v, t_0, grid)
        risk_v_est <- risk_overall_np_v_v2(dat_v_rd, Q_noS_n_v, omega_noS_n_v,
                                           t_0)
        # print("!!!!! DEBUGGING !!!!!")
        # print(paste("risk_v_est:", risk_v_est))
        # print(paste("risk_p_est:", risk_p_est))
        infl_fn_risk_v <- construct_infl_fn_risk_v_v2(dat_v_rd, Q_noS_n_v, # !!!!! Should still work
                                                   omega_noS_n_v, t_0, p_vacc) # !!!!! Should still work

      }

      # q_tilde_n <- memoise2(function(x,y,delta) {
      #   # !!!!! Replace with sapply
      #   num <- sum(apply(dat_v2_rd, 1, function(r) {
      #     r[["weights"]] * (omega_n(x,r[["s"]],y,delta)-Q_n(t_0,x,r[["s"]])) *
      #       f_n_srv(y, delta, x, r[["s"]]) * g_n(r[["s"]],x)
      #   }))
      #   den <- sum(apply(dat_v2_rd, 1, function(r) {
      #     r[["weights"]] * f_n_srv(y, delta, x, r[["s"]]) * g_n(r[["s"]],x)
      #   }))
      #   if (den==0) {
      #     return(0)
      #   } else {
      #     return(num/den)
      #   }
      # })

      # risk_v_est <- risk_overall_np_v(dat_v_rd, g_n, Q_n, omega_n, f_n_srv,
      #                                 q_tilde_n, t_0)
      # infl_fn_risk_v <- construct_infl_fn_risk_v(dat_v_rd, Q_n, g_n, omega_n,
      #                                            q_tilde_n, t_0, p_vacc)

      nde_est <- r_Mn_edge_est/risk_p_est

      sigma2_nde_est <- mean(apply(dat_rd, 1, function(r) {
        x <- as.numeric(r[1:dim_x])
        if_p <- infl_fn_risk_p(r[["a"]], r[["delta"]], r[["y"]], x)
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])
        return(((1/risk_p_est)*if_edge-(r_Mn_edge_est/risk_p_est^2)*if_p)^2)
      }), na.rm=T) # !!!!! Test getting rid of na.rm=T

      nde_se <- sqrt(sigma2_nde_est/n)
      nde_lo <- exp(log(nde_est)-1.96*(1/nde_est)*nde_se)
      nde_up <- exp(log(nde_est)+1.96*(1/nde_est)*nde_se)
      res[nrow(res)+1,] <- list("NDE", nde_est, nde_se, nde_lo, nde_up)

    }

    # Calculate NIE
    if (nie) {

      nie_est <- risk_v_est/r_Mn_edge_est

      sigma2_nie_est <- mean(apply(dat_rd, 1, function(r) {
        x <- as.numeric(r[1:dim_x])
        # if_v <- infl_fn_risk_v(a=r[["a"]], z=r[["z"]], weight=r[["weights"]],
        #                        s=r[["s"]], x=x, y=r[["y"]],
        #                        delta=r[["delta"]])
        if_v <- infl_fn_risk_v(r[["a"]], r[["delta"]], r[["y"]], x)
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])
        return(((1/r_Mn_edge_est)*if_v-(risk_v_est/r_Mn_edge_est^2)*if_edge)^2)
      }), na.rm=T) # !!!!! Test getting rid of na.rm=T

      nie_se <- sqrt(sigma2_nie_est/n)
      nie_lo <- exp(log(nie_est)-1.96*(1/nie_est)*nie_se)
      nie_up <- exp(log(nie_est)+1.96*(1/nie_est)*nie_se)
      res[nrow(res)+1,] <- list("NIE", nie_est, nie_se, nie_lo, nie_up)

    }

    # Calculate PM
    if (pm) {

      pm_est <- 1 - (log(r_Mn_edge_est/risk_p_est) / log(risk_v_est/risk_p_est))

      sigma2_pm_est <- mean(apply(dat_rd, 1, function(r) {
        rr <- risk_v_est/risk_p_est
        c_1 <- log(r_Mn_edge_est/risk_p_est) / risk_v_est
        c_2 <- log(risk_v_est/r_Mn_edge_est) / risk_p_est
        c_3 <- (-1*log(rr)) / r_Mn_edge_est
        x <- as.numeric(r[1:dim_x])
        # if_v <- infl_fn_risk_v(a=r[["a"]], z=r[["z"]], weight=r[["weights"]],
        #                        s=r[["s"]], x=x, y=r[["y"]],
        #                        delta=r[["delta"]])
        if_v <- infl_fn_risk_v(r[["a"]], r[["delta"]], r[["y"]], x)
        if_p <- infl_fn_risk_p(r[["a"]], r[["delta"]], r[["y"]], x)
        if_edge <- infl_fn_edge_2(a=r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])
        return((1/(log(rr))^2*(c_1*if_v+c_2*if_p+c_3*if_edge))^2)
      }), na.rm=T) # !!!!! Test getting rid of na.rm=T

      pm_se <- sqrt(sigma2_pm_est/n)
      pm_lo <- pm_est-1.96*pm_se
      pm_up <- pm_est+1.96*pm_se
      res[nrow(res)+1,] <- list("PM", pm_est, pm_se, pm_lo, pm_up)

    }

    if (scale=="VE") {

      # NDE
      res[res$effect=="NDE","est"] <- 1 - res[res$effect=="NDE","est"]
      res_lo_old <- res[res$effect=="NDE","ci_lower"]
      res_up_old <- res[res$effect=="NDE","ci_upper"]
      res[res$effect=="NDE","ci_lower"] <- 1 - res_up_old
      res[res$effect=="NDE","ci_upper"] <- 1 - res_lo_old

      # NIE
      res[res$effect=="NIE","est"] <- 1 - res[res$effect=="NIE","est"]
      res_lo_old <- res[res$effect=="NIE","ci_lower"]
      res_up_old <- res[res$effect=="NIE","ci_upper"]
      res[res$effect=="NIE","ci_lower"] <- 1 - res_up_old
      res[res$effect=="NIE","ci_upper"] <- 1 - res_lo_old

    }

  }

  return(res)

}
