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

    # Alias variables
    dat_copy <- dat # Consider keeping this and changing old `dat` variable
    dat_orig <- dat$v

    # Set params
    .default_params <- params_med_np()
    for (i in c(1:length(.default_params))) {
      p_name <- names(.default_params)[i]
      if (is.null(params_np[[p_name]])) { params_np[[p_name]] <- .default_params[[i]] }
    }
    p <- params_np

    # Rescale S to lie in [0,1] and create rounded data object
    s_min <- min(dat_orig$s, na.rm=T)
    s_max <- max(dat_orig$s, na.rm=T)
    s_shift <- -1 * s_min
    s_scale <- 1/(s_max-s_min)
    dat_orig$s <- (dat_orig$s+s_shift)*s_scale
    grid <- create_grid(dat_orig, p$grid_size, t_0)
    dat_orig_rounded <- round_dat(dat_orig, grid, p$grid_size)

    # Prepare precomputation values for conditional survival estimator
    x_distinct <- dplyr::distinct(dat_orig_rounded$x)
    x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
    vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index, s=grid$s)
    vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
    vals <- list(
      t = vals_pre$t,
      x = subset(vals_pre, select=names(dat_orig_rounded$x)),
      s = vals_pre$s
    )

    # Create phase-two data object (unrounded)
    dat <- ss(dat_orig_rounded, which(dat_orig_rounded$z==1)) # !!!!! uses the `dat` name again; change

    # Fit conditional survival estimator
    print("check 1")
    srvSL <- construct_Q_n(p$surv_type, dat, vals)
    print("check 2")
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens

    # Use rounded data objects moving forward
    dat_orig <- dat_orig_rounded

    # Compute various nuisance functions
    omega_n <- construct_omega_n(Q_n, Qc_n, t_0, grid)
    f_sIx_n <- construct_f_sIx_n(dat, type=p$density_type, k=p$density_bins, z1=F)
    f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    f_n_srv <- construct_f_n_srv(Q_n, Qc_n, grid)

    # Compute edge-corrected estimator and standard error
    n_orig <- attr(dat_orig, "n_orig")
    p_n <- (1/n_orig) * sum(dat$weights * In(dat$s!=0))
    print(paste0("Percent of mass at edge: ", round(100*(1-p_n),1), "%")) # !!!!! TEMP
    g_sn <- construct_g_sn(dat, f_n_srv, g_n, p_n)
    r_Mn_edge_est <- r_Mn_edge(dat_orig, g_sn, g_n, p_n, Q_n, omega_n, t_0)
    r_Mn_edge_est <- min(max(r_Mn_edge_est, 0), 1)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n, g_n,
                                                     r_Mn_edge_est, p_n, t_0)
    dat_orig_df <- as_df(dat_orig)
    dim_x <- attr(dat_orig, "dim_x")
    sigma2_edge_est <- mean(apply(dat_orig_df, 1, function(r) {
      (infl_fn_r_Mn_edge(r[["z"]], r[["weights"]], r[["s"]],
                         as.numeric(r[1:dim_x]), r[["y"]], r[["delta"]]))^2
    }))
    se_edge_est <- sqrt(sigma2_edge_est/n_orig)

    # Create results object
    res <- data.frame(
      "effect" = character(),
      "est" = double(),
      "se" = double(),
      "ci_lower" = double(),
      "ci_upper" = double()
    )

    # Create necessary objects and functions
    {

      ov <- est_overall(dat=dat_copy, t_0=t_0, method="KM", ve=F)
      risk_p <- ov[ov$group=="placebo","est"]
      risk_v <- ov[ov$group=="vaccine","est"]
      se_p <- ov[ov$group=="placebo","se"] # This equals sd_p/n_orig
      se_v <- ov[ov$group=="vaccine","se"] # This equals sd_v/n_orig

      # !!!!! Temporary until data is restructured
      # dat_v_copy <- dat_copy$v
      dat_orig_rounded_p <- round_dat(dat_copy$p, grid, p$grid_size)
      dat_v_copy <- dat_orig # !!!!! Changed; now using rounded data values
      dat_p_copy <- dat_orig_rounded_p
      num_v <- length(dat_v_copy$z)
      num_p <- length(dat_p_copy$z)
      PA1 <- num_v/(num_v+num_p)
      dat_v_copy$a <- rep(1, num_v)
      dat_p_copy$a <- rep(0, num_p)
      dat_combined <- cbind(
        rbind(dat_v_copy$x, dat_p_copy$x),
        a = c(dat_v_copy$a,dat_p_copy$a),
        y = c(dat_v_copy$y,dat_p_copy$y),
        delta = c(dat_v_copy$delta,dat_p_copy$delta),
        s = c(dat_v_copy$s,dat_p_copy$s),
        weights = c(dat_v_copy$weights,dat_p_copy$weights),
        strata = c(dat_v_copy$strata,dat_p_copy$strata),
        z = c(dat_v_copy$z,dat_p_copy$z)
      )

      infl_fn_km_v <- construct_km_infl_fn(dat_combined, which="vaccine") # OLD
      infl_fn_km_p <- construct_km_infl_fn(dat_combined, which="placebo") # OLD

      infl_fn_edge_2 <- function(a,z,weights,s,x,y,delta) {
        if (a==0) {
          return(0)
        } else {
          return((1/PA1)*infl_fn_r_Mn_edge(z,weights,s,x,y,delta))
        }
      }

    }

    # Calculate NDE
    if (nde) {

      # !!!! TEMP RESTRUCTURE
      if (T) {

        # !!!!! TEMP; new risk estimators
        x_distinct <- dplyr::distinct(dat_p_copy$x)
        x_distinct <- cbind("x_index"=c(1:nrow(x_distinct)), x_distinct)
        vals_pre <- expand.grid(t=grid$y, x_index=x_distinct$x_index)
        vals_pre <- dplyr::inner_join(vals_pre, x_distinct, by="x_index")
        vals <- list(
          t = vals_pre$t,
          x = subset(vals_pre, select=names(dat_p_copy$x))
        )
        print("check 3")
        srvSL_p <- construct_Q_noS_n(type="survSL", dat_p_copy, vals)
        print("check 4")
        Q_noS_n <- srvSL_p$srv
        Qc_noS_n <- srvSL_p$cens
        omega_noS_n <- construct_omega_noS_n(Q_noS_n, Qc_noS_n, t_0, grid)
        print("check 5 (before r_Mn_p3)")
        print(Sys.time())
        r_Mn_p3 <- risk_overall_np_p(dat_p_copy, dim_x, Q_noS_n,
                                     omega_noS_n, t_0)
        prob <- 1 - mean(dat_combined$a)
        infl_fn_risk_p <- construct_infl_fn_risk_p(dat_p_copy, Q_noS_n, # NEW
                                                   omega_noS_n, dim_x, t_0,
                                                   prob) # NEW
        print("check 6 (after r_Mn_p3; before r_Mn_v3)")
        print(Sys.time())

        dat_orig_df_short <- dat_orig_df[dat_orig_df$z==1,]
        q_tilde_n <- memoise2(function(x,y,delta) {
          # !!!!! Replace with sapply
          num <- sum(apply(dat_orig_df_short, 1, function(r) {
            r[["weights"]] * (omega_n(x,r[["s"]],y,delta)-Q_n(t_0,x,r[["s"]])) *
              f_n_srv(y, delta, x, r[["s"]]) * g_n(r[["s"]],x)
          }))
          den <- sum(apply(dat_orig_df_short, 1, function(r) {
            r[["weights"]] * f_n_srv(y, delta, x, r[["s"]]) * g_n(r[["s"]],x)
          }))
          if (den==0) {
            return(0)
          } else {
            return(num/den)
          }
        })

        r_Mn_v3 <- risk_overall_np_v(dat_orig, g_n, Q_n, omega_n, f_n_srv,
                                     q_tilde_n, t_0)
        infl_fn_risk_v <- construct_infl_fn_risk_v(dat_orig, Q_n, g_n, omega_n,
                                                   q_tilde_n, t_0, 1-prob)
        print("check 7 (after r_Mn_v3)")
        print(Sys.time())
        print("r_Mn_v3")
        print(r_Mn_v3)
        print("r_Mn_p3")
        print(r_Mn_p3)

        # !!!!! Special checks
        if (T) {

          sigma2_if_v_v3 <- mean(apply(dat_combined, 1, function(r) {
            x <- as.numeric(r[1:dim_x])
            if_v <- infl_fn_risk_v(a=r[["a"]], z=r[["z"]], weight=r[["weights"]],
                                   s=r[["s"]], x=x, y=r[["y"]],
                                   delta=r[["delta"]])
            return(if_v^2)
          }), na.rm=T)
          sigma2_if_p_v3 <- mean(apply(dat_combined, 1, function(r) {
            x <- as.numeric(r[1:dim_x])
            if_p <- infl_fn_risk_p(r[["a"]], r[["delta"]], r[["y"]], x)
            return(if_p^2)
          }), na.rm=T)
          se_v_v3 <- sqrt(sigma2_if_v_v3/(num_v+num_p))
          se_p_v3 <- sqrt(sigma2_if_p_v3/(num_v+num_p))
          res[nrow(res)+1,] <- list("Risk v3 V", r_Mn_v3, se_v_v3, 999, 999)
          res[nrow(res)+1,] <- list("Risk v3 P", r_Mn_p3, se_p_v3, 999, 999)

        }


      }

      nde_est <- r_Mn_edge_est/risk_p
      nde_se <- sqrt(r_Mn_edge_est^2*se_p^2/risk_p^4 + se_edge_est^2/risk_p^2)
      nde_lo <- exp(log(nde_est)-1.96*(1/nde_est)*nde_se)
      nde_up <- exp(log(nde_est)+1.96*(1/nde_est)*nde_se)
      res[nrow(res)+1,] <- list("NDE", nde_est, nde_se, nde_lo, nde_up)

      # !!!!! TEMP; for comparison (but maybe keep this one instead for consistency)
      print("check 8")
      sigma2_nde_est <- mean(apply(dat_combined, 1, function(r) {

        if_p <- infl_fn_km_p(A_i=r[["a"]],Delta_i=r[["delta"]], Y_i=r[["y"]],
                             t=t_0)
        if_edge <- infl_fn_edge_2(a=r[["a"]], r[["z"]], r[["weights"]],
                                  r[["s"]], as.numeric(r[1:dim_x]), r[["y"]],
                                  r[["delta"]])

        return(((1/risk_p)*if_edge-(r_Mn_edge_est/risk_p^2)*if_p)^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      print("check 9")
      se_nde_est <- sqrt(sigma2_nde_est/(num_v+num_p))
      res[nrow(res)+1,] <- list("NDE2", nde_est, se_nde_est, 999, 999)

      # !!!!! TEMP; return additional values
      ov2 <- est_overall(dat=dat_copy, t_0=t_0, method="Cox", ve=F)
      risk_p2 <- ov2[ov2$group=="placebo","est"]
      risk_se2 <- ov2[ov2$group=="placebo","se"]
      res[nrow(res)+1,] <- list("Risk_p", risk_p, se_p, 999, 999)
      res[nrow(res)+1,] <- list("Risk_p (Cox)", risk_p2, risk_se2, 999, 999)
      res[nrow(res)+1,] <- list("Risk_v", risk_v, se_v, 999, 999)
      res[nrow(res)+1,] <- list("CR(0)", r_Mn_edge_est, 999, 999, 999)
      print("check 10")

      sigma2_nde_est_v3 <- mean(apply(dat_combined, 1, function(r) {

        x <- as.numeric(r[1:dim_x])
        if_p <- infl_fn_risk_p(r[["a"]], r[["delta"]], r[["y"]], x)
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])

        return(((1/r_Mn_p3)*if_edge-(r_Mn_edge_est/r_Mn_p3^2)*if_p)^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      print("check 11")
      nde_est_v3 <- r_Mn_edge_est/r_Mn_p3
      se_nde_est_v3 <- sqrt(sigma2_nde_est_v3/(num_v+num_p))
      nde_lo_v3 <- exp(log(nde_est_v3)-1.96*(1/nde_est_v3)*se_nde_est_v3)
      nde_up_v3 <- exp(log(nde_est_v3)+1.96*(1/nde_est_v3)*se_nde_est_v3)
      res[nrow(res)+1,] <- list("NDE v3", nde_est_v3, se_nde_est_v3,
                                nde_lo_v3, nde_up_v3) # !!!!!
      print("check 12")

    }

    # Calculate NIE
    if (nie) {

      nie_est <- risk_v/r_Mn_edge_est
      sigma2_nie_est <- mean(apply(dat_combined, 1, function(r) {

        x <- as.numeric(r[1:dim_x])
        if_v <- infl_fn_km_v(A_i=r[["a"]], Delta_i=r[["delta"]], Y_i=r[["y"]], t=t_0)
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])

        return(((1/r_Mn_edge_est)*if_v-(risk_v/r_Mn_edge_est^2)*if_edge)^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      nie_se <- sqrt(sigma2_nie_est/(num_v+num_p))
      nie_lo <- exp(log(nie_est)-1.96*(1/nie_est)*nie_se)
      nie_up <- exp(log(nie_est)+1.96*(1/nie_est)*nie_se)
      res[nrow(res)+1,] <- list("NIE", nie_est, nie_se, nie_lo, nie_up)
      print("check 13")

      sigma2_nie_est_v3 <- mean(apply(dat_combined, 1, function(r) {

        x <- as.numeric(r[1:dim_x])
        if_v <- infl_fn_risk_v(a=r[["a"]], z=r[["z"]], weight=r[["weights"]],
                               s=r[["s"]], x=x, y=r[["y"]],
                               delta=r[["delta"]])
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])

        return(((1/r_Mn_edge_est)*if_v-(r_Mn_v3/r_Mn_edge_est^2)*if_edge)^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      print("check 14")
      nie_est_v3 <- r_Mn_v3/r_Mn_edge_est
      nie_se_v3 <- sqrt(sigma2_nie_est_v3/(num_v+num_p))
      nie_lo_v3 <- exp(log(nie_est_v3)-1.96*(1/nie_est_v3)*nie_se_v3)
      nie_up_v3 <- exp(log(nie_est_v3)+1.96*(1/nie_est_v3)*nie_se_v3)
      res[nrow(res)+1,] <- list("NIE v3", nie_est_v3, nie_se_v3, nie_lo_v3, nie_up_v3) # !!!!!

    }

    # Calculate PM
    if (pm) {

      pm_est <- 1 - (log(r_Mn_edge_est/risk_p) / log(risk_v/risk_p))
      print("check 15")
      sigma2_pm_est <- mean(apply(dat_combined, 1, function(r) {

        rr <- risk_v/risk_p
        c_1 <- log(r_Mn_edge_est/risk_p) / risk_v
        c_2 <- log(risk_v/r_Mn_edge_est) / risk_p
        c_3 <- (-1*log(rr)) / r_Mn_edge_est

        x <- as.numeric(r[1:dim_x])
        if_v <- infl_fn_km_v(A_i=r[["a"]], Delta_i=r[["delta"]], Y_i=r[["y"]], t=t_0)
        if_p <- infl_fn_km_p(A_i=r[["a"]], Delta_i=r[["delta"]], Y_i=r[["y"]], t=t_0)
        if_edge <- infl_fn_edge_2(r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])

        return((1/(log(rr))^2*(c_1*if_v+c_2*if_p+c_3*if_edge))^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      print("check 16")
      pm_se <- sqrt(sigma2_pm_est/(num_v+num_p))
      pm_lo <- pm_est-1.96*pm_se
      pm_up <- pm_est+1.96*pm_se
      res[nrow(res)+1,] <- list("PM", pm_est, pm_se, pm_lo, pm_up)

      pm_est_v3 <- 1 - (log(r_Mn_edge_est/r_Mn_p3) / log(r_Mn_v3/r_Mn_p3))
      print("check 17")
      sigma2_pm_est_v3 <- mean(apply(dat_combined, 1, function(r) {

        rr <- r_Mn_v3/r_Mn_p3
        c_1 <- log(r_Mn_edge_est/r_Mn_p3) / r_Mn_v3
        c_2 <- log(r_Mn_v3/r_Mn_edge_est) / r_Mn_p3
        c_3 <- (-1*log(rr)) / r_Mn_edge_est

        x <- as.numeric(r[1:dim_x])
        if_v <- infl_fn_risk_v(a=r[["a"]], z=r[["z"]], weight=r[["weights"]],
                               s=r[["s"]], x=x, y=r[["y"]],
                               delta=r[["delta"]])
        if_p <- infl_fn_risk_p(r[["a"]], r[["delta"]], r[["y"]], x)
        if_edge <- infl_fn_edge_2(a=r[["a"]], r[["z"]], r[["weights"]], r[["s"]], x, r[["y"]], r[["delta"]])

        return((1/(log(rr))^2*(c_1*if_v+c_2*if_p+c_3*if_edge))^2)

      }), na.rm=T) # na.rm included to avoid an issue where KM infl fn is NaN
      pm_se_v3 <- sqrt(sigma2_pm_est_v3/(num_v+num_p))
      pm_lo_v3 <- pm_est_v3-1.96*pm_se_v3
      pm_up_v3 <- pm_est_v3+1.96*pm_se_v3
      res[nrow(res)+1,] <- list("PM v3", pm_est_v3, pm_se_v3, pm_lo_v3, pm_up_v3) # !!!!!
      print("check 18")

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

      # !!!!! TEMP
      print("check 19")
      if (T) {
        res[res$effect=="NDE2","est"] <- 1 - res[res$effect=="NDE2","est"]
        res_lo_old <- res[res$effect=="NDE2","ci_lower"]
        res_up_old <- res[res$effect=="NDE2","ci_upper"]
        res[res$effect=="NDE2","ci_lower"] <- 1 - res_up_old
        res[res$effect=="NDE2","ci_upper"] <- 1 - res_lo_old

        res[res$effect=="NDE v3","est"] <- 1 - res[res$effect=="NDE v3","est"]
        res_lo_old <- res[res$effect=="NDE v3","ci_lower"]
        res_up_old <- res[res$effect=="NDE v3","ci_upper"]
        res[res$effect=="NDE v3","ci_lower"] <- 1 - res_up_old
        res[res$effect=="NDE v3","ci_upper"] <- 1 - res_lo_old

        res[res$effect=="NIE v3","est"] <- 1 - res[res$effect=="NIE v3","est"]
        res_lo_old <- res[res$effect=="NIE v3","ci_lower"]
        res_up_old <- res[res$effect=="NIE v3","ci_upper"]
        res[res$effect=="NIE v3","ci_lower"] <- 1 - res_up_old
        res[res$effect=="NIE v3","ci_upper"] <- 1 - res_lo_old

      }

    }

  }

  return(res)

}
