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
#' @param scale One of c("RR", "VE"). This determines whether estimates and CIs
#'     are computed on the risk ratio (RR) scale or the vaccine efficacy (VE)
#'     scale. The latter equals one minus the former.
#' @return A dataframe containing the following columns: \itemize{
#'     \item{\code{effect}: one of c("NDE", "NIE", "PM")}
#'     \item{\code{est}: point estimate}
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
    dat, type="NP", t_0, nde=TRUE, nie=TRUE, pm=TRUE, scale="RR"
    # ci_type="transformed", return_extras=FALSE,
    # params_cox=params_med_cox(), params_np=params_med_cox()
) {

  # !!!!! Need to refactor with est_ce functions
  # !!!!! Integrate params

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
    .default_params <- list(
      # surv_type = "survML-G",
      surv_type = "survSL",
      density_type = "binning",
      density_bins = 15,
      # deriv_type = "m-spline",
      # gamma_type = "Super Learner",
      q_n_type = "zero" # standard
    )

    # for (i in c(1:length(.default_params))) {
    #   p_name <- names(.default_params)[i]
    #   if (is.null(params[[p_name]])) { params[[p_name]] <- .default_params[[i]] }
    # }
    # p <- params
    p <- .default_params # !!!!! TEMP

    grid_size <- list(y=101, s=101, x=5) # !!!!! TEMP

    # Rescale S to lie in [0,1] and create rounded data object
    s_min <- min(dat_orig$s, na.rm=T)
    s_max <- max(dat_orig$s, na.rm=T)
    s_shift <- -1 * s_min
    s_scale <- 1/(s_max-s_min)
    dat_orig$s <- (dat_orig$s+s_shift)*s_scale
    grid <- create_grid(dat_orig, grid_size, t_0)
    dat_orig_rounded <- round_dat(dat_orig, grid, grid_size)

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
    dat <- ss(dat_orig_rounded, which(dat_orig_rounded$z==1)) # !!!!! uses the `dat` name again; change???

    # Fit conditional survival estimator
    srvSL <- construct_Q_n(p$surv_type, dat, vals)
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
      "ci_lower" = double(),
      "ci_upper" = double()
    )

    # Calculate NDE
    if (nde) {
      ov <- est_overall(dat=dat_copy, t_0=t_0, method="KM", ve=F)
      risk_p <- ov[ov$group=="placebo","est"]
      se_p <- ov[ov$group=="placebo","se"] # This equals sd_p/n_orig
      nde_est <- r_Mn_edge_est/risk_p
      nde_se <- sqrt(
        r_Mn_edge_est^2*se_p^2/risk_p^4 + se_edge_est^2/(n_orig*risk_p^2)
      )
      nde_lo <- exp(log(nde_est)-1.96*(1/nde_est)*nde_se)
      nde_up <- exp(log(nde_est)+1.96*(1/nde_est)*nde_se)
      res[nrow(res)+1,] <- c("NDE", nde_est, nde_lo, nde_up)
    }

    # Calculate PM
    if (pm) {
      pm_est <- 999 # !!!!! TO DO
      pm_lo <- 999 # !!!!! TO DO
      pm_up <- 999 # !!!!! TO DO
      res[nrow(res)+1,] <- c("PM", pm_est, pm_lo, pm_up)
    }

    # Calculate NIE
    if (nie) {
      nie_est <- 999 # !!!!! TO DO
      nie_lo <- 999 # !!!!! TO DO
      nie_up <- 999 # !!!!! TO DO
      res[nrow(res)+1,] <- c("NIE", nie_est, nie_lo, nie_up)
    }

    if (scale=="VE") {
      res$est <- 1 - res$est
      res_lo_old <- res$ci_lo
      res_up_old <- res$ci_up
      res$ci_upper <- 1 - res_lo_old
      res$ci_lower <- 1 - res_up_old
    }

  }

  return(res)

}
