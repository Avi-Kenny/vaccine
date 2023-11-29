
# Load data
set.seed(1)
data(hvtn505)
hvtn505$x_ch <- sample(letters[1:5], size=nrow(hvtn505), replace=T)
hvtn505$x_fac <- as.factor(hvtn505$x_ch)
set.seed(1)
hvtn505_sample <- hvtn505[sample(c(1:nrow(hvtn505)), size=1000),]
dat <- load_data(
  time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt", marker="IgG_V2",
  covariates=c("age","BMI","bhvrisk"), weights="wt", ph2="casecontrol",
  data=hvtn505
)

# Create a smaller dataset to run faster tests
dat_sample <- load_data(
  time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt", marker="IgG_V2",
  covariates=c("age","BMI","bhvrisk"), weights="wt", ph2="casecontrol",
  data=hvtn505_sample
)

# Make sure factor columns are handled correctly
dat_fac <- load_data(
  time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt", marker="IgG_V2",
  covariates=c("age","BMI","x_fac"), weights="wt", ph2="casecontrol",
  data=hvtn505
)

# Make sure character columns are handled correctly
dat_ch <- load_data(
  time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt", marker="IgG_V2",
  covariates=c("age","BMI","x_ch"), weights="wt", ph2="casecontrol",
  data=hvtn505
)

test_that("load_data", {
  expect_equal(class(dat), c("data.frame", "vaccine_dat"))
  expect_equal(attr(dat, "groups"), "both")
  expect_equal(attr(dat, "covariate_names"), c("age", "BMI", "bhvrisk"))
  expect_equal(attr(dat, "dim_x"), 3)
  expect_equal(attr(dat, "n"), 2302)
  expect_equal(attr(dat, "n_vacc"), 1161)
  expect_equal(attr(dat, "n_plac"), 1141)
  expect_equal(sum(dat[dat$a==1,"y"]), 391608)
  expect_equal(sum(dat[dat$a==0,"y"]), 380935)
  expect_equal(unique(dat$delta), c(0,1))
  expect_equal(sort(unique(dat$strata), na.last=T), c(1:8,NA))
  expect_equal(unique(dat$z), c(0,1))
})

test_that("load_data (character/factor columns)", {
  expect_equal(length(dat_fac), 14)
  expect_equal(length(dat_ch), 14)
  expect_equal(attr(dat_fac, "covariate_names"),
               c("age", "BMI", paste0("x_fac_",c(1:5))))
  expect_equal(attr(dat_ch, "covariate_names"),
               c("age", "BMI", paste0("x_ch_",c(1:5))))
  expect_equal(attr(dat_fac, "dim_x"), 7)
  expect_equal(attr(dat_ch, "dim_x"), 7)
  expect_equal(names(dat_fac)[1:7], paste0("x",c(1:7)))
  expect_equal(names(dat_ch)[1:7], paste0("x",c(1:7)))
  expect_equal(sort(unique(dat_fac$x3)), c(0,1))
  expect_equal(sort(unique(dat_fac$x7)), c(0,1))
  expect_equal(sort(unique(dat_ch$x3)), c(0,1))
  expect_equal(sort(unique(dat_ch$x7)), c(0,1))
})

ss <- summary_stats(dat, quietly=TRUE)

test_that("sumamry_stats", {
  expect_equal(class(ss), "list")
  expect_equal(ss$num_ph1_subj_v, 1161)
  expect_equal(ss$num_ph1_subj_p, 1141)
  expect_equal(ss$num_ph2_subj_v, 150)
  expect_equal(ss$num_ph2_subj_p, 39)
  expect_equal(ss$num_ph1_events_v, 27)
  expect_equal(ss$num_ph1_events_p, 21)
  expect_equal(ss$num_ph2_events_v, 25)
  expect_equal(ss$num_ph2_events_p, 19)
  expect_equal(ss$prop_ph1_events_v, 0.02326)
  expect_equal(ss$prop_ph1_events_p, 0.0184)
  expect_equal(ss$prop_ph2_events_v, 0.16667)
  expect_equal(ss$prop_ph2_events_p, 0.48718)
})

ests_o_KM <- est_overall(dat=dat, t_0=578, method="KM")
ests_risk_p <- ests_o_KM[ests_o_KM$stat=="risk"&ests_o_KM$group=="placebo",]
ests_risk_v <- ests_o_KM[ests_o_KM$stat=="risk"&ests_o_KM$group=="vaccine",]
ests_ve <- ests_o_KM[ests_o_KM$stat=="ve",]

test_that("est_overall (KM)", {
  expect_equal(class(ests_o_KM), c("data.frame", "vaccine_overall"))
  expect_equal(ests_risk_p$est, 0.02879861, tolerance=0.01)
  expect_equal(ests_risk_p$se, 0.006563785, tolerance=0.01)
  expect_equal(ests_risk_p$ci_lower, 0.0162236, tolerance=0.01)
  expect_equal(ests_risk_p$ci_upper, 0.04121288, tolerance=0.01)
  expect_equal(ests_risk_v$est, 0.04067009, tolerance=0.01)
  expect_equal(ests_risk_v$se, 0.008230842, tolerance=0.01)
  expect_equal(ests_risk_v$ci_lower, 0.02506853, tolerance=0.01)
  expect_equal(ests_risk_v$ci_upper, 0.05602199, tolerance=0.01)
  expect_equal(ests_ve$est, -0.4122241, tolerance=0.01)
  expect_equal(ests_ve$se, 0.4304518, tolerance=0.01)
  expect_equal(ests_ve$ci_lower, -1.5666, tolerance=0.01)
  expect_equal(ests_ve$ci_upper, 0.2229498, tolerance=0.01)
})

ests_o_Cox <- est_overall(dat=dat, t_0=578, method="Cox")
ests_risk_p <- ests_o_Cox[ests_o_Cox$stat=="risk"&ests_o_Cox$group=="placebo",]
ests_risk_v <- ests_o_Cox[ests_o_Cox$stat=="risk"&ests_o_Cox$group=="vaccine",]
ests_ve <- ests_o_Cox[ests_o_Cox$stat=="ve",]

test_that("est_overall (Cox)", {
  expect_equal(class(ests_o_Cox), c("data.frame", "vaccine_overall"))
  expect_equal(ests_risk_p$est, 0.02938706, tolerance=0.01)
  expect_equal(ests_risk_p$se, 0.006486545, tolerance=0.01)
  expect_equal(ests_risk_p$ci_lower, 0.0190193, tolerance=0.01)
  expect_equal(ests_risk_p$ci_upper, 0.04514638, tolerance=0.01)
  expect_equal(ests_risk_v$est, 0.04177642, tolerance=0.01)
  expect_equal(ests_risk_v$se, 0.008111679, tolerance=0.01)
  expect_equal(ests_risk_v$ci_lower, 0.02847302, tolerance=0.01)
  expect_equal(ests_risk_v$ci_upper, 0.06090588, tolerance=0.01)
  expect_equal(ests_ve$est, -0.4215925, tolerance=0.01)
  expect_equal(ests_ve$se, 0.4179152, tolerance=0.01)
  expect_equal(ests_ve$ci_lower, -1.529375, tolerance=0.01)
  expect_equal(ests_ve$ci_upper, 0.201018, tolerance=0.01)
})

test_that("est_overall (error handling)", {
  expect_error(
    est_overall(dat=list(), t_0=578, method="KM"),
    paste0("`dat` must be an object of class 'vaccine_dat' returned by load_da",
           "ta().")
  )
})

test_that("est_overall (error handling)", {
  expect_error(
    est_overall(dat=dat, t_0=578, method="hey"),
    "`method` must equal either 'Cox' or 'KM'."
  )
})

set.seed(1)
ests_med_NP <- est_med(dat=dat_sample, type="NP", t_0=578,
                       params_np=params_med_np(surv_type="Cox"))
ests_med_NP_nde <- ests_med_NP[ests_med_NP$effect=="NDE",]
ests_med_NP_nie <- ests_med_NP[ests_med_NP$effect=="NIE",]
ests_med_NP_pm <- ests_med_NP[ests_med_NP$effect=="PM",]

test_that("est_med", {
  expect_equal(class(ests_med_NP), "data.frame")
  expect_equal(ests_med_NP_nde$est, 10.18246, tolerance=0.01)
  expect_equal(ests_med_NP_nde$se, 11.89159, tolerance=0.01)
  expect_equal(ests_med_NP_nde$ci_lower, 1.032185, tolerance=0.01)
  expect_equal(ests_med_NP_nde$ci_upper, 100.4494, tolerance=0.01)
  expect_equal(ests_med_NP_nie$est, 0.154724, tolerance=0.01)
  expect_equal(ests_med_NP_nie$se, 0.1657846, tolerance=0.01)
  expect_equal(ests_med_NP_nie$ci_lower, 0.01894479, tolerance=0.01)
  expect_equal(ests_med_NP_nie$ci_upper, 1.263646, tolerance=0.01)
  expect_equal(ests_med_NP_pm$est, -4.105373, tolerance=0.01)
  expect_equal(ests_med_NP_pm$se, 5.058436, tolerance=0.01)
  expect_equal(ests_med_NP_pm$ci_lower, -14.01991, tolerance=0.01)
  expect_equal(ests_med_NP_pm$ci_upper, 5.809162, tolerance=0.01)
})

set.seed(1)
ests_cox <- est_ce(dat=dat, type="Cox", t_0=578, cve=T)

test_that("est_cox (CR)", {
  expect_equal(class(ests_cox), "vaccine_est")
  expect_equal(ests_cox$cr$s[1], 0, tolerance=0.01)
  expect_equal(ests_cox$cr$s[50], 1.15447037, tolerance=0.01)
  expect_equal(ests_cox$cr$est[1], 0.15486178, tolerance=0.01)
  expect_equal(ests_cox$cr$est[50], 0.08759399, tolerance=0.01)
  expect_equal(ests_cox$cr$se[1], 0.05514197, tolerance=0.01)
  expect_equal(ests_cox$cr$se[50], 0.01915460, tolerance=0.01)
  expect_equal(ests_cox$cr$ci_lower[1], 0.07427856, tolerance=0.01)
  expect_equal(ests_cox$cr$ci_lower[50], 0.05661917, tolerance=0.01)
  expect_equal(ests_cox$cr$ci_upper[1], 0.2950081, tolerance=0.01)
  expect_equal(ests_cox$cr$ci_upper[50], 0.1331231, tolerance=0.01)
})

test_that("est_cox (CVE)", {
  expect_equal(ests_cox$cve$s[1], 0, tolerance=0.01)
  expect_equal(ests_cox$cve$s[50], 1.15447037, tolerance=0.01)
  expect_equal(ests_cox$cve$est[1], -4.3774047, tolerance=0.01)
  expect_equal(ests_cox$cve$est[50], -2.0416048, tolerance=0.01)
  expect_equal(ests_cox$cve$se[1], 2.2734089, tolerance=0.01)
  expect_equal(ests_cox$cve$se[50], 0.9607152, tolerance=0.01)
  expect_equal(ests_cox$cve$ci_lower[1], -11.315226, tolerance=0.01)
  expect_equal(ests_cox$cve$ci_lower[50], -4.648935, tolerance=0.01)
  expect_equal(ests_cox$cve$ci_upper[1], -1.3480269019, tolerance=0.01)
  expect_equal(ests_cox$cve$ci_upper[50], -0.6377176338, tolerance=0.01)
})

set.seed(1)
ests_cox_spl <- est_ce(dat=dat_sample, type="Cox", t_0=578, cve=T,
                        params_cox=params_ce_cox(spline_df=3))

test_that("est_cox spline (CR)", {
  expect_equal(class(ests_cox_spl), "vaccine_est")
  expect_equal(ests_cox_spl$cr$s[1], 0, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$s[50], 1.144744, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$est[1], 0.2861402, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$est[50], 0.05516619, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$se[1], 0.1977738, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$se[50], 0.03191638, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$ci_lower[1], 0.05668296, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$ci_lower[50], 0.01727915, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$ci_upper[1], 0.7278059, tolerance=0.01)
  expect_equal(ests_cox_spl$cr$ci_upper[50], 0.162398, tolerance=0.01)
})

test_that("est_cox spline (CVE)", {
  expect_equal(ests_cox_spl$cve$s[1], 0, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$s[50], 1.144744, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$est[1], -12.76739, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$est[50], -1.654274, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$se[1], 10.91096, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$se[50], 1.848648, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$ci_lower[1], -64.08185, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$ci_lower[50], -9.39444, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$ci_upper[1], -1.912348, tolerance=0.01)
  expect_equal(ests_cox_spl$cve$ci_upper[50], 0.3222172, tolerance=0.01)
})

set.seed(1)
ests_np <- est_ce(dat=dat, type="NP", t_0=578, cve=T,
                  params_np=params_ce_np(surv_type="Cox"))

test_that("est_np (CR)", {
  expect_equal(class(ests_np), "vaccine_est")
  expect_equal(ests_np$cr$s[1], 0, tolerance=0.01)
  expect_equal(ests_np$cr$s[50], 1.15447, tolerance=0.01)
  expect_equal(ests_np$cr$est[1], 0.2692285, tolerance=0.01)
  expect_equal(ests_np$cr$est[50], 0.08646103, tolerance=0.01)
  # expect_equal(ests_np$cr$se[1], 999, tolerance=0.01)
  # expect_equal(ests_np$cr$se[50], 999, tolerance=0.01)
  # expect_equal(ests_np$cr$se[101], 999, tolerance=0.01)
  expect_equal(ests_np$cr$ci_lower[1], 0.161259, tolerance=0.01)
  expect_equal(ests_np$cr$ci_lower[50], 0.06596329, tolerance=0.01)
  expect_equal(ests_np$cr$ci_upper[1], 0.4155087, tolerance=0.01)
  expect_equal(ests_np$cr$ci_upper[50], 0.1151884, tolerance=0.01)
})

test_that("est_np (CVE)", {
  expect_equal(ests_np$cve$s[1], 0, tolerance=0.01)
  expect_equal(ests_np$cve$s[50], 1.15447, tolerance=0.01)
  expect_equal(ests_np$cve$est[1], -8.348665, tolerance=0.01)
  expect_equal(ests_np$cve$est[50], -2.002264, tolerance=0.01)
  expect_equal(ests_np$cve$se[1], 3.161112, tolerance=0.01)
  expect_equal(ests_np$cve$se[50], 0.8215806, tolerance=0.01)
  expect_equal(ests_np$cve$ci_lower[1], -17.13744, tolerance=0.01)
  expect_equal(ests_np$cve$ci_lower[50], -4.133193, tolerance=0.01)
  expect_equal(ests_np$cve$ci_upper[1], -3.837045, tolerance=0.01)
  expect_equal(ests_np$cve$ci_upper[50], -0.7789708, tolerance=0.01)
})

p <- plot_ce(ests_np)

test_that("plot_ce", {
  expect_equal(class(p), c("gg", "ggplot"))
})
