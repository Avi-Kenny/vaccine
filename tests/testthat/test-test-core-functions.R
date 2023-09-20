
t <- 1e-4
data(hvtn505)
set.seed(1)
hvtn505_sample <- hvtn505[sample(c(1:nrow(hvtn505)), size=500),]
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

test_that("load_data", {
  expect_equal(class(dat), "vaccine_dat")
  expect_equal(attr(dat, "groups"), "both")
  expect_equal(attr(dat, "covariate_names"), c("age", "BMI", "bhvrisk"))
  expect_equal(class(dat$v), "list")
  expect_equal(class(dat$p), "list")
  expect_equal(sum(dat$v$y), 391608)
  expect_equal(sum(dat$p$y), 380935)
  expect_equal(unique(dat$v$delta), c(0,1))
  expect_equal(unique(dat$v$strata), c(1,2))
  expect_equal(unique(dat$v$z), c(0,1))
  expect_equal(attr(dat$v, "n_orig"), 1161)
  expect_equal(attr(dat$v, "dim_x"), 3)
})

ests_o_KM <- est_overall(dat=dat, t_0=578, method="KM")
ests_risk_p <- ests_o_KM[ests_o_KM$stat=="risk"&ests_o_KM$group=="placebo",]
ests_risk_v <- ests_o_KM[ests_o_KM$stat=="risk"&ests_o_KM$group=="vaccine",]
ests_ve <- ests_o_KM[ests_o_KM$stat=="ve",]

test_that("est_overall (KM)", {
  expect_equal(class(ests_o_KM), c("data.frame", "vaccine_overall"))
  expect_equal(ests_risk_p$est, 0.02879861, tolerance=t)
  expect_equal(ests_risk_p$se, 0.006563785, tolerance=t)
  expect_equal(ests_risk_p$ci_lower, 0.0162236, tolerance=t)
  expect_equal(ests_risk_p$ci_upper, 0.04121288, tolerance=t)
  expect_equal(ests_risk_v$est, 0.04067009, tolerance=t)
  expect_equal(ests_risk_v$se, 0.008230842, tolerance=t)
  expect_equal(ests_risk_v$ci_lower, 0.02506853, tolerance=t)
  expect_equal(ests_risk_v$ci_upper, 0.05602199, tolerance=t)
  expect_equal(ests_ve$est, -0.4122241, tolerance=t)
  expect_equal(ests_ve$se, 0.4304518, tolerance=t)
  expect_equal(ests_ve$ci_lower, -1.5666, tolerance=t)
  expect_equal(ests_ve$ci_upper, 0.2229498, tolerance=t)
})

ests_o_Cox <- est_overall(dat=dat, t_0=578, method="Cox")
ests_risk_p <- ests_o_Cox[ests_o_Cox$stat=="risk"&ests_o_Cox$group=="placebo",]
ests_risk_v <- ests_o_Cox[ests_o_Cox$stat=="risk"&ests_o_Cox$group=="vaccine",]
ests_ve <- ests_o_Cox[ests_o_Cox$stat=="ve",]

test_that("est_overall (Cox)", {
  expect_equal(class(ests_o_Cox), c("data.frame", "vaccine_overall"))
  expect_equal(ests_risk_p$est, 0.02938706, tolerance=t)
  expect_equal(ests_risk_p$se, 0.006486545, tolerance=t)
  expect_equal(ests_risk_p$ci_lower, 0.0190193, tolerance=t)
  expect_equal(ests_risk_p$ci_upper, 0.04514638, tolerance=t)
  expect_equal(ests_risk_v$est, 0.04177642, tolerance=t)
  expect_equal(ests_risk_v$se, 0.008111679, tolerance=t)
  expect_equal(ests_risk_v$ci_lower, 0.02847302, tolerance=t)
  expect_equal(ests_risk_v$ci_upper, 0.06090588, tolerance=t)
  expect_equal(ests_ve$est, -0.4215925, tolerance=t)
  expect_equal(ests_ve$se, 0.4179152, tolerance=t)
  expect_equal(ests_ve$ci_lower, -1.529375, tolerance=t)
  expect_equal(ests_ve$ci_upper, 0.201018, tolerance=t)
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
  expect_equal(ests_med_NP_nde$est, 20.35041, tolerance=t)
  expect_equal(ests_med_NP_nde$se, 10.57582, tolerance=t)
  expect_equal(ests_med_NP_nde$ci_lower, 7.348651, tolerance=t)
  expect_equal(ests_med_NP_nde$ci_upper, 56.35583, tolerance=t)
  expect_equal(ests_med_NP_nie$est, 1.03442, tolerance=t)
  expect_equal(ests_med_NP_nie$se, 0.2177699, tolerance=t)
  expect_equal(ests_med_NP_nie$ci_lower, 0.6846922, tolerance=t)
  expect_equal(ests_med_NP_nie$ci_upper, 1.562782, tolerance=t)
  expect_equal(ests_med_NP_pm$est, 0.01110645, tolerance=t)
  expect_equal(ests_med_NP_pm$se, 0.068966, tolerance=t)
  expect_equal(ests_med_NP_pm$ci_lower, -0.1240669, tolerance=t)
  expect_equal(ests_med_NP_pm$ci_upper, 0.1462798, tolerance=t)
})

set.seed(1)
ests_cox <- est_ce(dat=dat, type="Cox", t_0=578, cve=T)

test_that("est_cox (CR)", {
  expect_equal(class(ests_cox), "vaccine_est")
  expect_equal(ests_cox$cr$s[1], 0, tolerance=t)
  expect_equal(ests_cox$cr$s[50], 1.15447037, tolerance=t)
  expect_equal(ests_cox$cr$est[1], 0.15486178, tolerance=t)
  expect_equal(ests_cox$cr$est[50], 0.08759399, tolerance=t)
  expect_equal(ests_cox$cr$se[1], 0.05514197, tolerance=t)
  expect_equal(ests_cox$cr$se[50], 0.01915460, tolerance=t)
  expect_equal(ests_cox$cr$ci_lower[1], 0.07427856, tolerance=t)
  expect_equal(ests_cox$cr$ci_lower[50], 0.05661917, tolerance=t)
  expect_equal(ests_cox$cr$ci_upper[1], 0.2950081, tolerance=t)
  expect_equal(ests_cox$cr$ci_upper[50], 0.1331231, tolerance=t)
})

test_that("est_cox (CVE)", {
  expect_equal(ests_cox$cve$s[1], 0, tolerance=t)
  expect_equal(ests_cox$cve$s[50], 1.15447037, tolerance=t)
  expect_equal(ests_cox$cve$est[1], -4.3774047, tolerance=t)
  expect_equal(ests_cox$cve$est[50], -2.0416048, tolerance=t)
  expect_equal(ests_cox$cve$se[1], 2.2734089, tolerance=t)
  expect_equal(ests_cox$cve$se[50], 0.9607152, tolerance=t)
  expect_equal(ests_cox$cve$ci_lower[1], -11.315226, tolerance=t)
  expect_equal(ests_cox$cve$ci_lower[50], -4.648935, tolerance=t)
  expect_equal(ests_cox$cve$ci_upper[1], -1.3480269019, tolerance=t)
  expect_equal(ests_cox$cve$ci_upper[50], -0.6377176338, tolerance=t)
})

set.seed(1)
ests_np <- est_ce(dat=dat, type="NP", t_0=578, cve=T,
                  params_np=params_ce_np(surv_type="Cox"))

test_that("est_np (CR)", {
  expect_equal(class(ests_cox), "vaccine_est")
  expect_equal(ests_np$cr$s[1], 0, tolerance=t)
  expect_equal(ests_np$cr$s[50], 1.15447, tolerance=t)
  expect_equal(ests_np$cr$est[1], 0.2692285, tolerance=t)
  expect_equal(ests_np$cr$est[50], 0.08646103, tolerance=t)
  # expect_equal(ests_np$cr$se[1], 999, tolerance=t)
  # expect_equal(ests_np$cr$se[50], 999, tolerance=t)
  # expect_equal(ests_np$cr$se[101], 999, tolerance=t)
  expect_equal(ests_np$cr$ci_lower[1], 0.161259, tolerance=t)
  expect_equal(ests_np$cr$ci_lower[50], 0.06596329, tolerance=t)
  expect_equal(ests_np$cr$ci_upper[1], 0.4155087, tolerance=t)
  expect_equal(ests_np$cr$ci_upper[50], 0.1151884, tolerance=t)
})

test_that("est_np (CVE)", {
  expect_equal(ests_np$cve$s[1], 0, tolerance=t)
  expect_equal(ests_np$cve$s[50], 1.15447, tolerance=t)
  expect_equal(ests_np$cve$est[1], -8.348665, tolerance=t)
  expect_equal(ests_np$cve$est[50], -2.002264, tolerance=t)
  expect_equal(ests_np$cve$se[1], 3.161112, tolerance=t)
  expect_equal(ests_np$cve$se[50], 0.8215806, tolerance=t)
  expect_equal(ests_np$cve$ci_lower[1], -17.13744, tolerance=t)
  expect_equal(ests_np$cve$ci_lower[50], -4.133193, tolerance=t)
  expect_equal(ests_np$cve$ci_upper[1], -3.837045, tolerance=t)
  expect_equal(ests_np$cve$ci_upper[50], -0.7789708, tolerance=t)
})
