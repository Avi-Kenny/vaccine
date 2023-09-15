
t <- 1e-4

data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)

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

ests <- est_overall(dat=dat, t_0=578, method="KM")
ests_risk_p <- ests[ests$stat=="risk" & ests$group=="placebo",]
ests_risk_v <- ests[ests$stat=="risk" & ests$group=="vaccine",]
ests_ve <- ests[ests$stat=="ve",]

test_that("est_overall (KM)", {
  expect_equal(class(ests), c("data.frame", "vaccine_overall"))
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

ests <- est_overall(dat=dat, t_0=578, method="Cox")
ests_risk_p <- ests[ests$stat=="risk" & ests$group=="placebo",]
ests_risk_v <- ests[ests$stat=="risk" & ests$group=="vaccine",]
ests_ve <- ests[ests$stat=="ve",]

test_that("est_overall (Cox)", {
  expect_equal(class(ests), c("data.frame", "vaccine_overall"))
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
