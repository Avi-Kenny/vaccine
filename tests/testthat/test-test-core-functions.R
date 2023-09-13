
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)

test_that("load_data", {
  expect_equal(class(dat), "vaccine_dat")
  expect_equal(attr(dat, "groups"), "both")
  expect_equal(attr(dat, "covariate_names"), c("age", "BMI", "bhvrisk"))

})

est_overall(dat=dat, t_0=578, method="KM")

test_that("est_overall", {
  expect_equal(2*2, 4)
})
