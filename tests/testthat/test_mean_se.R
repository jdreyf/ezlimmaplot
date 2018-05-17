context("mean se")

test_that("mean_se", {
  expect_error(mean_se(M))
  expect_error(mean_se(rep(NA, 3)))
})
