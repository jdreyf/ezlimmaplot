context("mean se")

test_that("mean_se", {
  expect_error(mean_se(M))
  expect_error(mean_se(rep(NA, 3)))

  x <- c(cumsum(1:4), NA)
  z <- round(mean_se(x))
  expect_equal(z, c(y=5, ymin=3, ymax=7))
})
