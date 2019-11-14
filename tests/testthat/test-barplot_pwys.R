context("barplot of significance per pathway")

test_that("direction properly specified", {
  expect_error(barplot_pwys(rc, prefix.v = "First3"))
  expect_silent(barplot_pwys(tab=rc, prefix.v = "First3", direction = "Up"))
})
