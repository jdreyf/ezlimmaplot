test_that("dotplot_pwy works", {
  expect_silent(dotplot_pwys(rc2, name = NA))
  expect_error(dotplot_pwys(rc2, name = NA, cut.sig = 10**-6))
  expect_error(dotplot_pwys(rc2, name = NA, type.sig = "FDR", cut.sig = 10**-6))
})
