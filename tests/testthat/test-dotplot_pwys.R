test_that("dotplot_pwy works", {
  expect_silent(dotplot_pwys(tab=rc2, name = "pwy", colorbar.title = "P<0.05"))
  expect_silent(dotplot_pwys(rc2, mixed="exclude", name="pwy"))
  expect_silent(dotplot_pwys(rc2, mixed="only", name="pwy"))
  rc.nomix <- rc2 |> dplyr::select(!contains("Mixed"))
  expect_silent(dotplot_pwys(rc.nomix, mixed="exclude", name="pwy"))
  unlink(x = "pwy_dotplot.pdf")

  expect_error(dotplot_pwys(rc2, name = NA, cut.sig = 10**-6))
  expect_error(dotplot_pwys(rc2, name = NA, type.sig = "FDR", cut.sig = 10**-6))
})
