test_that("dotplot_pwy works", {
  expect_silent(dotplot_pwys(tab=rc2, name = "pwy", colorbar.title = "P<0.05"))
  expect_silent(dotplot_pwys(rc2, mixed="exclude", name="pwy"))
  expect_silent(dotplot_pwys(rc2, mixed="only", name="pwy"))
  rc.nomix <- rc2 |> dplyr::select(!contains("Mixed"))
  expect_silent(dotplot_pwys(rc.nomix, mixed="exclude", name="pwy"))
  # make sure dotplot_pwys properly orders by p
  res <- rc2 |> dplyr::arrange(desc(First3.p)) |>
    dotplot_pwys(ntop=1)
  expect_equal(unique(res$data$Pwy), "pwy1")
  res2 <- rc2 |> dplyr::arrange(desc(First3.p)) |>
    dotplot_pwys(ntop=1, reorder.rows = FALSE)
  expect_equal(unique(res2$data$Pwy), "mother of all pathways")
  unlink(x = "pwy_dotplot.pdf")

  expect_error(dotplot_pwys(rc2, name = NA, cut.sig = 10**-6))
  expect_error(dotplot_pwys(rc2, name = NA, type.sig = "FDR", cut.sig = 10**-6))
})
