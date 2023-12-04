test_that("bubbleplot_pwy works", {
  expect_silent(bubbleplot_pwys(rc2, name = "pwy", cut.sig=1))
  rc.nomix <- rc2 |> dplyr::select(!contains("Mixed"))
  expect_silent(bubbleplot_pwys(rc.nomix, mixed="exclude", name="pwy", cut.sig=1))
  unlink(x = "pwy_bubbleplots.pdf")

  expect_error(bubbleplot_pwys(rc2, name = NA, cut.sig = 10**-6))
  expect_error(bubbleplot_pwys(rc2, name = NA, type.sig = "FDR", cut.sig = 10**-6))
})
