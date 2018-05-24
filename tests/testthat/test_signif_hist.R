context("signif_hist")

test_that("signif_hist", {
  df <- data.frame(p=1:4, FDR=2:5, p=0:3, check.names = FALSE)
  df2 <- data.frame(a.p=1:4, a.FDR=2:5, b.p=0:3)

  expect_equal(colnames(signif_hist(tab=res.df, plot=FALSE)), grep("p|FDR", colnames(res.df), value = TRUE))
  expect_error(signif_hist(tab=res.df, sep="_", plot=FALSE))

  expect_equal(signif_hist(tab=df[,1:2], plot=FALSE), df[,1:2])
  expect_error(signif_hist(tab=df, plot=FALSE))

  expect_error(signif_hist(tab=df2, plot=FALSE))
})
