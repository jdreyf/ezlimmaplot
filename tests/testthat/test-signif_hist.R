context("signif_hist")

df <- data.frame(p=1:4, FDR=2:5, p=0:3, check.names = FALSE)
df2 <- data.frame(a.p=1:4, a.FDR=2:5, b.p=0:3)

test_that("signif_hist testthat", {
  expect_equal(colnames(signif_hist(tab=res.df, plot=FALSE)), grep("p|FDR", colnames(res.df), value = TRUE))
  expect_error(signif_hist(tab=res.df, sep="_", plot=FALSE))

  expect_equal(signif_hist(tab=df[,1:2], sep=NA, plot=FALSE), df[,1:2])
  expect_error(signif_hist(tab=df, plot=FALSE))

  expect_error(signif_hist(tab=df2, plot=FALSE))
})

test_that("sh vdiffr", {
  res.ss <- res.df[,grep("^(First3|Last3)\\.", colnames(res.df))]
  sh <- function() res <- signif_hist(tab=res.ss, pi0 = TRUE, name=NA)
  expect_doppelganger(title="sh", fig=sh)
})

test_that("error messages fires if duplicate p-value column names, or no p-value columns", {
  res.df2 <- cbind(res.df, First3.p = res.df[,"First3.p"])
  expect_error(signif_hist(tab=res.df2, plot=FALSE))

  # The test on line 28 isn't needed because extract_prefix will catch the lack of p-val columns before you hit this line.
  # res.df3 <- res.df[,-c(3,7,11)] #removes all p-val columns
  # expect_error(signif_hist(tab=res.df3, plot=FALSE))
})
