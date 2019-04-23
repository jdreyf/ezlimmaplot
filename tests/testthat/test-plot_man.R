context("plot_man")

test_that("plot man", {
  pm <- function() plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  vdiffr::expect_doppelganger("g1", pm)
})
