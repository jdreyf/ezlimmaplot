context("plot_man")

test_that("plot man vdiffr", {
  pm <- function() plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  vdiffr::expect_doppelganger("g1", pm)
})

test_that("plot man non-vdiffr", {
  pm <- plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  # test pm$data
})
