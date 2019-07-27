context("plot_man")

test_that("plot man vdiffr", {
  pm <- function() plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  vdiffr::expect_doppelganger("g1", pm)
})

test_that("plot man non-vdiffr", {
  pm <- plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  # test pm$data: To check the location of the point
  expect_equal(pm$data["sample4", "Exprs"], M["gene1", "sample4"])
  # test pm$data
  expect_equal(as.character(pm$data["sample4", "Exposure"]), pheno["sample4", "grp"])
})
