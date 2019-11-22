context("plot_man")

test_that("plot man non-vdiffr", {
  pm <- plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  # test pm$data: To check the location of the point
  expect_equal(pm$data["sample4", "Exprs"], M["gene1", "sample4"])
  # test pm$data
  expect_equal(as.character(pm$data["sample4", "Exposure"]), pheno["sample4", "grp"])
  expect_equal(pm$labels$x, "Log2 abundance")
  expect_equal(pm$labels$y, "Covar")
  expect_equal(pm$labels$title, row.names(M[1,,drop=FALSE]))
  expect_equal(row.names(pm$data), colnames(M[1,,drop=FALSE]))
  expect_equal(pm$data$Outcome, pheno$covar_num)
  expect_equal(as.character(pm$data$Exposure), grp)

  #test for different labels and title
  pm <- plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, xlab = "Covar", main.v="MyTitle", name=NA)
  expect_equal(pm$labels$title, "MyTitle")
  expect_equal(pm$labels$x, "Covar")

  #test for E>2	and color
  pm <- plot_man(E=c("A1","A1","B2","B2","C3","C3"), M=M[1,,drop=FALSE], Y=pheno$covar_num, manual.color=c("red","orange","cyan"), name=NA)
  p.obj <- ggplot_build(pm)
  expect_equal(as.character(p.obj$data[[1]]$colour), c("red","red","orange","orange","cyan","cyan"))
  expect_equal(as.character(pm$data$Exposure), c("A1","A1","B2","B2","C3","C3"))

  #test for manual.color list not matching with factors in E
  expect_error(plot_man(E=c("A1","A1","B2","B2","C3","C3"), M=M[1,,drop=FALSE], Y=pheno$covar_num, manual.color=c("red","orange"), name=NA))
})

test_that("plot man vdiffr", {
  pm <- function() plot_man(E=grp, M=M[1,,drop=FALSE], Y=pheno$covar_num, ylab = "Covar", name=NA)
  vdiffr::expect_doppelganger("g1", pm)
})

