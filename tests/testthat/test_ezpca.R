context("ezpca")

test_that("ezpca", {
  #verified that this looks like multi.pca
  ezp <- function() ezpca(M, pheno, color="grp", name=NA)
  ezpl <- function() ezpca(M, pheno, color="grp", name=NA, labels=TRUE)
  ezps <- function() ezpca(M, pheno, color="grp", name=NA, all.size = 3)

  #addins ...
  vdiffr::expect_doppelganger("pca", ezp)
  vdiffr::expect_doppelganger("pca-labels", ezpl)
  vdiffr::expect_doppelganger("pca-size", ezps)
})
