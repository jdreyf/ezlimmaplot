context("ezpca")

test_that("ezpca", {
  #verified that this looks like multi.pca
  ezp <- function() ezpca(M, pheno, color="grp", name=NA)
  ezpl <- function() ezpca(M, pheno, color="grp", name=NA, labels=TRUE)
  ezps <- function() ezpca(M, pheno, color="grp", name=NA, all.size = 3)

  #addins ...
  vdiffr::expect_doppelganger(title="pca", fig=ezp)
  vdiffr::expect_doppelganger(title="pca-labels", fig=ezpl)
  vdiffr::expect_doppelganger(title="pca-size", fig=ezps)
})
