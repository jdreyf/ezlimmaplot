context("ezpca")

pheno2 <- data.frame(pheno, tissue=rep(c("muscle", "liver"), 3))

test_that("ezpca", {
  #verified that this looks like multi.pca
  ezp <- function() ezpca(M, pheno, shape="grp", name=NA, manual.shape = 1:2)
  ezpl <- function() ezpca(M, pheno, color="grp", name=NA, labels=TRUE, manual.color=c("red", "blue"))
  ezps <- function() ezpca(M, pheno, color="grp", name=NA, all.size = 3, rm.leg.title=TRUE)
  ezp2 <- function() ezpca(M, pheno2, color="grp", name=NA, facet = ". ~ tissue")

  #addins ...
  vdiffr::expect_doppelganger(title="pca", fig=ezp)
  vdiffr::expect_doppelganger(title="pca-labels", fig=ezpl)
  vdiffr::expect_doppelganger(title="pca-size", fig=ezps)
  vdiffr::expect_doppelganger(title="pca2", fig=ezp2)
})
