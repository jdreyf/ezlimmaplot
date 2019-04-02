context("ezpca")

test_that("ezpca", {
  #verified that this looks like multi.pca
  ezp <- function() ezpca(M, pheno, shape="grp", name=NA, manual.shape = 1:2)
  ezpl <- function() ezpca(M, pheno, color="grp", name=NA, labels=TRUE, manual.color=c("red", "blue"))
  ezps <- function() ezpca(M, pheno, color="grp", name=NA, all.size = 3, rm.leg.title=TRUE)
  ezp2 <- function() ezpca(M, pheno2, color="grp", name=NA, facet = ". ~ tissue")
  ezp3 <- function() ezpca(M, pheno, shape="tissue", color="covar_num", name=NA)

  ezp.tit <- function() ezpca(M, pheno, shape="grp", name=NA, title="main")
  ezp.tit2 <- function() ezpca(M, pheno, shape="grp", name=NA, title=NULL)
  ezp.tit3 <- function() ezpca(M, pheno, shape="grp", name=NA, title=NA)
  ezp.stit <- function() ezpca(M, pheno, shape="grp", name=NA, title="main", subtitle="submain")
  ezp.stit2 <- function() ezpca(M, pheno, shape="grp", name=NA, title=NA, subtitle="submain")
  ezp.stit3 <- function() ezpca(M, pheno, shape="grp", name=NA, title="", subtitle="submain")

  pheno.df2 <- data.frame(pheno, "covar num"=pheno$covar_num, check.names = FALSE)
  ezp.df <- function() ezpca(object=M, pheno.df=pheno.df2, color="`covar num`", name=NA)

  #addins ...
  vdiffr::expect_doppelganger(title="pca", fig=ezp)
  vdiffr::expect_doppelganger(title="pca-labels", fig=ezpl)
  vdiffr::expect_doppelganger(title="pca-size", fig=ezps)
  vdiffr::expect_doppelganger(title="pca2", fig=ezp2)
  vdiffr::expect_doppelganger(title="pca3", fig=ezp3)

  vdiffr::expect_doppelganger(title="pca.tit", fig=ezp.tit)
  vdiffr::expect_doppelganger(title="pca.tit2", fig=ezp.tit2)
  vdiffr::expect_doppelganger(title="pca.tit3", fig=ezp.tit3)
  vdiffr::expect_doppelganger(title="pca.stit", fig=ezp.stit)
  vdiffr::expect_doppelganger(title="pca.stit2", fig=ezp.stit2)
  vdiffr::expect_doppelganger(title="pca.stit3", fig=ezp.stit3)

  vdiffr::expect_doppelganger(title="pca.df", fig=ezp.df)
})

test_that("ezpca without pheno", {
  expect_silent(ezp <- ezpca(M, labels = TRUE))
})
