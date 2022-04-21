context("ezpca")

test_that("ezpca non-vdiffr", {
  #verified that this looks like multi.pca
  ezp <- ezpca(M, pheno, shape="grp", name=NA, manual.shape = 1:2, labels = TRUE, plot=FALSE)$data
  expect_lte(ezp["sample1", "PC1"], -2)
  expect_equal(ezp["sample1", "grp"], "First3")

  pheno.df2 <- data.frame(pheno, "covar num"=pheno$covar_num + 1, check.names = FALSE)
  ezp.df <- ezpca(object=M, pheno.df=pheno.df2, color="`covar num`", name=NA, plot=FALSE)$data
  expect_equal(ezp.df$covar_num, pheno.df2[rownames(ezp.df), "covar_num"])
})

test_that("ezpca without pheno", {
  expect_silent(ezp <- ezpca(M, labels = TRUE))
})

test_that("ezpca w/ ellipses", {
  M2 <- matrix(rnorm(100*12, sd=0.3), nrow=100, ncol=12)
  dimnames(M2) <- list(paste0("gene", 1:nrow(M2)), paste0("sample", 1:ncol(M2)))
  grp2 <- rep(c("First3", "Last3"), each=ncol(M2)/2)
  tissue2 <- rep(c("muscle", "liver"), times=ncol(M2)/2)
  covar_num2 <- rnorm(n = ncol(M2))
  pheno2 <- data.frame(row.names = colnames(M2), sample=colnames(M2), grp2, tissue2, covar_num2, stringsAsFactors = FALSE)
  M2[1, 1:(ncol(M2)/2)] <- M2[1, 1:(ncol(M2)/2)] + 2

  expect_silent(ezp <- ezpca(M2, pheno2, shape="grp2", name=NA, color="grp2", manual.shape = 1:2, labels = TRUE,
                             ellipses = TRUE, plot=FALSE))
})




test_that("ezpca vdiffr", {
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
