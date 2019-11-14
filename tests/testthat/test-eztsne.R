context("ez tsne")

test_that("ez tsne vdiffr", {
   ez.tsne <- function() eztsne(M, pheno, shape="grp", name=NA, manual.shape = 1:2)	
   eztsnel <- function() eztsne(M, pheno, color="grp", name=NA, labels=TRUE, manual.color=c("red", "blue"))
   eztsnes <- function() eztsne(M, pheno, color="grp", name=NA, all.size = 3, rm.leg.title=TRUE)
   eztsne2 <- function() eztsne(M, pheno2, color="grp", name=NA, facet = ". ~ tissue")
   eztsne3 <- function() eztsne(M, pheno, shape="tissue", color="covar_num", name=NA)

   eztsne.tit <- function() eztsne(M, pheno, shape="grp", name=NA, title="main")
   eztsne.tit2 <- function() eztsne(M, pheno, shape="grp", name=NA, title=NULL)
   eztsne.tit3 <- function() eztsne(M, pheno, shape="grp", name=NA, title=NA)
   eztsne.stit <- function() eztsne(M, pheno, shape="grp", name=NA, title="main", subtitle="submain")
   eztsne.stit2 <- function() eztsne(M, pheno, shape="grp", name=NA, title=NA, subtitle="submain")
   eztsne.stit3 <- function() eztsne(M, pheno, shape="grp", name=NA, title="", subtitle="submain")
  
  vdiffr::expect_doppelganger(title="tsne", fig=ez.tsne)
  vdiffr::expect_doppelganger(title="tsne-labels", fig=eztsnel)
  vdiffr::expect_doppelganger(title="tsne-size", fig=eztsnes)
  vdiffr::expect_doppelganger(title="tsne2", fig=eztsne2)
  vdiffr::expect_doppelganger(title="tsne3", fig=eztsne3)

  vdiffr::expect_doppelganger(title="tsne.tit", fig=eztsne.tit)
  vdiffr::expect_doppelganger(title="tsne.tit2", fig=eztsne.tit2)
  vdiffr::expect_doppelganger(title="tsne.tit3", fig=eztsne.tit3)
  vdiffr::expect_doppelganger(title="tsne.stit", fig=eztsne.stit)
  vdiffr::expect_doppelganger(title="tsne.stit2", fig=eztsne.stit2)
  vdiffr::expect_doppelganger(title="tsne.stit3", fig=eztsne.stit3)
})

test_that("ez tsne non-vdiffr", {
  eztsne.non.vdiffr <- eztsne(M, pheno, shape="grp", name=NA, manual.shape = 1:2, labels = TRUE, plot=FALSE)
  expect_lte(eztsne.non.vdiffr["sample1", "tSNE1"], 229)
  expect_equal(eztsne.non.vdiffr["sample1", "grp"], "First3")

  pheno.df2 <- data.frame(pheno, "covar.num"=pheno$covar_num + 1, check.names = FALSE)
  eztsne.df <- eztsne(object=M, pheno.df=pheno.df2, color="covar.num", name=NA, plot=FALSE)
  expect_equal(eztsne.df$covar_num, pheno.df2[rownames(eztsne.df), "covar_num"])
})

#test_that("ez tsne without pheno", {
#  expect_silent(eztsne(M, pheno ,labels = TRUE))
#})





