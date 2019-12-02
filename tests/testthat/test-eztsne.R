context("eztsne")

test_that("ez tsne non-vdiffr", {
  eztsne.non.vdiffr <- eztsne(M, pheno, shape="grp", manual.shape = 1:2, labels = TRUE, plot=FALSE)$data
  expect_lte(eztsne.non.vdiffr["sample1", "tSNE1"], 229)
  expect_equal(eztsne.non.vdiffr["sample1", "grp"], "First3")

  pheno.df2 <- data.frame(pheno, "covar.num"=pheno$covar_num + 1, check.names = FALSE)
  eztsne.df <- eztsne(object=M, pheno.df=pheno.df2, color="covar.num", plot=FALSE)$data
  expect_equal(eztsne.df$covar_num, pheno.df2[rownames(eztsne.df), "covar_num"])
})

test_that("ez tsne non vdiffr plot object", {
  ez.tsne <- eztsne(M, pheno, shape="grp", manual.shape = 1:2)
  ez.tsne.obj <- ggplot_build(ez.tsne)
  expect_equal(ez.tsne.obj$data[[1]]$shape, c(1,1,1,2,2,2))

  eztsnel <- eztsne(M, pheno, color="grp", labels=TRUE, manual.color=c("red", "blue"))
  eztsnel.obj <- ggplot_build(eztsnel)
  expect_equal(eztsnel.obj$data[[1]]$colour, c("red","red","red","blue","blue","blue"))

  eztsnes <- eztsne(M, pheno, color="grp", all.size = 3, rm.leg.title=TRUE)
  eztsnes.obj <- ggplot_build(eztsnes)
  expect_equal(eztsnes.obj$data[[1]]$size, rep(3,times=6))
  expect_equal(length(eztsnes$theme$legend.title), 0)

  eztsne2 <- eztsne(M, pheno2, color="grp", facet = ". ~ tissue")
  eztsne2.obj <- ggplot_build(eztsne2)
  expect_equal(as.character(eztsne2.obj$data[[1]]$PANEL), rep(c("2","1"), times=3))

  eztsne3 <- eztsne(M, pheno, shape="tissue", color="covar_num")
  eztsne3.obj <- ggplot_build(eztsne3)
  expect_equal(eztsne3.obj$data[[1]]$shape, rep(c(17,16), times=3))
  expect_equal(eztsne3.obj$data[[1]]$colour[1], "#1E4161")

  eztsne.tit <- eztsne(M, pheno, shape="grp", title="main")
  eztsne.tit.obj <- ggplot_build(eztsne.tit)
  expect_equal(eztsne.tit.obj$data[[1]]$shape, c(16,16,16,17,17,17))
  expect_equal(eztsne.tit$labels$title, "main")

  eztsne.tit2 <- eztsne(M, pheno, shape="grp", title=NULL)
  expect_null(eztsne.tit2$labels$title)

  eztsne.tit3 <- eztsne(M, pheno, shape="grp", title=NA)
  expect_true(is.na(eztsne.tit3$labels$title))

  eztsne.stit <- eztsne(M, pheno, shape="grp", title="main", subtitle="submain")
  expect_equal(eztsne.stit$labels$title, "main")
  expect_equal(eztsne.stit$labels$subtitle, "submain")

  eztsne.stit2 <- eztsne(M, pheno, shape="grp", title=NA, subtitle="submain")
  expect_true(is.na(eztsne.stit2$labels$title))
  expect_equal(eztsne.stit2$labels$subtitle, "submain")

  eztsne.stit3 <- eztsne(M, pheno, shape="grp", title="", subtitle="submain")
  expect_equal(eztsne.stit3$labels$title, c(""))
  expect_equal(eztsne.stit3$labels$subtitle, "submain")

  eztsnenl <- eztsne(M, pheno, shape="grp", labels=FALSE)
  expect_null(eztsnenl$labels$label)

  eztsnel <- eztsne(M, pheno, shape="grp", labels=TRUE)
  expect_equal(eztsnel$labels$label, "row_names")
})

#test_that("ez tsne without pheno", {
#  expect_silent(eztsne(M, pheno ,labels = TRUE))
#})

test_that("eztsne vdiffr", {
   ez.tsne <- function() eztsne(M, pheno, shape="grp", manual.shape = 1:2)
   eztsnel <- function() eztsne(M, pheno, color="grp", labels=TRUE, manual.color=c("red", "blue"))
   eztsnes <- function() eztsne(M, pheno, color="grp", all.size = 3, rm.leg.title=TRUE)
   eztsne2 <- function() eztsne(M, pheno2, color="grp", facet = ". ~ tissue")
   eztsne3 <- function() eztsne(M, pheno, shape="tissue", color="covar_num")

   eztsne.tit <- function() eztsne(M, pheno, shape="grp", title="main")
   eztsne.tit2 <- function() eztsne(M, pheno, shape="grp", title=NULL)
   eztsne.tit3 <- function() eztsne(M, pheno, shape="grp", title=NA)
   eztsne.stit <- function() eztsne(M, pheno, shape="grp", title="main", subtitle="submain")
   eztsne.stit2 <- function() eztsne(M, pheno, shape="grp", title=NA, subtitle="submain")
   eztsne.stit3 <- function() eztsne(M, pheno, shape="grp", title="", subtitle="submain")

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

teardown(unlink("tsne.pdf"))
