context("ezheat")

test_that("mat", {
  mat0 <- ezheat(M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1,
                 plot=FALSE)$mat
  expect_gt(max(mat0[1,1:2]), 0.75)
  expect_lt(min(mat0[1,4:6]), -0.5)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_equal(rownames(mat0)[11],  topgenes[11])
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- ezheat(M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA,
                                clip=1, only.labrows = TRUE, plot=FALSE, verbose=TRUE)$mat)
  expect_equal(nrow(mat1), sum(!is.na(res.df$Gene.Symbol[1:20])))
})

test_that("sc", {
  ezh.sc <- ezheat(M, pheno.df=pheno.df, name=NA, sc="z", plot=FALSE)$mat
  #same as pheatmap:::scale_mat(M, scale = "row")
  base.sc <- t(scale(t(M), center=TRUE, scale = TRUE))
  expect_equal(ezh.sc, base.sc)
})

test_that("vdiffr", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  ezh <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA)$gtable
  }
  vdiffr::expect_doppelganger(title="heat", fig=ezh)

  ezh2 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA, sc="z")$gtable
  }
  vdiffr::expect_doppelganger(title="heat2", fig=ezh2)

  #all avg+10 > 0
  ezh3 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(object=res.df[topgenes, 1:2]+10, labrows = lr, name=NA, sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")])$gtable
  }
  vdiffr::expect_doppelganger(title="heat3", fig=ezh3)

  ezh4 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE)$gtable
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat4", fig=ezh4)

  ezh5 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE, labcols=letters[1:6])$gtable
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat5", fig=ezh5)
})

test_that("vdiff low res", {
  # low res heat w/ clusters & w/o rowlabs
  cl.df <- data.frame(cl=LETTERS[as.numeric(sub("gene", "", rownames(M))) %% 5 + 1])
  rownames(cl.df) <- rownames(M)

  ezh6 <- function(){
    ret <- ezheat(M, pheno.df=pheno.df, labrows = "", name=NA, annotation_row = cl.df)
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat6", fig=ezh6)
})
