context("ezheat")

test_that("mat", {
  mat0 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1, plot=FALSE)$mat
  expect_gt(max(mat0[1,1:2]), 0.75)
  expect_lt(min(mat0[1,4:6]), -0.5)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_equal(rownames(mat0)[11],  topgenes[11])
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA,
                                clip=1, only.symbols = TRUE, plot=FALSE)$mat)
  expect_equal(nrow(mat1), sum(!is.na(res.df$Gene.Symbol[1:20])))
})

test_that("gtable", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  ezh <- function(){
    ret <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA)$gtable
  }
  vdiffr::expect_doppelganger(title="heat", fig=ezh)

  ezh2 <- function(){
    ret <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno2, name=NA, sc="z")$gtable
  }
  vdiffr::expect_doppelganger(title="heat2", fig=ezh2)

  #all avg+10 > 0
  ezh3 <- function(){
    ret <- ezheat(object=res.df[topgenes, 1:2]+10, symbols = res.df[topgenes, "Gene.Symbol"], name=NA, sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")])$gtable
  }
  vdiffr::expect_doppelganger(title="heat3", fig=ezh3)

  ezh4 <- function(){
    ret <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno2, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE)$gtable
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat4", fig=ezh4)
})
