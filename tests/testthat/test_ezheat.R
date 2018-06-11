context("ezheat")

pheno.df <- data.frame(row.names = rownames(pheno), grp=pheno$grp)
topgenes <- rownames(res.df)[1:20]

test_that("ezheat", {
  mat0 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1, plot=FALSE)
  expect_gt(max(mat0[1,1:2]), 0.75)
  expect_lt(min(mat0[1,4:6]), -0.5)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_equal(rownames(mat0)[11],  topgenes[11])
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA,
                                clip=1, only.symbols = TRUE, plot=FALSE))
  expect_equal(nrow(mat1), sum(!is.na(res.df$Gene.Symbol[1:20])))
})

# test_that("ezh viz", {
#   ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name="../figs/ezheat/heat")
#   #pheatmap creates two pages, but vdiffr supports only one
#   #viz checked heat; vdiffr crashes
#   ezh <- function() ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA)
#   # vdiffr::expect_doppelganger(title="heat", fig=ezh)
#
#   md5 <- tools::md5sum(files="../figs/ezheat/heat.pdf")
#   #as.char strips name from md5
#   expect_equal(as.character(md5), "acb389e70add03fb6f7629373d6716fd")
# })
