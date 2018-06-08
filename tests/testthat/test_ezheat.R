context("ezheat")

test_that("ezheat", {
  pheno.df <- data.frame(row.names = rownames(pheno), grp=pheno$grp)
  topgenes <- rownames(res.df)[1:20]

  mat0 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1)
  expect_gt(max(mat0[1,1:2]), 0.75)
  expect_lt(min(mat0[1,4:6]), -0.5)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_equal(rownames(mat0)[11],  topgenes[11])
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1,
                 only.symbols = TRUE))
  expect_equal(nrow(mat1), sum(!is.na(res.df$Gene.Symbol[1:20])))

  #pheatmap creates two pages, but vdiffr supports only one
  #viz checked heat; vdiffr crashes
  # ezh <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA)
  # vdiffr::expect_doppelganger("heat", ezh)
})
