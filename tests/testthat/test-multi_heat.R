context("multi_heat")

test_that("mat", {
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  # warning prune_mat n rows < 50
  expect_warning(mat0 <- multi_heat(tab=res.ss, object=M[topgenes,], labrows = lr, pheno.df=pheno.df,
                                    name=NA, clip=1, plot=FALSE)[[1]]$mat)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- multi_heat(tab=res.ss, M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"],
                                    pheno.df=pheno.df, name=NA, clip=1, only.labrows = TRUE, plot=FALSE,
                                    verbose=TRUE, ntop=10)[[1]]$mat)
})

# file empty :-(
# test_that("vdiffr", {
#   mh1 <- function(){
#     multi_heat(tab=res.ss, M, labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df,
#                       name="tmp", ntop=20)[[1]]$gtable
#   }
#   vdiffr::expect_doppelganger(title="mheat", fig=mh1)
# })
