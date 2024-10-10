context("multi_cor_heat")

test_that("non vdiffr multiple heatmap", {
  pheno2.num <- apply(pheno2, 2, FUN=\(xx) as.numeric(as.factor(xx))) |>
    data.frame()
  pheno2.num[1, 1] <- NA
  rownames(pheno2.num) <- rownames(pheno2)
  mc.res <- ezlimma::multi_cor(M, pheno2.num)
  expect_silent(mch <- multi_cor_heat(tab=mc.res, object = M, pheno.tab = pheno2.num, p.thresh = 0.1, plot=TRUE))
  expect_true("gene1" %in% rownames(mch$grp$mat))
})
