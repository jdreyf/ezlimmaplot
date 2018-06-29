context("multi_volcano")

test_that("vdiffr", {
  mvol <- multi_volcano(tab=res.df, name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                        ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol')
  vdiffr::expect_doppelganger(title="mvol1", fig=mvol[[1]])
  vdiffr::expect_doppelganger(title="mvol2", fig=mvol[[2]])

  mvol2 <- multi_volcano(tab=res.df, name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, type.sig="FDR", same.scale = TRUE,
                         cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol')
  vdiffr::expect_doppelganger(title="mvol3", fig=mvol2[[1]])
  vdiffr::expect_doppelganger(title="mvol4", fig=mvol2[[2]])
})
