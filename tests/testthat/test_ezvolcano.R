context("ezvolcano")

test_that("vdiffr", {
  ezvol <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol')
  vdiffr::expect_doppelganger(title="vol", fig=ezvol)

  #no lab column, so no labs
  ezvol2 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                                 type.sig = "FDR")
  vdiffr::expect_doppelganger(title="vol2", fig=ezvol2)

  ezvol3 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                 cut.sig=0.01, cut.color = "green", up.color = "red", x.bound = 3, y.bound=30,
                                 lab.col='Gene.Symbol')
  vdiffr::expect_doppelganger(title="vol3", fig=ezvol3)
})
