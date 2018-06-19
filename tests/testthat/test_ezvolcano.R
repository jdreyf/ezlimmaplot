context("ezvolcano")

test_that("vdiffr", {
  ezvol <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"))
  vdiffr::expect_doppelganger(title="vol", fig=ezvol)

  #no lab column, so no labs
  ezvol2 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                                 type.sig = "FDR")
  vdiffr::expect_doppelganger(title="vol2", fig=ezvol2)

  ezvol3 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                        cut.sig=0.01, cut.color = "black", up.ann.color = "red", x.bound = 3, y.bound=30)
  vdiffr::expect_doppelganger(title="vol3", fig=ezvol3)

  #check that no error with no cut or ann points
  ezvol4 <- function() ezvolcano(tab=res.df[-1,], comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
            cut.sig=0.01, cut.color = "blue", up.ann.color = "red", x.bound = 3, y.bound=30)
  vdiffr::expect_doppelganger(title="vol4", fig=ezvol4)

  #lab.col is NULL, but give pts to annotate
  ezvol5 <- function() suppressWarnings(ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                        lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
  vdiffr::expect_doppelganger(title="vol5", fig=ezvol5)

  expect_warning(ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                           lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
})
