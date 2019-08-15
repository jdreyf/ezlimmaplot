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
  ezvol5 <- function() suppressWarnings(ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1,
                        cut.lfc=1, lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
  vdiffr::expect_doppelganger(title="vol5", fig=ezvol5)

  #shape=1 -> empty dots
  ezvol6 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), shape = 1)
  vdiffr::expect_doppelganger(title="vol6", fig=ezvol6)

  expect_warning(ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                           lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), plot=FALSE))

  ezvol7 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), p05.line = TRUE)
  vdiffr::expect_doppelganger(title="vol7", fig=ezvol7)

  expect_warning(ezvolcano(tab=res.df, comparison = "First3", name=NA, type.sig="FDR", p05.line = TRUE, plot=FALSE))
})

test_that("non-vdiffr", {
  ezvol <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                     cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Plot layers match expectations
  expect_null(ezvol$layers[[1]]$geom$objname)
  expect_null(ezvol$layers[[1]]$stat$objname)
  expect_equal(ezvol$labels$label, "Gene.Symbol")
  #Scale range is NULL
  expect_null(ezvol$scales$scales[[1]]$range$range)

  ezvol2 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                      type.sig = "FDR", plot = FALSE)
  expect_equal(ezvol2$data["gene12","First3.logFC"], res.df["gene12","First3.logFC"])
  #Validating the location of geompoint
  expect_is(ezvol2$layers[[1]], "gg")
  expect_equal(ggplot2::ggplot_build(ezvol2)$data[[1]]$colour, "green")
  expect_null(ezvol2$layers[[1]]$geom$objname)
  expect_null(ezvol2$layers[[1]]$stat$objname)
  expect_null(ezvol2$labels$label)
  #Scale range is NULL
  expect_null(ezvol2$scales$scales[[1]]$range$range)

  ezvol3 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "black", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  expect_equal(ezvol3$data["gene25","First3.avg"], res.df["gene25","First3.avg"])

  expect_null(ezvol3$layers[[1]]$geom$objname)
  expect_null(ezvol3$layers[[1]]$stat$objname)
  expect_equal(ezvol3$labels$label, "Gene.Symbol")
  #Scale range is NULL
  expect_null(ezvol3$scales$scales[[1]]$range$range)

  ezvol4 <- ezvolcano(tab=res.df[-1,], comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "blue", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  expect_equal(ezvol4$data["gene25","First3.avg"], res.df["gene25","First3.avg"])

  expect_null(ezvol4$layers[[1]]$geom$objname)
  expect_null(ezvol4$layers[[1]]$stat$objname)
  expect_equal(ezvol4$labels$label, "Gene.Symbol")
  #Scale range is NULL
  expect_null(ezvol4$scales$scales[[1]]$range$range)


  ezvol5 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, plot = FALSE,
                      cut.lfc=1, lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"))

  expect_equal(ezvol5$data["gene25","First3.avg"], res.df["gene25","First3.avg"])

  expect_null(ezvol5$layers[[1]]$geom$objname)
  expect_null(ezvol5$layers[[1]]$stat$objname)
  expect_null(ezvol5$labels$label)
  #Scale range is NULL
  expect_null(ezvol5$scales$scales[[1]]$range$range)

  ezvol6 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), shape = 1)

  expect_equal(ezvol6$data["gene25","First3.avg"], res.df["gene25","First3.avg"])

  expect_null(ezvol6$layers[[1]]$geom$objname)
  expect_null(ezvol6$layers[[1]]$stat$objname)
  expect_equal(ezvol6$labels$label, "Gene.Symbol")
  #Scale range is NULL
  expect_null(ezvol6$scales$scales[[1]]$range$range)

  ezvol7 <- ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), p05.line = TRUE)
  expect_equal(ezvol7$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_null(ezvol7$layers[[1]]$geom$objname)
  expect_null(ezvol7$layers[[1]]$stat$objname)
  expect_equal(ezvol7$labels$label, "Gene.Symbol")
  #Scale range is NULL
  expect_null(ezvol7$scales$scales[[1]]$range$range)
})
