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
  #Validating geom and stat arguments for the layers
  expect_null(ezvol$layers[[1]]$geom$objname)
  expect_null(ezvol$layers[[1]]$stat$objname)
  #Validating the label
  expect_equal(ezvol$labels$label, "Gene.Symbol")
  #Validating the value of the outlier
  expect_equal(ezvol$data["gene1","nlg10sig"],boxplot.stats(ezvol$data$nlg10sig)$out[1])
  #Validating the limits of the scale of the ggplot
  expect_lte(ezvol$scales$scales[[1]]$limits[2], 3)
  expect_gte(ezvol$scales$scales[[1]]$limits[1], -3)
  #Validating the color of the point
  expect_equal(ggplot2::ggplot_build(ezvol)$data[[1]]$colour[[1]], "black")

  expect_equal(ggplot2::ggplot_build(ezvol)$data[[2]]$y[1], -log10(ezvol$data["gene1","First3.p"]))

  ezvol2 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                      type.sig = "FDR", plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol2$data["gene12","First3.logFC"], res.df["gene12","First3.logFC"])
  #Validating the location of geompoint
  expect_is(ezvol2$layers[[1]], "gg")
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol2)$data[[1]]$colour, "green")
  expect_equal(ggplot2::ggplot_build(ezvol2)$data[[1]]$y, -log10(ezvol2$data["gene1","First3.FDR"]))

  ezvol3 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "black", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol3$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_equal(ggplot2::ggplot_build(ezvol3)$data[[1]]$y, -log10(ezvol3$data["gene1","First3.p"]))
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol3)$data[[1]]$colour, "red")
  ezvol4 <- ezvolcano(tab=res.df[-1,], comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "blue", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol4$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating point D
  #Validating geom and stat arguments for the layers
  expect_equal(ggplot2::ggplot_build(ezvol4)$data[[1]]$y, -log10(ezvol4$data["gene59","First3.p"]))
  expect_null(ezvol4$layers[[1]]$geom$objname)
  expect_null(ezvol4$layers[[1]]$stat$objname)
  #Validating the label
  expect_equal(ezvol4$labels$label, "Gene.Symbol")
  #Validating the limits of the scale of the ggplot
  expect_lte(ezvol4$scales$scales[[1]]$limits[2], 3)
  expect_gte(ezvol4$scales$scales[[1]]$limits[1], -3)
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol4)$data[[1]]$colour, "black")

  expect_warning(ezvol5 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, plot = FALSE,
                                     cut.lfc=1, lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
  expect_equal(ezvol5$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol5)$data[[1]]$colour, "green")
  # no label
  expect_null(ggplot2::ggplot_build(ezvol5)$data[[1]]$label)
  expect_equal(ggplot2::ggplot_build(ezvol5)$data[[1]]$y, -log10(ezvol5$data["gene1","First3.p"]))

  ezvol6 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), shape = 1)
  #Validating the data point
  expect_equal(ezvol6$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating the label column
  expect_equal(ezvol6$labels$label, "Gene.Symbol")
  #Validating the shape of the plot
  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[1]]$shape[1], 19)
  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[3]]$shape[1], 1)

  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[1]]$y[[1]], -log10(ezvol6$data["gene1","First3.p"]))

  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[1]]$y[[2]], -log10(ezvol6$data["gene25","First3.p"]))

  ezvol7 <- ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), p05.line = TRUE)
  #Validating point A
  expect_gte(ggplot2::ggplot_build(ezvol7)$data[[2]]$y[1],29)
  #Validating the intercept point crossed by the straight line
  expect_equal(ggplot2::ggplot_build(ezvol7)$data[[4]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(ezvol7)$data[[1]]$y[[2]], -log10(ezvol7$data["gene25","First3.p"]))
  expect_equal(ggplot2::ggplot_build(ezvol7)$data[[1]]$y[[1]], -log10(ezvol7$data["gene1","First3.p"]))
})
