context("ezvolcano")

test_that("non-vdiffr", {
  set.seed(0)
  ezvol <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                     cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating the label
  expect_equal(ezvol$labels$title, "First3")
  #Validating the value of the outlier
  expect_equal(ezvol$data["gene1","nlg10sig"], -log10(res.df["gene1","First3.p"]))
  #Validating the limits of the scale of the ggplot
  expect_equal(ezvol$scales$scales[[1]]$limits[2], max(res.df$First3.avg))
  expect_equal(ezvol$scales$scales[[1]]$limits[1], -max(res.df$First3.avg))
  #Validating the color of the point
  expect_equal(ggplot2::ggplot_build(ezvol)$data[[3]]$colour[[1]], "black")
  expect_equal(ggplot2::ggplot_build(ezvol)$data[[3]]$y[1], -log10(ezvol$data["gene1", "First3.p"]))

  ezvol2 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                      type.sig = "FDR", plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol2$data["gene12","First3.logFC"], res.df["gene12","First3.logFC"])
  #Validating the location of geompoint
  expect_is(ezvol2$layers[[1]], "gg")
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol2)$data[[3]]$colour, "green")
  expect_equal(ggplot2::ggplot_build(ezvol2)$data[[3]]$y[1], -log10(ezvol2$data["gene1","First3.FDR"]))

  ezvol3 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "black", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol3$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_equal(ggplot2::ggplot_build(ezvol3)$data[[3]]$y, -log10(ezvol3$data["gene1","First3.p"]))
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol3)$data[[3]]$colour, "red")

  ezvol4 <- ezvolcano(tab=res.df[-1,], comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "blue", up.ann.color = "red", x.bound = 3, y.bound=30, plot = FALSE)
  #Validating the location of the specific point
  expect_equal(ezvol4$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating point D
  #Validating geom and stat arguments for the layers
  expect_equal(ggplot2::ggplot_build(ezvol4)$data[[3]]$y,  -log10(res.df["gene59", "First3.p"]))

  #Validating the label
  expect_equal(ezvol4$labels$title, "First3")
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol4)$data[[3]]$colour, "black")

  expect_warning(ezvol5 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, plot = FALSE,
                                     cut.lfc=1, lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
  expect_equal(ezvol5$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build(ezvol5)$data[[3]]$colour, "green")
  # no label
  expect_null(ggplot2::ggplot_build(ezvol5)$data[[3]]$label)
  expect_equal(ggplot2::ggplot_build(ezvol5)$data[[3]]$y, -log10(ezvol5$data["gene1","First3.p"]))

  ezvol6 <- ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), shape = 1)
  #Validating the data point
  expect_equal(ezvol6$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  #Validating the label column
  expect_equal(ezvol6$labels$title, "First3")
  #Validating the shape of the plot
  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[3]]$shape[1], 19)
  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[5]]$shape[1], 1)

  expect_equal(ggplot2::ggplot_build(ezvol6)$data[[3]]$y[1], -log10(ezvol6$data["gene1","First3.p"]))

  # expect_true(any(ggplot2::ggplot_build(ezvol6)$data[[3]]$y[1] == -log10(ezvol6$data["gene25","First3.p"])))

  expect_warning(ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                           lab.col=NULL, cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), plot=FALSE))

  ezvol7 <- ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 0.05)
  #Validating point A
  expect_gte(ggplot2::ggplot_build(ezvol7)$data[[3]]$y[1], 29)
  #Validating the intercept point crossed by the straight line
  expect_equal(ggplot2::ggplot_build(ezvol7)$data[[6]]$yintercept[1], -log10(0.05))
  # expect_true(any(ggplot2::ggplot_build(ezvol7)$data[[3]]$y == -log10(ezvol7$data["gene25", "First3.p"])))
  expect_equal(ggplot2::ggplot_build(ezvol7)$data[[3]]$y[1], -log10(ezvol7$data["gene1", "First3.p"]))
  # expect_equal(ezvol7$layers[[4]]$aes_params$linetype, 2)
  expect_equal(ggplot_build(ezvol7)$data[[6]]$linetype, "dashed")
  # also test color of A, which should be green, but isn't

  ezvol8 <- ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = c(0.05, 0.001))
  expect_equal(ggplot_build(ezvol8)$data[[6]]$linetype[1], "dashed")
  expect_equal(ggplot_build(ezvol8)$data[[6]]$linetype[2], "dotted")
  expect_equal(ggplot2::ggplot_build(ezvol8)$data[[6]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(ezvol8)$data[[6]]$yintercept[2], -log10(0.001))

  ezvol9 <- ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 0.05,
                      lab.col = "Gene.Symbol", up.ann.color = "blue")
  ev9.ann.mat <- ggplot2::ggplot_build(ezvol9)$data[[4]]
  expect_equal(ev9.ann.mat[ev9.ann.mat$label == "A", "colour"], "blue")
  expect_equal(ev9.ann.mat[ev9.ann.mat$label == "B", "colour"], "blue")

  expect_error(ezvolcano(tab=res.df, comparison = "First3", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 10**(-1*1:6)))

  ezvol10 <- ezvolcano(tab=res.df, comparison = "First3", plot=FALSE, alpha=1)
  expect_equal(min(ggplot2::ggplot_build(ezvol10)$data[[3]]$alpha), 1)
  expect_equal(max(ggplot2::ggplot_build(ezvol10)$data[[3]]$alpha), 1)

  ezvol11 <- ezvolcano(tab=res.df, comparison = "First3", plot=FALSE, type.sig="FDR")
  expect_equal(ezvol11$labels$y, expression("-" * log[10] ~ FDR))
})

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
                                                  cut.lfc=1, lab.col=NULL, cut.sig=0.01, cut.color = "green",
                                                  ann.rnames=c("gene1", "gene25")))
  vdiffr::expect_doppelganger(title="vol5", fig=ezvol5)

  #shape=1 -> empty dots
  ezvol6 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                 cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), shape = 1)
  vdiffr::expect_doppelganger(title="vol6", fig=ezvol6)

  ezvol7 <- function() ezvolcano(tab=res.df, comparison = "First3", name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                                 cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"))
  vdiffr::expect_doppelganger(title="vol7", fig=ezvol7)
})

