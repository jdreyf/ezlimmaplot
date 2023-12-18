context("multi_volcano")

test_that("missing logFC or p-value/FDR columns", {
 set.seed(0)
 cols <-  grep(paste0('\\', ".", 'logFC$'), colnames(res.df))
 res.df2 <- res.df[,-cols]
 expect_error(multi_volcano(tab=res.df2, name="tmp", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, type.sig="FDR", same.scale = TRUE,
                            cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol'), plot=FALSE)
 res.df3 <- res.df[,-c(grep(paste0('\\', ".", 'FDR'), colnames(res.df)))]
 expect_error(multi_volcano(tab=res.df3, name="tmp", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, type.sig="FDR", same.scale = TRUE,
                            cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol'), plot=FALSE)
})

test_that("non vdiffr",{
  # error that ntop.sig > 0 but lab.col is null.
  expect_error(multi_volcano(tab=res.df, lab.col=NULL, name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                     cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), plot = FALSE))
  mvl <- multi_volcano(tab=res.df, name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=0,
                       cut.sig=6.5*10**(-3), cut.color = "green", ann.rnames=c("gene1", "gene25"), plot = FALSE)
  #Validating the location of the specific point
  expect_equal(mvl$First3$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_equal(mvl$Last3$data["gene25","Last3.avg"], res.df["gene25","Last3.avg"])
  expect_equal(mvl$Last3vsFirst3$data["gene25","Last3vsFirst3.logFC"], res.df["gene25","Last3vsFirst3.logFC"])
  #Validating the title
  expect_equal(mvl$First3$labels$title, "First3")
  expect_equal(mvl$Last3$labels$title, "Last3")
  expect_equal(mvl$Last3vsFirst3$labels$title, "Last3vsFirst3")
  #Validating the value of the outlier
  expect_equal(mvl$First3$data["gene1","nlg10sig"], -log10(res.df["gene1","First3.p"]))
  expect_equal(mvl$Last3$data["gene1","nlg10sig"], -log10(res.df["gene1","Last3.p"]))
  expect_gte(mvl$Last3vsFirst3$data["gene1","nlg10sig"], -log10(res.df["gene1","Last3vsFirst3.p"]))
  #Validating the limits of the scale of the ggplot
  expect_equal(mvl$First3$scales$scales[[1]]$limits[2], max(res.df$First3.avg))
  expect_equal(mvl$First3$scales$scales[[1]]$limits[1], -max(res.df$First3.avg))
  expect_equal(mvl$Last3$scales$scales[[1]]$limits[2], -min(res.df$Last3.avg))
  expect_equal(mvl$Last3$scales$scales[[1]]$limits[1], min(res.df$Last3.avg))
  expect_equal(mvl$Last3vsFirst3$scales$scales[[1]]$limits[2], -min(res.df$Last3vsFirst3.logFC))
  expect_equal(mvl$Last3vsFirst3$scales$scales[[1]]$limits[1], min(res.df$Last3vsFirst3.logFC))
  #Validating the color of the point
  expect_equal(ggplot2::ggplot_build(mvl$First3)$data[[4]]$colour[4], "green")
  expect_equal(ggplot2::ggplot_build(mvl$First3)$data[[4]]$y[1], -log10(mvl$First3$data["gene1", "First3.p"]))
  expect_equal(ggplot2::ggplot_build(mvl$Last3)$data[[4]]$colour[[1]], "black")
  expect_equal(ggplot2::ggplot_build(mvl$Last3)$data[[4]]$y[1], -log10(mvl$Last3$data["gene1", "Last3.p"]))
  expect_equal(ggplot2::ggplot_build(mvl$Last3vsFirst3)$data[[4]]$colour[3], "green")
  expect_equal(ggplot2::ggplot_build(mvl$Last3vsFirst3)$data[[4]]$y[1], -log10(mvl$Last3vsFirst3$data["gene1", "Last3vsFirst3.p"]))

  mvl2 <-  multi_volcano(tab=res.df, name=NA, cut.lfc=1, cut.sig=0.01, cut.color = "green", up.ann.color = "green", down.ann.color = "green",
                      type.sig = "FDR", plot = FALSE)
  #Validating the location of the specific point
  expect_equal( mvl2$First3$data["gene12","First3.logFC"], res.df["gene12","First3.logFC"])
  expect_equal( mvl2$Last3$data["gene12","Last3.logFC"], res.df["gene12","Last3.logFC"])
  expect_equal( mvl2$Last3vsFirst3$data["gene12","Last3vsFirst3.logFC"], res.df["gene12","Last3vsFirst3.logFC"])
  #Validating the location of geompoint
  expect_is( mvl2$First3$layers[[1]], "gg")
  expect_is( mvl2$Last3$layers[[1]], "gg")
  expect_is( mvl2$Last3vsFirst3$layers[[1]], "gg")
  #Validating the colour of the point
  expect_equal(ggplot2::ggplot_build( mvl2$First)$data[[4]]$colour[1], "green")
  expect_equal(ggplot2::ggplot_build( mvl2$First)$data[[4]]$y[1], -log10(mvl2$First3$data["gene1","First3.FDR"]))
  expect_equal(ggplot2::ggplot_build( mvl2$Last3)$data[[4]]$colour[1], "black")
  expect_equal(ggplot2::ggplot_build( mvl2$Last3)$data[[4]]$y[1], -log10(mvl2$Last3$data["gene1","Last3.FDR"]))
  expect_equal(ggplot2::ggplot_build( mvl2$Last3vsFirst3)$data[[4]]$colour[1], "green")
  expect_equal(ggplot2::ggplot_build(mvl2$Last3vsFirst3)$data[[4]]$y[1],
                       -log10(mvl2$Last3vsFirst3$data["gene1","Last3vsFirst3.FDR"]))

  mvlc3 <- multi_volcano(tab=res.df, name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1,
                      cut.sig=0.01, cut.color = "black", up.ann.color = "red", plot = FALSE)
  #Validating the location of the specific point
  expect_equal(mvlc3$First$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_equal(ggplot2::ggplot_build(mvlc3$First)$data[[4]]$y[1], -log10(mvlc3$First$data["gene1","First3.p"]))
  expect_equal(mvlc3$Last3$data["gene25","Last3.avg"], res.df["gene25","Last3.avg"])
  expect_equal(ggplot2::ggplot_build(mvlc3$Last3)$data[[4]]$y[1], -log10(mvlc3$Last3$data["gene1","Last3.p"]))
  expect_equal(mvlc3$Last3vsFirst3$data["gene25","Last3vsFirst3.avg"], res.df["gene25","Last3vsFirst3.avg"])
  expect_equal(ggplot2::ggplot_build(mvlc3$Last3vsFirst3)$data[[4]]$y[1],
                       -log10(mvlc3$Last3vsFirst3$data["gene1","Last3vsFirst3.p"]))
  #Validating the colour of the point as black as col.lab is NULL, so it doesnt colors data
  expect_equal(ggplot2::ggplot_build(mvlc3$First)$data[[4]]$colour[1], "red")
  expect_equal(mvlc3$Last3$data["gene1", "color.point"], "black")
  expect_equal(mvlc3$Last3vsFirst3$data["gene1", "color.point"], "black")

  mvlc4 <- multi_volcano(tab=res.df[-1,], name=NA, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, lab.col = "Gene.Symbol",
                      cut.sig=0.01, cut.color = "blue", up.ann.color = "red", plot = FALSE)
  #Validating the location of the specific point
  expect_equal(mvlc4$First3$data["gene25","First3.avg"], res.df["gene25","First3.avg"])
  expect_equal(mvlc4$Last3$data["gene25","Last3.avg"], res.df["gene25","Last3.avg"])
  expect_equal(mvlc4$Last3vsFirst3$data["gene25","Last3vsFirst3.avg"], res.df["gene25","Last3vsFirst3.avg"])
  #Validating point D
  #Validating geom and stat arguments for the layers
  expect_true(any(ggplot2::ggplot_build(mvlc4$First3)$data[[4]]$y == -log10(res.df["gene59", "First3.p"])))
  #Validating point C
  #Validating geom and stat arguments for the layers
  expect_equal(ggplot2::ggplot_build(mvlc4$Last3)$data[[5]]$label[2], res.df["gene74", "Gene.Symbol"])
  #Validating point D
  #Validating geom and stat arguments for the layers
  expect_equal(ggplot2::ggplot_build(mvlc4$Last3vsFirst3)$data[[4]]$y[3],
                      -log10(mvlc4$Last3vsFirst3$data["gene59", "Last3vsFirst3.p"]))
   #Validating the label
  expect_equal(mvlc4$First3$labels$title, "First3")
  expect_equal(mvlc4$Last3$labels$title, "Last3")
  expect_equal(mvlc4$Last3vsFirst3$labels$title, "Last3vsFirst3")
  #Validating the colour of the point
  # expect_equal(ggplot2::ggplot_build(mvlc4$First3)$data[[4]]$colour[1], "black")
  # expect_equal(ggplot2::ggplot_build(mvlc4$Last3)$data[[4]]$colour, "black")
  # expect_equal(ggplot2::ggplot_build(mvlc4$Last3vsFirst3)$data[[4]]$colour, "red")

  mvlc5 <- multi_volcano(tab=res.df, name=NA, ntop.sig = 1, ntop.lfc = 1, plot = FALSE, cut.lfc=1,
                                     cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"))
  #Validating the colour of the point
  # expect_equal(ggplot2::ggplot_build(mvlc5$First3)$data[[3]]$colour, "green")
  # expect_equal(ggplot2::ggplot_build(mvlc5$Last3)$data[[3]]$colour[1], "black")
  # expect_equal(ggplot2::ggplot_build(mvlc5$Last3vsFirst3)$data[[3]]$colour, "green")

  expect_warning(mvlc6 <- multi_volcano(tab=res.df, name="Rplots", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, lab.col = "Gene.Symbol",
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25")))
  expect_equal(ggplot2::ggplot_build(mvlc6$First3)$data[[5]]$label[1], res.df["gene1","Gene.Symbol"])
  expect_equal(ggplot2::ggplot_build(mvlc6$Last3)$data[[5]]$label[2], res.df["gene25","Gene.Symbol"])
  expect_equal(ggplot2::ggplot_build(mvlc6$Last3)$data[[5]]$label[3], res.df["gene74","Gene.Symbol"])
  expect_equal(ggplot2::ggplot_build(mvlc6$Last3vsFirst3)$data[[5]]$label[1], res.df["gene1","Gene.Symbol"])

  mvlc7 <- multi_volcano(tab=res.df, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE, lab.col = "Gene.Symbol",
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 0.05)
  #Validating points
  expect_gte(ggplot2::ggplot_build(mvlc7$First3)$data[[4]]$y[1], 29)
  expect_equal(ggplot2::ggplot_build(mvlc7$Last3)$data[[5]]$label[2], res.df["gene25","Gene.Symbol"])
  expect_equal(ggplot2::ggplot_build(mvlc7$Last3)$data[[5]]$label[3], res.df["gene74","Gene.Symbol"])
  expect_equal(ggplot2::ggplot_build(mvlc7$Last3vsFirst3)$data[[5]]$label[1], res.df["gene1","Gene.Symbol"])

  #Validating the intercept point crossed by the straight line
  expect_equal(ggplot2::ggplot_build(mvlc7$First3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc7$First3)$data[[4]]$y[1], -log10(mvlc7$First3$data["gene1", "First3.p"]))
  expect_equal(ggplot_build(mvlc7$First3)$data[[1]]$linetype, "dashed")

  expect_equal(ggplot2::ggplot_build(mvlc7$Last3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc7$Last3)$data[[4]]$y[1], -log10(mvlc7$First3$data["gene1", "Last3.p"]))
  expect_equal(ggplot_build(mvlc7$Last3)$data[[1]]$linetype, "dashed")

  expect_equal(ggplot2::ggplot_build(mvlc7$Last3vsFirst3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc7$Last3vsFirst3)$data[[4]]$y[1],
                        -log10(mvlc7$Last3vsFirst3$data["gene1", "Last3vsFirst3.p"]))
  expect_equal(ggplot_build(mvlc7$Last3vsFirst3)$data[[1]]$linetype, "dashed")

  #Validating the line label and ylab
  expect_equal(mvlc7$First3$guides$linetype$title, "p")
  expect_equal(mvlc7$Last3$guides$linetype$title, "p")
  expect_equal(mvlc7$Last3vsFirst3$guides$linetype$title, "p")
  expect_equal(mvlc7$First3$labels$y, expression("-" * log[10] ~ p * "-" * value))
  expect_equal(mvlc7$Last3$labels$y, expression("-" * log[10] ~ p * "-" * value))
  expect_equal(mvlc7$Last3vsFirst3$labels$y, expression("-" * log[10] ~ p * "-" * value))

  mvlc8 <- multi_volcano(tab=res.df, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE, lab.col = "Gene.Symbol",
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = c(0.05, 0.001))
  expect_equal(ggplot_build(mvlc8$First3)$data[[1]]$linetype[1], "dashed")
  expect_equal(ggplot_build(mvlc8$First3)$data[[1]]$linetype[2], "dotted")
  expect_equal(ggplot2::ggplot_build(mvlc8$First3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc8$First3)$data[[1]]$yintercept[2], -log10(0.001))

  expect_equal(ggplot_build(mvlc8$Last3)$data[[1]]$linetype[1], "dashed")
  expect_equal(ggplot_build(mvlc8$Last3)$data[[1]]$linetype[2], "dotted")
  expect_equal(ggplot2::ggplot_build(mvlc8$Last3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc8$Last3)$data[[1]]$yintercept[2], -log10(0.001))

  expect_equal(ggplot_build(mvlc8$Last3vsFirst3)$data[[1]]$linetype[1], "dashed")
  expect_equal(ggplot_build(mvlc8$Last3vsFirst3)$data[[1]]$linetype[2], "dotted")
  expect_equal(ggplot2::ggplot_build(mvlc8$Last3vsFirst3)$data[[1]]$yintercept[1], -log10(0.05))
  expect_equal(ggplot2::ggplot_build(mvlc8$Last3vsFirst3)$data[[1]]$yintercept[2], -log10(0.001))

  mvlc9 <- multi_volcano(tab=res.df, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE,
                      cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 0.05,
                      lab.col = "Gene.Symbol", up.ann.color = "blue")
  mvlc9.ann.mat.First3 <- ggplot2::ggplot_build(mvlc9$First3)$data[[5]]
  expect_equal(mvlc9.ann.mat.First3[mvlc9.ann.mat.First3$label == "A", "colour"], "blue")
  expect_equal(mvlc9.ann.mat.First3[mvlc9.ann.mat.First3$label == "B", "colour"], "blue")

  mvlc9.ann.mat.Last3 <- ggplot2::ggplot_build(mvlc9$Last3)$data[[5]]
  expect_equal(mvlc9.ann.mat.Last3[mvlc9.ann.mat.Last3$label == "A", "colour"], "blue")
  expect_equal(mvlc9.ann.mat.Last3[mvlc9.ann.mat.Last3$label == "B", "colour"], "blue")

  ev9.ann.mat.Last3vsFirst3 <- ggplot2::ggplot_build(mvlc9$Last3vsFirst3)$data[[5]]
  expect_equal(ev9.ann.mat.Last3vsFirst3[ev9.ann.mat.Last3vsFirst3$label == "B", "colour"], "blue")

  mvlc10 <- multi_volcano(tab=res.df, ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, name=NA, plot=FALSE, lab.col = "Gene.Symbol",
                      type.sig="FDR", cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lines.sig = 0.05)
  expect_equal(mvlc10$First3$guides$linetype$title, "FDR")
  expect_equal(mvlc10$Last3$guides$linetype$title, "FDR")
  expect_equal(mvlc10$Last3vsFirst3$guides$linetype$title, "FDR")
  expect_equal(mvlc10$First3$labels$y, expression("-" * log[10] ~ FDR))
  expect_equal(mvlc10$Last3$labels$y, expression("-" * log[10] ~ FDR))
  expect_equal(mvlc10$Last3vsFirst3$labels$y, expression("-" * log[10] ~ FDR))

  mvlc11 <- multi_volcano(tab=res.df, same.scale=TRUE, plot=FALSE)
  expect_equal(mvlc11$First3$plot_env$x.bound, mvlc11$Last3$plot_env$x.bound)
  expect_equal(mvlc11$First3$plot_env$x.bound, mvlc11$Last3vsFirst3$plot_env$x.bound)
  expect_equal(mvlc11$First3$plot_env$y.bound, mvlc11$Last3$plot_env$y.bound)
  expect_equal(mvlc11$First3$plot_env$y.bound, mvlc11$Last3vsFirst3$plot_env$y.bound)

  mvlc12 <- multi_volcano(tab=res.df, plot=FALSE, alpha=1)
  expect_equal(min(ggplot2::ggplot_build(mvlc12$First3)$data[[3]]$alpha), 1)
  expect_equal(max(ggplot2::ggplot_build(mvlc12$First3)$data[[3]]$alpha), 1)
  expect_equal(min(ggplot2::ggplot_build(mvlc12$Last3)$data[[3]]$alpha), 1)
  expect_equal(max(ggplot2::ggplot_build(mvlc12$Last3)$data[[3]]$alpha), 1)
  expect_equal(min(ggplot2::ggplot_build(mvlc12$Last3vsFirst3)$data[[3]]$alpha), 1)
  expect_equal(max(ggplot2::ggplot_build(mvlc12$Last3vsFirst3)$data[[3]]$alpha), 1)

  mvlc13 <- multi_volcano(tab=res.df, plot=FALSE, type.sig="FDR")
  expect_equal(mvlc13$First3$labels$y, expression("-" * log[10] ~ FDR))
  expect_equal(mvlc13$Last3$labels$y, expression("-" * log[10] ~ FDR))
  expect_equal(mvlc13$Last3vsFirst3$labels$y, expression("-" * log[10] ~ FDR))
})

test_that("raster", {
  expect_warning(mvol.norast <- multi_volcano(tab=res.df, name="tmp_norast", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                        ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol'))
  expect_warning(mvol.rast <- multi_volcano(tab=res.df, name="tmp_rast", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                             ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol', raster = TRUE))
  expect_gt(file.size("tmp_norast.pdf"), file.size("tmp_rast.pdf"))
  unlink("tmp_rast.pdf")
})

test_that("vdiffr", {
  expect_warning(mvol <- multi_volcano(tab=res.df, name="tmp", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, cut.sig=0.01, cut.color = "green",
                        ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol'))
  vdiffr::expect_doppelganger(title="mvol1", fig=mvol[[1]])
  vdiffr::expect_doppelganger(title="mvol2", fig=mvol[[2]])

  mvol2 <- multi_volcano(tab=res.df, name="tmp", ntop.sig = 1, ntop.lfc = 1, cut.lfc=1, type.sig="FDR", same.scale = TRUE,
                         cut.sig=0.01, cut.color = "green", ann.rnames=c("gene1", "gene25"), lab.col='Gene.Symbol')
  vdiffr::expect_doppelganger(title="mvol3", fig=mvol2[[1]])
  vdiffr::expect_doppelganger(title="mvol4", fig=mvol2[[2]])
})
