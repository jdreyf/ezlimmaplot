context("multi_covar_pca")

test_that("multi_covar_pca agrees with ezpca", {
  ezp <- ezpca(M, pheno, name=NA, shape="grp", color="covar_num", manual.shape = 1:2, plot=FALSE)$data
  mcp <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num"),
                         manual.shape = 1:2, plot=FALSE)
  expect_equal(ezp, mcp[["covar_num"]]$data)
})

test_that("non-vdiffr test", {

  #basic sanity test
  mcp <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num","sample"),
                         plot=FALSE)
  expect_equal(mcp$tissue$labels$shape, "grp")
  expect_equal(mcp$covar_num$labels$shape, "grp")
  expect_equal(mcp$sample$labels$shape, "grp")
  expect_equal(mcp$tissue$labels$colour, "tissue")
  expect_equal(mcp$covar_num$labels$colour, "covar_num")
  expect_equal(mcp$sample$labels$colour, "sample")

  #test for manual color
  mcp1 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num"),
                         manual.color = c("red","cyan") , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp1$tissue)$data[[1]]$colour[1], "cyan")
  expect_equal(ggplot2::ggplot_build(mcp1$tissue)$data[[1]]$colour[2], "red")
  #excpect error , manual color provided for continous variable
  expect_error(plot(mcp1$covar_num))

  #test for manual shape
  mcp2 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num"),
                         manual.shape = 1:2 , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp2$tissue)$data[[1]]$shape[1], 1)
  expect_equal(ggplot2::ggplot_build(mcp2$tissue)$data[[1]]$shape[6], 2)
  expect_equal(ggplot2::ggplot_build(mcp2$covar_num)$data[[1]]$shape[1], 1)
  expect_equal(ggplot2::ggplot_build(mcp2$covar_num)$data[[1]]$shape[6], 2)

  #test for lebels being displayed as row.names of pheno
  mcp3 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num"),
                         labels = TRUE , plot=FALSE)
  expect_equal(mcp3$tissue$data$row_names, row.names(pheno))
  expect_equal(mcp3$covar_num$data$row_names, row.names(pheno))


  #test for labels being displayed as row.names of pheno
  mcp4 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="grp", covars = c("tissue", "covar_num"),
                         rm.leg.title = TRUE , plot=FALSE)
  expect_null(mcp4$tissue$theme$legend.title[1][[1]])
  expect_null(mcp4$tissue$theme$legend.title[2][[1]])
  expect_null(mcp4$covar_num$theme$legend.title[1][[1]])
  expect_null(mcp4$covar_num$theme$legend.title[2][[1]])

  #test for facet and labels in facet
  mcp5 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="sample", covars = c("tissue", "covar_num"),
                         facet = "grp" , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp5$tissue)$layout$layout$grp[1], "First3")
  expect_equal(ggplot2::ggplot_build(mcp5$tissue)$layout$layout$grp[2], "Last3")
  expect_equal(ggplot2::ggplot_build(mcp5$covar_num)$layout$layout$grp[1], "First3")
  expect_equal(ggplot2::ggplot_build(mcp5$covar_num)$layout$layout$grp[2], "Last3")
  expect_equal(mcp5$tissue$labels$shape, "sample")
  expect_equal(mcp$tissue$labels$colour, "tissue")

  #test for all.size variable size
  mcp6 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="sample", covars = c("tissue", "covar_num"),
                         all.size = 1:6 , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp6$tissue)$data[[1]]$size, 1:6)
  expect_equal(ggplot2::ggplot_build(mcp6$covar_num)$data[[1]]$size, 1:6)

  #test for all.size constant size
  mcp7 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="sample", covars = c("tissue", "covar_num"),
                         all.size = 8 , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp7$tissue)$data[[1]]$size[1], 8)
  expect_equal(ggplot2::ggplot_build(mcp7$covar_num)$data[[1]]$size[3], 8)

  #test for all.size alpha
  mcp8 <- multi_covar_pca(M, pheno, name="covar_pca", grp.var="sample", covars = c("tissue", "covar_num"),
                         alpha = 0.8 , plot=FALSE)
  expect_equal(ggplot2::ggplot_build(mcp8$tissue)$data[[1]]$alpha[1], 0.8)
  expect_equal(ggplot2::ggplot_build(mcp8$covar_num)$data[[1]]$alpha[5], 0.8)
})
