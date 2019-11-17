context("barplot of significance per pathway")

test_that("direction properly specified", {
  expect_error(barplot_pwys(rc))
  expect_error(barplot_pwys(rc, prefix.v = "First3"))
  expect_error(barplot_pwys(tab=rc, direction = "Up"))

  expect_silent(barplot_pwys(tab=rc, prefix.v = "First3", direction = "Up"))
  expect_silent(barplot_pwys(tab=rc, prefix.v = "Last3", direction = "Up"))
  expect_silent(barplot_pwys(tab=rc, prefix.v = "Last3vsFirst3", direction = "Down"))
})

test_that("barplot object non-vdiffr", {
  bar_plot <- barplot_pwys(tab=rc, prefix.v = "First3", direction = "Up")
  expect_equal(bar_plot$data["pwy1", "Direction"] , rc["pwy1", "First3.Direction"])
  expect_equal(bar_plot$data["pwy1", "neglog10p"] , log10(rc["pwy1", "First3.p"])*-1)
  expect_equal(bar_plot$layers[[1]]$aes_params$fill, "red")
  expect_equal(bar_plot$layers[[1]]$geom_params$width , 0.7)
  expect_equal(bar_plot$layers[[2]]$aes_params$hjust , -0.1)

  bar_plot_blue <- barplot_pwys(tab=rc, prefix.v = "Last3vsFirst3", direction = "Down")
  expect_equal(bar_plot_blue$layers[[1]]$aes_params$fill ,"blue")
  expect_equal(bar_plot_blue$layers[[2]]$aes_params$hjust , 1.1)
})

